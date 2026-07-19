#!/usr/bin/env bash
# Builds Vina-GPU+ (https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0) for
# GPU-accelerated docking and installs it into the active conda/mamba
# environment. Linux-only: dd_docking always falls back to the CPU `vina`
# package on macOS/Windows (see dd_docking/gpu_backend.py), so there is
# nothing to build there and this script just exits cleanly.
#
# Usage:
#   mamba activate dd_docking
#   bash scripts/build_vina_gpu.sh
set -euo pipefail

VINA_GPU_COMMIT=64bdb0927e0ffa839ff0820a3255eacbfc01c128
BOOST_VERSION=1.77.0
BOOST_VERSION_U=1_77_0

log() { echo "[build_vina_gpu] $*" >&2; }

if [[ "$(uname -s)" != "Linux" ]]; then
  log "Not on Linux ($(uname -s)) -- dd_docking uses the CPU vina backend here, nothing to build. Exiting."
  exit 0
fi

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  log "CONDA_PREFIX is not set -- activate the dd_docking environment first (mamba activate dd_docking)."
  exit 1
fi

if ! command -v clinfo >/dev/null 2>&1 && ! command -v nvidia-smi >/dev/null 2>&1; then
  log "WARNING: neither clinfo nor nvidia-smi found. Vina-GPU+ needs a working OpenCL runtime; continuing anyway, but the build or the binary may fail."
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THIRD_PARTY="$REPO_ROOT/third_party"
INSTALL_DIR="$CONDA_PREFIX/share/dd_docking/vina-gpu"
mkdir -p "$THIRD_PARTY" "$INSTALL_DIR"

# --- 1. Boost source (Vina-GPU+'s Makefile compiles a few Boost.Thread
#        translation units directly from source rather than only linking a
#        prebuilt libboost_thread, so a real Boost source tree is needed --
#        conda's boost-cpp package only ships headers + prebuilt libs).
BOOST_DIR="$THIRD_PARTY/boost_${BOOST_VERSION_U}"
if [[ ! -f "$BOOST_DIR/stage/lib/libboost_program_options.so" && ! -f "$BOOST_DIR/stage/lib/libboost_program_options.a" ]]; then
  if [[ ! -d "$BOOST_DIR" ]]; then
    log "Downloading Boost ${BOOST_VERSION} source..."
    curl -fL --progress-bar \
      "https://archives.boost.org/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION_U}.tar.gz" \
      -o "$THIRD_PARTY/boost_${BOOST_VERSION_U}.tar.gz"
    log "Extracting Boost source..."
    tar -xzf "$THIRD_PARTY/boost_${BOOST_VERSION_U}.tar.gz" -C "$THIRD_PARTY"
    rm -f "$THIRD_PARTY/boost_${BOOST_VERSION_U}.tar.gz"
  fi
  log "Bootstrapping and building boost_program_options/system/filesystem (this takes a few minutes)..."
  (
    cd "$BOOST_DIR"
    ./bootstrap.sh --with-libraries=program_options,system,filesystem,thread
    ./b2 --with-program_options --with-system --with-filesystem --with-thread -j"$(nproc)"
  )
else
  log "Boost ${BOOST_VERSION} already built in $BOOST_DIR, skipping."
fi

# --- 2. Vina-GPU-2.0 source, pinned commit.
#        The repo tracks its prebuilt Kernel*_Opt.bin/exe files via git-lfs,
#        so git-lfs must be installed (mamba install -c conda-forge git-lfs)
#        and registered before cloning.
if ! command -v git-lfs >/dev/null 2>&1; then
  log "ERROR: git-lfs not found. Install it first: mamba install -n dd_docking -c conda-forge git-lfs"
  exit 1
fi
git lfs install --skip-repo

VINA_GPU_SRC="$THIRD_PARTY/Vina-GPU-2.0"
if [[ ! -d "$VINA_GPU_SRC" ]]; then
  log "Cloning Vina-GPU-2.0 @ $VINA_GPU_COMMIT..."
  git clone https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0.git "$VINA_GPU_SRC"
  (cd "$VINA_GPU_SRC" && git checkout "$VINA_GPU_COMMIT")
else
  log "Vina-GPU-2.0 source already present in $VINA_GPU_SRC, skipping clone."
fi

BUILD_DIR="$VINA_GPU_SRC/Vina-GPU+"

# --- 3. Patch upstream OpenCL-C bugs in the kernel source. These are real
#        type errors (a redundant address-of that turns a row pointer into
#        a pointer-to-array, two pointers missing their __global address
#        space qualifier, and a call to get_global_linear_id() which is
#        OpenCL >=2.0 only) that NVIDIA's lenient OpenCL-3.0 compile path
#        happens to tolerate -- and then either segfaults while compiling
#        the kernel, or fails to reload the resulting program binary with
#        CL_BUILD_PROGRAM_FAILURE (reproduced on this project's GTX 1660
#        Ti / driver 560.35.05, whose OpenCL C compiler is natively 1.2 --
#        see `clinfo`'s "Device OpenCL C Version"). Fixing them is what
#        makes rigid-receptor GPU docking actually work here; see
#        DeltaGroupNJUPT/Vina-GPU-2.0 issues #1 and #26 for the same
#        symptom reported by others. Idempotent: skipped if already patched.
KERNEL1="$BUILD_DIR/OpenCL/src/kernels/kernel1.cl"
QUASI_NEWTON="$BUILD_DIR/OpenCL/src/kernels/quasi_newton.cpp"
KERNEL2="$BUILD_DIR/OpenCL/src/kernels/kernel2.cl"
if ! grep -q "a_coords\[3\]" "$KERNEL1"; then
  log "Patching upstream OpenCL-C kernel bugs..."
  sed -i 's/address = &(ar->relation\[temp\]);/address = ar->relation[temp];/' "$KERNEL1"
  sed -i 's/const atom_cl\* a = &pa->atoms\[i\];/const __global atom_cl* a = \&pa->atoms[i];/' "$KERNEL1"
  sed -i '/const float r2 = vec_distance_sqr(a->coords, probe_coords);/i\
		const float a_coords[3] = { a->coords[0], a->coords[1], a->coords[2] };' "$KERNEL1"
  sed -i 's/vec_distance_sqr(a->coords, probe_coords)/vec_distance_sqr(a_coords, probe_coords)/' "$KERNEL1"
  sed -i -E 's/const[[:space:]]+mis_cl\*([[:space:]]+)mis,/const __global mis_cl*\1mis,/' "$QUASI_NEWTON"
  sed -i 's/int gl = get_global_linear_id();/int gl = gy * gs + gx; \/\/ get_global_linear_id() needs OpenCL >=2.0/' "$KERNEL2"
else
  log "Kernel source already patched, skipping."
fi

# --- 4. Point the Makefile at our Boost build and the system OpenCL
#        (NVIDIA's OpenCL headers/libs live under the CUDA toolkit install;
#        AMD users should instead point OPENCL_LIB_PATH at their driver's
#        OpenCL install and set GPU_PLATFORM=-DAMD_PLATFORM below). Default
#        to OpenCL C 1.2 -- the patches above make the kernel build cleanly
#        under 1.2/2.0/3.0, and 1.2 is the most portable (it's this
#        project's dev GPU's *actual* OpenCL C compiler level per `clinfo`,
#        even though the platform itself advertises "OpenCL 3.0"). Override
#        with DD_DOCKING_OPENCL_VERSION=-DOPENCL_3_0 (or _2_0) if your
#        device's OpenCL C compiler genuinely supports it and you want it.
OPENCL_LIB_PATH="${DD_DOCKING_OPENCL_PATH:-/usr/local/cuda}"
GPU_PLATFORM="${DD_DOCKING_GPU_PLATFORM:--DNVIDIA_PLATFORM}"
OPENCL_VERSION="${DD_DOCKING_OPENCL_VERSION:--DOPENCL_1_2}"

if [[ ! -f "$OPENCL_LIB_PATH/include/CL/cl.h" ]]; then
  log "ERROR: no OpenCL headers found at $OPENCL_LIB_PATH/include/CL/cl.h. Set DD_DOCKING_OPENCL_PATH to your CUDA toolkit / OpenCL SDK root."
  exit 1
fi

sed -i \
  -e "s|^BOOST_LIB_PATH=.*|BOOST_LIB_PATH=$BOOST_DIR|" \
  -e "s|^OPENCL_LIB_PATH=.*|OPENCL_LIB_PATH=$OPENCL_LIB_PATH|" \
  -e "s|^OPENCL_VERSION=.*|OPENCL_VERSION=$OPENCL_VERSION|" \
  -e "s|^GPU_PLATFORM=.*|GPU_PLATFORM=$GPU_PLATFORM|" \
  "$BUILD_DIR/Makefile"

# --- 5. Build, in two passes per upstream's documented workflow:
#        `make source` first compiles the OpenCL kernels from source and
#        produces Kernel1_Opt.bin/Kernel2_Opt.bin; then plain `make`
#        rebuilds the binary to *load* those .bin files at startup instead
#        of recompiling the kernel every run (much faster startup).
#
#        Note: Vina-GPU+ only supports *rigid*-receptor docking -- its
#        main_procedure_cl.cpp asserts `m.num_other_pairs() == 0`, and
#        that count is nonzero whenever any flexible side chain is in
#        play (ligand-flex/flex-flex/flex-inflex pairs all count as
#        "other"). dd_docking's ensemble docking always uses flexible
#        side chains, so `--backend gpu/auto` still uses CPU QuickVina2
#        for those tasks; only a member prepared with no flexible
#        residues at all can actually run on the GPU here.
log "Building Vina-GPU+ (ulimit -s 8192, as required by upstream)..."
(
  cd "$BUILD_DIR"
  ulimit -s 8192
  make clean || true
  make source
  make clean
  make
)

BIN_NAME="Vina-GPU"
if [[ ! -f "$BUILD_DIR/$BIN_NAME" || ! -f "$BUILD_DIR/Kernel1_Opt.bin" || ! -f "$BUILD_DIR/Kernel2_Opt.bin" ]]; then
  log "ERROR: build finished but $BUILD_DIR/$BIN_NAME and/or the Kernel*_Opt.bin files were not produced."
  exit 1
fi

# --- 6. Install binary + kernel binaries together: Vina-GPU+ resolves
#        Kernel1_Opt.bin/Kernel2_Opt.bin relative to its *working directory*
#        at runtime, so they must stay next to each other, and dd_docking's
#        gpu_backend.py always runs the binary with this directory as cwd.
#        Also copy the Boost shared libs it's dynamically linked against
#        (gpu_backend.py points LD_LIBRARY_PATH at this same directory), so
#        the binary keeps working even after third_party/ is removed.
cp -f "$BUILD_DIR/$BIN_NAME" "$BUILD_DIR/Kernel1_Opt.bin" "$BUILD_DIR/Kernel2_Opt.bin" "$INSTALL_DIR/"
cp -fP "$BOOST_DIR"/stage/lib/libboost_program_options.so* \
       "$BOOST_DIR"/stage/lib/libboost_system.so* \
       "$BOOST_DIR"/stage/lib/libboost_filesystem.so* \
       "$INSTALL_DIR/"
chmod +x "$INSTALL_DIR/$BIN_NAME"

log "Installed to $INSTALL_DIR/$BIN_NAME (+ kernel binaries)."
log "Verify with: DD_DOCKING_VINA_GPU_DIR=$INSTALL_DIR python -c \"from dd_docking import gpu_backend; print(gpu_backend.gpu_binary_available())\""
