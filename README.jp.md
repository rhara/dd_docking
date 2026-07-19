# dd_docking — フレキシブルポケットのアンサンブルドッキング：MDで検証されたinduced-fit、GPUで高速に

バーチャルスクリーニング向けの、induced-fit様のポケットダイナミクスを備えた
アンサンブルドッキングツールキット。特定の標的やリガンドセットに紐付かない、
再利用可能なパッケージとして設計されている。

- **アンサンブル準備（`dd_docking-prep`）**: 複数の受容体コンフォメーション
  （PDB）をPDBFixerで構造修復し、各共結晶リガンドからドッキングボックスと
  ポケット近傍のフレキシブル残基を決定し、Meekoでrigid/flexのPDBQTペアを
  生成する。
- **アンサンブルドッキング（`dd_docking-dock`）**: `.smi` ライブラリ内の各
  分子を、フレキシブル側鎖付きでアンサンブルの全メンバーに対してドッキング
  し、メンバー間で最良（最低）のアフィニティによってランキングする
  （best-of-ensemble）。
- **MDベースの精密化（`dd_docking-refine`）**: 上位ヒットのみを短時間の
  暗黙溶媒MDシミュレーション（OpenMM、GBn2 + hydrogen mass repartitioning）
  で緩和し、リガンド重原子RMSDトラジェクトリがポケット内に留まるかどうかで
  再ランキングする——軽量なinduced-fitの健全性チェックである。

## インストール

meeko, pdbfixer, openmm, openmmforcefields, openff-toolkit, mdtraj
が必要（`pyproject.toml`の`dependencies`に明記済み）。加えてCPUドッキング
用に`qvina2` CLIバイナリ（conda-forgeの`qvina`パッケージ——QuickVina2、
AutoDock Vina 1.1.2を高速化したフォーク）がPATH上に必要である。これは
Pythonのimportではないため`pyproject.toml`の依存関係には書けない。
これらはconda-forge経由でのインストールが最善である。

```bash
mamba create -n dd_docking python=3.12 -c conda-forge \
    rdkit numpy pandas qvina meeko pdbfixer openmm openmmforcefields openff-toolkit mdtraj
mamba activate dd_docking
cd dd_docking
pip install -e .
```

これにより3つのコンソールコマンドがインストールされる: `dd_docking-prep`,
`dd_docking-dock`, `dd_docking-refine`。

### オプション: GPUアクセラレーションドッキング（Linux限定）

`dd_docking-dock`は、CPU版QuickVina2の代わりに
[Vina-GPU+](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0)を使用できる
（ビルド済みで、かつタスクのボックスサイズがOpenCLカーネルの制限に収まる
場合）。これはLinux限定の機能であり、macOS/Windows（またはバイナリを
ビルドしていないLinux）では、`dd_docking-dock`は自動的にCPU版QuickVina2に
フォールバックする——コードの変更は不要で、純粋にランタイムのフォール
バックである（後述の`--backend`を参照）。

```bash
mamba activate dd_docking
bash scripts/build_vina_gpu.sh
```

このスクリプトはVina-GPU+のソースを（固定コミットで）クローン・ビルドし
（`third_party/`配下、本リポジトリにはコミットしない）、生成された
バイナリとカーネルファイルを`$CONDA_PREFIX/share/dd_docking/vina-gpu/`に
インストールする。動作するOpenCLランタイムを持つNVIDIA/AMD GPUが必要
（`clinfo`で確認可能）——NVIDIAの場合は通常ドライバ/CUDA toolkitの
インストールに含まれる。AMD向け（`GPU_PLATFORM`）やOpenCLパスのカスタム
指定（`DD_DOCKING_OPENCL_PATH`）は`scripts/build_vina_gpu.sh`を参照。

Vina-GPU+カーネルの既知の制約: ドッキングボックスは各辺30Å未満である
必要がある。`dd_docking`はアンサンブルメンバーごとにこれをチェックし、
ボックスが大きすぎるメンバーについては自動的にCPUにフォールバックする
——後述のドッキングセクションの`--backend`を参照。

一部のGPU/ドライバの組み合わせでは、ビルド自体が成功していても
ランタイムでOpenCLカーネルのビルドに失敗する未解決の上流バグに遭遇する
（`CL_BUILD_PROGRAM_FAILURE`、またはカーネルコンパイル中のクラッシュ
——upstreamの[Issue #1](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0/issues/1)
と[Issue #26](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0/issues/26)
参照。本プロジェクトのGTX 1660 Tiでも再現した）。`dd_docking`はGPU側の
ドッキングタスク失敗を検知すると、1回だけ警告を出しつつ自動的にCPUで
再試行するため、結果は常に正しく得られる——ただし、ハードウェアによっては
GPUアクセラレーション自体が使えない場合がある。

## 使い方

### 1. アンサンブル準備（`dd_docking-prep`）

受容体コンフォメーションごとに `--member ID PDB_FILE COCRYSTAL_LIGAND_CODE`
を1回ずつ渡す。各メンバーはPDBFixerで個別に修復され、ドッキングボックスと
ポケット残基（デフォルトカットオフ5 Å）は共結晶リガンドの座標から導出
される。

```bash
dd_docking-prep \
  --member 6w63 data/raw/6W63.pdb X77 \
  --member 7l11 data/raw/7L11.pdb XF1 \
  --member 7l10 data/raw/7L10.pdb XEY \
  -o data/ensemble
```

出力（このリポジトリの `data/` から再現可能な実測値）:

```
6w63: 27 flexible residues -> data/ensemble/6w63_rigid.pdbqt
7l11: 24 flexible residues -> data/ensemble/7l11_rigid.pdbqt
7l10: 20 flexible residues -> data/ensemble/7l10_rigid.pdbqt

[done] 3 member(s) -> data/ensemble/manifest.json
```

`data/ensemble/` には、メンバーごとに `<id>_fixed.pdb`（PDBFixer適用後、
水素なし）、`<id>_rigid.pdbqt` / `<id>_flex.pdbqt`（Meeko出力）、そして
`manifest.json`（ボックス中心/サイズ、フレキシブル残基リスト、メンバーごとの
ファイルパス — 後続コマンドが読み込む）が生成される。

主なオプション:

| オプション | デフォルト | 説明 |
|---|---|---|
| `--chain` | `A` | 抽出する受容体鎖 |
| `--pocket-cutoff` | `5.0` | フレキシブル残基選定に用いる、共結晶リガンドからの距離（Å） |
| `--max-flexres` | `10` | フレキシブル残基数の上限（リガンドに近い順；`<=0`で上限なし） |
| `--box-padding` | `5.0` | 共結晶リガンドの範囲に対するドッキングボックスのパディング（Å） |
| `--charge-model` | `gasteiger` | Meeko受容体部分電荷モデル（`gasteiger`/`espaloma`/`zero`） |

`--max-flexres` のデフォルト値は経験に基づいて選ばれている。リジッド
ドッキング（事前計算済みのグリッドマップ）と異なり、フレキシブル側鎖の
1つ1つがVinaの探索空間に実際の回転自由度を追加する。本リポジトリの
Mpro試験データ（6W63/7L11/7L10、共結晶リガンドから5 Åのカットオフ）
では20〜27残基がカットオフ内に入り、上限なしで1リガンドを1メンバーに
対してドッキングすると1時間以内に終わらなかった。10（近い順）に制限
することで実用的な実行時間になった。「念のため」上限を外したり
カットオフを広げたりすると、容易にこの罠にはまる。

Python API:

```python
from dd_docking import prepare_ensemble

ensemble = prepare_ensemble(
    [("6w63", "data/raw/6W63.pdb", "X77"),
     ("7l11", "data/raw/7L11.pdb", "XF1"),
     ("7l10", "data/raw/7L10.pdb", "XEY")],
    "data/ensemble",
)
for m in ensemble:
    print(m.member_id, len(m.flexres), m.rigid_pdbqt)
```

### 2. アンサンブルドッキング（`dd_docking-dock`）

準備済みアンサンブルの全メンバーに対し、`.smi` ライブラリ（
`SMILES<TAB>ID` 形式、1行1分子）をドッキングする。

```bash
dd_docking-dock data/ensemble data/ligands.smi \
  -o data/screen --exhaustiveness 4 --n-poses 1 --n-jobs 4
```

作業は `(ligand, ensemble member)` のタスクに分割され（上記の例では
8リガンド × 3メンバー = 24タスク）、`--n-jobs`で並列化される（デフォルトは
逐次実行；`<=0`で全CPUコアを使用）。各リガンドが全メンバーに対して
ドッキングされるたびに進捗が1行出力される（本リポジトリのデータ、
16コアMac、`--n-jobs 4`での実測、天然化合物3種＋承認薬5種で約8分）:

```
[parallel] using 4 worker processes for 24 tasks
[1] X77_native_6w63  best_member=6w63  affinity=-6.837
[2] XF1_native_7l11  best_member=6w63  affinity=-7.697
[3] aspirin  best_member=7l11  affinity=-5.210
[4] XEY_native_7l10  best_member=7l11  affinity=-8.058
[5] ibuprofen  best_member=7l10  affinity=-5.282
[6] naproxen  best_member=7l10  affinity=-6.187
[7] acetaminophen  best_member=6w63  affinity=-4.376
[8] metformin  best_member=7l11  affinity=-5.539

[done] 8 ligand(s) ranked -> data/screen/ranked_results.csv
```

最終ランキング（`ranked_results.csv`、アフィニティ昇順）: 3つの共結晶
リガンド（セルフドッキング）が上位3位を占め、5つの無関係な承認薬は
それより下位にランキングされる（意図した弁別性能）。`best_member` 列と
メンバーごとの `affinity[6w63/7l11/7l10]` 列は行ごとに異なっており、
アンサンブルの各コンフォメーションが実際に異なる結果を生んでいることを
確認できる（アンサンブルドッキングが意図通り機能している証拠）:

| rank | ligand_id | best_member | best_affinity | affinity[6w63] | affinity[7l11] | affinity[7l10] |
|---|---|---|---|---|---|---|
| 1 | XEY_native_7l10 | 7l11 | -8.058 | -6.817 | -8.058 | -7.711 |
| 2 | XF1_native_7l11 | 6w63 | -7.697 | -7.697 | -7.414 | -7.418 |
| 3 | X77_native_6w63 | 6w63 | -6.837 | -6.837 | -6.773 | -6.554 |
| 4 | naproxen | 7l10 | -6.187 | -5.102 | -6.020 | -6.187 |
| 5 | metformin | 7l11 | -5.539 | -5.310 | -5.539 | -5.349 |
| 6 | ibuprofen | 7l10 | -5.282 | -4.388 | -4.739 | -5.282 |
| 7 | aspirin | 7l11 | -5.210 | -4.623 | -5.210 | -4.109 |
| 8 | acetaminophen | 6w63 | -4.376 | -4.376 | -4.022 | -4.355 |

出力:

- `data/screen/ranked_results.csv` — 列は `rank, ligand_id, smiles,
  best_member, best_affinity, receptor_pdb, pose_pdbqt,
  affinity[<member_id>]...`。`receptor_pdb` / `pose_pdbqt` は
  `dd_docking-refine` にそのまま渡される。
- `data/screen/top_hits.sdf` — 各ヒットの最良ポーズを `affinity` プロパティ
  付きで格納。このファイルは受容体PDB（`ranked_results.csv`の
  `receptor_pdb` 列）と一緒に蛋白-リガンドビューアへ直接読み込むことができる。
- `data/screen/ranked_results_poses/` — ヒットごとのPDBQTファイル
  （`dd_docking-refine` の入力）。

主なオプション:

| オプション | デフォルト | 説明 |
|---|---|---|
| `--exhaustiveness` | `16` | Vinaのexhaustiveness |
| `--n-poses` | `5` | メンバーごとに保持するポーズ数 |
| `--seed` | `0` | 乱数シード（埋め込みおよびドッキング） |
| `--top-n` | 全件 | 上位N件の結果のみ保持 |
| `--n-jobs` | `1` | `(ligand, member)` タスクごとの並列ワーカー数（`<=0`で全コア使用） |
| `--backend` | `auto` | ドッキングエンジン: `auto`はビルド済みかつメンバーのボックスがVina-GPU+の30Å未満制限に収まる場合にそれを使用、それ以外はCPU版QuickVina2。`cpu`は常にCPU版QuickVina2（全OSで動作）。`gpu`はVina-GPU+を優先し、使えないメンバーではCPUに警告付きでフォールバックする（[GPUアクセラレーションドッキング](#オプション-gpuアクセラレーションドッキングlinux限定)を参照） |
| `--no-progress` | - | 進捗ログを抑制 |

**`--n-jobs` の挙動とCPU割り当て**: `--n-jobs 1`（デフォルト）は逐次実行
され、各ドッキングタスクに内部的に利用可能な全コアを与える（`cpu=0`）。
それ以外の `--n-jobs` の値では、コアがワーカープロセス間で均等に分割
される（`cpu_count // n_jobs`、最小1）。フレキシブル側鎖ドッキングは
リジッドドッキングよりもはるかに重く、Vina自身のマルチスレッド化の恩恵を
大きく受けるため、軽量なリジッドドッキングでよく使われるパターン——
`--cpu 1` に固定してジョブ側を並列化する——はここでは大きく裏目に出る
（実測: 同条件で `cpu=1` 固定は1タスクあたり約850秒、複数コアを割り当てた
場合は1タスクあたり約150秒で、約5.6倍）。ワーカー数を増やしすぎると各
ワーカーのコア取り分が減り、これも遅くなる原因になり得るため、ライブラリ
規模とアンサンブル規模に応じて控えめな値（おおむねコア数の1/4〜1/2程度）
から始めるとよい。

Python API:

```python
from dd_docking import ligand_prep, screen_ensemble
from dd_docking.ensemble import load_manifest

ensemble = load_manifest("data/ensemble/manifest.json")
ligands = ligand_prep.read_smi("data/ligands.smi")
hits = screen_ensemble(ensemble, ligands, n_jobs=8,
                       out_csv="data/screen/ranked_results.csv",
                       out_sdf="data/screen/top_hits.sdf")
for hit in hits[:5]:
    print(hit.ligand_id, hit.best_member, hit.best_affinity, hit.member_affinities)
```

### 3. MD精密化と再ランキング（`dd_docking-refine`）

`dd_docking-dock` の `ranked_results.csv` から上位ヒットのみを、短時間の
暗黙溶媒MDシミュレーション（GBn2 + hydrogen mass repartitioning、4 fs）で
緩和・再評価する。これを（全件ではなく）上位N件に限定することで、計算
コストを実用的な範囲に抑えている。

```bash
dd_docking-refine data/screen/ranked_results.csv \
  -o data/screen/refine --top-n 2 --prod-ps 20 --equil-ps 5
```

実測出力（`--prod-ps` はここでは簡易テスト用に短縮している。本番では
より長い値を使用すること）:

```
[MD XEY_native_7l10] implicit solvent setup failed, retrying in vacuum: ''
[MD XEY_native_7l10] vacuum  rmsd_mean=1.90  rmsd_final=2.46  stable=True
[MD XF1_native_7l11] implicit  rmsd_mean=2.67  rmsd_final=3.32  stable=True

[done] 2 hit(s) refined -> data/screen/refine/md_rescore.csv
```

1件目のヒットではGBn2暗黙溶媒系のセットアップが失敗し（
`openmmforcefields` のフォースフィールドテンプレート生成における、一過性と
思われる不具合）、自動的に真空中MDにフォールバックした（`--vacuum` を
渡さなくてもこのフォールバックは自動的に発生する。想定内の挙動である）。
両ヒットとも安定性基準（RMSD平均 < 3 Å かつ最終値 < 4 Å）を満たし、
`stable=True` と判定された。

`md_rescore.csv` は、まず `stable`（リガンド重原子RMSD平均 < 3 Å かつ
最終値 < 4 Å）でソートし、その中でさらにVinaアフィニティ昇順でソート
する。列にはこのほか `rmsd_mean` / `rmsd_final` / `rmsd_max` /
`e_min_kcal` / `e_final_kcal` / `implicit`（GBn2暗黙溶媒を使用したか、
真空中にフォールバックしたか）が含まれる。各ヒットの複合体構造とトラジェ
クトリは `data/screen/refine/<rank>_<ligand_id>/`（`complex.pdb` +
`prod.dcd`）に保存される。

主なオプション:

| オプション | デフォルト | 説明 |
|---|---|---|
| `--top-n` | `5` | MDを実行する上位ヒット数 |
| `--prod-ps` | `100.0` | プロダクションランの長さ（ps） |
| `--equil-ps` | `20.0` | 平衡化の長さ（ps） |
| `--vacuum` | - | GBn2暗黙溶媒の代わりに真空中でMDを実行（暗黙溶媒のセットアップが失敗した場合の自動フォールバックでもある） |
| `--platform` | `auto` | OpenMMのプラットフォーム: `auto` はCUDA→OpenCLの順で優先し、どちらも使えない場合はCPUにフォールバックする。`CPU`/`CUDA`/`OpenCL`/`Reference` を指定して明示的に強制することも可能（そのプラットフォームが使えない場合はエラーになる） |

`--platform` はMDステップ（`dd_docking-refine`）にのみ適用される。
ドッキング自体のGPUオプションは`dd_docking-dock --backend`である
（[GPUアクセラレーションドッキング](#オプション-gpuアクセラレーションドッキングlinux限定)を参照）。

Python API:

```python
from dd_docking import refine_top_hits

result = refine_top_hits("data/screen/ranked_results.csv", "data/screen/refine",
                         top_n=3, prod_ps=100.0)
print(result[["name", "stable", "rmsd_mean", "best_affinity"]])
```

### エンドツーエンドパイプライン（Python API）

3つのステップをまとめて実行するには（各ステップは上記のように個別に
呼び出すこともできる）:

```python
from dd_docking.pipeline import run_ensemble_docking

df = run_ensemble_docking(
    "data/ensemble", "data/ligands.smi", "data/screen",
    exhaustiveness=8, n_jobs=8, refine=True, refine_top_n=3,
)
```

## 検証済みの挙動

- `data/ensemble/` にある3つのアンサンブルメンバー（SARS-CoV-2 Mpro、
  PDB 6W63/7L11/7L10 — それぞれ異なる阻害剤と結合した、実際に異なる
  ポケットコンフォメーション）のrigid/flex PDBQTファイルは、
  `qvina2 --receptor rigid.pdbqt --flex flex.pdbqt ...` で問題なく
  ドッキングできる。
- 各メンバーの共結晶リガンドをセルフドッキングすると、アフィニティに
  よって無関係な分子（承認薬）から弁別される（`data/ligands.smi` 参照）。
- `receptor_prep.py` の `regularize_carboxylate_geometry` は、PDBFixerの
  `addMissingAtoms()` の癖（新たに追加されたカルボキシラートのパートナー
  酸素——鎖末端のバックボーンOXT、あるいはAsp/GluのOD2/OE2——が兄弟酸素に
  対して化学的にありえない角度に配置される）を修正する。この不具合は
  未ブロックの末端を持つ実際のPDBエントリにおいて、以前はMeekoの
  `mk_prepare_receptor.py` を酸素価数エラーでクラッシュさせていた。
  PDB 4EQC（PAK1キナーゼドメイン）でエンドツーエンドの修正を確認済み:
  `dd_docking-prep` は以前失敗していたところで成功するようになり
  （3件の欠陥を検出・是正）、`dd_docking-dock` はNaringinを修復済み
  受容体に対して正常にドッキングする（最良アフィニティ -8.646 kcal/mol、
  フレキシブル側鎖使用）。

## モジュール構成（`dd_docking/`）

| モジュール | 役割 |
|---|---|
| `receptor_prep.py` | PDB取得/鎖単離/残基クリーンアップ（TER挿入、CYXリネーム）、PDBFixerベースの修復、PDBFixer適用後のカルボキシラート幾何是正（`regularize_carboxylate_geometry`） |
| `pocket.py` | ドッキングボックス計算、距離ベースのフレキシブル残基検出、Meeko flexres文字列フォーマット |
| `ensemble.py` | 複数コンフォメーションにまたがるreceptor_prep + pocket + Meeko PDBQT生成のバッチ処理、`manifest.json` としての保存/読み込み |
| `ligand_prep.py` | `.smi` 読み込み、SMILES -> 3D（ETKDGv3+MMFF）-> Meekoリガンド PDBQT |
| `docking.py` | QuickVina2（`qvina2` CLI）経由のCPUドッキング、フレキシブル受容体対応（`dock_ligand`） |
| `gpu_backend.py` | オプションのVina-GPU+バックエンド（Linux限定）: バックエンド選択（`resolve_backend`）、subprocess経由のドッキング（`dock_ligand_gpu`） |
| `screening.py` | 全リガンド × 全メンバーの並列ドッキング、best-of-ensembleランキング、CSV/SDF出力 |
| `refine_md.py` | 上位ヒットの短時間暗黙溶媒MD緩和、RMSD安定性評価、再ランキング |
| `io.py` | PDBQT -> RDKit変換、ポーズSDF出力、結果CSV I/O |
| `parallel.py` | `parallel_map` — 独立したタスクを `ProcessPoolExecutor` 上で並列化 |
| `progress.py` | `DockProgress` / `RefineProgress` — 完了タスクごとに1行出力 |
| `pipeline.py` | prep -> docking -> MD精密化を連鎖させる高レベル関数群 |
| `cli.py` | `dd_docking-prep` / `dd_docking-dock` / `dd_docking-refine` コマンド |

## 制約と拡張の可能性

- フレキシブル残基は共結晶リガンドからの距離のみで選ばれる。induced fit
  への関与が既知の特定残基を固定したい場合は、`pocket.find_pocket_residues`
  を経由せず、flexres文字列を直接 `ensemble.prepare_ensemble_member` に
  渡すこと。
- MD精密化のデフォルトは短時間・CPUのみ・暗黙溶媒（GBn2）シミュレーション
  である。明示溶媒やGPUベースの自由エネルギー計算はスコープ外
  （必要であれば `refine_md.py` を独立して拡張すること）。
- `screening.py` は全ての `(ligand, member)` の組み合わせにわたって並列化
  するため、Vinaのグリッドマップはタスクごとに再計算される。大規模な
  ライブラリでは、各メンバーのグリッドマップを一度だけ計算して再利用する
  （メンバー単位で並列化する）方が最適化として有効だろう。

## ライセンス

MIT — [LICENSE](LICENSE) を参照。
