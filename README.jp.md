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

### オプション: GPUアクセラレーションドッキング（Linux限定、リジッド受容体限定）

`dd_docking-dock`は、CPU版QuickVina2の代わりに
[Vina-GPU+](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0)を使用できる
（ビルド済みで、対象メンバーにフレキシブル側鎖が無く、かつタスクの
ボックスサイズがOpenCLカーネルの制限に収まる場合）。これはLinux限定の
機能であり、macOS/Windows（またはバイナリをビルドしていないLinux、
フレキシブル残基を持つメンバー、大きすぎるボックス）では、
`dd_docking-dock`は自動的にCPU版QuickVina2にフォールバックする——コードの
変更は不要で、純粋にランタイムのフォールバックである（後述の
`--backend`を参照）。

```bash
mamba activate dd_docking
bash scripts/build_vina_gpu.sh
```

このスクリプトはVina-GPU+のソースを（固定コミットで）クローンし
（`third_party/`配下、本リポジトリにはコミットしない）、**カーネル
ソースに実在するOpenCL Cの不具合をいくつかパッチしてから**（後述）
ビルドし、生成されたバイナリとカーネルファイルを
`$CONDA_PREFIX/share/dd_docking/vina-gpu/`にインストールする。動作する
OpenCLランタイムを持つNVIDIA/AMD GPUが必要（`clinfo`で確認可能）——
NVIDIAの場合は通常ドライバ/CUDA toolkitのインストールに含まれる。
AMD向け（`GPU_PLATFORM`）、OpenCLパスのカスタム指定
（`DD_DOCKING_OPENCL_PATH`）、OpenCL Cバージョンの指定
（`DD_DOCKING_OPENCL_VERSION`）は`scripts/build_vina_gpu.sh`を参照。

**修正不可能な根本的制約: Vina-GPU+はリジッド受容体しかサポートしない。**
`main_procedure_cl.cpp`はドッキング前に`m.num_other_pairs() == 0`という
アサーションを持つが、この値はligand-flex・flex-flex・flex-inflexの
相互作用が1つでもあれば非ゼロになる——つまりフレキシブル側鎖が1つでも
あれば即座に満たされない。`dd_docking`のアンサンブルドッキングは常に
フレキシブル側鎖を使うため、`--backend gpu`/`auto`はそれらのタスクでは
常にCPU版QuickVina2を使う（`gpu`を明示指定した場合は1回だけ警告する）。
フレキシブル残基が0のメンバーのみが実際にGPUで動作する。

以前のこのREADMEでは、Vina-GPU+は本プロジェクトのGTX 1660 Tiでは単純に
動作しない（`CL_BUILD_PROGRAM_FAILURE`、あるいはカーネルコンパイル中の
クラッシュ——upstreamの[Issue #1](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0/issues/1)
と[Issue #26](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0/issues/26)
と同じ症状）と報告していた。実際の根本原因: `clinfo`はこのデバイスの
実際のOpenCL Cコンパイラレベルが1.2であることを示している（*プラット
フォーム*自体はOpenCL 3.0対応と表示するにもかかわらず）。Vina-GPU+の
カーネルソースには実際の型エラーが複数存在する（行ポインタを誤った
ポインタ型に変換する冗長な`&`、`__global`アドレス空間修飾子が欠けている
ポインタ2箇所、そしてOpenCL 2.0以降専用の`get_global_linear_id()`の呼び出し）。
これらは厳格な1.2/2.0コンパイルでは正しく拒否されるが、NVIDIAのより緩い
OpenCL 3.0コンパイルパスでは（プログラムバイナリの段階でセグフォルトや
失敗を起こすまでの間だけ）通ってしまう。`build_vina_gpu.sh`は今これらを
ビルド前にパッチし、デフォルトを`-DOPENCL_1_2`にしており、この環境では
正しくビルド・動作する。修正後の実測では、GPU/subprocess起動コストを
償却できるだけの作業量があれば、リジッド受容体のGPUドッキングは実際に
CPUより速くなる——8リガンドを`--exhaustiveness 32`でドッキングした場合:
**GPUで16.2秒、CPUで35.2秒**（約2.2倍）。ただし単発の軽いタスクでは逆に
遅くなる（1リガンド・`--exhaustiveness 8`: GPU 3.2秒 vs CPU 2.3秒）ため、
小さなジョブでの高速化は期待しないこと。お使いのGPU/ドライバのOpenCL C
コンパイラが本当に2.0/3.0に対応しているなら、
`DD_DOCKING_OPENCL_VERSION=-DOPENCL_3_0`でより速いカーネルがビルドできる
かもしれない——これらのパッチはそのバージョン選択とは独立した正しさの
修正なので、どちらを選んでも適用される。

## 使い方

### 1. アンサンブル準備（`dd_docking-prep`）

受容体コンフォメーションごとに `--member ID PDB_FILE COCRYSTAL_LIGAND_CODE`
を1回ずつ渡す。各メンバーはPDBFixerで個別に修復され、ドッキングボックスと
ポケット残基（デフォルトカットオフ5 Å）は共結晶リガンドの座標から導出
される。

```bash
dd_docking-prep \
  --member 3ert data/raw/3ERT.pdb OHT \
  --member 1xpc data/raw/1XPC.pdb AIT \
  --member 1yim data/raw/1YIM.pdb CM4 \
  -o data/ensemble
```

出力（このリポジトリの `data/` から再現可能な実測値）:

```
3ert: 10 flexible residues -> data/ensemble/3ert_rigid.pdbqt
1xpc: 10 flexible residues -> data/ensemble/1xpc_rigid.pdbqt
1yim: 10 flexible residues -> data/ensemble/1yim_rigid.pdbqt

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
ERα試験データ（3ERT/1XPC/1YIM、共結晶リガンドから5 Åのカットオフ）
では22〜26残基がカットオフ内に入り、上限なしで1リガンドを1メンバーに
対してドッキングすると1時間以内に終わらなかった。10（近い順）に制限
することで実用的な実行時間になった。「念のため」上限を外したり
カットオフを広げたりすると、容易にこの罠にはまる。

ドッキングボックス自体は、共結晶リガンド*だけでなく*選ばれたフレキシブル
残基の原子も（`--box-padding`を加えた上で）カバーするように計算される。
リガンドの周囲だけでボックスを決めると、フレキシブル側鎖の可動原子が
探索空間の外に出てしまい、Vinaが「no conformations completely within the
search space」というエラーを出し、そのメンバーに対するすべてのリガンドで
ポーズが0件になることがある。ERαのフレキシブル残基は伸長したポケットの
周囲に比較的広く分布しているため、既定のパディングでもこの3メンバーは
ボックスが各辺25〜33Å程度になる（`1xpc`/`1yim`はVina-GPU+の別枠の30Å
OpenCLカーネル制限にも実際に引っかかっているが、ここでは無関係——3メンバー
とも全てフレキシブル側鎖を持つため、`--backend auto/gpu`はこのデータセット
では常にCPU版QuickVina2を使う。[GPUアクセラレーションドッキング](#オプション-gpuアクセラレーションドッキングlinux限定リジッド受容体限定)を参照）。

Python API:

```python
from dd_docking import prepare_ensemble

ensemble = prepare_ensemble(
    [("3ert", "data/raw/3ERT.pdb", "OHT"),
     ("1xpc", "data/raw/1XPC.pdb", "AIT"),
     ("1yim", "data/raw/1YIM.pdb", "CM4")],
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

ここでの`data/ligands.smi`は、各メンバーの共結晶リガンド3つ（セルフ
ドッキング確認用）＋ERα実活性化合物3つ＋物性が近いdecoy化合物9つ
（いずれも[DUD-E](https://dude.docking.org/targets/esr1/)から取得。後述の
[テストデータセット（ERα, DUD-Eデコイ）](#テストデータセットerα-dud-eデコイ)
を参照）から成る——無関係な承認薬を「明らかに違う」decoyとして数個混ぜる
より、はるかに厳しく現実的なスクリーニングになる。

作業は `(ligand, ensemble member)` のタスクに分割され（上記の例では
15リガンド × 3メンバー = 45タスク）、`--n-jobs`で並列化される（デフォルトは
逐次実行；`<=0`で全CPUコアを使用）。各リガンドが全メンバーに対して
ドッキングされるたびに進捗が1行出力される（本リポジトリのデータ、
16コアLinuxマシン、`--n-jobs 4`での実測、約29分）:

```
[parallel] using 4 worker processes for 45 tasks
[1] OHT_native_3ert  best_member=1xpc  affinity=-11.300
[2] AIT_native_1xpc  best_member=3ert  affinity=-11.400
[3] CM4_native_1yim  best_member=3ert  affinity=-11.700
[4] CHEMBL215428_active  best_member=1xpc  affinity=-11.700
[5] CHEMBL31127_active  best_member=1xpc  affinity=-11.300
[6] CHEMBL385358_active  best_member=1xpc  affinity=-12.100
[7] C48469840_decoy  best_member=1xpc  affinity=-10.600
[8] C14244874_decoy  best_member=3ert  affinity=-11.200
[9] C66896427_decoy  best_member=1yim  affinity=-11.100
[10] C27846680_decoy  best_member=1xpc  affinity=-11.300
[11] C37195230_decoy  best_member=1xpc  affinity=-8.600
[12] C37085981_decoy  best_member=1xpc  affinity=-9.000
[13] C09495369_decoy  best_member=3ert  affinity=-12.300
[14] C12776766_decoy  best_member=1xpc  affinity=-11.300
[15] C36904163_decoy  best_member=1xpc  affinity=-8.300

[done] 15 ligand(s) ranked -> data/screen/ranked_results.csv
```

最終ランキング（`ranked_results.csv`、アフィニティ昇順）。無関係な承認薬
とは違い、DUD-Eの物性が近いdecoy（実活性化合物と分子量・回転可能結合数・
電荷が近いが2次元構造は異なる）は意図的に難しい弁別テストであり、結果も
それを正直に反映している：1位はdecoy（`C09495369_decoy`）で、上位12行の
中に他のいくつかのdecoyも実活性化合物・天然リガンドと混在している。これは
バグではない——低exhaustivenessの単発AutoDock系ドッキングだけでは、物性が
近いdecoyを真の活性化合物から完全に分離できないことが多いという、よく
知られた現実であり、まさにこれが`dd_docking-refine`のMDベースinduced-fit
再ランキングが第二のフィルタとして存在する理由である（後述）。一方で
ドッキングが明確に達成していること：最も弱く、最も構造的に異なる3つの
decoy（順位13〜15、アフィニティ-9.0〜-8.3）は、他のすべて（-10.6以上）
からきれいに分離されている。また`best_member`列とメンバーごとの
`affinity[3ert/1xpc/1yim]`列は行ごとに異なっており、アンサンブルの各
コンフォメーションが実際に異なる結果を生んでいることを確認できる
（アンサンブルドッキングが意図通り機能している証拠）:

| rank | ligand_id | best_member | best_affinity | affinity[3ert] | affinity[1xpc] | affinity[1yim] |
|---|---|---|---|---|---|---|
| 1 | C09495369_decoy | 3ert | -12.3 | -12.3 | -12.2 | -11.4 |
| 2 | CHEMBL385358_active | 1xpc | -12.1 | -11.5 | -12.1 | -11.6 |
| 3 | CM4_native_1yim | 3ert | -11.7 | -11.7 | -11.7 | -10.8 |
| 4 | CHEMBL215428_active | 1xpc | -11.7 | -11.4 | -11.7 | -8.1 |
| 5 | AIT_native_1xpc | 3ert | -11.4 | -11.4 | -11.3 | -10.8 |
| 6 | OHT_native_3ert | 1xpc | -11.3 | -10.6 | -11.3 | -10.8 |
| 7 | CHEMBL31127_active | 1xpc | -11.3 | -10.7 | -11.3 | -10.5 |
| 8 | C27846680_decoy | 1xpc | -11.3 | -10.7 | -11.3 | -10.7 |
| 9 | C12776766_decoy | 1xpc | -11.3 | -10.5 | -11.3 | -10.1 |
| 10 | C14244874_decoy | 3ert | -11.2 | -11.2 | -11.0 | -10.4 |
| 11 | C66896427_decoy | 1yim | -11.1 | -11.0 | -10.2 | -11.1 |
| 12 | C48469840_decoy | 1xpc | -10.6 | -9.4 | -10.6 | -10.4 |
| 13 | C37085981_decoy | 1xpc | -9.0 | -8.2 | -9.0 | -9.0 |
| 14 | C37195230_decoy | 1xpc | -8.6 | -8.4 | -8.6 | -8.3 |
| 15 | C36904163_decoy | 1xpc | -8.3 | -7.6 | -8.3 | -8.1 |

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
| `--backend` | `auto` | ドッキングエンジン: `auto`はビルド済みかつメンバーのボックスがVina-GPU+の30Å未満制限に収まる場合にそれを使用、それ以外はCPU版QuickVina2。`cpu`は常にCPU版QuickVina2（全OSで動作）。`gpu`はVina-GPU+を優先し、使えないメンバーではCPUに警告付きでフォールバックする（[GPUアクセラレーションドッキング](#オプション-gpuアクセラレーションドッキングlinux限定リジッド受容体限定)を参照） |
| `--no-progress` | - | 進捗ログを抑制 |

**`--n-jobs` の挙動とCPU割り当て**: `--n-jobs 1`（デフォルト）は逐次実行
され、各ドッキングタスクに内部的に利用可能な全コアを与える（`cpu=0`）。
それ以外の `--n-jobs` の値では、コアがワーカープロセス間で均等に分割
される（`cpu_count // n_jobs`、最小1）。フレキシブル側鎖ドッキングは
リジッドドッキングよりもはるかに重く、QuickVina2自身のマルチスレッド化の
恩恵を大きく受けるため、軽量なリジッドドッキングでよく使われるパターン——
`--cpu 1` に固定してジョブ側を並列化する——はここでは大きく裏目に出る。
本リポジトリのERαデータ（1リガンド・1メンバー、`--exhaustiveness 16`、
16コアLinuxマシン）で実測した、`--cpu`ごとの単一タスクの実行時間:

| `--cpu` | 実行時間 | `--cpu 1`比の高速化 |
|---|---|---|
| 1 | 827.5秒 | 1.0倍 |
| 2 | 432.4秒 | 1.9倍 |
| 4 | 223.2秒 | 3.7倍 |
| 8 | 118.8秒 | 7.0倍 |
| 16 | 101.5秒 | 8.2倍 |

8コアを超えると収益が急激に減少する（8→16コアで得られる高速化はわずか
約17%）ため、多数の`(ligand, member)`タスクをまとめて処理するバッチでは、
余ったコアを1タスクに全部割り当てるより、いくつかの`--n-jobs`ワーカーに
分けて割り当てる（それぞれが`cpu_count // n_jobs`で複数コアを保持する）
方がマシンをうまく使える傾向にある。まずは控えめに（`--cpu`を1タスクあたり
4〜8程度、つまり`--n-jobs`を`cpu_count // 6`程度）始め、ライブラリ規模と
アンサンブル規模に応じて調整するとよい。

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
  -o data/screen/refine --top-n 3 --prod-ps 20 --equil-ps 5
```

実測出力（`--prod-ps` はここでは簡易テスト用に短縮している。本番では
より長い値を使用すること）:

```
[MD C09495369_decoy] implicit  platform=OpenCL  rmsd_mean=1.41  rmsd_final=1.35  stable=True
[MD CHEMBL385358_active] implicit  platform=OpenCL  rmsd_mean=0.94  rmsd_final=1.79  stable=True
[MD CM4_native_1yim] implicit  platform=OpenCL  rmsd_mean=1.48  rmsd_final=1.89  stable=True

[done] 3 hit(s) refined -> data/screen/refine/md_rescore.csv
```

上位3ヒット全て——ドッキングで全ての実活性化合物を上回った例のdecoy
（`C09495369_decoy`）も含めて——ここではMDの安定性基準（RMSD平均 < 3 Å
かつ最終値 < 4 Å）を満たし、GBn2暗黙溶媒系のセットアップも3件すべて成功し
（`--platform`はOpenCLを選択）、RMSDは3件とも2Å未満に収まった。これを
「見栄えの良い」結果に書き換えるのではなく、あえてそのまま示している。
これは物性が近いdecoyを使うことの正直で想定内の結果であり、バグではない。
DUD-Eの「decoy」はChEMBLに活性が報告されていないというだけの推定不活性
であって、不活性の証明ではない——一部は実際に、真のリガンドと同じように
互換性のあるサブポケットに安定して収まってしまう。まさにこれが、
ドッキング分野でproperty-matched decoyセットが本当に難しいベンチマークと
見なされている理由であり、このパイプラインのどの段階（ドッキングもMDも）
単独では最終的な答えとして信頼すべきではない理由でもある。

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
（[GPUアクセラレーションドッキング](#オプション-gpuアクセラレーションドッキングlinux限定リジッド受容体限定)を参照）。

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

## テストデータセット（ERα, DUD-Eデコイ）

`data/raw/` と `data/ligands.smi` には、特定の標的に縛られない、小規模だが
本格的なサンプルが入っている——複数の共結晶コンフォメーションとDUD-E
（または類似の）decoyセットを持つ任意の蛋白に入れ替え可能である:

- **アンサンブルメンバー**: ヒトエストロゲン受容体α（ERα）のリガンド結合
  ドメイン3構造。induced-fitの教科書的な例（アゴニストとSERMでヘリックス
  12が再配置され、ポケットの形が変わる）——
  [3ERT](https://www.rcsb.org/structure/3ERT)（4-ヒドロキシタモキシフェン）、
  [1XPC](https://www.rcsb.org/structure/1XPC)と
  [1YIM](https://www.rcsb.org/structure/1YIM)（異なる2種のSERM骨格）。
  いずれも野生型・単一鎖であり、リガンド結合ドメイン内に内部欠損（次節が
  重要視する理由を参照）がないことをRCSB APIで確認した上で選定した。
- **リガンドライブラリ**: 各メンバー自身の共結晶リガンド（セルフドッキング
  確認用）＋ERαの実活性化合物3つ＋物性が近いdecoy化合物9つ、いずれも
  [DUD-Eの`esr1`ターゲット](https://dude.docking.org/targets/esr1/)から
  取得し、RDKitで有効かつ3D埋め込み可能なSMILESのみを事前フィルタした
  （再現性のため固定シードでサブセット抽出）。DUD-Eのdecoyは、実活性化合物
  と分子量・回転可能結合数・電荷は似ているが2次元構造（トポロジー）が
  異なるように選ばれており、無関係な日常的な分子（例: 承認薬）よりも
  はるかに難しい弁別テストになる——承認薬はどんなドッキング手法でも
  真のリガンドから簡単に分離できてしまう。

## 検証済みの挙動

- `data/ensemble/` にある3つのアンサンブルメンバーのrigid/flex PDBQTファイル
  は`qvina2 --receptor rigid.pdbqt --flex flex.pdbqt ...`で問題なくドッキング
  できる。ボックス（リガンドの範囲＋全フレキシブル残基の原子＋パディング、
  `pocket.compute_box`の`extra_coords`を参照）はすべてのフレキシブル側鎖の
  可動原子を含んでおり、これは`dd_docking-dock`が全ての`(ligand, member)`
  タスクで実際のポーズ・アフィニティを返すこと（Vinaの「no conformations
  completely within the search space」失敗ではなく）で確認済みである。
- 各メンバー自身の共結晶リガンド、実際のERα活性化合物、DUD-Eの物性が近い
  decoyを一緒にドッキングすると、最も弱く最も構造的に異なるdecoyは下位に
  きれいに分離される一方、ドッキング単独ではより上位で全てのdecoyを全ての
  真の活性化合物/天然リガンドから完全には分離できず、この特定のケースでは
  MDベースのinduced-fit再ランキングもそれを解決しない（上記の実例を参照）。
  これは物性が近いベンチマークの正直に報告された想定内の挙動であり、
  隠すべき弁別性能の失敗ではない。
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
| `pocket.py` | ドッキングボックス計算（リガンド範囲＋`extra_coords`によるフレキシブル残基原子の任意追加）、距離ベースのフレキシブル残基検出、Meeko flexres文字列フォーマット |
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
