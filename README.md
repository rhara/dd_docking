# docking

アンサンブルドッキング + 誘導適合（インデュースドフィット）風のポケット
動的挙動を取り入れたバーチャルスクリーニングツールキット。特定の標的・
リガンド群に依存しない、任意のターゲットに使い回せる再利用可能パッケージ
として設計している（`rocslib` / `omegalib` / `plviewerlib` と同じ方針）。

- **アンサンブル準備 (`dockinglib-prep`)**: 複数の受容体コンフォメーション
  （PDB）をそれぞれ PDBFixer で構造修正し、共結晶リガンドからドッキング
  ボックスとポケット近傍のフレキシブル残基を自動決定して、Meeko で
  rigid/flex の PDBQT ペアを生成する。
- **アンサンブルドッキング (`dockinglib-dock`)**: `.smi` ライブラリの各分子
  を、アンサンブルの全メンバーに対してフレキシブル側鎖ドッキングし、
  メンバー横断で最良（最も低い）affinity によりランキングする
  (best-of-ensemble)。
- **MD後処理 (`dockinglib-refine`)**: 上位ヒットのみ、OpenMM による短時間の
  陰溶媒 MD（GBn2 + 水素質量再配分）で複合体を緩和し、リガンド重原子RMSD
  の時系列から「ポーズがポケットに留まるか」を評価して再ランキングする
  （誘導適合の簡易的な後段検証）。

## セットアップ

vina / meeko / pdbfixer / openmm / openmmforcefields / openff-toolkit /
mdtraj が入った conda env `mpro` を使う（新規 env 作成は不要）。

```bash
cd ~/work/docking
conda run -n mpro pip install -e .
```

インストールすると `dockinglib-prep` / `dockinglib-dock` / `dockinglib-refine`
の3コンソールコマンドが `~/miniforge3/envs/mpro/bin/` に入る。`conda activate
mpro` 済みなら PATH が通っているのでそのまま呼べる（以下の例はすべて
`conda run -n mpro <command>` の形で書くが、activate 済みならその接頭辞は不要）。

## 使い方

### 1. アンサンブル準備 (`dockinglib-prep`)

`--member ID PDBファイル 共結晶リガンドの3文字コード` を受容体コンフォ
メーションの数だけ繰り返し指定する。各メンバーは個別に PDBFixer で修正
され、共結晶リガンド座標からドッキングボックスとポケット残基（既定
カットオフ5Å以内）が決まる。

```bash
conda run -n mpro dockinglib-prep \
  --member 6w63 data/raw/6W63.pdb X77 \
  --member 7l11 data/raw/7L11.pdb XF1 \
  --member 7l10 data/raw/7L10.pdb XEY \
  -o data/ensemble
```

出力（実測、このリポジトリの `data/` で再現可能）:

```
6w63: 27 flexible residues -> data/ensemble/6w63_rigid.pdbqt
7l11: 24 flexible residues -> data/ensemble/7l11_rigid.pdbqt
7l10: 20 flexible residues -> data/ensemble/7l10_rigid.pdbqt

[done] 3 member(s) -> data/ensemble/manifest.json
```

`data/ensemble/` にメンバーごとの `<id>_fixed.pdb`（PDBFixer後、水素なし）、
`<id>_rigid.pdbqt` / `<id>_flex.pdbqt`（Meeko生成）、`manifest.json`（各
メンバーの box 中心/サイズ・フレキシブル残基一覧・ファイルパスをまとめた
もの、後続コマンドが読む）が書き出される。

主なオプション:

| オプション | 既定値 | 内容 |
|---|---|---|
| `--chain` | `A` | 抽出する受容体チェイン |
| `--pocket-cutoff` | `5.0` | 共結晶リガンドからのフレキシブル残基判定距離 (Å) |
| `--max-flexres` | `10` | フレキシブル残基の上限（リガンドに近い順に採用。`0`以下で上限なし） |
| `--box-padding` | `5.0` | 共結晶リガンド範囲に対するドッキングボックスの余白 (Å) |
| `--charge-model` | `gasteiger` | Meeko 受容体部分電荷モデル (`gasteiger`/`espaloma`/`zero`) |

`--max-flexres` は実測に基づく既定値: リジッドドッキングと違い、フレキシブル
側鎖は1残基ごとに本物の回転自由度がVinaの探索空間に追加される（グリッド
マップで済むリジッド部分とは違う）。このリポジトリのMpro実データ（6W63/
7L11/7L10、共結晶リガンドからの距離カットオフ5Å）ではカットオフ内に
20〜27残基もヒットし、上限なしで1リガンド×1メンバーのドッキングが小一時間
経っても終わらなかった。上限10（リガンドに近い順）にしたところ現実的な
時間に収まった。「なるべく広く効かせたい」という直感で上限を外したり
カットオフを広げたりすると、この罠に容易にはまるので注意。

Python API:

```python
from dockinglib import prepare_ensemble

ensemble = prepare_ensemble(
    [("6w63", "data/raw/6W63.pdb", "X77"),
     ("7l11", "data/raw/7L11.pdb", "XF1"),
     ("7l10", "data/raw/7L10.pdb", "XEY")],
    "data/ensemble",
)
for m in ensemble:
    print(m.member_id, len(m.flexres), m.rigid_pdbqt)
```

### 2. アンサンブルドッキング (`dockinglib-dock`)

`.smi`（`SMILES<TAB>ID` 形式、1行1分子）のライブラリを、準備済みアンサ
ンブルの全メンバーに対してドッキングする。

```bash
conda run -n mpro dockinglib-dock data/ensemble data/ligands.smi \
  -o data/screen --exhaustiveness 8 --n-poses 3 --n-jobs 8
```

処理は `(リガンド, アンサンブルメンバー)` の全組み合わせ（上の例なら
8リガンド×3メンバー=24タスク）に分解され、`--n-jobs` で並列化できる
（既定は逐次実行。`<=0` で全CPUコア使用）。1リガンドの全メンバーへの
ドッキングが完了するたびに1行進捗が出る:

```
[1] aspirin  best_member=6w63  affinity=-6.421
[2] X77_native_6w63  best_member=6w63  affinity=-9.883
...

[done] 8 ligand(s) ranked -> data/screen/ranked_results.csv
```

出力:

- `data/screen/ranked_results.csv` — `rank, ligand_id, smiles, best_member,
  best_affinity, receptor_pdb, pose_pdbqt, affinity[<member_id>]...` の列。
  `receptor_pdb` / `pose_pdbqt` は `dockinglib-refine` がそのまま読む。
- `data/screen/top_hits.sdf` — 各ヒットの最良ポーズ（`affinity` プロパティ
  付き）。**このプロパティ名は既存の蛋白-リガンドビューア
  (`~/work/viewer`, `plviewerlib`) がスコアとして自動認識する名前**なので、
  受容体 PDB (`ranked_results.csv` の `receptor_pdb` 列) と合わせてビュー
  アにそのままロードできる。
- `data/screen/ranked_results_poses/` — ヒットごとの個別 PDBQT ファイル
  （`dockinglib-refine` の入力）。

主なオプション:

| オプション | 既定値 | 内容 |
|---|---|---|
| `--exhaustiveness` | `16` | Vina exhaustiveness |
| `--n-poses` | `5` | メンバーごとに保持するポーズ数 |
| `--seed` | `0` | 乱数シード（埋め込み・ドッキング共通） |
| `--top-n` | 全件 | 上位N件のみ結果に残す |
| `--n-jobs` | `1` | 並列ワーカー数（`(リガンド,メンバー)` タスク単位、`<=0`で全コア） |
| `--no-progress` | - | 進捗ログを消す |

**`--n-jobs` の挙動とCPU割り当て**: `--n-jobs 1`（既定）は逐次実行で、各
ドッキングタスクにVinaが内部的に使えるコアを全部渡す（`cpu=0`）。
`--n-jobs`をそれ以外にすると、ワーカープロセス数ぶんコアを均等割り
（`(CPUコア数) // ワーカー数`、最低1）してVinaに渡す。フレキシブル側鎖
ドッキングはリジッドドッキングよりずっと重く、Vina自身のマルチスレッド化
が効きやすいため、`--cpu 1`固定でジョブだけ並列化する定石（軽いリジッド
ドッキング向け、`.archives/mpro/.archives/vspipe/dock.py`の方式）をそのまま
持ち込むと大幅に遅くなる（実測: 同条件でタスクあたり cpu=1固定は約850秒、
cpu=数コアなら約150秒 — 約5.6倍の差）。ワーカー数を上げすぎるとコア当たりの
割り当てが減って逆に遅くなるので、リガンド数・アンサンブルサイズに対して
`--n-jobs`は控えめに（コア数の1/4〜1/2程度)から試すとよい。

Python API:

```python
from dockinglib import ligand_prep, screen_ensemble
from dockinglib.ensemble import load_manifest

ensemble = load_manifest("data/ensemble/manifest.json")
ligands = ligand_prep.read_smi("data/ligands.smi")
hits = screen_ensemble(ensemble, ligands, n_jobs=8,
                       out_csv="data/screen/ranked_results.csv",
                       out_sdf="data/screen/top_hits.sdf")
for hit in hits[:5]:
    print(hit.ligand_id, hit.best_member, hit.best_affinity, hit.member_affinities)
```

### 3. MD後処理・再ランキング (`dockinglib-refine`)

`dockinglib-dock` が書き出した `ranked_results.csv` の上位ヒットだけを、
短時間の陰溶媒MD（GBn2 + 水素質量再配分・4fs）で緩和・再評価する。全ヒット
ではなく上位N件に限定することで、計算コストを現実的に抑えている。

```bash
conda run -n mpro dockinglib-refine data/screen/ranked_results.csv \
  -o data/screen/refine --top-n 3 --prod-ps 100 --equil-ps 20
```

```
[MD X77_native_6w63] implicit  rmsd_mean=0.87  rmsd_final=1.02  stable=True
[MD ibuprofen] implicit  rmsd_mean=4.13  rmsd_final=5.88  stable=False
...

[done] 3 hit(s) refined -> data/screen/refine/md_rescore.csv
```

`md_rescore.csv` は `stable`（リガンド重原子RMSD平均<3Å かつ終端<4Å）を
最優先し、その中で Vina affinity 昇順に再ソートされる。列には
`rmsd_mean` / `rmsd_final` / `rmsd_max` / `e_min_kcal` / `e_final_kcal` /
`implicit`（GBn2陰溶媒か真空フォールバックか）も含まれる。各ヒットの
複合体構造・軌道は `data/screen/refine/<rank>_<ligand_id>/`
（`complex.pdb` + `prod.dcd`）に保存される。

主なオプション:

| オプション | 既定値 | 内容 |
|---|---|---|
| `--top-n` | `5` | MD対象にする上位ヒット数 |
| `--prod-ps` | `100.0` | 本計算の長さ (ps) |
| `--equil-ps` | `20.0` | 平衡化の長さ (ps) |
| `--vacuum` | - | GBn2陰溶媒の代わりに真空でMD（陰溶媒系の構築に失敗した場合は自動でこのモードにフォールバックする） |

Python API:

```python
from dockinglib import refine_top_hits

result = refine_top_hits("data/screen/ranked_results.csv", "data/screen/refine",
                         top_n=3, prod_ps=100.0)
print(result[["name", "stable", "rmsd_mean", "best_affinity"]])
```

### 一気通貫パイプライン (Python API)

3ステップをまとめて実行したい場合は `dockinglib.pipeline` を使う（各
ステップは上記の通り個別にも呼べる）:

```python
from dockinglib.pipeline import run_ensemble_docking

df = run_ensemble_docking(
    "data/ensemble", "data/ligands.smi", "data/screen",
    exhaustiveness=8, n_jobs=8, refine=True, refine_top_n=3,
)
```

## 検証済みの動作

- `data/ensemble/` の3メンバー（SARS-CoV-2 Mpro、6W63/7L11/7L10 — 異なる
  阻害剤が結合した実際に異なるポケットコンフォメーション）で、rigid/flex
  PDBQT が `vina.Vina().set_receptor(rigid, flex)` に問題なくロードでき、
  グリッドマップも計算できることを確認済み。
- 各メンバーの共結晶リガンドをそれぞれ自己ドッキングし、承認薬など無関係
  な分子と affinity で識別できることを確認済み（`data/ligands.smi`）。

## モジュール構成 (`dockinglib/`)

| module | 役割 |
|---|---|
| `receptor_prep.py` | PDB取得・チェイン分離・残基整形（TER挿入・CYX化）・PDBFixer による構造修正 |
| `pocket.py` | ドッキングボックス算出、共結晶リガンド距離によるフレキシブル残基検出・Meeko flexres 文字列整形 |
| `ensemble.py` | 複数コンフォメーションを一括で `receptor_prep` + `pocket` + Meeko PDBQT化し、`manifest.json` として保存/読込 |
| `ligand_prep.py` | `.smi` 読み込み、SMILES→3D(ETKDGv3+MMFF)→Meeko 配位子PDBQT化 |
| `docking.py` | Vina のフレキシブル受容体対応薄ラッパー (`make_vina` / `dock_ligand`) |
| `screening.py` | 全リガンド×全メンバーの並列ドッキング、best-of-ensembleランキング、結果CSV/SDF出力 |
| `refine_md.py` | 上位ヒットのみ短時間陰溶媒MDで緩和・RMSD安定性評価・再ランキング |
| `io.py` | PDBQT→RDKit変換、ポーズSDF書き出し、結果CSV入出力 |
| `parallel.py` | 独立タスクを `ProcessPoolExecutor` で並列化する `parallel_map` |
| `progress.py` | 完了ごとに1行進捗を表示する `DockProgress` / `RefineProgress` |
| `pipeline.py` | 準備→ドッキング→MD後処理を一気通貫で呼べる高水準関数 |
| `cli.py` | `dockinglib-prep` / `dockinglib-dock` / `dockinglib-refine` コマンド |

## 制限・拡張の余地

- フレキシブル残基は共結晶リガンドからの距離のみで決めている。既知の
  誘導適合が起きやすい特定残基を明示指定したい場合は `pocket.py` の
  `find_pocket_residues` を使わず、`ensemble.prepare_ensemble_member` に
  直接 flexres 文字列を渡す形に拡張しやすい。
- MD後処理は既定で陰溶媒(GBn2)・CPUのみを想定した短時間シミュレーション
  であり、明示溶媒・GPUでの本格的な自由エネルギー計算などは対象外
  （必要になれば `refine_md.py` を独立に拡張できる）。
- `screening.py` は `(リガンド, メンバー)` の全組み合わせをタスクとして
  並列化するため、Vinaのグリッドマップはタスクごとに再計算される。
  大規模ライブラリでは、メンバーごとにグリッドマップを1回だけ計算して
  使い回す方式（メンバー単位の並列化）への変更が有効な最適化になる。
