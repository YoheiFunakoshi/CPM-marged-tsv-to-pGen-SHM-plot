# CPM merged AIRR TSV to pGen-SHM plot

CPM社のUMI集約後FASTAをIgBLAST `-outfmt 19`で解析した
`.umi_exact.igblast.airr.tsv`から、pGen、SHM、pGen-SHM図を作成するWindows GUIです。

実装の基本構造と図の考え方は
[RG版](https://github.com/YoheiFunakoshi/RG-marged-tsv-to-pGen-SHM-plot)
を踏襲しています。ただし、RGのreadとCPMのexact UMI familyを同じ単位として扱わないよう、
CPM固有の入力検証・重み付け・ファイル名・QC表示を追加しています。

## 最重要：解析単位

- 主解析単位は **exact UMI family** です。CPM AIRR TSVの1行を1 familyとして数えます。
- `sequence_id`の`|UMI=XXXXXXXXXXXX|DUPCOUNT=N`を必須とします。
- `DUPCOUNT`はそのfamilyを支持する **uncollapsed merged read数** です。
- exact UMI family数とDUPCOUNT supporting read数は別単位です。合算・直接比較しません。
- GUIの既定はfamily単位です。DUPCOUNT重み付き図は明示的にチェックした場合だけ副解析として作成します。
- UMIが12塩基のA/C/G/Tでない、DUPCOUNTが正整数でない、または注釈が欠ける入力は停止します。黙って1として補完しません。

このアプリはExcelを入力にしません。入力は元のIgBLAST AIRR TSVです。Excelは解析後の行データ出力の1つです。

## pGenの定義

- 入力AIRR TSVにはpGen列はありません。
- productive/canonical/IGHフィルター後の`junction_aa`ごとに、OLGA 1.3系の
  `human_B_heavy`モデルで`compute_aa_CDR3_pgen(junction_aa)`を計算します。
- 同一`junction_aa`のpGenはキャッシュできます。
- `pGen = 0`は行データには残し、対数軸の図からだけ除外します。
- 図のx軸は`log10(pGen)`です。

pGenは再構成配列の生成確率であり、観測family数やDUPCOUNTとは別概念です。

## SHMの定義

SHMはAIRR列`v_identity`から計算します。

- `v_identity`が0〜1スケールの場合：`SHM(%) = (1 - v_identity) × 100`
- `v_identity`が0〜100スケールの場合：`SHM(%) = 100 - v_identity`
- 負値は0に丸めます。

これはIgBLASTが報告したVアラインメント範囲に対するidentity由来のSHM proxyです。
V領域全長の厳密な塩基置換数ではなく、短いVアラインメントの影響を受けます。
必要な場合はGUIの`Min V alignment length`で最低長を設定します。

## 採用条件

既定はRG版と同じ探索的フィルターです。GUIの`Junction AA filter`で次を選べます。

- `RG reference`（既定）：空欄、`*`、`X`、非アルファベットを除外。先頭C、末尾F/W、長さ5〜40は必須にしません。
- `CPM conservative`：上記に加えて、長さ5〜40、先頭C、末尾F/Wを必須にします。CPM集計Excelの`Productive canonical`定義との整合確認用です。

モード名と除外理由はQC・run log・run conditions JSONへ記録します。

1. `locus`が存在し、値が入っている場合は`IGH`のみ採用
2. `locus`空欄は採用し、件数をQCへ記録
3. `productive == T`
4. `vj_in_frame`列があれば`T`
5. `stop_codon`列があれば`F`
6. `junction`はA/C/G/Tのみ
7. `junction_aa`は選択したfilter modeで判定
8. `v_identity`欠損を除外
9. 指定時はVアラインメント長が最低値以上

`productive`、`vj_in_frame`、`stop_codon`の判定は大文字小文字を無視し、
`T/TRUE/1/Y/YES`をtrueとして扱います。

## GUIの使い方

1. 初回だけ`setup_env.bat`をダブルクリックします。
2. `run_gui.bat`をダブルクリックします。
3. `AIRR TSV`にCPMの`.umi_exact.igblast.airr.tsv`を選びます。
4. 出力フォルダーとsample名を確認します。
5. 必要なら`Also create DUPCOUNT supporting-read-weighted plots`を選びます。
6. `Check setup`でOLGAモデルと依存関係を確認します。
7. `Run pGen + SHM + pGen-SHM plot`を押します。

処理はバックグラウンドスレッドで実行し、進捗・エラー原因をGUIログへ表示します。
元TSVや既存Excelは上書きしません。

## 主な出力

常に作成する主解析出力：

- `*_qc_summary.tsv`：除外理由、family数、supporting read数など
- `*_pgen_bins.tsv`：unique AA、exact UMI family、DUPCOUNT supportを別列で保存
- `*_pgen_bins_unique_junction_aa.png`
- `*_pgen_bins_exact_umi_family.png`
- `*_shm_hist_exact_umi_family.tsv/.png`
- `*_pgen_shm_rows.xlsx`：採用された1 exact UMI family/行
- `*_pgen_shm_points.tsv`
- `*_pgen_shm_beta1_unique_junction_points.tsv`
- `*_pgen_shm_roi_summary.tsv`
- `*_pgen_shm_scatter_exact_umi_family.png`
- `*_pgen_shm_kde_exact_umi_family.png`
- `*_pgen_shm_kde_exact_umi_family_log_density.png`
- `*_pgen_shm_kde_beta1_unique_junction_unweighted.png`
- `*_pgen_shm_kde_beta1_exact_umi_family.png`
- `*_run_conditions.json`
- `*_run_log.txt`
- `pgen_cache.tsv`

DUPCOUNT副解析を選択した場合だけ、ファイル名に
`dupcount_supporting_reads`を含むヒストグラム・散布図・KDEを追加します。
無単位の`weighted`という名前は使いません。

## 行データの重要列

`*_pgen_shm_rows.xlsx`と`*_pgen_shm_points.tsv`には次を保存します。

- `sequence_id`, `umi`
- `umi_family_count`（常に1）
- `supporting_read_count`（DUPCOUNT）
- `junction`, `junction_aa`
- `v_identity`, `shm`
- `pgen`, `log10_pgen`
- `locus`, `productive`, `v_call`, `j_call`
- 同じx-y座標に重なるfamily数とsupporting read数

## 最終産物の見方

CPM版では、最終産物を **主解析** と **補助解析** に分けて見ます。
最初に確認する中心データは `*_pgen_shm_rows.xlsx` です。これはフィルター後に採用されたAIRR各行を、1 exact UMI familyとして残したrow-levelデータです。図だけではなく、このExcel/TSVを解析本体として扱います。

### 主解析として見るもの

- `*_pgen_shm_rows.xlsx`：最終データ本体。1行 = 1 exact UMI familyです。
- `*_pgen_shm_points.tsv`：同じrow-levelデータをTSVで保存したものです。
- `*_pgen_shm_kde_exact_umi_family.png`：CPM版の主図候補です。exact UMI family数を単位としてpGen-SHM分布をKDE表示します。
- `*_pgen_shm_scatter_exact_umi_family.png`：KDEで見えにくい点の存在確認に使います。
- `*_qc_summary.tsv` と `*_run_conditions.json`：採用条件、除外理由、解析単位、DUPCOUNT副解析の有無を確認します。

`*_pgen_shm_kde_exact_umi_family.png`では、x軸が `log10(pGen)`、y軸がSHMです。
1点の元になるのは、フィルター後のAIRR 1行、すなわち1 exact UMI familyです。
同じx-y座標に複数familyが重なる場合は、その座標の密度にfamily数が反映されます。
したがって、CPM版の主図は「read数」ではなく「UMI family数」を反映した図です。

### 補助的に見るもの

- `*_pgen_shm_kde_exact_umi_family_log_density.png`：低密度集団を見落とさないための確認図です。通常KDEで薄く見える集団を確認する目的で使います。
- `*_pgen_bins.tsv`：pGen binごとに、unique junction AA、exact UMI family、DUPCOUNT supporting readを別列で保存します。
- `*_shm_hist_exact_umi_family.tsv/.png`：SHM分布をexact UMI family単位で確認します。
- `*_pgen_shm_roi_summary.tsv`：低SHM・高pGenなどの領域を図の色ではなく数値で確認します。

### beta1互換出力

`beta1`を含むファイルは、旧解析・前任者法との比較用です。現在の主解析ではありません。

- `*_pgen_shm_kde_beta1_unique_junction_unweighted.png`：unique junctionを1点として扱う旧法比較用の図です。
- `*_pgen_shm_kde_beta1_exact_umi_family.png`：beta1形式の座標にexact UMI family数を反映した比較図です。

beta1互換出力は、過去図との見え方を比較するために残しています。最終判断では、まずrow-levelの `*_pgen_shm_rows.xlsx` と `*_pgen_shm_kde_exact_umi_family.png` を見ます。

### DUPCOUNT重み付き副解析

DUPCOUNTは、1 exact UMI familyを支えていた元merged read数です。直感的には「UMIでまとめる前のread数の名残」です。
GUIで `Also create DUPCOUNT supporting-read-weighted plots` を選んだ場合だけ、`dupcount_supporting_reads` を含むファイルが追加されます。

DUPCOUNT重み付き図では、点の単位はUMI familyのままですが、KDEやヒストグラムの重みとして `supporting_read_count`、つまりDUPCOUNTを使います。
これはUMIなしデータのread数重み付きに近い見方ですが、PCR増幅やシーケンス深度の影響を受けやすいため、CPM版では副解析として扱います。

解釈の目安は次の通りです。

- exact UMI family主解析でもDUPCOUNT副解析でも同じ集団が強い：family数としてもread支持としても強い集団です。
- exact UMI family主解析では弱いがDUPCOUNT副解析で強い：少数familyが多く読まれた可能性があり、増幅・samplingの影響を考慮します。
- DUPCOUNT副解析だけを主結論にしない：UMIがあるCPMでは、主解析単位はexact UMI familyです。

短くまとめると、CPM版では **最終データ本体は `*_pgen_shm_rows.xlsx`、主図は `*_pgen_shm_kde_exact_umi_family.png`、DUPCOUNT重み付き図はread数寄りの補助確認** です。

## CLI

```powershell
.\.venv\Scripts\python.exe .\cpm_airr_pgen_shm_plot.py `
  --input "sample.umi_exact.igblast.airr.tsv" `
  --outdir "result" `
  --sample "SAMPLE" `
  --pgen-workers 6
```

DUPCOUNT副解析も作る場合：

```powershell
.\.venv\Scripts\python.exe .\cpm_airr_pgen_shm_plot.py `
  --input "sample.umi_exact.igblast.airr.tsv" `
  --outdir "result" `
  --include-supporting-read-outputs
```

CPM集計Excelのconservative canonical条件に合わせる場合は
`--canonical-mode cpm_conservative`を追加します。既定は`rg_reference`です。

## 依存関係

- Python 3.11推奨
- numpy
- matplotlib
- scipy
- olga
- openpyxl

## RGとの比較上の注意

- RG primary count：merged read
- CPM primary count：exact UMI family
- CPM secondary count：DUPCOUNT supporting merged read

同じ元PBMC由来でも別アリコート、別ライブラリー、別シーケンスランです。
RG read、CPM family、CPM DUPCOUNT supportを同一分母にした差・比・相関は自動作成しません。
pGen-SHM分布を比較する場合は、単位を揃えた並列図として解釈してください。

## テスト

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

テストはUMI/DUPCOUNT解析、異常入力拒否、productive/canonical/IGHフィルター、
family/supporting read保存則、SHM計算を確認します。
