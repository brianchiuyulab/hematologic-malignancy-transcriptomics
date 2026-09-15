血癌 transcriptomics 交接包｜主分析統一版 2026-09-15

現在主分析固定：within-sample CNV top 15% + 完整 KEGG BH FDR < 0.05。
先看 01_Summary 的 14 頁 PPT；圖在 02_Figures，CSV 在 03_Tables。

一、研究目的與流程
比較 T-ALL relapse 與 STAT5B 過度活化是否有同方向的 transcriptomic programs。
依序呈現 single-cell 完整 KEGG、bulk 完整 KEGG、同方向交集，最後檢視既有
proteasome / mitophagy 假說。Mitophagy 不是唯一或最強的共同 pathway。

二、固定主分析設定
Single-cell：GSE262271，3 位病人的 diagnosis / relapse 配對。
沿用已核對的原始 QC。在每個 sample 的 T-lineage candidates 中，保留 inferCNV
score 最高 15%，共 1,926 cells，最小 sample 107 cells。
原始 counts 以 sample 聚合，edgeR 模型 ~ patient_pair + condition_simple。
正 NES 表示 relapse 較高。GSEA 使用全部可排名基因，不另設 DEG FC 門檻。

Bulk：GSE218858，mouse Rag2-/- DN，STAT5B N642H 與 control 各 3 個獨立樣本。
DESeq2 模型 ~ condition，正 NES 表示 N642H 較高。

兩組均採各自完整 KEGG collection 的 BH FDR < 0.05。
交集要求：兩邊都測到、兩邊 FDR 都 < 0.05、NES 同號。
此輪只統一報告門檻，沒有改 QC、gene ranking、原始 P value 或 FDR 計算。

三、主分析結果
Single-cell：61 條上調、2 條下調，完整測試 353 條。
Bulk：107 條上調、6 條下調，完整測試 300 條。
兩邊共同測試的 275 條中：40 條共同上調、0 條共同下調。
Venn 的圈內數量限於共同測試範圍，因此不等於各自完整 catalog 的總數。

Proteasome：single-cell FDR 0.008379；bulk FDR 0.001078，均上調。
Mitophagy：single-cell FDR 0.013197；bulk FDR 0.028234，均上調。
General autophagy 不顯著：single-cell FDR 0.9831；bulk FDR 0.4308。
上述 FDR 已校正各自完整 KEGG collection，不是只校正兩個 target。

Fig03、04：各自顯著 pathway 概覽。Table04、05 保留全部測試結果。
Table11：只列 FDR < 0.05 的主分析顯著結果，適合先讀。
Fig05 / Table07：40 條同方向交集。
Fig06：全部交集熱圖。Fig06B：按兩邊較大的 FDR 排序之較強共同證據。
Fig07：target 結果，含不顯著的 autophagy，避免只呈現正結果。

四、為什麼用 top 15%，q95 如何支持？
本報告以 sample 內相同百分位規則，分析相對 CNV score 較高的 T-lineage core。
867 顆 non-T reference 的 q95 = 0.0269455。全體 12,828 顆 T-lineage candidates
中，1,978 顆（15.419%）超過此固定分數，可作為高 CNV 細胞比例的描述性對照。

但固定分數 > 0.0269455 和每個 sample 最高 15% 不是同一個選法：
top 15% 的 1,926 顆中，有 1,407 顆（73.1%）高於 reference q95。
因此 q95 可以作支持性資料，不能說所有 top15 cells 都超過正常上限，也不能
單憑這個比例唯一推導出 15%，或聲稱 top15 是與結果無關、事前指定的門檻。
Top15 是 sensitivity 檢查後固定的主分析設定。這項分析歷史在 Methods 揭露；
完整 KEGG FDR < 0.05 並沒有額外校正選擇分析條件的影響。
每個 sample 的實際 cutoff、細胞數與重疊比例：TableS03、TableS06；分布：FigS02。

五、補充分析與已核對事項
原始 global top25、within-sample top25、固定 reference q95 留作 sensitivity。
Global top25 在 FDR < 0.05 與 bulk 共同上調 13 條；兩個 target 的 FDR 均 0.1262。
補充數字不混入主分析 Fig05 或 Table07。
所有 10,950 個原始 single-cell、17,991 個 bulk gene-level 統計已重現到浮點精度。
本次採固定 seed 的高精度 GSEA；R 與 package 版本已記錄。
先前 FDR < 0.10 的 54 條交集不再作為目前主分析數字。

六、Annotation 與生物意義
Fig01A / 01B 提供廣義 lineage 與 cluster marker 證據。
CNV-high potential malignant T cells 是操作性定義，尚未完成獨立 malignant 身分驗證。
HTO singlet 過濾不能完全排除同一個 sample 內的 doublets。
Leading-edge 平均表達在兩個 target 皆為 3/3 病人上調，是描述性方向檢查。
Gene-level FDR 與 pathway FDR 分開標註，不能把 pathway 顯著說成每個基因顯著。
共同富集不證明 STAT5 的粒線體定位、mitophagy flux、proteasome 依賴性或 BTZ 抗藥性。

七、資料位置
01_Summary：目前唯一 PPT。
02_Figures：PNG 閱讀、PDF 排版；FigS 為補充。
03_Tables：可用 Excel 開啟，先看 Table11 與 Table07。
04_Code_Availability：程式、方法說明、版本與核對記錄。
05_Data_Locations：Raw data / 原始和更新版 R objects 的本機位置。
Raw data 與大型 R objects 保留在原始專案，不上傳 GitHub。
GitHub：https://github.com/brianchiuyulab/hematologic-malignancy-transcriptomics
