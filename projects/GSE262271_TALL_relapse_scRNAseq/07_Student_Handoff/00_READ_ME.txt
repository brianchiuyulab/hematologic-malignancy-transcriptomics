血癌 transcriptomics 交接包｜2026-09-15 核對版

先看 01_Summary 的 PPT，再看 02_Figures。03_Tables 是可直接用 Excel 開啟的結果。
不要再使用舊的 8-slide handoff、舊 Fig07 或舊 Venn 數字代表目前分析。

一、為什麼做這個分析？
我們想知道 T-ALL relapse 的 transcriptome，是否與 STAT5B 過度活化模型呈現共同的 pathway 變化，再據此評估 proteasome / mitophagy 機制假說。
閱讀順序是：兩組 dataset 的完整 KEGG GSEA → 同方向交集 → 再回到既有的 target 假說。
Mitophagy 不是唯一、也不是最強的共同 hit。選擇深入查看它，來自既有的生物學問題。

二、兩組 dataset
GSE262271：human T-ALL scRNA-seq，3 位病人的 diagnosis / relapse 配對。
GSE218858：mouse Rag2-/- DN bulk RNA-seq，STAT5B N642H 與 control 各 3 個獨立樣本。
Single-cell 以每個 patient sample 聚合原始 counts，用 paired edgeR；bulk 用 DESeq2 獨立組比較。
正 NES 分別表示 relapse 較高、STAT5B N642H 較高。

三、目前完整 KEGG 結果（FDR < 0.10）
Single-cell exploratory within-sample top 15%：83 條上調、5 條下調。
Bulk STAT5B N642H：128 條上調、8 條下調。
兩邊都測到的 275 條 pathway 中，共同上調 54 條，共同下調 0 條。
原始 global top 25% 保留作對照：47 條上調、1 條下調，與 bulk 共同上調 23 條。
每組的總數使用各自完整 catalog；Venn 僅限兩組都測到的 pathway，因此圈內加總不等於上述全部總數。
完整名單在 Table04、Table05；同方向交集在 Table07。圖中的 disease 名稱不代表病人真的患有該疾病，要看重疊的基因功能。

四、Proteasome / Mitophagy：方向沒有反轉，都是上調
                             Proteasome FDR     Mitophagy FDR
scRNA-seq global top 25%          0.126157           0.126157
scRNA-seq within-sample top 15%   0.008379           0.013197
STAT5B N642H bulk RNA-seq         0.001078           0.028234
這些都是完整 KEGG collection 的 BH FDR，不是只校正兩個 target。
General autophagy 未達 FDR < 0.10：top 15% FDR 0.9831，bulk FDR 0.4308。
Top 15% 的兩個 target leading-edge 平均表達都是 3/3 病人上升。但基因是由同一份 GSEA 選出，這是描述性方向檢查，不是額外獨立驗證。
Heatmap 的 gene-level FDR 與 pathway-level FDR 分開標示，沒有把 pathway 星號貼成每一個基因都顯著。

五、舊圖為何看起來矛盾？
1. 舊 cross-dataset 圖用的是 global top 25%，後來「兩條都顯著」用的是 within-sample top 15%。舊圖缺少清楚的 population 定義。
2. 舊 GSEA 接近門檻的數字有 Monte Carlo 精度問題。本次固定 seed、提高抽樣精度，重新跑全部 pathways，沒有選最小 P 的 seed。
3. 已用原來的 edgeR 4.4.2 成功重現全部 10,950 個 single-cell gene-level 結果；bulk 17,991 genes 的 DESeq2 結果也重現。電腦的另一套 R library 有不同版本，現已鎖定版本，避免交接時載錯。
4. 將原本儲存的 ranking 直接做高精度 GSEA，與本次從 counts 重建 ranking 的結果一致。因此這次 top 25% FDR 變動不是偷偷改了 disease contrast 或移除某位病人。
5. 舊 top-15% sensitivity 限制 pathway size ≤ 500，實際測試 352 條，但欄位寫成 353。本次與原始完整 KEGG 分析統一為至少 5 個可測基因、不設上限，實測 353 條並直接核對 BH。Top 15% 新舊數字差異同時包含抽樣精度及這個測試範圍差異。

六、15% 與 annotation 的必要界線
867 顆 non-T reference 的 empirical q95 = 0.0269455。1,978 / 12,828（15.419%）T-lineage candidates 高於這個值。
「固定 CNV score > 0.0269455」與「每個 sample 各取最高 15%」不同。後者保留 1,926 cells，且不同 sample 的 cutoff 不一樣。
因此 15.4% 的相近比例，不能讓 top 15% 變成 unbiased，也不能證明 95% specificity。部分 sample 的 top-15% cutoff 低於 pooled reference q95。
Top 15% 是看過 sensitivity 結果後採用的 exploratory core。完整 KEGG FDR 沒有額外校正挑選分析條件的影響。
新版 annotation 使用 marker 支持的 T-lineage / B / myeloid / pDC / erythroid / cytotoxic T-NK 廣義分類。HSC-like、quiescent、ER-stress 舊稱呼不再當作已驗證 T-ALL subtype。
CNV-high 是另外的操作性篩選，不等於已確診 malignant。投稿前仍需要獨立 reference mapping、同樣本 doublet 檢查、mutation / clonotype / DNA-CNV 或其他可用的正交證據。

七、怎麼解讀生物意義？
共同 pathway enrichment 支持相似的 transcriptional programs，不證明 STAT5 進入粒線體、mitophagy flux 增加、proteasome 依賴性或 bortezomib resistance。
兩組資料都不是 bortezomib responder vs resistant 的直接比較。

八、檔案
01_Summary：14 張 PPT，含中文交接說明。
02_Figures：PNG 供閱讀，PDF 供排版。主圖先看 Fig03、04、05、06、07。
03_Tables：有 00_TABLE_GUIDE.txt 說明用途。
04_Code_Availability：R scripts、原始 bulk pipeline、版本紀錄、核對表與重現說明。
05_Data_Locations：GEO accession、Raw data / R objects 的實際路徑。
Raw archive 和大型 R objects 沒有刪除，也不上傳 GitHub。

GitHub：https://github.com/brianchiuyulab/hematologic-malignancy-transcriptomics
本報告是可追溯的探索型分析交接，不聲稱已完成所有 reviewer 要求的驗證。
