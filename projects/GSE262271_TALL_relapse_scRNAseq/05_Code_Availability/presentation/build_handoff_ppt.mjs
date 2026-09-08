import fs from "node:fs/promises";
import path from "node:path";
import { pathToFileURL } from "node:url";
import { Presentation, PresentationFile } from "@oai/artifact-tool";

const { SKILL_DIR, TMP_DIR, PROJECT_ROOT, RUNTIME_PYTHON } = process.env;
if (![SKILL_DIR, TMP_DIR, PROJECT_ROOT, RUNTIME_PYTHON].every((p) => path.isAbsolute(p ?? ""))) {
  throw new Error("SKILL_DIR, TMP_DIR, PROJECT_ROOT and RUNTIME_PYTHON must be absolute paths");
}

const { finalizePresentation } = await import(
  pathToFileURL(path.join(SKILL_DIR, "container_tools/artifact_tool_utils.mjs")).href
);

await fs.mkdir(TMP_DIR, { recursive: true });
const FINAL_PPTX = path.join(
  TMP_DIR,
  "GSE262271_TALL_relapse_scRNAseq_HANDOFF_PUBLICATION_FINAL.pptx",
);
const FIG = path.join(PROJECT_ROOT, "01_Figures");

const W = 1280;
const H = 720;
const FONT = "Microsoft JhengHei";
const NAVY = "#14213D";
const BLUE = "#0072B2";
const ORANGE = "#D55E00";
const GOLD = "#E69F00";
const GREEN = "#009E73";
const RED = "#B2182B";
const INK = "#20242A";
const MUTED = "#59636E";
const LIGHT = "#EEF1F4";
const PALE = "#F7F7F4";
const WHITE = "#FFFFFF";

const presentation = Presentation.create({ slideSize: { width: W, height: H } });

function shape(slide, left, top, width, height, fill = "none", lineFill = "none", radius = "rect") {
  return slide.shapes.add({
    geometry: radius,
    position: { left, top, width, height },
    fill,
    line: { fill: lineFill, width: lineFill === "none" ? 0 : 1 },
  });
}

function textBox(slide, text, left, top, width, height, opts = {}) {
  const box = shape(slide, left, top, width, height, opts.fill ?? "none", opts.line ?? "none", opts.geometry ?? "rect");
  box.text = text;
  box.text.style = {
    typeface: FONT,
    fontSize: opts.fontSize ?? 24,
    bold: opts.bold ?? false,
    color: opts.color ?? INK,
    autoFit: opts.autoFit ?? "shrinkText",
  };
  return box;
}

function addBase(slide, title, section, page) {
  slide.background.fill = WHITE;
  shape(slide, 0, 0, W, 8, ORANGE);
  textBox(slide, section.toUpperCase(), 48, 22, 280, 24, { fontSize: 13, bold: true, color: ORANGE });
  textBox(slide, title, 48, 45, 1180, 48, { fontSize: 30, bold: true, color: NAVY });
  shape(slide, 48, 686, 1184, 1.2, "#D7DCE1");
  textBox(slide, "GSE262271 paired T-ALL relapse analysis", 48, 692, 520, 17, { fontSize: 10, color: MUTED });
  textBox(slide, String(page).padStart(2, "0"), 1185, 692, 46, 17, { fontSize: 10, bold: true, color: MUTED });
}

async function addImage(slide, imagePath, position, alt) {
  const bytes = await fs.readFile(imagePath);
  return slide.images.add({
    blob: bytes,
    contentType: "image/png",
    alt,
    fit: "contain",
    position,
  });
}

function callout(slide, title, body, left, top, width, height, accent = BLUE) {
  shape(slide, left, top, width, height, PALE, "#D8DEE5", "roundRect");
  shape(slide, left, top, 7, height, accent, accent, "roundRect");
  textBox(slide, title, left + 20, top + 14, width - 34, 26, { fontSize: 17, bold: true, color: NAVY });
  textBox(slide, body, left + 20, top + 47, width - 34, height - 60, { fontSize: 16, color: INK });
}

function setNotes(slide, text) {
  slide.speakerNotes.textFrame.setText(text);
}

// Slide 1
{
  const s = presentation.slides.add();
  s.background.fill = NAVY;
  shape(s, 0, 0, 18, H, ORANGE);
  textBox(s, "HEMATOLOGIC MALIGNANCY PROJECT", 72, 68, 780, 28, { fontSize: 15, bold: true, color: GOLD });
  textBox(s, "T-ALL relapse single-cell transcriptomics", 72, 125, 1000, 115, { fontSize: 48, bold: true, color: WHITE });
  textBox(s, "GSE262271 paired diagnosis and relapse analysis", 74, 250, 900, 45, { fontSize: 26, color: "#D9E1EA" });
  shape(s, 74, 332, 1080, 2.5, ORANGE);
  textBox(s, "研究問題", 74, 365, 180, 32, { fontSize: 19, bold: true, color: GOLD });
  textBox(s, "Relapse 的 potential malignant T cells 是否出現與 STAT5、粒線體品質控制及 proteostasis 一致的 pathway-level program？", 74, 405, 1080, 92, { fontSize: 26, bold: true, color: WHITE });
  textBox(s, "Publication figure audit, reproducible code and lab handoff", 74, 595, 780, 32, { fontSize: 17, color: "#BFC9D5" });
  textBox(s, "September 2026", 980, 595, 176, 32, { fontSize: 17, bold: true, color: GOLD });
  setNotes(s, "這份簡報是血癌 project 的 GSE262271 單細胞模組交接版。所有結果以 association 解讀，不宣稱 mitochondrial STAT5 的直接因果作用。");
}

// Slide 2
{
  const s = presentation.slides.add();
  addBase(s, "為什麼做這個分析", "Biological rationale", 2);
  textBox(s, "工作假說", 52, 116, 200, 30, { fontSize: 18, bold: true, color: ORANGE });
  const items = [
    ["1", "STAT5 activity", "構成性訊號可能改變 leukemia cell state"],
    ["2", "Mitochondrial program", "OXPHOS、stress response 與 organelle quality control"],
    ["3", "Mitophagy and proteostasis", "清除受損粒線體，並改變 proteasome demand"],
    ["4", "Relapse biology", "可能形成適應、依賴或治療抗性表型"],
  ];
  items.forEach((d, i) => {
    const x = 54 + i * 300;
    shape(s, x + 18, 176, 230, 4, i === 3 ? ORANGE : BLUE);
    shape(s, x, 202, 52, 52, i === 3 ? ORANGE : NAVY, "none", "ellipse");
    textBox(s, d[0], x + 13, 211, 28, 30, { fontSize: 19, bold: true, color: WHITE });
    textBox(s, d[1], x + 64, 200, 210, 34, { fontSize: 19, bold: true, color: NAVY });
    textBox(s, d[2], x + 64, 240, 210, 78, { fontSize: 15, color: INK });
  });
  callout(
    s,
    "這組資料能回答什麼",
    "Paired transcriptomic association、pathway direction、跨資料集 concordance。",
    70, 390, 520, 155, BLUE,
  );
  callout(
    s,
    "這組資料不能單獨證明什麼",
    "STAT5 的粒線體定位、mitophagy flux、proteasome dependency、bortezomib resistance 的直接因果。",
    635, 390, 575, 155, ORANGE,
  );
  textBox(s, "核心原則：先把 pathway-level evidence 與 gene-level evidence 分開，再談機制。", 82, 585, 1100, 40, { fontSize: 21, bold: true, color: NAVY });
  setNotes(s, "用四個層次交代研究動機。最重要的是不要把 GSEA 的 pathway-level 顯著直接寫成單一基因顯著或直接因果。");
}

// Slide 3
{
  const s = presentation.slides.add();
  addBase(s, "配對設計是整個統計分析的關鍵", "Study design", 3);
  const table = s.tables.add({
    rows: 4,
    columns: 4,
    left: 54,
    top: 132,
    width: 690,
    height: 260,
    columnWidths: [120, 175, 175, 220],
    values: [
      ["Patient", "Diagnosis", "Relapse", "Dx / Rel primary cells"],
      ["P1", "M104", "M127", "140 / 971"],
      ["P2", "M143", "M148", "712 / 269"],
      ["P3", "M187", "M187r", "456 / 659"],
    ],
  });
  table.styleOptions = { headerRow: true, bandedRows: true };
  table.borders.assign({ style: "solid", fill: "#CBD2D9", width: 1 });
  table.cells.block({ row: 0, column: 0, rowCount: 1, columnCount: 4 }).assign({
    fill: NAVY,
    textStyle: { typeface: FONT, fontSize: 16, bold: true, color: WHITE },
  });
  table.cells.block({ row: 1, column: 0, rowCount: 3, columnCount: 4 }).assign({
    textStyle: { typeface: FONT, fontSize: 16, color: INK },
  });
  callout(s, "Biological replicate", "病人才是 biological replicate，cell 不是獨立 replicate。", 790, 132, 420, 110, BLUE);
  callout(s, "Paired pseudobulk model", "每個 sample 先聚合 counts，再以 edgeR design ~ patient_pair + condition 比較 relapse 與 diagnosis。", 790, 263, 420, 150, ORANGE);
  callout(s, "GSEA ranking", "Genome-wide rank = sign(logFC) × sqrt(F)。配對資訊已經在上游 model 中，不是 GSEA 演算法自己判斷 pair。", 790, 434, 420, 150, GREEN);
  textBox(s, "N = 3 pairs，統計力有限，因此跨病人方向一致性與獨立資料集驗證格外重要。", 55, 445, 690, 90, { fontSize: 22, bold: true, color: RED });
  setNotes(s, "這張要特別講清楚 paired design。GSEA 使用的 rank 來自 paired edgeR model，pair 並不是 fgsea 自動辨識的。");
}

// Slide 4
{
  const s = presentation.slides.add();
  addBase(s, "Primary population 只保留 CNV-high top 25%", "Cell definition", 4);
  await addImage(s, path.join(FIG, "01_Main_Figures", "Fig02_primary_population_UMAP.png"), { left: 48, top: 104, width: 780, height: 560 }, "UMAP showing the primary potential malignant T-cell population");
  callout(s, "Operational definition", "Candidate malignant T-lineage cells 中，inferCNV burden 最高的 top 25%。", 850, 135, 360, 126, ORANGE);
  callout(s, "Primary analysis", "後續 paired pseudobulk DEG、KEGG、GO BP、STAT5 與 hypothesis-focused analyses 都使用這 3,207 顆細胞。", 850, 286, 360, 165, BLUE);
  callout(s, "Interpretation boundary", "這是 enriched malignant-like population，不等同臨床上已驗證的 malignant cell call。", 850, 478, 360, 130, RED);
  setNotes(s, "主分析只使用 top 25%。Top 50%、全 candidate T-lineage、cycling 和 quiescent subsets 只出現在 sensitivity figures。");
}

// Slide 5
{
  const s = presentation.slides.add();
  addBase(s, "Relapse 並沒有一致增加 potential malignant T-cell proportion", "Cell abundance", 5);
  await addImage(s, path.join(FIG, "01_Main_Figures", "Fig01_cohort_cell_counts_and_paired_proportions.png"), { left: 38, top: 102, width: 1204, height: 548 }, "Candidate T-cell counts and paired potential malignant T-cell proportions");
  textBox(s, "結論：P1 上升，P2 與 P3 下降。Exact paired Wilcoxon P = 0.75，因此不能宣稱 relapse 的 malignant-like fraction 普遍變多。", 80, 628, 1120, 42, { fontSize: 18, bold: true, color: NAVY, fill: "#FFF3EB" });
  setNotes(s, "先回答 abundance 問題。P1 的增加很大，但另外兩位病人下降，因此 group-level paired test 不顯著。");
}

// Slide 6
{
  const s = presentation.slides.add();
  addBase(s, "inferCNV 提供的是 relative CNV burden，不是最終腫瘤診斷", "Cell definition QC", 6);
  await addImage(s, path.join(FIG, "01_Main_Figures", "Fig03_inferCNV_nonT_reference_heatmap.png"), { left: 45, top: 108, width: 900, height: 545 }, "inferCNV heatmap using non-T reference cells");
  callout(s, "Reference", "B cell、myeloid、erythroid、pDC 等 non-T populations。", 960, 145, 270, 118, BLUE);
  callout(s, "Observation", "Candidate T-lineage cells 呈現相對較強的 chromosome-scale expression deviations。", 960, 285, 270, 142, ORANGE);
  callout(s, "Selection", "以 cell-level burden score 的 75th percentile 定義 top 25% primary population。", 960, 450, 270, 136, GREEN);
  setNotes(s, "原始 inferCNV heatmap 是軟體輸出，保留作方法 QC。文章主敘事仍應以 operational burden definition 描述，不要說成確定的 DNA CNV。");
}

// Slide 7
{
  const s = presentation.slides.add();
  addBase(s, "Relapse 的主訊號是 biosynthesis、replication 與 oxidative program", "KEGG GSEA", 7);
  await addImage(s, path.join(FIG, "01_Main_Figures", "Fig04_KEGG_GSEA_key_pathways.png"), { left: 52, top: 104, width: 880, height: 565 }, "Focused KEGG GSEA dot plot");
  callout(s, "FDR < 0.05", "Ribosome biogenesis、OXPHOS、DNA replication、ribosome、spliceosome。", 948, 126, 282, 155, ORANGE);
  callout(s, "Exploratory FDR < 0.10", "Alanine, aspartate and glutamate metabolism、JAK-STAT、HIF-1、mitophagy。", 948, 302, 282, 175, BLUE);
  callout(s, "Not significant", "Proteasome q = 0.109，sphingolipid signaling q = 0.359，general autophagy q = 0.754。", 948, 499, 282, 145, RED);
  setNotes(s, "最重要的精確數字：KEGG mitophagy NES 1.44, FDR 0.0998。Proteasome NES 1.53, FDR 0.109，不能寫成過 0.10。");
}

// Slide 8
{
  const s = presentation.slides.add();
  addBase(s, "GO BP 支持 translation 與 mitochondrial respiration，但不支持 broad autophagy", "GO BP GSEA", 8);
  await addImage(s, path.join(FIG, "01_Main_Figures", "Fig05_GO_BP_GSEA_curated_programs.png"), { left: 50, top: 103, width: 930, height: 567 }, "Curated nonredundant GO biological-process GSEA");
  callout(s, "一致的正向訊號", "Ribonucleoprotein biogenesis\nCytoplasmic translation\nElectron transport\nDNA replication", 995, 150, 240, 205, ORANGE);
  callout(s, "重要負結果", "GO mitophagy\nAutophagy\nProteasomal protein catabolism\n皆未達 FDR threshold", 995, 390, 240, 165, BLUE);
  textBox(s, "完整 4,746 個 tested GO terms 保留在 CSV，main figure 只呈現非冗餘且與研究問題直接相關的 terms。", 995, 575, 235, 78, { fontSize: 14, color: MUTED });
  setNotes(s, "GO BP 與 KEGG 的 gene-set definitions 不同。KEGG mitophagy 接近 0.10，不代表 GO mitophagy 也應顯著。");
}

// Slide 9
{
  const s = presentation.slides.add();
  addBase(s, "Mitophagy 是 pathway-level exploratory signal，proteasome 是 near-threshold trend", "Hypothesis focus", 9);
  await addImage(s, path.join(FIG, "01_Main_Figures", "Fig06_mitophagy_proteasome_leading_edge_heatmap.png"), { left: 38, top: 103, width: 1210, height: 560 }, "Paired leading-edge heatmap with pathway-level GSEA statistics");
  setNotes(s, "這張直接處理先前最容易誤讀的地方。Heatmap 是 sample pattern，不能把 pathway FDR 貼到單一 gene 上。");
}

// Slide 10
{
  const s = presentation.slides.add();
  addBase(s, "三位病人的方向一致性支持 mitophagy 訊號，但證據仍屬探索性", "Paired consistency", 10);
  const table = s.tables.add({
    rows: 4,
    columns: 6,
    left: 54,
    top: 145,
    width: 1170,
    height: 300,
    columnWidths: [225, 140, 160, 160, 150, 335],
    values: [
      ["Pathway", "NES", "FDR", "3 of 3", "2 of 3", "Interpretation"],
      ["KEGG Mitophagy", "1.44", "0.0998", "14 genes up", "16 genes up", "30 of 31 leading-edge genes were up in at least 2 pairs"],
      ["KEGG Autophagy animal", "-0.93", "0.754", "16 genes down", "24 genes down", "Direction is mostly down, but pathway is not significant"],
      ["KEGG Autophagy other", "-1.34", "0.299", "4 genes down", "7 genes down", "Direction is mostly down, but pathway is not significant"],
    ],
  });
  table.styleOptions = { headerRow: true, bandedRows: true };
  table.borders.assign({ style: "solid", fill: "#CBD2D9", width: 1 });
  table.cells.block({ row: 0, column: 0, rowCount: 1, columnCount: 6 }).assign({
    fill: NAVY,
    textStyle: { typeface: FONT, fontSize: 15, bold: true, color: WHITE },
  });
  table.cells.block({ row: 1, column: 0, rowCount: 3, columnCount: 6 }).assign({
    textStyle: { typeface: FONT, fontSize: 14, color: INK },
  });
  callout(s, "統計邏輯", "基因不需要全部 3 of 3 同方向才能進入 model。Paired model 估計平均 within-patient effect，方向一致性是可信度審查，不是額外硬門檻。", 105, 495, 1070, 125, BLUE);
  setNotes(s, "Mitophagy leading edge 中 14 genes 是 3/3 上升，16 genes 是 2/3 上升，只有 1 gene 是 2/3 下降。這支持 pathway direction，但樣本數仍只有三對。");
}

// Slide 11
{
  const s = presentation.slides.add();
  addBase(s, "STAT5A nominally increases, but STAT5 downstream gene sets are not significant", "STAT5 readout", 11);
  await addImage(s, path.join(FIG, "02_Supplementary_Figures", "FigS03_STAT5_paired_expression_and_downstream_GSEA.png"), { left: 44, top: 105, width: 1195, height: 555 }, "Paired STAT5 expression and downstream gene-set GSEA");
  textBox(s, "解讀：STAT5A log2FC +1.32，nominal P = 0.0173，但 gene-level FDR = 1.00。三個 custom STAT5 downstream sets 的 FDR 都是 0.458。", 80, 626, 1120, 44, { fontSize: 17, bold: true, color: NAVY, fill: "#FFF3EB" });
  setNotes(s, "不要把 STAT5A nominal P 當成 corrected significance。更不能把 cross-dataset concordance 寫成 STAT5 已經在這個 patient dataset 被證明活化。");
}

// Slide 12
{
  const s = presentation.slides.add();
  addBase(s, "STAT5B N642H bulk RNA-seq 提供部分同方向 validation", "Cross-dataset validation", 12);
  await addImage(s, path.join(FIG, "01_Main_Figures", "Fig07_cross_dataset_GSEA_concordance.png"), { left: 38, top: 102, width: 1210, height: 560 }, "Cross-dataset GSEA concordance between relapse scRNA-seq and STAT5B N642H bulk RNA-seq");
  setNotes(s, "Bulk dataset 是 STAT5B N642H model，不是同一病人 cohort。共同 pathway direction 支持 association，但 DNA replication 的反向結果提醒我們背景差異很大。");
}

// Slide 13
{
  const s = presentation.slides.add();
  addBase(s, "FDR < 0.10 時，共有 28 個同方向上升 KEGG pathways", "Cross-dataset overlap", 13);
  await addImage(s, path.join(FIG, "01_Main_Figures", "Fig08_cross_dataset_KEGG_overlap_venn.png"), { left: 48, top: 110, width: 1180, height: 535 }, "Venn diagrams of same-direction KEGG pathway overlap");
  textBox(s, "Up overlap = 28，down overlap = 0。這個結果描述 pathway set 的交集，不代表 28 條路徑都與 mitochondrial STAT5 有直接因果關係。", 90, 622, 1100, 48, { fontSize: 18, bold: true, color: NAVY, fill: "#FFF3EB" });
  setNotes(s, "同方向 overlap 的 threshold 是探索性 FDR < 0.10。要避免把 pathways 數量當成獨立證據，因為 KEGG gene sets 之間有大量基因重疊。");
}

// Slide 14
{
  const s = presentation.slides.add();
  addBase(s, "目前資料最合理的生物學解讀", "Interpretation", 14);
  textBox(s, "Relapse-associated leukemia state", 70, 130, 560, 40, { fontSize: 27, bold: true, color: NAVY });
  const leftItems = [
    [ORANGE, "Biosynthetic demand", "Ribosome、ribosome biogenesis、spliceosome 上升"],
    [BLUE, "Replication program", "DNA replication 上升，可能反映 proliferative state"],
    [GREEN, "Mitochondrial respiration", "OXPHOS 與 electron transport 穩定上升"],
    [GOLD, "Quality-control adaptation", "KEGG mitophagy 接近 FDR 0.10，proteasome 接近但未過門檻"],
  ];
  leftItems.forEach((d, i) => {
    const y = 198 + i * 93;
    shape(s, 76, y, 18, 18, d[0], d[0], "ellipse");
    textBox(s, d[1], 110, y - 5, 300, 28, { fontSize: 19, bold: true, color: NAVY });
    textBox(s, d[2], 110, y + 27, 500, 46, { fontSize: 16, color: INK });
  });
  callout(s, "最保守、也最可發表的說法", "Relapse-associated potential malignant T cells show a coordinated oxidative and biosynthetic transcriptional program with exploratory mitophagy enrichment and a near-threshold proteasome trend.", 690, 145, 500, 205, ORANGE);
  callout(s, "機制假說仍可保留", "STAT5 activation 可能與 mitochondrial quality control 及 proteostasis program 相連，但目前只能作為 mechanistic hypothesis。", 690, 382, 500, 155, BLUE);
  callout(s, "不能直接寫的結論", "Mitochondrial STAT5 activates mitophagy, causes bortezomib resistance, or creates proteasome addiction。", 690, 558, 500, 95, RED);
  setNotes(s, "這張提供 paper discussion 可用的核心句子。用 coordinated program、exploratory enrichment、near-threshold trend，避免 activation、dependency、resistance 等過度因果字眼。");
}

// Slide 15
{
  const s = presentation.slides.add();
  addBase(s, "下一步實驗要把 association 分解成 localization、flux、dependency 與 resistance", "Validation plan", 15);
  const rows = [
    ["1", "Localization", "Mitochondrial fractionation 或 imaging", "STAT5 是否進入粒線體"],
    ["2", "Mitophagy flux", "mt-Keima、MitoTracker turnover、lysosomal blockade", "不是只看 transcript，而是直接測 flux"],
    ["3", "Perturbation", "STAT5 inhibition 或 N642H induction", "STAT5 是否改變 mitophagy 與 OXPHOS"],
    ["4", "Proteasome dependency", "Bortezomib dose response、activity probe、UPR markers", "上升是適應、依賴，還是 resistance"],
    ["5", "Combination", "STAT5 plus mitophagy or proteasome perturbation", "是否有 selective vulnerability"],
  ];
  rows.forEach((d, i) => {
    const y = 130 + i * 99;
    shape(s, 58, y, 52, 52, i < 3 ? BLUE : ORANGE, "none", "ellipse");
    textBox(s, d[0], 74, y + 9, 22, 28, { fontSize: 18, bold: true, color: WHITE });
    textBox(s, d[1], 130, y - 2, 240, 27, { fontSize: 19, bold: true, color: NAVY });
    textBox(s, d[2], 385, y - 2, 410, 55, { fontSize: 16, color: INK });
    textBox(s, d[3], 825, y - 2, 370, 55, { fontSize: 16, color: MUTED });
    shape(s, 130, y + 68, 1065, 1, "#E0E4E8");
  });
  textBox(s, "優先順序建議：先證明 STAT5 localization 與 mitophagy flux，再測 bortezomib dependency。", 160, 642, 960, 32, { fontSize: 20, bold: true, color: RED });
  setNotes(s, "交接給學員時，先讓她知道哪些問題可用 transcriptomics 回答，哪些一定要做 functional validation。Priority 是 localization、flux、dependency。");
}

// Slide 16
{
  const s = presentation.slides.add();
  addBase(s, "交接與 code availability", "Repository map", 16);
  const table = s.tables.add({
    rows: 7,
    columns: 3,
    left: 50,
    top: 126,
    width: 760,
    height: 420,
    columnWidths: [220, 220, 320],
    values: [
      ["Folder", "User-facing", "Contents"],
      ["01_Figures", "Yes", "8 main, 5 supplementary, 3 sensitivity figures"],
      ["02_Tables", "Yes", "Key results, full results, sensitivity, method inputs"],
      ["03_Raw_Data", "Local only", "GEO archive and source metadata"],
      ["04_R_Objects", "Local only", "Final Seurat and inferCNV objects"],
      ["05_Code_Availability", "Yes", "Scripts 01 to 16, wrappers, environment records"],
      ["06_Handoff_PPT", "Yes", "This reviewed handoff deck"],
    ],
  });
  table.styleOptions = { headerRow: true, bandedRows: true };
  table.borders.assign({ style: "solid", fill: "#CBD2D9", width: 1 });
  table.cells.block({ row: 0, column: 0, rowCount: 1, columnCount: 3 }).assign({
    fill: NAVY,
    textStyle: { typeface: FONT, fontSize: 15, bold: true, color: WHITE },
  });
  table.cells.block({ row: 1, column: 0, rowCount: 6, columnCount: 3 }).assign({
    textStyle: { typeface: FONT, fontSize: 14, color: INK },
  });
  callout(s, "Canonical renderer", "R/16_publication_figure_suite.R\n（inside 05_Code_Availability）", 845, 130, 360, 105, BLUE);
  callout(s, "Statistical design", "Paired patient-level pseudobulk edgeR. GSEA rank uses sign(logFC) × sqrt(F).", 845, 260, 360, 125, ORANGE);
  callout(s, "GitHub scope", "Hematologic malignancy project repository。GSE262271 是 single-cell analysis module，raw data 與 R objects 不上傳。", 845, 410, 360, 145, GREEN);
  textBox(s, "閱讀順序：START_HERE → PPT → main figures → key CSV → code。", 92, 590, 1070, 44, { fontSize: 22, bold: true, color: NAVY });
  setNotes(s, "GitHub 版本只放 code、final figures、compact CSV、README 與 PPT。Raw data 與大型 RDS 保留在本機，避免 repo 膨脹與資料授權問題。");
}

const requirements = {
  explicitTotalSlideCount: 16,
  requiredNativeTableOwnerSlides: [3, 10, 16],
  requiredNativeChartOwnerSlides: [],
};
const fontPolicy = {
  basis: "design",
  families: [FONT],
  scriptFonts: { ea: FONT },
};
const stagingDir = path.join(TMP_DIR, ".codex-finalizer");
await fs.mkdir(stagingDir, { recursive: true });
await fs.mkdir(path.join(path.dirname(TMP_DIR), ".codex-presentation-validation"), { recursive: true });
await fs.mkdir(path.dirname(FINAL_PPTX), { recursive: true });
const candidatePath = path.join(stagingDir, "GSE262271_handoff_candidate.pptx");
await (await PresentationFile.exportPptx(presentation)).save(candidatePath);

const result = await finalizePresentation({
  ...requirements,
  workspaceDir: path.dirname(TMP_DIR),
  candidatePath,
  finalPath: FINAL_PPTX,
  pythonExecutable: RUNTIME_PYTHON,
  integrityValidatorPath: path.join(SKILL_DIR, "container_tools/inspect_presentation_package_integrity.py"),
  layoutValidatorPath: path.join(SKILL_DIR, "container_tools/inspect_presentation_layout_geometry.py"),
  layoutArgs: [
    "--expected-slide-size-emu", "12192000,6858000",
    "--validate-bullet-geometry",
    "--validate-heading-fit",
    "--require-native-table-slide", "3",
    "--require-native-table-slide", "10",
    "--require-native-table-slide", "16",
  ],
  requiredNativeTableOwnerSlides: requirements.requiredNativeTableOwnerSlides,
  fontPolicy,
  verifyArtifactToolImport: true,
  receiptPath: path.join(path.dirname(TMP_DIR), ".codex-presentation-validation", "GSE262271_handoff_final.validation.json"),
});

console.log(JSON.stringify({ finalPath: FINAL_PPTX, result }, null, 2));
