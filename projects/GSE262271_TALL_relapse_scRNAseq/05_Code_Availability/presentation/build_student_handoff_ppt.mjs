import fs from "node:fs/promises";
import path from "node:path";
import { pathToFileURL } from "node:url";
import { Presentation, PresentationFile } from "@oai/artifact-tool";

const { SKILL_DIR, TMP_DIR, PROJECT_ROOT, FINAL_PPTX, RUNTIME_PYTHON } = process.env;
if (![SKILL_DIR, TMP_DIR, PROJECT_ROOT, FINAL_PPTX, RUNTIME_PYTHON].every((p) => path.isAbsolute(p ?? ""))) {
  throw new Error("SKILL_DIR, TMP_DIR, PROJECT_ROOT, FINAL_PPTX and RUNTIME_PYTHON must be absolute paths");
}

const { finalizePresentation } = await import(
  pathToFileURL(path.join(SKILL_DIR, "container_tools/artifact_tool_utils.mjs")).href,
);

await fs.mkdir(TMP_DIR, { recursive: true });
await fs.mkdir(path.dirname(FINAL_PPTX), { recursive: true });

const W = 1280;
const H = 720;
const FONT = "Microsoft JhengHei";
const NAVY = "#15263C";
const BLUE = "#2878A8";
const ORANGE = "#D55E00";
const TEAL = "#148F77";
const RED = "#B23A48";
const INK = "#252B33";
const MUTED = "#5D6875";
const PALE = "#F5F7F8";
const LIGHT_BLUE = "#EAF3F8";
const LIGHT_ORANGE = "#FFF0E7";
const WHITE = "#FFFFFF";

const FIG = path.join(PROJECT_ROOT, "01_Figures");
const presentation = Presentation.create({ slideSize: { width: W, height: H } });

function rect(slide, left, top, width, height, fill = "none", lineFill = "none", geometry = "rect") {
  return slide.shapes.add({
    geometry,
    position: { left, top, width, height },
    fill,
    line: { fill: lineFill, width: lineFill === "none" ? 0 : 1 },
  });
}

function textBox(slide, text, left, top, width, height, opts = {}) {
  const box = rect(slide, left, top, width, height, opts.fill ?? "none", opts.line ?? "none", opts.geometry ?? "rect");
  box.text = text;
  box.text.style = {
    typeface: FONT,
    fontSize: opts.fontSize ?? 20,
    bold: opts.bold ?? false,
    color: opts.color ?? INK,
    autoFit: opts.autoFit ?? "shrinkText",
  };
  return box;
}

function base(slide, title, section, page) {
  slide.background.fill = WHITE;
  rect(slide, 0, 0, W, 8, ORANGE);
  textBox(slide, section.toUpperCase(), 50, 20, 360, 22, { fontSize: 12, bold: true, color: ORANGE });
  textBox(slide, title, 50, 45, 1175, 45, { fontSize: 29, bold: true, color: NAVY });
  rect(slide, 50, 683, 1175, 1, "#D7DDE2");
  textBox(slide, "GSE262271 paired scRNA-seq and GSE218858 STAT5B N642H bulk RNA-seq", 50, 690, 920, 16, { fontSize: 9, color: MUTED });
  textBox(slide, String(page).padStart(2, "0"), 1184, 690, 40, 16, { fontSize: 9, bold: true, color: MUTED });
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

function note(slide, text) {
  slide.speakerNotes.textFrame.setText(text);
}

function styleTable(table, rows, cols, bodySize = 15) {
  table.styleOptions = { headerRow: true, bandedRows: true };
  table.borders.assign({ style: "solid", fill: "#CAD2D9", width: 1 });
  table.cells.block({ row: 0, column: 0, rowCount: 1, columnCount: cols }).assign({
    fill: NAVY,
    textStyle: { typeface: FONT, fontSize: 15, bold: true, color: WHITE },
  });
  table.cells.block({ row: 1, column: 0, rowCount: rows - 1, columnCount: cols }).assign({
    textStyle: { typeface: FONT, fontSize: bodySize, color: INK },
  });
}

// 1. Cover
{
  const s = presentation.slides.add();
  s.background.fill = NAVY;
  rect(s, 0, 0, 16, H, ORANGE);
  textBox(s, "HEMATOLOGIC MALIGNANCY PROJECT", 72, 62, 760, 25, { fontSize: 14, bold: true, color: "#F4B66B" });
  textBox(s, "T-ALL relapse transcriptomics", 72, 123, 1040, 65, { fontSize: 46, bold: true, color: WHITE });
  textBox(s, "Single-cell discovery and STAT5B bulk validation", 74, 205, 970, 43, { fontSize: 27, color: "#DCE5EC" });
  rect(s, 74, 292, 1080, 2, ORANGE);
  textBox(s, "交接重點", 74, 332, 150, 30, { fontSize: 18, bold: true, color: "#F4B66B" });
  textBox(s, "Relapse-associated potential malignant T cells show an oxidative and proteostasis program. Mitophagy and proteasome signals strengthen in a high-confidence CNV core.", 74, 375, 1080, 120, { fontSize: 27, bold: true, color: WHITE });
  textBox(s, "8-slide student handoff | September 2026", 74, 606, 680, 28, { fontSize: 16, color: "#BFCAD4" });
  note(s, "這是精簡交接版。先分清兩組資料，再分清 frozen primary top 25% 與 exploratory high-confidence top 15%。");
}

// 2. Datasets and statistical design
{
  const s = presentation.slides.add();
  base(s, "兩組關鍵資料回答不同層次的問題", "Datasets", 2);
  const t = s.tables.add({
    rows: 3,
    columns: 4,
    left: 54,
    top: 125,
    width: 1170,
    height: 270,
    columnWidths: [180, 330, 300, 360],
    values: [
      ["Dataset", "Design", "核心比較", "在專案中的角色"],
      ["GSE262271", "Human T-ALL scRNA-seq\n3 patients, paired diagnosis and relapse", "Potential malignant T cells\nRelapse versus diagnosis", "疾病復發訊號與 CNV-high population 的主要分析"],
      ["GSE218858", "Mouse bulk RNA-seq\nSTAT5B N642H perturbation model", "STAT5B N642H versus control", "檢查 pathway direction 是否與 STAT5B hyperactivation 一致"],
    ],
  });
  styleTable(t, 3, 4, 15);
  textBox(s, "Single-cell 的統計單位是 patient sample，不是 cell", 64, 438, 710, 34, { fontSize: 22, bold: true, color: NAVY });
  textBox(s, "每個 sample 先做 pseudobulk，再用 edgeR model  ~ patient_pair + condition。GSEA 使用這個 paired model 產生的 genome-wide ranking。", 64, 484, 1120, 75, { fontSize: 19, color: INK, fill: LIGHT_BLUE });
  textBox(s, "跨資料集一致只支持 association，不證明 STAT5 直接活化 mitophagy 或造成 bortezomib resistance。", 64, 590, 1120, 44, { fontSize: 18, bold: true, color: RED });
  note(s, "GSE262271 是 paired human disease dataset。GSE218858 是 mouse STAT5B N642H bulk model，屬獨立 pathway-level validation。");
}

// 3. CNV cutoff and robustness
{
  const s = presentation.slides.add();
  base(s, "Top 15% 有 reference anchor，但仍須標示為 exploratory", "CNV definition", 3);
  await addImage(
    s,
    path.join(FIG, "03_Sensitivity_Checks", "FigSens10_CNV_fraction_definition_and_pathway_robustness.png"),
    { left: 40, top: 100, width: 1200, height: 545 },
    "CNV burden distributions, top-15 percent cutoffs, pathway robustness and retained cell counts",
  );
  textBox(s, "non-T reference q95 = 0.0269；15.4% candidate T-lineage cells 超過此正常參考上限。10% 至 50% 的 NES 方向大致一致。", 66, 640, 1148, 34, { fontSize: 16, bold: true, color: NAVY, fill: LIGHT_ORANGE });
  note(s, "Outcome-independent anchor 指的是 cutoff 來自 non-T reference CNV distribution，而不是 relapse-versus-diagnosis pathway P value。Top 15% 的精確比例仍來自敏感度分析，所以必須明示 exploratory。");
}

// 4. Frozen primary versus high-confidence core
{
  const s = presentation.slides.add();
  base(s, "Frozen primary 與 high-confidence core 不能混成同一分析", "Analysis hierarchy", 4);
  const t = s.tables.add({
    rows: 4,
    columns: 6,
    left: 42,
    top: 125,
    width: 1195,
    height: 360,
    columnWidths: [220, 160, 160, 180, 180, 295],
    values: [
      ["Analysis", "Cells", "Min/sample", "Mitophagy", "Proteasome", "正確定位"],
      ["Frozen primary\nGlobal top 25%", "3,207", "140", "NES 1.44\nFDR 0.100", "NES 1.53\nFDR 0.109", "保留原始分析與完整 provenance"],
      ["Fixed secondary\nWithin-sample top 25%", "3,210", "179", "NES 1.48\nFDR 0.081", "NES 1.60\nFDR 0.068", "樣本內平衡，探索型 FDR < 0.10"],
      ["High-confidence core\nWithin-sample top 15%", "1,926", "107", "NES 1.69\nFDR 0.015", "NES 1.89\nFDR 0.0059", "Reference-anchored exploratory sensitivity"],
    ],
  });
  styleTable(t, 4, 6, 14);
  textBox(s, "15% 不可因為 P value 最小而被稱為 unbiased", 65, 526, 700, 32, { fontSize: 22, bold: true, color: RED });
  textBox(s, "可辯護的理由是 non-T q95 提供與 outcome 無關的生物學 anchor，加上每個 sample 仍有至少 107 顆細胞，且其他合理 cutoff 的方向相符。", 65, 574, 1120, 65, { fontSize: 18, color: INK, fill: PALE });
  note(s, "Top 15% 的 pathway FDR 不能反過來當 cutoff 選擇依據。報告時要同時展示 primary、secondary、top 15 sensitivity。");
}

// 5. Primary GSEA
{
  const s = presentation.slides.add();
  base(s, "Primary GSEA 指向 oxidative、biosynthetic 與 replication programs", "Single-cell result", 5);
  await addImage(
    s,
    path.join(FIG, "01_Main_Figures", "Fig04_KEGG_GSEA_key_pathways.png"),
    { left: 40, top: 100, width: 1200, height: 545 },
    "Paired-pseudobulk KEGG GSEA in potential malignant T cells",
  );
  textBox(s, "Frozen top 25%：KEGG mitophagy FDR 0.0998；proteasome FDR 0.109；general autophagy 不顯著。", 74, 642, 1135, 32, { fontSize: 17, bold: true, color: NAVY, fill: LIGHT_ORANGE });
  note(s, "主要顯著 pathways 包括 oxidative phosphorylation、ribosome、ribosome biogenesis、DNA replication 與 spliceosome。Mitophagy 為探索型，proteasome 接近 0.10。");
}

// 6. Paired direction in top-15 core
{
  const s = presentation.slides.add();
  base(s, "Top 15% core 中兩個 target 都是 3/3 patients 上升", "Paired target check", 6);
  const t = s.tables.add({
    rows: 7,
    columns: 5,
    left: 60,
    top: 125,
    width: 1160,
    height: 385,
    columnWidths: [225, 120, 230, 230, 280],
    values: [
      ["Pathway", "Patient", "Diagnosis score", "Relapse score", "Paired change"],
      ["Proteasome", "P1", "6.794", "6.927", "+0.134"],
      ["Proteasome", "P2", "6.258", "6.971", "+0.713"],
      ["Proteasome", "P3", "6.305", "6.883", "+0.578"],
      ["Mitophagy", "P1", "6.246", "6.944", "+0.698"],
      ["Mitophagy", "P2", "5.927", "6.125", "+0.197"],
      ["Mitophagy", "P3", "5.565", "6.194", "+0.629"],
    ],
  });
  styleTable(t, 7, 5, 15);
  textBox(s, "Full KEGG correction", 76, 542, 260, 25, { fontSize: 17, bold: true, color: ORANGE });
  textBox(s, "Proteasome NES 1.89, FDR 0.0059\nMitophagy NES 1.69, FDR 0.015", 76, 570, 1070, 60, { fontSize: 21, bold: true, color: NAVY });
  textBox(s, "N = 3 pairs。方向一致性增加可信度，但不消除 small-sample uncertainty。", 76, 638, 1070, 25, { fontSize: 16, color: RED });
  note(s, "這張呈現 within-patient paired direction，不把 cells 當 replicate。Scores 是 leading-edge mean expression；FDR 是 pathway-level GSEA FDR。");
}

// 7. Cross-dataset evidence
{
  const s = presentation.slides.add();
  base(s, "STAT5B N642H bulk dataset 支持同方向 pathway association", "Cross-dataset validation", 7);
  await addImage(
    s,
    path.join(FIG, "01_Main_Figures", "Fig07_cross_dataset_GSEA_concordance.png"),
    { left: 42, top: 108, width: 745, height: 510 },
    "Cross-dataset pathway concordance",
  );
  await addImage(
    s,
    path.join(FIG, "01_Main_Figures", "Fig08_cross_dataset_KEGG_overlap_venn.png"),
    { left: 800, top: 126, width: 430, height: 470 },
    "Same-direction KEGG overlap at FDR below 0.10",
  );
  textBox(s, "Both up: OXPHOS, mitophagy, proteasome, ribosome", 74, 620, 700, 28, { fontSize: 18, bold: true, color: NAVY });
  textBox(s, "Shared KEGG at FDR < 0.10: 28 up, 0 down", 815, 620, 400, 28, { fontSize: 18, bold: true, color: ORANGE });
  note(s, "Cross-dataset overlap 只表示方向與 pathway 層級相符。物種、bulk versus single-cell 與 experimental context 不同，不能當作 direct causal validation。");
}

// 8. Handoff map
{
  const s = presentation.slides.add();
  base(s, "交接閱讀順序與 code availability", "Handoff", 8);
  const t = s.tables.add({
    rows: 6,
    columns: 3,
    left: 54,
    top: 120,
    width: 1170,
    height: 360,
    columnWidths: [170, 360, 640],
    values: [
      ["順序", "位置", "先理解什麼"],
      ["1", "00_START_HERE_請先看.txt", "兩組 dataset、主要結論、primary 與 top 15% 的關係"],
      ["2", "01_Summary", "本簡報，8 張 slide 完成全案導覽"],
      ["3", "02_Key_Figures", "只保留 6 張主圖和 1 張 cutoff sensitivity 圖"],
      ["4", "03_Key_Tables", "重點 GSEA、paired direction、cross-dataset overlap"],
      ["5", "04_Code_Availability", "GitHub link、script map、可重現性與限制"],
    ],
  });
  styleTable(t, 6, 3, 15);
  textBox(s, "GitHub", 64, 520, 120, 28, { fontSize: 20, bold: true, color: ORANGE });
  textBox(s, "github.com/brianchiuyulab/hematologic-malignancy-transcriptomics", 180, 520, 950, 28, { fontSize: 19, bold: true, color: NAVY });
  textBox(s, "Raw data 與大型 R objects 留在本機，不上傳 GitHub；交接包提供絕對路徑與 accession。", 64, 573, 1110, 35, { fontSize: 18, color: INK, fill: LIGHT_BLUE });
  textBox(s, "最先看的三張圖：CNV cutoff robustness、single-cell KEGG GSEA、cross-dataset concordance。", 64, 628, 1110, 32, { fontSize: 18, bold: true, color: RED });
  note(s, "交接包是 reader-facing subset。完整專案與全部 sensitivity tables 仍保留在原 module，GitHub 只排除 large raw and R object files。");
}

const candidatePath = path.join(TMP_DIR, "candidate_student_handoff.pptx");
await (await PresentationFile.exportPptx(presentation)).save(candidatePath);

await finalizePresentation({
  workspaceDir: path.dirname(PROJECT_ROOT),
  candidatePath,
  finalPath: FINAL_PPTX,
  pythonExecutable: RUNTIME_PYTHON,
  integrityValidatorPath: path.join(SKILL_DIR, "container_tools/inspect_presentation_package_integrity.py"),
  layoutValidatorPath: path.join(SKILL_DIR, "container_tools/inspect_presentation_layout_geometry.py"),
  layoutArgs: [
    "--expected-slide-size-emu", "12192000,6858000",
    "--validate-heading-fit",
    "--require-native-table-slide", "2",
    "--require-native-table-slide", "4",
    "--require-native-table-slide", "6",
    "--require-native-table-slide", "8",
  ],
  explicitTotalSlideCount: 8,
  requiredNativeTableOwnerSlides: [2, 4, 6, 8],
  requiredNativeChartOwnerSlides: [],
  fontPolicy: { basis: "design", families: [FONT] },
  verifyArtifactToolImport: true,
  receiptPath: path.join(TMP_DIR, `${path.basename(FINAL_PPTX)}.validation.json`),
});

console.log(FINAL_PPTX);
