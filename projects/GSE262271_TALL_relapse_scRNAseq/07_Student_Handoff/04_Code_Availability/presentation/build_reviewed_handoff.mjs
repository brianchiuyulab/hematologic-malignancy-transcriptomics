import fs from 'node:fs/promises';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
import {Presentation,PresentationFile} from '@oai/artifact-tool';
const {SKILL_DIR,TMP_DIR,REVIEWED_ROOT,FINAL_PPTX,RUNTIME_PYTHON,WORKSPACE_DIR}=process.env;
for(const p of [SKILL_DIR,TMP_DIR,REVIEWED_ROOT,FINAL_PPTX,RUNTIME_PYTHON,WORKSPACE_DIR])if(!path.isAbsolute(p??''))throw Error('All runtime paths must be absolute');
const {finalizePresentation}=await import(pathToFileURL(path.join(SKILL_DIR,'container_tools/artifact_tool_utils.mjs')).href);
const summary=JSON.parse(await fs.readFile(path.join(REVIEWED_ROOT,'04_Code_Availability','summary_for_slides.json'),'utf8'));
const pr=Presentation.create({slideSize:{width:1280,height:720}});
let page=0;
const FONT='Microsoft JhengHei',INK='#142F43',MUTED='#4D5B66',BLUE='#0072B2',ORANGE='#D55E00';
function txt(s,t,x,y,w,h,size=24,bold=false,color=INK){const a=s.shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});a.text=t;a.text.style={typeface:FONT,fontSize:size,bold,color,autoFit:'shrinkText'};return a;}
function slide(title,section){const s=pr.slides.add();s.background.fill='#FFFFFF';txt(s,section,50,20,1100,22,13,true,ORANGE);txt(s,title,50,55,1170,50,32,true);txt(s,'GSE262271 / GSE218858  |  Reviewed 15 September 2026',50,687,1080,18,11,false,MUTED);txt(s,String(++page).padStart(2,'0'),1190,687,40,18,11);return s;}
function notes(s,t){s.speakerNotes.textFrame.setText(t);}
async function image(s,n,x=45,y=115,w=1190,h=535){s.images.add({blob:await fs.readFile(path.join(REVIEWED_ROOT,'02_Figures',n+'.png')),contentType:'image/png',alt:n,fit:'contain',position:{left:x,top:y,width:w,height:h}});}
function table(s,values,widths,y=135,h=350,size=18){const t=s.tables.add({rows:values.length,columns:values[0].length,left:55,top:y,width:1170,height:h,columnWidths:widths,values});t.rows[0].height=44;for(let i=1;i<values.length;i++)t.rows[i].height=(h-44)/(values.length-1);t.styleOptions={headerRow:true,bandedRows:false};t.borders.assign({style:'solid',fill:'#D3DDE3',width:.6});t.cells.block({row:0,column:0,rowCount:1,columnCount:values[0].length}).assign({fill:INK,textStyle:{typeface:FONT,fontSize:size,bold:true,color:'#FFFFFF'}});t.cells.block({row:1,column:0,rowCount:values.length-1,columnCount:values[0].length}).assign({textStyle:{typeface:FONT,fontSize:size,color:INK}});return t;}
const q=x=>x<.001?x.toExponential(2):x.toFixed(4);
const core=summary.counts.find(x=>x.Dataset==='scRNAseq_within15'),bulk=summary.counts.find(x=>x.Dataset==='bulk_STAT5B'),overlap=summary.intersection.find(x=>x.analysis==='within15');
// 1
{
 const s=slide('T-ALL relapse 與 STAT5B 活化','RESEARCH QUESTION');
 txt(s,'兩組 transcriptome 有哪些共同的 pathway 變化？',60,175,1120,90,38,true);
 txt(s,'先比較完整 KEGG enrichment 與同方向交集，再檢視 proteasome、mitophagy 和相關生物機制。',62,300,1110,105,28);
 txt(s,'GSE262271：3 位病人 diagnosis / relapse 配對\nGSE218858：STAT5B N642H 與對照各 3 個 bulk samples',62,465,1100,100,24);
 txt(s,'主分析固定 within-sample CNV top 15%，pathway 與交集皆採完整 KEGG FDR < 0.05。',62,609,1110,44,20,false,MUTED);
 notes(s,'研究目的由先前 STAT5 粒線體、mitophagy、proteasome 假說延伸。跨物種、不同實驗情境的 pathway 同方向只屬支持性關聯。不可稱為 bortezomib sensitivity、直接因果或獨立人類臨床驗證。');
}
// 2
{
 const s=slide('資料集與比較方向','DATASETS');
 table(s,[['Dataset','設計','模型與正方向'],['GSE262271\nHuman scRNA-seq','3 paired patients\nP1: M104 / M127\nP2: M143 / M148\nP3: M187 / M187r','Patient-sample pseudobulk\nedgeR: patient + condition\nPositive = relapse higher'],['GSE218858\nMouse bulk RNA-seq','Rag2-/- DN stage\nSTAT5B N642H: n = 3\nControl: n = 3\nIndependent biological samples','DESeq2: condition\nR2S5b_DN vs R2b_DN\nPositive = N642H higher']],[240,440,490],130,395,20);
 txt(s,'Single-cell 的 replicate 是病人 sample，不是每一顆 cell。兩组資料不可直接合併成同一模型。',65,565,1120,70,24,true);
 notes(s,'Sources: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262271 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218858 . Bulk genotype uses GEO characteristics rather than inconsistent sample-title shorthand. Full selected sample IDs in bulk_sample_design.csv.');
}
// 3
{
 const s=slide('Cell lineage annotation','CELL IDENTITY');
 await image(s,'Fig01A_cell_lineage_UMAP',40,113,810,548);
 txt(s,'先保留 marker 支持的廣義分類',870,145,350,70,25,true);
 txt(s,'T-lineage candidate\nB cell / myeloid / pDC\nErythroid\nCytotoxic T/NK',875,250,330,155,22);
 txt(s,'HSC-like、quiescent 或 ER-stress 是待驗證的細胞狀態，不直接當成確定的 T-ALL subtype。',875,460,330,165,22,false,ORANGE);
 notes(s,'Reviewed labels are stored separately. Original Seurat annotations retained for traceability. Broad lineage assignments are supported by cluster markers but not independent reference mapping, flow cytometry or mutation/clonotype validation.');
}
// 4
{
 const s=slide('每個 cluster 的 marker 證據','ANNOTATION EVIDENCE');
 await image(s,'Fig01B_annotation_marker_dotplot',30,118,1220,490);
 txt(s,'大小 = 表達細胞比例；顏色 = 同一基因跨 cluster 的平均表達。先辨 lineage，再另外判斷 CNV-high。',60,622,1160,47,20);
 notes(s,'Figure 01B contains all 20 clusters and canonical marker panels. These descriptive scores are not donor-level significance tests. HTO singlet filtering does not rule out all same-sample doublets. Erythroid and cytotoxic T/NK labels remain coarse.');
}
// 5
{
 const s=slide('CNV 篩選與 paired pseudobulk','ANALYSIS DEFINITIONS');
 await image(s,'Fig02_paired_sample_cell_counts',40,112,800,460);
 txt(s,'主分析：within-sample top 15%\n1,926 cells，最小 sample 107',855,150,370,105,22,true);
 txt(s,'Non-T reference q95 = 0.0269\n全體候選 T 細胞有 15.4% 高於此值',855,305,370,110,22);
 txt(s,'~ patient_pair + condition\n以 condition coefficient 排名\nGSEA: sign(logFC) × √F',855,455,370,125,21);
 txt(s,'每個 sample 固定取 15%，無法用保留比例來證明 relapse 的 malignant population 增加。',60,605,1150,60,22,false,ORANGE);
 notes(s,'867 non-T reference cells yield empirical q95=0.0269455, and 1978/12828=15.419% candidates exceed it. The absolute cutoff is not equivalent to selecting 15% within every sample. The exploratory 15% was considered after outcome inspection, so it is not outcome-independent or unbiased.');
}
// 6
{
 const s=slide('Single-cell：完整 KEGG enrichment 概覽','DISCOVERY');
 await image(s,'Fig03B_scRNAseq_slide_overview',35,120,990,550);
 txt(s,`FDR < 0.05\n\n上調 ${core.up_FDR005} 條\n下調 ${core.down_FDR005} 條`,1025,190,220,185,24,true);
 txt(s,'主圖依 FDR 排序。完整結果在 Table04，顯著清單在 Table11。',1025,430,220,170,21);
 notes(s,'Figure03B uses exploratory within-sample top15, original stored QC, paired edgeR ranking, fixed seed 20260915 and full KEGG BH. Up to eight pathways per direction are displayed on this slide. Figure03 displays up to fifteen per direction and Table04 contains the complete catalog. Disease-labeled gene sets often share metabolic/proteostasis genes and do not diagnose those diseases.');
}
// 7
{
 const s=slide('Bulk：STAT5B N642H 的完整 KEGG enrichment','STAT5B PERTURBATION');
 await image(s,'Fig04B_bulk_slide_overview',35,120,990,550);
 txt(s,`FDR < 0.05\n\n上調 ${bulk.up_FDR005} 條\n下調 ${bulk.down_FDR005} 條`,1025,190,220,185,24,true);
 txt(s,'排名來自全部通過基因過濾的 DESeq2 Wald statistics。完整結果在 Table05。',1025,430,220,185,21);
 notes(s,'R2S5b_DN vs R2b_DN, independent n3 vs3. Enrichr KEGG2019 Mouse catalog differs from the frozen human KEGG catalog. Each complete catalog has its own BH denominator. Refit from deposited raw count matrix.');
}
// 8
{
 const s=slide('同方向 pathway 交集','CROSS-DATASET OVERLAP');
 await image(s,'Fig05_direction_matched_intersection',35,140,920,470);
 txt(s,`Top 15% 與 bulk\n\n共同上調 ${overlap.shared_up} 條\n共同下調 ${overlap.shared_down} 條`,970,173,270,195,25,true);
 txt(s,'只比較兩邊都測到的 pathway。\n兩邊都需 FDR < 0.05，且 NES 同號。\n\n其他 cutoff 與原始 top 25% 放在補充分析。',970,405,270,220,21);
 notes(s,'Figure05 restricts the universe to common tested canonical pathway names. It does not call pathways missing from one catalog non-significant. Circle area is schematic. Table06 has all membership; Table07 lists all40 primary shared pathways. TableS05 retains global25 sensitivity. Primary display uses within15 as requested, while cutoff selection remains exploratory.');
}
// 9
{
 const s=slide('交集中的 pathway 與 NES','SHARED PROGRAMS');
 await image(s,'Fig06B_shared_strongest',40,112,1200,550);
 notes(s,'All same-direction intersections are provided in Table07 and full Figure06. This slide shows at most15 pathways ranked by the worse FDR across the two datasets, not a target-selected list. Each tile contains NES and whole-collection significance.');
}
// 10
{
 const s=slide('Proteasome 與 mitophagy 的結果核對','TARGET PATHWAYS');
 await image(s,'Fig07_target_pathway_audit',35,143,1210,445);
 txt(s,'Proteasome、mitophagy 均上調，兩組資料的完整 KEGG FDR 都小於 0.05。',65,608,1150,59,22,true,ORANGE);
 notes(s,'Numeric results are read from regenerated tables, not hard-coded historical values. q means whole KEGG BH, not nominal P or two-target correction. Within15 is exploratory because cutoff was selected after multiverse results.');
}
// 11
{
 const s=slide('三位病人的變化方向','PAIRED DIRECTION');
 await image(s,'Fig08_paired_leading_edge_direction',40,125,900,520);
 txt(s,'每條線是一位病人',970,200,260,65,25,true);
 txt(s,'Leading-edge genes 由同一份 GSEA 選出。\n\n這裡只描述三位病人的方向，不作額外的獨立驗證。',970,300,260,230,22);
 notes(s,'Three patients limit inference. A gene need not rise in all three pairs for fitting a paired model. Gene-level DEG FDR and gene-set GSEA FDR answer different hypotheses.');
}
// 12
{
 const s=slide('CNV top 15% 設定與 reference 對照','CUTOFF INTERPRETATION');
 txt(s,'每個 sample 保留 CNV score 最高的 15%，作為相對高 CNV 的 T-lineage core',62,145,1150,72,28,true);
 table(s,[['規則','怎麼選','本次保留'],['主分析：sample top 15%','每個 sample 各自選 CNV score 最高 15%','1,926 cells\n最小 sample 107'],['Non-T reference q95','每顆 cell 都與同一數值 0.0269 比較','1,978 cells\n最小 sample 66']],[270,575,325],245,210,21);
 txt(s,'Reference q95 提供比例尺度對照。Top 15% 中有 1,407 / 1,926 cells（73.1%）高於 q95，兩種規則選到的細胞並不相同。',65,490,1120,87,23);
 txt(s,'Top 15% 在 sensitivity 檢查後固定為主分析。FDR < 0.05 是 pathway 門檻，沒有額外校正選擇分析條件的影響。',65,593,1120,66,21,false,ORANGE);
 notes(s,'inferCNV score is RNA-based deviation from reference, not measured DNA copy number, a probability of cancer, or a patient-independent validated specificity threshold. Reference sample/cell-type composition is unbalanced. Official workflow: https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV');
}
// 13
{
 const s=slide('生物意義與待驗證問題','INTERPRETATION');
 txt(s,'共同富集支持：復發與 STAT5B 活化可能伴隨相似的蛋白質穩態與代謝程式。',65,150,1120,105,29,true);
 txt(s,'Proteasome / mitophagy 的 transcript enrichment 可以用來形成後續假說。它尚未證明蛋白質降解活性、mitophagy flux 或 STAT5 的粒線體定位。',65,310,1120,135,26);
 txt(s,'這兩組資料也沒有直接比較 bortezomib responder 與 resistant samples，因此不能從交集直接推論 proteasome inhibitor resistance。',65,503,1120,125,25,false,ORANGE);
 notes(s,'These are limits of the actual contrasts. General autophagy is not interchangeable with mitophagy. Inspect shared leading-edge genes and full FDR tables before choosing individual targets.');
}
// 14
{
 const s=slide('交接檔案與可重現性','HANDOFF');
 table(s,[['資料','從哪裡看'],['PPT / 說明','01_Summary；根目錄 00_READ_ME.txt'],['Figures','02_Figures：PNG 閱讀、PDF 排版'],['Tables','03_Tables：完整 KEGG、同方向交集、target FDR、cell counts'],['Code / 稽核記錄','04_Code_Availability：R scripts 22–25、原始 bulk pipeline、版本與限制'],['Raw data / R objects','05_Data_Locations：原始 accession 與本機路徑，不重複搬大型檔案']],[270,900],126,380,21);
 txt(s,'github.com/brianchiuyulab/hematologic-malignancy-transcriptomics',65,548,1130,45,23,true,BLUE);
 txt(s,'投稿前仍需補強 malignant-cell 的獨立驗證。此包提供可追溯的探索型分析，不標榜已完成審稿驗證。',65,611,1120,59,22,false,ORANGE);
 notes(s,'Code22 recomputes count-based models and high-precision GSEA. Code23 renders reader figures/tables. The audit reports historical discrepancies rather than silently replacing provenance. No raw data, patient identifiers or large R objects are uploaded.');
}
await fs.mkdir(TMP_DIR,{recursive:true});await fs.mkdir(path.dirname(FINAL_PPTX),{recursive:true});
const candidatePath=path.join(TMP_DIR,'candidate_reviewed.pptx');await (await PresentationFile.exportPptx(pr)).save(candidatePath);
await finalizePresentation({workspaceDir:WORKSPACE_DIR,candidatePath,finalPath:FINAL_PPTX,pythonExecutable:RUNTIME_PYTHON,integrityValidatorPath:path.join(SKILL_DIR,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(SKILL_DIR,'container_tools/inspect_presentation_layout_geometry.py'),layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit','--require-native-table-slide','2','--require-native-table-slide','12','--require-native-table-slide','14'],explicitTotalSlideCount:14,requiredNativeTableOwnerSlides:[2,12,14],requiredNativeChartOwnerSlides:[],fontPolicy:{basis:'design',families:[FONT]},verifyArtifactToolImport:true,receiptPath:path.join(TMP_DIR,'reviewed.validation.json')});
console.log(FINAL_PPTX);
