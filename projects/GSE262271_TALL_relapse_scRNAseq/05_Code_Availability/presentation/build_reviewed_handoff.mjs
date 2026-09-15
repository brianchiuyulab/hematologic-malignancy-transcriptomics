import fs from 'node:fs/promises';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
import {Presentation,PresentationFile} from '@oai/artifact-tool';
const {SKILL_DIR,TMP_DIR,REVIEWED_ROOT,FINAL_PPTX,RUNTIME_PYTHON,WORKSPACE_DIR}=process.env;
for(const p of [SKILL_DIR,TMP_DIR,REVIEWED_ROOT,FINAL_PPTX,RUNTIME_PYTHON,WORKSPACE_DIR])if(!path.isAbsolute(p??''))throw Error('All runtime paths must be absolute');
const {finalizePresentation}=await import(pathToFileURL(path.join(SKILL_DIR,'container_tools/artifact_tool_utils.mjs')).href);
const summary=JSON.parse(await fs.readFile(path.join(REVIEWED_ROOT,'04_Code_Availability','summary_for_slides.json'),'utf8'));
const pr=Presentation.create({slideSize:{width:1280,height:720}});
const FONT='Microsoft JhengHei',INK='#202020',MUTED='#555555';let page=0;
function txt(s,t,x,y,w,h,size=24,bold=false){const a=s.shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});a.text=t;a.text.style={typeface:FONT,fontSize:size,bold,color:INK,autoFit:'none'};return a;}
function slide(title){const s=pr.slides.add();s.background.fill='#FFFFFF';txt(s,title,55,35,1170,62,34,true);page++;return s;}
function notes(s,t){s.speakerNotes.textFrame.setText(t);}
async function image(s,n,x=45,y=108,w=1190,h=554){s.images.add({blob:await fs.readFile(path.join(REVIEWED_ROOT,'02_Figures',n+'.png')),contentType:'image/png',alt:n,fit:'contain',position:{left:x,top:y,width:w,height:h}});}
function table(s,values,widths,y=150,h=360,size=23){const t=s.tables.add({rows:values.length,columns:values[0].length,left:60,top:y,width:1160,height:h,columnWidths:widths,values});t.rows[0].height=52;for(let i=1;i<values.length;i++)t.rows[i].height=(h-52)/(values.length-1);t.styleOptions={headerRow:true,bandedRows:false};t.borders.assign({style:'solid',fill:'#D5D5D5',width:.55});t.cells.block({row:0,column:0,rowCount:1,columnCount:values[0].length}).assign({fill:'#F0F0F0',textStyle:{typeface:FONT,fontSize:size,bold:true,color:INK}});t.cells.block({row:1,column:0,rowCount:values.length-1,columnCount:values[0].length}).assign({textStyle:{typeface:FONT,fontSize:size,color:INK}});return t;}
const core=summary.counts.find(x=>x.Dataset==='scRNAseq_within15'),bulk=summary.counts.find(x=>x.Dataset==='bulk_STAT5B'),overlap=summary.intersection.find(x=>x.analysis==='within15');
if(core.up_FDR005!==61||bulk.up_FDR005!==107||overlap.shared_up!==40||overlap.shared_down!==0)throw Error('Unexpected analysis summary');
const methods='Human GSE262271: 3 paired T-ALL patients. Potential malignant T cells are candidate T-lineage cells selected by within-sample CNV top15 under the stored original QC. Patient-sample pseudobulk raw counts are analysed with edgeR ~patient_pair + condition_simple. Rank = sign(logFC)*sqrt(F). Mouse GSE218858: STAT5B N642H vs control, independent n3 vs3, DESeq2 ~condition, ranked by Wald statistic. fgseaMultilevel uses fixed seed20260915, sampleSize501, nPermSimple50000, eps0. BH FDR is calculated separately over the full tested human KEGG353 and mouse KEGG300 collections. Cross-dataset overlap uses275 common tested canonical pathway names.';
const limits='Top15 was selected after exploratory sensitivity analysis. Within-catalog BH FDR does not adjust for selection across analysis settings. Pooled reference q95 is a different absolute threshold, not an outcome-independent derivation of within-sample15%. RNA-based CNV and broad marker annotations do not independently validate malignancy. GSEA measures transcriptional enrichment, not mitophagy flux, proteasome activity or drug response.';
{
 const s=slide('T-ALL relapse 與 STAT5B 活化');
 txt(s,'共同的轉錄路徑',65,195,1120,95,46,true);
 txt(s,'配對 single-cell RNA-seq 與 STAT5B N642H bulk RNA-seq',67,320,1110,75,27);
 txt(s,'GSE262271  /  GSE218858',67,470,1110,50,26);
 txt(s,'研究問題：復發與 STAT5B 活化是否伴隨相同的細胞程式？',67,572,1110,65,25);
 notes(s,'This handoff addresses pathway concordance between relapse-associated T-ALL transcription and constitutively active STAT5B. The biological hypothesis concerns mitochondrial and protein quality control. The comparison first examines complete KEGG enrichment and its direction-matched overlap. '+limits);
}
{
 const s=slide('資料集與分析設計');
 table(s,[['資料集','生物樣本','GSEA 排名來源'],['GSE262271\nHuman T-ALL scRNA-seq','3 位病人\nDiagnosis / relapse 配對','Paired pseudobulk\nedgeR: patient + condition'],['GSE218858\nMouse bulk RNA-seq','STAT5B N642H：n = 3\nControl：n = 3','DESeq2: condition\nWald statistic']],[340,390,430],155,330,24);
 txt(s,'正 NES：relapse 較高（single-cell）或 N642H 較高（bulk）',65,555,1140,70,24);
 notes(s,methods+' Sample pairs: P1 M104/M127, P2 M143/M148, P3 M187/M187r. GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262271 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218858 .');
}
{
 const s=slide('Cell type annotation');await image(s,'Fig01A_cell_lineage_UMAP',65,104,1150,557);
 notes(s,'Figure01A. UMAP of13953 cells. The six colours map to reviewed broad lineages, with the same palette used in the marker plot. T-lineage candidate12828, B459, myeloid284, pDC27, erythroid97, cytotoxic T/NK258. Coordinates and cell identities are unchanged. Broad annotation is supported by cluster markers. Candidate T lineage and CNV-high selection are separate steps. '+limits);
}
{
 const s=slide('Cell type markers');await image(s,'Fig01B_annotation_marker_dotplot',40,165,1200,435);
 notes(s,'Figure01B. Six reviewed broad lineages and the original canonical marker panel. Mean log-expression and percent expressing cells are pooled using cluster cell counts as weights, then mean expression is z-scored within gene across the six lineages. Each row label and its colour strip use the same cell type mapping as Figure01A. All20 individual cluster profiles remain in FigureS07. Exact pooled summaries are stored in annotation_lineage_marker_summary.csv; original cluster data are in annotation_marker_evidence.csv. This descriptive summary does not reassign any cells.');
}
{
 const s=slide('Potential malignant T cells');
 txt(s,'Within-sample CNV top 15%   ·   1,926 cells',65,112,1130,45,25);
 await image(s,'Fig02_paired_sample_cell_counts',85,175,1110,450);
 notes(s,'Figure02. Selected cell counts: P1 diagnosis434, relapse449; P2 diagnosis342, relapse295; P3 diagnosis107, relapse299. Each bar is one patient-sample pseudobulk. Fixed within-sample selection fractions do not estimate the prevalence of malignant cells. '+methods+' '+limits);
}
{
 const s=slide('Single-cell KEGG enrichment');
 txt(s,`${core.up_FDR005} 上調，${core.down_FDR005} 下調   (FDR < 0.05)`,65,102,1130,45,24);
 await image(s,'Fig03B_scRNAseq_slide_overview',40,155,1200,505);
 notes(s,'Figure03B. Top8 significant pathways per direction, ordered within direction by full-catalog FDR. Only2 downregulated pathways pass FDR<0.05, so both are shown. Complete353-pathway results are in Table04. Stars use unrounded full-KEGG FDR. Disease-related pathway names refer to gene sets, not diagnoses. '+methods);
}
{
 const s=slide('STAT5B N642H KEGG enrichment');
 txt(s,`${bulk.up_FDR005} 上調，${bulk.down_FDR005} 下調   (FDR < 0.05)`,65,102,1130,45,24);
 await image(s,'Fig04B_bulk_slide_overview',40,155,1200,510);
 notes(s,'Figure04B. Top8 significant pathways per direction, ordered within direction by full-catalog FDR. All6 downregulated pathways are shown. Complete300-pathway results are in Table05. Human and mouse KEGG collections differ. '+methods);
}
{
 const s=slide('同方向 KEGG pathway 交集');
 txt(s,'兩組均 FDR < 0.05，且 NES 同方向',65,112,1130,45,25);
 await image(s,'Fig05_direction_matched_intersection',65,190,1150,407);
 txt(s,'共同上調 40 條，共同下調 0 條',65,624,1130,42,25);
 notes(s,'Figure05. The overlap universe is275 common tested canonical KEGG pathway names. Up: sc-only8, shared40, bulk-only55. Down: sc-only2, shared0, bulk-only6. Circle areas are schematic, not proportional to counts. Pathways absent from one catalog are excluded from the overlap universe. Complete membership Table06;40 shared pathways Table07.');
}
{
 const s=slide('共同上調 pathways：1–20');await image(s,'Fig06B_shared_pathways_ranks01_20',40,111,1200,550);
 notes(s,'Figure06B. First20 of the complete40 shared pathways. Order is increasing max(sc_FDR,bulk_FDR), alphabetical pathway ties. This is a display sort, not a combined FDR. Proteasome is row14. X=NES, stars use unrounded full-KEGG BH FDR, all labels plain. Exact source Table07. Positive NES means higher in relapse or N642H respectively.');
}
{
 const s=slide('共同上調 pathways：21–40');await image(s,'Fig06C_shared_pathways_ranks21_40',40,111,1200,550);
 notes(s,'Figure06C. Remaining20 pathways, no omitted intersection members. Mitophagy is row28. Exact names and values are in Table07. Abbreviated labels retain the full names in the source table.');
}
{
 const s=slide('Proteostasis 與 mitochondrial pathways');await image(s,'Fig07_target_pathway_audit',45,148,1190,485);
 notes(s,'Figure07. Proteasome: sc NES1.88056919729581/FDR0.00837925505996343, bulk NES1.88731918419017/FDR0.00107780111297412. Mitophagy: sc NES1.6703493502465/FDR0.0131970176091592, bulk NES1.55764616223599/FDR0.028233910102446. General autophagy does not pass FDR0.05. Filled/open circles distinguish significant/non-significant results; dash means absent from the tested catalog. Table08 contains source values. '+limits);
}
{
 const s=slide('三位病人的配對變化');
 txt(s,'GSEA leading-edge genes 的平均表達',65,106,1130,45,25);
 await image(s,'Fig08_paired_leading_edge_direction',110,174,1060,480);
 notes(s,'Figure08. Mean pseudobulk log2CPM across pathway leading-edge genes from the same GSEA. Lines connect diagnosis and relapse for each of3 patients. Both pathway means increase in all3 pairs. The common y-axis allows comparable visual slopes. This descriptive display reuses selected genes and is not an independent test. Values in Table09.');
}
{
 const s=slide('CNV cell selection');
 table(s,[['選擇方式','門檻','保留細胞'],['主分析：within-sample top 15%','各 sample 的 CNV 第 85 百分位','1,926'],['Non-T reference q95','共同 CNV score > 0.0269','1,978']],[450,475,235],162,248,24);
 txt(s,'兩種選法重疊 1,407 顆細胞',65,456,1120,55,28);
 txt(s,'Top 15% 為敏感度分析後選定的探索型主分析',65,554,1120,60,24);
 notes(s,'Reference includes867 non-T cells. q95=0.026945532400026. Candidates12828, above q951978 (15.419%), top151926, overlap1407 (73.053% of top15). Minimum retained sample107 for top15,66 for reference-q95 selection. q95 is an empirical percentile of RNA-CNV scores, not a P value or probability of cancer. '+limits+' See TablesS02,S03,S06 and FigureS06 for sensitivity.');
}
{
 const s=slide('主要發現與後續問題');
 txt(s,'共同上調',65,148,1110,55,28,true);
 txt(s,'Oxidative phosphorylation、ribosome、proteasome、mitophagy',65,222,1110,85,28);
 txt(s,'研究假說',65,346,1110,55,28,true);
 txt(s,'STAT5B 活化與粒線體及蛋白質品質控制相關',65,420,1110,75,28);
 txt(s,'後續驗證：mitophagy flux、proteasome 活性與 BTZ 反應',65,571,1110,65,25);
 notes(s,'The shared pathways support transcriptional concordance across two different contexts and species. The mitochondrial-STAT5 mechanism is a research hypothesis. Neither dataset directly measures STAT5 mitochondrial localisation, autophagic flux, proteasome activity, bortezomib response or resistance. These functional links require new evidence. '+limits);
}
{
 const s=slide('交接與 code availability');
 table(s,[['內容','位置'],['精簡交接包','00_Quick_Start：PPT、5 張圖、3 個 CSV'],['完整圖表與方法','07_Student_Handoff'],['分析程式','05_Code_Availability：R scripts 22–27'],['Raw data / R objects','DATA_AND_OBJECT_LOCATIONS.txt']],[345,815],147,340,24);
 txt(s,'github.com/brianchiuyulab/hematologic-malignancy-transcriptomics',65,565,1130,60,23);
 notes(s,'Code22 refits count models and performs GSEA,23 generates tables,24 independently checks reproduction,25 generates supplementary evidence,26 draws the complete shared-pathway dots and27 applies the current figure style. run_reviewed_analysis.ps1 accepts project and bulk data roots. The original raw data and large R objects remain local. Figure design reference: https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/ . Scientific graphics are R-generated; PNG is embedded for presentation and vector PDF is supplied for figure editing.');
}
await fs.mkdir(TMP_DIR,{recursive:true});await fs.mkdir(path.dirname(FINAL_PPTX),{recursive:true});
const candidatePath=path.join(TMP_DIR,'candidate_reviewed.pptx');await(await PresentationFile.exportPptx(pr)).save(candidatePath);
await finalizePresentation({workspaceDir:WORKSPACE_DIR,candidatePath,finalPath:FINAL_PPTX,pythonExecutable:RUNTIME_PYTHON,integrityValidatorPath:path.join(SKILL_DIR,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(SKILL_DIR,'container_tools/inspect_presentation_layout_geometry.py'),layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit','--require-native-table-slide','2','--require-native-table-slide','13','--require-native-table-slide','15'],explicitTotalSlideCount:15,requiredNativeTableOwnerSlides:[2,13,15],requiredNativeChartOwnerSlides:[],fontPolicy:{basis:'design',families:[FONT]},verifyArtifactToolImport:true,receiptPath:path.join(TMP_DIR,'reviewed.validation.json')});
console.log(FINAL_PPTX);
