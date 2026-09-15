options(stringsAsFactors=FALSE, timeout=300)
try(Sys.setlocale('LC_CTYPE','English_United States.utf8'),silent=TRUE)
.libPaths(unique(c(Sys.getenv('GSE_REPRO_R_LIB', 'C:/Users/User/Documents/R/win-library/4.4'),.libPaths())))
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(edgeR);library(fgsea);library(dplyr);library(tidyr);library(AnnotationDbi);library(org.Hs.eg.db);library(org.Mm.eg.db)})
select<-dplyr::select
stopifnot(as.character(packageVersion('edgeR'))=='4.4.2',as.character(packageVersion('fgsea'))=='1.32.4')
reuse<-identical(Sys.getenv('AUDIT_RESUME'),'1')
args<-commandArgs(TRUE)
if(length(args)==0) args<-c('C:/Users/User/Desktop/Public dataset analysis/GEO dataset/2026/single cell/BTZ \u76f8\u95dc/GSE262271_TALL_relapse_scRNAseq','C:/Users/User/Desktop/Public dataset analysis/GEO dataset/2026/bulk/\u8840\u764c\u76f8\u95dc/R version','C:/Users/User/Documents/BTZ resistance/outputs/reviewer_audit/results')
project<-normalizePath(args[1],winslash='/',mustWork=TRUE)
bulk<-normalizePath(args[2],winslash='/',mustWork=TRUE)
out<-normalizePath(args[3],winslash='/',mustWork=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
wr<-function(x,name) write.csv(as.data.frame(x),file.path(out,paste0(name,'.csv')),row.names=FALSE,na='',fileEncoding='UTF-8')
msg<-function(...) {cat(format(Sys.time()),..., '\n');flush.console()}
cnv<-read.csv(file.path(project,'02_Tables/04_Method_Inputs/infercnv_nonTref_cell_cnv_burden.csv'))
kegg<-read.csv(file.path(project,'02_Tables/02_Full_Results/KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv'))
kegg<-kegg[!duplicated(kegg$pathway)&nzchar(kegg$linked_genes),]
sets<-setNames(lapply(strsplit(kegg$linked_genes,';',fixed=TRUE),unique),kegg$pathway)
ann<-kegg[,c('pathway','pathway_id','pathway_name')]
targets<-c('hsa03050','hsa04137','hsa04140','hsa04136')
msg('Loading final Seurat object')
obj<-readRDS(file.path(project,'04_R_Objects/GSE262271_seurat_qc_harmony_annotated.rds'))
obj<-JoinLayers(obj,assay='RNA')
counts<-GetAssayData(obj,assay='RNA',layer='counts')
norm<-GetAssayData(obj,assay='RNA',layer='data')
meta<-obj@meta.data;meta$cell<-rownames(meta)
stopifnot(!anyDuplicated(meta$cell),!anyDuplicated(cnv$cell),identical(colnames(counts),meta$cell))
meta$cnv_score<-cnv$infercnv_mean_abs_log2[match(meta$cell,cnv$cell)]
meta$patient_pair<-factor(meta$patient_pair,levels=c('P1','P2','P3'))
meta$condition_simple<-factor(meta$condition_simple,levels=c('diagnosis','relapse'))
type_map<-c(T_lineage_candidate_malignant='T-lineage candidate',B_cell_reference='B cell',Myeloid_reference='Myeloid',Erythroid_reference='Erythroid',pDC_reference='pDC',Cytotoxic_T_or_NK_reference_like='Cytotoxic T/NK')
meta$reviewed_cell_type<-unname(type_map[meta$broad_annotation])
stopifnot(!anyNA(meta$reviewed_cell_type))
umap<-as.data.frame(Embeddings(obj,'umap.harmony'));names(umap)<-c('UMAP1','UMAP2');umap$cell<-rownames(umap)
wr(left_join(umap,meta[,c('cell','sample_id','patient_pair','condition_simple','cluster_res0.6','reviewed_cell_type','cnv_score')],by='cell'),'annotation_umap')
wr(meta |> group_by(sample_id,patient_pair,condition_simple,reviewed_cell_type) |> summarise(cells=n(),median_features=median(nFeature_RNA),median_UMI=median(nCount_RNA),median_mt=median(percent.mt),.groups='drop'),'annotation_sample_counts')
wr(meta |> group_by(cluster_res0.6,reviewed_cell_type,major_annotation) |> summarise(cells=n(),samples=n_distinct(sample_id),median_features=median(nFeature_RNA),median_UMI=median(nCount_RNA),median_mt=median(percent.mt),.groups='drop'),'annotation_cluster_review')
panels<-list('T lineage'=c('CD3D','CD3E','TRAC','CD5','CD7'),'Immaturity'=c('DNTT','PTCRA','CD34','SOX4'),'B cell'=c('MS4A1','CD79A','CD79B'),'Myeloid'=c('LYZ','LST1','FCN1'),'pDC'=c('CLEC4C','LILRA4','TCF4'),'Erythroid'=c('HBB','ALAS2','AHSP'),'Cytotoxic'=c('NKG7','GNLY','KLRD1'),'Cell cycle'=c('MKI67','TOP2A','UBE2C'))
markers<-list();diffs<-list()
if(!reuse || !file.exists(file.path(out,'annotation_top30_markers_descriptive.csv'))){
for(cl in as.character(0:19)){
  ii<-which(meta$cluster_res0.6==cl);jj<-which(meta$cluster_res0.6!=cl)
  for(panel in names(panels))for(g in intersect(panels[[panel]],rownames(norm)))markers[[length(markers)+1]]<-data.frame(cluster=cl,cell_type=meta$reviewed_cell_type[ii[1]],panel=panel,gene=g,mean_log_expression=mean(norm[g,ii]),pct_expressing=100*mean(counts[g,ii]>0))
  avg<-Matrix::rowMeans(expm1(norm[,ii,drop=FALSE]));avg_other<-Matrix::rowMeans(expm1(norm[,jj,drop=FALSE]));pct<-Matrix::rowMeans(counts[,ii,drop=FALSE]>0);pct2<-Matrix::rowMeans(counts[,jj,drop=FALSE]>0)
  d<-data.frame(cluster=cl,gene=rownames(norm),descriptive_log2_ratio=log2((avg+1)/(avg_other+1)),pct_in=100*pct,pct_other=100*pct2)
  diffs[[cl]]<-d |> filter(pct_in>=25) |> arrange(desc(descriptive_log2_ratio)) |> head(30)
}
wr(bind_rows(markers),'annotation_marker_evidence');wr(bind_rows(diffs),'annotation_top30_markers_descriptive')
}
wr(meta |> count(HTO_classification.global),'HTO_doublet_QC')
candidate<-meta |> filter(broad_annotation=='T_lineage_candidate_malignant',is.finite(cnv_score))
reference<-cnv |> filter(grepl('^ref_',infercnv_group),is.finite(infercnv_mean_abs_log2))
q95<-unname(quantile(reference$infercnv_mean_abs_log2,.95))
wr(reference |> group_by(sample_id,broad_annotation) |> summarise(cells=n(),median=median(infercnv_mean_abs_log2),q95=quantile(infercnv_mean_abs_log2,.95),.groups='drop'),'reference_by_sample_celltype')
wr(data.frame(reference_cells=nrow(reference),q95=q95,candidate_cells=nrow(candidate),above_q95=sum(candidate$cnv_score>q95),fraction_above=mean(candidate$cnv_score>q95)),'reference_threshold')
selections<-list(global25=candidate |> filter(cnv_score>=quantile(cnv_score,.75)),within15=candidate |> group_by(sample_id) |> filter(cnv_score>=quantile(cnv_score,.85)) |> ungroup(),within25=candidate |> group_by(sample_id) |> filter(cnv_score>=quantile(cnv_score,.75)) |> ungroup(),reference95=candidate |> filter(cnv_score>q95))
sample_info<-candidate |> distinct(sample_id,patient_pair,condition_simple) |> arrange(patient_pair,condition_simple)
design<-model.matrix(~patient_pair+condition_simple,sample_info)
wr(cbind(sample_info,design),'paired_design_matrix')
selection_counts<-bind_rows(lapply(names(selections),function(id) sample_info |> left_join(selections[[id]] |> count(sample_id,name='cells'),by='sample_id') |> mutate(analysis=id,cells=replace_na(cells,0))))
wr(selection_counts,'selection_cell_counts')
wr(bind_rows(lapply(names(selections),function(id)selections[[id]] |> count(cluster_res0.6,sample_id) |> mutate(analysis=id))),'selection_cluster_composition')
quality<-list();ranks<-list();logcpms<-list();degs<-list()
for(id in names(selections)){
  sel<-selections[[id]];n<-selection_counts$cells[selection_counts$analysis==id]
  if(any(n<20)){msg(id,'not fit because <20 cells in a paired sample');next}
  msg('Paired edgeR',id)
  map<-sparseMatrix(i=seq_len(nrow(sel)),j=match(sel$sample_id,sample_info$sample_id),x=1,dims=c(nrow(sel),nrow(sample_info)))
  pb<-as.matrix(counts[,sel$cell,drop=FALSE]%*%map);colnames(pb)<-sample_info$sample_id
  y<-DGEList(pb);keep<-filterByExpr(y,design=design);y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE);test<-glmQLFTest(fit,coef='condition_simplerelapse')
  d<-topTags(test,n=Inf,sort.by='none')$table;d$gene<-rownames(d);d$rank_stat<-sign(d$logFC)*sqrt(pmax(d$F,0));d$analysis<-id
  stopifnot(max(abs(d$FDR-p.adjust(d$PValue,'BH')))<1e-12)
  r<-setNames(d$rank_stat,d$gene);r<-sort(r[is.finite(r)],decreasing=TRUE)
  ranks[[id]]<-r;degs[[id]]<-d;logcpms[[id]]<-cpm(y,log=TRUE,prior.count=2)
  quality[[id]]<-data.frame(analysis=id,cells=nrow(sel),minimum_cells=min(n),genes=length(r),residual_df=nrow(design)-qr(design)$rank)
  wr(d,paste0('DEG_',id))
}
wr(bind_rows(quality),'analysis_registry')
old<-read.csv(file.path(project,'02_Tables/02_Full_Results/pseudobulk_DEG_relapse_vs_diagnosis_all_genes.csv'))
cmp<-inner_join(old,degs$global25,by='gene',suffix=c('_old','_reviewed'))
wr(data.frame(genes_old=nrow(old),genes_reviewed=nrow(degs$global25),shared=nrow(cmp),max_abs_logFC=max(abs(cmp$logFC_old-cmp$logFC_reviewed)),max_abs_F=max(abs(cmp$F_old-cmp$F_reviewed)),max_abs_P=max(abs(cmp$PValue_old-cmp$PValue_reviewed))),'primary_DEG_reproduction')
fg_run<-function(r,gs,label,maxsize=Inf){
  cached<-file.path(out,paste0('KEGG_',if(label=='bulk_STAT5B')'bulk' else label,'.csv'))
  if(reuse && identical(Sys.getenv('AUDIT_GSEA_REUSE'),'1') && file.exists(cached)) {msg('Resuming completed GSEA',label);return(read.csv(cached))}
  n<-lengths(lapply(gs,intersect,names(r)));gs<-gs[n>=5 & n<=maxsize]
  msg('GSEA',label,length(gs),'sets')
  set.seed(20260915)
  z<-as.data.frame(fgseaMultilevel(gs,r,minSize=5,maxSize=if(is.finite(maxsize))maxsize else max(n),eps=0,sampleSize=501,nPermSimple=50000,nproc=1))
  stopifnot(all(is.finite(z$pval)),max(abs(z$padj-p.adjust(z$pval,'BH')))<1e-12)
  z$leadingEdge<-vapply(z$leadingEdge,paste,collapse=';',FUN.VALUE=character(1));z$analysis<-label;z$direction<-ifelse(z$NES>0,'Up','Down');z$tested_sets<-nrow(z);z
}
# First reproduce the exact frozen gene-set membership. A separate namespace
# sensitivity resolves only uniquely annotated aliases plus known MT symbols.
allgenes<-unique(unlist(sets));valid<-keys(org.Hs.eg.db,keytype='SYMBOL');unknown<-setdiff(allgenes,valid)
aliases<-suppressMessages(AnnotationDbi::select(org.Hs.eg.db,keys=intersect(unknown,keys(org.Hs.eg.db,keytype='ALIAS')),keytype='ALIAS',columns='SYMBOL')) |> filter(!is.na(SYMBOL)) |> distinct(ALIAS,SYMBOL) |> group_by(ALIAS) |> filter(n()==1) |> ungroup()
mapping<-setNames(allgenes,allgenes);mapping[aliases$ALIAS]<-aliases$SYMBOL
mt_map<-c(ATP6='MT-ATP6',ATP8='MT-ATP8',COX1='MT-CO1',COX2='MT-CO2',COX3='MT-CO3',CYTB='MT-CYB',ND1='MT-ND1',ND2='MT-ND2',ND3='MT-ND3',ND4='MT-ND4',ND4L='MT-ND4L',ND5='MT-ND5',ND6='MT-ND6');mapping[intersect(names(mt_map),allgenes)]<-mt_map[intersect(names(mt_map),allgenes)]
wr(data.frame(original=names(mapping),reviewed=unname(mapping),changed=names(mapping)!=mapping,reviewed_symbol_valid=unname(mapping)%in%valid),'gene_symbol_mapping_audit')
resolved_sets<-lapply(sets,function(g)unique(unname(mapping[g])))
res<-list()
for(id in names(ranks)){
  z<-fg_run(ranks[[id]],sets,id);if(!'pathway_id'%in%names(z))z<-left_join(z,ann,by='pathway');wr(z,paste0('KEGG_',id));res[[id]]<-z
}
resolved<-fg_run(ranks$within15,resolved_sets,'within15_symbol_sensitivity');if(!'pathway_id'%in%names(resolved))resolved<-left_join(resolved,ann,by='pathway');wr(resolved,'KEGG_within15_symbol_sensitivity')
allsc<-bind_rows(res);wr(allsc,'KEGG_single_cell_all_definitions')
wr(allsc |> filter(pathway_id%in%targets),'KEGG_target_summary')
pair<-bind_rows(lapply(names(res),function(id){z<-res[[id]] |> filter(pathway_id%in%targets);bind_rows(lapply(seq_len(nrow(z)),function(i){g<-intersect(strsplit(z$leadingEdge[i],';',fixed=TRUE)[[1]],rownames(logcpms[[id]]));scores<-colMeans(logcpms[[id]][g,,drop=FALSE]);s<-sample_info;s$score<-scores[s$sample_id];s$analysis<-id;s$pathway<-z$pathway_name[i];s |> select(analysis,pathway,patient_pair,condition_simple,score) |> pivot_wider(names_from=condition_simple,values_from=score) |> mutate(delta=relapse-diagnosis)}))}))
wr(pair,'target_paired_direction')
# Bulk data: re-fit six deposited count libraries using the original condition
# contrast, then compare gene-level output before recomputing its whole KEGG set.
msg('Re-fitting bulk DESeq2')
suppressPackageStartupMessages(library(DESeq2))
bm<-read.csv(file.path(bulk,'processed/R_sample_annotation.csv'),check.names=FALSE) |> filter(condition%in%c('R2b_DN','R2S5b_DN'))
bc<-data.table::fread(file.path(bulk,'processed/R_raw_count_matrix.csv.gz'),select=c('ensembl_id',bm$sample_id))
bc<-as.data.frame(bc);rownames(bc)<-bc$ensembl_id;bc$ensembl_id<-NULL;bc<-as.matrix(bc[,bm$sample_id]);storage.mode(bc)<-'integer'
bc<-bc[rowSums(bc)>=10 & rowSums(bc>=5)>=2,];rownames(bm)<-bm$sample_id;bm$condition<-factor(bm$condition,levels=c('R2b_DN','R2S5b_DN'))
dds<-DESeqDataSetFromMatrix(bc,bm,~condition);dds<-DESeq(dds,quiet=TRUE);bd<-as.data.frame(results(dds,contrast=c('condition','R2S5b_DN','R2b_DN'),alpha=.05));bd$ensembl_id<-rownames(bd)
oldb<-read.csv(file.path(bulk,'results/de_results/R_DESeq2_primary_R2S5b_DN_vs_R2b_DN.csv'))
cb<-inner_join(bd,oldb,by='ensembl_id',suffix=c('_reviewed','_old'))
wr(data.frame(genes=nrow(cb),max_abs_logFC=max(abs(cb$log2FoldChange_reviewed-cb$log2FoldChange_old),na.rm=TRUE),max_abs_stat=max(abs(cb$stat_reviewed-cb$stat_old),na.rm=TRUE),max_abs_P=max(abs(cb$pvalue_reviewed-cb$pvalue_old),na.rm=TRUE)),'bulk_DEG_reproduction')
ba<-read.csv(file.path(bulk,'processed/R_gene_annotation_orgMm.csv'))
bd<-left_join(bd,ba,by='ensembl_id') |> mutate(rank_metric=ifelse(is.na(stat),sign(log2FoldChange)*-log10(pmax(pvalue,1e-300)),stat))
br<-bd |> filter(!is.na(gene_symbol),gene_symbol!='',is.finite(rank_metric)) |> group_by(gene_symbol) |> slice_max(abs(rank_metric),n=1,with_ties=FALSE) |> ungroup() |> arrange(desc(rank_metric))
brank<-setNames(br$rank_metric,br$gene_symbol)
ln<-strsplit(readLines(file.path(bulk,'reference/KEGG_2019_Mouse.enrichr.gmt'),warn=FALSE),'\t',fixed=TRUE);bsets<-setNames(lapply(ln,function(x)unique(x[-c(1,2)])),vapply(ln,'[',character(1),1))
bg<-unique(unlist(bsets));upper<-setNames(names(brank)[!duplicated(toupper(names(brank)))],toupper(names(brank)[!duplicated(toupper(names(brank)))]));bmap<-setNames(names(brank)[match(bg,names(brank))],bg);missing<-is.na(bmap);bmap[missing]<-upper[toupper(names(bmap)[missing])]
ba2<-suppressMessages(AnnotationDbi::select(org.Mm.eg.db,keys=intersect(names(bmap)[is.na(bmap)],keys(org.Mm.eg.db,keytype='ALIAS')),keytype='ALIAS',columns='SYMBOL')) |> filter(!is.na(SYMBOL),SYMBOL%in%names(brank))
ba2<-ba2[!duplicated(ba2$ALIAS),];bmap[ba2$ALIAS]<-ba2$SYMBOL
bsets<-lapply(bsets,function(g)unique(na.omit(unname(bmap[g]))));bgsea<-fg_run(brank,bsets,'bulk_STAT5B',2000);bgsea$pathway_name<-bgsea$pathway;wr(bgsea,'KEGG_bulk');wr(bm[,c('sample_id','condition','genotype','stage','name.stage_replicate')],'bulk_sample_design')
wr(data.frame(pathway=names(bsets),genes=vapply(bsets,paste,collapse=';',FUN.VALUE=character(1))),'bulk_KEGG_sets_used')
saveRDS(list(results=res,bulk=bgsea,pair=pair,meta=meta,selections=lapply(selections,function(x)x$cell),logcpm=logcpms,sets=sets,ann=ann),file.path(out,'reviewed_cache.rds'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
msg('Complete');print(allsc |> filter(pathway_id%in%targets) |> select(analysis,pathway_name,NES,pval,padj));print(bgsea |> filter(grepl('Mitophagy|Proteasome|Autophagy',pathway)) |> select(pathway,NES,pval,padj))
