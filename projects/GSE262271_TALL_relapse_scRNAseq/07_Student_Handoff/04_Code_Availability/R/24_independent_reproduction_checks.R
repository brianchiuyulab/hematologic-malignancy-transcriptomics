options(stringsAsFactors=FALSE)
.libPaths(c(Sys.getenv('GSE_REPRO_R_LIB','C:/Users/User/Documents/R/win-library/4.4'),.libPaths()))
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(edgeR);library(fgsea);library(dplyr)})
args<-commandArgs(TRUE);project<-args[1];src<-args[2]
wr<-function(x,n)write.csv(x,file.path(src,paste0(n,'.csv')),row.names=FALSE)
ca<-readRDS(file.path(src,'reviewed_cache.rds'));obj<-JoinLayers(readRDS(file.path(project,'04_R_Objects/GSE262271_seurat_qc_harmony_annotated.rds')),assay='RNA');counts<-GetAssayData(obj,assay='RNA',layer='counts');meta<-ca$meta
sm<-meta |> distinct(sample_id,patient_pair,condition_simple) |> arrange(patient_pair,condition_simple);design<-model.matrix(~patient_pair+condition_simple,sm)
cells<-ca$selections$global25;cm<-meta[match(cells,meta$cell),]
pb<-do.call(cbind,lapply(sm$sample_id,function(s)Matrix::rowSums(counts[,cm$cell[cm$sample_id==s],drop=FALSE])));colnames(pb)<-sm$sample_id
map<-sparseMatrix(i=seq_along(cells),j=match(cm$sample_id,sm$sample_id),x=1,dims=c(length(cells),nrow(sm)))
matrix_difference<-max(abs(pb-as.matrix(counts[,cells,drop=FALSE]%*%map)))
fit_model<-function(x){y<-DGEList(counts=x,samples=as.data.frame(sm));keep<-filterByExpr(y,design=design);y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE);test<-glmQLFTest(fit,coef='condition_simplerelapse');z<-topTags(test,n=Inf,sort.by='none')$table;z$gene<-rownames(z);z}
z<-fit_model(as(pb,'dgCMatrix'));new<-read.csv(file.path(src,'DEG_global25.csv'));old<-read.csv(file.path(project,'02_Tables/02_Full_Results/pseudobulk_DEG_relapse_vs_diagnosis_all_genes.csv'))
a<-merge(z,new,by='gene');b<-merge(z,old,by='gene')
wr(data.frame(check=c('independent_aggregation_max_abs_count_difference','original_code_style_vs_reviewed_max_abs_F','original_code_style_vs_reviewed_max_abs_logFC','old_stored_table_vs_reviewed_max_abs_F','old_stored_table_vs_reviewed_max_abs_logCPM'),value=c(matrix_difference,max(abs(a$F.x-a$F.y)),max(abs(a$logFC.x-a$logFC.y)),max(abs(b$F.x-b$F.y)),max(abs(b$logCPM.x-b$logCPM.y)))),'independent_reproduction_checks')
# Recompute GSEA on the historical stored ranking to separate rank changes from
# Monte Carlo precision. Never replace a result by the minimum across seeds.
r<-sort(setNames(old$rank_stat,old$gene),decreasing=TRUE);set.seed(20260915)
g<-as.data.frame(fgseaMultilevel(ca$sets,r,minSize=5,maxSize=10000,eps=0,sampleSize=501,nPermSimple=50000,nproc=1));g$leadingEdge<-vapply(g$leadingEdge,paste,collapse=';',FUN.VALUE=character(1));g<-left_join(g,ca$ann,by='pathway');g$analysis<-'historical_rank_high_precision';wr(g,'KEGG_historical_rank_high_precision')
h<-read.csv(file.path(project,'02_Tables/02_Full_Results/KEGG_GSEA_all_pathways_relapse_vs_diagnosis.csv'));targets<-c('hsa03050','hsa04137');audit<-bind_rows(h |> filter(pathway_id%in%targets) |> transmute(source='Stored original figure',pathway_name,NES,pval,padj),g |> filter(pathway_id%in%targets) |> transmute(source='Stored ranking; high-precision GSEA',pathway_name,NES,pval,padj),ca$results$global25 |> filter(pathway_id%in%targets) |> transmute(source='Refit model; high-precision GSEA',pathway_name,NES,pval,padj));wr(audit,'historical_vs_reviewed_target_audit')
write.csv(meta[,c('cell','sample_id','patient_pair','condition_simple','cluster_res0.6','reviewed_cell_type','cnv_score')],file.path(src,'reviewed_cell_metadata.csv'),row.names=FALSE)
saveRDS(list(metadata=meta[,c('cell','sample_id','patient_pair','condition_simple','cluster_res0.6','reviewed_cell_type','cnv_score')],cell_selections=ca$selections,notes='Broad marker-reviewed lineage labels, not validated malignant classification.'),file.path(src,'reviewed_metadata_and_selections.rds'))
obj$reviewed_cell_type<-meta$reviewed_cell_type[match(colnames(obj),meta$cell)]
obj$reviewed_infercnv_score<-meta$cnv_score[match(colnames(obj),meta$cell)]
obj$potential_malignant_global25<-colnames(obj)%in%ca$selections$global25
obj$potential_malignant_within15_exploratory<-colnames(obj)%in%ca$selections$within15
Idents(obj)<-'reviewed_cell_type'
saveRDS(obj,file.path(src,'GSE262271_reviewed_lineage_annotation.rds'))
print(audit)
