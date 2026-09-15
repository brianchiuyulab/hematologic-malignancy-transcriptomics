options(stringsAsFactors=FALSE)
.libPaths(c(Sys.getenv('GSE_REPRO_R_LIB','C:/Users/User/Documents/R/win-library/4.4'),.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(patchwork);library(ComplexHeatmap);library(circlize)})
args<-commandArgs(TRUE);src<-args[1];out<-args[2]
fig<-file.path(out,'02_Figures');tab<-file.path(out,'03_Tables');ca<-readRDS(file.path(src,'reviewed_cache.rds'));m<-ca$meta
saveplot<-function(p,n,w,h){ggsave(file.path(fig,paste0(n,'.png')),p,width=w,height=h,dpi=320,bg='white',device=ragg::agg_png);ggsave(file.path(fig,paste0(n,'.pdf')),p,width=w,height=h,device=cairo_pdf,bg='white')}
theme_set(theme_classic(base_size=11,base_family='Arial')+theme(legend.position='bottom',strip.background=element_blank(),strip.text=element_text(face='bold'),plot.title=element_text(face='bold',size=13)))
ref<-m |> filter(reviewed_cell_type%in%c('B cell','Myeloid','Erythroid','pDC'));q95<-quantile(ref$cnv_score,.95)
can<-m |> filter(reviewed_cell_type=='T-lineage candidate');cuts<-can |> group_by(sample_id,patient_pair,condition_simple) |> summarise(cutoff=quantile(cnv_score,.85),.groups='drop');cuts$sample<-paste(cuts$patient_pair,ifelse(cuts$condition_simple=='diagnosis','Dx','Rel'))
can$sample<-paste(can$patient_pair,ifelse(can$condition_simple=='diagnosis','Dx','Rel'));dd<-bind_rows(can |> transmute(sample,cnv_score,Group='T-lineage candidate'),ref |> transmute(sample='Pooled reference',cnv_score,Group='Non-T reference'))
dd$sample<-factor(dd$sample,levels=c('Pooled reference','P1 Dx','P1 Rel','P2 Dx','P2 Rel','P3 Dx','P3 Rel'))
cuts$sample<-factor(cuts$sample,levels=levels(dd$sample))
p<-ggplot(dd,aes(cnv_score,color=Group))+stat_ecdf(geom='step',linewidth=.7)+geom_vline(xintercept=q95,linetype='dashed',linewidth=.45,color='#222222')+geom_vline(data=cuts,aes(xintercept=cutoff),color='#D55E00',linewidth=.55)+facet_wrap(~sample,ncol=4)+scale_color_manual(values=c('T-lineage candidate'='#0072B2','Non-T reference'='#555555'))+labs(title='CNV-score cutoffs across samples',subtitle='Dashed: pooled non-T q95; orange: within-sample top-15% threshold',x='Mean absolute log2 inferCNV expression ratio',y='Cumulative cell fraction',color=NULL,caption='The pooled absolute cutoff and within-sample percentile are different selection rules.')+theme(plot.caption=element_text(hjust=0))
saveplot(p,'FigS02_CNV_cutoff_definition',10,5.2)
write.csv(cuts,file.path(tab,'TableS03_CNV_thresholds_by_sample.csv'),row.names=FALSE)
# Gene heatmap: preserve patient pairing and keep gene-level significance distinct
# from the gene-set FDR. Leading-edge selection is descriptive.
de<-read.csv(file.path(src,'DEG_within15.csv'));tr<-ca$results$within15 |> filter(pathway_id%in%c('hsa03050','hsa04137'))
chosen<-bind_rows(lapply(seq_len(nrow(tr)),function(i){g<-strsplit(tr$leadingEdge[i],';',fixed=TRUE)[[1]];d<-de[match(g,de$gene),];d$Pathway<-if(tr$pathway_id[i]=='hsa03050')'Proteasome'else'Mitophagy';d |> arrange(desc(rank_stat)) |> slice_head(n=12)}))
allchosen<-bind_rows(lapply(seq_len(nrow(tr)),function(i){g<-strsplit(tr$leadingEdge[i],';',fixed=TRUE)[[1]];d<-de[match(g,de$gene),];d$Pathway<-if(tr$pathway_id[i]=='hsa03050')'Proteasome'else'Mitophagy';d}))
write.csv(allchosen |> select(Pathway,gene,logFC,PValue,FDR,rank_stat),file.path(tab,'TableS04_leading_edge_gene_statistics.csv'),row.names=FALSE)
sm<-m |> distinct(sample_id,patient_pair,condition_simple) |> arrange(patient_pair,condition_simple)
mat<-ca$logcpm$within15[chosen$gene,sm$sample_id];mat<-t(scale(t(mat)));mat[mat < -2]<- -2;mat[mat>2]<-2;rownames(mat)<-paste0(chosen$gene,ifelse(chosen$FDR<.05,' *',''));colnames(mat)<-ifelse(sm$condition_simple=='diagnosis','Dx','Rel')
ha<-HeatmapAnnotation(Condition=ifelse(sm$condition_simple=='diagnosis','Diagnosis','Relapse'),col=list(Condition=c(Diagnosis='#0072B2',Relapse='#D55E00')),show_annotation_name=FALSE,annotation_legend_param=list(Condition=list(nrow=1)))
ht<-Heatmap(mat,name='Row z-score',col=colorRamp2(c(-2,0,2),c('#2166AC','#FFFFFF','#B2182B')),cluster_rows=FALSE,cluster_columns=FALSE,row_split=factor(chosen$Pathway,levels=c('Mitophagy','Proteasome')),cluster_row_slices=FALSE,column_split=sm$patient_pair,column_title=c('P1','P2','P3'),top_annotation=ha,row_names_gp=grid::gpar(fontsize=10,fontface='italic'),column_names_gp=grid::gpar(fontsize=10),column_names_centered=TRUE,row_title_rot=0,row_title_gp=grid::gpar(fontsize=11,fontface='bold'),column_names_rot=0,row_names_side='left',rect_gp=grid::gpar(col='white',lwd=.7),row_gap=grid::unit(3,'mm'),column_gap=grid::unit(2,'mm'),heatmap_legend_param=list(direction='horizontal',legend_width=grid::unit(35,'mm')))
drawheat<-function(){draw(ht,heatmap_legend_side='bottom',annotation_legend_side='bottom',merge_legends=TRUE,padding=grid::unit(c(15,3,17,3),'mm'));grid::grid.text('GSEA leading-edge genes: within-sample top 15%',x=.02,y=.975,just='left',gp=grid::gpar(fontsize=13,fontface='bold'));grid::grid.text('Top 12 per pathway by paired rank; * gene-level BH FDR < 0.05',x=.02,y=.948,just='left',gp=grid::gpar(fontsize=10));grid::grid.text('Pathway FDR is reported separately in Fig07. Unmarked genes do not pass gene-level FDR < 0.05.',x=.02,y=.02,just='left',gp=grid::gpar(fontsize=9))}
ragg::agg_png(file.path(fig,'FigS03_paired_gene_heatmap.png'),width=8.2,height=8.5,units='in',res=320);drawheat();dev.off();cairo_pdf(file.path(fig,'FigS03_paired_gene_heatmap.pdf'),width=8.2,height=8.5);drawheat();dev.off()
write.csv(data.frame(genes_in_heatmap=nrow(chosen),genes_FDR_lt005=sum(chosen$FDR<.05),all_leading_edge_genes=nrow(allchosen)),file.path(src,'gene_heatmap_statistics.csv'),row.names=FALSE)
cat('Supporting figures complete\n')
