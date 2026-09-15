# Reader-facing scientific graphics. Reads frozen audit results; never refits models.
options(stringsAsFactors=FALSE)
.libPaths(c(Sys.getenv('GSE_REPRO_R_LIB','C:/Users/User/Documents/R/win-library/4.4'),.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(patchwork);library(ggrepel)})
args<-commandArgs(TRUE);stopifnot(length(args)==2)
src<-args[1];out<-args[2];fig<-file.path(out,'02_Figures');dir.create(fig,recursive=TRUE,showWarnings=FALSE)
rd<-function(n)read.csv(file.path(src,paste0(n,'.csv')),check.names=FALSE)
savep<-function(p,n,w=7.2,h=4.5){
 if(is.data.frame(p$data) && 'NES'%in%names(p$data)){
  b<-ggplot_build(p);i<-which(vapply(p$layers,function(l)inherits(l$geom,'GeomPoint'),logical(1)))[1]
  stopifnot(!is.na(i),isTRUE(all.equal(b$data[[i]]$x,p$data$NES,tolerance=0)),!anyNA(b$data[[i]]$x))
  if('star'%in%names(p$data))stopifnot(identical(p$data$star,stars(p$data$padj)))
 }
 ggsave(file.path(fig,paste0(n,'.png')),p,width=w,height=h,dpi=450,device=grDevices::png,type='cairo',bg='white');ggsave(file.path(fig,paste0(n,'.pdf')),p,width=w,height=h,device=cairo_pdf,bg='white')
}
base<-theme_classic(base_family='Arial',base_size=11)+theme(axis.text=element_text(color='#222222'),axis.title=element_text(size=11),axis.line=element_line(linewidth=.4),axis.ticks=element_line(linewidth=.3),strip.background=element_blank(),strip.text=element_text(size=11,face='plain'),legend.position='bottom',legend.title=element_text(size=10),legend.text=element_text(size=10),plot.title=element_text(size=12,face='plain'),plot.caption=element_text(size=9,hjust=0),plot.margin=margin(8,10,8,8),panel.grid=element_blank())
theme_set(base)
cols<-c('T-lineage candidate'='#79818C','B cell'='#0072B2','Myeloid'='#E69F00','pDC'='#CC79A7','Erythroid'='#D55E00','Cytotoxic T/NK'='#009E73')
datasets<-c('Relapse scRNA-seq'='#0072B2','STAT5B N642H bulk'='#D55E00')
stars<-function(x)ifelse(x<.001,'***',ifelse(x<.01,'**',ifelse(x<.05,'*','')))
keycaption<-'* FDR < 0.05     ** FDR < 0.01     *** FDR < 0.001'
canon<-function(s)sub(' - animal$','',s)
label<-function(s){s<-canon(s);s<-sub(' signaling pathway$',' signaling',s);s<-sub('Non-alcoholic fatty liver disease \\(NAFLD\\)','Non-alcoholic fatty liver disease',s);s<-sub('Cell adhesion molecules \\(CAMs\\)','Cell adhesion molecules',s);s<-sub('Chemical carcinogenesis - reactive oxygen species','Chemical carcinogenesis (ROS)',s);s<-sub('Pathways of neurodegeneration - multiple diseases','Neurodegeneration pathways',s);s}
u<-rd('annotation_umap');cl<-rd('annotation_cluster_review')
stopifnot(nrow(u)==13953,sum(cl$cells)==13953,setequal(unique(u$reviewed_cell_type),names(cols)))
u$reviewed_cell_type<-factor(u$reviewed_cell_type,levels=names(cols))
cent<-u |> filter(reviewed_cell_type!='T-lineage candidate') |> group_by(reviewed_cell_type) |> summarise(UMAP1=median(UMAP1),UMAP2=median(UMAP2),.groups='drop')
labelpos<-data.frame(reviewed_cell_type=c('B cell','Myeloid','pDC','Erythroid','Cytotoxic T/NK'),lx=c(1.4,.3,11.6,11.8,3.2),ly=c(-13.3,-17.6,-7.7,-10.8,-2.1))
cent<-left_join(cent,labelpos,by='reviewed_cell_type')
p<-ggplot(u,aes(UMAP1,UMAP2))+geom_point(aes(color=reviewed_cell_type),size=.38,alpha=.85,stroke=0)+
 geom_segment(data=cent,aes(x=lx,y=ly,xend=UMAP1,yend=UMAP2),color='#555555',linewidth=.25)+
 geom_label(data=cent,aes(x=lx,y=ly,label=reviewed_cell_type),size=3.3,color='#222222',fill='white',linewidth=0,label.padding=grid::unit(.12,'lines'))+
 scale_color_manual(values=cols,name='Cell type',guide=guide_legend(override.aes=list(size=3.5,alpha=1)))+coord_equal()+labs(x='UMAP 1',y='UMAP 2')+theme(legend.position='right',legend.key.height=grid::unit(7,'mm'))
savep(p,'Fig01A_cell_lineage_UMAP',8.8,5.7)
# Marker evidence: row colour strip links every cluster to the same UMAP key.
mk<-rd('annotation_marker_evidence');genes<-unique(mk$gene)
order_cl<-cl |> mutate(type=factor(reviewed_cell_type,levels=names(cols)),id=as.integer(cluster_res0.6)) |> arrange(type,id)
order_cl$y<-rev(seq_len(nrow(order_cl)))
mk$y<-order_cl$y[match(as.integer(mk$cluster),order_cl$id)];mk$x<-match(mk$gene,genes)
mk<-mk |> group_by(gene) |> mutate(z=as.numeric(scale(mean_log_expression))) |> ungroup()
pan<-mk |> group_by(panel) |> summarise(x=mean(range(x)),.groups='drop')
order_cl$type<-factor(order_cl$type,levels=names(cols))
p<-ggplot(mk,aes(x,y))+geom_tile(data=order_cl,aes(x=0,y=y,fill=type),inherit.aes=FALSE,width=.32,height=.86)+geom_point(aes(size=pct_expressing,color=z))+
 scale_fill_manual(values=cols,name='Cell type',guide=guide_legend(order=3,nrow=2,override.aes=list(alpha=1)))+
 scale_color_gradientn(colors=c('#E4E7EA','#6F95B3','#102C4B'),name='Mean expression\n(gene z-score)',guide=guide_colorbar(order=1,barwidth=grid::unit(25,'mm'),barheight=grid::unit(3,'mm')))+
 scale_size_area(max_size=4.5,breaks=c(25,50,75,100),name='Cells (%)',guide=guide_legend(order=2))+
 geom_text(data=pan,aes(x=x,y=21.2,label=panel),inherit.aes=FALSE,size=3.1)+
 scale_x_continuous(breaks=seq_along(genes),labels=genes,limits=c(-.3,length(genes)+.5),expand=expansion(mult=0))+
 scale_y_continuous(breaks=order_cl$y,labels=order_cl$id,limits=c(.3,22),expand=expansion(mult=0))+
 labs(x=NULL,y='Cluster')+theme(axis.text.x=element_text(angle=55,hjust=1,face='italic',size=9),axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.box='vertical',legend.spacing.y=grid::unit(1,'mm'))
savep(p,'FigS07_cluster_marker_dotplot',11.5,6.7)
# Pool cluster means with cell-count weights; retain the complete cluster figure above.
mk$cells<-cl$cells[match(as.integer(mk$cluster),as.integer(cl$cluster_res0.6))]
stopifnot(!anyNA(mk$cells),all(mk$cell_type==cl$reviewed_cell_type[match(as.integer(mk$cluster),as.integer(cl$cluster_res0.6))]))
agg<-mk |> group_by(cell_type,panel,gene) |> summarise(mean_log_expression=weighted.mean(mean_log_expression,cells),pct_expressing=weighted.mean(pct_expressing,cells),cells=sum(cells),.groups='drop') |> group_by(gene) |> mutate(z=as.numeric(scale(mean_log_expression))) |> ungroup()
agg$x<-match(agg$gene,genes);agg$y<-7-match(agg$cell_type,names(cols));agg$cell_type<-factor(agg$cell_type,levels=names(cols))
stopifnot(nrow(agg)==6*length(genes),all(agg$pct_expressing>=0 & agg$pct_expressing<=100))
rows<-data.frame(cell_type=factor(names(cols),levels=names(cols)),y=6:1)
p<-ggplot(agg,aes(x,y))+geom_tile(data=rows,aes(x=0,y=y,fill=cell_type),inherit.aes=FALSE,width=.32,height=.72)+geom_point(aes(size=pct_expressing,color=z))+
 scale_fill_manual(values=cols,guide='none')+scale_color_gradientn(colors=c('#E4E7EA','#6F95B3','#102C4B'),name='Mean expression\n(gene z-score)',guide=guide_colorbar(order=1,barwidth=grid::unit(27,'mm'),barheight=grid::unit(3,'mm')))+
 scale_size_area(max_size=5.6,breaks=c(25,50,75,100),name='Cells (%)',guide=guide_legend(order=2))+
 geom_text(data=pan,aes(x=x,y=7,label=panel),inherit.aes=FALSE,size=3.2)+
 scale_x_continuous(breaks=seq_along(genes),labels=genes,limits=c(-.3,length(genes)+.5),expand=expansion(mult=0))+
 scale_y_continuous(breaks=6:1,labels=names(cols),limits=c(.3,7.5),expand=expansion(mult=0))+labs(x=NULL,y=NULL)+
 theme(axis.text.x=element_text(angle=55,hjust=1,face='italic',size=10),axis.text.y=element_text(size=11),axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.box='horizontal')
savep(p,'Fig01B_annotation_marker_dotplot',11.5,4.2)
write.csv(agg |> select(cell_type,panel,gene,mean_log_expression,pct_expressing,cells,z),file.path(out,'04_Code_Availability','annotation_lineage_marker_summary.csv'),row.names=FALSE)
cts<-rd('selection_cell_counts') |> filter(analysis=='within15')
stopifnot(sum(cts$cells)==1926,min(cts$cells)==107)
cts$condition<-factor(cts$condition_simple,levels=c('diagnosis','relapse'),labels=c('Diagnosis','Relapse'))
p<-ggplot(cts,aes(condition,cells,fill=condition))+geom_col(width=.55)+geom_text(aes(label=cells),vjust=-.5,size=3.8)+facet_grid(.~patient_pair)+scale_fill_manual(values=c(Diagnosis='#0072B2',Relapse='#D55E00'),guide='none')+scale_y_continuous(limits=c(0,520),breaks=seq(0,500,100),expand=expansion(mult=0))+labs(x=NULL,y='Cells')+theme(panel.spacing.x=grid::unit(10,'mm'))
savep(p,'Fig02_paired_sample_cell_counts',8.2,4.3)
sc<-rd('KEGG_within15');bulk<-rd('KEGG_bulk')
overview<-function(x,dataset,n=8){
 z<-x |> filter(padj<.05) |> group_by(direction) |> arrange(padj,pathway_name,.by_group=TRUE) |> slice_head(n=n) |> ungroup() |> arrange(desc(NES>0),padj,pathway_name)
 z$label<-factor(label(z$pathway_name),levels=rev(label(z$pathway_name)));z$star<-stars(z$padj)
 ggplot(z,aes(NES,label))+geom_vline(xintercept=0,color='#999999',linewidth=.35)+geom_point(size=3.1,color=datasets[[dataset]])+geom_text(aes(label=star),nudge_x=.13,hjust=0,size=3.6)+
 scale_x_continuous(limits=c(-2.6,3.6),breaks=-2:3,expand=expansion(mult=0))+labs(x='Normalized enrichment score (NES)',y=NULL,caption=keycaption)+theme(axis.text.y=element_text(size=10),axis.line.y=element_blank(),axis.ticks.y=element_blank())
}
savep(overview(sc,'Relapse scRNA-seq'),'Fig03B_scRNAseq_slide_overview',9.5,5.2)
savep(overview(bulk,'STAT5B N642H bulk'),'Fig04B_bulk_slide_overview',9.5,5.7)
savep(overview(sc,'Relapse scRNA-seq',15),'Fig03_scRNAseq_all_pathway_overview',9.5,7.3)
savep(overview(bulk,'STAT5B N642H bulk',15),'Fig04_bulk_all_pathway_overview',9.5,8.1)
# Intersection membership uses the common tested universe, never absent pathways.
cx<-read.csv(file.path(out,'03_Tables','Table06_cross_dataset_pathway_membership.csv')) |> filter(analysis=='within15',tested_in_both)
stopifnot(nrow(cx)==275)
venn<-function(up=TRUE){
 a<-cx$key[cx$sc_significant & if(up)cx$sc_NES>0 else cx$sc_NES<0];b<-cx$key[cx$bulk_significant & if(up)cx$bulk_NES>0 else cx$bulk_NES<0]
 theta<-seq(0,2*pi,length.out=301);sep<-1.2
 cir<-bind_rows(data.frame(x=cos(theta),y=sin(theta),Dataset='Relapse scRNA-seq'),data.frame(x=sep+cos(theta),y=sin(theta),Dataset='STAT5B N642H bulk'))
 t<-data.frame(x=c(-.5,.6,1.7),y=0,label=c(length(setdiff(a,b)),length(intersect(a,b)),length(setdiff(b,a))))
 ggplot(cir,aes(x,y,group=Dataset))+geom_polygon(aes(fill=Dataset),alpha=.16,color=NA)+geom_path(aes(color=Dataset),linewidth=.65)+geom_text(data=t,aes(x,y,label=label),inherit.aes=FALSE,size=5.3)+
 scale_fill_manual(values=datasets,name=NULL)+scale_color_manual(values=datasets,name=NULL)+coord_equal(xlim=c(-1.15,2.35),ylim=c(-1.1,1.2))+labs(title=if(up)'Upregulated'else'Downregulated')+theme_void(base_family='Arial')+theme(plot.title=element_text(size=12,hjust=.5),legend.position='bottom')
}
p<-(venn(TRUE)+venn(FALSE))+plot_layout(guides='collect')&theme(legend.position='bottom')
savep(p,'Fig05_direction_matched_intersection',9,3.5)
targ<-bind_rows(sc |> mutate(Dataset='Relapse scRNA-seq'),bulk |> mutate(Dataset='STAT5B N642H bulk')) |> mutate(Pathway=canon(pathway_name)) |> filter(Pathway%in%c('Proteasome','Mitophagy','Autophagy','Autophagy - other','Oxidative phosphorylation','Ribosome'))
targ$Pathway<-factor(targ$Pathway,levels=rev(c('Proteasome','Mitophagy','Autophagy','Autophagy - other','Oxidative phosphorylation','Ribosome')));targ$Dataset<-factor(targ$Dataset,levels=names(datasets));targ$star<-stars(targ$padj)
p<-ggplot(targ,aes(NES,Pathway))+geom_vline(xintercept=0,color='#999999',linewidth=.35)+geom_point(aes(color=Dataset,shape=padj<.05),size=3.1,stroke=.9)+geom_text(aes(label=star),nudge_x=.15,hjust=0,size=3.7)+facet_grid(.~Dataset)+scale_color_manual(values=datasets,guide='none')+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),guide='none')+scale_x_continuous(limits=c(-1.7,3.6),breaks=-1:3)+labs(x='Normalized enrichment score (NES)',y=NULL,caption=paste0(keycaption,'\nOpen circle: FDR >= 0.05; dash: not tested'))+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),panel.spacing.x=grid::unit(7,'mm'))
missing<-data.frame(Dataset=factor('STAT5B N642H bulk',levels=names(datasets)),Pathway=factor('Autophagy - other',levels=levels(targ$Pathway)),NES=0)
p<-p+geom_text(data=missing,aes(label='\u2014'),size=4)
savep(p,'Fig07_target_pathway_audit',9.2,4.4)
pair<-rd('target_paired_direction') |> filter(analysis=='within15',pathway%in%c('Proteasome','Mitophagy - animal')) |> mutate(pathway=canon(pathway)) |> pivot_longer(c(diagnosis,relapse),names_to='condition',values_to='score')
pair$condition<-factor(pair$condition,levels=c('diagnosis','relapse'),labels=c('Diagnosis','Relapse'))
p<-ggplot(pair,aes(condition,score,group=patient_pair,color=patient_pair,shape=patient_pair))+geom_line(linewidth=.6)+geom_point(size=3)+facet_grid(.~pathway)+scale_color_manual(values=c(P1='#0072B2',P2='#D55E00',P3='#009E73'),name='Patient')+scale_shape_manual(values=c(P1=16,P2=17,P3=15),name='Patient')+scale_y_continuous(limits=c(5.4,7.1),breaks=c(5.5,6,6.5,7))+labs(x=NULL,y='Mean log2 CPM')+theme(panel.spacing.x=grid::unit(12,'mm'))
savep(p,'Fig08_paired_leading_edge_direction',7.6,4.7)
sens<-rd('KEGG_target_summary') |> mutate(Pathway=canon(pathway_name)) |> filter(analysis%in%c('within15','global25','within25','reference95'),Pathway%in%c('Proteasome','Mitophagy'))
sens<-bind_rows(sens,bulk |> mutate(Pathway=canon(pathway_name),analysis='bulk') |> filter(Pathway%in%c('Proteasome','Mitophagy')))
sens$definition<-factor(sens$analysis,levels=rev(c('within15','within25','global25','reference95','bulk')),labels=rev(c('Within-sample top 15%','Within-sample top 25%','Global top 25%','Reference q95','STAT5B N642H bulk')))
sens$Dataset<-ifelse(sens$analysis=='bulk','STAT5B N642H bulk','Relapse scRNA-seq');sens$star<-stars(sens$padj)
p<-ggplot(sens,aes(NES,definition))+geom_point(aes(color=Dataset,shape=padj<.05),size=3.1,stroke=.9)+geom_text(aes(label=star),nudge_x=.05,hjust=0,size=3.7)+facet_grid(.~Pathway)+scale_color_manual(values=datasets,name=NULL)+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),guide='none')+scale_x_continuous(limits=c(1.3,2.1),breaks=c(1.4,1.6,1.8,2))+labs(x='Normalized enrichment score (NES)',y=NULL,caption=paste0(keycaption,'\nOpen circle: FDR >= 0.05'))+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),panel.spacing.x=grid::unit(7,'mm'))
savep(p,'FigS06_target_definition_comparison',8.7,4.2)
# Store the exact colour mapping as reproducible metadata, not an extra reader file.
write.csv(data.frame(cell_type=names(cols),colour=unname(cols)),file.path(out,'04_Code_Availability','cell_type_palette.csv'),row.names=FALSE)
cat('PASS: data identities, 13953 annotated cells, 1926 selected cells and 275 common tested pathways.\n')
