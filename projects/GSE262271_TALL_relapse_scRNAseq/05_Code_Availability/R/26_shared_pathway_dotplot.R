# Display only: retains the original full-KEGG NES and FDR.
options(stringsAsFactors=FALSE)
if(nzchar(Sys.getenv('GSE_REPRO_R_LIB'))) .libPaths(c(Sys.getenv('GSE_REPRO_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(dplyr)})
args<-commandArgs(TRUE);stopifnot(length(args)==2)
d<-read.csv(args[1],check.names=FALSE)
out<-args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(nrow(d)==40,!anyDuplicated(d$key),all(d$sc_FDR<.05),all(d$bulk_FDR<.05),all(d$sc_NES>0),all(d$bulk_NES>0))
d<-d |> arrange(pmax(sc_FDR,bulk_FDR),Pathway)
d$label<-sub(' signaling pathway$',' signaling',d$Pathway)
map<-c('Human immunodeficiency virus 1 infection'='HIV-1 infection','Human T-cell leukemia virus 1 infection'='HTLV-1 infection','Kaposi sarcoma-associated herpesvirus infection'='KSHV infection','Human cytomegalovirus infection'='Cytomegalovirus infection')
ii<-match(d$Pathway,names(map));d$label[!is.na(ii)]<-unname(map[ii[!is.na(ii)]])
stopifnot(which(d$key=='proteasome')==14L,which(d$key=='mitophagy')==28L)
make_plot<-function(rows,full=FALSE){
 q<-d[rows,]
 z<-bind_rows(q |> transmute(Pathway,label,Dataset='Relapse scRNA-seq',NES=sc_NES,FDR=sc_FDR),q |> transmute(Pathway,label,Dataset='STAT5B N642H bulk',NES=bulk_NES,FDR=bulk_FDR))
 z$Dataset<-factor(z$Dataset,levels=c('Relapse scRNA-seq','STAT5B N642H bulk'))
 z$label<-factor(z$label,levels=rev(q$label))
 z$stars<-ifelse(z$FDR<.001,'***',ifelse(z$FDR<.01,'**','*'))
 ggplot(z,aes(NES,label))+geom_point(aes(color=Dataset),size=3.3)+
 geom_text(aes(label=stars),nudge_x=.085,hjust=0,vjust=.55,size=3.7,color='#252525')+
 facet_grid(.~Dataset)+scale_y_discrete(expand=expansion(add=.7))+
 scale_x_continuous(limits=c(1.3,3.13),breaks=c(1.5,2,2.5,3),expand=expansion(mult=0))+
 scale_color_manual(values=c('Relapse scRNA-seq'='#0072B2','STAT5B N642H bulk'='#D55E00'),guide='none')+
 labs(x='Normalized enrichment score (NES)',y=NULL,caption='* FDR < 0.05     ** FDR < 0.01     *** FDR < 0.001')+
 theme_classic(base_family='Arial',base_size=11)+
 theme(axis.text.y=element_text(size=if(full)10 else 10.7,face='plain',color='#202020',margin=margin(r=8)),axis.text.x=element_text(size=10,color='#202020'),axis.title.x=element_text(size=11,margin=margin(t=8)),axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.line.x=element_line(linewidth=.4),axis.ticks.x=element_line(linewidth=.3),panel.grid.major.y=element_line(color='#EDEDED',linewidth=.25),strip.background=element_blank(),strip.text=element_text(size=11,face='plain',margin=margin(b=10)),panel.spacing.x=grid::unit(6,'mm'),plot.caption=element_text(size=9,hjust=.5,margin=margin(t=8)),plot.margin=margin(10,10,8,5))
}
save_plot<-function(p,name,w,h){
 b<-ggplot_build(p)
 stopifnot(nrow(b$data[[1]])==nrow(p$data),isTRUE(all.equal(b$data[[1]]$x,p$data$NES,tolerance=0)),identical(b$data[[2]]$label,p$data$stars),!anyNA(b$data[[1]]$x))
 ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=320,bg='white',device='png',type='cairo')
 ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h,bg='white',device=cairo_pdf)
}
save_plot(make_plot(1:40,TRUE),'Fig06_shared_pathways_dotplot',8.5,10.5)
save_plot(make_plot(1:20),'Fig06B_shared_pathways_ranks01_20',11.5,5.8)
save_plot(make_plot(21:40),'Fig06C_shared_pathways_ranks21_40',11.5,5.8)
message('PASS: all NES positions and full-KEGG FDR stars match source; uniform plain labels.')
