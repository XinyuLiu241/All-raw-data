
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)


# library(GEOquery)
# library(limma)
# library(umap)
my_volcano=function(dat,p_cutoff=0.05,fc_cutoff=1,col=c("red","blue","black"),
                    ylab='-log10 (adj.PVal)',xlab='log2 (FoldChange)',leg.pos='right'){
  degs_dat=dat$DEG
  degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                              ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
  p=ggplot(degs_dat,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
    geom_point()+
    scale_color_manual(values=col)+#
    # geom_text_repel(
    #   data = tcga.diff$DEG[tcga.diff$DEG$adj.P.Val<p_fit & abs(tcga.diff$DEG$logFC)>fc_fit,],
    #   #aes(label = Gene),
    #   size = 3,
    #   segment.color = "black", show.legend = FALSE )+#
    theme_bw()+#
    theme(
      legend.title = element_blank(),#
      legend.position = leg.pos,
    )+
    ylab(ylab)+#
    xlab(xlab)+#
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),lty=3,col="black",lwd=0.5) +#|FoldChange|>2
    geom_hline(yintercept = -log10(p_cutoff),lty=3,col="black",lwd=0.5)#padj<0.05
  return(p)
}



#TCGA############
COAD.cli=read.delim('origin_datas/TCGA/COAD/Merge_COAD_clinical.txt')
head(COAD.cli)
COAD.cli=data.frame(Samples=paste0(COAD.cli$A0_Samples,'-01'),Age=COAD.cli$A17_Age,Gender=COAD.cli$A18_Sex,
                    T.stage=COAD.cli$A3_T,N.stage=COAD.cli$A4_N,M.stage=COAD.cli$A5_M,
                    Stage=COAD.cli$A6_Stage,
                    dataaset='COAD')
rownames(COAD.cli)=COAD.cli$Samples
COAD.exp=read.delim('origin_datas/TCGA/COAD/COAD_TPM.txt',row.names = 1,check.names = F)
table(substr(colnames(COAD.exp),14,15))
COAD.exp=COAD.exp[,which(substr(colnames(COAD.exp),14,15)%in%c('01','11'))]




READ.cli=read.delim('origin_datas/TCGA/READ/Merge_READ_clinical.txt')
head(READ.cli)
READ.cli=data.frame(Samples=paste0(READ.cli$A0_Samples,'-01'),Age=READ.cli$A17_Age,Gender=READ.cli$A18_Sex,
                    T.stage=READ.cli$A3_T,N.stage=READ.cli$A4_N,M.stage=READ.cli$A5_M,
                    Stage=READ.cli$A6_Stage,dataaset='READ')
rownames(READ.cli)=READ.cli$Samples
READ.exp=read.delim('origin_datas/TCGA/READ/READ_TPM.txt',row.names = 1,check.names = F)
table(substr(colnames(READ.exp),14,15))
READ.exp=READ.exp[,which(substr(colnames(READ.exp),14,15)%in%c('01','11'))]


######
tcga.cli=rbind(COAD.cli,READ.cli)
dim(tcga.cli)
table(tcga.cli$T.stage)
tcga.cli$T.stage=gsub('[ab]','',tcga.cli$T.stage)
tcga.cli$T.stage[tcga.cli$T.stage %in% c('','Tis')]=NA

table(tcga.cli$N.stage)
tcga.cli$N.stage=gsub('[abc]','',tcga.cli$N.stage)
tcga.cli$N.stage[tcga.cli$N.stage %in% c('','NX')]=NA

table(tcga.cli$M.stage)
tcga.cli$M.stage=gsub('[ab]','',tcga.cli$M.stage)
tcga.cli$M.stage[tcga.cli$M.stage %in% c('','MX')]=NA

table(tcga.cli$Stage)
tcga.cli$Stage=gsub('Stage ','',tcga.cli$Stage)
tcga.cli$Stage=gsub('[ABC]','',tcga.cli$Stage)
tcga.cli$Stage[tcga.cli$Stage %in% c('')]=NA

tcga.survival=read.delim('origin_datas/TCGA/survival_COADREAD_survival.txt')
colnames(tcga.survival)
tcga.survival=data.frame(Samples=tcga.survival$sample,
                         tcga.survival[,c("OS","OS.time","DSS","DSS.time","DFI","DFI.time","PFI","PFI.time")])

tcga.cli=merge(tcga.cli,tcga.survival,by='Samples')
tcga.cli$OS.time/365
tcga.cli=tcga.cli[tcga.cli$OS.time>30 ,]
tcga.cli=tcga.cli%>%drop_na(OS)
rownames(tcga.cli)=tcga.cli$Samples

table(rownames(COAD.exp)==rownames(READ.exp))
tcga.data=cbind(COAD.exp,READ.exp)

sample_T=colnames(tcga.data)[which(substr(colnames(tcga.data),14,15)=='01')]#肿瘤样本
sample_T=intersect(sample_T,tcga.cli$Samples)
sample_N=colnames(tcga.data)[which(substr(colnames(tcga.data),14,15)=='11')]#正常样本
tcga_type=data.frame(Samples=c(sample_T,sample_N),Type=rep(c('Tumor','Normal'),c(length(sample_T),length(sample_N))))
rownames(tcga_type)=tcga_type$Samples
table(tcga_type$Type)

tcga.cli=tcga.cli[sample_T,]
dim(tcga.cli)

tcga.data=tcga.data[,tcga_type$Samples]
range(tcga.data)
tcga.data=log2(tcga.data+1)

tcga.exp=tcga.data[,sample_T]
range(tcga.exp);dim(tcga.exp)


#01.#####
dir.create('results/01.DElncRNA')
genecode=read.delim
table(genecode$TYPE)
lncrna_genecode=genecode[which(genecode$TYPE=='lncRNA'),]
head(lncrna_genecode)


tcga.limma.lnc=mg_limma_DEG(exp = tcga.data[rownames(tcga.data)%in%lncrna_genecode$SYMBOL,tcga_type$Samples],
                            group = tcga_type$Type,ulab = 'Tumor',dlab = 'Normal')
tcga.limma.lnc$Summary
tcga.DElnc=tcga.limma.lnc$DEG[tcga.limma.lnc$DEG$adj.P.Val<0.05 & abs(tcga.limma.lnc$DEG$logFC)>log2(1.2),]
head(tcga.DElnc);dim(tcga.DElnc)
write.csv(tcga.DElnc,'results/01.DElncRNA/tcga_DElncRNA.csv')
tcga.DElnc[grep('CASC',rownames(tcga.DElnc)),]
tcga.DElnc['CASC15',]
fig1a=my_volcano(dat = tcga.limma.lnc,p_cutoff = 0.05,fc_cutoff = log2(1.2),col = c("orange","blue","grey"),)
fig1a

tcga.DElnc.filter=tcga.DElnc
tcga.DElnc.filter$type=ifelse(tcga.DElnc.filter$logFC>0,'up','down')
tcga.DElnc.filter <- tcga.DElnc.filter %>%
  mutate(gene = rownames(tcga.DElnc.filter))  %>%
  group_by(type) %>% slice_max(n =50, order_by = abs(logFC)) %>%  arrange(logFC)
head(tcga.DElnc.filter)

cli_anno=tcga_type[,'Type',drop=F]
fig1b=pheatmap(tcga.data[tcga.DElnc.filter$gene,rownames(cli_anno)],
               scale = 'row',name = 'Expression',  main="TOP50 DElncRNA in TCGA", # 
               color =  circlize::colorRamp2(c(-3, 0, 3), c('#009B9F', "white", '#C75DAA')),
               annotation_col = cli_anno,
               annotation_colors = list(Type=c(Tumor='#FB8072',Normal='#80B1D3')),
               cluster_cols = F, # 
               cluster_rows = T,
               show_rownames = F, #
               show_colnames = F)
library(ggplotify)
fig1b = as.ggplot(fig1b)
fig1b


fig1=mg_merge_plot(fig1a,fig1b,labels = c('A','B'))
ggsave('results/01.DElncRNA/Fig1.pdf',fig1,height = 6,width = 12)


upDElnc.enrich=read.delim('results/01.DElncRNA/up_DElncRNA_enrichment.txt')
head(upDElnc.enrich)
upDElnc.enrich.fit=upDElnc.enrich %>% slice_max(n = 10, order_by = Count)
upDElnc.enrich.fit
fig1c=ggplot(data=upDElnc.enrich.fit,aes(x=Count,y=reorder(Set,Count), color = -log10(FDR))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_continuous(type = "gradient")+
  scale_color_gradient(low = "blue", high = "red")+
  labs(x='Count',y='',title = 'Cancer Hallmark')+theme_bw()+
  theme(text = element_text(family = 'Times',size = 14))
fig1c

upDElnc.enrich2=read.delim('results/01.DElncRNA/up_DElncRNA_Tumor_Metastasis.txt')
head(upDElnc.enrich2)
# upDElnc.enrich2.fit=upDElnc.enrich2 %>% slice_max(n = 10, order_by = Count)

fig1d=ggplot(data=upDElnc.enrich2,aes(x=Count,y=reorder(Set,Count), color = -log10(FDR))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_continuous(type = "gradient")+
  scale_color_gradient(low = "blue", high = "red")+
  labs(x='Count',y='',title = 'Tumor metastasis')+theme_bw()+
  theme(text = element_text(family = 'Times',size = 14))
fig1d

fig1=mg_merge_plot(fig1a,fig1b,fig1c,fig1d,labels = LETTERS[1:4])
ggsave('results/01.DElncRNA/Fig1.pdf',fig1,height = 12,width = 12)




#02.CASC family#######
dir.create('results/02.Diagnosis_model')
CASC.family.genes=readxl::read_xls('origin_datas/Cancer Susceptibility family.xls')
CASC.family.genes=CASC.family.genes$Symbol
length(CASC.family.genes)
#14


CASC.DEGs=intersect(CASC.family.genes,rownames(tcga.DElnc))
CASC.DEGs
library(eulerr)
v=list(CASC.family.genes,rownames(tcga.DElnc))
names(v)=c('CASC related genes','TCGA DElncRNA')
fig2a=plot(venn(v),labels = list(col = "gray20", font = 2), 
           edges = list(col="gray60", lex=1),
           fills = list(fill = c("blue", "pink"), alpha = 0.6),
           quantities = list(cex=.8, col='gray20'))
fig2a
ggsave('results/02.Diagnosis_model/Fig2a.pdf',fig2a,height = 4.5,width = 5)

CASC.df=data.frame(Tissue=tcga_type$Type,t(tcga.data[CASC.DEGs,tcga_type$Samples]))
head(CASC.df)
CASC.df=melt(CASC.df)
table(CASC.df$Tissue)

fig2b=ggplot(CASC.df,aes(x=variable, y=value,fill=Tissue)) +
  geom_boxplot(width=0.3,outlier.colour = NA)+
  scale_fill_manual(values = c('blue','orange'))+
  theme_classic(base_size = 20)+ylab('Expression level')+xlab('')+
  ggpubr::stat_compare_means(aes(group=Tissue), label = 'p.signif', method ='wilcox.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'top',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
fig2b 
ggsave('results/02.Diagnosis_model/Fig2b.pdf',fig2b,height = 4.5,width = 5)

#######
hub.genes=CASC.DEGs
ML_dat=t(tcga.data[hub.genes,tcga_type$Samples])
ML_dat=cbind.data.frame(type=tcga_type$Type,ML_dat[tcga_type$Samples,])
ML_dat$type=as.factor(ML_dat$type)
head(ML_dat)
mat=as.matrix(ML_dat[,hub.genes])
group=as.factor(ML_dat$type)


#
library(pROC)
array.roc <- list()
for (i in 1:length(hub.genes)){
  roc1 <- roc(tcga_type$Type, t(tcga.data[hub.genes[i],tcga_type$Samples]))
  array.roc[[i]]=roc1
  names(array.roc)[i] <- paste0(hub.genes[i],' AUC=',round(roc1$auc[1],2))
}
fig2c=ggroc(array.roc,alpha = 1,size = 1)+
  geom_segment(aes(x = 1, y = 0, xend =0, yend = 1), color="darkgrey", linetype="dashed")+
  scale_colour_manual(name="lncRNA",values =RColorBrewer::brewer.pal(11,"Dark2"))+
  ggtitle("ROC Curves") +
  theme_bw()+ theme(panel.grid = element_blank(),legend.position = c(0.75,0.25),
                    text = element_text(family = 'Times'))
fig2c
ggsave('results/02.Diagnosis_model/Fig2c.pdf',fig2c,height = 4.5,width = 5)


library(ggbiplot)
tcga.pca <- prcomp(t(tcga.data[hub.genes,tcga_type$Samples]), scale=T)
fig2d <- ggbiplot(tcga.pca, scale=1, groups = tcga_type$Type,
                  ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_manual(values =c('blue','orange')) + 
  theme_bw() +xlab('PCA1') + ylab('PCA2') +
  theme(legend.direction = 'horizontal', legend.position = 'top',text = element_text(family = 'Times',size=14)) 
fig2d


fig2=mg_merge_plot(fig2a,fig2b,fig2c,fig2d,nrow=2,ncol=2,labels = LETTERS[1:4])
fig2
ggsave('results/02.Diagnosis_model/Fig2.pdf',fig2,height = 9,width = 10)


#03.RNA_Protein_Interaction#################
dir.create('results/03.RNA_Protein_Interaction')
RNA_Protein_interaction=read.delim('results/03.RNA_Protein_Interaction/RNA_Protein_Interaction.txt')
RNA_Protein_interaction=RNA_Protein_interaction[,c('Set','LncRNA')]
colnames(RNA_Protein_interaction)[1]='Protein'
head(RNA_Protein_interaction)

data_long <- RNA_Protein_interaction %>%  separate_rows(LncRNA, sep = ";")
data_long=as.data.frame(data_long)
write_tsv(data_long,file = 'results/03.RNA_Protein_Interaction/RNA_Protein_Interaction_long.tsv')

lncRNA=hub.genes
protein=RNA_Protein_interaction$Protein
nodes = tibble(gene=c(lncRNA,protein))
nodes$Type = ifelse(nodes$gene %in% protein,"protein","lncRNA")
nodes
write.table(nodes,'results/03.RNA_Protein_Interaction/nodes_anno.txt',quote = F,row.names = F,sep = '\t')


#04.############
dir.create('results/04.Immune_pathway')
tcga.est=read.delim('results/04.Immune_pathway/TCGA_ESTIMATE_score.txt',row.names = 1,check.names = F)
head(tcga.est)

hub.immune.df1=cbind(tcga.est[tcga.cli$Samples,1,drop=F],
                     t(tcga.exp[hub.genes,tcga.cli$Samples]))
head(hub.immune.df1)
hub.immune.df1=as.data.frame(hub.immune.df1)
hub.immune.df1=reshape::melt(hub.immune.df1,id='StromalScore')
fig4a=ggplot(hub.immune.df1, aes(x = value, y = StromalScore))+
  geom_smooth(aes(color = variable, fill = variable), method = "lm") +
  xlab('Hub genes Expression')+ylab('StromalScore')+
  ggpubr::stat_cor(aes(color = variable),method = 'spearman',label.x.npc = .6,label.y.npc = 1)+
  theme_bw()+
  theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
fig4a

hub.immune.df2=cbind(tcga.est[tcga.cli$Samples,2,drop=F],
                     t(tcga.exp[hub.genes,tcga.cli$Samples]))
head(hub.immune.df2)
hub.immune.df2=as.data.frame(hub.immune.df2)
hub.immune.df2=reshape::melt(hub.immune.df2,id='ImmuneScore')
fig4b=ggplot(hub.immune.df2, aes(x = value, y = ImmuneScore))+
  geom_smooth(aes(color = variable, fill = variable), method = "lm") +
  xlab('Hub genes Expression')+
  ggpubr::stat_cor(aes(color = variable),method = 'spearman',label.x.npc = .6,label.y.npc = 1)+
  theme_bw()+
  theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
fig4b


hub.immune.df3=cbind(tcga.est[tcga.cli$Samples,4,drop=F],
                     t(tcga.exp[hub.genes,tcga.cli$Samples]))
head(hub.immune.df3)
hub.immune.df3=as.data.frame(hub.immune.df3)
hub.immune.df3=reshape::melt(hub.immune.df3,id='TumorPurity')
fig4c=ggplot(hub.immune.df3, aes(x = value, y = TumorPurity))+
  geom_smooth(aes(color = variable, fill = variable), method = "lm") +
  xlab('Hub genes Expression')+
  ggpubr::stat_cor(aes(color = variable),method = 'spearman',label.x.npc = .6,label.y.npc = .5)+
  theme_bw()+
  theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
fig4c




load('results/04.Immune_pathway/tcga.hall.ssGSEA.RData')
dim(tcga.hall.ssGSEA)
rownames(tcga.hall.ssGSEA)=gsub('HALLMARK_','',rownames(tcga.hall.ssGSEA))
tcga.hall.ssGSEA[1:5,1:5]

hub.pathway.df2=cbind.data.frame(t(tcga.hall.ssGSEA[,tcga.cli$Samples]),t(tcga.exp[hub.genes,tcga.cli$Samples]))
head(hub.pathway.df2)
hub.pathway.df2=as.data.frame(hub.pathway.df2)
cor_res <- Hmisc::rcorr(as.matrix(hub.pathway.df2),type = 'spearman')
cor_res$P[is.na(cor_res$P)] <- 0

hallmark.anno=read.xlsx('HALLMARK_geneset.xlsx',check.names = F)
head(hallmark.anno)
hallmark.anno1=hallmark.anno[,'Process.category',drop=F]
rownames(hallmark.anno1)=hallmark.anno$Hallmark.name
head(hallmark.anno1)

cor_res$lab<-ifelse(cor_res$P<0.0001,'****',
                    ifelse(cor_res$P<0.001,'***',
                           ifelse(cor_res$P<0.01,'**',
                                  ifelse(cor_res$P<0.05,'*',''))))


anno.cols=pal_nejm()(8)
names(anno.cols)=unique(hallmark.anno1$Process.category)

color_anno=list(Process.category=anno.cols)

fig4d=pheatmap(mat = cor_res$r[rownames(hallmark.anno1),hub.genes],
               scale = 'none', main="",display_numbers = T,name = 'Correlation',
               number_format  = cor_res$lab[rownames(hallmark.anno1),hub.genes],  # 
               color =  colorRampPalette(c('deepskyblue', "grey", 'yellow'))(100),
               annotation_row = hallmark.anno1,
               annotation_colors = color_anno,
               cluster_cols = F,cluster_rows = F,
               show_rownames = T,show_colnames = T)

library(ggplotify)
fig4d = as.ggplot(fig4d)
fig4d


fig4=mg_merge_plot(mg_merge_plot(fig4a,fig4b,fig4c,nrow=3,ncol=1,labels = LETTERS[1:3],common.legend = T),
                   fig4d,nrow=1,ncol=2,labels = c('','D'),widths = c(1,1.2))
ggsave('results/04.Immune_pathway/Fig4.pdf',fig4,height = 15,width = 12)



#05.##########
dir.create('results/05.Clinical')
hub.genes
# "CASC15" "CASC16" "CASC8"  "CASC9"  "CASC19" "CASC18"
tcga.cli.merge=cbind.data.frame(t(tcga.exp[CASC.family.genes,tcga.cli$Samples]),tcga.cli)
head(tcga.cli.merge)
summary(coxph(formula=Surv(OS.time, OS)~CASC15,data=tcga.cli.merge))
summary(coxph(formula=Surv(DSS.time, DSS)~CASC15,data=tcga.cli.merge))
summary(coxph(formula=Surv(DFI.time, DFI)~CASC15,data=tcga.cli.merge))
summary(coxph(formula=Surv(PFI.time, PFI)~CASC15,data=tcga.cli.merge))


library(survival)
library(survminer)

km.list=list()
for (i in 1:6) {
  cutoff<-surv_cutpoint(tcga.cli.merge,time="PFI.time",event="PFI",variables=hub.genes[i])
  summary(cutoff)
  cutoff$cutpoint$cutpoint
  tcga.cli.merge$group <- ifelse(tcga.cli.merge[,hub.genes[i]] > cutoff$cutpoint$cutpoint, 'High', 'Low')
  km.list[[i]]=ggsurvplot(fit=survfit(Surv(PFI.time/365, PFI) ~ group,data = tcga.cli.merge),
                          data=tcga.cli.merge, surv.median.line = 'hv',
                          conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                          linetype = c("solid", "dashed","strata")[1],
                          palette = c('#756AB6','#EEC759'),legend = 'top',
                          legend.title = hub.genes[i],legend.labs=c('High','Low'))
  km.list[[i]]=mg_merge_plot( km.list[[i]]$plot, km.list[[i]]$table,heights = c(2.5,1),nrow=2,ncol=1)
  
}
length(km.list)
fig5a=mg_merge_plot(km.list,nrow=2,ncol=3)
fig5a


CASC15.df2=data.frame(tcga.cli,gene=as.numeric(tcga.exp['CASC15',tcga.cli$Samples]))
head(CASC15.df2)

fig5b=CASC15.df2%>%drop_na(Stage)%>%
  ggplot(aes(x=Stage, y=gene,fill=Stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  theme_classic(base_size = 20)+ylab('CASC15 Expression')+
  ggpubr::stat_compare_means(aes(group=Stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
fig5b


fig5c=CASC15.df2%>%drop_na(T.stage)%>%
  ggplot(aes(x=T.stage, y=gene,fill=T.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC15 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=T.stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
fig5c

fig5d=CASC15.df2%>%drop_na(N.stage)%>%
  ggplot(aes(x=N.stage, y=gene,fill=N.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC15 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=N.stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
fig5d

fig5e=CASC15.df2%>%drop_na(M.stage)%>%
  ggplot(aes(x=M.stage, y=gene,fill=M.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC15 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=M.stage), label = 'p.format', method ='wilcox.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
fig5e

fig5=mg_merge_plot(fig5a,
                   mg_merge_plot(fig5b,fig5c,fig5d,fig5e,ncol = 4,widths = c(1.2,1.2,1,0.8)),
                   nrow = 2,heights = c(2,1),labels = c('A','B'))
ggsave('results/05.Clinical/Fig5.pdf',fig5,height = 14,width = 14)


###########
head(tcga.cli.merge)
hub.genes

figs1a=tcga.cli.merge%>%drop_na(Stage)%>%
  ggplot(aes(x=Stage, y=CASC16,fill=Stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  theme_classic(base_size = 20)+ylab('CASC16 Expression')+
  ggpubr::stat_compare_means(aes(group=Stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs1b=tcga.cli.merge%>%drop_na(T.stage)%>%
  ggplot(aes(x=T.stage, y=CASC16,fill=T.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC16 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=T.stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs1c=tcga.cli.merge%>%drop_na(N.stage)%>%
  ggplot(aes(x=N.stage, y=CASC16,fill=N.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC16 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=N.stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs1d=tcga.cli.merge%>%drop_na(M.stage)%>%
  ggplot(aes(x=M.stage, y=CASC16,fill=M.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC16 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=M.stage), label = 'p.format', method ='wilcox.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))

figs2a=tcga.cli.merge%>%drop_na(Stage)%>%
  ggplot(aes(x=Stage, y=CASC8,fill=Stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  theme_classic(base_size = 20)+ylab('CASC8 Expression')+
  ggpubr::stat_compare_means(aes(group=Stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs2b=tcga.cli.merge%>%drop_na(T.stage)%>%
  ggplot(aes(x=T.stage, y=CASC8,fill=T.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC8 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=T.stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs2c=tcga.cli.merge%>%drop_na(N.stage)%>%
  ggplot(aes(x=N.stage, y=CASC8,fill=N.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC8 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=N.stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs2d=tcga.cli.merge%>%drop_na(M.stage)%>%
  ggplot(aes(x=M.stage, y=CASC8,fill=M.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC8 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=M.stage), label = 'p.format', method ='wilcox.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs3a=tcga.cli.merge%>%drop_na(Stage)%>%
  ggplot(aes(x=Stage, y=CASC9,fill=Stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  theme_classic(base_size = 20)+ylab('CASC9 Expression')+
  ggpubr::stat_compare_means(aes(group=Stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs3b=tcga.cli.merge%>%drop_na(T.stage)%>%
  ggplot(aes(x=T.stage, y=CASC9,fill=T.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC9 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=T.stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs3c=tcga.cli.merge%>%drop_na(N.stage)%>%
  ggplot(aes(x=N.stage, y=CASC9,fill=N.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC9 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=N.stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs3d=tcga.cli.merge%>%drop_na(M.stage)%>%
  ggplot(aes(x=M.stage, y=CASC9,fill=M.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC9 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=M.stage), label = 'p.format', method ='wilcox.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))



figs4a=tcga.cli.merge%>%drop_na(Stage)%>%
  ggplot(aes(x=Stage, y=CASC19,fill=Stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  theme_classic(base_size = 20)+ylab('CASC19 Expression')+
  ggpubr::stat_compare_means(aes(group=Stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs4b=tcga.cli.merge%>%drop_na(T.stage)%>%
  ggplot(aes(x=T.stage, y=CASC19,fill=T.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC19 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=T.stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs4c=tcga.cli.merge%>%drop_na(N.stage)%>%
  ggplot(aes(x=N.stage, y=CASC19,fill=N.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC19 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=N.stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs4d=tcga.cli.merge%>%drop_na(M.stage)%>%
  ggplot(aes(x=M.stage, y=CASC19,fill=M.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC19 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=M.stage), label = 'p.format', method ='wilcox.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))

figs5a=tcga.cli.merge%>%drop_na(Stage)%>%
  ggplot(aes(x=Stage, y=CASC18,fill=Stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  theme_classic(base_size = 20)+ylab('CASC18 Expression')+
  ggpubr::stat_compare_means(aes(group=Stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs5b=tcga.cli.merge%>%drop_na(T.stage)%>%
  ggplot(aes(x=T.stage, y=CASC18,fill=T.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC18 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=T.stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs5c=tcga.cli.merge%>%drop_na(N.stage)%>%
  ggplot(aes(x=N.stage, y=CASC18,fill=N.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC18 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=N.stage), label = 'p.format', method ='kruskal.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


figs5d=tcga.cli.merge%>%drop_na(M.stage)%>%
  ggplot(aes(x=M.stage, y=CASC18,fill=M.stage)) +
  geom_violin()+
  geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+ylab('CASC18 Expression')+
  scale_fill_manual(values = c("#FB8072","#80B1D3","#FDB462","#B3DE69"))+
  ggpubr::stat_compare_means(aes(group=M.stage), label = 'p.format', method ='wilcox.test')+
  theme(axis.text = element_text(color = 'black'),text = element_text(family = 'Times'),
        title = element_text(size = 12),legend.position = 'none',
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))



figs1=mg_merge_plot(figs1a,figs1b,figs1c,figs1d,ncol=4,widths = c(1.2,1.2,1,0.8))
figs2=mg_merge_plot(figs2a,figs2b,figs2c,figs2d,ncol=4,widths = c(1.2,1.2,1,0.8))
figs3=mg_merge_plot(figs3a,figs3b,figs3c,figs3d,ncol=4,widths = c(1.2,1.2,1,0.8))
figs4=mg_merge_plot(figs4a,figs4b,figs4c,figs4d,ncol=4,widths = c(1.2,1.2,1,0.8))
figs5=mg_merge_plot(figs5a,figs5b,figs5c,figs5d,ncol=4,widths = c(1.2,1.2,1,0.8))

fig_supple=mg_merge_plot(figs1,figs2,figs3,figs4,figs5,nrow=5,labels = LETTERS[1:5])
ggsave('results/05.Clinical/FigS1.pdf',fig_supple,height = 15,width = 12)


save.image(file = 'project.RData')




