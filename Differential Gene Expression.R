#####   Wilcox Test   #####
rna<-read.csv('TCGAbiolinks_RSEM_normalized_data.csv',row.names=1,check.names=F)
microbes<-read.csv('TCGA_RNA_cri3_raw_barcode.csv',row.names=1)
rna<-t(rna)[rownames(microbes),]

library(dplyr)
library(biomaRt)
library(DESeq2)

#remove non-protein-coding genes
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
genes <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name", "chromosome_name","transcript_biotype"), 
                        filters = c("transcript_biotype"),values = list("protein_coding"), mart = mart)
write.csv(genes,'biomart_protein_coding_genes.csv',row.names=F)

rna<-rna[,intersect(colnames(rna),genes$external_gene_name)]
rna<-rna[,which(unlist(apply(rna,2,function(x){sum(x>0)/length(x)}))>0.5)]

rna_var<-apply(rna,2,var)
rna<-rna[,names(rna_var[rna_var>quantile(rna_var)[2]])]

write.csv(rna,'rna_RSEM.csv',row.names=T)

#wilcox test
rna<-read.csv('TCGAbiolinks_combat_data_tumor.csv',row.names=1,check.names=F)
filter_rna<-read.delim('rna_RSEM.txt',row.names=1,check.names=F)
rna<-rna[colnames(filter_rna),]
rna<-t(rna)

library(DescTools)
info<-read.csv('sample_info_no_multibam_tumor_rf.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'
info<-info[rownames(rna),]
info$cms_merge<-as.factor(info$cms_merge)

cms21_pvalue<-c()
cms31_pvalue<-c()
cms41_pvalue<-c()
cmsmix1_pvalue<-c()
kw_pvalue<-c()
set.seed(1218)
for(i in 1:ncol(rna)){
  kw_pvalue[i]<-kruskal.test(x=rna[,i],g=info$cms_merge,control='CMS1')$p.value
  dun_pvalue<-DunnettTest(x=rna[,i],g=info$cms_merge,control='CMS1')$'CMS1'
  cms21_pvalue[i]<-dun_pvalue[1,4]
  cms31_pvalue[i]<-dun_pvalue[2,4]
  cms41_pvalue[i]<-dun_pvalue[3,4]
  cmsmix1_pvalue[i]<-dun_pvalue[4,4]
}

re<-data.frame(gene=colnames(rna),
               kw_pvalue=kw_pvalue,
               cms1_means=colMeans(rna[rownames(subset(info,info$cms_merge=='CMS1')),]),
               cms2_means=colMeans(rna[rownames(subset(info,info$cms_merge=='CMS2')),]),
               cms3_means=colMeans(rna[rownames(subset(info,info$cms_merge=='CMS3')),]),
               cms4_means=colMeans(rna[rownames(subset(info,info$cms_merge=='CMS4')),]),
               cmsmix_means=colMeans(rna[rownames(subset(info,info$cms_merge=='Mixed-CMS')),]),
               cms21_bias=colMeans(rna[rownames(subset(info,info$cms_merge=='CMS2')),])-colMeans(rna[rownames(subset(info,info$cms_merge=='CMS1')),]),
               cms31_bias=colMeans(rna[rownames(subset(info,info$cms_merge=='CMS3')),])-colMeans(rna[rownames(subset(info,info$cms_merge=='CMS1')),]),
               cms41_bias=colMeans(rna[rownames(subset(info,info$cms_merge=='CMS4')),])-colMeans(rna[rownames(subset(info,info$cms_merge=='CMS1')),]),
               cmsmix1=colMeans(rna[rownames(subset(info,info$cms_merge=='Mixed-CMS')),])-colMeans(rna[rownames(subset(info,info$cms_merge=='CMS1')),]),
               cms21_pvalue=cms21_pvalue,
               cms31_pvalue=cms31_pvalue,
               cms41_pvalue=cms41_pvalue,
               cmsmix1_pvalue=cmsmix1_pvalue)

re_diff<-subset(re,(re$kw_pvalue<0.05)&(re$cms21_pvalue<0.05)&(re$cms31_pvalue<0.05)
                &(re$cms41_pvalue<0.05)&(re$cmsmix1_pvalue<0.05))
head(re_diff)
dim(re_diff)

write.csv(re_diff,'rna_diff_vs_rf_dunt_f_CMS1vsOTHER.csv',quote=F)
