microbes_raw<-read.csv('TCGA_RNA_cri3_raw_barcode.csv',row.names=1)
microbes<-read.csv('TCGA_RNA_cri3_barcode.csv',row.names=1)

#prevlence>20%
microbes<-microbes[,which(unlist(apply(microbes_raw,2,function(x){sum(x>0)/length(x)}))>0.2)]

info<-read.csv('sample_info_no_multibam_tumor_rf.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'

info$cms_merge[which(info$cms_merge!='CMS1')]='OTHER'
wilcox_pvalue<-apply(t(microbes),1,function(x){wilcox.test(x~info[rownames(microbes),]$cms_merge,paired=FALSE)$p.value})
cmsmix_means<-rowMeans(t(microbes)[,rownames(subset(info,info$cms_merge=='CMS1'))])
other_means<-rowMeans(t(microbes)[,rownames(subset(info,info$cms_merge=='OTHER'))])
means_diff<-cmsmix_means-other_means
diff<-data.frame(cmsmix_means=cmsmix_means,other_means=other_means,means_diff=means_diff,pvalue=wilcox_pvalue)
rownames(diff)<-colnames(microbes)
diff<-subset(diff,diff$pvalue<0.05)
rcmsmix<-diff

write.csv(diff,'microbes_diff_vs_rf_CMS1vsOTHER.csv',quote=F)
