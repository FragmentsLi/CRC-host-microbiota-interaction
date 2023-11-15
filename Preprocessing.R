######   Normalization   ######

#Host gene expression matrix
data<-read.csv('TCGAbiolinks_RSEM_normalized_data.csv',row.names=1,check.names=F)
meta_RNA<-read.csv('meta_RNA_with_submitter_id_barcode.csv',row.names=1)
data<-data[,rownames(meta_RNA)]

#log2
data<-as.matrix(log2(data+1))

#ComBat
library(sva)
combat_data <- ComBat(dat = data, batch = meta_RNA$platform)

write.csv(combat_data,'TCGAbiolinks_combat_data_tumor.csv',row.names=T)

######   CMS classification   ######

data<-combat_data

library(clusterProfiler)
library(org.Hs.eg.db)
map<- bitr(rownames(data), fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)

data_id<-lapply(unique(map$ENTREZID),function(x){return(data[subset(map,map$ENTREZID==x)[1,]$SYMBOL,])})
data_id_matrix<-do.call("rbind", data_id)
rownames(data_id_matrix)<-unique(map$ENTREZID)

# CMSclassifier
library(CMSclassifier)
rna<-data_id_matrix
cms_label<-read.delim('cms_labels_public_all.txt')
cms_label<-subset(cms_label,cms_label$dataset=='tcga')
rownames(cms_label)<-cms_label$sample

library(stringr)
patient_id<-unlist(lapply(colnames(rna),function(x){substr(x,1,12)}))
cms_class<-unlist(lapply(patient_id,function(x)
  if(x %in% cms_label$sample)return(cms_label[x,5])
  else return('')))

#Random forest
Rfcms <- CMSclassifier::classifyCMS(rna,method="RF")[[3]]
cms_predict<-unlist(lapply(colnames(rna),function(x){Rfcms[x,'RF.predictedCMS']}))
cms_nearest<-unlist(lapply(colnames(rna),function(x){Rfcms[x,'RF.nearestCMS']}))
cms_merge<-unlist(lapply(patient_id,function(x)
  if(x %in% cms_label$sample)return(cms_label[x,5])
  else return(Rfcms[x,'RF.predictedCMS'])))
cms_merge[which(is.na(cms_merge))]='NOLBL'
  
info<-data.frame(barcode=barcode,patient_id=patient_id,cms_class=cms_class,
                 cms_predict=cms_predict,cms_nearest=cms_nearest,cms_merge=cms_merge)
write.csv(info,'TCGA_microbes/sample_info_no_multibam_tumor_rf.csv',row.names=F)                         

######   PCA   ######
data<-read.csv('TCGAbiolinks_RSEM_normalized_data.csv',row.names=1,check.names=F)
combat_data<-read.csv('TCGAbiolinks_combat_data_tumor.csv',row.names=1,check.names=F)
data<-data[,colnames(combat_data)]
meta_RNA<-read.csv('meta_RNA_with_submitter_id_barcode.csv',row.names=1)

#remove genes with expression equal to 0
combat_data<-combat_data[which(rowSums(data)!=0),]

#PCA
pca <- prcomp(t(combat_data), scale=TRUE)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)

                         
