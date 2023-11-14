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
