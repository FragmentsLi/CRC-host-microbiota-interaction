library(vegan)
genes<-read.csv('TCGAbiolinks_RSEM_normalized_data.csv',row.names=1,check.names=F)
microbes<-read.csv('TCGA_RNA_cri3_raw_barcode.csv',row.names=1)
genes<-t(genes)[rownames(microbes),]

otu_dis<-vegdist(microbes,method='bray')
rna_dis<-vegdist((genes+1),method='aitchison')

mds.o <- monoMDS(otu_dis)
mds.r <- monoMDS(rna_dis)
set.seed(123)
proc<-procrustes(X=mds.o,Y=mds.r,symmetric=TRUE)
summary(proc)
prot<-protest(X=mds.o,Y=mds.r,permutations=how(nperm=9999))
