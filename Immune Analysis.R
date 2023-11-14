#####   CIBERSORTx   #####
source('CIBERSORTx.R')
LM22.file<-'LM22.txt'

exp.file<-'/data3/data/JN/CRC/TCGAbiolinks/TCGAbiolinks_combat_data_tumor.txt'
exp.file<-'/data3/data/JN/CRC/TCGA_microbes/sparseCCA/rna_RSEM.csv'
result<-CIBERSORT(LM22.file,exp.file,perm=1000,QN=F)

write.csv(result,'TCGA_CIBERSORT.csv')

#####   ESTIMATE   #####
library(estimate)

## start analysis
#exp.file = "TCGAbiolinks_combat_data_tumor.txt"
exp.file = "TCGAbiolinks_combat_data.txt"
in.gct.file = "ESTIMATE_input_ALL.gct"
outputGCT(exp.file, in.gct.file)

## merge dataset with common genes
filterCommonGenes(input.f = exp.file,output.f = in.gct.file,id = "GeneSymbol")

## start estimate
#out.score.file = "ESTIMATE_score.gct"
out.score.file = "ESTIMATE_score_ALL.gct"
estimateScore(in.gct.file,out.score.file,platform = "illumina")
#plotPurity(scores=out.score.file)

## output score
ESTIMATE_score = read.table(out.score.file,skip = 2,header = T,row.names = 1,check.names=F)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[,2:ncol(ESTIMATE_score)]))
rownames(ESTIMATE_score)<-unlist(lapply(rownames(ESTIMATE_score),function(x){gsub('\\.','-',x)}))
ESTIMATE_score$Samples = rownames(ESTIMATE_score)
ESTIMATE_score = ESTIMATE_score[,c(ncol(ESTIMATE_score),2:ncol(ESTIMATE_score)-1)]
write.csv(ESTIMATE_score, "ESTIMATE/ESTIMATE_score_ALL.csv")
