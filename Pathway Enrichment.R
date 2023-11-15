library(clusterProfiler)
cms1<-read.csv('rna_diff_vs_rf_dunt_f_CMS1vsOTHER.csv',row.names=1)

library(R.utils)
library(DO.db)
library(org.Hs.eg.db)
R.utils::setOption("clusterProfiler.download.method","auto")

#to ENTREZ ID
map<- bitr(rownames(cms1), fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
cms1_kegg<-enrichKEGG(gene=map[,2],keyType = 'kegg',organism = 'hsa',
                    pvalueCutoff=1,pAdjustMethod = 'BH',qvalueCutoff = 1,
)
cms1_kegg <- setReadable(cmsmix_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.csv(subset(cmsmix_kegg,cmsmix_kegg$pvalue<0.05),'rna_KEGG_enrichment_CMS1.csv')

