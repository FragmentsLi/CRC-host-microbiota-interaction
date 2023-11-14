source('enrichment_lasso_disease_specific_genes.R')

genes <- read.csv("rna_diff_vs_rf_dunt_f_CMS1vsOTHER.csv",row.names=1)
background_genes  <- rownames(genes)
length(background_genes)

interest<-read.csv('LASSO_dunt_f_rb_CMS1.csv')
genes_of_interest <- unique(interest$gene)
length(genes_of_interest)

interest$FDR<-p.adjust(interest$pval,method='BH',n=length(interest$pval))
genes_of_interest <- unique(subset(interest,interest$FDR<0.1)$gene)
length(genes_of_interest)

#MSigDB
library(msigdbr)
msigdb_C2 <- msigdbr(species = "Homo sapiens", category = "C2")
msigdb_C2 <- as.data.frame(msigdb_C2)
table(msigdb_C2$gs_subcat)
## Only keep canonical pathways, or in other words, remove CGP (Chemical and genetic perturbations)
msigdb_C2_CP <- msigdb_C2[msigdb_C2$gs_subcat!= "CGP",]
table(msigdb_C2_CP$gs_subcat)
dim(msigdb_C2_CP)
length(unique(msigdb_C2_CP$gs_name))

## Only keep pathway DBs that we want to test gene sets against 
table(msigdb_C2_CP$gs_subcat)

## Pathway DB of interest
path_DB <- c("CP:KEGG","CP:PID", "CP:REACTOME")
msigdb_C2_CP <- msigdb_C2_CP[(msigdb_C2_CP$gs_subcat %in% path_DB), ]
table(msigdb_C2_CP$gs_subcat)

## unique pathways in DBs of interest
msigdb_C2_CP_unique_pathways <- unique(msigdb_C2_CP$gs_name)
length(msigdb_C2_CP_unique_pathways) 

pathway_DB <- msigdb_C2_CP
pathways <- msigdb_C2_CP_unique_pathways

##Call enrichment function 
CRC_pathways <- msigdb_enrichment(pathway_DB, pathways, background_genes, genes_of_interest) 
## convert the list to dataframe
CRC_pathways <- do.call(rbind.data.frame, CRC_pathways)
dim(CRC_pathways) 

##sort by pval
CRC_pathways <- CRC_pathways[order(CRC_pathways$p_val),]
## MHT correction (FDR)
CRC_pathways$p_adj <- p.adjust(CRC_pathways$p_val, method = "BH")

length(which(CRC_pathways$p_val < 0.05))

CRC_pathways <- CRC_pathways[CRC_pathways$p_val < 0.05,]

## write to file
write.csv(CRC_pathways, LASSO_CMS1_msigdb_FDR0.1.csv', row.names = F)
