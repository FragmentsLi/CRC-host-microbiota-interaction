#Step1：Import
source('lasso_tutorial.R')
genes <- load_gene_expr("rna_RSEM.txt")
microbes <- load_microbiome_abnd("microbes_vs_pre_0.2.txt")

info<-read.csv('sample_info_no_multibam_tumor_rf.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'

diff_genes<-read.csv('rna_diff_vs_rf_dunt_f_CMS1vsOTHER.csv',row.names=1)
diff_microbes<-read.csv('CMS/microbes_diff_vs_rf_CMS1vsOTHER.csv',row.names=1)

genes<-genes[rownames(subset(info,info$cms_merge=='CMS1')),rownames(diff_genes)]
microbes<-microbes[rownames(subset(info,info$cms_merge=='CMS1')),rownames(diff_microbes)]
dim(genes)
dim(microbes)

y<-genes#response
x<-microbes#predictors

#Step2：Fit LASSO model
## Extract expression of first gene in the matrix
i <- 1 ## replace with 2 or 3 to test other two genes
y_i <- y[,i]
gene_name <- colnames(y)[i]

stopifnot(class(y_i) == "numeric")

## Fit lasso model using LOOCV
fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i))
bestlambda <- fit.model$bestlambda
r.sqr <- fit.model$r.sqr

## Estimate sigma and beta init using the estimated LOOCV lambda.
## Sigma is the standard deviation of the error term or noise.
sigma.myfun <- estimate.sigma.loocv(x, y_i, bestlambda, tol=1e-4)
sigma <- sigma.myfun$sigmahat
beta <- as.vector(sigma.myfun$betahat)[-1] ## remove intercept term

## Perform inference using lasso projection method, also known as the de-sparsified Lasso,
## using an asymptotic gaussian approximation to the distribution of the estimator.
lasso.proj.fit <- lasso.proj(x, y_i, multiplecorr.method = "BH", betainit = beta, sigma = sigma, suppress.grouptesting = T)
## A few lines of log messages appear here along with a warning about substituting sigma value (i.e. standard deviation 
## of error term or noise) because we substituted value of sigma using our computation above.
## Warning message:
##   Overriding the error variance estimate with your own value.


## get 95% confidence interval for gene-taxa association
lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))

## prep lasso output dataframe
lasso.df <- data.frame(gene = rep(gene_name, length(lasso.proj.fit$pval)),
                       taxa = names(lasso.proj.fit$pval.corr),
                       r.sqr = r.sqr,
                       pval = lasso.proj.fit$pval,
                       ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
                       row.names=NULL)

## sort by p-value
lasso.df <- lasso.df[order(lasso.df$pval),]
head(lasso.df)

#Step 3: Stability selection

set.seed(0511)

## perform stability selection using glmnet lasso
stab.glmnet <- stabsel(x = x, y = y_i,
                       fitfun = glmnet.lasso, cutoff = 0.6,
                       PFER = 1)

taxa.selected <- names(stab.glmnet$selected)
if(length(taxa.selected) == 0) taxa.selected <-"None"

stabsel.df <- data.frame("gene" = gene_name, "taxa" = taxa.selected)
if(taxa.selected == "none"){
  stabsel.df$stability_selected = "no"
}else stabsel.df$stability_selected = "yes"

head(stabsel.df)

#Step 4: Merge output from steps 2. and 3. to get gene-taxa associations
overlap_lasso_stabsel <- merge(lasso.df,stabsel.df, by = c("gene","taxa"))
head(overlap_lasso_stabsel)

write.csv(overlap_lasso_stabsel,'LASSO_dunt_f_rb_FDR0.1_CMS1.csv')
