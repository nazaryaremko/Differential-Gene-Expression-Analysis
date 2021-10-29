#final code

library(GEOquery)
library(limma)
install.packages('umap')
library(umap)
install.packages('gProfileR')
library(gProfileR)

#gets the data
gse <- getGEO("GSE54839", GSEMatrix =TRUE, AnnotGPL=TRUE)
gse <- gse[[1]]
#gets the log of expression data
exprs(gse) <- log2(exprs(gse))

fvarLabels(gse) <- make.names(fvarLabels(gse))

#extracting phenotype data
sampleInfo <- pData(gse)
pData(gset)[,c("title","tissue:ch1")]

labels <- c('Cocaine', 'Control')
factor_label <- as.factor(c(rep(0,3), rep(1,3)))
boxplot(exprs(gse), boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=factor_label)
legend("topleft", labels, fill=palette(), bty="n")

#creates a design matrix where samples are separated based on
#whether they are treatment or control
#excludes intercept factor (~0)
design <- model.matrix(~0+sampleInfo$`disease state:ch1`)

#setting names of the design matrix
colnames(design) <- c("Cocaine","Control")
summary(exprs(gse))

## calculate median expression level
cutoff <- median(exprs(gse))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff

## Identify genes expressed in more than 2 samples
keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse <- gse[keep,]

#fits multiple linear models by weighted or generalized least squares
#read more
fit <- lmFit(exprs(gse), design)
head(fit$coefficients)

#expresses contrasts between a set of parameters as a numeric matrix
groups <- make.names(c("Cocaine","Control"))
cts <- paste(groups[1], groups[2], sep="-")
contrasts <- makeContrasts(Cocaine-Control, levels=design)
## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

#compute estimated coefficients and standard errors for a given set of contrasts
fit2 <- contrasts.fit(fit, contrasts)
#compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression 
#by empirical Bayes moderation of the standard errors towards a common value.
fit2 <- eBayes(fit2)
table_1 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
table_1 <- table_1[table_1$P.Value <0.05 , ]
IDs <- row.names(table_1)

#creating a list where I will store differentially expressed genes
dif_exp_genes <- c()
for (i in IDs){
  dif_exp_genes <- append(dif_exp_genes, fit2_website$genes[which(grepl(i, fit2_website$genes$ID)), 3])
}
dif_exp_genes


#table(decideTests(fit2))
#volcanoplot(fit2, coef=1, xlim=c(-2,2))

#GO analysis      
go_enrich_up <- gprofiler(as.vector(dif_exp_genes),organism='hsapiens')
go_enrich_down <- gprofiler(as.vector(dif_exp_genes),organism='hsapiens')

GOpval_oe_up <- go_enrich_up[ , c("term.id", "p.value")]
write.table(GOpval_oe_up, "/Users/nazaryaremko/Desktop/GOs_oe_up.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)

GOpval_oe_down <- go_enrich_down[ , c("term.id", "p.value")]
write.table(GOpval_oe_down, "/Users/nazaryaremko/Desktop/GOs_oe_down.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)

