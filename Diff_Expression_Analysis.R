#final code

library(GEOquery)
library(limma)
install.packages('umap')
library(umap)

#gets the data
#paramters:
#GSEMAtrix: A boolean telling GEOquery whether or not to use GSE Series Matrix files from GEO
#AnnotGPL: A boolean defaulting to FALSE as to whether or not to use the Annotation GPL information.
gse <- getGEO("GSE54839", GSEMatrix =TRUE, AnnotGPL=TRUE)
gse <- gse[[1]]

pData(gse) ## print the sample information
fData(gse) ## print the gene annotation
exprs(gse) ## print the expression data

#gets the log of expression data
exprs(gse) <- log2(exprs(gse))

#extracting phenotype data
sampleInfo <- pData(gse)
pData(gset)[,c("title","tissue:ch1")]

#creates names for the valriable labels from character vectors stored in gse
fvarLabels(gse) <- make.names(fvarLabels(gse))

#creating levels for the legend
labels <- c('Cocaine', 'Control')
factor_label <- as.factor(c(rep(0,3), rep(1,3)))
par(mar=c(8,5,2,3)+.1)
#parameters:
#boxwex: scale factor to be applied to all boxes.
#notch: if notch is TRUE, a notch is drawn in each side of the boxes.
#main: main title of teh boxplot
#outline: if outline is not true, the outliers are not drawn (as points whereas S+ uses lines).
#col: if col is non-null it is assumed to contain colors to be used to colour the bodies of the box plots.
boxplot(exprs(gse), boxwex=0.6, notch=T, xlab = 'Sample ID', ylab = 'Normalized Expression Value', main=title, outline=FALSE, col=factor_label)
#creates legend for the boxplot
#parameters
#fill: if specified, this argument will cause boxes filled with the specified colors (or shaded in the specified colors) to appear beside the legend text.
#bty: the type of box to be drawn around the legend
#bg: the background color for the legend box. (Note that this is only used if bty != "n".)
legend("topleft", labels, fill=palette(), bty="n", bg="black")

#creates a design matrix where samples are separated based on
#whether they are treatment or control
#excludes intercept factor (~0)
design <- model.matrix(~0+sampleInfo$`disease state:ch1`)

#setting names of the design matrix
colnames(design) <- c("Cocaine","Control")

#calculate median expression level
cutoff <- median(exprs(gse))
#TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff
#Identify genes expressed in more than 2 samples
keep <- rowSums(is_expressed) > 2

#check how many genes are removed / retained.
table(keep)

#subset to just those expressed genes
gse <- gse[keep,]

#fits multiple linear models by weighted or generalized least squares using design matrix
fit <- lmFit(exprs(gse), design)
#prints out the first couple of coefficients
head(fit$coefficients)

#expresses contrasts between a set of parameters as a numeric matrix
contrasts <- makeContrasts(Cocaine-Control, levels=design)

#compute estimated coefficients and standard errors for a given set of contrasts
fit2 <- contrasts.fit(fit, contrasts)

#compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression 
#by empirical Bayes moderation of the standard errors towards a common value.
fit2 <- eBayes(fit2)

#creates a table that contains each genes and various statistics such as fold change, p-value, adjusted p-value, etc
#parameters:
#adjust: specifies how we adjust for multiple corrections
#sort.by: character string specifying statistic to rank genes by
#number: number of entries to include 
table_1 <- topTable(fit2, adjust="fdr", sort.by="p", number=Inf)
#filtering out all genes which have an adjusted p-value higherthan 0.05
table_1 <- table_1[table_1$adj.P.Val <0.05, ]

IDs <- row.names(table_1)
#creating a list where I will store differentially expressed genes
dif_exp_genes <- c()
#here I am using a code from GEO website to create gene annotations based on their ID 
for (i in IDs){
  dif_exp_genes <- append(dif_exp_genes, fit2_website$genes[which(grepl(i, fit2_website$genes$ID)), 3])
}
dif_exp_genes

#adding a new column to my table
table_1$Gene_name <- dif_exp_genes

#adding a joined gene symbol - p-value column
for (i in 1:nrow(table_1)){
  table_1$Gene_name[i] <- str_remove_all(table_1$Gene_name[i], ' ')
  if (table_1$Gene_name_and_p_value[i] != '') {
    table_1$Gene_name_and_p_value[i] <- paste(table_1$Gene_name_and_p_value[i], table_1$adj.P.Val[i], sep=", ")
  }
}

#creating a final table
write.csv(table_1, "/Users/nazaryaremko/Desktop/diff_exp_genes.csv")



