# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gse_website <- getGEO("GSE54839", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gse_website) > 1) idx <- grep("GPL6947", attr(gse_website, "names")) else idx <- 1
gse_website <- gse_website[[idx]]

# make proper column names to match toptable 
fvarLabels(gse_website) <- make.names(fvarLabels(gse_website))

# group membership for all samples
gsms_website <- "000111000111000111000111000111000111000111000111000111000111"
sml_website <- strsplit(gsms_website, split="")[[1]]

# log2 transformation
ex_website <- exprs(gse_website)
qx_website <- as.numeric(quantile(ex_website, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC_website <- (qx_website[5] > 100) ||
  (qx_website[6]-qx_website[1] > 50 && qx_website[2] > 0)
if (LogC_website) { ex_website[which(ex_website <= 0)] <- NaN
exprs(gse_website) <- log2(ex_website) }

# assign samples to groups and set up design matrix
gs_website <- factor(sml_website)
groups_website <- make.names(c("control","cocaine"))
levels(gs_website) <- groups_website
gse_website$group_website <- gs_website
design_website  <- model.matrix(~group_website + 0, gse_website)
colnames(design_website) <- levels(gs_website)

fit_website <- lmFit(gse_website, design_website)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts_website <- paste(groups_website[1], groups_website[2], sep="-")
cont.matrix_website <- makeContrasts(contrasts=cts_website, levels=design_website)
fit2_website <- contrasts.fit(fit_website, cont.matrix_website)


          
# compute statistics and table of top significant genes
fit2_website <- eBayes(fit2_website, 0.01)
tT_website <- topTable(fit2_website, adjust="fdr", sort.by="B", number=Inf)

tT_website <- subset(tT_website, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT_website, file=stdout(), row.names=F, sep="\t")

as.vector(tT_website$Gene.symbol[tT_website$Cocaine...Control > 0])
          
# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2_website <- topTable(fit2_website, adjust="fdr", sort.by="B", number=Inf)
hist(tT2_website$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT_website <- decideTests(fit2_website, adjust.method="fdr", p.value=0.05)

genes_human_website <- tT2_website[tT2_website$P.Value<0.05,]$Gene.symbol
options(max.print=100000)
genes_human_website[2501:5350]

