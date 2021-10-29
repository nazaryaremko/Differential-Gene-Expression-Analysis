library(dplyr)
#gets the data
gse_mouse <- getGEO("GSE15774", GSEMatrix =TRUE, AnnotGPL=TRUE)
gse_mouse <- gse_mouse[[1]]
#gets the log of expression data

#selecting the right samples
gsms_mouse <- paste0("XXXXXXXXX1XXXXX0XXX0XXX1XXXXXXXXXXX0XXX1XXXXXXXXX0",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXX",
               "XXXXXXX0")
sml_mouse <- strsplit(gsms_mouse, split="")[[1]]

# filter out excluded samples (marked as "X")
sel_mouse <- which(sml_mouse != "X")
sml_mouse <- sml_mouse[sel_mouse]
gse_mouse <- gse_mouse[ ,sel_mouse]

exprs(gse_mouse) <- log2(exprs(gse_mouse))
fvarLabels(gse_mouse) <- make.names(fvarLabels(gse_mouse))

#extracting phenotype data
sampleInfo_mouse <- pData(gse_mouse)

#creates a design matrix where samples are separated based on
#whether they are treatment or control

#PROBABLY NEED TO GET RID OF EVERYTHING THAT IS NOT COCAINE OR SALINE
#THEN CONTINUE

gs_mouse <- factor(sml_mouse)
groups_mouse <- make.names(c("Control","Cocaine"))
levels(gs_mouse) <- groups_mouse
gse_mouse$group <- gs_mouse
design_mouse <- model.matrix(~group + 0, gse_mouse)
colnames(design_mouse) <- levels(gs_mouse)

#setting names of the design matrix

## calculate median expression level
cutoff_mouse <- median(exprs(gse_mouse))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed_mouse <- exprs(gse_mouse) > cutoff_mouse

## Identify genes expressed in more than 2 samples
keep_mouse <- rowSums(is_expressed_mouse) > 2

## check how many genes are removed / retained.

## subset to just those expressed genes
gse_mouse <- gse_mouse[keep_mouse,]

#fits multiple linear models by weighted or generalized least squares
#read more
fit_mouse <- lmFit(exprs(gse_mouse), design_mouse)
head(fit_mouse$coefficients)

#expresses contrasts between a set of parameters as a numeric matrix
cts_mouse <- paste(groups_mouse[1], groups_mouse[2], sep="-")
contrasts <- makeContrasts(contrasts_mouse=cts_mouse, levels=design_mouse)
## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

#compute estimated coefficients and standard errors for a given set of contrasts
fit2_mouse <- contrasts.fit(fit_mouse, contrasts)
#compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression 
#by empirical Bayes moderation of the standard errors towards a common value.
fit2_mouse <- eBayes(fit2_mouse)
table_1_mouse <- topTable(fit2_mouse, adjust="fdr", genelist=fit2_mouse$genes, sort.by="B", number=Inf)
write.table(table_1_mouse, file=stdout(), row.names=F, sep="\t")
table_1_mouse <- table_1_mouse[table_1_mouse$P.Value <0.05 , ]
IDs_mouse <- row.names(table_1_mouse)
dif_exp_genes_mouse <- c()

for (i in IDs_mouse){
  dif_exp_genes_mouse <- append(dif_exp_genes_mouse, tT2_mouse_website[which(grepl(i, tT2_mouse_website$ID)), 3])
}

dif_exp_genes_mouse
table(decideTests(fit2_mouse))

go_enrich_up <- gprofiler(as.vector(dif_exp_genes_mouse),organism='mmusculus')
go_enrich_down <- gprofiler(as.vector(dif_exp_genes_mouse),organism='mmusculus')
write.table(go_enrich_up, "GO_BM-MOs_YS-Macs_up_genes.csv")
write.table(go_enrich_down, "GO_BM-MOs_YS-Macs_down_genes.csv")

GOpval_oe_up_mouse <- go_enrich_up[ , c("term.id", "p.value")]
write.table(GOpval_oe_up_mouse, "/Users/nazaryaremko/Desktop/GOs_oe_up_mouse.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)

GOpval_oe_down_mouse <- go_enrich_down[ , c("term.id", "p.value")]
write.table(GOpval_oe_down_mouse, "/Users/nazaryaremko/Desktop/GOs_oe_down_mouse.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)


