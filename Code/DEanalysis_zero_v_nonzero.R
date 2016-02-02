#### DE analysis, similar to prostate_analysis_edgeR.R, except that groups are defined 
#### by SULT2B1b zero vs. nonzero

#source("https://bioconductor.org/biocLite.R")
library("edgeR")

setwd("/Users/fayezheng/Dropbox/Faye/Github/Ratliff_Prostate")

###---------- Getting SULT2B1b expression ----------###

d.full <- read.csv("GeneCountMatrix.AllBatches_withannot_.csv")
rownames(d.full) <- d.full[,1]
d.full <- d.full[,-c(1,2)]
samples.full <- names(d.full)

# Remove the control and repeat samples (for now)
ctrl.idx <- which(substr(samples.full, nchar(samples.full)-3+1, nchar(samples.full))=="trl", arr.ind=TRUE)
repeat.idx <- which(substr(samples.full, nchar(samples.full)-3+1, nchar(samples.full))=="eat", arr.ind=TRUE)
remove.idx <- c(ctrl.idx, repeat.idx)
ctrls <- d.full[, ctrl.idx]
repeats <- d.full[, repeat.idx]

d <- d.full[, -remove.idx] # 36135 genes, 399 samples

# Define groups by SULT2B1b
sult2b1b <- as.numeric(d[rownames(d)=="ENSG00000088002",])
sult2b1b.zero <- sult2b1b==0
sult2b1b.grp <- ifelse(sult2b1b.zero==TRUE, "Zero", "Nonzero")
table(sult2b1b.grp) # 352 zero vs. 47 nonzero

###---------- Getting filtered dataset ----------###

# Read in the dataset, already filtered to remove low-count genes
d <- read.table("prostatedata.txt")
batch <- substr(names(d),2,2)

#sult2b1b.grp <- ifelse(substr(names(data),1,1)=="C","Control","Knockdown")

# Note that indicator=1 for "Zero" group
design <- model.matrix(~batch+sult2b1b.grp) 
rownames(design) <- colnames(d)

###---------- Normalize, estimate disps, and BCV plot ----------###

# Calculate normalization factors to scale the library 
d <- DGEList(counts=d,group=sult2b1b.grp)
d <- calcNormFactors(d)

# Estimate dispersions
d <- estimateGLMCommonDisp(d, design, verbose=TRUE) 
# Disp = 3.74311, BCV = 1.9347
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)

# MDS plot: labeled/colored by batch
pdf("MDSplot_batch.pdf")
plotMDS(d, main="Multidimensional Scaling Plot of Batches", xlab="LogFC Distance, Dim 1", ylab="LogFC Distance, Dim 2", col=ifelse(batch=="1", "blue", ifelse(batch=="2", "black", "red")), labels=batch, cex=0.75)
legend(x=1.5, y=2.8, legend=c("Batch 1", "Batch 2", "Batch 3"), col=c("blue", "black", "red"), pch = 20, cex=0.75)
dev.off()

# MDS plot: labeled/colored by group
pdf("MDSplot_group.pdf")
plotMDS(d, main="Multidimensional Scaling Plot of Treatments", xlab="LogFC Distance, Dim 1", ylab="LogFC Distance, Dim 2", col=ifelse(sult2b1b.grp=="Zero", "blue", "red"), labels=sult2b1b.grp, cex=0.75)
legend(x=1.5, y=2.8, legend=unique(sult2b1b.grp), col=c("blue", "red"), pch = 20)
dev.off()

# Plot tagwise dispersions against log2-CPM
pdf("BCVplot.pdf")
plotBCV(d) 
dev.off()

###---------- Test batch effect ----------###

# Fit genewise GLMs
#fit <- glmFit(d, design)

# Conduct LRTs
#lrt <- glmLRT(fit, coef=2:3)
#FDR <- p.adjust(lrt$table$PValue, method="BH")
#sum(FDR < 0.5) #3880 genes DE between batches; it's good that we controlled for this

###---------- Test treatment effect ----------###

# Fit genewise GLMs
fit <- glmFit(d, design)
lrt <- glmLRT(fit)

# Total number of DE genes at FDR<0.05
summary(de <- decideTestsDGE(lrt, p=0.05))

# Table of DE results, ordered by FDR
table <- topTags(lrt, n=dim(d)[1])$table
# Table of DE results, for only DE genes
table.de <- table[table$FDR<0.05,] 

# Write to file
#write.table(table.de, "DEresults.txt", row.names=TRUE)
#write.csv(table.de, "DEresults.csv", row.names=TRUE)
write.table(table, "DEresults_full_zerovnonzero.txt", row.names=TRUE)
write.csv(table, "DEresults_full_zerovnonzero.csv", row.names=TRUE)

# Smearplot displays the log-fold changes with the DE genes highlighted:
pdf("smearplot.pdf")
DEnames <- rownames(d.filt)[as.logical(de)]
plotSmear(lrt, de.tags=DEnames)
abline(h = c(-2, 2), col = "blue")
dev.off()

# Gene ontology analysis
#go <- goana(lrt, species="Hs")
#topGO(go,ont="BP",sort="Up",n=30)

