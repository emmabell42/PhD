################################################################################
##
## This script loads Affymetrix microarray data into R, log2 transforms, quantile
## normalises, performs linear regression, and calls differentially expressed
## genes between two conditions. Input variables are capitalised.
##
################################################################################
#
########################################
#
## Load the relevant packages into R
#
########################################
library(limma)
########################################
#
## Load the data into R for processing
## as a lumi batch object
#
########################################
targets <- readTargets("TARGETS.txt")
rg <- read.maimages(targets,source="agilent")
########################################
#
## Background correct, log2 transform,
## and quantile normalise the probe
## intensity values
#
########################################
bgcorrected <- backgroundCorrect(rg,method="normexp")
ma <- normalizeWithinArrays(bgcorrected,method="loess")
ave <- avereps(ma, ID=ma$genes$ProbeName)
########################################
plot(hclust(dist(t(norm))))
boxplot(log(normexprs,base=2)
########################################
#
## Fit the probe intensities to a linear
## model
#
########################################
design <- array(0,dim=c(X,Y))
colnames(design) <- c("COLNAMES")
design[X,Y] <- 1
cont <- makeContrasts(CONTRAST="CONDITION1-CONDITION2",levels=design)
fit <- lmFit(norm,design)
fit <- contrasts.fit(fit,cont)
fit <- eBayes(fit)
########################################
#
## Annotate the probes with the names of
## their target genes
#
########################################
library(annotate)
library(ARRAY)
fit$genes$Symbol <- getSYMBOL(fit$genes$ID,"ARRAY")
########################################
#
## Output a table statistically
## evaluating differential expression
## between CONDITION1 and CONDITION2 
#
########################################
hits <- topTable(fit,number=nrow(norm))
