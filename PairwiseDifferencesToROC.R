#!/usr/bin/env Rscript

# clean environment
rm(list=ls())

# skip lines related to installation of libraries because there are supposed to be already installed
skip_instalation <- scan(what="character")
# install libraries
install.packages("ape")
install.packages("data.table")
install.packages("spaa")
install.packages("pROC")

rm(skip_instalation)

# load packages avoiding warning messages
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(spaa))
suppressPackageStartupMessages(library(pROC))

# get arguments
args = commandArgs(trailingOnly=TRUE)

# test if there are two arguments: if not, return an error
if (length(args)!=2) {
  stop("Please, provide a first cvs file of profiles and a second csv file of controls\n
       USAGE: Rscript PairwiseDifferencesToROC.R Profiles.csv Controls.csv", call.=FALSE)
}
# read first input argument (ape)
## read dataframe of profiles (i.e. Profiles.csv)
## S stands for sample (i.e. rows): n = 12
## L stands for locus (i.e. columns): n = 15
## A stands for allele (i.e. data): n= 180
dfp = read.table(args[1], dec = ".", header=TRUE, sep = ",", quote = "")
## make sure that each variable of the dataframe is characters
dfp = data.frame(lapply(dfp, as.character))
## transpose dataframe
tdfp <- transpose(dfp, keep.names = "locus", make.names = "sample")

## calculate pairwise differences between all vectors (i.e. independently of the sample amount)
### multiple for-loops adding results in an empty vector
v <- integer() # create empty vector v
for (i in tdfp[, 2:ncol(tdfp)]) { # first for loop from the column 2 of the dataframe
  for (j in tdfp[, 2:ncol(tdfp)]){ # second for loop from the column 2 of the dataframe
    #print(i) # for i checking
    #print(j) # for j checking
    output <- sum(!i == j) # count the number of FALSE (i.e. reverse (!) of paired vectors i versus j)
    v <- c(v, output) # over right the output vector into the empty vector v
  }
}
### create a matrix adding the vector v by row independently of the sample amount
mat <- matrix(v, nrow = (ncol(tdfp)-1) , ncol = (ncol(tdfp)-1), byrow = TRUE)
### add row and column names of the matrix
#### keep the names of sample identifiers into a vector
sample = dfp$sample
#### add sample names to matrix rows
rownames(mat) <- sample
#### add sample names to matrix rows
colnames(mat) <- sample
### export the matrix into a csv file
write.csv(mat,file="PairwiseMatrix.csv")

# transform the matrix of pairwise differences into a dataframe
## transform the matrix into a dist object
distobj = as.dist(mat)
## transform the dist object into a long dataframe of pairwise differences
dfl = dist2list(distobj)

# perform a Receiver Operating Characteristic (ROC) analysis from the dataframe of pairwise differences
## suppose that the samples S1, S2, S3, S4, S5, S6 and S7 are positive controls (PC) of an outbreak
## suppose that samples S8, S9, S10, S11 and S12 are negative controls (NC) of an outbreak
## complete the file Controls.csv
## prepare a dataframe of positive controls (PC)
### read second input argument (ape)
dfc <- read.table(args[2], dec = ".", header=TRUE, sep = ",", quote = "")
### subset PC
dfPC <- subset(dfc,dfc$control %in% c("PC"))
### subset NC
dfNC <- subset(dfc,dfc$control %in% c("NC"))

## derive variables "col" and "row" into variables "new col" and "new row" flagging positive (PC) and negative (NC) controls
### derivation of "col"
dfl$newcol <- ifelse(dfl$col %in% dfPC$sample, "PC",
                     ifelse(dfl$col %in% dfNC$sample, "NC",
                            "error"))
### derivation of "col"
dfl$newrow <- ifelse(dfl$row %in% dfPC$sample, "PC",
                     ifelse(dfl$row %in% dfNC$sample, "NC",
                            "error"))

## keep in mind that pairwise differences of outbreak related and unrelated samples must be compared and flag the elements below:
### the outbreak related pairs (PC versus PC)
### the extra outbreak pairs (NC versus NC)
### the outbreak unrelated pairs (PC versus NC OR NC versus PC)
dfl$status <- 
  ifelse((dfl$newcol == "PC") & (dfl$newrow == "PC"), "related",
         ifelse((dfl$newcol == "NC") & (dfl$newrow == "NC"), "extra",
                ifelse((dfl$newcol == "NC") & (dfl$newrow == "PC"), "unrelated",
                       ifelse((dfl$newcol == "PC") & (dfl$newrow == "NC"), "unrelated", "error"))))

### export the dataframe into a csv file
write.csv(dfl,file="Dataframe.csv")
### subset
dflsub <- subset(dfl,dfl$status %in% c("related","unrelated"))

## run a ROC analysis and plot
### reset graph devices
graphics.off()
### open a pdf file
pdf("ROC.pdf")
### run a analysis
ROC <- roc(dflsub$status, dflsub$value,
           percent=TRUE,
           ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
           plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
           print.auc=TRUE, show.thres=TRUE)
### add confidence intervals for sensitivity and specificity into the plot
sens.ci <- ci.se(ROC, specificities=seq(0, 100, 5))
spec.ci <- ci.sp(ROC, sensitivities=seq(0, 100, 5))
plot(sens.ci, type="shape", col="blue")
plot(spec.ci, type="shape", col="green")
### close the pdf file
dev.off()

## create MetricsROC.txt where ROC metrics will be saved
sink("MetricsROC.txt")
### print working directory path
cat("Working directory path:",getwd(),"\n")
### print arguments
cat("First argument:",args[1],"\n")
cat("Second argument:",args[2],"\n")
### area under the curve of raw data (i.e. 85.34%)
auc(ROC)
### variance (i.e. 11.34112)
var <- var(ROC)
cat("Variance:" ,round(var,2),"\n")
### close MetricsROC.txt
sink()

## extract sensitivity and specificity of all thresholds of pairwise differences
thresholds <- coords(ROC, ret=c("threshold", "sensitivity", "specificity"))
## export the thresholds into a csv file
write.csv(thresholds,file="Thresholds.csv")
## extract the threshold presenting the best combination of sensitivity and specificity
best <- coords(ROC, "best", ret=c("threshold", "sensitivity", "specificity"))
## export the best threshold into a csv file
write.csv(best,file="ThresholdBest.csv")