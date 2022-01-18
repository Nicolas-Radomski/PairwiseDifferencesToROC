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
install.packages("plyr")

rm(skip_instalation)

# load packages avoiding warning messages
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(spaa))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(plyr))

# get arguments
args = commandArgs(trailingOnly=TRUE)

# test if there are two arguments: if not, return an error
if (length(args)!=2) {
  stop("Please, provide a first cvs file of profiles and a second csv file of types\n
       USAGE: Rscript PairwiseDifferencesToROC.R Profiles.csv Types.csv", call.=FALSE)
}
## read dataframe of profiles (i.e. Profiles.csv)
## S stands for sample (i.e. rows): n = 12
## L stands for locus (i.e. columns): n = 15
## A stands for allele (i.e. data): n= 180
dfp = read.table(args[1], dec = ".", header=TRUE, sep = ",", quote = "")
## make sure that each variable of the dataframe is a character
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
#### add sample names to matrix columns
colnames(mat) <- sample
### export the matrix into a csv file
write.csv(mat,file="PairwiseMatrix.csv")

# transform the matrix of pairwise differences into a dataframe
## transform the matrix into a dist object
distobj = as.dist(mat)
## transform the dist object into a long dataframe of pairwise differences
dfl = dist2list(distobj)
#### rename variables
names(dfl)[names(dfl) == "col"] <- "FirstSample"
names(dfl)[names(dfl) == "row"] <- "SecondSample"
names(dfl)[names(dfl) == "value"] <- "Differences"

# perform a Receiver Operating Characteristic (ROC) analysis from the dataframe of pairwise differences
## suppose that the samples S1, S2, S3, S4, S5, S6 and S7 are positive controls (PC) of an outbreak
## suppose that the samples S8, S9, S10, S11 and S12 are negative controls (NC) of an outbreak
## suppose that the samples S13, S14, S15, S16, S17, S18, S19 and S20 are tested samples (TS)
## complete the file Types.csv
## prepare dataframes of positive controls (PC), negative controls (NC) and tested samples (TS)
### read the dataframe of samples types
dft <- read.table(args[2], dec = ".", header=TRUE, sep = ",", quote = "")
### subset PC
dfPC <- subset(dft,dft$Type %in% c("PC"))
### subset NC
dfNC <- subset(dft,dft$Type %in% c("NC"))
### subset TS
dfTS <- subset(dft,dft$Type %in% c("TS"))
## derive variables "FirstSample" and "SecondSample" into variables "FirstType" and "SecondType" flagging positive controls, (PC), negative controls (NC) and tested samples (TS)
### derivation of the "FirstSample" column
dfl$FirstType <- ifelse(dfl$FirstSample %in% dfPC$Sample, "PC",
                        ifelse(dfl$FirstSample %in% dfNC$Sample, "NC",
                               ifelse(dfl$FirstSample %in% dfTS$Sample, "TS",
                                      "error")))
### derivation of the "row" column
dfl$SecondType <- ifelse(dfl$SecondSample %in% dfPC$Sample, "PC",
                         ifelse(dfl$SecondSample %in% dfNC$Sample, "NC",
                                ifelse(dfl$SecondSample %in% dfTS$Sample, "TS",
                                       "error")))
## keep in mind that pairwise differences of outbreak related and unrelated samples must be compared and flag the elements below:
### the outbreak related pairs (PC versus PC)
### the extra outbreak pairs (NC versus NC)
### the outbreak unrelated pairs (PC versus NC OR NC versus PC)
### the test related pairs (TS versus TS OR TS versus PC OR PC versus TS OR TS versus NC OR NC versus TS
dfl$Status <- 
  ifelse((dfl$FirstType == "PC") & (dfl$SecondType == "PC"), "related",
         ifelse((dfl$FirstType == "NC") & (dfl$SecondType == "NC"), "extra",
                ifelse((dfl$FirstType == "NC") & (dfl$SecondType == "PC"), "unrelated",
                       ifelse((dfl$FirstType == "PC") & (dfl$SecondType == "NC"), "unrelated",
                              ifelse((dfl$FirstType == "TS") & (dfl$SecondType == "TS"), "test",
                                     ifelse((dfl$FirstType == "TS") & (dfl$SecondType == "PC"), "test",
                                            ifelse((dfl$FirstType == "PC") & (dfl$SecondType == "TS"), "test", 
                                                   ifelse((dfl$FirstType == "TS") & (dfl$SecondType == "NC"), "test", 
                                                          ifelse((dfl$FirstType == "NC") & (dfl$SecondType == "TS"), "test", 
                                                                 "error")))))))))
### export the dataframe into a csv file
write.csv(dfl,file="PairwiseDataframe.csv")
## subset dataframe for ROC analysis (i.e. related and unrelated status)
dflROC <- subset(dfl,dfl$Status %in% c("related","unrelated"))
## run a ROC analysis and plot
### reset graph devices
graphics.off()
### open a pdf file
pdf("ROC.pdf")
### run a analysis
ROC <- roc(dflROC$Status, dflROC$Differences,
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
## rename variables
names(thresholds)[names(thresholds) == "threshold"] <- "Threshold"
names(thresholds)[names(thresholds) == "sensitivity"] <- "Sensitivity"
names(thresholds)[names(thresholds) == "specificity"] <- "Specificity"
## export the thresholds into a csv file
write.csv(thresholds,file="Thresholds.csv")
## extract the threshold presenting the best combination of sensitivity and specificity
best <- coords(ROC, "best", ret=c("threshold", "sensitivity", "specificity"))
## rename variables
names(best)[names(best) == "threshold"] <- "Threshold"
names(best)[names(best) == "sensitivity"] <- "Sensitivity"
names(best)[names(best) == "specificity"] <- "Specificity"
## export the best threshold into a csv file
write.csv(best,file="ThresholdBest.csv")

# perform prediction on tested samples (TS)
## extract the best threshold from the ROC analysis
cutoff  <- best[1,1]
## extract pairwise differences of interest of the tested samples (i.e. TS versus PC from the test status)
### subset dataframe for prediction
dfltest <- subset(dfl,dfl$Status %in% c("test"))
### retain only forward pairs (i.e. without revers pairs)
dflTS <- subset(dfltest,dfltest$SecondType %in% c("PC"))
## flag lower and higher than cutoff
dflTS$Threshold <- 
  ifelse((dflTS$Differences < cutoff), "lower",
         ifelse((dflTS$Differences > cutoff), "higher",
                "error"))
## group by sample and count pairs lower than the threshold
dflower <- aggregate(Threshold ~ FirstSample, dflTS, function(x) sum(x == "lower"))
## rename variables
names(dflower)[names(dflower) == "FirstSample"] <- "TestedSample"
names(dflower)[names(dflower) == "Threshold"] <- "LowerThreshold"
## group by sample and count pairs higher than the threshold
dfhigher <- aggregate(Threshold ~ FirstSample, dflTS, function(x) sum(x == c("higher")))
## rename variables
names(dfhigher)[names(dfhigher) == "FirstSample"] <- "TestedSample"
names(dfhigher)[names(dfhigher) == "Threshold"] <- "HigherThreshold"
## joint dataframes
dfTS <- join(dflower, dfhigher, type = "inner")
## add proportion of lower
dfTS$ProportionLower <- round(((dfTS$LowerThreshold * 100) / (dfTS$LowerThreshold + dfTS$HigherThreshold)), digits = 1)
## add proportion of higher
dfTS$ProportionHigher <- round(((dfTS$HigherThreshold * 100) / (dfTS$LowerThreshold + dfTS$HigherThreshold)), digits = 1)
## add messages of prediction with 80%/20% (related/unrelated)
dfTS$Prediction <- 
  ifelse((dfTS$ProportionLower == 100) & (dfTS$ProportionHigher == 0), "probably related",
         ifelse((dfTS$ProportionLower == 0) & (dfTS$ProportionHigher == 100), "probably unrelated",
                ifelse((dfTS$ProportionLower > 80) & (dfTS$ProportionHigher < 20), "potentially related",
                       ifelse((dfTS$ProportionLower < 80) & (dfTS$ProportionHigher > 20), "potentially unrelated",
                              "ambiguous"))))
## export the predictions into a csv file
write.csv(dfTS,file="Predictions.csv")

# add a message
print("Developped by Nicolas Radomski on January 06 (2022) with the R version 4.1.2 (2021-11-01)")
