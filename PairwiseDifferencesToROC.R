#### Pairwise Differences to Receiving Operating Curve (ROC) Analysis ####

# clean environment
rm(list=ls())

# set working directory for Linux and Mac
setwd("/home/IZSNT/n.radomski/Documents/RstudioWorkingDirectory/EFSA-EOI-EFSA-SCIENCE-2020-01")

# install regular packages
install.packages("data.table")
install.packages("usedist")
install.packages("spaa")
install.packages("pROC")

# call library
library(data.table)
library(usedist)
library(spaa)
library(pROC)

# read a matrix as a dataframe
dfmat = read.table("PairwiseDifferences.csv", dec = ".", header=TRUE, sep = ",", quote = "")

# check class (i.e. data.frame)
class(dfmat)

# keep the labels of sample identifiers into a vector
labels = dfmat$sample

# remove the first column of sample identifiers from the dataframe
dfmat$sample <- NULL

# transform the dataframe as a matrix
mat = as.matrix(dfmat, rownames = TRUE)

# check class (i.e. matrix)
class(mat)

# transform the matrix as a dist object
distobj = as.dist(mat)

# check class (i.e. dist)
class(distobj)

# add the labels (i.e. sample identifiers) to the dist object
distobjlabels = dist_setNames(distobj, labels)

# transform the dist object as long dataframe
dflong = dist2list(distobjlabels)

# check class (i.e. data.frame)
class(dflong)

# suppose that the samples A, B, C, D and E are positive controls of an outbreak
# suppose that samples F, G, H, I and J are negative controls of an outbreak
# derive variables "col" and "row" into variables "new col" and "new row" flagging positive (PC) and negative (NC) controls
dflong$newcol <- 
ifelse((dflong$col == "A"), "PC",
ifelse((dflong$col == "B"), "PC",
ifelse((dflong$col == "C"), "PC",
ifelse((dflong$col == "D"), "PC",
ifelse((dflong$col == "E"), "PC", "NC")))))
dflong$newrow <- 
ifelse((dflong$row == "A"), "PC",
ifelse((dflong$row == "B"), "PC",
ifelse((dflong$row == "C"), "PC",
ifelse((dflong$row == "D"), "PC",
ifelse((dflong$row == "E"), "PC", "NC")))))

# keep in mind that pairwise differences of outbreak related and unrelated pairs must be compared and flag the elements below:
## the outbreak related pairs (PC versus PC)
## the extra outbreak pairs (NC versus NC)
## the outbreak unrelated pairs (PC versus NC OR NC versus PC)
dflong$status <- 
ifelse((dflong$newcol == "PC") & (dflong$newrow == "PC"), "related",
ifelse((dflong$newcol == "NC") & (dflong$newrow == "NC"), "extra",
ifelse((dflong$newcol == "NC") & (dflong$newrow == "PC"), "unrelated",
ifelse((dflong$newcol == "PC") & (dflong$newrow == "NC"), "unrelated", "error"))))

# subset the related and unrelated pairs
## check dimension ([1] 100   6)
dim(dflong)
## subset
dflongsub <- subset(dflong,dflong$status %in% c("related","unrelated"))
## check dimension ([1] 75  6)
dim(dflongsub)

# run a ROC analysis and plot
## reset graph devices
graphics.off()
## open a pdf file
pdf("ROC.pdf")
## run a analysis
ROC <- roc(dflongsub$status, dflongsub$value,
            percent=TRUE,
            ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
            plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
            print.auc=TRUE, show.thres=TRUE)
## add confidence intervals for sensitivity and specificity into the plot
sens.ci <- ci.se(ROC, specificities=seq(0, 100, 5))
spec.ci <- ci.sp(ROC, sensitivities=seq(0, 100, 5))
plot(sens.ci, type="shape", col="blue")
plot(spec.ci, type="shape", col="green")
## close the pdf file
dev.off()

# extract statistical features of the ROC analysis
## area under the curve of raw data (91.04%)
auc(ROC)
## area under the curve of smooth data (90.46%)
smooth(ROC)
# variance (10.48211)
var(ROC)

# extract sensitivity and specificity of all thresholds of pairwise differences
coords(ROC, ret=c("threshold", "sensitivity", "specificity"))
comment <- scan(what="character")
threshold sensitivity specificity
1       -Inf         100           0
2        1.0         100          20
3        3.0         100          28
4        4.5          96          52
5        5.5          88          68
6        6.5          88          76
7        9.0          88          84
8       11.5          84          84
9       13.5          76          84
10      16.0          72          92
11      17.5          68          92
12      18.5          64          92
13      19.5          60          92
14      21.5          60         100
15      23.5          56         100
16      24.5          48         100
17      25.5          44         100
18      30.0          40         100
19      39.5          36         100
20      45.5          28         100
21      50.0          24         100
22      59.0          20         100
23      64.5          16         100
24      71.5          12         100
25      81.0           4         100
26       Inf           0         100

rm(comment)

# extract the threshold presenting the best combination of sensitivity and specificity
coords(ROC, "best", ret=c("threshold", "sensitivity", "specificity"))
# threshold sensitivity specificity
# 9          88          84