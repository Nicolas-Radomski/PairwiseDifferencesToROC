#### profiles of microbial mutations to Receiver Operating Characteristic (ROC) analysis ####

# installation of packages

## clean environment
rm(list=ls())

## install regular packages
install.packages("data.table")
install.packages("spaa")
install.packages("pROC")
install.packages("plyr")

## call library
library(data.table)
library(spaa)
library(pROC)
library(plyr)

# transform profiles into a matrix of pairwise differences

## set working directory for Linux and Mac
setwd("/home/IZSNT/n.radomski/Documents/RstudioWorkingDirectory/PairwiseDifferencesToROC-20220118")

## read dataframe of profiles (i.e. Profiles.csv)
## S stands for sample (i.e. rows): n = 12
## L stands for locus (i.e. columns): n = 15
## A stands for allele (i.e. data): n= 180
dfp = read.table("Profiles.csv", dec = ".", header=TRUE, sep = ",", quote = "")

## make sure that each variable of the dataframe is a character
dfp = data.frame(lapply(dfp, as.character))

## check nature of variables (must be character for each variable)
str(dfp)

## check dimension (i.e. [1] 20 16)
dim(dfp)

## check 20 first lines
head(dfp, 20)
comment <- scan(what="character")
   sample  L1  L2  L3  L4  L5  L6  L7  L8   L9  L10 L11 L12 L13 L14 L15
1      S1 A20 A15 A55 A12 A30 A11 A24 A66  A12  A55 A66  A5 A86 A54 A47
2      S2 A20 A15 A55 A12 A30 A11 A24 A66  A12  A55 A66  A2 A87 A54 A47
3      S3 A20 A15 A55 A12 A30 A11 A24 A66  A12  A55 A66  A2 A86 A54 A47
4      S4 A20 A31 A55 A30 A30 A11 A55 A66  A55  A55 A66  A5 A87 A54 A47
5      S5 A20 A31 A55 A30 A30 A11 A55 A66  A55  A55 A66  A5 A98 A54 A47
6      S6 A10 A15 A10 A12 A30 A10 A24 A10  A12  A55 A66  A5 A98 A54 A47
7      S7 A41 A22 A41 A22 A22 A41 A27 A41  A27  A27 A66  A9 A86 A54 A47
8      S8 A41 A22 A41 A22 A22 A41 A27 A41  A27  A27 A66  A9 A86 A54 A47
9      S9 A41 A15 A41 A12 A30 A41 A24 A41  A12  A55 A66  A8 A97 A54 A47
10    S10 A50 A22 A50 A22 A55 A51 A27 A50  A27  A66 A66  A8 A97 A54 A47
11    S11 A10 A54 A15 A41 A65 A88 A75 A89 A420 A998 A66  A5 A86 A11 A10
12    S12 A10 A54 A98 A41 A65 A88 A75 A89 A420 A998 A66  A8 A86 A14  A1
13    S13 A20 A15 A55 A12 A30 A11 A24 A66  A12  A55 A66  A2 A87 A54 A47
14    S14 A41 A15 A41 A12 A30 A41 A24 A41  A12  A55 A66  A8 A97 A54 A47
15    S15 A20 A15 A55 A12 A30 A11 A24 A66  A12  A55 A66  A5 A86 A54 A47
16    S16 A10 A54 A15 A41 A65 A88 A75 A89 A420 A998 A66  A5 A86 A11 A10
17    S17 A20 A31 A55 A30 A30 A11 A55 A66  A55  A55 A66  A5 A87 A54 A47
18    S18 A41 A22 A41 A22 A22 A41 A27 A41  A27  A27 A66  A9 A86 A54 A47
19    S19 A41 A22 A41 A22 A22 A41 A27 A41  A27  A27 A66  A9 A86 A54 A47
20    S20 A10 A54 A98 A41 A65 A88 A75 A89 A420 A998 A66  A8 A86 A14  A1

rm(comment)

## transpose dataframe
tdfp <- transpose(dfp, keep.names = "locus", make.names = "sample")

## check nature of variables (must be character for each variable)
str(tdfp)

## check dimension (i.e. [1] 15 21)
dim(tdfp)

## check 20 first lines
head(tdfp, 20)
comment <- scan(what="character")
   locus  S1  S2  S3  S4  S5  S6  S7  S8  S9 S10  S11  S12 S13 S14 S15  S16 S17 S18 S19  S20
1     L1 A20 A20 A20 A20 A20 A10 A41 A41 A41 A50  A10  A10 A20 A41 A20  A10 A20 A41 A41  A10
2     L2 A15 A15 A15 A31 A31 A15 A22 A22 A15 A22  A54  A54 A15 A15 A15  A54 A31 A22 A22  A54
3     L3 A55 A55 A55 A55 A55 A10 A41 A41 A41 A50  A15  A98 A55 A41 A55  A15 A55 A41 A41  A98
4     L4 A12 A12 A12 A30 A30 A12 A22 A22 A12 A22  A41  A41 A12 A12 A12  A41 A30 A22 A22  A41
5     L5 A30 A30 A30 A30 A30 A30 A22 A22 A30 A55  A65  A65 A30 A30 A30  A65 A30 A22 A22  A65
6     L6 A11 A11 A11 A11 A11 A10 A41 A41 A41 A51  A88  A88 A11 A41 A11  A88 A11 A41 A41  A88
7     L7 A24 A24 A24 A55 A55 A24 A27 A27 A24 A27  A75  A75 A24 A24 A24  A75 A55 A27 A27  A75
8     L8 A66 A66 A66 A66 A66 A10 A41 A41 A41 A50  A89  A89 A66 A41 A66  A89 A66 A41 A41  A89
9     L9 A12 A12 A12 A55 A55 A12 A27 A27 A12 A27 A420 A420 A12 A12 A12 A420 A55 A27 A27 A420
10   L10 A55 A55 A55 A55 A55 A55 A27 A27 A55 A66 A998 A998 A55 A55 A55 A998 A55 A27 A27 A998
11   L11 A66 A66 A66 A66 A66 A66 A66 A66 A66 A66  A66  A66 A66 A66 A66  A66 A66 A66 A66  A66
12   L12  A5  A2  A2  A5  A5  A5  A9  A9  A8  A8   A5   A8  A2  A8  A5   A5  A5  A9  A9   A8
13   L13 A86 A87 A86 A87 A98 A98 A86 A86 A97 A97  A86  A86 A87 A97 A86  A86 A87 A86 A86  A86
14   L14 A54 A54 A54 A54 A54 A54 A54 A54 A54 A54  A11  A14 A54 A54 A54  A11 A54 A54 A54  A14
15   L15 A47 A47 A47 A47 A47 A47 A47 A47 A47 A47  A10   A1 A47 A47 A47  A10 A47 A47 A47   A1

rm(comment)

## calculate pairwise differences between two vectors (i.e. sum of reversed Boolean vector)
sum(!tdfp$S1 == tdfp$S1) # [1] 0
sum(!tdfp$S1 == tdfp$S2) # [1] 0
sum(!tdfp$S1 == tdfp$S3) # [1] 0
sum(!tdfp$S1 == tdfp$S4) # [1] 4
sum(!tdfp$S1 == tdfp$S5) # [1] 4
sum(!tdfp$S1 == tdfp$S6) # [1] 4
sum(!tdfp$S1 == tdfp$S7) # [1] 10
sum(!tdfp$S1 == tdfp$S8) # [1] 10
sum(!tdfp$S1 == tdfp$S9) # [1] 4
sum(!tdfp$S1 == tdfp$S10) # [1] 10
sum(!tdfp$S1 == tdfp$S11) # [1] 12
sum(!tdfp$S1 == tdfp$S12) # [1] 13
sum(!tdfp$S1 == tdfp$S13) # [1] 2
sum(!tdfp$S1 == tdfp$S14) # [1] 6
sum(!tdfp$S1 == tdfp$S15) # [1] 0
sum(!tdfp$S1 == tdfp$S16) # [1] 12
sum(!tdfp$S1 == tdfp$S17) # [1] 12
sum(!tdfp$S1 == tdfp$S18) # [1] 11
sum(!tdfp$S1 == tdfp$S19) # [1] 11
sum(!tdfp$S1 == tdfp$S20) # [1] 13

## calculate pairwise differences between one vector and all (i.e. dependently of the sample amount)
### S1 versus all
S1 <- integer() # empty output vector
for(i in tdfp[, 2:ncol(tdfp)]) { # for loop on sample colons
  output <- sum(!tdfp$S1 == i) # pairwise differences against S1 vector
  S1 <- c(S1, output) # add for loop output into the empty output vector 
}
print(S1) # [1]  0  2  1  5  5  5 11 11  6 12 12 13  2  6  0 12  5 11 11 13
### S2 versus all
S2 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S2 == i)
  S2 <- c(S2, output)
}
print(S2) # [1]  2  0  1  5  6  6 12 12  6 12 14 14  0  6  2 14  5 12 12 14
### S3 versus all
S3 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S3 == i)
  S3 <- c(S3, output)
}
print(S3) # [1]  1  1  0  6  6  6 11 11  6 12 13 13  1  6  1 13  6 11 11 13
### S4 versus all
S4 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S4 == i)
  S4 <- c(S4, output)
}
print(S4) # [1]  5  5  6  0  1  9 12 12 10 12 13 14  5 10  5 13  0 12 12 14
### S5 versus all
S5 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S5 == i)
  S5 <- c(S5, output)
}
print(S5) # [1]  5  6  6  1  0  8 12 12 10 12 13 14  6 10  5 13  1 12 12 14
### S6 versus all
S6 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S6 == i)
  S6 <- c(S6, output)
}
print(S6) # [1]  5  6  6  9  8  0 12 12  6 12 12 13  6  6  5 12  9 12 12 13
### S7 versus all
S7 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S7 == i)
  S7 <- c(S7, output)
}
print(S7) # [1] 11 12 11 12 12 12  0  0  8  8 13 13 12  8 11 13 12  0  0 13
### S8 versus all
S8 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S8 == i)
  S8 <- c(S8, output)
}
print(S8) # [1] 11 12 11 12 12 12  0  0  8  8 13 13 12  8 11 13 12  0  0 13
### S9 versus all
S9 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S9 == i)
  S9 <- c(S9, output)
}
print(S9) # [1]  6  6  6 10 10  6  8  8  0 10 14 13  6  0  6 14 10  8  8 13
### S10 versus all
S10 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S10 == i)
  S10 <- c(S10, output)
}
print(S10) # [1] 12 12 12 12 12 12  8  8 10  0 14 13 12 10 12 14 12  8  8 13
### S11 versus others
S11 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S11 == i)
  S11 <- c(S11, output)
}
print(S11) # [1] 12 14 13 13 13 12 13 13 14 14  0  4 14 14 12  0 13 13 13  4
### S12 versus all
S12 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S12 == i)
  S12 <- c(S12, output)
}
print(S12) # [1] 13 14 13 14 14 13 13 13 13 13  4  0 14 13 13  4 14 13 13  0
### S13 versus all
S13 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S13 == i)
  S13 <- c(S13, output)
}
print(S13) # [1]  2  0  1  5  6  6 12 12  6 12 14 14  0  6  2 14  5 12 12 14
### S14 versus all
S14 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S14 == i)
  S14 <- c(S14, output)
}
print(S14) # [1]  6  6  6 10 10  6  8  8  0 10 14 13  6  0  6 14 10  8  8 13
### S15 versus all
S15 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S15 == i)
  S15 <- c(S15, output)
}
print(S15) # [1]  0  2  1  5  5  5 11 11  6 12 12 13  2  6  0 12  5 11 11 13
### S16 versus all
S16 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S16 == i)
  S16 <- c(S16, output)
}
print(S16) # [1] 12 14 13 13 13 12 13 13 14 14  0  4 14 14 12  0 13 13 13  4
### S17 versus all
S17 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S17 == i)
  S17 <- c(S17, output)
}
print(S17) # [1]  5  5  6  0  1  9 12 12 10 12 13 14  5 10  5 13  0 12 12 14
### S18 versus all
S18 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S18 == i)
  S18 <- c(S18, output)
}
print(S18) # [1] 11 12 11 12 12 12  0  0  8  8 13 13 12  8 11 13 12  0  0 13
### S19 versus all
S19 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S19 == i)
  S19 <- c(S19, output)
}
print(S19) # [1] 11 12 11 12 12 12  0  0  8  8 13 13 12  8 11 13 12  0  0 13
### S20 versus all
S20 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S20 == i)
  S20 <- c(S20, output)
}
print(S20) # [1] 13 14 13 14 14 13 13 13 13 13  4  0 14 13 13  4 14 13 13  0

### combine vectors into dataframe
#### retrieve samples
sample <- dfp$sample
#### combine as dataframe
pwd <- data.frame(sample, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17, S18, S19, S20)
#### check pairwise differences
pwd
comment <- scan(what="character")
   sample S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19 S20
1      S1  0  2  1  5  5  5 11 11  6  12  12  13   2   6   0  12   5  11  11  13
2      S2  2  0  1  5  6  6 12 12  6  12  14  14   0   6   2  14   5  12  12  14
3      S3  1  1  0  6  6  6 11 11  6  12  13  13   1   6   1  13   6  11  11  13
4      S4  5  5  6  0  1  9 12 12 10  12  13  14   5  10   5  13   0  12  12  14
5      S5  5  6  6  1  0  8 12 12 10  12  13  14   6  10   5  13   1  12  12  14
6      S6  5  6  6  9  8  0 12 12  6  12  12  13   6   6   5  12   9  12  12  13
7      S7 11 12 11 12 12 12  0  0  8   8  13  13  12   8  11  13  12   0   0  13
8      S8 11 12 11 12 12 12  0  0  8   8  13  13  12   8  11  13  12   0   0  13
9      S9  6  6  6 10 10  6  8  8  0  10  14  13   6   0   6  14  10   8   8  13
10    S10 12 12 12 12 12 12  8  8 10   0  14  13  12  10  12  14  12   8   8  13
11    S11 12 14 13 13 13 12 13 13 14  14   0   4  14  14  12   0  13  13  13   4
12    S12 13 14 13 14 14 13 13 13 13  13   4   0  14  13  13   4  14  13  13   0
13    S13  2  0  1  5  6  6 12 12  6  12  14  14   0   6   2  14   5  12  12  14
14    S14  6  6  6 10 10  6  8  8  0  10  14  13   6   0   6  14  10   8   8  13
15    S15  0  2  1  5  5  5 11 11  6  12  12  13   2   6   0  12   5  11  11  13
16    S16 12 14 13 13 13 12 13 13 14  14   0   4  14  14  12   0  13  13  13   4
17    S17  5  5  6  0  1  9 12 12 10  12  13  14   5  10   5  13   0  12  12  14
18    S18 11 12 11 12 12 12  0  0  8   8  13  13  12   8  11  13  12   0   0  13
19    S19 11 12 11 12 12 12  0  0  8   8  13  13  12   8  11  13  12   0   0  13
20    S20 13 14 13 14 14 13 13 13 13  13   4   0  14  13  13   4  14  13  13   0

rm(comment)

### combine vectors into matrix
#### retrieve data
pwddata <- c(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17, S18, S19, S20)
#### add data into matrix
pwdmat <- matrix(pwddata,nrow=20,ncol=20,byrow=TRUE)
#### check pairwise differences
pwdmat
comment <- scan(what="character")
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20]
[1,]    0    2    1    5    5    5   11   11    6    12    12    13     2     6     0    12     5    11    11    13
[2,]    2    0    1    5    6    6   12   12    6    12    14    14     0     6     2    14     5    12    12    14
[3,]    1    1    0    6    6    6   11   11    6    12    13    13     1     6     1    13     6    11    11    13
[4,]    5    5    6    0    1    9   12   12   10    12    13    14     5    10     5    13     0    12    12    14
[5,]    5    6    6    1    0    8   12   12   10    12    13    14     6    10     5    13     1    12    12    14
[6,]    5    6    6    9    8    0   12   12    6    12    12    13     6     6     5    12     9    12    12    13
[7,]   11   12   11   12   12   12    0    0    8     8    13    13    12     8    11    13    12     0     0    13
[8,]   11   12   11   12   12   12    0    0    8     8    13    13    12     8    11    13    12     0     0    13
[9,]    6    6    6   10   10    6    8    8    0    10    14    13     6     0     6    14    10     8     8    13
[10,]   12   12   12   12   12   12    8    8   10     0    14    13    12    10    12    14    12     8     8    13
[11,]   12   14   13   13   13   12   13   13   14    14     0     4    14    14    12     0    13    13    13     4
[12,]   13   14   13   14   14   13   13   13   13    13     4     0    14    13    13     4    14    13    13     0
[13,]    2    0    1    5    6    6   12   12    6    12    14    14     0     6     2    14     5    12    12    14
[14,]    6    6    6   10   10    6    8    8    0    10    14    13     6     0     6    14    10     8     8    13
[15,]    0    2    1    5    5    5   11   11    6    12    12    13     2     6     0    12     5    11    11    13
[16,]   12   14   13   13   13   12   13   13   14    14     0     4    14    14    12     0    13    13    13     4
[17,]    5    5    6    0    1    9   12   12   10    12    13    14     5    10     5    13     0    12    12    14
[18,]   11   12   11   12   12   12    0    0    8     8    13    13    12     8    11    13    12     0     0    13
[19,]   11   12   11   12   12   12    0    0    8     8    13    13    12     8    11    13    12     0     0    13
[20,]   13   14   13   14   14   13   13   13   13    13     4     0    14    13    13     4    14    13    13     0

rm(comment)

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
### check the output vector v
print(v)
comment <- scan(what="character")
  [1]  0  2  1  5  5  5 11 11  6 12 12 13  2  6  0 12  5 11 11 13
 [21]  2  0  1  5  6  6 12 12  6 12 14 14  0  6  2 14  5 12 12 14
 [41]  1  1  0  6  6  6 11 11  6 12 13 13  1  6  1 13  6 11 11 13
 [61]  5  5  6  0  1  9 12 12 10 12 13 14  5 10  5 13  0 12 12 14
 [81]  5  6  6  1  0  8 12 12 10 12 13 14  6 10  5 13  1 12 12 14
[101]  5  6  6  9  8  0 12 12  6 12 12 13  6  6  5 12  9 12 12 13
[121] 11 12 11 12 12 12  0  0  8  8 13 13 12  8 11 13 12  0  0 13
[141] 11 12 11 12 12 12  0  0  8  8 13 13 12  8 11 13 12  0  0 13
[161]  6  6  6 10 10  6  8  8  0 10 14 13  6  0  6 14 10  8  8 13
[181] 12 12 12 12 12 12  8  8 10  0 14 13 12 10 12 14 12  8  8 13
[201] 12 14 13 13 13 12 13 13 14 14  0  4 14 14 12  0 13 13 13  4
[221] 13 14 13 14 14 13 13 13 13 13  4  0 14 13 13  4 14 13 13  0
[241]  2  0  1  5  6  6 12 12  6 12 14 14  0  6  2 14  5 12 12 14
[261]  6  6  6 10 10  6  8  8  0 10 14 13  6  0  6 14 10  8  8 13
[281]  0  2  1  5  5  5 11 11  6 12 12 13  2  6  0 12  5 11 11 13
[301] 12 14 13 13 13 12 13 13 14 14  0  4 14 14 12  0 13 13 13  4
[321]  5  5  6  0  1  9 12 12 10 12 13 14  5 10  5 13  0 12 12 14
[341] 11 12 11 12 12 12  0  0  8  8 13 13 12  8 11 13 12  0  0 13
[361] 11 12 11 12 12 12  0  0  8  8 13 13 12  8 11 13 12  0  0 13
[381] 13 14 13 14 14 13 13 13 13 13  4  0 14 13 13  4 14 13 13  0

rm(comment)

### check size of the vector v (i.e. [1] 144)
length(v)
### create a matrix adding the vector v by row independently of the sample amount
mat <- matrix(v, nrow = (ncol(tdfp)-1) , ncol = (ncol(tdfp)-1), byrow = TRUE)
### add row and column names of the matrix
#### keep the names of sample identifiers into a vector
sample = dfp$sample
#### add sample names to matrix rows
rownames(mat) <- sample
#### add sample names to matrix columns
colnames(mat) <- sample
### check class (i.e. [1] "matrix" "array")
class(mat)
### check pairwise differences
mat
comment <- scan(what="character")
    S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19 S20
S1   0  2  1  5  5  5 11 11  6  12  12  13   2   6   0  12   5  11  11  13
S2   2  0  1  5  6  6 12 12  6  12  14  14   0   6   2  14   5  12  12  14
S3   1  1  0  6  6  6 11 11  6  12  13  13   1   6   1  13   6  11  11  13
S4   5  5  6  0  1  9 12 12 10  12  13  14   5  10   5  13   0  12  12  14
S5   5  6  6  1  0  8 12 12 10  12  13  14   6  10   5  13   1  12  12  14
S6   5  6  6  9  8  0 12 12  6  12  12  13   6   6   5  12   9  12  12  13
S7  11 12 11 12 12 12  0  0  8   8  13  13  12   8  11  13  12   0   0  13
S8  11 12 11 12 12 12  0  0  8   8  13  13  12   8  11  13  12   0   0  13
S9   6  6  6 10 10  6  8  8  0  10  14  13   6   0   6  14  10   8   8  13
S10 12 12 12 12 12 12  8  8 10   0  14  13  12  10  12  14  12   8   8  13
S11 12 14 13 13 13 12 13 13 14  14   0   4  14  14  12   0  13  13  13   4
S12 13 14 13 14 14 13 13 13 13  13   4   0  14  13  13   4  14  13  13   0
S13  2  0  1  5  6  6 12 12  6  12  14  14   0   6   2  14   5  12  12  14
S14  6  6  6 10 10  6  8  8  0  10  14  13   6   0   6  14  10   8   8  13
S15  0  2  1  5  5  5 11 11  6  12  12  13   2   6   0  12   5  11  11  13
S16 12 14 13 13 13 12 13 13 14  14   0   4  14  14  12   0  13  13  13   4
S17  5  5  6  0  1  9 12 12 10  12  13  14   5  10   5  13   0  12  12  14
S18 11 12 11 12 12 12  0  0  8   8  13  13  12   8  11  13  12   0   0  13
S19 11 12 11 12 12 12  0  0  8   8  13  13  12   8  11  13  12   0   0  13
S20 13 14 13 14 14 13 13 13 13  13   4   0  14  13  13   4  14  13  13   0

rm(comment)

### export the matrix into a csv file
write.csv(mat,file="PairwiseMatrix.csv")

# transform the matrix of pairwise differences into a dataframe

## transform the matrix into a dist object
distobj = as.dist(mat)

## check class (i.e. [1] "dist")
class(distobj)

## check the dist object
distobj
comment <- scan(what="character")
    S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19
S2   2                                                                
S3   1  1                                                             
S4   5  5  6                                                          
S5   5  6  6  1                                                       
S6   5  6  6  9  8                                                    
S7  11 12 11 12 12 12                                                 
S8  11 12 11 12 12 12  0                                              
S9   6  6  6 10 10  6  8  8                                           
S10 12 12 12 12 12 12  8  8 10                                        
S11 12 14 13 13 13 12 13 13 14  14                                    
S12 13 14 13 14 14 13 13 13 13  13   4                                
S13  2  0  1  5  6  6 12 12  6  12  14  14                            
S14  6  6  6 10 10  6  8  8  0  10  14  13   6                        
S15  0  2  1  5  5  5 11 11  6  12  12  13   2   6                    
S16 12 14 13 13 13 12 13 13 14  14   0   4  14  14  12                
S17  5  5  6  0  1  9 12 12 10  12  13  14   5  10   5  13            
S18 11 12 11 12 12 12  0  0  8   8  13  13  12   8  11  13  12        
S19 11 12 11 12 12 12  0  0  8   8  13  13  12   8  11  13  12   0    
S20 13 14 13 14 14 13 13 13 13  13   4   0  14  13  13   4  14  13  13

rm(comment)

## transform the dist object into a long dataframe of pairwise differences
dfl = dist2list(distobj)

## check class (i.e. [1] "data.frame")
class(dfl)

## check dimension (i.e. [1] 400   3)
dim(dfl)

## check names of variables
colnames(dfl)

#### rename variables
names(dfl)[names(dfl) == "col"] <- "FirstSample"
names(dfl)[names(dfl) == "row"] <- "SecondSample"
names(dfl)[names(dfl) == "value"] <- "Differences"

## check the long dataframe
dfl
comment <- scan(what="character")
    FirstSample SecondSample Differences
1            S1           S1           0
2            S2           S1           2
3            S3           S1           1
4            S4           S1           5
5            S5           S1           5
6            S6           S1           5
7            S7           S1          11
8            S8           S1          11
9            S9           S1           6
10          S10           S1          12
11          S11           S1          12
12          S12           S1          13
13          S13           S1           2
14          S14           S1           6
15          S15           S1           0
16          S16           S1          12
17          S17           S1           5
18          S18           S1          11
19          S19           S1          11
20          S20           S1          13
21           S1           S2           2
22           S2           S2           0
23           S3           S2           1
24           S4           S2           5
25           S5           S2           6
26           S6           S2           6
27           S7           S2          12
28           S8           S2          12
29           S9           S2           6
30          S10           S2          12
31          S11           S2          14
32          S12           S2          14
33          S13           S2           0
34          S14           S2           6
35          S15           S2           2
36          S16           S2          14
37          S17           S2           5
38          S18           S2          12
39          S19           S2          12
40          S20           S2          14
41           S1           S3           1
42           S2           S3           1
43           S3           S3           0
44           S4           S3           6
45           S5           S3           6
46           S6           S3           6
47           S7           S3          11
48           S8           S3          11
49           S9           S3           6
50          S10           S3          12
51          S11           S3          13
52          S12           S3          13
53          S13           S3           1
54          S14           S3           6
55          S15           S3           1
56          S16           S3          13
57          S17           S3           6
58          S18           S3          11
59          S19           S3          11
60          S20           S3          13
61           S1           S4           5
62           S2           S4           5
63           S3           S4           6
64           S4           S4           0
65           S5           S4           1
66           S6           S4           9
67           S7           S4          12
68           S8           S4          12
69           S9           S4          10
70          S10           S4          12
71          S11           S4          13
72          S12           S4          14
73          S13           S4           5
74          S14           S4          10
75          S15           S4           5
76          S16           S4          13
77          S17           S4           0
78          S18           S4          12
79          S19           S4          12
80          S20           S4          14
81           S1           S5           5
82           S2           S5           6
83           S3           S5           6
84           S4           S5           1
85           S5           S5           0
86           S6           S5           8
87           S7           S5          12
88           S8           S5          12
89           S9           S5          10
90          S10           S5          12
91          S11           S5          13
92          S12           S5          14
93          S13           S5           6
94          S14           S5          10
95          S15           S5           5
96          S16           S5          13
97          S17           S5           1
98          S18           S5          12
99          S19           S5          12
100         S20           S5          14
101          S1           S6           5
102          S2           S6           6
103          S3           S6           6
104          S4           S6           9
105          S5           S6           8
106          S6           S6           0
107          S7           S6          12
108          S8           S6          12
109          S9           S6           6
110         S10           S6          12
111         S11           S6          12
112         S12           S6          13
113         S13           S6           6
114         S14           S6           6
115         S15           S6           5
116         S16           S6          12
117         S17           S6           9
118         S18           S6          12
119         S19           S6          12
120         S20           S6          13
121          S1           S7          11
122          S2           S7          12
123          S3           S7          11
124          S4           S7          12
125          S5           S7          12
126          S6           S7          12
127          S7           S7           0
128          S8           S7           0
129          S9           S7           8
130         S10           S7           8
131         S11           S7          13
132         S12           S7          13
133         S13           S7          12
134         S14           S7           8
135         S15           S7          11
136         S16           S7          13
137         S17           S7          12
138         S18           S7           0
139         S19           S7           0
140         S20           S7          13
141          S1           S8          11
142          S2           S8          12
143          S3           S8          11
144          S4           S8          12
145          S5           S8          12
146          S6           S8          12
147          S7           S8           0
148          S8           S8           0
149          S9           S8           8
150         S10           S8           8
151         S11           S8          13
152         S12           S8          13
153         S13           S8          12
154         S14           S8           8
155         S15           S8          11
156         S16           S8          13
157         S17           S8          12
158         S18           S8           0
159         S19           S8           0
160         S20           S8          13
161          S1           S9           6
162          S2           S9           6
163          S3           S9           6
164          S4           S9          10
165          S5           S9          10
166          S6           S9           6
167          S7           S9           8
168          S8           S9           8
169          S9           S9           0
170         S10           S9          10
171         S11           S9          14
172         S12           S9          13
173         S13           S9           6
174         S14           S9           0
175         S15           S9           6
176         S16           S9          14
177         S17           S9          10
178         S18           S9           8
179         S19           S9           8
180         S20           S9          13
181          S1          S10          12
182          S2          S10          12
183          S3          S10          12
184          S4          S10          12
185          S5          S10          12
186          S6          S10          12
187          S7          S10           8
188          S8          S10           8
189          S9          S10          10
190         S10          S10           0
191         S11          S10          14
192         S12          S10          13
193         S13          S10          12
194         S14          S10          10
195         S15          S10          12
196         S16          S10          14
197         S17          S10          12
198         S18          S10           8
199         S19          S10           8
200         S20          S10          13
201          S1          S11          12
202          S2          S11          14
203          S3          S11          13
204          S4          S11          13
205          S5          S11          13
206          S6          S11          12
207          S7          S11          13
208          S8          S11          13
209          S9          S11          14
210         S10          S11          14
211         S11          S11           0
212         S12          S11           4
213         S13          S11          14
214         S14          S11          14
215         S15          S11          12
216         S16          S11           0
217         S17          S11          13
218         S18          S11          13
219         S19          S11          13
220         S20          S11           4
221          S1          S12          13
222          S2          S12          14
223          S3          S12          13
224          S4          S12          14
225          S5          S12          14
226          S6          S12          13
227          S7          S12          13
228          S8          S12          13
229          S9          S12          13
230         S10          S12          13
231         S11          S12           4
232         S12          S12           0
233         S13          S12          14
234         S14          S12          13
235         S15          S12          13
236         S16          S12           4
237         S17          S12          14
238         S18          S12          13
239         S19          S12          13
240         S20          S12           0
241          S1          S13           2
242          S2          S13           0
243          S3          S13           1
244          S4          S13           5
245          S5          S13           6
246          S6          S13           6
247          S7          S13          12
248          S8          S13          12
249          S9          S13           6
250         S10          S13          12
251         S11          S13          14
252         S12          S13          14
253         S13          S13           0
254         S14          S13           6
255         S15          S13           2
256         S16          S13          14
257         S17          S13           5
258         S18          S13          12
259         S19          S13          12
260         S20          S13          14
261          S1          S14           6
262          S2          S14           6
263          S3          S14           6
264          S4          S14          10
265          S5          S14          10
266          S6          S14           6
267          S7          S14           8
268          S8          S14           8
269          S9          S14           0
270         S10          S14          10
271         S11          S14          14
272         S12          S14          13
273         S13          S14           6
274         S14          S14           0
275         S15          S14           6
276         S16          S14          14
277         S17          S14          10
278         S18          S14           8
279         S19          S14           8
280         S20          S14          13
281          S1          S15           0
282          S2          S15           2
283          S3          S15           1
284          S4          S15           5
285          S5          S15           5
286          S6          S15           5
287          S7          S15          11
288          S8          S15          11
289          S9          S15           6
290         S10          S15          12
291         S11          S15          12
292         S12          S15          13
293         S13          S15           2
294         S14          S15           6
295         S15          S15           0
296         S16          S15          12
297         S17          S15           5
298         S18          S15          11
299         S19          S15          11
300         S20          S15          13
301          S1          S16          12
302          S2          S16          14
303          S3          S16          13
304          S4          S16          13
305          S5          S16          13
306          S6          S16          12
307          S7          S16          13
308          S8          S16          13
309          S9          S16          14
310         S10          S16          14
311         S11          S16           0
312         S12          S16           4
313         S13          S16          14
314         S14          S16          14
315         S15          S16          12
316         S16          S16           0
317         S17          S16          13
318         S18          S16          13
319         S19          S16          13
320         S20          S16           4
321          S1          S17           5
322          S2          S17           5
323          S3          S17           6
324          S4          S17           0
325          S5          S17           1
326          S6          S17           9
327          S7          S17          12
328          S8          S17          12
329          S9          S17          10
330         S10          S17          12
331         S11          S17          13
332         S12          S17          14
333         S13          S17           5
[ reached 'max' / getOption("max.print") -- omitted 67 rows ]

rm(comment)

# perform a Receiver Operating Characteristic (ROC) analysis from the dataframe of pairwise differences
## suppose that the samples S1, S2, S3, S4, S5, S6 and S7 are positive controls (PC) of an outbreak
## suppose that the samples S8, S9, S10, S11 and S12 are negative controls (NC) of an outbreak
## suppose that the samples S13, S14, S15, S16, S17, S18, S19 and S20 are tested samples (TS)
## complete the file Types.csv as below
comment <- scan(what="character")
Sample	Type
S1	PC
S2	PC
S3	PC
S4	PC
S5	PC
S6	PC
S7	PC
S8	NC
S9	NC
S10	NC
S11	NC
S12	NC
S13	TS
S14	TS
S15	TS
S16	TS
S17	TS
S18	TS
S19	TS
S20	TS

rm(comment)

## prepare dataframes of positive controls (PC), negative controls (NC) and tested samples (TS)
### read the dataframe of samples types
dft <- read.table("Types.csv", dec = ".", header=TRUE, sep = ",", quote = "")
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
### derivation of the "SecondSample" column
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

## subset the related and unrelated pairs
### check dimension (i.e. [1] 400   6)
dim(dfl)
### check dataframe
dfl
comment <- scan(what="character")
FirstSample SecondSample Differences FirstType SecondType    Status
1            S1           S1           0        PC         PC   related
2            S2           S1           2        PC         PC   related
3            S3           S1           1        PC         PC   related
4            S4           S1           5        PC         PC   related
5            S5           S1           5        PC         PC   related
6            S6           S1           5        PC         PC   related
7            S7           S1          11        PC         PC   related
8            S8           S1          11        NC         PC unrelated
9            S9           S1           6        NC         PC unrelated
10          S10           S1          12        NC         PC unrelated
11          S11           S1          12        NC         PC unrelated
12          S12           S1          13        NC         PC unrelated
13          S13           S1           2        TS         PC      test
14          S14           S1           6        TS         PC      test
15          S15           S1           0        TS         PC      test
16          S16           S1          12        TS         PC      test
17          S17           S1           5        TS         PC      test
18          S18           S1          11        TS         PC      test
19          S19           S1          11        TS         PC      test
20          S20           S1          13        TS         PC      test
21           S1           S2           2        PC         PC   related
22           S2           S2           0        PC         PC   related
23           S3           S2           1        PC         PC   related
24           S4           S2           5        PC         PC   related
25           S5           S2           6        PC         PC   related
26           S6           S2           6        PC         PC   related
27           S7           S2          12        PC         PC   related
28           S8           S2          12        NC         PC unrelated
29           S9           S2           6        NC         PC unrelated
30          S10           S2          12        NC         PC unrelated
31          S11           S2          14        NC         PC unrelated
32          S12           S2          14        NC         PC unrelated
33          S13           S2           0        TS         PC      test
34          S14           S2           6        TS         PC      test
35          S15           S2           2        TS         PC      test
36          S16           S2          14        TS         PC      test
37          S17           S2           5        TS         PC      test
38          S18           S2          12        TS         PC      test
39          S19           S2          12        TS         PC      test
40          S20           S2          14        TS         PC      test
41           S1           S3           1        PC         PC   related
42           S2           S3           1        PC         PC   related
43           S3           S3           0        PC         PC   related
44           S4           S3           6        PC         PC   related
45           S5           S3           6        PC         PC   related
46           S6           S3           6        PC         PC   related
47           S7           S3          11        PC         PC   related
48           S8           S3          11        NC         PC unrelated
49           S9           S3           6        NC         PC unrelated
50          S10           S3          12        NC         PC unrelated
51          S11           S3          13        NC         PC unrelated
52          S12           S3          13        NC         PC unrelated
53          S13           S3           1        TS         PC      test
54          S14           S3           6        TS         PC      test
55          S15           S3           1        TS         PC      test
56          S16           S3          13        TS         PC      test
57          S17           S3           6        TS         PC      test
58          S18           S3          11        TS         PC      test
59          S19           S3          11        TS         PC      test
60          S20           S3          13        TS         PC      test
61           S1           S4           5        PC         PC   related
62           S2           S4           5        PC         PC   related
63           S3           S4           6        PC         PC   related
64           S4           S4           0        PC         PC   related
65           S5           S4           1        PC         PC   related
66           S6           S4           9        PC         PC   related
67           S7           S4          12        PC         PC   related
68           S8           S4          12        NC         PC unrelated
69           S9           S4          10        NC         PC unrelated
70          S10           S4          12        NC         PC unrelated
71          S11           S4          13        NC         PC unrelated
72          S12           S4          14        NC         PC unrelated
73          S13           S4           5        TS         PC      test
74          S14           S4          10        TS         PC      test
75          S15           S4           5        TS         PC      test
76          S16           S4          13        TS         PC      test
77          S17           S4           0        TS         PC      test
78          S18           S4          12        TS         PC      test
79          S19           S4          12        TS         PC      test
80          S20           S4          14        TS         PC      test
81           S1           S5           5        PC         PC   related
82           S2           S5           6        PC         PC   related
83           S3           S5           6        PC         PC   related
84           S4           S5           1        PC         PC   related
85           S5           S5           0        PC         PC   related
86           S6           S5           8        PC         PC   related
87           S7           S5          12        PC         PC   related
88           S8           S5          12        NC         PC unrelated
89           S9           S5          10        NC         PC unrelated
90          S10           S5          12        NC         PC unrelated
91          S11           S5          13        NC         PC unrelated
92          S12           S5          14        NC         PC unrelated
93          S13           S5           6        TS         PC      test
94          S14           S5          10        TS         PC      test
95          S15           S5           5        TS         PC      test
96          S16           S5          13        TS         PC      test
97          S17           S5           1        TS         PC      test
98          S18           S5          12        TS         PC      test
99          S19           S5          12        TS         PC      test
100         S20           S5          14        TS         PC      test
101          S1           S6           5        PC         PC   related
102          S2           S6           6        PC         PC   related
103          S3           S6           6        PC         PC   related
104          S4           S6           9        PC         PC   related
105          S5           S6           8        PC         PC   related
106          S6           S6           0        PC         PC   related
107          S7           S6          12        PC         PC   related
108          S8           S6          12        NC         PC unrelated
109          S9           S6           6        NC         PC unrelated
110         S10           S6          12        NC         PC unrelated
111         S11           S6          12        NC         PC unrelated
112         S12           S6          13        NC         PC unrelated
113         S13           S6           6        TS         PC      test
114         S14           S6           6        TS         PC      test
115         S15           S6           5        TS         PC      test
116         S16           S6          12        TS         PC      test
117         S17           S6           9        TS         PC      test
118         S18           S6          12        TS         PC      test
119         S19           S6          12        TS         PC      test
120         S20           S6          13        TS         PC      test
121          S1           S7          11        PC         PC   related
122          S2           S7          12        PC         PC   related
123          S3           S7          11        PC         PC   related
124          S4           S7          12        PC         PC   related
125          S5           S7          12        PC         PC   related
126          S6           S7          12        PC         PC   related
127          S7           S7           0        PC         PC   related
128          S8           S7           0        NC         PC unrelated
129          S9           S7           8        NC         PC unrelated
130         S10           S7           8        NC         PC unrelated
131         S11           S7          13        NC         PC unrelated
132         S12           S7          13        NC         PC unrelated
133         S13           S7          12        TS         PC      test
134         S14           S7           8        TS         PC      test
135         S15           S7          11        TS         PC      test
136         S16           S7          13        TS         PC      test
137         S17           S7          12        TS         PC      test
138         S18           S7           0        TS         PC      test
139         S19           S7           0        TS         PC      test
140         S20           S7          13        TS         PC      test
141          S1           S8          11        PC         NC unrelated
142          S2           S8          12        PC         NC unrelated
143          S3           S8          11        PC         NC unrelated
144          S4           S8          12        PC         NC unrelated
145          S5           S8          12        PC         NC unrelated
146          S6           S8          12        PC         NC unrelated
147          S7           S8           0        PC         NC unrelated
148          S8           S8           0        NC         NC     extra
149          S9           S8           8        NC         NC     extra
150         S10           S8           8        NC         NC     extra
151         S11           S8          13        NC         NC     extra
152         S12           S8          13        NC         NC     extra
153         S13           S8          12        TS         NC      test
154         S14           S8           8        TS         NC      test
155         S15           S8          11        TS         NC      test
156         S16           S8          13        TS         NC      test
157         S17           S8          12        TS         NC      test
158         S18           S8           0        TS         NC      test
159         S19           S8           0        TS         NC      test
160         S20           S8          13        TS         NC      test
161          S1           S9           6        PC         NC unrelated
162          S2           S9           6        PC         NC unrelated
163          S3           S9           6        PC         NC unrelated
164          S4           S9          10        PC         NC unrelated
165          S5           S9          10        PC         NC unrelated
166          S6           S9           6        PC         NC unrelated
[ reached 'max' / getOption("max.print") -- omitted 234 rows ]

rm(comment)

### export the dataframe into a csv file
write.csv(dfl,file="PairwiseDataframe.csv")

### subset dataframe for ROC analysis (i.e. related and unrelated status)
dflROC <- subset(dfl,dfl$Status %in% c("related","unrelated"))
## check dimension (i.e. [1] 119  6)
dim(dflROC)

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

## check dataframe
thresholds
comment <- scan(what="character")
   Threshold Sensitivity Specificity
1       -Inf   100.00000     0.00000
2        0.5    97.14286    14.28571
3        1.5    97.14286    26.53061
4        3.5    97.14286    30.61224
5        5.5    97.14286    46.93878
6        7.0    85.71429    67.34694
7        8.5    80.00000    71.42857
8        9.5    80.00000    75.51020
9       10.5    74.28571    75.51020
10      11.5    68.57143    83.67347
11      12.5    34.28571   100.00000
12      13.5    11.42857   100.00000
13       Inf     0.00000   100.00000

rm(comment)

## export the thresholds into a csv file
write.csv(thresholds,file="Thresholds.csv")

## extract the threshold presenting the best combination of sensitivity and specificity
best <- coords(ROC, "best", ret=c("threshold", "sensitivity", "specificity"))

## rename variables
names(best)[names(best) == "threshold"] <- "Threshold"
names(best)[names(best) == "sensitivity"] <- "Sensitivity"
names(best)[names(best) == "specificity"] <- "Specificity"

## check dataframe
best
comment <- scan(what="character")
threshold sensitivity specificity
1       9.5          80     75.5102

rm(comment)

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
### check dataframe
dfTS
comment <- scan(what="character")
TestedSample LowerThreshold HigherThreshold ProportionLower ProportionHigher            Prediction
1          S13              6               1            85.7             14.3   potentially related
2          S14              5               2            71.4             28.6 potentially unrelated
3          S15              6               1            85.7             14.3   potentially related
4          S16              0               7             0.0            100.0    probably unrelated
5          S17              6               1            85.7             14.3   potentially related
6          S18              1               6            14.3             85.7 potentially unrelated
7          S19              1               6            14.3             85.7 potentially unrelated
8          S20              0               7             0.0            100.0    probably unrelated

rm(comment)

## export the predictions into a csv file
write.csv(dfTS,file="Predictions.csv")

# add a message
print("Developped by Nicolas Radomski on January 06 (2022) with the R version 4.1.2 (2021-11-01)")
