#### profiles of microbial mutations to Receiver Operating Characteristic (ROC) analysis ####

# installation of packages

## clean environment
rm(list=ls())

## install regular packages
install.packages("data.table")
install.packages("spaa")
install.packages("pROC")

## call library
library(data.table)
library(spaa)
library(pROC)

# transform profiles into a matrix of pairwise differences

## set working directory for Linux and Mac
setwd("/PATH/TO/RstudioWorkingDirectory")

## read dataframe of profiles (i.e. Profiles.csv)
## S stands for sample (i.e. rows): n = 12
## L stands for locus (i.e. columns): n = 15
## A stands for allele (i.e. data): n= 180
dfp = read.table("Profiles.csv", dec = ".", header=TRUE, sep = ",", quote = "")

## make sure that each variable of the dataframe is a character
dfp = data.frame(lapply(dfp, as.character))

## check nature of variables (must be character for each variable)
str(dfp)

## check dimension (i.e. [1] 12 16)
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

rm(comment)

## transpose dataframe
tdfp <- transpose(dfp, keep.names = "locus", make.names = "sample")

## check nature of variables (must be character for each variable)
str(tdfp)

## check dimension (i.e. [1] 15 13)
dim(tdfp)

## check 20 first lines
head(tdfp, 20)
comment <- scan(what="character")
locus  S1  S2  S3  S4  S5  S6  S7  S8  S9 S10  S11  S12
1     L1 A20 A20 A20 A20 A20 A10 A41 A41 A41 A50  A10  A10
2     L2 A15 A15 A15 A31 A31 A15 A22 A22 A15 A22  A54  A54
3     L3 A55 A55 A55 A55 A55 A10 A41 A41 A41 A50  A15  A98
4     L4 A12 A12 A12 A30 A30 A12 A22 A22 A12 A22  A41  A41
5     L5 A30 A30 A30 A30 A30 A30 A22 A22 A30 A55  A65  A65
6     L6 A11 A11 A11 A11 A11 A10 A41 A41 A41 A51  A88  A88
7     L7 A24 A24 A24 A55 A55 A24 A27 A27 A24 A27  A75  A75
8     L8 A66 A66 A66 A66 A66 A10 A41 A41 A41 A50  A89  A89
9     L9 A12 A12 A12 A55 A55 A12 A27 A27 A12 A27 A420 A420
10   L10 A55 A55 A55 A55 A55 A55 A27 A27 A55 A66 A998 A998
11   L11 A66 A66 A66 A66 A66 A66 A66 A66 A66 A66  A66  A66
12   L12  A5  A2  A2  A5  A5  A5  A9  A9  A8  A8   A5   A8
13   L13 A86 A87 A86 A87 A98 A98 A86 A86 A97 A97  A86  A86
14   L14 A54 A54 A54 A54 A54 A54 A54 A54 A54 A54  A11  A14
15   L15 A47 A47 A47 A47 A47 A47 A47 A47 A47 A47  A10   A1

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

## calculate pairwise differences between one vector and all (i.e. dependently of the sample amount)
### S1 versus all
S1 <- integer() # empty output vector
for(i in tdfp[, 2:ncol(tdfp)]) { # for loop on sample colons
  output <- sum(!tdfp$S1 == i) # pairwise differences against S1 vector
  S1 <- c(S1, output) # add for loop output into the empty output vector 
}
print(S1) # [1]  0  2  1  5  5  5 11 11  6 12 12 13
### S2 versus all
S2 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S2 == i)
  S2 <- c(S2, output)
}
print(S2) # [1]  2  0  1  5  6  6 12 12  6 12 14 14
### S3 versus all
S3 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S3 == i)
  S3 <- c(S3, output)
}
print(S3) # [1]  1  1  0  6  6  6 11 11  6 12 13 13
### S4 versus all
S4 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S4 == i)
  S4 <- c(S4, output)
}
print(S4) # [1]  5  5  6  0  1  9 12 12 10 12 13 14
### S5 versus all
S5 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S5 == i)
  S5 <- c(S5, output)
}
print(S5) # [1]  5  6  6  1  0  8 12 12 10 12 13 14
### S6 versus all
S6 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S6 == i)
  S6 <- c(S6, output)
}
print(S6) # [1]  5  6  6  9  8  0 12 12  6 12 12 13
### S7 versus all
S7 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S7 == i)
  S7 <- c(S7, output)
}
print(S7) # [1] 11 12 11 12 12 12  0  0  8  8 13 13
### S8 versus all
S8 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S8 == i)
  S8 <- c(S8, output)
}
print(S8) # [1] 11 12 11 12 12 12  0  0  8  8 13 13
### S9 versus all
S9 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S9 == i)
  S9 <- c(S9, output)
}
print(S9) # [1]  6  6  6 10 10  6  8  8  0 10 14 13
### S10 versus all
S10 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S10 == i)
  S10 <- c(S10, output)
}
print(S10) # [1] 12 12 12 12 12 12  8  8 10  0 14 13
### S11 versus others
S11 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S11 == i)
  S11 <- c(S11, output)
}
print(S11) # [1] 12 14 13 13 13 12 13 13 14 14  0  4
### S12 versus all
S12 <- integer()
for(i in tdfp[, 2:ncol(tdfp)]) {
  output <- sum(!tdfp$S12 == i)
  S12 <- c(S12, output)
}
print(S12) # [1] 13 14 13 14 14 13 13 13 13 13  4  0

### combine vectors into dataframe
#### retrieve samples
sample <- dfp$sample
#### combine as dataframe
pwd <- data.frame(sample, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12)
#### check pairwise differences
pwd
comment <- scan(what="character")
   sample S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12
1      S1  0  2  1  5  5  5 11 11  6  12  12  13
2      S2  2  0  1  5  6  6 12 12  6  12  14  14
3      S3  1  1  0  6  6  6 11 11  6  12  13  13
4      S4  5  5  6  0  1  9 12 12 10  12  13  14
5      S5  5  6  6  1  0  8 12 12 10  12  13  14
6      S6  5  6  6  9  8  0 12 12  6  12  12  13
7      S7 11 12 11 12 12 12  0  0  8   8  13  13
8      S8 11 12 11 12 12 12  0  0  8   8  13  13
9      S9  6  6  6 10 10  6  8  8  0  10  14  13
10    S10 12 12 12 12 12 12  8  8 10   0  14  13
11    S11 12 14 13 13 13 12 13 13 14  14   0   4
12    S12 13 14 13 14 14 13 13 13 13  13   4   0

rm(comment)

### combine vectors into matrix
#### retrieve data
pwddata <- c(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12)
#### add data into matrix
pwdmat <- matrix(pwddata,nrow=12,ncol=12,byrow=TRUE)
#### check pairwise differences
pwdmat
comment <- scan(what="character")
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
 [1,]    0    2    1    5    5    5   11   11    6    12    12    13
 [2,]    2    0    1    5    6    6   12   12    6    12    14    14
 [3,]    1    1    0    6    6    6   11   11    6    12    13    13
 [4,]    5    5    6    0    1    9   12   12   10    12    13    14
 [5,]    5    6    6    1    0    8   12   12   10    12    13    14
 [6,]    5    6    6    9    8    0   12   12    6    12    12    13
 [7,]   11   12   11   12   12   12    0    0    8     8    13    13
 [8,]   11   12   11   12   12   12    0    0    8     8    13    13
 [9,]    6    6    6   10   10    6    8    8    0    10    14    13
[10,]   12   12   12   12   12   12    8    8   10     0    14    13
[11,]   12   14   13   13   13   12   13   13   14    14     0     4
[12,]   13   14   13   14   14   13   13   13   13    13     4     0

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
  [1]  0  2  1  5  5  5 11 11  6 12 12 13
 [13]  2  0  1  5  6  6 12 12  6 12 14 14
 [25]  1  1  0  6  6  6 11 11  6 12 13 13
 [37]  5  5  6  0  1  9 12 12 10 12 13 14
 [49]  5  6  6  1  0  8 12 12 10 12 13 14
 [61]  5  6  6  9  8  0 12 12  6 12 12 13
 [73] 11 12 11 12 12 12  0  0  8  8 13 13
 [85] 11 12 11 12 12 12  0  0  8  8 13 13
 [97]  6  6  6 10 10  6  8  8  0 10 14 13
[109] 12 12 12 12 12 12  8  8 10  0 14 13
[121] 12 14 13 13 13 12 13 13 14 14  0  4
[133] 13 14 13 14 14 13 13 13 13 13  4  0

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
    S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12
S1   0  2  1  5  5  5 11 11  6  12  12  13
S2   2  0  1  5  6  6 12 12  6  12  14  14
S3   1  1  0  6  6  6 11 11  6  12  13  13
S4   5  5  6  0  1  9 12 12 10  12  13  14
S5   5  6  6  1  0  8 12 12 10  12  13  14
S6   5  6  6  9  8  0 12 12  6  12  12  13
S7  11 12 11 12 12 12  0  0  8   8  13  13
S8  11 12 11 12 12 12  0  0  8   8  13  13
S9   6  6  6 10 10  6  8  8  0  10  14  13
S10 12 12 12 12 12 12  8  8 10   0  14  13
S11 12 14 13 13 13 12 13 13 14  14   0   4
S12 13 14 13 14 14 13 13 13 13  13   4   0

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
    S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11
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

rm(comment)

## transform the dist object into a long dataframe of pairwise differences
dfl = dist2list(distobj)

## check class (i.e. [1] "data.frame")
class(dfl)

## check dimension (i.e. [1] 144   3)
dim(dfl)

## check the long dataframe
dfl
comment <- scan(what="character")
    col row value
1    S1  S1     0
2    S2  S1     2
3    S3  S1     1
4    S4  S1     5
5    S5  S1     5
6    S6  S1     5
7    S7  S1    11
8    S8  S1    11
9    S9  S1     6
10  S10  S1    12
11  S11  S1    12
12  S12  S1    13
13   S1  S2     2
14   S2  S2     0
15   S3  S2     1
16   S4  S2     5
17   S5  S2     6
18   S6  S2     6
19   S7  S2    12
20   S8  S2    12
21   S9  S2     6
22  S10  S2    12
23  S11  S2    14
24  S12  S2    14
25   S1  S3     1
26   S2  S3     1
27   S3  S3     0
28   S4  S3     6
29   S5  S3     6
30   S6  S3     6
31   S7  S3    11
32   S8  S3    11
33   S9  S3     6
34  S10  S3    12
35  S11  S3    13
36  S12  S3    13
37   S1  S4     5
38   S2  S4     5
39   S3  S4     6
40   S4  S4     0
41   S5  S4     1
42   S6  S4     9
43   S7  S4    12
44   S8  S4    12
45   S9  S4    10
46  S10  S4    12
47  S11  S4    13
48  S12  S4    14
49   S1  S5     5
50   S2  S5     6
51   S3  S5     6
52   S4  S5     1
53   S5  S5     0
54   S6  S5     8
55   S7  S5    12
56   S8  S5    12
57   S9  S5    10
58  S10  S5    12
59  S11  S5    13
60  S12  S5    14
61   S1  S6     5
62   S2  S6     6
63   S3  S6     6
64   S4  S6     9
65   S5  S6     8
66   S6  S6     0
67   S7  S6    12
68   S8  S6    12
69   S9  S6     6
70  S10  S6    12
71  S11  S6    12
72  S12  S6    13
73   S1  S7    11
74   S2  S7    12
75   S3  S7    11
76   S4  S7    12
77   S5  S7    12
78   S6  S7    12
79   S7  S7     0
80   S8  S7     0
81   S9  S7     8
82  S10  S7     8
83  S11  S7    13
84  S12  S7    13
85   S1  S8    11
86   S2  S8    12
87   S3  S8    11
88   S4  S8    12
89   S5  S8    12
90   S6  S8    12
91   S7  S8     0
92   S8  S8     0
93   S9  S8     8
94  S10  S8     8
95  S11  S8    13
96  S12  S8    13
97   S1  S9     6
98   S2  S9     6
99   S3  S9     6
100  S4  S9    10
101  S5  S9    10
102  S6  S9     6
103  S7  S9     8
104  S8  S9     8
105  S9  S9     0
106 S10  S9    10
107 S11  S9    14
108 S12  S9    13
109  S1 S10    12
110  S2 S10    12
111  S3 S10    12
112  S4 S10    12
113  S5 S10    12
114  S6 S10    12
115  S7 S10     8
116  S8 S10     8
117  S9 S10    10
118 S10 S10     0
119 S11 S10    14
120 S12 S10    13
121  S1 S11    12
122  S2 S11    14
123  S3 S11    13
124  S4 S11    13
125  S5 S11    13
126  S6 S11    12
127  S7 S11    13
128  S8 S11    13
129  S9 S11    14
130 S10 S11    14
131 S11 S11     0
132 S12 S11     4
133  S1 S12    13
134  S2 S12    14
135  S3 S12    13
136  S4 S12    14
137  S5 S12    14
138  S6 S12    13
139  S7 S12    13
140  S8 S12    13
141  S9 S12    13
142 S10 S12    13
143 S11 S12     4
144 S12 S12     0

rm(comment)

# perform a Receiver Operating Characteristic (ROC) analysis from the dataframe of pairwise differences

## suppose that the samples S1, S2, S3, S4, S5, S6 and S7 are positive controls (PC) of an outbreak
## suppose that samples S8, S9, S10, S11 and S12 are negative controls (NC) of an outbreak
## complete the file Controls.csv as below
comment <- scan(what="character")
sample	control
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

rm(comment)

## prepare dataframes of positive (PC) and negative (NC) controls
### read the dataframe of controls
dfc <- read.table("Controls.csv", dec = ".", header=TRUE, sep = ",", quote = "")
### subset PC
dfPC <- subset(dfc,dfc$control %in% c("PC"))
### subset NC
dfNC <- subset(dfc,dfc$control %in% c("NC"))

## derive variables "col" and "row" into variables "new col" and "new row" flagging positive (PC) and negative (NC) controls
### derivation of the "col" column
dfl$newcol <- ifelse(dfl$col %in% dfPC$sample, "PC",
                     ifelse(dfl$col %in% dfNC$sample, "NC",
                            "error"))
### derivation of the "row" column
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

## subset the related and unrelated pairs
### check dimension (i.e. [1] 144   6)
dim(dfl)
### check dataframe
dfl
comment <- scan(what="character")
col row value newcol newrow    status
1    S1  S1     0     PC     PC   related
2    S2  S1     2     PC     PC   related
3    S3  S1     1     PC     PC   related
4    S4  S1     5     PC     PC   related
5    S5  S1     5     PC     PC   related
6    S6  S1     5     PC     PC   related
7    S7  S1    11     PC     PC   related
8    S8  S1    11     NC     PC unrelated
9    S9  S1     6     NC     PC unrelated
10  S10  S1    12     NC     PC unrelated
11  S11  S1    12     NC     PC unrelated
12  S12  S1    13     NC     PC unrelated
13   S1  S2     2     PC     PC   related
14   S2  S2     0     PC     PC   related
15   S3  S2     1     PC     PC   related
16   S4  S2     5     PC     PC   related
17   S5  S2     6     PC     PC   related
18   S6  S2     6     PC     PC   related
19   S7  S2    12     PC     PC   related
20   S8  S2    12     NC     PC unrelated
21   S9  S2     6     NC     PC unrelated
22  S10  S2    12     NC     PC unrelated
23  S11  S2    14     NC     PC unrelated
24  S12  S2    14     NC     PC unrelated
25   S1  S3     1     PC     PC   related
26   S2  S3     1     PC     PC   related
27   S3  S3     0     PC     PC   related
28   S4  S3     6     PC     PC   related
29   S5  S3     6     PC     PC   related
30   S6  S3     6     PC     PC   related
31   S7  S3    11     PC     PC   related
32   S8  S3    11     NC     PC unrelated
33   S9  S3     6     NC     PC unrelated
34  S10  S3    12     NC     PC unrelated
35  S11  S3    13     NC     PC unrelated
36  S12  S3    13     NC     PC unrelated
37   S1  S4     5     PC     PC   related
38   S2  S4     5     PC     PC   related
39   S3  S4     6     PC     PC   related
40   S4  S4     0     PC     PC   related
41   S5  S4     1     PC     PC   related
42   S6  S4     9     PC     PC   related
43   S7  S4    12     PC     PC   related
44   S8  S4    12     NC     PC unrelated
45   S9  S4    10     NC     PC unrelated
46  S10  S4    12     NC     PC unrelated
47  S11  S4    13     NC     PC unrelated
48  S12  S4    14     NC     PC unrelated
49   S1  S5     5     PC     PC   related
50   S2  S5     6     PC     PC   related
51   S3  S5     6     PC     PC   related
52   S4  S5     1     PC     PC   related
53   S5  S5     0     PC     PC   related
54   S6  S5     8     PC     PC   related
55   S7  S5    12     PC     PC   related
56   S8  S5    12     NC     PC unrelated
57   S9  S5    10     NC     PC unrelated
58  S10  S5    12     NC     PC unrelated
59  S11  S5    13     NC     PC unrelated
60  S12  S5    14     NC     PC unrelated
61   S1  S6     5     PC     PC   related
62   S2  S6     6     PC     PC   related
63   S3  S6     6     PC     PC   related
64   S4  S6     9     PC     PC   related
65   S5  S6     8     PC     PC   related
66   S6  S6     0     PC     PC   related
67   S7  S6    12     PC     PC   related
68   S8  S6    12     NC     PC unrelated
69   S9  S6     6     NC     PC unrelated
70  S10  S6    12     NC     PC unrelated
71  S11  S6    12     NC     PC unrelated
72  S12  S6    13     NC     PC unrelated
73   S1  S7    11     PC     PC   related
74   S2  S7    12     PC     PC   related
75   S3  S7    11     PC     PC   related
76   S4  S7    12     PC     PC   related
77   S5  S7    12     PC     PC   related
78   S6  S7    12     PC     PC   related
79   S7  S7     0     PC     PC   related
80   S8  S7     0     NC     PC unrelated
81   S9  S7     8     NC     PC unrelated
82  S10  S7     8     NC     PC unrelated
83  S11  S7    13     NC     PC unrelated
84  S12  S7    13     NC     PC unrelated
85   S1  S8    11     PC     NC unrelated
86   S2  S8    12     PC     NC unrelated
87   S3  S8    11     PC     NC unrelated
88   S4  S8    12     PC     NC unrelated
89   S5  S8    12     PC     NC unrelated
90   S6  S8    12     PC     NC unrelated
91   S7  S8     0     PC     NC unrelated
92   S8  S8     0     NC     NC     extra
93   S9  S8     8     NC     NC     extra
94  S10  S8     8     NC     NC     extra
95  S11  S8    13     NC     NC     extra
96  S12  S8    13     NC     NC     extra
97   S1  S9     6     PC     NC unrelated
98   S2  S9     6     PC     NC unrelated
99   S3  S9     6     PC     NC unrelated
100  S4  S9    10     PC     NC unrelated
101  S5  S9    10     PC     NC unrelated
102  S6  S9     6     PC     NC unrelated
103  S7  S9     8     PC     NC unrelated
104  S8  S9     8     NC     NC     extra
105  S9  S9     0     NC     NC     extra
106 S10  S9    10     NC     NC     extra
107 S11  S9    14     NC     NC     extra
108 S12  S9    13     NC     NC     extra
109  S1 S10    12     PC     NC unrelated
110  S2 S10    12     PC     NC unrelated
111  S3 S10    12     PC     NC unrelated
112  S4 S10    12     PC     NC unrelated
113  S5 S10    12     PC     NC unrelated
114  S6 S10    12     PC     NC unrelated
115  S7 S10     8     PC     NC unrelated
116  S8 S10     8     NC     NC     extra
117  S9 S10    10     NC     NC     extra
118 S10 S10     0     NC     NC     extra
119 S11 S10    14     NC     NC     extra
120 S12 S10    13     NC     NC     extra
121  S1 S11    12     PC     NC unrelated
122  S2 S11    14     PC     NC unrelated
123  S3 S11    13     PC     NC unrelated
124  S4 S11    13     PC     NC unrelated
125  S5 S11    13     PC     NC unrelated
126  S6 S11    12     PC     NC unrelated
127  S7 S11    13     PC     NC unrelated
128  S8 S11    13     NC     NC     extra
129  S9 S11    14     NC     NC     extra
130 S10 S11    14     NC     NC     extra
131 S11 S11     0     NC     NC     extra
132 S12 S11     4     NC     NC     extra
133  S1 S12    13     PC     NC unrelated
134  S2 S12    14     PC     NC unrelated
135  S3 S12    13     PC     NC unrelated
136  S4 S12    14     PC     NC unrelated
137  S5 S12    14     PC     NC unrelated
138  S6 S12    13     PC     NC unrelated
139  S7 S12    13     PC     NC unrelated
140  S8 S12    13     NC     NC     extra
141  S9 S12    13     NC     NC     extra
142 S10 S12    13     NC     NC     extra
143 S11 S12     4     NC     NC     extra
144 S12 S12     0     NC     NC     extra

rm(comment)

### export the dataframe into a csv file
write.csv(dfl,file="PairwiseDataframe.csv")

### subset
dflsub <- subset(dfl,dfl$status %in% c("related","unrelated"))
## check dimension (i.e. [1] 95  6)
dim(dflsub)

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
### area under the curve of raw data (i.e. 85.34%)
auc(ROC)
### variance (i.e. 11.34112)
var <- var(ROC)
cat("Variance:" ,round(var,2),"\n")
### close MetricsROC.txt
sink()

## extract sensitivity and specificity of all thresholds of pairwise differences
thresholds <- coords(ROC, ret=c("threshold", "sensitivity", "specificity"))
thresholds
comment <- scan(what="character")
   threshold sensitivity specificity
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
best
comment <- scan(what="character")
threshold sensitivity specificity
1       9.5          80     75.5102

rm(comment)

## export the best threshold into a csv file
write.csv(best,file="ThresholdBest.csv")

# add a message
print("Developped by Nicolas Radomski on January 06 (2022) with the R version 4.1.2 (2021-11-01)")
