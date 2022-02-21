# Usage
The R scripts PairwiseDifferences2ROC.R (detailed algorithm with Rstudio) and PairwiseDifferencesToROC.R (automatic algorithm with Rscript) aim at deriving profiles of microbial mutations (e.g. cg/wgMLST, genes, SNPs, InDels, kmers) (Profiles.csv) into a matrix (PairwiseMatrix.csv) and a dataframe (PairwiseDataframe.csv) of pairwise differences with the objective to perform:
- a Receiver Operating Characteristic (ROC) analysis (ROC.pdf and MetricsROC.txt) assessing the threshold of pairwise differences (Thresholds.csv) presenting the best combination of sensitivity and specificity (ThresholdBest.csv) to distinguish between samples related (i.e. positive controls PC in Types.csv) and unrelated (i.e. negative controls NC in Types.csv) to a studied outbreak
- proportions of pairwise differences (between tested and positive samples, i.e. tested samples TS in Types.csv) lower or higher than the best ROC threshold (i.e. Proportions.csv).
# Input
## Profiles of microbial mutations (i.e. Profiles.csv)
- S stands for sample
- L stands for locus
- A stands for allele
- 0 stands dor missing data
```
sample  L1  L2  L3  L4  L5  L6  L7  L8   L9  L10 L11 L12 L13 L14 L15 L16 L17 L18  L19 L20
    S1 G20 G15 G55 G12 G30 G11 G24 G66  G12  G55 G66  G5 G86 G54 G47   0 G63 G76 G100 G32
    S2 G20 G15 G55 G12 G30 G11 G24 G66  G12  G55 G66  G2 G87 G54 G47 G23   0 G76 G100 G32
    S3 G20 G15 G55 G12 G30 G11 G24 G66  G12  G55 G66  G2 G86 G54 G47 G23 G63   0 G100 G32
    S4 G20 G31 G55 G30 G30 G11 G55 G66  G55  G55 G66  G5 G87 G54 G47 G23 G63 G76    0 G32
    S5 G20 G31 G55 G30 G30 G11 G55 G66  G55  G55 G66  G5 G98 G54 G47 G23 G63 G76 G100   0
    S6 G10 G15 G10 G12 G30 G10 G24 G10  G12  G55 G66  G5 G98 G54 G47 G23 G63 G76 G100 G32
    S7 G41 G22 G41 G22 G22 G41 G27 G41  G27  G27 G66  G9 G86 G54 G47 G23 G63 G76 G100 G32
    S8 G41 G22 G41 G22 G22 G41 G27 G41  G27  G27 G66  G9 G86 G54 G47 G23 G63 G76 G100 G32
    S9 G41 G15 G41 G12 G30 G41 G24 G41  G12  G55 G66  G8 G97 G54 G47 G23 G63 G76 G100 G32
   S10 G50 G22 G50 G22 G55 G51 G27 G50  G27  G66 G66  G8 G97 G54 G47 G23 G63 G76 G100 G32
   S11 G10 G54 G15 G41 G65 G88 G75 G89 G420 G998 G66  G5 G86 G11 G10 G23 G63 G76 G100 G32
   S12 G10 G54 G98 G41 G65 G88 G75 G89 G420 G998 G66  G8 G86 G14  G1 G23 G63 G76 G100 G32
   S13 G20 G15 G55 G12 G30 G11 G24 G66  G12  G55 G66  G2 G87 G54 G47 G23 G63 G76 G100 G32
   S14 G41 G15 G41 G12 G30 G41 G24 G41  G12  G55 G66  G8 G97 G54 G47 G23 G63 G76 G100 G32
   S15 G20 G15 G55 G12 G30 G11 G24 G66  G12  G55 G66  G5 G86 G54 G47 G23 G63 G76 G100 G32
   S16 G10 G54 G15 G41 G65 G88 G75 G89 G420 G998 G66  G5 G86 G11 G10 G23 G63 G76 G100 G32
   S17 G20 G31 G55 G30 G30 G11 G55 G66  G55  G55 G66  G5 G87 G54 G47 G23 G63 G76 G100 G32
   S18 G41 G22 G41 G22 G22 G41 G27 G41  G27  G27 G66  G9 G86 G54 G47 G23 G63 G76 G100 G32
   S19 G41 G22 G41 G22 G22 G41 G27 G41  G27  G27 G66  G9 G86 G54 G47 G23 G63 G76 G100 G32
   S20 G10 G54 G98 G41 G65 G88 G75 G89 G420 G998 G66  G8 G86 G14  G1 G23 G63 G76 G100 G32
```
## Positive (PC) and negative (NC) controls of the outbreak, as well as tested samples (TS) (i.e. Types.csv)
```
Sample	Type
    S1    PC
    S2    PC
    S3    PC
    S4    PC
    S5    PC
    S6    PC
    S7    PC
    S8    NC
    S9    NC
   S10    NC
   S11    NC
   S12    NC
   S13    TS
   S14    TS
   S15    TS
   S16    TS
   S17    TS
   S18    TS
   S19    TS
   S20    TS
```

# Output
## Matrix of pairwise differences of microbial mutations (i.e. PairwiseMatrix.csv)
```
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
```
## Dataframe of pairwise differences of microbial mutations (i.e. PairwiseDataframe.csv)
```
FirstSample SecondSample Differences FirstType SecondType    Status
         S1           S2           2        PC         PC   related
         S1           S3           1        PC         PC   related
         S1           S4           5        PC         PC   related
         S1           S5           5        PC         PC   related
         S1           S6           5        PC         PC   related
         S1           S7          11        PC         PC   related
         S1           S8          11        PC         NC unrelated
         S1           S9           6        PC         NC unrelated
         S1          S10          12        PC         NC unrelated
         S1          S11          12        PC         NC unrelated
         S1          S12          13        PC         NC unrelated
         S1          S13           2        PC         TS      test
         S1          S14           6        PC         TS      test
         S1          S15           0        PC         TS      test
         S1          S16          12        PC         TS      test
         S1          S17           5        PC         TS      test
         S1          S18          11        PC         TS      test
         S1          S19          11        PC         TS      test
         S1          S20          13        PC         TS      test
...
```
## ROC analysis and related metrics (i.e. ROC.pdf and MetricsROC.txt)
```
Working directory path: /home/IZSNT/n.radomski/Documents/RstudioWorkingDirectory/PairwiseDifferencesToROC-20220217 
First argument: Profiles.csv 
Second argument: Types.csv 
Area under the curve: 83.13%
Variance: 28.63 
```
## Sentitivity and specificity of thresholds (i.e. Thresholds.csv)
```
Threshold Sensitivity Specificity
      Inf   100.00000     0.00000
     13.5   100.00000    11.42857
     12.5   100.00000    34.28571
     11.5    80.95238    68.57143
     10.5    71.42857    74.28571
      9.5    71.42857    80.00000
      8.5    66.66667    80.00000
      7.0    61.90476    85.71429
      5.5    38.09524    97.14286
      3.5    19.04762    97.14286
      1.5    14.28571    97.14286
      0.5     0.00000    97.14286
     -Inf     0.00000   100.00000
```
## Sentitivity and specificity of the best threshold (i.e. ThresholdBest.csv)
```
Threshold Sensitivity Specificity
      9.5    71.42857          80
```

## Proportions of pairwise differences (between tested and positive samples) lower or higher than the threshold (i.e. Proportions.csv)
```
TestedSample LowerThreshold HigherThreshold ProportionLower ProportionHigher
         S13              5               1            83.3             16.7
         S14              4               2            66.7             33.3
         S15              5               1            83.3             16.7
         S16              0               6             0.0            100.0
         S17              5               1            83.3             16.7
         S18              1               5            16.7             83.3
         S19              1               5            16.7             83.3
         S20              0               5             0.0            100.0
```

# Install R and Rstudio
## 1/ Install R (from configured sources)
```
sudo apt update
sudo apt install r-base
R --version
```
## 2/ Install RStudio (from dowloaded rstudio-2021.09.1-372-amd64.deb)
```
sudo apt install gdebi-core
sudo gdebi /home/Downloads/rstudio-2021.09.1-372-amd64.deb
rstudio --version
```
# Update R version from 3.x to 4.x
## 1/ Check the current R version
```
R --version
```
## 2/ Add key and application repository
```
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
```
## 3/ Install the R version 4.x
```
sudo apt install r-base
```
## 4/ Check the current R version
```
R --version
```
# Update R 3.x packages to R 4.x
## 1/ Update installed packages from the R console
```
update.packages(ask = FALSE, checkBuilt = TRUE)
```
## 2/ Update missing packages from the R console
```
old_packages <- installed.packages(lib.loc = "/home/IZSNT/n.radomski/R/x86_64-pc-linux-gnu-library/3.6/")
new_packages <- installed.packages()
missing_packages <- as.data.frame(old_packages[!old_packages[, "Package"] %in% new_packages[, "Package"],])
install.packages(missing_packages$Package)
```
# Install R dependencies and launch with R
## Install dependencies from the R console
The R scripts PairwiseDifferences2ROC.R (detailed algorithm with Rstudio) and PairwiseDifferencesToROC.R (automatic algorithm with Rscript) were prepared and tested with R version 4.1.2 and RStudio 2021.09.1.
```
install.packages("benchmarkme")
install.packages("data.table")
install.packages("spaa")
install.packages("pROC")
install.packages("plyr")
```
## Launch each command from Rstudio (i.e. PairwiseDifferences2ROC.R detailed algorithm)
```
git clone https://github.com/Nicolas-Radomski/PairwiseDifferencesToROC.git
cd PairwiseDifferencesToROC
rstudio PairwiseDifferences2ROC.R
```
## Launch the whole script from Rscript (i.e. PairwiseDifferencesToROC.R automatic algorithm)
```
git clone https://github.com/Nicolas-Radomski/PairwiseDifferencesToROC.git
cd PairwiseDifferencesToROC
Rscript PairwiseDifferencesToROC.R Profiles.csv Types.csv
```
# Install Docker image and launch with Docker
## 1/ Pull Docker image from Docker Hub
```
docker pull nicolasradomski/pairwisedifferencestoroc
```
## 2/ Launch with Docker and different paired-trees
```
docker run --name nicolas --rm -v /home/data:/data -v /home/output:/output nicolasradomski/pairwisedifferencestoroc:latest sh -c 'Rscript code/PairwiseDifferencesToROC.R data/Profiles.csv data/Types.csv' > output/std.log 2>&1

```
# Illustration
![ROC figure](https://github.com/Nicolas-Radomski/PairwiseDifferencesToROC/blob/main/illustration.png)
# References
- Jinlong Zhang J. Package 'spaa' - The Comprehensive R Archive Network. 2016, cran.r-project.org, Version 0.2.2
- Robin X., Turck N., Hainard A., Tiberti N., Lisacek F., Sanchez J.C., Müller M., Siegert S., Doering M. and Billings Z. pROC: Display and Analyze ROC Curves. 2021, cran.r-project.org, Version 1.18.0
- Robin, X., Turck, N., Hainard, A. et al. pROC: an open-source package for R and S+ to analyze and compare ROC curves. 2011, BMC Bioinformatics, 12(77): 1-8
- Docker Hub: https://hub.docker.com/r/nicolasradomski/pairwisedifferencestoroc
# Acknowledgment
Adriano Di Pasquale for our discussions about algorithmic approaches
# Author
Nicolas Radomski
