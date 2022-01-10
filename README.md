# Usage
The R scripts PairwiseDifferences2ROC.R (detailed algorithm with Rstudio) and PairwiseDifferencesToROC.R (automatic algorithm with Rscript) aim at deriving profiles of microbial mutations (e.g. cg/wgMLST, genes, SNPs, InDels, kmers) (Profiles.csv) into a matrix (PairwiseMatrix.csv) and a dataframe (PairwiseDataframe.csv) of pairwise differences with the objective to perform a Receiver Operating Characteristic (ROC) analysis (ROC.pdf and MetricsROC.txt) assessing the threshold of pairwise differences (Thresholds.csv) presenting the best combination of sensitivity and specificity (ThresholdBest.csv) to distinguish between samples related (i.e. positive controls PC in Controls.csv) and unrelated (i.e. negative controls NC in Controls.csv) to a studied outbreak.
# Input
## Profiles of microbial mutations (i.e. Profiles.csv)
- S stands for sample
- L stands for locus
- A stands for allele
```
sample  L1  L2  L3  L4  L5  L6  L7  L8   L9  L10 L11 L12 L13 L14 L15
    S1 A20 A15 A55 A12 A30 A11 A24 A66  A12  A55 A66  A5 A86 A54 A47
    S2 A20 A15 A55 A12 A30 A11 A24 A66  A12  A55 A66  A2 A87 A54 A47
    S3 A20 A15 A55 A12 A30 A11 A24 A66  A12  A55 A66  A2 A86 A54 A47
    S4 A20 A31 A55 A30 A30 A11 A55 A66  A55  A55 A66  A5 A87 A54 A47
    S5 A20 A31 A55 A30 A30 A11 A55 A66  A55  A55 A66  A5 A98 A54 A47
    S6 A10 A15 A10 A12 A30 A10 A24 A10  A12  A55 A66  A5 A98 A54 A47
    S7 A41 A22 A41 A22 A22 A41 A27 A41  A27  A27 A66  A9 A86 A54 A47
    S8 A41 A22 A41 A22 A22 A41 A27 A41  A27  A27 A66  A9 A86 A54 A47
    S9 A41 A15 A41 A12 A30 A41 A24 A41  A12  A55 A66  A8 A97 A54 A47
   S10 A50 A22 A50 A22 A55 A51 A27 A50  A27  A66 A66  A8 A97 A54 A47
   S11 A10 A54 A15 A41 A65 A88 A75 A89 A420 A998 A66  A5 A86 A11 A10
   S12 A10 A54 A98 A41 A65 A88 A75 A89 A420 A998 A66  A8 A86 A14  A1
```
## Positive (PC) and negative (NC) controls of the outbreak (i.e. Controls.csv)
```
sample	control
    S1	     PC
    S2	     PC
    S3	     PC
    S4	     PC
    S5	     PC
    S6	     PC
    S7	     PC
    S8	     NC
    S9	     NC
   S10	     NC
   S11	     NC
   S12	     NC
```

# Output
## Matrix of pairwise differences of microbial mutations (i.e. PairwiseMatrix.csv)
```
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
```
## Dataframe of pairwise differences of microbial mutations (i.e. PairwiseDataframe.csv)
```
col row value newcol newrow    status
 S1  S1     0     PC     PC   related
 S2  S1     2     PC     PC   related
 S3  S1     1     PC     PC   related
 S4  S1     5     PC     PC   related
 S5  S1     5     PC     PC   related
 S6  S1     5     PC     PC   related
 S7  S1    11     PC     PC   related
 S8  S1    11     NC     PC unrelated
 S9  S1     6     NC     PC unrelated
S10  S1    12     NC     PC unrelated
S11  S1    12     NC     PC unrelated
S12  S1    13     NC     PC unrelated
...
```
## ROC analysis and related metrics (i.e. ROC.pdf and MetricsROC.txt)
```
Working directory path: /home/IZSNT/n.radomski/Downloads/PairwiseDifferencesToROC-main 
First argument: Profiles.csv 
Second argument: Controls.csv 
Area under the curve: 85.34%
Variance: 11.34
```
## Sentitivity and specificity of thresholds (i.e. Thresholds.csv)
```
threshold sensitivity specificity
     -Inf   100.00000     0.00000
      0.5    97.14286    14.28571
      1.5    97.14286    26.53061
      3.5    97.14286    30.61224
      5.5    97.14286    46.93878
      7.0    85.71429    67.34694
      8.5    80.00000    71.42857
      9.5    80.00000    75.51020
     10.5    74.28571    75.51020
     11.5    68.57143    83.67347
     12.5    34.28571   100.00000
     13.5    11.42857   100.00000
      Inf     0.00000   100.00000
```
## Sentitivity and specificity of the best threshold (i.e. ThresholdBest.csv)
```
threshold sensitivity specificity
      9.5          80     75.5102
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
# Start
## Install dependencies from the R console
The R scripts PairwiseDifferences2ROC.R (detailed algorithm with Rstudio) and PairwiseDifferencesToROC.R (automatic algorithm with Rscript) were prepared and tested with R version 4.1.2 and RStudio 2021.09.1.
```
install.packages("ape")
install.packages("data.table")
install.packages("spaa")
install.packages("pROC")
```
## Launch each command from Rstudio (i.e. detailed algorithm)
```
git clone https://github.com/Nicolas-Radomski/PairwiseDifferencesToROC.git
cd PairwiseDifferencesToROC
rstudio PairwiseDifferences2ROC.R
```
## Launch the whole script from Rscript (i.e. automatic algorithm)
```
git clone https://github.com/Nicolas-Radomski/PairwiseDifferencesToROC.git
cd PairwiseDifferencesToROC
Rscript PairwiseDifferencesToROC.R Profiles.csv Controls.csv
```
# Illustration
![ROC figure](https://github.com/Nicolas-Radomski/PairwiseDifferencesToROC/blob/main/illustration.png)
# Reference
- Jinlong Zhang J. Package 'spaa' - The Comprehensive R Archive Network. 2016, cran.r-project.org, Version 0.2.2
- Robin X., Turck N., Hainard A., Tiberti N., Lisacek F., Sanchez J.C., Müller M., Siegert S., Doering M. and Billings Z. pROC: Display and Analyze ROC Curves. 2021, cran.r-project.org, Version 1.18.0
- Robin, X., Turck, N., Hainard, A. et al. pROC: an open-source package for R and S+ to analyze and compare ROC curves. 2011, BMC Bioinformatics, 12(77): 1-8
# Acknowledgment
Adriano Di Pasquale for our discussions about algorithmic approaches
# Author
Nicolas Radomski
