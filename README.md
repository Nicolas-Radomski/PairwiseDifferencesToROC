# Usage
The R script PairwiseDifferencesToROC.R aims at deriving pairwise differences of microbial mutations (e.g. cg/wgMLST, genes, SNPs, InDels, kmers) into a long dataframe with the objective to perform receiving operating curve (ROC) analysis assessing the threshold of pairwise differences presenting the best combination of sensitivity and specificy to distinguish between samples related and unrelated to a studied outbreak.
# Input paiwise differences (i.e. PairwiseDifferences.csv)
```
    A  B  C  D  E  F  G  H  I  J
1   0  2  4  5  6  5 45 84 65 12
2   2  0  7  4  5 15 24 23  5 11
3   4  7  0  4 15 17 18 19 24 26
4   5  4  4  0 20 45 64 78 25 34
5   6  5 15 20  0 12 78 46 54  4
6   5 15 17 45 12  0 30 67  8 16
7  45 24 18 64 78 30  0 50 24 12
8  84 23 19 78 46 67 50  0 60 50
9  65  5 24 25 54  8 24 60  0 70
10 12 11 26 34  4 16 12 50 70  0
```
# Dependencies
The R script PairwiseDifferencesToROC.R was prepared and tested with R version 4.1.2 and RStudio 2021.09.1.
- library(data.table)
- library(usedist)
- library(spaa)
- library(pROC)
# Install R and Rstudio
## 1/ Install R (from configured sources)
```
sudo apt update
sudo apt install r-base
R --version
```
## 2/ Install RStudio (from dowloaded rstudio-1.3.1093-amd64.deb)
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
```
rstudio PairwiseDifferencesToROC
```
# Illustration
![ROC figure](https://github.com/Nicolas-Radomski/PairwiseDifferencesToROC/blob/main/illustration.png)
# Reference
- Bittinger K. Package ‘usedist’. 2020, cran.r-project.org, Version 0.4.0
- Jinlong Zhang J. Package 'spaa' - The Comprehensive R Archive Network. 2016, cran.r-project.org, Version 0.2.2
- Robin X., Turck N., Hainard A., Tiberti N., Lisacek F., Sanchez J.C., Müller M., Siegert S., Doering M. and Billings Z. pROC: Display and Analyze ROC Curves. 2021, cran.r-project.org, Version 1.18.0
# Acknowledgment
My former colleagues with whom I learned a lot about R
# Author
Nicolas Radomski
