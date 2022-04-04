# R-package irrNA: Coefficients of Interrater Reliability - Generalized for Randomly Incomplete Datasets
## Indication
irrNA provides coefficients of interrater reliability that are generalized to cope with randomly incomplete (i.e. unbalanced) datasets without any imputation of missing values or any (row-wise or column-wise) omissions of actually available data. Applied to complete (balanced) datasets, these generalizations yield the same results as the common procedures, namely the Intraclass Correlation according to McGraw & Wong (1996) and the Coefficient of Concordance according to Kendall & Babington Smith (1939).

## Usage
To get startet you could type the following line-commands into the R-console:
```
install.packages(irrNA)
library(irrNA)
example(iccNA)
example(icc_corr)
example(kendallNA)
```
For metrically scaled data please use iccNA; kendallNA is for ordinally scaled data.
For further information please refer to https://cran.r-project.org/web/packages/irrNA/irrNA.pdf

## Repository
You can find the latest version of irrNA on all CRAN-mirrors or at https://CRAN.R-project.org/package=irrNA .