
### R/GWASRenl: An R package for multi-(environment/variate/treatment) GWAS analysis

**Authors**: Julian Taylor

This is the public facing GitHub repository version of the R package GWASReml.

**R/GWASReml** is a multi-(environment/variate/treatment) genome wide
association analysis R package that allows users to flexibly scan markers
within one or more chromosomes. The package allows various markers scanning
approaches including the well established approach outlined in Malosetti et
al. (2013). The main analysis function uses ASReml-R V4 for its core linear mixed modelling. To use full functionality of the package users will require a valid license for ASReml-R V4 and this can be obtained from [https://www.vsni.co.uk/software/asreml-r](https://www.vsni.co.uk/software/asreml-r). 

To install the package from GitHub you will need to do the following: 

1. Install the [devtools](https://cran.r-project.org/package=devtools) package. Do this by invoking R and then typing


```r
install.packages("devtools")
```

2. Install GWASreml using 


```r
devtools::install_github("DrJ001/GWASReml")
```

#### Getting Started

For a quick but complete introduction of the functionality of the package please
see the extensive help files of the package.

#### References

Malosetti,M., Ribaut,J.-M., and van Eeuwijk,F. A. (2013) The statistical analysis of multi-environment data: modeling genotype-by-environment interaction and its genetic basis. *Frontiers in physiology*, **4**, 44.

