# MGC

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [License](./LICENSE)
- [Issues](https://github.com/neurodata/mgc/issues)
- [Usage](#usage)
- [Pseudocode](#pseudocode)

# Overview

In modern scientific discovery, it is becoming increasingly critical to uncover whether one property of a dataset is related to another. The `MGC` (pronounced *magic*), or Multiscale Generalized Correlation, provides a framework for investigation into the relationships between properties of a dataset and the underlying geometries of the relationships, all while requiring sample sizes feasible in real data scenarios.

# Repo Contents

- [R](./R): `R` package code.
- [docs](./docs): package documentation.
- [man](./man): package manual for help in R session.
- [tests](./tests): `R` unit tests written using the `testthat` package.
- [vignettes](./vignettes): `R` vignettes for R session html help pages.

This package is supported for *Linux*, *Mac OSX*, and *Windows* operating systems.

# System Requirements

Before setting up the MGC package, users should have R version 3.4.0 or higher, and several packages set up from CRAN.

## Installing R version 3.4.2 on Ubuntu

the latest version of R can be installed by adding the latest repository to `apt`:

```
sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -
sudo apt-get update
sudo apt-get install r-base r-base-dev
```

## Package dependencies

### Required

From inside an `R` session, we can install the `MGC` dependencies with the following command:

```
install.packages(c('ggplot2', 'reshape2', 'Rmisc', 'devtools'))  # set up package dependencies
```

### Optional

If the user wishes to use the interactive demos (`R` vignettes) or run the tests, then the following packages must also be set up:

```
install.packages(c('testthat', 'knitr', 'rmarkdown'))
```

# Installation Guide

## No Vignettes

From the same `R` session, type:

```
require(devtools)
install_github('neurodata/mgc')  # install MGC
require(MGC)  # source the package now that it is set up
```

The package should take approximately 4 seconds on a standard computer (4 cores @ 3.3 GHz/core, 16 GB RAM with internet at least 25 Mbps). 

## Vignettes

If the user wishes to use the interactive demos (`R` vignettes) or run the tests, from the same `R` session, type:

```
install_github('neurodata/mgc', build_vignettes=TRUE)  # install MGC with the vignettes
require(MGC)  # source the package now that it is set up
vignette("MGC Statistic", package="MGC")  # view one of the basic vignettes
```

The package should take approximately 8 seconds with vignettes on a standard computer (4 cores @ 3.3 GHz/core, 16 GB RAM with internet at least 25 Mbps). 

# Demo

## Simulated Dataset

The `MGC` demo can be run as follows:

```
require(MGC)
set.seed(12345)
mgc.sample(MGC::test_xunif, MGC::test_ylin)$statMGC  # test with a linear relationship between x and y
```

the x data provided is by sampling 100 times from a uniform distribution on the interval [-1, 1], and the y data is formed by adding normally distributed error with variance 0.2 (indicating a linear relationship).

and is expected  to produce the following result exactly approximately *instantaneously*:

```
0.891153
```

a more interactive demo is provided in the package vignette (if installed):

```
vignette("MGC Statistic", package="MGC")
```


## Real Dataset

In the below demo, we show the result of `MGC` to determine the relationship between the first (sepal length) and third (petal length) dimensions of the `iris` dataset:

```
set.seed(12345)
res <- mgc.sample(iris[,1], iris[,2])
mgc.plot.plot_matrix(res$localCorr, title="MGC Corr Map, Sepal Length and Petal Length",
    xlab="Sepal Length Neighbors", ylab="Petal Length Neighbors", legend.name = "statMGC")
```

![image](https://user-images.githubusercontent.com/8883547/32354581-bb214476-c000-11e7-995e-716e00fd8863.png)

viewing the corr map above we see that the relationship betweel Sepal and Petal Length is somewhat linear.

# Usage

Instructions for usage of the `MGC` statistic can be invoked from the `R` terminal window:

```
help(MGC)
```

which produces the `man` page:

```
The main function that computes the MGC measure between two datasets: It first computes all local correlations, then use the maximal statistic among all local correlations based on thresholding.

Description

The main function that computes the MGC measure between two datasets: It first computes all local correlations, then use the maximal statistic among all local correlations based on thresholding.

Usage

mgc.sample(A, B, option = "mgc")
Arguments

+ A	is interpreted as:
a [nxn] distance matrix, A is a square matrix with zeros on diagonal
a [nxd] data matrix, Otherwise
+ B	is interpreted as:
a [nxn] distance matrix, B is a square matrix with zeros on diagonal
a [nxd] data matrix, Otherwise
+ option='mgc'	is a string that specifies which global correlation to build up-on.
    + 'mgc' use the MGC global correlation.
    + 'dcor' use the dcor global correlation.
    + 'mantel' use the mantel global correlation.
    + 'rank' use the rank global correlation.

and Returns:
+ statMGC is the sample MGC statistic within [-1,1];

+ localCorr consists of all local correlations by double matrix index;

+ optimalScale the estimated optimal scale in matrix single index.

Author(s)

C. Shen
```

The discriminability statistic can be invoked as follows:

```
help(discr.discr)
```

and shows:

```
Description

A function for computing the discriminability from a distance matrix and a set of associated labels.

Usage

discr.discr(X, ids, thresh = 0, verbose = FALSE)
Arguments

+ X	is interpreted as:
    + a [n x n] distance matrix if X is a square matrix with zeros on diagonal
    + a [n x d] data matrix, Otherwise
+ ids	[n] is a vector containing the labels for our n observations.
+ remove_outliers=TRUE	boolean indicating whether to ignore observations with rdf below a certain cutoff.
+ thresh=0	the threshold below which to ignore observations. If thresh > 0, ignores observations where the rdf is < thresh in the discriminability computation.
+ verbose=FALSE	a boolean indicating whether to:
    + TRUE print output to console
    + FALSE Do not print output to console

Returns
+ discr the discriminability statistic.

Author(s)

Eric Bridgeford and Gregory Kiar
```

# Pseudocode

Pseudocode for the meethods employed in the `MGC` package can be found on the [arXiv - MGC](https://arxiv.org/abs/1609.05148) in Appendix C (starting on page 30). 
