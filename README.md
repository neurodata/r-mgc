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

From inside an `R` session, we can install the `MGC` dependencies with the following command:

```
install.packages(c('ggplot2', 'reshape2', 'Rmisc'))  # set up package dependencies
```

# Installation Guide

From inside an `R` session, we must first install and source `devtools`:

```
install.packages('devtools')  # set up devtools
require(devtools)  # source it
install_github('neurodata/mgc')  # install MGC
require(MGC)  # source the package now that it is set up
```

The package should take approximately 8 seconds on a standard computer (4 cores @ 3.3 GHz/core, 16 GB RAM with internet at least 25 Mbps). 

# Demo

The `MGC` demo can be run as follows:

```
mgc.run_demo()  # run the demo
```

and is expected  to produce the following result (results will not be exact as the data is randomly generated) in about 10 seconds:

```
Under Linear Dependency of 20 observations, MGC Statistic, p-value, and optimal scales (matrix single index) are:  1 , 0 , 156 173 174 175 176 191 192 193 194 195 196 210 211 212 213 214 215 216 230 231 233 234 235 236 249 250 251 252 256 269 270 271 272 274 275 276 289 290 291 292 294 296 308 309 310 311 312 313 314 315 
Under Quadratic Dependency of 20 observations, MGC Statistic, p-value, and optimal scales (matrix single index) are:  0.56 , 0 , 276 
```

The `discriminability` demo can be run as follows:

```
discr.run_demo()
```


and is expected to produce the following result, again in about 10 seconds:

```
[1] "Parameters of first random variable: mu=0.00, sd=0.20"
[1] "Parameters of second random variable: mu=2.00, sd=0.20"
[1] "Discriminability statistic for 2 observations per random variable: 1.00"
```



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
