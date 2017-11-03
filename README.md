# MGC

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo and Usage](#demo-usage)
    + [MGC](#mgc-usage)
    + [Discriminability](#discriminability-usage)
- [License](./LICENSE)
- [Issues](https://github.com/neurodata/mgc/issues)
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

# Demo and Usage

## MGC Usage

### Simulated Dataset

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

### Real Dataset

In the below demo, we show the result of `MGC` to determine the relationship between the first (sepal length) and third (petal length) dimensions of the `iris` dataset, which should run in about 2 seconds:

```
set.seed(12345)
res <- mgc.sample(iris[,1], iris[,3])
mgc.plot.plot_matrix(res$localCorr, title="MGC Corr Map, Sepal Length and Petal Length",
    xlab="Sepal Length Neighbors", ylab="Petal Length Neighbors", legend.name = "statMGC")
print(res$statMGC)
```

![image](https://user-images.githubusercontent.com/8883547/32355967-7de64590-c008-11e7-9c3b-e24470fdbdaa.png)

viewing the corr map above we see that the relationship betweel Sepal and Petal Length is somewhat linear.


a more interactive demo is provided in the package vignette (if installed):

```
vignette("MGC Statistic", package="MGC")
```

### Usage

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

## Discriminability Usage

### Simulated Dataset

Here, we assume that we have 5 independent sources of a measurement, and take 10 measurements from each source. Each measurement source *i = 1:n* is a random variable where measurements taken at the source for the *j*'th dimension is normally distributed with mean *i*.

```
require(MGC)
require(latex2exp)

nsrc <- 5
nobs <- 10
d <- 20
src_id <- array(1:nsrc)
set.seed(12345)
labs <- sample(rep(src_id, nobs))
dat <- t(sapply(labs, function(lab) rnorm(d, mean=lab, sd=0.5)))
discr.discr(dat, labs)  # expect high discriminability since measurements taken at a source have the same mean and sd of only 0.5
```

which should show approximately *instantaneously*:

```
0.9983889
```

### Real Dataset

Below, we show how discriminability might be used on real data, by demonstrating its usage on the first $4$ dimensions of the `iris` dataset, to determine the relationship between the flower species and the distances between the different dimensions of the iris dataset (sepal width/length and petal width/length):

```{r, fig.width=6, fig.height=4}
Dx <- as.matrix(dist(iris[sort(as.vector(iris$Species), index=TRUE)$ix,]))

mgc.plot.plot_matrix(Dx, xlab="Sorted by Species", ylab="Sorted by Species", title="Distance Matrix for Iris Dataset", legend.name="dist(x, y)"))
```

which produces the following distance matrix. The block structure indicates that the flowers from a particular species have more similar properties (as defined by the euclidian distance) to the flowers in the same species than to flowers of different species:

![image](https://user-images.githubusercontent.com/8883547/32355935-466bf308-c008-11e7-980a-64ac3a1b22da.png)

this is confirmed by the high discriminability statistic which should occur approximately *instantaneously*:

```
discr.discr(iris[,c(1,2,3,4)], as.vector(iris$Species))
```

which should show:

```
0.9320476
```

### Usage

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
