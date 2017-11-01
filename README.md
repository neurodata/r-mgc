# MGC

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [License](./LICENSE)
- [Issues](https://github.com/neurodata/mgc/issues)

# Overview

In modern scientific discovery, it is becoming increasingly critical to uncover whether one property of a dataset is related to another. The `MGC` (pronounced *magic*), or Multiscale Generalized Correlation, provides a framework for investigation into the relationships between properties of a dataset and the underlying geometries of the relationships, all while requiring sample sizes feasible in real data scenarios.

# Repo Contents

- [./R]: `R` package code.
- [./docs]: package documentation.
- [./man]: package manual for help in R session.
- [./tests]: `R` unit tests written using the `testthat` package.
- [./vignettes]: `R` vignettes for R session html help pages.

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

# Demo

The `MGC` demo can be run as follows:

```
mgc.run_demo()  # run the demo
```

and is expected  to produce the following result (results will not be exact as the data is randomly generated):

```
Under Linear Dependency of 20 observations, MGC Statistic, p-value, and optimal scales (matrix single index) are:  1 , 0 , 156 173 174 175 176 191 192 193 194 195 196 210 211 212 213 214 215 216 230 231 233 234 235 236 249 250 251 252 256 269 270 271 272 274 275 276 289 290 291 292 294 296 308 309 310 311 312 313 314 315 
Under Quadratic Dependency of 20 observations, MGC Statistic, p-value, and optimal scales (matrix single index) are:  0.56 , 0 , 276 
```

The `discriminability` demo can be run as follows:

```
discr.run_demo()
```

and is expected to produce the following result:

```
[1] "Parameters of first random variable: mu=0.00, sd=0.20"
[1] "Parameters of second random variable: mu=2.00, sd=0.20"
[1] "Discriminability statistic for 2 observations per random variable: 1.00"
```

