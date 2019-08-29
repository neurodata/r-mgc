This R repo is a development branch, the actively developed repo is in Python at https://github.com/neurodata/mgcpy.

# Multiscale Graph Correlation (MGC)

[![CRAN Status Badge](http://www.r-pkg.org/badges/version/mgc)](http://cran.r-project.org/web/packages/mgc)
[![eLife Shield](https://img.shields.io/badge/eLife-41690-green.svg?style=flat)](https://elifesciences.org/articles/41690)
[![Downloads badge](https://cranlogs.r-pkg.org/badges/mgc)](https://cranlogs.r-pkg.org/badges/mgc)
[![Zenodo Shield](https://zenodo.org/badge/103285533.svg)](https://zenodo.org/badge/latestdoi/103285533)
[![Travis](https://travis-ci.org/neurodata/r-mgc.svg?branch=master)](https://travis-ci.org/neurodata/r-mgc)

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Instructions for Use](#instructions-for-use)
- [License](./LICENSE)
- [Issues](https://github.com/neurodata/mgc/issues)
- [Citation](#citation)
- [Reproducibility](#reproducibility)

# Overview

In modern scientific discovery, it is becoming increasingly critical to uncover whether one property of a dataset is related to another. The `MGC` (pronounced *magic*), or Multiscale Graph Correlation, provides a framework for investigation into the relationships between properties of a dataset and the underlying geometries of the relationships, all while requiring sample sizes feasible in real data scenarios.

# Repo Contents

- [R](./R): `R` package code.
- [docs](./docs): package documentation.
- [man](./man): package manual for help in R session.
- [tests](./tests): `R` unit tests written using the `testthat` package.
- [vignettes](./vignettes): `R` vignettes for R session html help pages.


# System Requirements

## Hardware Requirements

The `MGC` package requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB  
CPU: 4+ cores, 3.3+ GHz/core

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, 4 cores@3.3 GHz) and internet of speed 25 Mbps.

## Software Requirements

### OS Requirements

This package is supported for *Linux* operating systems. The package has been tested on the following systems:

Linux: Ubuntu 16.04  
Mac OSX:  
Windows:  

Before setting up the `MGC` package, users should have `R` version 3.4.0 or higher, and several packages set up from CRAN.

#### Installing R version 3.4.2 on Ubuntu 16.04

the latest version of R can be installed by adding the latest repository to `apt`:

```
sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -
sudo apt-get update
sudo apt-get install r-base r-base-dev
```

which should install in about 20 seconds.

#### Package dependencies

Users should install the following packages prior to installing `mgc`, from an `R` terminal:

```
install.packages(c('ggplot2', 'reshape2', 'Rmisc', 'devtools', 'testthat', 'knitr', 'rmarkdown', 'latex2exp', 'MASS'))
```

which will install in about 80 seconds on a recommended machine.

#### Package Versions

The `mgc` package functions with all packages in their latest versions as they appear on `CRAN` on October 15, 2017. Users can check [CRAN snapshot](https://mran.microsoft.com/timemachine/) for details. The versions of software are, specifically:
```
ggplot2: 2.2.1
reshape2: 1.4.2
Rmisc: 1.5
devtools: 1.13.3
testthat: 0.2.0
knitr: 1.17
rmarkdown: 1.6
latex2exp: 0.4.0
MASS: 7.3
```

If you are having an issue that you believe to be tied to software versioning issues, please drop us an [Issue](https://github.com/neurodata/mgc/issues). 

# Installation Guide

From an `R` session, type:

```
require(devtools)
install_github('neurodata/mgc', build_vignettes=TRUE)  # install mgc with the vignettes
require(mgc)  # source the package now that it is set up
vignette("MGC", package="mgc")  # view one of the basic vignettes
```

The package should take approximately 20 seconds to install with vignettes on a recommended computer. 


### Instructions for Use

Please see the vignettes for help using the package:

```
vignette("MGC", package="mgc")
vignette("Discriminability", package="mgc")
vignette("simulations", package="mgc")
```


# Pseudocode

Pseudocode for the methods employed in the `mgc` package can be found on the [arXiv - MGC](https://arxiv.org/abs/1609.05148) in Appendix C (starting on page 30).

# Citation

For citing code or the paper, please use the citations found in [citation.bib](./citation.bib).

# Reproducibility

All the code to reproduce any figures from https://arxiv.org/abs/1609.05148 is available here https://github.com/neurodata/mgc-paper
