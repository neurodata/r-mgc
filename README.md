This R repo is a development branch, the actively developed repo is in Python at https://github.com/neurodata/hyppo.

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

Linux: Ubuntu 20.04, 18.04 
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
install_github('neurodata/r-mgc', build_vignettes=TRUE)  # install mgc with the vignettes
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

### MGC

All the code to reproduce any figures from https://arxiv.org/abs/1609.05148 is available here https://github.com/neurodata/mgc-paper.

### Discriminability

Here, we describe how to reproduce the manuscript figures from the discriminability paper. To begin, clone this repository locally:

```
git clone https://github.com/neurodata/r-mgc.git
```

We assume that the directory `r-mgc` placed locally on the system is `<package_root>`. Note that all figures were stylized using [Adobe Photoshop](www.adobe.com/Photoshop) prior to submission.

- [Figure 1: Mini Sims Figure](https://github.com/neurodata/r-mgc/tree/master/docs/discriminability/paper/simulations/dummy_sims.Rmd) This figure demonstrates the behavior of discriminability, Fingerprinting, ICC/I2C2, and Kernel methods under a range of basic simulation settings in 1 dimension.

- [Figure 2: Multisim Figure](https://github.com/neurodata/r-mgc/tree/master/docs/discriminability/paper/simulations) This figure demonstrates the behavior of discriminability, ICC, and I2C2 under a variety of simulation benchmark settings. To execute the script with fresh data:

```
setwd('<package_root>/docs/discriminability/paper/simulations')
source('shared_scripts.R`)
```

Note: the scripts will automatically multithread, however, the simulation benchmarks take quite a while to execute (1.5 days on a 96 core machine with 1 TB of RAM).

Using the included [bound](https://github.com/neurodata/r-mgc/blob/master/docs/discriminability/paper/data/sims/discr_sims_bound.rds), [one sample](https://github.com/neurodata/r-mgc/blob/master/docs/discriminability/paper/data/sims/discr_sims_os.rds), and [two sample](https://github.com/neurodata/r-mgc/blob/master/docs/discriminability/paper/data/sims/discr_sims_ts.rds) data, you can proceed to duplicate the figure by opening the R notebook [simulation plots](https://github.com/neurodata/r-mgc/blob/master/docs/discriminability/paper/simulations/multisim_figure.Rmd), and executing the script.

- [Figure 3: 64 pipelines figure](https://github.com/neurodata/r-mgc/tree/master/docs/discriminability/paper/64pipes_fig). To regenerate the source data for this portion of the manuscript, users can use the following two scripts from an R terminal:

```
setwd('<package_root>/docs/discriminability/paper/discr_computation')
# edit lines 17 and 18, and lines 210 and 211, and set to your local path where
# preprocessed brains are located
source('./real_data_driver.R')  # runs the discriminability calculations
# edit lines 17 and 18, and lines 108 and 109, to the location of the 
# preprocessed brains
source('./realdat_perm_testing.R')  # runs the two sample testing
```
Again, the scripts will multithread, but can be expected to take approximately 3 days on a 96 core, 1 TB RAM machine.

To regenerate Figure 2 from the manuscript, users can execute the [64 Pipelines Figure](https://github.com/neurodata/r-mgc/blob/master/docs/discriminability/paper/64pipes_fig/real_data.Rmd) notebook.

- [Figure 4: Marginalized Options Comparison](https://github.com/neurodata/r-mgc/blob/master/docs/discriminability/paper/multi_modal_opts') Users can regenerate the figure by using the notebook [Multi Modal Opts](https://github.com/neurodata/r-mgc/blob/master/docs/discriminability/paper/multi_modal_opts/multi_modal_opts.Rmd).

- [Figure 5: Effect Size Investigation](https://github.com/neurodata/r-mgc/edit/master/docs/discriminability/paper/dcor_fig) Users can reproduce the data collected with:

```
setwd('<package_root>/docs/discriminability/paper/dcor_fig')
source('./dep_wt_driver.R')
```

Results can be expected to take 2 days on a 96 core, 1 TB machine.

To reproduce the figure, users can use the [Effect Size Investigation](https://github.com/neurodata/r-mgc/blob/master/docs/discriminability/paper/dcor_fig/dcor_bypipe_exps.Rmd) notebook.
