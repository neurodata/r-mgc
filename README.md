# mgc
R package for MGC code


## Installation:

Set up `devtools` package for R:

```
sudo apt-get -y install libcurl4-openssl-dev libssl-dev libxml2-dev
```

Set up the `fmriutils` package:

Next, we can install the R package dependencies directly. From an R terminal (or a GUI such as Rstudio), type:
```
install.packages(c('devtools', 'Rcpp', 'RSpectra', 'stringr', 'ggplot2', 'reshape2', 'abind', 'plyr'))
require(devtools)
install_github('igraph/rigraph')  # installs latest version of igraph
install_github('neurodata/fmriutils')
require('fmriutils')
```
Install the `MGC` package:

Still from your R terminal (or GUI such as RStudio), type:

```
install_github('ebridge2/mgc')
```

## Demo

### Discriminability

Check out our demo tutorial for [discriminability-brain analysis]().

### MGC
Into your R console, type the following:

```
library(MGC)
run_demo()  # runs the MGC demo
```
which takes < 10 seconds to run, which shows the actual MGC and p-value for testing on linear and quadratic data.
