# mgc
R package for MGC code


## Installation:

Set up `devtools` package for R:

```
sudo apt-get -y install libcurl4-openssl-dev libssl-dev libxml2-dev
```

Optionally, set up the `fmriutils` package (required for several tutorials) following the [fmriutils-repo](https://github.com/neurodata/fmriutils) install instructions.

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
