library(igraph)
library(lattice)
library(graphclass) #github jesusdaniel/graphclass (only needed for the plots)
library(ggplot2)
library(multiRDPG)
library(Matrix)
source("R/mase.R")

source("R/parametric-bootstrap.R")
source("R/sample_graphs.R")


# read getElbows function from github
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/youngser/gmmase/master/R/getElbows.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
