# Multiple Adjacency Spectral Embedding


The *multiple adjacency spectral embedding* (MASE) performs a joint embedding of graphs based on a common invariant subspace assumption (see [Arroyo et al. (2019)](https://arxiv.org/abs/1906.10026)).

# Overview

Given a sample of *m* graphs <img src="https://latex.codecogs.com/gif.latex?A_1,&space;\ldots,&space;A_m" title="A_1, \ldots, A_m" />, each with *n* aligned vertices, MASE obtains a joint embedding of their adjacency matrices. The embedding of the graphs is based on the *common subspace independent-edge* (COSIE) random graph model, which assumes that the expected adjacency of all the graphs share the same  invariant subspace of dimension <img src="https://latex.codecogs.com/gif.latex?d" title="d" /> . Under this model, the expected adjacency of a graph <img src="https://latex.codecogs.com/gif.latex?\hat{A}_i" title="\hat{A}_i" /> is given by

<img src="https://latex.codecogs.com/svg.latex?\Large&space;E[A_i]=VR_iV^T," title="\Large E[A_i]=VR_iV^T," />

where <img src="https://latex.codecogs.com/gif.latex?V\in\mathbb{R}^{n\times&space;d}" title="V\in\mathbb{R}^{n\times d}" />
 is a matrix with orthogonal columns that represents the basis of the common invariant subspace of the model, and <img src="https://latex.codecogs.com/gif.latex?R_1,&space;\ldots,&space;R_m\in\mathbb{R}^{d\times&space;d}" title="R_1, \ldots, R_m\in\mathbb{R}^{d\times d}" /> are individual *score matrices* that represent each graph.
 
To estimate the parameters of the embedding, MASE calculates a separate *adjacency spectral embedding*  (ASE) for each graph, which consists in computing the eigendecomposition of each adjacency matrix, with a possible eigenvalue scaling, followed by a joint singular value decomposition of the concatenated ASEs. This procedure obtains an estimated common invariant subspace matrix <img src="https://latex.codecogs.com/gif.latex?\hat{V}" title="\hat{V}" /> (depicted in the picture below). After that, each adjacency matrix <img src="https://latex.codecogs.com/gif.latex?\hat{A}_i" title="\hat{A}_i" />
is projected back to the eigenspace to obtain an estimated score matrix  <img src="https://latex.codecogs.com/gif.latex?\hat{R}_i" title="\hat{R}_i" />

![mase embedding](https://raw.githubusercontent.com/jesusdaniel/mase/master/img/mase1.jpg)

The COSIE model and the MASE algorithm can be deployed for a number of subsequent network inference tasks, including
+ Graph hypothesis testing
+ Multiple graph dimensionality reduction
+ Joint community detection on multiple networks
+ Graph classification

# R Code
This repository contains R code that implements MASE and other auxiliary functions. To use MASE, download all the content from the [R/](https://github.com/jesusdaniel/mase/tree/master/R) folder. The file [MASE-hyptest-example.R](https://github.com/jesusdaniel/mase/blob/master/MASE-hyptest-example.R) contains an example of MASE for two-sample graph hypothesis testing.

# References

Arroyo, J., Athreya, A., Cape, J., Chen, G., Priebe, C.E., Vogelstein, J.T., "Inference for multiple heterogeneous networks with a common invariant subspace", 
[![arXiv shield](https://img.shields.io/badge/arXiv-1906.10026-red.svg?style=flat)](https://arxiv.org/abs/1906.10026)
