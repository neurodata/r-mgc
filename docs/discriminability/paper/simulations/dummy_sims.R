## Dummy Sims

require(ggplot2)
require(reshape2)
require(grid)
require(gridExtra)
require(mgc)
require(ICC)
require(I2C2)
require(cowplot)
require(scales)
require(latex2exp)
require(MASS)
source('./shared_scripts.R')
require(dplyr)

n = 10

data <- lapply(c("(i) Discriminable", "(ii) Offset", "(iii) Outlier"), function (sim) {
  X.d <- matrix(NaN, nrow=n*2, ncol=1)
  if (sim == "(ii) Offset") {
    centers <- seq(-1, 1, length.out=n)
    ids <- 1:n
    dist <- 2/n
    X.d[,1] <- c(centers + dist, centers - dist)
    Y <- c(ids, ids)
    Z <- c(rep(1, n), rep(2, n))
  }
  else if (sim == "(iii) Outlier") {
    centers <- seq(-1, 1, length.out=n)
    ids <- 1:n
    dist <- .5/n
    X.d[,1] <- c(centers + dist, centers - dist)
    rpts <- c(3, 17)
    Y <- c(ids, ids)
    set.seed(1234)
    X.d[rpts,1] <- X.d[rpts,1] + 8
    Z <- c(rep(1, n), rep(2, n))
  }
  else if (sim == "(i) Discriminable") {
    centers <- seq(-1, 1, length.out=n)
    ids <- 1:n
    dist <- .5/n
    X.d[,1] <- c(centers + dist, centers - dist)
    Z <- c(rep(1, n), rep(2, n))
    Y <- c(ids, ids)
  }
  else if (sim == "(iv) Unreliable") {
    centers <- rep(0, length.out=n)
    ids <- 1:n
    dist <- .5/n
    X.d[,1] <- c(centers + dist, centers - dist)
    Z <- c(rep(1, n), rep(2, n))
    Y <- c(ids, rev(ids))
  }
  return(data.frame(X=X.d, Y=Y, Z=Z, Simulation=sim))
})
names(data) <- c("(i) Discriminable", "(ii) Offset", "(iii) Outlier")

results <- lapply(names(data), function(sim) {
  dat <- data[[sim]]
  icc <- icc.os(dat$X, dat$Y)
  discr <- discr.os(dat$X, dat$Y, is.dist=FALSE)
  fpi <- fpi.os(matrix(dat$X, ncol=1), dat$Y, dat$Z, dist.xfm=dist.x, is.sim=FALSE)
  disco <- disco.os(dat$X, dat$Y, is.dist=FALSE)
  data.frame(Statistic=c(icc, discr, fpi, disco), Simulation=sim,
             Algorithm=c("ICC", "Discr.", "Finger.", "Kernel"))
})

var.ord <- c(sapply(1:10, function(i) rep(i, 2)))
Dmtx_results <- lapply(names(data), function(sim) {
  dat <- data[[sim]]
  new_ord <- sort(dat$Y, index.return=TRUE)$ix
  Dmtx <- as.matrix(dist(matrix(dat$X[new_ord], ncol=1))) %>%
    melt() %>%
    dplyr::rename(Distance=value) %>%
    dplyr::mutate(Simulation = sim) %>%
    dplyr::mutate(Individual1=var.ord[Var1], Individual2=var.ord[Var2])
})

do.call(rbind, data) %>%
  group_by(Simulation) %>%
  mutate(X=(X - min(X))/(max(X) - min(X))) %>%
  ungroup() %>%
  mutate(Simulation=factor(Simulation,
                           levels=c("(i) Discriminable", "(ii) Offset", "(iii) Outlier"), ordered=TRUE),
         Z=factor(Z)) %>%
  filter(Simulation == "(ii) Offset") +
  ggplot(aes(x=X - min(X)/max(X), y=as.factor(Y),
             fill=as.factor(Y), shape=Z)) +
  geom_point(color='white', shape=21) +
  #scale_shape_manual(values=c('1'=22, '2'=23), guide=FALSE) +
  xlab("Value") +
  ylab("Individual ID") +
  ggtitle("(ii) Offset Simulation") +
  scale_x_continuous(breaks=c(0, 0.5, 1), labels=c(-1, 0, 1)) +
  scale_fill_discrete(guide=FALSE) +
  scale_y_discrete(breaks=c(1, 10), labels=c(1, 10)) +
  theme_bw(base_size=18)

