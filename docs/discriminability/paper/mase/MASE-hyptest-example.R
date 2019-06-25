# Load packages and MASE source code
source("R/loadAll.R")

# Generate a pair of mixed-membership SBM graphs with same community structure
m = 2 # sample size
n = 300 # number of vertices
d=3 # common subspace dimension
p = 0.4 # within community edge probability
q = 0.1 # between community edge probability

library(MCMCpack)
# Sample mixed-memberships
Z <- rdirichlet(n, alpha = rep(0.1, d))


b = 500 #bootstrap samples
numsims = 100

# Vector of differences between graphs 
epsilons = seq(0, 0.1, 0.01)

p.values.list <- list()

i = 1

# Equal distribution on the two graphs ----------------------------------------------
epsilon = 0
B1 = diag(3)*(p-q) + q
B2 = B1
B2[1, 1] = B2[1, 1] + epsilon

# Edge probabilities of the graphs
P1 = Z %*% B1 %*% t(Z)
P2 = Z %*% B2 %*% t(Z)
Ps.true = list(P1, P2)
# Obtain COSIE parameters from true matrices
mase.true = mase(Ps.true, d)
V = mase.true$V
R1 = mase.true$R[[1]]
R2 = mase.true$R[[2]]

# Sample two graphs
set.seed(1989)
A1 = sample_from_P(P1)
A2 = sample_from_P(P2)

# MASE algorithm
mase.embed = mase(list(A1, A2), d = 3)
mase.embed$R

MASE.test.stat = norm(mase.embed$R[[1]] - mase.embed$R[[2]], "F")

# Calculate null distribution of the test
mase.tstat.null.true1 <- mase.test.statistic.distr(V = V, R = R1, B = b)
mase.tstat.null.true2 <- mase.test.statistic.distr(V = V, R = R2, B = b)

# Obtain conservative p-value 
mase.pval.equal.distr = max(p.values(MASE.test.stat, mase.tstat.null.true1),
                p.values(MASE.test.stat, mase.tstat.null.true2))
print(mase.pval.equal.distr)


# Different distribution on the two graphs -------------------------------------------
epsilon = 0.05
B1 = diag(3)*(p-q) + q
B2 = B1
B2[1, 1] = B2[1, 1] + epsilon

# Edge probabilities of the graphs
P1 = Z %*% B1 %*% t(Z)
P2 = Z %*% B2 %*% t(Z)
Ps.true = list(P1, P2)
# Obtain COSIE parameters from true matrices
mase.true = mase(Ps.true, d)
V = mase.true$V
R1 = mase.true$R[[1]]
R2 = mase.true$R[[2]]

# Sample two graphs
A1 = sample_from_P(P1)
A2 = sample_from_P(P2)

# MASE algorithm
mase.embed = mase(list(A1, A2), d = 3)
mase.embed$R

MASE.test.stat = norm(mase.embed$R[[1]] - mase.embed$R[[2]], "F")

# Calculate null distribution of the test
mase.tstat.null.true1 <- mase.test.statistic.distr(V = V, R = R1, B = b)
mase.tstat.null.true2 <- mase.test.statistic.distr(V = V, R = R2, B = b)

# Obtain conservative p-value 
mase.pval.diff.distr = max(p.values(MASE.test.stat, mase.tstat.null.true1),
                p.values(MASE.test.stat, mase.tstat.null.true2))
print(mase.pval.diff.distr)
