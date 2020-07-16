## Testing for Fingerprint Index
# Data perfectly indexable
require(abind)

X.base <- matrix(0, nrow = 20, ncol = 10)
X.base[1:10,] <- cbind(rep(1, 10), rep(-1, 10), rep(1, 10), rep(-1, 10),
                  rep(1, 10), rep(-1, 10), rep(1, 10), rep(-1, 10),
                  rep(1, 10), rep(-1, 10))
X.base[11:20,] <- X.base[1:10,]
X <- X.base + matrix(0.1*rnorm(200), nrow=20)
Y <- c(rep("i1", 10), rep("i2", 10))
Z <- rep(sapply(1:10, function(i) sprintf("ses%s", i)), 2)

# should be near 0.5
fpi.os(X, Y, Z)

X.base[11:20,] <- -X.base[1:10,]
X <- X.base + matrix(0.1*rnorm(200), nrow=20)
# should be near 1
fpi.os(X, Y, Z)

# should be somewhere in btwn
X <- X.base + matrix(2*rnorm(200), nrow=20)
fpi.os(X, Y, Z)

## 3 Subject Test
X.base <- matrix(0, nrow = 30, ncol = 10)
X.base[1:10,] <- cbind(rep(1, 10), rep(-1, 10), rep(1, 10), rep(-1, 10),
                       rep(1, 10), rep(-1, 10), rep(1, 10), rep(-1, 10),
                       rep(1, 10), rep(-1, 10))
X.base[11:20,] <- X.base[1:10,]
X.base[21:30,] <- X.base[1:10,]
X <- X.base + matrix(0.1*rnorm(300), nrow=30)
Y <- c(rep("i1", 10), rep("i2", 10), rep("i3", 10))
Z <- rep(sapply(1:10, function(i) sprintf("ses%s", i)), 3)
# should be near 0.33
fpi.os(X, Y, Z)

X.base[11:20,] <- -X.base[1:10,]
X.base[21:30,] <- cbind(rep(1, 10), rep(1, 10), rep(-1, 10), rep(1, 10),
                        rep(1, 10), rep(-1, 10), rep(1, 10), rep(1, 10),
                        rep(-1, 10), rep(1, 10))
X <- X.base + matrix(0.1*rnorm(300), nrow=30)
# should be near 1
fpi.os(X, Y, Z)

X <- X.base + matrix(2*rnorm(300), nrow=30)
# should be somewhere in btwn
fpi.os(X, Y, Z)

