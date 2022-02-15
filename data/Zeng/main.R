source("../../symulacje/2-models.R")
source("../../symulacje/3-exp_matrix.R")
library(RcppEigen)
library(ggplot2)
library(bigstep)
replaceNA <- function(X) {
  ## replace missing values with the mean
  ## (function uses replaceOneNA)
  na_pos <- apply(X, 1, function(x) which(is.na(x)))
  for (i in seq_along(na_pos)) {
    if (length(na_pos[[i]]) == 0) next
    X[i, na_pos[[i]]] <- sapply(na_pos[[i]], replaceOneNA, x = X[i, ])
  }
  return(X)
}
# data
d1 <- read.table("bm4zb.out", skip = 9, row.names = 1)
d2 <- read.table("bm6zb.out", row.names = 1)
d3 <- read.table("bs4zb.out", skip = 9, row.names = 1)
d4 <- read.table("bs6zb.out", row.names = 1)

# X
#X <- as.matrix(rbind(d1[, 1:45], d2[, 1:45], d3[, 1:45], d4[, 1:45]))
#X <- as.matrix(d1[, 1:45])
X <- as.matrix(rbind(d1[, 1:45], d2[, 1:45]))
X <- X[, -(1:6)] # odrzucam 1. chromosom X
X[X == 9] <- NA
X <- replaceNA(X) - 1
# group <- rep(c("bm4", "bm6", "bs4", "bs6"), c(nrow(d1), nrow(d2), nrow(d3), nrow(d4)))
#X <- X[group %in% c("bm4", "bm6"), ]
#X <- X[group %in% c("bs4", "bs6"), ]
n <- nrow(X)
M <- ncol(X)

# y
y <- as.matrix(c(d1[, 46], d2[, 46], d3[, 46], d4[, 46]), ncol = 1)
#y <- as.matrix(d1[, 46], ncol = 1)
y <- as.matrix(c(d1[, 46], d2[, 46]), ncol = 1)
y <- scale(y)
#y <- y + 1 * X[, 38]
q <- 1
#y <- y[group %in% c("bm4", "bm6"), , drop = FALSE]
#y <- y[group %in% c("bs4", "bs6"), , drop = FALSE]

# map
markers <- c(6, 16, 23)[-1] # bez 1. chromosomu
chr <- length(markers)
len <- c(0, 3.60, 10.60, 9.20, 17.20, 18.70, 0, 6.98, 10.10, 4.94, 6.51, 6.19,
         20.46, 12.78, 3.90, 4.55, 7.49, 30.02, 16.85, 4.34, 3.71, 7.03, 0, 4.99,
         9.34, 6.97, 7.44, 14.46, 6.79, 3.55, 6.32, 11.86, 4.58, 6.85, 6.35,
         11.79, 12.88, 9.15, 3.30, 7.98, 13.09, 10.04, 3.70, 9.79, 3.43)[-(1:6)]
len_cum <- unlist(tapply(len, rep(1:chr, markers), cumsum), use.names = FALSE)

# pseudomarkers
by <- 2
res <- pseudomarkers(X, markers, len_cum, by)
X <- res$X
M <- ncol(X)
markers <- res$markers

source("../../symulacje/4-help_functions.R")

findDuplicate2 <- function(X) {
  X[X < -0.7] <- -1
  X[X > 0.7] <- 1
  X[X > -0.3 & X < 0.3] <- 0
  dupl.ind <- which(duplicated(X, MARGIN = 2)) # ktore kolumny sie powtarzaja
  dupl.repl <- which(duplicated(X, MARGIN = 2, fromLast = TRUE)) # z tymi sie powtarzaja
  dupl <- cbind(dupl.ind, dupl.repl)
  return(dupl)
}
dupl <- findDuplicate2(X)

#X <- scale(X) # problem

# ----------------------------------------------------------------------
# simple model
pb <- txtProgressBar(min = 0, max = q, style = 3)
for(i in 1:q) {
  t.simple[i, ] <- simpleModel(X, y[, i])
  crit.simple[i] <- calculateCrit(t.simple[i, ], markers)
  find.simple[[i]] <- which(abs(t.simple[i, ]) > crit.simple[i])

  res <- simpleModelForward(X, y[, i], markers, dupl = NULL)
  t.simple.forward[i, ] <- res$t
  find.simple.forward[[i]] <- res$find
  crit.simple.forward[i] <- res$crit
  setTxtProgressBar(pb, i)
}
close(pb)

plotTrait(abs(t.simple[1, ]), markers, crit.simple[1])
# pdf("zeng_simple_forward.pdf", w = 5, h = 3)
plotTrait(abs(t.simple.forward[1, ]), markers, crit.simple.forward[1])
# dev.off()

# ----------------------------------------------------------------------
# CIM model
t.CIM <- matrix(0, q, M)
crit.CIM <- numeric(q)
find.CIM <- list()
pb <- txtProgressBar(min = 0, max = q, style = 3)
for(i in 1:q) {
  t.CIM[i, ] <- CIMModel(X, y[, i], window = 10)
  crit.CIM[i] <- calculateCrit(t.CIM[i, ], markers)
  find.CIM[[i]] <- which(abs(t.CIM[i, ]) > crit.CIM[i])
  setTxtProgressBar(pb, i)
}
close(pb)

# pdf("zeng_CIM.pdf", w = 5, h = 3)
plotTrait(abs(t.CIM[1, ]), markers, crit.CIM[1])
# dev.off()

# ----------------------------------------------------------------------
# mixed model
pb <- txtProgressBar(min = 0, max = q, style = 3)
for(i in 1:q) {
  res <- mixedModel(X, y[, i], withoutTau = FALSE)
  t.mix[i, ] <- res$t
  crit.mix[i] <- calculateCrit(t.mix[i, ], markers)
  find.mix[[i]] <- which(abs(t.mix[i, ]) > crit.mix[i])
  tau.est[i] <- res$tau.est
  sigma.est[i] <- res$sigma.est

  res <- mixedModelForward(X, y[, i], markers = markers, dupl = NULL, withoutTau = FALSE)
  t.mix.forward[i, ] <- res$t
  find.mix.forward[[i]] <- res$find
  tau.est.forward[i] <- res$tau.est
  sigma.est.forward[i] <- res$sigma.est
  ignoreTau[i] <- res$ignoreTau
  ignoreS[i] <- res$ignoreS
  crit.mix.forward[i] <- res$crit

  # res <- mixedModel(X, y[, i], sigma = 0.213, tau = 0.0023)
  # t.mix.known[i, ] <- res$t
  # find.mix.known[[i]] <- which(abs(t.mix.known[i, ]) > calculateCrit(t.mix.known[i, ], markers))
  #
  # res <- mixedModelForward(X, y[, i], markers = markers, dupl = dupl, sigma = 0.213, tau = 0.0023)
  # t.mix.forward.known[i, ] <- res$t
  # find.mix.forward.known[[i]] <- res$find
  setTxtProgressBar(pb, i)
}
close(pb)

plotTrait(abs(t.mix[1, ]), markers, crit.mix[1])
# pdf("zeng_mixed_forward.pdf", w = 5, h = 3)
plotTrait(abs(t.mix.forward[1, ]), markers, crit.mix.forward[1])
# dev.off()

# H2 = VarG/VarY = VarG (standaryzacja) = VarQTL + VarPoly = VarY - sig^2 = 1 - sig^2

# H2 efektu genetycznego dla mixed (QTL + poly)
1 - sigma.est.forward^2
1 - sigma.est.forward^2 - 0.041 # var(Xm[, 1:4] %*% beta[1:4]) w mixedModelForward

# dla simple forward
1 - 0.2735095 # fastLmPure(Xm, y)$s^2 w simpleModelForward
