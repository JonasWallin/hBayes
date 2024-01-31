simulateDesignMatrix <- function(chr, m, n, d, seed = NA) {
  ## generate a design matrix
  # chr -- number of chromosomes
  # m -- number of markers on each chromosome
  # n -- number of observations
  # d -- distance between markers (cM)
  # seed -- seed for a random number generator
  if (length(m) == 1) m <- rep(m, chr)
  m_cum <- c(0, cumsum(m))
  if(is.na(seed)==F)
    set.seed(seed)
  r <- (1 - exp(-2*d/100)) / 2
  X <- matrix(0, n, sum(m))
  for(i in 1:chr) {
    x <- 2 * sample(2, n, rep = TRUE) - 3 # szybsze od c(-1, 1)
    x <- matrix(rep(x, m[i]), ncol = m[i])
    z <- 2 * sample(2, n*m[i], rep = TRUE, prob = c(r, 1 - r)) - 3
    z <- cumprod(z)
    z <- t(matrix(z, ncol = n))
    X[, (m_cum[i]+1):m_cum[i+1]] <- (x*z + 1)/2
  }
  return(2 * X - 1) # -1/1 coding
}

simulateTraits <- function(X, 
                           q, 
                           mu, 
                           tau, 
                           cis = NULL, 
                           qtl = c(0, 1, 1),
                           qtl2 = c(0, 1, 1), 
                           qtl3 = c(0, 1, 1),
                           mixedd = 'normal') {
  ## generate vector of q traits for a given design matrix (X) with polygenci effects
  ## and up to 3 true QTLs
  # mu, tau -- polygenic effect (mean and sd)
  # cis -- c(mean, sd)
  # qtl -- c(position: marker, trait; beta)
  M <- ncol(X)
  if(mixedd == "normal"){
    beta.matrix <- matrix(rnorm(q * M, mu, tau), ncol = q, byrow = TRUE)
  }else if(mixedd  == "Laplace"){
    beta.matrix <- matrix(rep(mu,M) + sqrt(rgamma(q*M, 1)) * rnorm(q * M, 0, tau), ncol = q, byrow = TRUE)
  }
  if (!is.null(cis)) {
    cis.matrix <- diag(rnorm(M, cis[1], cis[2]))
  } else cis.matrix <- 0
  qtl.matrix <- matrix(0, M, q)
  qtl.matrix[qtl[2], qtl[3:length(qtl)]] <- qtl[1]
  qtl2.matrix <- matrix(0, M, q)
  qtl2.matrix[qtl2[2], qtl2[3:length(qtl2)]] <- qtl2[1]
  qtl3.matrix <- matrix(0, M, q)
  qtl3.matrix[qtl3[2], qtl3[3:length(qtl3)]] <- qtl3[1]
  e.matrix <- replicate(q, rnorm(n, 0, sigma))
  y <- X %*% (beta.matrix + cis.matrix + qtl.matrix + qtl2.matrix + qtl3.matrix) + e.matrix
  return(y)
}

findDuplicate <- function(X) {
  ## find which columns are duplicated
  dupl.ind <- which(duplicated(X, MARGIN = 2)) # these columns are duplicated...
  dupl.repl <- which(duplicated(X, MARGIN = 2, fromLast = TRUE)) # ...with these columns
  dupl <- cbind(dupl.ind, dupl.repl)
  return(dupl)
}

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

replaceOneNA <- function(na, x) {
  ## replace one missing value
  if (!is.na(x[na])) stop("Value is not NA.")
  n <- length(x)
  left <- which(!is.na(x[na:1]))[1] - 1
  right <- which(!is.na(x[na:n]))[1] - 1
  if (is.na(left) & is.na(right)) stop("All values are NA.")
  if (is.na(left)) left <- Inf else if (is.na(right)) right <- Inf
  if (left < right) {
    replace <- x[na - left]
  } else if (left > right) {
    replace <- x[na + right]
  } else {
    replace <- (x[na - left] + x[na + right])/2
  }
  return(replace)
}

pseudovalue <- function(m1, m2, d1, d2){
  ## calculate expected value of marker m between m1 and m2 (conditional)
  ## for m in {-1, 0, 1}
  r1 <- (1 - exp(-2*d1/100)) / 2
  r2 <- (1 - exp(-2*d2/100)) / 2
  r <- r1 + r2 - 2*r1*r2
  num <- (1 - r1)^(1 + m1) * r1^(1 - m1) * (1 - r2)^(1 + m2) * r2^(1 - m2) -
    (1 - r1)^(1 - m1) * r1^(1 + m1) * (1 - r2)^(1 - m2) * r2^(1 + m2)
  den <- r^abs(m1 - m2) * (1 - r)^(2 - abs(m1 - m2))
  return(num/den)
}

pseudomarkers <- function(X, markers, len_cum, by = 2){
  ## put pseudomarkers between the first and last markers
  # markers -- vector, number of markers on each chromosome
  # len_cum -- distance from the begining of chromosome
  n <- nrow(X)
  chr <- length(markers)
  markers_cum <- c(0, cumsum(markers))
  # first marker on a given chromosome must be at the distance of 0:
  len_cum <- len_cum - rep(len_cum[markers_cum[-(chr + 1)] + 1], markers)
  chr_len <- len_cum[cumsum(markers)] # lengths of chromosomes
  markers <- c(0, cumsum(markers))
  new_markers <- numeric(chr)
  XX <- NULL
  l <- 1
  for(i in 1:chr) {
    new_markers[i] <- ceiling(chr_len[i] / by)
    XX <- cbind(XX, matrix(0, n, new_markers[i]))
    for(j in 1:new_markers[i]) {
      m <- sum(len_cum[(markers[i] + 1) : markers[i + 1]] <= (j - 1)*by) + markers[i]
      m1 <- X[, m]
      m2 <- X[, m + 1]
      d1 <- (j - 1)*by - len_cum[m]
      d2 <- len_cum[m + 1] - (j - 1)*by
      XX[, l] <- pseudovalue(m1, m2, d1, d2)
      l <- l + 1
    }
  }
  return(list(X = XX, markers = new_markers))
}
