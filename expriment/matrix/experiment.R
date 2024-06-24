p <- 24
n	<- 40
Gamma <- smp.rnd.orth.mat(p)
Lambda  <- sort(rgamma(p,shape = 10,scale = 10),dec = T)
Sigma <- Gamma %*% diag(Lambda) %*% t(Gamma)      # General covariance matrix
