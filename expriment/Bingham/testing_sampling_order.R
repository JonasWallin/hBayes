#' Comparing samplers ordering
#'
#'
graphics.off()
values <- (10:1)^2
n <- 100
eigA <- list(values = 1/values,
     vectors = diag(nrow=length(values),ncol=length(values)))
eigB <-  list(values = -values*n/2,
             vectors = diag(nrow=length(values),ncol=length(values)))

p <- length(values)
n.mcmc <- 1000 * p
res.iter <- rbing.iter(n.mcmc/p, NULL, NULL,eigA = eigA, eigB = eigB,ret = 2)
res      <- rbing(n.mcmc, NULL, NULL,eigA = eigA, eigB = eigB,ret = 2)
plot(abs(res.iter$Es[1,1,]))
points(abs(res$Es[1,1,]), col='red')
acf(abs(res.iter$Es[5,5,seq(1,n.mcmc,by=p)]),10)
acf(abs(res$Es[5,5,seq(1,n.mcmc,by=p)]),10)
