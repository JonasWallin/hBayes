library(NPBayes)
graphics.off()
c <- -0.1
sim <- 4000
x.grid <- seq(pi/2,pi,length.out=100)
x.grid.corse <- seq(pi/2,pi,length.out=50)

# TODO: start timing (not needed!!)
# TODO: Do 5x5 gibbs sampler
# TODO: Do acf and speed analysis




b <- c(1,0.1)
a <- c(1,3.1)
E <- rbing.2diagmatrix(a, b)
A <- E%*%diag(a)%*%t(E)
E.A <- rbing.2(A, b)
E.A.rstiefel<-rbing.Op(-A,diag(b))
