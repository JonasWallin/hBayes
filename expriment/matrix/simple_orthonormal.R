graphics.off()
#library(NPBayes)
library(rstiefel)
library(far)
library(MASS)
n.mcmc <- 2000
n	<- 1

set.seed(2)
sample.rstiefel = F
Lambda  <- diag(c(6.5,1,0.001))#diag(c(100,rep(10,p/2),rep(1,p/2)))
p <- dim(Lambda)[1]
E <- diag(p)
E2 <- matrix(c(0.59857780 , 0.80105315,
               -0.00429662,  0.80016093,
               -0.59814860 ,-0.04428022,
               -0.03804083 , 0.02306717, -0.99900991),
             3,3)
E2 <- orthonormalization(E2)
#E2 <- E
E3 <- diag(p)
E.vec <- array(NA,dim=c(p,p,n.mcmc))
E2.vec <- array(NA,dim=c(p,p,n.mcmc))
Z2.vec <-array(NA,dim=c(2,n.mcmc))
ind2.vec <- array(NA,dim=c(2,n.mcmc))
E3.vec <- array(NA,dim=c(p,p,n.mcmc))


Lambda_hatE <- list(values = 1/diag(Lambda),
                    vectors = diag(nrow=p,ncol=p))
XXtE <- list(values = -diag(Lambda)*n/2,
             vectors = diag(nrow=p,ncol=p))
lik.1 <-lik.2 <- lik.3 <- rep(0,n.mcmc)
R.XtX <- diag(sqrt(diag(n*Lambda)))
X <- matrix(0,nrow=n,ncol=p)
res.E.o <- oracle.metrop.sampler(E2,
                                 X=X,
                                 diag(Lambda),
                                 10,
                                 target.acc = 0.23,
                                 R = R.XtX)
E <- res.E.o$Es[,,10]
for(i in 1:n.mcmc){
    cat('*')
    res.E.o <- oracle.metrop.sampler(E,
                                   X=X,
                                   diag(Lambda),
                                   10,
                                   target.acc = 0.23,
                                   MH.objs = res.E.o$MH.objs,
                                   R = R.XtX)
    E <- res.E.o$Es[,,10]
    E.vec[,,i] <- E
    #oracle.gamma.R.loglik(R.XtX,  mcmc.gamma[,,1],mcmc.lambda[1, ], n)
    lik.1[i] <- oracle.gamma.R.loglik(R.XtX,  E,diag(Lambda), n)

    res.E <- rbing(1,
                        A = NULL,
                        B = NULL,
                        eigA=XXtE ,
                        eigB=Lambda_hatE,
                       E0=E2,
                        EtAE=NULL,
                        ret = 2)
    E2 <- res.E$Es[,,1]
    Z2.vec[,i] <- res.E$Z
    ind2.vec[,i] <- res.E$inds
    #E2 <-  rbing.2diagmatrix(1/diag(Lambda),-diag(Lambda)*n/2) #res.E$Es[,,p]
    lik.2[i] <- oracle.gamma.R.loglik(R.XtX,  res.E$Es[,,1],diag(Lambda), n)
    E2.vec[,,i] <- E2

    if(sample.rstiefel){
        E3 <-rbing.matrix.gibbs(diag(1/diag(Lambda)),-0.5*n*Lambda, E3)
        lik.3[i] <- oracle.gamma.R.loglik(R.XtX,  E3,diag(Lambda), n)
        E3.vec[,,i] <- E3
    }
}

cat('\n')
# Run sampler
#res <- Eigen.sampler(x.smp,n.mcmc=n.mcmc,n.E.mcmc = 1e1,L = 6, sample.eigenvalues=F)

#res_old <- Eigen.sampler.old(x.smp,n.mcmc=n.mcmc,n.E.mcmc = 1e1,L = 6, sample.eigenvalues=F)

# Display eigenvalue estimation

 if(sample.rstiefel){
    par(mfcol = c(1,3))
     plot(lik.1,type='l')
     plot(lik.2,type='l')
     plot(lik.3,type='l')
 }else{
     par(mfcol = c(1,2))
     plot(lik.1,type='l')
     plot(lik.2,type='l')
 }
#hist.seq <- seq(min(lik.1,lik.2),max(lik.1,lik.2),length.out=40)
#res.hist <- hist(lik.1,breaks=hist.seq)

#hist((E2.vec[2,p,]),breaks = seq(-1,1,length.out=20))
#hist(lik.2,breaks=hist.seq)
#hist((E.vec[2,p,]),breaks = seq(-1,1,length.out=20))
