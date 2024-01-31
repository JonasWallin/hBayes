##
# student exam score
#
##
library(NPBayes)
library(devtools)
graphics.off()
#devtools::load_all()
mcmc.samples = 20000
X <- as.matrix(read.table('../../data/student_score.txt',header=T))
X <- scale(X, center=T,scale=F)
n <- dim(X)[1]
E <-  eigen(t(X)%*%X/n)
set.seed(1)
res <- oracle.metrop.sampler(E$vectors,X, E$values, mcmc.samples)
set.seed(1)
res_old <- oracle.metrop.sampler.old(E$vectors,X, E$values, mcmc.samples)
plot(cumsum(res$acc[,2])/(1:length(res$acc[,2])),ylim=c(0,1),type='l')
for(i in 3:5){
    lines(cumsum(res$acc[,i])/(1:length(res$acc[,i])))
}
for(i in 2:5){
    cat('sigma[',i,'] = ', round(res$MH.objs[[i]]$sigma,3),' n = ',res$MH.objs[[i]]$n.SA,'\n')
}
print(max(abs(res$loglik-res_old$loglik)))


