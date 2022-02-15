library(NPBayes)
L	<- 6
I	<- 2^L

log.odds <- function(p) {
    log(p / (1 - p))
}
inv.log.odds <- function(theta)	{
    exp(theta) / (1 + exp(theta))
}

tree.parent.lst	<- as.list(1:L)
for (l in 1:L)
{
    if (l == 1)	tree.parent.lst[[l]]	<- as.list(0)
    if (l > 1)	tree.parent.lst[[l]]	<- as.list(rep(1:(2^(l - 2)),each = 2))
    names(tree.parent.lst[[l]])	<- paste("Level", l, " -- Node", 1:(2^(l - 1)))
}


tree.child.node.lst		<- as.list(1:L)
for (l in 1:L)
{
    tree.child.node.lst[[l]]	<- as.list(1:(2^(l - 1)))
    names(tree.child.node.lst[[l]])		<- paste("Level", l, " --  Node", 1:(2^(l - 1)))
    for (node in 1:(2^(l - 1)))
    {
        tree.child.node.lst[[l]][[node]]	<- c(2*(node - 1) + 1, 2*(node - 1) + 2)
    }
}

tree.a.ind.lst <- as.list(1:L)
for (l in 1:L)
{
    tree.a.ind.lst[[l]]	<- as.list(1:(2^(l - 1)))
    names(tree.a.ind.lst[[l]]) <- paste("Level", l, " -- Node", 1:(2^(l - 1)))
    for (node in 1:(2^(l - 1)))
    {
        p.ind <- (1:(2^(L - l))) + (node - 1)*2^(L - l + 1)
        q.ind <- p.ind + 2^(L - l)
        tree.a.ind.lst[[l]][[node]]	<- list(p.ind = p.ind, q.ind = q.ind)
    }
}


#	This a matrix of indices that is used to efficiently sample distributions

prod.ind.mat	<- array(dim = c(L, I))
for (l in 1:L)
{
    p.ind			<- (2^(l - 1)):((2^l) - 1)
    q.ind			<- p.ind + 2^L - 1
    p.q				<- c(rbind(p.ind,q.ind))
    prod.ind.mat[l,] <- rep(p.q, each = 2^(L - l))
}


############################################################################

#	2.	Application of eBayes with Gibbs sampler for hierarchical beta prior

############################################################################

beta.dist.gbbs.sampler <- function(X.mat,beta.initial,y.vec,n.gbbs)
{

    a.vec	<- seq(-24,24,length = I + 1)
    a.med	<- (a.vec[-1] + a.vec[-(I + 1)])/2

    beta.vals	<- seq(-24,24,length = I*20 + 1)
    num.vals <- length(beta.vals)
    beta.vals.a.ind	<- c(rep(1:I, each = 20), I)

    K			<- dim(X.mat)[2]
    n			<- dim(X.mat)[1]
    pi.gbbs		<- array(dim = c(n.gbbs,I))
    beta.gbbs	<- array(dim = c(n.gbbs,K))
    delta.gbbs	<- array(dim = c(n.gbbs,K))
    loglik.tmp	<- rep(NA, num.vals)

    num.y.1		<- sum(y.vec == 1)
    ind.y.1		<- which(rep(y.vec, I + 1) == 1)

    for (g in 1:n.gbbs)
    {
        if (g == 1)	#	Initialize beta.gbbs, pi.gbbs
        {
            beta.gbbs[g,]	<- beta.initial
            pi.gbbs[g,]		<- (c(sum(beta.gbbs[g,] < a.vec[2]),diff(rank(c(a.vec[2:I],100,beta.gbbs[g,]))[1:I] -
                                                                      (1:I))) + 100/I) / 900
            lin.pred.tmp	<- X.mat %*% beta.gbbs[g,]
        }

        #	Gibbs step 1: for k = 1 ... K and l = 1 ... L, sample  delta.k[g] conditionally on  x.vec[k], pi.gbbs


        if (g > 1)	beta.gbbs[g,]	<- beta.gbbs[g - 1,]

        beta.prior	<- pi.gbbs[g,beta.vals.a.ind]
        nvec		<- rep(0,I)

        for (k in 1:K)
        {
            if (k %% 100 == 1)	print(paste(g," -- ", k))

            lin.pred.tmp			<- lin.pred.tmp - X.mat[,k]*beta.gbbs[g,k]
            mult.beta.k.lin.pred	<- rep(lin.pred.tmp, I + 1) + rep(X.mat[,k], I + 1)*rep(a.vec, each = n)

            mult.beta.k.loglik.num	<- rbind(rep(1, num.y.1)) %*% array(mult.beta.k.lin.pred[ind.y.1], dim = c(num.y.1,I + 1))
            mult.beta.k.loglik.denom <- rbind(rep(1, n)) %*% array(log(1 + exp(mult.beta.k.lin.pred)), dim = c(n, I + 1))
            mult.beta.k.loglik <- mult.beta.k.loglik.num - mult.beta.k.loglik.denom
            mult.beta.k.loglik <- mult.beta.k.loglik - max(mult.beta.k.loglik)
            beta.vals.loglik <- approx(a.vec,mult.beta.k.loglik,xout = beta.vals)$y
            beta.pst <- exp(beta.vals.loglik) * beta.prior
            beta.ind <- sample(1:num.vals,1,prob = beta.pst / sum(beta.pst))

            delta.gbbs[g, k] <- beta.vals.a.ind[beta.ind]
            nvec[delta.gbbs[g, k]]	<- nvec[delta.gbbs[g, k]] + 1
            beta.gbbs[g, k] <- beta.vals[beta.ind]
            lin.pred.tmp <- lin.pred.tmp  + X.mat[, k] * beta.gbbs[g, k]
        }



        #	Gibbs step 2: for each node independently sample node prob conditionally on n.vec and then hierarchically construct pi.smp

        if (g < n.gbbs)
        {

            pp.vec		<- rep(NA, 2^L - 1)
            node.ind	<- 1

            for (l in 1:L)
                for (node in 1:(2^(l - 1)))
                {
                    n.p	<- sum(nvec[tree.a.ind.lst[[l]][[node]]$p.ind])
                    n.q	<- sum(nvec[tree.a.ind.lst[[l]][[node]]$q.ind])
                    pp.vec[node.ind] <- rbeta(1, n.p + 1, n.q + 1)
                    node.ind <- node.ind + 1
                }

            pq.vec <- c(pp.vec,1-pp.vec)
            pq.mat <- array(pq.vec[prod.ind.mat], dim=c(L,I))
            pi.gbbs[g + 1, ]		<- apply(pq.mat, 2, prod)
        }
    }
    return(list(a.vec      = a.vec,
                pi.gibbs   = pi.gbbs,
                beta.gibbs = beta.gbbs))
}

n.gibbs <- 300
n <- 200
kappa <- 0.1 #p/n
gamma2 <- 5 #gamma^2
eps <- 1/2 #czesc zmiennych niezerowych

p <- round(kappa*n) #liczba zmiennych
k <- round(eps*p) #liczba zmiennych z niezerowymi wspolczynnikami
gamma <- sqrt(gamma2) #gamma
amplitude <- sqrt(gamma2)/(sqrt(kappa*eps)) #wartosci niezerowych wspolczynnikow

beta.1 <- c(rep(amplitude,times = k), rep(0,times = p - k)) #wektor wspolczynnikow, p-k z nich to 0, a k to 'amplitude'
X <- matrix(rnorm(n*p)/sqrt(n), n, p)
Y		<- rbinom(n,1,inv.log.odds(X %*% beta.1))
fit.1	<- glm(Y ~ X + 0,family = binomial) #change later
res.Yekutieli <- gibbs.beta.Yekutieli(n.gibbs,
                                         Y,
                                         X,
                                         c(-24,24),
                                         L,
                                      fit.1$coefficients)




aa1	<- beta.dist.gbbs.sampler(X, fit.1$coef, Y, n.gibbs)
