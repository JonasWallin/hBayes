#' eigenvalue, vector sampler
#' @param X -       (n x p)
#' @param n.mcmc.   (int) number of mcmc samples
#' @param n.E.mcmc. (int) number of Eigenvector samples in each loop
#' @param L.        (int) tree size
Eigen.sampler <- function(X,
                                  n.mcmc=5e4,
                                  n.E.mcmc = 1e1,
                                  L = 6,
                          sample.eigenvalues=T,
                          gamma=NULL){

    return(Eigen.sampler.bingham(X, n.mcmc, n.E.mcmc, L,
                                    sample.eigenvalues = sample.eigenvalues,
                                    gamma = gamma,
                                 shuffle=FALSE))
}



#' eigensampler using bingham
#' @param X -       (n x p)
#' @param n.mcmc.   (int) number of mcmc samples
#' @param n.E.mcmc. (int) number of Eigenvector samples in each loop
#' @param L.        (int) tree size
Eigen.sampler.bingham <- function(X,
                                  n.mcmc=5e4,
                                  n.E.mcmc = 1e1,
                                  L = 6,
                                  sample.eigenvalues = T,
                                  gamma = NULL,
                                  shuffle=FALSE){

    #fixing data
    n <- dim(X)[1]
    p <- dim(X)[2]
    XtX <- t(X) %*% X
    E <-  eigen(XtX/n)

    R.XtX <- chol(XtX)


    tree.data <- build.tree(L)
    I = tree.data$I

    #mcmc values
    mcmc.lambda <- array(dim = c(n.mcmc, p))
    mcmc.pi <- array(dim = c(n.mcmc, I))
    mcmc.gamma <- array(dim = c(p, p, n.mcmc))



    # Initialization
    init_params_Lambda <- initLambdaSampling(E, n, I)

    mcmc.lambda[1,] = E$values
    mcmc.pi[1,] = rep(1/I,I)
    if(is.null(gamma)==F){
        mcmc.gamma[,,1] =  t(E$vectors)%*%gamma
    }else{
        mcmc.gamma[,,1] = diag(nrow=p,ncol=p)
    }


    Lambda_hatE <- list(values = 1/E$values,
                        vectors = diag(nrow=length(E$values),ncol=length(E$values)))
    XXtE <- list(values = -E$values*n/2,
                 vectors = diag(nrow=length(E$values),ncol=length(E$values)))

    res.E <- NULL
    for(i in 2:n.mcmc){
        # Example lambda_hat calculation (to be done per iteration in your actual Gibbs sampling)
        lambda_hat <- rowSums(mcmc.gamma[,,i-1]^2%*%diag(n*E$values))

        # Sampling lambda in the Gibbs sampler
        # \pi(lambda|X) \propto lambda^{-n/2} exp( - \sum lambda_hat_i/2 * lambda_i^{-1} )*pi(\lambda)
        mcmc.lambda[i, ] <- sampleLambda(lambda_hat, n, mcmc.pi[i-1, ],init_params_Lambda)

        #sample prior
        lambda.ind <- findInterval(log(mcmc.lambda[i,]), init_params_Lambda$a_vec, left.open = T)
        nvec <- table(factor(lambda.ind, levels = 1:I))
        mcmc.pi[i, ] <- exp(polya.gibs(nvec, tree.data))     # Generate eigenvalue prior distribution

        # sample eigenvectors
        # exp(-0.5 \Sigma E D E^T) = exp(tr(-0.5  (E^T E_S) D_S (E^T E_S)^T D ))
        ind.lambda <- order(mcmc.lambda[i, ], decreasing = T)
        mcmc.lambda[i, ] <- mcmc.lambda[i,ind.lambda]
        if(sample.eigenvalues){
            Lambda_hatE$values <- 1/mcmc.lambda[i, ]
            if(shuffle==FALSE){
                res.E <- rbing.iter(n.E.mcmc,
                               A = NULL,
                               B = NULL,
                               eigA=Lambda_hatE,
                               eigB=XXtE,
                               E0=mcmc.gamma[,ind.lambda,i-1],
                               EtAE=NULL,
                               ret = 2)
            }else{
                res.E <- rbing(n.E.mcmc,
                                    A = NULL,
                                    B = NULL,
                                    eigA=Lambda_hatE,
                                    eigB=XXtE,
                                    E0=mcmc.gamma[,ind.lambda,i-1],
                                    EtAE=NULL,
                                    ret = 2)

            }
            mcmc.gamma[,,i] <- res.E$Es[,,n.E.mcmc]
        }else{
            mcmc.gamma[,,i] <- mcmc.gamma[,ind.lambda,i-1]
        }
    }
    res_Egamma <- apply(mcmc.gamma,3,function(x) E$vectors%*%x) #transform to correct scale
    mcmc.gamma <- array(res_Egamma, dim = dim(mcmc.gamma))
    return(list(lambda = mcmc.lambda,
                pi = mcmc.pi,
                gamma = mcmc.gamma,
                gamma.obj =res.E,
                lambda.param = init_params_Lambda
    ))
}



#' first version of the eigensampler
#' @param X -       (n x p)
#' @param n.mcmc.   (int) number of mcmc samples
#' @param n.E.mcmc. (int) number of Eigenvector samples in each loop
#' @param L.        (int) tree size
Eigen.sampler.old <- function(X,
                          n.mcmc=5e4,
                          n.E.mcmc = 1e1,
                          L = 6,
                          sample.eigenvalues = T,
                          gamma = NULL){

    #fixing data
    n <- dim(X)[1]
    p <- dim(X)[2]
    XtX <- t(X) %*% X
    E <-  eigen(XtX/n)

    R.XtX <- chol(XtX)


    tree.data <- build.tree(L)
    I = tree.data$I

    #mcmc values
    mcmc.lambda <- array(dim = c(n.mcmc, p))
    mcmc.pi <- array(dim = c(n.mcmc, I))
    mcmc.gamma <- array(dim = c(p, p, n.mcmc))
    mcmc.loglik <- rep(NA, n.mcmc)



    # Initialization
    init_params_Lambda <- initLambdaSampling(E, n, I)

    mcmc.lambda[1,] = E$values
    mcmc.pi[1,] = rep(1/I,I)
    if(is.null(gamma)==F){
        mcmc.gamma[,,1] = gamma
    }else{
        mcmc.gamma[,,1] = E$vectors
    }
    mcmc.loglik[1]  <- oracle.gamma.R.loglik(R.XtX,  mcmc.gamma[,,1],mcmc.lambda[1, ], n)
    if(sample.eigenvalues){
        res.E <- oracle.metrop.sampler(E$vectors,X, E$values, n.E.mcmc)
    }else{
        res.E <- NULL
    }

    for(i in 2:n.mcmc){
        # Example lambda_hat calculation (to be done per iteration in your actual Gibbs sampling)
        lambda_hat <- colSums((R.XtX%*% mcmc.gamma[,,i-1])^2)

        # Sampling lambda in the Gibbs sampler
        # \pi(lambda|X) \propto lambda^{-n/2} exp( - \sum lambda_hat_i/2 * lambda_i^{-1} )*pi(\lambda)
        mcmc.lambda[i, ] <- sampleLambda(lambda_hat, n,  mcmc.pi[i-1, ],init_params_Lambda)

        #sample prior
        lambda.ind <- findInterval(log(mcmc.lambda[i,]), init_params_Lambda$a_vec, left.open = T)
        nvec <- table(factor(lambda.ind, levels = 1:I))
        mcmc.pi[i, ] <- exp(polya.gibs(nvec, tree.data))     # Generate eigenvalue prior distribution

        # sample eigenvectors
        # exp(-0.5 \Sigma E D E^T) = exp(tr(-0.5  (E^T E_S) D_S (E^T E_S)^T D ))
        ind.lambda <- order(mcmc.lambda[i, ], decreasing = T)
        mcmc.lambda[i, ] <- mcmc.lambda[i,ind.lambda]
        if(sample.eigenvalues){
            res.E <- oracle.metrop.sampler(mcmc.gamma[,ind.lambda,i-1],
                                           X,
                                           mcmc.lambda[i, ],
                                           n.E.mcmc,
                                           target.acc = 0.23,
                                           MH.objs = res.E$MH.objs,
                                           R = R.XtX)
            mcmc.gamma[,,i] <- res.E$Es[,,n.E.mcmc]
        }else{
            mcmc.gamma[,,i] <- mcmc.gamma[,ind.lambda,i-1]
        }
    }

    return(list(lambda = mcmc.lambda,
           pi = mcmc.pi,
           gamma = mcmc.gamma,
           gamma.obj =res.E
           ))
}
