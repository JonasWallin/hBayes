#' Eigenvalue and Eigenvector Sampler using MCMC assuming Bingham Distribution
#'
#'
#' @param X Matrix (n x p), where n is the number of samples, and p is the number of dimensions.
#' @param n.mcmc Number of MCMC samples to draw (default: 50,00).
#' @param n.E.mcmc Number of eigenvector samples per MCMC iteration (default: 10).
#' @param L Tree size, used for eigenvalue sampling (default: 6).
#' @param sample.eigenvalues Logical, whether to sample eigenvalues (default: TRUE).
#' @param sample.eigenvector Logical, whether to sample eigenvectors (default: TRUE).
#' @param gamma Optional initialization for eigenvector sampling.
#' @param lambda Optional initialization for eigenvalue sampling.
#' @param verbose Logical, whether to print progress updates (default: FALSE).
#' @return A list with sampled eigenvalues (lambda), eigenvector prior (pi),
#'         sampled eigenvectors (gamma), and additional sampling info.
Eigen.sampler <- function(X,
                          n.mcmc = 5e3,
                          n.E.mcmc = 1e1,
                          L = 6,
                          sample.eigenvalues = TRUE,
                          sample.eigenvector = TRUE,
                          gamma = NULL,
                          lambda = NULL,
                          verbose = FALSE) {
    # Calls the Bingham-based sampler
    return(Eigen.sampler.bingham(X, n.mcmc, n.E.mcmc, L,
                                 sample.eigenvalues = sample.eigenvalues,
                                 sample.eigenvector = sample.eigenvector,
                                 lambda = lambda,
                                 gamma = gamma,
                                 shuffle = FALSE, verbose = verbose))
}
#' Eigenvalue and Eigenvector Sampler using Bingham Distribution
#'
#' This function implements the MCMC-based sampling for eigenvalues and eigenvectors using
#' the Bingham distribution. It is designed to handle large datasets by iterating over
#' the matrix of interest.
#'
#' @param X Matrix (n x p), where n is the number of samples, and p is the number of dimensions.
#' @param n.mcmc Number of MCMC samples to draw (default: 50,000).
#' @param n.E.mcmc Number of eigenvector samples per MCMC iteration (default: 10).
#' @param L Tree size, used for eigenvalue sampling (default: 6).
#' @param sample.eigenvalues Logical, whether to sample eigenvalues (default: TRUE).
#' @param sample.eigenvector Logical, whether to sample eigenvectors (default: TRUE).
#' @param gamma Optional initialization for eigenvector sampling.
#' @param lambda Optional initialization for eigenvalue sampling.
#' @param shuffle Logical, whether to shuffle indices of the eigenvalues (default: FALSE).
#' @param verbose Logical, whether to print progress updates (default: FALSE).
#' @return A list containing MCMC samples of eigenvalues, prior distributions,
#'         sampled eigenvectors, and additional parameters.
Eigen.sampler.bingham <- function(X,
                                  n.mcmc = 5e4,
                                  n.E.mcmc = 1e1,
                                  L = 6,
                                  sample.eigenvalues = TRUE,
                                  sample.eigenvector = TRUE,
                                  gamma = NULL,
                                  lambda = NULL,
                                  shuffle = FALSE,
                                  verbose = FALSE) {

    # Matrix dimensions
    n <- dim(X)[1]
    p <- dim(X)[2]

    # Compute XtX and its eigen decomposition
    XtX <- t(X) %*% X
    E <- eigen(XtX / n)  # Eigen decomposition of XtX/n

    # Cholesky decomposition of XtX
    R.XtX <- chol(XtX)

    # Build tree structure for eigenvalue sampling
    tree.data <- build.tree(L)
    I <- tree.data$I

    # Initialize MCMC storage arrays
    mcmc.lambda <- array(dim = c(n.mcmc, p))  # Eigenvalues
    mcmc.pi <- array(dim = c(n.mcmc, I))      # Prior for eigenvalues
    mcmc.gamma <- array(dim = c(p, p, n.mcmc))  # Eigenvectors

    # Initial values for eigenvalues and eigenvectors
    init_params_Lambda <- initLambdaSampling(E, n, I)

    if (is.null(lambda)) {
        mcmc.lambda[1, ] <- E$values
    } else {
        mcmc.lambda[1, ] <- lambda
    }

    mcmc.pi[1, ] <- rep(1 / I, I)

    if (!is.null(gamma)) {
        mcmc.gamma[, , 1] <- t(E$vectors) %*% gamma
    } else {
        mcmc.gamma[, , 1] <- diag(p)
    }

    Lambda_hatE <- list(values = 1 / E$values,
                        vectors = diag(nrow = length(E$values), ncol = length(E$values)))
    XXtE <- list(values = -E$values * n / 2,
                 vectors = diag(nrow = length(E$values), ncol = length(E$values)))

    res.E <- NULL
    for (i in 2:n.mcmc) {
        if (verbose && i %% 500 == 0) {
            cat('sample =', i, 'of', n.mcmc, '\n')
        }

        # Sample eigenvalues if specified
        if (sample.eigenvalues) {
            lambda_hat <- colSums(diag(n * E$values) %*% mcmc.gamma[, , i - 1]^2)
            mcmc.lambda[i, ] <- sampleLambda(lambda_hat, n, mcmc.pi[i - 1, ], init_params_Lambda)
        } else {
            mcmc.lambda[i, ] <- mcmc.lambda[i - 1, ]
        }

        # Sample prior for eigenvalues
        lambda.ind <- findInterval(log(mcmc.lambda[i, ]), init_params_Lambda$a_vec, left.open = TRUE)
        nvec <- table(factor(lambda.ind, levels = 1:I))
        mcmc.pi[i, ] <- exp(polya.gibs(nvec, tree.data))  # Generate eigenvalue prior distribution

        # Sample eigenvectors if specified
        ind.lambda <- order(mcmc.lambda[i, ], decreasing = TRUE)
        mcmc.lambda[i, ] <- mcmc.lambda[i, ind.lambda]

        if (sample.eigenvector) {
            Lambda_hatE$values <- 1 / mcmc.lambda[i, ]
            res.E <- rbing.iter(n.E.mcmc, eigA = XXtE, eigB = Lambda_hatE,
                                E0 = mcmc.gamma[, ind.lambda, i - 1], ret = 2)
            mcmc.gamma[, , i] <- res.E$Es[, , n.E.mcmc]
        } else {
            mcmc.gamma[, , i] <- mcmc.gamma[, ind.lambda, i - 1]
        }
    }

    # Transform the sampled eigenvectors back to the correct scale
    res_Egamma <- apply(mcmc.gamma, 3, function(x) E$vectors %*% x)
    mcmc.gamma <- array(res_Egamma, dim = dim(mcmc.gamma))

    return(list(lambda = mcmc.lambda,
                pi = mcmc.pi,
                gamma = mcmc.gamma,
                gamma.obj = res.E,
                lambda.param = init_params_Lambda))
}



#' Eigenvalue and Eigenvector Sampler using Bingham Distribution usng Gibbs sampling
#'
#'
#' @param X Matrix (n x p), where n is the number of samples, and p is the number of dimensions.
#' @param n.mcmc Number of MCMC samples to draw (default: 50,00).
#' @param n.E.mcmc Number of eigenvector samples per MCMC iteration (default: 10).
#' @param L Tree size, used for eigenvalue sampling (default: 6).
#' @param sample.eigenvalues Logical, whether to sample eigenvalues (default: TRUE).
#' @param sample.eigenvector Logical, whether to sample eigenvectors (default: TRUE).
#' @param gamma Optional initialization for eigenvector sampling.
#' @param lambda Optional initialization for eigenvalue sampling.
#' @param shuffle Logical, whether to shuffle indices of the eigenvalues when sampling (default: FALSE).
#' @param verbose Logical, whether to print progress updates (default: FALSE).
#' @return A list containing MCMC samples of eigenvalues, prior distributions,
#'         sampled eigenvectors, and additional parameters.
Eigen.sampler.bingham <- function(X,
                                  n.mcmc=5e3,
                                  n.E.mcmc = 1e1,
                                  L = 6,
                                  sample.eigenvalues = T,
                                  sample.eigenvector=T,
                                  gamma = NULL,
                                  lambda = NULL,
                                  shuffle=FALSE,
                                  verbose=FALSE){

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

    if(is.null(lambda)){
        mcmc.lambda[1,] = E$values
    }else{
        mcmc.lambda[1,] = lambda
    }

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
        if(verbose){
            if(i%%500==0)
                cat('sample = ',i,' of ',n.mcmc,'\n')
        }
        if(sample.eigenvalues){
        # Example lambda_hat calculation (to be done per iteration in your actual Gibbs sampling)
            lambda_hat <- colSums(diag(n*E$values)%*%mcmc.gamma[,,i-1]^2)
            # Sampling lambda in the Gibbs sampler
            # \pi(lambda|X) \propto lambda^{-n/2} exp( - \sum lambda_hat_i/2 * lambda_i^{-1} )*pi(\lambda)
            mcmc.lambda[i, ] <- sampleLambda(lambda_hat, n, mcmc.pi[i-1, ],init_params_Lambda)
        }else{
            mcmc.lambda[i, ] <- mcmc.lambda[i-1, ]
        }

        #E_ <-  E$vectors%*%mcmc.gamma[,,i-1]
        #E_n <- mcmc.gamma[,,i-1]
        #sample prior
        lambda.ind <- findInterval(log(mcmc.lambda[i,]), init_params_Lambda$a_vec, left.open = T)
        nvec <- table(factor(lambda.ind, levels = 1:I))
        mcmc.pi[i, ] <- exp(polya.gibs(nvec, tree.data))     # Generate eigenvalue prior distribution

        # sample eigenvectors
        # exp(-0.5 \Sigma E D E^T) = exp(tr(-0.5  (E^T E_S) D_S (E^T E_S)^T D ))
        ind.lambda <- order(mcmc.lambda[i, ], decreasing = T)
        mcmc.lambda[i, ] <- mcmc.lambda[i,ind.lambda]
        if(sample.eigenvector){
            Lambda_hatE$values <- 1/mcmc.lambda[i, ]
            if(shuffle==FALSE){
                res.E <- rbing.iter(n.E.mcmc,
                               A = NULL,
                               B = NULL,
                               eigA=XXtE,
                               eigB=Lambda_hatE,
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



#' Eigenvalue and Eigenvector Sampler (First Version) Using Householder Decomposition
#'
#' This function performs MCMC sampling of eigenvalues and eigenvectors of a matrix using Householder decomposition.
#' It is designed to handle large datasets and complex eigenstructure by iterating over the matrix of interest.
#' The function samples both eigenvalues and eigenvectors using a Metropolis-Hastings algorithm and computes
#' the posterior distributions over these quantities.
#'
#' @param X Matrix (n x p), where n is the number of samples and p is the number of dimensions.
#'          This matrix is assumed to be the data from which eigenvalues and eigenvectors are sampled.
#' @param n.mcmc Number of MCMC samples to generate (default: 50,000).
#' @param n.E.mcmc Number of eigenvector samples to generate in each MCMC loop (default: 10).
#' @param L Integer representing the tree size used for the eigenvalue prior sampling (default: 6).
#' @param sample.eigenvalues Logical flag to indicate if eigenvalues should be sampled (default: TRUE).
#' @param sample.eigenvector Logical flag to indicate if eigenvectors should be sampled (default: TRUE).
#' @param lambda Optional initialization for eigenvalue sampling. If NULL, the eigenvalues of X'X/n are used.
#' @param gamma Optional initialization for eigenvector sampling. If NULL, the eigenvectors of X'X/n are used.
#' @param verbose Logical flag to control printing of progress updates during MCMC sampling (default: FALSE).
#' @return A list containing:
#'   - lambda: MCMC samples of eigenvalues (matrix of size n.mcmc x p),
#'   - pi: MCMC samples of prior distribution over eigenvalues (matrix of size n.mcmc x I),
#'   - gamma: MCMC samples of eigenvectors (array of size p x p x n.mcmc),
#'   - gamma.obj: Results from the Metropolis-Hastings sampling of eigenvectors,
#'   - lambda.param: Parameters used for lambda sampling initialization.

#'
Eigen.sampler.old <- function(X,
                          n.mcmc=5e4,
                          n.E.mcmc = 1e1,
                          L = 6,
                          sample.eigenvalues = T,
                          sample.eigenvector=T,
                          lambda = NULL,
                          gamma = NULL,
                          verbose=FALSE){

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
    if(is.null(lambda)){
    mcmc.lambda[1,] = E$values
    }else{
        mcmc.lambda[1,] = lambda
    }

    mcmc.pi[1,] = rep(1/I,I)
    if(is.null(gamma)==F){
        mcmc.gamma[,,1] = gamma
    }else{
        mcmc.gamma[,,1] = E$vectors
    }
    mcmc.loglik[1]  <- oracle.gamma.R.loglik(R.XtX,  mcmc.gamma[,,1],mcmc.lambda[1, ], n)
    if(sample.eigenvector){
        res.E <- oracle.metrop.sampler(E$vectors,X, E$values, n.E.mcmc)
    }else{
        res.E <- NULL
    }
    ind.lambda <- 1:p
    for(i in 2:n.mcmc){
        if(verbose){
            if(i%%500==0)
                cat('sample = ',i,' of ',n.mcmc,'\n')
        }
        if(sample.eigenvalues){
            # Example lambda_hat calculation (to be done per iteration in your actual Gibbs sampling)
            lambda_hat <- colSums((R.XtX%*% mcmc.gamma[,,i-1])^2)

            # Sampling lambda in the Gibbs sampler
            # \pi(lambda|X) \propto lambda^{-n/2} exp( - \sum lambda_hat_i/2 * lambda_i^{-1} )*pi(\lambda)
            mcmc.lambda[i, ] <- sampleLambda(lambda_hat, n,  mcmc.pi[i-1, ],init_params_Lambda)
        }else{
            mcmc.lambda[i, ] <- mcmc.lambda[i-1, ]
        }
        #sample prior
        lambda.ind <- findInterval(log(mcmc.lambda[i,]), init_params_Lambda$a_vec, left.open = T)
        nvec <- table(factor(lambda.ind, levels = 1:I))
        mcmc.pi[i, ] <- exp(polya.gibs(nvec, tree.data))     # Generate eigenvalue prior distribution

        # sample eigenvectors
        # exp(-0.5 \Sigma E D E^T) = exp(tr(-0.5  (E^T E_S) D_S (E^T E_S)^T D ))
        ind.lambda <- order(mcmc.lambda[i, ], decreasing = T)
        mcmc.lambda[i, ] <- mcmc.lambda[i,ind.lambda]
        if(sample.eigenvector){
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
           gamma.obj =res.E,
           lambda.param = init_params_Lambda
           ))
}
