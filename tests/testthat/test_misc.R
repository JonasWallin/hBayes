test_that("C++ and R implementations match", {
    # Set up some example data
    set.seed(123)
    Gamma <- matrix(rnorm(16), nrow = 4)
    v_p <- c(1, 2, 3)

    # Run the C++ implementation
    result_cpp <- house_reflection_mult_cpp(Gamma, v_p, right = TRUE)

    # Run the R implementation
    result_R <- house.reflctn.mult(Gamma, v_p, right = TRUE)

    # Run the C++ implementation
    result_cpp_left <- house_reflection_mult_cpp(Gamma, v_p, right = FALSE)

    # Run the R implementation
    result_R_left <- house.reflctn.mult(Gamma, v_p, right = FALSE)

    house_reflection_mult_inplace_cpp(Gamma, v_p, right = FALSE)
    # Compare the results
    expect_equal(result_cpp, result_R)
    expect_equal(result_cpp_left, result_R_left)
    expect_equal(Gamma, result_R_left)
})
test_that("Compare likelihood estimation", {
    library(mvtnorm)
    set.seed(123)
    p <- 4
    n <- 10
    #Gamma <- matrix(rnorm(p^2), nrow = p)
    S <- solve(toeplitz(c(2,1,rep(0,p-2))))
    X	<- rmvnorm(n,rep(0,p),S)
    E <-  eigen(S)
    Gamma <- E$vectors

    XtX <- t(X)%*%X
    Sigma <- (Gamma)%*%diag(E$values)%*%t(Gamma)
    R<-chol(XtX)
    #V <- R%*%Gamma #p^2 (complexity)
    #t(v)%*%v - t(Gamma)%*%XtX%*%Gamma
    #colSums(V^2) - diag(t(Gamma)%*%XtX%*%Gamma)
    #sum(diag(XtX%*%solve(Sigma)) - colSums(V^2)/E$values)
    #-0.5*sum(colSums(V^2)/E$values) - 0.5 *n* sum(log(E$values)) - 0.5*n*p*log(2*pi)
    oracle.gamma.x.loglik.new(R, Gamma, E$values, n)


    expect_equal(sum(dmvnorm(X,sigma = Sigma,log=T)), oracle.gamma.x.loglik(X, (Gamma), E$values))
    expect_equal(oracle.gamma.x.loglik.new(R, Gamma, E$values, n), oracle.gamma.x.loglik(X, (Gamma), E$values))

})

