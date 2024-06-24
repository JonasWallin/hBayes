test_that("test 2x2 bingham for diagonal", {
    set.seed(1)
    sim <- 10000
    b <- c(2,1)
    a <- c(2,1)
    E <- array(dim=c(dim = c(2,2,sim)))
    for(i in 1:sim){
        E[,,i] <-   rbing.2diagmatrix(-a, b)
    }

    V1 <- round(apply(E[1,,], 1, var),2)
    expect_equal(V1, c(0.38,0.62)) #taking from bingham sampling in rstiefel

    b <- c(2,.1)
    a <- c(.1,1)
    E <- array(dim=c(dim = c(2,2,sim)))
    for(i in 1:sim)
        E[,,i] <-   rbing.2diagmatrix(-a, b)

    V1 <- round(apply(E[1,,], 1, var),2)
    expect_equal(V1, c(0.7,0.3)) #taking from bingham sampling in rstiefel

})

test_that("test 2x2 bingham for non-diagonal", {
    set.seed(2)

    sim <- 10000
    b <- c(2,1)
    a <- c(2,1)
    E <- rbing.2diagmatrix(a, b)
    A <- E%*%diag(a)%*%t(E)
    E <- array(dim=c(dim = c(2,2,sim)))
    for(i in 1:sim)
        E[,,i] <-   rbing.2diagmatrix(-A, b)

    V1 <- round(apply(E[1,,], 1, var),2)
    expect_equal(V1, c(0.26,0.74)) #taking from bingham sampling in rstiefel
})

test_that("test 2x2 eigenvalue ", {
    set.seed(1)

    for(i in 1:sim){
        A <- matrix(rnorm(4),2,2)
        A <- t(A)%*%A
        EA <- eigen_2x2(A)
        EAm <- eigen_2x2(-A)
        expect_equal(A,EA$vectors%*%diag(EA$values)%*%t(EA$vectors))
        expect_equal(-A,EAm$vectors%*%diag(EAm$values)%*%t(EAm$vectors))

        EA_true <- eigen(A)
        expect_equal(EA$values,EA_true$values)

        EA$vectors <- EA$vectors*(EA_true$vectors[1,1]/EA$vectors[1,1])
        expect_equal(EA$vectors,EA_true$vectors)

    }
})
