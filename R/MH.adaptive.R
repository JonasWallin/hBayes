

MH.init <- function(sigma, accept.target = 0.234){

    adaptive.obj <- list(sigma = sigma, # sampling prob
                         accept.target = accept.target, #target acceptance rate
                         adap.num = 50, #how often to adapt
                         accept.count = 0, #count acceptence
                         count = 0, #how many samples
                         c_p = 0.5, # scaling factor
                         n.SA  = 1 #how many updates step
                         )
    return(adaptive.obj)
}

MH.adaptive <- function(acc, adaptive.obj){

    if(adaptive.obj$count %% adaptive.obj$adap.num  == 0){
        acc          <-  adaptive.obj$accept.count / adaptive.obj$adap.num
        w0 <- (1/adaptive.obj$n.SA)^adaptive.obj$c_p
        adaptive.obj$sigma <- adaptive.obj$sigma * exp(w0 *  (acc - adaptive.obj$accept.target))
        adaptive.obj$n.SA <- adaptive.obj$n.SA + 1
        adaptive.obj$accept.count <- 0
        adaptive.obj$count <- 1
    }else{
        adaptive.obj$accept.count <- adaptive.obj$accept.count + acc
        adaptive.obj$count <- adaptive.obj$count + 1
    }
    return(adaptive.obj)
}
