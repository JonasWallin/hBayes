#' @param x.sim simulations to get posterior conf int
plot.graph <- function(x,
                       window.size,
                       markers,
                       beta_true= NULL,
                       name='',
                       x.sim = NULL){
    x.window <- rep(0,length(x))
    beta_true.window = rep(0, length(x))
    marker.tot <- c(0, cumsum(markers))
    CI <- matrix(0, ncol=2, nrow=length(x))
    for(i in 2:length(marker.tot)){
        x.window[(1+marker.tot[i-1]):marker.tot[i]] <- zoo::rollsum(x[(1+marker.tot[i-1]):marker.tot[i]],
                                                                    window.size,
                                                                    na.pad=T,
                                                                    fill="extend")

        if(!is.null(beta_true))
            beta_true.window[(1+marker.tot[i-1]):marker.tot[i]] <- zoo::rollsum(beta_true[(1+marker.tot[i-1]):marker.tot[i]],
                                                                                window.size,
                                                                                na.pad=T,
                                                                                fill="extend")

        if(!is.null(x.sim)){
            CI.window = t(apply(x.sim[,(1+marker.tot[i-1]):marker.tot[i]], 1,zoo::rollsum,k= window.size,
                                na.pad=T,
                                fill="extend"))
            CI[(1+marker.tot[i-1]):marker.tot[i],] = t(apply(CI.window, 2, quantile, probs=c(0.05,0.95)))
        }
    }
    df <- data.frame(locus = 1:length(x.window), beta = x.window)
    fig <- ggplot(df)
    if(!is.null(beta_true))

        fig <- fig + geom_line(data=data.frame(locus = 1:length(x.window), beta = beta_true.window),
                               aes(x=locus, y= beta),
                               colour ="red",
                               alpha = 0.5)
    fig <- fig + geom_line(aes(x=locus, y=beta))
    fig <- fig +  geom_vline(xintercept = marker.tot,
                             colour ="blue",
                             linetype = "dashed",
                             alpha = 0.5)
    fig <- fig + ylab(TeX("$\\beta$"))+ theme(axis.text.x  = element_text(size=12),
                                              axis.text.y  = element_text(size=12),
                                              axis.title.y = element_text(angle=0,vjust=0.5,size=22,face="bold"),
                                              axis.title.x = element_text(size=20,face="bold"))
    if(!is.null(x.sim)){
        fig <- fig  + geom_ribbon(data = data.frame(locus=1:length(x),
                                                    low = CI[,1]
                                                    ,upp = CI[,2]),
                                  aes(x=locus,ymin = low, ymax=upp), alpha=0.2, fill='blue')
    }
    fig <- fig +  labs(title=paste(name))
    fig <- fig + coord_cartesian(ylim=c(-0.3, 0.1))
    #fig <- fig +  labs(title=paste(name,', RMSE = ',round(sqrt(mean((beta_true.window-x.window)^2))/sd(beta_true.window),3) , ", window size = ", window.size,sep=""))
    fig <- fig+ theme(plot.title = element_text(hjust = 0.5, size=20))

    return(list(fig=fig, res = beta_true.window-x.window, beta_true = beta_true.window ) )
}

smooth.beta <- function(x, beta_true, markers, window.size=10){
    x.window <- x
    beta.smooth <- beta_true
    marker.tot <- c(0, cumsum(markers))
    for(i in 2:length(marker.tot)){
        x.window[(1+marker.tot[i-1]):marker.tot[i]] <- zoo::rollsum(x[(1+marker.tot[i-1]):marker.tot[i]],
                                                                    window.size,
                                                                    na.pad=T,
                                                                    fill="extend")
        beta.smooth[(1+marker.tot[i-1]):marker.tot[i]] <- zoo::rollsum(beta_true[(1+marker.tot[i-1]):marker.tot[i]],
                                                                       window.size,
                                                                       na.pad=T,
                                                                       fill="extend")
    }
    return(list(x=x.window, RMSE = sqrt(mean( (x.window-beta.smooth)^2 ))/sd(beta.smooth)))
}

scores <- function(XB, X, beta, intercept, beta.true, Xbeta=NULL){
    if(is.null(Xbeta)){
        MSE_XB <- sqrt(mean((XB- X%*%beta - intercept)^2))/sd(XB)
    }else{
        MSE_XB <- sqrt(mean((XB- Xbeta)^2))/sd(XB)
    }
    MSE_beta <- sqrt(mean((beta.true  -beta )^2))/sd(beta.true)
    res.smooth <- smooth.beta(beta, beta.true, markers, window.size=10)

    return(c(MSE_XB,
             MSE_beta,
             res.smooth$RMSE,
             beta[qtl.pos],
             res.smooth$x[qtl.pos]
             ))
}

