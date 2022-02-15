
###
# Polya tree for resturant sampling
#
#
###



build.tree <- function(L){

    I	<- 2^L

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
        q.ind			<- p.ind + I - 1
        p.q				<- c(rbind(p.ind,q.ind))
        prod.ind.mat[l,] <- rep(p.q, each = 2^(L - l))
    }
    return(list( L = L,
                 I = I,
                 tree.a.ind.lst =  tree.a.ind.lst,
                 prod.ind.mat = prod.ind.mat
    ))
}

#'
#' Gibbs sampling for polya data
#' generates pi (log)
#' @param nvec  - (vec) number of elements in each bin
#' @param tree.data list
#'        contains relevant data need for sampling pi
#'
polya.gibs <- function(nvec, tree.data){


    L <- tree.data$L
    I <- tree.data$I
    pp.vec		<- rep(NA, I - 1)
    node.ind	<- 1
    for (l in 1:L){
        for (node in 1:(2^(l - 1)))
        {
            n.p	<- sum(nvec[tree.data$tree.a.ind.lst[[l]][[node]]$p.ind])
            n.q	<- sum(nvec[tree.data$tree.a.ind.lst[[l]][[node]]$q.ind])
            pp.vec[node.ind] <- rbeta(1, n.p + 1, n.q + 1)
            node.ind <- node.ind + 1
        }
    }
    pq.vec <- c(pp.vec,1-pp.vec)
    pq.mat <- array(log(pq.vec[tree.data$prod.ind.mat]), dim=c(L,I))
    pi.log		<- apply(pq.mat, 2, sum)

    return(pi.log)
}
