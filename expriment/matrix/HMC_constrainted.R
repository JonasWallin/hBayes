library(tmg)
# Set number of samples
n=15000;
p = 4^2
#Define precision matrix and linear term
M = solve(toeplitz(c(2,1,rep(0,p-2))))
r = rep(0,p)

# Set initial point for the Markov chain
initial = rep(0,p)



constr = list()
for(i in 1:sqrt(p)){
    A1 = matrix(0,p,p)
    diag(A1) <- c(rep(0,sqrt(p)*(i-1)),rep(1,sqrt(p)),rep(0,sqrt(p)*(sqrt(p)-i)))
    B1 = rep(0,p)
    C1 = 1
    constr[[i]] = list(A1,B1,C1)
}



# Sample and plot
samples = rtmg(n, M, r, initial, f=NULL, g =NULL,q =  constr);
plot(samples[,1], pch=".")
