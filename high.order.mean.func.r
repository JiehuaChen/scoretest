###Davies Formula for the case where the mean function is unknown
### with high order correction
### Jiehua Chen Nov. 26th.
### For exponential family
### P(Z > b)
###calculating $xi$, when n is large, xi should be close to zero.

davis.third <- function(rho, dist, b, X, k){
    n <- dim(dist)[1]
    #the hat matrix
    H <- diag(1, n) - X%*% solve(t(X) %*% X) %*% t(X)
    V <- rho^(dist)
    V.H <- V %*% H
    H.V.H <- H%*% V %*% H
    W.2 <- sum(diag(H.V.H)^2)
    V.H.2 <- V.H %*% V.H
    tr.V.H.2 <- sum(diag(V.H.2))
    G.denom <- sqrt(2*tr.V.H.2 + k*W.2)
    G <- H.V.H/G.denom
    tr.G.3 <- sum(diag(G %*% G %*% G))
    xi <- (-1 +  sqrt(1 + 16 * b * tr.G.3))/(8 *tr.G.3)
    var.z.xi <- 1 + 2 *(b - xi)/xi
    exp.phi <- exp(-2/3 * xi * b + 1/6 * xi^2)
    p.z.xi <- (var.z.xi)^(-1/2)* exp.phi
    return(c(p.z.xi, xi, var.z.xi, exp.phi))
}


### var(diff(Z))

var.diffz <- function(rho, dist, b, X, k){
    n <- dim(dist)[1]
    #the hat matrix
    H <- diag(1, n) - X%*% solve(t(X) %*% X) %*% t(X)
    V <- rho^(dist)
    V.H <- V %*% H
    H.V.H <- H%*% V %*% H
    W.2 <- sum(diag(H.V.H)^2)
    V.H.2 <- V.H %*% V.H
    tr.V.H.2 <- sum(diag(V.H.2))
    G.denom <- sqrt(2*tr.V.H.2 + k*W.2)
    G <- H.V.H/G.denom


    # calculate G.derivative
    V.diff <- dist*rho^(dist-1)
    H.V.diff.H <- H %*% V.diff %*% H
    tr.V.diff.H.V.H <- sum(diag(V.diff%*%H %*% V.H))
    sum.W.diff <- sum(diag(H.V.H)*diag(H.V.diff.H))

    G.diff.1 <- H.V.diff.H/G.denom
    G.diff.2 <- H.V.H*(2*tr.V.diff.H.V.H + k*sum.W.diff)/G.denom^(3)
    G.diff <- G.diff.1 - G.diff.2
    ##first order
    var.diff.1 <- 2 * sum(diag(G.diff %*% G.diff))


    var.diff.2 <- sum(diag(G%*%G.diff%*% G.diff))
    ##calculate xi
    tr.G.3 <- sum(diag(G %*% G %*% G))
    xi <- (-1 +  sqrt(1 + 16 * b * tr.G.3))/(8 *tr.G.3)
    var.diff <- var.diff.1 +  8 * xi * var.diff.2
    var.diff.sr <- sqrt(var.diff)
    var.diff.no <- sqrt(var.diff.1)
    return(c(var.diff.no, var.diff.sr, var.diff.2, xi))
}



## Equation (12)
davis.prob.high.order <- function(l, u, b, dist, X,k){
    f <- 1/2*( davis.third(u, dist, b, X,k)[[1]] * var.diffz(u, dist, b, X,k)[[2]] + davis.third(l, dist, b,X,k)[[1]] * var.diffz(l, dist, b, X,k)[[2]])
    for ( i in 1:99){
        f <- f +  davis.third((l+(u -l)/100 * i), dist, b, X,k)[[1]] * var.diffz((l+(u -l)/100 * i), dist, b, X,k)[[2]] 
    }
    p.l <- davis.third(l, dist, b, X, k)
    var.z.l <- p.l[[3]]
    xi.l <- p.l[[2]]
    ##P(Z(L) > b)
    p.l.z <- p.l[[4]]*(1-pnorm(0, mean=-var.z.l*xi.l, sd=sqrt(var.z.l)))*exp(var.z.l *xi.l^2/2)
    f.2 <- p.l.z+ f * (u-l)/100 *1/(2*pi)
    return(c(p.l.z, f.2)) 
}  

