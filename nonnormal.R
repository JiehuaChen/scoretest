#nonnormal 

#calculating xi for nonnormal distribuiton
xi.cal <- function(m6,m3, m2, m4, rho, dist, b, X){
    #m6 is the 6th moment, # m3 is the 3rd moment # k is the kurtosis
    n <- dim(dist)[1]
    #the hat matrix
    H <- diag(1, n) - X%*% solve(t(X) %*% X) %*% t(X)
    V <- rho^(dist)
    V.H <- V %*% H
    H.V.H <- H%*% V %*% H
    W.2 <- sum(diag(H.V.H)^2)
    V.H.2 <- V.H %*% V.H
    tr.V.H.2 <- sum(diag(V.H.2))
    G.denom <- m2*sqrt(2*tr.V.H.2 + k*W.2)
    G <- H.V.H/G.denom
    mem3 <-8*m2^3*sum(diag(G%*%G%*%G))+sum(diag(G)^3)*(m6- 15*m2^3-3*k*m2^3 - 10*m3^2)+(4*sum(G^3) + 6*sum(diag(diag(G))%*%G%*%diag(diag(G))))*m3^2 
    mem3 <- mem3 + 12*k*m2^3*(sum( G^2 %*% diag(diag(G)))- sum(diag(G)^3))
    xi <- (-1 +  sqrt(1 + 2* b * mem3))/(mem3)
    return(xi)
}

davis.third.nonnormal <- function(b, xi){
    var.z.xi <- 1 + 2 *(b - xi)/xi
    exp.phi <- exp(-2/3 * xi * b + 1/6 * xi^2)
    p.z.xi <- (var.z.xi)^(-1/2)* exp.phi
    return(c(p.z.xi, xi, var.z.xi, exp.phi))
}

var.diffz.nonnormal <- function(rho, dist, b, X, k, xi,m2){
    n <- dim(dist)[1]
    #the hat matrix
    H <- diag(1, n) - X%*% solve(t(X) %*% X) %*% t(X)
    V <- rho^(dist)
    V.H <- V %*% H
    H.V.H <- H%*% V %*% H
    W.2 <- sum(diag(H.V.H)^2)
    V.H.2 <- V.H %*% V.H
    tr.V.H.2 <- sum(diag(V.H.2))
    G.denom <- m2*sqrt(2*tr.V.H.2 + k*W.2)
    G <- H.V.H/G.denom
    # calculate G.derivative
    V.diff <- (dist * rho^(dist - 1))
    H.V.diff.H <- H %*% V.diff %*% H
    tr.V.diff.H.V.H <- sum(diag(V.diff%*%H %*% V.H))
    sum.W.diff <- sum(diag(H.V.H)*diag(H.V.diff.H))

    G.diff.1 <- H.V.diff.H/G.denom
    G.diff.2 <- m2^2*H.V.H*(2*tr.V.diff.H.V.H + k*sum.W.diff)/G.denom^(3)
    G.diff <- G.diff.1 - G.diff.2
    ##first order
    var.diff.1 <- 2 * sum(diag(G.diff %*% G.diff)) 
    var.diff.2 <- sum(diag(G%*%G.diff%*% G.diff))
    var.diff <- var.diff.1 +  8 * xi * var.diff.2
    var.diff.sr <- sqrt(var.diff)
    var.diff.no <- sqrt(var.diff.1)
    return(c(var.diff.no, var.diff.sr, var.diff.2, xi))
}


##the integration
davis.prob.high.order.nonnormal <- function(l, u, b, dist, X, k, m6, m3, m2){
    xiu <- xi.cal(m6,m3,m2, k, u, dist, b, X)
    xil <- xi.cal(m6,m3,m2, k, l, dist, b, X)
    f <- 1/2*(davis.third.nonnormal(b, xiu)[[1]] * var.diffz.nonnormal(u, dist, b, X,k, xiu,m2)[[2]] + davis.third.nonnormal(b, xil)[[1]] * var.diffz.nonnormal(l, dist, b, X,k, xil, m2)[[2]])
    for ( i in 1:99){
        xi <- xi.cal(m6,m3, m2, k,(l+(u -l)/100 * i) , dist, b, X)
        f <- f +  davis.third.nonnormal(b, xi)[[1]] * var.diffz.nonnormal((l+(u -l)/100 * i), dist, b, X,k, xi,m2)[[2]] 
    }
    p.l <- davis.third.nonnormal(b, xil)
    var.z.l <- p.l[[3]]
    xi.l <- p.l[[2]]
    ##P(Z(L) > b) 
    p.l.z <- p.l[[4]]*(1-pnorm(0, mean=-var.z.l*xi.l, sd=sqrt(var.z.l)))*exp(var.z.l *xi.l^2/2)
    f.2 <- p.l.z+ f * (u-l)/100 *1/(2*pi)
    return(c(p.l.z, f.2)) 
}  

##the integration of no correction
davis.prob.no.correction.nonnormal <- function(l, u, b, dist, X, k, m6, m3, m2){
    xiu <- xi.cal(m6,m3,m2, k, u, dist, b, X)
    xil <- xi.cal(m6,m3,m2, k, l, dist, b, X)
    f <- 1/2*( var.diffz.nonnormal(u, dist, b, X,k, xiu,m2)[[1]] +  var.diffz.nonnormal(l, dist, b, X,k, xil, m2)[[1]])
    for ( i in 1:99){
        xi <- xi.cal(m6,m3, m2, k,(l+(u -l)/100 * i) , dist, b, X)
        f <- f +  var.diffz.nonnormal((l+(u -l)/100 * i), dist, b, X,k, xi,m2)[[1]] 

    }
    f <- pnorm(-b) + f * (u-l)/100 *1/(2*pi) *exp(-b^2/2)
    return(f) 
}  

for(i in 1:1000){

    f <- function(x, l, u, dist, X, k, m6, m3,p){
        res <- (davis.prob.high.order(l, u, x, dist, X, k, m6, m3))[2] -p)^2
    return(res)
    }
    c.value[i] <- nlm(f, l, u, dist, X, k, m6, m3,p)
}
