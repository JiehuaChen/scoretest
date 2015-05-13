### Maximum Score Statistic for various spatial covariance models
### Jiehua Chen jc3288@columbia.edu


# maximum score statistics with exponential spatial correlation
davis.real.data.normknown.exp <-function(dist, d, X){
    rao.score <- 0
    n <- dim(dist)[1]
    z <- rnorm(n, 0, 1)
    H <- diag(1, n) - X%*% solve(t(X) %*% X) %*% t(X)
    sigmasq <- sum((z-mean(z))^2)/(n-p)
    for (i in seq(d, 1-d, by= 0.01)) {
        V <- i^(dist) 
        rao.score2 <- (t(z)%*% H %*% V %*% H %*% z - sum(diag(H %*% V))) / sqrt(2*sum(diag(V %*% H %*% V%*%H)))
        if (rao.score2 > rao.score){
            rao.score <- rao.score2
        }
    }
    return(rao.score)
}

# maximum score statistics with gaussian spatial correlation
davis.real.data.normknown.gaussian <-function(dist, d, X){
    rao.score <- 0
    n <- dim(dist)[1]
    z <- rnorm(n, 0, 1)
    H <- diag(1, n) - X%*% solve(t(X) %*% X) %*% t(X)
    sigmasq <- sum((z-mean(z))^2)/(n-p)
    for (i in seq(d, 1-d, by= 0.01)) {
        V <- i^(dist^2) 
        rao.score2 <- (t(z)%*% H %*% V %*% H %*% z - sum(diag(H %*% V))) / sqrt(2*sum(diag(V %*% H %*% V%*%H)))
        if (rao.score2 > rao.score){
            rao.score <- rao.score2
        }
    }
    return(rao.score)
}

# maximum score statistics with spherical spatial correlation
davis.real.data.normknown.sph <-function(dist, d, X){
    rao.score <- 0
    n <- dim(dist)[1]
    z <- rnorm(n, 0, 1)
    H <- diag(1, n) - X%*% solve(t(X) %*% X) %*% t(X)
    sigmasq <- sum((z-mean(z))^2)/(n-p)
    for (i in seq(d, 8-d, by= 0.5)) {
        V <- ifelse(dist>i, 0, 1-3/2*dist/i + 1/2*(dist/i)^3)
        rao.score2 <- (t(z)%*% H %*% V %*% H %*% z - sum(diag(H %*% V))) / sqrt(2*sum(diag(V %*% H %*% V%*%H)))
        if (rao.score2 > rao.score){
            rao.score <- rao.score2
        }
    }
    return(rao.score)
}

# maximum score statistics with exponential spatial correlation with estimated variance
davis.real.data.normest.exp <-function(dist, d, X){
    rao.score <- 0
    n <- dim(dist)[1]
    z <- rnorm(n, 0, 1)
    H <- diag(1, n) - X%*% solve(t(X) %*% X) %*% t(X)
    p <- dim(X)[2]
    sigmasq <- sum((z-mean(z))^2)/(n-p)
    k<- (1/sigmasq)^2*(sum((z-mean(z))^4)/(n-p))-3
    for (i in seq(d, 1-d, by= 0.01)) {
        V <- i^(dist) 
        rao.score2 <- (t(z)%*% H %*% V %*% H %*% z/sigmasq - sum(diag(H %*% V))) / sqrt(2*sum(diag(V %*% H %*% V%*%H))+ k*sum((diag(H%*%V%*%H)^2)))
        if (rao.score2 > rao.score){
            rao.score <- rao.score2
        }
    }
    return(rao.score)
}



## two dim with the locations of the real dataset
## t distribution residuals with exponential spatial covariances
## maximum score statistic with estimated kurtosis
davis.real.data.t.exp <- function(dist, d, X){
    rao.score <- 0
    n <- dim(dist)[1]
    z <- rt(n, 8)
    H <- diag(1, n) - X%*% solve(t(X) %*% X) %*% t(X)
    p<- dim(X)[2]
    sigmasq <- sum((z-mean(z))^2)/(n-p)
    k<- (1/sigmasq)^2*(sum((z-mean(z))^4)/(n-p))-3
    for (i in seq(d, 1-d, by= 0.01)) {
        V <- i^(dist) 
        rao.score2 <- (t(z)%*% H %*% V %*% H %*% z/(sigmasq) - sum(diag(H %*% V))) / sqrt(2*sum(diag(V %*% H %*% V%*%H)) + k*sum((diag(H%*%V%*%H)^2)))
        if (rao.score2 > rao.score){
            rao.score <- rao.score2
        }
    }
    return(rao.score)
}


## two dim with the locations of the real dataset
## chi-square distribution residuals with exponential spatial covariances
## maximum score statistic with estimated kurtosis 

davis.real.data.chisq.exp <- function(dist, d, X){
    rao.score <- 0
    n <- dim(dist)[1]
    z <- rchisq(n, 1)
    z <- z-1
    H <- diag(1, n) - X%*% solve(t(X) %*% X) %*% t(X)
    p<- dim(X)[2]
    sigmasq <- sum((z-mean(z))^2)/(n-p)
    k<- (1/sigmasq)^2*(sum((z-mean(z))^4)/(n-p))-3
    for (i in seq(d, 1-d, by= 0.01)) {
        V <- i^(dist) 
        rao.score2 <- (t(z)%*% H %*% V %*% H %*% z/sigmasq - sum(diag(H %*% V))) / sqrt(2*sum(diag(V %*% H %*% V%*%H)) + k*sum((diag(H%*%V%*%H)^2)))
        if (rao.score2 > rao.score){
            rao.score <- rao.score2
        }
    }
    return(rao.score)
}




