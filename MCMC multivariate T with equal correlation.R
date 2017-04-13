##############################################################################################################
    #Gibbs Sampling for equal correlation
##############################################################################################################

t.f <- function(x, mu, sigma, nu){

    gamma((nu + 1) / 2) / gamma(nu / 2) / sqrt(nu * pi) / sqrt(sigma) * 
    (1 + nu ^ (-1) * (x - mu) ^ 2 / sigma) ^ (-(nu + 1) / 2)

}

t.F <- function(x, mu, sigma, nu) {

    integrate(t.f, lower = -Inf, upper = x, mu = mu, sigma = sigma, nu = nu)

}

root.t.F <- function(p, x, mu, sigma, nu) {

    p - t.F(x, mu, sigma, nu)$value

}
    
MVT.Gibbs.Sampling <- function(sample.amt, nu, Sigma, x.init = 0, burn = 100){

    n <- dim(Sigma)[1]

    x <- matrix(x.init, nrow = sample.amt + burn + 1, ncol = n)

    Sigma11 <- Sigma[1, 1]
    Sigma12 <- t(as.matrix(Sigma[1, -1]))
    Sigma21 <- as.matrix(Sigma[-1, 1])
    Sigma22 <- Sigma[-1, -1]

    Sigma.x <- Sigma11 - Sigma12 %*% solve(Sigma22) %*% Sigma21

    for (iter in 2:(sample.amt + burn + 1)) {
        
        for (i in 1:n) {
               
            if (i + 1 > n) {
            
                x.star <- c(x[iter, 0:(i - 1)])
            
            } else {
            
                x.star <- c(x[iter - 1, (i + 1):n], x[iter, 0:(i - 1)])
                
            }

            mu.x <- Sigma12 %*% solve(Sigma22) %*% x.star

            u.p <- runif(1, 0, 1)

            x.t <- qt(u.p, nu, mu.x)
            
            x[iter, i] <- uniroot(
                            root.t.F, 
                            interval = c(x.t - 3, x.t + 3), 
                            p = u.p, 
                            mu = mu.x, 
                            sigma = Sigma.x, 
                            nu = nu
                        )$root
        
        }

    }
    
    x[-c(1:(burn + 1)), ]

}

##############################################################################################################
    #Example
##############################################################################################################
#library(mvtnorm)
#
#Sigma <- matrix(c(1, -0.5, -0.5, 1), ncol = 2, byrow = TRUE)
#
##debug(MVT.Gibbs.Sampling)
##debug(root.t.F)
##debug(t.f)
##debug(t.F)
#
#
#a <- MVT.Gibbs.Sampling(10000, 100, Sigma, burn = 10000)
#b <- rmvt(10000, sigma = Sigma, 100)
#
#max.ab <- max(c(a, b))
#min.ab <- min(c(a, b))
#
#par(mfrow=c(1,2))
#plot(c(max.ab, min.ab), c(max.ab, min.ab), type = 'n', main = 'a', xlab = 'x1', ylab = 'x2') 
#points(a[,1], a[, 2])
#
#plot(c(max.ab, min.ab), c(max.ab, min.ab), type = 'n', main = 'b', xlab = 'x1', ylab = 'x2') 
#points(b[,1], b[, 2])
#
#k <- 0
#p.a <- sum(a[, 1] < k & a[, 2] < k) / 10000
#p.b <- sum(b[, 1] < k & b[, 2] < k) / 10000
#p.c <- pmvt(upper = rep(k, 2), sigma = Sigma, df = 5)