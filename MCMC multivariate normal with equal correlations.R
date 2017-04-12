##############################################################################################################
    #Gibbs Sampling for equal correlation
##############################################################################################################

MVN.Gibbs.Sampling <- function(sample.amt, Sigma, x.init = 0, burn = 100){

    n <- dim(Sigma)[1]

    x <- matrix(x.init, nrow = sample.amt + burn + 1, ncol = n)
    
    TRS <- diag(n)

    Sigma11 <- Sigma[1, 1]
    Sigma12 <- t(as.matrix(Sigma[1, -1]))
    Sigma21 <- as.matrix(Sigma[-1, 1])
    Sigma22 <- Sigma[-1, -1]

    Sigma.x <- Sigma11 - Sigma12 %*% solve(Sigma22) %*% Sigma21

    for (iter in 2:(sample.amt + burn + 1)) {
        
        for (i in 1:n) {
        
            TRS.star <- rbind(TRS[i, ], TRS[-i, ])
    
        
            if (i + 1 > n) {
            
                x.star <- c(x[iter, 0:(i - 1)])
            
            } else {
            
                x.star <- c(x[iter - 1, (i + 1):n], x[iter, 0:(i - 1)])
                
            }

            mu.x <- Sigma12 %*% solve(Sigma22) %*% x.star
            
            x[iter, i] <- rnorm(1, mu.x, sqrt(Sigma.x))
        
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
#
#a <- MVN.Gibbs.Sampling(10000, Sigma, burn = 1000)
#b <- rmvnorm(10000, sigma = Sigma)
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
#p.b <- pmvnorm(upper = rep(k, 2), sigma = Sigma)