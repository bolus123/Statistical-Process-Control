####################################################################################################################################################   
    #Based on Champ and Jones(2004), Rose and Does(1995)
####################################################################################################################################################
library(mvtnorm)
library(adehabitatLT)
####################################################################################################################################################

source('https://raw.githubusercontent.com/bolus123/Statistical-Process-Control/master/MKLswitch.R')

####################################################################################################################################################
    #parts of MCMC
####################################################################################################################################################
P.empirical <- function(samp, ll, sample.amt, tails = '2-sided') {
                                                              #The purpose of this function is
    n <- dim(samp)[2]                                         #to obtain empirical probability

    upper <- ll
    
    if (tails == '2-sided'){
    
        lower <- -ll
    
    } else {
       
        lower <- rep(-Inf, n)
    
    }

    P.ind.matrix <- matrix(NA, nrow = sample.amt, ncol = n)
    P.ind <- lapply(1:n, function(i) samp[, i] < upper[i] & samp[, i] > lower[i])
    
    for (i in 1:n) {
    
        P.ind.matrix[, i] <- P.ind[[i]]
    
    }
    
    sum(rowSums(P.ind.matrix) == n) / sample.amt

}

root.P.empirical <- function(p, L, samp, nlength, sample.amt, tails = '2-sided'){
                                                             #The purpose of this function is
    ll <- rep(L, nlength)                                    #to obtain appropriate L.
    
    p - P.empirical(samp, ll, sample.amt, tails)             #L is the root of this function

}


####################################################################################################################################################
    #Multivariate T based on MCMC
####################################################################################################################################################

t.f <- function(x, mu, sigma, nu){
                                                                              #The purpose of this function
    gamma((nu + 1) / 2) / gamma(nu / 2) / sqrt(nu * pi) / sqrt(sigma) *       #is to calculate pdf of t distribution
    (1 + nu ^ (-1) * (x - mu) ^ 2 / sigma) ^ (-(nu + 1) / 2)                  #based on Peng Li(2016)

}

t.F <- function(x, mu, sigma, nu) {                                            #The purpose of this function is 
                                                                               #to calculate cdf of t distribution
    integrate(t.f, lower = -Inf, upper = x, mu = mu, sigma = sigma, nu = nu)   #
                                                                               #
}

root.t.F <- function(p, x, mu, sigma, nu) {                                    #The purpose of this function is
                                                                               #obtain appropriate quantile 
    p - t.F(x, mu, sigma, nu)$value                                            #of t distribution
                                                                               #
}
    
MVT.Gibbs.Sampling <- function(sample.amt, nu, Sigma, x.init = 0, burn = 1000){

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
            
            p.t.1 <- root.t.F(u.p, x.t, mu.x, Sigma.x, nu)
            
            for (j in 1:1000){

                x.t.alt <- x.t + 0.1 * j * (-1) ^ j
            
                p.t.2 <- root.t.F(u.p, x.t.alt, mu.x, Sigma.x, nu)
                
                if (p.t.1 * p.t.2 < 0) break 
                
            }
            
            x[iter, i] <- uniroot(
                            root.t.F, 
                            interval = c(x.t, x.t.alt), 
                            p = u.p, 
                            mu = mu.x, 
                            sigma = Sigma.x, 
                            nu = nu
                        )$root
        
        }

    }
    
    x[-c(1:(burn + 1)), ]

}




MVT.F.Gibbs.Sampling <- function(ll, sample.amt, nu, Sigma, tails = '2-sided', x.init = 0, burn = 1000){

    n <- dim(Sigma)[1]
    
    samp <- MVT.Gibbs.Sampling(sample.amt, nu, Sigma, x.init = x.init, burn = burn)

    P.empirical(samp, ll, sample.amt, tails)
    

}



MVT.Q.Gibbs.Sampling <- function(p, sample.amt, nu, Sigma, tails = '2-sided', x.init = 0, burn = 1000, search.interval = c(1, 5)){

    n <- dim(Sigma)[1]
    
    samp <- MVT.Gibbs.Sampling(sample.amt, nu, Sigma, x.init = x.init, burn = burn)
    
    uniroot(
        root.P.empirical, 
        interval = search.interval, 
        samp = samp, 
        p = p, 
        nlength = n, 
        sample.amt = sample.amt, 
        tails = tails
    )$root

}


####################################################################################################################################################
    #Multivariate Normal based on MCMC
####################################################################################################################################################


MVN.Gibbs.Sampling <- function(sample.amt, Sigma, x.init = 0, burn = 1000){

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
            
            x[iter, i] <- rnorm(1, mu.x, sqrt(Sigma.x))
        
        }

    }
    
    x[-c(1:(burn + 1)), ]

}

MVN.F.Gibbs.Sampling <- function(ll, sample.amt, Sigma, tails = '2-sided', x.init = 0, burn = 1000){

    n <- dim(Sigma)[1]
    
    samp <- MVN.Gibbs.Sampling(sample.amt, Sigma, x.init = x.init, burn = burn)
    
    P.empirical(samp, ll, sample.amt, tails)
    
    

}

MVN.Q.Gibbs.Sampling <- function(p, sample.amt, Sigma, tails = '2-sided', x.init = 0, burn = 1000, search.interval = c(1, 5)){

    n <- dim(Sigma)[1]
    
    samp <- MVN.Gibbs.Sampling(sample.amt, Sigma, x.init = x.init, burn = burn)
    
    uniroot(
        root.P.empirical, 
        interval = search.interval, 
        samp = samp, 
        p = p, 
        nlength = n, 
        sample.amt = sample.amt, 
        tails = tails
    )$root

}

####################################################################################################################################################
    #parts of getting L
####################################################################################################################################################

c4.f <- function(nu) sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi)             #c4.function

corr.f <- function(m, off.diag = - 1 / (m - 1)){
                                                                                   #correlation matrix
    crr <- diag(m)
    crr[which(crr == 0)] <- off.diag

    crr
    
}

####################################################################################################################################################
    #get L by multivariate T
####################################################################################################################################################

get.L.mvt <- function(
                 m 
                 ,nu
                 ,FAP = 0.1
                 #,ARL = 370
                 ,Phase1 = TRUE
                 ,off.diag = NULL
                 ,alternative = '2-sided'
                 ,maxiter = 10000
                 #,MCMC = FALSE
                 #,MCMC.search.interval = c(1, 5)
                 #,MCMC.maxsim = 10000
                 #,MCMC.burn = 1000

){
                                                                              #The purpose of this function is 
    MCMC <- FALSE                                                             #to obtain L and K based on
                                                                              #multivariate T.
                                                                              #MCMC part is not available now.

    #if (off.diag == NULL) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
    
    corr.P <- corr.f(m = m, off.diag = off.diag)                              #get correlation matrix with equal correlations

    pu <- 1 - FAP
    
    if (MCMC == TRUE) {
        #MVN.Q.Gibbs.Sampling(
        #    pu, 
        #    MCMC.maxsim, 
        #    corr.P, 
        #    tails = alternative, 
        #    burn = MCMC.burn, 
        #    search.interval = MCMC.search.interval
        #)
    
        
    
    } else {
 
        L <- ifelse(
                alternative == '2-sided',
                qmvt(pu, df = nu, sigma = corr.P, tail = 'both.tails', maxiter = maxiter)$quantile,
                qmvt(pu, df = nu, sigma = corr.P, maxiter = maxiter)$quantile
            )
                                                      #get L by multivariate T
        

    }
    
    
    K <- L * c4.f(nu) * sqrt((m - 1) / m)             #get K
   
    list(L = L, K = K)

}

####################################################################################################################################################
    #get L by multivariate Normal
####################################################################################################################################################

#root.mvn.F <- function(K, m, nu, Y, sigma, pu, alternative = '2-sided') {
#                                                            #The purpose of this function is
#    s <- length(Y)                                          #to obtain appropriate K
#                                                            #by multivariate normal
#    L <- K / sqrt((m - 1) / m * nu) * Y * c4.f(nu)          #
#
#    pp <- lapply(
#            1:s,
#            function(i){
#            
#                LL <- rep(L[i], m)
#                ifelse(
#                    alternative == '2-sided',
#                    pmvnorm(lower = -LL, upper = LL, sigma = sigma),
#                    pmvnorm(lower = rep(-Inf, m), upper = LL, sigma = sigma)
#                )
#            
#            } 
#    
#    )
#    
#    pp <- mean(unlist(pp))
#    
#    pu - pp
#
#
#}


joint.pdf.mvn.chisq <- function(Y, K, m, nu, sigma, alternative = '2-sided') {

    s <- length(Y)

    L <- K / sqrt((m - 1) / m * nu) * Y * c4.f(nu)
    
    dpp <- lapply(
            1:s,
            function(i){

                LL <- rep(L[i], m)
                
                ifelse(
                    alternative == '2-sided',
                    pmvnorm(lower = -LL, upper = LL, sigma = sigma),
                    pmvnorm(lower = rep(-Inf, m), upper = LL, sigma = sigma)
                )
            
            } 
    
    )    

    dpp <- unlist(dpp)
    
    dpp * dchi(Y, df = nu)


}



root.mvn.F <- function(K, m, nu, sigma, pu, alternative = '2-sided') {
                                                            #The purpose of this function is
    #s <- length(Y)                                          #to obtain appropriate L
    #                                                        #by multivariate normal
    #L <- K / sqrt((m - 1) / m * nu) * Y * c4.f(nu)          #
    #
    #pp <- lapply(
    #        1:s,
    #        function(i){
    #        
    #            LL <- rep(L[i], m)
    #            ifelse(
    #                alternative == '2-sided',
    #                pmvnorm(lower = -LL, upper = LL, sigma = sigma),
    #                pmvnorm(lower = rep(-Inf, m), upper = LL, sigma = sigma)
    #            )
    #        
    #        } 
    #
    #)
    #
    #pp <- mean(unlist(pp))
    
    pp <- integrate(joint.pdf.mvn.chisq, lower = 0, upper = Inf, K = K, m = m, nu = nu, sigma = sigma, alternative = alternative)$value
    
    pu - pp


}

#sigma = corr.f(20)
#joint.pdf.mvn.chisq(1, 3, 20, 80, sigma, alternative = '2-sided')
#root.mvn.F(3, 20, 80, sigma, 0.7)

get.L.mvn <- function(
                 m
                 ,nu
                 ,FAP = 0.1
                 #,ARL = 370
                 ,Phase1 = TRUE
                 ,off.diag = NULL
                 ,alternative = '2-sided'
                 #,maxsim = 10000
                 ,maxiter = 10000
                 #,MCMC = FALSE
                 #,MCMC.search.interval = c(1, 5)
                 #,MCMC.maxsim = 10000
                 #,MCMC.burn = 1000

){                                                          #The purpose of this function is 
                                                            #to obtain L and K based on
    MCMC <- FALSE                                           #multivariate normal.
                                                            #MCMC part is not available now.
    if (is.null(off.diag)) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
    
    
    corr.P <- corr.f(m = m, off.diag = off.diag)  

    pu <- 1 - FAP
    
    #Y <- sqrt(rchisq(maxsim, df = nu))
    
    if (MCMC == TRUE) {
    
        #MVN.Q.Gibbs.Sampling(
        #    pu, 
        #    MCMC.maxsim, 
        #    corr.P, 
        #    tails = alternative, 
        #    burn = MCMC.burn, 
        #    search.interval = MCMC.search.interval
        #)
    
        
    
    } else {

        #L <- ifelse(
        #        alternative == '2-sided',
        #        qmvnorm(pu, sigma = corr.P, tail = 'both.tails', maxiter = maxiter)$quantile,
        #        qmvnorm(pu, sigma = corr.P, maxiter = maxiter)$quantile
        #    )

        K <- uniroot(
                root.mvn.F,
                interval = c(1, 7),
                m = m,
                nu = nu,
                #Y = Y,
                sigma = corr.P,
                pu = pu,
                alternative = alternative,
                maxiter = maxiter
        )$root

    }
    
    
    #K <- mean(L / Y * sqrt((m - 1) / m * nu)) #/ c4.f(nu))
    
    L <- K / c4.f(nu) * sqrt(m / (m - 1))
    
    list(L = L, K = K)
    

}

#get.L.mvn(20, 80)

####################################################################################################################################################

get.L <- function(
            m
            ,nu
            ,FAP = 0.1
            #,ARL = 370
            #,Phase1 = TRUE
            ,off.diag = NULL
            ,alternative = '2-sided'
            ,maxiter = 10000
            ,method = 'direct'
            #,indirect.maxsim = 10000
            #,MCMC = FALSE
            #,MCMC.search.interval = c(1, 5)
            #,MCMC.maxsim = 10000
            #,MCMC.burn = 1000

){                                                  #The purpose of this function is to obtain L and K
                                                    #by multivariate T or multivariate normal.
                                                    #Multivariate normal is so time-consuming
                                                    #that I do not recommend.
    Phase1 <- TRUE
    
    
    if (is.null(off.diag)) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))


    if (method == 'direct') {                       #using multivariate T to obtain L and K
    
        get.L.mvt(
            m = m
            ,nu = nu
            ,FAP = FAP
            ,Phase1 = Phase1
            ,off.diag = off.diag
            ,alternative = alternative
            ,maxiter = maxiter
        )
    
    } else {                                       #using multivariate normal to obtain L and K
    
        get.L.mvn(
            m = m
            ,nu = nu
            ,FAP = FAP
            ,Phase1 = Phase1
            ,off.diag = off.diag
            ,alternative = alternative
            #,maxsim = indirect.maxsim
            ,maxiter = maxiter
        )
    
    }


}


get.L(
            20
            ,80
            ,FAP = 0.1
            ,Phase1 = TRUE
            ,off.diag = NULL
            ,alternative = '2-sided'
            ,maxiter = 10000
            ,method = 'direct'
)
get.L(
            20
            ,80
            ,FAP = 0.1
            ,Phase1 = TRUE
            ,off.diag = NULL
            ,alternative = '2-sided'
            ,maxiter = 10000
            ,method = 'indirect'
)

####################################################################################################################################################
    #Example
####################################################################################################################################################
#a <- get.L(
#            m = 15
#            ,nu = 30
#            ,FAP = 0.1
#            ,Phase1 = TRUE
#            ,alternative = '2-sided'
#            ,maxiter = 10000
#            ,method = 'direct'
#
#
#)   
#
#
#b <- get.L(
#            m = 15
#            ,nu = 30
#            ,FAP = 0.1
#            ,Phase1 = TRUE
#            ,alternative = '2-sided'
#            ,maxiter = 10000
#            ,method = 'indirect'
#            ,indirect.maxsim = 10000
#
#)
