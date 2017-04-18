####################################################################################################################################################   
    #Based on Champ and Jones(2004), Rose and Does(1995)
####################################################################################################################################################
library(mvtnorm)
####################################################################################################################################################


if (sum(c("RevoUtils", "RevoUtilsMath") %in% rownames(installed.packages())) == 2){ #Maximize the computation performance
                                                                                    #based on intel MKL library
    require(RevoUtils)                                                              #In this case, you may need packages, RevoUtils and RevoUtilsMath
    require(RevoUtilsMath)                                                          #but they are just optional libraries
    setMKLthreads(getMKLthreads())

}
####################################################################################################################################################
    #parts of MCMC
####################################################################################################################################################
P.empirical <- function(samp, ll, sample.amt, tails = '2-sided') {

    n <- dim(samp)[2]

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

    ll <- rep(L, nlength)
    
    p - P.empirical(samp, ll, sample.amt, tails)

}


####################################################################################################################################################
    #Multivariate T based on MCMC
####################################################################################################################################################

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
                 ,off.diag = ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
                 ,alternative = '2-sided'
                 ,maxiter = 10000
                 #,MCMC = FALSE
                 #,MCMC.search.interval = c(1, 5)
                 #,MCMC.maxsim = 10000
                 #,MCMC.burn = 1000

){

    MCMC <- FALSE

    corr.P <- corr.f(m = m, off.diag = off.diag)  
    #pu <- ifelse(alternative == '2-sided', 1 - (1 - (1 - 1 / ARL) ^ (1 / m)) / 2, (1 - 1 / ARL) ^ (1 / m))

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
 
        

    }
    
    
    K <- L * c4.f(nu) * sqrt((m - 1) / m)
   
    list(L = L, K = K)

}

####################################################################################################################################################
    #get L by multivariate Normal
####################################################################################################################################################

root.mvn.F <- function(K, m, nu, Y, sigma, pu, alternative = '2-sided') {

    s <- length(Y)

    L <- K / sqrt((m - 1) / m * nu) * Y * c4.f(nu)

    pp <- lapply(
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
    
    pp <- mean(unlist(pp))
    
    pu - pp


}


get.L.mvn <- function(
                 m
                 ,nu
                 ,FAP = 0.1
                 #,ARL = 370
                 ,Phase1 = TRUE
                 ,off.diag = ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
                 ,alternative = '2-sided'
                 ,maxsim = 10000
                 ,maxiter = 10000
                 #,MCMC = FALSE
                 #,MCMC.search.interval = c(1, 5)
                 #,MCMC.maxsim = 10000
                 #,MCMC.burn = 1000

){

    MCMC <- FALSE
    
    corr.P <- corr.f(m = m, off.diag = off.diag)  
    #pu <- ifelse(alternative == '2-sided', 1 - (1 - (1 - 1 / ARL) ^ (1 / m)) / 2, (1 - 1 / ARL) ^ (1 / m))

    pu <- 1 - FAP
    
    Y <- sqrt(rchisq(maxsim, df = nu))
    
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
                Y = Y,
                sigma = corr.P,
                pu = pu,
                alternative = alternative
        )$root

    }
    
    
    #K <- mean(L / Y * sqrt((m - 1) / m * nu)) #/ c4.f(nu))
    
    L <- K / c4.f(nu) * sqrt(m / (m - 1))
    
    list(L = L, K = K)
    

}

####################################################################################################################################################

get.L <- function(
            m
            ,nu
            ,FAP = 0.1
            #,ARL = 370
            ,Phase1 = TRUE
            ,off.diag = ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
            ,alternative = c('2-sided')
            ,maxsim = 10000
            ,maxiter = 10000
            ,method = c('direct')
            #,MCMC = FALSE
            #,MCMC.search.interval = c(1, 5)
            #,MCMC.maxsim = 10000
            #,MCMC.burn = 1000

){

    if (method == 'direct') {
    
        get.L.mvt(
            m = m
            ,nu = nu
            ,FAP = FAP
            ,Phase1 = Phase1
            ,off.diag = off.diag
            ,alternative = alternative
            ,maxsim = maxsim
            ,maxiter = maxiter
            ,method = method
        )
    
    } else {
    
        get.L.mvn(
            m = m
            ,nu = nu
            ,FAP = FAP
            ,Phase1 = Phase1
            ,off.diag = off.diag
            ,alternative = alternative
            ,maxsim = maxsim
            ,maxiter = maxiter
            ,method = method
        )
    
    }


}





###test
#
#get.c <- function(L, m, nu, Phase1 = TRUE) ifelse(Phase1 == TRUE, L * sqrt((m - 1) / m), L * sqrt((m + 1) / m))
#
#l.R.D <- rep(NA, 3)
#l.R.D[1] <- get.L(30, 29, 1 - (1 - 0.0029)^30, Phase1 = FALSE, alternative = '1-sided')
#l.R.D[2] <- get.L(60, 59, 1 - (1 - 0.0020)^60, Phase1 = FALSE, alternative = '1-sided')
#l.R.D[3] <- get.L(100, 99, 1 - (1 - 0.0018)^100, Phase1 = FALSE, alternative = '1-sided')
#get.c(l.R.D[1], 30, 29, Phase1 = FALSE)
#get.c(l.R.D[2], 60, 59, Phase1 = FALSE)
#get.c(l.R.D[3], 100, 99, Phase1 = FALSE)
#
#get.k <- function(L, m, nu, Phase1 = TRUE) ifelse(Phase1 == TRUE, L * c4.f(nu) * sqrt((m - 1) / m), L * c4.f(nu) * c4.f(nu) * sqrt((m - 1) / m))
#
#l.R.D <- rep(NA, 3)
#l.C.J[1] <- get.L(5, 20, 0.05, Phase1 = TRUE, alternative = '2-sided')
#l.C.J[2] <- get.L(10, 40, 0.05, Phase1 = TRUE, alternative = '2-sided')
#l.C.J[3] <- get.L(15, 60, 0.05, Phase1 = TRUE, alternative = '2-sided')
#get.k(l.C.J[1], 5, 20, Phase1 = TRUE)
#get.k(l.C.J[2], 10, 40, Phase1 = TRUE)
#get.k(l.C.J[3], 15, 60, Phase1 = TRUE)
#
#
#
#m.seq <- c(
#    3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30,
#    50, 75, 100, 125, 150, 200, 250, 300
#)
#
#nu.seq <- c(
#    3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30,
#    50, 75, 100, 125, 150, 200, 250, 300
#)
#
#alpha.seq <- c(
#    0.1, 0.05
#)
#
#record.phase1 <- matrix(NA, ncol = 4, 
#        nrow = length(m.seq) * length(nu.seq) * length(alpha.seq)
#)
#
#i <- 0
#for (m in m.seq){
#    for (nu in nu.seq){
#        for (alpha in alpha.seq){
#        
#            #if (nu < m) next
#        
#            i <- i + 1
#            record.phase1[i, 1] <- alpha
#            record.phase1[i, 2] <- m
#            record.phase1[i, 3] <- nu
#            record.phase1[i, 4] <- get.L(m, nu, alpha, Phase1 = TRUE, alternative = '2-sided')
#        
#        }
#    }
#}
#
#rho.seq <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
#
#record.phase2 <- matrix(NA, ncol = 5, 
#        nrow = length(m.seq) * length(n.seq) * length(alpha.seq) * length(rho.seq)
#)
#
#i <- 0
#for (rho in rho.seq){ 
#    for (m in m.seq){
#        for (n in n.seq){
#            for (alpha in alpha.seq){
#                i <- i + 1
#                record.phase2[i, 1] <- alpha
#                record.phase2[i, 2] <- m
#                record.phase2[i, 3] <- n
#                record.phase2[i, 4] <- rho
#                record.phase2[i, 5] <- get.L(m, m * (n - 1), alpha, Phase1 = FALSE, alternative = '2-sided')
#            
#            }
#        }
#    }
#}
#
#
#####################################################################################################################################################
#    #get alpha by multivariate normal(the conditioning-unconditioning method)
#####################################################################################################################################################
#get.alpha.cond.uncond <- function(
#                                C, 
#                                m, 
#                                nu, 
#                                Phase1 = TRUE, 
#                                off.diag = ifelse(Phase1 == TRUE, - 1 /(m - 1), 0), 
#                                alternative = '2-sided',
#                                max.sim = 10000, 
#                                cores = 3
#){                                                                                      #Based on the conditioning-unconditioning method
#                                                                                        #
#    require(parallel)
#
#    Y <- rchisq(max.sim, df = nu)                                                       #Generate random values from Chi Square
#    
#    sqrt.Y <- sqrt(Y)                                                                   #get random values chi with nu degrees of freedom
#                                                                                        
#    c4 <- c4.f(nu)                                                                      #define c4
#         
#    corr.P <- corr.f(m, off.diag = off.diag)                                            #define correlation matrix
#    
#    if (Phase1 == TRUE){                                                                #define the upper bounds for integration
#                                                                                        #
#        l <- C / c4 * sqrt(m / (m - 1)) / sqrt(nu) * sqrt.Y                             #
#                                                                                        #
#    } else {                                                                            #
#                                                                                        #
#        l <- C / c4 * sqrt(m / (m + 1)) / sqrt(nu) * sqrt.Y                             #
#                                                                                        #
#    }                                                                                   #
#                                                          
#    cl <- makeCluster(cores)
#
#    clusterExport(cl, c('l', 'm', 'alternative', 'corr.P'), envir = environment())
#    clusterEvalQ(cl, library(mvtnorm))     
#    
#    p <- clusterApply(
#            cl,
#            1:max.sim, 
#            function(x){
#            
#                pp = ifelse(                                                                    #get the integration of multivariate normal
#                        alternative == '2-sided',                                               #
#                        pmvnorm(lower = -rep(l[x], m), upper = rep(l[x], m), corr = corr.P),    #
#                        pmvnorm(lower = -Inf, upper = rep(l[x], m), corr = corr.P)              #
#                )                                                                               #
#                
#                1 - pp
#            }
#    
#    )
#    
#    stopCluster(cl)
#    
#    mean(unlist(p))
#                                                                                       
#}     
#
##debug(get.alpha.cond.uncond)
#get.alpha.cond.uncond(
#                3, 
#                30, 
#                29, 
#                Phase1 = FALSE, 
#                off.diag = 0, 
#                alternative = '1-sided',
#                max.sim = 10000, 
#                cores = 3
#)
#
#    
#####################################################################################################################################################
#
