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
    #get L by multivariate T
####################################################################################################################################################

c4.f <- function(nu) sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi)             #c4.function

corr.f <- function(m, off.diag = - 1 / (m - 1)){
                                                                                   #correlation matrix
    crr <- diag(m)
    crr[which(crr == 0)] <- off.diag

    crr
    
}

get.F <- function(                                                                 #get the integration of multivariate T
                    L, 
                    m, 
                    nu, 
                    alpha = 0.05, 
                    Phase1 = TRUE, 
                    off.diag = ifelse(Phase1 == TRUE, - 1 /(m - 1), 0), 
                    alternative = '2-sided'
) {     

    corr.P <- corr.f(m = m, off.diag = off.diag)                                  #get correlation matrix
 
    ll <- rep(L, m)                                                               #define bounds
    
    p <- 1 - alpha                                                                #coverage
    
    ifelse(
            alternative == '2-sided',
            p - pmvt(lower = -ll, upper = ll, df = nu, corr = corr.P),
            p - pmvt(lower = rep(-Inf, m), upper = ll, df = nu, corr = corr.P)
    )
    
 
}


get.L <- function(                                                                #get the limits of multivariate T 
                    m, 
                    nu, 
                    alpha = 0.05, 
                    Phase1 = TRUE,
                    off.diag = ifelse(Phase1 == TRUE, - 1 /(m - 1), 0),
                    alternative = '2-sided',
                    search.interval = c(1, 20),
                    tol = 1e-6,
                    maxiter = 100000
) {

    L <- uniroot(get.F, m = m, nu = nu, alpha = alpha, Phase1 = Phase1,           #find the root of integration of multivariate T
        off.diag = off.diag, alternative = alternative, interval = search.interval, 
        tol = tol, maxiter = maxiter)
    

    L$root

}

###test

get.c <- function(L, m, nu, Phase1 = TRUE) ifelse(Phase1 == TRUE, L * sqrt((m - 1) / m), L * sqrt((m + 1) / m))

l.R.D <- rep(NA, 3)
l.R.D[1] <- get.L(30, 29, 1 - (1 - 0.0029)^30, Phase1 = FALSE, alternative = '1-sided')
l.R.D[2] <- get.L(60, 59, 1 - (1 - 0.0020)^60, Phase1 = FALSE, alternative = '1-sided')
l.R.D[3] <- get.L(100, 99, 1 - (1 - 0.0018)^100, Phase1 = FALSE, alternative = '1-sided')
get.c(l.R.D[1], 30, 29, Phase1 = FALSE)
get.c(l.R.D[2], 60, 59, Phase1 = FALSE)
get.c(l.R.D[3], 100, 99, Phase1 = FALSE)

get.k <- function(L, m, nu, Phase1 = TRUE) ifelse(Phase1 == TRUE, L * c4.f(nu) * sqrt((m - 1) / m), L * c4.f(nu) * c4.f(nu) * sqrt((m - 1) / m))

l.R.D <- rep(NA, 3)
l.C.J[1] <- get.L(5, 20, 0.05, Phase1 = TRUE, alternative = '2-sided')
l.C.J[2] <- get.L(10, 40, 0.05, Phase1 = TRUE, alternative = '2-sided')
l.C.J[3] <- get.L(15, 60, 0.05, Phase1 = TRUE, alternative = '2-sided')
get.k(l.C.J[1], 5, 20, Phase1 = TRUE)
get.k(l.C.J[2], 10, 40, Phase1 = TRUE)
get.k(l.C.J[3], 15, 60, Phase1 = TRUE)



m.seq <- c(
    3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30,
    50, 75, 100, 125, 150, 200, 250, 300
)

nu.seq <- c(
    3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30,
    50, 75, 100, 125, 150, 200, 250, 300
)

alpha.seq <- c(
    0.1, 0.05
)

record.phase1 <- matrix(NA, ncol = 4, 
        nrow = length(m.seq) * length(nu.seq) * length(alpha.seq)
)

i <- 0
for (m in m.seq){
    for (nu in nu.seq){
        for (alpha in alpha.seq){
        
            #if (nu < m) next
        
            i <- i + 1
            record.phase1[i, 1] <- alpha
            record.phase1[i, 2] <- m
            record.phase1[i, 3] <- nu
            record.phase1[i, 4] <- get.L(m, nu, alpha, Phase1 = TRUE, alternative = '2-sided')
        
        }
    }
}

rho.seq <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

record.phase2 <- matrix(NA, ncol = 5, 
        nrow = length(m.seq) * length(n.seq) * length(alpha.seq) * length(rho.seq)
)

i <- 0
for (rho in rho.seq){ 
    for (m in m.seq){
        for (n in n.seq){
            for (alpha in alpha.seq){
                i <- i + 1
                record.phase2[i, 1] <- alpha
                record.phase2[i, 2] <- m
                record.phase2[i, 3] <- n
                record.phase2[i, 4] <- rho
                record.phase2[i, 5] <- get.L(m, m * (n - 1), alpha, Phase1 = FALSE, alternative = '2-sided')
            
            }
        }
    }
}


####################################################################################################################################################
    #get alpha by multivariate normal(the conditioning-unconditioning method)
####################################################################################################################################################
get.alpha.cond.uncond <- function(
                                C, 
                                m, 
                                nu, 
                                Phase1 = TRUE, 
                                off.diag = ifelse(Phase1 == TRUE, - 1 /(m - 1), 0), 
                                alternative = '2-sided',
                                max.sim = 10000, 
                                cores = 3
){                                                                                      #Based on the conditioning-unconditioning method
                                                                                        #
    require(parallel)

    Y <- rchisq(max.sim, df = nu)                                                       #Generate random values from Chi Square
    
    sqrt.Y <- sqrt(Y)                                                                   #get random values chi with nu degrees of freedom
                                                                                        
    c4 <- c4.f(nu)                                                                      #define c4
         
    corr.P <- corr.f(m, off.diag = off.diag)                                            #define correlation matrix
    
    if (Phase1 == TRUE){                                                                #define the upper bounds for integration
                                                                                        #
        l <- C / c4 * sqrt(m / (m - 1)) / sqrt(nu) * sqrt.Y                             #
                                                                                        #
    } else {                                                                            #
                                                                                        #
        l <- C / c4 * sqrt(m / (m + 1)) / sqrt(nu) * sqrt.Y                             #
                                                                                        #
    }                                                                                   #
                                                          
    cl <- makeCluster(cores)

    clusterExport(cl, c('l', 'm', 'alternative', 'corr.P'), envir = environment())
    clusterEvalQ(cl, library(mvtnorm))     
    
    p <- clusterApply(
            cl,
            1:max.sim, 
            function(x){
            
                pp = ifelse(                                                                    #get the integration of multivariate normal
                        alternative == '2-sided',                                               #
                        pmvnorm(lower = -rep(l[x], m), upper = rep(l[x], m), corr = corr.P),    #
                        pmvnorm(lower = -Inf, upper = rep(l[x], m), corr = corr.P)              #
                )                                                                               #
                
                1 - pp
            }
    
    )
    
    stopCluster(cl)
    
    mean(unlist(p))
                                                                                       
}     

#debug(get.alpha.cond.uncond)
get.alpha.cond.uncond(
                3, 
                30, 
                29, 
                Phase1 = FALSE, 
                off.diag = 0, 
                alternative = '1-sided',
                max.sim = 10000, 
                cores = 3
)

    
####################################################################################################################################################
