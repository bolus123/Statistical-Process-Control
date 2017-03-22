####################################################################################################################################################
    #Based on Champ and Jones(2004)
####################################################################################################################################################
                                                                        #This code needs package mvtnorm to obtain 
library(mvtnorm)                                                        #percentile of multivariate t distribution.
                                                                        #
####################################################################################################################################################

####################################################################################################################################################
    #get table III
####################################################################################################################################################

get.Kv.exact <- function(Kv, m, n, alpha, alternative = '2-sided'){     #The purpose of this function is 
                                                                        #to get the exact root Kv by solving the equation
                                                                        #P(-l(Kv) < Tv,i < l(Kv)) = 1 - alpha

    nu <- m * (n - 1)                                                   #Because I try to obtain Kv, 
                                                                        #degree of freedom nu = m * (n - 1)
                                                                          
        
    #c4 <- sqrt(2 / (nu)) * gamma((nu + 1) / 2) / gamma(( nu ) / 2)     #According to var(sigma0.hat.v)
                                                                        #c4 should be rely on nu and table II
    
    c4 <- sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi)             #gamma((nu + 1) / 2) / gamma(( nu ) / 2) =
                                                                        #1 / beta(nu / 2, 1 / 2) * sqrt(pi)
    
    
    c <- 1 / c4                                                         #Because I try to obtain Kv, c = 1 / c4
                                                                        #
                                                                        #
    corr.P <- diag(m)                                                   #According to the content of bottom of page 501
    corr.P[which(corr.P == 0)] <- - 1 / (m - 1)                         #Tv,i's are mutually dependent
                                                                        #And the correlation matrix is with the off-diagonal
                                                                        #elements of pij = - 1 / (m - 1)
                                                                        
    p <- (1 - alpha)                                                    #According to the top of page 503,
                                                                        #p is the integration of joint density of Tv,i

    
    l <- Kv * c * sqrt(m / (m - 1))                                     #According to the bottom of page 501
                                                                        #
    l <- rep(l, m)                                                      #According to the top of page 503,
                                                                        #each Tv,i's should be between -l and l
                                                                        
                                                                        
    p - pmvt(lower = -l, upper = l, df = nu, corr = corr.P, type = "Kshirsagar")[1]                                                                   
                                                                        #When (1 - alpha) - P(-l(Kv) < Tv,i < l(Kv)) = 0
                                                                        #I can get the root by solving this equation
    
    
}


get.Kv.approx <- function(m, n, alpha, alternative = '2-sided'){        #The purpose of this function is 
                                                                        #to get the approximate Kv

    nu <- m * (n - 1)                                                   #Because I try to obtain Kv, 
                                                                        #degree of freedom nu = m * (n - 1)
                                                                               
    #c4 <- sqrt(2 / (nu)) * gamma((nu + 1) / 2) / gamma(( nu ) / 2)     #According to var(sigma0.hat.v)
                                                                        #c4 should be rely on nu and table II
    
    c4 <- sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi)             #gamma((nu + 1) / 2) / gamma(( nu ) / 2) =
                                                                        #1 / beta(nu / 2, 1 / 2) * sqrt(pi)
    
    
    c <- 1 / c4                                                         #Because I try to obtain Kv, c = 1 / c4
                                                                     
    p <- 1 - (1 - (1 - alpha) ^ (1 / m)) / 2                            #According to the equation (8)
   
    l <- qt(p, df = nu)                                                 #l is the quantile of p
                                                                        #
    Kv <- l / c * sqrt((m - 1) / m)                                     #According to (7)
    
    return(Kv)
    
    
}                                                               
                                                                        
                                                                        
####################################################################################################################################################
    #check table III
####################################################################################################################################################
#alpha.seq <- c(0.1)
#m.seq <- 4:15
#n.seq <- 3:10
#   
#Kv.table <- matrix(NA, ncol = 5, nrow = length(m.seq) * length(n.seq) * length(alpha.seq))
#
#i <- 0
#for (alpha in alpha.seq){
#    for (m in m.seq){
#        for (n in n.seq){
#            i <- i + 1
#               
#            Kv.table[i, 1] <- alpha
#            Kv.table[i, 2] <- m
#            Kv.table[i, 3] <- n
#            Kv.table[i, 4] <- uniroot(get.Kv.exact, m = m, n = n, alpha = alpha, interval = c(1, 20), maxiter = 10000)$root
#            Kv.table[i, 5] <- get.Kv.approx(m, n, alpha)
#                
#        }
#    }
#}

                                                                        
####################################################################################################################################################
    #extend table III
####################################################################################################################################################
alpha.seq <- c(0.1, 0.05)
m.seq <- c(5:15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250, 300)
n.seq <- c(3, 5, 10)
   
Kv.table <- matrix(NA, ncol = 5, nrow = length(m.seq) * length(n.seq) * length(alpha.seq))

i <- 0
for (alpha in alpha.seq){
    for (m in m.seq){
        for (n in n.seq){
            i <- i + 1
               
            Kv.table[i, 1] <- alpha
            Kv.table[i, 2] <- m
            Kv.table[i, 3] <- n
            Kv.table[i, 4] <- uniroot(get.Kv.exact, m = m, n = n, alpha = alpha, interval = c(1, 20), maxiter = 10000)$root
            Kv.table[i, 5] <- get.Kv.approx(m, n, alpha)
                
        }
    }
}
 