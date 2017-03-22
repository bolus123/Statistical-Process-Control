####################################################################################################################################################
                                                                        #This code needs package mvtnorm to obtain 
library(mvtnorm)                                                        #percentile of multivariate t distribution.
                                                                        #
library(LaplacesDemon)                                                  #This library is for obtain random values
                                                                        #from double exponential distribution
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
    Kv <- l / c * sqrt((m - 1) / m)                                      
    
    return(Kv)
    
    
}                                                               


####################################################################################################################################################
    #simulation
####################################################################################################################################################

getData <- function(m, n, datatype = 'norm', mean = 0, sd = 1, df = 5, ncp = 0, location = 0, scale = 1, min = -1, max = 1, shape = 1, rate = 1){
                                                                                                     #main purpose of this function is
                                                                                                     #to get data from a specific distribution
    if (datatype == 'norm'){                                                                         #
        matrix(rnorm(m * n, mean = mean, sd = sd), nrow = m, ncol = n)                               #get data from normal distribution
    } else if (datatype == 't') {                                                                    #
        matrix(rt(m * n, df = df, ncp = ncp), nrow = m, ncol = n)                                    #get data from t distribution
    } else if (datatype == 'doubexp') {                                                              #
        matrix(rlaplace(m * n, location = location, scale = scale), nrow = m, ncol = n)              #get data from double exponential distribution
    } else if (datatype == 'unif') {                                                                 #
        matrix(runif(m * n, min = min, max = max), nrow = m, ncol = n)                               #get data from uniform distribution
    } else if (datatype == 'logis'){                                                                 #
        matrix(rlogis(m * n, location = location, scale = scale), nrow = m, ncol = n)                #get data from logistic distribution
    } else if (datatype == 'gamma') {                                                                #
        matrix(rgamma(m * n, shape = shape, rate = rate), nrow = m, ncol = n)                        #get data from gamma distribution
    }                                                                                                #
    
}



sim <- function(m, n, alpha = 0.0027, getKvtype = 'approx', datatype = 'norm', mean = 0, sd = 1, df = 5, ncp = 0, location = 0, scale = 1, min = 0, max = 1, shape = 1, rate = 1, max.sim = 10000){

    nu <- m * (n - 1)                                                                                               #degree of freedom
    c4 <- sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi)                                                         #get c4
    c <- 1 / c4                                                                                                     #get c
                                                                                                                    
    T.max <- rep(NA, max.sim)
    
    for (sim in 1:max.sim){
    
        x <- getData(m, n, datatype = datatype, mean = mean, sd = sd, df = df, ncp = ncp, location = location, scale = scale, min = min, max = max, shape = shape, rate = rate)
                                                                                                                    #get data from a specific distribution
        sigma.hat <- sqrt(mean(diag(var(t(x))))) / c4                                                               #calculate sigma hat
        x.bar <- rowMeans(x)                                                                                        #calculate x bar
        x.bar.bar <- mean(x)                                                                                        #calculate x bar bar
        
        T <- c * sqrt(m / (m - 1)) * (x.bar - x.bar.bar) / sigma.hat * sqrt(n)                                      #According to (6)
                                                                                                                    #calculate charting statistics 
        T.max[sim] <- max(abs(T))
    
    }

    T.l <- quantile(T.max, 1 - alpha)
    
    T.l / c * sqrt((m - 1) / m)                                                                                     #According to (7)
    
}                                                                                                                   

set.seed(12345)
#m <- 10
#n <- 5
#alpha <- 0.05
#
#sim(m = m, n = n, alpha = alpha, getKvtype = 'exact', datatype = 'norm', max.sim = 100000)
#uniroot(get.Kv.exact, m = m, n = n, alpha = alpha, interval = c(1, 20), maxiter = 10000)$root
#
#m <- 10
#n <- 5
#alpha <- 0.1
#sim(m = m, n = n, alpha = alpha, getKvtype = 'exact', datatype = 'norm', max.sim = 100000)
#uniroot(get.Kv.exact, m = m, n = n, alpha = alpha, interval = c(1, 20), maxiter = 10000)$root


#max.sim <- 1000000
#
#alpha.seq <- c(0.05, 0.1)
#m.seq <- c(3,4,5,6,7,8,9, 10, 15,20,25,30,50,75,100,125,150,200,250,300)
#n.seq <- c(3,5,7,10)
#
#record <- as.data.frame(matrix(NA, nrow = length(alpha.seq) * length(m.seq) * length(n.seq), ncol = 4))
#
#i <- 0
#for (alpha in alpha.seq){
#
#    for (m in m.seq){
#    
#        for (n in n.seq){
#        
#            i <- i + 1
#            record[i, 1] <- alpha
#            record[i, 2] <- m
#            record[i, 3] <- n
#            record[i, 4] <- sim(m = m, n = n, alpha = alpha, getKvtype = 'exact', datatype = 'norm', max.sim = max.sim)
#            record[i, 5] <- uniroot(get.Kv.exact, m = m, n = n, alpha = alpha, interval = c(1, 20), maxiter = 10000)$root
#
#        }
#    
#    }
#
#}
