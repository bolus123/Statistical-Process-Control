####################################################################################################################################################   
    #Based on Champ and Jones(2004), Rose and Does(1995)
####################################################################################################################################################

if (sum(c("RevoUtils", "RevoUtilsMath") %in% rownames(installed.packages())) == 2){ #Maximize the computation performance
                                                                                    #based on intel MKL library
    require(RevoUtils)                                                              #In this case, you may need packages, RevoUtils and RevoUtilsMath
    require(RevoUtilsMath)                                                          #but they are just optional libraries
    setMKLthreads(getMKLthreads())

}

####################################################################################################################################################
    #get c analytically
####################################################################################################################################################

get.c <- function(k, signal.rate, phase1 = FALSE){                                  #based on Champ and Jones' approximate Kv
                                                                                    #
    p <- 1 - signal.rate                                                            #
                                                                                    #
    l <- qt(p, df = k - 1)                                                          #
    
    c4 <- sqrt(2 / (k - 1)) * 1 / beta((k - 1) / 2, 1 / 2) * sqrt(pi)
    
                                                                                    #
    Kv <- ifelse(phase1 == FALSE, l * sqrt((k + 1) / k), l * sqrt((k - 1) / k)) / c4#
                                                                                    #
    return(Kv)                                                                      #
   
}   
 
####################################################################################################################################################
    #Check Table 3 on R&D
####################################################################################################################################################

sim.signal.rate <- function(k, c, max.sim = 100000){                                   #Based on the conditioning-unconditioning method
                                                                                       #
    Z <- rnorm(max.sim)                                                                #Generate random values from Z
    Chi.sq <- rchisq(max.sim, df = k - 1)                                              #Generate random values from Chi Square
                                                                                       #
    c4 <- sqrt(2 / (k - 1)) * 1 / beta((k - 1) / 2, 1 / 2) * sqrt(pi)                  #
                                                                                       #
    1 - mean(pnorm(Z / sqrt(k) + c / c4 / sqrt(k - 1) * sqrt(Chi.sq)))                 #Do a 1-sided test like R&D
                                                                                       #
}                                                                                      #

####################################################################################################################################################

    get.c(30, 0.0029)                                                                      #Values are from table 3
    get.c(60, 0.0020)                                                                      #Results should be closed to 3
    get.c(100, 0.0018)                                                                     #

get.c(30, 0.0029, phase1 = TRUE)  
get.c(60, 0.0020, phase1 = TRUE)  
get.c(100, 0.0018, phase1 = TRUE) 

####################################################################################################################################################

max.sim <- 10000000                                                                    #Values are from table 3
                                                                                       #
sim.signal.rate(30, 3, max.sim = max.sim)                                              #should be closed to 0.0029
sim.signal.rate(60, 3, max.sim = max.sim)                                              #should be closed to 0.0020
sim.signal.rate(100, 3, max.sim = max.sim)                                             #should be closed to 0.0018

####################################################################################################################################################
    #Check the value from analytical method
####################################################################################################################################################

sim.signal.rate(30, get.c(30, 0.0029), max.sim = max.sim)                              #should be closed to 0.0029
sim.signal.rate(60, get.c(60, 0.0020), max.sim = max.sim)                              #should be closed to 0.0020
sim.signal.rate(100, get.c(100, 0.0018), max.sim = max.sim)                            #should be closed to 0.0018

####################################################################################################################################################


