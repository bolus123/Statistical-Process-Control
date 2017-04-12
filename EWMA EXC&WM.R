##########################################################################
    #No extra package is needed
##########################################################################

Z.f <- function(x, lambda, Z0 = 0){                                  #Basic structure of EWMA chart
                                                                     #Z.i = lambda * X + (1 - lambda) * Z.i-1
    tt <- length(x)                                                  #
                                                                     #
    z <- rep(Z0, tt + 1)                                             #
                                                                     #
    for (i in 2:(tt + 1)){                                           #
                                                                     #
        z[i] <- lambda * x[i - 1] + (1 - lambda) * z[i - 1]          #
                                                                     #
    }                                                                #
                                                                     #
    z[-1]                                                            #
                                                                     #
}                                                                    #

##########################################################################

EX.f <- function(x, y, percentile){
                                                                     #The purpose of this function is
    if (is.vector(y)) y <- t(as.matrix(y))                           #to get amount of Y exceed X(r)

    m <- length(x)
    n <- dim(y)[2]
    t <- dim(y)[1]
                                                                     
    r <- ceiling((m + 1) * percentile)                               #Get specific rank corresponding to percentile

    rr <- rank(x)                                                    #rank of all X
    
    x.r <- x[which.min(abs(rr - r))[1]]                              #Find the X(r), but rank may have ties
    
    rowSums(y > x.r)


}



EWMA.EX.Chart <- function(x, y, percentile, lambda = 0.1, L = 3, tt = NULL, graph = TRUE, Z0 = 0){
                                                                                   #This purpose of this function is
    if (is.vector(y)) y <- t(as.matrix(y))                                         #to draw EWMA Exceedance chart
                                                                                   #
    m <- length(x)                                                                 #
    n <- dim(y)[2]                                                                 #
    if (is.null(tt)){                                                              #When t = NULL, t will be obtain from the very begining
                                                                                   #When t != NULL, t will be fixed for every charting staitics
        t <- c(1:dim(y)[1])                                                        #
                                                                                   #
    } else {
    
        t <- rep(tt, dim(y)[1])
    
    }

    r <- ceiling((m + 1) * percentile)
    a <- r / (m + 1)
    
    X <- EX.f(x, y, percentile)                                                     #EX statistics
    
    Z <- Z.f(X, lambda, Z0 = Z0)                                                    #Get charting statistics
    
    CL <- n * (1 - a) * (1 - (1 - lambda) ^ t)                                      #Get CL
    V <- sqrt(n * a * (1 - a) / (m + 2) * (n * (1 - (1 - lambda) ^ t)^2 + lambda * (m + 1) / (2 - lambda) * (1 - (1 - lambda) ^ (2 * t))))
    
    UCL <- CL + L * V                                                               #Get UCL
     
    LCL <- CL - L * V                                                               #Get LCL
    
    if (graph == TRUE){
        plot(c(1, dim(y)[1]), c(3 / 4 * min(c(Z, LCL)), 5 / 4 * max(c(Z, UCL))), type = 'n', xlab = 't', ylab = 'Z')
        points(c(1:dim(y)[1]), UCL, type = 'l', lty = 2)
        points(c(1:dim(y)[1]), CL, type = 'l')
        points(c(1:dim(y)[1]), LCL, type = 'l', lty = 2)
        points(c(1:dim(y)[1]), Z, type = 'l')
    }
    return(list(CL = CL, LCL = LCL, UCL = UCL, Z = Z))
    
}

#x <- rnorm(100)
#y <- matrix(rnorm(100), ncol = 5)

#debug(EWMA.EX.Chart)

#EWMA.EX.Chart(x, y, 0.5, lambda = 0.1, L = 3, tt = 100000)


##########################################################################################

WM.f <- function(x, y){
                                                                            #The purpose of this function is
    if (is.vector(y)) y <- t(as.matrix(y))                                  #to get MW statistics
                                                                            
    m <- length(x)                                                          
    n <- dim(y)[2]                                                          
    t <- dim(y)[1]                                                          
                                                                            
    X <- x                                                                  
    Y <- y                                                                  
                                                                            
    T <- lapply(                                                            
            1:t,                                                            
            function(x) {
                wilcox.test(X, Y[x, ], alternative = 'greater')$statistic   #Obtain WM by willcox.test
            }
    )
    
    unlist(T)

}


EWMA.WM.Chart <- function(x, y, lambda = 0.1, L = 3, tt = NULL, graph = TRUE, 
Z0 = ifelse(is.vector(y), dim(t(as.matrix(y)))[2] * (length(x) + dim(t(as.matrix(y)))[2] + 1) / 2, dim(y)[2] * (length(x) + dim(y)[2] + 1) / 2)){
                                                                            #The purpose of this function is
    if (is.vector(y)) y <- t(as.matrix(y))                                  #to draw EWMA WM chart
                                                                            #
    m <- length(x)                                                          #
    n <- dim(y)[2]                                                          #
    if (is.null(tt)){                                                       #
                                                                            #
        t <- c(1:dim(y)[1])                                                 #
                                                                            #
    } else {                                                                #
                                                                            #
        t <- rep(tt, dim(y)[1])
    
    }
    
    X <- WM.f(x, y)                                                         #Get WM statitics
    
    Z <- Z.f(X, lambda, Z0 = Z0)                                            #Get charting statistics
    
    CL <- n * (m + n + 1) / 2                                               #Get CL
    sigma2.W <- m * n * (m + n + 1) / 12
    V <- sqrt(sigma2.W * (lambda / (2 - lambda) * (1 - (1 - lambda) ^ (2 * t))))
    
    UCL <- CL + L * V                                                       #Get UCL
     
    LCL <- CL - L * V                                                       #Get LCL
    
    if (graph == TRUE){
        plot(c(1, dim(y)[1]), c(3 / 4 * min(c(Z, LCL)), 5 / 4 * max(c(Z, UCL))), type = 'n', xlab = 't', ylab = 'Z')
        points(c(1:dim(y)[1]), UCL, type = 'l', lty = 2)
        points(c(1:dim(y)[1]), rep(CL, dim(y)[1]), type = 'l')
        points(c(1:dim(y)[1]), LCL, type = 'l', lty = 2)
        points(c(1:dim(y)[1]), Z, type = 'l')
    }
    return(list(CL = CL, LCL = LCL, UCL = UCL, Z = Z))
    
}

#x <- rnorm(100)
#y <- matrix(rnorm(10000), ncol = 5)

#debug(EWMA.WM.Chart)

#EWMA.WM.Chart(x, y, lambda = 0.1, L = 3)


##########################################################################################

#get.data <- function(n, subgroup.num = 10000, type = c('norm', 't', 'exp', 'chisq'), shift = 0){
#                                                                                                        #Generate data
#    if (type == 'norm'){                                                                                #
#                                                                                                        #
#        matrix(rnorm(n * subgroup.num), ncol = n, nrow = subgroup.num) + shift                          #
#                                                                                                        #Get normal(0,1) + shift
#    } else if (type == 't') {                                                                           #
#                                                                                                        #
#        matrix(rt(n * subgroup.num, 5), ncol = n, nrow = subgroup.num) + shift                          #
#                                                                                                        #Get t(1) + shift
#    } else if (type == 'exp'){                                                                          #
#                                                                                                        #
#        matrix(rexp(n * subgroup.num, 1) - log(2), ncol = n, nrow = subgroup.num) + shift               #
#                                                                                                        #Get exp(1) - median + shift
#    } else if (type == 'chisq') {                                                                       #median = log(2)
#                                                                                                        #
#        matrix(rchisq(n * subgroup.num, 1) - (1 - 2 / 9) ^ 3, ncol = n, nrow = subgroup.num) + shift    #
#                                                                                                        #Get chisq(1) - median + shift
#    }                                                                                                   #median = (1 - 2 / 9) ^ 3
#                                                                                                        #
#                                                                                                        #
#}                                                                                                       #
#                                                                                                        #
#
#
#RL.EX <- function(m, n, lambda, percentile, L, tt = NULL, subgroup.num = 10000, x.type = 'norm', y.type = 'norm', shift = 0, burn = tt){ 
#    rl <- 0                                                                                             #The purpose of this function is
#    Z0 <- 0                                                                                             #obtain RL
#                                                                                                        #Here burn means data points before burn
#    iter <- 0                                                                                           #will be dropped
#                                                                                                        #
#    repeat {                                                                                            #
#                                                                                                        #
#        iter <- iter + 1                                                                                #
#                                                                                                        #
#        x <- get.data(m, 1, x.type)
#        y <- get.data(n, subgroup.num, y.type, shift)
#
#        EX <- EWMA.EX.Chart(x, y, percentile, lambda, L, tt, graph = FALSE, Z0)
#
#        if (iter == 1 && length(burn) > 0){
#        
#            Z <- EX$Z[-c(1:burn)]
#            LCL <- EX$LCL[-c(1:burn)]
#            UCL <- EX$UCL[-c(1:burn)]
#        
#        } else {
#            Z <- EX$Z
#            LCL <- EX$LCL
#            UCL <- EX$UCL
#        }
#        
#
#        
#        if (any(Z - UCL > 0) || any(Z - LCL < 0)){
#        
#            rl <- rl + min(c(which(Z - UCL > 0), which(Z - LCL < 0)))
#            
#            break
#        
#        } else {
#        
#            rl <- ifelse(iter == 1 && length(burn) > 0, rl + subgroup.num - burn, rl + subgroup.num)
#            Z0 <- tail(Z, 1)
#        
#        }
#        
#    } 
#    return(rl)
#}
#
#RL.WM <- function(m, n, lambda, L, tt = NULL, subgroup.num = 10000, x.type = 'norm', y.type = 'norm', shift = 0, burn = tt){ 
#    rl <- 0                                                                                            #The purpose of this function is
#    Z0 <- n * (m + n + 1) / 2                                                                          #obtain RL
#                                                                                                       #Here burn means data points before burn
#    iter <- 0                                                                                          #will be dropped
#                                                                                                       #
#    repeat {                                                                                           #
#        
#        iter <- iter + 1
#    
#    
#        x <- get.data(m, 1, x.type)
#        y <- get.data(n, subgroup.num, y.type, shift)
#        WM <- EWMA.WM.Chart(x, y, lambda, L, tt, graph = FALSE, Z0)
#
#        if (iter == 1 && length(burn) > 0){
#        
#            Z <- WM$Z[-c(1:burn)]
#            LCL <- WM$LCL[-c(1:burn)]
#            UCL <- WM$UCL[-c(1:burn)]
#        
#        } else {
#            Z <- WM$Z
#            LCL <- WM$LCL
#            UCL <- WM$UCL
#        }
#        
#        
#        if (any(Z - UCL > 0) || any(Z - LCL < 0)){
#        
#            rl <- rl + min(c(which(Z - UCL > 0), which(Z - LCL < 0)))
#            
#            break
#        
#        } else {
#        
#            rl <- ifelse(iter == 1 && length(burn) > 0, rl + subgroup.num - burn, rl + subgroup.num)
#            Z0 <- tail(Z, 1)
#        
#        }
#        
#    }
#    return(rl)
#}
#
#
##debug(RL.WM)
#
#ARL.EX <- function(m, n, lambda, percentile, L, tt = NULL, subgroup.num = 10000, x.type = 'norm', y.type = 'norm', shift = 0, burn = tt, maxsim = 100){
#                                                                                                          #The purpose of this function is
#    ll <- rep(NA, maxsim)                                                                                 #to get RL distributuion and its statistics
#                                                                                                          #
#    for (i in 1:maxsim){                                                                                  #
#                                                                                                          #
#        ll[i] <- RL.EX(m, n, lambda, percentile, L, tt, subgroup.num, x.type, y.type, shift, burn)        #
#                                                                                                          #
#    }                                                                                                     #
#                                                                                                          #
#    return(list(RL = ll, ARL = mean(ll), SDRL = sd(ll), Qs = quantile(ll, c(0.05, 0.25, 0.5, 0.75, 0.95))))     #
#
#}
#
#ARL.WM <- function(m, n, lambda, L, tt = NULL, subgroup.num = 10000, x.type = 'norm', y.type = 'norm', shift = 0, burn = tt, maxsim = 100){
#                                                                                                          #The purpose of this function is
#    ll <- rep(NA, maxsim)                                                                                 #to get RL distributuion and its statistics
#                                                                                                          #
#    for (i in 1:maxsim){                                                                                  #
#                                                                                                          #
#        ll[i] <- RL.WM(m, n, lambda, L, tt, subgroup.num, x.type, y.type, shift, burn)                    #
#                                                                                                          #
#    }                                                                                                     #
#                                                                                                          #
#    return(list(RL = ll, ARL = mean(ll), SDRL = sd(ll), Qs = quantile(ll, c(0.05, 0.25, 0.5, 0.75, 0.95))))     #
#
#}
#

#EWMA.EX <- ARL.EX(100, 20, 0.1, 0.5, 1.5, maxsim = 1000)
#save(EWMA.EX, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/EX.Rdata')
#
#EWMA.WM <- ARL.WM(100, 20, 0.1, 2.7, maxsim = 1000)
#save(EWMA.WM, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/WM.Rdata')





#m = 100
#n = 5
#
#CL = n * (m + n + 1) / 2 
#
#UCL = 310.5
#sigma2 = m * n * (m + n + 1) / 12
#lambda = 0.1
#L = (UCL - CL) / sqrt(sigma2 * lambda / (2 - lambda))
#
#debug(RL.WM)

#EWMA.WM <- ARL.WM(100, 5, 0.1, 2.984289, tt = 1000000, maxsim = 10000, burn = NULL)
#save(EWMA.WM, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/WM.Rdata')



##########################################################################################
    #Example
##########################################################################################

#X <- c(
#74.030,	74.002,	74.019,	73.992,	74.008,
#73.995,	73.992,	74.001,	74.011,	74.004,
#73.998,	74.024,	74.021,	74.005,	74.002,
#74.002,	73.996,	73.993,	74.015,	74.009,
#73.992,	74.007,	74.015,	73.989,	74.014,
#74.009,	73.994,	73.997,	73.985,	73.993,
#73.995,	74.006,	73.994,	74.000,	74.005,
#73.985,	74.003,	73.993,	74.015,	73.988,
#74.008,	73.995,	74.009,	74.005,	74.004,
#73.998,	74.000,	73.990,	74.007,	73.995,
#73.994,	73.998,	73.994,	73.995,	73.990,
#74.004,	74.000,	74.007,	74.000,	73.996,
#73.983,	74.002,	73.998,	73.997,	74.012,
#74.006,	73.967,	73.994,	74.000,	73.984,
#74.012,	74.014,	73.998,	73.999,	74.007,
#74.000,	73.984,	74.005,	73.998,	73.996,
#73.994,	74.012,	73.986,	74.005,	74.007,
#74.006,	74.010,	74.018,	74.003,	74.000,
#73.984,	74.002,	74.003,	74.005,	73.997,
#74.000,	74.010,	74.013,	74.020,	74.003,
#73.982,	74.001,	74.015,	74.005,	73.996,
#74.004,	73.999,	73.990,	74.006,	74.009,
#74.010,	73.989,	73.990,	74.009,	74.014,
#74.015,	74.008,	73.993,	74.000,	74.010,
#73.982,	73.984,	73.995,	74.017,	74.013
#)
#
#Y <- matrix(c(
#74.012,	74.015,	74.030,	73.986,	74.000,
#73.995,	74.010,	73.990,	74.015,	74.001,
#73.987,	73.999,	73.985,	74.000,	73.990,
#74.008,	74.010,	74.003,	73.991,	74.006,
#74.003,	74.000,	74.001,	73.986,	73.997,
#73.994,	74.003,	74.015,	74.020,	74.004,
#74.008,	74.002,	74.018,	73.995,	74.005,
#74.001,	74.004,	73.990,	73.996,	73.998,
#74.015,	74.000,	74.016,	74.025,	74.000,
#74.030,	74.005,	74.000,	74.016,	74.012,
#74.001,	73.990,	73.995,	74.010,	74.024,
#74.015,	74.020,	74.024,	74.005,	74.019,
#74.035,	74.010,	74.012,	74.015,	74.026,
#74.017,	74.013,	74.036,	74.025,	74.026,
#74.010,	74.005,	74.029,	74.000,	74.020
#
#
#), ncol = 5, byrow = TRUE)
#
#debug(EX.f)
#
#EWMA.EX.Chart(X, Y, 0.5, 0.1, 2.31, NULL, graph = TRUE, 2.5)
#
#
#
#
#
#
##########################################################################################












##########################################################################################
