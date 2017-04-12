##########################################################################
    #No extra package is needed
##########################################################################

Z.f <- function(x, lambda, Z0 = 0){                                #Basic structure of EWMA chart
                                                                   #Z.i = lambda * X + (1 - lambda) * Z.i-1
    tt <- length(x)                                                #
                                                                   #
    z <- rep(Z0, tt + 1)                                           #
                                                                   #
    for (i in 2:(tt + 1)){                                         #
                                                                   #
        z[i] <- lambda * x[i - 1] + (1 - lambda) * z[i - 1]        #
                                                                   #
    }                                                              #
                                                                   #
    z[-1]                                                          #
                                                                   #
}                                                                  #


##########################################################################

EWMA.N.Chart <- function(x, med = 0, lambda = 0.1, L = 3, graph = TRUE){        #Get EWMA chart based on normal distribution
                                                                                #
    t <- 1:dim(x)[1]                                                            #
    n <- dim(x)[2]                                                              #
                                                                                #
    X <- rowMeans(x)                                                            #Get x.bar.i
                                                                                #
    Z <- Z.f(X, lambda)                                                         #Get Z based on x.bar.i
                                                                                #
    UCL <- L * 1 * sqrt(lambda / (2 - lambda) * (1 - (1 - lambda) ^ (2 * t)))   #Get limits
    LCL <- -UCL                                                                 #
    
    if (graph == TRUE){                                                         #When graph is TRUE, it will output a chart visually
        plot(c(1, max(t)), c(5 / 4 * min(c(Z, LCL)), 5 / 4 * max(c(Z, UCL))), type = 'n', xlab = 't', ylab = 'Z')
        points(t, UCL, type = 'l', lty = 2)
        points(t, rep(0, max(t)), type = 'l')
        points(t, LCL, type = 'l', lty = 2)
        points(t, Z, type = 'l')
    }
    return(list(CL = 0, LCL = LCL, UCL = UCL, Z = Z))                           
    
}



##########################################################################

SN.f <- function(x, med){                                                  #Get SN
                                                                           #
    n <- dim(x)[2]                                                         #
                                                                           #
    2 * rowSums(x > med) - n                                               #
                                                                           #
}                                                                          #

EWMA.SN.Chart <- function(x, med, lambda = 0.1, L = 3, graph = TRUE){      #Get EWMA chart based on SN
                                                                           #
    t <- 1:dim(x)[1]                                                       #
    n <- dim(x)[2]                                                         #
                                                                           #
    SN <- SN.f(x, med)                                                     #Get SN
                                                                           #
    Z <- Z.f(SN, lambda)                                                   #Get Z based on SN
    
    UCL <- L * sqrt(n) * sqrt(lambda / (2 - lambda) * (1 - (1 - lambda) ^ (2 * t))) #Get limits
    LCL <- -UCL                                                                     #
    
    if (graph == TRUE){                                                    #When graph is TRUE, it will output a chart visually
        plot(c(1, max(t)), c(5 / 4 * min(c(Z, LCL)), 5 / 4 * max(c(Z, UCL))), type = 'n', xlab = 't', ylab = 'Z')
        points(t, UCL, type = 'l', lty = 2)
        points(t, rep(0, max(t)), type = 'l')
        points(t, LCL, type = 'l', lty = 2)
        points(t, Z, type = 'l')
    }
    return(list(CL = 0, LCL = LCL, UCL = UCL, Z = Z))
    
}



##########################################################################

SR.f <- function(x, med){                                                      #Get SR
                                                                               #
    X <- x                                                                     #
                                                                               #
    t <- dim(x)[1]                                                             #
    n <- dim(x)[2]                                                             #
                                                                               #
    T <- lapply(                                                               #
            1:t,                                                               #
            function(x) {                                                      #
                wilcox.test(X[x, ], mu = med, alternative = 'greater')$statistic  #wilcox.test will give us the W stat
            }
    )

    
    2 * unlist(T) - n * (n + 1) / 2

}


EWMA.SR.chart <- function(x, med, lambda, L, graph = TRUE){                     #Get EWMA chart based on SR
                                                                                #
    t <- 1:dim(x)[1]                                                            #
    n <- dim(x)[2]                                                              #
                                                                                #
    SR <- SR.f(x, med)                                                          #Get SR
                                                                                #
    Z <- Z.f(SR, lambda)                                                        #Get Z based on SR
                                                                                #
    UCL <- L * sqrt(n * (n + 1) * (2 * n + 1) / 6) * sqrt(lambda / (2 - lambda) * (1 - (1 - lambda) ^ (2 * t))) #Get limits
    LCL <- -UCL                                                                                                 #
    
    if (graph == TRUE){                                                         #When graph is TRUE, it will output a chart visually
        plot(c(1, max(t)), c(5 / 4 * min(c(Z, LCL)), 5 / 4 * max(c(Z, UCL))), type = 'n', xlab = 't', ylab = 'Z')
        points(t, UCL, type = 'l', lty = 2)
        points(t, rep(0, max(t)), type = 'l')
        points(t, LCL, type = 'l', lty = 2)
        points(t, Z, type = 'l')
    }
    return(list(CL = 0, LCL = LCL, UCL = UCL, Z = Z))
}

##########################################################################
#
#get.data <- function(n, subgroup.num = 10000, type = c('norm', 't', 'exp', 'chisq'), shift = 0){        #Generate data
#                                                                                                        #
#    if (type == 'norm'){                                                                                #
#                                                                                                        #
#        matrix(rnorm(n * subgroup.num), ncol = n, nrow = subgroup.num) + shift                          #Get normal(0,1) + shift
#                                                                                                        #
#    } else if (type == 't') {                                                                           #
#                                                                                                        #
#        matrix(rt(n * subgroup.num, 5), ncol = n, nrow = subgroup.num) + shift                          #Get t(1) + shift
#                                                                                                        #
#    } else if (type == 'exp'){                                                                          #
#                                                                                                        #
#        matrix(rexp(n * subgroup.num, 1) - log(2), ncol = n, nrow = subgroup.num) + shift               #Get exp(1) - median + shift
#                                                                                                        #median = log(2)
#    } else if (type == 'chisq') {                                                                       #
#                                                                                                        #
#        matrix(rchisq(n * subgroup.num, 1) - (1 - 2 / 9) ^ 3, ncol = n, nrow = subgroup.num) + shift    #Get chisq(1) - median + shift
#                                                                                                        #median = (1 - 2 / 9) ^ 3
#    }                                                                                                   #
#                                                                                                        #
#                                                                                                        #
#}                                                                                                       #
#
#RL.SN <- function(n, med, lambda, UCL, LCL, subgroup.num = 10000, type = 'norm', shift = 0){            #Get Run-length distribution by simulations
#    rl <- 0                                                                                             #it is based on EWMA SN
#    Z0 <- 0                                                                                             #
#    repeat {                                                                                            #
#        x <- get.data(n, subgroup.num, type, shift)                                                     #This code will calculate charting statistics
#        SN <- SN.f(x, med)                                                                              #subgroup.num times in the meanwhile
#        Z <- Z.f(SN, lambda, Z0)                                                                        #
#        if (any(Z > UCL) || any(Z < LCL)){                                                              #When one of them is out of control
#                                                                                                        #it will find the first out-of-control point
#            rl <- rl + min(c(which(Z > UCL), which(Z < LCL)))                                           #and run-length will be the run-length at the last iteration
#                                                                                                        #plus the length of the first out-of-control point
#            break                                                                                       #
#                                                                                                        #
#        } else {                                                                                        #When none of them is out of control,
#                                                                                                        #it will only calculate run-length
#            rl <- rl + subgroup.num                                                                     #run-length will be the run-length at the last iteration
#            Z0 <- tail(Z, 1)                                                                            #plus subgroup.num
#                                                                                                        #
#        }                                                                                               #
#                                                                                                        #
#    }                                                                                                   #
#    return(rl)                                                                                          #
#}                                                                                                       #
#
#RL.SR <- function(n, med, lambda, UCL, LCL, subgroup.num = 10000, type = 'norm', shift = 0){            #Get Run-length distribution by simulations
#    rl <- 0                                                                                             #it is based on EWMA SR
#    Z0 <- 0                                                                                             #
#    repeat {                                                                                            #
#        x <- get.data(n, subgroup.num, type, shift)                                                     #This code will calculate charting statistics
#        SR <- SR.f(x, med)                                                                              #subgroup.num times in the meanwhile
#        Z <- Z.f(SR, lambda, Z0)                                                                        #
#        if (any(Z > UCL) || any(Z < LCL)){                                                              #When one of them is out of control
#                                                                                                        #it will find the first out-of-control point
#            rl <- rl + min(c(which(Z > UCL), which(Z < LCL)))                                           #and run-length will be the run-length at the last iteration
#                                                                                                        #plus the length of the first out-of-control point
#            break                                                                                       #
#                                                                                                        #
#        } else {                                                                                        #When none of them is out of control,
#                                                                                                        #it will only calculate run-length
#            rl <- rl + subgroup.num                                                                     #run-length will be the run-length at the last iteration
#            Z0 <- tail(Z, 1)                                                                            #plus subgroup.num
#                                                                                                        #
#        }                                                                                               #
#                                                                                                        #
#    }                                                                                                   #
#    return(rl)                                                                                          #
#}                                                                                                       #
#
#ARL.SN <- function(t, n, med, lambda, L, subgroup.num = 10000, type = 'norm', shift = 0, maxiter = 100){#Get stat from run-length distribution
#    ll <- rep(0, maxiter)                                                                               #based on EWMA SN
#    UCL <- L * sqrt(n) * sqrt(lambda / (2 - lambda) * (1 - (1 - lambda) ^ (2 * t)))                     #
#    LCL <- -UCL                                                                                         #
#    for (i in 1:maxiter){                                                                               #
#        ll[i] <- RL.SN(n, med, lambda, UCL, LCL, subgroup.num = subgroup.num, type = type, shift)       #Record RL
#    }                                                                                                   #
#    return(list(RL=ll, ARL=mean(ll), SDRL=sd(ll), med. = median(ll)))                                   #
#}                                                                                                       #
#
#ARL.SR <- function(t, n, med, lambda, L, subgroup.num = 10000, type = 'norm', shift = 0, maxiter = 100){         #Get stat from run-length distribution
#    ll <- rep(0, maxiter)                                                                                        #based on EWMA SR
#    UCL <- L * sqrt(n * (n + 1) * (2 * n + 1) / 6) * sqrt(lambda / (2 - lambda) * (1 - (1 - lambda) ^ (2 * t)))  #
#    LCL <- -UCL                                                                                                  #
#    for (i in 1:maxiter){                                                                                        #
#        ll[i] <- RL.SR(n, med, lambda, UCL, LCL, subgroup.num = subgroup.num, type = type, shift)                #Record RL
#    }                                                                                                            #
#    return(list(RL=ll, ARL=mean(ll), SDRL=sd(ll), med. = median(ll)))                                            #
#}                                                                                                                #
#
###########################################################################
#    #simulate run-length distribution
###########################################################################
#
#tt <- 1000000   #This will give the run-length in steady-state
#n <- 20
#med <- 0
#
#
#EWMA.SN.N <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000)
#EWMA.SN.T <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000)
#EWMA.SN.EXP <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000)
#EWMA.SN.CHISQ <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000)
#
#save(EWMA.SN.N, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.Rdata')
#save(EWMA.SN.T, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.Rdata')
#save(EWMA.SN.EXP, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.Rdata')
#save(EWMA.SN.CHISQ, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.Rdata')
#
#EWMA.SR.N <- ARL.SR(tt, n, 0, 0.1, 2.7, subgroup.num = 1000, maxiter = 10000)
#EWMA.SR.T <- ARL.SR(tt, n, 0, 0.1, 2.7, subgroup.num = 1000, type = 't', maxiter = 10000)
#EWMA.SR.EXP <- ARL.SR(tt, n, 0, 0.1, 2.7, subgroup.num = 1000, type = 'exp', maxiter = 10000)
#EWMA.SR.CHISQ <- ARL.SR(tt, n, 0, 0.1, 2.7, subgroup.num = 1000, type = 'chisq', maxiter = 10000)
#
#save(EWMA.SR.N, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.Rdata')
#save(EWMA.SR.T, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.Rdata')
#save(EWMA.SR.EXP, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.Rdata')
#save(EWMA.SR.CHISQ, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.Rdata')
#
#EWMA.SN.N.1 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000, shift = 0.1)
#EWMA.SN.T.1 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000, shift = 0.1)
#EWMA.SN.EXP.1 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000, shift = 0.1)
#EWMA.SN.CHISQ.1 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000, shift = 0.1)
#
#EWMA.SN.N.2 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000, shift = 0.25)
#EWMA.SN.T.2 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000, shift = 0.25)
#EWMA.SN.EXP.2<- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000, shift = 0.25)
#EWMA.SN.CHISQ.2 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000, shift = 0.25)
#
#EWMA.SN.N.3 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000, shift = 0.5)
#EWMA.SN.T.3 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000, shift = 0.5)
#EWMA.SN.EXP.3 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000, shift = 0.5)
#EWMA.SN.CHISQ.3 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000, shift = 0.5)
#
#EWMA.SN.N.4 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000, shift = 0.75)
#EWMA.SN.T.4 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000, shift = 0.75)
#EWMA.SN.EXP.4 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000, shift = 0.75)
#EWMA.SN.CHISQ.4 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000, shift = 0.75)
#
#EWMA.SN.N.5 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000, shift = 1)
#EWMA.SN.T.5 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000, shift = 1)
#EWMA.SN.EXP.5 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000, shift = 1)
#EWMA.SN.CHISQ.5 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000, shift = 1)
#
#save(EWMA.SN.N.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.1.Rdata')
#save(EWMA.SN.T.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.1.Rdata')
#save(EWMA.SN.EXP.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.1.Rdata')
#save(EWMA.SN.CHISQ.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.1.Rdata')
#
#save(EWMA.SN.N.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.2.Rdata')
#save(EWMA.SN.T.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.2.Rdata')
#save(EWMA.SN.EXP.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.2.Rdata')
#save(EWMA.SN.CHISQ.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.2.Rdata')
#
#save(EWMA.SN.N.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.3.Rdata')
#save(EWMA.SN.T.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.3.Rdata')
#save(EWMA.SN.EXP.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.3.Rdata')
#save(EWMA.SN.CHISQ.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.3.Rdata')
#
#save(EWMA.SN.N.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.4.Rdata')
#save(EWMA.SN.T.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.4.Rdata')
#save(EWMA.SN.EXP.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.4.Rdata')
#save(EWMA.SN.CHISQ.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.4.Rdata')
#
#save(EWMA.SN.N.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.5.Rdata')
#save(EWMA.SN.T.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.5.Rdata')
#save(EWMA.SN.EXP.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.5.Rdata')
#save(EWMA.SN.CHISQ.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.5.Rdata')
#
#EWMA.SR.N.1 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000, shift = 0.1)
#EWMA.SR.T.1 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000, shift = 0.1)
#EWMA.SR.EXP.1 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000, shift = 0.1)
#EWMA.SR.CHISQ.1 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000, shift = 0.1)
#
#EWMA.SR.N.2 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000, shift = 0.25)
#EWMA.SR.T.2 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000, shift = 0.25)
#EWMA.SR.EXP.2<- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000, shift = 0.25)
#EWMA.SR.CHISQ.2 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000, shift = 0.25)
#
#EWMA.SR.N.3 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000, shift = 0.5)
#EWMA.SR.T.3 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000, shift = 0.5)
#EWMA.SR.EXP.3 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000, shift = 0.5)
#EWMA.SR.CHISQ.3 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000, shift = 0.5)
#
#EWMA.SR.N.4 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000, shift = 0.75)
#EWMA.SR.T.4 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000, shift = 0.75)
#EWMA.SR.EXP.4 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000, shift = 0.75)
#EWMA.SR.CHISQ.4 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000, shift = 0.75)
#
#EWMA.SR.N.5 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000, shift = 1)
#EWMA.SR.T.5 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000, shift = 1)
#EWMA.SR.EXP.5 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000, shift = 1)
#EWMA.SR.CHISQ.5 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000, shift = 1)
#
#save(EWMA.SR.N.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.1.Rdata')
#save(EWMA.SR.T.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.1.Rdata')
#save(EWMA.SR.EXP.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.1.Rdata')
#save(EWMA.SR.CHISQ.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.1.Rdata')
#
#save(EWMA.SR.N.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.2.Rdata')
#save(EWMA.SR.T.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.2.Rdata')
#save(EWMA.SR.EXP.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.2.Rdata')
#save(EWMA.SR.CHISQ.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.2.Rdata')
#
#save(EWMA.SR.N.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.3.Rdata')
#save(EWMA.SR.T.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.3.Rdata')
#save(EWMA.SR.EXP.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.3.Rdata')
#save(EWMA.SR.CHISQ.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.3.Rdata')
#
#save(EWMA.SR.N.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.4.Rdata')
#save(EWMA.SR.T.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.4.Rdata')
#save(EWMA.SR.EXP.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.4.Rdata')
#save(EWMA.SR.CHISQ.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.4.Rdata')
#
#save(EWMA.SR.N.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.5.Rdata')
#save(EWMA.SR.T.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.5.Rdata')
#save(EWMA.SR.EXP.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.5.Rdata')
#save(EWMA.SR.CHISQ.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.5.Rdata')
#
#EWMA.SN.N.6 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000, shift = 0.05)
#EWMA.SN.T.6 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000, shift = 0.05)
#EWMA.SN.EXP.6 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000, shift = 0.05)
#EWMA.SN.CHISQ.6 <- ARL.SN(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000, shift = 0.05)
#
#save(EWMA.SN.N.6, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.6.Rdata')
#save(EWMA.SN.T.6, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.6.Rdata')
#save(EWMA.SN.EXP.6, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.6.Rdata')
#save(EWMA.SN.CHISQ.6, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.6.Rdata')
#
#
#
#
#EWMA.SR.N.6 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, maxiter = 10000, shift = 0.05)
#EWMA.SR.T.6 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 't', maxiter = 10000, shift = 0.05)
#EWMA.SR.EXP.6 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'exp', maxiter = 10000, shift = 0.05)
#EWMA.SR.CHISQ.6 <- ARL.SR(tt, n, 0, 0.1, 2.71, subgroup.num = 1000, type = 'chisq', maxiter = 10000, shift = 0.05)
#
#save(EWMA.SR.N.6, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.6.Rdata')
#save(EWMA.SR.T.6, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.6.Rdata')
#save(EWMA.SR.EXP.6, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.6.Rdata')
#save(EWMA.SR.CHISQ.6, file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.6.Rdata')
#
###########################################################################
#    #get other quantiles
###########################################################################
#
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.1.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.2.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.3.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.4.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.5.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.N.6.Rdata')
#
#quantile(EWMA.SN.N$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.N.1$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.N.2$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.N.3$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.N.4$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.N.5$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.N.6$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#
#
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.1.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.2.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.3.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.4.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.5.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.T.6.Rdata')
#
#quantile(EWMA.SN.T$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.T.1$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.T.2$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.T.3$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.T.4$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.T.5$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.T.6$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#
#
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.1.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.2.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.3.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.4.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.5.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.EXP.6.Rdata')
#
#quantile(EWMA.SN.EXP$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.EXP.1$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.EXP.2$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.EXP.3$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.EXP.4$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.EXP.5$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.EXP.6$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.1.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.2.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.3.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.4.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.5.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SN.CHISQ.6.Rdata')
#
#quantile(EWMA.SN.CHISQ$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.CHISQ.1$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.CHISQ.2$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.CHISQ.3$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.CHISQ.4$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.CHISQ.5$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SN.CHISQ.6$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#
#
#
#
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.1.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.2.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.3.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.4.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.5.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.N.6.Rdata')
#
#quantile(EWMA.SR.N$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.N.1$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.N.2$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.N.3$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.N.4$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.N.5$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.N.6$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#
#
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.1.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.2.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.3.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.4.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.5.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.T.6.Rdata')
#
#quantile(EWMA.SR.T$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.T.1$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.T.2$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.T.3$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.T.4$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.T.5$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.T.6$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#
#
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.1.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.2.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.3.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.4.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.5.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.EXP.6.Rdata')
#
#quantile(EWMA.SR.EXP$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.EXP.1$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.EXP.2$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.EXP.3$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.EXP.4$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.EXP.5$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.EXP.6$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.1.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.2.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.3.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.4.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.5.Rdata')
#load(file = 'C:/Users/bolus/Desktop/st697 chakraborti/nonparametric/ewma/SR.CHISQ.6.Rdata')
#
#quantile(EWMA.SR.CHISQ$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.CHISQ.1$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.CHISQ.2$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.CHISQ.3$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.CHISQ.4$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.CHISQ.5$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#quantile(EWMA.SR.CHISQ.6$RL, c(0.05, 0.25, 0.5, 0.75, 0.95))
#
