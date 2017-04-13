####################################################################################################################################################
    #For running this code, you may need to run FP.R first 
####################################################################################################################################################

source(file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.R') #This means I run FP.R. You can put FP.R at any folder.
                                                   
####################################################################################################################################################

get.data <- function(n, subgroup.num = 10000, type = c('norm', 't', 'exp', 'chisq'), shift = 0){
                                                                             #The purpose of this function is 
                                                                             #to generate different samples.
                                                                             #All samples have 0 mediam
                                                                             #shift is to decide their locations
    if (type == 'norm'){
    
        matrix(rnorm(n * subgroup.num), ncol = n, nrow = subgroup.num) + shift
    
    } else if (type == 't') {
    
        matrix(rt(n * subgroup.num, 5), ncol = n, nrow = subgroup.num) + shift
    
    } else if (type == 'exp'){
    
        matrix(rexp(n * subgroup.num, 1) - log(2), ncol = n, nrow = subgroup.num) + shift
    
    } else if (type == 'chisq') {
    
        matrix(rchisq(n * subgroup.num, 1) - (1 - 2 / 9) ^ 3, ncol = n, nrow = subgroup.num) + shift
    
    }


}


Uncond.RL <- function(m, n, shift, ARL = 370, chart = 't', xtype = 'norm', ytype = 'norm', subgroup.num = 1000, U.modified.pmf = NULL, cores = 1, maxiter = 10000){ 
                                                                             #The purpose of this function is 
                                                                             #to get run lengths for different charts
                                                                             #
    alpha <- 1/ ARL
    rl <- 0   
    
    if (chart == 'FP' && is.null(U.modified.pmf)){
        U.modified.pmf <- get.U.Modi.Dist(m, n, cores = cores, maxiter = maxiter)
    }
    
    repeat { 
        X <- get.data(m, subgroup.num, xtype)
        Y <- get.data(n, subgroup.num, ytype, shift)
        
        
        if (chart == 't'){
        
            p.value <- unlist(lapply(
                        1:subgroup.num, 
                        function(x) t.test(X[x, ], Y[x, ])$p.value           #Run 2-sample T test
            ))
        
        } else if (chart == 'WM'){
        
            p.value <- unlist(lapply(
                        1:subgroup.num, 
                        function(x) wilcox.test(X[x, ], Y[x, ])$p.value      #Run 2-sample WM test
            ))
        
        } else if (chart == 'FP'){
        
            p.value <- unlist(lapply(
                        1:subgroup.num, 
                        function(x) {                                        #Because I already boil down FP modification
                            U.modified.stat <- get.U.Modi.Stat(X[x, ], Y[x, ])#this part is to get p.value from FP modification
                            sum(U.modified.pmf[which(abs(U.modified.pmf[, 1]) > abs(U.modified.stat)), 2])
                        
                        }
            ))
        }
        
        if (any(p.value < alpha)){
        
            rl <- rl + min(which(p.value < alpha))
            
            break
        
        } else {
        
            rl <- rl + subgroup.num
        
        }

            
    } 
    
    return(rl)
}

Uncond.ARL <- function(m, n, shift, ARL = 370, chart = 't', xtype = 'norm', ytype = 'norm', subgroup.num = 1000, cores = 1, maxiter = 10000, maxsim = 100){ 
                                                                                 #The purpose of this function is
                                                                                 #to get the unconditional run-length distribution
                                                                                 #for different charts by simulations
                                                                                 #
    RL <- rep(NA, maxsim)
    
    U.modified.pmf <- NULL                                                        #If you do not have a specific
                                                                                  #the U modified distribution
    if (chart == 'FP'){                                                           #it will get it
        U.modified.pmf <- get.U.Modi.Dist(m, n, cores = cores, maxiter = maxiter) #by simulations
    }                                                                             #
    
    for (sim in 1:maxsim){
    
        RL[sim] <- Uncond.RL(m, n, shift, ARL, chart, xtype, ytype, subgroup.num, U.modified.pmf, cores, maxiter)
    
    }

    list(RL = RL, ARL = mean(RL), SDRL = sd(RL), Qs = quantile(RL, c(0.025, 0.25, 0.5, 0.75, 0.975)))

}

T.2samp.0 <- Uncond.ARL(100, 20, 0)
T.2samp.1 <- Uncond.ARL(100, 20, 0.1)
T.2samp.2 <- Uncond.ARL(100, 20, 0.25)
T.2samp.3 <- Uncond.ARL(100, 20, 0.5)
T.2samp.4 <- Uncond.ARL(100, 20, 0.75)
T.2samp.5 <- Uncond.ARL(100, 20, 1)

save(T.2samp.0, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.0.Rdata')
save(T.2samp.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.1.Rdata')
save(T.2samp.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.2.Rdata')
save(T.2samp.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.3.Rdata')
save(T.2samp.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.4.Rdata')
save(T.2samp.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.5.Rdata')

WM.2samp.0 <- Uncond.ARL(100, 20, 0, chart = 'WM')
WM.2samp.1 <- Uncond.ARL(100, 20, 0.1, chart = 'WM')
WM.2samp.2 <- Uncond.ARL(100, 20, 0.25, chart = 'WM')
WM.2samp.3 <- Uncond.ARL(100, 20, 0.5, chart = 'WM')
WM.2samp.4 <- Uncond.ARL(100, 20, 0.75, chart = 'WM')
WM.2samp.5 <- Uncond.ARL(100, 20, 1, chart = 'WM')

save(WM.2samp.0, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.0.Rdata')
save(WM.2samp.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.1.Rdata')
save(WM.2samp.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.2.Rdata')
save(WM.2samp.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.3.Rdata')
save(WM.2samp.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.4.Rdata')
save(WM.2samp.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.5.Rdata')


FP.2samp.0 <- Uncond.ARL(100, 20, 0, chart = 'FP')
FP.2samp.1 <- Uncond.ARL(100, 20, 0.1, chart = 'FP')
FP.2samp.2 <- Uncond.ARL(100, 20, 0.25, chart = 'FP')
FP.2samp.3 <- Uncond.ARL(100, 20, 0.5, chart = 'FP')
FP.2samp.4 <- Uncond.ARL(100, 20, 0.75, chart = 'FP')
FP.2samp.5 <- Uncond.ARL(100, 20, 1, chart = 'FP')

save(FP.2samp.0, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.0.Rdata')
save(FP.2samp.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.1.Rdata')
save(FP.2samp.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.2.Rdata')
save(FP.2samp.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.3.Rdata')
save(FP.2samp.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.4.Rdata')
save(FP.2samp.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.5.Rdata')







####################################################################################################################################################

T.2samp.0 <- Uncond.ARL(100, 20, 0, ytype = 't')
T.2samp.1 <- Uncond.ARL(100, 20, 0.1, ytype = 't')
T.2samp.2 <- Uncond.ARL(100, 20, 0.25, ytype = 't')
T.2samp.3 <- Uncond.ARL(100, 20, 0.5, ytype = 't')
T.2samp.4 <- Uncond.ARL(100, 20, 0.75, ytype = 't')
T.2samp.5 <- Uncond.ARL(100, 20, 1, ytype = 't')

save(T.2samp.0, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.0.y.t.Rdata')
save(T.2samp.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.1.y.t.Rdata')
save(T.2samp.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.2.y.t.Rdata')
save(T.2samp.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.3.y.t.Rdata')
save(T.2samp.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.4.y.t.Rdata')
save(T.2samp.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.5.y.t.Rdata')

WM.2samp.0 <- Uncond.ARL(100, 20, 0, chart = 'WM', ytype = 't')
WM.2samp.1 <- Uncond.ARL(100, 20, 0.1, chart = 'WM', ytype = 't')
WM.2samp.2 <- Uncond.ARL(100, 20, 0.25, chart = 'WM', ytype = 't')
WM.2samp.3 <- Uncond.ARL(100, 20, 0.5, chart = 'WM', ytype = 't')
WM.2samp.4 <- Uncond.ARL(100, 20, 0.75, chart = 'WM', ytype = 't')
WM.2samp.5 <- Uncond.ARL(100, 20, 1, chart = 'WM', ytype = 't')

save(WM.2samp.0, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.0.y.t.Rdata')
save(WM.2samp.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.1.y.t.Rdata')
save(WM.2samp.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.2.y.t.Rdata')
save(WM.2samp.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.3.y.t.Rdata')
save(WM.2samp.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.4.y.t.Rdata')
save(WM.2samp.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.5.y.t.Rdata')


FP.2samp.0 <- Uncond.ARL(100, 20, 0, chart = 'FP', ytype = 't')
FP.2samp.1 <- Uncond.ARL(100, 20, 0.1, chart = 'FP', ytype = 't')
FP.2samp.2 <- Uncond.ARL(100, 20, 0.25, chart = 'FP', ytype = 't')
FP.2samp.3 <- Uncond.ARL(100, 20, 0.5, chart = 'FP', ytype = 't')
FP.2samp.4 <- Uncond.ARL(100, 20, 0.75, chart = 'FP', ytype = 't')
FP.2samp.5 <- Uncond.ARL(100, 20, 1, chart = 'FP', ytype = 't')

save(FP.2samp.0, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.0.y.t.Rdata')
save(FP.2samp.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.1.y.t.Rdata')
save(FP.2samp.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.2.y.t.Rdata')
save(FP.2samp.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.3.y.t.Rdata')
save(FP.2samp.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.4.y.t.Rdata')
save(FP.2samp.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.5.y.t.Rdata')

####################################################################################################################################################

T.2samp.0 <- Uncond.ARL(100, 20, 0, ytype = 'exp')
T.2samp.1 <- Uncond.ARL(100, 20, 0.1, ytype = 'exp')
T.2samp.2 <- Uncond.ARL(100, 20, 0.25, ytype = 'exp')
T.2samp.3 <- Uncond.ARL(100, 20, 0.5, ytype = 'exp')
T.2samp.4 <- Uncond.ARL(100, 20, 0.75, ytype = 'exp')
T.2samp.5 <- Uncond.ARL(100, 20, 1, ytype = 'exp')

save(T.2samp.0, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.0.y.exp.Rdata')
save(T.2samp.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.1.y.exp.Rdata')
save(T.2samp.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.2.y.exp.Rdata')
save(T.2samp.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.3.y.exp.Rdata')
save(T.2samp.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.4.y.exp.Rdata')
save(T.2samp.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/T.2samp.5.y.exp.Rdata')

WM.2samp.0 <- Uncond.ARL(100, 20, 0, chart = 'WM', ytype = 'exp')
WM.2samp.1 <- Uncond.ARL(100, 20, 0.1, chart = 'WM', ytype = 'exp')
WM.2samp.2 <- Uncond.ARL(100, 20, 0.25, chart = 'WM', ytype = 'exp')
WM.2samp.3 <- Uncond.ARL(100, 20, 0.5, chart = 'WM', ytype = 'exp')
WM.2samp.4 <- Uncond.ARL(100, 20, 0.75, chart = 'WM', ytype = 'exp')
WM.2samp.5 <- Uncond.ARL(100, 20, 1, chart = 'WM', ytype = 'exp')

save(WM.2samp.0, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.0.y.exp.Rdata')
save(WM.2samp.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.1.y.exp.Rdata')
save(WM.2samp.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.2.y.exp.Rdata')
save(WM.2samp.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.3.y.exp.Rdata')
save(WM.2samp.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.4.y.exp.Rdata')
save(WM.2samp.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/WM.2samp.5.y.exp.Rdata')


FP.2samp.0 <- Uncond.ARL(100, 20, 0, chart = 'FP', ytype = 'exp')
FP.2samp.1 <- Uncond.ARL(100, 20, 0.1, chart = 'FP', ytype = 'exp')
FP.2samp.2 <- Uncond.ARL(100, 20, 0.25, chart = 'FP', ytype = 'exp')
FP.2samp.3 <- Uncond.ARL(100, 20, 0.5, chart = 'FP', ytype = 'exp')
FP.2samp.4 <- Uncond.ARL(100, 20, 0.75, chart = 'FP', ytype = 'exp')
FP.2samp.5 <- Uncond.ARL(100, 20, 1, chart = 'FP', ytype = 'exp')

save(FP.2samp.0, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.0.y.exp.Rdata')
save(FP.2samp.1, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.1.y.exp.Rdata')
save(FP.2samp.2, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.2.y.exp.Rdata')
save(FP.2samp.3, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.3.y.exp.Rdata')
save(FP.2samp.4, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.4.y.exp.Rdata')
save(FP.2samp.5, file = 'C:/Users/bolus/Desktop/st697 chakraborti/JSM/FP.2samp.5.y.exp.Rdata')

####################################################################################################################################################

####################################################################################################################################################

####################################################################################################################################################

####################################################################################################################################################

####################################################################################################################################################
