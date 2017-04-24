####################################################################################################################################################

source('https://raw.githubusercontent.com/bolus123/Statistical-Process-Control/master/FP.R')
                        #Getting the function of calculating modified U statistics
                        
#source('https://raw.githubusercontent.com/bolus123/Statistical-Process-Control/Nonparametric-Shewhart/precedence%20chart.R')
#                        #Getting the function of precedance test
                        
####################################################################################################################################################
    #Simulate data
####################################################################################################################################################
get.data <- function(
                n
                ,subgroup.amt = 1000 
                ,type = c('norm', 't', 'exp', 'chisq', 'unif', 'doubexp')
                ,shift = 0
                ,scale = 1
){
                                                                                                       #The purpose of this function is to get
	require(LaplacesDemon)                                                                             #data from different distributions with
                                                                                                       #same median which is equal to their shifts
    if (type == 'norm'){                                                                               #
                                                                                                       #
        (matrix(rnorm(n * subgroup.amt), ncol = n, nrow = subgroup.amt) + shift) / scale                         #Get data from normal(0, 1)
                                                                                                       #
    } else if (type == 't.1') {                                                                          #
                                                                                                       #
        (matrix(rt(n * subgroup.amt, 1), ncol = n, nrow = subgroup.amt) + shift) / scale                          #Get data from t(1)
                                                                                                       #
    } else if (type == 't.5' || type == 't') {                                                                          #
                                                                                                       #
        (matrix(rt(n * subgroup.amt, 5), ncol = n, nrow = subgroup.amt) + shift) / scale                          #Get data from t(5)
                                                                                                       #
    } else if (type == 'exp'){                                                                         #
                                                                                                       #
        (matrix(rexp(n * subgroup.amt, 1) - log(2), ncol = n, nrow = subgroup.amt) + shift) / scale               #Get data from exp(1) - median
                                                                                                       #
    } else if (type == 'chisq') {                                                                      #
                                                                                                       #
        (matrix(rchisq(n * subgroup.amt, 1) - (1 - 2 / 9) ^ 3, ncol = n, nrow = subgroup.amt) + shift) / scale    #Get data from chisq(1) - median
                                                                                                       #
    } else if (type == 'unif'){                                                                        #
	                                                                                                   #
        (matrix(runif(n * subgroup.amt, -1, 1), ncol = n, nrow = subgroup.amt) + shift) / scale                   #Get data from unif(-1, 1)
	                                                                                                   #
    } else if (type == 'doubexp'){                                                                     #
	                                                                                                   #
        (matrix(rlaplace(n * subgroup.amt), ncol = n, nrow = subgroup.amt) + shift) / scale                       #Get data from double exp(0, 1)
	                                                                                                   #
    }                                                                                                  #


}
####################################################################################################################################################
    #Contrust a Chart
####################################################################################################################################################
SPC.Chart <- function(CS, mu, L, V, graph = TRUE){
                                                 #The purpose of this function is to draw a SPC chart
                                                 #And this chart will always be LCL < CS < UCL
                                                 #CS, mu and V can be a vector.
                                                 #L must be a constant
    n <- length(CS)

    CL <- mu                                     #
    LCL <- mu - L * V                            #LCL = mean - L * V
    UCL <- mu + L * V                            #UCL = mean + L * V

    if (length(CL) == 1) CL <- rep(CL, n)
    if (length(LCL) == 1) LCL <- rep(LCL, n)
    if (length(UCL) == 1) UCL <- rep(UCL, n)

    if (graph == TRUE){
        plot(c(1, n), c(3 / 4 * min(c(CS, LCL)), 5 / 4 * max(c(CS, UCL))), type = 'n', xlab = 't', ylab = 'CS')
        points(c(1:n), UCL, type = 'l', lty = 2)
        points(c(1:n), CL, type = 'l')
        points(c(1:n), LCL, type = 'l', lty = 2)
        points(c(1:n), CS, type = 'l')
    }

    return(list(CS = CS, CL = CL, LCL = LCL, UCL = UCL))

}


####################################################################################################################################################
    #Get WM statistics
####################################################################################################################################################
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

WM.stat <- function(m, n){

    mu <- n * (m + n + 1) / 2 - m * (m + 1) / 2
    V <- m * n * (m + n + 1) / 12

    list(mu = mu, V = V)
    
}

####################################################################################################################################################
    #RL simulation
####################################################################################################################################################

check.RL <- function(CS, LCL = NULL, UCL = NULL){
                                                    #The purpose of this function is to check RL at each iteration
                                                    #CS, LCL and UCL must have same lengths or CS have n length
    n <- length(CS)                                 #and LCL and UCL have 1 length
                                                    #
                                                    #When LCL or UCL is length 1, this function will recognize 
                                                    #this chart having fixed limits
                                                    #
    if (is.null(UCL) && !is.null(LCL)){             #If UCL is null, it will try to obtain RL by LCL < CS < Inf
                                                    #
        if (n != length(LCL)) {                     #
                                                    #
                                                    #
            if (length(LCL) == 1){                  #
                                                    #
                LCL <- rep(LCL, n)                  #
                                                    #
            } else {                                #
                                                    #
                cat('Error2')
            
            }
        
        
        }
        
        UCL <- rep(Inf, n)
 
    
        
    
    } else if (!is.null(UCL) && is.null(LCL)) {     #If LCL is null, it will try to obtain RL by -Inf < CS < UCL
    
    
        if (n != length(UCL)) {
        
        
            if (length(UCL) == 1){
        
                UCL <- rep(UCL, n)
        
            } else {
            
                cat('Error2')
            
            }
        
        
        }
        
        LCL <- rep(-Inf, n)
    
    
    } else if (!is.null(UCL) && !is.null(LCL))  {  #If both of UCL and LCL are not null, 
                                                   #it will try to obtain RL by LCL < CS < UCL
    
        if (n != length(UCL)) {
        
        
            if (length(UCL) == 1){
        
                UCL <- rep(UCL, n)
        
            } else {
            
                cat('Error2')
            
            }
        
        
        }
        
        
        if (n != length(LCL)) {
        
        
            if (length(LCL) == 1){
        
                LCL <- rep(LCL, n)
        
            } else {
            
                cat('Error2')
            
            }
        
        
        }
        
    
    
    } else {
    
        cat('Error1')
    
    }
        
    if (any(CS - UCL > 0) || any(CS - LCL < 0)){   #OOC is an indicator to show whether this chart is out of control
    
        list(RL = min(c(which(CS - UCL > 0), which(CS - LCL < 0))), OOC = TRUE)
    
    } else {
    
        list(RL = n, OOC = FALSE)
    
    }
        
        
}


get.RL <- function(
                m
                , n = NULL
                , Chart = 'WM'
                , xtype = 'norm'
                , ytype = 'norm'
                , L = 3
                , shift = 0
                , scale = 1
                , subgroup.amt = 1000
){
                                                                                                       #The purpose of this function is
    rl <- 0                                                                                            #to get each RL
    
    x <- get.data(m, 1, type = xtype)                                                                  #generate x
    
    repeat{
    
        if (is.null(n)) {
        
        
        } else {
                                                           
            y <- get.data(n, subgroup.amt = subgroup.amt, type = ytype, shift = shift, scale = scale)                 #generate y
        
            if (Chart == 'WM') {                                                                       #do WM chart
                                                                                                       #
                CS <- WM.f(x, y)                                                                       #get charting staitsitcs
                stat <- WM.stat(m, n)                                                                  #get mean and variance
                Ch <- SPC.Chart(CS, stat$mu, L, sqrt(stat$V), graph = FALSE)                           #get LCL and UCL
            
            } else if (Chart == 'T'){                                                                  #do T chart
                                                                                                       #
                CS <- rowMeans(y)                                                                      #get charting staitsitcs
                mu <- mean(x)                                                                          #get mean
                V <- var(as.vector(x))                                                                            #get variance
                
                Ch <- SPC.Chart(CS, mu, L, sqrt(V), graph = FALSE)                                     #get LCL and UCL
            
            
            } else if (Chart == 'FP'){                                                                 #do WM chart with FP replacement
                             
                CS <- lapply(
                        1:subgroup.amt,
                        function(ii) get.U.Modi.Stat(as.vector(x), y[ii, ]) 
                )
                
                CS <- unlist(CS)
                
                Ch <- SPC.Chart(CS, 0, L, 1, graph = FALSE)                                            #get LCL and UCL
                                                                                                       #with 0 mean and 1 variance by definition
            }
        
        }
        
        
        ll <- check.RL(CS, Ch$LCL, Ch$UCL)                                                             #check whether this iteration is out of control
                                                                                                       #If so, OOC = TRUE.
        rl <- rl + ll$RL
    
        if (ll$OOC == TRUE) {                                                                          #If OOC = TRUE, stop this iteration
                                                                                                       #
            break                                                                                      #
                                                                                                       #
        }                                                                                              #
    
    
    }
    
    rl


}



RL.stat.sim <- function(
                m, n = NULL, 
                Chart = 'WM', 
                xtype = 'norm', 
                ytype = 'norm', 
                L = 3, 
                shift = 0, 
                scale = 1,
                subgroup.amt = 1000, 
                maxsim = 10000,
                cores = 1
){
   
    require(parallel)
                                                                                                        #The purpose of this function is
    #RES <- rep(NA, maxsim)                                                                             #to get RL distribution by simulations
     
    cl <- makeCluster(cores)
    
    clusterExport(
            cl, 
            c('m', 'n', 'Chart', 'xtype', 'ytype', 'L', 'shift', 'subgroup.amt'), 
            envir = environment()
    )
    
    clusterExport(
            cl,
            ls(envir = .GlobalEnv)[sapply(ls(envir = .GlobalEnv), function(x) is.function(get(x)))]
    )

    

    RES <- parLapply(
                cl,
                1:maxsim,
                function(sim) get.RL(m, n, Chart, xtype, ytype, L, shift, scale, subgroup.amt)
    
    
    )
    
    #for (sim in 1:maxsim) {                                                                             #
    #                                                                                                    #
    #    RES[sim] <- get.RL(m, n, Chart, xtype, ytype, L, shift, subgroup.amt)                           #repeatly get RL until sim > maxsim
    #                                                                                                    #
    #}                                                                                                   #
                                                                                                        #
    RES <- unlist(RES)

     
    stopCluster(cl)
         #
    list(
        RL = RES, 
        ARL = mean(RES), 
        SDRL = sd(RES), 
        QRL = quantile(RES, c(0.5, 0.25, 0.5, 0.75, 0.95))
    )                                                                                                   #Output RL and its statistics
    

}


####################################################################################################################################################
    #Search L for a specific ARL
####################################################################################################################################################
bisec.RLsim <- function(
                    ARL, 
                    m, 
                    n = NULL, 
                    Chart = 'WM', 
                    xtype = 'norm', 
                    ytype = 'norm', 
                    shift = 0, 
                    scale = 1,
                    subgroup.amt = 1000, 
                    maxsim = 10000,
                    L.lower.init = 1, 
                    L.upper.init = 5,
                    maxiter = 1000, 
                    tol = 1e-6,
                    cores = 1,
                    double.check = TRUE
){
                                                                                #The purpose of this function is to search
    Given.ARL <- ARL                                                            #L by bisection method
                    
    L.lower <- L.lower.init
    L.upper <- L.upper.init
    
    for (iter in 1:maxiter){
    
        ##cat('iteration:', iter, '\n')
    
        L.mid <- (L.lower + L.upper) / 2
	    
        cat('iteration:', iter, ', L:', L.mid,  '\n')
	    
        Sim.RLD <- RL.stat.sim(m, n, Chart, xtype, ytype, L.mid, shift, scale, subgroup.amt, maxsim, cores)
        
        Sim.ARL <- Sim.RLD$ARL
        
        cat('iteration:', iter, ', ARL:', Sim.ARL,  '\n')
        
        if (abs(Given.ARL - Sim.ARL) < tol) {
        
            if (double.check == TRUE) {
            
                check.RLD <- RL.stat.sim(m, n, Chart, xtype, ytype, L.mid, shift, scale, subgroup.amt, maxsim, cores)
                check.ARL <- check.RLD$ARL
                
                cat('Double Check ARL:', check.ARL,  '\n')
                
                if (abs(Given.ARL - check.ARL) < tol ) {
                
                    break
                    
                } else {
                
                    L.lower <- L.lower.init
                    L.upper <- L.upper.init
                    
                    next
                    
                
                }
                
            } else {
            
                break
            
            }
            
        }
        
        if (Given.ARL - Sim.ARL > 0){
        
            L.lower <- L.mid
        
        } else if (Given.ARL - Sim.ARL < 0) {
        
            L.upper <- L.mid
        
        }
    
    }
    
    #if 

    list(L = L.mid, RL = Sim.RLD)

}

################################################
    #Example
################################################
#start.time <- Sys.time()
#
#a <- bisec.RLsim(
#        ARL = 370, 
#        m = 100, 
#        n = 5, 
#        Chart = 'WM', 
#        xtype = 'norm', 
#        ytype = 'norm', 
#        shift = 0, 
#        subgroup.amt = 1000, 
#        maxsim = 1000,
#        L.lower.init = 1, 
#        L.upper.init = 5, 
#        maxiter = 1000, 
#        tol = 5
#    )
#
#end.time <- Sys.time()
