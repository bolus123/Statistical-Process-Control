if (sum(c("RevoUtils", "RevoUtilsMath") %in% rownames(installed.packages())) == 2){ #Maximize the computation performance
                                                                                    #based on intel MKL library
    require(RevoUtils)                                                              #In this case, you may need packages, RevoUtils and RevoUtilsMath
    require(RevoUtilsMath)                                                          #but they are just optional libraries
    setMKLthreads(getMKLthreads())

}