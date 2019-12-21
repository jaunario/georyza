# Author: Jorrel Khalil S. Aunario
# International Rice Research Institute
# Date :  25 May 2011
# Version 0.0.1
# Licence GPL v3

modis.compute <- function(modis, funlist ,datatype="numeric"){
	# DEPRECATE
    result <- do.call(datatype,as.list(0))
    for (i in 1:length(funlist)){
        argnames <- names(formals(get(funlist[i])))
        if (sum(argnames %in% colnames(modis))!=length(argnames)) {
            assign(toupper(funlist[i]), rep(NA,nrow(modis)))
        } else {
            arglist <- list()
            for (arg in argnames){
                arglist[[arg]] <- modis[,arg] 
            }
            assign(toupper(funlist[i]), do.call(funlist[i],arglist))
        }
        result <- as.data.frame(cbind(result, get(toupper(funlist[i]))))        
    }
    colnames(result) <- funlist
    return(result)
}
