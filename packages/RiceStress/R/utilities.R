# TODO: Add comment
# 
# Author: jaunario
###############################################################################

periodWthVar <- function(wth, var, pstart, pend){
	wth <- wth[order(wth$date),]
	vardat <- vector()
	for (i in 1:length(pstart)){         
		vardat <- cbind(vardat,wth[wth$date>=pstart[i] & wth$date<=pend[i],var])    
	}
	return(vardat)
}
