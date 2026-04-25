# TODO: Add comment
# 
# Author: jaunario
###############################################################################

#getWthXY <- function (dbcon, tablename, xy, years = "all", verbose = FALSE) {
#	
#	vars <- list(c("tmin","tmax"))
#	names(vars) <- tablename
#	if (length(years)>1){
#		w <- fetch.daily(xy, srcvars=vars, connection=dbcon, stdate=as.Date(paste(min(years),"-1-1",sep="")), endate=as.Date(paste(max(years)+1,"-12-31",sep="")))[[1]]
#	} else if (length(years)==1 & years == "all") {
#		w <- fetch.daily(xy, srcvars=vars, connection=dbcon)[[1]]
#	}	
#	w <- w@w
#	year <- yearFromDate(w$date)
#	doy <- doyFromDate(w$date)
#	#vars <- colnames(w)[4:ncol(w)]
#	w <- cbind(w$date, year, doy, w[, -1])
#	colnames(w) <- c("date", "year", "doy", unlist(vars))
#	#if (nrow(w) > 0) {
#	#    tavg <- which(colnames(w) %in% c("t2m", "tavg"))
#	#    if (length(tavg) == 0) {
#	w$tavg <- (w$tmin + w$tmax)/2
#	#        tavg <- which(colnames(w) == "tavg")
#	#    }
#	#}
#	return(w)
#}
#
#
#
stress.cold <- function(con, table, xy, sowdoy, dae, years, ...){
	# INITIALIZE OUTPUT VECTOR AND NAMES
		out <- vector()
		name <- vector()
	# END OF OUTPUT INITIALIZATION
	
	# SET THRESHOLDS
		# Seedling stage thresholds
		tminLL_seed <- 10    #threshold for cold stress + 1.4 adjustment (Bai et al., 2010)
		tavgLL_seed <- 15    #threshold for cold stress + (2.8 + 1.4) / 2 adjustment 
		
		# Reproductive stage thresholds
		tminLL_pnin <- 17    # As per email from aglaborte 21 Oct 2013 
		tminLL_flwr <- 15    # As per email from aglaborte 21 Oct 2013
	# END OF THRESHOLDS
		
	# GET WEATHER DATA FROM DATABASE
		z <- withRetry(getWthXY(dbcon=con, tablename=table, xy=xy, years=years, ...)) # query database
		# check if query successfull
		if (length(z)==0){ 
			show.message("Failed to connect to server. Rerun later.\n",eol="\r")
			return(out)
		}
		# check if resultset is not empty
		if (nrow(z) <= 0) { 
			show.message("No weather data.\n", eol="\r")
			return(out)
		}
	# END OF GET WEATHER
	
	
	# Derive stage dates from planting and dae
		sowdate <- dateFromDoy(sowdoy, years) # Sowing/planting date
		seeddate <- sowdate+21 # End of seedling stage date
		
		harvdate <- sowdate+dae # Harvest date
		flwrdate <- harvdate-30 # Flowering date
		pnindate <- harvdate-65 # Panicle initiation date
			
	
	# Get dates of the period that covers the stage of interest
		tminseed <-  periodWthVar(z, "tmin", sowdate, seeddate)    
		tavgseed <-  periodWthVar(z, "tavg", sowdate, seeddate)

		tminflwr <-  periodWthVar(z, "tmin", flwrdate-4, flwrdate+4)
		tminpnin <-  periodWthVar(z, "tmin", pnindate+14-5, pnindate+14+5)
		

	# Check which dates exceed the thresholds
		tmindiff_seed <- -(tminseed-tminLL_seed)
		tmindiff_seed[tmindiff_seed<=0] <- NA
	
		tavgdiff_seed <- -(tavgseed-tavgLL_seed)
		tavgdiff_seed[tavgdiff_seed<=0] <- NA
		
		tmindiff_flwr <- -(tminflwr-tminLL_flwr)
		tmindiff_flwr[tmindiff_flwr<=0] <- NA
		
		tmindiff_pnin <- -(tminpnin-tminLL_pnin)
		tmindiff_pnin[tmindiff_pnin<=0] <- NA
		
		
	# TODO check original computation of statistics
	tnsavg <- colMeans(tminseed, na.rm=TRUE)
	out <- c(out, mean(tnsavg, na.rm=TRUE), sd(tnsavg, na.rm=TRUE))
	name <- c(name, paste("tmin_seed", c("avg","sd"), sep="_"))
	
	tnfavg <- colMeans(tminflwr, na.rm=TRUE)
	out <- c(out, mean(tnfavg, na.rm=TRUE), sd(tnfavg, na.rm=TRUE))
	name <- c(name, paste("tmin_flwr", c("avg","sd"), sep="_"))
	
	tnpavg <- colMeans(tminpnin, na.rm=TRUE)
	out <- c(out, mean(tnpavg, na.rm=TRUE), sd(tnpavg, na.rm=TRUE))
	name <- c(name, paste("tmin_pnin", c("avg","sd"), sep="_"))
	
	tnavg <- colMeans(rbind(tminseed, tminpnin, tminflwr), na.rm=TRUE)
	out <- c(out, mean(tnavg, na.rm=TRUE), sd(tnavg, na.rm=TRUE))
	name <- c(name, paste("tmin_all", c("avg","sd"), sep="_"))
	
	tndsavg <- colMeans(tmindiff_seed, na.rm=TRUE)
	out <- c(out, mean(tndsavg, na.rm=TRUE), sd(tndsavg, na.rm=TRUE))
	name <- c(name, paste("tmindiff_seed", c("avg","sd"), sep="_"))
	
	tndfavg <- colMeans(tmindiff_flwr, na.rm=TRUE)
	out <- c(out, mean(tndfavg, na.rm=TRUE), sd(tndfavg, na.rm=TRUE))
	name <- c(name, paste("tmindiff_flwr", c("avg","sd"), sep="_"))
	
	tndpavg <- colMeans(tmindiff_pnin, na.rm=TRUE)
	out <- c(out, mean(tndpavg, na.rm=TRUE), sd(tndpavg, na.rm=TRUE))
	name <- c(name, paste("tmindiff_pnin", c("avg","sd"), sep="_"))
	
	tndavg <- colMeans(rbind(tndsavg, tndfavg, tndpavg), na.rm=TRUE)
	out <- c(out, mean(tndavg, na.rm=TRUE), sd(tndavg,na.rm=TRUE))
	name <- c(name, paste("tmindiff_all", c("avg","sd"), sep="_"))
	
	stressseed <- !is.na(tmindiff_seed) | !is.na(tavgdiff_seed)
	stressseedC <- colSums(stressseed)
	out <- c(out, stressseedC)
	name <- c(name, paste("Cold_seed", years, sep="_"))
	
	out <- c(out, mean(stressseedC, na.rm=TRUE), sd(stressseedC, na.rm=TRUE))
	name <- c(name, paste("ColdDays_seed", c("avg","sd"), sep="_"))
	
	stressflwr <- !is.na(tmindiff_flwr)
	stressflwrC <- colSums(stressflwr)
	out <- c(out, stressflwrC)
	name <- c(name, paste("Cold_flwr", years, sep="_"))
	out <- c(out, mean(stressflwrC, na.rm=TRUE), sd(stressflwrC, na.rm=TRUE))
	name <- c(name, paste("ColdDays_flwr", c("avg","sd"), sep="_"))
	
	stresspnin <- !is.na(tmindiff_pnin)
	stresspninC <- colSums(stresspnin)
	out <- c(out, stresspninC)
	name <- c(name, paste("Cold_pnin", years, sep="_"))
	out <- c(out, mean(stresspninC, na.rm=TRUE), sd(stresspninC, na.rm=TRUE))
	name <- c(name, paste("ColdDays_pnin", c("avg","sd"), sep="_"))
	
	stressC <- colSums(rbind(stressseedC,stressflwrC,stresspninC), na.rm=TRUE)
	out <- c(out, stressC)
	name <- c(name, paste("Cold_all", years, sep="_"))
	
	out <- c(out, mean(stressC, na.rm=TRUE), sd(stressC, na.rm=TRUE))
	name <- c(name, paste("ColdDays_all", c("avg","sd"), sep="_"))
	
	LL_seed <- apply(stressseed,2,maxConsecutive)
	stressS <- LL_seed >=5
	#LL_repr <- apply(stressrepr,2,maxConsecutive)
	#stressR <- LL_repr >=10
	stressR <- stressflwrC>=3 | stresspninC>=3
	
	out <- c(out, sum(stressS))
	name <- c(name, "ColdYr_seed") 
	out <- c(out, sum(stressR)) 
	name <- c(name, "ColdYr_repr")
	out <- c(out,sum(stressS|stressR)) 
	name <- c(name, "ColdYr_all")
	
	names(out) <- name
	return(out)    
	
}
