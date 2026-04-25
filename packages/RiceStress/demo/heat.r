# TODO: Add comment
# 
# Author: jaunario
###############################################################################
library(RiceStress)

#focalDisagg <- function(objraster, maskraster, rm.zeros=TRUE){
#	rs1 <- disaggregate(objraster, fact=6)
#	rs1 <- mask(rs1, maskraster)
#	rs1 <- focal(rs1, fun=mean, w=c(5,5), na.rm=TRUE)
#	rs1 <- focal(rs1, fun=mean, w=c(5,5), na.rm=TRUE)
#	rs1 <- mask(rs1, maskraster)
#	if(rm.zeros) rs1[values(rs1)<1] <- NA 
#	rs1 <- disaggregate(rs1, fact=3)
#	return(rs1)    
#}



stress.nighttimeHeat <- function(con, xy, sowdoy,dae,year, consecutive=TRUE, ...){
	#tminLL <- 25   #threshold for nightime heat stress                    
	tminLL <- 25    #threshold for nightime heat stress + 1.4 adjustment (Bai et al., 2010)                    
	ndays <- 15    # # days that temperature is above the threshold 
	
	# Output vector and names
	out <- vector()
	name <- vector()
	
	sowdate <- dateFromDoy(sowdoy, years)
	harvdate <- sowdate+dae
	seeddate <- sowdate+21
	flwrdate <- harvdate-30
	dstart <- flwrdate-65                               # PI to maturity
	dend <- harvdate
	# GET WEATHER DATA FROM DATABASE
	z <- 
	z <- withRetry(getWthXY(con, xy=xy, years=year, ...)) # query database
	
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
	
	tminseed <-  periodWthVar(z, "tmin", dstart, dend)    
	
	tmindiff_seed <- tminseed-tminLL
	tmindiff_seed[tmindiff_seed<=0] <- NA
	
	tnavg <- colMeans(tminseed, na.rm=TRUE) # Mean tmin per year - NO OUTPUT
	out <- c(out, mean(tnavg, na.rm=TRUE), sd(tnavg, na.rm=TRUE)) # average annual tmin and sd
	name <- c(name, paste("tmin_all", c("avg","sd"), sep="_"))
	
	tndsavg <- colMeans(tmindiff_seed, na.rm=TRUE) # Mean tmin diff from limit - NO OUTPUT
	out <- c(out, mean(tndsavg, na.rm=TRUE), sd(tndsavg, na.rm=TRUE)) # average annual anomaly of tmin from limit and sd
	name <- c(name, paste("tmindiff_all", c("avg","sd"), sep="_"))
	
	potstress <- !is.na(tmindiff_seed)
	potStressCnt <- colSums(potstress)
	potStressCnt[potStressCnt<ndays] <- 0
	out <- c(out, potStressCnt)
	name <- c(name, paste("HotN_all", years, sep="_"))
	
	out <- c(out, mean(potStressCnt, na.rm=TRUE), sd(potStressCnt, na.rm=TRUE))
	name <- c(name, paste("HotN_all", c("avg","sd"), sep="_"))
	
	if (consecutive) LL_seed <- apply(potstress,2,maxConsecutive) else LL_seed <- potStressCnt
	stress <- LL_seed >=ndays
	
	out <- c(out, sum(stress))
	name <- c(name, "HotNYr_all") 
	
	names(out) <- name
	odbcCloseAll()
	return(out)    
}

stress.daytimeHeat <- function(cell,sowdoy,dae,years, consecutive=TRUE, ...){
	#tmaxLL <- 35   #threshold for day time heat stress
	tmaxLL <- 35    #threshold for day time heat stress + 2.8 adjustment (Bai, 2010)
	ndays <- 10    # # days that temperature is above the threshold 
	
	# Output vector and names
	out <- vector()
	name <- vector()
	
	sowdate <- dateFromDoy(sowdoy, years)
	harvdate <- sowdate+dae
	seeddate <- sowdate+21
	flwrdate <- harvdate-30
	dstart <- flwrdate-15
	dend <- flwrdate+21
	# GET WEATHER DATA FROM DATABASE
	xy <- xyFromCell(raster(),cell)
	wthvars <- list(nasa_correction=c("tmax", "tmin", "tdew"))
	wth <- fetch.daily(xy=xy, srcvars=wthvars, connection=con)[[1]]
	z <- wth@w
	#z <- withRetry(getWthXY(cell=cell,...)) # query database
	
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
	
	tmax <-  periodWthVar(z, "tmax", dstart, dend)    
	
	tmaxdiff <- tmax-tmaxLL
	tmaxdiff[tmaxdiff<=0] <- NA
	
	tnavg <- colMeans(tmax, na.rm=TRUE) # Mean tmax per year - NO OUTPUT
	out <- c(out, mean(tnavg, na.rm=TRUE), sd(tnavg, na.rm=TRUE)) # average annual tmax and sd
	name <- c(name, paste("tmax_all", c("avg","sd"), sep="_"))
	
	tndsavg <- colMeans(tmaxdiff, na.rm=TRUE) # Mean tmax diff from limit - NO OUTPUT
	out <- c(out, mean(tndsavg, na.rm=TRUE), sd(tndsavg, na.rm=TRUE)) # average annual anomaly of tmax from limit and sd
	name <- c(name, paste("tmaxdiff_all", c("avg","sd"), sep="_"))
	
	potstress <- !is.na(tmaxdiff) # Check if tmax exceeds thresholds
	potStressCnt <- colSums(potstress) # Count of number of days that tmax exceeds threshold
	potStressCnt[potStressCnt<ndays] <- 0 # Disregard years less than ndays
	out <- c(out, potStressCnt)
	name <- c(name, paste("HotD_all", years, sep="_"))
	
	out <- c(out, mean(potStressCnt, na.rm=TRUE), sd(potStressCnt, na.rm=TRUE))
	name <- c(name, paste("HotD_all", c("avg","sd"), sep="_"))
	
	if (consecutive) LL <- apply(potstress,2,maxConsecutive) else LL <- potStressCnt
	stress <- LL >=ndays
	
	out <- c(out, sum(stress))
	name <- c(name, "HotDYr_all") 
	
	names(out) <- name
	#odbcCloseAll()
	return(out)    
}

thermal.stress <- function(planting, dae, years, type, climdsn, climset, maskraster=NA, outdir=getwd(), verbose=TRUE, ...){
	sowdoy <- getValues(planting)
	daevals <- getValues(dae)
	cells <- which(!is.na(sowdoy) & sowdoy>0)
	stdcell <- cellsFromExtent(raster(),planting)[cells]
	
	stressdata <- vector()
	rem <- vector()    
	for (i in 1:length(cells)){
		if ( is.na(dae[cells[i]])|is.na(planting[cells[i]])){
			rem <- c(rem,i)
			next
		} 
		if (verbose) show.message(type, ": Analyzing weather on cell ", sprintf("%05d",cells[i]), ".", eol="\r")
		if(type=="cold") {
			stressrec <- stress.cold(stdcell[i], sowdoy[cells[i]], daevals[cells[i]], years, database=climdsn, table=climset)
		} else if (type=="nighttimeheat"){
			stressrec <- stress.nighttimeHeat(stdcell[i], sowdoy[cells[i]], daevals[cells[i]], years, consecutive=FALSE, database=climdsn, table=climset)
		} else if (type=="daytimeheat"){
			stressrec <- stress.daytimeHeat(stdcell[i], sowdoy[cells[i]], daevals[cells[i]], years, consecutive=FALSE, database=climdsn, table=climset)
		}
		
		if (length(stressrec)>0) {
			stressdata <- rbind(stressdata,stressrec,deparse.level=0)
		} else {
			if (verbose) show.message(sprintf("%05d",cells[i]), ": No valid output ", eol="\n")
			rem <- c(rem,i)
		}
	}
	if (length(rem)>0) cells <- cells[-rem]
	rownames(stressdata) <- cells 
	stress.map(stressdata, planting, maskraster, writeto=outdir)
	return(TRUE)
}

inpath <- "D:/Google Drive/Projects/Heat Stress Models/Heat Stress Models/Inputs"

fd1  <- list.files(path=inpath,pattern="PLANT_PK[1-3]_60.tif$", full.names=TRUE)
fde1  <- list.files(path=inpath,pattern="HARV_PK[1-3]_60.tif$", full.names=TRUE)

m <- raster(paste(inpath, "rice_1deg_new2.tif", sep="/"))
projection(m) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

planting <- raster(fd1[1])
NAvalue(planting) <- 65535
projection(planting) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
m <- crop(m, planting)
dend1 <-  raster(fde1[1])
NAvalue(dend1) <- 65535
projection(dend1) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
dend1 <- mask(dend1, m)
dend1 <-calc(dend1, fun=function(x){return(ifelse(x > 0, x, NA) )} )
dae<- overlay(planting, dend1, fun=function(x,y){return(ifelse(y > x, y - x, y - x + 365) )} )
dae  <- calc(dae, fun=function(x){return(ifelse(x < 70, 70, x) )} )                    #adjust too low maturity
dae  <- calc(dae, fun=function(x){return(ifelse(x > 200, 200, x) )} )                    #adjust too low maturity


planting  <- calc(planting, fun=function(x){return(ifelse(x > 0, x, NA) )} )

ricemask <- raster(paste(inpath, "rice_1deg_new2.tif", sep="/"))
ricemask <- crop(ricemask, planting)
projection(ricemask) <- projection(planting)


#years <- 1983:2008
years <- 2012
climdsn <- "climate_SRV3A"
climset <- "nasa_correction"

