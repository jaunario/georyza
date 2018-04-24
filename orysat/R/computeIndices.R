# Authors: Federico Filipponi, Francesco Nutini, Giacinto Manfron, Mirco Boschetti
# Consiglio Nazionale delle Ricerche (CNR-IREA) - Italy
# Date :  Dec 2012
# Version 1.0
# Licence GPL v3
# Maintainer: Mirco Boschetti <boschetti.m@irea.cnr.it>, Federico Filipponi <federico.filipponi@gmail.com>
# Description: Compute evi, ndfi, noise indices from MODIS data required for phenorice computation

computeIndices <- function(modis, tile, y, writeto=NULL, verbose=TRUE, asbrick=FALSE, ...) {
	
	# check if output folder exists
	if (is.character(writeto)){
		noisedir <- paste(writeto, y, "noise", sep="/")
		if (file.exists(noisedir)){
			noisefiles <- dir(path=noisedir, pattern=modis@acqdate)
		} else {
			force.directories(noisedir, recursive=TRUE)	
			noisefiles <- vector()
		}
		
		vegdir <- paste(writeto, y, "veg", sep="/")
		if (file.exists(noisedir)){
			vegfiles <- dir(path=vegdir, pattern=modis@acqdate)
		} else {
			force.directories(vegdir, recursive=TRUE)	
			vegfiles <- vector()
		}
	}
	
	if (length(noisefiles)==1 & length(vegfiles)==2){
		if (verbose) message(modis@acqdate, ": Indices previously computed.\r")
		
	} else {
		# duplicate modis data structure to store indices and cloud
		noisem <- modis.data(modis)
		vegind <- modis.data(modis)
		
		# use band 12 to create raster mask and cloud weights from snow, cloud and land coverage information (expressed in bit)
		#if (!(require(RemoteSensing))) {stop("You need to install the RemoteSensing package to use this function")}
		
		# extract cloud presence quality flag information
		iscloud <- modis.sqa500f(modis@imgvals$state_500m)
		iscloud[is.na(iscloud)] <- 1
		
		# extract snow presence quality flag information
		issnow <- modis.sqa500h(modis@imgvals$state_500m)
		issnow[is.na(issnow)] <- 0
		cloud <- as.vector(as.numeric((iscloud+issnow)>0))
		
		# extract cloud quality flag information
		cloudw <- modis.sqa500a(modis@imgvals$state_500m)
		cloudw[is.na(cloudw)] <- 1
		cloudw <- replace(cloudw,cloudw==0,3)
		
		# extract cloud quality information from blue
		cloudblue <- ifelse(modis@imgvals$b03>1800, 1, (ifelse(modis@imgvals$b03<500, 3, (ifelse(modis@imgvals$b03<1800 & modis@imgvals$b03>1100, 2, 0)))))
		blu <- ifelse(modis@imgvals$b03>=1800, 1, (ifelse(modis@imgvals$b03<=1100, 3, (ifelse(modis@imgvals$b03<1800 & modis@imgvals$b03>1100, 2, 0)))))
		noisev <- as.vector(as.numeric(ifelse(cloud==1, 1, ifelse(cloudblue==0, cloudw, ifelse(cloudblue<=cloudw, cloudblue, cloudw)))))
		
		# store cloud results in noise modis data stack
		noisem@imgvals <- data.frame(cbind(noisev))
		names(noisem@imgvals) <- c("noise")
		
		## compute vegetation indices
		# compute EVI index
		evi <- (2.5 * ( ( modis@imgvals$b02 - modis@imgvals$b01 ) / ( modis@imgvals$b02 + ( 6.0 * modis@imgvals$b01 ) - ( 7.5 * modis@imgvals$b03 )  + 10000 ) ) *10000 )
		evi[is.infinite(evi)] <- -10000
		evi[is.na(evi)] <- -10000
		evi[evi < -10000] <- -10000
		evi[evi > 10000] <- 10000
		evi <- as.vector(as.numeric(as.integer(evi)))
		
		## compute NDFI index
		ndfi <- as.vector(rep(NA, length(modis@imgvals$b01)))
		ndfi <- (modis@imgvals$b01 - modis@imgvals$b07) / (modis@imgvals$b01 + modis@imgvals$b07)
		ndfi[is.infinite(ndfi)] <- 0
		ndfi[is.na(ndfi)] <- 0
		ndfi[ndfi < -1] <- -1
		ndfi[ndfi > 1] <- 1
		ndfi <- as.vector(as.numeric(as.integer(ndfi*10000)))
		# Reference: Mirco Boschetti: Boschetti, M., Nelson, A., Nutini, F. rancesc., 2012. Temporal analysis of a robust spectral index for flood detection: a contribution to mapping and monitoring lowland rice. Remote Sensing of Environment under revision.
		
		# store vegetation indices in vegindices modis data stack
		vegind@imgvals <- data.frame(cbind(evi,ndfi))
		
		# write noise to disk
		# check if output folder exists
		if (is.character(writeto)){
			# save noise raster to disk
			if (verbose) message(noisem@acqdate, ": Writing noise raster to disk.\r")
			nbrick <- modis.brick(noisem, intlayers=1:ncol(noisem@imgvals), writeto=noisedir, intNA=-32768, fltNA=-32768, options="COMPRESS=LZW", overwrite=TRUE,...)
			rm(nbrick)
			gc(verbose=FALSE)  
			# save vegetation indices raster to disk
			if (verbose) message(vegind@acqdate, ": Writing 'evi' 'ndfi' indices to disk.\r")
			vbrick <- modis.brick(vegind, intlayers=1:ncol(vegind@imgvals), writeto=vegdir, intNA=-32768, fltNA=-32768, options="COMPRESS=LZW", overwrite=TRUE,...)
			rm(vbrick)
			gc(verbose=FALSE)  
		}
		
	}  
  return(0)
}
