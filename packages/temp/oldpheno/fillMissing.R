# Author: Federico Filipponi
# Consiglio Nazionale delle Ricerche (CNR-IREA) - Italy
# Date :  Dec 2012
# Version 1.0
# Licence GPL v3
# Maintainer: Mirco Boschetti <boschetti.m@irea.cnr.it>, Federico Filipponi <federico.filipponi@gmail.com>
# Description: Create evi, ndfi, noise indices with neutral values for filling MODIS missing date required for phenorice computation

fabricate.missing <- function(acqdate, layers, writeto=".", value=-10000, verbose=TRUE) {
	success <- FALSE
	
	# check if output folder exists
	if (!file.exists(writeto)){
		stop("Output folder ", writeto, " not found.")
	}
	
	# check values and layers correspondence
	if (length(layers)!=length(value)){	
		stop("Length of layers and values don't match.")
	}
		
	mfiles <- modisFiles(path=writeto, full.names=TRUE, pattern=layers[1])
	
	if (nrow(mfiles)==0){
		stop("No MODIS images found.")
	}
	
	mdata <- modis.data(x=mfiles$filename[1], nodata=TRUE)
	mdata@imgvals <- data.frame(matrix(data=value,ncol=length(layers),nrow=5760000, byrow=TRUE))
	colnames(mdata@imgvals) <- layers
	mdata@proddate <- format(Sys.time(),"%Y%j%H%M%S")
	for (i in acqdate) {
		success <- FALSE
		mdata@acqdate <- i
		# save noise raster to disk
		if (verbose) message(mdata@acqdate, ": Writing layers (", paste(layers, collapse=", "),") raster to disk.", eol="\r")
		nbrick <- modis.brick(mdata, intlayers=1:ncol(mdata@imgvals), writeto=writeto, intNA=-32768, fltNA=-32768, options="COMPRESS=LZW", overwrite=TRUE)
		rm(nbrick)
		gc(verbose=FALSE,reset=TRUE)
			
	}
	success <- TRUE
	rm(mdata)
	gc(verbose=FALSE,reset=TRUE)
	return(success)
}
