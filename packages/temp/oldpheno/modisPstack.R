# Author: Federico Filipponi on a previous code (modisStack.R in RiceMap package) created bySonia Asilo, Ritsuko Fuchiyama, Robert J. Hijmans, Yann Chemin, Angelo Carlo Pacheco, Jorrel Khalil S. Aunario, Andrew Nelson (International Rice Research Institute)
# Consiglio Nazionale delle Ricerche (CNR-IREA) - Italy
# Date :  Dec 2012
# Version 1.0
# Licence GPL v3
# Maintainer: Mirco Boschetti <boschetti.m@irea.cnr.it>, Federico Filipponi <federico.filipponi@gmail.com>
# Decription: Assemble MODIS layers in a modis.data stack


modisPstack <- function(modfiles, modisdate, bands=c("b01", "b02", "b03", "b07", "state_500m"), verbose=TRUE){
		
  # get only files with acqdate=modisdate
	files <- modfiles[modfiles$acqdate==modisdate,] ###@ query file list using required date
	if(nrow(files)<1){
		bdata <- NULL
	} else {
		if (verbose) show.message(modisdate, ": Reading MODIS images. ", eol="\r")
		pstack <- suppressMessages(stack(files$filename)) ###@ create a raster stack using selected date
		NAvalue(pstack) <- -28672
		stkvals <- as.data.frame(values(pstack))
		colnames(stkvals) <- files$band
		
		# check if required modis for computation have been imported
		bandscheck <- sum(as.numeric(bands %in% names(stkvals)))
		if (bandscheck<5) {
			stop("Required MODIS bands for computation are not present in TIF folder")
		}
		
		# Create modis.data object for mask
		bdata <- new("modis.data") ###@ create an empty "modis.data" object, which function is defined in "modisdata.R"
		bdata@product <- files$product[1] ###@ assign metadata to mdata object retrieved from "files" object
		bdata@acqdate <-  files$acqdate[1]
		bdata@zone <- files$zone[1]
		bdata@version <- files$version[1]
		bdata@proddate <- files$proddate[1]
		bdata@projection <- projection(pstack)
		bdata@extent <- extent(pstack)
		bdata@ncols <- ncol(pstack)
		bdata@nrows <- nrow(pstack)
		
		# fill bdata imgvals slot with band values
		bdata@imgvals <- as.data.frame(stkvals)
		
		# remove unuseful data from workspace
		rm(stkvals,pstack)
		
	}
       
	return(bdata)
}
