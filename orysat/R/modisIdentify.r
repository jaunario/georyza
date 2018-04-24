# Author: Sonia Asilo, Jorrel Khalil S. Aunario
# IRRI
# License GPL3
# Version 2, March 2009

modis.identify <- function(modis, what=c("flooded1", "persistentwater", "drought", "forest", "shrub"), writeto=NA, verbose=TRUE){	
	# check if 64-bit 
	is64 <- version$arch=="x86_64"
	
	if (verbose) message(modis@acqdate, ": Delineating ", paste(what, collapse=", "))
	
	features <- modis.data(modis)
	if (is64) {
		features@imgvals <- modis.compute(modis@imgvals,funlist=what)
	} else {
		for (i in 1:length(what)){
			rbands <- getRequiredBands(what[i])
			input <- as.data.frame(modis@imgvals[,rbands])
			colnames(input) <- rbands
			if (i==1) {
				features@imgvals <- modis.compute(input,funlist=what[i], datatype="logical")
			} else {
				features@imgvals[,what[i]] <- modis.compute(input,funlist=what[i],datatype="logical")[,what[i]]
			}
			rm(input)
			gc(verbose=FALSE)
		}
	}
	
	if (is.character(writeto)){
		if(!file.exists(writeto)){
			dir.create(writeto, recursive=TRUE)	
		} 
		
		outdir <- normalizePath(writeto)
				
		if (verbose) message(modis@acqdate, ": Writing features to disk.")
		modis.brick(features, process="identify", intlayers=1:ncol(modis@imgvals), writeto=outdir, options="COMPRESS=LZW", overwrite=TRUE)		
	}
	if (verbose) message(modis@acqdate, ": ------------------ DONE IDENTIFYING -------------------")
	rm(modis)
	gc(verbose=FALSE)
	return(features)	
}

modis.composite <- function(composite, features=NULL){
	# check if 64-bit 
	is64 <- version$arch=="x86_64"
	
	if(is.null(features)){
		composite@acqdate <- paste(substr(composite@acqdate,1,5),"000",sep="")
		composite@proddate <- ""
		composite@imgvals$nimgs <- 1 
	} else {
		features@imgvals$nimgs <- 1
		if(sum(colnames(features@imgvals)==colnames(composite@imgvals))!=ncol(composite@imgvals)) stop("Invalid features column")
				
		for (i in 1:ncol(composite@imgvals)){
			x <- composite@imgvals[,i]
			y <- features@imgvals[,i]
			nas <- is.na(x) & is.na(y)
			x[is.na(x)] <- 0
			y[is.na(y)] <- 0
			composite@imgvals[,i] <- x+y
			composite@imgvals[nas,i] <- NA
		}
	}
	return(composite)
}