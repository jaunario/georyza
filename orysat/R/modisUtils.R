# Author: Jorrel Khalil Aunario 
# International Rice Research Institute
# Date : 21 May 2010
# Version 0,1
# Licence GPL v3
# DEPRECATE

extract.subdataset <- function(dataset.file, subdataset.idx=0, savepath=".", skip.existing=TRUE){
	dataset.specs <- strsplit(dataset.file, ".")
	gdal.obj <- GDAL.open(dataset.file) 
}

run.MRT <- function(hdffile, outdir=getwd(), MRT_HOME=Sys.getenv("MRT_HOME"), rm.hdf=FALSE, res.files=TRUE, spectral_subset=c(1,1,1,1,0,1,1,0,0,0,0,1,0), output_projection="SIN", resampling_type="NEAREST_NEIGHBOR", OPP="6371007.181 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0",options=vector(),...){
	# DEPRECATE: use extract.subdataset (general version)
	success <- FALSE
	
	if(!file.exists(outdir)){
		dir.create(outdir, recursive=TRUE)	
	} 
	#Check existing TIFF images related to hdffile. 
	xoutput <- dir(outdir, pattern=sub(".hdf","",basename(hdffile)), ...)
	
	# Skip if exists.
	if (length(xoutput)<sum(spectral_subset)){
		
		if(!is.character(hdffile)) {
			message(hdffile," is not a valid HDF file name character string?", appendLF=TRUE)
			return(FALSE)
		}	
		
		if (MRT_HOME=="") {
			message("MRT not installed. Download it here (https://lpdaac.usgs.gov/lpdaac/tools/modis_reprojection_tool)", appendLF=TRUE)
		} else {
			MRT <- paste(MRT_HOME,"bin", sep="/")
			
			filename <- paste(MRT, "/modisconfig.prm", sep="")
			mrtconfig <- c(paste('INPUT_FILENAME = ', hdffile, sep=""), 
					paste('SPECTRAL_SUBSET = ( ', paste(spectral_subset, collapse=" "),' )', sep=""),
					paste('OUTPUT_FILENAME = ', outdir,"/", sub(".hdf","",basename(hdffile)),'.tif', sep=""), 
					paste('RESAMPLING_TYPE =', resampling_type), 
					paste('OUTPUT_PROJECTION_TYPE =', toupper(output_projection)),
					paste('OUTPUT_PROJECTION_PARAMETERS = (', OPP,')'),
					options)
			writeLines(mrtconfig,filename)
			success <- system(paste(MRT, '/resample -p ', MRT, '/modisconfig.prm', sep=""))
			if (success==0) {
				success <- TRUE
				xoutput <- dir(outdir, pattern=sub(".hdf","",basename(hdffile)), ...)
			} else success <- FALSE 
			if (rm.hdf) unlink(hdffile)
		}
		
	} else success <- TRUE
	
	if (res.files){
		success <- xoutput
	}
	return(success)
}

fillwith.adjmean <- function(x,...){
	center <- x[5]
	if(!is.na(center)) {
		adjmean <- center
	} else {
		adjmean <- mean(x, ...)  
		if (is.nan(adjmean)) adjmean <- NA
	}
	return(adjmean)
}

inventory.modis <- function(modisfiles, sep="\\.", modisinfo=c('product', 'acqdate', 'zone', 'version', 'band'), file.ext="tif") {
	
	#if (format=="GTiff") info <- sub(".hdf","",basename(modisfiles))	
	info <- basename(modisfiles)
	
	# Remove Extension
	info <- gsub(paste(".", file.ext ,sep=""),"",info)
	
	x <- unlist(strsplit(info, sep))
	m <- data.frame(matrix(x, ncol=length(modisinfo), byrow=TRUE), stringsAsFactors=FALSE)
	if (ncol(m) != length(modisinfo)) { 
		message("Non-standard filenames found ")
		return(vector()) 
	}
	colnames(m) <- modisinfo
	m$year <- as.numeric(substr(m$acqdate,2,5))
	m$doy <- as.numeric(substr(m$acqdate,6,8))
	m$filename <- modisfiles
	m <- m[,c("filename", "year", "doy", modisinfo)]
	if ("band" %in% modisinfo) m$band <- as.character(sub("sur_refl_","",m$band))
	return(m)
}

getRequiredBands <- function(funlist, byfun=FALSE){
	if (byfun) argnames <- list() else argnames <- NULL  
	for (i in 1:length(funlist)){
		if (byfun) argnames[[funlist[i]]] <- names(formals(get(funlist[i]))) else argnames <- c(argnames,names(formals(get(funlist[i]))))
	}
	if (byfun) names(argnames) <- funlist else argnames <- unique(argnames)
	return(argnames)
}

raster2SGDF <- function(baseraster, vals=NULL){
	if (!is.null(vals)) {
		baseraster <- setValues(baseraster, vals)
	}
	baseraster <- as(baseraster, 'SpatialGridDataFrame')
	return(baseraster)
}

formatExt <- function(myformat){
    ext <- rep(NA, length(myformat))
    ext[tolower(myformat)=="raster"] <- "grd"
    ext[tolower(myformat)=="gtiff"] <- "tif" 
    ext[tolower(myformat)=="hdf"] <- "hdf"
    return(ext)    
}

extFormat <- function(filext){
	formt <- rep(NA,length(filext))
	formt[tolower(filext)=="tif"] <- "GTiff"
	formt[tolower(filext)=="grd"] <- "raster"
	return(formt)
}

bandnames <- function(bandnum, ref="ricemap"){
	bands <- as.data.frame(cbind(c("red", "nir1", "blue", "green", "nir2", "swir1", "swir2"), c("red", "nir", "blue", "green", NA, "swir1", "swir2")),stringsAsFactors=FALSE)
	colnames(bands) <- c("default","ricemap")
	return(bands[bandnum,ref])
} 

bandnumber <- function(bandname, ref="ricemap", asString=TRUE){
	bands <- cbind(c("red", "nir1", "blue", "green", "nir2", "swir1", "swir2"), c("red", "nir", "blue", "green", NA, "swir1", "swir2"), c("red", "nir", "blue", "green", "swir1", "swir2","swir3"))
	colnames(bands) <- c("default","ricemap","web")
	if (asString) result <- paste("b",gsub(" ",0,format(which(bands[,ref] %in% bandname),width=2)),sep="") else result <- which(bands[,ref] %in% bandname)
	return(result)
} 

validFolders <- function(styear=2000, enyear=as.numeric(format(Sys.Date(),"%Y"))){
	valid <- vector()
	for (y in styear:enyear){
		st <- ifelse(y==2000,paste(y,"2","18",sep="-"),paste(y,"1","1",sep="-"))
		valid <- c(valid,format(seq(from=as.Date(st), to=as.Date(paste(y,"12","31",sep="-")), by=8),"%Y.%m.%d"))
	}
	return(valid)
}

subFolderFromDoy <- function(doy,year){
	valid <- validFolders()
	dirdate <- format(as.Date(paste(doy, year),"%j %Y"),"%Y.%m.%d")
	dirdate[dirdate %in% valid] 
	return(dirdate[dirdate %in% valid])
}
