# Author: Federico Filipponi
# Consiglio Nazionale delle Ricerche (CNR-IREA) - Italy
# Date :  Dec 2012
# Version 1.0
# Licence GPL v3
# Maintainer: Mirco Boschetti <boschetti.m@irea.cnr.it>, Federico Filipponi <federico.filipponi@gmail.com>
# Description: Assemble Phenorice results in original tile extent using GDAL libraries

subsetStack <- function(mpath, subpath, season=3, tile=tile, startYear=startYear){
  
  # set input path
  subsetdir <- paste(subpath, "/",sep="")
  product <- "MOD09A1"
  phenolist <- c("flood","SoS","peak","EoS")
  productlist <- c("phenoricemap","cropseasons","riceseasons","seasonofrice","ricemap.s1","ricemap.s2","ricemap.s3")
  
  # set output path
  outputlayerdir <- paste(mpath, startYear, "phenorice", sep="/")
  force.directories(outputlayerdir)
  
  stackinfo <- paste(product, startYear, tile, sep="_")
  
  FWTOOLS_HOME <- Sys.getenv("FWTOOLS_HOME")
  
  if (FWTOOLS_HOME=="") {
    show.message("FWTOOLS not installed. Download it here (http://fwtools.maptools.org/)", appendLF=TRUE)
  } else {
    if (Sys.info()["sysname"]=="Windows"){
      gdalbuildvrt <- shQuote(paste(FWTOOLS_HOME,"bin", "gdalbuildvrt.exe", sep="/"))
      gdaltranslate <- shQuote(paste(FWTOOLS_HOME,"bin", "gdal_translate.exe", sep="/"))
    } else {
      gdalbuildvrt <- paste(FWTOOLS_HOME,"bin_safe", "gdalbuildvrt", sep="/")
      gdaltranslate <- paste(FWTOOLS_HOME,"bin_safe", "gdal_translate", sep="/")
    } 
	# start cicle for ricemap output assembling
	for(p in 1:length(productlist)){
		message("Create output file for ", productlist[p], " results", appendLF=TRUE)
		system(paste(gdalbuildvrt, 
					paste(subpath, "/", stackinfo, ".", productlist[p], ".vrt", sep=""),
					paste(subpath, "/*", productlist[p], ".tif", sep="")),invisible=TRUE)
		
		system(paste(gdaltranslate, "-of GTiff -ot Byte -co COMPRESS=LZW -a_nodata 255", 
					paste(subpath, "/", stackinfo, ".", productlist[p], ".vrt", sep=""),
					paste(outputlayerdir, "/", stackinfo, "_", sub("\\.","_",productlist[p]), ".tif", sep="")),invisible=TRUE)
	}
	
	# start cicle for phenological output assembling
	message("seasons=", season, appendLF=TRUE)
	for (s in 1:season){
		message("season=", s, appendLF=TRUE)
		for(p in 1:length(phenolist)){
			message("Create output file for ", phenolist[p], " results", appendLF=TRUE)
			system(paste(gdalbuildvrt,  
						paste(subpath, "/", stackinfo, ".", phenolist[p], ".A.", s, ".vrt", sep=""),
						paste(subpath, "/*", phenolist[p], ".A.", s, ".tif", sep="")),invisible=TRUE)
			system(paste(gdalbuildvrt, 
						paste(subpath, "/", stackinfo, ".", phenolist[p], ".B.", s, ".vrt", sep=""),
						paste(subpath, "/*", phenolist[p], ".B.", s, ".tif", sep="")),invisible=TRUE)
			system(paste(gdalbuildvrt, 
						paste(subpath, "/", stackinfo, ".", phenolist[p], ".C.", s, ".vrt", sep=""),
						paste(subpath, "/*", phenolist[p], ".C.", s, ".tif", sep="")),invisible=TRUE)
			system(paste(gdalbuildvrt, "-separate", 
							paste(subpath, "/", stackinfo, ".", phenolist[p], ".", s, ".vrt", sep=""), 
							paste(subpath, "/", stackinfo, ".", phenolist[p], ".A.", s, ".vrt", sep=""),
							paste(subpath, "/", stackinfo, ".", phenolist[p], ".B.", s, ".vrt", sep=""),
							paste(subpath, "/", stackinfo, ".", phenolist[p], ".C.", s, ".vrt", sep="")),invisible=TRUE)
			system(paste(gdaltranslate, "-of GTiff -ot Int32 -co COMPRESS=LZW -a_nodata 0", 
							paste(subpath, "/", stackinfo, ".", phenolist[p], ".", s, ".vrt", sep=""),
							paste(outputlayerdir, "/", stackinfo, "_", phenolist[p], "_", s, ".tif", sep="")),invisible=TRUE)
		}
	}
	
  }
  
}

subsetEsStack <- function(mpath, subpath, season=3, tile=tile, startYear=startYear, verbose=FALSE){
  
	# set output path
	outputlayerdires <- paste(mpath, startYear, "stack", sep="/")
	force.directories(outputlayerdires)
	
	# set input path
	stackinfo <- paste("MOD09A1", startYear, tile, sep=".")
	eslist <- c("flood","delta","peak")
	
	# create variables for assembling evismooth
	evismoothdir <- paste(subpath, "/evismooth", sep="")
	esdirspath <- list.files(path=evismoothdir,full.names=TRUE)
	
	FWTOOLS_HOME <- Sys.getenv("FWTOOLS_HOME")
	  
	if (FWTOOLS_HOME=="") {
		message("FWTOOLS not installed. Download it here (http://fwtools.maptools.org/)", appendLF=TRUE)
	} else {
		if (Sys.info()["sysname"]=="Windows"){
			gdalbuildvrt <- shQuote(paste(FWTOOLS_HOME,"bin", "gdalbuildvrt.exe", sep="/"))
			gdaltranslate <- shQuote(paste(FWTOOLS_HOME,"bin", "gdal_translate.exe", sep="/"))
		} else {
			gdalbuildvrt <- paste(FWTOOLS_HOME,"bin_safe", "gdalbuildvrt", sep="/")
			gdaltranslate <- paste(FWTOOLS_HOME,"bin_safe", "gdal_translate", sep="/")
	    }
		
		message("Create output file for smoothed EVI", appendLF=TRUE)
		if (length(esdirspath==8)){
			for(s in 1:8){
				system(paste(gdalbuildvrt, "-separate",
								paste(subpath, "/", stackinfo, ".", s, ".es.vrt", sep=""), 
								paste(esdirspath[s], "/*es.tif", sep="")))
			}
			
			# create evismooth output in TIF format for 77 doys
			system(paste(gdalbuildvrt,
							paste(subpath, "/", stackinfo, ".evi.smooth.vrt", sep=""),
							paste(subpath, "/*.es.vrt", sep="")))
			system(paste(gdaltranslate, "-of GTiff -ot Int16 -co COMPRESS=LZW -a_nodata -32768",
							paste(subpath, "/", stackinfo, ".evi.smooth.vrt", sep=""),
							paste(outputlayerdires, "/", stackinfo, ".evi.smooth.tif", sep="")))
			
		} else {
			if (verbose) {message("EVI smooth calculation failed for at least one subtile. Intermediate results are stored in: ", subpath, appendLF=TRUE)}
		}
		
		# create flood, peak and delta using evismoothed values
		for(s in 1:season){
			message("season=", s, appendLF=TRUE)
			for(p in 1:length(eslist)){
				message("Create output file for ", eslist[p], " results", appendLF=TRUE)
				system(paste(gdalbuildvrt, 
								paste(subpath, "/", stackinfo, ".", eslist[p], ".evi.", s, ".vrt", sep=""),
								paste(subpath, "/*", eslist[p], ".evi.", s, ".tif", sep="")))
				system(paste(gdaltranslate, "-of GTiff -ot Int16 -co COMPRESS=LZW -a_nodata -32768",
								paste(subpath, "/", stackinfo, ".", eslist[p], ".evi.", s, ".vrt", sep=""),
								paste(outputlayerdires, "/", stackinfo, ".", eslist[p], ".evi.", s, ".tif", sep="")))
			}
		}
	}
}
