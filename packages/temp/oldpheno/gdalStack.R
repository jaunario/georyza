# Author: Federico Filipponi
# Editor: Jorrel Khalil Aunario
# Consiglio Nazionale delle Ricerche (CNR-IREA) - Italy
# Date :  Dec 2012
# Version 1.2
# Licence GPL v3
# Maintainer: Mirco Boschetti <boschetti.m@irea.cnr.it>, Federico Filipponi <federico.filipponi@gmail.com>
# Description: Create multitemporal stacks from indices required for phenorice computing using GDAL libraries

# gdalbuildvrt is not available on FWTools-2.0.6 (current stable release for linux)
# TODO: check FWTOOLS version if at least 3.0
# Correct Installed FWTOOLS bin folder is "bin_safe"

# normalizePath ruins links which have embedded permissions

tiftostack <- function(index, inputdir, outputstackdir, stackinfo, ot=Int16 , FWTOOLS_HOME=Sys.getenv("FWTOOLS_HOME"), subset=TRUE) {
  
  # normalize directories path
  # outputstackdir <- normalizePath(outputstackdir, winslash="/")
  #inputdir <- normalizePath(inputdir, winslash="/")
	if (file.exists(outputstackdir)) {
		indexfiles <- dir(path=outputstackdir, pattern=paste(index, ".stack.vrt$", sep="")) 
	} else {
		force.directories(outputstackdir,recursive=TRUE)
		indexfiles <- vector()
	}
	
	if (length(indexfiles)==0) {
		success <- FALSE
		if (FWTOOLS_HOME=="") {
			show.message("FWTOOLS not installed. Download it here (http://fwtools.maptools.org/)", eol="\n")
		} else {
			if (Sys.info()["sysname"]=="Windows"){
				gdalbuildvrt <- shQuote(paste(FWTOOLS_HOME,"bin", "gdalbuildvrt.exe", sep="/"))
				gdaltranslate <- shQuote(paste(FWTOOLS_HOME,"bin", "gdal_translate.exe", sep="/"))
			} else {
				gdalbuildvrt <- paste(FWTOOLS_HOME,"bin_safe", "gdalbuildvrt", sep="/")
				gdaltranslate <- paste(FWTOOLS_HOME,"bin_safe", "gdal_translate", sep="/")
			} 
			
			system(paste(gdalbuildvrt, " -separate ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", inputdir, "/*", index, ".tif", sep=""))
			
			if (subset) {
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " -srcwin 0 0 600 800 ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack_sub1.tif", sep=""))
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " -srcwin 600 0 600 800 ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack_sub2.tif", sep=""))
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " -srcwin 1200 0 600 800 ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack_sub3.tif", sep=""))
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " -srcwin 1800 0 600 800 ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack_sub4.tif", sep=""))
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " -srcwin 0 800 600 800 ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack_sub5.tif", sep=""))
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " -srcwin 600 800 600 800 ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack_sub6.tif", sep=""))
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " -srcwin 1200 800 600 800 ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack_sub7.tif", sep=""))
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " -srcwin 1800 800 600 800 ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack_sub8.tif", sep=""))
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " -srcwin 0 1600 600 800 ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack_sub9.tif", sep=""))
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " -srcwin 600 1600 600 800 ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack_sub10.tif", sep=""))
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " -srcwin 1200 1600 600 800 ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack_sub11.tif", sep=""))
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " -srcwin 1800 1600 600 800 ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack_sub12.tif", sep=""))
			} else {
				system(paste(gdaltranslate, " -of GTiff -ot ", ot, " ", outputstackdir, "/", stackinfo, ".", index, ".stack.vrt ", outputstackdir, "/", stackinfo, ".", index, ".stack.tif", sep=""))
			}
		}	
	} else {
		message(index, " multitemporal stack previously processed.", appendLF=TRUE)
	}
  
}

lsttostack <- function(mpath=Sys.getenv("MODIS_LOCAL_DIR"), startYear=startYear, tile=tile, verbose=TRUE, subset=TRUE, lst=lststacklist, FWTOOLS_HOME=Sys.getenv("FWTOOLS_HOME")) {
  
  # set folders for stack creation
  outputstackdir <- paste(mpath, startYear, "stack", sep="/")
  if (file.exists(outputstackdir)) {
	  lstfiles <- dir(path=outputstackdir, pattern="lst.stack.vrt$") 
  } else {
	  force.directories(outputstackdir,recursive=TRUE)
	  lstfiles <- vector()
  }
  
  if (length(lstfiles)==0) { 
	  product <- c("MOD11A2")
	  stackinfo <- paste(product, startYear, tile, sep=".")
	  if (verbose) message("Creating LST multitemporal stack for year ", startYear)
	  
	  # normalize directories path REDUNDANT
	  # lststackdir <- normalizePath(outputstackdir, winslash="/")
	  
	  success <- FALSE
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
		  
		  system(paste(gdalbuildvrt, " -allow_projection_difference -separate ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", lst, sep=""))
		  
		  if (subset) {
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 -srcwin 0 0 600 800 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack_sub1.tif", sep=""))
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 -srcwin 600 0 600 800 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack_sub2.tif", sep=""))
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 -srcwin 1200 0 600 800 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack_sub3.tif", sep=""))
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 -srcwin 1800 0 600 800 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack_sub4.tif", sep=""))
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 -srcwin 0 800 600 800 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack_sub5.tif", sep=""))
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 -srcwin 600 800 600 800 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack_sub6.tif", sep=""))
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 -srcwin 1200 800 600 800 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack_sub7.tif", sep=""))
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 -srcwin 1800 800 600 800 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack_sub8.tif", sep=""))
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 -srcwin 0 1600 600 800 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack_sub9.tif", sep=""))
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 -srcwin 600 1600 600 800 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack_sub10.tif", sep=""))
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 -srcwin 1200 1600 600 800 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack_sub11.tif", sep=""))
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 -srcwin 1800 1600 600 800 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack_sub12.tif", sep=""))
		  } else {
			  system(paste(gdaltranslate, " -of GTiff -ot UInt16 ", outputstackdir, "/", stackinfo, ".lst.stack.vrt ", outputstackdir, "/", stackinfo, ".lst.stack.tif", sep=""))
		  }
	  }  
  } else {
	  message("LST multitemporal stack for year ", startYear, " previously processed.", appendLF=TRUE)
  }	  
}
