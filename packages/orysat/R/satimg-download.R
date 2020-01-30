# Author: Jorrel Khalil S. Aunario
# IRRI
# License GPL3
# Version 1, August 2011

DL.FAST  = 0 # If already exists, skip
DL.SMART = 1 # If already exists, check integrity. If integrity passed, skip else delete then download again
DL.RENEW = 2 # If exists redownload

check.integrity <- function(modis.hdf, xml){
	cksumver <- Sys.which("cksum")	
	if (cksumver==""){
		cksum <- file.info(modis.hdf)$size
		chk <- xml[grep("FileSize>",xml)]
		idx <- unlist(gregexpr("[[:digit:]]", chk))
		chk <- as.numeric(substr(chk, min(idx), max(idx)))
	} else {
		cksum <- system(paste("cksum", modis.hdf), intern=TRUE)
		cksum <- unlist(strsplit(cksum[length(cksum)], " "))[1]
		chk <- xml[grep("Checksum>",xml)]						
		idx <- unlist(gregexpr("[[:digit:]]", chk))
		chk <- substr(chk, min(idx), max(idx))
		if (grepl("\\|", chk)) chk <- unlist(strsplit(cksum[length(chk)], "\\|"))[1]
	}
	return(cksum==chk)
}

download.modis <- function(tile, years, userpwd, doy=NULL, product="MOD09A1", prod.ver=6, savepath=NULL, modis.site="http://e4ftl01.cr.usgs.gov/", dl.mode=DL.SMART, integrity=TRUE, skip.exists=TRUE, verbose=TRUE, ...){
	# Check time, stop if server downtime
	
	
	
	# Initialize required objects
	if (is.null(savepath)){
		dirs <- apply(expand.grid(product, tile, stringsAsFactors = FALSE), FUN = paste, MARGIN = 1, collapse="/")
		lapply(dirs, FUN=dir.create,recursive=TRUE)
	} else if (!file.exists(savepath)){
		dir.create(savepath,recursive=TRUE)
	}
	
	if (dl.mode==DL.FAST) {
		skip.exists <- TRUE
		integrity <- FALSE
	} else if (dl.mode==DL.SMART) {
		integrity <- TRUE
	} else if (dl.mode==DL.RENEW) {
		skip.exists <- FALSE			
	}
	
	result <- vector()
	curl_handle <- new_handle(userpwd=userpwd)
	
	for (pr in product){
		prod.info <- modprods[grep(pr,modprods$ShortName),]
		if(!grepl("day", prod.info$Temporal.Granularity)){
			warning("Unsupported product ", product,". Kindly contact the developer.")
			next
		}
		
		prod.site <- paste(modis.site, "MO", switch(prod.info$Platform, Aqua="LA", Terra="LT", Combined="TA"), "/", paste(pr,sprintf(paste("%03d",sep=""),prod.ver),sep="."),sep="")
		tim.gran <- paste("t",gsub(" ", "", prod.info$Temporal.Granularity),sep="")
		validdoys <- switch(tim.gran, t4day=seq(from=1,to=365, by=4), t8day=seq(from=1,to=365, by=8), t16day=ifelse(rep(prod.info$Platform,23)=="Aqua", seq(from=9,to=365, by=16), seq(from=1,to=365, by=16)))
		
		if (is.null(doy)){
			doy <- validdoys				
		} else {
			doy <- doy[doy %in% validdoys]
		}
			
		for (yr in years) {
			for (d in doy){
				date.site <- paste(prod.site, "/",	format(as.Date(paste(yr, d), "%Y %j"), "%Y.%m.%d"), "/", sep="")				
				if (verbose) message("Acquiring file list in ", date.site, appendLF=TRUE)
				date.page <- readLines(curl(date.site))
				# TODO: If date page does not exist (Error 404) skip 
				acqdate <- paste("A",yr,sprintf("%03d",d),sep="")
				
				for (tl in tile){
					if (is.null(savepath)) tile.dir <- paste(pr,tl, sep="/") else tile.dir <- savepath
					tile.page <- date.page[grep(paste(pr, acqdate, tl, sep="."),date.page)]
					tile.page <- tile.page[-grep("BROWSE",tile.page)]
					link.st <- regexpr(paste(">",product,".*./",sep=""), tile.page)
					link.en <- regexpr("</a>", tile.page)
					tilefiles <- substr(tile.page, link.st+1,link.en-1)
					tile.site <- paste(date.site, tilefiles, sep="")
					
					if (length(tilefiles)==2){ 
						# extract filenames from html
						hdffile <- tilefiles[1]
						xmlfile <- tilefiles[2]
						
						if (file.exists(paste(tile.dir, hdffile, sep="/")) & skip.exists) {
							if (verbose) message(hdffile, " exists locally.", appendLF=TRUE)
							result <- c(result,paste(savepath,hdffile,sep="/"))
							next
							# File already present in local savepath	
						} 
						
						if (integrity) {							
							#download.file(tile.site[2], destfile=xmlfile, ...)
							xml <- rawToChar(curl_fetch_memory(tile.site[2],handle=curl_handle)$content)
							xml <-  unlist(strsplit(xml,"\n"))
							#xml <- unlist(strsplit(getURL(paste(product.site, xmlfile, sep="")),"\n"))
						}
						
						if (file.exists(paste("./", tile.dir, "/", hdffile, sep="")) & integrity){
							if (verbose) {
								message(hdffile, " exists locally.", appendLF=TRUE)
								message("Checking integrity...", appendLF=FALSE)
							}
							
							if(check.integrity(modis.hdf=paste(tile.dir, hdffile, sep="/"),xml=xml)) {
								message(" SUCCESS!", appendLF=TRUE)
								result <- c(result,paste(savepath,hdffile,sep="/"))
								next
							} else {
								#if (verbose) message("Downloading ", product.site, hdffile, appendLF=TRUE)		
								
								#hdf <- download.file(paste(product.site, hdffile, sep=""), destfile=paste(savepath,hdffile, sep="/"), method='internal', mode='wb',quiet=!verbose)								
								message("FAILED. Removing old file.", appendLF=TRUE)
								file.remove(paste(tile.dir, hdffile, sep="/")) 
							}
						} else if (file.exists(paste(tile.dir, hdffile, sep="/")) & !skip.exists){
							message("Removing old file.", appendLF=TRUE)
							file.remove(paste(savepath,hdffile, sep="/"))
						}
						
						
						# File not yet downloaded - attempt to get it!
						if (verbose) message("Downloading ", hdffile, appendLF=TRUE)
						hdf <- try(curl_fetch_disk(tile.site[1], path=paste(tile.dir, hdffile, sep="/"), handle=curl_handle))
						if(class(hdf)=="try-error"){
							break
						}
						#hdf <- download.file(tile.site[1], destfile=paste(tile.dir, hdffile, sep="/"), ...)
						
						# check integrity
						if (integrity){
							if (verbose) message("Checking integrity...", appendLF=FALSE)
							# Verify successful download
							if(check.integrity(modis.hdf=paste(tile.dir, hdffile, sep="/"),xml=xml)) {
								message(" SUCCESS!", appendLF=TRUE)
								result <- c(result,paste(tile.dir, hdffile, sep="/"))
								next
							} else {
								message("FAILED. Try to redownload later.", appendLF=TRUE)
								file.remove(paste(tile.dir, hdffile, sep="/")) 
							}
						}					
					} else {
						message(tile, " not found in ", prod.site, appendLF=TRUE) 
					}
					
				}			
			}
		}
	}
	
	return(result)
}


# TODO: Functions below may be deprecated
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

