# Author: Federico Filipponi, Francesco Nutini
# Editor: Jorrel Khalil Aunario
# Consiglio Nazionale delle Ricerche (CNR-IREA) - Italy
# Date :  Gen 2013
# Version 1.2
# Licence GPL v3
# Maintainer: Mirco Boschetti <boschetti.m@irea.cnr.it>, Federico Filipponi <federico.filipponi@gmail.com>
# Description: Run phenorice algorithm on indices stack and save results to disk

# require libraries

doyFrom2000 <- function(x){
	x <- as.Date(as.character(x),"%Y%j")
	return(as.numeric(x-as.Date("2000-1-1")+1))
}

PhenoRiceRun <- function(sub, startYear, tile, mpath=Sys.getenv("MODIS_LOCAL_DIR"), evismoothopt=FALSE, verbose=TRUE, debug=FALSE){
	
	if(debug){
    	secdat <- vector()
    	errdat <- vector()    
    	subst <- Sys.time()
  	}
  
  # Get Static Mask 
  subs <- paste("_sub", sub, sep="")
  maskfile <- dir(paste(mpath,"masks",sep="/"), pattern=paste(tile, subs, ".tif$", sep=""), full.names=TRUE)
  # import mask file if exists in the mask path
  if(length(maskfile)==1){      
  	# Import mask values
	maskraster <- raster(maskfile)		
    maskvector <- getValues(maskraster)
  } else {
    if(verbose) message("Mask file found for tile ",tile, " subtile ", sub, " is ", length(maskfile))
    maskraster <- vector()
  } 
  
  # sum mask values
  if(maskfile==""){
    continua=1
  } else {
    masksum <- sum(maskvector, na.rm=TRUE)
    continua <- ifelse(masksum==0, 0, 1)
  }
  
  # skip PhenoRice process if all pixels are masked
  if(continua==0){
    maxseason <- 0
  } else {
    
	  # there are pixel to be processed. PhenoRice process start  
	    
	  # set stack path
	  stackdir <- paste(mpath, tile, "Pheno", startYear, "stack", sep="/")
	  evistackdir <- paste(stackdir, paste("MOD09A1", ".", startYear, ".", tile, ".evi.stack", subs, ".tif", sep=""), sep="/")
	  ndfistackdir <- paste(stackdir, paste("MOD09A1", ".", startYear, ".", tile, ".ndfi.stack", subs,".tif", sep=""), sep="/")
	  noisestackdir <- paste(stackdir, paste("MOD09A1", ".", startYear, ".", tile, ".noise.stack", subs, ".tif", sep=""), sep="/")
	  lststackdir <- paste(stackdir, paste("MOD11A2", ".", startYear, ".", tile, ".lst.stack", subs, ".tif", sep=""), sep="/")
	  
	  # create layer stack index for each stack
	  evistack <- stack(evistackdir)
	  ndfistack <- stack(ndfistackdir)
	  noisestack <- stack(noisestackdir)
	  
	  
	  # set parameters using raster stack
	  t <- nlayers(evistack)
	  
	  # Create modis.data object for raster output
	  mdata <- new("modis.data") ###@ create an empty "modis.data" object, which function is defined in "modisdata.R"
	  mdata@product <- "MOD09A1" ###@ assign metadata to mdata object retrieved from "files" object
	  mdata@acqdate <-  as.character(startYear)
	  mdata@zone <- tile
	  mdata@version <- "005"
	  mdata@proddate <- as.character(startYear)
	  mdata@projection <- projection(evistack)
	  mdata@extent <- extent(evistack)
	  mdata@ncols <- ncol(evistack)
	  mdata@nrows <- nrow(evistack)
	  
	  # import tile cubes
	  if(verbose) message("Importing indices for year ", startYear)
	  evicubo <- getValues(evistack)
	  ndficubo <- getValues(ndfistack)
	  noisecubo <- getValues(noisestack)
	  rm(evistack,ndfistack,noisestack)
	  gc(verbose=FALSE,reset=TRUE)
	  
	  # replace NA with numeric values
	  ndficubo[is.na(ndficubo)] <- 0
	  evicubo[is.na(evicubo)] <- -10000
	  noisecubo[is.na(noisecubo)] <- 1
	  noisecubo[which(noisecubo==0)] <- 1
	  
	  # set pixel number to create empty objects for storing results
	  if(maskfile!=""){
		if (length(maskvector)!=length(evicubo[,1])) {
	      maskfile <- ""
	      npixel <- length(evicubo[,1])
	    } else {
	      npixel <- length(maskvector)
	    }
	  } else {
	    npixel <- length(evicubo[,1])
	  }
	  
	  # import lst stack if present in stack folder
	  if (file.exists(lststackdir)) {
	    lstcubo <- getValues(stack(lststackdir))
	        
	    # update ndfi index stack using lst values
	    lstcubo <- !(lstcubo>0 & lstcubo<14408) # value 14408 corresponds to land surface temperature of 15 ?C
		
		ndficubo <- lstcubo*ndficubo
	  }
	  
	  # remove LST stack
	  rm(lstcubo)
	  gc(verbose=FALSE,reset=TRUE)
	    
	  # create empty arrays for storing results
	  # TODO: Do not allocate. Just use 'c' (combine function)
	  ricemap <- as.vector(as.integer(rep(0, npixel)))
	  season <- as.vector(as.integer(rep(0, npixel)))
	  riceseason <- as.vector(as.integer(rep(0, npixel)))
	  riceseason1 <- as.vector(as.integer(rep(0, npixel)))
	  riceseason2 <- as.vector(as.integer(rep(0, npixel)))
	  riceseason3 <- as.vector(as.integer(rep(0, npixel)))
	  seasonofrice <- as.vector(as.integer(rep(0, npixel)))	  
	  startseason <- array(0, dim=c(npixel,3))
	  endseason <- array(0, dim=c(npixel,3))
	  absolutemin <- array(0, dim=c(npixel,3))
	  absolutemax <- array(0, dim=c(npixel,3))
	  
	  if (evismoothopt) {
	    # create empty array for storing smoothed evi
	    minievi <- array(0, dim=c(npixel,3))
	    massievi <- array(0, dim=c(npixel,3))
	    deltaevi <- array(0, dim=c(npixel,3))
	    evismooth <- array(0, dim=c(npixel,t))
		colnames(minievi) <- c(paste(subs, ".flood.evi.", 1:3, sep=""))
		colnames(massievi) <- c(paste(subs, ".peak.evi.", 1:3, sep=""))
		colnames(deltaevi) <- c(paste(subs, ".delta.evi.", 1:3, sep=""))
		
	    # duplicate mdata for evi smooth output
	    edata <- mdata
	  }
	
	  gc(verbose=FALSE,reset=TRUE)
	  
	  # create empty vectors for mask update
	  if(verbose) message("Creating mask for tile ", tile)
	
	  ncb <- noisecubo>1
	  emn <- evicubo*ncb>4500
	  maskeviforest <- !rowSums(emn)>=(rowSums(ncb)/2)
	  maskevi <- !apply(evicubo, 1, max)<2000
	  #maskndfi <- !apply(ndficubo, 1, max)<0
	  
	  #nmn <- ndficubo*ncb>0
	  #maskndfi_3 <- !rowSums(nmn)>=(rowSums(ncb)/1.1)
	  
	  rm(ncb, emn)
	  gc(verbose=FALSE,reset=TRUE)
	   
	  if(maskfile!=""){
	    maskdef <- as.integer(maskvector*maskevi*maskeviforest)
		rm(maskvector)
	  } else {
	    maskdef <- maskevi*maskeviforest
	  }
	  
	  # clean R workspace performing garbage collection
		rm(maskeviforest,maskevi,maskfile)
		gc(verbose=FALSE,reset=TRUE)
	  
	  ## calculate phenorice parameters
	  
	  # extract time series one pixel per time
	  # rice identification stage (time it - start)
	  if (debug) {
		  stage.st <- Sys.time()
		  if (verbose) message("Computing phenorice for ", sum(maskdef), " pixel. Started at: ", stage.st)
	  }
	  analyze <- which(maskdef==1)
	  gc(verbose=FALSE,reset=TRUE)
	  
	  for (p in analyze) {
	      
	      # run phenorice calculation
		  phenoriceout <- try(PhenoRice(evi=evicubo[p,], ndfi=ndficubo[p,], noise=noisecubo[p,], evismoothopt=evismoothopt))
		  
		  if(class(phenoriceout)=="try-error"){
			errdat <- rbind(errdat, c(as.character(sub), as.character(p), phenoriceout, format(p.st, "%c"), format(Sys.time(), "%c")))  
			next
		  }
		  
	      # skip phenorice output storage in result matrix if rice is not present 
	      if (phenoriceout[[1]]$ricemap[1]==0) {
	        next
	      } else {
	        # write phenorice outputs in result matrix
	        ricemap[p] <- round(phenoriceout[[1]]$ricemap[1],0)
	        season[p] <- round(phenoriceout[[1]]$season[1],0)
			# FROM 1.3.7
			riceseason[p] <- round(phenoriceout[[1]]$riceseason[1],0)
			riceseason1[p] <- round(phenoriceout[[1]]$riceseason1[1],0)
			riceseason2[p] <- round(phenoriceout[[1]]$riceseason2[1],0)
			riceseason3[p] <- round(phenoriceout[[1]]$riceseason3[1],0)
			seasonofrice[p] <- phenoriceout[[1]]$ricedet[1]
			
	        startseason[p,1:3] <- phenoriceout[[1]]$start
	        endseason[p,1:3] <- phenoriceout[[1]]$end
	        absolutemin[p,1:3] <- phenoriceout[[1]]$min
	        absolutemax[p,1:3] <- phenoriceout[[1]]$max
	        
	        if (evismoothopt) {
	          minievi[p,1:3] <- phenoriceout[[2]]$mini_evi
	          massievi[p,1:3] <- phenoriceout[[2]]$massi_evi
	          deltaevi[p,1:3] <- phenoriceout[[2]]$delta_evi
	          evismooth[p,] <- phenoriceout[[3]]
	        }
	      }
	  	rm(phenoriceout)
		gc(verbose=FALSE, reset=TRUE)
  }
  # rice identification (log it - end)  
  if (debug){
	  stage.en <- Sys.time()
	  secdat <- rbind(secdat, c(as.character(sub), "identification", format(stage.st, "%c"), format(stage.en, "%c")))
  } 
  
  # CLEANUP
  rm(evicubo,ndficubo,noisecubo)
  gc(verbose=FALSE,reset=TRUE)
  
  # Output formatting and storage (time it - start)
  if (debug) stage.st <- Sys.time()
  # evaluate maximum number of season
  maxseason <- max(season)
  
  
  # store ricemap and season results in data frame
  ricemap <- data.frame(ricemap,season,riceseason,seasonofrice,riceseason1,riceseason2,riceseason3)
  colnames(ricemap) <- c(paste(subs,".phenoricemap", sep=""),paste(subs, ".cropseasons", sep=""),paste(subs, ".riceseasons", sep=""),paste(subs, ".seasonofrice", sep=""),paste(subs, ".ricemap.s1", sep=""),paste(subs, ".ricemap.s2", sep=""),paste(subs, ".ricemap.s3", sep=""))
  
  rm(season,riceseason,seasonofrice,riceseason1,riceseason2,riceseason3)
  gc(verbose=FALSE, reset=TRUE)
  
  # convert positions to YYYYDoY
  absoluteminA <- ifelse(absolutemin==0, NA, ifelse(absolutemin<20, (startYear-1)*1000+(absolutemin+26)*8+1, ifelse(absolutemin>65, (startYear+1)*1000+(absolutemin-66)*8+1, (startYear)*1000+(absolutemin-20)*8+1)))
  absoluteminA[is.na(absoluteminA)] <- 0
  colnames(absoluteminA) <- c(paste(subs, ".flood.A.", 1:3, sep=""))
  # convert positions to MODIS composite number
  absoluteminB <- ifelse(absolutemin==0, NA, ifelse(absolutemin<20, absolutemin+27, ifelse(absolutemin>65, absolutemin+1935, absolutemin+981)))
  absoluteminB[is.na(absoluteminB)] <- 0
  colnames(absoluteminB) <- c(paste(subs, ".flood.B.", 1:3, sep=""))
  # convert positions to progressive Julian date starting 1 January 2000
  absoluteminC <- ifelse(absoluteminA=="NA", NA, ifelse(absoluteminA==0, NA, as.numeric(doyFrom2000(absoluteminA))))
  absoluteminC[is.na(absoluteminC)] <- 0
  colnames(absoluteminC) <- c(paste(subs, ".flood.C.", 1:3, sep=""))

  # convert positions to YYYYDoY
  startseasonA <- ifelse(startseason==0, NA, ifelse(startseason<20, (startYear-1)*1000+(startseason+26)*8+1, ifelse(startseason>65, (startYear+1)*1000+(startseason-66)*8+1, (startYear)*1000+(startseason-20)*8+1)))
  startseasonA[is.na(startseasonA)] <- 0
  colnames(startseasonA) <- c(paste(subs, ".SoS.A.", 1:3, sep=""))
  # convert positions to MODIS composite number
  startseasonB <- ifelse(startseason==0, NA, ifelse(startseason<20, startseason+27, ifelse(startseason>65, startseason+1935, startseason+981)))
  startseasonB[is.na(startseasonB)] <- 0
  colnames(startseasonB) <- c(paste(subs, ".SoS.B.", 1:3, sep=""))
  # convert positions to progressive Julian date starting 1 January 2000
  startseasonC <- ifelse(startseasonA=="NA", NA, ifelse(startseasonA==0, NA, as.numeric(doyFrom2000(startseasonA))))
  startseasonC[is.na(startseasonC)] <- 0
  colnames(startseasonC) <- c(paste(subs, ".SoS.C.", 1:3, sep=""))

	
  # convert positions to YYYYDoY
  absolutemaxA <- ifelse(absolutemax==0, NA, ifelse(absolutemax<20, (startYear-1)*1000+(absolutemax+26)*8+1, ifelse(absolutemax>65, (startYear+1)*1000+(absolutemax-66)*8+1, (startYear)*1000+(absolutemax-20)*8+1)))
  absolutemaxA[is.na(absolutemaxA)] <- 0
  colnames(absolutemaxA) <- c(paste(subs, ".peak.A.", 1:3, sep=""))
  # convert positions to MODIS composite number
  absolutemaxB <- ifelse(absolutemax==0, NA, ifelse(absolutemax<20, absolutemax+27, ifelse(absolutemax>65, absolutemax+1935, absolutemax+981)))
  absolutemaxB[is.na(absolutemaxB)] <- 0
  colnames(absolutemaxB) <- c(paste(subs, ".peak.B.", 1:3, sep=""))
  # convert positions to progressive Julian date starting 1 January 2000
  absolutemaxC <- ifelse(absolutemaxA=="NA", NA, ifelse(absolutemaxA==0, NA, as.numeric(doyFrom2000(absolutemaxA))))
  absolutemaxC[is.na(absolutemaxC)] <- 0
  colnames(absolutemaxC) <- c(paste(subs, ".peak.C.", 1:3, sep=""))

  
  # convert positions to YYYYDoY
  endseasonA <- ifelse(endseason==0, NA, ifelse(endseason<20, (startYear-1)*1000+(endseason+26)*8+1, ifelse(endseason>65, (startYear+1)*1000+(endseason-66)*8+1, (startYear)*1000+(endseason-20)*8+1)))
  endseasonA[is.na(endseasonA)] <- 0
  colnames(endseasonA) <- c(paste(subs, ".EoS.A.", 1:3, sep=""))
  # convert positions to MODIS composite number
  endseasonB <- ifelse(endseason==0, NA, ifelse(endseason<20, endseason+27, ifelse(endseason>65, endseason+1935, endseason+981)))
  endseasonB[is.na(endseasonB)] <- 0
  colnames(endseasonB) <- c(paste(subs, ".EoS.B.", 1:3, sep=""))
  # convert positions to progressive Julian date starting 1 January 2000
  endseasonC <- ifelse(endseasonA=="NA", NA, ifelse(endseasonA==0, NA, as.numeric(doyFrom2000(endseasonA))))
  endseasonC[is.na(endseasonC)] <- 0
  colnames(endseasonC) <- c(paste(subs, ".EoS.C.", 1:3, sep=""))

  
  # convert positions to YYYYDoY
  endseason <- ifelse(endseason==0, NA, ifelse(endseason<23, (startYear-1)*1000+(endseason+23)*8+1, ifelse(endseason>68, (startYear+1)*1000+(endseason-69)*8+1, (startYear)*1000+(endseason-23)*8+1)))
  colnames(endseason) <- c(paste(subs, ".EoS.", 1:3, sep=""))

  rm(absolutemin,startseason,absolutemax,endseason)
  gc(verbose=FALSE,reset=TRUE)
  
  if (evismoothopt) {
    # store phenorice results in modis.data object
	mdata@imgvals <- data.frame(ricemap,absoluteminA,absoluteminB,absoluteminC,startseasonA,startseasonB,startseasonC,absolutemaxA,absolutemaxB,absolutemaxC,endseasonA,endseasonB,endseasonC,minievi,massievi,deltaevi)
	rm(ricemap,absoluteminA,absoluteminB,absoluteminC,startseasonA,startseasonB,startseasonC,absolutemaxA,absolutemaxB,absolutemaxC,endseasonA,endseasonB,endseasonC,minievi,massievi,deltaevi)
  } else {
	mdata@imgvals <- data.frame(ricemap,absoluteminA,absoluteminB,absoluteminC,startseasonA,startseasonB,startseasonC,absolutemaxA,absolutemaxB,absolutemaxC,endseasonA,endseasonB,endseasonC)
	rm(ricemap,absoluteminA,absoluteminB,absoluteminC,startseasonA,startseasonB,startseasonC,absolutemaxA,absolutemaxB,absolutemaxC,endseasonA,endseasonB,endseasonC)
  }
  gc(verbose=FALSE,reset=TRUE)
  
  # Output formatting and storing (log it - end)  
  if (debug){
	  stage.en <- Sys.time()
	  secdat <- rbind(secdat, c(as.character(sub), "format and store", format(stage.st, "%c"), format(stage.en, "%c")))
  } 
  
  # Write to disk (time it - start)
  if (debug) stage.st <- Sys.time()
  
  # save phenorice results to disk
  outdir <- stackdir <- paste(mpath, tile, "Pheno", startYear, "subset", sep="/")
  force.directories(outdir,recursive=TRUE)
  modisricemap <- modis.brick(mdata, process="", intlayers=1:ncol(mdata@imgvals),writeto=outdir, options="COMPRESS=LZW", overwrite=TRUE)
  rm(mdata,modisricemap)
  gc(verbose=FALSE,reset=TRUE)
  
  if (evismoothopt) {
    # save evismooth results to disk
    edata@imgvals <- as.data.frame(evismooth)
    colnames(edata@imgvals) <- paste(".evi.smooth.",sprintf("%02d",1:77),sep="") 			
    outesdir <- paste(outdir, "evismooth", paste("evi_smooth", subs, sep=""), sep="/")
    force.directories(outesdir,recursive=TRUE)
    modisricemap <- modis.brick(edata, process="es", intlayers=1:ncol(edata@imgvals), writeto=outesdir, fltNA=-32768, options="COMPRESS=LZW", overwrite=TRUE)
	rm(edata,modisricemap)
	gc(verbose=FALSE,reset=TRUE)	
  }  
  gc(verbose=FALSE,reset=TRUE)  
  }
  
  # Write to disk (log it - end)  
  if (debug){
	  fname <- paste(paste(tile, sub, startYear, format(subst, "%H-%M"), format(Sys.time(), "%H-%M"), Sys.time()-subst, length(analyze), "v2", sep="_"),".txt", sep="")
	  stage.en <- Sys.time()
	  secdat <- rbind(secdat, c(as.character(sub), "disk save", format(stage.st, "%c"), format(stage.en, "%c")))
	  write.csv(secdat, fname, row.names=FALSE)
	  if(length(errdat)!=0){
		  errfname <- paste(paste(tile, sub, startYear, format(subst, "%H-%M"), format(Sys.time(), "%H-%M"), Sys.time()-subst, length(analyze), "err", sep="_"),".txt", sep="")
		  write.csv(errdat, errfname, row.names=FALSE)
	  }	  	  
  }
  return(maxseason) 
}

#mtrace(PhenoRiceRun)
#tmp <- try(PhenoRiceRun(sub=1, startYear=year, tile=tile, evismoothopt=smoothen,debug=TRUE),silent=TRUE)
#go(127)


