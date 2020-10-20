# Title: Tools for processing MODIS HDF files
# Description:
# Author: Jorrel Khalil Aunario 
# International Rice Research Institute
# Date : 21 Jan 2019
# Version 0.2
# Licence GPL v3
# 

PROJ.MODIS <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

modis.acqdoys <- function(){
  
}

modis.productinfo <- function(product){
  # TODO: Create Function that utilizes online LPDAAC Catalog
  modprods <- read.csv(system.file("satproducts/modis.products.ref.csv", package="orysat"), stringsAsFactors=FALSE)
  return(modprods[modprods$ShortName==product,])
}

modis.readHDF <- function(hdffile, layer=1, verbose=TRUE, ...){
  if(!exists("pygdal")| !exists("pyosr")) {
    # #stop("Run python start")
    # use_python("~/../AppData/Local/r-miniconda/python")
    # list.libs=list(pygdal ="gdal", pyosr="osr")
    # for(i in 1:length(list.libs)){
    #   assign(names(list.libs)[i], import(list.libs[[i]]))
    # }
    pylibs.start(...)
  }  
  
  info.hdffile <- unlist(strsplit(basename(hdffile),"\\."))
  m <- new("modis.data")
  m@product=info.hdffile[1]
  m@acqdate=info.hdffile[2]
  m@zone=info.hdffile[3]
  m@version=info.hdffile[4] 
  m@proddate=info.hdffile[5]
  m@cached <- FALSE
  m@imgvals <- data.frame()
  obj.gdal <- pygdal$Open(hdffile, pygdal$GA_ReadOnly)
  obj.subdatasetlist <- obj.gdal$GetSubDatasets()
  layernames <- vector()
  if(length(layer)<1 || is.na(layer)|is.null(layer)) layer <- 1:length(obj.subdatasetlist)
  for(i in layer){
    sdsinfo <- unlist(strsplit(obj.subdatasetlist[[i]][[1]],":"))
    layernames <- c(layernames, sdsinfo[length(sdsinfo)])
    obj.sds <- pygdal$Open(obj.subdatasetlist[[i]][[1]], pygdal$GA_ReadOnly)
    fillvalue <- as.numeric(obj.sds$GetMetadataItem("_FillValue"))
    scalefactor <- as.numeric(obj.sds$GetMetadataItem("scale_factor"))
    if(length(scalefactor)==0) scalefactor <- 1
    if(log10(scalefactor)>0) scalefactor <- 1/scalefactor
    
    if(m@projection==""){
      proj.crs <- pyosr$SpatialReference()
      proj.crs$ImportFromWkt(obj.sds$GetProjection())

      obj.origin <- obj.sds$GetGeoTransform()
      minX <- obj.origin[[1]]
      resX <- obj.origin[[2]]
      maxX <- minX + resX * obj.sds$RasterXSize
      
      maxY <- obj.origin[[4]]
      resY <- obj.origin[[6]]
      minY <- maxY + resY * obj.sds$RasterYSize
      m@extent <- extent(minX, maxX, minY, maxY)
      m@projection <- proj.crs$ExportToProj4()
      m@ncols <- obj.sds$RasterXSize
      m@nrows <- obj.sds$RasterYSize
      rm(proj.crs, minX, maxX, resX, minY, maxY, resY)
    }
    values <- as.numeric(t(obj.sds$ReadAsArray()))*scalefactor
    values[values==fillvalue] <- NA 
    if(nrow(m@imgvals)==0){ 

      m@imgvals <- data.frame(values)
      
    } else m@imgvals[,layernames[length(layernames)]] <- values
    colnames(m@imgvals) <- layernames
    rm(obj.sds)
    gc(reset=TRUE)
  }
  return(m)
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

modis.compute <- function(modis, funlist){
  result <- layers <- vector()
  for (i in 1:length(funlist)){
    argnames <- names(formals(get(funlist[i])))
    if (sum(argnames %in% colnames(modis@imgvals))!=length(argnames)) {
      warning("Missing arguments. Did not run ", funlist[i])
    } else {
      layers <- c(layers, toupper(funlist[i]))
      result <- cbind(result, do.call(funlist[i],as.list(modis@imgvals[,argnames])))
    }
  }
  colnames(result) <- layers
  return(as.data.frame(result))
}

modis.brick <- function(mdata, process=NULL, intlayers=NULL, writeto=NULL, intNA=-15, fltNA=-9999.0, format="GTiff", skipx=FALSE, ...){
  if(!file.exists(writeto)) dir.create(writeto, recursive=TRUE)
  mraster <- raster(modis@extent, ncols=modis@ncols, nrows=modis@nrows, crs=modis@projection)
  
  if(is.character(writeto)){                
    fname <- gsub("\\.\\.", "\\.", paste(modis@product, modis@acqdate, modis@zone, modis@version, modis@proddate, colnames(modis@imgvals), process, formatExt(format), sep="."))        
    fname <- as.character(paste(writeto,fname,sep="/"))          
  } else mbrick <- brick(mraster)
  
  for(i in 1:ncol(modis@imgvals)){
    if(skipx & file.exists(fname[i])) next
    mraster <- setValues(mraster,modis@imgvals[,i])
    
    if(!is.null(intlayers) & is.numeric(intlayers)){
      dataType(mraster) <- ifelse(i %in% intlayers,"INT1U","FLT4S")                
    } else if(!is.null(intlayers)){
      message("Ignoring intlayers. Should be 1:ncol(modis@imgvals) instead of", paste(intlayers, collapse=","))
    }
    
    if(is.character(writeto)) {            
      if (dataType(mraster)== "INT1U"){
        writeRaster(mraster,filename=fname[i], format=format, NAflag=intNA, ...)
      } else {
        writeRaster(mraster,filename=fname[i], format=format, NAflag=fltNA, ...)
      }
    } else  mbrick <- addLayer(mbrick,mraster)
  }
  if(is.character(writeto)) mbrick <- TRUE else mbrick@layernames <- colnames(modis@imgvals)
  return(mbrick)
}
