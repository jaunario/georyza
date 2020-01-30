# Author: Federico Filipponi
# Editor: Jorrel Khalil Aunario
# Consiglio Nazionale delle Ricerche (CNR-IREA) - Italy
# Date :  Dec 2012
# Version 1.0
# Licence GPL v3
# Maintainer: Mirco Boschetti <boschetti.m@irea.cnr.it>, Federico Filipponi <federico.filipponi@gmail.com>
# Description: Stack creation for evi, ndfi, noise indices needed in phenorice calculation using GDAL libraries

# normalizePath ruins links which have embedded permissions

assembleStack <- function(mpath=getwd(), startYear=startYear, tile=tile, verbose=TRUE, subset=TRUE) {
  
  # set folders for stack creation
  stackdir <- paste(mpath, startYear, "stack", sep="/")
  force.directories(stackdir,recursive=TRUE)
  # set vegetation indices directory
  vegindicesdir <- paste(mpath, startYear, "veg", sep="/")
  # set cloud directory
  noisedir <- paste(mpath, startYear, "noise", sep="/")
  # print message error if there are not 77 tif images in veg and noise folder
  if (length(list.files(path=vegindicesdir, pattern=c("evi.tif"), full.names=FALSE))!=77){
    stop("Folder ", vegindicesdir, " is not containing 77 evi images")
  }
  if (length(list.files(path=vegindicesdir, pattern=c("ndfi.tif"), full.names=FALSE))!=77){
    stop("Folder ", vegindicesdir, " is not containing 77 ndfi images")
    
  }
  if (length(list.files(path=noisedir, pattern=c("noise.tif"), full.names=FALSE))!=77){
    stop("Folder ", noisedir, " is not containing 77 noise images")

  }
  # set TIF directory
  #lstdir <- normalizePath(paste(mpath, tile, "TIF", sep="/"), winslash="/", mustWork=TRUE)

  # create output names path
  stackinfo <- paste("MOD09A1", startYear, tile, sep=".")
  
  # create multitemporal raster stack in ENVI bsq format
  
  
  if (!subset) {
    if (verbose) message("Creating evi index multitemporal stack for year ", startYear)
    evistack <- tiftostack(index="evi", vegindicesdir, outputstackdir=stackdir, stackinfo=stackinfo, ot="Int16", subset=FALSE)
    if (verbose) message("Creating ndfi index multitemporal stack for year ", startYear)
    ndfistack <- tiftostack(index="ndfi", vegindicesdir, outputstackdir=stackdir, stackinfo=stackinfo, ot="Int16", subset=FALSE)
    if (verbose) message("Creating MODIS noise multitemporal stack for year ", startYear)
    noisestack <- tiftostack(index="noise", noisedir, outputstackdir=stackdir, stackinfo=stackinfo, ot="Byte", subset=FALSE)
  } else {
    if (verbose) message("Creating evi index multitemporal stack for year ", startYear)
    evistack <- tiftostack(index="evi", vegindicesdir, outputstackdir=stackdir, stackinfo=stackinfo, ot="Int16", subset=TRUE)
    if (verbose) message("Creating ndfi index multitemporal stack for year ", startYear)
    ndfistack <- tiftostack(index="ndfi", vegindicesdir, outputstackdir=stackdir, stackinfo=stackinfo, ot="Int16", subset=TRUE)
    if (verbose) message("Creating MODIS noise multitemporal stack for year ", startYear)
    noisestack <- tiftostack(index="noise", noisedir, outputstackdir=stackdir, stackinfo=stackinfo, ot="Byte", subset=TRUE)
  }
}
