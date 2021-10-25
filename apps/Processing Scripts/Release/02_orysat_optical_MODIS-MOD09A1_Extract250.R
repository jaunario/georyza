# Basic Pre-Processing for MODIS MOD09A1 optical data rice mapping
# Use of MODIS Land COver product, MCD12Q1 allows the masking of features that are expected to be non-rice. 
# This reduces the processing time significantly as this will exclude the masked areas from the analysis  

# DIRECTO/RY SETTINGS
MODIS_HOME = "F:/MODIS"                           # Local directory where MODIS HDF files were saved
WORKSPACE  = "E:/WORKSPACE/FLAGSHIP"              # A google dive folder for saving crop.curve outputs and rice maps
PROCESS_DIR = "E:/WORKSPACE/Demeter Initiative"   # Local directory for saving intermediate files

# MODIS SETTINGS
TILE        = "h25v07"
PRODUCTS    =  "MOD13Q1" #"MOD09A1"
PROD_VER    = 6

# TODO: Automate based on products using mod.prods for now then using lpdaac.search later

# 500m Products
# LAYERS        = c(1:7,12,13)
# LAYER_NAMES   = c("red", "nir", "blue", "green", "nir2","swir", "swir2", "state_500m", "julianday")

# 250m products
# LAYERS        = c(1:7,11,12)
# LAYER_NAMES   = c("ndvi", "evi", "qaflag", "red", "nir", "blue","mir", "doy", "reliability")

SQA_FILL     = c("cloud", "internalCloud")

# METHOD SETTINGS
YEAR       = 2018

# PROCESS OPTIONS
#INDICES       = c("ndvi", "evi", "lswi", "ndsi", "mndwi", "ndfi")
#INDICES       = c("lswi", "ndsi", "mndwi", "ndfi")
# LAYERS_TOFILL = 1:6


# SYSTEM SETTINGS
SETTINGS_FILE = "00_orysat_RiceMap-Flagship_SETTINGS-Template.R"
INDICES_DIR   = paste0(WORKSPACE, "/250m/INDICES/", TILE)
SKIP_EXISTING = FALSE

# PROCESS PROPER
proc.st <- Sys.time()
message("ORYSAT: MODIS Preprocessing start.")

if(SETTINGS_FILE!="") source(SETTINGS_FILE)
message("ORYSAT: Preparing R environment.")

setwd(WORKSPACE)
if(SETTINGS_FILE!="" && file.exists(SETTINGS_FILE)) source(SETTINGS_FILE)
if(!dir.exists(INDICES_DIR)) dir.create(INDICES_DIR, recursive = TRUE)

library(orysat)
library(manipulateR)
pylibs.start(python="~/../miniconda3", required=TRUE)

modprods <- read.csv(system.file("satPRODUCTS/modis.PRODUCTS.ref.csv", package="orysat"), stringsAsFactors=FALSE)

# List required images for the specified year and method
# TODO: Make a function to create this using a function
info.prod <- modprods[modprods$ShortName==PRODUCTS,]
ACQDOYS <-  switch(info.prod$Temporal.Granularity, Yearly=1, Daily=1:365, "16 day"={if(info.prod$Platform=="Terra") seq(1,361,by=16) else seq(9,361,by=16)}, "8 day"=seq(1,361,by=8),NA)
required.acqdates <- paste0("A",YEAR,sprintf("%03g", ACQDOYS))
  

message("ORYSAT: Creating inventory required MODIS files in ", MODIS_HOME)
filelist.modis <- dir(paste(MODIS_HOME, paste0("v", PROD_VER),PRODUCTS, TILE, sep="/"), pattern = paste0(PRODUCTS, ".*.", TILE, ".*.hdf$"), recursive = TRUE, full.names = TRUE)
inv.hdffiles <-  inventory.modis(filelist.modis, modisinfo = c("product", "acqdate", "zone", "version", "proddate"), file.ext = "hdf")
inv.hdffiles <- inv.hdffiles[inv.hdffiles$acqdate %in% required.acqdates,]

timest.INDICES <- Sys.time()
for(i in 1:nrow(inv.hdffiles)){
  st <- Sys.time()
  message("ORYSAT | Extract-Compute,", inv.hdffiles$acqdate[i], ": Reading ", basename(inv.hdffiles$filename[i]))
  mdata.mod09 <- modis.readHDF(inv.hdffiles$filename[i], layer = LAYERS)
  if(nrow(mdata.mod09@imgvals)==0) stop("not read") #mdata.mod09 <- modis.readHDF(inv.hdffiles$filename[i], layer = LAYERS)
  colnames(mdata.mod09@imgvals) <- LAYER_NAMES
  
  if(!exists("baseraster")){
    baseraster <- raster(mdata.mod09@extent, ncol=mdata.mod09@ncols, nrow=mdata.mod09@nrows, crs=mdata.mod09@projection)
  } 
  
  message("ORYSAT | Extract-Compute,", inv.hdffiles$acqdate[i],": Computing INDICES.")
  # Compute INDICES required for the analysis
  mdata.INDICES <- modis.data(mdata.mod09)
  mdata.INDICES@imgvals <- modis.compute(mdata.mod09, funlist = INDICES)
  mdata.INDICES@imgvals$actualdate <- as.vector(dateFromDoy(mdata.mod09@imgvals$julianday, year=as.numeric(substr(mdata.mod09@acqdate,2,5))))
  mdata.INDICES@imgvals$state_flags <- mdata.mod09@imgvals$state_500m
  
  # SQA-based filters, TO be removed and later tried to be filled 
  for(j in 1:length(SQA_FILL)){
    m <- SQA_FILL[j]
    bool.fill <- eval(parse(text=paste0("stateflags.",m, "(mdata.mod09@imgvals$state_500m)")))
    if(j==1) tofill <- bool.fill else tofill <- tofill | bool.fill
  }
  
  water.mask <- stateflags.water(mdata.mod09@imgvals$state_500m, exclude = c(2,3)) 
  tofill <- tofill | water.mask   
  # Fill filtered pixels
  message("ORYSAT | Extract-Compute,", inv.hdffiles$acqdate[i],": Removing data from cloud pixels.")
  mdata.INDICES@imgvals[which(tofill), LAYERS_TOFILL] <- NA
  
  fname.base <- paste(mdata.INDICES@product, mdata.INDICES@zone, mdata.INDICES@acqdate, sep=".")
  for(j in 1:ncol(mdata.INDICES@imgvals)){
    fname <- paste0(INDICES_DIR,"/", fname.base, ".", colnames(mdata.INDICES@imgvals)[j],".tif")
    if(SKIP_EXISTING & file.exists(fname)) next
    
    thislayer <- baseraster
    values(thislayer) <- mdata.INDICES@imgvals[,j]
    
    
    if(!colnames(mdata.INDICES@imgvals)[j] %in% c("actualdate", "state_flags")){
      #if(maxValue(thislayer)>3.2767 | minValue(thislayer) < -3.2767) stop("CHANGE DATATYPE")
      thislayer <- round(thislayer,4)*10000
    }
    
    message("ORYSAT | Extract-Compute,", inv.hdffiles$acqdate[i],": Writing ", colnames(mdata.INDICES@imgvals)[j], " to disk.")
    
    if(colnames(mdata.INDICES@imgvals)[j]=="state_flags"){
      thislayer <- writeRaster(thislayer, 
                               filename = fname,
                               options = "COMPRESS=LZW", datatype="INT2U",
                               overwrite=TRUE)    
    } else {
      thislayer <- writeRaster(thislayer, 
                               filename = fname,
                               options = "COMPRESS=LZW", datatype="INT2S",
                               overwrite=TRUE)
    }
    rm(thislayer)
  }
  rm(mdata.INDICES,mdata.mod09)
  gc(reset = TRUE)
  
  en <- Sys.time()
  aproctime <- en-st
  message("ORYSAT | Extract-Compute,", inv.hdffiles$acqdate[i], ": Done. (", round(aproctime,2), " ", attr(aproctime, "unit"), ", ", i, " of ", nrow(inv.hdffiles), "). ")

}

timeen.INDICES <- Sys.time()
timeen.INDICES-timest.INDICES
#toremove <- ls()[!ls() %in% c("MODIS_HOME", "required.acqdates", "INDICES_DIR", "METHOD", "PRODUCTS", "TILE", "YEAR")]
rm(list=ls())
.rs.restartR()
#source('~/GitHub/georyza/apps/Processing Scripts/orysat_Xiao.R')
