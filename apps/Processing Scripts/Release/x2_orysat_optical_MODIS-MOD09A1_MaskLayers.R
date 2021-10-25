# Basic Pre-Processing for MODIS MOD09A1 optical data rice mapping
# Use of MODIS Land COver product, MCD12Q1 allows the masking of features that are expected to be non-rice. 
# This reduces the processing time significantly as this will exclude the masked areas from the analysis  

# DIRECTORY SETTINGS
MODIS_HOME = "F:/MODIS"

WORKSPACE  = "E:/WORKSPACE/FLAGSHIP"
INDICES_DIR = paste0(WORKSPACE, "/indices")
IDXFILL_DIR = paste0(WORKSPACE, "/Filled")
CACHE_DIR    = paste0(WORKSPACE, "/cache")
   
# MODIS-SPECIFIC SETTINGS
TILE        = "h26v06"
PRODUCTS    = "MOD09A1"

LAYER_NAMES   = c("evi", "lswi", "ndvi", "mndwi", "ndfi")

# METHOD SETTINGS
YEAR       = 2018

# PROCESS OPTIONS
USE_MODISLC = TRUE
MODISLC_PROD = "MCD12Q1"
LC_LAYER     = "LC_Type1"
LC_FILTER    = c(1:5, # Forests
                 8,   # Savannas
                 13)  # Built-up

SAVE_MASKS    = FALSE
DO_SPATFILL   = TRUE # Perform spatial-filling using focal
XIAO_SNOWFILTER = FALSE

SETTINGS_FILE = ""

if(SETTINGS_FILE!="") source(SETTINGS_FILE)

message("ORYSAT: MODIS Preprocessing start.")
proc.st <- Sys.time()

message("ORYSAT: Preparing R environment.")


setwd(WORKSPACE)
if(!dir.exists(paste0(IDXFILL_DIR,"/", TILE))) dir.create(paste0(IDXFILL_DIR,"/", TILE), recursive = TRUE)
#if(SAVE_MASKS & !dir.exists(MASK_DIR)) dir.create(MASK_DIR)  

library(orysat)
library(manipulateR)
pylibs.start(python="~/../miniconda3", required=TRUE)

# List required images for the specified year and method
# TODO: Make a function to create this using a function
ACQDOYS <- seq(from=1,to=361, by=8)
required.acqdates <- paste("A",c(paste(YEAR-1,sprintf("%03g", ACQDOYS)[(length(ACQDOYS)-7):length(ACQDOYS)], sep=""),  #PREVIOUS YEAR
                                 paste(YEAR,sprintf("%03g", ACQDOYS), sep=""),           # CURRENT YEAR
                                 paste(YEAR+1,sprintf("%03g",ACQDOYS)[1:8], sep="")),   # SUCCEEDING YEAR
                           sep="")
date.acqdates <- as.Date(required.acqdates, "A%Y%j")

message("ORYSAT: Creating inventory required MODIS files in ", MODIS_HOME)
filelist.indices <- dir(paste(INDICES_DIR, TILE, sep="/"), pattern = paste0(PRODUCTS, ".", TILE, ".*.tif$"), full.names = TRUE)
inv.idxfiles <-  inventory.modis(filelist.indices, modisinfo = c("product", TILE, "acqdate", "layer"), file.ext = "tif")
inv.idxfiles <- inv.idxfiles[inv.idxfiles$acqdate %in% required.acqdates,]

date.acqdates <- date.acqdates[required.acqdates %in% unique(inv.idxfiles$acqdate)]
required.acqdates <- required.acqdates[required.acqdates %in% unique(inv.idxfiles$acqdate)]
YOI <- grep(YEAR, required.acqdates)

# CHeck if a snowing tile
shp.tileinfo <- shapefile(dir(MODIS_HOME, pattern="modis.*.shp", full.names = TRUE))
ADD_SNOW <- shp.tileinfo$has_snow[shp.tileinfo$tile==TILE]==1
rm(shp.tileinfo)

message("ORYSAT - Mask and Fill: Creating mask using DEM and SLOPE." )
rst.dem <- raster(dir(paste0(MODIS_HOME, "/dem"), pattern=paste0(TILE,".*.DEM.tif$"), full.names = TRUE, recursive = TRUE))
rst.slp <- raster(dir(paste0(MODIS_HOME, "/dem"), pattern=paste0(TILE,".*.SLOPE.tif$"), full.names = TRUE, recursive=TRUE))

rst.demmask <- (rst.dem<=2000 | rst.slp<=2)
idx.tomask <- which(rst.demmask[]==0)
rm(rst.dem, rst.slp, rst.demmask)
gc(reset=TRUE)

if(USE_MODISLC) {
  message("ORYSAT - Mask and Fill: Extracting MODIS Land Cover data." )
  mdata.lc <- modis.readHDF(dir(MODIS_HOME, pattern=paste(paste0("MCD12Q1.A",YEAR,"001.",TILE),"hdf",sep=".*."), recursive = TRUE, full.names = TRUE))
  lc_mask <- which(mdata.lc@imgvals[,LC_LAYER] %in% LC_FILTER)
  rm(mdata.lc)
}

timest.indices <- Sys.time()


message("ORYSAT - Mask and Fill: PROPER START." )
dat.timelog <- data.frame(index=toupper(LAYER_NAMES), st=Sys.time(), en=Sys.time())
# TODO: CREATE THIS ADDITIONAL MASK LAYERS USING STATE FLAGS FOR MORE FLEXIBILITY
# USE THIS MASKS ON 03-MaskAndFill stage

for(i in 1:length(LAYER_NAMES)){
  
  if(ADD_SNOW) {
    message("ORYSAT - Mask and Fill, ", LAYER_NAMES[i], ": FINDING Snow covered pixels.")
    stk.sf <- raster::stack(inv.idxfiles$filename[inv.idxfiles$layer=="state_flags"])
    mat.sf <- raster::values(stk.sf)
    mat.sf <- mat.sf[-idx.tomask, YOI]
    is.snow <- stateflags.snow(mat.sf)
    
    if(length(is.snow)>0) {
      is.snow <- matrix(is.snow, ncol = ncol(mat.sf))
      is.snow <- t(is.snow)
      cells <- 1:ncell(stk.sf)
      is.snow <- cbind(cells[!cells %in% idx.tomask], is.snow)
      saveRDS(is.snow, file=paste0(CACHE_DIR, "/MAT_BOOL_", paste(TILE, YEAR, "SNOW", "rds", sep=".")))
      
    }
    rm(stk.sf, mat.sf)
    gc(reset = TRUE)
  }
  
  dat.timelog$st[i] <- Sys.time()
  message("ORYSAT - Mask and Fill, ", LAYER_NAMES[i], ": Loading raster files.")
  stk.idx <- raster::stack(inv.idxfiles$filename[inv.idxfiles$layer==toupper(LAYER_NAMES[i])])
  
#  message("ORYSAT - Mask and Fill,", LAYER_NAMES[i], ": Removing pixels unsuitable for rice cultivation.")
#  stk.idx <- mask(stk.idx, rst.demmask)
  
  message("ORYSAT - Mask and Fill, ", LAYER_NAMES[i], ": Extracting ", LAYER_NAMES[i], "." )
  mat.idx <- raster::values(stk.idx)
  mat.idx[idx.tomask,] <- NA
  if(USE_MODISLC) mat.idx[lc_mask,] <- NA 
  # Count NA's within the year of data being filled
  if(!exists("idx.toprocess")){
    message("ORYSAT - Mask and Fill, ", LAYER_NAMES[i], ": Identifying candidate pixels for filling." )
    
    #can.interp <- function(x){return(sum(diff(x)>1)>0)} # If true, data has gaps (missing data in between valid data)
    
    mat.nna <- !is.na(mat.idx)  # Determine if value is missing
    nna.count <- rowSums(mat.nna[,YOI]) # Count how many are valid values within YOI
    idx.toprocess <- which(nna.count>8)
    #mat.nna <- mat.nna[idx.toprocess,]
    mat.dates <- mat.dates[idx.toprocess,]
    # nna.idx <- apply(mat.nna,1,which) # Get Indices of not NA 
    # rm(mat.nna)
    # data.patches <- lapply(nna.idx, consecutive.groups)
    # rm(nna.idx)
    # 
    # cnt.patch <- sapply(data.patches, length)
    # rm(data.patches)
    # #dif.idx <- lapply(nna.idx, diff)  # Determine if there are NAs in between not NA data points 
    # gc(reset=TRUE)
    # 
    #idx.toapprox  <- sapply(nna.idx, FUN=can.interp)
    saveRDS(idx.toprocess, file=paste0(CACHE_DIR, "/NUM_INDEX_", paste(TILE, YEAR, "IDX_TOPROC", "rds", sep=".")))
    
    #rst.masked[-idx.toprocess] <- 2 # Not enough data for interpolation and smoothing
    #rst.masked[nna.count<=1] <- 1 # Oceans, Deep water, High-elevation and High-Slope pixels
    #if(USE_MODISLC) rst.masked[lc_mask] <- 3 # Forests, Savannas, and Built-up Areas
    
    rm(nna.count, data.patches, cnt.patch)
    gc(reset = TRUE)
    
    mat.dates <- mat.dates[which((nna.count==length(YOI)) | ((cnt.patch>=1) & (nna.count>10))),]
    #mat.dates[rowSums(is.na(mat.dates))==ncol(mat.dates),] <- date.acqdates
    mat.dates <- t(mat.dates)
    mat.dates <- as.data.frame(mat.dates)
    
    
    
    message("ORYSAT - Mask and Fill: Extracting dates" )
    stk.dates <- raster::stack(inv.idxfiles$filename[inv.idxfiles$layer=="actualdate"])
    mat.dates <- raster::values(stk.dates)
    mat.dates[idx.tomask,] <- NA
    if(USE_MODISLC) mat.dates[lc_mask,] <- NA
    
  }
  message("ORYSAT - Mask and Fill, ", LAYER_NAMES[i], ": Preparing data for interpolation.")
  
  mat.idx <- mat.idx[idx.toprocess,]
  mat.idx <- t(mat.idx)
  mat.idx <- as.data.frame(mat.idx)
  message("ORYSAT - Mask and Fill, ", LAYER_NAMES[i], ": Iterpolating.")
  fill.idx <- mapply(approx_int, y=mat.idx, x=mat.dates, xout=data.frame(xout=date.acqdates))
  
  message("ORYSAT - Mask and Fill, ", LAYER_NAMES[i], ": Smoothing time-series.")
  smooth.idx <- apply(fill.idx, 2, safe.sg, n=11)
  smooth.idx <- smooth.idx[YOI,]

  message("ORYSAT - Mask and Fill, ", LAYER_NAMES[i], ": Freeing up memory.")
  rm(mat.idx, fill.idx)
  gc(reset=TRUE)
  smooth.idx <- round(smooth.idx,0)
  smooth.idx <- matrix(as.integer(smooth.idx), nrow = length(YOI), ncol=length(idx.toprocess))
  saveRDS(smooth.idx, file=paste0(CACHE_DIR, "/MAT_SMOOTH_", paste(TILE, YEAR, toupper(LAYER_NAMES[i]), "rds", sep=".")))
  
  for (j in 1:nrow(smooth.idx)){
    fname <- paste0(IDXFILL_DIR,"/", TILE, "/", paste(PRODUCTS,TILE,required.acqdates[YOI[j]],toupper(LAYER_NAMES[i]),"tif",sep="."))
    #if(file.exists(fname)) next
    message("ORYSAT - Mask and Fill, ", LAYER_NAMES[i], ": Saving filled ", required.acqdates[YOI[j]], " to disk.")
    rst.idx <- raster(rst.masked)
    rst.idx[idx.toprocess] <- smooth.idx[j,]
    rst.idx <- writeRaster(rst.idx,filename = fname, options="COMPRESS=LZW", overwrite=TRUE)
  }

  rm(smooth.idx)
  dat.timelog$en[i] <- Sys.time()

}

#toremove <- ls()[!ls() %in% c("MODIS_HOME", "required.acqdates", "INDICES_DIR", "METHOD", "PRODUCTS", "TILE", "YEAR")]
rm(list=ls())
.rs.restartR()
#source('~/GitHub/georyza/apps/Processing Scripts/orysat_Xiao.R')
