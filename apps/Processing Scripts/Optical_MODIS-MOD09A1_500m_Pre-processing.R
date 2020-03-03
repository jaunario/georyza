# Basic Pre-Processing for MODIS MOD09A1 optical data rice mapping
# Use of MODIS Land COver product, MCD12Q1 allows the masking of features that are expected to be non-rice. 
# This reduces the processing time significantly as this will exclude the masked areas from the analysis  


# DIRECTORY SETTINGS
MODIS_HOME = "C:/Data/MODIS"

WORKSPACE  = "C:/WORKSPACE/Orysat"
INDICES_DIR = paste0(WORKSPACE, "/indices")
MASK_DIR    = paste0(WORKSPACE, "/mask")
   
# MODIS-SPECIFIC SETTINGS
TILE        = "h26v06"
PRODUCTS    = "MOD09A1"
PROD_VER    = 6

LAYERS        = c(1:7,12,13)
LAYER_NAMES   = c("red", "nir", "blue", "green", "nir2","swir", "swir1", "state_500m", "julianday")
LAYERS_TOFILL = 1:4

# METHOD SETTINGS
METHOD     = "xiao-v1"
YEAR       = 2015

# PROCESS OPTIONS
SAVE_MASKS    = TRUE
DO_SPATFILL   = TRUE # Perform spatial-filling using focal
SETTINGS_FILE = ""


if(SETTINGS_FILE!="") source(SETTINGS_FILE)

focal.mean <- function(x){
  if(sum(is.na(x))==length(x)) ans <- NA else ans <- mean(x, na.rm=TRUE)
  return(ans)
}

message("ORYSAT: RiceMapper ", METHOD, " start.")
proc.st <- Sys.time()

message("ORYSAT: Preparing R environment.")


setwd(WORKSPACE)
if(!dir.exists(INDICES_DIR)) dir.create(INDICES_DIR)
if(SAVE_MASKS & !dir.exists(MASK_DIR)) dir.create(MASK_DIR)  

library(orysat)
library(manipulateR)
pylibs.start(python="~/../miniconda3", required=TRUE)

# List required images for the specified year and method
# TODO: Make a function to create this using a function
ACQDOYS <- seq(from=1,to=361, by=8)
required.acqdates <- paste("A",c(paste(YEAR-1,sprintf("%03g", ACQDOYS)[(length(ACQDOYS)-7):length(ACQDOYS)], sep=""),  #PREVIOUS YEAR
                                 paste(YEAR,sprintf("%03g", ACQDOYS), sep=""),           # CURRENT YEAR
                                 paste(YEAR+1,sprintf("%03g",ACQDOYS)[1:11], sep="")),   # SUCCEEDING YEAR
                           sep="")

indices <- c("ndvi", "evi", "lswi", "ndsi", "mndwi")

message("ORYSAT: Creating inventory required MODIS files in ", MODIS_HOME)
filelist.modis <- dir(paste(MODIS_HOME, paste0("v", PROD_VER),PRODUCTS, TILE, sep="/"), pattern = paste0(PRODUCTS, ".*.", TILE, ".*.hdf$"), recursive = TRUE, full.names = TRUE)
inv.hdffiles <-  inventory.modis(filelist.modis, modisinfo = c("product", "acqdate", "zone", "version", "proddate"), file.ext = "hdf")
inv.hdffiles <- inv.hdffiles[inv.hdffiles$acqdate %in% required.acqdates,]

rst.dem <- raster(dir(paste0(MODIS_HOME, "/dem"), pattern=TILE, full.names = TRUE, recursive=TRUE))
rst.slp <- raster(dir(paste0(MODIS_HOME, "/slp"), pattern=TILE, full.names = TRUE, recursive=TRUE))

rst.demmask <- (rst.dem<=2000 | rst.slp<=2)
rst.demmask[rst.demmask[]==0] <- NA

timest.indices <- Sys.time()
for(i in 1:nrow(inv.hdffiles)){
  st <- Sys.time()
  message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i], ": Reading ", basename(inv.hdffiles$filename[i]))
  mdata.mod09 <- modis.readHDF(inv.hdffiles$filename[i], layer = LAYERS)
  if(nrow(mdata.mod09@imgvals)==0) stop("not read") #mdata.mod09 <- modis.readHDF(inv.hdffiles$filename[i], layer = LAYERS)
  colnames(mdata.mod09@imgvals) <- LAYER_NAMES
  
  if(!exists("baseraster")){
    baseraster <- raster(mdata.mod09@extent, ncol=mdata.mod09@ncols, nrow=mdata.mod09@nrows, crs=mdata.mod09@projection)
  } 
  
  message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i],": Computing indices.")
  # Compute indices required for the analysis
  mdata.indices <- modis.data(mdata.mod09)
  mdata.indices@imgvals <- modis.compute(mdata.mod09, funlist = indices)
  mdata.indices@imgvals$julianday <- mdata.mod09@imgvals$julianday
  
  # Includes cloud, internal cloud, blue-based cloud
  fillable.cloud <- (stateflags.cloud(mdata.mod09@imgvals$state_500m) + stateflags.internalCloud(mdata.mod09@imgvals$state_500m))*0.5 +  xiaoflags.cloud(mdata.mod09@imgvals$blue, scale=1)*2  
  
  # Fill suspected cloudy and snow pixels
  message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i],": Removing data from cloud pixels.")
  mdata.indices@imgvals[which(fillable.cloud>.5), LAYERS_TOFILL] <- NA
  
  flag.snow <- xiaoflags.snow(mdata.indices@imgvals$NDSI,mdata.mod09@imgvals$nir) + stateflags.snow(mdata.mod09@imgvals$state_500m)   
  
  fname.base <- paste(mdata.indices@product, mdata.indices@zone, mdata.indices@acqdate, sep=".")
  for(j in 1:ncol(mdata.indices@imgvals)){
    fname <- paste0(INDICES_DIR,"/", fname.base, ".", colnames(mdata.indices@imgvals)[j],".tif")
    if(file.exists(fname)) next
    thislayer <- baseraster
    
    if(colnames(mdata.indices@imgvals)[j]=="julianday"){
      actualdates <- as.vector(dateFromDoy(mdata.indices@imgvals[,j], year=as.numeric(substr(mdata.mod09@acqdate,2,5))))
      values(thislayer) <- actualdates
      colnames(mdata.indices@imgvals)[j] <- "actualdate"
    } else {
      values(thislayer) <- mdata.indices@imgvals[,j]  
      
      if(DO_SPATFILL){
        message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i],": Trying to fill ", indices[j])
        thislayer <- focal(thislayer, matrix(rep(1,9),ncol=3), fun=focal.mean, NAonly=TRUE)
      }
      
      message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i],": Flagging snow.")
      thislayer[which(flag.snow>0 & !is.na(thislayer[]))] <- -3.2767
      
      if(maxValue(thislayer)>3.2767 | minValue(thislayer) < -3.2767) stop("CHANGE DATATYPE")
      thislayer <- round(thislayer,4)*10000
    }
    
    thislayer <- mask(thislayer,rst.demmask)
    message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i],": Writing ", colnames(mdata.indices@imgvals)[j], " to disk.")
    thislayer <- writeRaster(thislayer, 
                               filename = fname,
                               options = c("COMPRESS=LZW"), datatype="INT2S")    
    rm(thislayer)
  }
  rm(mdata.indices,mdata.mod09)
  gc(reset = TRUE)
  
  if(SAVE_MASKS){
    message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i], ": Saving mask.")
    MASK_DIR <- "./mask"
    if(!dir.exists(MASK_DIR)) dir.create(MASK_DIR)
    
    rst.mask <- (rst.dem>2000 | rst.slp>2)            # High elevation/ steep slope
    rst.mask[fillable.cloud>0.5 & rst.mask[]==0] <- 2 # Cloud 
    rst.mask[flag.snow>0 & rst.mask[]==0] <- 3        # Snow
    fname <- paste0(MASK_DIR,"/", fname.base, ".PREPROCMASK.tif")
    if(!file.exists(fname)){
      rst.mask <- writeRaster(rst.mask, filename = fname,
                         options = c("COMPRESS=LZW"), overwrite=TRUE, datatype="INT1U")
    }
  }  
  rm(rst.mask)
  gc(reset=TRUE)
  en <- Sys.time()
  aproctime <- en-st
  message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i], ": Done. (", round(aproctime,2), " ", attr(aproctime, "unit"), ", ", i, " of ", nrow(inv.hdffiles), "). ")

}

timeen.indices <- Sys.time()
timeen.indices-timest.indices
toremove <- ls()[!ls() %in% c("MODIS_HOME", "required.acqdates", "INDICES_DIR", "METHOD", "PRODUCTS", "TILE", "YEAR")]
rm(list=toremove)
#source('~/GitHub/georyza/apps/Processing Scripts/orysat_Xiao.R')
