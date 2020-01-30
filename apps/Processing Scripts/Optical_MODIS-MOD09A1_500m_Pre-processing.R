# Basic Pre-Processing for MODIS MOD09A1 optical data rice mapping

# DIRECTORY SETTINGS
MODIS_HOME = "C:/Data/MODIS/v006"
WORKSPACE  = "C:/WORKSPACE/Orysat"
INDICES_DIR = paste0(WORKSPACE, "/indices-filled")
MASK_DIR    = paste0(WORKSPACE, "/mask")
   
# MODIS-SPECIFIC SETTINGS
TILE        = "h26v06"
PRODUCTS    = "MOD09A1"
LCLU_PROD   = "MCD12Q1"

LAYERS        = c(1:6,12,13)
LAYER_NAMES   = c("red", "nir", "blue", "green", "nir2","swir","state_500m", "julianday")
LAYERS_TOFILL = 1:3

# METHOD SETTINGS
METHOD     = "xiao-v1"
YEAR       = 2018

# PROCESS OPTIONS
SAVE_MASKS = FALSE
DO_FILL    = TRUE # Perform spatial-filling using focal
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
python.start()

# List required images for the specified year and method
# TODO: Make a function to create this using a function
ACQDOYS <- seq(from=1,to=361, by=8)
required.acqdates <- paste("A",c(paste(YEAR-1,sprintf("%03g", ACQDOYS)[(length(ACQDOYS)-7):length(ACQDOYS)], sep=""),  #PREVIOUS YEAR
                                 paste(YEAR,sprintf("%03g", ACQDOYS), sep=""),           # CURRENT YEAR
                                 paste(YEAR+1,sprintf("%03g",ACQDOYS)[1:11], sep="")),   # SUCCEEDING YEAR
                           sep="")

indices <- c("ndvi", "evi", "lswi", "ndsi")

if(LCLU_PROD!=""){
  mdata.LCLU <- modis.readHDF(dir(paste0(MODIS_HOME,"/", LCLU_PROD),pattern = paste(LCLU_PROD, TILE, "hdf$", sep=".*."), recursive = TRUE, full.names = TRUE))
  baseraster <- raster(mdata.LCLU@extent, ncol=mdata.LCLU@ncols, nrow=mdata.LCLU@nrows, crs=mdata.LCLU@projection)
  maskable.features <-  which(mdata.LCLU@imgvals$LC_Type1 %in% c(1:7, 13, 15, 17))
  pixels.cropland <- which(!mdata.LCLU@imgvals$LC_Type1 %in% c(1:7, 13, 15, 17))
} 

message("ORYSAT: Creating inventory required MODIS files in ", MODIS_HOME)
filelist.modis <- dir(paste(MODIS_HOME, PRODUCTS, TILE, sep="/"), pattern = paste0(PRODUCTS, ".*.", TILE, ".*.hdf$"), recursive = TRUE, full.names = TRUE)
inv.hdffiles <-  inventory.modis(filelist.modis, modisinfo = c("product", "acqdate", "zone", "version", "proddate"), file.ext = "hdf")
inv.hdffiles <- inv.hdffiles[inv.hdffiles$acqdate %in% required.acqdates,]


ind.st <- Sys.time()
for(i in 1:nrow(inv.hdffiles)){
  st <- Sys.time()
  message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i], ": Reading ", basename(inv.hdffiles$filename[i]))
  mdata.mod09 <- modis.readHDF(inv.hdffiles$filename[i], layer = LAYERS)
  colnames(mdata.mod09@imgvals) <- LAYER_NAMES
  if(exists("pixels.cropland")){
    mdata.mod09@imgvals <- mdata.mod09@imgvals[pixels.cropland,]
  } else {
    pixels.cropland <- 1:nrow(mdata@imgvals)
  }

  message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i],": Computing indices.")
  # Compute indices required for the analysis
  mdata.indices <- modis.data(mdata.mod09)
  mdata.indices@imgvals <- modis.compute(mdata.mod09, funlist = indices)
  
  # Includes cloud, internal cloud, blue-based cloud and snow
  fillable.cloud <- (stateflags.cloud(mdata.mod09@imgvals$state_500m) + stateflags.internalCloud(mdata.mod09@imgvals$state_500m))*0.5 +  xiaoflags.cloud(mdata.mod09@imgvals$blue, scale=1)*2  
  fillable.snow <- xiaoflags.snow(mdata.indices@imgvals$ndsi,mdata.mod09@imgvals$nir) + stateflags.snow(mdata.mod09@imgvals$state_500m)   
  
  tofill <- unique(which(fillable.cloud>.5), which(fillable.snow>0))
  
  # Fill suspected cloudy and snow pixels
  if(length(tofill)>0){

    message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i],": Removing data from cloudy and snowy pixels.")
    mdata.indices@imgvals[tofill, LAYERS_TOFILL] <- NA
    
    if(DO_FILL){
      for(j in 1:length(LAYERS_TOFILL)){
        message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i],": Trying to fill ", indices[j])
        thislayer <- raster(baseraster)
        thislayer[pixels.cropland] <- mdata.indices@imgvals[,j]
        thislayer <- focal(thislayer, matrix(rep(1,9),ncol=3), fun=focal.mean, NAonly=TRUE)
        mdata.indices@imgvals[tofill,j] <- thislayer[pixels.cropland[tofill]]
        rm(thislayer)
      }
    }
  }
  mdata.indices@imgvals <- mdata.indices@imgvals[,LAYERS_TOFILL]
  mdata.indices@imgvals$julianday <- mdata.mod09@imgvals$julianday
  gc(reset = TRUE)
  
  # Save indices as GeoTiff raster files,
  fname.base <- paste(mdata.indices@product, mdata.indices@zone, mdata.indices@acqdate, sep=".")
  for(j in 1:ncol(mdata.indices@imgvals)){
    message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i],": Writing ", colnames(mdata.indices@imgvals)[j], " to disk.") 
    thislayer <- raster(baseraster)
    if(colnames(mdata.indices@imgvals)[j]!="julianday"){
      thislayer[pixels.cropland] <- round(mdata.indices@imgvals[,j],4)*10000
    } else thislayer[pixels.cropland] <- mdata.indices@imgvals[,j]
    thislayer <- writeRaster(thislayer, 
      filename = paste0(INDICES_DIR,"/", fname.base, ".", colnames(mdata.indices@imgvals)[j],".tif"),
      options = c("COMPRESS=LZW"), overwrite=TRUE, datatype="INT2S")
    rm(thislayer)
  }
  
  if(SAVE_MASKS){
    message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i], ": Saving cloud layer.")
    MASK_DIR <- "./mask"
    if(!dir.exists(MASK_DIR)) dir.create(MASK_DIR)
    
    cld <- raster(baseraster)
    cld[pixels.cropland] <- fillable.cloud
    cld[cld[]<=0.5] <- NA
    cld[cld[]>0.5] <- 1
    cld <- writeRaster(cld,
      filename = paste0(MASK_DIR,"/", fname.base, ".CLOUD.tif"),
      options = c("COMPRESS=LZW"), overwrite=TRUE, datatype="INT1U")
    
    if(length(fillable.snow)>0){
      snw <- raster(baseraster)
      snw[pixels.cropland] <- fillable.snow
      snw[snw[]<=0] <- NA
      snw[snw[]>0] <- 1
      snw <- writeRaster(snw,
                         filename = paste0(MASK_DIR,"/", fname.base, ".SNOW.tif"),
                         options = c("COMPRESS=LZW"), overwrite=TRUE, datatype="INT1U")
      
    }
    
  }  
  rm(mdata.indices, mdata.mod09)
  gc(reset=TRUE)
  en <- Sys.time()
  aproctime <- en-st
  message("ORYSAT-", METHOD,",", inv.hdffiles$acqdate[i], ": Done. (", round(aproctime,2), " ", attr(aproctime, "unit"), ", ", i, " of ", nrow(inv.hdffiles), "). ")

}
ind.en <- Sys.time()



# xx[] <- mdata.LCLU@imgvals$LC_Type1
# rr <- gg <- bb <- raster(baseraster)
# rr[pixels.cropland] <- mdata.mod09@imgvals$red
# gg[pixels.cropland] <- mdata.mod09@imgvals$green
# bb[pixels.cropland] <- mdata.mod09@imgvals$blue
# rgb <- stack(rr,gg,bb)
# 
# plotRGB(rgb, stretch="lin")
# zz <- raster(xx)
# zz[pixels.cropland] <- 1
# 
# plotRGB(rgb, stretch="lin")
# plot(cld, add=T)
# 
# others <- raster(baseraster)
# others[pixels.cropland] <- mdata.indices@imgvals$NDVI
# #others[!others[]] <- NA
# #others[others[]>=1] <- 1
# 
# plotRGB(rgb, stretch="lin")
# plot(others)
