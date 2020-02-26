# Title: IRRI-GIS Implementation of Xiao (2006) Rice Mapping Algorithm

#     -4 - Forest/Thick Vegetation 
#     -5 - Shrub
#     -6 - persistent water
# library(future.apply)
library(orysat)
library(manipulateR)

APPROX_BY_JULIANDAY = TRUE

maxlen.condition <- function(x){
  x <- which(x)
  return(maxConsecutive(x))
}

approx_int <- function(...){
  result <- approx(...)$y
  result <- round(result)
  return(result)
}

# Extract data
message("ORYSAT-", METHOD, ": Listing files.")
filelist.index <- dir(INDICES_DIR, pattern = paste0(PRODUCTS, ".", TILE, ".*.tif$"), full.names = TRUE)
inv.idxfiles <-  inventory.modis(filelist.index, modisinfo = c("product", "zone", "acqdate", "band"), file.ext = "tif")
inv.idxfiles <- inv.idxfiles[inv.idxfiles$acqdate %in% required.acqdates,]
inv.idxfiles <- inv.idxfiles[order(inv.idxfiles$band, inv.idxfiles$acqdate),]

# NAs in actualdate layers are water pixels
message("ORYSAT-", METHOD, ": Extracting Image date.")
stk.jday <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="actualdate"])
mat.jday <- raster::values(stk.jday) 
na.count <- rowSums(is.na(mat.jday))
if (exists("baseraster")) maskraster <- baseraster else maskraster <- raster(stk.jday)
maskraster[na.count==65] <- -8 #Ocean


save(mat.jday, file = "mat_jday.Rdata")
rm(mat.jday, na.count, stk.jday)
gc(reset=TRUE)

timest.extract <- Sys.time()
message("ORYSAT-", METHOD, ": Extracting NDVI.")
stk.ndvi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="NDVI"])
mat.ndvi <- raster::values(stk.ndvi)

# For masking out pixels with too many consecutive NAs (cloud or missing)
message("ORYSAT-", METHOD, ": Identifying Clouds.")
na.maxlen <- apply(is.na(mat.ndvi), 1, maxlen.condition)
maskraster[na.maxlen>8] <- -9 # Cloud
rm(na.maxlen, stk.ndvi)
# For masking out pixels with too many consecutive Snow or too many snow instances
# in case snw.maxlen is not good enough, the 2 conditions below can be explored
# It is possible that there's only snow at one time but rice can be planted some other time within the year
#
# snw.count <- rowSums(mat.ndvi == -32767) 
# nsnw.maxlen <- apply(mat.ndvi != -32767, 1, maxlen.condition)
# valid.maxlen <- apply(mat.ndvi != -3.2767, 2, maxlen.condition)
# Initial snow masking

pixels.toprocess <- which(is.na(maskraster[]))
mat.ndvi <- mat.ndvi[pixels.toprocess,]
mat.ndvi <- apply(mat.ndvi, 1, approx_int, x=1:ncol(mat.ndvi), xout=1:ncol(mat.ndvi)) # Output is a dataframe of time-series ndvi values where pixels are at the columns


message("ORYSAT-", METHOD, ": Identifying snow.")
snw.maxlen <- apply(mat.ndvi == -32767, 2, maxlen.condition)
maskraster[pixels.toprocess[snw.maxlen>3]] <- -7 # Snow

message("ORYSAT-", METHOD, ": Identifying forests.")
save(mat.ndvi, file="mat_ndvi.Rdata")
mat.ndvi <- mat.ndvi[grepl(YEAR, required.acqdates),]
forest.count <- colSums(mat.ndvi>=7000)
maskraster[pixels.toprocess[forest.count>=20]] <- -4 # Forests

save(mat.ndvi, file=paste0("mat_ndvi", YEAR, ".Rdata"))
rm(snw.maxlen, forest.count, mat.ndvi)


timest.lswi <- Sys.time()
message("ORYSAT-", METHOD, ": Extracting LSWI.")
stk.lswi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="LSWI"])
mat.lswi <- raster::values(stk.lswi)
mat.lswi <- mat.lswi[pixels.toprocess,]
mat.lswi <- apply(mat.lswi, 1, approx_int, x=1:ncol(mat.lswi), xout=1:ncol(mat.lswi))


# Caching full lswi series
save(mat.lswi, file="mat_lswi.Rdata")

message("ORYSAT-", METHOD, ": Identifying shrubs.")
mat.lswi <- mat.lswi[grepl(YEAR, required.acqdates),]
shrb.count <- colSums(mat.lswi<1000)
maskraster[pixels.toprocess[shrb.count==0]] <- -5
rm(shrb.count)

load(paste0("mat_ndvi", YEAR, ".Rdata"))

message("ORYSAT-", METHOD, ": Identifying permanent water.")
nh2o <- mat.ndvi>=1000 | mat.ndvi >= mat.lswi
nh2o.count <- apply(nh2o, 2, maxlen.condition)
maskraster[pixels.toprocess[nh2o.count<15]] <- -6
rm(nh2o, nh2o.count, mat.lswi, mat.ndvi)
gc(reset=TRUE)
potential.rice <- is.na(maskraster[pixels.toprocess])

load("mat_jday.Rdata")
mat.jday <- mat.jday[pixels.toprocess[potential.rice],]
mat.jday <- t(mat.jday)
mat.jday <- as.data.frame(mat.jday)

message("ORYSAT-", METHOD, ": Extracting EVI.")
stk.evi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="EVI"])
mat.evi <- raster::values(stk.evi)
mat.evi <- mat.evi[pixels.toprocess[potential.rice],]
mat.evi <- apply(mat.evi, 1, approx_int, x=1:ncol(mat.evi), xout=1:ncol(mat.evi))
mat.evi <- mat.evi/10000
mat.evi <- as.data.frame(mat.evi)

message("ORYSAT-", METHOD, ": Filling EVI Time Series.")
if(APPROX_BY_JULIANDAY){
  fill.evi <- mapply(approx, y=mat.evi, x=mat.jday,xout=data.frame(xout=as.Date(required.acqdates, "A%Y%j")))
  fill.evi <- unlist(fill.evi[seq(2, length(fill.evi), by=2)])
  fill.evi <- matrix(fill.evi, nrow=length(required.acqdates))  
  fill.evi <- round(fill.evi,4)
  fill.evi <- as.data.frame(fill.evi)
  rm(mat.evi)
} 

load("mat_ndvi.Rdata")
mat.ndvi <- mat.ndvi[, potential.rice]
mat.ndvi <- mat.ndvi/10000
mat.ndvi <- as.data.frame(mat.ndvi)

message("ORYSAT-", METHOD, ": Filling NDVI Time Series.")
if(APPROX_BY_JULIANDAY){
  fill.ndvi <- mapply(approx, y=mat.ndvi, x=mat.jday,xout=data.frame(xout=as.Date(required.acqdates, "A%Y%j")))
  fill.ndvi <- unlist(fill.ndvi[seq(2, length(fill.ndvi), by=2)])
  fill.ndvi <- matrix(fill.ndvi, nrow=length(required.acqdates))  
  fill.ndvi <- round(fill.ndvi,4)
  fill.ndvi <- as.data.frame(fill.ndvi)
  rm(mat.ndvi)
} 

load("mat_lswi.Rdata")
mat.lswi <- mat.lswi[, potential.rice]
mat.lswi <- mat.lswi/10000
mat.lswi <- as.data.frame(mat.lswi)

message("ORYSAT-", METHOD, ": Filling LSWI Time Series.")
if(APPROX_BY_JULIANDAY){
  fill.lswi <- mapply(approx, y=mat.lswi, x=mat.jday,xout=data.frame(xout=as.Date(required.acqdates, "A%Y%j")))
  fill.lswi <- unlist(fill.lswi[seq(2, length(fill.lswi), by=2)])
  fill.lswi <- matrix(fill.lswi, nrow=length(required.acqdates))  
  fill.lswi <- round(fill.lswi,4)
  fill.lswi <- as.data.frame(fill.lswi)
  rm(mat.lswi)
} 

rm(list=grep("stk.", ls(), value = TRUE))
gc(reset=TRUE)

riceraster <- raster(maskraster)
riceraster[pixels.toprocess[potential.rice]] <- mapply(FUN=rice.Xiao_v1, evi=fill.evi, lswi=fill.lswi, ndvi=fill.ndvi) #, evi.ricemax=0, evi.halfricemax=0) 

message("ORYSAT-", METHOD, ": Writing to disk.")
maskraster <- writeRaster(maskraster, paste0(paste("MASK", METHOD, TILE, YEAR, sep="_"),".tif"), datatype="INT4S", overwrite=TRUE)
riceraster <- writeRaster(riceraster, paste0(paste("RiceMap", METHOD, TILE, YEAR, ifelse(APPROX_BY_JULIANDAY,"ABJ","ACQ"), sep="_"),".tif"), datatype="INT4S", overwrite=TRUE)
timeen.proc <- timeen.rice <- Sys.time()
timedur.rice <- timeen.proc-timest.extract

message("ORYSAT-", METHOD, ": Done. (", round(timedur.rice,2), " ", attr(timedur.rice, "unit"), ") ")
