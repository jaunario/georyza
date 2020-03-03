# Title: IRRI-GIS Implementation of Xiao (2006) Rice Mapping Algorithm

#     -4 - Forest/Thick Vegetation 
#     -5 - Shrub
#     -6 - persistent water
# library(future.apply)
APPROX_BY_JULIANDAY = TRUE
CONSMISSING_THRES   = 5
library(orysat)
library(manipulateR)


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

timest.extract <- Sys.time()
message("ORYSAT-", METHOD, ": Extracting NDVI.")
stk.ndvi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="NDVI"])
mat.ndvi <- raster::values(stk.ndvi)

# For masking out pixels with too many consecutive NAs (cloud or missing)
message("ORYSAT-", METHOD, ": Identifying pixels with sufficient data.")
na.maxlen <- apply(is.na(mat.ndvi), 1, maxlen.condition)

pixels.toprocess <- which(na.maxlen<=CONSMISSING_THRES)

rst.lc <- raster(stk.ndvi) 
rst.lc[na.maxlen>CONSMISSING_THRES & na.maxlen<65] <- -11 # Cloud

mat.ndvi <- mat.ndvi[pixels.toprocess,]
mat.ndvi <- t(mat.ndvi)
mat.ndvi <- as.data.frame(mat.ndvi)

message("ORYSAT-", METHOD, ": Identifying forests.")
forest.count <- colSums(mat.ndvi[grepl(YEAR, required.acqdates),]>=7000, na.rm = TRUE)
rst.lc[pixels.toprocess[forest.count>=20]] <- -4 # Forests

message("ORYSAT-", METHOD, ": Identifying snow.")
# For masking out pixels with too many consecutive Snow or too many snow instances
# in case snw.maxlen is not good enough, the 2 conditions below can be explored
# It is possible that there's only snow at one time but rice can be planted some other time within the year
#
# snw.count <- rowSums(mat.ndvi == -32767) 
nsnw.maxlen <- apply(mat.ndvi != -32767, 2, maxlen.condition)
snw.maxlen <- apply(mat.ndvi == -32767, 2, maxlen.condition)
# valid.maxlen <- apply(mat.ndvi != -3.2767, 2, maxlen.condition)
# Initial snow masking
rst.lc[pixels.toprocess[snw.maxlen>3 & nsnw.maxlen<10]] <- -7 # Snow


save(mat.ndvi, file=paste0("mat_ndvi", YEAR, ".Rdata"))
timest.lswi <- Sys.time()
message("ORYSAT-", METHOD, ": Extracting LSWI.")
stk.lswi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="LSWI"])
mat.lswi <- raster::values(stk.lswi)
mat.lswi <- mat.lswi[pixels.toprocess,]
mat.lswi <- t(mat.lswi)
mat.lswi <- as.data.frame(mat.lswi)

message("ORYSAT-", METHOD, ": Identifying shrubs.")
shrb.count <- colSums(mat.lswi[grepl(YEAR, required.acqdates),]<1000, na.rm = TRUE)
rst.lc[pixels.toprocess[shrb.count==0 & is.na(rst.lc[pixels.toprocess])]] <- -5

message("ORYSAT-", METHOD, ": Identifying permanent water.")
nh2o <- mat.ndvi>=1000 | mat.ndvi >= mat.lswi
nh2o.maxlen <- apply(nh2o, 2, maxlen.condition)
rst.lc[pixels.toprocess[nh2o.maxlen<15 & is.na(rst.lc[pixels.toprocess])]] <- -6
rm(nh2o, nh2o.maxlen, shrb.count, forest.count, na.maxlen)
gc(reset=TRUE)
potential.rice <- is.na(rst.lc[pixels.toprocess])

mat.ndvi <- mat.ndvi[, potential.rice]
mat.lswi <- mat.lswi[, potential.rice]

if(APPROX_BY_JULIANDAY){
  message("ORYSAT-", METHOD, ": Extracting pixel acqdate date.")
  stk.jday <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="actualdate"])
  mat.jday <- raster::values(stk.jday) 
  mat.jday <- mat.jday[pixels.toprocess[potential.rice],]
  mat.jday <- t(mat.jday)
  mat.jday <- as.data.frame(mat.jday)
}

if(APPROX_BY_JULIANDAY){
  fill.ndvi <- mapply(approx, y=mat.ndvi, x=mat.jday, xout=data.frame(xout=as.Date(required.acqdates, "A%Y%j")))
  fill.ndvi <- unlist(fill.ndvi[seq(2, length(fill.ndvi), by=2)])
  fill.ndvi <- matrix(fill.ndvi, nrow=length(required.acqdates))
  
} else {
  # Linear interpolation of time series
  fill.ndvi <- apply(mat.ndvi, 2, approx_int, x=1:ncol(mat.ndvi), xout=1:ncol(mat.ndvi)) # Output is a dataframe of time-series ndvi values where pixels are at the columns
}
rm(mat.ndvi)
fill.ndvi <- fill.ndvi/10000
fill.ndvi <- round(fill.ndvi,4)
fill.ndvi <- as.data.frame(fill.ndvi)

if(APPROX_BY_JULIANDAY){
  fill.lswi <- mapply(approx, y=mat.lswi, x=mat.jday, xout=data.frame(xout=as.Date(required.acqdates, "A%Y%j")))
  fill.lswi <- unlist(fill.lswi[seq(2, length(fill.lswi), by=2)])
  fill.lswi <- matrix(fill.lswi, nrow=length(required.acqdates))
  
} else {
  # Linear interpolation of time series
  fill.lswi <- apply(mat.lswi, 1, approx_int, x=1:ncol(mat.lswi), xout=1:ncol(mat.lswi)) # Output is a dataframe of time-series lswi values where pixels are at the columns
}
rm(mat.lswi)
fill.lswi <- fill.lswi/10000
fill.lswi <- round(fill.lswi,4)
fill.lswi <- as.data.frame(fill.lswi)

message("ORYSAT-", METHOD, ": Extracting EVI.")
stk.evi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="EVI"])
mat.evi <- raster::values(stk.evi)
mat.evi <- mat.evi[pixels.toprocess[potential.rice],]
mat.evi <- t(mat.evi)
mat.evi <- as.data.frame(mat.evi)

if(APPROX_BY_JULIANDAY){
  fill.evi <- mapply(approx, y=mat.evi, x=mat.jday, xout=data.frame(xout=as.Date(required.acqdates, "A%Y%j")))
  fill.evi <- unlist(fill.evi[seq(2, length(fill.evi), by=2)])
  fill.evi <- matrix(fill.evi, nrow=length(required.acqdates))
  
} else {
  # Linear interpolation of time series
  fill.evi <- apply(mat.evi, 1, approx_int, x=1:ncol(mat.evi), xout=1:ncol(mat.evi)) # Output is a dataframe of time-series evi values where pixels are at the columns
}
rm(mat.evi)
fill.evi <- fill.evi/10000
fill.evi <- round(fill.evi,4)
fill.evi <- as.data.frame(fill.evi)

rst.rice <- raster(rst.lc)
rst.rice[pixels.toprocess[potential.rice]] <- mapply(FUN=rice.Xiao_v1, evi=fill.evi, lswi=fill.lswi, ndvi=fill.ndvi) #, evi.ricemax=0, evi.halfricemax=0) 

message("ORYSAT-", METHOD, ": Writing to disk.")
rst.lc <- writeRaster(rst.lc, paste0(paste("LANDCOVER", METHOD, TILE, YEAR, sep="_"),".tif"), datatype="INT4S", overwrite=TRUE)
rst.rice <- writeRaster(rst.rice, paste0(paste("RiceMapFin", METHOD, TILE, YEAR, ifelse(APPROX_BY_JULIANDAY,"ABJ","ACQ"), sep="_"),".tif"), datatype="INT4S", overwrite=TRUE)
timeen.proc <- timeen.rice <- Sys.time()
timedur.rice <- timeen.proc-timest.extract

message("ORYSAT-", METHOD, ": Done. (", round(timedur.rice,2), " ", attr(timedur.rice, "unit"), ") ")
