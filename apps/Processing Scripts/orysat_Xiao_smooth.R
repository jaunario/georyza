# Title: IRRI-GIS Implementation of Xiao (2006) Rice Mapping Algorithm

#     -4 - Forest/Thick Vegetation 
#     -5 - Shrub
#     -6 - persistent water
# library(future.apply)
WORKSPACE  = "E:/WORKSPACE/Orysat"
INDICES_DIR = paste0(WORKSPACE, "/mosaic/bgd/indices")

APPROX_BY_JULIANDAY = TRUE
CONSMISSING_THRES   = 5

METHOD     = "Xiao-V1"
YEAR       = 2016
TILE       = "BGD"
PRODUCTS    = "SMOOTH"

APPROX_BY_JULIANDAY = TRUE
CONSMISSING_THRES   = 5

library(orysat)
library(manipulateR)

maxlen.condition <- function(x){
  x <- which(x)
  return(maxConsecutive(x))
}

setwd(WORKSPACE)

ACQDOYS <- seq(from=1,to=361, by=8)
required.acqdates <- paste("A",c(paste(YEAR-1,sprintf("%03g", ACQDOYS)[(length(ACQDOYS)-7):length(ACQDOYS)], sep=""),  #PREVIOUS YEAR
                                 paste(YEAR,sprintf("%03g", ACQDOYS), sep=""),           # CURRENT YEAR
                                 paste(YEAR+1,sprintf("%03g",ACQDOYS)[1:11], sep="")),   # SUCCEEDING YEAR
                           sep="")

# Extract data
message("ORYSAT-", METHOD, ": Listing files.")
filelist.index <- dir(INDICES_DIR, pattern = paste0(PRODUCTS, ".", TILE, ".*.tif$"), full.names = TRUE)
inv.idxfiles <-  inventory.modis(filelist.index, modisinfo = c("product", "zone", "acqdate", "band", "res", "proj"), file.ext = "tif", sep="_")
inv.idxfiles <- inv.idxfiles[inv.idxfiles$acqdate %in% required.acqdates,]
inv.idxfiles <- inv.idxfiles[order(inv.idxfiles$band, inv.idxfiles$acqdate),]

timest.extract <- Sys.time()
message("ORYSAT-", METHOD, ": Extracting NDVI.")
stk.ndvi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="NDVI"])
mat.ndvi <- raster::values(stk.ndvi)

# For masking out pixels with too many consecutive NAs (cloud or missing)
message("ORYSAT-", METHOD, ": Identifying pixels with sufficient data.")
pixels.toprocess <- which(rowSums(!is.na(mat.ndvi))>0)
mat.ndvi <- mat.ndvi[pixels.toprocess,]
mat.ndvi <- t(mat.ndvi)
mat.ndvi <- as.data.frame(mat.ndvi)

timest.lswi <- Sys.time()
message("ORYSAT-", METHOD, ": Extracting LSWI.")
stk.lswi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="LSWI"])
mat.lswi <- raster::values(stk.lswi)
mat.lswi <- mat.lswi[pixels.toprocess,]
mat.lswi <- t(mat.lswi)
mat.lswi <- as.data.frame(mat.lswi)


message("ORYSAT-", METHOD, ": Extracting EVI.")
stk.evi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="EVI"])
mat.evi <- raster::values(stk.evi)
mat.evi <- mat.evi[pixels.toprocess,]
mat.evi <- t(mat.evi)
mat.evi <- as.data.frame(mat.evi)
mat.evi <- mat.evi/10000
mat.ndvi <- mat.ndvi/10000
mat.lswi <- mat.lswi/10000

rst.rice <- raster(stk.evi)
rst.rice[pixels.toprocess] <- mapply(FUN=rice.Xiao_v1, evi=mat.evi, lswi=mat.lswi, ndvi=mat.ndvi) #, evi.ricemax=0, evi.halfricemax=0) 


vals.rice <- values(rst.rice)
cells.rice <- which(!is.na(vals.rice) & vals.rice>0)
vals.rice <- vals.rice[cells.rice]

season1 <- as.numeric(substr(vals.rice, 2,3))
season2 <- as.numeric(substr(vals.rice,5,6))
season3 <- as.numeric(substr(vals.rice,8,9))
dat.planting <- cbind(rep(cells.rice,3),c(season1,season2,season3))
dat.planting <- na.omit(dat.planting)
colnames(dat.planting) <- c("cell","acq.index")
dat.planting <- as.data.frame(dat.planting)
dat.planting$acqdate <- required.acqdates[dat.planting$acq.index]
dat.planting$date <- as.Date(dat.planting$acqdate, "A%Y%j")
dat.planting$month <- manipulateR::monthFromDate(dat.planting$date)

seasons <- list(aman=6:9, boro=c(1:2,10:12), aus=3:5) # v2

for(ss in 1:length(seasons)){
  dat.season <- dat.planting[dat.planting$month %in% seasons[[ss]],]
  dat.season$doy <- sapply(dat.season$date, manipulateR::doyFromDate)
  rst.season <- raster(rst.rice)
  rst.season[dat.season$cell] <- dat.season$doy
  writeRaster(rst.season, filename = paste(paste0("RiceMap_", ifelse(OTHER_ID!="",paste0(OTHER_ID, "_"), ""), METHOD), COUNTRY, y, names(seasons)[ss], "SIN", sep="_"), format="GTiff", overwrite=TRUE)
  
  rst.seasonwgs <- projectRaster(rst.season, crs = projection(shp.adm), method="ngb")
  rst.seasonwgs <- crop(rst.seasonwgs, manipulateR::cellListExtent(rst.seasonwgs,which(!is.na(rst.seasonwgs[]))))
  
  writeRaster(rst.seasonwgs, filename = paste(paste0("RiceMap_", ifelse(OTHER_ID!="",paste0(OTHER_ID, "_"), ""), METHOD), COUNTRY, y, names(seasons)[ss], "WGS", sep="_"), format="GTiff", overwrite=TRUE)
  rm(rst.season, rst.seasonwgs)
}

message("ORYSAT-", METHOD, ": Writing to disk.")
rst.rice <- writeRaster(rst.rice, paste0(paste("RiceMap", METHOD, TILE, YEAR, ifelse(APPROX_BY_JULIANDAY,"ABJ","ACQ"), sep="_"),".tif"), datatype="INT4S", overwrite=TRUE)
timeen.proc <- timeen.rice <- Sys.time()
timedur.rice <- timeen.proc-timest.extract

message("ORYSAT-", METHOD, ": Done. (", round(timedur.rice,2), " ", attr(timedur.rice, "unit"), ") ")
rm(list=ls())
.rs.restartR()
