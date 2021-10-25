# Title: IRRI-GIS Implementation of Xiao (2006) Rice Mapping Algorithm
#     -1 - Forest/Thick Vegetation 
#     -2 - Built-up
#     -3 - Ocean and High-elevation
#     -4 - Cloud
#     -5 - Snow
#     -6 - Persistent water
# library(future.apply)

# DIRECTORY SETTINGS
MODIS_HOME = "F:/MODIS/v6"
WORKSPACE  = "E:/WORKSPACE/Orysat"
INDICES_DIR = paste0(WORKSPACE, "/indices")

APPROX_BY_JULIANDAY = TRUE
CONSMISSING_THRES   = 5
MISSING_THRES   = 50


METHOD     = "ModXiao"
YEAR       = 2017
TILE       = "h25v06"
PRODUCTS    = "MOD09A1"

library(orysat)
library(manipulateR)

maxlen.condition <- function(x){
  if(sum(is.na(x))==length(x)) {
    result <- NA
  } else {
    x <- which(x)
    result <- maxConsecutive(x)
  }
  return(result)
}

setwd(WORKSPACE)

pylibs.start(python="~/../miniconda3", required=TRUE)
mdata.lc <- modis.readHDF(dir(MODIS_HOME, pattern=paste(paste0("MCD12Q1.A",YEAR,"001.",TILE),"hdf",sep=".*."), recursive = TRUE, full.names = TRUE))

rst.lc <- raster(mdata.lc@extent, ncol=mdata.lc@ncols, nrow=mdata.lc@nrows, crs=mdata.lc@projection)
rst.lc[mdata.lc@imgvals$LC_Type1<=5 | mdata.lc@imgvals$LC_Type1==8] <- -1 # FORESTS AND OTHER VEGETATION
rst.lc[mdata.lc@imgvals$LC_Type1==13] <- -2 # Built-up

ACQDOYS <- seq(from=1,to=361, by=8)
required.acqdates <- paste("A",c(paste(YEAR-1,sprintf("%03g", ACQDOYS)[(length(ACQDOYS)-7):length(ACQDOYS)], sep=""),  #PREVIOUS YEAR
                                 paste(YEAR,sprintf("%03g", ACQDOYS), sep=""),           # CURRENT YEAR
                                 paste(YEAR+1,sprintf("%03g",ACQDOYS)[1:11], sep="")),   # SUCCEEDING YEAR
                           sep="")
date.acqdates <- as.Date(required.acqdates, "A%Y%j")

# Extract data
message("ORYSAT-", METHOD, ": Listing files.")
filelist.index <- dir(INDICES_DIR, pattern = paste0(PRODUCTS, ".", TILE, ".*.tif$"), full.names = TRUE)
inv.idxfiles <-  inventory.modis(filelist.index, modisinfo = c("product", "zone", "acqdate", "band"), file.ext = "tif")
inv.idxfiles <- inv.idxfiles[inv.idxfiles$acqdate %in% required.acqdates,]
inv.idxfiles <- inv.idxfiles[order(inv.idxfiles$band, inv.idxfiles$acqdate),]

timest.extract <- Sys.time()

message("ORYSAT-", METHOD, ": Extracting MNDWI.")
stk.mndwi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="MNDWI"])
mat.mndwi <- raster::values(stk.mndwi)

rst.lc[rowSums(is.na(mat.mndwi))==65] <- -3 # Ocean and Mountains 
pixels.toprocess <- which(is.na(rst.lc[]))

mat.mndwi <- mat.mndwi[pixels.toprocess,]

message("ORYSAT-", METHOD, ": Extracting pixel acqdate date.")
stk.jday <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="actualdate"])
mat.jday <- raster::values(stk.jday)
mat.jday <- mat.jday[pixels.toprocess,]
mat.jday <- t(mat.jday)
mat.jday <- as.data.frame(mat.jday)

mat.mndwi <- t(mat.mndwi)
mat.mndwi <- as.data.frame(mat.mndwi)

message("ORYSAT-", METHOD, ": Adjusting mndwi to equally-spaced acqdate date.")
fill.mndwi <- mapply(approx_int, y=mat.mndwi, x=mat.jday, xout=list(xout=date.acqdates))
# smooth.mndwi <- apply(fill.mndwi, 2, safe.sg, n=11) Smoothing water index seem problematic

water <- fill.mndwi >= -2280

nwater.maxlen <- apply(!water, 2, maxlen.condition)
#water.maxlen <- apply(water, 1, maxlen.condition)

rst.lc[pixels.toprocess[nwater.maxlen<10]] <- -6 # Permanent water
#rst.lc[water.maxlen>=45] <- -5 # Permanent water

pixels.toprocess <- pixels.toprocess[nwater.maxlen>=10]
water <- water[,nwater.maxlen>=10]
water <- as.data.frame(water)
mat.jday <- mat.jday[,nwater.maxlen>=10]

rm(mat.mndwi, fill.mndwi, smooth.mndwi, nwater.maxlen)
gc(reset=TRUE)

message("ORYSAT-", METHOD, ": Extracting EVI.")
stk.evi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="EVI"])
mat.evi <- raster::values(stk.evi)
mat.evi <- mat.evi[pixels.toprocess,]

# TODO: If non-snow area, snow pixels are NA, else period of snow is non-rice 
mat.evi[mat.evi == -32767] <- NA 

# For masking out pixels with too many consecutive NAs (cloud or missing)
message("ORYSAT-", METHOD, ": Identifying pixels with sufficient data.")
#na.maxlen <- apply(is.na(mat.evi), 1, maxlen.condition)
na.count <- rowSums(is.na(mat.evi)) 

rst.lc[pixels.toprocess[which(na.count>(ncol(mat.evi)*0.8))]] <- -7 # Cloud

# rst.lc[na.maxlen>CONSMISSING_THRES & na.maxlen<65] <- -4 # Cloud

pixels.toprocess <- pixels.toprocess[which(!(na.count>(ncol(mat.evi)*0.8)))]
mat.evi  <- mat.evi[which(!(na.count>(ncol(mat.evi)*0.8))),]
mat.jday <- mat.jday[,which(!(na.count>(ncol(mat.evi)*0.8)))]
water <- water[,which(!(na.count>(ncol(mat.evi)*0.8)))]

#water <- as.data.frame(water)

mat.evi <- t(mat.evi)
mat.evi <- as.data.frame(mat.evi)
#mat.evi <- mat.evi[,nwater.maxlen>=10]

fill.evi <- mapply(approx_int, y=mat.evi, x=mat.jday, xout=data.frame(xout=date.acqdates))
#smooth.evi <- apply(fill.evi, 2, safe.sg, n=11)



rice <- mapply(FUN=rice.modxiao, ts.vi=as.data.frame(fill.evi), ts.flood=water)



# message("ORYSAT-", METHOD, ": Identifying snow.")
# For masking out pixels with too many consecutive Snow or too many snow instances
# in case snw.maxlen is not good enough, the 2 conditions below can be explored
# It is possible that there's only snow at one time but rice can be planted some other time within the year
#
# snw.count <- rowSums(mat.evi == -32767)
# nsnw.maxlen <- apply(mat.evi != -32767, 1, maxlen.condition)
# snw.maxlen <- apply(mat.evi == -32767, 1, maxlen.condition)
# rst.lc[snw.maxlen>3] <- -5 # Snow


# LOOKING AT GT POINTS ONLY (REMOVE LATER)
#shp.gt <- shapefile("BGD_GT_2010.shp")

# gt.evi <- mat.evi[which(pixels.toprocess %in% shp.gt$h26v06_), ]
# gt.jday <- mat.jday[which(pixels.toprocess %in% shp.gt$h26v06_), ]
# gt.mndwi <- mat.mndwi[which(pixels.toprocess %in% shp.gt$h26v06_), ]
# 
# mat.evi <- t(gt.evi)
# mat.mndwi <- t(gt.mndwi)
# mat.jday <- t(gt.jday)

# END OF GT TRIAL

#rice <- mapply(polregIT.rice, pix.evi = smooth.evi, evi.date = list(evi.date=date.acqdates), dates.covered = list(date=date.acqdates), alpha=0.5)

# rm(mat.evi)
# fill.evi <- fill.evi/10000
# fill.evi <- round(fill.evi,4)
# fill.evi <- as.data.frame(fill.evi)
# smooth.evi <- as.data.frame(smooth.evi)
# 
# rm(mat.mndwi)
# 
# 
# water <- fill.mndwi >= -0.2280 # Optimum threshold for MNDWI to detact surface water according to Boschetti 2014
# 
# water <- water[, nwater.maxlen>=10]
# water <- as.data.frame(water)


#xx <- apply(mat.ndvi[grepl(YEAR, required.acqdates),]-7000,2, ttest.p)
#rst.lc2[pixels.toprocess[xx>0.05]] <- -4 # Forests


rm(mdata.lc, stk.mndwi, fill.mndwi, na.maxlen, snw.maxlen, nwater.maxlen)
# 
# message("ORYSAT-", METHOD, ": Identifying shrubs.")
# shrb.count <- colSums(mat.lswi[grepl(YEAR, required.acqdates),]<1000, na.rm = TRUE)
# rst.lc[pixels.toprocess[shrb.count==0 & is.na(rst.lc[pixels.toprocess])]] <- -5
# 
# message("ORYSAT-", METHOD, ": Identifying permanent water.")
# nh2o <- mat.ndvi>=1000 | mat.ndvi >= mat.lswi
# nh2o.maxlen <- apply(nh2o, 2, maxlen.condition)
# rst.lc[pixels.toprocess[nh2o.maxlen<15 & is.na(rst.lc[pixels.toprocess])]] <- -6
# rm(nh2o, nh2o.maxlen, shrb.count, forest.count, na.maxlen)
# gc(reset=TRUE)


rst.rice <- raster(rst.lc)
rst.rice[pixels.toprocess] <- rice # mapply(FUN=rice.modxiao, ts.vi=fill.evi, ts.flood=water) #, evi.ricemax=0, evi.halfricemax=0) 


message("ORYSAT-", METHOD, ": Writing to disk.")
rst.lc <- writeRaster(rst.lc, paste0(paste("LANDCOVER", METHOD, TILE, YEAR, sep="_"),".tif"), datatype="INT4S", overwrite=TRUE)
rst.rice <- writeRaster(rst.rice, paste0(paste("RiceMap", METHOD, TILE, YEAR, ifelse(APPROX_BY_JULIANDAY,"ABJ","ACQ"), sep="_"),".tif"), datatype="INT4S", overwrite=TRUE)
timeen.proc <- timeen.rice <- Sys.time()
timedur.rice <- timeen.proc-timest.extract

rm(pixels.toprocess, mat.jday, rst.lc, rst.rice, stk.evi, stk.jday, water, fill.evi)

message("ORYSAT-", METHOD, ": Done. (", round(timedur.rice,2), " ", attr(timedur.rice, "unit"), ") ")
