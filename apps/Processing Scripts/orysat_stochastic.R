# Title: IRRI-GIS Implementation of Xiao (2006) Rice Mapping Algorithm
#     -1 - Forest/Thick Vegetation 
#     -2 - Built-up
#     -3 - Ocean and High-elevation
#     -4 - Cloud
#     -5 - Snow
#     -6 - Persistent water
# library(future.apply)

# DIRECTORY SETTINGS
MODIS_HOME  = "F:/MODIS/v6"
WORKSPACE   = "E:/WORKSPACE/Orysat"
INDICES_DIR = paste0(WORKSPACE, "/indices")
OUTPUT_DIR  = "AUTO"


MISSING_THRES   = 50

METHOD     = "IterRegr"
VERSION    = 0.1 
YEAR       = 2017
TILE       = "h26v06"
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


approx_int <- function(...){
  result <- approx(...)$y
  result <- as.integer(round(result))
  return(result)
}

setwd(WORKSPACE)

if(OUTPUT_DIR=="AUTO") {
  OUTPUT_DIR <- paste0("./",paste(METHOD, VERSION, PRODUCTS, sep="_"))
  if(!dir.exists(OUTPUT_DIR)) force.directories(OUTPUT_DIR)
}

pylibs.start(python="~/../miniconda3", required=TRUE)
mdata.lc <- modis.readHDF(dir(MODIS_HOME, pattern=paste(paste0("MCD12Q1.A",YEAR,"001.",TILE),"hdf",sep=".*."), recursive = TRUE, full.names = TRUE))

rst.lc <- raster(mdata.lc@extent, ncol=mdata.lc@ncols, nrow=mdata.lc@nrows, crs=mdata.lc@projection)
rst.lc[mdata.lc@imgvals$LC_Type1<=5 | mdata.lc@imgvals$LC_Type1==8] <- -1 # FORESTS AND OTHER VEGETATION
rst.lc[mdata.lc@imgvals$LC_Type1==13] <- -2 # Built-up

ACQDOYS <- seq(from=1,to=361, by=8)
required.acqdates <- paste("A",c(paste(YEAR-1,sprintf("%03g", ACQDOYS)[(length(ACQDOYS)-10):length(ACQDOYS)], sep=""),  #PREVIOUS YEAR
                                 paste(YEAR,sprintf("%03g", ACQDOYS), sep=""),           # CURRENT YEAR
                                 paste(YEAR+1,sprintf("%03g",ACQDOYS)[1:12], sep="")),   # SUCCEEDING YEAR
                           sep="")
date.acqdates <- as.Date(required.acqdates, "A%Y%j")

# Extract data
message("ORYSAT-", METHOD, ": Listing files.")
filelist.index <- dir(INDICES_DIR, pattern = paste0(PRODUCTS, ".", TILE, ".*.tif$"), full.names = TRUE)
inv.idxfiles <-  inventory.modis(filelist.index, modisinfo = c("product", "zone", "acqdate", "band"), file.ext = "tif")
inv.idxfiles <- inv.idxfiles[inv.idxfiles$acqdate %in% required.acqdates,]
inv.idxfiles <- inv.idxfiles[order(inv.idxfiles$band, inv.idxfiles$acqdate),]

timest.extract <- Sys.time()

# Water Dynamics Analysis
message("ORYSAT-", METHOD, ": Extracting MNDWI.")
stk.mndwi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="MNDWI"])
mat.mndwi <- raster::values(stk.mndwi)

# TODO: If non-snow area, snow pixels are NA, else period of snow is non-rice 
mat.mndwi[mat.mndwi == -32767] <- NA 
na.count <- rowSums(is.na(mat.mndwi)) 
rst.lc[na.count == 69] <- -3 # Ocean and Mountains 
rst.lc[na.count > MISSING_THRES & na.count < 69] <- -4 # Insufficent Data

pixels.toprocess <- which(is.na(values(rst.lc)))
mat.mndwi <- mat.mndwi[pixels.toprocess, ]
mat.mndwi <- t(mat.mndwi)
mat.mndwi <- as.data.frame(mat.mndwi)

message("ORYSAT-", METHOD, ": Extracting pixel acqdate date.")
stk.jday <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="actualdate"])
mat.jday <- raster::values(stk.jday) 
mat.jday <- mat.jday[pixels.toprocess,]
mat.jday <- t(mat.jday)
mat.jday <- as.data.frame(mat.jday)

fill.mndwi <- mapply(approx_int, y=mat.mndwi, x=mat.jday, xout=list(xout=date.acqdates))
#smooth.mndwi <- apply(fill.mndwi, 2, safe.sg, n=11)

water <- fill.mndwi >= -2280 # Optimum threshold for MNDWI to detact surface water according to Boschetti 2014
nwater.maxlen <- apply(!water, 2, maxlen.condition)

rst.lc[pixels.toprocess[nwater.maxlen<7]] <- -5 # Permanent water
pixels.toprocess <- pixels.toprocess[nwater.maxlen>=7]
mat.jday <- mat.jday[,nwater.maxlen>=7]

rm(stk.mndwi, mat.mndwi, fill.mndwi, smooth.mndwi, na.count, nwater.maxlen)
gc(reset = TRUE)
message("ORYSAT-", METHOD, ": Extracting EVI.")
stk.evi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="EVI"])
mat.evi <- raster::values(stk.evi)
mat.evi <- mat.evi[pixels.toprocess,]
mat.evi[mat.evi == -32767] <- NA 
mat.evi <- t(mat.evi)
mat.evi <- as.data.frame(mat.evi)

# Check low variability
# means.pixel <- colMeans(mat.evi, na.rm = TRUE)
# sd.pixel <- mapply(sd, mat.evi, na.rm=TRUE)
# min.pixel <- mapply(min, mat.evi, na.rm=TRUE)
# max.pixel <- mapply(max, mat.evi, na.rm=TRUE)
# 
# rst.sd <- raster(rst.lc)
# rst.sd[pixels.toprocess] <- means.pixel>3000 & sd.pixel<800 #(max.pixel - min.pixel)<2000
# plot(rst.sd)
# 
fill.evi <- mapply(approx_int, y=mat.evi, x=mat.jday, xout=data.frame(xout=date.acqdates))
smooth.evi <- apply(fill.evi, 2, safe.sg, n=11)

# fill.evi <- fill.evi[seq(2, length(fill.evi), by=2)]
smooth.evi <- round(unlist(smooth.evi))
smooth.evi <- as.integer(smooth.evi)
smooth.evi <- matrix(smooth.evi, nrow=length(required.acqdates))
colnames(smooth.evi) <- paste0("p", pixels.toprocess)
rm(stk.evi, mat.evi, fill.evi) #, fill.evi
gc(reset=TRUE)
CHUNK_SIZE = 100000
timest.rice <- Sys.time()
if(length(pixels.toprocess)>1000000){
  chunk <- ceiling(length(pixels.toprocess)/CHUNK_SIZE)
  chunks.st <- seq(1,length(pixels.toprocess), by = CHUNK_SIZE)
  chunks.en <- c(chunks.st[-1]-1, length(pixels.toprocess))

  
  for(i in 1:length(chunks.st)){
    message("chunk ", i, " of ",length(chunks.st))
    #rst.rice[pixels.toprocess[i]] <- rice.ts(ts.vi=fill.evi[,i],ts.date = mat.jday[,i], ts.flood = water[,i])
    thischunk  <- apply(smooth.evi[,chunks.st[i]:chunks.en[i]], 2, rice.itreg, evi.date=date.acqdates, eos.evires=0.4)
    if(i==1) rice <- thischunk else rice <- c(rice, thischunk)
    rm(thischunk)
  }

} else{
  rice <- apply(smooth.evi, 2, rice.itreg, evi.date=date.acqdates, eos.evires=0.4)
}
timeen.rice <- Sys.time()
save(rice, file=paste0("itreg_", TILE,"_w-15_eprop-p4_a-p05.Rdata"))

rice <- mapply(labelResult, x=rice, label=as.list(pixels.toprocess))
rice <- do.call(rbind,rice)
# rst.intensity <- raster(rst.lc)
# rst.intensity[pixels.toprocess] <- sapply(rice, length)/4

# toverify <- which(rst.intensity[]>3)
# toverify.adj <- match(toverify,pixels.toprocess) 
# FOR DEBUGGING
# rice <- vector()
# for(i in 1:length(pixels.toprocess)){
#   thispix <- rice.itreg(pix.evi = smooth.evi[,i], evi.date=date.acqdates, eos.evires=0.4)
#   message(pixels.toprocess[i], ":", length(thispix)/4, " seasons of rice (", i, " of ", length(pixels.toprocess),").")
#   if (length(thispix)>0)  rice <- rbind(rice, thispix)
# }

colnames(rice) <- c("cell", "sos", "eos", "intercept", "b1", "b2", "rsq")
rice <- as.data.frame(rice)

rice$sosdate <- as.Date(rice$sos, "1970-1-1")
rice$eosdate <- as.Date(rice$eos, "1970-1-1")
rice$eosyear <- yearFromDate(rice$eosdate)
rice$sosmonth <- monthFromDate(rice$sosdate)

rice <- rice[rice$eosyear==YEAR, ]

quarters <- list(1:3, 4:6, 7:9, 10:12)
for(i in 1:4){
  thisQ <- rice[rice$sosmonth %in% quarters[[i]],]
  thisQ <- thisQ[order(thisQ$cell, thisQ$rsq),]
  thisQ$rsq <- round(thisQ$rsq,4)*10000
  basefname <- paste("RiceMap", METHOD, TILE, YEAR, paste0("Q",i), sep="_")
  for (j in 2:7){
    rst.riceresults <- raster(rst.lc)
    rst.riceresults[thisQ$cell] <- thisQ[,j]
    rst.riceresults <- writeRaster(rst.riceresults, filename = paste0(OUTPUT_DIR, "/", basefname, "_", toupper(colnames(thisQ)[j]), ".tif"), datatype="INT4S", overwrite=TRUE)
  }
}

# colnames(rice) <- c("cell", "sos", "eos", "datapoints", "rsq")
# rice <- as.data.frame(rice)
# rice$sosdate <- as.Date(rice$sos, "1970-1-1")
# rice$eosdate <- as.Date(rice$eos, "1970-1-1")
# rice$eosyear <- yearFromDate(rice$eosdate)
# rice$sosmonth <- monthFromDate(rice$sosdate)
# 
# rice <- rice[rice$eosyear==YEAR, ]
# 
# quarters <- list(1:3, 4:6, 7:9, 10:12)
# for(i in 1:4){
#   thisQ <- rice[rice$sosmonth %in% quarters[[i]],]
#   thisQ <- thisQ[order(thisQ$cell, thisQ$rsq),]
#   thisQ$rsq <- round(thisQ$rsq,4)*10000
#   basefname <- paste("RiceMap", METHOD, TILE, YEAR, paste0("Q",i), sep="_")
#   for (j in c(2,3,5)){
#     rst.riceresults <- raster(rst.lc)
#     rst.riceresults[thisQ$cell] <- thisQ[,j]
#     rst.riceresults <- writeRaster(rst.riceresults, filename = paste0(OUTPUT_DIR, "/", basefname, "_", toupper(colnames(thisQ)[j]), ".tif"), datatype="INT4S", overwrite=TRUE)
#   }
# }

# TODO: M ove this part on actual analysis rather than nasking them outright as some pixel could have enough data depending on the period
# On Final Rice Map indicate pixels with insufficent data
# For masking out pixels with too many consecutive NAs (cloud or missing)
# message("ORYSAT-", METHOD, ": Identifying pixels with sufficient data.")
# na.maxlen <- apply(is.na(mat.mndwi), 1, maxlen.condition)
# rst.lc[na.maxlen>CONSMISSING_THRES & na.maxlen<65] <- -4 # Cloud

# TODO: Snow masking should be confined to regions where there is actual snow (Arctic regions)
# message("ORYSAT-", METHOD, ": Identifying snow.")
# For masking out pixels with too many consecutive Snow or too many snow instances
# in case snw.maxlen is not good enough, the 2 conditions below can be explored
# It is possible that there's only snow at one time but rice can be planted some other time within the year
#
# snw.count <- rowSums(mat.ndvi == -32767) 
# nsnw.maxlen <- apply(mat.mndwi != -32767, 1, maxlen.condition)
# snw.maxlen <- apply(mat.mndwi == -32767, 1, maxlen.condition)

# Initial snow masking
# rst.lc[snw.maxlen>3] <- -5 # Snow


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




message("ORYSAT-", METHOD, ": Writing to disk.")
rst.lc <- writeRaster(rst.lc, paste0(OUTPUT_DIR, "/", paste("LANDCOVER", METHOD, TILE, YEAR, sep="_"),".tif"), datatype="INT4S", overwrite=TRUE)
timeen.proc <- timeen.rice <- Sys.time()
timedur.rice <- timeen.proc-timest.extract


message("ORYSAT-", METHOD, ": Done. (", round(timedur.rice,2), " ", attr(timedur.rice, "unit"), ") ")
rm(list=ls())

.rs.restartR()
