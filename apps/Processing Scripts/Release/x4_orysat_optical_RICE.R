# Title: IRRI-GIS Implementation of Xiao (2006) Rice Mapping Algorithm
#     -1 - Forest/Thick Vegetation 
#     -2 - Built-up
#     -3 - Ocean and High-elevation
#     -4 - Cloud
#     -5 - Snow
#     -6 - Persistent water
# library(future.apply)

# DIRECTORY SETTINGS
MODIS_HOME  = "F:/MODIS"
WORKSPACE   = "E:/WORKSPACE/FLAGSHIP"
INDICES_DIR = paste0(WORKSPACE, "/Filled")
CACHE_DIR    = paste0(WORKSPACE, "/cache")
OUTPUT_DIR  = "AUTO"


MISSING_THRES   = 50

METHOD     = "VWIDynamics"
VERSION    = 0.1 
YEAR       = 2007
TILE       = "h25v06"
PRODUCTS    = "MOD09A1"

# TODO: Have default layer-parameter mapping but allow custom mapping for experimental methods
LAYERPARAM_MAPPING  = list(vi="EVI", wi="MNDWI")

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

#pylibs.start(python="~/../miniconda3", required=TRUE)
#mdata.lc <- modis.readHDF(dir(MODIS_HOME, pattern=paste(paste0("MCD12Q1.A",YEAR,"001.",TILE),"hdf",sep=".*."), recursive = TRUE, full.names = TRUE))

# rst.lc <- raster(mdata.lc@extent, ncol=mdata.lc@ncols, nrow=mdata.lc@nrows, crs=mdata.lc@projection)
# rst.lc[mdata.lc@imgvals$LC_Type1<=5 | mdata.lc@imgvals$LC_Type1==8] <- -1 # FORESTS AND OTHER VEGETATION
# rst.lc[mdata.lc@imgvals$LC_Type1==13] <- -2 # Built-up

ACQDOYS <- seq(from=1,to=361, by=8)
required.acqdates <- paste("A",c(paste(YEAR-1,sprintf("%03g", ACQDOYS)[(length(ACQDOYS)-7):length(ACQDOYS)], sep=""),  #PREVIOUS YEAR
                                 paste(YEAR,sprintf("%03g", ACQDOYS), sep=""),           # CURRENT YEAR
                                 paste(YEAR+1,sprintf("%03g",ACQDOYS)[1:4], sep="")),   # SUCCEEDING YEAR
                           sep="")
date.acqdates <- as.Date(required.acqdates, "A%Y%j")

# Extract data
message("ORYSAT-", METHOD, ": Listing files.")
filelist.index <- dir(paste0(INDICES_DIR,"/", TILE), pattern = paste0(PRODUCTS, ".", TILE, ".*.tif$"), full.names = TRUE, recursive = TRUE)
inv.idxfiles <-  inventory.modis(filelist.index, modisinfo = c("product", "zone", "acqdate", "band"), file.ext = "tif")
inv.idxfiles <- inv.idxfiles[inv.idxfiles$acqdate %in% required.acqdates,]
inv.idxfiles <- inv.idxfiles[order(inv.idxfiles$band, inv.idxfiles$acqdate),]

timest.extract <- Sys.time()

# message("ORYSAT - Mask and Fill: Creating mask using DEM and SLOPE." )
# rst.dem <- raster(dir(paste0(MODIS_HOME, "/dem"), pattern=paste0(TILE,".*.DEM.tif$"), full.names = TRUE, recursive = TRUE))
# rst.slp <- raster(dir(paste0(MODIS_HOME, "/dem"), pattern=paste0(TILE,".*.SLOPE.tif$"), full.names = TRUE, recursive=TRUE))
# 
# rst.demmask <- (rst.dem>2000 | rst.slp>2)
# idx.tomask <- which(rst.demmask[]==1)
# rm(rst.dem, rst.slp, rst.demmask)
# gc(reset=TRUE)


pixels.toprocess <- readRDS(file=paste0(CACHE_DIR, "/NUM_INDEX_", paste(TILE, YEAR, "IDX_TOPROC", "rds", sep=".")))
# pixels.toprocess <- pixels.toprocess[which(!pixels.toprocess %in% idx.tomask )]

message("ORYSAT-", METHOD, ": Extracting EVI.")
stk.evi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="EVI"])
mat.evi <- raster::values(stk.evi)
mat.evi <- mat.evi[pixels.toprocess,]
mat.evi <- t(mat.evi)
mat.evi <- as.data.frame(mat.evi)
avg.evi <- rowMeans(mat.evi, na.rm = TRUE)
sd.evi  <- apply(mat.evi, 1, sd)

message("ORYSAT-", METHOD, ": Extracting MNDWI.")
stk.mndwi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="MNDWI"])
mat.mndwi <- raster::values(stk.mndwi)
mat.mndwi <- mat.mndwi[pixels.toprocess, ]
mat.mndwi <- t(mat.mndwi)
mat.mndwi <- as.data.frame(mat.mndwi)

rice <- mapply(rice.VWIdynamics, vi=mat.evi, wi=mat.mndwi)
rm(mat.evi, mat.mndwi)
gc(reset=T)
# CHUNK_SIZE = 100000
# timest.rice <- Sys.time()
# if(length(pixels.toprocess)>1000000){
#   chunk <- ceiling(length(pixels.toprocess)/CHUNK_SIZE)
#   chunks.st <- seq(1,length(pixels.toprocess), by = CHUNK_SIZE)
#   chunks.en <- c(chunks.st[-1]-1, length(pixels.toprocess))
# 
#   
#   for(i in 1:length(chunks.st)){
#     message("chunk ", i, " of ",length(chunks.st))
#     #rst.rice[pixels.toprocess[i]] <- rice.ts(ts.vi=fill.evi[,i],ts.date = mat.jday[,i], ts.flood = water[,i])
#     thischunk  <- mapply(rice.VWIdynamics, vi=mat.evi[,chunks.st[i]:chunks.en[i]], wi=mat.mndwi[,chunks.st[i]:chunks.en[i]])
#     if(i==1) rice <- thischunk else rice <- c(rice, thischunk)
#     rm(thischunk)
#   }
# 
# } else{
#   rice <- mapply(rice.VWIdynamics, vi=mat.evi, wi=mat.mndwi)
# }
timeen.rice <- Sys.time()

rice <- mapply(labelResult, x=rice, label=as.list(pixels.toprocess))
rice <- do.call(rbind,rice)
rice <- as.data.frame(rice)
colnames(rice)[1] <- "cell"
saveRDS(rice, file=paste0(CACHE_DIR, "/DF_RICE-", METHOD, "_", paste(TILE, YEAR, "VWDProds", "rds", sep=".")))

#colnames(rice) <- c("cell", "sos", "eos", "intercept", "b1", "b2", "rsq")

rice$sosdate <- date.acqdates[rice$sos]
rice$eosdate <- date.acqdates[rice$eos]
rice$flwrdate <- date.acqdates[rice$flwr]

rice$eosyear <- yearFromDate(rice$eosdate)
rice$sosmonth <- monthFromDate(rice$sosdate)

rice <- rice[rice$eosyear==YEAR, ]
rice$crop.period <- rice$eosdate-rice$sosdate 
rice <- subset(rice, mean.vi>=2000 & sd.vi>=350 & sd.wi>= 500 & crop.period>=90)

quarters <- list(1:3, 4:6, 7:9, 10:12)
for(i in 1:4){
  thisQ <- rice[rice$sosmonth %in% quarters[[i]],]
  thisQ <- thisQ[,c("cell","sosdate", "eosdate")]
  colnames(thisQ) <- c("cell","sos", "eos")
  
  #thisQ <- thisQ[order(thisQ$cell, thisQ$rsq),]
  #thisQ$rsq <- round(thisQ$rsq,4)*10000
  basefname <- paste("RiceMap", METHOD, TILE, YEAR, paste0("Q",i), sep="_")
  for (j in 2:3){
    rst.riceresults <- raster(stk.evi)
    rst.riceresults[thisQ$cell] <- as.numeric(format(as.Date(thisQ[,j]),"%j"))
    rst.riceresults <- writeRaster(rst.riceresults, filename = paste0(OUTPUT_DIR, "/", basefname, "_", toupper(colnames(thisQ)[j]), ".tif"), datatype="INT4S", overwrite=TRUE)
  }
}

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
# 
# colnames(rice) <- c("cell", "sos", "eos", "intercept", "b1", "b2", "rsq")
# rice <- as.data.frame(rice)
# 
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
#   for (j in 2:7){
#     rst.riceresults <- raster(rst.lc)
#     rst.riceresults[thisQ$cell] <- thisQ[,j]
#     rst.riceresults <- writeRaster(rst.riceresults, filename = paste0(OUTPUT_DIR, "/", basefname, "_", toupper(colnames(thisQ)[j]), ".tif"), datatype="INT4S", overwrite=TRUE)
#   }
# }

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




# message("ORYSAT-", METHOD, ": Writing to disk.")
# rst.lc <- writeRaster(rst.lc, paste0(OUTPUT_DIR, "/", paste("LANDCOVER", METHOD, TILE, YEAR, sep="_"),".tif"), datatype="INT4S", overwrite=TRUE)
# timeen.proc <- timeen.rice <- Sys.time()
# timedur.rice <- timeen.proc-timest.extract


message("ORYSAT-", METHOD, ": Done. (", round(timedur.rice,2), " ", attr(timedur.rice, "unit"), ") ")
rm(list=ls())

.rs.restartR()
