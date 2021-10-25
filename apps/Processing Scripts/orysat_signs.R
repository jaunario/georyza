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
WORKSPACE   = "E:/WORKSPACE/Orysat/mosaic/BGD/"
INDICES_DIR = paste0(WORKSPACE, "/indices")
OUTPUT_DIR  = "AUTO"

MISSING_THRES   = 50

METHOD     = "VWIDynamics"
VERSION    = 0.1 
YEAR       = 2007
TILE       = "h26v06"
PRODUCTS    = "MOD09A1"

CHUNKSIZE = 200000 # In-terms of pixels


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

water <- fill.mndwi >= -2280 # Optimum threshold for MNDWI to detact surface water according to Boschetti 2014
nwater.maxlen <- apply(!water, 2, maxlen.condition)

rst.lc[pixels.toprocess[nwater.maxlen<7]] <- -5 # Permanent water
pixels.toprocess <- pixels.toprocess[nwater.maxlen>=7]
pixels.toprocess <- cbind(pixels.toprocess, rowFromCell(rst.lc, pixels.toprocess), colFromCell(rst.lc, pixels.toprocess))

save(pixels.toprocess, rst.lc, date.acqdates, file=paste0(METHOD,"_", TILE, "_250m_", YEAR, ".rdata"))
mat.jday <- mat.jday[,nwater.maxlen>=7]
fill.mndwi <- fill.mndwi[,nwater.maxlen>=7]
smooth.mndwi <- apply(fill.mndwi, 2, safe.sg, n=11)
smooth.mndwi <- round(smooth.mndwi,0)
smooth.mndwi <- as.integer(smooth.mndwi)
smooth.mndwi <- matrix(smooth.mndwi, nrow=length(required.acqdates))

rm(mat.mndwi, fill.mndwi, na.count, water, nwater.maxlen)
gc(reset = TRUE)

message("ORYSAT-", METHOD, ": Extracting EVI.")
stk.evi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="EVI"])
mat.evi <- raster::values(stk.evi)
mat.evi <- mat.evi[pixels.toprocess[,1],]
mat.evi[mat.evi == -32767] <- NA 
mat.evi <- t(mat.evi)
mat.evi <- as.data.frame(mat.evi)

fill.evi <- mapply(approx_int, y=mat.evi, x=mat.jday, xout=list(xout=date.acqdates))

smooth.evi <- apply(fill.evi, 2, safe.sg, n=11)
rm(mat.evi, fill.evi, mat.jday)
gc(reset=TRUE)

if(nrow(pixels.toprocess)>2000000){
  batch_nrows <- nrow(stk.jday)/4 
} else {
  batch_nrows <- nrow(stk.jday)
}

batch.rowst <- seq(1, nrow(stk.evi), by=batch_nrows)
batch.rowen <- pmin(batch.rowst+batch_nrows-1, nrow(stk.evi))

for(i in 1:length(batch.rowen)){
  message("ORYSAT-", METHOD, " BATCH-", i, ": Reading raster files.")
  mat.jday <- raster::values(x=stk.jday, batch.rowst[i], batch_nrows) 
  mat.evi <- raster::values(x=stk.evi, batch.rowst[i], batch_nrows) 
  
  batch.pixels <- pixels.toprocess[pixels.toprocess[,2] >= batch.rowst[i] & pixels.toprocess[,2] <= batch.rowen[i],] 
  if(i>1){
    adj.pix <- batch.pixels[,1]-cellFromRowCol(stk.evi, row = batch.rowen[i-1], col = ncol(stk.evi))
  } else {
    adj.pix <- batch.pixels[,1]
  }
  mat.jday <- mat.jday[adj.pix,]
  mat.evi <- mat.evi[adj.pix,]
  
  chunk.st <- seq(1, length(adj.pix), by=CHUNKSIZE)
  chunk.en <- pmin(chunk.st-1+CHUNKSIZE, length(adj.pix))
  
  for (j in 1:length(chunk.st)){
    chunktime.st <- Sys.time()
    message("ORYSAT-", METHOD, " BATCH-", i, " CHUNK-", j, ": Extracting EVI.")
    chunk.evi <- mat.evi[chunk.st[j]:chunk.en[j],]
    chunk.evi <- t(chunk.evi)
    chunk.evi <- as.data.frame(chunk.evi)
    
    message("ORYSAT-", METHOD, " BATCH-", i, " CHUNK-", j, ": Extracting actual ACQDATES.")
    chunk.jday <- mat.jday[chunk.st[j]:chunk.en[j],]
    chunk.jday <- t(chunk.jday)
    chunk.jday <- as.data.frame(chunk.jday)
    
    message("ORYSAT-", METHOD, " BATCH-", i, " CHUNK-", j, ": Smoothing the time-series.")
    fill.evi <- mapply(approx_int, y=chunk.evi, x=chunk.jday, xout=data.frame(xout=date.acqdates))
    smooth.evi <- apply(fill.evi, 2, safe.sg, n=11)
    colnames(smooth.evi) <- paste0("p", batch.pixels[chunk.st[j]:chunk.en[j],1])
    #saveRDS(smooth.evi, file=paste0("./cache/smooth-evi_batch-", i,"_chunk-",j,".rds"))  
    thischunk  <- apply(smooth.evi, 2, rice.trendanalysis)
    thischunk <- mapply(labelResult, x=thischunk, label=as.list(batch.pixels[chunk.st[j]:chunk.en[j],1]))
    thischunk <- do.call(rbind,thischunk)
    
    if(i==1 & j==1) rice <- thischunk else rice <- rbind(rice, thischunk)
    
    rm(chunk.evi, chunk.jday, fill.evi,smooth.evi)
    gc(reset=TRUE)
    chunktime.en <- Sys.time()
    chunktime.elapsed <- chunktime.en-chunktime.st
    message("ORYSAT-", METHOD, " BATCH-", i, " CHUNK-", j, ": Done. ", format(chunktime.elapsed, unit="min"))
  }
}

# smooth.diff <- apply(smooth.evi, 2, diff)
# idx.maxevi <- apply(smooth.evi, 2, which.max)
# pr1 <- which(idx.maxevi>=8 & idx.maxevi<53)
# 
# 
# # Analyze Maturity Signature... Consistent decline
# post.maxevi <- mapply('[', as.data.frame(smooth.evi[,pr1]), as.data.frame(sapply(idx.maxevi[pr1], seq, length.out=10)))
# post.diff <- apply(post.maxevi, 2, diff)
# post.sign <- apply(post.diff, 2, sign)
# 
# post.consdesc <- colSums(post.sign[1:5,])/-5 
# post.wcmin <- apply(post.diff, 2, which.min)
# 
# # Remove pixels with (idx.maxevi-13) < 0
# psv <- which(idx.maxevi<13)
# pre.maxevi <- mapply('[', as.data.frame(smooth.evi[,pr1]), as.data.frame(sapply(idx.maxevi[pr1]-13, seq, length.out=13)))
# pre.diff <- apply(pre.maxevi, 2, diff)
# pre.sign <- apply(pre.diff, 2, sign)
# 
# 
# 
# xx <- raster(rst.lc)
# xx[pixels.toprocess] <- 0
# xx[pixels.toprocess[pr1]] <- post.consdesc
# xx[pixels.toprocess[post.consdesc]] <- 1
# xx[pixels.toprocess[psv]] <- 2
# 
# #(idx.maxevi+10)>length(date.acqdates))
# pixels.ncv <- pixels.toprocess
# 
# pre.maxevi <- smooth.evi[(idx.maxevi-13):idx.maxevi,]
# fill.evi <- fill.evi[seq(2, length(fill.evi), by=2)]
smooth.evi <- round(smooth.evi,0)
smooth.evi <- as.integer(smooth.evi)
smooth.evi <- matrix(smooth.evi, nrow=length(required.acqdates))
colnames(smooth.evi) <- paste0("p", pixels.toprocess[,1])
rm(stk.evi, mat.evi, fill.evi) #, fill.evi
gc(reset=TRUE)
CHUNK_SIZE = 100000
timest.rice <- Sys.time()
if(length(pixels.toprocess)>2000000){
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
  smooth.evi <- t(smooth.evi)
  smooth.evi <- as.data.frame(smooth.evi)
  
  smooth.mndwi <- t(smooth.mndwi)
  smooth.mndwi <- as.data.frame(smooth.mndwi)
  st <- Sys.time()
  rice <- apply(smooth.evi, 2, rice.VIdynamics)
  en <- Sys.time()
  en-st
                #, evi.date=date.acqdates, eos.evires=0.4)
}
timeen.rice <- Sys.time()
save(rice, file=paste0("trend_", TILE,"_w-15_eprop-p4_a-p05.Rdata"))

smooth.evi <- as.data.frame(smooth.evi)
smooth.mndwi <- as.data.frame(smooth.mndwi)
st <- Sys.time()
rice <- mapply(rice.VWIdynamics, vi = smooth.evi, wi = smooth.mndwi)
en <- Sys.time()
en-st


rice <- mapply(labelResult, x=rice, label=as.list(pixels.toprocess[,1]))

#
rice <- mapply(labelResult, x=rice, label=as.list(pixels.toprocess[match(gt.data$cells.gt,pixels.toprocess)]))

rice <- do.call(rbind,rice)
rice <- as.data.frame(rice)
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

rice <- subset(rice, mean.vi>=2000 & sd.vi>=350 & sd.wi>= 500 & crop.period>=90)
#colnames(rice) <- c("cell", "sos", "eos", "intercept", "b1", "b2", "rsq")
colnames(rice)[1] <- "cell"
rice <- as.data.frame(rice)

rice$sosdate <- date.acqdates[rice$sos]
rice$eosdate <- date.acqdates[rice$eos]
rice$flwrdate <- date.acqdates[rice$flwr]

rice$eosyear <- yearFromDate(rice$eosdate)
rice$sosmonth <- monthFromDate(rice$sosdate)

rice <- rice[rice$eosyear==YEAR, ]
rice$crop.period <- rice$eosdate-rice$sosdate 
quarters <- list(1:3, 4:6, 7:9, 10:12)
for(i in 1:4){
  thisQ <- rice[rice$sosmonth %in% quarters[[i]],]
  thisQ <- thisQ[,c("cell","sosdate", "eosdate")]
  colnames(thisQ) <- c("cell","sos", "eos")
  
  #thisQ <- thisQ[order(thisQ$cell, thisQ$rsq),]
  #thisQ$rsq <- round(thisQ$rsq,4)*10000
  basefname <- paste("RiceMap", METHOD, TILE, YEAR, paste0("Q",i), sep="_")
  for (j in 2:3){
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
