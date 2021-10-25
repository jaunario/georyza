# Title: IRRI-GIS Implementation of Xiao (2006) Rice Mapping Algorithm
#     -1 - Forest/Thick Vegetation 
#     -2 - Built-up
#     -3 - Ocean and High-elevation
#     -4 - Cloud
#     -5 - Snow
#     -6 - Persistent water
# library(future.apply)

# DIRECTORY SETTINGS
MODIS_HOME  = "D:/MODIS/v6"
WORKSPACE   = "E:/WORKSPACE/Orysat/250m"
SMOOTH_DIR = paste0(WORKSPACE, "/cache")
OUTPUT_DIR  = "AUTO"


MISSING_THRES  = 50

METHOD         = "IterRegr"
VERSION        = 0.1
YEAR           = 2010
TILE           = "h25v06"
PRODUCTS       = "MD13Q1"
SKIPBY         = 5
ASSIGNED_GROUP = 5

library(orysat)
library(manipulateR)

setwd(WORKSPACE)

if(OUTPUT_DIR=="AUTO") {
  OUTPUT_DIR <- paste0("./",paste(METHOD, VERSION, PRODUCTS, sep="_"))
  if(!dir.exists(OUTPUT_DIR)) force.directories(OUTPUT_DIR)
}


ACQDOYS <- seq(from=1,to=361, by=8)
required.acqdates <- paste("A",c(paste(YEAR-1,sprintf("%03g", ACQDOYS)[(length(ACQDOYS)-10):length(ACQDOYS)], sep=""),  #PREVIOUS YEAR
                                 paste(YEAR,sprintf("%03g", ACQDOYS), sep=""),           # CURRENT YEAR
                                 paste(YEAR+1,sprintf("%03g",ACQDOYS)[1:12], sep="")),   # SUCCEEDING YEAR
                           sep="")
date.acqdates <- as.Date(required.acqdates, "A%Y%j")

# Extract data
timest.proc <- Sys.time()

message("ORYSAT-", METHOD, ": Listing files.")
filelist.smoothevi <- dir(SMOOTH_DIR, pattern = paste0("smooth-evi_.*.rds$"), full.names = TRUE)
chunks_toproces <- seq(ASSIGNED_GROUP, ceiling(length(filelist.smoothevi)/SKIPBY)*SKIPBY, by=SKIPBY)
chunks_toproces <- pmin(chunks_toproces,length(filelist.smoothevi))

for (i in chunks_toproces){
  if(file.exists(paste0(paste("RICE", METHOD, TILE, ASSIGNED_GROUP, i, sep="_"), ".rds"))) next
  rice <- vector()
  saveRDS(rice, file=paste0(paste("RICE", METHOD, TILE, ASSIGNED_GROUP, i, sep="_"), ".rds"))
  chunktime.st <- Sys.time()
  message("ORYSAT-", METHOD, " CHUNK-", i, ": Analyzing temporal signature.")
  smooth.evi <- readRDS(filelist.smoothevi[i])
  thischunk  <- apply(smooth.evi, 2, rice.itreg, evi.date=date.acqdates, eos.evires=0.4)
  thischunk <- mapply(labelResult, x=thischunk, label=as.list(as.integer(sub("p", "", names(thischunk)))))
  
  rice <- do.call(rbind,thischunk) 
  saveRDS(rice, file=paste0(paste("RICE", METHOD, TILE, ASSIGNED_GROUP, i, sep="_"), ".rds"))
  
  message("ORYSAT-", METHOD, " CHUNK-", i, ": Cleaning-up.")
  rm(rice, smooth.evi, thischunk) #, fill.evi
  gc(reset=TRUE)
  chunktime.en <- Sys.time()
  chunktime.elapsed <- chunktime.en-chunktime.st
  message("ORYSAT-", METHOD, " CHUNK-", i, ": Done. ", format(chunktime.elapsed, unit="min"))
}

load("itregr_h25v06_250m_2010.rdata")
#
# COMBINE ALL RICE OUTPUTS FROM OTHER INSTANCES
files.rice <- dir(pattern = paste0(paste("RICE", METHOD, TILE, sep="_"), ".*.rds$"))
for(i in 1:length(files.rice)){
  rice <- readRDS(files.rice[i])
  if(i==1) all.rice <- rice else all.rice <- rbind(all.rice, rice)
}
rice <- unique(all.rice)
rm(all.rice)
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

# # TODO: Move this part on actual analysis rather than asking them outright as some pixel could have enough data depending on the period
# # On Final Rice Map indicate pixels with insufficent data
# # For masking out pixels with too many consecutive NAs (cloud or missing)
# # message("ORYSAT-", METHOD, ": Identifying pixels with sufficient data.")
# # na.maxlen <- apply(is.na(mat.mndwi), 1, maxlen.condition)
# # rst.lc[na.maxlen>CONSMISSING_THRES & na.maxlen<65] <- -4 # Cloud
# 
# # TODO: Snow masking should be confined to regions where there is actual snow (Arctic regions)
# # message("ORYSAT-", METHOD, ": Identifying snow.")
# # For masking out pixels with too many consecutive Snow or too many snow instances
# # in case snw.maxlen is not good enough, the 2 conditions below can be explored
# # It is possible that there's only snow at one time but rice can be planted some other time within the year
# #
# # snw.count <- rowSums(mat.ndvi == -32767) 
# # nsnw.maxlen <- apply(mat.mndwi != -32767, 1, maxlen.condition)
# # snw.maxlen <- apply(mat.mndwi == -32767, 1, maxlen.condition)
# 
# # Initial snow masking
# # rst.lc[snw.maxlen>3] <- -5 # Snow
# 
# 
# # message("ORYSAT-", METHOD, ": Identifying shrubs.")
# # shrb.count <- colSums(mat.lswi[grepl(YEAR, required.acqdates),]<1000, na.rm = TRUE)
# # rst.lc[pixels.toprocess[shrb.count==0 & is.na(rst.lc[pixels.toprocess])]] <- -5
# # 
# # message("ORYSAT-", METHOD, ": Identifying permanent water.")
# # nh2o <- mat.ndvi>=1000 | mat.ndvi >= mat.lswi
# # nh2o.maxlen <- apply(nh2o, 2, maxlen.condition)
# # rst.lc[pixels.toprocess[nh2o.maxlen<15 & is.na(rst.lc[pixels.toprocess])]] <- -6
# # rm(nh2o, nh2o.maxlen, shrb.count, forest.count, na.maxlen)
# # gc(reset=TRUE)
# 
# 
# 
# 
# message("ORYSAT-", METHOD, ": Writing to disk.")
# rst.lc <- writeRaster(rst.lc, paste0(OUTPUT_DIR, "/", paste("LANDCOVER", METHOD, TILE, YEAR, sep="_"),".tif"), datatype="INT4S", overwrite=TRUE)
timeen.rice <- Sys.time()
timedur.rice <- timeen.rice-timest.proc


message("ORYSAT-", METHOD, ": Done. (", round(timedur.rice,2), " ", attr(timedur.rice, "unit"), ") ")
rm(list=ls())

.rs.restartR()
