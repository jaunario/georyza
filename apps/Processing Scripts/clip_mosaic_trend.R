library(raster)
library(manipulateR)

WORKSPACE  = "E:/WORKSPACE/FLAGSHIP/VWIDynamics_0.1_MOD09A1"
OUTPUT_DIR = "AUTO"
YEARS = 2007
COUNTRY  = "BGD" 
METHOD = "VWIDynamics"
OTHER_ID = ""

ACQDOYS <- seq(from=1,to=361, by=8)

#seasons <- list(aman=6:9, boro=c(1,10:12), aus=2:5) # v1
seasons <- list(aman=6:9, boro=c(1:2,10:12), aus=3:5) # v2
#seasons <- list(aman=7:10, boro=c(1:3,11:12), aus=4:6) # v3

setwd(WORKSPACE)
if(OUTPUT_DIR=="AUTO") {
  OUTPUT_DIR <- paste0("./",paste("MOSAIC", COUNTRY, "500m", sep="_"))
  if(!dir.exists(OUTPUT_DIR)) force.directories(OUTPUT_DIR)
}

files.rice <- dir(pattern=paste0("RiceMap_", ifelse(OTHER_ID!="",paste0(OTHER_ID, "_"), ""), METHOD, "_.*.tif"), ignore.case = TRUE)
dat.files <- do.call(rbind, sapply(basename(files.rice),strsplit, split="_", USE.NAMES = FALSE))
colnames(dat.files) <- c("Prod", "Method", "Tile", "Year", "Quarter", "Variable")
dat.files <- data.frame(dat.files, stringsAsFactors = FALSE)
dat.files$Year <- as.integer(dat.files$Year)
dat.files$Variable <- sub(".tif", "", dat.files$Variable)
dat.files$filename <- files.rice

shp.adm <- getData(country=COUNTRY, level=1)

for (y in YEARS){
  
  for (q in 1:4){
    
    for (v in unique(dat.files$Variable)){
      message(v, ": Merging")
      tomerge <- subset(dat.files, Year==y & Quarter==paste0("Q",q) & Variable==v)
      
      for(ff in tomerge$filename){
        rst.ricemapyear <- raster(ff)
        if(!exists("rst.mosaic")) rst.mosaic <- rst.ricemapyear else rst.mosaic <- merge(rst.mosaic, rst.ricemapyear) 
      }
      
      if(!exists("shp.admsinu")) {
        shp.admsinu <- spTransform(shp.adm, CRSobj = CRS(projection(rst.ricemapyear)))
        rst.adm <- rasterize(shp.admsinu,rst.mosaic)
        rst.adm <- crop(rst.adm, shp.admsinu)
      }                                                   
      
      rst.mosaic <- crop(rst.mosaic, rst.adm)
      rst.mosaic[rst.mosaic[]==0] <- NA
      rst.mosaic <- mask(rst.mosaic, rst.adm)
      #rst.mosaic <- crop(rst.mosaic, manipulateR::cellListExtent(rst.mosaic,which(!is.na(rst.mosaic[]))))
      
      vals.rice <- values(rst.mosaic)
      cells.rice <- which(!is.na(vals.rice))
      vals.rice <- vals.rice[cells.rice]
      if (v %in% c("EOS","SOS")){
        doy.rice <- vals.rice
        vals.rice <- dadateFromDoy(vals.rice, y)
        mo.rice <- monthFromDate(vals.rice)
        doy.rice <- vals.rice
        #as.integer(julian(vals.rice,as.Date(paste0(y,"-1-1"))-1))
        dat.v <- data.frame(cells.rice, doy.rice, mo.rice)
        colnames(dat.v) <- c("cell", paste0(v,"_doy"), paste0(v,"_mo"))
        rm(mo.rice, doy.rice)
      } else {
        dat.v  <- data.frame(cell=cells.rice, vals.rice)
        colnames(dat.v) <- c("cell", v)
      }
      
      if(!exists("dat.q")) dat.q <- dat.v else dat.q <- plyr::join(dat.q, dat.v, by="cell")
      rm(dat.v,rst.mosaic)
    }
    gc(reset=TRUE)
    
    dat.qs1 <- dat.q[dat.q$SOS_mo %in% seasons[[1]],]
    dat.qs2 <- dat.q[dat.q$SOS_mo %in% seasons[[2]],]
    dat.qs3 <- dat.q[dat.q$SOS_mo %in% seasons[[3]],]
    
    if(!exists("dat.s1")) dat.s1 <- dat.qs1 else dat.s1 <- rbind(dat.s1, dat.qs1) 
    if(!exists("dat.s2")) dat.s2 <- dat.qs2 else dat.s2 <- rbind(dat.s2, dat.qs2) 
    if(!exists("dat.s3")) dat.s3 <- dat.qs3 else dat.s3 <- rbind(dat.s3, dat.qs3) 
    rm(dat.q, dat.qs1, dat.qs2, dat.qs3)
  }
  #writeRaster(rst.mosaicwgs, filename = paste(paste0("RiceMap_", ifelse(OTHER_ID!="",paste0(OTHER_ID, "_"), ""), METHOD), COUNTRY, y, "SIN", sep="_"), format="GTiff", overwrite=TRUE)
  
  dat.seasons <- list(AMAN=dat.s1,BORO=dat.s2, AUS=dat.s3)
  rm(list = paste0("dat.s",1:3))
  for (i in 1:length(dat.seasons)){
    dat.s <- dat.seasons[[i]]
    rst.sos <- raster(rst.adm)
    rst.sos[dat.s$cell] <- dat.s$SOS_doy
    rst.sos <- writeRaster(rst.sos, filename=paste0(OUTPUT_DIR,"/", paste("RiceMap",METHOD,COUNTRY, "SOS", y, names(dat.seasons)[i], "SIN", sep="_"),".tif"), datatype="INT2S")
    rst.sos <- projectRaster(rst.sos, crs = projection(shp.adm))
    rst.sos <- crop(rst.sos, shp.adm)
    rst.sos <- writeRaster(rst.sos, filename=paste0(OUTPUT_DIR,"/", paste("RiceMap",METHOD,COUNTRY, "SOS", y, names(dat.seasons)[i], "WGS84", sep="_"),".tif"), datatype="INT2S")
    rst.eos <- raster(rst.adm)
    rst.eos[dat.s$cell] <- dat.s$EOS_doy
    rst.eos <- writeRaster(rst.eos, filename=paste0(OUTPUT_DIR,"/", paste("RiceMap",METHOD,COUNTRY, "EOS", y, names(dat.seasons)[i], "SIN", sep="_"),".tif"), datatype="INT2S")
    rst.eos <- projectRaster(rst.eos, crs = projection(shp.adm))
    rst.eos <- crop(rst.eos, shp.adm)
    rst.eos <- writeRaster(rst.eos, filename=paste0(OUTPUT_DIR,"/", paste("RiceMap",METHOD,COUNTRY, "EOS", y, names(dat.seasons)[i], "WGS84", sep="_"),".tif"), datatype="INT2S")
    # rst.rsq <- raster(rst.adm)
    # rst.rsq[dat.s$cell] <- dat.s$RSQ
    # rst.rsq <- writeRaster(rst.rsq, filename=paste0(OUTPUT_DIR,"/", paste("RiceMap",METHOD,COUNTRY, "RSQ", y, names(dat.seasons)[i],sep="_"),".tif"), datatype="INT2S")
  }
  
  # plot(rst.sos)
  # plot(rst.sos==6)
  # 
  # 
  # dat.planting <- cbind(rep(cells.rice,3),c(season1,season2,season3))
  # dat.planting <- na.omit(dat.planting)
  # colnames(dat.planting) <- c("cell","acq.index")
  # dat.planting <- as.data.frame(dat.planting)
  # dat.planting$acqdate <- required.acqdates[dat.planting$acq.index]
  # dat.planting$date <- as.Date(dat.planting$acqdate, "A%Y%j")
  # dat.planting$month <- manipulateR::monthFromDate(dat.planting$date)
  # 
  # for(ss in 1:length(seasons)){
  #   dat.season <- dat.planting[dat.planting$month %in% seasons[[ss]],]
  #   dat.season$doy <- manipulateR::doyFromDate(dat.season$date)
  #   rst.season <- raster(rst.mosaic)
  #   rst.season[dat.season$cell] <- dat.season$doy
  #   writeRaster(rst.season, filename = paste(paste0("RiceMap_", ifelse(OTHER_ID!="",paste0(OTHER_ID, "_"), ""), METHOD), COUNTRY, y, names(seasons)[ss], "SIN", sep="_"), format="GTiff", overwrite=TRUE)
  # 
  #   rst.seasonwgs <- projectRaster(rst.season, crs = projection(shp.adm), method="ngb")
  #   rst.seasonwgs <- crop(rst.seasonwgs, manipulateR::cellListExtent(rst.seasonwgs,which(!is.na(rst.seasonwgs[]))))
  #   
  #   writeRaster(rst.seasonwgs, filename = paste(paste0("RiceMap_", ifelse(OTHER_ID!="",paste0(OTHER_ID, "_"), ""), METHOD), COUNTRY, y, names(seasons)[ss], "WGS", sep="_"), format="GTiff", overwrite=TRUE)
  #   rm(rst.season, rst.seasonwgs)
  # }
  # rm(rst.mosaic)
}
