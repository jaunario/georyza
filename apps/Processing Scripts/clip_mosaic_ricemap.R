library(raster)

WORKSPACE  = "E:/WORKSPACE/Orysat"
YEARS = c(2007, 2016, 2017)
COUNTRY  = "BGD" 
METHOD = "xiao-v1"
OTHER_ID = ""
OUTPUT_DIR = "../SPIA/LATEST_RICE"

ACQDOYS <- seq(from=1,to=361, by=8)

#seasons <- list(aman=6:9, boro=c(1,10:12), aus=2:5) # v1
seasons <- list(aman=6:9, boro=c(1:2,10:12), aus=3:5) # v2
#seasons <- list(aman=7:10, boro=c(1:3,11:12), aus=4:6) # v3

setwd(WORKSPACE)

if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if(!exists("rst.adm")){
  if(!file.exists("BGD_ADM.tif")){
    #shp.admbuff <- rgeos::gBuffer(shp.adm, width = 0.05)
    shp.admsin <- spTransform(shp.adm, CRSobj = CRS(PROJ.MODIS))
    
    rst.adm <- crop(rst.mosaic, shp.admsin)
    rst.adm <- extend(rst.adm, c(2,2))
    rst.adm <- rasterize(shp.admsin, rst.adm)
    rst.adm <- writeRaster(rst.adm, "BGD_ADM.tif", datatype="INT2S")
  } else rst.adm <- raster("BGD_ADM.tif")
} 

if(!exists("rst.utm")){
  if(!file.exists("BGD_UTM495.tif")){
    #shp.admbuff <- rgeos::gBuffer(shp.adm, width = 0.05)
    rst.utm <- projectRaster(rst.adm, crs = PROJ.UTM45N, method="ngb")
    xx <- raster(rst.utm)
    res(xx) <- 495
    rst.utm <- resample(rst.utm, xx)
    rst.utm <- crop(rst.utm, extentFromCells(rst.utm, which(!is.na(rst.utm[]))))
    rst.utm <- writeRaster(rst.utm, "BGD_UTM495.tif", datatype="INT2S",overwrite=TRUE)
  } else rst.utm <- raster("BGD_UTM495.tif")
} 

#shp.adm <- getData(country=COUNTRY, level=1)
for (y in YEARS){
  YEAR <- y
  required.acqdates <- paste("A",c(paste(YEAR-1,sprintf("%03g", ACQDOYS)[(length(ACQDOYS)-7):length(ACQDOYS)], sep=""),  #PREVIOUS YEAR
                                   paste(YEAR,sprintf("%03g", ACQDOYS), sep=""),           # CURRENT YEAR
                                   paste(YEAR+1,sprintf("%03g",ACQDOYS)[1:11], sep="")),   # SUCCEEDING YEAR
                             sep="")
  files.rice <- dir(pattern=paste0("RiceMap_", ifelse(OTHER_ID!="",paste0(OTHER_ID, "_"), ""), METHOD, ".*.",y), ignore.case = TRUE)
  for(ff in files.rice){
    rst.ricemapyear <- raster(ff)
    if(!exists("rst.mosaic")) rst.mosaic <- rst.ricemapyear else rst.mosaic <- merge(rst.mosaic, rst.ricemapyear) 
  }
  
  
  rst.mosaic <- crop(rst.mosaic, rst.adm)
  rst.mosaic <- mask(rst.mosaic, rst.adm)
  
  vals.rice <- values(rst.mosaic)
  cells.rice <- which(!is.na(vals.rice))
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
  

  for(ss in 1:length(seasons)){
    dat.season <- dat.planting[dat.planting$month %in% seasons[[ss]],]
    dat.season$doy <- sapply(dat.season$date, manipulateR::doyFromDate)
    message(names(season)[ss],": Creating Rice Map.")
    rst.season <- raster(rst.mosaic)
    rst.season[dat.season$cell] <- dat.season$doy
    writeRaster(rst.season, filename = paste(paste0(OUTPUT_DIR, "/RiceMap_", ifelse(OTHER_ID!="",paste0(OTHER_ID, "_"), ""), METHOD), COUNTRY, y, names(seasons)[ss], "SIN", sep="_"), format="GTiff", overwrite=TRUE)
    
    message(names(season)[ss],": Reprojecting to WGS84.")
    rst.seasonwgs <- projectRaster(rst.season, crs = projection(shp.adm), method="ngb")
    rst.seasonwgs <- crop(rst.seasonwgs,extentFromCells(rst.seasonwgs,which(!is.na(rst.seasonwgs[]))))
    
    writeRaster(rst.seasonwgs, filename = paste(paste0(OUTPUT_DIR, "/RiceMap_", ifelse(OTHER_ID!="",paste0(OTHER_ID, "_"), ""), METHOD), COUNTRY, y, names(seasons)[ss], "WGS", sep="_"), format="GTiff", overwrite=TRUE)
    
    message(names(season)[ss],": Reprojecting to UTM.")
    rst.seasonutm <- projectRaster(rst.season, crs = PROJ.UTM45N, method="ngb")
    message(names(season)[ss],": Resampling to 495m.")
    rst.seasonutm <- resample(rst.seasonutm, rst.utm, method="ngb")
    message(names(season)[ss],": Writing 495m UTM to disk.")
    rst.seasonutm <- writeRaster(rst.seasonutm, filename = paste0(paste0(OUTPUT_DIR, "/RiceMap_", ifelse(OTHER_ID!="",paste0(OTHER_ID, "_"), ""), METHOD), paste(COUNTRY, y, names(seasons)[ss], "495m_UTM45N", sep="_"), ".tif"), format="GTiff", overwrite=TRUE)
    rm(rst.season, rst.seasonwgs, rst.seasonutm)
  
  }
  rm(rst.mosaic)
}
