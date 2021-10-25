library(raster)

WORKSPACE  = "E:/WORKSPACE/Orysat"
INPUT_DIR  = "../FLAGSHIP/Filled/"
YEARS     = c(2006:2008, 2015:2018)
COUNTRY  = "BGD" 

INDEX_NAME = "LSWI"

ACQDOYS <- seq(from=1,to=361, by=8)
PROJ.MODIS <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
PROJ.UTM45N <- "+proj=utm +zone=45 +datum=WGS84 +units=m +no_defs"

setwd(WORKSPACE)
OUTPUT_DIR = paste("./mosaic", COUNTRY, 'indices_utm', sep="/")
if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
shp.adm <- getData(country=COUNTRY, level=0)

for (y in YEARS){
  YEAR <- y
  required.acqdates <- paste(YEAR,sprintf("%03g", ACQDOYS), sep="")           # CURRENT YEAR
                       
  files.tif <- dir(path=INPUT_DIR, pattern=paste0(YEAR, ".*.", INDEX_NAME, ".tif$"), ignore.case = TRUE, recursive = TRUE, full.names = TRUE)
  for(acq in required.acqdates){
    acq.files <- grep(acq, files.tif)
    message(acq,": Mosaicing.")
    
    for (ff in acq.files){
      rst.tile <- raster(files.tif[ff])
      if(!exists("rst.mosaic")) rst.mosaic <- rst.tile else rst.mosaic <- raster::merge(rst.mosaic, rst.tile) 
    }
    
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
    
    rst.mosaic <- crop(rst.mosaic, rst.adm)
    

    message(acq,": Masking.")
    rst.mosaic <- mask(rst.mosaic, rst.adm)
    
    message(acq,": Reprojecting.")
    #rst.mosaic <- projectRaster(rst.mosaic, crs = projection(shp.adm), method="ngb")
    rst.mosaic <- projectRaster(rst.mosaic, crs = PROJ.UTM45N, method="ngb")
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
    rst.mosaic <- resample(rst.mosaic, rst.utm, method="ngb")
    message(acq,": Writing to disk.")
    rst.mosaic <- writeRaster(rst.mosaic, filename = paste0(OUTPUT_DIR, "/", paste("SMOOTH", COUNTRY, paste0("A",acq), INDEX_NAME, "495m_UTM45N", sep="_"), ".tif"), format="GTiff", overwrite=TRUE)
    
    message(acq,": Done.")
    rm(rst.mosaic)
  }
  
}

rm(list=ls())

.rs.restartR()
