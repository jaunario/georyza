library(raster)

WORKSPACE = "E:/WORKSPACE/Orysat"

AOI_COUNTRY = "Bangladesh"
FAOSTAT_FILE = "C:/Users/jaunario/Downloads/FAOSTAT_data_3-24-2020.csv"
LEGACY_RICEMAP_FILE = "C:/Data/Raster/rice_modis.tif"
RICEMAP_FNAME = "RiceMap_Xiao-v1_BGD" # "Ricemap_NOWATERCOUNT_BGD"

YEARS = 2010

setwd(WORKSPACE)

dat.faorice <- read.csv(FAOSTAT_FILE, stringsAsFactors = FALSE)

dat.aoiricearea <- dat.faorice[dat.faorice$Area==AOI_COUNTRY & dat.faorice$Element=="Area harvested", ]


# validation with FAOSTAT Rice Area Harvested
dat.statscompare <- data.frame(year=numeric(0), ricemap_area=numeric(0), faostat_area=numeric(0))
for(h in 1:length(YEARS)){
  YEAR <- YEARS[h]
  files.seasonalrice <- dir(pattern=paste0(RICEMAP_FNAME, "_", YEAR, ".*.tif"), full.names = TRUE)
  val.total <- 0
  for(i in 1:length(files.seasonalrice)){
    rst.thisseason <- raster(files.seasonalrice[i])
    #val.total <- val.total + sum(!is.na(rst.thisseason[]))*463.3127^2/10000
    rst.area <- raster::area(rst.thisseason)
    rst.area <- mask(rst.area, rst.thisseason)
    val.total <- val.total + sum(rst.area[], na.rm = TRUE)*100
  }
  dat.statscompare[h,] <- rbind(YEAR, val.total, dat.aoiricearea$Value[dat.aoiricearea$Year==YEAR])
}

#shp.adm <- raster::getData(country="BGD", level=2)
shp.adm <- shapefile("APY_BGD.shp")
dat.districArea <- data.frame(year=numeric(0), ricemap_area=numeric(0), faostat_area=numeric(0))
dat.distRA <- shp.adm@data
for(h in 1:length(YEARS)){
  YEAR <- YEARS[h]
  files.seasonalrice <- dir(pattern=paste0(RICEMAP_FNAME, "_", YEAR, ".*.tif"), full.names = TRUE, recursive=FALSE)
  files.seasonalrice <- dir(path="G:/My Drive/Projects/SPIA/SPIA_Flood/2_Process/JA_RiceMaps/By Season",pattern=paste0(RICEMAP_FNAME, "_", YEAR, ".*.tif"), full.names = TRUE, recursive=TRUE)
  val.total <- 0
  for(i in 1:length(files.seasonalrice)){
    rst.thisseason <- raster(files.seasonalrice[i])
    #val.total <- val.total + sum(!is.na(rst.thisseason[]))*463.3127^2/10000
    rst.area <- raster::area(rst.thisseason)
    if(!exists("lst.dstta")){
      lst.dstta <- extract(rst.area, shp.adm)
      dat.distRA[,paste0("LANDAREA")] <- sapply(lst.dstta, sum, na.rm=TRUE)*100
    }
    rst.area <- mask(rst.area, rst.thisseason)
    lst.distra <- extract(rst.area,shp.adm)
    
    dat.distRA[,paste0("RA",YEAR,"AM_v2fin")] <- sapply(lst.distra, sum, na.rm=TRUE)*100
  }
  dat.statscompare[h,] <- rbind(YEAR, val.total, dat.aoiricearea$Value[dat.aoiricearea$Year==YEAR])
}

write.csv(dat.distRA[,c(1:4, 42, 46:47)], "./BGD_DISTRICT_RiceArea_Aman_2016v2fin.csv")
pracma::rmserr(as.numeric(dat.distRA$All_Ha_16_), dat.distRA$RA2016AM_v2fin, summary = TRUE)

# Accuracy assessment
# 2010
shp.points <- shapefile("C:/WORKSPACE/SPIA/Validation/BGD_BrahmaputraRegion_2017_AMAN_RnR_SRSP_VGrid.shp")

rst.rice <- raster("G:/My Drive/Projects/SPIA/SPIA_Flood/2_Process/JA_RiceMaps/By Season/Ricemap_BGD_2016_aman.tif")
rst.rice <- raster("E:/WORKSPACE/Orysat/BGD/ModXIaoV2-S2/RiceMap_ModXiao_BGD_2016_aman.tif")
rst.rice <- raster("E:/WORKSPACE/Orysat/Trend_0.1_MOD09A1/MOSAIC_BGD_500m/RiceMap_Trend_BGD_EOS_2017_AMAN_WGS84.tif")
rst.rice <- projectRaster(rst.rice, crs=projection(shp.points))
shp.points <- spTransform(shp.points, projection(rst.rice))
vp.points <- rgeos::gCentroid(shp.points, byid = TRUE)

pts.modra <- data.frame(cell=cellFromXY(rst.rice, vp.points))
pts.modra$value <- rst.rice[pts.modra$cell]
pts.modra$rice <- !is.na(pts.modra$value)
pts.modra$gtRnR <- as.numeric(shp.points@data$RnR)
pts.modra <- unique(pts.modra)
#pts.modra <- pts.modra[!duplicated(pts.modra[,1]),]
table(pts.modra$rice,pts.modra$gtRnR)
table(pts.modra$rice,pts.modra$gtRnR)/nrow(pts.modra)*100
all.yrs  <- cbind(2010, pts.modra)


# 2015
shp.points <- shapefile("G:/My Drive/Projects/SPIA/SPIA_Flood/1_Source/RMS/Shapefile/RMS_2016_Bangladesh_Aman_2015_PlotLevel.shp")
rst.rice <- raster("./BGD/Xiao_v1-S1/Ricemap_BGD_2015_aman.tif")
pts.modra <- extract(rst.rice,shp.points, cellnumbers=TRUE)
pts.modra <- cbind(unique(pts.modra),1)
table(pts.modra[,3]==1,!is.na(pts.modra[,2]))
table(pts.modra[,3]==1,!is.na(pts.modra[,2]))/nrow(pts.modra)*100
all.yrs <- rbind(all.yrs,cbind(2015,pts.modra))

# 2016
shp.points <- shapefile("G:/My Drive/Projects/SPIA/SPIA_Flood/1_Source/RMS/Shapefile/RMS_2017_Bangladesh_Aman_2016.shp")
rst.rice <- raster("./BGD/Xiao_v1-S1/Ricemap_BGD_2016_aman.tif")
pts.modra <- extract(rst.rice,shp.points, cellnumbers=TRUE)
pts.modra <- cbind(unique(pts.modra),1)
table(pts.modra[,3]==1,!is.na(pts.modra[,2]))
table(pts.modra[,3]==1,!is.na(pts.modra[,2]))/nrow(pts.modra)*100
all.yrs <- rbind(all.yrs,cbind(2016,pts.modra))


shp.points$rs.rice <- extract(rst.rice, shp.points)
shp.ptsin$rs.rice2 <- extract(rst.season, shp.ptsin)

dat.evi <- extract(stk.evi, shp.ptsinu)
colnames(dat.evi) <- required.acqdates


dat.acc <- data.frame(ground=shp.ptsin$rice, modis=as.numeric(!is.na(shp.ptsin$rs.rice2)))
dat.acc$commission <- dat.acc$ground==0 & !is.na(dat.acc$modis)
dat.acc$ommission <- dat.acc$ground==1 & is.na(dat.acc$modis)
table(dat.acc$ground==1, dat.acc$modis==1)
