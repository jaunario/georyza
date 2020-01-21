library(orynfect)
library(GloCR)
library(raster)

WEATHER_DIR = "files/PHL"
RICE_SOS    = "files/CropCal/WORLD_PLANT_PK1_5.tif"
YEARS       = 2010:2018

files.wth <- dir(WEATHER_DIR, pattern=".csv", full.names = TRUE)

rst.riceplant <- raster(RICE_SOS)
norice <- vector()

dat.infection <- vector()
for (i in 1:length(files.wth)){
  info.wth <- unlist(strsplit(basename(files.wth[i]),"_"))
  xy <- c(as.numeric(info.wth[3]),as.numeric(info.wth[4]))
  cell <- cellFromXY(rst.riceplant,xy)
  doy.emergence <- rst.riceplant[cell]
  
  if(is.na(doy.emergence)| cell %in% norice) {
    norice <- unique(c(norice,cell))
    next
  }

  dat.wth <- read.csv(files.wth[i], stringsAsFactors = FALSE)

  dat.wth$date <- as.Date(dat.wth$date)
  colnames(dat.wth) <- c("date", "tmax", "tmin", "prec", "srad", "rhmax","rhmin", "wind.morningmax", "wind.daymax", "wind.avg")
  dat.wth$tavg <- (dat.wth$tmax + dat.wth$tmin)/2

  wth <- data.frame(date=dat.wth$date,tavg=dat.wth$tavg,rhmax=dat.wth$rhmax,prec=dat.wth$prec)

  for(yr in YEARS){
    message("Cell-",cell, "_Year-",yr)

    audpc <- leafBlast(wth=wth, crop.estabdate=as.character(dateFromDoy(doy.emergence,yr)))

    print(audpc)
  }
  
}