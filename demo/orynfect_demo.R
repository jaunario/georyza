library(orynfect)
library(GloCR)

WEATHER_DIR = "G:/My Drive/Shared Data/AWHERE/daily/PHL"
RICE_SOS    = "C:/Data/Raster/CropCal/WORLD_PLANT_PK1_5.tif"
YEARS       = 2010:2018

infect <- function(disease=leafBlast, summary.fun=sum, field="severity", ...){
  infection <- disease(...)@d
  return(summary.fun(infection[,field]))
}

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

  wth <- new("weather")
  wth@w <- dat.wth
  
  for(yr in YEARS){
    message("Cell-",cell, "_Year-",yr)
    infection <- data.frame(cell, year=yr, season="main", audpc=infect(wth=wth, emergence=dateFromDoy(doy.emergence,yr)))
    dat.infection <- rbind(dat.infection, infection)
  }  
  
}

dat.infection <- cbind(dat.infection, xyFromCell(rst.riceplant, dat.infection$cell)) 
dat.infection$audpc.class <- cut(dat.infection$audpc,breaks=seq(0,450,length.out = 9))
ggplot()+ geom_tile(data = dat.infection, aes(x=x,y=y, fill=audpc.class )) +facet_wrap(~year) + scale_fill_brewer(palette = "Reds") + theme_dark() + theme(legend.position = "bottom")

# shp.phl <- getData(country="PHL", level=0)
# if(exists("stk.infection")) rm(stk.infection)
# for(yr in YEARS){
#   message("Year-",yr)
#   dat.yrinect <- dat.infection[dat.infection$year==yr,]
#   rst.infection <- raster(rst.riceplant)
#   rst.infection[dat.yrinect$cell] <- dat.yrinect$audpc
#   rst.infection <- crop(rst.infection, shp.phl)
#   dat.gginfect <- as.data.frame(stk.infection, xy=TRUE)
#   
#   #if(!exists("stk.infection")) stk.infection <- raster::stack(rst.infection) else stk.infection <- addLayer(stk.infection, rst.infection)
# #  plot(rst.infection, col=colorRampPalette(c("white", "#aa0000"))(25), main=yr, zlim=c(0,400))
# #  plot(shp.phl, add=TRUE, border="#aaaaaa")
# #  Sys.sleep(1)
# }  
# names(stk.infection) <- paste("Year", YEARS)
# levelplot(stk.infection, layout=c(3,3), par.settings=BuRdTheme)
# st <- Sys.time()
# yy <- leafBlast(wth, emergence = wth@w$date[2])@d
# Sys.time()-st
