# Title: Mapping potential heat-stress affected rice areas
# Author: Jorrel Khalil Aunario, Alice Laborte
# Date Created: March 18, 2019
# Version: 0.101
# Description: 
# 
#


# SETTINGS
WORKSPACE  = paste(normalizePath(Sys.getenv("GDRIVE"),winslash="/"), "Projects/Heat Stress Models", sep = "/")
INPUT_DIR  = "Models/Inputs"
OUTPUT_DIR = "Projects/Heat Stress Models/Output/DaytimeHeat/NP2.30arcmin"

# WEATHER DATASET SETTINGS
WEATHERSET = "NP2" # NASA-POWER v2
WEATHERSCHEMA = "power30m_asia"

# Version Info and File Naming
VERSION = 0.101
RUNNO   = 1
RUNDATE = Sys.Date()
OUTDEST  = "Auto"

# Geographics Settings
AOI_LEVEL = "region"
OUTRES   = 0.5

# Crop-related Paramters
CROP_LOWMAT =  FALSE # If TRUE, adjusts harvest doy to reflect low-maturity variety
RICE_SEASON = "PK3"

#INPUT 
FILE_RICEMASK = paste(INPUT_DIR, "rice_10min_new.tif", sep="/")
FILE_ADMRASTER = system.file("data/WORLD_ID_0_30.tif",package="manipulateR")

AOI_NAME = "Asia"
YEAR_ST  = 1983
YEAR_EN  = 2017

library(RiceStress)

setwd(WORKSPACE)

if (OUTDEST=="Auto"){
  OUTPUT_DIR <- paste0("./Output/DTH_v", VERSION, "_", WEATHERSET, "_", format(RUNDATE, "%Y%m%d"))
  force.directories(OUTPUT_DIR, recursive=TRUE)
}

# TODO: save settings
years <- YEAR_ST:YEAR_EN
#TODO: yrs may not be needed 
yrs <- paste(YEAR_ST,YEAR_EN, sep="-")



FILE_PLANTDOY = paste(INPUT_DIR, "/WORLD_PLANT_", RICE_SEASON, "_", OUTRES*60, ".tif", sep="")
FILE_HARVDOY  = paste(INPUT_DIR, "/WORLD_HARV_", RICE_SEASON, "_", OUTRES*60, ".tif", sep="")

if(!file.exists(FILE_PLANTDOY)) stop("Planting day raster not found")
if(!file.exists(FILE_HARVDOY)) stop("Harvest day raster not found")

rst_mask <- raster(FILE_RICEMASK)
if(res(rst_mask)[1]<OUTRES) {
  tempraster <- raster()
  res(tempraster) <- OUTRES
  rst_mask <- raster::resample(rst_mask,tempraster, method="ngb")
  rm(tempraster)
}

rst_plantdoy <- raster(FILE_PLANTDOY)
rst_plantdoy <- mask(rst_plantdoy,rst_mask) 

rst_harvdoy <- raster(FILE_HARVDOY)
rst_harvdoy <-  mask(rst_harvdoy, rst_mask)

# TODO: Double check calculations
if (CROP_LOWMAT){ # Adjust harvest doy to low-maturity variety
	rst_harvdoy <- calc(rst_harvdoy, fun=function(x){return(ifelse(x > 0, x, NA) )})
	rst_harvdoy <- overlay(rst_plantdoy, rst_harvdoy, fun=function(x,y){return(ifelse(y > x, y - x, y - x + 365) )} )
	rst_harvdoy  <- calc(rst_harvdoy, fun=function(x){return(ifelse(x < 70, 70, x) )} )
	rst_harvdoy  <- calc(rst_harvdoy, fun=function(x){return(ifelse(x > 200, 200, x) )} )
	# floweringdate <- dateFromDoy(rst_plantdoy[poi[i]]+dae[poi[i]]-30,y)
	rst_harvdoy <- rst_plantdoy+rst_harvdoy
}

#
rst_adm <- raster(FILE_ADMRASTER)
poi <- switch(AOI_LEVEL, region=admAttribCells(filter=AOI_NAME, admattrib="UNREGION2", refraster=FILE_ADMRASTER), country= admAttribCells(filter=AOI_NAME, admattrib="NAME_ENGLI", refraster=FILE_ADMRASTER))
poi <- poi[which(!is.na(rst_mask[poi]))]
poi <- poi[which(!is.na(rst_harvdoy[poi]))]

wthvars <- list(c("tmax", "tmin", "tdew"))
names(wthvars) <- WEATHERSET
con <- odbcConnect("climate_SRV3A")
stress.dat <- vector()
for ( i in 1:length(poi)){
	show.message("Analyzing pixel ", poi[i], "\r", EL=TRUE, appendLF=FALSE)
	if(is.na(rst_plantdoy[poi[i]])|rst_plantdoy[poi[i]]<=0) next
  query <- paste0("SELECT wdate date, ", paste0(wthvars[[WEATHERSET]],"/100 ", wthvars[[WEATHERSET]],collapse = ", "), " FROM ", WEATHERSCHEMA, ".c", poi[i], " WHERE wdate BETWEEN '", YEAR_ST, "-1-1' AND '", YEAR_EN, "-12-31'")
  dat.wth <- sqlQuery(con, query)
	xy <- xyFromCell(rst_adm,poi[i])
	
	#wth <- fetch.daily(xy=xy, srcvars=wthvars, connection=con)[[1]]
  # tmin should be <= tdew
	dat.wth$tdew[dat.wth$tdew<dat.wth$tmin] <- dat.wth$tmin[dat.wth$tdew<dat.wth$tmin]
	wth <- new("weather")
	wth@lon <- xy[1,"x"]
	wth@lat <- xy[1,"y"]
	wth@w   <- dat.wth
	
	for (y in years){
		floweringdate <- dateFromDoy(rst_harvdoy[poi[i]]-30,ifelse(rst_harvdoy[poi[i]]<rst_plantdoy[poi[i]],y+1,y))
		if(is.na(floweringdate)) next		
		if((floweringdate-15)<as.Date(paste0(YEAR_ST,"-1-1"))|(floweringdate+21)>as.Date(paste0(YEAR_EN,"-12-31"))) {
			show.message(poi[i], ": Required data exceeds coverage of database -", format(floweringdate,"%x"), appendLF=TRUE)
			next
		}
		stress.dat <- rbind(stress.dat, cbind(poi[i], y,ricestress.dayheat(wth=wth,flwrdate=floweringdate, sf.method=NA, tmax.threshold=35)))		
	}
}
con <- odbcClose(con)
colnames(stress.dat)[1:2] <- c("cell", "year")
stress.dat <- as.data.frame(stress.dat)
write.csv(stress.dat, paste0(OUTPUT_DIR,"/stressdat_", RICE_SEASON, "_", AOI_NAME, "_", YEAR_ST, "-", YEAR_EN, ".csv"))
#stress.dat <- read.csv(paste0(OUTPUT_DIR,"/stressdat_", RICE_SEASON, "_", AOI_NAME, "_", YEAR_ST, "-", YEAR_EN, ".csv"), stringsAsFactors = FALSE)
#stress.dat <- stress.dat[,-1]
ricemask_orig <- raster(FILE_RICEMASK)

# Mapping stress data

for (y in years){		
	thisyr <- stress.dat[stress.dat$year==y,]
	for (x in 3:ncol(thisyr)){
	  if(x<12) next
		base <- raster(rst_harvdoy)
		base[thisyr$cell] <- thisyr[,x]
		#base <- smooth.raster(base, ricemask_orig)
		base.crop <- crop(base,cellListExtent(rst_harvdoy, poi,adjust=2))
		dt <- ifelse(grepl("sum",colnames(stress.dat)[x], ignore.case=TRUE)|grepl("stress",colnames(stress.dat)[x], ignore.case=TRUE),"INT4S","FLT4S")
		base.crop <- writeRaster(base.crop, filename=paste(OUTPUT_DIR, "/", paste(toupper(colnames(thisyr)[x]), RICE_SEASON, AOI_NAME, y, sep="_"),".tif", sep=""), datatype=dt, format="GTiff", options="COMPRESS=LZW",overwrite=TRUE)
	}				
} 


for (x in 6:ncol(stress.dat)){	
	fun <- ifelse(grepl("stress",colnames(stress.dat)[x], ignore.case=TRUE), "sum", "mean")
	overall <- aggregate(x=stress.dat[,x], by=list(cell=stress.dat$cell), FUN=fun, na.rm=TRUE)
	if(sum(!is.na(overall[,2]))==0|sum(overall[,2],na.rm=TRUE)==0) next
	base <- raster(rst_harvdoy)
	base[overall$cell] <- overall[,2]
	base.crop <- crop(base,cellListExtent(rst_harvdoy, poi, adjust=2))	
	base.crop <- writeRaster(base.crop, filename=paste(OUTPUT_DIR, "/", paste(toupper(colnames(stress.dat)[x]), RICE_SEASON, AOI_NAME, fun, YEAR_ST, YEAR_EN ,sep="_"),".tif", sep=""), format="GTiff", options="COMPRESS=LZW",overwrite=TRUE, datatype="INT4S")
	base.smooth <- smooth.raster(base, ricemask_orig)	
	base.smooth <- crop(base.smooth,cellListExtent(rst_harvdoy, poi, adjust=2))	
	base.smooth <- writeRaster(base.smooth, filename=paste(OUTPUT_DIR, "/", paste(toupper(colnames(stress.dat)[x]), RICE_SEASON, AOI_NAME, fun, "smooth", YEAR_ST, YEAR_EN ,sep="_"),".tif", sep=""), format="GTiff", options="COMPRESS=LZW",overwrite=TRUE, datatype="INT4S")
	
}


stress.tmax <- raster("Asia_1983-2012_Stress.Based.on.TMax.tif")
plot(stress.tmax, col=brewer.pal(9,"Reds"),zlim=c(1,30))
plot(world, add=TRUE, border="light grey")

stress.tpan <- raster("Asia_1983-2012_Stress.Based.on.MaxTPan.tif")
plot(stress.tpan, col=brewer.pal(9,"Reds"),zlim=c(1,30))
plot(world, add=TRUE, border="light grey")
stress.sf <- raster("Asia_1983-2012_Stress.Based.on.Spikelet.Fertility.tif")
plot(stress.sf, col=brewer.pal(9,"Reds"),zlim=c(1,30))
plot(world, add=TRUE, border="light grey")


stress.map <- function(con, plantraster, harvraster, model=ricestress.dayheat,...){
	con <- odbcConnect("climate_SRV3A")
} 

# TODO: Investigate code below... Probably for debugging
for (y in years){
	floweringdate <- dateFromDoy(rst_plantdoy[11109]+dae[11109]-30,y)
	if(is.na(floweringdate)) next		
	if(floweringdate-15<as.Date("1983-1-1")|floweringdate-15>as.Date("2014-12-31")) {
		show.message(poi[i], ": Required data exceeds coverage of database -", format(floweringdate,"%x"), appendLF=TRUE)
		next
	}
	stress.dat <- rbind(stress.dat, cbind(poi[i], y,ricestress.dayheat(wth=wth@w,flwrdate=floweringdate, method="Jagadish")))		
}