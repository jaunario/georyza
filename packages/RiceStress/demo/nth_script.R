library(RiceStress)

# SETTINGS
low.maturity <- FALSE 
adm.level <- "region"
inpath <- paste(normalizePath(Sys.getenv("GDRIVE"),winslash="/"), "Projects/Heat Stress Models/Models/Inputs", sep = "/")
WORKSPACE = paste(normalizePath(Sys.getenv("GDRIVE"),winslash="/"), "Projects/Heat Stress Models/NightTime/30min", sep = "/")
force.directories(WORKSPACE, recursive=TRUE)
setwd(WORKSPACE)

#INPUT 

coi <- "WORLD"
#coi <- "Pakistan"
years <- 1987:2016
yrs <- paste(min(years),max(years), sep="-")

m <- raster(paste(inpath, "rice_1d_new2.tif", sep="/"))
ricemask_orig <- raster(paste(inpath, "rice_10min_new.tif", sep="/"))
m <- aggregate(ricemask_orig, fact=res(praster)[1]/res(ricemask_orig)[1])
praster <- raster(paste(inpath, "/WORLD_PLANT_PK1_30.tif", sep=""))
praster <- mask(praster,m) 

hraster <- raster(paste(inpath, "/WORLD_HARV_PK1_30.tif", sep=""))
hraster <-  mask(hraster, m)
if (low.maturity){
	hraster <- calc(hraster, fun=function(x){return(ifelse(x > 0, x, NA) )})
	hraster <- overlay(praster, hraster, fun=function(x,y){return(ifelse(y > x, y - x, y - x + 365) )} )
	hraster  <- calc(hraster, fun=function(x){return(ifelse(x < 70, 70, x) )} )
	hraster  <- calc(hraster, fun=function(x){return(ifelse(x > 200, 200, x) )} )
	# floweringdate <- dateFromDoy(praster[poi[i]]+dae[poi[i]]-30,y)
	hraster <- praster+hraster
}
if (!is.na(adm.level)) {
  poi <- switch(adm.level, 
                region=admAttribCells(filter=coi, admattrib="UNREGION2"), 
                country=admAttribCells(filter=coi, admattrib="NAME_ENGLI"))
  poi <- poi[which((!is.na(hraster[poi])|hraster[poi]>0))]  
} else {
  poi <- which(!is.na(praster[]))
}

outdir <- WORKSPACE #paste(normalizePath(Sys.getenv("GDRIVE"),winslash="/") , "Projects/Heat Stress Models/Output/uncor.35.direct.sol",coi,sep="/")
#force.directories(outdir, recursive=TRUE)
#TODO: Fix TS.daily in GloCR Fetch
TS.daily =TRUE
wthvars <- list(nasacomp_1d=c("tmax", "tmin", "tdew"))
con <- odbcConnect("climate_SRV3A")
stress.dat <- vector()
fail <- vector()
for ( i in 1229:length(poi)){
	show.message("Analyzing pixel ", poi[i], "\r", EL=TRUE, appendLF=FALSE)
	if(is.na(praster[poi[i]])|praster[poi[i]]<=0) next
	xy <- xyFromCell(praster,poi[i])
	wth <- withRetry(get.power(xy[,1],xy[,2], vars=c("T2MDEW", "T2M_MIN", "T2M_MAX"),stdate = as.Date(paste(min(years),"-1-1",sep="")), endate=as.Date(paste(max(years),"-12-31",sep="")), np.site="https://power.larc.nasa.gov/cgi-bin/v1/DataAccess.py?request=execute&identifier=SinglePoint"))
	if(class(wth)!="weather") {
	  fail <- c(fail, poi[i])
	  next
	}
	colnames(wth@w) <- c("date", "tdew", "tmax", "tmin")
	
	#wth <- fetch.daily(xy=xy, srcvars=wthvars, connection=con)
	#wth@w$tdew[wth@w$tdew<wth@w$tmin] <- wth@w$tmin[wth@w$tdew<wth@w$tmin] 
	
	for (y in years){
	  period.dates <- ricestage.period(praster[poi[i]], hraster[poi[i]],year = y, stage = rp.repr)
		#floweringdate <- dateFromDoy(hraster[poi[i]]-30,y)
		#if(is.na(floweringdate)) next		
		# if(floweringdate-15<as.Date("1983-1-1")|floweringdate-15>as.Date("2014-12-31")) {
		# 	show.message(poi[i], ": Required data exceeds coverage of database -", format(floweringdate,"%x"), appendLF=TRUE)
		# 	next
		# }
		stress.dat <- rbind(stress.dat, data.frame(cell=poi[i], year=y,t(ricestress.nighttimeHeat(wth=wth,dstart=period.dates$start[1], dend = period.dates$end[1], consecutive = FALSE))))		
	}
}
con <- odbcClose(con)
#colnames(stress.dat)[1:2] <- c("cell", "year")
stress.dat <- as.data.frame(stress.dat)
write.csv(stress.dat, paste(outdir,"out_notconsecutive.csv",sep="/"))
#ricemask_orig <- raster(paste(inpath, "rice_10min_new.tif", sep="/"))

# Mapping stress data

for (y in years){		
	thisyr <- stress.dat[stress.dat$year==y,]
	for (x in 3:ncol(thisyr)){
		base <- raster(hraster)
		base[thisyr$cell] <- thisyr[,x]
		base <- smooth.raster(base, ricemask_orig)
		base <- crop(base,cellListExtent(poi,hraster, adjust=2))
		dt <- ifelse(grepl("hotdays",colnames(stress.dat)[x], ignore.case=TRUE)|grepl("stress",colnames(stress.dat)[x], ignore.case=TRUE),"INT4S","FLT4S")
		base <- writeRaster(base, filename=paste(outdir, "/", paste(colnames(thisyr)[x],coi,y,sep="_"),".tif", sep=""), datatype=dt, format="GTiff", options="COMPRESS=LZW",overwrite=TRUE)
	}				
} 


for (x in 3:ncol(stress.dat)){	
  fun <- ifelse(grepl("stress",colnames(stress.dat)[x], ignore.case=TRUE)|grepl("cnt",colnames(stress.dat)[x], ignore.case=TRUE), "sum", "mean")
  overall <- aggregate(x=stress.dat[,x], by=list(cell=stress.dat$cell), FUN=fun)
  overall <- aggregate(stressed~cell, stress.dat, FUN=fun)
  if(sum(!is.na(overall[,2]))==0|sum(overall[,2],na.rm=TRUE)==0) next
  dt <- ifelse(fun=="sum","INT4S","FLT4S")
  
  base <- raster(hraster)
  base[overall$cell] <- overall[,2]
  b1d.crop <- crop(base,cellListExtent(poi,hraster, adjust=2))
  b1d.crop <- writeRaster(b1d.crop, filename=paste(outdir, "/", paste(colnames(thisyr)[x],coi,yrs,"1deg",fun,sep="_"),".tif", sep=""), datatype=dt, format="GTiff", options="COMPRESS=LZW",overwrite=TRUE)
  base <- smooth.raster(base, ricemask_orig)	
  base <- crop(base,cellListExtent(poi, hraster, adjust=2))	
  base <- writeRaster(base, filename=paste(outdir, "/", paste(colnames(stress.dat)[x],coi,yrs,"smooth",fun,sep="_"),".tif", sep=""), format="GTiff", options="COMPRESS=LZW",overwrite=TRUE, datatype=dt)
  if(fun=="mean"){
    overall <- aggregate(x=stress.dat[,x], by=list(cell=stress.dat$cell), FUN="sd")
    
    base <- raster(hraster)
    base[overall$cell] <- overall[,2]
    b1d.crop <- crop(base,cellListExtent(poi,hraster, adjust=2))
    b1d.crop <- writeRaster(b1d.crop, filename=paste(outdir, "/", paste(colnames(thisyr)[x],coi,yrs,"1deg","sd",sep="_"),".tif", sep=""), datatype=dt, format="GTiff", options="COMPRESS=LZW",overwrite=TRUE)
    base <- smooth.raster(base, ricemask_orig)	
    base <- crop(base,cellListExtent(poi, hraster, adjust=2))	
    base <- writeRaster(base, filename=paste(outdir, "/", paste(colnames(stress.dat)[x],coi,yrs,"smooth","sd",sep="_"),".tif", sep=""), format="GTiff", options="COMPRESS=LZW",overwrite=TRUE, datatype=dt)
    
  } else {
    overall[,2] <- round(overall[,2]/length(years)*100,0)
    
    base <- raster(hraster)
    base[overall$cell] <- overall[,2]
    b1d.crop <- crop(base,cellListExtent(poi,hraster, adjust=2))
    b1d.crop <- writeRaster(b1d.crop, filename=paste(outdir, "/", paste(colnames(thisyr)[x],coi,yrs,"1deg","pct",sep="_"),".tif", sep=""), datatype=dt, format="GTiff", options="COMPRESS=LZW",overwrite=TRUE)
    base <- smooth.raster(base, ricemask_orig)	
    base <- crop(base,cellListExtent(poi, hraster, adjust=2))	
    base <- writeRaster(base, filename=paste(outdir, "/", paste(colnames(stress.dat)[x],coi,yrs,"smooth","pct",sep="_"),".tif", sep=""), format="GTiff", options="COMPRESS=LZW",overwrite=TRUE, datatype=dt)
    
  }
  
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


for (y in years){
	floweringdate <- dateFromDoy(praster[11109]+dae[11109]-30,y)
	if(is.na(floweringdate)) next		
	if(floweringdate-15<as.Date("1983-1-1")|floweringdate-15>as.Date("2014-12-31")) {
		show.message(poi[i], ": Required data exceeds coverage of database -", format(floweringdate,"%x"), appendLF=TRUE)
		next
	}
	stress.dat <- rbind(stress.dat, cbind(poi[i], y,ricestress.dayheat(wth=wth@w,flwrdate=floweringdate, method="Jagadish")))		
}