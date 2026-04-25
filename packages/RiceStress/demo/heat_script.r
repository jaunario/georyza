
# SETTINGS
low.maturity =  FALSE 
adm.level = "country"
inpath <- paste(normalizePath(Sys.getenv("GDRIVE"),winslash="/"), "Projects/Heat Stress Models/Models/Inputs", sep = "/")
setwd(paste(normalizePath(Sys.getenv("GDRIVE"),winslash="/"), "Projects/Eclipse/neon2/RiceStress/data", sep = "/"))

library(RiceStress)

#INPUT 

coi <- "Asia"
#coi <- "Pakistan"
years <- 1981:2018

yrs <- paste(min(years),max(years), sep="-")

m <- raster(paste(inpath, "rice_1d_new2.tif", sep="/"))
praster <- raster("WORLD_PLANT_PK1_60.tif")
praster <- mask(praster,m) 

hraster <- raster("WORLD_HARV_PK1_60.tif")
hraster <-  mask(hraster, m)
if (low.maturity){
	hraster <- calc(hraster, fun=function(x){return(ifelse(x > 0, x, NA) )})
	hraster <- overlay(praster, hraster, fun=function(x,y){return(ifelse(y > x, y - x, y - x + 365) )} )
	hraster  <- calc(hraster, fun=function(x){return(ifelse(x < 70, 70, x) )} )
	hraster  <- calc(hraster, fun=function(x){return(ifelse(x > 200, 200, x) )} )
	# floweringdate <- dateFromDoy(praster[poi[i]]+dae[poi[i]]-30,y)
	hraster <- praster+hraster
}
poi <- switch(adm.level, region=admAttribCells(filter=coi, admattrib="UNREGION2"), country= admAttribCells(filter=coi, admattrib="NAME_ENGLI"))
poi <- poi[which((!is.na(hraster[poi])|hraster[poi]>0))]


outdir <- paste( normalizePath(Sys.getenv("GDRIVE"),winslash="/") , "Projects/Heat Stress Models/Output/uncor.35.direct.sol",coi,sep="/")
force.directories(outdir, recursive=TRUE)

wthvars <- list(nasa_correction=c("tmax", "tmin", "tdew"))
con <- odbcConnect("climate_SRV3A")
stress.dat <- vector()
for ( i in 1:length(poi)){
	show.message("Analyzing pixel ", poi[i], "\r", EL=TRUE, appendLF=FALSE)
	if(is.na(praster[poi[i]])|praster[poi[i]]<=0) next
	xy <- xyFromCell(raster(),poi[i])
	wth <- fetch.daily(xy=xy, srcvars=wthvars, connection=con)[[1]]
	wth@w$tdew[wth@w$tdew<wth@w$tmin] <- wth@w$tmin[wth@w$tdew<wth@w$tmin]
	
	for (y in years){
		floweringdate <- dateFromDoy(hraster[poi[i]]-30,y)
		if(is.na(floweringdate)) next		
		if(floweringdate-15<as.Date("1983-1-1")|floweringdate-15>as.Date("2014-12-31")) {
			show.message(poi[i], ": Required data exceeds coverage of database -", format(floweringdate,"%x"), appendLF=TRUE)
			next
		}
		stress.dat <- rbind(stress.dat, cbind(poi[i], y,ricestress.dayheat(wth=wth,flwrdate=floweringdate, method="Jagadish", tmax.threshold=35)))		
	}
}
con <- odbcClose(con)
colnames(stress.dat)[1:2] <- c("cell", "year")
stress.dat <- as.data.frame(stress.dat)
write.csv(stress.dat, paste(outdir,"out.csv",sep="/"))
ricemask_orig <- raster(paste(inpath, "rice_10min_new.tif", sep="/"))

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


for (x in 6:ncol(stress.dat)){	
	fun <- ifelse(grepl("stress",colnames(stress.dat)[x], ignore.case=TRUE), "sum", "mean")
	overall <- aggregate(x=stress.dat[,x], by=list(cell=stress.dat$cell), FUN=fun)
	if(sum(!is.na(overall[,2]))==0|sum(overall[,2],na.rm=TRUE)==0) next
	base <- raster(hraster)
	base[overall$cell] <- overall[,2]
	base <- smooth.raster(base, ricemask_orig)	
	base <- crop(base,cellListExtent(poi, hraster, adjust=2))	
	base <- writeRaster(base, filename=paste(outdir, "/", paste(colnames(stress.dat)[x],coi,yrs,sep="_"),".tif", sep=""), format="GTiff", options="COMPRESS=LZW",overwrite=TRUE, datatype="FLT4S")
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