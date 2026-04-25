library(RiceStress)

# SETTINGS
low.maturity <- FALSE 
adm.level <- "region"
inpath <- paste(normalizePath(Sys.getenv("GDRIVE"),winslash="/"), "Projects/Heat Stress Models/Models/Inputs", sep = "/")
setwd(paste(normalizePath(Sys.getenv("GDRIVE"),winslash="/"), "Projects/Eclipse/neon2/RiceStress/data", sep = "/"))

#INPUT 
coi <- "Asia"
#coi <- "Pakistan"
years <- 1983:2016
yrs <- paste(min(years),max(years), sep="-")

wthvars <- list(nasacomp_1d=c("tmin"))

m <- raster(paste(inpath, "rice_1d_new2.tif", sep="/"))
ricemask_orig <- raster(paste(inpath, "rice_10min_new.tif", sep="/"))

con <- odbcConnect("climate_SRV3A")

for(h in 2:2){
	outdir <- paste( normalizePath(Sys.getenv("GDRIVE"),winslash="/") , "Projects/Heat Stress Models/Output/nighttime_2018005",coi,paste("Season", h, sep=""),sep="/")
	force.directories(outdir, recursive=TRUE)
	praster <- raster(paste("WORLD_PLANT_PK",h,"_60.tif", sep=""))
	praster <- mask(praster,m) 
	
	hraster <- raster(paste("WORLD_HARV_PK",h,"_60.tif", sep=""))
	hraster <-  mask(hraster, m)
	if (low.maturity){
		hraster <- calc(hraster, fun=function(x){return(ifelse(x > 0, x, NA) )})
		hraster <- overlay(praster, hraster, fun=function(x,y){return(ifelse(y > x, y - x, y - x + 365) )} )
		hraster  <- calc(hraster, fun=function(x){return(ifelse(x < 70, 70, x) )} )
		hraster  <- calc(hraster, fun=function(x){return(ifelse(x > 200, 200, x) )} )
		# floweringdate <- dateFromDoy(praster[poi[i]]+dae[poi[i]]-30,y)
		hraster <- praster+hraster
	}
	poi <- switch(adm.level, region=admAttribCells(filter=coi, admattrib="UNREGION2", reftable=system.file("data/WORLD_GADM_reference.csv", package="RiceStress"), refraster=system.file("data/WORLD_ID_0_60.tif", package="RiceStress")), country= admAttribCells(filter=coi, admattrib="NAME_ENGLI", reftable=system.file("data/WORLD_GADM_reference.csv", package="RiceStress"), refraster=system.file("data/WORLD_ID_0_60.tif", package="RiceStress")))
	poi <- poi[which((!is.na(hraster[poi])|hraster[poi]>0))]
	
	
	stress.dat <- vector()
	
	
	for ( i in 1:length(poi)){
		show.message("Analyzing pixel ", poi[i], "\r", EL=TRUE, appendLF=FALSE)
		if(is.na(praster[poi[i]])|praster[poi[i]]<=0) next
		xy <- xyFromCell(raster(),poi[i])
		wth <- fetch.daily(xy=xy, srcvars=wthvars, connection=con)
		#wth@w$tdew[wth@w$tdew<wth@w$tmin] <- wth@w$tmin[wth@w$tdew<wth@w$tmin]
		cell <- poi[i]
		for (y in years){
			#floweringdate <- dateFromDoy(hraster[poi[i]]-30,y)
			#if(is.na(floweringdate)) next		
			#if(floweringdate-15<as.Date("1983-1-1")|floweringdate-15>as.Date("2014-12-31")) {
			#	show.message(poi[i], ": Required data exceeds coverage of database -", format(floweringdate,"%x"), appendLF=TRUE)
			#	next
			#}
			stress.dat <- rbind(stress.dat, cbind(cell, y,t(ricestress.nighttimeHeat(wth, praster[poi[i]], hraster[poi[i]], year=y, consecutive=FALSE))))		
		}
	}
	
	colnames(stress.dat)[2] <- "year"
	stress.dat <- as.data.frame(stress.dat)
	write.csv(stress.dat, paste(outdir,"/stress-out_bypixel_byyear_season", h, ".csv",sep=""))
	
# Mapping stress data
	
	for (y in years){		
		thisyr <- stress.dat[stress.dat$year==y,]
		for (x in 3:ncol(thisyr)){
			dt <- ifelse(grepl("cnt",colnames(stress.dat)[x], ignore.case=TRUE)|grepl("stress",colnames(stress.dat)[x], ignore.case=TRUE),"INT4S","FLT4S")
			base <- raster(hraster)
			base[thisyr$cell] <- thisyr[,x]
			b1d.crop <- crop(base,cellListExtent(poi,hraster, adjust=2))
			b1d.crop <- writeRaster(b1d.crop, filename=paste(outdir, "/", paste(colnames(thisyr)[x],coi,y,"1deg",sep="_"),".tif", sep=""), datatype=dt, format="GTiff", options="COMPRESS=LZW",overwrite=TRUE)
			base <- smooth.raster(base, ricemask_orig)
			base <- crop(base,cellListExtent(poi,hraster, adjust=2))
			base <- writeRaster(base, filename=paste(outdir, "/", paste(colnames(thisyr)[x],coi,y,"smooth",sep="_"),".tif", sep=""), datatype=dt, format="GTiff", options="COMPRESS=LZW",overwrite=TRUE)
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

con <- odbcClose(con)


#stress.tmax <- raster("Asia_1983-2012_Stress.Based.on.TMax.tif")
#plot(stress.tmax, col=brewer.pal(9,"Reds"),zlim=c(1,30))
#plot(world, add=TRUE, border="light grey")
#
#stress.tpan <- raster("Asia_1983-2012_Stress.Based.on.MaxTPan.tif")
#plot(stress.tpan, col=brewer.pal(9,"Reds"),zlim=c(1,30))
#plot(world, add=TRUE, border="light grey")
#stress.sf <- raster("Asia_1983-2012_Stress.Based.on.Spikelet.Fertility.tif")
#plot(stress.sf, col=brewer.pal(9,"Reds"),zlim=c(1,30))
#plot(world, add=TRUE, border="light grey")
#
#
#stress.map <- function(con, plantraster, harvraster, model=ricestress.dayheat,...){
#	con <- odbcConnect("climate_SRV3A")
#} 
#
#
#for (y in years){
#	floweringdate <- dateFromDoy(praster[11109]+dae[11109]-30,y)
#	if(is.na(floweringdate)) next		
#	if(floweringdate-15<as.Date("1983-1-1")|floweringdate-15>as.Date("2014-12-31")) {
#		show.message(poi[i], ": Required data exceeds coverage of database -", format(floweringdate,"%x"), appendLF=TRUE)
#		next
#	}
#	stress.dat <- rbind(stress.dat, cbind(poi[i], y,ricestress.dayheat(wth=wth@w,flwrdate=floweringdate, method="Jagadish")))		
#}