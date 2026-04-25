source("D:/Google Drive/Projects/kepler_workspace/RiceStress/R/cold.r")

# INPUTS
target.res <- 0.25
#countries <- "WORLD"
#countries <- c("Ethiopia", "Madagascar","Burundi", "Rwanda")
countries <- "Indonesia"

# INPUT SETTINGS
inpath <- "D:/Google Drive/Projects/Heat Stress Models/Models/Inputs"
ctryref <- read.csv("D:/Google Drive/Shared Data/GADM/WORLD_GADM_reference.csv")
m <- raster(paste(inpath, "rice_1deg_new2.tif", sep="/"))
projection(m) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

ricemask_orig <- raster(paste(inpath, "rice_10min_new.tif", sep="/"))

years <- 1983:2012
climdsn <- "climate_SRV3A"
climset <- "nasa_correction"
type <- "cold"
filled <- FALSE

#fd1  <- list.files(path=inpath,pattern="PLANT_PK[1-3]_60.tif$", full.names=TRUE)
#fde1  <- list.files(path=inpath,pattern="HARV_PK[1-3]_60.tif$", full.names=TRUE)

#fd1  <- list.files(path=inpath,pattern="PLANT_ST[1-3]_60.tif$", full.names=TRUE)
#fde1  <- list.files(path=inpath,pattern="HARV_ST[1-3]_60.tif$", full.names=TRUE)

#fd1  <- list.files(path=inpath,pattern="PLANT_END[1-3]_60.tif$", full.names=TRUE)
#fde1  <- list.files(path=inpath,pattern="HARV_END[1-3]_60.tif$", full.names=TRUE)

for (i in 1:length(countries)){
  country <- countries[i]
  if(filled){
    ricecal.dir <- paste("E:/Rice Calendar", country, "filled", target.res*60, "latest",sep="/")
  } else {
    ricecal.dir <- paste("E:/Rice Calendar", country, target.res*60, "latest",sep="/")  
  }
  
  plantfiles <- dir(ricecal.dir, pattern=paste("PLANT.*.",target.res*60,".tif$",sep=""), full.names=TRUE)
  harvfiles  <- dir(ricecal.dir, pattern=paste("HARV.*.",target.res*60,".tif$",sep=""), full.names=TRUE)
  filemeta <- matrix(unlist(strsplit(sub(".tif","",basename(plantfiles)),"_")),ncol=4, byrow=TRUE)


  #adm <- raster(dir(path="D:/Google Drive/Shared Data/GADM",pattern=paste(country,".*.",target.res*60,".tif$",sep=""),full.names=TRUE))
  
  #ctrymask <- adm
  #ctrymask[] <- !is.na(ctrymask[])
  #ctrymask[ctrymask[]!=1] <-NA
  #ctryxy <- xyFromCell(ctrymask,which(!is.na(ctrymask[])))
  #ctryext <-extent(min(ctryxy[,"x"]),max(ctryxy[,"x"]),min(ctryxy[,"y"]),max(ctryxy[,"y"]))
  
  for (j in 1:length(plantfiles)){
    planting <- raster(plantfiles[j])
    #planting <- mask(planting, ctrymask)
    #planting <- mask(planting, m)
    #planting <- crop(planting,ctryext)
    
    if (length(which(!is.na(planting[])))==0) {
        show.message(country,": Not planting for season ",filemeta[j,3],appendLF=TRUE)
        next
    }
    planting  <- calc(planting, fun=function(x){return(ifelse(x > 0, x, NA) )} )
    
    dend1 <-  raster(harvfiles[j])
    #dend1 <- mask(dend1, ctrymask)
    #dend1 <- mask(dend1, m)
    #dend1 <- crop(dend1,ctryext)
    dend1 <-calc(dend1, fun=function(x){return(ifelse(x > 0, x, NA) )} )
    dae  <- overlay(planting, dend1, fun=function(x,y){return(ifelse(y > x, y - x, y - x + 365) )} )
    dae  <- calc(dae, fun=function(x){return(ifelse(x < 70, 70, x) )} )                    #adjust too low maturity
    dae  <- calc(dae, fun=function(x){return(ifelse(x > 200, 200, x) )} )                    #adjust too low maturity

    #ricemask <- crop(ricemask_orig, planting)
    #projection(ricemask) <- projection(planting)
    
    procstart <- Sys.time() 
    outpath1 <- paste("D:/Google Drive/Projects/Heat Stress Models/Output", country, type, filemeta[j,3],sep="/")
    #outpath1 <- "c:/projects/jrb/cold/peak/crop1"
    #thermal.stress(planting, dae, years, type, climdsn,climset, maskraster=ricemask, outdir=outpath1)
    thermal.stress(planting, dae, years, type, climdsn,climset, outdir=outpath1) #NO RICE MASK
    procend <- Sys.time()    
    stresslog <- c(paste("Process Start:", procstart),paste("Process End:",procend),paste("Process time:", round(difftime(procend,procstart, unit="mins"),2),"mins.\n"))
    writeLines(stresslog, paste(outpath1,"cold_log.txt",sep="/"))

  }
    
}

#m <- crop(m, planting)
#years <- 1983:2008


planting <- raster(fd1[2])
NAvalue(planting) <- 65535
projection(planting) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
planting <- mask(planting, m)

dend1 <-  raster(fde1[2])
NAvalue(dend1) <- 65535
projection(dend1) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
dend1 <- mask(dend1, m)

dend1 <-calc(dend1, fun=function(x){return(ifelse(x > 0, x, NA) )} )
dae<- overlay(planting, dend1, fun=function(x,y){return(ifelse(y > x, y - x, y - x + 365) )} )
dae  <- calc(dae, fun=function(x){return(ifelse(x < 70, 70, x) )} )                    #adjust too low maturity
dae  <- calc(dae, fun=function(x){return(ifelse(x > 200, 200, x) )} )                    #adjust too low maturity

planting  <- calc(planting, fun=function(x){return(ifelse(x > 0, x, NA) )} )

procstart <- Sys.time() 
type <- "cold"
#outpath1 <- paste("E:/stress","test","cold",sep="/")
outpath1 <- "D:/Google Drive/Projects/Heat Stress Models/END2"
#outpath1 <- "c:/projects/jrb/cold/peak/crop1"
thermal.stress(planting, dae, years, type, climdsn,climset, maskraster=ricemask, outdir=outpath1)
procend <- Sys.time()    
stresslog <- c(paste("Process Start:", procstart),paste("Process End:",procend),paste("Process time:", round(difftime(procend,procstart, unit="mins"),2),"mins.\n"))
writeLines(stresslog, paste(outpath1,"cold_log.txt",sep="/"))


planting <- raster(fd1[3])
NAvalue(planting) <- 65535
projection(planting) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
planting <- mask(planting, m)

dend1 <-  raster(fde1[3])
NAvalue(dend1) <- 65535
projection(dend1) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
dend1 <- mask(dend1, m)

dend1 <-calc(dend1, fun=function(x){return(ifelse(x > 0, x, NA) )} )
dae<- overlay(planting, dend1, fun=function(x,y){return(ifelse(y > x, y - x, y - x + 365) )} )
dae  <- calc(dae, fun=function(x){return(ifelse(x < 70, 70, x) )} )                    #adjust too low maturity
dae  <- calc(dae, fun=function(x){return(ifelse(x > 200, 200, x) )} )                    #adjust too low maturity

planting  <- calc(planting, fun=function(x){return(ifelse(x > 0, x, NA) )} )

procstart <- Sys.time() 
type <- "cold"
#outpath1 <- paste("E:/stress","test","cold",sep="/")
outpath1 <- "D:/Google Drive/Projects/Heat Stress Models/END3"
#outpath1 <- "c:/projects/jrb/cold/peak/crop1"
thermal.stress(planting, dae, years, type, climdsn,climset, maskraster=ricemask, outdir=outpath1)
procend <- Sys.time()    
stresslog <- c(paste("Process Start:", procstart),paste("Process End:",procend),paste("Process time:", round(difftime(procend,procstart, unit="mins"),2),"mins.\n"))
writeLines(stresslog, paste(outpath1,"cold_log.txt",sep="/"))

