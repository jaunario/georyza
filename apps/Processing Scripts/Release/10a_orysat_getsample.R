# DIRECTORY SETTINGS
WORKSPACE   = "E:/WORKSPACE/Demeter Initiative/Analysis"
MODIS_DIR  = "E:/WORKSPACE/Demeter Initiative/500m/Indices"
MASK_DIR = "./OUTPUT/mask_highres/"
CACHE_DIR    = paste0(WORKSPACE, "/cache")
OUTPUT_DIR  = "AUTO"

VERSION    = 0.1 
PRODUCTS    = list(m500="MOD09A1", m250="MXD13Q1")
# CONSOLIDATE MONTHLY RICE AREAS ON EACH TILE

# TODO: Have default layer-parameter mapping but allow custom mapping for experimental methods
TARGET_RES=500
TMPLATE_PROD <- PRODUCTS[paste0("m", TARGET_RES)][[1]]

#PROJ.SINU = "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext"
PROJ.SINU = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
PROJ.CEA  = "+proj=cea +lon_0=100 +lat_ts=0 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

binIndextoString <- function(x, length=100){
  bins <- rep(0,length)
  bins[x] <- 1
  bins <- rev(bins)
  return(paste(bins, collapse = ""))
}

binStringSum <- function(x){
  result <- strsplit(x,"")[[1]]
  result <- as.integer(result)
  return(sum(result))
}

countUnique <- function(x){
  x <- na.omit(x)
  if(length(x)==0) result <- NA else result <- length(unique(x))
  return(result)
}

library(orysat)
library(terra)
#library(parallel)


block.apply <- function(x, baseraster=NULL, blocks = 25, fun = NULL, extend.to.base=FALSE, retain.blanks=TRUE, ...){
  # TODO: Find the most efficient block size
  if(is.null(baseraster)) rst.base <- x else rst.base <- baseraster
  if (is.na(blocks) || blocks=="auto"){
    blocks <- 25   
  }
  if(length(blocks)==1){
    blk.ncol <- ceiling(ncol(rst.base)/sqrt(blocks))
    blk.nrow <- ceiling(nrow(rst.base)/sqrt(blocks))
  
    blk.colst <- seq(1, ncol(rst.base), by=blk.ncol)
    blk.colen <- pmin(blk.colst+blk.ncol - 1,ncol(rst.base))
    blk.rowst <- seq(1, nrow(rst.base), by=blk.nrow)
    blk.rowen <- pmin(blk.rowst+blk.nrow - 1,nrow(rst.base))
    
    blk.xmin <- xFromCol(rst.base, blk.colst)-xres(rst.base)/2
    blk.xmax <- xFromCol(rst.base, blk.colen)+xres(rst.base)/2
    
    blk.ymax <- yFromRow(rst.base, blk.rowst)+yres(rst.base)/2
    blk.ymin <- yFromRow(rst.base, blk.rowen)-yres(rst.base)/2
    
    blk.idx <- expand.grid(1:sqrt(blocks), 1:sqrt(blocks))
    
    blk.coords <- data.frame(xmin=blk.xmin[blk.idx[,2]], xmax=blk.xmax[blk.idx[,2]], ymin=blk.ymin[blk.idx[,1]], ymax=blk.ymax[blk.idx[,1]])
    
    # Check if blocks intersects/overlaps with x
    blk.extents <- apply(t(blk.coords),2, extent)
    x.extent <- extent(as.vector(ext(x)))
    blk.intersect <- sapply(blk.extents, raster::intersect, y=x.extent)
    blk.extents <- apply(blk.coords, 1, ext)
  }
  # TODO: Parallel block operation
  #cores = 2
  #registerDoParallel(2)
  #minmax(x)[1]:minmax(x)[2]
  
  xx <- list()
  for(i in 1:blocks) {
    message("blk ", i)
    rst.blkbase <- crop(rst.base, blk.extents[[i]])
    
    if(!is.null(blk.intersect[[i]])) {
      aa <- crop(x, rst.blkbase)
      if((sum(is.na(minmax(aa)))<2 && minmax(aa)[1]<minmax(aa)[2]) && sum(minmax(aa)>0)>0) {
        aa <- resample(aa, rst.blkbase, method="ngb")
        if(!is.null(fun)) aa <- fun(aa, ...)
      } else aa <- rst.blkbase 
      
      #subset.wrapper(rst = x, subset.ext = blk.extents[[i]], baseraster = rst.base,  ...)  
    } else if(retain.blanks){
      aa <- rst.blkbase
    } else aa <- NULL
    
    if(!is.null(aa)) {
      xx <- append(xx, list(aa))
    }
  }

  if(length(xx)>1) {
    xx <- do.call(merge, xx) 
  } else if (length(xx)==1) {
      xx <- xx[[1]]
  } else xx <- NULL
  if(!is.null(xx) & extend.to.base) xx <- extend(xx, rst.base)
  return(xx)
} 


setwd(WORKSPACE)

if(OUTPUT_DIR=="AUTO") {
  OUTPUT_DIR = paste0("./OUTPUT/",paste0("m", TARGET_RES, "_sample"))
  if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
}

filenames.raster <- dir("./OUTPUT/mask_highres", recursive=TRUE, full.names = TRUE)
fname.attr <- sapply(basename(filenames.raster), strsplit, split="_", USE.NAMES = FALSE)
info.attr <- do.call(rbind, fname.attr)  
info.attr <- data.frame(filename=filenames.raster, layername=info.attr[,2], resolution=info.attr[,1], tile=sub(".tif", "", info.attr[,3]), stringsAsFactors = FALSE)

# years <- unique(info.fnames$year)
# years <- 2018:2020
# seasons <- unique(info.fnames$season)
if(!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR)
#pylibs.start()
#pylibs.start(python="~/../miniconda3")

# shp.modis <- shapefile(MODIS_GRIDFILE)
# projection(shp.modis) <- PROJ.MODIS 
#i=j=k=l=1

done <- vector()

# COMPUTE SAMPLE SIZE
assumed.u <- 100  # Assumed p
effect.size <- 0.02 # 0.02=small, 0.15=medium, and 0.35 large effect sizes  
siglevel <- 0.05
power.model <- 0.95

power.calc <- pwr::pwr.f2.test(u=assumed.u, f2=effect.size, sig.level = siglevel, power = power.model)
samplesize <- round(round(power.calc$v)/10^(nchar(round(power.calc$v))-1))*10^(nchar(round(power.calc$v))-1)
sample.validation <- samplesize*0.25

# Creating Monthly Rice Area Maps
samplesize.tile <- data.frame(tile=info.attr$tile, npop=0, nsamples=0, stringsAsFactors = FALSE)
cells.pop <- list()
for (i in 1:length(info.attr$tile)){
  tmplate_file <- dir(paste0(MODIS_DIR, "/", info.attr$tile[i]), pattern = paste(TMPLATE_PROD, info.attr$tile[i],sep="."), full.names = TRUE)[1]
  rst.base <- rast(tmplate_file)
  rst.base <- rast(rst.base)
  
  rst.mask <- rast(info.attr$filename[i])
  agg.fact <- round(xres(rst.base)/xres(rst.mask))
  rst.mask <- aggregate(rst.mask, fact=agg.fact, fun=sum)
  rst.mask <- rst.mask == (agg.fact^2)
  vals.mask <- values(rst.mask)
  cells.pop[[info.attr$tile[i]]] <- which(vals.mask==1)
  rm(rst.base, rst.mask, vals.mask)
  gc(reset = TRUE)
  tmpFiles(remove = TRUE)

}
samplesize <- ceiling(samplesize+samplesize*0.25)
samplesize.tile$npop <- sapply(cells.pop, length)
samplesize.tile$nsamples <- ceiling(samplesize*samplesize.tile$npop/sum(samplesize.tile$npop))

  
cells.samp <- mapply(sample, cells.pop, as.list(samplesize.tile$nsamples))
sample.validation <- mapply(sample, cells.samp, as.list(ceiling(samplesize.tile$nsamples*0.2)))
sample.training <- mapply("[", cells.samp, lapply(mapply(match, cells.samp, sample.validation),is.na))
gc(reset=TRUE)

# Extract Rice Area percentage (RAP) by month and year
files.ra <- dir("./OUTPUT/By_tile-year-month", pattern=".tif$", full.names = TRUE)
info.filesra <- sapply(basename(files.ra), strsplit, "_")
info.filesra <- do.call(rbind, info.filesra)
info.filesra[,5] <- sub(".tif", "", info.filesra[,5])
info.filesra <- data.frame(filename=files.ra, tile=info.filesra[,3], year=as.integer(info.filesra[,4]), month=as.integer(info.filesra[,5]))

for ( i in 1:length(info.attr$tile)){
  tile <- info.attr$tile[i]
  this.tile <- info.filesra[info.filesra$tile==tile,]
  dat.train <- dat.validate <- vector()
  for(j in 1:nrow(this.tile)){
    rst.ra <- rast(this.tile$filename[j])
    rst.ra <- aggregate(rst.ra, fact=agg.fact, fun=sum)
    rst.ra <- rst.ra/(agg.fact^2)
    dat.train <- cbind(dat.train, rst.ra[sample.training[[i]]]) 
    dat.validate <- cbind(dat.validate, rst.ra[sample.validation[[i]]]) 
    rm(rst.ra)
    gc(reset=TRUE)
    tmpFiles(remove = TRUE)
  }
  colnames(dat.train) <- colnames(dat.validate) <- sprintf("Y%04d_M%02d", this.tile$year,this.tile$month)
  dat.train <- data.frame(cell=sample.training[[i]], dat.train)
  dat.validate <- data.frame(cell=sample.validation[[i]], dat.validate)
  write.csv(dat.train, file=paste0("./sample/training/", paste("RAP", tile, TARGET_RES, sep="_"), ".csv"), row.names = FALSE)
  write.csv(dat.validate, file=paste0("./sample/validation/", paste("RAP", tile, TARGET_RES, sep="_"), ".csv"), row.names = FALSE)
  rm(dat.train, dat.validate)
  gc(reset=TRUE)
}

# Extract sample Training Data
INDICES = c("evi","lswi", "mndwi", "ndfi", "ndvi")
YEARS = 2017:2020
for ( i in 1:length(info.attr$tile)){
  tile <- info.attr$tile[i]
  #this.tile <- info.filesra[info.filesra$tile==tile,]
  dat.train <- read.csv(file=paste0("./sample/training/", paste("RAP", tile, TARGET_RES, sep="_"), ".csv"))
  for(j in 1:length(INDICES)){
    index <- INDICES[j]
    message(index,"-", tile, ": Making inventory.")
    files.index <- dir(paste0(MODIS_DIR,"/", tile), pattern=paste0(index,".tif$"), ignore.case = TRUE, full.names = TRUE)
    info.fromfname <- sapply(sub(".tif", "", basename(files.index)), strsplit, split="\\.")
    info.fromfname <- do.call(rbind, info.fromfname)
    inv.index <- data.frame(filename=files.index, acqdate=info.fromfname[,3], year=substr(info.fromfname[,3],2,5))
    inv.index <- inv.index[inv.index$year %in% YEARS,]
    
    message(index,"-", tile, ": loading rasters")
    stk.index <- rast(inv.index$filename)
    
    message(index, "-", tile, ": extracting values of sample ", tile)
    dat.index <- stk.index[dat.train$cell]
    dat.index <- cbind(dat.train$cell, dat.index)
    colnames(dat.index) <- c("cell", inv.index$acqdate)
    message(index,"-", tile, ": Saving index data to csv.")
    write.csv(dat.index, paste0(OUTPUT_DIR, "/", paste(tile,index,min(YEARS),max(YEARS), sep="_"),".csv"), row.names = FALSE)
    message(index,"-", tile, ": cleaning up.")
    rm(files.index, info.fromfname, inv.index, stk.index, dat.index)
    gc(reset=TRUE)
    tmpFiles(remove=TRUE)
  }
}

       for (l in 1:length(tiles.aoi)){
        if(!exists("rst.base")){
          if (!is.na(tmplate_file)) {
            rst.base <- rast(tmplate_file)
            rst.base <- rast(rst.base)
            rst.base <- disaggregate(rst.base, fact=20)
          } else {
            stop("Template file = ", tmplate_file)
          }
        } 
        
        #rst.reprojtile <- crop(rst.reproj, rst.base)
        message(basename(this.season$filename[k]), "-", tiles.aoi[l], ": Cropping to MODIS tile.")
        rst.reprojtile <- block.apply(rst.reproj, rst.base, extend.to.base = TRUE)
        names(rst.reprojtile) <- paste(tiles.aoi[l],"SAR-based","SOS-Month", years[i], sep = "_")
        if(sum(minmax(rst.reprojtile)==0, na.rm=TRUE)>1) {
          rm(rst.reprojtile, rst.base)
          gc(reset=TRUE)
          next
        }
        
        if(i==1 & j==1){
          message(basename(this.season$filename[k]), ": Mask for ", tiles.aoi[l])
          mask_file <- paste0("./OUTPUT/mask_highres/", paste("MASK", paste0(ceiling(xres(rst.reproj)),"m"), tiles.aoi[l], sep="_") ,".tif")
          rst.masktile <- !is.na(rst.reprojtile)
          
          if(!minmax(rst.masktile)[2]>0) {
            rm(rst.masktile)
            gc(reset=TRUE)
            next
          }
          
          message(basename(this.season$filename[k]), ": Converting NAs to 0.")
          rst.masktile <- classify(rst.masktile, rcl=t(c(NA,0)))
          names(rst.masktile) <- paste(tiles.aoi[l],"SAR-based","MASK", sep = "_")
          
          if(file.exists(mask_file)){
            message(basename(this.season$filename[k]), "-", tiles.aoi[l], ": Combining mask from other files.")
            rst.old <- rast(mask_file)
            rst.masktile <- sum(rst.masktile, rst.old, na.rm=TRUE)
            if(minmax(rst.masktile)[2]>1) rst.masktile <- rst.masktile>0
            rm(rst.old)
            #rst.monthtile <- extend(rst.monthtile, rst.old)
            #rst.old <- extend(rst.old, rst.monthtile)
            #rst.old <- classify(rst.old, t(c(NA,0)))
          } 
          message(basename(this.season$filename[k]), ": Saving mask to disk. (", basename(mask_file), ")")
          rst.masktile <- terra::writeRaster(rst.masktile, mask_file, overwrite=TRUE, wopt=list(datatype="INT2S", gdal="COMPRESS=LZW"))
          rm(rst.masktile)
          gc(reset=TRUE)
        }
        
        rst.reprojtile <- classify(rst.reprojtile, rcl=t(c(0,NA)))
        
        
        for(m in 2:12){
       
          if(sum(is.na(minmax(rst.reprojtile)))==2 || (m < minmax(rst.reprojtile)[1] | m > minmax(rst.reprojtile)[2])) next
          
          output_file <- paste0("./OUTPUT/By_tile-year-month/", paste("RiceSAR", paste0(ceiling(xres(rst.reprojtile)),"m"), tiles.aoi[l], years[i], sprintf("%02d", m),sep="_") ,".tif")
          
          message(basename(this.season$filename[k]), "-", tiles.aoi[l], ": Isolating rice with SOS falling on ", month.abb[m], ".")
          rst.monthtile <- rst.reprojtile==m #block.apply(rst.reproj, rst.base, extend.to.base = FALSE, fun=isolate.value, value=m, rst.old=rst.old)
      
          if(!minmax(rst.monthtile)[2]>0) {
            rm(rst.monthtile)
            gc(reset=TRUE)
            next
          }
          
          rst.monthtile <- classify(rst.monthtile, t(c(NA,0)))
          if(file.exists(output_file)){
            message(basename(this.season$filename[k]), "-", tiles.aoi[l], ": Combining data from other files.")
            rst.old <- rast(output_file)
            rst.monthtile <- sum(rst.monthtile, rst.old, na.rm=TRUE)
            if(minmax(rst.monthtile)[2]>1) rst.monthtile <- rst.monthtile>0
            rm(rst.old)
            #rst.monthtile <- extend(rst.monthtile, rst.old)
            #rst.old <- extend(rst.old, rst.monthtile)
            #rst.old <- classify(rst.old, t(c(NA,0)))
          } 
    
          #rst.monthtile <- block.apply(rst.monthtile, rst.base, extend.to.base = FALSE, fun=isolate.value, value=1, rst.old=rst.old)
          
          message(basename(this.season$filename[k]), "-", tiles.aoi[l], ": Writing to disk")
          rst.monthtile <- terra::writeRaster(rst.monthtile, output_file, overwrite=TRUE, wopt=list(datatype="INT2S", gdal="COMPRESS=LZW"))
          
          rm(rst.monthtile)
          gc(reset=TRUE)
        }
        message(basename(this.season$filename[k]), ": Cleaning up.")
        rm(rst.base, rst.reprojtile)
        gc(reset=TRUE)
      }
      rm(rst.reproj)
      tmpFiles(remove = TRUE)
      done <- c(done, this.season$filename[k])
    } 
  }
}



setwd(WORKSPACE)

if(OUTPUT_DIR=="AUTO") {
  OUTPUT_DIR = paste0("./OUTPUT/",paste0("m", TARGET_RES))
  if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
}

filenames.raster <- dir("./OUTPUT/By_tile-year-month", pattern=".tif$", recursive=TRUE, full.names = TRUE)
fname.attr <- sapply(basename(filenames.raster), strsplit, split="_", USE.NAMES = FALSE)
fname.attr <- do.call(rbind, fname.attr)
colnames(fname.attr) <- c("src", "res", "tile", "year", "month")
fname.attr <- data.frame(filename=filenames.raster, fname.attr, stringsAsFactors = FALSE)
fname.attr$year <- as.integer(fname.attr$year)
fname.attr$month <- as.integer(sub(".tif", "", fname.attr$month))

tiles <- unique(fname.attr$tile)
years <- unique(fname.attr$year)
for (i in 1:length(tiles)){
  tmplate_file <- dir(paste0(MODIS_DIR, "/", tiles[i]), pattern = paste(TMPLATE_PROD, tiles[i], sep="."), full.names = TRUE)[1]

  for(y in  years){
    stk.yr <- rast(fname.attr$filename[fname.attr$year==y & fname.attr$tile==tiles[i]])
    rast.yr <- sum(stk.yr, na.rm=TRUE)
  }
  #shp.tile <- shp.modis[shp.modis$tile==tiles.aoi[l],]
  # Create/update mask
  maskoutput_file <- paste0(OUTPUT_DIR,"/../mask/", paste("MASK", TARGET_RES, tiles.aoi[l], years[i], sep="_") ,".tif")
  if(!file.exists(maskoutput_file)){
    tmplate_file <- dir(paste0(MODIS_DIR, "/", tiles.aoi[l]), pattern = paste(TMPLATE_PROD, tiles.aoi[l],paste0("A",years[i]),sep="."), full.names = TRUE)[1]
    if(!exists("rst.base") && !is.na(tmplate_file)) {
      rst.base <- rast(tmplate_file)
    } else {
      stop("Template file = ", tmplate_file)
    }
  } else {
    rst.base <- rast(maskoutput_file)
  }
  #pct.area <- xres(rst.masked)^2/xres(rst.base)^2
  rst.masktile <- crop(rst.mask, rst.base)
  
  #rst.monthtile <- rst.maskedtile==m
  message(basename(this.season$filename[k]), ": Getting exploration area covered by tile ", tiles.aoi[l])
  rst.masktile <- terra::aggregate(rst.masktile,fact=round(xres(rst.base)/xres(rst.masktile),0), fun=sum, na.rm=TRUE)
  rst.masktile <- terra::resample(rst.masktile, rst.base, method="ngb")
  rst.masktile <- classify(rst.masktile, rcl=cbind(1, Inf, 1), othersNA=TRUE)
  rst.masktile <- !is.na(rst.masktile)
  
  if(file.exists(maskoutput_file)) {
    message(basename(this.season$filename[k]), ": Combining other areas in the mask")
    rst.masktile <- rst.masktile+rst.base
    rst.masktile <- rst.masktile>0
  }
  
  rst.masktile <- terra::writeRaster(rst.masktile, maskoutput_file, overwrite=TRUE)
  rst.reprojtile <- crop(rst.reproj, rst.base)
  
  rm(rst.base)
  
  pct.area <- xres(rst.reprojtile)^2/xres(rst.masktile)^2
  
  for(m in 1:12){
    # If month not in rst.masked, skip 
    if(m < minmax(rst.reproj)[1] | m > minmax(rst.masked)[2]) next
    output_file <- paste0(OUTPUT_DIR, "/", paste("RiceArea", "PCT", TARGET_RES, tiles.aoi[l], years[i], sprintf("m%02d",m),sep="_") ,".tif")
    
    if(!file.exists(output_file)){
      rst.old <- rast(rst.masktile)
      terra::values(rst.old) <- 0
    } else {
      rst.old <- rast(output_file)
    }
    
    message(basename(this.season$filename[k]), ": Isolating rice with SOS falling on ", month.abb[m], " in ", tiles.aoi[l])
    rst.monthtile <- rst.reprojtile==m
    
    message(basename(this.season$filename[k]), ": Calculating PCT AREA ", month.abb[m], " in ", tiles.aoi[l])
    
    rst.pct <- terra::aggregate(rst.monthtile,fact=round(xres(rst.old)/xres(rst.monthtile),0), fun=sum, na.rm=TRUE)
    rst.pct <- terra::resample(rst.pct, rst.old, method="ngb")
    rst.pct <- rst.pct*pct.area
    rst.pct <- classify(rst.pct, rcl=matrix(c(NA,0), ncol=2)) # Fill NA with 0
    
    #rst.pct <- terra::aggregate(rst.pct, fact=10, fun=sum, na.rm=TRUE)
    if(file.exists(output_file)) {
      rst.pct <- rst.pct+rst.old
    }
    if(minmax(rst.pct)[2]>1.5) warning("Pixels exceeding 150% rice. ", basename(fname.reproj))
    rst.pct <- classify(rst.pct, rcl=matrix(c(1,Inf,1), ncol=3)) # Fill NA with 0
    message(basename(this.season$filename[k]), ": Writing to disk ", basename(output_file), ".")
    names(rst.pct) <- month.abb[m]
    rst.pct <- terra::writeRaster(rst.pct, output_file, overwrite=TRUE)
    
    # if(1 %in% minmax(rst.monthtile)){
    #   rst.maskedtile <- mask(rst.maskedtile, rst.monthtile)
    # } else {
    #   next
    # }
    rm(rst.old, rst.pct, rst.monthtile)
  }
  
  #rst.modal <- terra::aggregate(rst.masked,fact=round(xres(rst.disaggcrop)/xres(rst.reproj),0), fun=modal, na.rm=TRUE)
  
  # rst.count <- terra::aggregate(rst.masked,fact=round(xres(rst.base)/xres(rst.monthtile),0), fun=countUnique)
  # rst.count <- terra::zonal(rst.maskedtile, rst.base, fun=sum, na.rm=TRUE)
  # #rst.max <- terra::aggregate(rst.masked,fact=round(xres(rst.disaggcrop)/xres(rst.reproj),0), fun=max)
  # rst.disaggcrop <- terra::crop(rst.disagg,rst.masked)
  # rst.avg <- terra::aggregate(rst.masked,fact=round(xres(rst.disaggcrop)/xres(rst.reproj),0), fun=mean, na.rm=TRUE)
  # #rst.sd <- terra::aggregate(rst.masked,fact=round(xres(rst.disaggcrop)/xres(rst.reproj),0), fun=sd, na.rm=TRUE)
  # rst.temp <- terra::resample(rst.avg,rst.disagg, method="ngb")
  # rst.temp2 <- rst.temp>0
  # rst.rice <- terra::aggregate(rst.temp2, fact=10, fun=sum, na.rm=TRUE)
  # rst.rice <- terra::resample(rst.rice, rst.base,method="ngb")
  
}





# TODO: Separate sos by month
blksize <- blockSize(rst.reproj)
for (m in unique(rat.sos$month)){
  month.ids <- unique(rat.sos$ID[rat.sos$month==m])
  fname.month <- paste0(CACHE_DIR,"/", paste(this.season$region[k], this.season$year[k], this.season$season[k], month.abb[m], sep="_"), ".tif")
  if(!file.exists(fname.month)) {
    rst.month <- raster(rst.reproj)
  } else {
    next
    rst.month <- raster(fname.month)
  }
  for(r in 1:blksize$n){
    #yFromRow(rst.reproj,blksize$row[2]), yFromRow(rst.masked,blksize$row[r])
    message("Month ", m, "- Block ", r, ": Extracting values.")
    vals <- values(rst.reproj,row=blksize$row[r], nrows=blksize$nrows[r])
    cells.block <- seq(min(cellFromRow(rst.reproj, blksize$row[r])), max(cellFromRow(rst.reproj, blksize$row[r]+blksize$nrows[r]-1)))
    message("Month ", m, "- Block ", r, ": Filtering values within month.")
    cells.block <- cells.block[vals %in% month.ids]
    if(length(cells.block)>0) {
      message("Month ", m, "- Block ", r, ": inserting values to monthly raster.")
      rst.month[cells.block] <- vals[vals %in% month.ids]
    }
    rm(vals, cells.block)
    gc(reset=TRUE)
  }
  message("Month ", m, ": saving to disk.")
  rst.month <- writeRaster(rst.month,fname.month, datatype="INT1S", format="GTiff", options="COMPRESS=LZW", overwrite=TRUE)
  message("Month ", m, ": ", basename(fname.month), " saved. Cleaning up.")
  rm(rst.month)
  unlink(dir(paste0(tempdir(),"/raster"), full.names=TRUE))
  message("Month ", m, ": Done.")
}

shp.aoi <- raster::crop(shp.modis, rst.reproj)
tiles.aoi <- shp.aoi$tile
TMPLATE_PROD <- PRODUCTS[paste0("m", TARGET_RES)][[1]]
for (l in 1:length(tiles.aoi)){
  output_file <- paste0(OUTPUT_DIR, "/", paste(TMPLATE_PROD,tiles.aoi[l],years[i],seasons[j], sep="_"),".tif")
  if(!file.exists(output_file)){
    tmplate_file <- dir(paste0(MODIS_DIR,"/", tiles.aoi[l]), pattern = paste(TMPLATE_PROD, tiles.aoi[l],paste0("A",years[i]),sep="."), full.names = TRUE)[1]
    if(!is.na(tmplate_file)) {
      rst.base <- rast(tmplate_file)
      rst.disagg <- terra::disaggregate(rst.base, fact=10)
    }
    
    #rst.modal <- terra::aggregate(rst.masked,fact=round(xres(rst.disaggcrop)/xres(rst.reproj),0), fun=modal, na.rm=TRUE)
    rst.count <- terra::aggregate(rst.riceog,fact=round(xres(rst.base)/xres(rst.riceog),0), fun=sum, na.rm=TRUE)
    rst.count <- terra::aggregate(rst.masked,fact=round(xres(rst.base)/xres(rst.reproj),0), fun=countUnique)
    rst.count <- resample(rst.count,rst.base, method="ngb")
    #rst.max <- terra::aggregate(rst.masked,fact=round(xres(rst.disaggcrop)/xres(rst.reproj),0), fun=max)
    rst.disaggcrop <- terra::crop(rst.disagg,rst.masked)
    rst.avg <- terra::aggregate(rst.masked,fact=round(xres(rst.disaggcrop)/xres(rst.reproj),0), fun=mean, na.rm=TRUE)
    #rst.sd <- terra::aggregate(rst.masked,fact=round(xres(rst.disaggcrop)/xres(rst.reproj),0), fun=sd, na.rm=TRUE)
    rst.temp <- terra::resample(rst.avg,rst.disagg, method="ngb")
    rst.temp2 <- rst.temp>0
    rst.rice <- terra::aggregate(rst.temp2, fact=10, fun=sum, na.rm=TRUE)
    rst.rice <- terra::resample(rst.rice, rst.base,method="ngb")
    rst.rice <- terra::writeRaster(rst.rice, output_file)
  } 
  
}

#rst.reproj[rst.reproj[]<=1] <- NA 
#rst.file <- rast(fname.reproj)
rst.tmp <- raster(this.season$filename[k])

#blx <- blockSize(rst.tmp)
#l=1
#blk.extent <- extentFromCells(rst.tmp,cellFromRow(rst.tmp,blx$row[l]:(blx$row[l]+blx$nrow[l]-1)))
#rst.blk <- crop(rst.tmp,blk.extent)
#rst.blk <- projectRaster(rst.blk, crs=PROJ.MODIS)
#rm(rst.cells)
#gc(reset = TRUE)
#projection(rst.tmp) <- PROJ.CEA





#rst.modis <- raster(paste0("./indices/",TILE,paste("/MOD09A1", TILE, "A2019001", "EVI", "tif", sep = ".")))

rst.modis <- raster(paste0("./250m/indices/",TILE,paste("/MXD13Q1", TILE, "A2018001", "ndvi", "tif", sep = ".")))
rst.base <- raster::disaggregate(rst.modis, fact=10)

values(rst.base) <- matrix(rep(matrix(1:100, ncol = 10, byrow = TRUE),ncell(rst.modis)),ncol=ncol(rst.base))

shp.validator <- shapefile("./Validation/rice_20200831.shp")
#rst.rice <- raster(rst.base)

if(file.exists(paste0(VALIDATION_DIR, "/",paste("Rice_Raw", TILE, 2020, "GISTDA", "THA", sep = "_") ,".rds"))){
  dat.valid <- readRDS(file=paste0(VALIDATION_DIR, "/",paste("Rice_Raw", TILE, 2020, "GISTDA", "THA", sep = "_") ,".rds"))
} else {
  for (i in nrow(shp.validator):1){
    message("Polygon ", i, ": Start.")
    shp.singlepoly <- shp.validator[i,]
    shp.singlepoly <- spTransform(shp.singlepoly, CRSobj = projection(rst.base))
    
    rst.poly <- try(crop(rst.modis, shp.singlepoly,snap="out"), silent = TRUE)
    if(class(rst.poly)=="try-error") next
    #stop()
    # Create raster of subpixel_ids
    rst.tmp <- disaggregate(rst.poly, fact=10)
    subcell.row <- matrix(rep(matrix(1:100, ncol=10, byrow = T), ncol(rst.poly)), ncol=ncol(rst.tmp))
    values(rst.tmp) <- rep(t(subcell.row), nrow(rst.poly))
    
    message("Polygon ", i, ": Getting rice pixels.")
    
    rice.cells <- cellFromPolygon(rst.tmp, shp.singlepoly)[[1]]
    # xx <- raster(rst.tmp) 
    # xx[rice.cells] <- rst.tmp[rice.cells]

    message("Polygon ", i, ": ", length(rice.cells), " pixels found.")
    if(length(rice.cells)==0) next
    rice.xy <- xyFromCell(rst.tmp, rice.cells)
    dat.rice_base <- data.frame(cell.modis=cellFromXY(rst.modis, rice.xy), cell.base=rice.cells, subcell.modis=rst.tmp[rice.cells])
    rice.modcell <- unique(dat.rice_base$cell.modis)
    
    for(j in 1:length(rice.modcell)){
      dat.rice_modcell <- data.frame(cell=rice.modcell[j],subcell.rice=binIndextoString(dat.rice_base$subcell.modis[dat.rice_base$cell.modis==rice.modcell[j]]), stringsAsFactors = FALSE)
      if(j==1) dat.rice_mod <- dat.rice_modcell else dat.rice_mod <- rbind(dat.rice_mod, dat.rice_modcell)
      rm(dat.rice_modcell)
    }
    rm(dat.rice_base, rst.poly,rst.tmp)
    # 
    # rice.cells <- rice.cells[which(rst.poly[]==1)]
    # if(length(rice.cells)==0) next
    # if(!exists("rst.rice")) {
    #   rst.rice <- rst.poly
    # } else {
    #   rst.rice <- raster::merge(rst.rice, rst.poly)
    #   rst.rice <- raster::mosaic(rst.rice, rst.poly, fun=sum)
    # }
    # 
    
    message("Polygon ", i, ": Adding to tile data.")
    if(!exists("dat.valid")) {
      dat.valid <- dat.rice_mod
    } else {
      to.update <- match(dat.rice_mod$cell, dat.valid$cell)
      idx.update <- which(!is.na(to.update))
      if(length(idx.update)>0){
        for(cc in idx.update){
          subcell.new <- as.numeric(unlist(strsplit(dat.rice_mod$subcell.rice[cc], "")))
          subcell.old <- as.numeric(unlist(strsplit(dat.valid$subcell.rice[to.update[cc]], "")))
          subcell.old[which(subcell.new==1)] <- 1
          dat.valid$subcell.rice[to.update[cc]] <- paste(subcell.old, collapse = "")
        }
      }
      dat.valid <- rbind(dat.valid, dat.rice_mod[is.na(to.update),])
      
    }
    
    message("Polygon ", i, ": Done.")
  }
  dat.valid$ra_pct <- sapply(dat.valid$subcell.rice, binStringSum)
  saveRDS(dat.valid, file=paste0(VALIDATION_DIR, "/",paste("Rice_Raw_", TILE, "250m", 2020, "GISTDA", "THA", sep = "_") ,".rds"))
  xx <- raster(rst.modis)
  xx[dat.valid$cell] <- dat.valid$ra_pct
}

# if (file.exists(paste0(VALIDATION_DIR, "/",paste("Rice_SUMMARY", TILE, 2020, "GISTDA", "THA", sep = "_") ,".rds"))){
#   dat.summary <- readRDS(file=paste0(VALIDATION_DIR, "/",paste("Rice_SUMMARY", TILE, 2020, "GISTDA", "THA", sep = "_") ,".rds"))
# } else {
#   dat.summary <- aggregate(rice ~ cell_50m + x + y, dat.valid, sum)
#   dat.summary$cell500m <- cellFromXY(rst.modis, dat.summary[,c("x","y")])
#   dat.summary$ra_cnt <- 1 
#   
#   saveRDS(dat.summary, file=paste0(VALIDATION_DIR, "/",paste("Rice_SUMMARY", TILE, 2020, "GISTDA", "THA", sep = "_") ,".rds"))
# }




#dat.modis <- aggregate(ra_cnt ~ cell500m, dat.summary, sum) 
rst.vrice <- raster(rst.modis)
rst.vrice[dat.valid$cell] <- dat.valid$ra_pct

writeRaster(rst.vrice, paste0(VALIDATION_DIR, "/", paste("RICEHOMOGENEITY",TILE , "250m", "THA", "GISTDA", "082020", "tif", sep=".")), overwrite=TRUE)


`# rst.cropper <- projectRaster(rst.modis, crs = projection(shp.validator))
# shp.validator <- crop(shp.validator,rst.cropper)
# shp.validator <- spTransform(shp.validator, CRSobj = projection(rst.modis))
