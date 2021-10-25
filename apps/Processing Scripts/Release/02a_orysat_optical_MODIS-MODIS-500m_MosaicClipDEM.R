#
#  Title: ASTER DEM Mosaicer-Clipper
#  Description: This script stitches ASTER DEM tiles for a specific region. The required parameters are
#  DEM_HOME = Directory where all ASTER DEM zip files are located
#  AOI      = Area of interest which can be a shapefile, raster or extent vector (xmin, xmax, ymin, ymax) 
#  AOI_TYPE = Indicate if a file is shapefile, raster or extent
#  

# System Settings
OUTPUT_DIR  = "F:/MODIS/dem/250m/"
WORKSPACE   = "E:/WORKSPACE/FLAGSHIP"
OUTPUT_PROJ = "auto"
OUTFILE_PREFIX = "h26v06"
# DEM Settings
DEM_HOME       = "F:/ASTER2/"
DEM_COMPRESSED = TRUE
COMPRESS_TYPE  = "zip"


# Area of Interest Settings
AOI        = "E:/WORKSPACE/Demeter Initiative/250m/Indices/h26v06/evi.250m.h26v06.A2014001.tif"
AOI_TYPE   = "raster"
AOI_BUFFER =  5


setwd(WORKSPACE)
library(raster)
library(terra)
library(manipulateR)
DEM_PROJ       = crs(rast())
# Create database of DEM files. Make sure the DEM_HOME a certain collection of data (e.g. ASTER2 or SRTM) 
# and of the same file format (e.g. tif or zip) 

if(!file.exists("ASTER2_FileInventory.rds")){
  DEM_FILES <- dir(DEM_HOME, pattern=ifelse(DEM_COMPRESSED, COMPRESS_TYPE, "tif"), recursive = TRUE)
  dat.uncompress <- sapply(paste0(DEM_HOME,DEM_FILES), unzip, list=TRUE, simplify = FALSE)
  dat.uncompress <- mapply(labelResult, dat.uncompress, label=as.list(DEM_FILES), stringsAsFactors=FALSE, SIMPLIFY = FALSE)
  dat.uncompress <- do.call(rbind, dat.uncompress)
  dat.uncompress <-  dat.uncompress[-grep("num", dat.uncompress$Name),]
  dat.uncompress$FILETYPE <- substr(dat.uncompress$Name, nchar(dat.uncompress$Name)-2, nchar(dat.uncompress$Name))
  dat.uncompress <- dat.uncompress[dat.uncompress$FILETYPE %in% c("zip","tif"),]
  dat.uncompress$Coords <- substr(basename(dat.uncompress$Name), 9, 15)
  dat.uncompress$ymin <- as.numeric(substr(dat.uncompress$Coords, 2, 3))
  dat.uncompress$ymin[grepl("S", dat.uncompress$Coords)] <- dat.uncompress$ymin[grepl("S", dat.uncompress$Coords)] * -1
  dat.uncompress$xmin <- as.numeric(substr(dat.uncompress$Coords, 5, 7))
  dat.uncompress$xmin[grepl("W", dat.uncompress$Coords)] <- dat.uncompress$xmin[grepl("W", dat.uncompress$Coords)] * -1
  saveRDS(dat.uncompress, file = "ASTER2_FileInventory.rds")  
} else dat.uncompress <- readRDS("ASTER2_FileInventory.rds")

aoi <- switch(AOI_TYPE, raster=rast(AOI), shpfile=shapefile(AOI), extent=extent(AOI))
if(AOI_BUFFER>0) {
  aoi.extended <- extend(aoi, rep(xres(aoi)*1.5,2))
  aoi.extended <- classify(aoi.extended, t(c(NA, 0)))
  aoi.extent <- extent(raster(aoi))
}

if(AOI_TYPE!="extent" & crs(aoi)!=DEM_PROJ){
  aoi_demmatch  <- switch(AOI_TYPE, raster=project(aoi.extended, rast()), shpfile=spTransform(aoi, CRSobj = DEM_PROJ))
  aoi_demmatch <- crop(aoi_demmatch, extentFromCells(aoi_demmatch, which(!is.na(aoi_demmatch[])))) 
  aoi_demmatch <- !is.na(aoi_demmatch)
} else aoi_demmatch <- aoi

if(OUTPUT_PROJ=="auto") {
  OUTPUT_PROJ <- projection(raster(aoi))
}
to.mosaic <- subset(dat.uncompress, xmin>=floor(xmin(aoi_demmatch)) & xmin<=ceiling(xmax(aoi_demmatch)) & ymin>=floor(ymin(aoi_demmatch)) & ymin<=ceiling(ymax(aoi_demmatch)))
to.mosaic$cell <- cellFromXY(aoi_demmatch, as.matrix(to.mosaic[,c("xmin", "ymin")]+0.5))
to.mosaic <- to.mosaic[!is.na(to.mosaic$cell),]
# xx <- raster(aoi_demmatch)
# xx[to.mosaic$cell] <- 1:nrow(to.mosaic)
to.mosaic$include <- aoi_demmatch[to.mosaic$cell]
to.mosaic <- to.mosaic[to.mosaic$include==1, ]
to.mosaic <- to.mosaic[order(to.mosaic$xmin,to.mosaic$ymin), ]

if(exists("aoi.dem")) rm(aoi.dem)
vals.aoidem <- done <- vector()
gc(reset=TRUE)

for(zp in unique(to.mosaic$label)){
  message(zp, ": Extracting files.")
  this.zip <- to.mosaic[to.mosaic$label==zp,] 
  #this.zip$cell <- cellFromXY(aoi_demmatch, as.matrix(this.zip[,c("xmin", "ymin")]+0.5))
  unzip(paste0(DEM_HOME, zp), files = this.zip$Name, junkpaths = TRUE)
  
  for(i in 1:nrow(this.zip)){
    # if(is.na(aoi_demmatch[this.zip$cell[i]])){
    #   if(file.exists(this.zip$Name[i])) file.remove(this.zip$Name[i])
    #   next
    # }
    if(this.zip$Name[i] %in% done) next
    if(this.zip$FILETYPE[i]=="zip"){
      lst.files <- unzip(this.zip$Name[i], list = TRUE)
      lst.files <- grep("dem", lst.files$Name, value = TRUE)
      unzip(this.zip$Name[i], files = lst.files, junkpaths = TRUE)
      rst.dem <- raster(basename(lst.files))
      #file.remove(this.zip$Name[i])
    } else rst.dem <- raster(this.zip$Name[i])
    
    message(zp, ": ", basename(lst.files), " reprojecting.")
    rst.dem <- projectRaster(rst.dem, crs = OUTPUT_PROJ)
    dem.extent <- extent(rst.dem)
    rst.dem <- rasterToRast(rst.dem)
    crs(rst.dem) <- crs(aoi) 
    if(is.null(intersect(dem.extent, aoi.extent))) {
      rm(rst.dem)
      file.remove(basename(lst.files))
      next
    }
    #stop("Finally")
    
    rst.dem <- aggregate(rst.dem, fact=ceiling(xres(aoi)/xres(rst.dem)))
    
    rst.piece <- crop(aoi, rst.dem)
    rst.dem <- crop(rst.dem,rst.piece)
    
    message(zp, ": ", basename(lst.files), " resampling.")
    rst.demrs <- resample(rst.dem, rst.piece)
    
    
    message(zp, ": ", basename(lst.files), " merging.")
    
    vals.dem <- values(rst.demrs)
    cells.dem <- which(!is.na(vals.dem))
    if(length(cells.dem)==0) next
    xy.dem <- xyFromCell(rst.demrs, cells.dem)
    aoi.cells <- cellFromXY(aoi, xy.dem)
    aoi.cells <- as.integer(aoi.cells)
    aoi.dem <-vals.dem[cells.dem]*1000
    aoi.dem <-round(aoi.dem,0)
    aoi.dem <-as.integer(aoi.dem)
    vals.thisastertile <- data.frame(cell=aoi.cells, dem=aoi.dem)
    
    vals.aoidem <- rbind(vals.aoidem, vals.thisastertile)
    vals.aoidem <- aggregate(dem ~ cell, vals.aoidem, FUN=mean, na.rm=TRUE)
    # if(!exists("vals.aoidem")){
    #   
    # } else {
    #   #aoi.dem <- merge(aoi.dem, rst.demrs)
    #   aoi.dem <- extend(aoi.dem, rst.demrs)
    #   
    #   aoi.dem <- 
    # }
    # 
    done <- c(done,this.zip$Name[i])
    message(zp, ": ", basename(lst.files), " cleaning.")
    file.remove(basename(lst.files))
    rm(rst.dem, rst.demrs, rst.piece, vals.dem, cells.dem, xy.dem, aoi.cells, aoi.dem, vals.thisastertile)
    gc(reset=TRUE)
  }
}

aoi.dem <- rast(aoi)
vals <- rep(NA, ncell(aoi.dem))
vals[vals.aoidem$cell] <- vals.aoidem$dem
values(aoi.dem) <- vals
aoi.dem <-  focal(aoi.dem, na.only=TRUE, fun="mean")


#aoi <- crop(aoi, og_extent)
#aoi.dem <- crop(aoi.dem, aoi)
#aoi.dem <- resample(aoi.dem, aoi, method="ngb")
aoi.dem <- writeRaster(aoi.dem, filename = paste0(OUTPUT_DIR, OUTFILE_PREFIX, "_ASTER2DEM.tif"), wopt=list(datatype="INT4S", gdal="COMPRESS=LZW"), overwrite=TRUE)
aoi.dem <- aoi.dem/1000
aoi.slope <- slope(aoi.dem, unit="degrees")
aoi.slope <- round(aoi.slope*100)
aoi.slope <- writeRaster(aoi.slope, filename = paste0(OUTPUT_DIR, OUTFILE_PREFIX, "_ASTER2SLOPE.tif"), wopt=list(datatype="INT2U", gdal="COMPRESS=LZW"), overwrite=TRUE)
