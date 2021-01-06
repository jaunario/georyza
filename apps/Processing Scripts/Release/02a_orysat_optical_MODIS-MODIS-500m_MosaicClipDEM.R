#
#  Title: ASTER DEM Mosaicer-Clipper
#  Description: This script stitches ASTER DEM tiles for a specific region. The required parameters are
#  DEM_HOME = Directory where all ASTER DEM zip files are located
#  AOI      = Area of interest which can be a shapefile, raster or extent vector (xmin, xmax, ymin, ymax) 
#  AOI_TYPE = Indicate if a file is shapefile, raster or extent
#  

# System Settings
OUTPUT_DIR  = "F:/MODIS/dem/500m/"
WORKSPACE   = "E:/WORKSPACE/FLAGSHIP"
OUTPUT_PROJ = "auto"
OUTFILE_PREFIX = "h29v07"
# DEM Settings
DEM_HOME       = "F:/ASTER2/"
DEM_COMPRESSED = TRUE
COMPRESS_TYPE  = "zip"
DEM_PROJ       = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# Area of Interest Settings
AOI        = "E:/WORKSPACE/FLAGSHIP/indices/h27v07/MOD09A1.h27v07.A2016001.EVI.tif"
AOI_TYPE   = "raster"
AOI_BUFFER =  5


setwd(WORKSPACE)
library(raster)
# Create database of DEM files. Make sure the DEM_HOME a certain collection of data (e.g. ASTER2 or SRTM) 
# and of the same file format (e.g. tif or zip) 
DEM_FILES <- dir(DEM_HOME, pattern=ifelse(DEM_COMPRESSED, COMPRESS_TYPE, "tif"), recursive = TRUE)


if(!file.exists("ASTER2_FileInventory.rds")){
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

aoi <- switch(AOI_TYPE, raster=raster(AOI), shpfile=shapefile(AOI), extent=extent(AOI))
if(AOI_BUFFER>0) {
  og_extent <- extent(aoi)
  aoi <- extend(aoi, rep(AOI_BUFFER,2))
}
if(AOI_TYPE!="extent" & projection(aoi)!=DEM_PROJ){
  aoi_demmatch  <- switch(AOI_TYPE, raster=projectRaster(aoi, crs=DEM_PROJ), shpfile=spTransform(aoi, CRSobj = DEM_PROJ))
   
} else aoi_demmatch <- aoi
if(OUTPUT_PROJ=="auto") OUTPUT_PROJ <- projection(aoi)

to.mosaic <- subset(dat.uncompress, xmin>=floor(xmin(aoi_demmatch)) & xmin<ceiling(xmax(aoi_demmatch)) & ymin>=floor(ymin(aoi_demmatch)) & ymin<ceiling(ymax(aoi_demmatch)))
to.mosaic <- to.mosaic[order(to.mosaic$label), ]

if(exists("aoi.dem")) rm(aoi.dem)
for(zp in unique(to.mosaic$label)){
  message(zp, ": Extracting files.")
  this.zip <- to.mosaic[to.mosaic$label==zp,] 
  unzip(paste0(DEM_HOME, zp), files = this.zip$Name, junkpaths = TRUE)
  
  for(i in 1:nrow(this.zip)){
    if(this.zip$FILETYPE[i]=="zip"){
      lst.files <- unzip(this.zip$Name[i], list = TRUE)
      lst.files <- grep("dem", lst.files$Name, value = TRUE)
      unzip(this.zip$Name[i], files = lst.files, junkpaths = TRUE)
      rst.dem <- raster(basename(lst.files))
      file.remove(this.zip$Name[i])
    } else rst.dem <- raster(this.zip$Name[i])
    message(zp, ": ", basename(lst.files), " reprojecting.")
    rst.dem <- projectRaster(rst.dem, crs = OUTPUT_PROJ)
    if(is.null(intersect(extent(aoi),extent(rst.dem)))) {
      file.remove(basename(lst.files))
      next
    }
    rst.dem <- raster::aggregate(rst.dem, fact=ceiling(res(aoi)[1]/res(rst.dem)[1]))
    
    rst.piece <- crop(aoi, rst.dem)
    rst.dem <- crop(rst.dem,rst.piece)
    
    message(zp, ": ", basename(lst.files), " resampling.")
    rst.demrs <- raster::resample(rst.dem, rst.piece)
    
    
    message(zp, ": ", basename(lst.files), " merging.")
    if(!exists("aoi.dem")){
      aoi.dem <- rst.demrs
    } else {
      aoi.dem <- merge(aoi.dem, rst.demrs)
    }
    message(zp, ": ", basename(lst.files), " cleaning.")
    file.remove(basename(lst.files))
  }
  
}
aoi <- crop(aoi, og_extent)
aoi.dem <- crop(aoi.dem, aoi)
aoi.dem <- resample(aoi.dem, aoi, method="ngb")
aoi.dem <- writeRaster(aoi.dem, filename = paste0(OUTPUT_DIR, OUTFILE_PREFIX, "_ASTER2DEM.tif"))
aoi.slope <- terrain(aoi.dem, unit="degrees")
aoi.slope <- writeRaster(aoi.slope, filename = paste0(OUTPUT_DIR, OUTFILE_PREFIX, "_ASTER2SLOPE.tif"))
