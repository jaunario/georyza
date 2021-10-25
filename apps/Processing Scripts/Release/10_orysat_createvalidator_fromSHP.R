# DIRECTORY SETTINGS
MODIS_HOME  = "F:/MODIS"
WORKSPACE   = "E:/WORKSPACE/FLAGSHIP"
VALIDATION_DIR = paste0(WORKSPACE, "/Validation")
CACHE_DIR    = paste0(WORKSPACE, "/cache")
OUTPUT_DIR  = "AUTO"


MISSING_THRES   = 50

METHOD     = "VWIDynamics"
VERSION    = 0.1 
YEAR       = 2019
TILE       = "h27v07"
PRODUCTS    = "MOD09A1"

# TODO: Have default layer-parameter mapping but allow custom mapping for experimental methods
LAYERPARAM_MAPPING  = list(vi="EVI", wi="MNDWI")

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

library(orysat)

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


# rst.cropper <- projectRaster(rst.modis, crs = projection(shp.validator))
# shp.validator <- crop(shp.validator,rst.cropper)
# shp.validator <- spTransform(shp.validator, CRSobj = projection(rst.modis))
