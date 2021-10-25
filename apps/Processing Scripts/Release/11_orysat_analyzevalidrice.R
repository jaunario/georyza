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

library(ggplot2)
library(orysat)
pixels.toprocess <- readRDS(file=paste0(CACHE_DIR, "/NUM_INDEX_", paste(TILE, YEAR, "IDX_TOPROC", "rds", sep=".")))
crops <- readRDS(file=paste0(CACHE_DIR, "/", paste("CROP", METHOD,TILE, YEAR, sep = "_") ,".rds"))


rst.vrice <- raster(paste0(VALIDATION_DIR, "/", paste("RICEHOMOGENEITY",TILE , "50m", "THA", "GISTDA", "082020", "tif", sep=".")))
pixels.rice <- which(rst.vrice[]>50)

dat.rice <- data.frame(cell=pixels.rice, ra_pct=rst.vrice[pixels.rice])

ws.crop <- crops[crops$sos>27,]

dat.joined <- plyr::join(dat.rice, ws.crop, by="cell")
plot(dat.joined$ra_pct, dat.joined$crop.duration)

p <- ggplot(data = dat.joined, aes(x=ra_pct, y=crop.duration))

xx <- crops[crops$cell %in% pixels.rice,]
