# NOTE: For windows, make sure you have wget executable and add it to path.
DOWNLOAD_DIR = "F:/MODIS/v6" # Path to where you want to store downloaded images
TILES        =  "h28v07" #c("h26v06", "h25v06") # c("h30v08", "h30v07", "h29v08", "h29v07") # specify TILES to download using modis tile code, seprated by comma
PRODUCTS     = "MOD09A1" # c("MCD12Q1") #, "MYD13Q1") #, "MYD16A2") # specify PRODUCTS to download using product short name, separated by comma
YEARS        = 2020
USGS.USERPWD = "geospatial.irri:8&-yFm*v2R-s?_!"

library(curl)
library(orysat)
library(manipulateR)

setwd(DOWNLOAD_DIR) 

modprods <- read.csv(system.file("satPRODUCTS/modis.PRODUCTS.ref.csv", package="orysat"), stringsAsFactors=FALSE)

for (i in 1:length(TILES)){
  for (j in 1:length(PRODUCTS)) {
    doy <-  seq(1,366, by=8)
    for (l in 1:length(YEARS)){
      for (k in 1:length(doy)) {
        withRetry(download.modis(tile=TILES[i],
                                 product=PRODUCTS[j],
                                 years=YEARS[l], userpwd=USGS.USERPWD,
                                 doy = doy[k],
                                 savedir=paste(PRODUCTS[j], TILES[i],sep="/"),
                                 prod.ver=6 ,
                                 skip.exists=TRUE,
                                 modis.site = "https://e4ftl01.cr.usgs.gov/"),
                  delay = 5)
      }
    }
  }  
}
