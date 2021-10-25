# NOTE: For windows, make sure you have wget executable and add it to path.
MODIS_HOME   = "F:/MODIS"   # Path to where you want to store downloaded images
TILES        = "h26v06"  #c("h25v07","h25v06", "h26v06", "h26v07")     #c("h25v06","h26v06", "h27v07", "h28v07", "h29v07") #c("h26v06", "h25v06") # c("h30v08", "h30v07", "h29v08", "h29v07") # specify TILES to download using modis tile code, seprated by comma
PRODUCTS     = "MCD12Q1" # c("MOD13Q1", "MYD13Q1") #, "MYD16A2") # specify PRODUCTS to download using product short name, separated by comma
PROD_VERSION = 6
YEARS        = 2002:2006
USGS.USERPWD = "geospatial.irri:8&-yFm*v2R-s?_!"
USGS.USERPWD = "jaunario:Ragnarok09"

library(curl)
library(orysat)
library(manipulateR)

modprods <- read.csv(system.file("satPRODUCTS/modis.PRODUCTS.ref.csv", package="orysat"), stringsAsFactors=FALSE)

arg.quote <- function(x){
  return(paste0("'", x, "'"))
}

for (i in 1:length(TILES)){
  for (j in 1:length(PRODUCTS)) {
    prod.info <- modis.productinfo(product = PRODUCTS[j])
    prod.site <- paste("https://e4ftl01.cr.usgs.gov/", "MO", switch(prod.info$Platform, Aqua="LA", Terra="LT", Combined="TA"), "/", paste(PRODUCTS[j],sprintf(paste("%03d",sep=""),PROD_VERSION),sep="."),sep="")
    prod.dateinv <- readLines(curl(prod.site))
    dateinv.stdate <- as.Date(strsplit(prod.dateinv[grep("Parent Directory", prod.dateinv)+1],">")[[1]][3], "%Y.%m.%d/</a")
    dateinv.endate <- as.Date(strsplit(prod.dateinv[length(prod.dateinv)-2],">")[[1]][3], "%Y.%m.%d/</a")
    prod.interval <- modprods$Temporal.Granularity[modprods$ShortName==PRODUCTS[j]]
    acq.by <- switch(prod.interval, Daily=1, Yearly=365, "4 day"=4, "8 day"=8, "16 day"=16, NA)
    if(is.na(acq.by)){
      message("Unsupported MODIS Product. Skipping.")
      next
    } else if (acq.by!=16){
      doy <-  seq(1,365, by=acq.by)
    } else {
      mod.platform <- modprods$Platform[modprods$ShortName==PRODUCTS[j]]
      if(mod.platform=="Terra"){
        doy <-  seq(1,366, by=acq.by)
      } else if(mod.platform=="Aqua"){
        doy <-  seq(9,366, by=acq.by)
      } else {
        message("Unsupported MODIS Product. Skipping.")
        next
      }
    }

    for (l in 1:length(YEARS)){
      action <- "start"
      # system(paste("RScript 00a_orysat_RiceMap-Flagship_SESSION.R", paste(
      #   paste0("MODE=", arg.quote("DOWNLOAD")), 
      #   paste0("TILE=",arg.quote(TILES[i])),
      #   paste0("PRODUCT=", arg.quote(PRODUCTS[j])),
      #   paste0("YEAR=",YEARS[l]),
      #   paste0("PROGRESS=",0),
      #   paste0("ACTION=",arg.quote(action)), collapse = " ")
      # ), intern = FALSE,  show.output.on.console = FALSE)
      # 
      success <- 0
      for (k in 1:length(doy)) {
        acqdate <- dateFromDoy(doy[k],YEARS[l])
        if(acqdate<dateinv.stdate | acqdate>dateinv.endate) next
        if(length(dir(paste(MODIS_HOME, paste0("v", PROD_VERSION), PRODUCTS[j], TILES[i],sep="/"),  pattern=sprintf("%s.A%04d%03d.%s.", PRODUCTS[j],YEARS[l],doy[k],TILES[i])))>0) next
        result <- withRetry(download.modis(tile=TILES[i],
                                 product=PRODUCTS[j],
                                 years=YEARS[l], userpwd=USGS.USERPWD,
                                 doy = doy[k],
                                 savepath=paste(MODIS_HOME, paste0("v", PROD_VERSION), PRODUCTS[j], TILES[i],sep="/"),
                                 prod.ver=PROD_VERSION,
                                 skip.exists=TRUE,
                                 modis.site = "https://e4ftl01.cr.usgs.gov/"),
                  delay = 5, retries = 3)
        
        if(file.exists(result)){
          if(action=="start") action <- "wait"
          success <- success + 1
          prog <- round(success/length(doy)*100,2)
        } else {
          action <- "rerun"
        }
        
        if(prog>=100){
          action <- "none"
        } 
        
        # if(!exists("GSHEET_ID") || GSHEET_ID==""){
        #     system(paste("RScript 00a_orysat_RiceMap-Flagship_SESSION.R", paste(
        #     paste0("MODE=", arg.quote("DOWNLOAD")), 
        #     paste0("TILE=",arg.quote(TILES[i])),
        #     paste0("PRODUCT=", arg.quote(PRODUCTS[j])),
        #     paste0("YEAR=",YEARS[l]),
        #     paste0("PROGRESS=",prog),
        #     paste0("ACTION=",arg.quote(action)), collapse = " ")
        #   ), intern = FALSE,  show.output.on.console = FALSE)
        # }
      }
      # if(success==0){
      #   action <- "debug"
      #   system(paste("RScript 00a_orysat_RiceMap-Flagship_SESSION.R", paste(
      #     paste0("MODE=", arg.quote("DOWNLOAD")), 
      #     paste0("TILE=",arg.quote(TILES[i])),
      #     paste0("PRODUCT=", arg.quote(PRODUCTS[j])),
      #     paste0("YEAR=",YEARS[l]),
      #     paste0("PROGRESS=",prog),
      #     paste0("ACTION=",arg.quote(action)), collapse = " ")
      #   ), intern = FALSE,  show.output.on.console = FALSE)
      # }      
    }
  }
}  

