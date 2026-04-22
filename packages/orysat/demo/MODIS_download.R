# NOTE: For windows, make sure you have wget executable and add it to path.
# TODO: Create a script that will transfer downloaded files to MODIS_STORE
# MODIS_STORE = "/media/jaunario/warehouse/raw/MODIS"
setwd("/media/jaunario/omnistore/raw/MODIS")
setwd("/media/jaunario/scratch/MODIS")
MODIS_HOME   = "."
TILE        = "h30v07" #, "h28v08")      # c("h25v06","h26v06", "h27v07", "h28v07", "h29v07") #c("h26v06", "h25v06") # c("h30v08", "h30v07", "h29v08", "h29v07") # specify TILES to download using modis tile code, seprated by comma
PRODUCTS     = c("MOD13Q1", "MYD13Q1") #, "MYD16A2") # specify PRODUCTS to download using product short name, separated by comma
PROD_VERSION = 61

# TODO: Update RemoteSensing Package or Integrate all in orysat
library(orysat)

YEARS = 2021:2023 

info.runner <- Sys.info() 
#runner_id <- ifelse(grepl("rpi-thrall", info.runner["nodename"]), as.numeric(sub("rpi-thrall", "", info.runner["nodename"])), 0) 
runner_id <- 0
prog <- 0
if(!file.exists(paste0(MODIS_HOME,"/modis_inventory_", sprintf("%03.1f",PROD_VERSION)  ,".rds"))){
  files.modis <- dir(paste0(MODIS_HOME, "/v", PROD_VERSION), pattern = ".hdf$", full.names = TRUE, recursive = TRUE)
  dat.modisinv <- orysat::inventory.modis(files.modis, modisinfo=c('product', 'acqdate', 'tile', 'version', 'proddate'), file.ext="hdf")
  dat.modisinv$action <- "store"
  saveRDS(dat.modisinv,file = paste0(MODIS_HOME,"/modis_inventory_", sprintf("%03.1f",PROD_VERSION)  ,".rds"))
} else {
  dat.modisinv <- readRDS(paste0(MODIS_HOME,"/modis_inventory_", sprintf("%03.1f",PROD_VERSION)  ,".rds"))
}

for (pr in PRODUCTS){

  prod.doy <- orysat::modis.acqdoys(pr)
  
  for (yr in YEARS){
    
    if(runner_id==0){
      runner_doy <- prod.doy
    } else {
      runner_doy <- prod.doy[seq(runner_id,length(doy), by=3)]
    }
    
    if(length(runner_doy)==1 & (manipulateR::yearFromDate(Sys.Date())-yr)<2) {
      message("MODIS Product, ", pr, " is probably not yet available for year ", yr)
      next
    }
    # status.action is variable for progress monitoring. Pass to network log (gsheet) or Endeavor System Cache
    status.action <- "start"
    success <- 0
    for (i in runner_doy) {
      done.recs <- subset(dat.modisinv, acqdate==sprintf("A%g%03g", yr, i) & product==pr & tile==TILE)
      if(nrow(done.recs)>0) next
      #stop()
      message(pr, ": Downloading ", paste0(yr, " -", prettyNum(i,width=3)))

      # result <- try(orysat::download.modis2(tile=TILE,
      #                                       product=pr,
      #                                       year=yr,
      #                                       doy = i,
      #                                       savepath=paste(MODIS_HOME, paste0("v", PROD_VERSION), pr, TILE, sep="/"),
      #                                       prod.ver=PROD_VERSION,
      #                                       skip.exists=TRUE,
      #                                       site.url = "https://e4ftl01.cr.usgs.gov/",
      #                                       method="wget", mode="wb")
      #                                       )

      result <- try(orysat::download.modis2(tile=TILE,
                                    product=pr,
                                    year=yr,
                                    doy = i,
                                    savepath=paste(MODIS_HOME, paste0("v", PROD_VERSION), pr, TILE, sep="/"),
                                    prod.ver=PROD_VERSION,
                                    skip.exists=TRUE,
                                    site.url = "https://e4ftl01.cr.usgs.gov/",
                                    method="curl", mode="wb",
                                    extra=c("-0", "-L", "-n",
                                            "-b ~/.urs_cookies",
                                            "-c ~/.urs_cookies",
                                            "--expect100-timeout 200",
                                            "--retry 10",
                                            "--retry-max-time 400")
                                    ))

      # result <- try(orysat::download.modis(tile=TILE,
      #                                       product=pr,
      #                                       year=yr,
      #                                       doy = i,
      #                                       savepath=paste(MODIS_HOME, paste0("v", PROD_VERSION), pr, TILE, sep="/"),
      #                                       prod.ver=PROD_VERSION,
      #                                       skip.exists=TRUE,
      #                                       site.url = "https://e4ftl01.cr.usgs.gov/",
      #                                       handle=usgs.handle))
      
      if(is.character(result) && file.exists(result)){
        success <- success + 1
        new.dl <- orysat::inventory.modis(result, modisinfo=c('product', 'acqdate', 'tile', 'version', 'proddate'), file.ext="hdf")
        new.dl$action <- "store"
        dat.modisinv <- rbind(dat.modisinv, new.dl)
        prog <- round(success/length(runner_doy)*100,2)
        message(pr, ": Downloaded ", prettyNum(prog,width=6), "%")
        #manipulateR::timer.message("Server cooldown ", delaytime = 150)
      } else {
        if(yr==manipulateR::yearFromDate(Sys.Date())) break
      }

      if(status.action=="start") status.action <- "wait"
      
      if(prog>=100){
        action <- "none"
      } 
    }  
  }
}

#CLEANUP
rm(action, success, doy, result, YEARS)
gc(reset=TRUE)
saveRDS(dat.modisinv, file = paste0(MODIS_HOME,"/modis_inventory_", sprintf("%03.1f",PROD_VERSION)  ,".rds"))

