# TODO: Add comment
# 
# Author: jaunario
###############################################################################


ricestress.dayheat <- function(wth, flwrdate, tmax.threshold=35, nstressdays.threshold=10, sf.method="Horie", sf.threshold=0.75, cool.effect=TRUE, summaries.only=TRUE, consecutive=FALSE, anthesis=seq(0, 23.9, by=0.1),tpeakfl=c(9.5,10,10.5,11,11.5), ...){
  stressed <- FALSE

  # Flowering window
  flwrstart <- flwrdate-15
  flwrend <- flwrdate+21
  dates <- t(c(flwrdate, flwrstart, flwrend))

  if(class(wth) == "weather") {
    x <- wth@lon
    y <- wth@lat
    wth <- wth@w
  } else {
    x <- 0
    y <- 0
  }
  weather <- subset(wth, date >= flwrstart & wth$date <= flwrend)

  if(!is.na(sf.method)){
  # Check data if weather parameters are available
    prereq <- c("tmax", "tmin","tdew")
    if (sum(prereq %in% colnames(weather))!=length(prereq)) {
      stop(paste(prereq[!(prereq%in%colnames(weather))], collapse=", "), "not found in weather data.")
    }

    for (i in 1:nrow(weather)){
      sol.time <- solar.time(lat=y, lon=x, date=weather$date[i])
      TAir <- sapply(anthesis, ttime, tmax=weather$tmax[i], tmin=weather$tmin[i], t.sunrise=sol.time$rise, t.sunset=sol.time$set, DL=sol.time$duration)
      TPan <- tpan(tdew=weather$tdew[i],tair=TAir, no.cooling=!cool.effect)
      FFL <- p.fltime(anthesis)
      SFPOP <- sum(FFL*spikelet.fertility(temp=TPan, method=sf.method))
      SFPEAK <- mean(spikelet.fertility(temp=tpan(tdew=weather$tdew[i],tair=TAir[which(anthesis %in% tpeakfl)], no.cooling=!cool.effect), method=sf.method))
      if(!exists("details")){
        details <- data.frame(date=weather$date[i], maxtpan=max(TPan), SFPOP, SFPEAK)
      } else {
        details <- rbind(details,data.frame(date=weather$date[i], maxtpan=max(TPan), SFPOP, SFPEAK))
      }
      rm(TAir,TPan,FFL,SFPOP,SFPEAK)
      gc(reset=TRUE, verbose=FALSE)
    }

    details <- merge(weather,details)
    #details$part.sterile.pop <- details$SFPOP<.75
    #details$high.sterile.pop <- details$SFPOP<.50
    #details$comp.sterile.pop <- details$SFPOP==0

    #details$part.sterile.pk <- details$SFPEAK<.75
    #details$high.sterile.pk <- details$SFPEAK<.50
    #details$comp.sterile.pk <- details$SFPEAK==0

    details$maxtpan.vthres <- details$maxtpan>=tmax.threshold
    details$maxtpan.xthres <- details$maxtpan-tmax.threshold
    details$maxtpan.xthres[details$maxtpan.xthres<0] <- NA

  } else {
    details <- weather
  }

  details$tmax.vthres <- details$tmax>=tmax.threshold 

  # STEP X: Original Threshold checking
  details$tmax.xthres <- details$tmax-tmax.threshold
  details$tmax.xthres[details$tmax.xthres<0] <- NA # Remove negative values



  # SUMMARY
  means <- c("tmax", "tmax.xthres", "maxtpan", "SFPOP", "SFPEAK") # vars to compute means and sd
  sums  <- c("tmax.vthres", "maxtpan.vthres") # vars to compute sums and maxConsecutive

  stress.summary <- c(apply(details[,means[means %in% colnames(details)]], 2, FUN = mean, na.rm=TRUE),
    apply(details[,means[means %in% colnames(details)]], 2, FUN = sd, na.rm=TRUE),
    sum(details[,"tmax.vthres"], na.rm=TRUE),
    maxConsecutive(details[,"tmax.vthres"]))

  names(stress.summary) <- c(paste(means[means %in% colnames(details)],"mean", sep="-"),
                             paste(means[means %in% colnames(details)],"sd", sep="-"),
                             "tmax.vthres-sum", "tmax.vthres-maxcons")

  tmax.stress <- ifelse(consecutive, stress.summary["tmax.vthres-maxcons"]>=nstressdays.threshold, stress.summary["tmax.vthres-sum"]>=nstressdays.threshold)
  og.names <- names(stress.summary)
  stress.summary <- c(stress.summary, tmax.stress)
  names(stress.summary) <- c(og.names,"tmax.stress")

  if(sum(grepl("SFPOP", names(stress.summary)))>0){
    og.names <- names(stress.summary)
    stress.summary <- c(stress.summary, sum(details[,"maxtpan.vthres"], na.rm=TRUE), maxConsecutive(details[,"maxtpan.vthres"]))
    names(stress.summary) <- c(og.names, "maxtpan.vthres-sum","maxtpan.vthres-maxcons")

    maxtpan.stress <- ifelse(consecutive, stress.summary["maxtpan.vthres-maxcons"]>=nstressdays.threshold, stress.summary["maxtpan.vthres-sum"]>=nstressdays.threshold)
    sfpop.stress <- stress.summary["SFPOP-mean"] < sf.threshold
    sfpeak.stress <- stress.summary["SFPEAK-mean"] < sf.threshold
    og.names <- names(stress.summary)
    stress.summary <- c(stress.summary, maxtpan.stress, sfpop.stress, sfpeak.stress)
    names(stress.summary) <- c(og.names,"maxtpan.stress","sfpop.stress","sfpeak.stress")
  }

  # stress.summary <- c(mean(details$tmax, na.rm=TRUE), # Average tmax
  #           sd(details$tmax, na.rm=TRUE), # Std. Dev tmax
  #           mean(details$tmax.xthres, na.rm=TRUE), # Average tmax excess from threshold
  #            sd(details$tmax.xthres, na.rm=TRUE), # Std. Dev tmax excess from threshold
  #           mean(details$maxtpan, na.rm=TRUE), # Average Max. Panicle Temperature
  #           sd(details$maxtpan, na.rm=TRUE), # Std. Dev Max. Panicle Temperature 
  #           mean(details$maxtpan.xthres, na.rm=TRUE), # Average maxtpan excess from threshold
  #           sd(details$maxtpan.xthres, na.rm=TRUE), # Std. Dev maxtpan excess from threshold
  #           sum(details$tmax.vthres), # Number of days tmax exceeding threshold - HOTDAYS based on TMax
  #           sum(details$maxtpan.vthres), # Number of days maxtpan exceeding threshold - HOTDAYS based on TPan
  #           maxConsecutive(details$tmax.vthres), # Max consecutive days tmax exceeded threshold
  #           maxConsecutive(details$maxtpan.vthres), # Max consecutive days maxtpan exceeded threshold
  #           mean(details$SFPOP, na.rm=TRUE), # Average spikelet fertility pop
  #           mean(details$SFPEAK, na.rm=TRUE), # Average spikelet fertility pk
  #           ifelse(consecutive, maxConsecutive(details$tmax.vthres)>=nstressdays.threshold, sum(details$tmax.vthres)>=nstressdays.threshold), # is stressed year if TMax Hotdays exceeds stressdays.threshold
  #           ifelse(consecutive, maxConsecutive(details$maxtpan.vthres)>=nstressdays.threshold, sum(details$maxtpan.vthres)>=nstressdays.threshold), # is stressed year if Max TPan Hotdays exceeds stressdays.threshold
  #           mean(details$SFPOP, na.rm=TRUE)<0.75,
  #           mean(details$SFPEAK, na.rm=TRUE)<0.75
  #           )
  # 
  # names(stress.summary) <- c("Average.TMax", 
  #               "StDev.TMax", 
  #               "Average.TMax.Excess.From.Threshold", 
  #               "StDev.TMax.Excess.From.Threshold", 
  #               "Average.MaxTPan", 
  #               "StDev.MaxTPan", 
  #               "Average.MaxTPan.Excess.From.Threshold", 
  #               "StDev.MaxTPan.Excess.From.Threshold",
  #               "Hotdays.based.on.TMax", 
  #               "Hotdays.based.on.MaxTpan", 
  #               "Max.Consecutive.Hotdays.based.on.TMax", 
  #               "Max.Consecutive.Hotdays.based.on.MaxTpan",
  #               "Average.spikelet.fertility.pop",
  #               "Average.spikelet.fertility.pk",
  #               "Stress.Based.on.TMax", "Stress.Based.on.MaxTPan", "Stress.Based.on.Spikelet.Fertility.Pop", "Stress.Based.on.Spikelet.Fertility.Pk"
  #               )

  if (summaries.only){
    out <- c(stress.summary)
    out <- data.frame(t(out))
    names <- c("flowering.date", "flwrwindow.start", "flwrwindow.end", colnames(out))
    out <- cbind(flwrdate, flwrstart, flwrend, out)
    colnames(out) <- names
  } else {
    out <- list(stress.summary, details)
  }
  return(out)
}

ricestress.nighttimeHeat <- function(wth, dstart, dend, consecutive=TRUE, tmin.threshold=25,...){ #sowdoy, harvdoy, year, 
  #tminLL <- 25    #threshold for nightime heat stress + 1.4 adjustment (Bai et al., 2010)                    
  tminLL <- tmin.threshold                       
  ndays <- 15    # days that temperature is above the threshold 
  
  # Output vector and names
  #sowdate <- dateFromDoy(sowdoy, year)
  #harvdate <- dateFromDoy(harvdoy, ifelse(sowdoy>harvdoy,year+1, year))
  #seeddate <- sowdate+21
  #flwrdate <- harvdate-30
  #dstart <- flwrdate-65                               # PI to maturity
  #dend <- harvdate

  z <- wth@w
  
  tminseed <-  periodWthVar(z, "tmin", dstart, dend)    
  
  tmindiff_seed <- tminseed-tminLL
  tmindiff_seed[tmindiff_seed<=0] <- NA # If tmin is less than tminLL 
  
  tnavg <- mean(tminseed, na.rm=TRUE) # Mean tmin for year specified
  tndavg <- mean(tmindiff_seed, na.rm=TRUE) # Mean tmin diff from limit
  
  potstress <- !is.na(tmindiff_seed)
  potStressCnt <- sum(potstress) # number of days tmin is above tmin_LL
  
  if (consecutive) LL_seed <- maxConsecutive(potstress) else LL_seed <- potStressCnt
  stress <- LL_seed >=ndays
  
  out <- c(tnavg, tndavg, potStressCnt, LL_seed, stress) 
  names(out) <- c("tmin.avg", "tminLL_excess.avg", "excess.cnt", "stress.cnt", "stressed")
  return(out)    
}

thermal.stress <- function(planting, dae, years, type, climdsn, climset, maskraster=NA, outdir=getwd(), verbose=TRUE, ...){
  sowdoy <- getValues(planting)
  daevals <- getValues(dae)
  cells <- which(!is.na(sowdoy) & sowdoy>0)
  xy <- xyFromCell(planting, cells)
  #stdcell <- cellsFromExtent(raster(),planting)[cells]
  
  stressdata <- vector()
  rem <- vector()    
  
  connection <- odbcConnect(climdsn) 
  for (i in 1:length(cells)){
    #if(cells[i]==25308) message("stop")
    if (verbose) show.message(type, ": Analyzing weather on cell ", sprintf("%05d",cells[i]), ".", eol="\r")
    if(type=="cold") {
      stressrec <- stress.cold(con=connection, table=climset, xy=t(xy[i,]), sowdoy=sowdoy[cells[i]], dae=daevals[cells[i]], years=years)
    } else if (type=="nighttimeheat"){
      stressrec <- stress.nighttimeHeat(con=connection, table=climset, xy=t(xy[i,]), sowdoy=sowdoy[cells[i]], dae=daevals[cells[i]], years=years, consecutive=FALSE)
    } else if (type=="daytimeheat"){
      stressrec <- stress.daytimeHeat(con=connection, table=climset, xy=t(xy[i,]), sowdoy=sowdoy[cells[i]], dae=daevals[cells[i]], years=years, consecutive=FALSE)
    }
    
    if (length(stressrec)>0) {
      stressdata <- rbind(stressdata,stressrec,deparse.level=0)
    } else {
      if (verbose) show.message(sprintf("%05d",cells[i]), ": No valid output ", eol="\n")
      rem <- c(rem,i)
    }
  }
  if (length(rem)>0) cells <- cells[-rem]
  rownames(stressdata) <- cells 
  stress.map(stressdata, planting, maskraster, writeto=outdir)
  return(TRUE)
}
