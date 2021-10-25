crop.extract <- function(vi, dates, vi.label="evi", interval=1, senescence.min=5, growth.min=6, seasonends.pctdiff=0.1, negligible.change=100, data.save=null, as.z = FALSE){
  #if(is.null(dates)) dates <- seq(1,366, by=8, lengt)
  crops <- NA
  if(!is.na(vi[1]) & (sum(is.na(vi))!=length(vi))){
    tile.h <- trunc(vi[1])
    tile.v <- (vi[1]-tile.h) * 100
    cell <- vi[2]
    
    vi <- vi[-(1:2)]
    
  
    #if(verbose) message(tile, "-", cell)
    if(as.z) vi <- (vi-mean(vi))/sd(vi) # Compute for z standard value
    
    dx <- diff(vi) # Find out which ones are increasing and decreasing
    
    crops <- 0
    
    trend.senescence <- consecutive.groups(which(dx <= 0)) # Find consecutively decreasing VI's indicating senescence
    trend.senescence <- lapply(trend.senescence, function(x) {return(c(x, x[length(x)] + 1))}) # Extends the index to the last value of the trend
    st.senescence <- sapply(trend.senescence, min)
    valid.ts <- which(!is.na(st.senescence))
    
    trend.growth <-consecutive.groups(which(dx > 0)) # Find consecutively increasing VI's indicating plant growth
    trend.growth <- lapply(trend.growth, function(x) {return(c(x, x[length(x)] + 1))}) # Extends the index to the last value of the trend
    st.growth <- sapply(trend.growth, min)
    valid.tg <- which(!is.na(st.growth))
    
    if(length(valid.tg)>1 & length(valid.ts)>1){
      trend.growth <- trend.growth[valid.tg]
      st.growth <- st.growth[valid.tg]
      en.growth <- sapply(trend.growth, max)
      
      trend.senescence <- trend.senescence[valid.ts]
      st.senescence <- st.senescence[valid.ts]
      en.senescence <- sapply(trend.senescence, max)
      rm(valid.ts, valid.tg)
    
      # Discard last growth series when it doesn't have a corresponding senescence series
      repeat{
        if(length(trend.growth)==0) break
        if (en.growth[length(trend.growth)] > en.senescence[length(trend.senescence)]) {
          trend.growth <- trend.growth[-length(trend.growth)]
          en.growth <- sapply(trend.growth, max)
        } else break
      }
      
      repeat{
        # Discard first senescence series when it doesn't have a corresponding growth series
        if(length(trend.senescence)==0) break
        if (st.growth[1] > st.senescence[1]) {
          trend.senescence <- trend.senescence[-1] 
          st.senescence <- sapply(trend.senescence, min)
        } else break
      }
      
      if (length(trend.senescence)>0 & length(trend.growth)>0) {
        #if(length(trend.senescence)!=length(trend.growth)) stop(cell,": Detected number of growth series do not match number of senescence series.")
        
        # Try to find interruptions in the curve, i.e. an observation of decrease 
        # after a supposed peak and then increased again to exhibit a second peak
        # This is rather easier done in reverse since the end-state (lowest point of senescence)
                
        # Get the VIs of the senescence series
        vi.senescence <- mapply("[", as.data.frame(vi), idx = trend.senescence, SIMPLIFY = FALSE)
        
        # Compute for the minimum VI per series
        min.senescence <- sapply(vi.senescence, min)
        #max.senescence <- sapply(vi.senescence, max)
        
        # Get the VIs of the growth series
        vi.growth <- mapply("[", as.data.frame(vi), idx = trend.growth, SIMPLIFY = FALSE)
        # Compute for the maximum VI per series
        max.growth <- sapply(vi.growth, max)
        # Compute for the minimum VI per series
        min.growth <- sapply(vi.growth, min)
        # Compute for the halfway-point VI per series
        md.growth <- min.growth+(max.growth-min.growth)/2 #sapply(vi.growth, median)
        
        # Incomplete crop cycle since the min.senescence of is lower than halfway-point 
        merge.potential <- which(md.growth <= min.senescence)
        
        if (length(merge.potential)>0){
          # Group merges with consecutive indexes
          merge.potential <- consecutive.groups(merge.potential)
          
          # TODO: Log disruptions disruption <- st.growth[merge.potential]
          merge.senesce <- mapply("[", list(trend.senescence), idx=merge.potential, SIMPLIFY = FALSE)
          merge.senesce <- lapply(merge.senesce, unlist)
          
          merge.min <- lapply(merge.potential,min)
          merge.max <- lapply(merge.potential,max)
          merge.max <- lapply(merge.max, "+", 1)
          merge.pot.growth <- mapply(c, as.list(merge.potential), merge.max, SIMPLIFY = FALSE)
          merge.growth <- mapply("[", list(trend.growth), idx=merge.pot.growth, SIMPLIFY = FALSE)
          merge.growth <- lapply(merge.growth, unlist)
          merged <- mapply(c, merge.growth, merge.senesce, SIMPLIFY = FALSE)
          merged <- lapply(merged, unique)
          merged <- lapply(merged, sort)
          
          trend.senescence <- trend.senescence[-unlist(merge.potential)]
          
          trend.growth[unlist(merge.min)] <- merged 
          trend.growth <- trend.growth[-unlist(mapply("[", merge.pot.growth, idx=list(-1)))]
          if(length(trend.growth)>length(trend.senescence)) trend.growth <- trend.growth[-length(trend.growth)]
        }
        
        len.senescence <- sapply(trend.senescence, length)
        len.growth <- sapply(trend.growth, length)
        
        # Remove those with growth trends shorter than minimum growth period and 
        # senescence trend shorter than minimum senescence period 
        crop.potential <- which(!((len.growth < growth.min) | (len.senescence < senescence.min)))
        
        if(length(crop.potential)>0){
          trend.growth <- trend.growth[crop.potential]
          trend.senescence <- trend.senescence[crop.potential]
          
          trend.crop <- mapply(c, trend.growth, trend.senescence, SIMPLIFY = FALSE)
          rm(trend.growth, trend.senescence)
          trend.crop <- lapply(trend.crop, unique)
          trend.crop <- lapply(trend.crop, sort)
          vi.crop <- mapply("[", as.data.frame(vi), idx = trend.crop, SIMPLIFY = FALSE)
          
          pos <- lapply(vi.crop, which.max)
          pos <- mapply("[", trend.crop, idx=pos, SIMPLIFY = FALSE)
          pos.evi <- lapply(vi.crop, max)
          evi.slope <- mapply(orysat::slope, x1=trend.crop, x2=pos, y1=vi.crop, y2=pos.evi, SIMPLIFY = FALSE)
          
          sos <- lapply(evi.slope, which.max)
          sos <- lapply(sos, "-", 2)
          sos <- lapply(sos, max, 1)
          sos.evi <- mapply("[", vi.crop, idx=sos)
          #sos.slope <- mapply("[", evi.slope, idx=sos)
          sos <- mapply("[", trend.crop, idx=sos)
          
          #FINDING EOS
          evi.amp <- mapply("-", pos.evi, sos.evi, SIMPLIFY = FALSE)
          evi.ampp5 <- lapply(evi.amp, "/",2)
          evi.mid <- mapply("+", sos.evi, evi.ampp5, SIMPLIFY = FALSE)
          
          filter.initevi <- mapply("<", vi.crop, evi.mid, SIMPLIFY = FALSE)
          filter.initevi <- lapply(filter.initevi, which)
          filtered.slope <- mapply("[", evi.slope, filter.initevi, SIMPLIFY = FALSE)
          
          eos <- lapply(filtered.slope, which.min)
          #eos.slope <- mapply("[", filtered.slope, idx=eos)
          eos <- mapply("[", filter.initevi, idx=eos)
          #eos.evi <- mapply("[", vi.crop, idx=eos)
          eos <- mapply("[", trend.crop, idx=eos)
          #crop.duration <- eos-sos+1

          actualcrop.idx <- mapply(":", as.list(sos), as.list(eos), SIMPLIFY = FALSE)
          vi.actualcrop  <- mapply("[", as.data.frame(vi), idx=actualcrop.idx, SIMPLIFY = FALSE)
          
          if(length(vi.actualcrop)>0) {
            # Save pixel data. Indicate Tile, Cell, and Start of season
            #if(class(sos.date)!="Date" && class(sos.date)=="numeric") sos.date <- as.Date(sos.date, "1970-1-1")
            mapply(saveRDS, vi.actualcrop, file=as.list(paste0(data.save, "/", format(dates[sos], "m%m/"), vi.label, ".", sprintf("h%02gv%02g.c%02i.st-%02g.en-%02g.rds", tile.h, tile.v, cell, sos, eos))))
            
            # for(i in 1:length(vi.actualcrop)){
            #   if(!is.null(dates)) sos.date <- dates[sos]
            #   if(class(sos.date)!="Date" && class(sos.date)=="numeric") sos.date <- as.Date(sos.date, "1970-1-1")
            #   
            #   #if(!dir.exists(paste0(data.save, "/",format(sos.date, "m%m")))) dir.create(paste0(data.save, "/",format(sos.date, "m%m")), recursive=TRUE)
            #   saveRDS(vi.actualcrop[[i]], file = paste0(data.save, "/", format(sos.date, "m%m/"), sprintf("h%02gv%02g.c%02i.st-%02g.en-%02g.rds", tile.h, tile.v, cell, sos[i], eos[i])))
            # }
            crops <- length(vi.actualcrop)
          }
          
          rm(trend.crop, sos.evi, evi.amp, evi.ampp5, evi.mid, filter.initevi, filtered.slope, vi.actualcrop, actualcrop.idx, eos, sos)
        }
        #if(!((length(st.growth)==1 && is.na(st.growth)) | (length(st.senescence)==1 && is.na(st.senescence)))){
      } 
    }
  }
  #gc(reset=TRUE)
  return(crops)
}