crop.extract <- function(vi, dates, vi.label="evi", interval=1, senescence.min=5, growth.min=6, seasonends.pctdiff=0.1, negligible.change=100, data.save=NULL, as.z = FALSE){
  # vi should include h and v tile value (vi[1]) and cell_id (vi[2]). This will be used when caching VI's of crop curve
  
  crops <- NA # Default output if nothing comes out from the analysis and filtering
  
  if(!is.na(vi[1]) & (sum(is.na(vi))!=length(vi))){
    tile.h <- trunc(vi[1])                # extract tile h-value
    tile.v <- round((vi[1]-tile.h) * 100) # extract tile v-value. round ensures Whole number
    cell <- vi[2]                         # extract cell ID  
    
    vi <- vi[-(1:2)]                      # Remove tile and cell values
    
    if(as.z) vi <- (vi-mean(vi))/sd(vi)   # Compute for z standard value
    
    dx <- diff(vi)  # Find out which ones are increasing and decreasing
    
    crops <- list()
    
    trend.senescence <- consecutive.groups(which(dx <= 0))                                     # Find consecutively decreasing VI's indicating senescence (sensence vectors)
    trend.senescence <- lapply(trend.senescence, function(x) {return(c(x, x[length(x)] + 1))}) # Extends the index to the last value of the trend
    st.senescence <- sapply(trend.senescence, min)                                             # Get the index of the first element of the senescence vectors
    valid.ts <- which(!is.na(st.senescence))                                                   # Find non-empty senescence vectors
    
    trend.growth <-consecutive.groups(which(dx > 0))                                           # Find consecutively increasing VI's indicating plant growth
    trend.growth <- lapply(trend.growth, function(x) {return(c(x, x[length(x)] + 1))})         # Extends the index to the last value of the trend
    st.growth <- sapply(trend.growth, min)                                                     # Get the index of the first element of the growth vectors
    valid.tg <- which(!is.na(st.growth))                                                       # Find non-empty growth vectors
    
    # if there are valid (non-empty) sensence and growth vectors do the following
    if(length(valid.tg)>1 & length(valid.ts)>1){ 
      trend.growth <- trend.growth[valid.tg]    # Remove empty growth vectors
      st.growth <- st.growth[valid.tg]          # Remove empty growth starts
      en.growth <- sapply(trend.growth, max)    # Extract index of last element of growth vectors
      
      trend.senescence <- trend.senescence[valid.ts]  # Remove empty senescence vectors
      st.senescence <- st.senescence[valid.ts]        # Remove empty senescence starts
      en.senescence <- sapply(trend.senescence, max)  # Extract index of last element of senescence vectors
      rm(valid.ts, valid.tg) # Clean up
      
      
      repeat{
        # Discard last growth series when it doesn't have a corresponding senescence series
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
        # Try to find interruptions in the curve, i.e. an observation of decrease 
        # after a supposed peak and then increased again to exhibit a second peak
        # I'm sure there's be a simpler way to do this, but it works so this is the way for now
        
        vi.growth <- mapply("[", as.data.frame(vi), idx = trend.growth, SIMPLIFY = FALSE)         # Get the VIs of the growth series
        max.growth <- sapply(vi.growth, max)                                                      # Compute for the maximum VI per series
        min.growth <- sapply(vi.growth, min)                                                      # Compute for the minimum VI per series
        md.growth <- min.growth+(max.growth-min.growth)/2                                         # Compute for the halfway-point VI per series
        
        vi.senescence <- mapply("[", as.data.frame(vi), idx = trend.senescence, SIMPLIFY = FALSE) # Get the VIs of the senescence series
        min.senescence <- sapply(vi.senescence, min)                                              # Compute for the minimum VI per series
        
        # Probably an incomplete crop cycle since the min.senescence is lower than halfway-point 
        merge.potential <- which(md.growth <= min.senescence)
        
        # if merge.potential is not empty, merge these potentials to growth vectors
        if (length(merge.potential)>0){
          # Try merging growth and senescence vectors which are potentially disruptions (due to some factor; e.g. weather, disease; that affected the signal)
          
          merge.potential <- consecutive.groups(merge.potential) # Group merge.potentials with consecutive indices 
          
          merge.senesce <- mapply("[", list(trend.senescence), idx=merge.potential, SIMPLIFY = FALSE) # get indices of senescence vectors using consecutive indices groups
          merge.senesce <- lapply(merge.senesce, unlist)         # Collapse each group into one vector
          
          merge.max <- lapply(merge.potential,max)  # Get the ending index of the merge vector
          merge.max <- lapply(merge.max, "+", 1)    # Include next index after merge.max 
          merge.pot.growth <- mapply(c, as.list(merge.potential), merge.max, SIMPLIFY = FALSE)    # Add merge.max to merge.potential
          merge.growth <- mapply("[", list(trend.growth), idx=merge.pot.growth, SIMPLIFY = FALSE) # Get indices of trend.growth 
          merge.growth <- lapply(merge.growth, unlist)                                            # Collapse merge.growth into 1 vector for each merge.growth group
          merged <- mapply(c, merge.growth, merge.senesce, SIMPLIFY = FALSE) # Combine merge.growth and merge.senscence
          merged <- lapply(merged, unique)                                   # Remove duplicate indices
          merged <- lapply(merged, sort)                                     # Sort the indices
          
          trend.senescence <- trend.senescence[-unlist(merge.potential)]     # Remove senescence vectors which are in merge.potential
          
          merge.min <- lapply(merge.potential,min)   # Get the starting index of the merge vector
          trend.growth[unlist(merge.min)] <- merged  # Append 'merged' to growth vectors
          trend.growth <- trend.growth[-unlist(mapply("[", merge.pot.growth, idx=list(-1)))] # Remove growth vectors which were part of the merge
          
          # if there are more growth vector than senescence vectors, remove the last growth vector
          if(length(trend.growth)>length(trend.senescence)) trend.growth <- trend.growth[-length(trend.growth)] 
        }
        
        len.senescence <- sapply(trend.senescence, length) # Compute for the length of the senescence vectors
        len.growth <- sapply(trend.growth, length)         # Compute for the length of the growth vectors
        
        # Remove those with growth trends shorter than minimum growth period and 
        # senescence trend shorter than minimum senescence period 
        crop.potential <- which(!((len.growth < growth.min) | (len.senescence < senescence.min)))
        
        # If crop.potential is not empty, compute vi metrics and save into cache
        if(length(crop.potential)>0){
          trend.growth <- trend.growth[crop.potential]            # Retain growth vectors which are potential crops
          trend.senescence <- trend.senescence[crop.potential]    # Retain senescence vectors which are potential crops
          
          trend.crop <- mapply(c, trend.growth, trend.senescence, SIMPLIFY = FALSE) # combine growth and senescence vectors into one crop vector
          rm(trend.growth, trend.senescence) # Clean up
          trend.crop <- lapply(trend.crop, unique) # Remove duplicate indices per crop vector
          trend.crop <- lapply(trend.crop, sort)   # Sort crop vectors to ensure VI's will be chronological
          vi.crop <- mapply("[", as.data.frame(vi), idx = trend.crop, SIMPLIFY = FALSE) # Extract VI values using crop vectors (indices)
          
          # Peak of Season (POS)- When the maximum VI was observed
          pos <- lapply(vi.crop, which.max) 
          pos <- mapply("[", trend.crop, idx=pos, SIMPLIFY = FALSE)
          
          # pos.evi = Maximum VI value in the crop vector
          pos.evi <- lapply(vi.crop, max)
          
          # Compute slope relative to the Maximum VI
          evi.slope <- mapply(orysat::slope, x1=trend.crop, x2=pos, y1=vi.crop, y2=pos.evi, SIMPLIFY = FALSE)
          
          # Start of season (SOS) = 2 acquisition day points before maximum slope (experimental)
          # For rice, SOS is simply the acquisition day point where agronomic flooding is detected
          # This is done with the idea that agronomic flooding should be followed by a fast plant growth 
          sos <- lapply(evi.slope, which.max)
          sos <- lapply(sos, "-", 2)
          sos <- lapply(sos, max, 1)
          sos <- mapply("[", trend.crop, idx=sos)
          
          # sos.evi =  VI value at the assumed start of the season
          sos.evi <- mapply("[", vi.crop, idx=sos)
          
          # VI amplitude (Difference of Peak VI and SOS VI)
          evi.amp <- mapply("-", pos.evi, sos.evi, SIMPLIFY = FALSE)
          # Half of the amplitude
          evi.ampp5 <- lapply(evi.amp, "/",2)
          # VI mid is the threshold that needs to be met, also to indicate fast vegetation growth, 
          # which should be attained after agronomic flooding, (One criterion to be considered rice)
          # Also used to estimate End of Season
          evi.mid <- mapply("+", sos.evi, evi.ampp5, SIMPLIFY = FALSE)
          
          # Logically one could assume that the VI at the end of season should be more or less the same as
          # VI at the start of the season, but this is not the case. Hence, the criterion for eos VI is to 
          # at least pass evi.mid
          filter.initevi <- mapply("<", vi.crop, evi.mid, SIMPLIFY = FALSE) # find VI's less than evi.mid
          filter.initevi <- lapply(filter.initevi, which)                   # Get the index 
          filtered.slope <- mapply("[", evi.slope, filter.initevi, SIMPLIFY = FALSE) # From the VI's that passed evi.mid, get the slope
          eos <- lapply(filtered.slope, which.min)     # Index of the steepest slope (min because slopes are negative)
          eos <- mapply("[", filter.initevi, idx=eos)  # Index wrt filter.initevi (VI below evi.mid) with the steepest slope
          eos <- mapply("[", trend.crop, idx=eos)      # Index wrt crop vector
          
          actualcrop.idx <- mapply(":", as.list(sos), as.list(eos), SIMPLIFY = FALSE) # sequence of indices from sos to eos
          vi.actualcrop  <- mapply("[", as.data.frame(vi), idx=actualcrop.idx, SIMPLIFY = FALSE) # trimmed VI's from sos to eos
          
          # If vi.actual crop is not empty
          if(length(vi.actualcrop)>0) {
            # Save pixel data. Indicate Tile, Cell, and Start of season
            if (!is.null(data.save)) { 
              # save to disk
              mapply(saveRDS, vi.actualcrop, file=as.list(paste0(data.save, "/", format(dates[sos], "m%m/"), vi.label, ".", sprintf("h%02gv%02g.c%02i.st-%02g.en-%02g.rds", tile.h, tile.v, cell, sos, eos))))
            } 
            names(vi.actualcrop) <- sprintf("st-%02g.en-%02g", sos, eos) # Indicate sos and eos in names of list
            # output to memory
            crops <- vi.actualcrop
          }
          # Clean up
          rm(trend.crop, sos.evi, evi.amp, evi.ampp5, evi.mid, filter.initevi, filtered.slope, vi.actualcrop, actualcrop.idx, eos, sos)
          crops <- crops[tmp.curvlen>1] # Remove empty crop vectors
          
          # Compute metrics 
          crop.duration <- lapply(crops, length)
          sos.vi        <- lapply(crops, "[", 1) # Get vi at start of season (index 1) for each crop vi vector
          mid1.vi       <- lapply(crops, "[", floor(pos/2))  # Get vi at mid-growth of season (index 1) for each crop vi vector
          pos.vi        <- lapply(crops, max)                # Get vi at start of season (index 1) for each crop vi vector
          mid2.vi       <- lapply(crops, "[", pos+floor((pos-eos)/2)) # Get vi at mid-harvest of season (index 1) for each crop vi vector
          eos.vi        <- lapply(crops, "[", crop.duration) # Get vi at end of season (index length of vi vector) for each crop vi vector
          pos           <- lapply(crops, which.max)
          auc           <- lapply(crops, curve.Area)
          up.slope      <- orysat::slope(x1=1, y1=sos.vi, x2=pos, y2=pos.vi)
          down.slope    <- orysat::slope(x1=pos, y1=pos.vi, x2=eos, y2=eos.vi)
          q1.slope      <- orysat::slope(x1=1, y1=sos.vi, x2=floor(pos/2), y2=pos.vi)
          q2.slope      <- orysat::slope(x1=floor(pos/2), y1=sos.vi, x2=floor(pos/2), y2=pos.vi)
          q3.slope      <- orysat::slope(x1=pos, y1=pos.vi, x2=pos+floor((pos-eos)/2), y2=mid2.vi)
          q4.slope      <- orysat::slope(x1=pos+floor((pos-eos)/2), y1=mid2.vi, x2=floor(pos/2), y2=eos.vi)
          vec.crops <- c(sos, pos, eos, sos.vi, xiao5.vi, mid1.vi, pos.vi, mid2.vi, eos.vi, auc, up.slope, down.slope, q1.slope, q2.slope, q3.slope, q4.slope)
        }
      } 
    }
  } 
  #gc(reset=TRUE)
  return(vec.crops)
}
