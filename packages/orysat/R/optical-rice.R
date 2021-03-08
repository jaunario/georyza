# Title: Rice detection using vegetation indices from Optical Sensors
# Description: The inputs of the functions are a time series of vegetation indices (EVI, NDVI) and water indices (LSWI, NDWI, MNDWI)
# Author: JAunario
###############################################################################

pheno.cropstart <- function(){
  
}

# This is a pixel level implementation of the original Rice Mapping algorithm used in the IRRI GSM (formerly GIS) Laboratory
# It uses EVI, NDVI and LSWI to first detect flooding (assumed start of season), then find an increase
# in the vegetation index signaling growth of the supposed rice plant. 
# This implementation aims to reduce disk space usage and processing time
# Input is a time series of the mentioned indices above 46 in the year of interest + 15 from previous
# Returned value are as follows:
#   rnr (rice-nonrice) = TRUE
#     0 - Non-Rice
#     1 - Rice
#
#   rnr = FALSE
#     0 - Fallow
#     100+idx - Single 
#     20000+idx2*100+idx1 - 2 Seasons of Rice
#     3000000+idx3*10000+idx2*100+idx1 - 3 Seasons of Rice


rice.Xiao_v1 <- function(evi, lswi, ndvi=NULL, evi.floodmax=0.25, evi.ricemax=-1, evi.halfricemax=-1, rnr.only=FALSE, data.interval=8, crop.duration=120){
  
  rice <- 0
  if(is.null(ndvi)) ndvi <- evi
  if(length(evi)!=length(lswi)|length(lswi)!=length(ndvi)) stop("Unequal number of data points")
  #evi.ricemax=-1
  #evi.halfricemax=-1
  evi.ricemax <- ifelse(evi.ricemax==0, 0.75, evi.ricemax) 
  evi.halfricemax <- ifelse(evi.halfricemax==0, 0.4, evi.halfricemax)
  
  pts.40d <- ceiling(40/data.interval) # No of data points to reach 40 days
  pts.cd <- ceiling(crop.duration/data.interval) # No of data points to reach crop duration
  
  # Expected latest possible start of season based on crop duration specified
  flooding.last <- length(evi)-pts.cd
  # Determine rice and non-rice properties, if necessary, on the pixel data
  
  # Find flooding. Normal lowland rice fields are initially flooded
  flooding <- which(flooded1(lswi=lswi[1:flooding.last], ndvi=ndvi[1:flooding.last], evi=evi[1:flooding.last]) &
                      evi[1:flooding.last] <= evi.floodmax & evi[1:flooding.last] >= -1)
  evi.flood <- evi[flooding]

  if (length(flooding)>0){
    # Get evi from start of flooding upto crop duration
    evi.cropdur <- matrix(evi[sapply(flooding, FUN=seq, length.out=pts.cd)],ncol=pts.cd, byrow=TRUE)
    #Compute for evi.rice max if < 0
    if(evi.ricemax < 0) evi.ricemax <- apply(evi.cropdur, 1, max, na.rm = TRUE)
    if(evi.halfricemax < 0)  evi.halfricemax <- evi.ricemax/2
    
    # find flooding if no increase within cropping duration
    to.retain <- which(evi.halfricemax>evi.flood & rowSums(evi.cropdur>evi.flood, na.rm=TRUE)>0)
    
    # If there is still valid flooding
    if(length(to.retain)>0){
      flooding <- flooding[to.retain]
      evi.ricemax <- evi.ricemax[to.retain]
      evi.halfricemax <- evi.halfricemax[to.retain]
      evi.cropdur <- matrix(evi.cropdur[to.retain,], nrow = length(to.retain))
      
      # Create matrix of EVI's 40 days after the flooding 
      evi.fldp40 <- matrix(evi.cropdur[,(1:pts.40d)+1], nrow = nrow(evi.cropdur))
      rice.potential <- which(rowSums(evi.fldp40>=evi.halfricemax, na.rm=TRUE)>1)
      
      # TODO : add option if too many rice check evi.ricemax should be after 40days
      # evi.after40 <- matrix(evi.cropdur[,(pts.40d+1):ncol(evi.cropdur)], nrow = nrow(evi.cropdur))
      
      # If there is rice.potential
      if(length(rice.potential)>0 & rnr.only){
        rice <- 1
      } else if(length(rice.potential)>0){
        if(length(rice.potential)>1){
          
          flooding <- flooding[rice.potential]
          
          # remove consecutive, use latest flooding
          flooding <- sapply(consecutive.groups(flooding), max)
          
          if(length(flooding)>1){
            #TODO: check this
            difs <- sapply(flooding, "-", flooding)
            difs[difs<0] <- NA 
            difs <- difs-pts.cd
            difs[difs == -15] <- NA
            flooding <- flooding[which(colSums(difs<0, na.rm=TRUE)==0)]
          } 
          
          rice <- sum(10^((1:length(flooding)-1)*3)*flooding+(10^(c(2,5,8)[1:length(flooding)])*(1:length(flooding))))
        }
      }
      if(length(flooding)==0) rice <- 0
    }
  }
  return(rice)
}

# TODO: My interpretation of the Xiao method, with some improvements
# rice.Xiao2006 <- function(evi, lswi, ndvi){
#   
#   if(length(evi)!=length(lswi)|length(lswi)!=length(ndvi)) stop("Unequal number of data points")
#   
#   # Expected latest possible start of season based on crop duration specified
#   flooding.last <- length(evi)-pts.cd
#   
#   # Find flooding. Normal lowland rice fields are initially flooded
#   flooding <- which(flooded1(lswi=lswi[1:flooding.last], ndvi=ndvi[1:flooding.last], evi=evi[1:flooding.last]) &
#                       evi[1:flooding.last] <= evi.floodmax & evi[1:flooding.last] >= -1)
#   
#   evi.flood <- evi[flooding]
#   
#   if (length(flooding)>0){
#     # Get evi from start of flooding upto crop duration
#     evi.cropdur <- matrix(evi[sapply(flooding, FUN=seq, length.out=pts.cd)],ncol=pts.cd, byrow=TRUE)
#     #Compute for evi.rice max if < 0
#     if(evi.ricemax < 0) evi.ricemax <- apply(evi.cropdur, 1, max, na.rm = TRUE)
#     if(evi.halfricemax < 0)  evi.halfricemax <- evi.ricemax/2
#     
#     # find flooding if no increase within cropping duration
#     to.retain <- which(evi.halfricemax>evi.flood & rowSums(evi.cropdur>evi.flood, na.rm=TRUE)>0)
#     
#     # If there is still valid flooding
#     if(length(to.retain)>0){
#       flooding <- flooding[to.retain]
#       evi.ricemax <- evi.ricemax[to.retain]
#       evi.halfricemax <- evi.halfricemax[to.retain]
#       evi.cropdur <- matrix(evi.cropdur[to.retain,], nrow = length(to.retain))
#       
#       
#       
#       # If there is rice.potential
#       if(length(rice.potential)>0 & rnr.only){
#         rice <- 1
#       } else if(length(rice.potential)>0){
#         if(length(rice.potential)>1){
#           
#           flooding <- flooding[rice.potential]
#           
#           # remove consecutive, use latest flooding
#           flooding <- sapply(consecutive.groups(flooding), max)
#           
#           if(length(flooding)>1){
#             #TODO: check this
#             difs <- sapply(flooding, "-", flooding)
#             difs[difs<0] <- NA 
#             difs <- difs-pts.cd
#             difs[difs == -15] <- NA
#             flooding <- flooding[which(colSums(difs<0, na.rm=TRUE)==0)]
#           } 
#           
#           rice <- sum(10^((1:length(flooding)-1)*3)*flooding+(10^(c(2,5,8)[1:length(flooding)])*(1:length(flooding))))
#         }
#       }
#       if(length(flooding)==0) rice <- 0
#     }
#   }
#   return(rice)
# }


rice.modxiao <- function(ts.vi, ts.flood=NULL, analysis.fun=xiaoflags.rice, data.interval=8, crop.duration=120, rnr.only=FALSE, ...){
# Adjustable Xiao-based rice detection algorithm.
# Input any vegetation index, any flood vector
# 
  pts.cd <- ceiling(crop.duration/data.interval) # No. of data points to reach crop duration
  
  # Expected latest possible start of season based on crop duration specified
  flooding.last <- length(ts.vi)-pts.cd
  # Determine rice and non-rice properties, if necessary, on the pixel data
  
  # Find flooding. Normal lowland rice fields are initially flooded
  flooding <- which(ts.flood[1:flooding.last])
  
  if (length(flooding)>0){
    
    flooding <- sapply(consecutive.groups(flooding), max) # For 
    #evi.flood <- ts.vi[flooding]
    
    # Get evi from start of flooding upto crop duration
    evi.cropdur <- matrix(ts.vi[sapply(flooding, FUN=seq, length.out=pts.cd)],ncol=pts.cd, byrow=TRUE)

    with.snow <- which(rowSums(evi.cropdur < -1, na.rm=TRUE)>0)
    if(length(with.snow)>0){
      evi.cropdur <- evi.cropdur[-with.snow,]
    }
    if(is.null(dim(evi.cropdur))){
      if(length(evi.cropdur)>0)rice.potential <- which(analysis.fun(evi.cropdur)) else rice.potential <- vector()
    } else if(nrow(evi.cropdur)>0){
      rice.potential <- which(apply(evi.cropdur, 1, analysis.fun))
    } else rice.potential <- vector()
    
    # If there is rice.potential
    if(length(rice.potential)>0){
      if (rnr.only) {
        rice <- 1
      } else {
        flooding <- flooding[rice.potential]
        if(length(flooding)>1){
          # remove consecutive, use latest flooding
          true.rice <- vector()
          repeat{
            this.flood <- flooding[1]
            if(length(true.rice)==0) true.rice <- this.flood else true.rice <- c(true.rice, this.flood)
            chk.overlap <- flooding-this.flood-15
            flooding <- flooding[chk.overlap>0]
            if(length(flooding)==0) break
          }
        } else true.rice <- flooding
        rice <- sum(10^((1:length(true.rice)-1)*3)*true.rice+(10^(c(2,5,8)[1:length(true.rice)])*(1:length(true.rice))))
      }
    } else {
      rice <- 0
    }
  } else rice <- 0   
  
  
  # rice <- 0
  # if (length(flooding)>0){
  #   
  #   flooding <- sapply(consecutive.groups(flooding), max) # For 
  #   #evi.flood <- ts.vi[flooding]
  #   
  #   # Get evi from start of flooding upto crop duration
  #   evi.cropdur <- matrix(ts.vi[sapply(flooding, FUN=seq, length.out=pts.cd)],ncol=pts.cd, byrow=TRUE)
  #   #Compute for evi.rice max if < 0
  #   rice.potential <- which(apply(evi.cropdur, 1, analysis.fun))
  #   
  #   # If there is rice.potential
  #   if(length(rice.potential)>0 & rnr.only){
  #       rice <- 1
  #   } else if(length(rice.potential)>0){
  #     flooding <- flooding[rice.potential]
  #     if(length(flooding)>1){
  #       # remove consecutive, use latest flooding
  #       true.rice <- vector()
  #       repeat{
  #         this.flood <- flooding[1]
  #         if(length(true.rice)==0) true.rice <- this.flood else true.rice <- c(true.rice, this.flood)
  #         chk.overlap <- flooding-this.flood-15
  #         flooding <- flooding[chk.overlap>0]
  #         if(length(flooding)==0) break
  #       }
  #     } else true.rice <- flooding
  #     rice <- sum(10^((1:length(true.rice)-1)*3)*true.rice+(10^(c(2,5,8)[1:length(true.rice)])*(1:length(true.rice))))
  #   }
  # }   
  return(rice)
}

#TODO: Could be good to convert into a more generic function that can detect "hills" 
#      in a time series
crop.signature <- function(vi, dates=null, interval=1, senescence.min=5, growth.min=6, seasonends.pctdiff=0.1, negligible.change=100){
  
  # stdvi <- (vi-mean(vi))/sd(vi)
  dx <- diff(vi) # Find out which ones are increasing and decreasing
  # dx <- diff(stdvi) # Find out which ones are increasing and decreasing
  trend.senescence <- consecutive.groups(which(dx <= 0)) # Find consecutively decreasing VI's indicating senescence
  trend.senescence <- lapply(trend.senescence, function(x) {return(c(x, x[length(x)] + 1))}) # Extends the index to the last value of the trend
  #len.senescence <- sapply(trend.senescence, length)
  st.senescence <- sapply(trend.senescence, min)
  en.senescence <- sapply(trend.senescence, max)
  
  trend.growth <-consecutive.groups(which(dx > 0)) # Find consecutively increasing VI's indicating plant growth
  trend.growth <- lapply(trend.growth, function(x) {return(c(x, x[length(x)] + 1))}) # Extends the index to the last value of the trend
  #len.growth <- sapply(trend.growth, length)
  st.growth <- sapply(trend.growth, min)
  en.growth <- sapply(trend.growth, max)
  
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
  crops <- data.frame()
  
  if (length(trend.senescence)>0 & length(trend.growth)>0) {
    if(length(trend.senescence)!=length(trend.growth)) stop()
    
    # Try to find interruptions in the curve, i.e. an observation of decrease 
    # after a supposed peak and then increased again to exhibit a second peak
    # This is rather easier done in reverse since the end-state (lowest point of senescence)
    # 
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
      trend.crop <- lapply(trend.crop, unique)
      trend.crop <- lapply(trend.crop, sort)
      vi.crop <- mapply("[", as.data.frame(vi), idx = trend.crop, SIMPLIFY = FALSE)
      
      pos <- lapply(vi.crop, which.max)
      pos <- mapply("[", trend.crop, idx=pos, SIMPLIFY = FALSE)
      pos.evi <- lapply(vi.crop, max)
      evi.slope <- mapply(slope, x1=trend.crop, x2=pos, y1=vi.crop, y2=pos.evi, SIMPLIFY = FALSE)
      
      sos <- lapply(evi.slope, which.max)
      sos <- lapply(sos, "-", 2)
      sos <- lapply(sos, max, 1)
      sos.evi <- mapply("[", vi.crop, idx=sos)
      sos.slope <- mapply("[", evi.slope, idx=sos)
      sos <- mapply("[", trend.crop, idx=sos)
      
      evi.amp <- mapply("-", pos.evi, sos.evi, SIMPLIFY = FALSE)
      evi.ampp5 <- lapply(evi.amp, "/",2)
      evi.mid <- mapply("+", sos.evi, evi.ampp5, SIMPLIFY = FALSE)
      
      filter.initevi <- mapply("<", vi.crop, evi.mid, SIMPLIFY = FALSE)
      filter.initevi <- lapply(filter.initevi, which)
      filtered.slope <- mapply("[", evi.slope, filter.initevi, SIMPLIFY = FALSE)
      
      
      eos <- lapply(filtered.slope, which.min)
      eos.slope <- mapply("[", filtered.slope, idx=eos)
      eos <- mapply("[", filter.initevi, idx=eos)
      eos.evi <- mapply("[", vi.crop, idx=eos)
      eos <- mapply("[", trend.crop, idx=eos)
      crop.duration <- eos-sos+1
      
      actualcrop.idx <- mapply(":", as.list(sos), as.list(eos), SIMPLIFY = FALSE)
      vi.actualcrop <- mapply("[", as.data.frame(vi), idx=actualcrop.idx, SIMPLIFY = FALSE)
      meanVI <- sapply(vi.actualcrop, mean)
      sdVI <- sapply(vi.actualcrop, sd)
      aucVI <- mapply(curve.Area, y=vi.actualcrop, x=actualcrop.idx)
      crops <- data.frame(sos, pos=unlist(pos), eos, sos.evi, pos.evi=unlist(pos.evi), eos.evi, sos.slope, eos.slope, meanVI, sdVI, aucVI, crop.duration, row.names = NULL)
      
      crops <- crops[crops$pos-crops$sos>=growth.min & crops$eos-crops$pos>=senescence.min,  ] 
    }
  }
  
  return(crops)
  # len.cp <- lapply(cropping.potential, length)
  # len.cp < growth.min + senescence.min
  # Find point of harvest where EVI decrease is maximum (Probably should be )  
  
  # Find min.growth closest to min.senescence
  # min.diff <- mapply("-", as.list(min.senescence), as.data.frame(min.growth), SIMPLIFY = FALSE)
  # min.diff <- lapply(min.diff, abs)
  # 
  # max.diff <- mapply("-", as.list(max.senescence), as.data.frame(max.growth), SIMPLIFY = FALSE)
  # max.diff <- lapply(max.diff, abs)
  # 
  # close.diff <- lapply(min.diff, which.min)
  # 
  # min.senescence<(max.growth/2)
  # disrupt.idx <- which((len.growth<growth.min & len.senescence>=senescence.min) | (len.growth>=growth.min & len.senescence<senescence.min))
  # Growth duration that is less than growth.min probably is not a crop unless 
  # the following senescence is only growth disruption or it is a part of longer duration crop
  # noncrop.growth <- which(len.growth<growth.min & len.senescence>=senescence.min)
  
  # min.growth >= min.senescence & len.growth >=
  
  # Potential disruption in growth if senescence duration (len) is less than senescence.min
  # disrupt.growth <- which(len.senescence<senescence.min)
  # acceptable.senescence
  
   
  
  # Merge noise desc to corresponding asc
  # trend.growth <- as.list(rep(NA, length))
  # if(sum(len.senescence<senescence.min)>0){
  #   trend.growth <- mapply(c, trend.growth[len.senescence<=noise.length], trend.senescence[len.senescence<=noise.length])
  # }
  # 
  #len.senescence <- sapply(trend.senescence, length) # Duration of decreasing VI's
  #len.growth <- sapply(trend.growth, length) # Duration of increasing VI's
  

    
  # pct.diff <- mapply("-", as.list(min.senescence), as.data.frame(min.growth), SIMPLIFY = FALSE)
  # pct.diff <- mapply("/", pct.diff, as.data.frame(min.growth), SIMPLIFY = FALSE)
  # pct.diff <- lapply(pct.diff, abs)
  # pct.diff <- lapply(pct.diff, "<=", 0.2)
  #(min.senescence[1]-min.growth)/min.growth
  # flr.growth <- min.growth-(max.growth-min.growth)*3/4
  # len.growth <- sapply(vi.growth, length)
  
  #len.growth <- sapply(vi.growth, length)
  # Try to combine ASC and DESC
  # If max of Gn is < Gn+1 and Length < growth.min
  #dx2 <- diff(max.growth)
  #dx2.idx <- which(c(FALSE,dx2<0) & len.growth<growth.min)
  
  #TODO: Probably good to compute other stats here
  # Check if min vi of senescence within range of growth
  #lvl1 <- lapply(pct.diff, "<=", 0.65)#pct.diff < 0.2 
  # lvl1 <- mapply("<", as.list(min.senescence), as.data.frame(md.growth)) & mapply(">", as.list(min.senescence), as.data.frame(flr.growth))  
  #lvl1 <- lvl1 & len.growth>=growth.min
  # Get index of lvl1=TRUE
  # lvl1.match <- sapply(as.data.frame(lvl1), which, simplify = FALSE)
  # Remove matches where potential match occur after trend.senescence group
  # lvl1.filter <- mapply("<=", lvl1.match, as.list(1:length(trend.senescence)), SIMPLIFY = FALSE)
  # lvl1.match <- mapply("[", lvl1.match, lvl1.filter, SIMPLIFY = FALSE)
  # 
  # len.match <- sapply(lvl1.match, length)
  # 
  # if(sum(len.match)>0){
  #   dt.pts <- mapply("[", as.data.frame(len.growth), idx=lvl1.match[len.match > 0], SIMPLIFY = FALSE)
  #   dt.pts <- lapply(dt.pts, rev)
  #   dt.pts <- lapply(dt.pts, cumsum)
  #   dt.pts <- mapply(">=", dt.pts, as.list(growth.min), SIMPLIFY = FALSE)
  #   dt.pts <- lapply(dt.pts, rev)
  #   dt.pts <- lapply(dt.pts, which)
  #   dt.pts.len <- lapply(dt.pts, length)
  #   dt.pts[dt.pts.len>0] <- lapply(dt.pts[dt.pts.len>0], max) 
  #   
  #   lvl1.match[len.match > 0] <- mapply("[", lvl1.match[len.match > 0], dt.pts)
  #   
  #   len.match <- sapply(lvl1.match, length)
  # }
  # lvl2.filter <- mapply("<", lvl1.match, as.list(1:length(trend.senescence)))
  # lvl1.ok <- mapply("%in%", lvl1.match, lvl1.match)
  # lvl1.inc <- mapply("=", lvl1.match, lvl1.ok)
  #
  # Check for multiple and overlapping matches
  
  # if(sum(len.match)>0){
  #   grps.tocombine <-  mapply(":" , lvl1.match[len.match > 0], as.list((1:length(trend.senescence))[len.match > 0]), SIMPLIFY = FALSE)
  #   asc.idx <- mapply("[", list(trend.growth), idx = grps.tocombine, SIMPLIFY = FALSE)
  #   asc.idx <- lapply(asc.idx, unlist)
  #   #asc.len <- sapply(asc.idx, length)
  #   
  #   desc.idx <- mapply("[", list(trend.senescence), idx = grps.tocombine, SIMPLIFY = FALSE)
  #   desc.idx <- lapply(desc.idx, unlist)
  #   #desc.len <- sapply(desc.idx, length)
  #   
  #   #asc.idx <- asc.idx[asc.len >= growth.min & desc.len >= senescence.min]
  #   #desc.idx <- desc.idx[asc.len >= growth.min & desc.len >= senescence.min]
  #   
  #   crop.idx <- mapply(c, asc.idx, desc.idx, SIMPLIFY = FALSE)
  #   crop.idx <- lapply(crop.idx, sort)
  #   crop.idx <- lapply(crop.idx, unique)
  #   
  #   # Check if cropping long enough 
  #   crop.len <- sapply(crop.idx, length)
  #   crop.idx <- crop.idx[crop.len>=(growth.min+senescence.min)]
  #   if(length(crop.idx)>0){
  #     vi.crop <-  mapply("[", as.data.frame(vi), idx = crop.idx, SIMPLIFY = FALSE)
  #     meanVI <- sapply(vi.crop, mean)
  #     sdVI <- sapply(vi.crop, sd)
  #     
  #     stdVI.crop <- mapply("-", vi.crop, as.list(meanVI), SIMPLIFY = FALSE) 
  #     stdVI.crop <- mapply("/", vi.crop, as.list(sdVI), SIMPLIFY = FALSE)
  #     
  #     crop.sigchange <- lapply(stdVI.crop, ">", 0.5)
  #     crop.sigchange <- lapply(crop.sigchange, which)
  #     sigchange.len <- sapply(crop.sigchange, length)
  #     
  #   } else sigchange.len <- 0
  #   
  #   if(sum(sigchange.len)>0){
  #     crop.idx <- crop.idx[sigchange.len>2]
  #     vi.crop <- vi.crop[sigchange.len>2]
  #     crop.start <- lapply(crop.sigchange[sigchange.len>2], min)
  #     crop.end <- lapply(crop.sigchange[sigchange.len>2], max)
  #     
  #     sos <- mapply("[", crop.idx, idx=crop.start)
  #     sosVI <- mapply("[", vi.crop, idx=crop.start)
  #     eos <- mapply("[", crop.idx, idx=crop.end)
  #     eosVI <- mapply("[", vi.crop, idx=crop.end)
  #     peakVI <- sapply(vi.crop, max)
  #     peak <- lapply(vi.crop, which.max)
  #     peak <- mapply("[", crop.idx, idx=peak)
  #     
  #     retain <- (peak-sos)>=growth.min & (eos-peak)>=senescence.min
      # crops <- crops[retain, ]
      # crops <- crops[!duplicated(crops$peak),]
      # crops <- crops[crops$crop.duration>=(growth.min+senescence.min),] 
                     #& crops$sdVI>=200,]      
  
  
  # vidiff.grsns <- mapply("-", as.list(min.senescence), as.data.frame(min.senescence), SIMPLIFY = FALSE)
  
  # find compatible min
  # vidiff.grsns <- mapply("-", as.list(min.senescence), as.data.frame(min.senescence), SIMPLIFY = FALSE)
  # vidiff.grsns <- sapply(vidiff.grsns, simplify = FALSE)
  
  # len.growth >= growth.min 
  # Remove len.senescence < senesce.min & len.senescence > senesce.max (senescence period should not be too short or too long)
  # valid.desc <- trend.senescence[len.senescence>=senescence.min & len.senescence<=senescence.max]
  
  # Algorithm on whether a crop starts and ends on a "non-crop state" 
  # where VI's are relatively close assuming that there is 
  # no plant or very little vegetation is detected on the ground
  
  # Compute for percent abs difference of the min VI's WRT min of growth period (vi.growth)
  # lvl1.test <- mapply("<",vi.senescence, min.growth+min.growth*seasonends.pctdiff) 
  #& min.growth+min.growth*seasonends.pctdiff
  
  # crop.idx <- mapply(c, trend.growth, trend.senescence)
  # crop.idx <- sapply(crop.idx, unique)
  
  # Investigate crop duration 
  # crop.duration <- sapply(crop.idx, length)
  
  
  # UNEXPECTED BREAKS IN INCREASE 
  # brk.asc <- which(len.senescence<=break.max)
  # brk.asc <- brk.asc[!brk.asc %in% c(1,length(trend.senescence))]
  # Break should not start at 1 or end in length(vi)

  
  # LEVEL 1 CROPPING CHECK start and end VI's should be below specified seasonends.pctdiff
  # lvl1.idx <- mapply(c, trend.growth[lvl1.test<=seasonends.pctdiff], trend.senescence[lvl1.test<=seasonends.pctdiff])

  # If minima adjacent crop periods is at least higher than seasonends.pctdiff, 
  # they could be 1 crop with several harvest like fruit or a disease could have affected the crop 
  # but was able to recover 

  # attempt to combine apparent incomplete crop seasons
  # to.combine <- which(lvl1.test > seasonends.pctdiff)
  #if(length(to.combine) >2) stop("To many crops to combine")
  # lvl2.idx <-  sort(do.call(c, c(trend.senescence[to.combine], trend.growth[to.combine])))
  
  # crops.idx <- 
  #
  
  # vi.crop <- mapply(c, vi.growth,vi.senescence)
  # mean.crop <- sapply(vi.crop, mean)
  # diff.crop <- sapply(vi.crop, diff)
  # absdiff <- sapply(diff.crop, abs)
  # Use quantile of absolute diffs to remove possible plateauing 
  
  
  # Investigate VI values
  # Get VI at 
  #   crop.start [i=1 from all trend.growth list items], 
  #   crop.end  [i=which.min dx[tred.desc]]
  #   crop.peak [i=1 from all trend.senescence list items]

  # VI at crop.end should be "near" VI crop.start
  # VI at crop.peak should be significantly different from VI at crop.start and VI at crop.end
  
  
  # MAX VI ( get first vi from each remaining trend.senescence)
  # TODO: Determine the following: 
  # CROP START
  # CROP END
  # GROWTH RATE
  # SENESCE RATE
  
}

rice.VWIdynamics <- function(vi, wi, flowering.period=7, latest.planting=50, max.evi = 7000, ripening.period=8){
  result <- vector()
  #steps
  # find maxima
  # Should be consistently increasing from x1(DEFAULT=13 to complete 100 days) acqdates before max evi, and consistently decrease x2(DEFAULT=10 maximum) acqdates after max evi)
  # test poly reg to get rsquare
  
  trend.growth <- consecutive.groups(which(dx>0))
  len.growth <- sapply(trend.growth, length)
  min.desc <- sapply(trend.growth, min)
  
  min.asc <- sapply(trend.growth, min)
  enoughtime_toflower <- len.growth>=flowering.period # Period to flowering (max.evi) is constantly increasing trend for at least 56 days (7 ACQDATES, more lenient than 8)
  within_timeframe <- min.asc<=latest.planting # Limits the analysis based on latest planting date
  pot.rice <- which(within_timeframe & enoughtime_toflower)
  
  if(length(pot.rice)>0){
    trend.growth <- trend.growth[pot.rice]    # Limit analysis to segments that passes potentiality
    min.asc <- min.asc[pot.rice]
    max.asc <- sapply(trend.growth, max)+1 # get the point of flowering EVI
    len.growth <- len.growth[pot.rice]
    maturity.period <- sapply(max.asc, seq, length.out=ripening.period)
    
    vi.ripening <- matrix(vi[maturity.period], ncol=ripening.period, byrow = TRUE)
    dx.ripening <- apply(vi.ripening, 1, diff) # Compute Decrements
    #trend.senescence <- apply(vi.ripening, 1, diff) # Compute Decrements
    cons.dec5 <- colSums(matrix(dx.ripening[1:5,], nrow=5)<0, na.rm = TRUE) # 
    is.rice <- cons.dec5>=4
    if(sum(is.rice)>0){
      trend.growth <- trend.growth[is.rice]
      min.asc <- min.asc[is.rice]
      max.asc <- max.asc[is.rice]
      
      len.growth <- len.growth[is.rice]
      dx.ripening <- matrix(dx.ripening[,is.rice], ncol = sum(is.rice))
      max_decrement <- apply(matrix(dx.ripening[-(1:3), ], ncol=sum(is.rice)), 2, which.min)+3 # Maximum decrement will indicate harvest date max_decrement <- max_decrement[is.rice]
      vi.flower <- vi[max.asc]
      eos <- max.asc + max_decrement
      sos <- min.asc
      idx <- mapply(seq, from=as.list(sos), to=as.list(eos), SIMPLIFY = FALSE)
      vi.season <- mapply("[", data.frame(vi), i=idx, SIMPLIFY = FALSE)
      vi.mean <- as.numeric(sapply(vi.season, mean, na.rm=TRUE, USE.NAMES = FALSE))
      vi.sd <- as.numeric(sapply(vi.season, sd, na.rm=TRUE, USE.NAMES = FALSE))
      
      wi.season <- mapply("[", data.frame(wi), i=idx, SIMPLIFY = FALSE)
      wi.mean <- as.numeric(sapply(wi.season, mean, na.rm=TRUE, USE.NAMES = FALSE))
      wi.sd <- as.numeric(sapply(wi.season, sd, na.rm=TRUE, USE.NAMES = FALSE))
      # for(i in 1:length(sos)){
      #   vi.mean <- c(vi.mean, mean(vi[sos[i]:eos[i]]))
      #   vi.sd <- c(vi.sd, sd(vi[sos[i]:eos[i]]))
      # }
      
      # Additional Points
      vi.halfflwr <- vi.flower/2 #vi.sos+(vi.flower-vi.sos)/2
      vi.5 <- matrix(vi[sapply(min.asc, seq, length.out=5)], ncol=5, byrow = TRUE)
      criterion.xiao  <- rowSums((vi.5 > vi.halfflwr))>0
      
      result <- data.frame(sos, eos, flwr=max.asc, vi.sos=vi[sos], vi.eos=vi[eos], vi.flwr=vi.flower, mean.vi=vi.mean, sd.vi=vi.sd, wi.sos=wi[sos], wi.eos=wi[eos], wi.flwr=wi[max.asc], mean.wi=wi.mean, sd.wi=wi.sd, criterion.xiao, met.riceflwrperiod=len.growth>(flowering.period-1))
      #result <- data.frame(sos, eos, flwr=max.asc, vi.sos=vi[sos], vi.eos=vi[eos], vi.flwr=vi.flower, mean.vi=vi.mean, sd.vi=vi.sd, criterion.xiao, met.riceflwrperiod=len.growth>(flowering.period-1))
    } 
  } 
  return(result)
}

rice.itreg <- function(pix.evi, evi.date, dates.covered=NULL, width = 15, increments = 5, eos.evires=0.4, alpha=0.05){
  # Remove NAs
  missing <- which(is.na(pix.evi))
  if(length(missing)>0){
    pix.evi <- pix.evi[-missing]
    evi.date <- evi.date[-missing]
  }
  
  if(is.null(dates.covered)) {
    dates.covered <- evi.date
  } 
  # TODO: dates.covered will be used to look at a specific window. else statement should remove data outside date.covered
  
  idx.starts <- seq(1,length(pix.evi)-width, by=increments)
  
  # Prepare segment data frames for lm
  segments.evi <- as.data.frame(matrix(pix.evi[sapply(idx.starts, FUN=seq, length.out=width)], ncol = length(idx.starts)))
  segments.dates <- as.data.frame(matrix(evi.date[sapply(idx.starts, FUN=seq, length.out=width)], ncol = length(idx.starts)))
  segments.list <- mapply(data.frame, x=segments.dates, y=segments.evi, SIMPLIFY = FALSE)
  
  # Quadratic regression per segment
  inc.lm <- lapply(segments.list, FUN=lm, formula=y~poly(x,2))
  inc.pvals <- sapply(inc.lm, lm.pvalue, lower.tail=FALSE)
  inc.coef1 <- sapply(inc.lm, lm.coeff)
  inc.coef2 <- sapply(inc.lm, lm.coeff, coef.id=2)
  inc.cfpv1 <- sapply(inc.lm, lm.coeff,pvalue=TRUE)
  inc.cfpv2 <- sapply(inc.lm, lm.coeff, coef.id=2, pvalue=TRUE)
  
  # Get segments with significant models, sig coef2 is negative 
  # sig. coef1 is positive proven need not be true  7-Jul-2020, on h25v06 for 2007
  
  potential.rice <- which(inc.pvals<=alpha & inc.coef2<0 & inc.cfpv2<=alpha)
  dat.rice <- vector() 
  
  if(length(potential.rice)>0){
    
    last.rice <- 0 #as.numeric(dates.covered[length(dates.covered)])
    for(i in 1:length(potential.rice)){
      if(last.rice>=dates.covered[idx.starts[potential.rice[i]]]) next
      newx <- data.frame(x=as.numeric(seq(dates.covered[idx.starts[potential.rice[i]]], dates.covered[idx.starts[potential.rice[i]]+width]+7, by="day")))
      new.evi <- predict(inc.lm[[potential.rice[i]]], newdata=newx)
      max.evi <- max(new.evi)
      eos.evi <- max.evi - ((max.evi-new.evi[1])*eos.evires)
      rice.sos <- newx$x[1]
      cnt.rec <- sum(!is.na(segments.evi[,potential.rice[i]]))
      intercept <- lm.coeff(inc.lm[[potential.rice[i]]], coef.id = 0)
      b1 <- lm.coeff(inc.lm[[potential.rice[i]]], coef.id = 1)
      b2 <- lm.coeff(inc.lm[[potential.rice[i]]], coef.id = 2)
      rsq <- lm.rsquare(inc.lm[[potential.rice[i]]])
      pkevi <- which(new.evi==max.evi)    
      # TODO: Verify max.evi as maturity
      post.pk <- which(new.evi[(pkevi+1):length(new.evi)]<=eos.evi)
      if(length(post.pk)>0) {
        post.pk <- min(pkevi+post.pk,length(newx$x))
        rice.eos <- newx$x[post.pk]
        dat.rice <- rbind(dat.rice, c(rice.sos, rice.eos, intercept, b1, b2, rsq))
        last.rice <- rice.eos
      } else next
    }   
  }
  
  return(dat.rice)
}


nonrice.regression <- function(ts.vi, vi.date){
  dat.vi <- data.frame(y=ts.vi, x=vi.date)
  dat.vi <- na.omit(dat.vi)
  # Perform linear regression to immediately check if theres a need for iterative regression
  lm.all <- lm(y~x, data = dat.vi)
  lm.allsum <- summary(lm.all)
  
  dat.rice <- c(lm.allsum$coefficients[1,1],lm.allsum$coefficients[2,1], lm.pvalue(lm.all), lm.allsum$r.squared)
  return(dat.rice)
}

rice.ts <- function(ts.vi, ts.date, ts.flood=NULL, analysis.fun=regression.rice, data.interval=8, crop.duration=160, rnr.only=FALSE, ...){
  # Adjustable Xiao-based rice detection algorithm.
  # Input any vegetation index, any flood vector
  # 
  pts.cd <- ceiling(crop.duration/data.interval) # No. of data points to reach crop duration
  
  # Expected latest possible start of season based on crop duration specified
  flooding.last <- length(ts.vi)-pts.cd
  # Determine rice and non-rice properties, if necessary, on the pixel data
  
  # Find flooding. Normal lowland rice fields are initially flooded
  flooding <- which(ts.flood[1:flooding.last])

  if (length(flooding)>0){
    
    flooding <- sapply(consecutive.groups(flooding), max) # For 
    #evi.flood <- ts.vi[flooding]
    
    # Get evi from start of flooding upto crop duration
    evi.cropdur <- matrix(ts.vi[sapply(flooding, FUN=seq, length.out=pts.cd)],ncol=pts.cd, byrow=TRUE)
    date.cropdur <- matrix(ts.date[sapply(flooding, FUN=seq, length.out=pts.cd)],ncol=pts.cd, byrow=TRUE)
    
    with.snow <- which(rowSums(evi.cropdur < -1, na.rm=TRUE)>0)
    if(length(with.snow)>0){
      evi.cropdur <- evi.cropdur[-with.snow,]
      date.cropdur <- date.cropdur[-with.snow,]
    }
    if(is.null(dim(evi.cropdur))){
      if(length(evi.cropdur)>0)rice.potential <- which(analysis.fun(y=evi.cropdur, x=date.cropdur)) else rice.potential <- vector()
    } else if(nrow(evi.cropdur)>0){
      evi.cropdur <- t(evi.cropdur)
      evi.cropdur <- as.data.frame(evi.cropdur)
      date.cropdur <- t(date.cropdur)
      date.cropdur <- as.data.frame(date.cropdur)
      
      #Compute for evi.rice max if < 0
      rice.potential <- which(unlist(mapply(analysis.fun, y=evi.cropdur, x=date.cropdur)))
    } else rice.potential <- vector()
        
    # If there is rice.potential
    if(length(rice.potential)>0){
      if (rnr.only) {
        rice <- 1
      } else {
        flooding <- flooding[rice.potential]
        if(length(flooding)>1){
          # remove consecutive, use latest flooding
          true.rice <- vector()
          repeat{
            this.flood <- flooding[1]
            if(length(true.rice)==0) true.rice <- this.flood else true.rice <- c(true.rice, this.flood)
            chk.overlap <- flooding-this.flood-15
            flooding <- flooding[chk.overlap>0]
            if(length(flooding)==0) break
          }
        } else true.rice <- flooding
        rice <- sum(10^((1:length(true.rice)-1)*3)*true.rice+(10^(c(2,5,8)[1:length(true.rice)])*(1:length(true.rice))))
      }
    } else {
      rice <- 0
    }
  } else rice <- 0   
  return(rice)
}


# TODO: Quadratic Regression Based Method
optical.rice <- function(data, season_start=1, season_end=1){
# General algorithm for finding rice involves 2 steps.
#   1. Find flooding in the field 
#   2. See if vegetation follows a sudden increase
#
# A more elaborate analysis is to look at the signature from the flooding until 
# the expected end of cropping, characterized by a decrease in vegetation index
    
  	
}



