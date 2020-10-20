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

rice.VWIdynamics <- function(vi, wi, flowering.period=7, latest.planting=50, max.evi = 7000, ripening.period=8){
  result <- vector()
  #steps
  # find maxima
  # Should be consistently increasing from x1(DEFAULT=13 to complete 100 days) acqdates before max evi, and consistently decrease x2(DEFAULT=10 maximum) acqdates after max evi)
  # test poly reg to get rsquare
  dx <- diff(vi)
  trend.asc <- consecutive.groups(which(dx>0))
  len.asc <- sapply(trend.asc, length)
  min.asc <- sapply(trend.asc, min)
  enoughtime_toflower <- len.asc>=flowering.period # Period to flowering (max.evi) is constantly increasing trend for at least 56 days (7 ACQDATES, more lenient than 8)
  within_timeframe <- min.asc<=latest.planting # Limits the analysis based on latest planting date
  pot.rice <- which(within_timeframe & enoughtime_toflower)
  
  if(length(pot.rice)>0){
    trend.asc <- trend.asc[pot.rice]    # Limit analysis to segments that passes potentiality
    min.asc <- min.asc[pot.rice]
    max.asc <- sapply(trend.asc, max)+1 # get the point of flowering EVI
    len.asc <- len.asc[pot.rice]
    maturity.period <- sapply(max.asc, seq, length.out=ripening.period)
    
    vi.ripening <- matrix(vi[maturity.period], ncol=ripening.period, byrow = TRUE)
    dx.ripening <- apply(vi.ripening, 1, diff) # Compute Decrements
    #trend.desc <- apply(vi.ripening, 1, diff) # Compute Decrements
    cons.dec5 <- colSums(matrix(dx.ripening[1:5,], nrow=5)<0, na.rm = TRUE) # 
    is.rice <- cons.dec5>=4
    if(sum(is.rice)>0){
      trend.asc <- trend.asc[is.rice]
      min.asc <- min.asc[is.rice]
      max.asc <- max.asc[is.rice]
      
      len.asc <- len.asc[is.rice]
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
      
      result <- data.frame(sos, eos, flwr=max.asc, vi.sos=vi[sos], vi.eos=vi[eos], vi.flwr=vi.flower, mean.vi=vi.mean, sd.vi=vi.sd, wi.sos=wi[sos], wi.eos=wi[eos], wi.flwr=wi[max.asc], mean.wi=wi.mean, sd.wi=wi.sd, criterion.xiao, met.riceflwrperiod=len.asc>(flowering.period-1))
      #result <- data.frame(sos, eos, flwr=max.asc, vi.sos=vi[sos], vi.eos=vi[eos], vi.flwr=vi.flower, mean.vi=vi.mean, sd.vi=vi.sd, criterion.xiao, met.riceflwrperiod=len.asc>(flowering.period-1))
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



