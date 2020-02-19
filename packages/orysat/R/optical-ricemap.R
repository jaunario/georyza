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

# TODO: My interpretation of the Xiao method, with some improvements 
rice.Xiao_v1 <- function(evi, lswi, ndvi=NULL, out.rice=TRUE, flood.fun=flooded1, evi.floodmax=0.25, evi.ricemax=-1, evi.halfricemax=-1, rnr.only=FALSE, data.interval=8, crop.duration=120){
  
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
  flooding <- which(flood.fun(lswi=lswi[1:flooding.last], ndvi=ndvi[1:flooding.last], evi=evi[1:flooding.last]) &
                      evi[1:flooding.last] <= evi.floodmax & evi[1:flooding.last] > -3)
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
      
      # TODO : if too many rice check evi.ricemax should be after 40days
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

# TODO: Quadratic Regression Based Method




optical.rice <- function(vi, h20i, season_start=1, season_end=1){
# General algorithm for finding rice involves 2 steps.
#   1. Find flooding in the field 
#   2. See if vegetation follows a sudden increase
#
# A more elaborate analysis is to look at the signature from the flooding until 
# the expected end of cropping, characterized by a decrease in vegetation index
    
  	
}



