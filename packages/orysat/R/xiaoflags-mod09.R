# Author: Sonia Asilo, Ritsuko Fuchiyama, Robert J. Hijmans, Yann Chemin, Angelo Carlo Pacheco, Jorrel Khalil S. Aunario 
# Description: This file contains MOD09A1-specific quality flags
# International Rice Research Institute
# Date : 21 May 2010
# Version 0,1
# Licence GPL v3


# Blue band-based cloud mask
xiaoflags.cloud <- function(b03, scale=0.0001){    
  return((b03*scale) >= 0.2)
}

#second snow mask
xiaoflags.snow <- function(ndsi, nir) {    
  #	res <- 
  #	res[(nir > 0.11) & (ndsi > 0.40)] <- 0
  #	res[is.na(res)] <- -15
  return((nir > 0.11) & (ndsi > 0.40))
}

xiaoflags.snowadj <- function(ndsi, green, nir) {
  # TODO: Reference?
  #	res <- 
  #	res[(nir > 0.10) & (green > 0.10) & (ndsi >= 0.40)] <- 0
  #	res[is.na(res)] <- -15
  return((nir > 0.10) & (green > 0.10) & (ndsi >= 0.40))
}


xiaoflags.forest <- function(ndvi){
  return(sum(ndvi>=0.7, na.rm=TRUE)>=20)
}


xiaoflags.shrub <- function(lswi){
  return(sum(lswi<0.1, na.rm=TRUE)==0)
}


xiaoflags.persistentwater <- function(ndvi, lswi){
  return(sum((ndvi < 0.10) & (ndvi < lswi), na.rm=TRUE) >= 10)
}


xiaoflags.flooded <- function(evi, ndvi, lswi) {    
  #	res <- 
  #	res[(nir > 0.11) & (ndsi > 0.40)] <- 0
  #	res[is.na(res)] <- -15
  return(!((nir > 0.11) & (ndsi > 0.40)))
}

xiaoflags.rice <- function(evi, evi.ricemax=NULL, evi.halfricemax=NULL, data.interval=8, ricemax.post40=TRUE){
  pts.40d <- ceiling(40/data.interval) # No of data points to reach 40 days
  
  if(is.null(evi.ricemax)){
    if(ricemax.post40)  evi.ricemax <- max(evi[(pts.40d+1):length(evi)],na.rm = TRUE) else evi.ricemax <- max(evi,na.rm = TRUE)
  }

  if(is.null(evi.halfricemax)){
    evi.halfricemax <- evi.ricemax/2 # evi.ricemax-((evi.ricemax-evi[1])/2)
  }
  
  if(ricemax.post40){
    result <- sum(evi[1:pts.40d]>=evi.halfricemax)>0 & sum(evi[1:pts.40d]>=evi.ricemax)==0
  } else {
    result <- sum(evi[1:pts.40d]>=evi.halfricemax)>0
  } 
  return(result)
}

modis.mask <- function(modvals, masks){
	#DEPRECATE
    masks <- as.matrix(masks)
    m <- rowSums(masks, na.rm=TRUE)
    mm <- which(m<ncol(masks)) 
    if (is.null(ncol(modvals))) modvals[mm] <- NA else {
        for (i in 1:ncol(modvals)) modvals[mm,i] <- NA
    }
    return(modvals)
}
