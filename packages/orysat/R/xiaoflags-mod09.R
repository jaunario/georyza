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
  return(!((nir > 0.11) & (ndsi > 0.40)))
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
