# Author: Sonia Asilo, Robert Hijmans, Jorrel Khalil S. Aunario, Andy Nelson
# IRRI
# License GPL3
# Version 2, March 2009

# TODO: verify if this is NDFI
# New NDWI using RED band and SWIR2 in MODIS reflectance product
ndwi7 <- function(red, swir2){ 
  n7 <- (red - swir2) / (red + swir2)
  return(n7)
}

dvel <- function(evi, lswi){
# DVEL: Difference in the Values of EVI and LSWI
  return(evi-lswi)
}

flood.bythreshold <- function(index.h2o, low.thres= -inf, up.thres=inf, open.circle=TRUE){
  if(open.circle) {
    result <- index.h2o>low.thres & index.h2o < up.thres 
  } else {
    result <- index.h2o >= low.thres & index.h2o <= up.thres 
  }
  return(result)
}

flooded1 <- function(lswi, ndvi, evi) { 
#Xiao X., Boles S., Liu J., Zhuang D., Frokling S., Li C., Salas W., Moore III B. (2005). 
#Mapping paddy rice agriculture in southern China using multi-temporal MODIS images. 
#Remote Sensing of Environment 95:480-492.
	res <- (lswi+0.05 >= evi) | (lswi+0.05 >= ndvi) 
    return(res)
}

flooded2 <- function(evi, lswi){
#Yan-Er Yan, Zu-Tao Ouyangm Hai-Qiang Guo, Shu-Song Jin, Bin Zhao (2010)
#Detecting the spatiotemporal changes of tidal flood in the estuarine wetland by using MODIS time series data
#Journal of Hydrology 384:156-163
  flood <- rep(NA,length(evi))
  dvl <- dvel(evi,lswi)
  flood[evi>0.2] <- 0
  flood[!(evi<=0.2 & dvl<=0.05)] <- 0
  wrp <- which(evi<=0.2 & dvl<=0.05) #water-related points
  flood[wrp[evi[wrp]<=0.1]] <- 1
  flood[wrp[lswi[wrp]>0 & lswi[wrp]<0.2]] <- 0
  
  return(flood)
}

flooded3 <- function(evi,lswi){
  flood <- rep(NA, length(evi))
  dvl <- dvel(evi,lswi)
  flood[evi>0.3] <- 0
  flood[!(evi<=0.3 & dvl<=0.05)] <- 0
  wrp <- which((evi<=0.3 & dvl<=0.05)|(evi<=0.05 & lswi<=0)) #water-related points
  flood[wrp[evi[wrp]<=0.1]] <- 1
  return(flood)     
}

flooded4 <- function(ndwi7){
	return(ndwi7 > 0)	
}

persistentwater <- function(ndvi,lswi){ 
#Xiao X., Boles S., Liu J., Zhuang D., Frokling S., Li C., Salas W., Moore III B. (2005). 
#Mapping paddy rice agriculture in southern China using multi-temporal MODIS images. 
#Remote Sensing of Environment 95:480-492.
	return((ndvi < 0.10) & (ndvi < lswi))
}
