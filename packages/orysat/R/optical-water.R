# Author: Sonia Asilo, Robert Hijmans, Jorrel Khalil S. Aunario, Andy Nelson, Yann Chemin, Aileen Maunahan
# IRRI
# License GPL3
# Version 2, March 2009

# TODO: verify if this is NDFI
# New NDWI using RED band and SWIR2 in MODIS reflectance product

ndfi <- function(red, swir2){ 
  result <- normalizedDifference(red,swir2) 
  result[result < -1] <- -1
  result[result > 1] <- 1
  return(result)
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

nddi <- function(ndvi, ndwi) {
  # NDDI: Normalized Difference Drought Index {
  result<- normalizedDifference(ndvi, ndwi)
  result[is.infinite(result)] <- NA
  result[result < 0] <- 0

  result[result > 2] <- 2
  return(result)
}

drought <- function(ndvi, ndwi) {
  # DROUGHT where drought = 1, no drought=0
  res <- ((ndvi < 0.5 & ndwi < 0.3)*2) + ((ndvi > 0.6 & ndwi > 0.4)*1) - 1
  return(!res)
}

ndwi <- function(green, nir) {
  # NDWI: Normalized Difference Water Index
  # Stuart McFeeters. 1996. The Use of Normalized Difference Water Index in the Delination of
  #	Open Water Features. International Journal of Remote Sensing 27(14):3025-3033
  result<- normalizedDifference(green, nir)
  result[is.infinite(result)] <- NA
  result[result < -1] <- -1
  result[result > 1] <- 1
  return(result)
}

mndwi <- function(green, swir) {
  # MNDWI: Modified Normalized Difference Water Index
  # Hanqui XU. 2006. Modification of Normalized Difference Water Index to Enhance Open Water
  #	Features om Remotely Sensed Imagery. International Journal of Remote Sensing 17(7):1425-1432
  result<- normalizedDifference(green,swir)
  result[is.infinite(result)] <- NA
  result[result < -1] <- -1
  result[result > 1] <- 1
  return(result)
}

lswi<-function(nir, swir) {
  #LSWI: Land Surface Water Index
  result <- normalizedDifference(nir , swir)
  result[is.infinite(result)] <- NA
  result[result < -1] <- -1
  result[result > 1] <- 1
  return(result)
}


water<-function(ndvi, albedo) {
  #water: generic water mapping tool
  return( (ndvi < 0.1) & (albedo < 0.1) )
}


waterModis<-function(ndvi, band7) {
  #water.modis: Terra-MODIS water mapping tool
  #Xiao X., Boles S., Liu J., Zhuang D., Frokling S., Li C., Salas W., Moore III B. (2005). 
  #Mapping paddy rice agriculture in southern China using multi-temporal MODIS images. 
  #Remote Sensing of Environment 95:480-492.
  #
  #Roy D.P., Jin Y., Lewis P.E., Justice C.O. (2005). 
  #Prototyping a global algorithm for systematic fire-affected
  #area mapping using MODIS time series data. 
  #Remote Sensing of Environment 97:137-162.
  result<- ndvi * band7
  return( (ndvi < 0.1) & (band7 < 0.04) )	
  
}

