# Author: Sonia Asilo, Robert Hijmans, Jorrel Khalil S. Aunario, Andy Nelson, Yann Chemin, Aileen Maunahan
# IRRI
# License GPL3
# Version 2, March 2009

# TODO: verify if this is NDFI
# New NDWI using RED band and SWIR2 in MODIS reflectance product

#' @title Normalized Difference Frequency Index (NDFI)
#' @description Computes the Normalized Difference Frequency Index (or similar) using RED and SWIR2 bands.
#' @param red Vector of RED band reflectance values.
#' @param swir2 Vector of SWIR2 band reflectance values.
#' @return Normalized difference value clipped between -1 and 1.
#' @export
ndfi <- function(red, swir2) {
  result <- normalizedDifference(red, swir2)
  result[result < -1] <- -1
  result[result > 1] <- 1
  result
}

#' @title Difference in EVI and LSWI (DVEL)
#' @description Computes the difference between EVI and LSWI.
#' @param evi Vector of EVI values.
#' @param lswi Vector of LSWI values.
#' @return Difference (evi - lswi).
#' @export
dvel <- function(evi, lswi) {
  evi - lswi
}

#' @title Flood Detection by Threshold
#' @description Detects flooding based on index thresholds.
#' @param index.h2o Vector of water index values.
#' @param low.thres Lower threshold value.
#' @param up.thres Upper threshold value.
#' @param open.circle Logical; if TRUE, uses exclusive range (x > low & x < up), else inclusive.
#' @return Logical vector indicating flooded pixels.
#' @export
flood.bythreshold <- function(index.h2o, low.thres = -inf, up.thres = inf, open.circle = TRUE) {
  if (open.circle) {
    result <- index.h2o > low.thres & index.h2o < up.thres
  } else {
    result <- index.h2o >= low.thres & index.h2o <= up.thres
  }
  result
}

#' @title Xiao Flood Detection Algorithm (flooded1)
#' @description Implements the flood detection algorithm from Xiao et al. (2005).
#' @param lswi Vector of LSWI values.
#' @param ndvi Vector of NDVI values.
#' @param evi Vector of EVI values.
#' @return Logical vector indicating flooded pixels.
#' @references Xiao X., et al. (2005). Mapping paddy rice agriculture in southern China using multi-temporal MODIS images. RSE 95:480-492.
#' @export
flooded1 <- function(lswi, ndvi, evi) {
  res <- (lswi + 0.05 >= evi) | (lswi + 0.05 >= ndvi)
  res
}

#' @title Tidal Flood Detection Algorithm (flooded2)
#' @description Implements the tidal flood detection algorithm from Yan et al. (2010).
#' @param evi Vector of EVI values.
#' @param lswi Vector of LSWI values.
#' @return Numeric vector (1 for flood, 0 for no flood).
#' @references Yan-Er Yan, et al. (2010). Detecting the spatiotemporal changes of tidal flood in the estuarine wetland by using MODIS time series data. Journal of Hydrology 384:156-163.
#' @export
flooded2 <- function(evi, lswi) {
  flood <- rep(NA, length(evi))
  dvl <- dvel(evi, lswi)
  flood[evi > 0.2] <- 0
  flood[!(evi <= 0.2 & dvl <= 0.05)] <- 0
  wrp <- which(evi <= 0.2 & dvl <= 0.05) # water-related points
  flood[wrp[evi[wrp] <= 0.1]] <- 1
  flood[wrp[lswi[wrp] > 0 & lswi[wrp] < 0.2]] <- 0

  flood
}

#' @title Flood Detection Algorithm (flooded3)
#' @description Alternative flood detection algorithm based on EVI and LSWI thresholds.
#' @param evi Vector of EVI values.
#' @param lswi Vector of LSWI values.
#' @return Numeric vector (1 for flood, 0 for no flood).
#' @export
flooded3 <- function(evi, lswi) {
  flood <- rep(NA, length(evi))
  dvl <- dvel(evi, lswi)
  flood[evi > 0.3] <- 0
  flood[!(evi <= 0.3 & dvl <= 0.05)] <- 0
  wrp <- which((evi <= 0.3 & dvl <= 0.05) | (evi <= 0.05 & lswi <= 0)) # water-related points
  flood[wrp[evi[wrp] <= 0.1]] <- 1
  flood
}

#' @title Flood Detection by NDWI7 (flooded4)
#' @description Detects flooding if NDWI (using band 7) is greater than 0.
#' @param ndwi7 Vector of NDWI values using band 7.
#' @return Logical vector indicating flooded pixels.
#' @export
flooded4 <- function(ndwi7) {
  ndwi7 > 0
}

#' @title Persistent Water Detection
#' @description Detects persistent water based on NDVI and LSWI.
#' @param ndvi Vector of NDVI values.
#' @param lswi Vector of LSWI values.
#' @return Logical vector indicating persistent water.
#' @references Xiao X., et al. (2005). Mapping paddy rice agriculture in southern China using multi-temporal MODIS images. RSE 95:480-492.
#' @export
persistentwater <- function(ndvi, lswi) {
  (ndvi < 0.10) & (ndvi < lswi)
}

#' @title Normalized Difference Drought Index (NDDI)
#' @description Computes the NDDI using NDVI and NDWI.
#' @param ndvi Vector of NDVI values.
#' @param ndwi Vector of NDWI values.
#' @return Normalized difference value clipped between 0 and 2.
#' @export
nddi <- function(ndvi, ndwi) {
  result <- normalizedDifference(ndvi, ndwi)
  result[is.infinite(result)] <- NA
  result[result < 0] <- 0

  result[result > 2] <- 2
  result
}

#' @title Drought Detection
#' @description Detects drought based on NDVI and NDWI thresholds.
#' @param ndvi Vector of NDVI values.
#' @param ndwi Vector of NDWI values.
#' @return Logical vector (TRUE for drought).
#' @export
drought <- function(ndvi, ndwi) {
  res <- ((ndvi < 0.5 & ndwi < 0.3) * 2) + ((ndvi > 0.6 & ndwi > 0.4) * 1) - 1
  !res
}

#' @title Normalized Difference Water Index (NDWI)
#' @description Computes the NDWI using GREEN and NIR bands.
#' @param green Vector of GREEN band reflectance values.
#' @param nir Vector of NIR band reflectance values.
#' @return Normalized difference value clipped between -1 and 1.
#' @references McFeeters, S. (1996). The use of Normalized Difference Water Index in the delineation of open water features. IJRS 17(7):1425-1432.
#' @export
ndwi <- function(green, nir) {
  result <- normalizedDifference(green, nir)
  result[is.infinite(result)] <- NA
  result[result < -1] <- -1
  result[result > 1] <- 1
  result
}

#' @title Modified Normalized Difference Water Index (MNDWI)
#' @description Computes the MNDWI using GREEN and SWIR bands.
#' @param green Vector of GREEN band reflectance values.
#' @param swir Vector of SWIR band reflectance values.
#' @return Normalized difference value clipped between -1 and 1.
#' @references Xu, H. (2006). Modification of Normalized Difference Water Index to enhance open water features on remotely sensed imagery. IJRS 27(14):3025-3033.
#' @export
mndwi <- function(green, swir) {
  result <- normalizedDifference(green, swir)
  result[is.infinite(result)] <- NA
  result[result < -1] <- -1
  result[result > 1] <- 1
  result
}

#' @title Land Surface Water Index (LSWI)
#' @description Computes the LSWI using NIR and SWIR bands.
#' @param nir Vector of NIR band reflectance values.
#' @param swir Vector of SWIR band reflectance values.
#' @return Normalized difference value clipped between -1 and 1.
#' @export
lswi <- function(nir, swir) {
  result <- normalizedDifference(nir, swir)
  result[is.infinite(result)] <- NA
  result[result < -1] <- -1
  result[result > 1] <- 1
  result
}


#' @title Generic Water Mapping
#' @description Detects water based on NDVI and Albedo thresholds.
#' @param ndvi Vector of NDVI values.
#' @param albedo Vector of albedo values.
#' @return Logical vector indicating water.
#' @export
water <- function(ndvi, albedo) {
  (ndvi < 0.1) & (albedo < 0.1)
}


#' @title MODIS Water Mapping Tool
#' @description Detects water in MODIS data based on NDVI and Band 7 thresholds.
#' @param ndvi Vector of NDVI values.
#' @param band7 Vector of Band 7 reflectance values.
#' @return Logical vector indicating water.
#' @references Xiao X., et al. (2005). Mapping paddy rice agriculture in southern China using multi-temporal MODIS images. RSE 95:480-492.
#' @references Roy D.P., et al. (2005). Prototyping a global algorithm for systematic fire-affected area mapping using MODIS time series data. RSE 97:137-162.
#' @export
waterModis <- function(ndvi, band7) {
  result <- ndvi * band7
  (ndvi < 0.1) & (band7 < 0.04)
}
