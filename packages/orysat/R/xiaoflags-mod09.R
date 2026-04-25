# Author: Sonia Asilo, Ritsuko Fuchiyama, Robert J. Hijmans, Yann Chemin, Angelo Carlo Pacheco, Jorrel Khalil S. Aunario
# Description: This file contains MOD09A1-specific quality flags
# International Rice Research Institute
# Date : 21 May 2010
# Version 0,1
# Licence GPL v3


# Blue band-based cloud mask
#' @title Xiao Cloud Mask (Blue Band)
#' @description Detects clouds based on the blue band (B03) reflectance threshold.
#' @param b03 Vector of blue band reflectance values.
#' @param scale Scale factor for the reflectance values (default: 0.0001).
#' @return Logical vector indicating clouds.
#' @export
xiaoflags.cloud <- function(b03, scale = 0.0001) {
  (b03 * scale) >= 0.2
}

# second snow mask
#' @title Xiao Snow Mask
#' @description Detects snow based on NDSI and NIR reflectance thresholds.
#' @param ndsi Vector of NDSI values.
#' @param nir Vector of NIR reflectance values.
#' @return Logical vector indicating snow.
#' @export
xiaoflags.snow <- function(ndsi, nir) {
  (nir > 0.11) & (ndsi > 0.40)
}

#' @title Xiao Adjusted Snow Mask
#' @description Detects snow based on NDSI, GREEN, and NIR reflectance thresholds.
#' @param ndsi Vector of NDSI values.
#' @param green Vector of GREEN reflectance values.
#' @param nir Vector of NIR reflectance values.
#' @return Logical vector indicating snow.
#' @export
xiaoflags.snowadj <- function(ndsi, green, nir) {
  (nir > 0.10) & (green > 0.10) & (ndsi >= 0.40)
}


#' @title Xiao Forest Check
#' @description Checks if a pixel belongs to a forest based on NDVI history.
#' @param ndvi Vector of NDVI values over time.
#' @return Logical; TRUE if the pixel is considered forest (at least 20 acquisitions with NDVI >= 0.7).
#' @export
xiaoflags.forest <- function(ndvi) {
  sum(ndvi >= 0.7, na.rm = TRUE) >= 20
}


#' @title Xiao Shrub Check
#' @description Checks if a pixel belongs to a shrub based on LSWI history.
#' @param lswi Vector of LSWI values over time.
#' @return Logical; TRUE if the pixel is considered shrub (no LSWI < 0.1 detected).
#' @export
xiaoflags.shrub <- function(lswi) {
  sum(lswi < 0.1, na.rm = TRUE) == 0
}


#' @title Xiao Persistent Water Check
#' @description Checks if a pixel belongs to persistent water based on NDVI and LSWI history.
#' @param ndvi Vector of NDVI values over time.
#' @param lswi Vector of LSWI values over time.
#' @return Logical; TRUE if the pixel is considered persistent water (at least 10 acquisitions meeting the water criterion).
#' @export
xiaoflags.persistentwater <- function(ndvi, lswi) {
  sum((ndvi < 0.10) & (ndvi < lswi), na.rm = TRUE) >= 10
}


#' @title Xiao Flooded Flag (Place-holder)
#' @description Place-holder for Xiao-based flood detection.
#' @param evi Vector of EVI values.
#' @param ndvi Vector of NDVI values.
#' @param lswi Vector of LSWI values.
#' @return Logical vector (inverse of snow detection logic in current implementation).
#' @export
xiaoflags.flooded <- function(evi, ndvi, lswi) {
  !((nir > 0.11) & (ndsi > 0.40))
}

#' @title Xiao Rice Flag
#' @description Checks if a pixel matches rice growth criteria based on EVI dynamics after flooding.
#' @param evi Vector of EVI values for a cropping season.
#' @param evi.ricemax Optional maximum EVI for rice.
#' @param evi.halfricemax Optional half-maximum EVI for rice.
#' @param data.interval Data interval in days.
#' @param ricemax.post40 Logical; if TRUE, uses max EVI after 40 days of flooding.
#' @return Logical; TRUE if criteria are met.
#' @export
xiaoflags.rice <- function(evi, evi.ricemax = NULL, evi.halfricemax = NULL, data.interval = 8, ricemax.post40 = TRUE) {
  pts.40d <- ceiling(40 / data.interval) # No of data points to reach 40 days

  if (is.null(evi.ricemax)) {
    if (ricemax.post40) evi.ricemax <- max(evi[(pts.40d + 1):length(evi)], na.rm = TRUE) else evi.ricemax <- max(evi, na.rm = TRUE)
  }

  if (is.null(evi.halfricemax)) {
    evi.halfricemax <- evi.ricemax / 2 # evi.ricemax-((evi.ricemax-evi[1])/2)
  }

  if (ricemax.post40) {
    result <- sum(evi[1:pts.40d] >= evi.halfricemax) > 0 & sum(evi[1:pts.40d] >= evi.ricemax) == 0
  } else {
    result <- sum(evi[1:pts.40d] >= evi.halfricemax) > 0
  }
  result
}

#' @title Mask MODIS Values
#' @description Applies masks to MODIS data values.
#' @param modvals Data frame or matrix of MODIS values.
#' @param masks Logical matrix of masks.
#' @return Masked MODIS values (with NA where mask is FALSE).
#' @note This function is potentially slated for deprecation.
#' @export
modis.mask <- function(modvals, masks) {
  # DEPRECATE
  masks <- as.matrix(masks)
  m <- rowSums(masks, na.rm = TRUE)
  mm <- which(m < ncol(masks))
  if (is.null(ncol(modvals))) {
    modvals[mm] <- NA
  } else {
    for (i in seq_len(ncol(modvals))) modvals[mm, i] <- NA
  }
  modvals
}
