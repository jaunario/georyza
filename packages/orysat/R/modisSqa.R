# Author: Yann Chemin
# IRRI
# License GPL3
# Veon 1, October 2008

# MODIS State Quality Assessment Extractor
# Makes Human-readable images of State Quality Assessment binary bits from MOD09A 500m products (Vermote et al., 2008), this is the proper place to find cloud related information.
# Vermote E.F., Kotchenova S.Y., Ray J.P. MODIS Surface Reflectance User's Guide. Veon 1.2. June 2008. MODIS Land Surface Reflectance Science Computing Facility. Homepage

#' @title Extract Cloud State
#' @description Extracts Cloud State unsigned int bits[0-1] from MODIS 500m State Quality Assessment.
#' @param pixel Unsigned integer representing the SQA pixel.
#' @return Integer (0-3): 0 (clear), 1 (cloudy), 2 (mixed), 3 (not set, assumed clear).
#' @export
modis.sqa500a <- function(pixel) {
  bitAnd(pixel, 3)
}

#' @title Extract Cloud Shadow
#' @description Extracts Cloud Shadow unsigned int bit[2] from MODIS 500m State Quality Assessment.
#' @param pixel Unsigned integer representing the SQA pixel.
#' @return Integer (0-1): 0 (yes), 1 (no).
#' @export
modis.sqa500b <- function(pixel) {
  pixel <- bitShiftR(pixel, 2)
  bitAnd(pixel, 1)
}

#' @title Extract Land/Water Flag
#' @description Extracts Land/Water Flag unsigned int bits[3-5] from MODIS 500m State Quality Assessment.
#' @param pixel Unsigned integer representing the SQA pixel.
#' @return Integer (0-7): 0 (Shallow ocean), 1 (Land), 2 (Ocean coastlines/lake shorelines), 3 (Shallow inland water), 4 (Ephemeral water), 5 (Deep inland water), 6 (Continental/moderate ocean), 7 (Deep ocean).
#' @export
modis.sqa500c <- function(pixel) {
  pixel <- bitShiftR(pixel, 3)
  bitAnd(pixel, 7)
}

#' @title Extract Aerosol Quantity
#' @description Extracts Aerosol Quantity unsigned int bits[6-7] from MODIS 500m State Quality Assessment.
#' @param pixel Unsigned integer representing the SQA pixel.
#' @return Integer (0-3): 0 (Climatology), 1 (Low), 2 (Average), 3 (High).
#' @export
modis.sqa500d <- function(pixel) {
  pixel <- bitShiftR(pixel, 6)
  bitAnd(pixel, 3)
}

#' @title Extract Cirrus Detected
#' @description Extracts Cirrus Detected unsigned int bits[8-9] from MODIS 500m State Quality Assessment.
#' @param pixel Unsigned integer representing the SQA pixel.
#' @return Integer (0-3): 0 (None), 1 (Small), 2 (Average), 3 (High).
#' @export
modis.sqa500e <- function(pixel) {
  pixel <- bitShiftR(pixel, 8)
  bitAnd(pixel, 3)
}

#' @title Extract Internal Cloud Algorithm Flag
#' @description Extracts Internal Cloud Algorithm Flag unsigned int bit[10] from MODIS 500m State Quality Assessment.
#' @param pixel Unsigned integer representing the SQA pixel.
#' @return Integer (0-1): 0 (Cloud), 1 (No cloud).
#' @export
modis.sqa500f <- function(pixel) {
  pixel <- bitShiftR(pixel, 10)
  bitAnd(pixel, 1)
}

#' @title Extract Internal Fire Algorithm Flag
#' @description Extracts Internal Fire Algorithm Flag unsigned int bit[11] from MODIS 500m State Quality Assessment.
#' @param pixel Unsigned integer representing the SQA pixel.
#' @return Integer (0-1): 0 (Fire), 1 (No fire).
#' @export
modis.sqa500g <- function(pixel) {
  pixel <- bitShiftR(pixel, 11)
  bitAnd(pixel, 1)
}

#' @title Extract MOD35 Snow/Ice Flag
#' @description Extracts MOD35 Snow/Ice Flag unsigned int bit[12] from MODIS 500m State Quality Assessment.
#' @param pixel Unsigned integer representing the SQA pixel.
#' @param bandno Optional band number (not used in current implementation).
#' @return Integer (0-1): 0 (Yes, snow), 1 (No).
#' @export
modis.sqa500h <- function(pixel, bandno) {
  pixel <- bitShiftR(pixel, 12)
  bitAnd(pixel, 1)
}

#' @title Extract Pixel Adjacent to Cloud Flag
#' @description Extracts Pixel Adjacent to Cloud flag unsigned int bit[13] from MODIS 500m State Quality Assessment.
#' @param pixel Unsigned integer representing the SQA pixel.
#' @return Integer (0-1): 0 (Yes), 1 (No).
#' @export
modis.sqa500i <- function(pixel) {
  pixel <- bitShiftR(pixel, 13)
  bitAnd(pixel, 1)
}

#' @title Extract BRDF Correction Flag
#' @description Extracts BRDF correction performed flag unsigned int bit[14] from MODIS 500m State Quality Assessment.
#' @param pixel Unsigned integer representing the SQA pixel.
#' @return Integer (0-1): 0 (Yes), 1 (No).
#' @export
modis.sqa500j <- function(pixel) {
  pixel <- bitShiftR(pixel, 14)
  bitAnd(pixel, 1)
}

#' @title Extract Internal Snow Mask
#' @description Extracts Internal Snow Mask unsigned int bit[15] from MODIS 500m State Quality Assessment.
#' @param pixel Unsigned integer representing the SQA pixel.
#' @return Integer (0-1): 0 (Snow), 1 (No snow).
#' @export
modis.sqa500k <- function(pixel) {
  pixel <- bitShiftR(pixel, 15)
  bitAnd(pixel, 1)
}
