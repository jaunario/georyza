# Author: Sonia Asilo, Ritsuko Fuchiyama, Robert J. Hijmans, Yann Chemin, Angelo Carlo Pacheco, Jorrel Khalil S. Aunario
# Description: This file contains MOD09A1-specific flags
# International Rice Research Institute
# Date : 21 May 2010
# Version 0,1
# Licence GPL v3

# Cloud Mask

# Internal cloud  algorithm flag
#' @title Internal Cloud Flag
#' @description Checks the internal cloud algorithm flag from MODIS 500m state.
#' @param state_500m MODIS 500m state quality assessment value.
#' @return Logical; TRUE if no cloud is detected.
#' @export
stateflags.internalCloud <- function(state_500m) {
  modis.sqa500f(state_500m) == 1
}

# Cloud bit mask
#' @title Cloud Flag
#' @description Checks the cloud state flag from MODIS 500m state.
#' @param state_500m MODIS 500m state quality assessment value.
#' @return Logical; TRUE if cloudy.
#' @export
stateflags.cloud <- function(state_500m) {
  modis.sqa500a(state_500m) == 1
}

# Cloud shadow bit mask
#' @title Cloud Shadow Flag
#' @description Checks the cloud shadow flag from MODIS 500m state.
#' @param state_500m MODIS 500m state quality assessment value.
#' @return Logical; TRUE if no shadow is detected.
#' @export
stateflags.cloudShadow <- function(state_500m) {
  modis.sqa500b(state_500m) > 0
}

# Water mask
#' @title Water Flag
#' @description Checks the land/water flag from MODIS 500m state.
#' @param state_500m MODIS 500m state quality assessment value.
#' @param exclude Vector of classes to exclude (default: 2 - Ocean coastlines/lake shorelines).
#' @return Logical; TRUE if water (not in exclude list).
#' @export
stateflags.water <- function(state_500m, exclude = 2) {
  exclude <- c(1, exclude)
  !modis.sqa500c(state_500m) %in% exclude
}

# Internal Snow mask
#' @title Snow Flag
#' @description Checks the internal snow mask from MODIS 500m state.
#' @param state_500m MODIS 500m state quality assessment value.
#' @return Logical; TRUE if no snow is detected.
#' @export
stateflags.snow <- function(state_500m) {
  modis.sqa500k(state_500m) == 1
}
