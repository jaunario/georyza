# Author: Sonia Asilo, Ritsuko Fuchiyama, Robert J. Hijmans, Yann Chemin, Angelo Carlo Pacheco, Jorrel Khalil S. Aunario 
# Description: This file contains MOD09A1-specific flags
# International Rice Research Institute
# Date : 21 May 2010
# Version 0,1
# Licence GPL v3

# Cloud Mask

# Internal cloud  algorithm flag
stateflags.internalCloud <- function (state_500m) {
  # state_500m <- modis.sqa500f(state_500m)
  # state_500m <- state_500m + 1
  # state_500m[state_500m > 1] <- 0
  return(modis.sqa500f(state_500m)==1)
}

# Cloud bit mask
stateflags.cloud <- function (state_500m) {
  # cloudbit <- modis.sqa500a(state_500m)
  # cloudbit <- state_500m + 1
  # state_500m[state_500m > 1] <- 0
  return(modis.sqa500a(state_500m)==1)
}

# Cloud shadow bit mask
stateflags.cloudShadow <- function (state_500m) {
  # cldshbit <- modis.sqa500b(state_500m)
  # cldshbit <- state_500m + 1
  # cldshbit[state_500m > 1] <- 0
  return(modis.sqa500b(state_500m)>0)
}

# Water mask
stateflags.water <- function(state_500m, exclude = 2){
  exclude <- c(1, exclude)
  return(!modis.sqa500c(state_500m) %in% exclude)
}

#Internal Snow mask
stateflags.snow <- function(state_500m){    
  return(modis.sqa500k(state_500m)==1)
}

