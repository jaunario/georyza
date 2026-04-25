# Author: jaunario
#
#
###############################################################################

ricestage.period <- function(planting.doy, harvest.doy, year = 1990, stage = rp.whole) {
  if (stage == rp.whole) {
    st <- dateFromDoy(planting.doy, year)
    en <- dateFromDoy(harvest.doy, ifelse(rep(planting.doy > harvest.doy, length(year)), year + 1, year))
  } else if (stage == rp.seed) {
    st <- dateFromDoy(planting.doy, year)
    en <- st + 21
  } else if (stage == rp.pini) {
    harvdate <- dateFromDoy(harvest.doy, ifelse(rep(planting.doy>harvest.doy, length(year)), year + 1, year))
    st <- harvdate - 65 + 14 - 5
    en <- harvdate - 65 + 14 + 5
  } else if(stage == rp.flwr){
    harvdate <- dateFromDoy(harvest.doy, ifelse(rep(planting.doy > harvest.doy, length(year)), year + 1, year))
    st <- harvdate - 30 - 4
    en <- harvdate - 30 + 4
  } else if (stage == rp.mat) {
    harvdate <- dateFromDoy(harvest.doy, ifelse(rep(planting.doy > harvest.doy, length(year)), year + 1, year))
    st <- harvdate - 30 + 4
    en <- harvdate
  } else if(stage == rp.veg) {
    harvdate <- dateFromDoy(harvest.doy, ifelse(rep(planting.doy > harvest.doy, length(year)), year + 1, year))
    st <- dateFromDoy(planting.doy, year) + 21
    en <- dateFromDoy(harvest.doy, ifelse(rep(planting.doy > harvest.doy, length(year)), year + 1, year)) - 65 + 14 - 5
  } else if(stage == rp.repr){
    harvdate <- dateFromDoy(harvest.doy, ifelse(rep(planting.doy > harvest.doy, length(year)), year + 1, year))
    st <- harvdate - 65 + 14 - 5
    en <- harvdate - 30 + 4
  } else {
    stop("Unknown Rice Stage")
  }
  return(list(start = st, end = en))
}
