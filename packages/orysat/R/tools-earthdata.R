# Title: Tools for using NASA Earthdata
# Author: Jorrel Khalil S. Aunario
# IRRI
# License GPL3
# Version 1, April 2025

APPEARS_URL <- "https://appeears.earthdatacloud.nasa.gov/api/"

DATAPOOL_PLATFORMS <- data.frame(
  platform = c("ASTER GDEM", "Aqua MODIS", "Combined MODIS", "Terra MODIS", "Suomi NPP VIIRS", "ECOSTRESS", "MEaSUREs"),
  subdirectory = c("ASTT", "MOLA", "MOTA", "MOLT", "VIIRS", "ECOSTRESS", "MEASURES"),
  stringsAsFactors = FALSE
)

#' @title Get Earthdata Product Information
#' @description Retrieves product information from the NASA AppEEARS API.
#' @param product Optional product name filter (ProductAndVersion).
#' @param update Logical; if TRUE, updates the local JSON reference file from the API.
#' @return A data frame or list containing product information.
#' @export
earthdata.product <- function(product = NULL, update = TRUE) {
  if (update || !file.exists(system.file("satproducts/earthdata.products.ref.json", package = "orysat"))) {
    prod.info <- httr::GET(paste0(APPEARS_URL, "product"))
    prod.info <- jsonlite::toJSON(httr::content(prod.info), auto_unbox = TRUE)
    writeLines(prod.info, con = file.path(system.file("satproducts", package = "orysat"), "earthdata.products.ref.json"))
  } else {
    prod.info <- readLines(system.file("satproducts/earthdata.products.ref.json", package = "orysat"))
  }
  prod.info <- jsonlite::fromJSON(prod.info)
  # prod.info <- setNames(split(prod.info, seq_len(nrow(prod.info))), prod.info$ProductAndVersion)
  if (!is.null(product)) prod.info <- prod.info[grep(product, prod.info$ProductAndVersion), ]
  prod.info
}

#' @title Configure Earthdata Credentials
#' @description Creates a .netrc file with NASA Earthdata credentials.
#' @param netrc.file Path to the .netrc file (default: ".netrc").
#' @param user Earthdata username.
#' @param password Earthdata password.
#' @return Logical; TRUE if successful.
#' @note On Windows, the .netrc file usually works best in the working directory. On Linux, it should be in the home directory.
#' @export
earthdata.config <- function(netrc.file = ".netrc", user, password) {
  if (!file.exists(netrc.file)) {
    success <- try(writeLines(c(
      "machine urs.earthdata.nasa.gov",
      paste("login", user),
      paste("password", password)
    ),
    con = netrc.file
    ), silent = TRUE)
  } else {
    success <- TRUE
  }
  class(success) != "try-error"
}
