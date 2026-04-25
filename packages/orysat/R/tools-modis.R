# Title: Tools for processing MODIS HDF files
# Description:
# Author: Jorrel Khalil Aunario
# International Rice Research Institute
# Date : 21 Jan 2019
# Version 0.2
# Licence GPL v3
#

PROJ.MODIS <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

# lpdaac.search

#' @title Get MODIS Product Information
#' @description Retrieves information about a MODIS product from a local reference file.
#' @param product Short name of the MODIS product (e.g., "MOD09A1").
#' @param online Logical; if TRUE, attempts to fetch information online (not yet implemented).
#' @return A data frame containing product information.
#' @export
modis.productinfo <- function(product, online = FALSE) {
  # TODO: Create Function that utilizes online LPDAAC Catalog
  if (online) {} else {
    modprods <- read.csv(system.file("satproducts/modis.products.ref.csv", package = "orysat"), stringsAsFactors = FALSE)
  }
  modprods[modprods$ShortName == product, ]
}

#' @title Get MODIS Acquisition DOYs
#' @description Returns the standard acquisition Days of Year (DOYs) for a given MODIS product.
#' @param product Short name of the MODIS product.
#' @return A numeric vector of DOYs.
#' @export
modis.acqdoys <- function(product) {
  prod.info <- modis.productinfo(product)
  acq.by <- switch(prod.info$Temporal.Granularity,
    Daily = 1,
    Yearly = 365,
    "4 day" = 4,
    "8 day" = 8,
    "16 day" = 16,
    NA
  )
  if (is.na(acq.by)) {
    message("Unsupported MODIS Product. Skipping.")
    doy <- vector()
  } else if (acq.by == 16) {
    if (prod.info$Platform == "Terra") {
      doy <- seq(1, 366, by = acq.by)
    } else if (prod.info$Platform == "Aqua") {
      doy <- seq(9, 366, by = acq.by)
    } else {
      message("Unsupported MODIS Product. Skipping.")
      doy <- vector()
    }
  } else {
    doy <- seq(1, 365, by = acq.by)
  }
  doy
}

#' @title Read MODIS HDF File
#' @description Reads a MODIS HDF file and extracts specified layers into a `modis.data` object.
#' @param hdffile Path to the HDF file.
#' @param layer Indices of layers to read (default: all).
#' @param verbose Logical; if TRUE, prints progress messages.
#' @param ... Additional arguments passed to `pylibs.start`.
#' @return A `modis.data` object.
#' @export
modis.readHDF <- function(hdffile, layer = 1, verbose = TRUE, ...) {
  if (!exists("pygdal") || !exists("pyosr")) {
    pylibs.start(...)
  }

  info.hdffile <- unlist(strsplit(basename(hdffile), "\\."))
  m <- new("modis.data")
  m@product <- info.hdffile[1]
  m@acqdate <- info.hdffile[2]
  m@zone <- info.hdffile[3]
  m@version <- info.hdffile[4]
  m@proddate <- info.hdffile[5]
  m@cached <- FALSE
  m@imgvals <- data.frame()
  obj.gdal <- pygdal$Open(hdffile, pygdal$GA_ReadOnly)
  obj.subdatasetlist <- obj.gdal$GetSubDatasets()
  layernames <- vector()
  if (length(layer) < 1 || is.na(layer) || is.null(layer)) layer <- seq_along(obj.subdatasetlist)
  for (i in layer) {
    sdsinfo <- unlist(strsplit(obj.subdatasetlist[[i]][[1]], ":"))
    layernames <- c(layernames, sdsinfo[length(sdsinfo)])
    obj.sds <- pygdal$Open(obj.subdatasetlist[[i]][[1]], pygdal$GA_ReadOnly)
    fillvalue <- as.numeric(obj.sds$GetMetadataItem("_FillValue"))
    scalefactor <- as.numeric(obj.sds$GetMetadataItem("scale_factor"))
    if (length(scalefactor) == 0) scalefactor <- 1
    if (log10(scalefactor) > 0) scalefactor <- 1 / scalefactor

    if (m@projection == "") {
      proj.crs <- pyosr$SpatialReference()
      # proj.crs$ImportFromWkt(obj.sds$GetProjection())

      obj.origin <- obj.sds$GetGeoTransform()
      minX <- obj.origin[[1]]
      resX <- obj.origin[[2]]
      maxX <- minX + resX * obj.sds$RasterXSize

      maxY <- obj.origin[[4]]
      resY <- obj.origin[[6]]
      minY <- maxY + resY * obj.sds$RasterYSize
      m@extent <- extent(minX, maxX, minY, maxY)
      m@projection <- obj.sds$GetProjection()
      m@ncols <- obj.sds$RasterXSize
      m@nrows <- obj.sds$RasterYSize
      rm(proj.crs, minX, maxX, resX, minY, maxY, resY)
    }
    values <- as.numeric(t(obj.sds$ReadAsArray())) * scalefactor
    values[values == fillvalue] <- NA
    if (nrow(m@imgvals) == 0) {
      m@imgvals <- data.frame(values)
    } else {
      m@imgvals[, layernames[length(layernames)]] <- values
    }
    colnames(m@imgvals) <- layernames
    rm(obj.sds)
    gc(reset = TRUE)
  }
  m
}

#' @title Fill with Adjacent Mean
#' @description Fills missing values in a vector with the mean of its neighbors.
#' @param x A numeric vector (typically length 9 for a 3x3 window).
#' @param ... Additional arguments passed to `mean`.
#' @return The filled value (center of the vector if not NA, else the mean).
#' @export
fillwith.adjmean <- function(x, ...) {
  center <- x[5]
  if (!is.na(center)) {
    adjmean <- center
  } else {
    adjmean <- mean(x, ...)
    if (is.nan(adjmean)) adjmean <- NA
  }
  adjmean
}

#' @title Inventory MODIS Files
#' @description Parses MODIS filenames to extract metadata like product, acquisition date, and zone.
#' @param modisfiles Vector of MODIS file paths.
#' @param sep Separator used in filenames (default: ".").
#' @param modisinfo Names of the metadata fields to extract.
#' @param file.ext File extension (default: "tif").
#' @return A data frame with extracted metadata and filenames.
#' @export
inventory.modis <- function(modisfiles, sep = "\\.", modisinfo = c("product", "acqdate", "zone", "version", "band"), file.ext = "tif") {
  # if (format=="GTiff") info <- sub(".hdf","",basename(modisfiles))
  info <- basename(modisfiles)

  # Remove Extension
  info <- gsub(paste(".", file.ext, sep = ""), "", info)

  x <- unlist(strsplit(info, sep))
  m <- data.frame(matrix(x, ncol = length(modisinfo), byrow = TRUE), stringsAsFactors = FALSE)
  if (ncol(m) != length(modisinfo)) {
    message("Non-standard filenames found ")
    return(data.frame())
  }
  colnames(m) <- modisinfo
  m$year <- as.numeric(substr(m$acqdate, 2, 5))
  m$doy <- as.numeric(substr(m$acqdate, 6, 8))
  m$filename <- modisfiles
  m <- m[, c("filename", "year", "doy", modisinfo)]
  if ("band" %in% modisinfo) m$band <- as.character(sub("sur_refl_", "", m$band))
  m
}

#' @title Get Required Bands for Functions
#' @description Extracts the argument names (assumed to be bands) for a list of functions.
#' @param funlist Vector of function names.
#' @param byfun Logical; if TRUE, returns a list of arguments per function.
#' @return A vector or list of required band names.
#' @export
getRequiredBands <- function(funlist, byfun = FALSE) {
  if (byfun) argnames <- list() else argnames <- NULL
  for (i in seq_along(funlist)) {
    if (byfun) argnames[[funlist[i]]] <- names(formals(get(funlist[i]))) else argnames <- c(argnames, names(formals(get(funlist[i]))))
  }
  if (byfun) names(argnames) <- funlist else argnames <- unique(argnames)
  argnames
}

#' @title Format to Extension
#' @description Converts a raster format name to its standard file extension.
#' @param myformat Raster format name (e.g., "raster", "GTiff").
#' @return A string representing the file extension.
#' @export
formatExt <- function(myformat) {
  ext <- rep(NA, length(myformat))
  ext[tolower(myformat) == "raster"] <- "grd"
  ext[tolower(myformat) == "gtiff"] <- "tif"
  ext[tolower(myformat) == "hdf"] <- "hdf"
  ext
}

#' @title Extension to Format
#' @description Converts a file extension to its standard raster format name.
#' @param filext File extension (e.g., "tif", "grd").
#' @return A string representing the raster format name.
#' @export
extFormat <- function(filext) {
  formt <- rep(NA, length(filext))
  formt[tolower(filext) == "tif"] <- "GTiff"
  formt[tolower(filext) == "grd"] <- "raster"
  formt
}

#' @title Get Band Names
#' @description Returns the name of a band given its number and reference set.
#' @param bandnum Band number.
#' @param ref Reference set (default: "ricemap").
#' @return Band name string.
#' @export
bandnames <- function(bandnum, ref = "ricemap") {
  bands <- as.data.frame(cbind(c("red", "nir1", "blue", "green", "nir2", "swir1", "swir2"), c("red", "nir", "blue", "green", NA, "swir1", "swir2")), stringsAsFactors = FALSE)
  colnames(bands) <- c("default", "ricemap")
  bands[bandnum, ref]
}

#' @title Get Band Number
#' @description Returns the number of a band given its name and reference set.
#' @param bandname Band name.
#' @param ref Reference set (default: "ricemap").
#' @param asString Logical; if TRUE, returns band number with "b" prefix.
#' @return Band number (numeric or string).
#' @export
bandnumber <- function(bandname, ref = "ricemap", asString = TRUE) {
  bands <- cbind(c("red", "nir1", "blue", "green", "nir2", "swir1", "swir2"), c("red", "nir", "blue", "green", NA, "swir1", "swir2"), c("red", "nir", "blue", "green", "swir1", "swir2", "swir3"))
  colnames(bands) <- c("default", "ricemap", "web")
  if (asString) result <- paste("b", gsub(" ", 0, format(which(bands[, ref] %in% bandname), width = 2)), sep = "") else result <- which(bands[, ref] %in% bandname)
  result
}

#' @title Compute MODIS Indices
#' @description Computes multiple indices or functions on a `modis.data` object.
#' @param modis A `modis.data` object.
#' @param funlist Vector of function names to apply.
#' @return A data frame containing the results of each function as columns.
#' @export
modis.compute <- function(modis, funlist) {
  result <- layers <- vector()
  for (i in seq_along(funlist)) {
    argnames <- names(formals(get(funlist[i])))
    if (sum(argnames %in% colnames(modis@imgvals)) != length(argnames)) {
      warning("Missing arguments. Did not run ", funlist[i])
    } else {
      layers <- c(layers, toupper(funlist[i]))
      result <- cbind(result, do.call(funlist[i], as.list(modis@imgvals[, argnames])))
    }
  }
  colnames(result) <- layers
  as.data.frame(result)
}

#' @title Convert MODIS Data to Raster Brick
#' @description Converts `modis.data` image values into a RasterBrick or writes them to files.
#' @param mdata A `modis.data` object.
#' @param process Optional label for the process.
#' @param intlayers Indices of layers that should be saved as integers.
#' @param writeto Directory path to write rasters to. If NULL, returns a RasterBrick.
#' @param intNA NA flag for integer layers.
#' @param fltNA NA flag for float layers.
#' @param format Raster format (default: "GTiff").
#' @param skipx Logical; if TRUE, skips existing files.
#' @param ... Additional arguments passed to `writeRaster`.
#' @return A RasterBrick or TRUE if written to files.
#' @export
modis.brick <- function(mdata, process = NULL, intlayers = NULL, writeto = NULL, intNA = -15, fltNA = -9999.0, format = "GTiff", skipx = FALSE, ...) {
  if (!file.exists(writeto)) dir.create(writeto, recursive = TRUE)
  mraster <- raster(modis@extent, ncols = modis@ncols, nrows = modis@nrows, crs = modis@projection)

  if (is.character(writeto)) {
    fname <- gsub("\\.\\.", "\\.", paste(modis@product, modis@acqdate, modis@zone, modis@version, modis@proddate, colnames(modis@imgvals), process, formatExt(format), sep = "."))
    fname <- as.character(paste(writeto, fname, sep = "/"))
  } else {
    mbrick <- brick(mraster)
  }

  for (i in seq_len(ncol(modis@imgvals))) {
    if (skipx && file.exists(fname[i])) next
    mraster <- setValues(mraster, modis@imgvals[, i])

    if (!is.null(intlayers) && is.numeric(intlayers)) {
      dataType(mraster) <- ifelse(i %in% intlayers, "INT1U", "FLT4S")
    } else if (!is.null(intlayers)) {
      message("Ignoring intlayers. Should be seq_len(ncol(modis@imgvals)) instead of", paste(intlayers, collapse = ","))
    }

    if (is.character(writeto)) {
      if (dataType(mraster) == "INT1U") {
        writeRaster(mraster, filename = fname[i], format = format, NAflag = intNA, ...)
      } else {
        writeRaster(mraster, filename = fname[i], format = format, NAflag = fltNA, ...)
      }
    } else {
      mbrick <- addLayer(mbrick, mraster)
    }
  }
  if (is.character(writeto)) mbrick <- TRUE else mbrick@layernames <- colnames(modis@imgvals)
  mbrick
}
