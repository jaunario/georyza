#' @title Initialize Python Libraries
#' @description Configures reticulate to use a specific Python version and imports specified libraries into the global environment.
#' @param list.libs Named list of Python libraries to import (default: gdal and osr).
#' @param ... Additional arguments passed to `reticulate::use_python`.
#' @return Logical; status of initialization (currently always returns FALSE).
#' @export
pylibs.start <- function(list.libs = list(pygdal = "gdal", pyosr = "osr"), ...) {
  pystatus <- FALSE
  use_python(...)
  if (length(list.libs) > 0) {
    for (i in seq_along(list.libs)) {
      assign(names(list.libs)[i], import(list.libs[[i]]), envir = .GlobalEnv)
    }
  }
  pystatus
}
