#' @title Calculate Slope
#' @description Computes the slope between two points (x1, y1) and (x2, y2).
#' @param y1 Y-coordinate of the first point.
#' @param y2 Y-coordinate of the second point.
#' @param x1 X-coordinate of the first point.
#' @param x2 X-coordinate of the second point.
#' @return The calculated slope.
#' @export
slope <- function(y1, y2, x1, x2) {
  (y2 - y1) / (x2 - x1)
}

#' @title Calculate Area Under Curve (AUC)
#' @description Computes the area under a curve using the trapezoidal rule.
#' @param x Vector of X-coordinates.
#' @param y Vector of Y-coordinates.
#' @return The calculated area.
#' @export
curve.Area <- function(x, y) {
  h <- diff(x)
  w <- zoo::rollsum(y, k = 2)
  sum(h * w / 2)
}
