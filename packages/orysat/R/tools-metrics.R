slope <- function(y1, y2, x1, x2) {
  (y2 - y1) / (x2 - x1)
}

curve.Area <- function(x, y) {
  h <- diff(x)
  w <- zoo::rollsum(y, k = 2)
  sum(h * w / 2)
}
