
slope <- function(y1, y2, x1, x2){
  return((y2-y1)/(x2-x1))
}

curve.Area <- function(x, y){
  h <- diff(x)
  w <- zoo::rollsum(y,k=2)
  return(sum(h * w/2))
}
