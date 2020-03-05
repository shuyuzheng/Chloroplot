lighten <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- grDevices::col2rgb(color)
  col <- col + (255 - col)*factor
  col <- grDevices::rgb(t(col), maxColorValue=255)
  col
}
darken <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- grDevices::col2rgb(color)
  col <- col - col*factor
  col <- grDevices::rgb(t(col), maxColorValue=255)
  col
}

MildColor <- function(color){
  out <- t(grDevices::col2rgb(color))
  if (sum(out < 127) > 1){
    out <- lighten(color, 0.5)
  } else if(sum(out >= 127) > 1) {
    out <- darken(color, 0.5)
  } else {
    out <- grDevices::rgb(0, 0, 0)
  }
  return(out)
}
CompColor <- function(color){
  out <- t(grDevices::col2rgb(color))
  if (all(out < 127)){
    out <- "white"
  } else {
    out <- "black"
  }
  return(out)
}

Transparent <- function(color_v, degree = 0.2){
  c <- apply(grDevices::col2rgb(color_v), 2,
             function(x) grDevices::rgb(t(x), alpha = round(degree * 255),
                             maxColorValue = 255))
  return(c)
}
