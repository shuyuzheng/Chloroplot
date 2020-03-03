lighten <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col + (255 - col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
darken <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col - col*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

MildColor <- function(color){
  out <- t(col2rgb(color))
  if (sum(out < 127) > 1){
    out <- lighten(color, 0.5)
  } else if(sum(out >= 127) > 1) {
    out <- darken(color, 0.5)
  } else {
    out <- rgb(0, 0, 0)
  }
  return(out)
}
CompColor <- function(color){
  out <- t(col2rgb(color))
  if (all(out < 127)){
    out <- "white"
  } else {
    out <- "black"
  }
  return(out)
}

unTransparent <- function(color_v){
  c <- apply(col2rgb(color_v), 2,
             function(x) rgb(t(x), maxColorValue = 255))
}
