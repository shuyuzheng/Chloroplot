CompColor <- function(color){
  out <- t(col2rgb(color))
  if (all(out > 127)){
    out <- rgb(0, 0, 0, maxColorValue = 255)
  } else {
    out <- rgb(255, 255, 255, maxColorValue = 255)
  }
  return(out)
}

