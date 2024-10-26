
########### Class ##############

density_object <- function(x, h = NULL) {
  structure(
    list(
      x = x,
      h = h,
      cv_bw = cv_bw_l_fast(x),
      amise_bw = AMISE_bw(x)
    ),
    class = "density_object"
  )
}

kern_dens <- function(x) {
  UseMethod("kern_dens")
}

kern_dens <- function(object, h = "CV", binned = FALSE) {
  if(is.null(h)){
    h <- object$cv_bw
  }
  else if(h == "CV"){
    h <- object$cv_bw
  }
  else if(h == "AMISE"){
    h <- object$amise_bw
  }
  if(binned){
      kern_dens_bin(object$x, h)
    }
  kern_dens1(object$x, h)
}
