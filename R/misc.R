# --------------
# Miscellanous
# --------------
# TODO: create function for plotting H and G
# #' @title Plotting the function H or G
# #' @export
# plot.H <- function(x, n=100, ...) {
#   # x is an object of class 'H' (see 'nfcm_mle')
#   xx <- seq(1e-3, 1-1e-3, length.out = n)
#   z <- matrix()
#   
#   persp(z, ...)
# }

#' @title Plotting the function H or G approximated by splines
#' @export
plot.H.spline <- function(x, n=100, ...) {
  # x is an object of class 'H' (see 'nfcm_mle')
  xx <- seq(1e-3, 1-1e-3, length.out = n)
  control <- x$control
  
  # basis functions (at observations)
  control$x <- xx
  phi <- do.call(paste0(x$type,"Spline"), control)
  phi_par <- phi %*% x$par
  z <- matrix(ncol = n, nrow = n) 
  for(i in 1:n) {
    for(j in 1:n) {
      z[i,j] <- phi_par[i,] %*% phi[j,]
    }
  }
  # z <- matrix(rowSums(matrix(rep(phi, nrow(phi)),ncol=ncol(phi),byrow=TRUE) %*% x$par * matrix(rep(phi, each=nrow(phi)), ncol=ncol(phi), byrow=TRUE)), ncol = nrow(phi))
  persp(z, ...)
}

#' @rdname plot.H.spline
#' @export
plot.G.spline <- function(x, n=100, ...) {
  # x is an object of class 'G' (see 'nfcm_mle')
  xx <- seq(1e-3, 1-1e-3, length.out = n)
  control <- x$control
  
  # basis functions (at observations)
  control$x <- xx
  phi <- do.call(paste0(x$type,"Spline"), control)
  phi_par <- phi %*% x$par
  z <- matrix(ncol = n, nrow = n) 
  for(i in 1:n) {
    for(j in 1:n) {
      z[i,j] <- phi_par[i,] %*% phi[j,]
    }
  }
  # z <- matrix(rowSums(matrix(rep(phi, nrow(phi)),ncol=ncol(phi),byrow=TRUE) %*% x$par * matrix(rep(phi, each=nrow(phi)), ncol=ncol(phi), byrow=TRUE)), ncol = nrow(phi))
  persp(z, ...)
}

