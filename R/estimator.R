#' @title Linear program to estimate spline coefficients
#' @description 
#' Estimate the spline coefficients of the tensor product spline
#' by minimising
#' \deqn{\sum_{i=1}^N|x_i-\phi(u)^TX\psi(y_i)|,}
#' under monotonicity constraint (nonnegative partial derivatives).
#' @param u vector of latent uniform(0,1) variable.
#' @param v uniform(0,1) vector of error.
#' @param w uniform(0,1) vector of observation.
#' @param G if \code{TRUE} estimate function \code{G} (\code{\link{G}}), 
#' otherwise estimate function \code{H} (\code{\link{H}}) (see 'Details').
#' @param type specify spline basis, either \code{"b"} (default), 
#' \code{"c"}, \code{"i"} or \code{"m"};
#' @param splines_control control (see \code{\link{splines.control}}).
#' @return A matrix of spline coefficients
#' @importFrom slam as.simple_triplet_matrix
#' @importFrom Rglpk Rglpk_solve_LP
#' @importFrom stats deriv
lp_fit_Spline <- function(u, v, w, G=TRUE, type="b", splines_control=splines.control()){
  # u - common latent variable from unif(0,1)
  # v - unif(0,1) errors
  # w - observations from unif(0,1)
  # splines_control - list of control for Spline()
  # Note: G(u,v) is suppose to be monotonic increasing in both arguments. 
  # If dim(u) < dim(v), it means that dim(u) = n and dim(v) = n * N.
  
  # verification
  splines_control <- do.call("splines.control",splines_control)
  if(as.integer(splines_control$derivs)!=0) warning("The order of derivative is not 0")
  if(!is.logical(G)) stop("'G' must be logical")
  if(!is.numeric(u)) stop("'u' must be numeric")
  if(!is.numeric(v)) stop("'v' must be numeric")
  if(!is.numeric(w)) stop("'w' must be numeric")
  if(any(u<0 | u>1)) stop("'u' must be between 0 and 1")
  if(any(v<0 | v>1)) stop("'v' must be between 0 and 1")
  if(any(w<0 | w>1)) stop("'w' must be between 0 and 1")
  if(!is.null(dim(u))) stop("'u' should be an object of dimension NULL")
  if(!is.null(dim(v))) stop("'v' should be an object of dimension NULL")
  if(!is.null(dim(w))) stop("'w' should be an object of dimension NULL")
  if(length(w) != length(v)) stop("'w' and 'v' lengths should match")
  if(length(v) %% length(u) != 0) stop("incorrect dimensions for 'u' and 'v'")
  
  # computation
  N <- length(w)
  n <- length(u)
  d <- N / n
  if(d>1) u <- rep(u,d)
  k <- length(splines_control$knots) + 1L
  
  # basis functions (at observations)
  splines_control$x <- u
  phi <- do.call(paste0(type,"Spline"), splines_control)
  if(G) splines_control$x <- w
  if(!G) splines_control$x <- v
  psi <- do.call(paste0(type,"Spline"), splines_control)
  
  # partial derivatives (at the knots)
  splines_control$x <- c(splines_control$knots, 1)
  psi0 <- do.call(paste0(type,"Spline"), splines_control)
  d_psi0 <- deriv(psi0)

  # linear problem is of the form:
  # min C^T\xi under constraint A\xi <= b
  
  # Constraints:
  p <- ncol(phi)
  p2 <- p * p
  C <- c(rep(0,p2),rep(1,N)) 
  pp <- matrix(rep(phi,p), nrow = N) * matrix(rep(t(psi),each=p),nrow = N, byrow=TRUE)
  
  # "regular" constrain
  A1 <- cbind(rbind(pp,-pp)[c(rbind(seq(N),seq(N)+N)),],
              diag(N) %x% rbind(-1,-1))
  if(G) b1 <- cbind(c(rbind(v,-v)))
  if(!G) b1 <- cbind(c(rbind(w,-w)))
  
  # monotonic constrain on partial derivatives
  partial_d_phi <- partial_d_psi <- matrix(nrow = k * k, ncol = p2)
  ind <- 0L
  for(i in seq_len(k)){
    for(j in seq_len(k)){
      ind <- ind + 1L
      partial_d_phi[ind,] <- c(tcrossprod(d_psi0[i,],psi0[j,]))
      partial_d_psi[ind,] <- c(tcrossprod(psi0[i,],d_psi0[j,]))
    }
  }
  A2 <- cbind(
    rbind(-partial_d_phi,-partial_d_psi),
    matrix(0, nrow = 2 * k * k, ncol = N)
  )
  b2 <- cbind(rep(0,2 * k * k))
  
  # equality constraint for I-splines
  # (coefficients sum up to 1)
  if(type=="i"){
    A2 <- rbind(A2, c(rep(1,p2),rep(0,N)))
    b2 <- rbind(b2, 1.0)
  }
  
  # overall constrains
  # every coefficient lies in (0,1) (by default they are assumed to be positive)
  bounds <- list(upper = list(ind=seq_len(p2), val = rep(1.0, p2)))
  A <- rbind(A1,A2)
  b <- rbind(b1,b2)
  f_dir <- rep("<=",nrow(A))
  if(type == "i") f_dir[length(f_dir)] <- "=="
  A <- as.simple_triplet_matrix(A)
  C <- as.simple_triplet_matrix(C)
  fit_lp <- Rglpk_solve_LP(obj = C, mat = A, dir = f_dir, rhs = b, bounds = bounds)
  matrix(fit_lp$solution[seq_len(p2)],nrow=p)
}