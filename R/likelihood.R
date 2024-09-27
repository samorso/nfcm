#' @title Control for spline
#' @description Control for the spline basis, currently from
#' the \code{\link{splines2}} package.
#' @param df degree of freedom (can be \code{NULL});
#' @param knots intenal breakpoints that define the splines;
#' @param degree the degree of the piecewise polynomial;
#' @param intercept boolean;
#' @param periodic boolean;
#' @param Boundary.knots boundary points at which to anchor the splines;
#' @param derivs the order of the derivative the spline;
#' @param integral boolean.
#' @seealso \code{\link[splines2]{bSpline}}, \code{\link[splines2]{mSpline}},
#' \code{\link[splines2]{cSpline}}, \code{\link[splines2]{iSpline}}
#' @export
splines.control <- function(df = NULL, knots = seq(.2,.8,by=.2), degree = 3L, 
                            intercept = FALSE, Boundary.knots = c(0,1),
                            periodic = FALSE, derivs = 0L, integral = FALSE){
  if(!is.null(df)) if(as.integer(df)<0) stop("'df' must be positive")
  if(min(knots)<0 || max(knots)>1) warning("'knots' defined outside (0,1)")
  if(as.integer(degree)<0) stop("'degree' must be positive")
  if(!is.logical(intercept)) stop("'intercept' must be a boolean")
  if(any(Boundary.knots>1) || any(Boundary.knots<0)) stop("'Boundary.knots' must be in (0,1)")
  if(!is.logical(periodic)) stop("'periodic' must be a boolean")
  if(as.integer(derivs)<0) stop("'derivs' must be positive")
  if(!is.logical(integral)) stop("'integral' must be a boolean")
  
  list(df=df, knots=knots, degree=degree, intercept=intercept, Boundary.knots=Boundary.knots, periodic=periodic, derivs=derivs, integral=integral)
}

# --------------
# negative log-likelihood (bivariate case)
# --------------
#' @title Negative approximate log-likelihood of Spline Based Nonlinear Factor Copula Model
#' @description 
#' Compute an upper bound of the negative log-likelihood for a bivariate model. 
#' It is defined 
#' as \deqn{\frac{1}{n}\frac{1}{N}\sum_{i=1}^n\sum_{j=1}^N\Phi^\top\Lambda \psi'(w_{ij}),}
#' where \eqn{\psi'} is the first derivative of a spline basis \eqn{\psi},
#' \eqn{\Lambda} is a matrix of spline coefficients and \eqn{\Phi} is 
#' defined as \deqn{\int_0^1\psi(u)du.} It is assumed that the variables have a
#' positive quadrant dependence. 
#' @param lambda vector of spline coefficients (\eqn{\Lambda} vectorized).
#' @param w uniform(0,1) \eqn{n\times N} \code{matrix} of observations.
#' @param type specify spline basis, either \code{"b"} (default), 
#' \code{"c"}, \code{"i"} or \code{"m"};
#' @param splines_control control (see \code{\link{splines.control}}).
#' @return The approximate negative log-likelihood divided by \eqn{n} and \eqn{N}.
#' @importFrom splines2 bSpline cSpline iSpline mSpline
#' @export
nfcm_nll <- function(lambda, w, type="b", splines_control=splines.control()){
  # lambda            - a square matrix (not symmetric) vectorized of spline coefficients
  # w                 - a matrix of samples
  # type              - spline basis
  # splines_control   - splines control for spline
  splines_control <- do.call("splines.control",splines_control)
  
  # verification
  if(!type %in% c("b","c","i","m")) stop("Spline basis not supported")
  k <- splines_control$df
  if(is.null(k)) k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
  if(!is.matrix(w)) stop("'w' is not a matrix")
  if(sqrt(length(lambda)) != k) stop("'lambda' has not the correct dimension")
  
  # compute C := \Lambda^T \Phi^T
  Lambda <- matrix(lambda, nrow = k, ncol = k)
  if(as.integer(splines_control$derivs)!=0L) splines_control$derivs <- 0L
  if(!isTRUE(splines_control$integral)) splines_control$integral <- TRUE
  splines_control$x <- 1
  Phi <- do.call(paste0(type,"Spline"), splines_control)
  splines_control$integral <- FALSE
  C <- Phi %*% Lambda
  
  # Compute the tensor products
  splines_control$derivs <- 1L
  n <- nrow(w)
  N <- ncol(w)
  tp_spline <- matrix(nrow = n, ncol = N)
  for(j in 1:N){
    splines_control$x <- w[,j]
    psip <- do.call(paste0(type,"Spline"), splines_control)
    tp_spline[,j] <- tcrossprod(psip, C)
  }
  
  # compute the negative log-likelihood
  -mean(log(tp_spline))
}

#' @rdname nfcm_nll
#' @return The gradient of the approximate negative log-likelihood divided by \eqn{n} and \eqn{N}.
#' @export
nfcm_grad_nll <- function(lambda, w, type="b",splines_control=splines.control()){
  # lambda            - a square matrix (not symmetric) vectorized of spline coefficients
  # w                 - a matrix of samples
  # type              - spline basis
  # splines_control   - splines control for spline
  splines_control <- do.call("splines.control",splines_control)
  
  # verification
  if(!type %in% c("b","c","i","m")) stop("Spline basis not supported")
  k <- splines_control$df
  if(is.null(k)) k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
  if(!is.matrix(w)) stop("'w' is not a matrix")
  if(sqrt(length(lambda)) != k) stop("'lambda' has not the correct dimension")
  
  # compute C := \Lambda^T \Phi^T
  Lambda <- matrix(lambda, nrow = k, ncol = k)
  if(as.integer(splines_control$derivs)!=0L) splines_control$derivs <- 0L
  if(!isTRUE(splines_control$integral)) splines_control$integral <- TRUE
  splines_control$x <- 1
  Phi <- do.call(paste0(type,"Spline"), splines_control)
  splines_control$integral <- FALSE
  C <- Phi %*% Lambda
  
  # Compute the tensor products
  splines_control$derivs <- 1L
  n <- nrow(w)
  N <- ncol(w)
  A <- matrix(0, nrow = k, ncol = k)
  for(j in 1:N){
    splines_control$x <- w[,j]
    psip <- do.call(paste0(type,"Spline"), splines_control)
    tp_spline <- tcrossprod(psip, C)
    for(i in 1:n){
      mat <- crossprod(Phi,psip[i,]) / tp_spline[i]
      A <- A + mat
    }
  }
  -c(A) / n / N
}

#' @rdname nfcm_nll
#' @return Akaike information criterion
#' @export
nfcm_aic <- function(lambda, w, type="b", splines_control=splines.control()){
  k <- length(lambda)
  n <- nrow(w)
  nfcm_nll(lambda=lambda,w=w,type=type,splines_control=splines_control) + 2.0 * k / n
}

#' @rdname nfcm_nll
#' @return Bayesian information criterion
#' @export
nfcm_bic <- function(lambda, w, type="b", splines_control=splines.control()){
  k <- length(lambda)
  n <- nrow(w)
  nfcm_nll(lambda=lambda,w=w,type=type,splines_control=splines_control) + log(n) * k / n
}