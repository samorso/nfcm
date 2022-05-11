
#' @title Control for spline
#' @description Control for the spline basis, currently from
#' the \code{\link{splines2}} package.
#' @param df degree of freedom (can be \code{NULL});
#' @param knots intenal breakpoints that define the splines;
#' @param degree the degree of the piecewise polynomial;
#' @param intercept boolean;
#' @param Boundary.knots boundary points at which to anchor the splines;
#' @param derivs the order of the derivative the spline.
#' @seealso \code{\link[splines2]{bSpline}}, \code{\link[splines2]{mSpline}},
#' \code{\link[splines2]{cSpline}}, \code{\link[splines2]{iSpline}}
#' @export
splines.control <- function(df = NULL, knots = seq(.2,.8,by=.2), degree = 3L, 
                            intercept = FALSE, Boundary.knots = c(0,1),
                            derivs = 0L){
  if(!is.null(df)) if(as.integer(df)<0) stop("'df' must be positive")
  if(min(knots)<0 || max(knots)>1) warning("'knots' defined outside (0,1)")
  if(as.integer(degree)<0) stop("'degree' must be positive")
  if(!is.logical(intercept)) stop("'intercept' must be a boolean")
  if(as.integer(derivs)<0) stop("'derivs' must be positive")
  
  list(df=df, knots=knots, degree=degree, intercept=intercept, Boundary.knots=Boundary.knots, derivs=derivs)
}

#' @title "P" matrix
#' @description Compute the "P" matrix, which is defined as 
#' \deqn{\int_0^1\phi(u)\phi(u)^T d u},
#' where \eqn{\phi} is a spline basis.
#' @param type specify spline basis, either \code{"b"} (default), 
#' \code{"c"}, \code{"i"} or \code{"m"};
#' @param splines_control control (see \code{\link{splines.control}}).
#' @details Integral is approximated using \code{\link[stats]{integrate}}.
#' @importFrom stats integrate
#' @importFrom splines2 bSpline cSpline iSpline mSpline
#' @export
P <- function(type = "b", splines_control = splines.control()){
  if(!type %in% c("b","c","i","m")) stop("Spline basis not supported")
  # splines.control 
  splines_control <- do.call("splines.control",splines_control)
  if(splines_control$derivs!=0) warning("Order of derivative is not 0") 
  
  # function to integrate
  f <- function(x,type,i,j,control){
    control$x <- x
    crossprod(do.call(paste0(type,"Spline"),control))[i,j]
  }
  
  # compute P
  k <- splines_control$df
  if(is.null(k)) k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
  mat <- matrix(nrow = k, ncol = k)
  for(i in 1:k){
    for(j in i:k){
      mat[i,j] <- mat[j,i] <- integrate(Vectorize(f, vectorize.args = "x"), 
                                        type = type, lower = 0, upper = 1, 
                                        i = i, j = j, control = splines_control)$value
    }
  }
  
  mat
}

# --------------
# negative log-likelihood (bivariate case)
# --------------
#' @title Negative log-likelihood of Spline Based Nonlinear Factor Copula Model
#' @description 
#' Compute the negative log-likelihood for a bivariate model. It is defined
#' as \deqn{c(w_1,w_2) = [\psi'(w_1)]^T \lambda^T P \lambda \psi'(w_2),}
#' where \eqn{\psi'} is the first derivative of a spline basis \eqn{\psi},
#' \eqn{\lambda} is a matrix of spline coefficients and \eqn{P} is 
#' some convenient matrix (see \link{P}). The copula models
#' positive quadrant dependent variables. If one has \eqn{N}-dimensional
#' variable, the variable is transformed to a 2-dimensional matrix (see 'Details').
#' @details 
#' If one inputs a matrix of \eqn{w_{i,1},\dots,w_{i,N}, i=1,\dots,n,} data points
#' throught \code{w1}, the data is internally transformed to a 
#' 2 columns matrix. If \eqn{N} is odd, the last column is dropped.
#' @param x vector spline coefficients (\eqn{\lambda} vectorized).
#' @param w1 uniform(0,1) vector of observations, or \code{matrix}
#' if \code{w2} is \code{NULL}.
#' @param w2 uniform(0,1) vector of observations, can be set to \code{NULL} (by default)
#' if \code{w1} is a \code{matrix}.
#' @param P the "P" matrix (see \code{\link{P}}), if \code{NULL} (by default) it is 
#' computed when calling the function.
#' @param type specify spline basis, either \code{"b"} (default), 
#' \code{"c"}, \code{"i"} or \code{"m"};
#' @param splines_control control (see \code{\link{splines.control}}).
#' @return The negative log-likelihood divided by \eqn{n}.
#' @export
nfcm_nll <- function(x,w1,w2=NULL,P=NULL,type="b",splines_control=splines.control()){
  # x         - a square matrix (not symmetric) vectorized of spline coefficients
  # w1        - sample 1 or a matrix of samples
  # w2        - sample 2 (optional)
  # P         - a square symmetric matrix of weights
  # type      - spline basis
  # splines_control   - splines control for spline
  splines_control <- do.call("splines.control",splines_control)
  
  # verification
  if(!type %in% c("b","c","i","m")) stop("Spline basis not supported")
  k <- splines_control$df
  if(is.null(k)) k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
  if(!is.null(w2)) if(length(w1) != length(w2)) stop("'w1' and 'w2' sizes mismatch")
  if(is.null(w2)) if(!is.matrix(w1)) stop("'w1' is not a matrix")
  if(!is.null(P)) if(nrow(P) != k || ncol(P) != k) stop("'P' has not the correct dimension")
  if(sqrt(length(x)) != k) stop("'x' has not the correct dimension")
  if(!is.null(P)) if(!identical(P,t(P))) stop("'P' is not symmetric")
  
  # rearrange observations (when 'w2' is NULL)
  if(is.null(w2)) {
    if(ncol(w1)==2){
      w2 <- w1[,2]
      w1 <- w1[,1]
    } else {
      warning("sample rearranged as a bivariate case")
      if(ncol(w1) %% 2 != 0) {
        warning("one column lost in conversion")
        tmp <- matrix(c(w1[,-ncol(w1)]), ncol = 2)
      } else {
        tmp <- matrix(c(w1), ncol = 2)
      }
      w1 <- tmp[,1]
      w2 <- tmp[,2]
    }
  }
  
  # compute C = x^T P x
  if(is.null(P)) P <- P(type, splines_control)
  x <- matrix(x, nrow = k, ncol = k)
  C <- crossprod(x, P) %*% x
  
  # generate 'a' and 'b'
  if(as.integer(splines_control$derivs)!=1L) splines_control$derivs <- 1L
  splines_control$x <- w1
  a <- do.call(paste0(type,"Spline"), splines_control)
  splines_control$x <- w2
  b <- do.call(paste0(type,"Spline"), splines_control)
  
  # compute the negative log-likelihood
  -mean(log(diag(tcrossprod(a %*% C, b))))
}

#' @rdname nfcm_nll
#' @return The gradient of the negative log-likelihood divided by \eqn{n}.
#' @export
nfcm_grad_nll <- function(x,w1,w2=NULL,P=NULL,type="b",splines_control=splines.control()){
  # x         - a square matrix (not symmetric) vectorized of spline coefficients
  # w1        - sample 1 or a matrix of samples
  # w2        - sample 2 (optional)
  # P         - a square symmetric matrix of weights
  # type      - spline basis
  # splines_control   - splines control for spline
  splines_control <- do.call("splines.control",splines_control)
  
  # verification
  if(!type %in% c("b","c","i","m")) stop("Spline basis not supported")
  k <- splines_control$df
  if(is.null(k)) k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
  if(!is.null(w2)) if(length(w1) != length(w2)) stop("'w1' and 'w2' sizes mismatch")
  if(is.null(w2)) if(!is.matrix(w1)) stop("'w1' is not a matrix")
  if(!is.null(P)) if(nrow(P) != k || ncol(P) != k) stop("'P' has not the correct dimension")
  if(sqrt(length(x)) != k) stop("'x' has not the correct dimension")
  if(!is.null(P)) if(!identical(P,t(P))) stop("'P' is not symmetric")
  
  # rearrange observations (when 'w2' is NULL)
  if(is.null(w2)) {
    if(ncol(w1)==2){
      w2 <- w1[,2]
      w1 <- w1[,1]
    } else {
      warning("sample rearranged as a bivariate case")
      if(ncol(w1) %% 2 != 0) {
        warning("one column lost in conversion")
        tmp <- matrix(c(w1[,-ncol(w1)]), ncol = 2)
      } else {
        tmp <- matrix(c(w1), ncol = 2)
      }
      w1 <- tmp[,1]
      w2 <- tmp[,2]
    }
  }
  
  # compute C = x^T P x
  x <- matrix(x, nrow = k, ncol = k)
  D <- crossprod(x, P)
  C <- D %*% x
  
  # generate 'a' and 'b'
  if(as.integer(splines_control$derivs)!=1L) splines_control$derivs <- 1L
  splines_control$x <- w1
  a <- do.call(paste0(type,"Spline"), splines_control)
  splines_control$x <- w2
  b <- do.call(paste0(type,"Spline"), splines_control)
  
  # compute the gradient of the negative log-likelihood
  tmp <- diag(tcrossprod(a %*% C, b))
  B <- matrix(0, nrow = k, ncol = k)
  n <- length(w1)
  for(i in seq_len(n)){
    mat <- tcrossprod(a[i,],b[i,]) / tmp[i]
    B <- B + mat + t(mat)
  }
  
  - c(crossprod(D, B)) / n
}

#' @rdname nfcm_nll
#' @return Akaike information criterion
#' @export
nfcm_aic <- function(x,w1,w2=NULL,P=NULL,type="b",splines_control=splines.control()){
  k <- length(x)
  if(is.null(w2)) n <- nrow(w1) else n <- length(w2)
  nfcm_nll(x=x,w1=w1,w2=w2,P=P,type=type,splines_control=splines_control) + 2.0 * k / n
}

#' @rdname nfcm_nll
#' @return Bayesian information criterion
#' @export
nfcm_bic <- function(x,w1,w2=NULL,P=NULL,type="b",splines_control=splines.control()){
  k <- length(x)
  if(is.null(w2)) n <- nrow(w1) else n <- length(w2)
  nfcm_nll(x=x,w1=w1,w2=w2,P=P,type=type,splines_control=splines_control) + log(n) * k / n
}