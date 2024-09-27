# --------------
# LP estimator (assuming everything is observable)
# --------------
#' @title Linear program to estimate spline coefficients
#' @description 
#' Estimate the spline coefficients of the tensor product spline
#' by minimizing
#' \deqn{\sum_{i=1}^N|x_i-\phi(u)^TX\psi(y_i)|,}
#' under monotonicity constraint (non-negative partial derivatives).
#' @param u vector of latent uniform(0,1) variable.
#' @param v uniform(0,1) vector of error.
#' @param w uniform(0,1) vector of observation.
#' @param G if \code{TRUE} estimate function \code{G} (\code{\link{G}}), 
#' otherwise estimate function \code{H} (\code{\link{H}}) (see 'Details').
#' @param type specify spline basis, either \code{"b"} (default), 
#' \code{"c"}, \code{"i"} or \code{"m"};
#' @param splines_control control (see \code{\link{splines.control}}).
#' @return A matrix of spline coefficients
#' @importFrom slam as.simple_triplet_matrix simple_triplet_matrix
#' @importFrom Rglpk Rglpk_solve_LP
#' @importFrom stats deriv
#' @export
lp_fit_Spline <- function(u, v, w, G=TRUE, type="b", splines_control=splines.control()){
  # u - common latent variable from unif(0,1)
  # v - unif(0,1) errors
  # w - observations from unif(0,1)
  # splines_control - list of control for Spline()
  # Note: G(u,v) is suppose monotonic non-decreasing in both arguments. 
  # If dim(u) < dim(v), it means that dim(u) = n and dim(v) = n * N.
  # linear problem is of the form:
  # min C^T\xi under constraint A\xi <= b
  # where \xi = (vec(splines coefficients), vec(data))^T
  
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
  
  # Constraints:
  p <- ncol(phi)
  p2 <- p * p
  # C <- c(rep(0,p2),rep(1,N)) 
  C <- simple_triplet_matrix(i = p2 + seq_len(N), j = rep(1,N), v = rep(1,N), nrow = p2 + N, ncol = 1)
  pp <- matrix(rep(phi,p), nrow = N) * matrix(rep(t(psi),each=p), nrow = N, byrow=TRUE)
  
  # "regular" constrain
  A1 <- cbind(rbind(pp,-pp)[c(rbind(seq(N),seq(N)+N)),],
              diag(N) %x% rbind(-1,-1))
  if(G) b1 <- cbind(c(rbind(v,-v)))
  if(!G) b1 <- cbind(c(rbind(w,-w)))
  
  # specific constrain:
  # monotonic constrain for B-splines (increasing in both direction)
  if(type=="b"){
    ind <- matrix(1:p2, nrow = p)
    tmp <- matrix(0, nrow = 2 * p * (p - 1), ncol = p2)
    it <- 0L
    for (i in seq_len(p)) {
      for (j in seq_len(p - 1)) {
        it <- it + 1L
        tmp[it, ind[i, j]] <- 1
        tmp[it, ind[i, j + 1]] <- -1
        tmp[it + p * (p - 1), ind[j, i]] <- 1
        tmp[it + p * (p - 1), ind[j + 1, i]] <- -1
      }
    }
    A2 <- cbind(tmp,matrix(0, nrow = 2 * p * (p - 1), ncol = N))
    b2 <- cbind(rep(0, 2 * p * (p - 1)))
  }
  
  # specific constrain:
  # equality constrain for I-splines and C-splines
  # (coefficients sum up to 1)
  if(type=="i" || type=="c"){
    # A2 <- rbind(A2, c(rep(1,p2),rep(0,N)))
    # b2 <- rbind(b2, 1.0)
    A2 <- c(rep(1,p2),rep(0,N))
    b2 <- 1.0
  }
  
  # overall constrains
  # every coefficient lies in R (by default they are assumed to be positive)
  # additionally, they are assumed to lie in [0,1]
  bounds <- list(upper = list(ind=seq_len(p2), val = rep(1.0, p2)))
  # if(type=="b"){ bounds <- list(upper = list(ind=seq_len(p2), val = rep(1.0, p2))) }else {
  #   bounds <- list(lower = list(ind=seq_len(p2), val = rep(-Inf, p2)))
  # }
  if(type != "m" ) A <- rbind(A1,A2) else A <- A1
  if(type != "m") b <- rbind(b1,b2) else b <- b1
  f_dir <- rep("<=",nrow(A))
  if(type == "i" || type == "c") f_dir[length(f_dir)] <- "=="
  A <- as.simple_triplet_matrix(A)
  # C <- as.simple_triplet_matrix(C)
  fit_lp <- Rglpk_solve_LP(obj = C, mat = A, dir = f_dir, rhs = b, bounds = bounds)
  res <- list(fit=fit_lp, par = matrix(fit_lp$solution[seq_len(p2)],nrow=p), control=splines_control, type=type)
  if(G) structure(res, class = "G.spline") else structure(res, class = "H.spline")
}

# --------------
# Maximum likelihood (bivariate) estimator
# --------------
#' @title Control for fitting the MLE
#' @description Control for the optimization routine
#' for the MLE (see \code{\link{nfcm_mle}}). Currently the
#' optimization routine is a sequential quadratic programming
#' (SQP) algorithm \code{\link[nloptr]{slsqp}} from the \code{nloptr} 
#' package, and thus the control are the same as \code{\link[nloptr]{nl.opts}}
#' @param stopval stop minimization at this value
#' @param xtol_rel tolerance for step
#' @param maxeval max number of function evaluations (default is 30000)
#' @param ftol_rel relative tolerance
#' @param ftol_abs absolute tolerance
#' @param check_derivatives verification
#' @seealso \code{\link[nloptr]{nl.opts}}, \code{\link[nloptr]{slsqp}}
#' @export
mle.control <- function(stopval = -Inf, xtol_rel = 1e-6, maxeval = 30000,
                        ftol_rel = 0.0, ftol_abs = 0.0, check_derivatives = FALSE){
  if(!is.numeric(stopval)) stop("'stopval' must be numeric")
  if(!is.numeric(xtol_rel)) stop("'xtol_rel' must be numeric")
  if(as.integer(maxeval)<=0) stop("'maxeval' must be positive")
  if(!is.numeric(ftol_rel)) stop("'ftol_rel' must be numeric")
  if(!is.numeric(ftol_abs)) stop("'ftol_abs' must be numeric")
  if(!is.logical(check_derivatives)) stop("'check_derivatives' must be a boolean")
  
  list(stopval=stopval, xtol_rel=xtol_rel, maxeval=maxeval, ftol_rel=ftol_rel, ftol_abs=ftol_abs, check_derivatives=check_derivatives)
}


#' @title Maximum likelihood estimator of Spline Based Nonlinear Factor Copula Model
#' @description 
#' Minimize the negative log-likelihood for a bivariate model. The copula models assume
#' positive quadrant dependent variables. If one has \eqn{N}-dimensional
#' variable, the variable is transformed to a 2-dimensional matrix (see 'Details').
#' @param lambda vector of starting values for spline coefficients (\eqn{\Lambda} vectorized).
#' @param w uniform(0,1) \eqn{n\times N} \code{matrix} of observations.
#' @param type specify spline basis, either \code{"b"} (default), 
#' \code{"c"}, \code{"i"} or \code{"m"};
#' @param splines_control control (see \code{\link{splines.control}}).
#' @param mle_control control to pass to the optimization routine (see 'Details').
#' @return An object from the optimization routine (see 'Details')
#' @importFrom nloptr slsqp
#' @seealso \code{\link[nloptr]{slsqp}}
#' @export
nfcm_mle <- function(lambda, w, type="b", splines_control=splines.control(), mle_control=mle.control()){
  # lambda            - a square matrix (not symmetric) vectorized of spline coefficients
  # w                 - a matrix of samples
  # type              - spline basis
  # splines_control   - splines control for spline
  # mle_control   - control for optimization routine
  splines_control <- do.call("splines.control",splines_control)
  mle_control <- do.call("mle.control", mle_control)
  
  # verification
  if(!type %in% c("b","c","i","m")) stop("Spline basis not supported")
  k <- splines_control$df
  if(is.null(k)) k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
  if(!type %in% c("b","c","i","m")) stop("Spline basis not supported")
  k <- splines_control$df
  if(is.null(k)) k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
  if(!is.matrix(w)) stop("'w' is not a matrix")
  if(sqrt(length(lambda)) != k) stop("'lambda' has not the correct dimension")
  if(any(lambda<0) || any(lambda>1)) stop("'lambda' must be between 0 and 1")
  
  # inequality constraint
  # inequality constraint (for B-spline)
  hin <- hinjac <- NULL
  if(type == "b") {
    hin <- function(x) {
      p2 <- length(x)
      p <- sqrt(p2)
      ind <- matrix(seq_len(p2), nrow = p)
      tmp <- matrix(0, nrow = 2 * p * (p - 1), ncol = p2)
      it <- 0L
      for (i in seq_len(p)) {
        for (j in seq_len(p - 1)) {
          it <- it + 1L
          tmp[it, ind[i, j]] <- 1
          tmp[it, ind[i, j + 1]] <- -1
          tmp[it + p * (p - 1), ind[j, i]] <- 1
          tmp[it + p * (p - 1), ind[j + 1, i]] <- -1
        }
      }
      tmp %*% x
    }
    
    hinjac <- function(x) {
      p2 <- length(x)
      p <- sqrt(p2)
      ind <- matrix(seq_len(p2), nrow = p)
      tmp <- matrix(0, nrow = 2 * p * (p - 1), ncol = p2)
      it <- 0L
      for (i in seq_len(p)) {
        for (j in seq_len(p - 1)) {
          it <- it + 1L
          tmp[it, ind[i, j]] <- 1
          tmp[it, ind[i, j + 1]] <- -1
          tmp[it + p * (p - 1), ind[j, i]] <- 1
          tmp[it + p * (p - 1), ind[j + 1, i]] <- -1
        }
      }
      tmp
    }
  }
  
  # equality constraint (for I-spline and C-spline)
  heq <- heqjac <- NULL
  if(type == "i" || type == "c"){
    heq <- function(x){
      sum(x) - 1.0
    }
    heqjac <- function(x){
      rep(1.0, length(x))
    }
  }
  
  # fit the MLE
  fit <- slsqp(x0 = lambda, fn = nfcm_nll, gr = nfcm_grad_nll, 
               lower = rep(0.0, length(lambda)), upper = rep(1.0, length(lambda)),
               hin = hin, hinjac = hinjac, heq = heq, heqjac = heqjac, control = mle_control,
               deprecatedBehavior = FALSE,
               w = w, type = type, splines_control = splines_control)
  res <- list(fit=fit, par = matrix(fit$par, nrow=k), control=splines_control, type=type)
  structure(res, class = "G.spline")
}

# --------------
# Prediction of the latent variable
# --------------
#' @title Predict latent variable
#' @description 
#' Estimate the latent variable \eqn{u}.
#' @param lambda vector spline coefficients (\eqn{\Lambda} vectorized).
#' @param w uniform(0,1) \code{matrix} of observations.
#' @param type specify spline basis, either \code{"b"} (default), 
#' \code{"c"}, \code{"i"} or \code{"m"};
#' @param splines_control control (see \code{\link{splines.control}}).
#' @param method method to estimate the latent variable, either \code{"map"} (default) or \code{"blup"}.
#' @return Estimates for \eqn{u}, return a vector if \code{method="map"} and a matrix if \code{method="blup"}.
#' @importFrom stats optimize
#' @export
predict_u <- function(lambda,w,type="b",splines_control=splines.control(), method="map"){
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
  if(!method %in% c("map","blup")) stop("Method not supported")
  
  # set coefficients as matrix
  Lambda <- matrix(lambda, ncol = k)

  # generate psi
  derivs <- as.integer(splines_control$derivs)
  if(derivs!=1L) splines_control$derivs <- 1L
  n <- nrow(w)
  N <- ncol(w)
  psip <- array(dim = c(n,k,N))
  for(j in 1:N){
    splines_control$x <- w[,j]
    psip[,,j] <- do.call(paste0(type,"Spline"), splines_control)
  }
  
  # reset derivs to 0
  splines_control$derivs <- 0L
  
  
  # MAP estimator
  if(method == "map"){
    # predictions
    u <- rep(NA_real_, n)
    
    # conditional distribution of u given w1, w2
    f <- function(u,Lambda,psip,type="b",splines_control){
      splines_control$x <- u
      phi <- do.call(paste0(type,"Spline"), splines_control)
      -prod(phi %*% Lambda %*% psip)
    }
    
    for(i in seq_len(n)){
      u[i] <- optimize(f, interval = c(0,1), Lambda=Lambda, psip = psip[i,,], 
                       type = type, splines_control = splines_control)$minimum
    }
  }
  
  # BLUP estimator
  if(method == "blup"){
    # predictions
    u <- matrix(nrow=n, ncol=2)
    
    # compute integrals
    q_vec <- q_vec(type = type, splines_control = splines_control)
    if(!isTRUE(splines_control$integral)) splines_control$integral <- TRUE
    splines_control$x <- 1
    Phi <- do.call(paste0(type,"Spline"), splines_control)
    C1 <- q_vec %*% Lambda
    C2 <- Phi %*% Lambda
    
    for(i in seq_len(n)) {
      u[i,1] <- prod(C1 %*% psip[i,,]) / prod((C2 %*% psip[i,,])^(1 / N))
      u[i,2] <- prod(C1 %*% psip[i,,])^(1/N) / prod(C2 %*% psip[i,,])
    }
        
  }
  
  u
}


# --------------
# Integral of the spline basis
# --------------
# Approximation of the integral of \int_0^1 u\phi(u)du
#' @importFrom stats integrate
q_vec <- function(type = "b", splines_control = splines.control()){
  if(!type %in% c("b","c","i","m")) stop("Spline basis not supported")
  # splines.control 
  if(!is.null(splines_control$x)) splines_control$x <- NULL
  splines_control <- do.call("splines.control",splines_control)
  if(splines_control$derivs!=0) warning("Order of derivative is not 0") 
  
  # function to integrate
  f <- function(x,type,i,control){
    control$x <- x
    x * do.call(paste0(type,"Spline"),control)[i]
  }
  
  # compute p_vec
  k <- splines_control$df
  if(is.null(k)) k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
  vec <- rep(NA_real_, k)
  for(i in 1:k){
    vec[i] <- integrate(Vectorize(f, vectorize.args = "x"), lower = 0, upper = 1, type = type, i = i, control = splines_control)$value
  }
  
  vec
}
