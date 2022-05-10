# --------------
# library
# --------------
library(splines2)

# --------------
# Miscelleaneous
# --------------
# control for splines
# see ? splines2::iSpline
splines.control <- function(df = NULL, knots = seq(.2,.8,by=.2), degree = 3L, 
                            intercept = TRUE, Boundary.knots = c(0,1),
                            derivs = 0L){
  list(df=df, knots=knots, degree=degree, intercept=intercept, Boundary.knots=Boundary.knots, derivs=derivs)
}

# compute the matrix "P" (integral of crossprod of splines)
p_fct <- function(type = "i", splines_control = list()){
  if(!type %in% c("b","i")) stop("Support only 'b' or 'i' spline basis")
  # splines.control 
  splines_control <- do.call("splines.control",splines_control)

  # function to integrate
  if(type=="i")
    f <- function(x,i,j,control) crossprod(iSpline(x, df = NULL,
                                                   knots=control$knots,
                                                   degree=control$degree,
                                                   intercept=control$intercept,
                                                   Boundary.knots=control$Boundary.knots,
                                                   derivs = 0L))[i,j]
    
  if(type=="b")
    f <- function(x,i,j,control) crossprod(bSpline(x, df = NULL,
                                                   knots=control$knots,
                                                   degree=control$degree,
                                                   intercept=control$intercept,
                                                   Boundary.knots=control$Boundary.knots,
                                                   derivs = 0L))[i,j]

  # compute P
  k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
  P <- matrix(nrow = k, ncol = k)
  for(i in 1:k){
    for(j in i:k){
      P[i,j] <- P[j,i] <- integrate(Vectorize(f, vectorize.args = "x"), lower = 0, upper = 1, i = i, j = j, 
                                    control = splines_control)$value
    }
  }
  
  P
}

# --------------
# negative log-likelihood (bivariate case)
# --------------
neg_log_likelihood <- function(x,w1,w2,P,type="i",splines_control=list()){
  # x         - a square matrix (not symmetric) vectorized of spline coefficients
  # w1        - sample 1
  # w2        - sample 2
  # P         - a square symmetric matrix of weights
  # type      - spline basis
  # splines_control   - splines control for spline
  splines_control <- do.call("splines.control",splines_control)
  
  # verification
  if(!type %in% c("b","i")) stop("Support only 'b' or 'i' spline basis")
  k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
  if(length(w1) != length(w2)) stop("sample size should match")
  # if(!is.matrix(x)) stop("'x' should be a matrix")
  if(nrow(P) != k || ncol(P) != k) stop("'P' has not the correct dimension")
  # if(nrow(x) != k || ncol(x) != k) stop("'x' has not the correct dimension")
  if(sqrt(length(x)) != k) stop("'x' has not the correct dimension")
  if(!identical(P,t(P))) stop("'P' is not symmetric")
  
  # compute C = x^T P x
  x <- matrix(x, nrow = k, ncol = k)
  C <- crossprod(x, P) %*% x
  
  # generate 'a' and 'b'
  if(type == "i"){
    a <- iSpline(w1, knots=splines_control$knots, degree=splines_control$degree, 
                 intercept=splines_control$intercept, 
                 Boundary.knots=splines_control$Boundary.knots, derivs=1L)
    b <- iSpline(w2, knots=splines_control$knots, degree=splines_control$degree, 
                 intercept=splines_control$intercept, 
                 Boundary.knots=splines_control$Boundary.knots, derivs=1L)
  }
  if(type == "b"){
    a <- bSpline(w1, knots=splines_control$knots, degree=splines_control$degree, 
                 intercept=splines_control$intercept, 
                 Boundary.knots=splines_control$Boundary.knots, derivs=1L)
    b <- bSpline(w2, knots=splines_control$knots, degree=splines_control$degree, 
                 intercept=splines_control$intercept, 
                 Boundary.knots=splines_control$Boundary.knots, derivs=1L)
  }
  
  # compute the negative log-likelihood
  -mean(log(diag(tcrossprod(a %*% C, b))))
}

# gradient
grad_neg_log_likelihood <- function(x,w1,w2,P,type="i",splines_control=list()){
  # x         - a square matrix (not symmetric) vectorized of spline coefficients
  # w1        - sample 1
  # w2        - sample 2
  # P         - a square symmetric matrix of weights
  # type      - spline basis
  # splines_control   - splines_control for spline
  splines_control <- do.call("splines.control",splines_control)
  
  # verification
  if(!type %in% c("b","i")) stop("Support only 'b' or 'i' spline basis")
  k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
  if(length(w1) != length(w2)) stop("sample size should match")
  # if(!is.matrix(x)) stop("'x' should be a matrix")
  if(nrow(P) != k || ncol(P) != k) stop("'P' has not the correct dimension")
  # if(nrow(x) != k || ncol(x) != k) stop("'x' has not the correct dimension")
  if(sqrt(length(x)) != k) stop("'x' has not the correct dimension")
  if(!identical(P,t(P))) stop("'P' is not symmetric")

  # compute C = x^T P x
  x <- matrix(x, nrow = k, ncol = k)
  D <- crossprod(x, P)
  C <- D %*% x
  
  # generate 'a' and 'b'
  if(type=="i"){
    a <- iSpline(w1, knots=splines_control$knots, degree=splines_control$degree, 
                 intercept=splines_control$intercept, 
                 Boundary.knots=splines_control$Boundary.knots, derivs=1L)
    b <- iSpline(w2, knots=splines_control$knots, degree=splines_control$degree, 
                 intercept=splines_control$intercept, 
                 Boundary.knots=splines_control$Boundary.knots, derivs=1L)
  }
  if(type=="b"){
    a <- bSpline(w1, knots=splines_control$knots, degree=splines_control$degree, 
                 intercept=splines_control$intercept, 
                 Boundary.knots=splines_control$Boundary.knots, derivs=1L)
    b <- bSpline(w2, knots=splines_control$knots, degree=splines_control$degree, 
                 intercept=splines_control$intercept, 
                 Boundary.knots=splines_control$Boundary.knots, derivs=1L)
  }
  
  
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

# test
# set.seed(1)
# w1 <- runif(100)
# w2 <- runif(100)
# x <- runif(8 * 8)
# x <- x / sum(x)
# P <- p_fct()
# neg_log_likelihood(x,w1,w2,P)
# grad_neg_log_likelihood(x,w1,w2,P)
# 
# tol <- 1e-7
# grad <- matrix(nrow=8,ncol=8)
# for(i in 1:8){
#   for(j in 1:8){
#     eps <- matrix(0,nrow=8,ncol=8)
#     eps[i,j] <- tol
#     grad[i,j] <- (neg_log_likelihood(x+c(eps),w1,w2,P) - neg_log_likelihood(x-c(eps),w1,w2,P)) / tol / 2.0
#   }
# }
# grad
# abs(c(grad) - grad_neg_log_likelihood(x,w1,w2,P))

# --------------
# MAP estimator (bivariate case)
# --------------
map <- function(x,w1,w2,P,type="i",splines_control=list()){
  # x         - a square matrix (not symmetric) of spline coefficients
  # w1        - sample 1
  # w2        - sample 2
  # P         - a square symmetric matrix of weights
  # type      - spline basis
  # splines_control   - splines control for spline
  splines_control <- do.call("splines.control",splines_control)
  
  # verification
  if(!type %in% c("b","i")) stop("Support only 'b' or 'i' spline basis")
  k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
  n <- length(w1)
  if(n != length(w2)) stop("sample size should match")
  if(!is.matrix(x)) stop("'x' should be a matrix")
  if(nrow(P) != k || ncol(P) != k) stop("'P' has not the correct dimension")
  # if(nrow(x) != k || ncol(x) != k) stop("'x' has not the correct dimension")
  if(length(x) != length(P)) stop("the dimension of 'x' and 'P' should match")
  if(!identical(P,t(P))) stop("'P' is not symmetric")
  
  # compute C = x^T P x
  C <- crossprod(x, P) %*% x
  
  # generate 'a' and 'b'
  if(type == "i"){
    a <- iSpline(w1, knots=splines_control$knots, degree=splines_control$degree, 
                 intercept=splines_control$intercept, 
                 Boundary.knots=splines_control$Boundary.knots, derivs=1L)
    b <- iSpline(w2, knots=splines_control$knots, degree=splines_control$degree, 
                 intercept=splines_control$intercept, 
                 Boundary.knots=splines_control$Boundary.knots, derivs=1L)
  }
  if(type == "b"){
    a <- bSpline(w1, knots=splines_control$knots, degree=splines_control$degree, 
                 intercept=splines_control$intercept, 
                 Boundary.knots=splines_control$Boundary.knots, derivs=1L)
    b <- bSpline(w2, knots=splines_control$knots, degree=splines_control$degree, 
                 intercept=splines_control$intercept, 
                 Boundary.knots=splines_control$Boundary.knots, derivs=1L)
  }
  y <- diag(tcrossprod(a %*% C, b))
  
  # conditional distribution of u given w1, w2
  f <- function(u,x,a,b,y,type="i",splines_control){
    if(type == "i")
      phi2 <- crossprod(iSpline(u, knots=splines_control$knots,
                                degree=splines_control$degree, 
                                intercept=splines_control$intercept, 
                                Boundary.knots=splines_control$Boundary.knots, 
                                derivs=0L))
    if(type == "b")
      phi2 <- crossprod(bSpline(u, knots=splines_control$knots,
                                degree=splines_control$degree, 
                                intercept=splines_control$intercept, 
                                Boundary.knots=splines_control$Boundary.knots, 
                                derivs=0L))
    -tcrossprod(a %*% crossprod(x, phi2) %*% x, b) / y
  }
  
  # MAP estimator
  u <- rep(NA_real_, n)
  for(i in seq_len(n)){
    u[i] <- optimize(f, interval = c(0,1), x=x, a = a[i,], b = b[i,], y = y[i], 
                     type = type, splines_control = splines_control)$minimum
  }
  u
}
