#---------------------------------
# Estimation of spline based nonlinear factor copula model
# Knots backward selection
#---------------------------------
library(nfcm)
spline_type <- "b"
assumed_model <- "clayton"
n <- 1000 # sample size
param <- as.list(0.5)
names(param) <- "alpha"
d <- 2 # number of dimension (min 2)
splines_control <- list(degree = 3L, Boundary.knots =  c(0,1), intercept=FALSE)

# data
set.seed(7322)
u <- runif(n)
v <- matrix(runif(n * d), ncol = 2)
w <- H(u, v, copula = assumed_model, param = param)

# initial values
kn0 <- seq(.05,.95,by=.05)
splines_control$knots <- kn0
x0 <- lp_fit_Spline(u = u, v = c(v), w = c(w), type = spline_type, splines_control = splines_control)
bic0 <- nfcm_bic(c(x0), w1 = w[,1], w2 = w[,2], type = spline_type, splines_control = splines_control)
bic1 <- bic0 - 1L
ind <- NULL
step <- 0L

while(bic1 < bic0){
  if(step>0) bic0 <- bic1
  bic <- rep(NA_real_, length(kn0))
  for(i in seq_along(kn0)){
    splines_control$knots <- kn0[-i]
    x0 <- lp_fit_Spline(u = u, v = c(v), w = c(w), type = spline_type, splines_control = splines_control)
    bic[i] <- nfcm_bic(c(x0), w1 = w[,1], w2 = w[,2], type = spline_type, splines_control = splines_control)
  }
  new_ind <- which.min(bic)
  ind <- c(ind, new_ind)
  if(min(bic, na.rm=TRUE)>=bic0) break
  bic1 <- bic[new_ind]
  kn0 <- kn0[-new_ind]
  step <- step + 1L
  cat(bic1,"\n")
}

# kn0 <- c(0.05, 0.65) after 17 steps


# initial values
kn0 <- seq(.05,.95,by=.05)
splines_control$knots <- kn0
x0 <- lp_fit_Spline(u = u, v = c(v), w = c(w), type = spline_type, splines_control = splines_control)
aic0 <- nfcm_aic(c(x0), w1 = w[,1], w2 = w[,2], type = spline_type, splines_control = splines_control)
aic1 <- aic0 - 1L
ind <- NULL
step <- 0L

while(aic1 < aic0){
  if(step>0) aic0 <- aic1
  aic <- rep(NA_real_, length(kn0))
  for(i in seq_along(kn0)){
    splines_control$knots <- kn0[-i]
    x0 <- lp_fit_Spline(u = u, v = c(v), w = c(w), type = spline_type, splines_control = splines_control)
    aic[i] <- nfcm_aic(c(x0), w1 = w[,1], w2 = w[,2], type = spline_type, splines_control = splines_control)
  }
  new_ind <- which.min(aic)
  ind <- c(ind, new_ind)
  if(min(aic, na.rm=TRUE)>=aic0) break
  aic1 <- aic[new_ind]
  kn0 <- kn0[-new_ind]
  step <- step + 1L
  cat(aic1,"\n")
}

# kn0 <- c(0.05, 0.25, 0.65, 0.85) after 15 steps
