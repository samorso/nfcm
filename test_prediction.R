#---------------------------------
# Estimation of spline based nonlinear factor copula model
#---------------------------------
library(nfcm)
library(splines2)

#------------------------
# Simulation setting
#------------------------
spline_type <- "b"
assumed_model <- "normal"
n <- 1000
param <- list(0.5)
d <- 2 # number of dimension (min 2)
splines_control <- list(degree = 3L, Boundary.knots =  c(0,1), intercept=FALSE)
kn <- c(.1,.3,.7,.9)
splines_control$knots <- kn
k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
P <- P(type = spline_type, splines_control = splines_control)
MC <- 100

# seeds for random generator
set.seed(round(1e3 * param[[1]]))
random_number1 <- sample.int(1e6,1)
set.seed(643 * n)
random_number2 <- sample.int(1e6,1)
set.seed(7483 * d)
random_number3 <- sample.int(1e6,1)
set.seed(6543 + random_number1 + random_number2 + random_number3) # random hand typing algorithm
se1 <- sample.int(1e7, MC)
se2 <- sample.int(1e7, MC)
res <- list(
  starting_value = matrix(nrow = MC, ncol = k * k),
  spline_coef_G = matrix(nrow = MC, ncol = k * k),
  oracle_G = matrix(nrow = MC, ncol = k * k),
  time = matrix(nrow = MC, ncol = 4),
  u_blup = matrix(nrow = MC, ncol = n),
  u_map = matrix(nrow = MC, ncol = n),
  mse = matrix(nrow = MC, ncol = 2)
)

#------------------------
# Simulation
#------------------------
for(m in 1:MC){
  # Sample 
  set.seed(se1[m])
  u0 <- runif(n)
  v0 <- matrix(runif(n * d), ncol = d)
  w0 <- H(u0, v0, copula = assumed_model, param = param)
  # transform to bivariate case
  w1 <- matrix(c(w0),ncol=2)[,1]
  w2 <- matrix(c(w0),ncol=2)[,2]
  
  ### Estimation of G
  # Oracle for G
  fit_oracle_G <- NULL
  t1 <- Sys.time()
  try(fit_oracle_G <- lp_fit_Spline(u = u0, v = c(v0), w = c(w0), type = spline_type, splines_control = splines_control), silent = TRUE)
  t2 <- Sys.time()
  if(!is.null(fit_oracle_G)){
    res$oracle_G[m,] <- fit_oracle_G$par
    res$time[m,1] <- difftime(t2, t1, units = "secs")
  }
  
  # Find starting values (based on a normal copula)
  est_cor <- cor(w1, w2)
  set.seed(se2[m])
  u <- runif(n)
  v <- matrix(runif(n * d), ncol = 2)
  w <- H(u, v, copula = "normal", param = list(corr = est_cor))
  fit_sv <- NULL
  t1 <- Sys.time()
  try(fit_sv <- lp_fit_Spline(u = u, v = c(v), w = c(w), type = spline_type, splines_control = splines_control), silent = TRUE)
  t2 <- Sys.time()
  starting_value <- c(fit_sv$par)
  
  
  # Maximum (bivariate) likelihood estimator
  fit_G <- NULL
  t1 <- Sys.time()
  try(fit_G <- nfcm_mle(x = starting_value, w1 = w1, w2 = w2, P = P, type = spline_type, splines_control = splines_control), silent = TRUE)
  t2 <- Sys.time()
  if(is.null(fit_G)) next
  res$time[m,2] <- difftime(t2, t1, units = "secs")
  res$spline_coef_G[m,] <- fit_G$par
  
  ### Prediction of U
  # MAP
  t1 <- Sys.time()
  res$u_map[m,] <- predict_u(x = c(fit_G$par), w = w0, type = spline_type, splines_control = splines_control, method = "map")
  t2 <- Sys.time()
  res$time[m,3] <- difftime(t2, t1, units = "secs")
  res$mse[m,1] <- mean((res$u_map[m,] - u0)^2)
  
  # BLUP
  t1 <- Sys.time()
  res$u_blup[m,] <- predict_u(x = c(fit_G$par), w = w0, type = spline_type, splines_control = splines_control, method = "blup")
  t2 <- Sys.time()
  res$time[m,4] <- difftime(t2, t1, units = "secs")
  res$mse[m,2] <- mean((res$u_blup[m,] - u0)^2)
  
  cat(m, "\n")
}