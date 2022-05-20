#---------------------------------
# Estimation of spline based nonlinear factor copula model
#---------------------------------
library(nfcm)
library(splines2)

#------------------------
# Simulation setting
#------------------------
spline_type <- as.character(Sys.getenv("SPLINE"))
assumed_model <- as.character(Sys.getenv("MODEL"))
n <- as.integer(Sys.getenv("N")) # sample size
param <- as.list(na.omit(as.double(Sys.getenv(paste0("PARAM",1:5)))))
names(param) <- na.omit(as.character(Sys.getenv(paste0("NPARAM",1:5),unset=NA)))
d <- as.integer(Sys.getenv("D")) # number of dimension (min 2)
splines_control <- list(degree = 2L, Boundary.knots =  c(0,1), intercept=FALSE)
kn <- c(.1,.3,.7,.9)
splines_control$knots <- kn
k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
P <- P(type = spline_type, splines_control = splines_control)
MC <- 1e3

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
  spline_coef_H = matrix(nrow = MC, ncol = k * k),
  oracle_G = matrix(nrow = MC, ncol = k * k),
  oracle_H = matrix(nrow = MC, ncol = k * k),
  convergence = rep(NA_integer_, MC),
  iteration = rep(NA_integer_, MC),
  time = matrix(nrow = MC, ncol = 5),
  nll = matrix(nrow = MC, ncol = 3),
  mise = matrix(nrow = MC, ncol = 5),
  rmse = matrix(nrow = MC, ncol = 2)
)

# MISE
eps <- 1e-3
xx <- seq(eps,1-eps,length.out=1000)
splines_control$x <- xx
psi <- do.call(paste0(spline_type,"Spline"), splines_control)
splines_control$x <- NULL
z <- matrix(G(rep(xx,length(xx)),rep(xx,each=length(xx)), copula = assumed_model, param = param), nc = length(xx))
y <- matrix(H(1-rep(xx,length(xx)),rep(xx,each=length(xx)), copula = assumed_model, param = param), nc = length(xx))

#------------------------
# slurm specs
#------------------------
n_array <- 1000
ind <- matrix(seq_len(MC), nr=n_array, byr=TRUE)
id_slurm <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#------------------------
# Simulation
#------------------------
for(m in na.omit(ind[id_slurm,])){
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
    res$oracle_G[m,] <- fit_oracle_G
    res$time[m,1] <- difftime(t2, t1, units = "secs")
    res$nll[m,1] <- nfcm_nll(c(fit_oracle_G), w1 = w1, w2 = w2, P = P, type = spline_type, splines_control = splines_control)
    res$mise[m,1] <- mean((z - psi %*% matrix(fit_oracle_G, ncol = k) %*% t(psi))^2)
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
  if(is.null(fit_sv)) next
  x0 <- c(fit_sv)
  x0[x0<0] <- 0.0
  x0[x0>1] <- 1.0
  res$starting_value[m,] <- fit_sv
  res$nll[m,2] <- nfcm_nll(x0, w1 = w1, w2 = w2, P = P, type = spline_type, splines_control = splines_control)
  res$mise[m,2] <- mean((z - psi %*% matrix(fit_sv, ncol = k) %*% t(psi))^2)
  res$time[m,2] <- difftime(t2, t1, units = "secs")
  
  # Maximum (bivariate) likelihood estimator
  fit_G <- NULL
  t1 <- Sys.time()
  try(fit_G <- nfcm_mle(x = x0, w1 = w1, w2 = w2, P = P, type = spline_type, splines_control = splines_control), silent = TRUE)
  t2 <- Sys.time()
  if(is.null(fit_G)) next
  res$time[m,3] <- difftime(t2, t1, units = "secs")
  res$spline_coef_G[m,] <- fit_G$par
  res$convergence[m] <- fit_G$convergence
  res$iteration[m] <- fit_G$iter
  res$nll[m,3] <- nfcm_nll(c(fit_G$par), w1 = w1, w2 = w2, P = P, type = spline_type, splines_control = splines_control)
  res$mise[m,3] <- mean((z - psi %*% matrix(fit_G$par, ncol = k) %*% t(psi))^2)
  
  ### Estimation of H
  # Oracle for H
  fit_oracle_H <- NULL
  t1 <- Sys.time()
  try(fit_oracle_H <- lp_fit_Spline(u = 1.0 - u0, v = c(v0), w = c(w0), G=FALSE, type = spline_type, splines_control = splines_control), silent = TRUE)
  t2 <- Sys.time()
  if(!is.null(fit_oracle_H)){
    res$oracle_H[m,] <- fit_oracle_H
    res$time[m,4] <- difftime(t2, t1, units = "secs")
    res$mise[m,4] <- mean((y - psi %*% matrix(fit_oracle_H, ncol = k) %*% t(psi))^2)
  }
  
  # MAP
  u_hat <- map(x = c(fit_G$par), w = w0, type = spline_type, splines_control = splines_control)
  res$rmse[m,1] <- sqrt(sum((u0-u_hat)^2))
  splines_control$x <- u_hat
  C <- do.call(paste0(spline_type,"Spline"), splines_control) %*% matrix(fit_G$par, ncol = k)
  v_hat <- apply(w0, 2, function(x, control) {
    control$x <- x
    rowSums(do.call(paste0(spline_type,"Spline"), control) * C)
  }, control = splines_control)
  res$rmse[m,2] <- sqrt(sum((v0-v_hat)^2))
  splines_control$x <- NULL
  fit_H <- NULL
  t1 <- Sys.time()
  try(fit_H <- lp_fit_Spline(u = 1.0 - u_hat, v = c(v_hat), w = c(w0), G=FALSE, type = spline_type, splines_control = splines_control), silent = TRUE)
  t2 <- Sys.time()
  if(!is.null(fit_H)){
    res$spline_coef_H[m,] <- fit_H
    res$time[m,5] <- difftime(t2, t1, units = "secs")
    res$mise[m,5] <- mean((y - psi %*% matrix(fit_H, ncol = k) %*% t(psi))^2)
  }
  
  # Save results
  save(res, file=paste0("tmp/",assumed_model,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}