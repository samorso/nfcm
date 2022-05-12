#---------------------------------
# Estimation of spline based nonlinear factor copula model
#---------------------------------
library(nfcm)
library(splines2)
library(nloptr)

#------------------------
# Simulation setting
#------------------------
spline_type <- as.character(Sys.getenv("SPLINE"))
assumed_model <- as.character(Sys.getenv("MODEL"))
n <- as.integer(Sys.getenv("N")) # sample size
param <- as.list(na.omit(as.double(Sys.getenv(paste0("PARAM",1:5)))))
names(param) <- na.omit(as.character(Sys.getenv(paste0("NPARAM",1:5),unset=NA)))
d <- as.integer(Sys.getenv("D")) # number of dimension (min 2)
splines_control <- list(degree = 3L, Boundary.knots =  c(0,1), intercept=FALSE)
kn <- c(0.05, 0.25, 0.65, 0.85)
splines_control$knots <- kn
k <- length(splines_control$knots) + splines_control$degree + as.integer(splines_control$intercept)
P <- P(type = spline_type, splines_control = splines_control)
MC <- 1e3

# inequality constraint
hin <- function(x){
  tol <- 1e-10
  knots <- c(0.05, 0.25, 0.65, 0.85)
  xx <- c(knots, 1.0)
  lambda <- matrix(x, nrow = 7)
  psi <- bSpline(xx, knots = knots,
                 degree = 3L,
                 intercept = FALSE,
                 Boundary.knots = c(0,1))
  d_psi <- deriv(psi)
  c(c(tcrossprod(psi %*% lambda, d_psi)), c(tcrossprod(d_psi %*% lambda, psi))) + tol
}
hinjac <- function(x){
  knots <- c(0.05, 0.25, 0.65, 0.85)
  xx <- c(knots, 1.0)
  psi <- bSpline(xx, knots = knots,
                 degree = 3L,
                 intercept = FALSE,
                 Boundary.knots = c(0,1))
  d_psi <- deriv(psi)
  mat <- matrix(nrow = length(xx) * length(xx) * 2, ncol = length(x))
  id <- 0L
  for(i in seq_along(xx)){
    for(j in seq_along(xx)){
      id <- id + 1L
      mat[id,] <- c(tcrossprod(psi[i,],d_psi[j,]))
    }
  }
  for(i in seq_along(xx)){
    for(j in seq_along(xx)){
      id <- id + 1L
      mat[id,] <- c(tcrossprod(d_psi[i,],psi[j,]))
    }
  }
  mat
}

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
  spline_coef = matrix(nrow = MC, ncol = k * k),
  oracle = matrix(nrow = MC, ncol = k * k),
  convergence = rep(NA_integer_, MC),
  iteration = rep(NA_integer_, MC),
  time = matrix(nrow = MC, ncol = 3),
  nll = matrix(nrow = MC, ncol = 3)
)

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
  
  # Oracle
  fit_oracle <- NULL
  t1 <- Sys.time()
  try(fit_oracle <- lp_fit_Spline(u = u0, v = c(v0), w = c(w0), type = spline_type, splines_control = splines_control), silent = TRUE)
  t2 <- Sys.time()
  if(!is.null(fit_oracle)){
    res$oracle[m,] <- fit_oracle
    res$time[m,1] <- difftime(t2, t1, units = "secs")
  }
  
  # Find starting values (based on a normal copula)
  w1 <- matrix(c(w0),ncol=2)[,1]
  w2 <- matrix(c(w0),ncol=2)[,2]
  est_cor <- cor(w1, w2)
  set.seed(se2[m])
  u <- runif(n)
  v <- matrix(runif(n * d), ncol = 2)
  w <- H(u, v, copula = "normal", param = list(corr = est_cor))
  fit_sv <- NULL
  try(fit_sv <- lp_fit_Spline(u = u, v = c(v), w = c(w), type = spline_type, splines_control = splines_control), silent = TRUE)
  if(is.null(fit_sv)) next
  res$starting_value[m,] <- fit_sv
  res$time[m,2] <- difftime(t2, t1, units = "secs")
  
  # starting values
  x0 <- c(fit_sv)
  x0[x0<0] <- 0.0
  x0[x0>1] <- 1.0
  
  t1 <- Sys.time()
  fit <- slsqp(x0 = x0, fn = nfcm_nll, gr = nfcm_grad_nll, 
                lower = rep(0.0, length(x0)), upper = rep(1.0, length(x0)),
                hin = hin, hinjac = hinjac, control = list(maxeval = 20000),
                w1 = w1, w2 = w2, P = P, type = spline_type, splines_control = splines_control)
  t2 <- Sys.time()
  res$time[m,3] <- difftime(t2, t1, units = "secs")
  res$spline_coef[m,] <- fit$par
  res$convergence[m] <- fit$convergence
  res$iteration[m] <- fit$iter
  if(!is.null(fit_oracle)) res$nll[m,1] <- nfcm_nll(c(fit_oracle), w1 = w1, w2 = w2, P = P, type = spline_type, splines_control = splines_control)
  res$nll[m,2] <- nfcm_nll(x0, w1 = w1, w2 = w2, P = P, type = spline_type, splines_control = splines_control)
  res$nll[m,3] <- nfcm_nll(c(fit$par), w1 = w1, w2 = w2, P = P, type = spline_type, splines_control = splines_control)
  
  
  # Save results
  save(res, file=paste0("tmp/",assumed_model,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}