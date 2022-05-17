#---------------------------------
# Estimation of spline based nonlinear factor copula model
#---------------------------------
library(nfcm)
library(splines2)

# setting
spline_type <- "b"
assumed_model <- "clayton"
n <- 1000 # sample size
param <- list(alpha=1)
d <- 2 # number of dimension (min 2)
splines_control <- list(degree = 2L, Boundary.knots =  c(0,1), intercept=FALSE)
kn <- c(.1,.3,.7,.9)
splines_control$knots <- kn
MC <- 1e3
k <- 7

# load data
load(file=paste0(".simu/data/",assumed_model,"_",spline_type,"spline","_n_",n,"_d_",d,"_param_",param[[1]],".rds"))

# MISE
xx <- seq(0.01,1-0.01,by=.01)
# xx <- seq(1e-4,1-1e-4,length.out=1000)
splines_control$x <- xx
psi <- do.call(paste0(spline_type,"Spline"), splines_control)
z <- matrix(G(rep(xx,length(xx)),rep(xx,each=length(xx)), copula = assumed_model, param = param), nc = length(xx))
mise <- matrix(nrow = MC, nc = 3)
for(i in seq_len(MC)){
  z_hat <- psi %*% matrix(res$starting_value[i,], ncol = k) %*% t(psi)
  mise[i,1] <- mean((z - z_hat)^2)
  z_hat <- psi %*% matrix(res$spline_coef[i,], ncol = k) %*% t(psi)
  mise[i,2] <- mean((z - z_hat)^2)
  z_oracle <- psi %*% matrix(res$oracle[i,], ncol = k) %*% t(psi)
  mise[i,3] <- mean((z - z_oracle)^2)
}
colMeans(mise)
