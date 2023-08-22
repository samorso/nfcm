#---------------------------------
# Analysis of spline based nonlinear factor copula model
#---------------------------------
# different dependence
mise <- matrix(nr=9,nc=5)
colnames(mise) <- c("Oracle G", "Initial G", "MLE G", "Oracle H", "MLE H")
rownames(mise) <- paste0(rep(c("Clayton_", "Gumbel_", "Normal_"), each=3), rep(c(.2,.5,.8),3))
nll <- matrix(nr=9,nc=3)
colnames(nll) <- c("Oracle G", "Initial G", "MLE G")
rownames(nll) <- paste0(rep(c("Clayton_", "Gumbel_", "Normal_"), each=3), rep(c(.2,.5,.8),3))

load(".simu/data/clayton_bspline_n_250_d_20_param_0.5.rds")
mise[1,] <- colMeans(res$mise, na.rm=T)
nll[1,] <- colMeans(res$nll, na.rm=T)
load(".simu/data/clayton_bspline_n_250_d_20_param_2.rds")
mise[2,] <- colMeans(res$mise, na.rm=T)
nll[2,] <- colMeans(res$nll, na.rm=T)
load(".simu/data/clayton_bspline_n_250_d_20_param_5.rds")
mise[3,] <- colMeans(res$mise, na.rm=T)
nll[3,] <- colMeans(res$nll, na.rm=T)
load(".simu/data/gumbel_bspline_n_250_d_20_param_1.25.rds")
mise[4,] <- colMeans(res$mise, na.rm=T)
nll[4,] <- colMeans(res$nll, na.rm=T)
load(".simu/data/gumbel_bspline_n_250_d_20_param_2.rds")
mise[5,] <- colMeans(res$mise, na.rm=T)
nll[5,] <- colMeans(res$nll, na.rm=T)
load(".simu/data/gumbel_bspline_n_250_d_20_param_5.rds")
mise[6,] <- colMeans(res$mise, na.rm=T)
nll[6,] <- colMeans(res$nll, na.rm=T)
load(".simu/data/normal_bspline_n_250_d_20_param_0.3.rds")
mise[7,] <- colMeans(res$mise, na.rm=T)
nll[7,] <- colMeans(res$nll, na.rm=T)
load(".simu/data/normal_bspline_n_250_d_20_param_0.7.rds")
mise[8,] <- colMeans(res$mise, na.rm=T)
nll[8,] <- colMeans(res$nll, na.rm=T)
load(".simu/data/normal_bspline_n_250_d_20_param_0.95.rds")
mise[9,] <- colMeans(res$mise, na.rm=T)
nll[9,] <- colMeans(res$nll, na.rm=T)
mise
nll

# impact of dimension on estimator of u
d <- c(2,4,8,16,32,64)
MC <- 1e3
rmse <- matrix(nr = MC, ncol = length(d))
mise <- rep(NA_real_, length(d))
for(i in seq_along(d)){
  load(paste0(".simu/data/clayton_bspline_n_250_d_",d[i],"_param_2.rds"))
  rmse[,i] <- res$rmse[,1]
  mise[i] <- mean(res$mise[,5], na.rm=T)
}
boxplot(rmse)
mise
plot(x=d, y=colMeans(rmse,na.rm=T), type="b", pch = 16, xlab = "dimension", ylab = "Mean RMSE of U")

# impact of n on estimation of G
n <- c(250,500,1000,2000)
MC <- 1e3
mise <- nll <- matrix(nr = MC, ncol = length(n))
for(i in seq_along(n)){
  load(paste0(".simu/data/gumbel_bspline_n_",n[i],"_d_10_param_2.rds"))
  mise[,i] <- res$mise[,5]
  nll[,i] <- res$nll[,3]
}
boxplot(nll)
mise
plot(x=n, y=colMeans(mise,na.rm=T), type="b", pch = 16, xlab = "sample size", ylab = "MISE")
