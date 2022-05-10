# --------------
# library
# --------------
library(nloptr)
source("R/loglikelihood.R")
source("R/lp_estimator.R")

# --------------
# DGP
# --------------
H_Clayton <- function(u,v,alpha) (0.1e1 - log(v) / qgamma(1-u,shape=alpha))^(-alpha)
G_Clayton <- function(u,w,alpha) exp((1.0 - w^(-1/alpha)) * qgamma(1-u, shape = alpha))
H_Normal <- function(u,v,cor) pnorm(sqrt(cor) * qnorm(1-u) + sqrt(1-cor) * qnorm(v))
G_Normal <- function(u,w,cor) pnorm(qnorm(w) / sqrt(1 - cor) - sqrt(cor / (1 - cor)) * qnrom(1 - u))
alpha <- .5 # true parameter for Clayton H
n <- 1000 # sample size
d <- 2 # dimension 
set.seed(125)
u <- runif(n)
v <- matrix(runif(n*d),nc=d)
w <- apply(v,2,function(x)H_Clayton(u,x,alpha))

# --------------
# Fit comparison iSpline vs bSpline
# --------------
# c(.2,.5,.8), c(.04,.08,.2,.5,.8,.9,.95), c(.02,.04,.06,.08,seq(.1,.9,by=.1),.92,.94,.96,.98)
knots <- c(.04,.08,.2,.5,.8,.9,.95)
bsplines_control <- list(knots = knots, degree = 3L, Boundary.knots =  c(0,1), intercept=FALSE)
isplines_control <- list(knots = knots, degree = 3L, Boundary.knots =  c(0,1), intercept=FALSE)
# "ideal" fit
fit_ispline <- lp_fit_iSpline(u,c(v),c(w),isplines_control)
fit_bspline <- lp_fit_bSpline(u,c(v),c(w),bsplines_control)
# approximative fit with Normal copula
set.seed(1)
u2 <- runif(n)
v2 <- matrix(runif(n*d),nc=d)
w2 <- apply(v2,2,function(x)H_Normal(u2,x,cor(w)[1,2]))
fit_bspline2 <- lp_fit_bSpline(u2,c(v2),c(w2),bsplines_control)

# check partial derivatives
psi <- iSpline(c(isplines_control$knots,1),
               intercept = isplines_control$intercept,
               knots=isplines_control$knots,
               degree=isplines_control$degree,
               Boundary.knots=isplines_control$Boundary.knots)
d_psi <- deriv(psi)
psi %*% fit_ispline %*% t(d_psi)
d_psi %*% fit_ispline %*% t(psi)

psi <- bSpline(bsplines_control$knots,
               intercept = bsplines_control$intercept,
               knots=bsplines_control$knots,
               degree=bsplines_control$degree,
               Boundary.knots=bsplines_control$Boundary.knots)
d_psi <- bSpline(bsplines_control$knots,
               intercept = bsplines_control$intercept,
               knots=bsplines_control$knots,
               degree=bsplines_control$degree,
               Boundary.knots=bsplines_control$Boundary.knots,
               derivs = 1L)
psi %*% fit_bspline %*% t(d_psi)
d_psi %*% fit_bspline %*% t(psi)

# check the fit
angle <- 250
xx <- seq(0.01,1-0.01,by=.01)
z <- matrix(G_Clayton(rep(xx,length(xx)),rep(xx,each=length(xx)),alpha),nc=length(xx))
persp(z, theta=angle,xlab="u",ylab="w",zlab="G",main="True function")
# isplines_control <- splines.control()
a <- iSpline(xx,knots=isplines_control$knots,
             degree=isplines_control$degree,intercept=isplines_control$intercept,
             Boundary.knots=isplines_control$Boundary.knots)
b <- iSpline(rev(xx),knots=isplines_control$knots,
             degree=isplines_control$degree,intercept=isplines_control$intercept,
             Boundary.knots=isplines_control$Boundary.knots)
persp(a %*% fit_ispline %*% t(a),theta=angle,xlab="u",ylab="w",zlab="G",main="iSpline")
a2 <- bSpline(xx,knots=bsplines_control$knots,
             degree=bsplines_control$degree,
             Boundary.knots=bsplines_control$Boundary.knots)
b2 <- bSpline(rev(xx),knots=bsplines_control$knots,
             degree=bsplines_control$degree,
             Boundary.knots=bsplines_control$Boundary.knots)
# persp(b2 %*% fit_bspline %*% t(a2),theta=-25,xlab="u",ylab="w",zlab="G",main="bSpline")
persp(a2 %*% fit_bspline2 %*% t(a2),theta=angle,xlab="u",ylab="w",zlab="G",main="bSpline2")

persp(a2 %*% fit_bspline %*% t(a2),theta=0,xlab="u",ylab="w",zlab="G",main="bSpline2")
plot(u,w[,1])
plot(u,w[,2])

# --------------
# Estimation (iSplines)
# --------------

# equality constraint
heq <- function(x){
  sum(x) - 1.0
}
heqjac <- function(x){
  rep(1.0, length(x))
}

# inequality constraint
hin <- function(x){
  # xx <- seq(.05,.95,by=.05)
  xx <- seq(.1,.9,by=.1)
  lambda <- matrix(x, nrow = 6)
  psi <- iSpline(xx, knots = c(0.2,0.5,0.8),
                 degree = 3L,
                 intercept = FALSE,
                 Boundary.knots = c(0,1))
  d_psi <- iSpline(xx, knots = c(0.2,0.5,0.8), 
                   intercept = FALSE,
                   degree = 3L,
                   Boundary.knots = c(0,1), 
                   derivs = 1L)
  c(c(tcrossprod(psi %*% lambda, d_psi)), c(tcrossprod(d_psi %*% lambda, psi)))
}
hinjac <- function(x){
  # xx <- seq(.05,.95,by=.05)
  xx <- seq(.1,.9,by=.1)
  psi <- iSpline(xx, knots = c(0.2,0.5,0.8),
                 degree = 3L,
                 intercept = FALSE,
                 Boundary.knots = c(0,1))
  d_psi <- iSpline(xx, knots = c(0.2,0.5,0.8), 
                   intercept = FALSE,
                   degree = 3L,
                   Boundary.knots = c(0,1), 
                   derivs = 1L)
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

# starting values
# set.seed(1)
# x0 <- runif(8 * 8)
# x0 <- x0 / sum(x0)
x0 <- c(fit_ispline)

# P matrix
P <- p_fct(splines_control=isplines_control)

# verification
# check.derivatives(.x = runif(length(P)), func = neg_log_likelihood, 
#                   func_grad = grad_neg_log_likelihood, w1 = runif(10),
#                   w2 = runif(10), P = P)

fit <- slsqp(x0 = x0, fn = neg_log_likelihood, gr = grad_neg_log_likelihood, 
             lower = rep(0.0, length(x0)), upper = rep(1.0, length(x0)),
             heq = heq, heqjac = heqjac, nl.info = TRUE,
             hin = hin, hinjac = hinjac, control = list(maxeval = 5000),
             w1 = w[,1], w2 = w[,2], P = P, splines_control = isplines_control) 

lambda_hat <- matrix(fit$par, nrow = nrow(fit_ispline))

neg_log_likelihood(x0, w[,1], w[,2], P, splines_control = isplines_control)
neg_log_likelihood(fit$par, w[,1], w[,2], P, splines_control = isplines_control)

# check the fit
persp(z, theta=-25,xlab="u",ylab="w",zlab="G")
persp(b %*% fit_ispline %*% t(a),theta=-25,xlab="u",ylab="w",zlab="G",main="iSpline")
persp(b %*% lambda_hat %*% t(a),theta=-25,xlab="u",ylab="w",zlab="G",main="iSpline")

# --------------
# Estimation (bSplines)
# --------------

# inequality constraint
hin <- function(x){
  tol <- 1e-10
  knots <- c(.04,.08,.2,.5,.8,.9,.95)
  xx <- c(knots, 1.0)
  lambda <- matrix(x, nrow = 10)
  psi <- bSpline(xx, knots = knots,
                 degree = 3L,
                 intercept = FALSE,
                 Boundary.knots = c(0,1))
  d_psi <- deriv(psi)
  c(c(tcrossprod(psi %*% lambda, d_psi)), c(tcrossprod(d_psi %*% lambda, psi))) + tol
}
hinjac <- function(x){
  knots <- c(.04,.08,.2,.5,.8,.9,.95)
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

# starting values
x0 <- c(fit_bspline2)
x0[x0<0] <- 0.0
x0[x0>1] <- 1.0

# P matrix
P <- p_fct(type="b",bsplines_control)

fit2 <- slsqp(x0 = x0, fn = neg_log_likelihood, gr = grad_neg_log_likelihood, 
             lower = rep(0.0, length(x0)), upper = rep(1.0, length(x0)),
             hin = hin, hinjac = hinjac, nl.info = TRUE, control = list(maxeval = 10000),
             w1 = w[,1], w2 = w[,2], P = P, type = "b", splines_control = bsplines_control)

lambda_hat2 <- matrix(fit2$par, nrow = nrow(fit_bspline))

neg_log_likelihood(x0, w[,1], w[,2], P, type = "b", splines_control = bsplines_control)
neg_log_likelihood(c(fit_bspline), w[,1], w[,2], P, type = "b", splines_control = bsplines_control)
neg_log_likelihood(fit2$par, w[,1], w[,2], P, type = "b", splines_control = bsplines_control)

# check the fit
angle <- 240
persp(z, theta=angle,xlab="u",ylab="w",zlab="G",main="True function")
persp(a2 %*% fit_bspline %*% t(a2),theta=angle,xlab="u",ylab="w",zlab="G",main="bSpline")
persp(a2 %*% lambda_hat2 %*% t(a2),theta=angle,xlab="u",ylab="w",zlab="G",main="bSpline")

library(plotly)
fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~z, opacity = 0.5, colorscale = "red")
# fig <- fig %>% add_surface(z = ~a2 %*% lambda_hat2 %*% t(a2), opacity = 1, colorscale = "cyan", autocolorscale = FALSE)
fig <- fig %>% add_surface(z = ~a2 %*% fit_bspline %*% t(a2), opacity = 1, colorscale = "cyan", autocolorscale = FALSE)
fig


z <- a2 %*% lambda_hat2 %*% t(a2)
persp(z,xlab="u",ylab="w",zlab="G",main="bSpline", theta = 25, phi = 15, d = 1, border = "gray60")

# MAP
u_hat <- map(lambda_hat2, w[,1], w[,2], P, type = "b", bsplines_control)
plot(u,u_hat)
abline(coef=c(0,1),lwd=2)

# knots backward selection
knots <- seq(.1,.9,by=.1)
bsplines_control <- list(knots = knots, degree = 3L, Boundary.knots =  c(0,1), intercept=FALSE)
fit_bspline <- lp_fit_bSpline(u,c(v),c(w),bsplines_control)
P <- p_fct(type="b",bsplines_control)
# AIC
2 * neg_log_likelihood(c(fit_bspline), w[,1], w[,2], P, type = "b", splines_control = bsplines_control) + 2.0 * (3 + length(knots)) / n
# BIC
2 * neg_log_likelihood(c(fit_bspline), w[,1], w[,2], P, type = "b", splines_control = bsplines_control) + log(n) * (3 + length(knots)) / n

# Step 1: 4, knots <- knots[-4] (both AIC and BIC)
# Step 2: (AIC stops), 5, knots <- knots[-5]
# Step 3: 6, knots <- knots[-6]
# Step 4: stop for BIC
for(i in seq_along(knots)){
  knots2 <- knots[-i]
  bsplines_control <- list(knots = knots2, degree = 3L, Boundary.knots =  c(0,1), intercept=FALSE)
  fit_bspline <- lp_fit_bSpline(u,c(v),c(w),bsplines_control)
  P <- p_fct(type="b",bsplines_control)
  # AIC
  cat(c(
  2 * neg_log_likelihood(c(fit_bspline), w[,1], w[,2], P, type = "b", splines_control = bsplines_control) + 2.0 * (3 + length(knots)) / n,
  # BIC
  2 * neg_log_likelihood(c(fit_bspline), w[,1], w[,2], P, type = "b", splines_control = bsplines_control) + log(n) * (3 + length(knots)) / n),"\n")
}
