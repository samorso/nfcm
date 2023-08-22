library(splines2)

x <- seq.int(0, 1, 0.01)
knots <- c(0.3, 0.5, 0.6)

## quadratic B-splines
bsMat <- bSpline(x, knots = knots, degree = 3, intercept = TRUE)
Jn <- length(knots) + 3 + 1

op <- par(mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
matplot(x, bsMat, type = "l", ylab = "Cubic B-splines")
abline(v = knots, lty = 2, col = "gray")
y <- bsMat %*% seq(0,1,length=Jn)
y <- bsMat %*% c(0,.3,.5,.6,.7,.8,1)
y <- bsMat %*% c(0,.3,.1,.38,.5,.8,1)
lines(x,y,lwd=2)

d1Mat <- deriv(bSpline(c(0,knots,1), knots = knots, degree = 3, intercept = TRUE))
d1Mat %*% c(0,.3,.1,.4,.5,.8,1)

deriv(bsMat) %*% c(0,0,.3,.1,.5,.35,.8,1)
