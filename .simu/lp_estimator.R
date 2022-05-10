library("splines2")
library("Rglpk")
library("slam")

# Function for estimating spline coefficients based on Linear Programming
lp_fit_iSpline <- function(u,v,w,control){
  # u - common latent variable from unif(0,1)
  # v - unif(0,1) errors
  # w - observations from unif(0,1)
  # control - list of control for iSpline()
  # Note: for the moment G(u,w) is suppose to be monotonic increasing in
  # both arguments. If dim(u) < dim(w), it means that dim(u) = n and dim(w) = n * d = N.
  
  N <- length(w)
  n <- length(u)
  d <- N / n
  if(d>1) u <- rep(u,d)
  # u <- 1.0 - u # TODO: check direction
  
  # basis functions
  phi <- iSpline(u, knots = control$knots, degree = control$degree,
                 Boundary.knots = control$Boundary.knots, 
                 intercept = control$intercept)
  psi <- iSpline(w, knots = control$knots, degree = control$degree,
                 Boundary.knots = control$Boundary.knots,
                 intercept = control$intercept)
  
  # linear problem is of the form:
  # min C^T\xi under constraint A\xi <= b
  
  # Constraints:
  p <- ncol(phi)
  p2 <- p * p
  C <- c(rep(0,p2),rep(1,N))
  pp <- matrix(rep(phi,p),nr=N) * matrix(rep(t(psi),each=p),nr=N,byrow=T)
  
  # "regular" constrain
  A1 <- cbind(rbind(pp,-pp)[c(rbind(seq(N),seq(N)+N)),],
             diag(N) %x% rbind(-1,-1))
  b1 <- cbind(c(rbind(v,-v)))
  
  # equality constrain
  A2 <- c(rep(1,p2),rep(0,N))
  b2 <- 1.0
  
  # overall constrains
  # all coefficients in (0,1)
  bounds <- list(upper = list(ind=seq_len(p2),val=rep(1,p2)))
  A <- rbind(A1,A2)
  b <- rbind(b1,b2)
  f_dir <- c(rep("<=",nrow(A1)), "==")
  A <- as.simple_triplet_matrix(A)
  C <- as.simple_triplet_matrix(C)
  fit_lp <- Rglpk_solve_LP(obj = C, mat = A, dir = f_dir, rhs = b, bounds = bounds)
  matrix(fit_lp$solution[seq_len(p2)],nrow = p)
}

# Function for estimating spline coefficients based on Linear Programming
# (no constraint such as w_11==0 and w_pp==1)
lp_fit_bSpline <- function(u,v,w,control){
  # u - common latent variable from unif(0,1)
  # v - unif(0,1) errors
  # w - observations from unif(0,1)
  # control - list of control for bSpline()
  # Note: for the moment G(u,v) is suppose to be monotonic increasing in
  # both arguments. If dim(u) < dim(v), it means that dim(u) = n and dim(v) = n * d = N.
  
  N <- length(w)
  n <- length(u)
  d <- N / n
  if(d>1) u <- rep(u,d)
  # u <- 1.0 - u # TODO: check direction
  k <- length(control$knots) + 1L
  
  # basis functions
  # at true observations (we will need to change that since we do not observe u0,v0)
  phi <- bSpline(u, knots = control$knots, degree = control$degree,
                 Boundary.knots = control$Boundary.knots,
                 intercept = control$intercept)
  psi <- bSpline(w, knots = control$knots, degree = control$degree,
                 Boundary.knots = control$Boundary.knots,
                 intercept = control$intercept)
  
  # partial derivatives (at the knots)
  psi0 <- bSpline(c(control$knots,1), knots = control$knots, degree = control$degree,
                  Boundary.knots = control$Boundary.knots,
                  intercept = control$intercept)
  d_psi0 <- deriv(psi0)
  # d_phi <- deriv(phi)
  # d_psi <- deriv(psi)
  
  # linear problem is of the form:
  # min C^T\xi under constraint A\xi <= b
  
  # Constraints:
  p <- ncol(phi)
  p2 <- p * p
  C <- c(rep(0,p2),rep(1,N)) 
  pp <- matrix(rep(phi,p),nr=N) * matrix(rep(t(psi),each=p),nr=N,byrow=T)

    # "regular" constrain
  A1 <- cbind(rbind(pp,-pp)[c(rbind(seq(N),seq(N)+N)),],
              diag(N) %x% rbind(-1,-1))
  b1 <- cbind(c(rbind(v,-v)))
  # additional monotonic constrain on partial derivatives
  partial_d_phi <- partial_d_psi <- matrix(nrow = k * k, ncol = p2)
  ind <- 0L
  for(i in seq_len(k)){
    for(j in seq_len(k)){
      ind <- ind + 1L
      partial_d_phi[ind,] <- c(tcrossprod(d_psi0[i,],psi0[j,]))
      partial_d_psi[ind,] <- c(tcrossprod(psi0[i,],d_psi0[j,]))
    }
  }
  A2 <- cbind(
    rbind(-partial_d_phi,-partial_d_psi),
    matrix(0, nrow = 2 * k * k, ncol = N)
  )
  b2 <- cbind(rep(0,2 * k * k))
  
  # overall constrains
  # all coefficients in (0,1)
  bounds <- list(upper = list(ind=seq_len(p2),val=rep(1,p2)))
  A <- rbind(A1,A2)
  b <- rbind(b1,b2)
  f_dir <- rep("<=",nrow(A))
  A <- as.simple_triplet_matrix(A)
  C <- as.simple_triplet_matrix(C)
  fit_lp <- Rglpk_solve_LP(obj = C, mat = A, dir = f_dir, rhs = b, bounds = bounds)
  matrix(fit_lp$solution[seq_len(p2)],nrow=p)
}