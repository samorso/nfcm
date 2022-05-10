# Copula H and G functions

#' @title Copula H and G functions
#' @param u latent uniform(0,1) variable.
#' @param v uniform(0,1) errors.
#' @param copula model (see 'Details').
#' @param param \code{list} of parameters to pass to the model.
#' @details 
#' Currently we support the following \code{model}: 
#' "\code{normal}", "\code{student}", "\code{skewt}" (skew-t distribution),
#' "\code{ghyp}" (generalized hyperbolic distribution), "\code{clayton}", "\code{gumbel}".
#' The list of \code{param} depends on the model:
#' \itemize{
#'    \item \code{normal}: either correlation \code{corr} (between -1 and 1) or 
#'    variances of \eqn{Z} \code{var_z}\eqn{>0} and \eqn{\epsilon} \code{var_e}\eqn{>0};
#'    \item \code{student}: degrees of freedom \code{nu}\eqn{>0};
#'    \item \code{skewt}: degrees of freedom \code{nu}\eqn{>0}, 
#'    the variance of \eqn{\epsilon} \code{var_e}\eqn{>0} and 
#'    the parameter of asymetry \code{lambda} (betweem -1 and 1);
#'    \item \code{ghyp}: \code{lambda}, \code{chi} (positive), \code{psi} (positive),
#'    the variance of \eqn{\epsilon} \code{var_e}\eqn{>0} and \code{gamma};
#'    \item \code{clayton}: \code{alpha} (positive), the parameter of dependence;
#'    \item \code{gumbel}: \code{alpha}\eqn{>1}, the parameter of dependence.
#' }
#' @return A vector or matrix of the same dimension as \code{v} for \code{H}
#' or as \code{w} for \code{G}.
#' @importFrom stats pnorm qnorm pt qchisq qgamma
#' @importFrom sgt psgt
#' @importFrom ghyp ghyp pghyp qgig
#' @importFrom stabledist qstable
#' @export
H <- function(u, v, copula = "normal", param = list()){
  # general verifications
  if(!is.numeric(u)) stop("'u' must be numeric")
  if(!is.numeric(v)) stop("'v' must be numeric")
  if(any(u<0 | u>1)) stop("'u' must be between 0 and 1")
  if(any(v<0 | v>1)) stop("'v' must be between 0 and 1")
  if(length(v) %% length(u) != 0) stop("incorrect dimensions for 'u' and 'v'")
  if(!copula%in%c("normal","student","skewt","ghyp","clayton","gumbel")) 
    stop(sprintf("%s is not supported"), copula)
  if(!is.list(param)) stop("'param' must be provided as a list")
  if(!all(sapply(param,is.numeric))) stop("elements of 'param' must be numeric")
  
  # specific implementation
  # normal copula
  if(copula == "normal"){
    if(any(length(param)>2 || length(param)==0)) stop("'param' must be of length 1 (correlation) or 2 (variances)")
    
    # standard normal: param is correlation
    if(length(param)==1){
      corr <- param[[1]]
      if(corr<=-1 || corr>=1) stop("correlation must be between -1 and 1")
      return(pnorm(sqrt(corr) * qnorm(u) + sqrt(1 - corr) * qnorm(v)))
    }
    
    # centered normal: param is variances
    if(length(param)==2){
      if(any(param<=0)) stop("'param' must be positive")
      if(!all(names(param)%in%c("var_z","var_e"))) stop("'param' name error")
      return(pnorm((sqrt(param$var_z) * qnorm(u) + sqrt(param$var_e) * qnorm(v))) / sqrt(sum(unlist(param))))
    }
  }
  
  # student copula
  if(copula == "student"){
    if(length(param)!=1) stop("'param' must be of length 1: nu or degrees of freedom")
    if(param<=0) stop("'param' must be positive")
    if(names(param)!="nu") stop("'param' name error")
    return(pt(qnorm(v) * sqrt(param$nu) / sqrt(qchisq(u, df = param$nu)), df = param$nu))
  }
  
  # skew-t copula
  if(copula == "skewt"){
    if(length(param)!=3) stop("'param' must be of length 3: lambda, nu and var_e")
    if(!all(names(param)%in%c("lambda","nu","var_e"))) stop("'param' name error")
    if(param$lambda< -1 || param$lambda>1) stop("lambda must be between -1 and 1")
    if(param$nu<=0) stop("nu must be positive")
    if(param$var_e<=0) stop("var_e must be positive")
    return(psgt(param$lambda * qchisq(u, df = param$nu) / param$nu + sqrt(param$var_e) * qnorm(v) * sqrt(param$nu) / sqrt(qchisq(u, df = param$nu)),
                lambda = param$lambda, q = param$nu, sigma = param$var_e, var.adj = FALSE))
  }
  
  # generalized hyperbolic copula
  if(copula == "ghyp"){
    if(length(param)!=5) stop("'param' must be of length 5: lambda, chi, psi, var_e and gamma")
    if(!all(names(param)%in%c("lambda","chi","psi","var_e","gamma"))) stop("'param' name error")
    if(param$chi<=0) stop("chi must be positive")
    if(param$psi<=0) stop("psi must be positive")
    if(param$var_e<=0) stop("var_e must be positive")
    ghyp_obj <- ghyp(lambda = param$lambda, chi = param$chi, psi = param$psi, sigma = param$var_e, gamma = param$gamma)
    xx <- qgig(u, lambda = param$lambda, chi = param$chi, psi = param$psi)
    return(pghyp(param$lambda * xx + sqrt(param$var_e) * qnorm(v) * sqrt(xx), object = ghyp_obj))
  }
  
  # clayton copula
  if(copula == "clayton"){
    if(length(param)!=1) stop("'param' must be of length 1: alpha")
    if(names(param)!="alpha") stop("'param' name error")
    if(param$alpha<=0) stop("'alpha' must be positive")
    return(1.0 / (1.0 - log(v) / qgamma(u, shape = param$alpha))^param$alpha)
  }
  
  # gumbel copula
  if(copula == "gumbel"){
    if(length(param)!=1) stop("'param' must be of length 1: alpha")
    if(names(param)!="alpha") stop("'param' name error")
    if(param$alpha<=1) stop("'alpha' must be greater than 1")
    gamma <- cos(pi / 2.0 / param$alpha)^param$alpha
    return(exp(- (- log(v) / qstable(u, alpha = param$alpha, beta = 1, gamma = gamma))^(1.0 / param$alpha)))
  }
}

#' @rdname H
#' @param w uniform(0,1) observable variables
#' @importFrom stats qt
#' @importFrom sgt qsgt
#' @importFrom ghyp qghyp
#' @export
G <- function(u, w, copula = "normal", param = list()){
  # general verifications
  if(!is.numeric(u)) stop("'u' must be numeric")
  if(!is.numeric(w)) stop("'w' must be numeric")
  if(any(u<0 | u>1)) stop("'u' must be between 0 and 1")
  if(any(w<0 | w>1)) stop("'w' must be between 0 and 1")
  if(length(w) %% length(u) != 0) stop("incorrect dimensions for 'u' and 'v'")
  if(!copula%in%c("normal","student","skewt","ghyp","clayton","gumbel")) 
    stop(sprintf("%s is not supported"), copula)
  if(!is.list(param)) stop("'param' must be provided as a list")
  if(!all(sapply(param,is.numeric))) stop("elements of 'param' must be numeric")
  
  # specific implementation
  # normal copula
  if(copula == "normal"){
    if(any(length(param)>2 || length(param)==0)) stop("'param' must be of length 1 (correlation) or 2 (variances)")
    
    # standard normal: param is correlation
    if(length(param)==1){
      corr <- param[[1]]
      if(corr<=-1 || corr>=1) stop("correlation must be between -1 and 1")
      return(pnorm(qnorm(w) / sqrt(1.0 - corr) - sqrt(corr / (1.0 - corr)) * qnorm(u)))
    }
    
    # centered normal: param is variances
    if(length(param)==2){
      if(any(param<=0)) stop("'param' must be positive")
      if(!all(names(param)%in%c("var_z","var_e"))) stop("'param' name error")
      return(pnorm(qnorm(w) * sqrt(sum(unlist(param))) / sqrt(param$var_e) - sqrt(param$var_z) * qnorm(u) / sqrt(param$var_e)))
    }
  }
  
  # student copula
  if(copula == "student"){
    if(length(param)!=1) stop("'param' must be of length 1: nu or degrees of freedom")
    if(any(param<=0)) stop("'param' must be positive")
    if(names(param)!="nu") stop("'param' name error")
    return(pnorm(qt(w, df = param$nu) * sqrt(qchisq(u, df = param$nu)) / sqrt(param$nu)))
  }
  
  # skew-t copula
  if(copula == "skewt"){
    if(length(param)!=3) stop("'param' must be of length 3: lambda, nu and var_e")
    if(!all(names(param)%in%c("lambda","nu","var_e"))) stop("'param' name error")
    if(param$lambda< -1 || param$lambda>1) stop("lambda must be between -1 and 1")
    if(param$nu<=0) stop("nu must be positive")
    if(param$var_e<=0) stop("var_e must be positive")
    return(pnorm((qsgt(w, lambda = param$lambda, q = param$nu, sigma = param$var_e, var.adj = FALSE) - 
                    param$lambda * qchisq(u, df = param$nu) / param$nu) / 
                   sqrt(param$var_e) / sqrt(param$nu) * sqrt(qchisq(u, df = param$nu))))
  }
  
  # generalized hyperbolic copula
  if(copula == "ghyp"){
    if(length(param)!=5) stop("'param' must be of length 5: lambda, chi, psi, var_e and gamma")
    if(!all(names(param)%in%c("lambda","chi","psi","var_e","gamma"))) stop("'param' name error")
    if(param$chi<=0) stop("chi must be positive")
    if(param$psi<=0) stop("psi must be positive")
    if(param$var_e<=0) stop("var_e must be positive")
    ghyp_obj <- ghyp(lambda = param$lambda, chi = param$chi, psi = param$psi, sigma = param$var_e, gamma = param$gamma)
    xx <- qgig(u, lambda = param$lambda, chi = param$chi, psi = param$psi)
    return(pnorm((qghyp(w, object = ghyp_obj) - param$lambda * xx) / sqrt(xx) / sqrt(param$var_e)))
  }
  
  # clayton copula
  if(copula == "clayton"){
    if(length(param)!=1) stop("'param' must be of length 1: alpha")
    if(names(param)!="alpha") stop("'param' name error")
    if(param$alpha<=0) stop("'alpha' must be positive")
    return(exp((1.0 - w^(-1.0 / param$alpha)) * qgamma(u, shape = param$alpha)))
  }
  
  # gumbel copula
  if(copula == "gumbel"){
    if(length(param)!=1) stop("'param' must be of length 1: alpha")
    if(names(param)!="alpha") stop("'param' name error")
    if(param$alpha<=1) stop("'alpha' must be greater than 1")
    gamma <- cos(pi / 2.0 / param$alpha)^param$alpha
    return(exp(-qstable(u, alpha = param$alpha, beta = 1, gamma = gamma) * (-log(w))^param$alpha))
  }
}