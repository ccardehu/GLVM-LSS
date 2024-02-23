#' Fit a Generalized latent variable model for location, scale, and shape parameters
#'
#' @param data Matrix \code{(n x p)} of observed variables. Missing data should be coded as \code{NA}.
#' @param family String vector of size \code{p} with family for each observed variable.
#' @param mu.eq Formula for location parameter e.g.,  \code{mu.eq = ~ Z1}.
#' @param sg.eq (optional) Formula for scale parameter.
#' @param nu.eq (optional) Formula for shape parameter.
#' @param ta.eq (optional) Formula for shape parameter.
#' @param control List of control parameters (see 'Details').
#' @param ... Further arguments to be passed to \code{control}.
#'
#' @return A list with components:
#' \itemize{
#'  \item \code{b}: An array with the estimated intercept and slopes.
#'  \item \code{R}: Estimated factor correlation matrix (if \code{control$corLV = T}).
#'  \item \code{ll}: Marginal log-likelihood evaluated at \code{b} and \code{Rz}.
#'  \item \code{seb}: Estimated standard errors for \code{b}.
#'  \item \code{seR}: Estimated standard errors for \code{R} (if \code{control$corLV = T}).
#' }
#' @details Test
#' @author Camilo Cárdenas-Hurtado (\email{c.a.cardenas-hurtado@@lse.ac.uk}).
#' @export
glvmlss <- function(data, family,
                    mu.eq = ~ Z1, sg.eq = NULL,
                    nu.eq = NULL, ta.eq = NULL,
                    control = list(), ...){

  Y <- as.matrix(data)
  control <- pr_control(control, ...)
  form <- pr_form(mu.eq, sg.eq, nu.eq, ta.eq)
  q <- pr_q(form)
  bcube <- pr_b(Y,family,form,control,q)
  dnam <- dimnames(bcube)
  if(is.null(control$start.R) & control$corLV == F) control$start.R <- diag(q)
  if(is.null(control$start.R) & control$corLV == T){
    control$start.R <- matrix(0.05, q, q)
    diag(control$start.R) <- 1
  }
  out <- glvmlss_rcpp(Yr = Y,famr = family,q = q,b = bcube[], control = control)
  dimnames(out$b) <- dimnames(out$seb) <- dnam
  return(out)
}

#' Simulate data for a Generalized latent variable model for location, scale, and shape parameters
#'
#' @param n Number of simulated entries.
#' @param family String vector of size \code{p} with family for each observed variable.
#' @param mu.eq Formula for location parameter e.g.,  \code{mu.eq = ~ Z1}.
#' @param sg.eq (optional) Formula for scale parameter.
#' @param nu.eq (optional) Formula for shape parameter.
#' @param ta.eq (optional) Formula for shape parameter.
#' @param control List of control parameters (see 'Details').
#' @param ... Further arguments to be passed to \code{control}.
#'
#' @return A list with components:
#' \itemize{
#'  \item \code{Y}: A matrix \code{(n x p)} of observed variables, each column from a distribution specified in \code{family}.
#'  \item \code{Z}: Simulated latent variables (includes intercept).
#'  \item \code{b}: Original coefficients.
#'  \item \code{R}: Original factor correlation matrix.
#' }
#' @details Test
#' @author Camilo Cárdenas-Hurtado (\email{c.a.cardenas-hurtado@@lse.ac.uk}).
#' @export
glvmlss_sim <- function(n, family,
                        mu.eq = ~ Z1, sg.eq = NULL,
                        nu.eq = NULL, ta.eq = NULL,
                        control = list(), ...){

  control <- pr_control(control, ...)
  form <- pr_form(mu.eq, sg.eq, nu.eq, ta.eq)
  q <- pr_q(form)
  if(is.null(control$start.R)) control$start.R <- diag(q)
  bcube <- pr_bsim(control,q,length(family),form)
  dnam <- dimnames(bcube)
  out <- glvmlss_rcpp_sim(n = n,famr = family,q = q,b = bcube[], control = control)
  dimnames(out$b) <- dnam
  colnames(out$Y) <- rownames(bcube)
  colnames(out$Z) <- colnames(bcube)
  return(out)
}
