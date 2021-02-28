
# https://biostatmatt.com/archives/2754

library('statmod')
library('fastGHQuad')

## Calculate nodes and weights for Gauss-Hermite quadrature
## with weights corresponding to standard normal distribution
gH <- function(n) {
  ## Compute Gauss-Hermite quadrature nodes and weights
  # pts <- gauss.quad(n, kind="hermite")
  pts <- list("nodes" = gaussHermiteData(n)$x, "weights" = gaussHermiteData(n)$w)
  
  ## By default, the weight function for Hermite quadrature is 
  ## w(x) = exp(-x^2). Need to adjust weights so that 
  ## w(x) = (2*pi)^(1/2)*exp(-1/2*x^2)
  pts$weights <- pts$weights * (2*pi)^(-1/2)*exp(pts$nodes^2/2)
  #pts$nodes <- pts$nodes * sqrt(2)
  #pts$weights <- pts$weights * 1/(sqrt(pi))
  
  ## rename 'nodes' to conform with code below
  pts$points <- pts$nodes
  pts$nodes <- NULL
  return(pts)
}


# compute multivariate Gaussian quadrature points
# n     - number of points each dimension before pruning
# mu    - mean vector
# sigma - covariance matrix
# prune - NULL - no pruning; [0-1] - fraction to prune
# # # # # #

mvgH <- function(n, mu, sigma, prune = NULL, formula = "~ Z1 + Z2") {
  
  # evaluated @ q. = 2; n = 10; mu = rep(0,q.); sigma = diag(q.); formula = form
  
  nl <- unique(unlist(lapply(1:length(formula), function(i) all.vars(as.formula(formula[[i]])))))
  nl <- nl[grep("Z", nl, fixed = T)]
  if(missing(mu) && missing(sigma)){
   mu <- rep(0,length(nl))
   sigma <- diag(length(nl))
  }
  
  if(!all(dim(sigma) == length(mu))) stop("mu and sigma have nonconformable dimensions")
  if(length(mu) > 0) dm  <- length(mu) else dm <- 1
  if(missing(n)) n <- round(100^(1/dm))
  gh  <- gH(n)
  # idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh$points[idx],nrow(idx),dm)
  wts <- as.matrix(apply(matrix(gh$weights[idx],nrow(idx),dm), 1, prod))
  
  # Pruning
  # ~~~~~~~
   
  if(!is.null(prune)) {
   qwt <- quantile(wts, probs=prune)
   pts <- pts[wts > qwt,]
   wts <- wts[wts > qwt]
  }
  
  # Rotating, scaling, and translate points
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  eig <- eigen(sigma) 
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  
  colnames(pts) <- nl

  out <- NULL
  for(i in names(formula)){
    out[[i]] <- as.data.frame(model.matrix(as.formula(formula[[i]]), as.data.frame(pts)))
    colnames(out[[i]])[1] <- "Int"
  }
  
  return(list(points = pts, weights = wts, out = out, n = n))
}


# gr <- mvgH(n = 50, mu = rep(0,q.), sigma = diag(q.), formula = form)