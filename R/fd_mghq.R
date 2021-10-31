
##############################################################
## SPLVM: Multivariate Gaussian Hermite Quadrature function ##
##############################################################

mvghQ <- function(n, mu, sigma, formula = "~ Z1 + Z2") {
  
# Goal: Computes mvgHQ (points & weights)
# Input : n (number of qp), mu (mean of LV), sigma (LV covariance matrix),
#         prune (reducing points / NOT TESTED YET),
#         formula (list of formulas for measurement eqs. for each parameter).
# Output: List with points, weights and design matrices with ghQs.
# Testing: n = 15; mu = c(0,0); sigma = diag(2); prune = NULL;
#          formula = form.

nl <- unique(unlist(lapply(1:length(formula), function(i) all.vars(as.formula(formula[[i]])))))
nl <- nl[grep("Z", nl, fixed = T)]
if(missing(mu) && missing(sigma)){
 mu <- rep(0,length(nl))
 sigma <- diag(length(nl))
}
  
if(!all(dim(sigma) == length(mu))) stop("mu and sigma have nonconformable dimensions (in ghQ)")
if(length(mu) > 0) dm  <- length(mu) else dm <- 1
if(missing(n)) n <- round(100^(1/dm))
gh  <- list("points" = fastGHQuad::gaussHermiteData(n)$x, "weights" = fastGHQuad::gaussHermiteData(n)$w)
gh$weights <-  gh$weights * (2*pi)^(-1/2)*exp(gh$points^2/2)
idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
pts <- matrix(gh$points[idx],nrow(idx),dm)
wts <- as.matrix(apply(matrix(gh$weights[idx],nrow(idx),dm), 1, prod))
# Rotating if mu & sigma != NULL
eig <- eigen(sigma) 
rot <- eig$vectors %*% diag(sqrt(eig$values))
pts <- t(rot %*% t(pts) + mu)
colnames(pts) <- nl
out <- NULL
for(i in names(formula)){
 out[[i]] <- as.data.frame(model.matrix(as.formula(formula[[i]]), as.data.frame(pts)))
}
return(list(points = pts, weights = wts, out = out, n = n))
}