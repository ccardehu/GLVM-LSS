
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Functions for conditional, marginal, and complete log-likelihoods #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Vector function f(Y_m|z_m) for each Y_m

mfyz <- function(z,Y,Z,b,fam){ # should be evaluated at Y, fam, for i = 1:p; \forall Z = gr$out[z,] (rows)
  
  # To evaluate:
  # Y = simR$Y; Z = gr$out; b = bold; fam = fam; z = 10
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  fyz <- sapply(1:ncol(Y), FUN = function(i) dFun(i,z, fam = fam[i], Y = Y[,i], Z = Z, b = b))
  fyz <- exp(rowSums(fyz))
  return(drop(fyz))
}

# Vector function for f(Y_m) = integral f(Y_m|z_m)

mfy <- function(Y,b,gr,fam){ # should be evaluated at Y = Y; Z = gr$out; b = borg
  
  # To evaluate:
  # Y = simR$Y; b = borg; fam = fam;
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  fy <- sapply(1:nrow(gr$points), function(z){mfyz(z,Y,gr$out,b,fam)})%*%gr$weights
  return(drop(fy))
}

# Log-likelihood function (to be used with optim)

loglik <- function(B, ghQ, beta, loadmt){
  
  # To evaluate: B = unlist(borg); ghQ = gr; beta = borg; loadmt = loadmt
  
  B <- coefmod(bet = B, beta = beta, gr = ghQ, loadmt = loadmt)
  tmp <- sum(log(mfy(Y,B,ghQ,fam)))
  return(tmp)
}

# Log-likelihood function (to be used with trust)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Functions for computing Score and Hessian #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Step (3): Compute score vectors (with quadrature) for each \alpha_i (function of j and theta)

bmsc <- function(i,Y,b,gr,fam,ec){
  
  # To evaluate:
  # i = 1; Y = simR$Y; b = bold; gr = gr; fam = fam; ec = EC; z = 10
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  ed <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1mu)
  bmsc <- colSums(gr$out$mu*(diag(crossprod(ec,ed))*drop(gr$weights)))
  return(as.matrix(bmsc))
}

bssc <- function(i,Y,b,gr,fam,ec){
  
  # To evaluate:
  # i = 1; Y = simR$Y; b = bold; fam = fam; emfy = mfy(Y,b,gr,fam); z = 10
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  ed <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1sg)
  bssc <- colSums(gr$out$sigma*(diag(crossprod(ec,ed))*drop(gr$weights)))
  return(as.matrix(bssc))
}

bmhe <- function(i,Y,b,gr,fam,ec){
  
  # To evaluate:
  # i = 1; Y = simR$Y; b = bold; fam = fam; emfy = mfy(Y,b,gr,fam); z = 10
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  ed <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d2mu)
  t1 <- diag(crossprod(ec,ed))*drop(gr$weights)
  eh <- lapply(1:nrow(gr$points), function(z) crossprod(as.matrix(gr$out$mu[z,]))*t1[z])
  eh <- Reduce("+",eh)
  return(as.matrix(eh))
}

bshe <- function(i,Y,b,gr,fam,ec){
  
  # To evaluate:
  # i = 1; Y = simR$Y; b = bold; fam = fam; emfy = mfy(Y,b,gr,fam); z = 45
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  ed <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d2sg)
  t1 <- diag(crossprod(ec,ed))*drop(gr$weights)
  eh <- lapply(1:nrow(gr$points), function(z) crossprod(as.matrix(gr$out$sigma[z,]))*t1[z])
  eh <- Reduce("+",eh)
  return(as.matrix(eh))
}

