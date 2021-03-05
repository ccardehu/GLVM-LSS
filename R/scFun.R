
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

mfy <- function(Y,b,gr,fam){
  
 # To evaluate:
 # Y = simR$Y; b = borg; fam = fam;
 if(!is.matrix(Y)) Y <- as.matrix(Y)
 fy <- sapply(1:nrow(gr$points), function(z){mfyz(z,Y,gr$out,b,fam)})%*%gr$weights
 return(drop(fy))
}

# Log-likelihood function (to be used with optim)

loglikFun <- function(B, Y, beta, ghQ, fam){
 
 # To evaluate: B = unlist(borg); Y = Y; ghQ = gr; beta = borg; fam = fam
 B <- coefmod(bet = B, beta = beta)
 tmp <- sum(log(mfy(Y,B,ghQ,fam)))
 return(tmp)
}

# Score function (to numerically compare outputs)

ScoreFun <- function(B, Y, beta, ghQ, fam){
 # To evaluate: B = unlist(ex1$b); ghQ = ex1$gr; beta = ex1$b; loadmt = ex1$loadmt
 if(!is.matrix(Y)) Y <- as.matrix(Y)
 B <- coefmod(bet = B, beta = beta)
 efy <- mfy(Y,B,ghQ,fam)
 EC <- sapply(1:nrow(ghQ$points), function(z) mfyz(z,Y,ghQ$out,B,fam))/efy
 Sm <- NULL
 parY <- unique(unlist(lapply(1:length(fam),function(i) pFun(fam[i]))))
 pC <- vector(mode = "list", length = length(parY)); names(pC) <- parY
 for(i in parY){
  for(j in 1:ncol(Y)){ if(i %in% pFun(fam[j])) { pC[[i]] <- append(pC[[i]],j) } }
  if(i == "mu") Sm[[i]] <- c(t(sapply(pC[[i]], function(r) bmsc(r,Y,B,ghQ,fam,EC))))
  if(i == "sigma") Sm[[i]] <- c(t(sapply(pC[[i]], function(r) bssc(r,Y,B,ghQ,fam,EC))))
  if(i == "tau") Sm[[i]] <- c(t(sapply(pC[[i]], function(r) btsc(r,Y,B,ghQ,fam,EC))))
  if(i == "nu") Sm[[i]] <- c(t(sapply(pC[[i]], function(r) bnsc(r,Y,B,ghQ,fam,EC))))
 }
 return(unname(unlist(Sm)))
}

# Hessian function (to numerically compare outputs)

HessFun <- function(B, Y, beta, ghQ, fam){
 # To evaluate: B = unlist(ex1$b); ghQ = ex1$gr; beta = ex1$b; loadmt = ex1$loadmt
 if(!is.matrix(Y)) Y <- as.matrix(Y)
 B <- coefmod(bet = B, beta = beta)
 efy <- mfy(Y,B,ghQ,fam)
 EC <- sapply(1:nrow(ghQ$points), function(z) mfyz(z,Y,ghQ$out,B,fam))/efy
 Hm <- matrix()
 idxmu <- 1:ncol(ghQ$out$mu)

 for(ii in 1:ncol(Y)){
  HH <- bmhe(ii,Y,B,ghQ,fam,EC)
  if("sigma" %in% pFun(fam[ii])){
   HH <- as.matrix(Matrix::bdiag(HH, bshe(ii,Y,B,ghQ,fam,EC)))
   idxsg <- c(1:ncol(HH))[1:ncol(HH) %!in% idxmu]
   HH[idxmu,idxsg] <- mshe(ii,Y,B,ghQ,fam,EC)
   HH[idxsg,idxmu] <- t(mshe(ii,Y,B,ghQ,fam,EC))
  }
  if("tau" %in% pFun(fam[ii])){
   HH <- as.matrix(Matrix::bdiag(HH, bthe(ii,Y,B,ghQ,fam,EC)))
   idxta <- c(1:ncol(HH))[1:ncol(HH) %!in% c(idxmu,idxsg)]
   HH[idxmu,idxta] <- mthe(ii,Y,B,ghQ,fam,EC)
   HH[idxta,idxmu] <- t(mthe(ii,Y,B,ghQ,fam,EC))
   HH[idxsg,idxta] <- sthe(ii,Y,B,ghQ,fam,EC)
   HH[idxta,idxsg] <- t(sthe(ii,Y,B,ghQ,fam,EC))
  }
  if("nu" %in% pFun(fam[ii])){
   HH <- as.matrix(Matrix::bdiag(HH, bnhe(ii,Y,B,ghQ,fam,EC)))
   idxnu <- c(1:ncol(HH))[1:ncol(HH) %!in% c(idxmu,idxsg,idxta)]
   HH[idxmu,idxnu] <- mnhe(ii,Y,B,ghQ,fam,EC)
   HH[idxnu,idxmu] <- t(mnhe(ii,Y,B,ghQ,fam,EC))
   HH[idxsg,idxnu] <- snhe(ii,Y,B,ghQ,fam,EC)
   HH[idxnu,idxsg] <- t(snhe(ii,Y,B,ghQ,fam,EC))
   HH[idxta,idxnu] <- tnhe(ii,Y,B,ghQ,fam,EC)
   HH[idxnu,idxta] <- t(tnhe(ii,Y,B,ghQ,fam,EC))
  }
  Hm <- Matrix::bdiag(Hm,HH)
}
 return(as.matrix(Hm)[-1,-1])
}

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
  bmsc <- colSums(gr$out$mu*(colSums(ec*ed)*drop(gr$weights)))
  return(as.matrix(bmsc))
}

bssc <- function(i,Y,b,gr,fam,ec){
  
  # To evaluate:
  # i = 1; Y = simR$Y; b = bold; fam = fam; emfy = mfy(Y,b,gr,fam); z = 10
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  ed <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1sg)
  bssc <- colSums(gr$out$sigma*(colSums(ec*ed)*drop(gr$weights)))
  return(as.matrix(bssc))
}

bmhe <- function(i,Y,b,gr,fam,ec){
  
  # To evaluate:
  # i = 1; Y = simR$Y; b = bold; fam = fam; emfy = mfy(Y,b,gr,fam); z = 10
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  ed <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d2mu)
  t1 <- colSums(ec*ed)*drop(gr$weights)
  hm <- array(0,dim = c(ncol(gr$out$mu),ncol(gr$out$mu),length(gr$weights)))
  for(z in 1:nrow(gr$points)){ hm[,,z] <- crossprod(as.matrix(gr$out$mu[z,]))*t1[z] }
  hm <- rowSums(hm, dims = 2)
  return(as.matrix(hm))
}

bshe <- function(i,Y,b,gr,fam,ec){
  
  # To evaluate:
  # i = 1; Y = simR$Y; b = bold; fam = fam; emfy = mfy(Y,b,gr,fam); z = 45
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  ed <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d2sg)
  t1 <- colSums(ec*ed)*drop(gr$weights)
  hm <- array(0,dim = c(ncol(gr$out$sigma),ncol(gr$out$sigma),length(gr$weights)))
  for(z in 1:nrow(gr$points)){ hm[,,z] <- crossprod(as.matrix(gr$out$sigma[z,]))*t1[z] }
  hm <- rowSums(hm, dims = 2)
  return(as.matrix(hm))
}

mshe <- function(i,Y,b,gr,fam,ec){

  # To evaluate:
  # i = 1; Y = simR$Y; b = bold; fam = fam; ec = EC
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  ed <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$dcms)
  t1 <- colSums(ec*ed)*drop(gr$weights)
  hm <- array(0,dim = c(ncol(gr$out$mu),ncol(gr$out$sigma),length(gr$weights)))
  for(z in 1:nrow(gr$points)){ hm[,,z] <- crossprod(as.matrix(gr$out$mu[z,]), t(as.matrix(gr$out$sigma[z,])))*t1[z] }
  hm <- rowSums(hm, dims = 2)
  return(as.matrix(hm))  
}

 
# rbenchmark::benchmark(
#   "diag" = {diag(crossprod(ec,ed))},
#   "mat" = {colSums(ec*ed)},
#   replications = 10000,
#   columns = c("test", "replications", "elapsed",
#               "relative", "user.self", "sys.self")
# )
