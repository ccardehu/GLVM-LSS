
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
 return(c(fyz))
}

# Vector function for f(Y_m) = integral f(Y_m|z_m)

mfy <- function(Y,b,gr,fam){
  
 # To evaluate:
 # Y = simR$Y; b = borg; fam = fam;
 if(!is.matrix(Y)) Y <- as.matrix(Y)
 fy <- sapply(1:nrow(gr$points), function(z){mfyz(z,Y,gr$out,b,fam)})%*%gr$weights
 return(c(fy))
}

# Log-likelihood function (to be used with optim)

loglikFun <- function(B, Y, beta, ghQ, fam){
 
 # To evaluate: B = unlist(borg); Y = Y; ghQ = gr; beta = borg; fam = fam
 B <- coefmod2(bet = B, beta = beta)
 tmp <- sum(log(mfy(Y,B,ghQ,fam)))
 return(tmp)
}

llkFun <- function(B, Y, beta, ghQ, fam){
  # To evaluate: B = decoefmod(borg); Y = Y; ghQ = ex1$gr; beta = borg; fam = fam
  B <- coefmod2(bet = B, beta = beta)
  tmp <- sum(log(mfy(Y,B,ghQ,fam)))
  tmp2 <- ScHesFun(decoefmod(B), Y, beta, ghQ, fam)
  return(list(value = tmp, gradient = tmp2$score, hessian = tmp2$hessian))
}

zsc <- function(mod){
 EC <- sapply(1:nrow(mod$gr$points), function(z) mfyz(z,mod$Y,mod$gr$out,mod$b,mod$fam)*mod$gr$weights[z])/mfy(mod$Y,mod$b,mod$gr,mod$fam)
 zscore <- EC%*%mod$gr$points[,,drop=F]
 zscore <- as.data.frame(zscore); names(zscore) <- colnames(mod$gr$points)
 return(zscore)
}

# Score function (to numerically compare outputs)

ScHesFun <- function(B, Y, beta, ghQ, fam){
 # To evaluate: B = decoefmod(ex1$b); ghQ = ex1$gr; beta = ex1$b; loadmt = ex1$loadmt
 if(!is.matrix(Y)) Y <- as.matrix(Y)
 di <- sum(B != 0)
 B <- coefmod2(bet = B, beta = beta)
 efy <- mfy(Y,B,ghQ,fam)
 EC <- sapply(1:nrow(ghQ$points), function(z) mfyz(z,Y,ghQ$out,B,fam))/efy
 Sm <- NULL
 idx <- idp <- rep(0,ncol(Y))
 Hm <- matrix(0,di,di)
 for(i in 1:ncol(Y)){
  tmpx <- bsc(i,Y,B,ghQ,fam,EC); idx[i] <- length(tmpx)
  Sm <- c(Sm,tmpx); idp[i] <- length(Sm)
  Hr <- NULL
  for(j in i:ncol(Y)){
   Hr <- rbind(Hr,bhe(i,j,Y,B,ghQ,fam,EC))
  }
  Hm[(idp[i]-idx[i]+1):nrow(Hm),(idp[i]-idx[i]+1):idp[i]] <- Hr
  Hm[(idp[i]-idx[i]+1):idp[i],(idp[i]-idx[i]+1):ncol(Hm)] <- t(Hr)
 }
 return(list("score" = Sm, "hessian" = Hm))
}

iScHesFun <- function(i, B, Y, beta, ghQ, fam){
# To evaluate: B = decoefmod(ex1$b); ghQ = ex1$gr; beta = ex1$b; loadmt = ex1$loadmt
if(!is.matrix(Y)) Y <- as.matrix(Y)
# di <- sum(B != 0)
B <- coefmod2(bet = B, beta = beta)
efy <- mfy(Y,B,ghQ,fam)
EC <- sapply(1:nrow(ghQ$points), function(z) mfyz(z,Y,ghQ$out,B,fam))/efy
sc <- bsc(i,Y,B,ghQ,fam,EC)
he <- bhe(i,i,Y,B,ghQ,fam,EC)
idx <- ((i)*length(sc)):(length(sc))
return(list("score" = sc, "hessian" = he))
}

ScoreFun <- function(B, Y, beta, ghQ, fam){
 # To evaluate: B = decoefmod(ex1$b); ghQ = ex1$gr; beta = ex1$b; loadmt = ex1$loadmt
 if(!is.matrix(Y)) Y <- as.matrix(Y)
 B <- coefmod2(bet = B, beta = beta)
 efy <- mfy(Y,B,ghQ,fam)
 EC <- sapply(1:nrow(ghQ$points), function(z) mfyz(z,Y,ghQ$out,B,fam))/efy
 Sm <- NULL
 for(i in 1:ncol(Y)){
 Sm <- c(Sm,bsc(i,Y,B,ghQ,fam,EC))
 }
 return(unname(unlist(Sm)))
}

# Hessian function (to numerically compare outputs)

HessFun <- function(B, Y, beta, ghQ, fam){
 # To evaluate: B = decoefmod(ex1$b); ghQ = ex1$gr; beta = ex1$b; loadmt = ex1$loadmt
 if(!is.matrix(Y)) Y <- as.matrix(Y)
 B <- coefmod2(bet = B, beta = beta)
 efy <- mfy(Y,B,ghQ,fam)
 EC <- sapply(1:nrow(ghQ$points), function(z) mfyz(z,Y,ghQ$out,B,fam))/efy
 Hm <- NULL
 for(i in 1:ncol(Y)){
  Hr <- NULL
  for(j in 1:ncol(Y)){
   Hr <- rbind(Hr,bhe(i,j,Y,B,ghQ,fam,EC))
  }
  Hm <- cbind(Hm,Hr)
 }
 return(as.matrix(Hm))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Functions for computing Score and Hessian #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Step (3): Compute score vectors (with quadrature) for each \alpha_i (function of j and theta)

bsc <- function(i,Y,b,gr,fam,ec){
# To evaluate:
# i = 1; Y = simR$Y; b = bold; gr = gr; fam = fam; ec = EC; z = 10
# ec <- sapply(1:nrow(gr$points), function(z) mfyz(z,Y,gr$out,b,fam))/mfy(Y,b,gr,fam)
 if(!is.matrix(Y)) Y <- as.matrix(Y)
 wmu <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1mu)*ec)*c(gr$weights)
 if("sigma" %in%  pFun(fam[i])) wsg <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1sg)*ec)*c(gr$weights)
 if("tau" %in%  pFun(fam[i])) wta <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1ta)*ec)*c(gr$weights)
 if("nu" %in%  pFun(fam[i])) wnu <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1nu)*ec)*c(gr$weights)
 ps <- gr$out$mu[,which(b$mu[i,] != 0)] * wmu # change here
 if("sigma" %in% pFun(fam[i])) ps <- cbind(ps,gr$out$sigma[,which(b$sigma[i,] != 0)] * wsg) # change here
 if("tau" %in% pFun(fam[i])) ps <- cbind(ps,gr$out$tau[,which(b$tau[i,] != 0)] * wta) # change here
 if("nu" %in% pFun(fam[i])) ps <- cbind(ps,gr$out$nu[,which(b$nu[i,] != 0)] * wnu) # change here
 return(as.matrix(colSums(unname(ps))))
}

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
  for(z in 1:nrow(gr$points)){ hm[,,z] <- t(as.matrix(gr$out$mu[z,,drop=F])) %*% as.matrix(gr$out$sigma[z,,drop=F])*t1[z] }
  hm <- rowSums(hm, dims = 2)
  return(as.matrix(hm))  
}

bhe1 <- function(i,j,Y,b,gr,fam,ec){
# To evaluate:
# i = 1; Y = simR$Y; b = ex1$b; gr = ex1$gr; fam = fam; z = 1
# ec <- sapply(1:nrow(gr$points), function(z) mfyz(z,Y,gr$out,b,fam))/mfy(Y,b,gr,fam)
if(!is.matrix(Y)) Y <- as.matrix(Y)
if(i == j){
 wmu <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d2mu)*ec)*c(gr$weights)
 if("sigma" %in%  pFun(fam[i])){
  wsg <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d2sg)*ec)*c(gr$weights)
  wms <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$dcms)*ec)*c(gr$weights)
 } 
 if("tau" %in%  pFun(fam[i])){
  wta <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d2ta)*ec)*c(gr$weights)
  wmt <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$dcmt)*ec)*c(gr$weights)
  wst <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$dcst)*ec)*c(gr$weights)
 }
 if("nu" %in%  pFun(fam[i])){
  wnu <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d2nu)*ec)*c(gr$weights)
  wmn <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$dcmn)*ec)*c(gr$weights)
  wsn <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$dcsn)*ec)*c(gr$weights)
  wtn <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$dctn)*ec)*c(gr$weights)
 }
 hm <- array(0,dim = c(ncol(gr$out$mu[,which(b$mu[i,] != 0),drop=F]),ncol(gr$out$mu[,which(b$mu[i,] != 0),drop=F]),length(gr$weights)))
 for(z in 1:nrow(gr$points)){ hm[,,z] <- crossprod(as.matrix(gr$out$mu[z,which(b$mu[i,] != 0),drop=F]))*wmu[z] }
 He <- hm <- rowSums(hm, dims = 2)
 if("sigma" %in% pFun(fam[i])){
  hs <- array(0,dim = c(ncol(gr$out$sigma[,which(b$sigma[i,] != 0),drop=F]),ncol(gr$out$sigma[,which(b$sigma[i,] != 0),drop=F]),length(gr$weights)))
  hms <- array(0,dim = c(ncol(gr$out$mu[,which(b$mu[i,] != 0)]),ncol(gr$out$sigma[,which(b$sigma[i,] != 0),drop=F]),length(gr$weights)))
  for(z in 1:nrow(gr$points)){
   hs[,,z] <- crossprod(as.matrix(gr$out$sigma[z,which(b$sigma[i,] != 0),drop=F]))*wsg[z]
   hms[,,z] <- t(as.matrix(gr$out$mu[z,which(b$mu[i,] != 0),drop=F])) %*% as.matrix(gr$out$sigma[z,which(b$sigma[i,] != 0),drop=F])*wms[z]
  }
  hs <- rowSums(hs, dims = 2)
  hms <- rowSums(hms, dims = 2)
  He <- cbind(rbind(hm,t(hms)),rbind(hms,hs))
 }
 if("tau" %in% pFun(fam[i])){
  ht <- array(0,dim = c(ncol(gr$out$tau[,which(b$tau[i,] != 0),drop=F]),ncol(gr$out$tau[,which(b$tau[i,] != 0),drop=F]),length(gr$weights)))
  hmt <- array(0,dim = c(ncol(gr$out$mu[,which(b$mu[i,] != 0),drop=F]),ncol(gr$out$tau[,which(b$tau[i,] != 0),drop=F]),length(gr$weights)))
  hst <- array(0,dim = c(ncol(gr$out$sigma[,which(b$sigma[i,] != 0),drop=F]),ncol(gr$out$tau[,which(b$tau[i,] != 0),drop=F]),length(gr$weights)))
  for(z in 1:nrow(gr$points)){
   ht[,,z] <- crossprod(as.matrix(gr$out$tau[z,which(b$tau[i,] != 0),drop=F]))*wta[z]
   hmt[,,z] <- t(as.matrix(gr$out$mu[z,which(b$mu[i,] != 0),drop=F])) %*% as.matrix(gr$out$tau[z,which(b$tau[i,] != 0),drop=F])*wmt[z]
   hst[,,z] <- t(as.matrix(gr$out$sigma[z,which(b$sigma[i,] != 0),drop=F])) %*% as.matrix(gr$out$tau[z,which(b$tau[i,] != 0),drop=F])*wst[z]
  }
  ht <- rowSums(ht, dims = 2)
  hmt <- rowSums(hmt, dims = 2)
  hst <- rowSums(hst, dims = 2)
  He <- cbind(rbind(He,cbind(t(hmt),t(hst))),rbind(hmt,hst,ht))
 }
 if("nu" %in% pFun(fam[i])){
  hn <- array(0,dim = c(ncol(gr$out$nu[,which(b$nu[i,] != 0),drop=F]),ncol(gr$out$nu[,which(b$nu[i,] != 0),drop=F]),length(gr$weights)))
  hmn <- array(0,dim = c(ncol(gr$out$mu[,which(b$mu[i,] != 0),drop=F]),ncol(gr$out$nu[,which(b$nu[i,] != 0),drop=F]),length(gr$weights)))
  hsn <- array(0,dim = c(ncol(gr$out$sigma[,which(b$sigma[i,] != 0),drop=F]),ncol(gr$out$nu[,which(b$nu[i,] != 0),drop=F]),length(gr$weights)))
  htn <- array(0,dim = c(ncol(gr$out$tau[,which(b$tau[i,] != 0),drop=F]),ncol(gr$out$nu[,which(b$nu[i,] != 0),drop=F]),length(gr$weights)))
  for(z in 1:nrow(gr$points)){
   hn[,,z] <- crossprod(as.matrix(gr$out$nu[z,which(b$nu[i,] != 0),drop=F]))*wnu[z]
   hmn[,,z] <- t(as.matrix(gr$out$mu[z,which(b$mu[i,] != 0),drop=F])) %*% as.matrix(gr$out$nu[z,which(b$nu[i,] != 0),drop=F])*wmn[z]
   hsn[,,z] <- t(as.matrix(gr$out$sigma[z,which(b$sigma[i,] != 0),drop=F])) %*% as.matrix(gr$out$nu[z,which(b$nu[i,] != 0),drop=F])*wsn[z]
   htn[,,z] <- t(as.matrix(gr$out$tau[z,which(b$tau[i,] != 0),drop=F])) %*% as.matrix(gr$out$nu[z,which(b$nu[i,] != 0),drop=F])*wtn[z]
  }
  hn <- rowSums(hn, dims = 2)
  hmn <- rowSums(hmn, dims = 2)
  hsn <- rowSums(hsn, dims = 2)
  htn <- rowSums(htn, dims = 2)
  He <- cbind(rbind(He,cbind(t(hmn),t(hsn),t(htn))),rbind(hmn,hsn,htn,hn))
 }
} else {
 dimi <- ncol(gr$out$mu[,which(b$mu[i,] != 0),drop=F])
 dimj <- ncol(gr$out$mu[,which(b$mu[j,] != 0),drop=F])
 if("sigma" %in% pFun(fam[i])) dimi <- dimi + ncol(gr$out$sigma[,which(b$sigma[i,] != 0),drop=F])
 if("sigma" %in% pFun(fam[j])) dimj <- dimj + ncol(gr$out$sigma[,which(b$sigma[j,] != 0),drop=F])
 if("tau" %in% pFun(fam[i])) dimi <- dimi + ncol(gr$out$tau[,which(b$tau[i,] != 0),drop=F])
 if("tau" %in% pFun(fam[j])) dimj <- dimj + ncol(gr$out$tau[,which(b$tau[j,] != 0),drop=F])
 if("nu" %in% pFun(fam[i])) dimi <- dimi + ncol(gr$out$nu[,which(b$nu[i,] != 0),drop=F])
 if("nu" %in% pFun(fam[j])) dimj <- dimj + ncol(gr$out$nu[,which(b$nu[j,] != 0),drop=F])
 He <- matrix(0,dimi,dimj)
}
 return(t(He))
}

bhe2 <- function(i,j,Y,b,gr,fam,ec){
 # To evaluate:
 # i = 1; Y = simR$Y; b = ex2$b; gr = ex2$gr; fam = fam; ec = EC; z = 10
 if(!is.matrix(Y)) Y <- as.matrix(Y) 
 mi <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1mu); listi <- "mi"
 mj <- sapply(1:nrow(gr$points), function(z) dvFun(j,z,fam[j],Y[,j],gr$out,b)$d1mu); listj <- "mj"
 ri<- unname(as.matrix(gr$out$mu[,which(b$mu[i,] != 0)]))
 rj <- unname(as.matrix(gr$out$mu[,which(b$mu[j,] != 0)]))
 imi <- 1:ncol(ri); idxi <- "imi"
 imj <- 1:ncol(rj); idxj <- "imj"
 for(k in c("i","j")){
  if("sigma" %in% pFun(fam[get(k)])){
   assign(paste0("s",k), sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1sg))
   assign(paste0("list",k), append(get(paste0("list",k)), paste0("s",k)))
   assign(paste0("r",k), cbind(get(paste0("r",k)),unname(as.matrix(gr$out$sigma[,which(b$sigma[get(k),] != 0)]))))
   assign(paste0("is",k), c(1:ncol(get(paste0("r",k))))[c(1:ncol(get(paste0("r",k)))) %!in% get(paste0("im",k))])
   assign(paste0("idx",k), append(get(paste0("idx",k)), paste0("is",k)))
  }
  if("tau" %in% pFun(fam[get(k)])){
   assign(paste0("t",k), sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1ta))
   assign(paste0("list",k), append(get(paste0("list",k)), paste0("t",k)))
   assign(paste0("r",k), cbind(get(paste0("r",k)),unname(as.matrix(gr$out$tau[,which(b$tau[get(k),] != 0)]))))
   assign(paste0("it",k), c(1:ncol(get(paste0("r",k))))[c(1:ncol(get(paste0("r",k)))) %!in% c(get(paste0("im",k)),get(paste0("is",k)))])
   assign(paste0("idx",k), append(get(paste0("idx",k)), paste0("it",k)))
  }
  if("nu" %in% pFun(fam[get(k)])){
   assign(paste0("n",k), sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1nu))
   assign(paste0("list",k), append(get(paste0("list",k)), paste0("n",k)))
   assign(paste0("r",k), cbind(get(paste0("r",k)),unname(as.matrix(gr$out$nu[,which(b$nu[get(k),] != 0)]))))
   assign(paste0("in",k), c(1:ncol(get(paste0("r",k))))[c(1:ncol(get(paste0("r",k)))) %!in% c(get(paste0("im",k)),get(paste0("is",k)), get(paste0("it",k)))])
   assign(paste0("idx",k), append(get(paste0("idx",k)), paste0("in",k)))
  }
 }
 idM <- unname(as.matrix(expand.grid(1:length(listi), 1:length(listj))))
 hm <- wij <- array(0,dim = c(ncol(ri),ncol(rj),length(gr$weights)))
 for(z in 1:nrow(gr$points)){
  for(zz in 1:nrow(idM)){
   tidi <- get(idxi[idM[zz,][1]])
   tidj <- get(idxj[idM[zz,][2]])
   wij[tidi,tidj,z] <- c(colSums(get(listi[idM[zz,][1]])*get(listj[idM[zz,][2]])*ec)*c(gr$weights))[z]
   hm[tidi,tidj,z] <- tcrossprod(ri[z,tidi],rj[z,tidj]) * wij[tidi,tidj,z]
  }
 }
 return(t(rowSums(hm, dims = 2)))
}

bhe3 <- function(i,j,Y,b,gr,fam,ec){
 # To evaluate:
 # i = 1; Y = simR$Y; b = bold; gr = gr; fam = fam; ec = EC; z = 10
 if(!is.matrix(Y)) Y <- as.matrix(Y)
 mi <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1mu)*ec
 mj <- sapply(1:nrow(gr$points), function(z) dvFun(j,z,fam[j],Y[,j],gr$out,b)$d1mu)*ec
 dimi <- ncol(gr$out$mu[,which(b$mu[i,] != 0)])
 dimj <- ncol(gr$out$mu[,which(b$mu[j,] != 0)])
 for(k in c("i","j")){
  if("sigma" %in% pFun(fam[get(k)])){
   assign(paste0("s",k), sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1sg)*ec)
   assign(paste0("dim",k), get(paste0("dim",k)) + ncol(gr$out$sigma[,which(b$sigma[get(k),] != 0)]))
  }
  if("tau" %in% pFun(fam[get(k)])){
   assign(paste0("t",k), sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1ta)*ec)
   assign(paste0("dim",k), get(paste0("dim",k)) + ncol(gr$out$tau[,which(b$tau[get(k),] != 0)])) 
  }
  if("nu" %in% pFun(fam[get(k)])){
   assign(paste0("n",k), sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1nu)*ec)
   assign(paste0("dim",k), get(paste0("dim",k)) + ncol(gr$out$nu[,which(b$nu[get(k),] != 0)]))
  }
 }
 hm <- array(0,dim = c(dimi, dimj, nrow(Y)))
 for(m in 1:nrow(Y)){
  ri <- colSums(gr$out$mu[,which(b$mu[i,] != 0)]*mi[m,]*c(gr$weights))
  rj <- colSums(gr$out$mu[,which(b$mu[j,] != 0)]*mj[m,]*c(gr$weights))
  for(k in c("i","j")){
   if("sigma" %in% pFun(fam[get(k)])) assign(paste0("r",k), c(get(paste0("r",k)), colSums(gr$out$sigma[,which(b$sigma[get(k),] != 0)]*get(paste0("s",k))[m,]*c(gr$weights))))
   if("tau" %in% pFun(fam[get(k)])) assign(paste0("r",k), c(get(paste0("r",k)), colSums(gr$out$tau[,which(b$tau[get(k),] != 0)]*get(paste0("t",k))[m,]*c(gr$weights))))
   if("nu" %in% pFun(fam[get(k)])) assign(paste0("r",k), c(get(paste0("r",k)), colSums(gr$out$nu[,which(b$nu[get(k),] != 0)]*get(paste0("n",k))[m,]*c(gr$weights))))
  }
  hm[,,m] <- tcrossprod(ri,rj)
 }
 return(t(rowSums(hm, dims = 2)))
}

bhe3a <- function(i,j,Y,b,gr,fam,ec){
# To evaluate:
# i = 1; Y = Y; b = ex1$b; gr = ex1$gr; fam = fam
# ec <- sapply(1:nrow(gr$points), function(z) mfyz(z,Y,gr$out,b,fam))/mfy(Y,b,gr,fam)
  
if(!is.matrix(Y)) Y <- as.matrix(Y) 
mi <- t(t(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1mu)*ec)*c(gr$weights))
mj <- t(t(sapply(1:nrow(gr$points), function(z) dvFun(j,z,fam[j],Y[,j],gr$out,b)$d1mu)*ec)*c(gr$weights))
ri <- unname(as.matrix(gr$out$mu[,which(b$mu[i,] != 0)]))
rj <- unname(as.matrix(gr$out$mu[,which(b$mu[j,] != 0)]))
imi <- 1:ncol(ri)
imj <- 1:ncol(rj)
for(k in c("i","j")){
 if("sigma" %in% pFun(fam[get(k)])){
  assign(paste0("s",k), t(t(sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1sg)*ec)*c(gr$weights)))
  assign(paste0("r",k), cbind(get(paste0("r",k)),unname(as.matrix(gr$out$sigma[,which(b$sigma[get(k),] != 0)]))))
  assign(paste0("is",k), c(1:ncol(get(paste0("r",k))))[c(1:ncol(get(paste0("r",k)))) %!in% get(paste0("im",k))])
 }
 if("tau" %in% pFun(fam[get(k)])){
  assign(paste0("t",k), t(t(sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1ta)*ec)*c(gr$weights)))
  assign(paste0("r",k), cbind(get(paste0("r",k)),unname(as.matrix(gr$out$tau[,which(b$tau[get(k),] != 0)]))))
  assign(paste0("it",k), c(1:ncol(get(paste0("r",k))))[c(1:ncol(get(paste0("r",k)))) %!in% c(get(paste0("im",k)),get(paste0("is",k)))])
 }
 if("nu" %in% pFun(fam[get(k)])){
  assign(paste0("n",k), t(t(sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1nu)*ec)*c(gr$weights)))
  assign(paste0("r",k), cbind(get(paste0("r",k)),unname(as.matrix(gr$out$nu[,which(b$nu[get(k),] != 0)]))))
  assign(paste0("in",k), c(1:ncol(get(paste0("r",k))))[c(1:ncol(get(paste0("r",k)))) %!in% c(get(paste0("im",k)),get(paste0("is",k)), get(paste0("it",k)))])
 }
}
wi <- array(0,dim = c(nrow(gr$points),ncol(ri),nrow(Y)))
wj <- array(0,dim = c(nrow(gr$points),ncol(rj),nrow(Y)))
for(m in 1:nrow(Y)){
 wi[,imi,m] <- matrix(rep(mi[m,],length(imi)),nrow = ncol(mi))
 if("sigma" %in% pFun(fam[i])) wi[,isi,m] <- matrix(rep(si[m,],length(isi)),nrow = ncol(si))
 if("tau" %in% pFun(fam[i])) wi[,iti,m] <- matrix(rep(ti[m,],length(iti)),nrow = ncol(ti))
 if("nu" %in% pFun(fam[i])) wi[,ini,m] <- matrix(rep(ni[m,],length(ini)),nrow = ncol(ni))
 wj[,imj,m] <- matrix(rep(mj[m,],length(imj)),nrow = ncol(mj))
 if("sigma" %in% pFun(fam[j])) wj[,isj,m] <- matrix(rep(sj[m,],length(isj)),nrow = ncol(sj))
 if("tau" %in% pFun(fam[j])) wj[,itj,m] <- matrix(rep(tj[m,],length(itj)),nrow = ncol(tj))
 if("nu" %in% pFun(fam[j])) wj[,inj,m] <- matrix(rep(nj[m,],length(inj)),nrow = ncol(nj))
}
zs <- array(sapply(1:nrow(Y), function(m) tcrossprod(colSums(ri*wi[,,m]),colSums(rj*wj[,,m]))),
            dim = c(ncol(ri),ncol(rj),nrow(Y)))
return(t(rowSums(zs,dim = 2)))
}

bhe <- function(i,j,Y,b,gr,fam,ec){
  return(bhe1(i,j,Y,b,gr,fam,ec) + bhe2(i,j,Y,b,gr,fam,ec) - bhe3a(i,j,Y,b,gr,fam,ec))
}

# ec <- sapply(1:nrow(gr$points), function(z) mfyz(z,Y,gr$out,b,fam))/mfy(Y,b,gr,fam)

h1he <- function(i,Y,b,gr,fam,ec){
  
  # To evaluate:
  # efy <- mfy(Y,b,gr,fam)
  # EC <- sapply(1:nrow(gr$points), function(z) mfyz(z,Y,gr$out,b,fam))/efy 
  # i = 1; Y = simR$Y; b = bold; gr = gr; fam = fam; ec = EC; z = 10
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  ed <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1mu)
  dim <- ncol(gr$out$mu)
  tmpo <- gr$out$mu
  if("sigma" %in% pFun(fam[i])){ 
   ed <- ed*sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1sg)
   dim <- dim + ncol(gr$out$sigma)
   tmpo <- cbind(tmpo, gr$out$sigma) }
  if("tau" %in% pFun(fam[i])){ 
   ed <- ed*sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1ta)
   dim <- dim + ncol(gr$out$tau)
   tmpo <- cbind(tmpo, gr$out$tau) }
  if("nu" %in% pFun(fam[i])){
   ed <- ed*sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1nu)
   dim <- dim + ncol(gr$out$nu)
   tmpo <- cbind(tmpo, gr$out$nu) }
  t1 <- colSums(ec*ed)*drop(gr$weights)
  hm <- array(0,dim = c(dim,dim,length(gr$weights)))
  for(z in 1:nrow(gr$points)){ hm[,,z] <- crossprod(as.matrix(tmpo[z,]))*t1[z] }
  hm <- rowSums(hm, dims = 2)
  return(as.matrix(hm))
}

h2he <- function(i,Y,b,gr,fam,ec){
  
  # To evaluate:
  # efy <- mfy(Y,b,gr,fam)
  # EC <- sapply(1:nrow(gr$points), function(z) mfyz(z,Y,gr$out,b,fam))/efy 
  # i = 1; Y = simR$Y; b = bold; gr = gr; fam = fam; ec = EC; z = 10
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  ed <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1mu)
  w1 <- ec*ed*matrix(c(gr$weights),byrow = T,ncol = length(gr$weights), nrow = nrow(Y))
  dim <- ncol(gr$out$mu)
  tmpo <- gr$out$mu
  if("sigma" %in% pFun(fam[i])){ 
    ed <- ed*sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1sg)
    dim <- dim + ncol(gr$out$sigma)
    tmpo <- cbind(tmpo, gr$out$sigma) }
  if("tau" %in% pFun(fam[i])){ 
    ed <- ed*sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1ta)
    dim <- dim + ncol(gr$out$tau)
    tmpo <- cbind(tmpo, gr$out$tau) }
  if("nu" %in% pFun(fam[i])){
    ed <- ed*sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1nu)
    dim <- dim + ncol(gr$out$nu)
    tmpo <- cbind(tmpo, gr$out$nu) }
  
  t1 <- colSums(ec*ed)*drop(gr$weights)
  hm <- array(0,dim = c(dim,dim,length(gr$weights)))
  for(z in 1:nrow(gr$points)){ hm[,,z] <- crossprod(as.matrix(tmpo[z,]))*t1[z]^2 }
  
  # t1 <- rowSums(ec*ed*matrix(c(gr$weights),byrow = T,ncol = length(gr$weights), nrow = nrow(Y)))
  # hm <- array(0,dim = c(ncol(gr$out$sigma),ncol(gr$out$sigma),length(gr$weights)))
  # for(z in 1:nrow(gr$points)){ hm[,,z] <- crossprod(as.matrix(gr$out$sigma[z,]))*t1[z] }
  
  hm <- rowSums(hm, dims = 2)
  return(as.matrix(hm))
}

# rbenchmark::benchmark(
#   "bhe" = {bhe(1,2,Y,b,gr,fam,ec)},
#   "hess" = {hess(1,2,ghQ,b,fam,dvL,pD)},
#   replications = 100,
#   columns = c("test", "replications", "elapsed",
#               "relative", "user.self", "sys.self")
# )
