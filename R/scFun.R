
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
 B <- coefmod2(bet = B, beta = beta)
 tmp <- sum(log(mfy(Y,B,ghQ,fam)))
 return(tmp)
}

# Score function (to numerically compare outputs)

ScoreFun <- function(B, Y, beta, ghQ, fam){
 # To evaluate: B = unlist(ex1$b); ghQ = ex1$gr; beta = ex1$b; loadmt = ex1$loadmt
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
 # To evaluate: B = unlist(ex1$b); ghQ = ex1$gr; beta = ex1$b; loadmt = ex1$loadmt
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
 

 # for(ii in 1:ncol(Y)){
 #  HH <- bmhe(ii,Y,B,ghQ,fam,EC)
 #  bo <- B$mu
 #  if("sigma" %in% pFun(fam[ii])){
 #   HH <- as.matrix(Matrix::bdiag(HH, bshe(ii,Y,B,ghQ,fam,EC)))
 #   bo <- cbind(bo,B$sigma)
 #   idxsg <- c(1:ncol(HH))[1:ncol(HH) %!in% idxmu]
 #   HH[idxmu,idxsg] <- mshe(ii,Y,B,ghQ,fam,EC)
 #   HH[idxsg,idxmu] <- t(mshe(ii,Y,B,ghQ,fam,EC))
 #  }
 #  if("tau" %in% pFun(fam[ii])){
 #   HH <- as.matrix(Matrix::bdiag(HH, bthe(ii,Y,B,ghQ,fam,EC)))
 #   bo <- cbind(bo,B$tau)
 #   idxta <- c(1:ncol(HH))[1:ncol(HH) %!in% c(idxmu,idxsg)]
 #   HH[idxmu,idxta] <- mthe(ii,Y,B,ghQ,fam,EC)
 #   HH[idxta,idxmu] <- t(mthe(ii,Y,B,ghQ,fam,EC))
 #   HH[idxsg,idxta] <- sthe(ii,Y,B,ghQ,fam,EC)
 #   HH[idxta,idxsg] <- t(sthe(ii,Y,B,ghQ,fam,EC))
 #  }
 #  if("nu" %in% pFun(fam[ii])){
 #   HH <- as.matrix(Matrix::bdiag(HH, bnhe(ii,Y,B,ghQ,fam,EC)))
 #   bo <- cbind(bo,B$nu)
 #   idxnu <- c(1:ncol(HH))[1:ncol(HH) %!in% c(idxmu,idxsg,idxta)]
 #   HH[idxmu,idxnu] <- mnhe(ii,Y,B,ghQ,fam,EC)
 #   HH[idxnu,idxmu] <- t(mnhe(ii,Y,B,ghQ,fam,EC))
 #   HH[idxsg,idxnu] <- snhe(ii,Y,B,ghQ,fam,EC)
 #   HH[idxnu,idxsg] <- t(snhe(ii,Y,B,ghQ,fam,EC))
 #   HH[idxta,idxnu] <- tnhe(ii,Y,B,ghQ,fam,EC)
 #   HH[idxnu,idxta] <- t(tnhe(ii,Y,B,ghQ,fam,EC))
 #  }
 #  indx <- seq(from = (ii-1)*ncol(HH)+1, length.out = ncol(HH))
 #  HH <- HH + h1he(ii,Y,B,ghQ,fam,EC) - h2he(ii,Y,B,ghQ,fam,EC) # tcrossprod(ScoreFun(c(t(bo)), Y, B, ghQ, fam)[indx])
 #  Hm <- Matrix::bdiag(Hm,HH)
 # }
 return(as.matrix(Hm)) # return(as.matrix(Hm)[-1,-1])
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
 if(!is.matrix(Y)) Y <- as.matrix(Y)
 wmu <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1mu)*ec)*c(gr$weights)
 if("sigma" %in%  pFun(fam[i])) wsg <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1sg)*ec)*c(gr$weights)
 if("tau" %in%  pFun(fam[i])) wta <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1ta)*ec)*c(gr$weights)
 if("nu" %in%  pFun(fam[i])) wnu <- colSums(sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1nu)*ec)*c(gr$weights)
 ps <- gr$out$mu * wmu
 if("sigma" %in% pFun(fam[i])) ps <- cbind(ps,gr$out$sigma * wsg)
 if("tau" %in% pFun(fam[i])) ps <- cbind(ps,gr$out$tau * wta)
 if("nu" %in% pFun(fam[i])) ps <- cbind(ps,gr$out$nu * wnu)
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
# i = 1; Y = simR$Y; b = bold; gr = gr; fam = fam; ec = EC; z = 10
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
 hm <- array(0,dim = c(ncol(gr$out$mu),ncol(gr$out$mu),length(gr$weights)))
 for(z in 1:nrow(gr$points)){ hm[,,z] <- crossprod(as.matrix(gr$out$mu[z,]))*wmu[z] }
 He <- hm <- rowSums(hm, dims = 2)
 if("sigma" %in% pFun(fam[i])){
  hs <- array(0,dim = c(ncol(gr$out$sigma),ncol(gr$out$sigma),length(gr$weights)))
  hms <- array(0,dim = c(ncol(gr$out$mu),ncol(gr$out$sigma),length(gr$weights)))
  for(z in 1:nrow(gr$points)){
   hs[,,z] <- crossprod(as.matrix(gr$out$sigma[z,]))*wsg[z]
   hms[,,z] <- t(as.matrix(gr$out$mu[z,,drop=F])) %*% as.matrix(gr$out$sigma[z,,drop=F])*wms[z]
  }
  hs <- rowSums(hs, dims = 2)
  hms <- rowSums(hms, dims = 2)
  He <- cbind(rbind(hm,t(hms)),rbind(hms,hs))
 }
 if("tau" %in% pFun(fam[i])){
  ht <- array(0,dim = c(ncol(gr$out$tau),ncol(gr$out$tau),length(gr$weights)))
  hmt <- array(0,dim = c(ncol(gr$out$mu),ncol(gr$out$tau),length(gr$weights)))
  hst <- array(0,dim = c(ncol(gr$out$sigma),ncol(gr$out$tau),length(gr$weights)))
  for(z in 1:nrow(gr$points)){
   ht[,,z] <- crossprod(as.matrix(gr$out$tau[z,]))*wta[z]
   hmt[,,z] <- t(as.matrix(gr$out$mu[z,,drop=F])) %*% as.matrix(gr$out$tau[z,,drop=F])*wmt[z]
   hst[,,z] <- t(as.matrix(gr$out$sigma[z,,drop=F])) %*% as.matrix(gr$out$tau[z,,drop=F])*wst[z]
  }
  ht <- rowSums(ht, dims = 2)
  hmt <- rowSums(hmt, dims = 2)
  hst <- rowSums(hst, dims = 2)
  He <- cbind(rbind(He,cbind(t(hmt),t(hst))),rbind(hmt,hst,ht))
 }
 if("nu" %in% pFun(fam[i])){
  hn <- array(0,dim = c(ncol(gr$out$nu),ncol(gr$out$nu),length(gr$weights)))
  hmn <- array(0,dim = c(ncol(gr$out$mu),ncol(gr$out$nu),length(gr$weights)))
  hsn <- array(0,dim = c(ncol(gr$out$sigma),ncol(gr$out$nu),length(gr$weights)))
  htn <- array(0,dim = c(ncol(gr$out$tau),ncol(gr$out$nu),length(gr$weights)))
  for(z in 1:nrow(gr$points)){
   hn[,,z] <- crossprod(as.matrix(gr$out$nu[z,]))*wnu[z]
   hmn[,,z] <- t(as.matrix(gr$out$mu[z,,drop=F])) %*% as.matrix(gr$out$nu[z,,drop=F])*wmn[z]
   hsn[,,z] <- t(as.matrix(gr$out$sigma[z,,drop=F])) %*% as.matrix(gr$out$nu[z,,drop=F])*wsn[z]
   htn[,,z] <- t(as.matrix(gr$out$tau[z,,drop=F])) %*% as.matrix(gr$out$nu[z,,drop=F])*wtn[z]
  }
  hn <- rowSums(hn, dims = 2)
  hmn <- rowSums(hmn, dims = 2)
  hsn <- rowSums(hsn, dims = 2)
  htn <- rowSums(htn, dims = 2)
  He <- cbind(rbind(He,cbind(t(hmn),t(hsn),t(htn))),rbind(hmn,hsn,htn,hn))
 }
} else {
 dimi <- dimj <- ncol(gr$out$mu)
 if("sigma" %in% pFun(fam[i])) dimi <- dimi + ncol(gr$out$sigma)
 if("sigma" %in% pFun(fam[j])) dimj <- dimj + ncol(gr$out$sigma)
 if("tau" %in% pFun(fam[i])) dimi <- dimi + ncol(gr$out$tau)
 if("tau" %in% pFun(fam[j])) dimj <- dimj + ncol(gr$out$tau)
 if("nu" %in% pFun(fam[i])) dimi <- dimi + ncol(gr$out$nu)
 if("nu" %in% pFun(fam[j])) dimj <- dimj + ncol(gr$out$nu)
 He <- matrix(0,dimi,dimj)
}
 return(He)
}

bhe2 <- function(i,j,Y,b,gr,fam,ec){
 # To evaluate:
 # i = 1; Y = simR$Y; b = bold; gr = gr; fam = fam; ec = EC; z = 10
 if(!is.matrix(Y)) Y <- as.matrix(Y) 
 mi <- sapply(1:nrow(gr$points), function(z) dvFun(i,z,fam[i],Y[,i],gr$out,b)$d1mu); listi <- "mi"
 mj <- sapply(1:nrow(gr$points), function(z) dvFun(j,z,fam[j],Y[,j],gr$out,b)$d1mu); listj <- "mj"
 ri<- rj <- unname(as.matrix(gr$out$mu))
 imi <- imj <- 1:ncol(ri); idxi <- "imi"; idxj <- "imj"
 for(k in c("i","j")){
  if("sigma" %in% pFun(fam[get(k)])){
   assign(paste0("s",k), sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1sg))
   assign(paste0("list",k), append(get(paste0("list",k)), paste0("s",k)))
   assign(paste0("r",k), cbind(get(paste0("r",k)),unname(as.matrix(gr$out$sigma))))
   assign(paste0("is",k), c(1:ncol(get(paste0("r",k))))[c(1:ncol(get(paste0("r",k)))) %!in% get(paste0("im",k))])
   assign(paste0("idx",k), append(get(paste0("idx",k)), paste0("is",k)))
  }
  if("tau" %in% pFun(fam[get(k)])){
   assign(paste0("t",k), sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1ta))
   assign(paste0("list",k), append(get(paste0("list",k)), paste0("t",k)))
   assign(paste0("r",k), cbind(get(paste0("r",k)),unname(as.matrix(gr$out$tau))))
   assign(paste0("it",k), c(1:ncol(get(paste0("r",k))))[c(1:ncol(get(paste0("r",k)))) %!in% c(get(paste0("im",k)),get(paste0("is",k)))])
   assign(paste0("idx",k), append(get(paste0("idx",k)), paste0("it",k)))
  }
  if("nu" %in% pFun(fam[get(k)])){
   assign(paste0("n",k), sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1nu))
   assign(paste0("list",k), append(get(paste0("list",k)), paste0("n",k)))
   assign(paste0("r",k), cbind(get(paste0("r",k)),unname(as.matrix(gr$out$nu))))
   assign(paste0("it",k), c(1:ncol(get(paste0("r",k))))[c(1:ncol(get(paste0("r",k)))) %!in% c(get(paste0("im",k)),get(paste0("is",k)), get(paste0("it",k)))])
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
 dimi <- dimj <- ncol(gr$out$mu)
 for(k in c("i","j")){
  if("sigma" %in% pFun(fam[get(k)])){
   assign(paste0("s",k), sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1sg)*ec)
   assign(paste0("dim",k), get(paste0("dim",k)) + ncol(gr$out$sigma))
  }
  if("tau" %in% pFun(fam[get(k)])){
   assign(paste0("t",k), sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1ta)*ec)
   assign(paste0("dim",k), get(paste0("dim",k)) + ncol(gr$out$tau)) 
  }
  if("nu" %in% pFun(fam[get(k)])){
   assign(paste0("n",k), sapply(1:nrow(gr$points), function(z) dvFun(get(k),z,fam[get(k)],Y[,get(k)],gr$out,b)$d1nu)*ec)
   assign(paste0("dim",k), get(paste0("dim",k)) + ncol(gr$out$nu))
  }
 }
 hm <- array(0,dim = c(dimi, dimj, nrow(Y)))
 for(m in 1:nrow(Y)){
  ri <- colSums(gr$out$mu*mi[m,]*c(gr$weights))
  rj <- colSums(gr$out$mu*mj[m,]*c(gr$weights))
  for(k in c("i","j")){
   if("sigma" %in% pFun(fam[get(k)])) assign(paste0("r",k), c(get(paste0("r",k)), colSums(gr$out$sigma*get(paste0("s",k))[m,]*c(gr$weights))))
   if("tau" %in% pFun(fam[get(k)])) assign(paste0("r",k), c(get(paste0("r",k)), colSums(gr$out$tau*get(paste0("t",k))[m,]*c(gr$weights))))
   if("nu" %in% pFun(fam[get(k)])) assign(paste0("r",k), c(get(paste0("r",k)), colSums(gr$out$nu*get(paste0("n",k))[m,]*c(gr$weights))))
  }
  hm[,,m] <- tcrossprod(ri,rj)
 }
 return(t(rowSums(hm, dims = 2)))
}

bhe <- function(i,j,Y,b,gr,fam,ec){
  return(bhe1(i,j,Y,b,gr,fam,ec) + bhe2(i,j,Y,b,gr,fam,ec) - bhe3(i,j,Y,b,gr,fam,ec))
}


# ec <- sapply(1:nrow(gr$points), function(z) mfyz(z,Y,gr$out,b,fam))/mfy(Y,b,gr,fam)

# 
# 
# round(ana.hess1[1:3,1:3],3)
# round(num.hess1[1:9,1:9],3)
# round(bhe(1,1,Y,b,gr,fam,ec),3)

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
#   "diag" = {diag(crossprod(ec,ed))},
#   "mat" = {colSums(ec*ed)},
#   replications = 10000,
#   columns = c("test", "replications", "elapsed",
#               "relative", "user.self", "sys.self")
# )
