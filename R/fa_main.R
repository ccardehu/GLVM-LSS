
##############################################
## SPLVM: Matrices functions for estimation ##
##############################################

dY <- function(Y,ghQ,b,fam){

# Goal: To compute (conditional) density of items in (nxp) matrix Y
# Input : Y (item matrix), ghQ (GHQ object), b (loadings matrix), fam (distributions)
# Output: Array n x qp x p with f(y_p|z_qp)
# Testing: Y = simR$Y; ghQ = ghQ; b = borg; fam = fam

if(!is.matrix(Y)) Y <- as.matrix(Y)
fyz <- array(NA, dim = c(nrow(Y),nrow(ghQ$points),ncol(Y)))
Z <- ghQ$out
for(z in 1:nrow(ghQ$points)){
for(i in 1:ncol(Y)){
if(fam[i] == "normal"){
 mu = c(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
 sigma = c(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
 fyz[,z,i] <- dnorm(Y[,i], mu, sigma, log = T)
}
if(fam[i] == "lognormal"){
 mu = c(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
 sigma = exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,]))
 fyz[,z,i] <- dlnorm(Y[,i], mu, sigma, log = T)
}
if(fam[i] == "poisson"){
 mu = c(exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
 fyz[,z,i] <- dpois(Y[,i], mu, log = T)
}
if(fam[i] == "gamma"){
 mu = c(exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
 sigma = c(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
 fyz[,z,i] <- dgamma(Y[,i],shape = mu, scale = sigma, log = T)
}
if(fam[i] == "binomial"){
 mu = c(probs(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
 fyz[,z,i] <- dbinom(Y[,i],size = 1,prob = mu, log = T)
}
if(fam[i] == "ZIpoisson"){
 mu = c(exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
 sigma = c(probs(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
 dZIpoisson <- function(Y,mu,sigma,log = T){
  u <- as.numeric(Y == 0)
  lf <- u*c(log(sigma + (1-sigma)*exp(-mu))) + (1-u)*(c(log(1-sigma)) - c(mu) + Y*c(log(mu)) - lfactorial(Y))
  if(log == T) return(lf) else return(exp(lf))
 }
 fyz[,z,i] <- dZIpoisson(Y[,i], mu, sigma, log = T)
}
  # Add other distributions
} }

return(fyz)

}

dvY <- function(Y,ghQ,b,fam,information = "Fisher"){

# Goal: To compute first and second derivatives of log-likelihood Y
# Input : Y (item matrix), ghQ (GHQ object), b (loadings matrix), fam (distributions)
# Output: List of arrays of dims n x qp x p with d-logf(y_p|z_qp) / d-eta_p,phi
# Testing: Y = simR$Y; ghQ = ghQ; b = borg; fam = fam

info.F <- match.arg(information,c("Fisher","Hessian")) == "Fisher"
if(!is.matrix(Y)) Y <- as.matrix(Y)
Z <- ghQ$out
dvY <- vector(mode = "list", 14)  
names(dvY) <- c("d1mu","d1sg","d1ta","d1nu","d2mu","d2sg","d2ta","d2nu","dcms",
                "dcmt","dcmn","dcst","dcsn","dctn")
for(i in 1:length(dvY)){ 
  dvY[[i]] <- array(NA,dim = c(nrow(Y),nrow(ghQ$points),ncol(Y))) }

for(z in 1:nrow(ghQ$points)){
for(i in 1:ncol(Y)){
if(fam[i] == "normal"){ # structure(gamlss.dist::NO)
 mu = c(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
 sigma = c(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
 dvY$d1mu[,z,i] = (Y[,i] - mu)/(sigma^2) * 1
 dvY$d2mu[,z,i] = rep(-(1/sigma^2) * 1, length(Y[,i])) + 0
 dvY$d1sg[,z,i] = (((Y[,i] - mu)^2 - sigma^2)/(sigma^3)) * sigma
 if(info.F) dvY$d2sg[,z,i] = rep(-(2/(sigma^2)), length(Y[,i]))* sigma^2 else
  dvY$d2sg[,z,i] = (-3*(Y[,i]-mu)^2/sigma^4 + 1/sigma^2) * sigma^2 + dvY$d1sg[,z,i] # rep(-(2/(sigma^2)), length(Y))* sigma^2 # 
 if(info.F) dvY$dcms[,z,i] = rep(0, length(Y[,i])) else
  dvY$dcms[,z,i] = -2*(Y[,i]-mu)/sigma^3 * 1 * sigma # rep(0, length(Y)) # the expected cross derivative mu and sigma
}
if(fam[i] == "lognormal"){ # structure(gamlss.dist::LOGNO)
 mu = c(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
 sigma = c(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
 dvY$d1mu[,z,i] = (log(Y[,i]) - mu)/(sigma^2) * 1
 dvY$d2mu[,z,i] = rep(-(1/sigma^2) * 1, length(Y[,i])) + 0
 dvY$d1sg[,z,i] = (((log(Y[,i]) - mu)^2 - sigma^2)/(sigma^3)) * sigma
 if(info.F) dvY$d2sg[,z,i] = rep(-(2/(sigma^2)) * sigma^2, length(Y[,i])) else 
  dvY$d2sg[,z,i] = (-3*(log(Y[,i])-mu)^2/sigma^4 + 1/sigma^2) * sigma^2 + dvY$d1sg[,z,i] # rep(-(2/(sigma^2)) * sigma^2, length(Y)) # -3*(Y-mu)^2/sigma^4 + 1/sigma^2
 if(info.F) dvY$dcms[,z,i] = rep(0, length(Y[,i])) else
  dvY$dcms[,z,i] = -2*(log(Y[,i])-mu)/sigma^3 * 1 * sigma # rep(0, length(Y)) # the expected cross derivative mu and sigma
}
if(fam[i] == "poisson"){ # structure(gamlss.dist::PO)
 mu = c(exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
 dvY$d1mu[,z,i] = ((Y[,i] - mu)/mu) * mu
 if(info.F) dvY$d2mu[,z,i] = rep((-1/mu) * mu^2, length(Y[,i])) else
  dvY$d2mu[,z,i] = -Y[,i]/mu^2 * mu^2 + dvY$d1mu[,z,i] # rep((-1/mu) * mu^2, length(Y))
}
if(fam[i] == "gamma"){ # structure(gamlss.dist::GA)
 mu = c(exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
 sigma = c(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
 dvY$d1mu[,z,i] = (log(Y[,i]) - digamma(mu) - log(sigma)) * mu
 dvY$d2mu[,z,i] = rep(-trigamma(mu) * mu^2, length(Y[,i])) + dvY$d1mu[,z,i]
 dvY$d1sg[,z,i] = (Y[,i]/(sigma^2) - mu/sigma) * sigma
 if(info.F) dvY$d2sg[,z,i] = rep((-mu/(sigma)) * sigma^2, length(Y[,i])) else
  dvY$d2sg[,z,i] = ((-2*Y[,i]+mu*sigma)/sigma^3) * sigma^2 + dvY$d1sg[,z,i] # rep((-mu/(sigma)) * sigma^2, length(Y)) # 
 if(info.F) dvY$dcms[,z,i] = rep(0, length(Y[,i])) else
  dvY$dcms[,z,i] = rep(-1/sigma*mu*sigma,length(Y[,i])) # rep(0, length(Y))  # the expected cross derivative mu and sigma
}
if(fam[i] == "binomial"){ # structure(gamlss.dist::BI)
 mu = c(probs(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
 dvY$d1mu[,z,i] = (Y[,i] - mu)/(mu * (1 - mu)) * (mu*(1-mu))
 if(info.F) dvY$d2mu[,z,i] = rep(-(1/(mu * (1 - mu))) * (mu*(1-mu))^2, length(Y[,i])) else
  dvY$d2mu[,z,i] = -(mu^2 - 2*Y[,i]*mu + Y[,i])/((mu-1)*mu)^2 * (mu*(1-mu))^2 + (1-2*mu)*dvY$d1mu[,z,i]# rep(-(1/(mu * (1 - mu))) * (mu*(1-mu))^2, length(Y))
}
if(fam[i] == "ZIpoisson"){ # structure(gamlss.dist::ZIP)
 mu = c(exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
 sigma = c(probs(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
 u = ifelse(Y[,i] == 0, (1 + exp(-logit(sigma)-mu))^(-1),0) # as.numeric(Y == 0)
 dvY$d1mu[,z,i] = (1-u)*(Y[,i]/mu-1) * mu
 if(info.F) dvY$d2mu[,z,i] = (1-u)*(-1/mu^2) * mu^2 else 
  dvY$d2mu[,z,i] = (1-u)*(-Y[,i]/mu^2) * mu^2 + dvY$d1mu[,z,i] # expected value is -1/mu; original is (-Y/mu^2)
 dvY$d1sg[,z,i] = (u-sigma)/(sigma*(1-sigma)) * (sigma*(1-sigma))
 if(info.F) dvY$d2sg[,z,i] = -(u/(sigma^2) + (1-u)/((1-sigma)^2)) * (sigma*(1-sigma))^2 else  # CHECK THIS
  dvY$d2sg[,z,i] = -(sigma^2 - 2*u*sigma + u)/((sigma-1)*sigma)^2 * (sigma*(1-sigma))^2 + (1-2*sigma)*dvY$d1sg[,z,i] # -(u/(sigma^2) + (1-u)/((1-sigma)^2)) * (sigma*(1-sigma))^2
 dvY$dcms[,z,i] = rep(0, length(Y[,i]))
} } }
for(i in names(dvY)){ if(all(is.na(dvY[[i]]))) dvY[[i]] <- NULL } 
return(dvY)
}

sco <- function(i,ghQ,b,fam,dvL,pD){
  
# Goal: To compute Score vector for item i (eq. A1)
# Input : i (item), ghQ (GHQ object), b (loadings matrix),
#         fam (distributions), dvL (list of derivatives),
#         pD (posterior density, to be compute with joint with weights)
# Output: vector of length (total # of parameters for item i)
# Testing: i = 1; Y = simR$Y; ghQ = ghQ; b = borg; fam = fam; 
#          dvL = dvY(Y,ghQ,b,fam)
#          pD = exp(rowSums(dY(Y,ghQ,b,fam),dim = 2)) /
#          c(exp(rowSums(dY(Y,ghQ,b,fam),dim = 2))%*%ghQ$weights)

iN0 <- which(lb2mb(b)[i,] != 0)
Z <- ghQ$out
wmu <- colSums(dvL$d1mu[,,i]*pD)*c(ghQ$weights)
ps <- Z$mu * wmu
if("sigma" %in%  pFun(fam[i])){
 wsg <- colSums(dvL$d1sg[,,i]*pD)*c(ghQ$weights)
 ps <- cbind(ps,Z$sigma * wsg)
}
if("tau" %in%  pFun(fam[i])){
 wta <- colSums(dvL$d1ta[,,i]*pD)*c(ghQ$weights)
 ps <- cbind(ps,Z$tau * wta)
}
if("nu" %in%  pFun(fam[i])){
 wnu <- colSums(dvL$d1nu[,,i]*pD)*c(ghQ$weights)
 ps <- cbind(ps,Z$nu * wnu)
}
ps <- as.matrix(colSums(unname(ps[,iN0])))
return(list("score" = ps, "pos" = iN0))
}

hess <- function(i,j,ghQ,b,fam,dvL,pD){
  
# Goal: To compute Hessian matrix for items i & j (eq. A3)
# Input : i & j (items), ghQ (GHQ object), b (loadings matrix),
#         fam (distributions), dvL (list of derivatives),
#         pD (part of posterior density, to be jointly computed with weights)
# Output: vector of length (total # of parameters for item i)
# Testing: i = 1; Y = simR$Y; ghQ = ghQ; b = borg; fam = fam; 
#          dvL = dvY(Y,ghQ,b,fam,control$information)
#          pD = exp(rowSums(dY(Y,ghQ,b,fam),dim = 2)) /
#          c(exp(rowSums(dY(Y,ghQ,b,fam),dim = 2))%*%ghQ$weights)

Z <- ghQ$out

if(i == j){
iN0 <- jN0 <- which(lb2mb(b)[i,] != 0)
# For first Hessian (eq. A3, row 1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Zmu <- unname(as.matrix(Z$mu[,which(b$mu[i,] != 0),drop=F]))
wmu1 <- colSums(dvL$d2mu[,,i]*pD)*c(ghQ$weights)
hmu1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zmu[z,]) * wmu1[z]),
            dim = c(ncol(Zmu),ncol(Zmu),length(ghQ$weights)))
hess1 <- hmu1 <- rowSums(hmu1, dims = 2)
# For second Hessian (eq. A3, row 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wmu2 <- colSums(dvL$d1mu[,,i]*dvL$d1mu[,,i]*pD)
hmu2 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zmu[z,]) * 
                  c(ghQ$weights)[z] * wmu2[z]),
            dim = c(ncol(Zmu),ncol(Zmu),length(ghQ$weights)))
hess2 <- hmu2 <- rowSums(hmu2, dims = 2)
# For third Hessian (eq. A3, row 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wmu3 <- sweep(dvL$d1mu[,,i]*pD,MARGIN = 2,c(ghQ$weights),"*")
Zhe3 <- array(sapply(1:nrow(pD), function(m) Zmu*wmu3[m,]),
              dim = c(nrow(ghQ$points),ncol(Zmu),nrow(pD)))
if("sigma" %in%  pFun(fam[i])){
 # For first Hessian (eq. A3, row 1)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Zsg <- unname(as.matrix(Z$sigma[,which(b$sigma[i,] != 0),drop=F]))
 wsg1 <- colSums(dvL$d2sg[,,i]*pD)*c(ghQ$weights)
 wms1 <- colSums(dvL$dcms[,,i]*pD)*c(ghQ$weights)
 hsg1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zsg[z,]) * wsg1[z]),
            dim = c(ncol(Zsg),ncol(Zsg),length(ghQ$weights)))
 hms1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zmu[z,],Zsg[z,]) * wms1[z]),
              dim = c(ncol(Zmu),ncol(Zsg),length(ghQ$weights)))
 hsg1 <- rowSums(hsg1, dims = 2)
 hms1 <- rowSums(hms1, dims = 2)
 hess1 <- cbind(rbind(hmu1,t(hms1)),rbind(hms1,hsg1))
 # For second Hessian (eq. A3, row 2)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 wsg2 <- colSums(dvL$d1sg[,,i]*dvL$d1sg[,,i]*pD)
 wms2 <- colSums(dvL$d1mu[,,i]*dvL$d1sg[,,i]*pD)
 hsg2 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zsg[z,]) * 
                     c(ghQ$weights)[z] * wsg2[z]),
               dim = c(ncol(Zsg),ncol(Zsg),length(ghQ$weights)))
 hms2 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zmu[z,],Zsg[z,]) *
                     c(ghQ$weights)[z] * wms2[z]),
               dim = c(ncol(Zmu),ncol(Zsg),length(ghQ$weights)))
 hsg2 <- rowSums(hsg2, dims = 2)
 hms2 <- rowSums(hms2, dims = 2)
 hess2 <- cbind(rbind(hmu2,t(hms2)),rbind(hms2,hsg2))
 # For third Hessian (eq. A3, row 3)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 wsg3 <- sweep(dvL$d1sg[,,i]*pD,MARGIN = 2,c(ghQ$weights),"*")
 the3 <- array(sapply(1:nrow(pD), function(m) Zsg*wsg3[m,]),
              dim = c(nrow(ghQ$points),ncol(Zsg),nrow(pD)))
 Zhe3 <- abind::abind(Zhe3,the3,along = 2)
} 
if("tau" %in%  pFun(fam[i])){
 # For first Hessian (eq. A3, row 1)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Zta <- unname(as.matrix(Z$tau[,which(b$tau[i,] != 0),drop=F]))
 wta1 <- colSums(dvL$d2ta[,,i]*pD)*c(ghQ$weights)
 wmt1 <- colSums(dvL$dcmt[,,i]*pD)*c(ghQ$weights)
 wst1 <- colSums(dvL$dcst[,,i]*pD)*c(ghQ$weights)
 hta1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zta[z,]) * wta1[z]),
               dim = c(ncol(Zta),ncol(Zta),length(ghQ$weights)))  
 hmt1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zmu[z,],Zta[z,]) * wmt1[z]),
               dim = c(ncol(Zmu),ncol(Zta),length(ghQ$weights))) 
 hst1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zsg[z,],Zta[z,]) * wst1[z]),
               dim = c(ncol(Zsg),ncol(Zta),length(ghQ$weights)))   
 hta1 <- rowSums(hta1, dims = 2)
 hmt1 <- rowSums(hmt1, dims = 2)
 hst1 <- rowSums(hst1, dims = 2)
 hess1 <- cbind(rbind(hess1,cbind(t(hmt1),t(hst1))),rbind(hmt1,hst1,hta1))
 # For second Hessian (eq. A3, row 2)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 wta2 <- colSums(dvL$d1ta[,,i]*dvL$d1ta[,,i]*pD)
 wmt2 <- colSums(dvL$d1mu[,,i]*dvL$d1ta[,,i]*pD)
 wst2 <- colSums(dvL$d1sg[,,i]*dvL$d1ta[,,i]*pD)
 hta2 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zta[z,]) * 
                     c(ghQ$weights)[z] * wta2[z]),
               dim = c(ncol(Zta),ncol(Zta),length(ghQ$weights)))
 hmt2 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zmu[z,],Zta[z,]) *
                     c(ghQ$weights)[z] * wmt2[z]),
               dim = c(ncol(Zmu),ncol(Zta),length(ghQ$weights)))
 hst2 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zsg[z,],Zta[z,]) *
                     c(ghQ$weights)[z] * wst2[z]),
               dim = c(ncol(Zsg),ncol(Zta),length(ghQ$weights)))
 hta2 <- rowSums(hta2, dims = 2)
 hmt2 <- rowSums(hmt2, dims = 2)
 hst2 <- rowSums(hst2, dims = 2)
 hess2 <- cbind(rbind(hess2,cbind(t(hmt2),t(hst2))),rbind(hmt2,hst2,hta2))
 # For third Hessian (eq. A3, row 3)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 wta3 <- sweep(dvL$d1ta[,,i]*pD,MARGIN = 2,c(ghQ$weights),"*")
 the3 <- array(sapply(1:nrow(pD), function(m) Zta*wta3[m,]),
              dim = c(nrow(ghQ$points),ncol(Zta),nrow(pD)))
 Zhe3 <- abind::abind(Zhe3,the3,along = 2)
}
if("nu" %in%  pFun(fam[i])){
 # For first Hessian (eq. A3, row 1)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 wnu1 <- colSums(dvL$d2nu[,,i]*pD)*c(ghQ$weights)
 wmn1 <- colSums(dvL$dcmn[,,i]*pD)*c(ghQ$weights)
 wsn1 <- colSums(dvL$dcsn[,,i]*pD)*c(ghQ$weights)
 wtn1 <- colSums(dvL$dctn[,,i]*pD)*c(ghQ$weights)
 Znu <- unname(as.matrix(Z$nu[,which(b$nu[i,] != 0),drop=F]))
 hnu1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Znu[z,]) * wnu1[z]),
               dim = c(ncol(Znu),ncol(Znu),length(ghQ$weights)))
 hmn1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zmu[z,],Znu[z,]) * wmn1[z]),
               dim = c(ncol(Zmu),ncol(Znu),length(ghQ$weights)))
 hsn1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zsg[z,],Znu[z,]) * wsn1[z]),
               dim = c(ncol(Zsg),ncol(Znu),length(ghQ$weights)))
 htn1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zta[z,],Znu[z,]) * wtn1[z]),
               dim = c(ncol(Zta),ncol(Znu),length(ghQ$weights)))
 hnu1 <- rowSums(hnu1, dims = 2)
 hmn1 <- rowSums(hmn1, dims = 2)
 hsn1 <- rowSums(hsn1, dims = 2)
 htn1 <- rowSums(htn1, dims = 2)
 hess1 <- cbind(rbind(hess1,cbind(t(hmn1),t(hsn1),t(htn1))),rbind(hmn1,hsn1,htn1,hnu1))
 # For second Hessian (eq. A3, row 2)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 wnu2 <- colSums(dvL$d1nu[,,i]*dvL$d1nu[,,i]*pD)
 wmn2 <- colSums(dvL$d1mu[,,i]*dvL$d1nu[,,i]*pD)
 wsn2 <- colSums(dvL$d1sg[,,i]*dvL$d1nu[,,i]*pD)
 wtn2 <- colSums(dvL$d1ta[,,i]*dvL$d1nu[,,i]*pD)
 hnu2 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Znu[z,]) * 
                     c(ghQ$weights)[z] * wnu2[z]),
               dim = c(ncol(Znu),ncol(Znu),length(ghQ$weights)))
 hmn2 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zmu[z,],Znu[z,]) *
                     c(ghQ$weights)[z] * wmn2[z]),
               dim = c(ncol(Zmu),ncol(Znu),length(ghQ$weights)))
 hsn2 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zsg[z,],Znu[z,]) *
                     c(ghQ$weights)[z] * wsn2[z]),
               dim = c(ncol(Zsg),ncol(Znu),length(ghQ$weights)))
 htn2 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zta[z,],Znu[z,]) *
                     c(ghQ$weights)[z] * wtn2[z]),
               dim = c(ncol(Zta),ncol(Znu),length(ghQ$weights)))
 hnu2 <- rowSums(hnu2, dims = 2)
 hmn2 <- rowSums(hmn2, dims = 2)
 hsn2 <- rowSums(hsn2, dims = 2)
 htn2 <- rowSums(htn2, dims = 2)
 hess2 <- cbind(rbind(hess2,cbind(t(hmn2),t(hsn2),t(htn2))),rbind(hmt2,hst2,htn2,hta2))
 # For third Hessian (eq. A3, row 3)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 wnu3 <- sweep(dvL$d1nu[,,i]*pD,MARGIN = 2,c(ghQ$weights),"*")
 the3 <- array(sapply(1:nrow(pD), function(m) Znu*wnu3[m,]),
              dim = c(nrow(ghQ$points),ncol(Znu),nrow(pD)))
 Zhe3 <- abind::abind(Zhe3,the3,along = 2)
}
Zhe3 <- t(sapply(1:nrow(pD), function(m) colSums(Zhe3[,,m])))
hess3 <- matrix(rowSums(sapply(1:nrow(pD), function(m) tcrossprod(Zhe3[m,]))),
                ncol(Zhe3))
hess <- hess1 + hess2 - hess3
} else {
iN0 <- which(lb2mb(b)[i,] != 0)
jN0 <- which(lb2mb(b)[j,] != 0)
# For first Hessian (eq. A3, row 1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Zhe2i <- Zmui <- unname(as.matrix(Z$mu[,which(b$mu[i,] != 0)]))
Zhe2j <- Zmuj <- unname(as.matrix(Z$mu[,which(b$mu[j,] != 0)]))
# For second Hessian (eq. A3, row 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wmu2i <- dvL$d1mu[,,i]; wl2i <- "wmu2i"
wmu2j <- dvL$d1mu[,,j]; wl2j <- "wmu2j"
imui <- 1:ncol(Zmui); idxi <- "imui"
imuj <- 1:ncol(Zmuj); idxj <- "imuj"
# For third Hessian (eq. A3, row 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wmu3i <- sweep(dvL$d1mu[,,i]*pD,MARGIN = 2,c(ghQ$weights),"*")
Zhe3i <- array(sapply(1:nrow(pD), function(m) Zmui*wmu3i[m,]),
               dim = c(nrow(ghQ$points),ncol(Zmui),nrow(pD)))
wmu3j <- sweep(dvL$d1mu[,,j]*pD,MARGIN = 2,c(ghQ$weights),"*")
Zhe3j <- array(sapply(1:nrow(pD), function(m) Zmuj*wmu3j[m,]),
               dim = c(nrow(ghQ$points),ncol(Zmuj),nrow(pD)))
for(k in c("i","j")){
if("sigma" %in% pFun(fam[get(k)])){
 assign(paste0("Zsg",k), unname(as.matrix(Z$sigma[,which(b$sigma[get(k),] != 0)])))
 assign(paste0("wsg2",k), dvL$d1sg[,,get(k)] )
 assign(paste0("wl2",k), append(get(paste0("wl2",k)), paste0("wsg2",k)))
 assign(paste0("Zhe2",k), cbind(get(paste0("Zhe2",k)), get(paste0("Zsg",k))))
 assign(paste0("isg",k), (ncol(get(paste0("Zmu",k)))+1):ncol(get(paste0("Zhe2",k))))
 assign(paste0("idx",k), append(get(paste0("idx",k)), paste0("isg",k)))
 # For third Hessian (eq. A3, row 3)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 assign(paste0("wsg3",k), sweep(dvL$d1sg[,,get(k)]*pD,MARGIN = 2,c(ghQ$weights),"*"))
 assign(paste0("the3",k), array(sapply(1:nrow(pD), function(m) get(paste0("Zsg",k)) * get(paste0("wsg3",k))[m,]),
              dim = c(nrow(ghQ$points),ncol(get(paste0("Zsg",k))),nrow(pD))))
 assign(paste0("Zhe3",k), abind::abind(get(paste0("Zhe3",k)),get(paste0("the3",k)),along = 2))
}
if("tau" %in% pFun(fam[get(k)])){
 assign(paste0("Zta",k), unname(as.matrix(Z$tau[,which(b$tau[get(k),] != 0)])))
 assign(paste0("wta2",k), dvL$d1ta[,,get(k)] )
 assign(paste0("wl2",k), append(get(paste0("wl2",k)), paste0("wta2",k)))
 assign(paste0("Zhe2",k), cbind(get(paste0("Zhe2",k)), get(paste0("Zta",k))))
 assign(paste0("ita",k), (ncol(get(paste0("Zsg",k)))+1):ncol(get(paste0("Zhe2",k))))
 assign(paste0("idx",k), append(get(paste0("idx",k)), paste0("ita",k)))
 # For third Hessian (eq. A3, row 3)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 assign(paste0("wta3",k), sweep(dvL$d1ta[,,get(k)]*pD,MARGIN = 2,c(ghQ$weights),"*"))
 assign(paste0("the3",k), array(sapply(1:nrow(pD), function(m) get(paste0("Zta",k)) * get(paste0("wta3",k))[m,]),
              dim = c(nrow(ghQ$points),ncol(get(paste0("Zta",k))),nrow(pD))))
 assign(paste0("Zhe3",k), abind::abind(get(paste0("Zhe3",k)),get(paste0("the3",k)),along = 2))
}
if("nu" %in% pFun(fam[get(k)])){
 assign(paste0("Znu",k), unname(as.matrix(Z$nu[,which(b$nu[get(k),] != 0)])))
 assign(paste0("wnu2",k), dvL$d1nu[,,get(k)] )
 assign(paste0("wl2",k), append(get(paste0("wl2",k)), paste0("wnu2",k)))
 assign(paste0("Zhe2",k), cbind(get(paste0("Zhe2",k)), get(paste0("Znu",k))))
 assign(paste0("inu",k), (ncol(get(paste0("Zta",k)))+1):ncol(get(paste0("Zhe2",k))))
 assign(paste0("idx",k), append(get(paste0("idx",k)), paste0("inu",k)))
 # For third Hessian (eq. A3, row 3)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 assign(paste0("wnu3",k), sweep(dvL$d1nu[,,get(k)]*pD,MARGIN = 2,c(ghQ$weights),"*"))
 assign(paste0("the3",k), array(sapply(1:nrow(pD), function(m) get(paste0("Znu",k)) * get(paste0("wnu3",k))[m,]),
              dim = c(nrow(ghQ$points),ncol(get(paste0("Znu",k))),nrow(pD))))
 assign(paste0("Zhe3",k), abind::abind(get(paste0("Zhe3",k)),get(paste0("the3",k)),along = 2))
} }
# For second Hessian (eq. A3, row 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idM <- unname(as.matrix(expand.grid(1:length(wl2i), 1:length(wl2j))))
wij2 <- array(0,dim = c(ncol(Zhe2i),ncol(Zhe2j),nrow(ghQ$points)))
for(z in 1:nrow(ghQ$points)){ for(zz in 1:nrow(idM)){
 tidi <- get(idxi[idM[zz,][1]])
 tidj <- get(idxj[idM[zz,][2]])
 wij2[tidi,tidj,z] <- c(colSums(get(wl2i[idM[zz,][1]])*get(wl2j[idM[zz,][2]])*pD)*c(ghQ$weights))[z]
}}; rm(idM)
hess2 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zhe2i[z,],Zhe2j[z,])),
               dim = c(ncol(Zhe2i),ncol(Zhe2j),nrow(ghQ$points)))
hess2 <- t(rowSums(hess2*wij2,dims = 2))
# For third Hessian (eq. A3, row 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Zhe3i <- t(sapply(1:nrow(pD), function(m) colSums(Zhe3i[,,m])))
Zhe3j <- t(sapply(1:nrow(pD), function(m) colSums(Zhe3j[,,m])))
hess3 <- t(matrix(rowSums(sapply(1:nrow(pD), function(m) tcrossprod(Zhe3i[m,],Zhe3j[m,]))),
                nrow = ncol(Zhe3i), ncol = ncol(Zhe3j)))
hess1 <- matrix(0,nrow = nrow(hess2),ncol = ncol(hess2))
hess <- hess1 + hess2 - hess3
}
return(list("hessian" = hess, "posi" = iN0, "posj" = jN0, "fisher" = hess1))
}

upB <- function(bold,ghQ,fam,dvL,pD,full.hess = T,information = c("Hessian","Fisher")){

# Goal: To update factor loadings
# Input : bold (loadings), ghQ (GHQ object), fam (distributions),
#         dvL (list of derivatives), pD (posterior density)
#         full.hess (whether to update item by item or using full hessian)
# Output: bnew (new loadings)
# Testing: bold = borg; ghQ = ghQ; fam = fam; dvL = dvY(Y,ghQ,bold,fam)
#          pD = exp(rowSums(dY(Y,ghQ,bold,fam),dim = 2)) /
#          c(exp(rowSums(dY(Y,ghQ,bold,fam),dim = 2))%*%ghQ$weights)  
#          full.hess = T;
          
info = match.arg(information, c("Hessian","Fisher"))
if(info == "Fisher" && full.hess == T) full.hess <- F
b = lb2mb(bold)
if(full.hess == T){
# Update using Full Hessian
# ~~~~~~~~~~~~~~~~~~~~~~~~~
cb <- lb2cb(bold)
Sm <- NULL
i1 <- i2 <- rep(0,length(fam))
Hm <- matrix(0,sum(b != 0),sum(b != 0))
for(i in 1:length(fam)){
 A <- sco(i,ghQ,bold,fam,dvL,pD); i1[i] <- length(A$score)
 Sm <- c(Sm,A$score); i2[i]<- length(Sm) 
 Hr <- NULL
 for(j in i:length(fam)){
  B <- hess(i,j,ghQ,bold,fam,dvL,pD)
  Hr <- rbind(Hr,B$hessian)
 }
Hm[(i2[i]-i1[i]+1):nrow(Hm),(i2[i]-i1[i]+1):i2[i]] <- Hr
Hm[(i2[i]-i1[i]+1):i2[i],(i2[i]-i1[i]+1):ncol(Hm)] <- t(Hr) 
}
# Check PD
# ~~~~~~~~
# Hm <- m2pdm(-Hm)
# cb[cb != 0] <- cb[cb != 0] + c(Hm$inv.mat%*%matrix(Sm))
cb[cb != 0] <- cb[cb != 0] - c(solve(Hm)%*%matrix(Sm))
# cb[cb != 0] <- cb[cb != 0] + 0.0005*c(matrix(Sm))
bnew <- cb2lb(cb,bold)
} else {
# Update by item
# ~~~~~~~~~~~~~~
Hm <- matrix()
Sm <- NULL
for(i in 1:length(fam)){
 A <- sco(i,ghQ,bold,fam,dvL,pD)
 B <- hess(i,i,ghQ,bold,fam,dvL,pD)
 # Check PD
 # ~~~~~~~~
 # B <- m2pdm(-B$hessian)
 # b[i, A$pos] <- c(b[i,A$pos]) + c(B$inv.mat%*%(A$score))
 if(info == "Hessian") b[i, A$pos] <- c(b[i,A$pos]) - c(solve(B$hess)%*%(A$score))
 if(info == "Fisher") b[i, A$pos] <- c(b[i,A$pos]) - c(solve(B$fisher)%*%(A$score))
 # b[i, A$pos] <- c(b[i,A$pos]) + 0.0005*c(A$score)
 Sm <- c(Sm,A$score)
 Hm <- Matrix::bdiag(Hm,B$hess)
 } 
bnew <- mb2lb(b,bold)
Hm <- Hm[-1,-1] }
return(list(b = bnew, gradient = Sm, hessian = as.matrix(Hm)))
}

upB.pen <- function(bold,ghQ,fam,dvL,pD,full.hess = T,pen.idx,information = c("Hessian","Fisher"),
                    pml.control = list(type = "lasso", lambda = 1, w.alasso = NULL,a = NULL)){
                
# Goal: To update factor loadings (when Penalised EM)
# Input : bold (loadings), ghQ (GHQ object), fam (distributions),
#         dvL (list of derivatives), pD (posterior density)
#         full.hess (whether to update item by item or using full hessian)
# Output: bnew (new loadings)
# Testing: bold = borg; ghQ = ghQ; fam = fam; dvL = dvY(Y,ghQ,bold,fam,control$info)
#          pD = exp(rowSums(dY(Y,ghQ,bold,fam),dim = 2)) /
#          c(exp(rowSums(dY(Y,ghQ,bold,fam),dim = 2))%*%ghQ$weights)  
#          full.hess = T;
#          pml.control = list(type = "lasso", lambda = 1, w.alasso = NULL, gamma = 1, a = 3.7)

info = match.arg(information, c("Hessian","Fisher"))
if(info == "Fisher" && full.hess == T) full.hess <- F
b = lb2mb(bold)
if(full.hess == T){
# Update using Full Hessian
# ~~~~~~~~~~~~~~~~~~~~~~~~~
cb <- lb2cb(bold)
Sm <- NULL
i1 <- i2 <- rep(0,length(fam))
Hm <- matrix(0,sum(b != 0),sum(b != 0))
for(i in 1:length(fam)){
 A <- sco(i,ghQ,bold,fam,dvL,pD) ; i1[i] <- length(A$score)
 Sm <- c(Sm,A$score); i2[i]<- length(Sm) 
 Hr <- NULL
 for(j in i:length(fam)){
  B <- hess(i,j,ghQ,bold,fam,dvL,pD)
  Hr <- rbind(Hr,B$hessian)
 }
Hm[(i2[i]-i1[i]+1):nrow(Hm),(i2[i]-i1[i]+1):i2[i]] <- Hr
Hm[(i2[i]-i1[i]+1):i2[i],(i2[i]-i1[i]+1):ncol(Hm)] <- t(Hr) 
}
if(is.list(pml.control$w.alasso)) pml.control$w.alasso <- lb2mb(pml.control$w.alasso)[pen.idx]
if(length(cb[c(t(pen.idx))]) != 0){
pS <- nrow(pD)*penM(cb[c(t(pen.idx))],type = pml.control$type, lambda = pml.control$lambda, #  == (b != 0)
                    w.alasso = pml.control$w.alasso, a = pml.control$a)
pS. <- rep(0,length(cb)); pS.[c(t(pen.idx))] <- diag(pS)
pS <- diag(pS.); pS <- pS[cb != 0,cb != 0]; } else { pS <- diag(rep(0,sum(cb != 0))) }
Sm <- Sm - pS%*%cb[cb != 0]
Hm <- Hm - pS
# Check PD
# ~~~~~~~~
# Hm <- m2pdm(-Hm)
# cb[cb != 0] <- cb[cb != 0] + c(Hm$inv.mat%*%matrix(Sm))
cb[cb != 0] <- cb[cb != 0] - c(solve(Hm)%*%matrix(Sm))
cb[abs(cb) < 1*sqrt(.Machine$double.eps)] <- 0
bnew <- cb2lb(cb,bold)
} else {
# Update by item
# ~~~~~~~~~~~~~~
Hm <- matrix()
Sm <- NULL
for(i in 1:length(fam)){
if(is.list(pml.control$w.alasso)) pml.control$w.alasso <- lb2mb(pml.control$w.alasso)[i,pen.idx[i,]]
 if(length(c(b[i,pen.idx[i,]])) != 0){
 pS <- nrow(pD)*penM(c(b[i,pen.idx[i,]]),type = pml.control$type, lambda = pml.control$lambda, w.alasso = pml.control$w.alasso, a = pml.control$a);
 if(length(pS) != 1) pS <- diag(pS)
 pS. <- diag(rep(0,length(c(b[i,])))); diag(pS.)[c(t(pen.idx[i,]))] <- pS
 pS <- pS.; pS <- pS[b[i,] != 0, b[i,] != 0]; } else { pS <- diag(rep(0,sum(c(b[i,] != 0)))) }
 A <- sco(i,ghQ,bold,fam,dvL,pD)
 As <- A$score - pS%*%b[i,b[i,] != 0]
 if(info == "Hessian") B <- hess(i,i,ghQ,bold,fam,dvL,pD)$hessian - pS
 if(info == "Fisher") B <- hess(i,i,ghQ,bold,fam,dvL,pD)$fisher - pS
 # Check PD
 # ~~~~~~~~
 # B <- m2pdm(-B)
 # b[i, A$pos] <- c(b[i,A$pos]) + c(B$inv.mat%*%(As))
 b[i, A$pos] <- c(b[i,A$pos]) - c(solve(B)%*%(As))
 # if(info == "Hessian") b[i, A$pos] <- c(b[i,A$pos]) - c(solve(B$hess)%*%(As))
 # if(info == "Fisher") b[i, A$pos] <- c(b[i,A$pos]) - c(solve(B$fisher)%*%(As))
 b[i, abs(b[i,]) < 1*sqrt(.Machine$double.eps)] <- 0
 Sm <- c(Sm,As)
 Hm <- Matrix::bdiag(Hm,B)
 }
bnew <- mb2lb(b,bold)
Hm <- Hm[-1,-1]}
return(list(b = bnew, hessian = as.matrix(Hm)))
}

sche <- function(ghQ,b,fam,dvL,pD,info){

# Goal: To compute Score vector & full Hessian matrix for all items
# Input : ghQ (GHQ object), b (loadings matrix), fam (distributions),
#         dvL (list of derivatives),
#         pD (posterior density, to be compute with joint with weights)
# Output: list of gradient vector and Hessian matrix
# Testing: ghQ = gr; b = borg; fam = fam; dvL = dvY(Y,ghQ,b,fam)
#          pD = exp(rowSums(dY(Y,ghQ,b,fam),dim = 2)) /
#          c(exp(rowSums(dY(Y,ghQ,b,fam),dim = 2))%*%ghQ$weights)

info.F <- match.arg(info,c("Fisher","Hessian")) == "Fisher"
cb <- lb2cb(b)
Sm <- NULL
i1 <- i2 <- rep(0,length(fam))
Hm <- matrix(0,sum(cb != 0),sum(cb != 0))
for(i in 1:length(fam)){
 A <- sco(i,ghQ,b,fam,dvL,pD); i1[i] <- length(A$score)
 Sm <- c(Sm,A$score); i2[i]<- length(Sm) 
 Hr <- NULL
 for(j in i:length(fam)){
  B <- hess(i,j,ghQ,b,fam,dvL,pD)
  if(info.F) Hr <- rbind(Hr,B$fisher) else
  Hr <- rbind(Hr,B$hessian)
 }
Hm[(i2[i]-i1[i]+1):nrow(Hm),(i2[i]-i1[i]+1):i2[i]] <- Hr
Hm[(i2[i]-i1[i]+1):i2[i],(i2[i]-i1[i]+1):ncol(Hm)] <- t(Hr) 
}
return(list(gradient = Sm, hessian = Hm))
}

penM <- function(params, type = "lasso", lambda = 1, w.alasso = NULL,a = 3.7){ 

# Goal: To produce a Penalty matrix of parameters (params)
# Input : model parameters (params)
# Output: Penalty matrix
# Testing: params = lb2cb(borg); lambda = 1; gamma = 1; a = 3.7;
# 

eps = 1e-8 # sqrt(.Machine$double.eps) # protective tolerance level
if(type == "ridge"){
 A1 <-  lambda*rep(1, length(params))
}
if(type == "lasso"){
 A1 <- lambda*1/sqrt(params^2 + eps) 
}
if(type == "alasso"){
 if(is.null(a)) a = 2
 if( is.null(w.alasso) ) w.alasso <- 1
 w.al <- 1/abs(w.alasso)^a
 A1 <- lambda*w.al/sqrt(params^2 + eps)
}
if(type == "scad"){
 if(is.null(a)) a = 3.7
 theta <- abs(params) 
 f1 <- sapply(theta, function(theta) { max(a*lambda - theta, 0)/((a-1)*lambda + eps) })
 f.d <- ((theta <= lambda) + f1 * (theta > lambda))
 A1 <- lambda* f.d / ( sqrt(params^2 + eps) )
}
if(type == "mcp"){
 if(is.null(a)) a = 2.5
 theta <- abs(params) 
 f.d <- (lambda-theta/a)*(theta < lambda*a)
 A1 <- f.d / ( sqrt(params^2 + eps) )
}
if(length(A1) == 1) S <- matrix(A1) else S <- diag(A1)
return(S)
}
