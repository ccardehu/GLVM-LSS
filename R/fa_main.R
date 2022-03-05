
##############################################
## SPLVM: Matrices functions for estimation ##
##############################################

dY <- function(Y,ghQ,b,fam){

# Goal: To compute (conditional) density of items in (nxp) matrix Y
# Input : Y (item matrix), ghQ (GHQ object), b (loadings matrix), fam (distributions)
# Output: Array n x qp x p with f(y_p|z_qp)
# Testing: Y = simR$Y; ghQ = ghQ; b = borg; fam = fam
# C++ alternative: DY

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
if(fam[i] == "beta"){
 mu = c(probs(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
 sigma = c(probs(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
 alpha = mu*(1-sigma^2)/sigma^2
 beta = (1-mu)*(1-sigma^2)/sigma^2
 fyz[,z,i] <- dbeta(badj(Y[,i]),shape1 = alpha, shape2 = beta, log = T)
}
if(fam[i] == "gumbel"){
 mu = c(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
 sigma = c(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
 fyz[,z,i] <- gamlss.dist::dGU(Y[,i], mu, sigma, log = T)
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
  dvY$d2sg[,z,i] = (-3*(Y[,i]-mu)^2/sigma^4 + 1/sigma^2) * sigma^2 + dvY$d1sg[,z,i]
 if(info.F) dvY$dcms[,z,i] = rep(0, length(Y[,i])) else
  dvY$dcms[,z,i] = -2*(Y[,i]-mu)/sigma^3 * 1 * sigma
}
if(fam[i] == "lognormal"){ # structure(gamlss.dist::LOGNO)
 mu = c(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
 sigma = c(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
 dvY$d1mu[,z,i] = (log(Y[,i]) - mu)/(sigma^2) * 1
 dvY$d2mu[,z,i] = rep(-(1/sigma^2) * 1, length(Y[,i])) + 0
 dvY$d1sg[,z,i] = (((log(Y[,i]) - mu)^2 - sigma^2)/(sigma^3)) * sigma
 if(info.F) dvY$d2sg[,z,i] = rep(-(2/(sigma^2)) * sigma^2, length(Y[,i])) else 
  dvY$d2sg[,z,i] = (-3*(log(Y[,i])-mu)^2/sigma^4 + 1/sigma^2) * sigma^2 + dvY$d1sg[,z,i]
 if(info.F) dvY$dcms[,z,i] = rep(0, length(Y[,i])) else
  dvY$dcms[,z,i] = -2*(log(Y[,i])-mu)/sigma^3 * 1 * sigma
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
 if(info.F) dvY$d2mu[,z,i] = (1-u)*(-1/mu) * mu^2 else 
  dvY$d2mu[,z,i] = (1-u)*(-Y[,i]/mu^2) * mu^2 + dvY$d1mu[,z,i] # expected value is -1/mu; original is (-Y/mu^2)
 dvY$d1sg[,z,i] = (u-sigma)/(sigma*(1-sigma)) * (sigma*(1-sigma))
 if(info.F) dvY$d2sg[,z,i] = 1/(sigma*(sigma-1)) * (sigma*(1-sigma))^2 else  # CHECK THIS
  dvY$d2sg[,z,i] = -(sigma^2 - 2*u*sigma + u)/((sigma-1)*sigma)^2 * (sigma*(1-sigma))^2 + (1-2*sigma)*dvY$d1sg[,z,i] # -(u/(sigma^2) + (1-u)/((1-sigma)^2)) * (sigma*(1-sigma))^2
 dvY$dcms[,z,i] = rep(0, length(Y[,i]))
}
if(fam[i] == "beta"){ # structure(gamlss.dist::BEo)
 mu = c(probs(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
 sigma = c(probs(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
 a = mu*(1 - sigma^2)/(sigma^2)
 beta = (1-mu)*(1 - sigma^2)/(sigma^2)
 y. = badj(Y[,i])
 dvY$d1mu[,z,i] = (((1 - sigma^2)/(sigma^2))*(-digamma(a) + digamma(beta) + log(y.) - log(1 - y.)))*(mu*(1-mu))
 if(info.F) dvY$d2mu[,z,i] = rep(-(((1 - sigma^2)^2)/(sigma^4)) * (trigamma(a) + trigamma(beta))*(mu*(1-mu))^2, length(y.)) else 
  dvY$d2mu[,z,i] = rep(-(((1 - sigma^2)^2)/(sigma^4)) * (trigamma(a) + trigamma(beta))*(mu*(1-mu))^2, length(y.)) + (1-2*mu)*dvY$d1mu[,z,i]
 dvY$d1sg[,z,i] = -(2/(sigma^3)) * (mu * (-digamma(a) + digamma(a + beta) + log(y.)) + (1 - mu) * (-digamma(beta) + digamma(a + beta) + log(1 - y.)))*(sigma*(1-sigma))
 if(info.F) dvY$d2sg[,z,i] = rep(-(4/(sigma^6)) * ((mu^2) * trigamma(a) + ((1 - mu)^2) * trigamma(beta) - trigamma(a + beta)) * (sigma*(1-sigma))^2, length(y.)) else
  dvY$d2sg[,z,i] = rep(-(4/(sigma^6)) * ((mu^2) * trigamma(a) + ((1 - mu)^2) * trigamma(beta) - trigamma(a + beta)) * (sigma*(1-sigma))^2, length(y.)) + (1-2*sigma)*dvY$d1sg[,z,i]
 dvY$dcms[,z,i] = rep((2 * (1 - sigma^2)/(sigma^5)) * (mu * trigamma(a) - (1 - mu) * trigamma(beta))*(mu*(1-mu))*(sigma*(1-sigma)), length(y.))
}
if(fam[i] == "gumbel"){ # structure(gamlss.dist::GU)
 mu = c(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
 sigma = c(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
 dvY$d1mu[,z,i] = (exp((Y[,i] - mu)/sigma) - 1)/sigma * 1
 if(!info.F) dvY$d2mu[,z,i] = exp((Y[,i] - mu)/sigma)/sigma^2 + dvY$d1mu[,z,i] else
  dvY$d2mu[,z,i] = rep(-(1/sigma^2) * 1, length(Y[,i]))
 dvY$d1sg[,z,i] = (((Y[,i] - mu)/sigma^2) * (exp((Y[,i] - mu)/sigma) - 1) - (1/sigma)) * sigma
 if(!info.F) dvY$d2sg[,z,i] = ((-(Y[,i] - mu)^2/sigma^4) * exp((Y[,i] - mu)/sigma) * (-2*(Y[,i] - mu)/sigma^3) * (exp((Y[,i] - mu)/sigma) - 1) + (1/sigma^2)) * sigma^2 + dvY$d1sg[,z,i] else
  dvY$d2sg[,z,i] = rep(-1.82368/sigma^2, length(Y[,i]))* sigma^2
 if(!info.F) dvY$dcms[,z,i] =  (-(sigma + mu + Y[,i])*exp((Y[,i] - mu)/sigma)-sigma)/sigma^3 else
  dvY$dcms[,z,i] = rep(-0.422784/sigma^2, length(Y[,i])) * 1 * sigma # rep(0, length(Y)) # the expected cross derivative mu and sigma
} } }
for(i in names(dvY)){ if(all(is.na(dvY[[i]]))) dvY[[i]] <- NULL } 
return(dvY)
}

sco <- function(i,ghQ,b,rs,fam,dvL,pD){
  
# Goal: To compute Score vector for item i (eq. A1)
# Input : i (item), ghQ (GHQ object), b (loadings matrix), rs (restriction matrix),
#         fam (distributions), dvL (list of derivatives),
#         pD (posterior density, to be compute with joint with weights)
# Output: vector of length (total # of parameters for item i)
# Testing: i = 1; Y = simR$Y; ghQ = ghQ; b = borg; fam = fam; 
#          dvL = dvY(Y,ghQ,b,fam); rs = loadmt;
#          pD = exp(rowSums(dY(Y,ghQ,b,fam),dim = 2)) /
#          c(exp(rowSums(dY(Y,ghQ,b,fam),dim = 2))%*%ghQ$weights)

iN0 <- which(lb2mb(rs)[i,])
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
ps <- as.matrix(colSums(unname(ps[,iN0,drop = F])))
return(list("score" = ps, "pos" = iN0))
}

hess <- function(i,j,ghQ,b,rs,fam,dvL,pD,information = "Fisher"){
  
# Goal: To compute Hessian matrix for items i & j (eq. A3)
# Input : i & j (items), ghQ (GHQ object), b (loadings matrix),
#         fam (distributions), dvL (list of derivatives),
#         pD (part of posterior density, to be jointly computed with weights)
# Output: vector of length (total # of parameters for item i)
# Testing: i = 1; Y = simR$Y; ghQ = ghQ; b = borg; fam = fam; rs = loadmt;
#          dvL = dvY(Y,ghQ,b,fam,control$information)
#          pD = exp(rowSums(dY(Y,ghQ,b,fam),dim = 2)) /
#          c(exp(rowSums(dY(Y,ghQ,b,fam),dim = 2))%*%ghQ$weights)

info.F <- match.arg(information,c("Fisher","Hessian")) == "Fisher"
Z <- ghQ$out
if(!info.F){ # If Hessian
if(i == j){
iN0 <- jN0 <- which(lb2mb(rs)[i,])
# For first Hessian (eq. A3, row 1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Zmu <- unname(as.matrix(Z$mu[,rs$mu[i,],drop=F]))
wmu1 <- colSums(dvL$d2mu[,,i]*pD)*c(ghQ$weights)
hmu1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zmu[z,]) * wmu1[z]),
            dim = c(ncol(Zmu),ncol(Zmu),length(ghQ$weights)))
hess1 <- hmu1 <- rowSums(hmu1, dims = 2)
# For second Hessian (eq. A3, row 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wmu2 <- colSums(dvL$d1mu[,,i]*dvL$d1mu[,,i]*pD)*c(ghQ$weights)
hmu2 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zmu[z,]) * wmu2[z]),
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
 Zsg <- unname(as.matrix(Z$sigma[,rs$sigma[i,],drop=F]))
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
 Zta <- unname(as.matrix(Z$tau[,rs$tau[i,],drop=F]))
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
 Znu <- unname(as.matrix(Z$nu[,rs$nu[i,],drop=F]))
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
Zhe3 <- t(sapply(1:nrow(pD), function(m) colSums(Zhe3[,,m,drop = T])))
hess3 <- matrix(rowSums(sapply(1:nrow(pD), function(m) tcrossprod(Zhe3[m,]))),
                ncol(Zhe3))
hess <- hess1 + hess2 - hess3
} else {
iN0 <- which(lb2mb(rs)[i,])
jN0 <- which(lb2mb(rs)[j,])
# For first Hessian (eq. A3, row 1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Zhe2i <- Zmui <- unname(as.matrix(Z$mu[,rs$mu[i,]]))
Zhe2j <- Zmuj <- unname(as.matrix(Z$mu[,rs$mu[j,]]))
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
 assign(paste0("Zsg",k), unname(as.matrix(Z$sigma[,rs$sigma[get(k),]])))
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
 assign(paste0("Zta",k), unname(as.matrix(Z$tau[,rs$tau[get(k),]])))
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
 assign(paste0("Znu",k), unname(as.matrix(Z$nu[,rs$nu[get(k),]])))
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
return(list("hessian" = hess, "posi" = iN0, "posj" = jN0))} else { # If Fisher Information
if(i == j){
iN0 <- jN0 <- which(lb2mb(rs)[i,])
# For first Hessian (eq. A3, row 1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Zmu <- unname(as.matrix(Z$mu[,rs$mu[i,],drop=F]))
wmu1 <- colSums(dvL$d2mu[,,i]*pD)*c(ghQ$weights)
hmu1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zmu[z,]) * wmu1[z]),
            dim = c(ncol(Zmu),ncol(Zmu),length(ghQ$weights)))
hess1 <- hmu1 <- rowSums(hmu1, dims = 2)
if("sigma" %in%  pFun(fam[i])){
 # For first Hessian (eq. A3, row 1)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Zsg <- unname(as.matrix(Z$sigma[,rs$sigma[i,],drop=F]))
 wsg1 <- colSums(dvL$d2sg[,,i]*pD)*c(ghQ$weights)
 wms1 <- colSums(dvL$dcms[,,i]*pD)*c(ghQ$weights)
 hsg1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zsg[z,]) * wsg1[z]),
            dim = c(ncol(Zsg),ncol(Zsg),length(ghQ$weights)))
 hms1 <- array(sapply(1:nrow(ghQ$points), function(z) tcrossprod(Zmu[z,],Zsg[z,]) * wms1[z]),
              dim = c(ncol(Zmu),ncol(Zsg),length(ghQ$weights)))
 hsg1 <- rowSums(hsg1, dims = 2)
 hms1 <- rowSums(hms1, dims = 2)
 hess1 <- cbind(rbind(hmu1,t(hms1)),rbind(hms1,hsg1))
} 
if("tau" %in%  pFun(fam[i])){
 # For first Hessian (eq. A3, row 1)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Zta <- unname(as.matrix(Z$tau[,rs$tau[i,],drop=F]))
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
}
if("nu" %in%  pFun(fam[i])){
 # For first Hessian (eq. A3, row 1)
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 wnu1 <- colSums(dvL$d2nu[,,i]*pD)*c(ghQ$weights)
 wmn1 <- colSums(dvL$dcmn[,,i]*pD)*c(ghQ$weights)
 wsn1 <- colSums(dvL$dcsn[,,i]*pD)*c(ghQ$weights)
 wtn1 <- colSums(dvL$dctn[,,i]*pD)*c(ghQ$weights)
 Znu <- unname(as.matrix(Z$nu[,rs$nu[i,],drop=F]))
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
}
hess <- hess1
} else {
iN0 <- which(lb2mb(rs)[i,])
jN0 <- which(lb2mb(rs)[j,])
hess1 <- matrix(0,nrow = length(iN0),ncol = length(jN0))
hess <- hess1
}
return(list("hessian" = hess, "posi" = iN0, "posj" = jN0)) }
}

sche <- function(ghQ,b,rs,fam,dvL,pD,info,full.hess = F){

# Goal: To compute Score vector & full Hessian matrix for all items
# Input : ghQ (GHQ object), b (loadings matrix), fam (distributions),
#         dvL (list of derivatives),
#         pD (posterior density, to be compute with joint with weights)
# Output: list of gradient vector and Hessian matrix
# Testing: ghQ = ghQ; b = bold; fam = fam; dvL = dvY(Y,ghQ,b,fam); rs = loadmt;
#          pD = exp(rowSums(dY(Y,ghQ,b,fam),dim = 2)) /
#          c(exp(rowSums(dY(Y,ghQ,b,fam),dim = 2))%*%ghQ$weights); info = "Fisher"

cb <- lb2cb(b)
Sm <- array(NA,dim=dim(lb2mb(b)))
posM <- matrix(1:length(cb),nrow = nrow(Sm),byrow = T)
Hm <- matrix(0,length(cb),length(cb))
Hm[c(t(!lb2mb(rs))),] <- Hm[,c(t(!lb2mb(rs)))] <- NA
for(i in 1:length(fam)){
 A <- sco(i,ghQ,b,rs,fam,dvL,pD);
 Sm[i,A$pos] <- c(A$score)
 if(full.hess){
 for(j in i:length(fam)){
  B <- hess(i,j,ghQ,b,rs,fam,dvL,pD,match.arg(info,c("Fisher","Hessian")))
  Hm[posM[i,B$posi],posM[j,B$posj]] <- t(B$hessian)
  Hm[posM[j,B$posj],posM[i,B$posi]] <- B$hessian
 }
 } else {
  B <- hess(i,i,ghQ,b,rs,fam,dvL,pD,match.arg(info,c("Fisher","Hessian")))
  Hm[posM[i,B$posi],posM[i,B$posi]] <- t(B$hessian)
  Hm[posM[i,B$posi],posM[i,B$posi]] <- B$hessian
}}
Sm <- c(t(Sm)[t(lb2mb(rs))])
Hm <- Hm[c(t(lb2mb(rs))),c(t(lb2mb(rs)))]
# if(!m2pdm(-Hm)$is.PDM){ Hm <- -m2pdm(-Hm)$mat }
return(list(gradient = Sm, hessian = Hm))
}

upB <- function(bold,shObj,loadmt){

# Goal: To update factor loadings
# Input : bold (loadings), shObj (score-hessian object), loadmt (restri)
# Output: bnew (new loadings)
# Testing: bold = borg; 
#          shObj = sche(ghQ,bold,loadmt,fam,dvL,
#          pD,control$information,control$full.hess)

cb <- lb2cb(bold)
bu <- cb[t(lb2mb(loadmt))]
iM <- tryCatch({solve(shObj$hess)}, error = function(e){return(-m2pdm(-shObj$hess)$inv.mat)})
bn <- bu - c(iM%*%matrix(shObj$gradient))
cb[t(lb2mb(loadmt))] <- bn
bnew <- cb2lb(cb,bold)
adsol <- m2pdm(-shObj$hess)$is.PDM
return(list(b = bnew, gradient = shObj$gradient, hessian = shObj$hessian, adsol = adsol))
}

upB.pen <- function(bold,shObj,pmObj,loadmt){
                
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

cb <- lb2cb(bold)
bu <- cb[t(lb2mb(loadmt))]
M <- shObj$hess - pmObj$pM[seq_along(bu),seq_along(bu)]
iM <- tryCatch({solve(M)}, error = function(e){return(-m2pdm(-M)$inv.mat)})
bn <- bu - c(iM%*%matrix(shObj$gradient - pmObj$pM[seq_along(bu),seq_along(bu)]%*%bu))
# bn[abs(bn) < 1e5*sqrt(.Machine$double.eps)] <- 0
cb[t(lb2mb(loadmt))] <- bn
bnew <- cb2lb(cb,bold)
adsol <- m2pdm(-shObj$hess)$is.PDM
return(list(b = bnew,  gradient = shObj$gradient, hessian = M, adsol = adsol))
}

upSa <- function(ghQ,pD){ # pz,ghQ,pD,respz
 # if(is.null(respz)){ corlst <- which(lower.tri(pz),arr.ind = T) } else corlst <- which(lower.tri(pz),arr.ind = T)[-respz,]
 # V <- Reduce("+",lapply(1:nrow(pD), function(m){
 #  Reduce("+",lapply(1:nrow(ghQ$points), function(i) tcrossprod(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i])) }))
 # sc <- -nrow(pD)/2*solve(pz) + 0.5*solve(pz)%*%V%*%solve(pz)
 # sc <- 2*sc - diag(diag(sc))
 # sc <- c(sc[lower.tri(sc)],sc[upper.tri(sc)])
 # kR <- kronecker(solve(pz),solve(pz))
 # for(y in 1:nrow(corlst)){
 #  idr <- (((y-1)*nrow(pz))+nrow(pz)+1):((y)*nrow(pz)+nrow(pz))
 #  idc <- (((y-1)*nrow(pz))+nrow(pz)+1):((y)*nrow(pz)+nrow(pz))
 #  he <- nrow(pD)/2*(kR[id,id]) - 0.5*(solve(pz)%*%solve(pz)%*%V%*%solve(pz) - solve(pz)%*%V%*%solve(pz)%*%solve(pz)) 
 # }
 # 
 # # he <- 2*he - diag(diag(he))
 # # for(y in 1:nrow(corlst)){
 # #  pz[corlst[y,1],corlst[y,2]] <- pz[corlst[y,1],corlst[y,2]] - solve(he)[corlst[y,1],corlst[y,2]]*sc[corlst[y,1],corlst[y,2]]
 # # }
 # # pz[upper.tri(pz)] <- pz[lower.tri(pz)]
 pz <- (1/nrow(pD))*Reduce("+",lapply(1:nrow(pD), function(m){
  Reduce("+",lapply(1:nrow(ghQ$points), function(i) tcrossprod(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i]))}))
 pz <- diag(1/sqrt(diag(pz)))%*%pz%*%diag(1/sqrt(diag(pz)))
 return(pz)
}

upS <- function(pz,ghQ,pD,respz){
 if(is.null(respz)){ corlst <- which(lower.tri(pz),arr.ind = T) } else corlst <- which(lower.tri(pz),arr.ind = T)[-respz,]
 sv <- pz[lower.tri(pz)]
 gra <- matrix(0,nrow = length(sv))
 hess <- matrix(0, nrow = length(sv), ncol = length(sv))
 Iq <- diag(nrow(pz))
 for(y in 1:nrow(corlst)){
  Dy <- matrix(0,nrow = nrow(pz),ncol(pz)); Dy[corlst[y,1],corlst[y,2]] <- 1
  Dy <- Dy + t(Dy) - Dy%*%Dy
  Gy <- solve(pz)%*%Dy%*%solve(pz)
  # Fy <- solve(pz)%*%Dy%*%Gy + Gy%*%Dy%*%solve(pz)
  P1 <- -nrow(pD)/2*sum(diag((2*solve(pz) - solve(pz)*Iq)%*%Dy))
  P2 <- Reduce("+",lapply(1:nrow(pD), function(m){sum(diag(
        Gy%*%Reduce("+",lapply(1:nrow(ghQ$points), function(i) tcrossprod(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i])) )) })) 
  P3 <- Reduce("+",lapply(1:nrow(pD), function(m){
        v <- Reduce("+",lapply(1:nrow(ghQ$points), function(i) matrix(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i]));
        t(v)%*%Gy%*%v }))
  gra[y] <- P1+P2+P3
  for(yy in 1:nrow(corlst)){
   Dyy <- matrix(0,nrow = nrow(pz),ncol(pz)); Dyy[corlst[yy,1],corlst[yy,2]] <- 1
   Dyy <- Dyy + t(Dyy) - Dyy%*%Dyy
   alp <- sum(diag(-solve(pz)%*%solve(pz)%*%Dyy))
   Fy <- solve(pz)%*%Dy + Dy%*%solve(pz)
   H1 <- -nrow(pD)/2*sum(diag(alp*Dy))
   H2 <- 0.5*alp*Reduce("+",lapply(1:nrow(pD), function(m){sum(diag(
    Fy%*%Reduce("+",lapply(1:nrow(ghQ$points), function(i) tcrossprod(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i])) )) }))
   H3 <- 0.5*alp*Reduce("+",lapply(1:nrow(pD), function(m){
    v <- Reduce("+",lapply(1:nrow(ghQ$points), function(i) matrix(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i]));
    t(v)%*%Fy%*%v }))
   hess[yy,y] <- hess[y,yy] <- H1+H2+H3
  }
 }
 pz[lower.tri(pz)] <- matrix(sv) - solve(hess)%*%gra
 
  # 
  # H1 <- nrow(pD)/2*sum(diag(Gy%*%Dy))
  # H2 <- -0.5*Reduce("+",lapply(1:nrow(pD), function(m){sum(diag(
  # Fy%*%Reduce("+",lapply(1:nrow(ghQ$points), function(i) tcrossprod(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i])) )) }))
  # H3 <- -0.5*Reduce("+",lapply(1:nrow(pD), function(m){
  #   v <- Reduce("+",lapply(1:nrow(ghQ$points), function(i) matrix(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i]));
  #   t(v)%*%Fy%*%v }))
  # pz[corlst[y,1],corlst[y,2]] <- pz[corlst[y,1],corlst[y,2]] - c(P1+P2+P3)/c(H1+H2+H3)
  # 
  # gradient <- c(gradient, P1+P2+P3)
  # hessian <- c(hessian, H1+H2+H3)
 # }
 # if(length(hessian) == 1) hessian <- matrix(hessian) else hessian <- diag(hessian)
 pz[upper.tri(pz)] <- pz[lower.tri(pz)]
 return(list(pz = pz, gradient = c(gra), hessian = hess))
}

penM <- function(params, type = "lasso", id = pen.idx, rs = loadmt, lambda = 1,
                 w.alasso = NULL,a = NULL){  # , pz, clv

# Goal: To produce a Penalty matrix of parameters (params)
# Input : model parameters (params)
# Output: Penalty matrix
# Testing: params = lb2cb(bold); lambda = pml.control$lambda; id = pen.idx; pz = sigZ; clv = T;
#          w.alasso = pml.control$w.alasso; rs = loadmt2; type = pml.control$type; a = NULL
 
if(!is.null(w.alasso) && is.list(w.alasso)){ 
  # if(clv) w.alasso <- c(t(lb2mb(w.alasso))[t(id)],pz[lower.tri(pz)]) else 
  w.alasso <- c(t(lb2mb(w.alasso))[t(id)]) }
# if(clv){S. <- diag(0,length(c(params,pz[lower.tri(pz)])))} else 
  S. <- diag(0,length(c(params)))
# if(clv){param <- c(params[c(t(id))],pz[lower.tri(pz)])} else 
  param <- params[c(t(id))]
eps = 1e-7 # sqrt(.Machine$double.eps) # protective tolerance level
if(type == "ridge"){
 A1 <-  c(lambda)*rep(1, length(param))
}
if(type == "lasso"){
 A1 <- c(lambda)/sqrt(param^2 + eps) 
}
if(type == "alasso"){
 if(is.null(a)) a = 2
 if( is.null(w.alasso) ) w.alasso <- 1
 w.al <- abs(w.alasso)^a
 A1 <- c(lambda)/(w.al*sqrt(param^2 + eps))
}
if(type == "scad"){
 if(is.null(a)) a = 3.7
 theta <- abs(param)
 f1 <- sapply(theta, function(theta) { max(a*c(lambda) - theta, 0)/((a-1)*c(lambda) + eps) })
 f.d <- ((theta <= c(lambda)) + f1 * (theta > c(lambda)))
 A1 <- c(lambda)* f.d / ( sqrt(params^2 + eps) )
}
if(type == "mcp"){
 if(is.null(a)) a = 2.5
 theta <- abs(param) 
 f.d <- (c(lambda)-theta/a)*(theta < c(lambda)*a)
 A1 <- f.d / ( sqrt(param^2 + eps) )
}
if(length(A1) == 1) S <- matrix(A1) else S <- diag(A1)
#if(clv){diag(S.)[c(t(id),rep(T,nrow(S.)-length(t(id))))] <- A1} else 
diag(S.)[c(t(id))] <- A1
#if(clv){S. <- S.[c(t(lb2mb(rs)),rep(T,nrow(S.)-length(t(id)))),c(t(lb2mb(rs)),rep(T,nrow(S.)-length(t(id))))]} else 
S. <- S.[t(lb2mb(rs)),t(lb2mb(rs))]
return(list(full = S., red = S))
}

uplm <- function(b,loadmt){
  
# Goal: To update loadmt (with new zeroes in b)
# Input : b (parameters) and loadmt (restrictions)
# Output: loadmt
# Testing: b = bold;
 
lmt <- loadmt
for(i in names(b)){ lmt[[i]] <- (b[[i]] != 0) & loadmt[[i]] }
return(lmt)
}

op.lambda <- function(b,Y,idx,rs,pz,clv,pml,shObj,itl,tol = 1e-3){
  
# Goal: to update automatic lambda
# Input: lambda, updated-Beta (penalised) obj, penalty object, restriction matrix (rs)
# Output: lambda (updated)
# Testing: b = bnew; Y = Y; idx = pen.idx; rs = loadmt2; pml = pml.control; shObj = A3; itl = control$iter.lim

it <- 0
es <- 1
nY <- nrow(Y)
iT <- T
# step <- 1

pmObj <- lb2pM(b,Y,idx,rs,pml)
pM <- pmObj$pM/c(nY*pml$lambda)
sJ <- m2pdm(-shObj$hessian)$sq
Q <- qr.Q(qr(sJ))
R <- qr.R(qr(sJ))
B <- m2pdm(pmObj$pM)$sq
th <- t(lb2mb(b))[t(lb2mb(rs))] # t(lb2mb(b))[lb2mb(rs)]
svdRB <- svd(rbind(R,B))
U1 <- svdRB$u[1:nrow(R),]
V <- svdRB$v
D <- diag(svdRB$d)
iD <- m2pdm(D)$inv.mat
K = sJ%*%th + solve(t(sJ))%*%shObj$gradient
K1 = t(U1)%*%t(Q)%*%K
Zl = iD%*%t(V)%*%pM%*%V%*%iD
Cl = Zl%*%crossprod(U1)
A <- sJ%*%m2pdm(-shObj$hessian + pmObj$pM)$inv.mat%*%t(sJ)
MSEo <- c(crossprod(K-A%*%K)) + (pml$gamma)*2*sum(diag(A)) - sum(lb2mb(rs)) # /sum(lb2mb(rs)
   
while(it < itl && es > tol){

lambda1 <- pml$lambda

dtrA = -pml$lambda*nY*sum(diag(Cl))
d2trA = 2*(pml$lambda*nY)^2*sum(diag(Zl%*%Cl)) + dtrA
dP = 2*pml$lambda*nY*(t(K1)%*%Zl%*%K1 - t(K1)%*%Cl%*%K1)
d2P = 2*(pml$lambda*nY)^2*t(K1)%*%(Zl%*%Cl + Zl%*%Cl - Zl%*%Zl - Zl%*%Zl + Cl%*%Zl)%*%K1 + dP

dV = (dP + (pml$gamma)*2*dtrA)#/sum(lb2mb(rs))
d2V = (d2P + (pml$gamma)*2*d2trA)#/sum(lb2mb(rs))
i2dV <- if(iT){ tryCatch({solve(d2V)}, error = function(e){return(1/d2V)}) } else { i2dV/2 } # 0.1*1/it
lambda2 <- c(exp(log(pml$lambda) - i2dV%*%dV))
es <- abs(lambda2 - pml$lambda)
it = it + 1 
pml$lambda <- lambda2

# Test if MSE is reduced

pmObj <- lb2pM(b,Y,idx,rs,pml)
pM <- pmObj$pM/c(nY*pml$lambda)
B <- m2pdm(pmObj$pM)$sq
th <- t(lb2mb(b))[t(lb2mb(rs))]
svdRB <- svd(rbind(R,B))
U1 <- svdRB$u[1:nrow(R),]
V <- svdRB$v
D <- diag(svdRB$d)
iD <- m2pdm(D)$inv.mat
K = sJ%*%th + solve(t(sJ))%*%shObj$gradient
K1 = t(U1)%*%t(Q)%*%K
Zl = iD%*%t(V)%*%pM%*%V%*%iD
Cl = Zl%*%crossprod(U1)
A <- sJ%*%m2pdm(-shObj$hessian + pmObj$pM)$inv.mat%*%t(sJ)
MSEn <- c(crossprod(K-A%*%K)) + (pml$gamma)*2*sum(diag(A)) - sum(lb2mb(rs)) # /sum(lb2mb(rs)

if(MSEo < MSEn){

if(iT){ i2dV <- 0.1 }
iT <- F
it = it - 1   
pml$lambda <- lambda1

pmObj <- lb2pM(b,Y,idx,rs,pml)
pM <- pmObj$pM/c(nY*pml$lambda)
sJ <- m2pdm(-shObj$hessian)$sq
Q <- qr.Q(qr(sJ))
R <- qr.R(qr(sJ))
B <- m2pdm(pmObj$pM)$sq
th <- t(lb2mb(b))[t(lb2mb(rs))] # t(lb2mb(b))[lb2mb(rs)]
svdRB <- svd(rbind(R,B))
U1 <- svdRB$u[1:nrow(R),]
V <- svdRB$v
D <- diag(svdRB$d)
iD <- m2pdm(D)$inv.mat
K = sJ%*%th + solve(t(sJ))%*%shObj$gradient
K1 = t(U1)%*%t(Q)%*%K
Zl = iD%*%t(V)%*%pM%*%V%*%iD
Cl = Zl%*%crossprod(U1)
# A <- sJ%*%m2pdm(-shObj$hessian + pmObj$pM)$inv.mat%*%t(sJ)
# MSEo <- c(crossprod(K-A%*%K)/sum(lb2mb(rs)) + 2*sum(diag(A))/sum(lb2mb(rs)) - 1)

dtrA = -pml$lambda*nY*sum(diag(Cl))
d2trA = 2*(pml$lambda*nY)^2*sum(diag(Zl%*%Cl)) + dtrA
dP = 2*pml$lambda*nY*(t(K1)%*%Zl%*%K1 - t(K1)%*%Cl%*%K1)
d2P = 2*(pml$lambda*nY)^2*t(K1)%*%(Zl%*%Cl + Zl%*%Cl - Zl%*%Zl - Zl%*%Zl + Cl%*%Zl)%*%K1 + dP

dV = (dP + (pml$gamma)*2*dtrA)#/sum(lb2mb(rs))
d2V = (d2P + (pml$gamma)*2*d2trA)#/sum(lb2mb(rs))
#i2dV <- 1
lambda2 <- c(exp(log(pml$lambda) - i2dV%*%dV))
es <- abs(lambda2 - pml$lambda)
it = it + 1 
pml$lambda <- lambda2 

pmObj <- lb2pM(b,Y,idx,rs,pml)
pM <- pmObj$pM/c(nY*pml$lambda)
B <- m2pdm(pmObj$pM)$sq
th <- t(lb2mb(b))[t(lb2mb(rs))]
svdRB <- svd(rbind(R,B))
U1 <- svdRB$u[1:nrow(R),]
V <- svdRB$v
D <- diag(svdRB$d)
iD <- m2pdm(D)$inv.mat
K = sJ%*%th + solve(t(sJ))%*%shObj$gradient
K1 = t(U1)%*%t(Q)%*%K
Zl = iD%*%t(V)%*%pM%*%V%*%iD
Cl = Zl%*%crossprod(U1)
A <- sJ%*%m2pdm(-shObj$hessian + pmObj$pM)$inv.mat%*%t(sJ)
MSEn <- c(crossprod(K-A%*%K) + (pml$gamma)*2*sum(diag(A)) - sum(lb2mb(rs)))

}

MSEo <- MSEn

}

if(it == itl) it <- paste0("max ~ ",itl)
return(list(lambda = pml$lambda, iter = it, sse = MSEn))
}
