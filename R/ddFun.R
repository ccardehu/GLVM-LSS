
coefmod <- function(bet,beta){ # bet = unlist(borg); beta = borg; gr = gr
  coefM <- NULL
  tmpl <- NULL
  for(i in 1:length(beta)){
    if(i == 1) tmpl[[i]] <- 1:length(beta[[i]])
    else tmpl[[i]] <- seq(from = (length(beta[[i-1]]) + 1), length.out = length(beta[[i]]))
  }
  for(i in 1:length(beta)){
    coefM[[i]] <- matrix(bet[tmpl[[i]]],nrow = p, ncol = ncol(beta[[i]])) #*loadmt[[i]]
    colnames(coefM[[i]]) <- colnames(beta[[i]]); rownames(coefM[[i]]) <- rownames(beta[[i]])
  }
  names(coefM) <- names(beta)
  return(coefM)
}

coefmod2 <- function(bet,beta){ # bet = unlist(borg); beta = borg; gr = gr
 
 tmpl <- NULL
 for(i in 1:length(beta)){
  if(i == 1) tmpl[[i]] <- 1:ncol(beta[[i]])
  else tmpl[[i]] <- seq(from = tail(tmpl[[i-1]],1)+1, to = tail(tmpl[[i-1]],1) + ncol(beta[[i]]))
 }
 
 bet <- matrix(bet,byrow = T,nrow = nrow(beta[[1]]))
 coefM <- NULL
 for(i in 1:length(beta)){
  coefM[[i]] <- bet[,tmpl[[i]],drop=F]
  colnames(coefM[[i]]) <- colnames(beta[[i]]); rownames(coefM[[i]]) <- rownames(beta[[i]])
 }
 names(coefM) <- names(beta)
 return(coefM)
}

coefmod3 <- function(beta){ # bet = unlist(borg); beta = borg; gr = gr
r <- NULL
for(i in 1:length(beta)){ r <- cbind(r,beta[[i]]) }
return(unname(r))
}

decoefmod <- function(beta){
 tmp <- NULL
 for(i in names(beta)){tmp <- cbind(tmp,beta[[i]])}
 return(c(t(tmp)))
}

probs <- function(x){
  pr <- plogis(x)
  if (any(ind <- pr == 1)) pr[ind] <- 1 - sqrt(.Machine$double.eps)
  if (any(ind <- pr == 0)) pr[ind] <- sqrt(.Machine$double.eps)
  return(pr)
}

expit <- function(x) exp(x)/(1+exp(x)); logit <- function(x) log(x/(1-x))

# To use in fFun() below and outside

Zm <- function(lv = c("Z1"), Z. = Z){
# Z. = simR$Z$sigma; lv = "Z1"
  tmpr <- vector(mode = "list", length = length(lv)); names(tmpr) <- lv
  for(i in lv){ tmpr[[i]] <- grep(i, colnames(Z.), fixed = T) }
  rmv <- unname(unlist(tmpr))
  if(length(rmv) != 0) Z.[,-c(1,rmv)] <- 0 else Z.[,-1] <- 0
  return(Z.)
} # Returns Z matrix for plotting marginals with all but "lv" = 0

pFun <- function(fam){ # should be evaluated at i = 1:p; fam = fam[i]
 if(fam == "normal") pars <- c("mu","sigma")
 if(fam == "lognormal") pars <- c("mu","sigma")
 if(fam == "poisson") pars <- c("mu")
 if(fam == "gamma") pars <- c("mu", "sigma")
 if(fam == "binomial") pars <- c("mu")
 if(fam == "ZIpoisson") pars <- c("mu", "sigma")
 # Add other parameters for other distributions
 return(pars)
}

rFun <- function(n,i,fam,Z,b){ #should be evaluated at i = 1:p; fam = fam[i]; Z = simR$Z; b = borg
  
 if(fam == "normal"){
  mu = as.matrix(Z$mu)%*%matrix(b$mu[i,])
  sigma = exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))
  fyz <- rnorm(n, c(mu), c(sigma))
 }
 if(fam == "lognormal"){
  mu = as.matrix(Z$mu)%*%matrix(b$mu[i,])
  sigma = exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))
  fyz <- rlnorm(n, drop(mu), drop(sigma))
 }
 if(fam == "poisson"){
  mu = exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
  fyz <- rpois(n, c(mu)) 
 }
 if(fam == "gamma"){
  mu = exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
  sigma = exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))
  fyz <- rgamma(n,shape = mu, scale = sigma)
 }
 if(fam == "binomial"){
  mu = probs(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
  fyz <- rbinom(n,1,prob = mu)
 }  
 if(fam == "ZIpoisson"){
  mu = exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
  sigma = probs(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))
  fyz <- ifelse(rbinom(n,1,1-sigma) == 0, 0, rpois(n,mu))
 }
 # Add other distributions in SimFA::fod
 return(fyz)
}

gFun <- function(Z,b,fam,i){
  # To be evaluated at:
  # Z <- gr$out; b <- borg; fam = fam[i]
  gout <- vector(mode = "list", length = length(pFun(fam))); names(gout) <- pFun(fam)
  if(fam == "normal"){
   gout$mu = drop(unname(as.matrix(Z$mu)%*%matrix(b$mu[i,])))
   gout$sigma = drop(unname(exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))))
  }
  if(fam == "lognormal"){
   gout$mu = drop(unname(as.matrix(Z$mu)%*%matrix(b$mu[i,])))
   gout$sigma = drop(unname(exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))))
  }
  if(fam == "poisson"){ gout$mu = drop(unname(exp(as.matrix(Z$mu)%*%matrix(b$mu[i,])))) }
  if(fam == "gamma"){
    gout$mu = drop(unname(exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))))
    gout$sigma = drop(unname(exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))))
  }
  if(fam == "binomial"){ gout$mu = drop(unname(probs(as.matrix(Z$mu)%*%matrix(b$mu[i,])))) }
  if(fam == "ZIpoisson"){
    gout$mu = drop(unname(exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))))
    gout$sigma = drop(unname(probs(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))))
  }
 return(gout)
}

dFun2 <- function(Y,fam,Z,b,j){
dY <- Y
for(i in 1:ncol(Y)){
 if(fam[i] == "normal"){
   dY[,i] <- dnorm(Y[,i], gFun(Z,b,fam[i],i)$mu[j], gFun(Z,b,fam[i],i)$sigma[j], log = T)
  }
 if(fam[i] == "lognormal"){
   dY[,i] <- dlnorm(Y[,i], gFun(Z,b,fam[i],i)$mu[j], gFun(Z,b,fam[i],i)$sigma[j], log = T)
  }
 if(fam[i] == "poisson"){
   dY[,i] <- dpois(Y[,i], gFun(Z,b,fam[i],i)$mu[j], log = T) 
  }
 if(fam[i] == "gamma"){
   dY[,i] <- dgamma(Y[,i], shape = gFun(Z,b,fam[i],i)$mu[j], scale = gFun(Z,b,fam[i],i)$sigma[j], log = T)
  }
 if(fam[i] == "binomial"){
   dY[,i] <- dbinom(Y[,i],size = 1,prob = gFun(Z,b,fam[i],i)$mu[j], log = T)
  }
 if(fam[i] == "ZIpoisson"){
   dZIpoisson <- function(Y,mu,sigma, log = T){
    u <- as.numeric(Y == 0)
    lf <- u*c(log(sigma + (1-sigma)*exp(-mu))) + (1-u)*(c(log(1-sigma)) - c(mu) + Y*c(log(mu)) - lfactorial(Y))
    if(log == T) return(lf) else return(exp(lf))
   }
   dY[,i] <- dZIpoisson(Y[,i], gFun(Z,b,fam[i],i)$mu[j], gFun(Z,b,fam[i],i)$sigma[j], log = T)
  }
 # Add other distributions
}
return(dY)
}

dFun <- function(i,z,fam,Y,Z,b){ #should be evaluated at i = 1:p; z = 10; fam = fam[i]; Y = simR$Y[,i]; Z = gr$out; b = borg
 # if(!is.matrix(Z)) Z <- as.matrix(Z)
 # if(ncol(Z) == 1) Z <- t(Z)
  
 if(fam == "normal"){
  mu = as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])
  sigma = exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,]))
  fyz <- dnorm(Y, c(mu), c(sigma), log = T)
 }
 if(fam == "lognormal"){
  mu = as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])
  sigma = exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,]))
  fyz <- dlnorm(Y, c(mu), c(sigma), log = T)
 }
 if(fam == "poisson"){
  mu = exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
  fyz <- dpois(Y, c(mu), log = T) 
 }
 if(fam == "gamma"){
  mu = exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
  sigma = exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,]))
  fyz <- dgamma(Y,shape = c(mu), scale = c(sigma), log = T)
 }
 if(fam == "binomial"){
  mu = probs(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
  fyz <- dbinom(Y,size = 1,prob = c(mu), log = T)
 }
 if(fam == "ZIpoisson"){
  mu = exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
  sigma = probs(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,]))
  dZIpoisson <- function(Y,mu,sigma,log = T){
   u <- as.numeric(Y == 0)
   lf <- u*c(log(sigma + (1-sigma)*exp(-mu))) + (1-u)*(c(log(1-sigma)) - c(mu) + Y*c(log(mu)) - lfactorial(Y))
   if(log == T) return(lf) else return(exp(lf))
  }
  fyz <- dZIpoisson(Y, c(mu), c(sigma))
  #fyz <- gamlss.dist::dZIP(Y, drop(mu), drop(sigma), log = T)
 }
 # Add other distributions
 return(fyz)
}

dvFun <- function(i,z,fam,Y,Z,b){ # i = 1; z = 10; fam = fam[i]; Y = simR$Y[,i]; Z = gr$out; b = borg
 # if(!is.matrix(Z)) Z <- as.matrix(Z)
 # if(ncol(Z) == 1) Z <- t(Z)
  
 if(fam == "normal"){ # structure(gamlss.dist::NO)
  mu = c(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
  sigma = c(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
  d1mu = (Y - mu)/(sigma^2)  * 1
  d2mu = rep(-(1/sigma^2) * 1, length(Y))# + 0
  d1sg = (((Y - mu)^2 - sigma^2)/(sigma^3)) * sigma
  d2sg = (-3*(Y-mu)^2/sigma^4 + 1/sigma^2) * sigma^2 + d1sg # rep(-(2/(sigma^2)), length(Y))* sigma^2 # 
  dcms = -2*(Y-mu)/sigma^3 * 1 * sigma # rep(0, length(Y)) # the expected cross derivative mu and sigma
  rL <- list("d1mu" = d1mu, "d2mu" = d2mu, "d1sg" = d1sg, "d2sg" = d2sg, "dcms" = dcms)
 }
 if(fam == "lognormal"){ # structure(gamlss.dist::LOGNO)
  mu = drop(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
  sigma = drop(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
  d1mu = (log(Y) - mu)/(sigma^2)  * 1
  d2mu = rep(-(1/sigma^2) * 1, length(Y))
  d1sg = (((log(Y) - mu)^2 - sigma^2)/(sigma^3)) * sigma
  d2sg = (-3*(log(Y)-mu)^2/sigma^4 + 1/sigma^2) * sigma^2 + d1sg # rep(-(2/(sigma^2)) * sigma^2, length(Y)) # -3*(Y-mu)^2/sigma^4 + 1/sigma^2
  dcms = -2*(log(Y)-mu)/sigma^3 * 1 * sigma # rep(0, length(Y)) # the expected cross derivative mu and sigma
  rL <- list("d1mu" = d1mu, "d2mu" = d2mu, "d1sg" = d1sg, "d2sg" = d2sg, "dcms" = dcms)
 }
 if(fam == "poisson"){ # structure(gamlss.dist::PO)
    mu = exp(c(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
    d1mu = ((Y - mu)/mu) * mu
    d2mu = -Y/mu^2 * mu^2 + d1mu # rep((-1/mu) * mu^2, length(Y))
    rL <- list("d1mu" = d1mu, "d2mu" = d2mu)
  }
 if(fam == "gamma"){ # structure(gamlss.dist::GA)
    mu = c(exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
    sigma = c(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
    d1mu = (log(Y) - digamma(mu) - log(sigma)) * mu
    d2mu = rep(-trigamma(mu) * mu^2, length(Y)) + d1mu
    d1sg = (Y/(sigma^2) - mu/sigma) * sigma
    d2sg = ((-2*Y+mu*sigma)/sigma^3) * sigma^2 + d1sg # rep((-mu/(sigma)) * sigma^2, length(Y)) # 
    dcms = rep(-1/sigma*mu*sigma,length(Y)) # rep(0, length(Y))  # the expected cross derivative mu and sigma
    rL <- list("d1mu" = d1mu, "d2mu" = d2mu, "d1sg" = d1sg, "d2sg" = d2sg, "dcms" = dcms)
  }
 if(fam == "binomial"){ # structure(gamlss.dist::BI)
    mu = probs(c(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
    d1mu = (Y - mu)/(mu * (1 - mu)) * (mu*(1-mu))
    d2mu = -(mu^2 - 2*Y*mu + Y)/((mu-1)*mu)^2 * (mu*(1-mu))^2 + (1-2*mu)*d1mu# rep(-(1/(mu * (1 - mu))) * (mu*(1-mu))^2, length(Y))
    rL <- list("d1mu" = d1mu, "d2mu" = d2mu)
  }
 if(fam == "ZIpoisson"){ # structure(gamlss.dist::ZIP)
    mu = exp(c(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
    sigma = probs(c(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
    u = ifelse(Y == 0, (1 + exp(-logit(sigma)-mu))^(-1),0) # as.numeric(Y == 0)
    d1mu = (1-u)*(Y/mu-1) * mu
    d2mu = (1-u)*(-Y/mu^2) * mu^2 + d1mu # expected value is -1/mu; original is (-Y/mu^2)
    d1sg = (u-sigma)/(sigma*(1-sigma)) * (sigma*(1-sigma))
    d2sg = -(sigma^2 - 2*u*sigma + u)/((sigma-1)*sigma)^2 * (sigma*(1-sigma))^2 + (1-2*sigma)*d1sg # -(u/(sigma^2) + (1-u)/((1-sigma)^2)) * (sigma*(1-sigma))^2
    dcms = rep(0, length(Y))
    rL <- list("d1mu" = d1mu, "d2mu" = d2mu, "d1sg" = d1sg, "d2sg" = d2sg, "dcms" = dcms)
  }
  
 return(rL)
}

# Function to compute E(Y) and SD(Y) (or quantiles)

fFun <- function(i,fam,Z,b,qnt = c(0.2,0.4,0.6,0.8),forms,lvp){
  #should be evaluated at i = 1:p; fam = fam[i]; Z = simR$Z; b = borg; qnt = c(0.2,0.4,0.6,0.8) forms = ex1$formula
 if(missing(lvp)){ lvp = c(1); rtF = F } else rtF = T
 if(missing(qnt)){ qnt = c(0.2,0.8)}
 pars <- pFun(fam)
 lvar <- unique(unlist(lapply(pars, function(i) all.vars(forms[[i]]))))
 lvar <- grep("Z", lvar, fixed = T, value = T); qMM <- list()

 if(fam == "normal"){
  mu = drop(unname(as.matrix(Z$mu)%*%matrix(b$mu[i,])))
  sigma = drop(unname(exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))))
  EY = mu
  SY = sigma
  qM = matrix(nrow = length(mu), ncol = length(qnt)); colnames(qM) <- paste0("q",qnt*100)
  for(j in seq_along(qnt)){ qM[,j] <- qnorm(qnt[j],mu,sigma) }
   eM <- sM <- matrix(nrow = length(mu), ncol = length(lvar)); colnames(eM) <- colnames(sM) <- lvar
  for(k in lvar){
   muM = drop(unname(as.matrix(Zm(k,Z$mu))%*%matrix(b$mu[i,])))
   sgM = drop(unname(exp(as.matrix(Zm(k,Z$sigma))%*%matrix(b$sigma[i,]))))
   eM[,k] <- muM
   sM[,k] <- sgM
   qMM[[k]] <- matrix(nrow = length(mu), ncol = length(qnt)); colnames(qMM[[k]]) <- paste0("q",qnt*100)
   for(j in seq_along(qnt)){ qMM[[k]][,j] <- qnorm(qnt[j],muM,sgM) }
  }
  lvpl = paste0("Z",lvp)
  mu2M = drop(unname(as.matrix(Zm(lvpl,Z$mu))%*%matrix(b$mu[i,])))
  sg2M = drop(unname(exp(as.matrix(Zm(lvpl,Z$sigma))%*%matrix(b$sigma[i,]))))
  EY2M = mu2M
 }
 if(fam == "lognormal"){
  mu = drop(unname(as.matrix(Z$mu)%*%matrix(b$mu[i,])))
  sigma = drop(unname(exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))))
  EY = exp(mu + 0.5*sigma^2)
  SY = sqrt((exp(sigma^2)-1)*exp(2*mu+sigma^2))
  qM = matrix(nrow = length(mu), ncol = length(qnt)); colnames(qM) <- paste0("q",qnt*100)
  for(j in seq_along(qnt)){ qM[,j] <- qlnorm(qnt[j],mu,sigma) }
  eM <- sM <- matrix(nrow = length(mu), ncol = length(lvar)); colnames(eM) <- colnames(sM) <- lvar
  for(k in lvar){
   muM = drop(unname(as.matrix(Zm(k,Z$mu))%*%matrix(b$mu[i,])))
   sgM = drop(unname(exp(as.matrix(Zm(k,Z$sigma))%*%matrix(b$sigma[i,]))))
   eM[,k] <- muM
   sM[,k] <- sgM
   qMM[[k]] <- matrix(nrow = length(mu), ncol = length(qnt)); colnames(qMM[[k]]) <- paste0("q",qnt*100)
   for(j in seq_along(qnt)){ qMM[[k]][,j] <- qlnorm(qnt[j],muM,sgM) }
  }
  lvpl = paste0("Z",lvp)
  mu2M = drop(unname(as.matrix(Zm(lvpl,Z$mu))%*%matrix(b$mu[i,])))
  sg2M = drop(unname(exp(as.matrix(Zm(lvpl,Z$sigma))%*%matrix(b$sigma[i,]))))
  EY2M = mu2M
 }
 if(fam == "poisson"){
    mu = drop(unname(exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))))
    EY = mu
    SY = sqrt(mu)
    qM = matrix(nrow = length(mu), ncol = length(qnt)); colnames(qM) <- paste0("q",qnt*100)
    for(i in seq_along(qnt)){ qM[,i] <- qpois(qnt[i],mu) }
    eM <- sM <- matrix(nrow = length(mu), ncol = length(lvar)); colnames(eM) <- colnames(sM) <- lvar
    for(k in lvar){
     muM = drop(unname(exp(as.matrix(Zm(k,Z$mu))%*%matrix(b$mu[i,]))))
     eM[,k] <- muM
     sM[,k] <- sqrt(muM)
     qMM[[k]] <- matrix(nrow = length(mu), ncol = length(qnt)); colnames(qMM[[k]]) <- paste0("q",qnt*100)
     for(j in seq_along(qnt)){ qMM[[k]][,j] <- qpois(qnt[j],muM) }
    }
    lvpl = paste0("Z",lvp)
    mu2M = drop(unname(exp(as.matrix(Zm(lvpl,Z$mu))%*%matrix(b$mu[i,]))))
    EY2M = mu2M
  }
 if(fam == "gamma"){
    mu = drop(unname(exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))))
    sigma = drop(unname(exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))))
    EY = mu*sigma
    SY = sqrt(mu*sigma^2)
    qM = matrix(nrow = length(mu), ncol = length(qnt)); colnames(qM) <- paste0("q",qnt*100)
    for(i in seq_along(qnt)){ qM[,i] <- qgamma(qnt[i], shape = mu, scale = sigma) }
    eM <- sM <- matrix(nrow = length(mu), ncol = length(lvar)); colnames(eM) <- colnames(sM) <- lvar
    for(k in lvar){
     muM = drop(unname(exp(as.matrix(Zm(k,Z$mu))%*%matrix(b$mu[i,]))))
     sgM = drop(unname(exp(as.matrix(Zm(k,Z$sigma))%*%matrix(b$sigma[i,]))))
     eM[,k] <- muM*sgM
     sM[,k] <- sqrt(muM*sgM^2)
     qMM[[k]] <- matrix(nrow = length(mu), ncol = length(qnt)); colnames(qMM[[k]]) <- paste0("q",qnt*100)
     for(j in seq_along(qnt)){ qMM[[k]][,j] <- qgamma(qnt[j],  shape = muM, scale = sgM) }
    }
    lvpl = paste0("Z",lvp)
    mu2M = drop(unname(exp(as.matrix(Zm(lvpl,Z$mu))%*%matrix(b$mu[i,]))))
    sg2M = drop(unname(exp(as.matrix(Zm(lvpl,Z$sigma))%*%matrix(b$sigma[i,]))))
    EY2M = mu2M*sg2M
  }
 if(fam == "binomial"){
    mu = drop(unname(probs(as.matrix(Z$mu)%*%matrix(b$mu[i,]))))
    EY = mu
    SY = sqrt(mu*(1-mu))
    qM = matrix(nrow = length(mu), ncol = length(qnt)); colnames(qM) <- paste0("q",qnt*100)
    for(i in seq_along(qnt)){ qM[,i] <- qbinom(qnt[i], 1, mu) }
    eM <- sM <- matrix(nrow = length(mu), ncol = length(lvar)); colnames(eM) <- colnames(sM) <- lvar
    for(k in lvar){
     muM = drop(unname(probs(as.matrix(Zm(k,Z$mu))%*%matrix(b$mu[i,]))))
     eM[,k] <- muM
     sM[,k] <- sqrt(muM*(1-muM))
     qMM[[k]] <- matrix(nrow = length(mu), ncol = length(qnt)); colnames(qMM[[k]]) <- paste0("q",qnt*100)
     for(j in seq_along(qnt)){ qMM[[k]][,j] <- qbinom(qnt[j],1,muM) }
    }
    lvpl = paste0("Z",lvp)
    mu2M = drop(unname(probs(as.matrix(Zm(lvpl,Z$mu))%*%matrix(b$mu[i,]))))
    EY2M = mu2M
  }
 if(fam == "ZIpoisson"){
    mu = drop(unname(exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))))
    sigma = drop(unname(probs(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))))
    qZIpoisson <- function (p,mu,sigma){
      ly <- max(length(p), length(mu), length(sigma))
      p <- rep(p, length = ly); sigma <- rep(sigma, length = ly); mu <- rep(mu, length = ly)
      pnew <- ((p - sigma)/(1 - sigma)) - (1e-07)
      pnew <- ifelse(pnew > 0, pnew, 0)
      #q <- qpois(pnew, lambda = mu)
      q <- ifelse(p > sigma*(1-sigma)*dpois(0,mu), qpois(pnew, mu), 0)
      return(q)
    }
    # https://stats.stackexchange.com/questions/463844/is-there-an-equation-for-the-median-and-percentile-of-a-zero-inflated-poisson
    EY = (1-sigma)*mu
    SY = sqrt(mu*(1-sigma)*(1+mu*sigma))
    qM = matrix(nrow = length(mu), ncol = length(qnt)); colnames(qM) <- paste0("q",qnt*100)
    for(i in seq_along(qnt)){ qM[,i] <- qZIpoisson(qnt[i], mu, sigma) }
    eM <- sM <- matrix(nrow = length(mu), ncol = length(lvar)); colnames(eM) <- colnames(sM) <- lvar
    for(k in lvar){
     muM = drop(unname(exp(as.matrix(Zm(k,Z$mu))%*%matrix(b$mu[i,]))))
     sgM = drop(unname(probs(as.matrix(Zm(k,Z$sigma))%*%matrix(b$sigma[i,]))))
     eM[,k] <- (1-sgM)*muM
     sM[,k] <- sqrt(muM*(1-sgM)*(1+muM*sgM))
     qMM[[k]] <- matrix(nrow = length(mu), ncol = length(qnt)); colnames(qMM[[k]]) <- paste0("q",qnt*100)
     for(j in seq_along(qnt)){ qMM[[k]][,j] <- qZIpoisson(qnt[j],muM,sgM) }
    }
    lvpl = paste0("Z",lvp)
    mu2M = drop(unname(exp(as.matrix(Zm(lvpl,Z$mu))%*%matrix(b$mu[i,]))))
    sg2M = drop(unname(probs(as.matrix(Zm(lvpl,Z$sigma))%*%matrix(b$sigma[i,]))))
    EY2M = (1-sg2M)*mu2M
  }
 # Add other distributions in SimFA::fod
 if(rtF == T){ return(list(mean = EY,  sd = SY, quant = as.data.frame(qM), eM = eM, sM = sM, quantM = qMM, EY2M = EY2M))
 } else return(list(mean = EY,  sd = SY, quant = as.data.frame(qM), eM = eM, sM = sM, quantM = qMM))
}

f2Fun <- function(Y,fam,g,i){
  #For graphics: Y = mod$Y[,item]; fam = mod$fam[item]; g = gFun(mod$gr$out,mod$b,fam[[item]],item); i = cuts
  
 if(fam == "normal"){
  fyz <- dnorm(Y, g$mu[i], g$sigma[i])
 }
 if(fam == "lognormal"){
  fyz <- dlnorm(Y, g$mu[i], g$sigma[i])
 }
 if(fam == "poisson"){
  fyz <- dpois(ceiling(Y), g$mu[i]) 
 }
 if(fam == "gamma"){
  fyz <- dgamma(Y,shape = g$mu[i], scale = g$sigma[i])
 }
 if(fam == "binomial"){
  fyz <- dbinom(round(Y),1,prob = g$mu[i])
 }  
 if(fam == "ZIpoisson"){
  dZIpoisson <- function(Y,mu,sigma,log=T){
   u <- as.numeric(Y == 0)
   lf <- u*c(log(sigma + (1-sigma)*exp(-mu))) + (1-u)*(c(log(1-sigma)) - c(mu) + Y*c(log(mu)) - lfactorial(Y))
   if(log == T) return(lf) else return(exp(lf))
   }
  fyz <- dZIpoisson(ceiling(Y), g$mu[i], g$sigma[i], log = F)
 }
  
  # Add other distributions in SimFA::fod
  return(fyz)
}

p2Fun <- function(fam){ # should be evaluated at i = 1:p; fam = fam[i]
  if(fam == "normal") lims <- c(-Inf,Inf)
  if(fam == "lognormal") lims <- c(0  + .Machine$double.eps, Inf)
  if(fam == "poisson") lims <- c(0,Inf)
  if(fam == "gamma") lims <- c(0  + .Machine$double.eps, Inf)
  if(fam == "binomial") lims <- c(0,1)
  if(fam == "ZIpoisson") lims <- c(0,Inf)
  # Add other parameters for other distributions
  return(lims)
}

