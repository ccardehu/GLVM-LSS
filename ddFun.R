
coefmod <- function(bet,beta,gr,loadmt){ # bet = unlist(borg); beta = borg; gr = gr
  coefM <- NULL
  tmpl <- NULL
  for(i in 1:length(beta)){
    if(i == 1) tmpl[[i]] <- 1:length(beta[[i]])
    else tmpl[[i]] <- seq(from = (length(beta[[i-1]]) + 1), length.out = length(beta[[i]]))
  }
  for(i in 1:length(beta)){
    coefM[[i]] <- matrix(bet[tmpl[[i]]],nrow = p, ncol = ncol(beta[[i]]))*loadmt[[i]]
    colnames(coefM[[i]]) <- colnames(gr$out[[i]]); rownames(coefM[[i]]) <- colnames(Y)
  }
  names(coefM) <- names(beta)
  return(coefM)
}

probs <- function(x){
  pr <- plogis(x)
  if (any(ind <- pr == 1)) pr[ind] <- 1 - sqrt(.Machine$double.eps)
  if (any(ind <- pr == 0)) pr[ind] <- sqrt(.Machine$double.eps)
  return(pr)
}

expit <- function(x) exp(x)/(1+exp(x)); logit <- function(x) log(x/(1-x))

pFun <- function(i,fam){ # should be evaluated at i = 1:p; fam = fam[i]
  if(fam == "normal") pars <- c("mu","sigma")
  if(fam == "poisson") pars <- c("mu")
  if(fam == "gamma") pars <- c("mu", "sigma")
  if(fam == "binom") pars <- c("mu")
  if(fam == "ZIpoisson") pars <- c("mu", "sigma")
  # Add other parameters for other distributions
  return(pars)
}

rFun <- function(n,i,fam,Z,b){ #should be evaluated at i = 1:p; fam = fam[i]; Z = simR$Z; b = borg
  
  if(fam == "normal"){
    mu = as.matrix(Z$mu)%*%matrix(b$mu[i,])
    sigma = exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))
    fyz <- rnorm(n, drop(mu), drop(sigma))
  }
  if(fam == "poisson"){
    mu = exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
    fyz <- rpois(n, drop(mu)) 
  }
  if(fam == "gamma"){
    mu = exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
    sigma = exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))
    fyz <- rgamma(n,shape = mu, scale = sigma)
  }
  if(fam == "binom"){
    mu = probs(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
    fyz <- rbinom(n,1,prob = mu)
  }  
  if(fam == "ZIpoisson"){
    mu = exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
    sigma = probs(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))
    #u <- rbinom(n,1,1-sigma)
    #pois <- rpois(n,mu)
    #fyz <- u * pois
    fyz <- ifelse(rbinom(n,1,1-sigma) == 0, 0, rpois(n,mu))
    #fyz <- gamlss.dist::rZIP(n, drop(mu), drop(sigma))
  }
  
  # Add other distributions in SimFA::fod
  return(fyz)
}

dFun <- function(i,z,fam,Y,Z,b){ #should be evaluated at i = 1:p; z = 10; fam = fam[i]; Y = simR$Y[,i]; Z = gr$out; b = borg
  # if(!is.matrix(Z)) Z <- as.matrix(Z)
  # if(ncol(Z) == 1) Z <- t(Z)
  
  if(fam == "normal"){ # normalHT
    mu = as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])
    sigma = exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,]))
    fyz <- dnorm(Y, drop(mu), drop(sigma), log = T)
  }
  if(fam == "poisson"){
    mu = exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
    fyz <- dpois(Y, drop(mu), log = T) 
  }
  if(fam == "gamma"){
    mu = exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
    sigma = exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,]))
    fyz <- dgamma(Y,shape = drop(mu), scale = drop(sigma), log = T)
  }
  if(fam == "binom"){
    mu = probs(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
    fyz <- dbinom(Y,size = 1,prob = drop(mu), log = T)
  }
  if(fam == "ZIpoisson"){
    mu = exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
    sigma = probs(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,]))
    dZIpoisson <- function(Y,mu,sigma){
      u <- as.numeric(Y == 0)
      lf <- u*c(log(sigma + (1-sigma)*exp(-mu))) + (1-u)*(c(log(1-sigma)) - c(mu) + Y*c(log(mu)) - lfactorial(Y))
      return(lf)
    }
    fyz <- dZIpoisson(Y, drop(mu), drop(sigma))
    #fyz <- gamlss.dist::dZIP(Y, drop(mu), drop(sigma), log = T)
  }
  # Add other distributions
  return(fyz)
}

dvFun <- function(i,z,fam,Y,Z,b){ # i = 1; z = 10; fam = fam[i]; Y = simR$Y[,i]; Z = gr$out; b = borg
  # if(!is.matrix(Z)) Z <- as.matrix(Z)
  # if(ncol(Z) == 1) Z <- t(Z)
  
  if(fam == "normal"){
    mu = drop(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,]))
    sigma = drop(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
    d1mu = (Y - mu)/(sigma^2)  * 1
    d2mu = rep(-(1/sigma^2) * 1, length(Y))
    d1sg = (((Y - mu)^2 - sigma^2)/(sigma^3)) * sigma
    d2sg = rep(-(2/(sigma^2)) * sigma^2, length(Y))
    dcms = rep(0, length(Y)) # the expected cross derivative mu and sigma
    rL <- list("d1mu" = d1mu, "d2mu" = d2mu, "d1sg" = d1sg, "d2sg" = d2sg, "dc" = dcms)
  }
  if(fam == "poisson"){
    mu = exp(drop(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
    d1mu = ((Y - mu)/mu) * mu
    d2mu = rep((-1/mu) * mu^2, length(Y))
    rL <- list("d1mu" = d1mu, "d2mu" = d2mu)
  }
  if(fam == "gamma"){
    mu = drop(exp(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
    sigma = drop(exp(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
    d1mu = (log(Y) - digamma(mu) - log(sigma)) * mu
    d2mu = rep(-trigamma(mu) * mu^2, length(Y))
    d1sg = (Y/(sigma^2) - mu/sigma) * sigma
    d2sg = rep((-mu/(sigma)) * sigma^2, length(Y))
    dcms = rep(0, length(Y)) # the expected cross derivative mu and sigma
    rL <- list("d1mu" = d1mu, "d2mu" = d2mu, "d1sg" = d1sg, "d2sg" = d2sg, "dc" = dcms)
  }
  if(fam == "binom"){
    mu = probs(drop(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
    d1mu = (Y - mu)/(mu * (1 - mu)) * (mu*(1-mu))
    d2mu = rep(-(1/(mu * (1 - mu))) * (mu*(1-mu))^2, length(Y))
    rL <- list("d1mu" = d1mu, "d2mu" = d2mu)
  }
  if(fam == "ZIpoisson"){
    mu = exp(drop(as.matrix(Z$mu[z,])%*%matrix(b$mu[i,])))
    sigma = probs(drop(as.matrix(Z$sigma[z,])%*%matrix(b$sigma[i,])))
    u = ifelse(Y == 0, (1 + exp(-logit(sigma)-mu))^(-1),0) # as.numeric(Y == 0)
    d1mu = (1-u)*(Y/mu-1) * mu
    d2mu = (1-u)*(-1/mu) * mu^2 # expected value is -1/mu; original is (-Y/mu^2)
    d1sg = (u-sigma)/(sigma*(1-sigma)) * (sigma*(1-sigma))
    d2sg = -(u/(sigma^2) + (1-u)/((1-sigma)^2)) * (sigma*(1-sigma))^2
    dcms = rep(0, length(Y))
    rL <- list("d1mu" = d1mu, "d2mu" = d2mu, "d1sg" = d1sg, "d2sg" = d2sg, "dc" = dcms)
  }
  
  return(rL)
}

# Function to compute E(Y) and SD(Y) (or quantiles)

fFun <- function(i,fam,Z,b,qnt = c(0.2,0.4,0.6,0.8)){
  #should be evaluated at i = 1:p; fam = fam[i]; Z = simR$Z; b = borg; qnt = c(0.2,0.4,0.6,0.8)
  
  if(fam == "normal"){
    mu = drop(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
    sigma = drop(exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,])))
    EY = mu
    SY = sigma
    qM = matrix(nrow = length(mu), ncol = length(qnt)); colnames(qM) <- paste0("q",qnt*100)
    for(i in seq_along(qnt)){ qM[,i] <- qnorm(qnt[i],mu,sigma) }
  }
  if(fam == "poisson"){
    mu = drop(exp(as.matrix(Z$mu)%*%matrix(b$mu[i,])))
    EY = mu
    SY = sqrt(mu)
    qM = matrix(nrow = length(mu), ncol = length(qnt)); colnames(qM) <- paste0("q",qnt*100)
    for(i in seq_along(qnt)){ qM[,i] <- qpois(qnt[i],mu) }
  }
  if(fam == "gamma"){
    mu = drop(exp(as.matrix(Z$mu)%*%matrix(b$mu[i,])))
    sigma = drop(exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,])))
    EY = mu*sigma
    SY = sqrt(mu*sigma^2)
    qM = matrix(nrow = length(mu), ncol = length(qnt)); colnames(qM) <- paste0("q",qnt*100)
    for(i in seq_along(qnt)){ qM[,i] <- qgamma(qnt[i], scale = mu, scale = sigma) }
  }
  if(fam == "binom"){
    mu = drop(probs(as.matrix(Z$mu)%*%matrix(b$mu[i,])))
    EY = mu
    SY = sqrt(mu*(1-mu))
    qM = matrix(nrow = length(mu), ncol = length(qnt)); colnames(qM) <- paste0("q",qnt*100)
    for(i in seq_along(qnt)){ qM[,i] <- qbinom(qnt[i], 1, mu) }
  }
  if(fam == "ZIpoisson"){
    mu = drop(exp(as.matrix(Z$mu)%*%matrix(b$mu[i,])))
    sigma = drop(probs(as.matrix(Z$sigma)%*%matrix(b$sigma[i,])))
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
  }
  # Add other distributions in SimFA::fod
  return(list(mean = EY,  sd = SY, quant = as.data.frame(qM)))
}

