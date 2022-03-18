Normal <- function(mu.link = "identity", sg.link = "log"){
  
  sta.mu <- make.link(mu.link)
  sta.sg <- make.link(sg.link)
  
  dy <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    sapply(1:nrow(ghQ$points), function(r) dnorm(y,mu[r],sg[r],log = T))
  }
  
  dv1y <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    dvy <- list(mu = NULL, sg = NULL)
    dvy$mu = sapply(1:qp, function(r) (y-mu[r])/(sg[r]^2) ) #  * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))
    dvy$sg = sapply(1:qp, function(r) (((y - mu[r])^2 - sg[r]^2)/(sg[r]^3)) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])))
    return(dvy)
  }
  
  dv2y <- function(i,y,b,ghQ,info,dv1Y = NULL){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    dvy <- list(mu = NULL, sg = NULL)
    dvy$mu = sapply(1:qp, function(r) rep(-1/sg[r]^2 , length(y)) ) # * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2
    if(info == "Fisher"){
      dvy$sg = sapply(1:qp, function(r) rep(-2/(sg[r]^2) * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2 , length(y)) )
    } else {
      dvy$sg = sapply(1:qp, function(r){ (-3*(y-mu[r])^2/sg[r]^4 + 1/sg[r]^2) * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2 + dv1Y$sg[,r] })
    }
    return(dvy)
  }
  
  dvCy <- function(i,y,b,ghQ,info){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    dvy <- list(mu = list(sg = NULL))
    if(info == "Fisher"){
      dvy$mu$sg = sapply(1:qp, function(r) rep(0, length(y)) )
    } else {
      dvy$mu$sg = sapply(1:qp, function(r){ -2*(y-mu[r])/sg[r]^3 * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])) })
    }
    return(dvy)
  }
  
  sfun <- function(i,n,b,Z){
    mu = sta.mu$linkinv(c(as.matrix(Z$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(Z$sig)%*%b$sig[i,]))
    rnorm(n,mu,sg)
  }
  
  structure(list(family = "Normal", npar = 2, pars = c("mu","sigma"),
                 iuse = quote(NO()), dY = dy, sf = sfun,
                 dv1Y = dv1y, dv2Y = dv2y, dvCY = dvCy,
                 link.mu = sta.mu$name, linkf.mu = sta.mu$linkfun,
                 link.sg = sta.sg$name, linkf.sg = sta.sg$linkfun),
            class = "dist_glvmlss")
}

Binomial <- function(n = 1, mu.link = "logit"){
  
  sta.mu <- make.link(mu.link)
  
  dy <- function(i,y,b,ghQ){
    mu = c(sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,])))
    sapply(1:nrow(ghQ$points), function(r) dbinom(y,n,mu[r], log = T))
  }
  
  dv1y <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    qp = length(mu)
    dvy <- list(mu = NULL)
    dvy$mu = sapply(1:qp, function(r){ (y - mu[r])/(mu[r] * (1 - mu[r])) * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) } )
    return(dvy)
  }
  
  dv2y <- function(i,y,b,ghQ,info, dv1Y = NULL){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    qp = length(mu)
    dvy <- list(mu = NULL)
    if(info == "Fisher"){
      dvy$mu = sapply(1:qp, function(r) rep(-(1/(mu[r] * (1 - mu[r]))) * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2, length(y)))
    } else {
      dvy$mu = sapply(1:qp, function(r){ ( -(mu[r]^2 - 2*y*mu[r] + y)/((mu[r]-1)*mu[r])^2 * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2 
                                              + numDeriv::grad(function(x) sta.mu$mu.eta(sta.mu$linkfun(pmin(pmax(x,0),1))), mu[r])*dv1Y$mu[,r] ) } ) 
    }
    return(dvy)
  }
  
  dvCy <- function(i,y,b,ghQ,info){ NULL }
  
  sfun <- function(i,n.,b,Z){
    mu = sta.mu$linkinv(c(as.matrix(Z$mu)%*%b$mu[i,]))
    rbinom(n.,n,mu)
  }
  
  structure(list(family = "Binomial", npar = 1, pars = c("mu"),
                 iuse = quote(BI()), dY = dy, sf = sfun,
                 dv1Y = dv1y, dv2Y = dv2y, dvCY = dvCy,
                 link.mu = sta.mu$name, linkf.mu = sta.mu$linkfun),
            class = "dist_glvmlss")
}

ZIpoisson <- function(mu.link = "log", sg.link = "logit"){
  
  sta.mu <- make.link(mu.link)
  sta.sg <- make.link(sg.link)
  
  dy <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    dZIpoisson <- function(Y,mu,sigma,log = T){
      u <- as.numeric(Y == 0)
      lf <- u*c(log(sigma + (1-sigma)*exp(-mu))) + (1-u)*(c(log(1-sigma)) - c(mu) + Y*c(log(mu)) - lfactorial(Y))
      if(log == T) return(lf) else return(exp(lf)) }
    sapply(1:nrow(ghQ$points), function(r) dZIpoisson(y,mu[r],sg[r], log = T))
  }

  dv1y <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    u <- function(y,mu.,sg.){ifelse(y == 0, (1 + exp(-sta.sg$linkfun(sg.)-mu.))^(-1),0)}
    dvy <- list(mu = NULL, sg = NULL)
    dvy$mu = sapply(1:qp, function(r){ (1-u(y,mu[r],sg[r]))*(y/mu[r]-1) * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) } )
    dvy$sg = sapply(1:qp, function(r){ (u(y,mu[r],sg[r])-sg[r])/(sg[r]*(1-sg[r])) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])) } )  
    return(dvy)
  }
  
  dv2y <- function(i,y,b,ghQ,info, dv1Y = NULL){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    u <- function(y,mu.,sg.){ifelse(y == 0, (1 + exp(-sta.sg$linkfun(sg.)-mu.))^(-1),0)}
    dvy <- list(mu = NULL, sg = NULL)
    if(info == "Fisher"){
      dvy$mu = sapply(1:qp, function(r){ (1-u(y,mu[r],sg[r])) * (-1/mu[r]) * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2 } )
      dvy$sg = sapply(1:qp, function(r){ 1/(sg[r]*(sg[r]-1)) * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2 } )
    } else {
      dvy$mu = sapply(1:qp, function(r){ ( (1-u(y,mu[r],sg[r]))*(-y/mu[r]^2) * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2
                                            + dv1Y$mu[,r] ) } )
      dvy$sg = sapply(1:qp, function(r){ ( -(sg[r]^2 - 2*u(y,mu[r],sg[r])*sg[r] + u(y,mu[r],sg[r]))/((sg[r]-1)*sg[r])^2 * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2
                                            + numDeriv::grad(function(x) sta.sg$mu.eta(sta.sg$linkfun(pmin(pmax(x,0),1))), sg[r])*dv1Y$sg[,r] ) } ) }  
    return(dvy)
  }
  
  dvCy <- function(i,y,b,ghQ,info){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    u <- function(y,mu.,sg.){ifelse(y == 0, (1 + exp(-sta.sg$linkfun(sg.)-mu.))^(-1),0)}
    dvy <- list(mu = list(sg = NULL))
    dvy$mu$sg = sapply(1:qp, function(r) rep(0, length(y)) )
    return(dvy)
  }
  
  sfun <- function(i,n,b,Z){
    mu = sta.mu$linkinv(c(as.matrix(Z$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(Z$sig)%*%b$sig[i,]))
    ifelse(rbinom(n,1,1-sg) == 0, 0, rpois(n,mu))
  }
  
  structure(list(family = "ZI-poisson", npar = 2, pars = c("mu","sigma"),
                 iuse = quote(ZIP()), dY = dy, sf = sfun,
                 dv1Y = dv1y, dv2Y = dv2y, dvCY = dvCy,
                 link.mu = sta.mu$name, linkf.mu = sta.mu$linkfun,
                 link.sg = sta.sg$name, linkf.sg = sta.sg$linkfun),
            class = "dist_glvmlss")
}

Beta <- function(mu.link = "logit", sg.link = "logit"){
  
  sta.mu <- make.link(mu.link)
  sta.sg <- make.link(sg.link)
  
  dy <- function(i,y,b,ghQ){
    y <- y.(y)
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    alp = mu*(1-sg^2)/sg^2
    bet = (1-mu)*(1-sg^2)/sg^2
    sapply(1:nrow(ghQ$points), function(r) dbeta(y, shape1 = alp[r],shape2 = bet[r],log = T))
  }
  
  dv1y <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    alp = mu*(1-sg^2)/sg^2
    bet = (1-mu)*(1-sg^2)/sg^2
    y <- y.(y)
    dvy <- list(mu = NULL, sg = NULL)
    dvy$mu = sapply(1:qp, function(r){ (((1 - sg[r]^2)/(sg[r]^2))*(-digamma(alp[r]) + digamma(bet[r]) + log(y) - log(1 - y))) * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) })
    dvy$sg = sapply(1:qp, function(r){ -(2/(sg[r]^3)) * (mu[r] * (-digamma(alp[r]) + digamma(alp[r] + bet[r]) + log(y)) + (1 - mu[r]) * (-digamma(bet[r]) + digamma(alp[r] + bet[r]) + log(1 - y))) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])) } )
    return(dvy)
  }
  
  dv2y <- function(i,y,b,ghQ,info, dv1Y = NULL){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    alp = mu*(1-sg^2)/sg^2
    bet = (1-mu)*(1-sg^2)/sg^2
    y <- y.(y)
    dvy <- list(mu = NULL, sg = NULL)
    if(info == "Fisher"){
      dvy$mu = sapply(1:qp, function(r){ rep(-(((1 - sg[r]^2)^2)/(sg[r]^4)) * (trigamma(alp[r]) + trigamma(bet[r])) * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2, length(y)) } )
      dvy$sg = sapply(1:qp, function(r){ rep(-(4/(sg[r]^6)) * ((mu[r]^2) * trigamma(alp[r]) + ((1 - mu[r])^2) * trigamma(bet[r]) - trigamma(alp[r] + bet[r])) * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2, length(y)) } )
    } else {
      dvy$mu = sapply(1:qp, function(r){( rep(-(((1 - sg[r]^2)^2)/(sg[r]^4)) * (trigamma(alp[r]) + trigamma(bet[r])) * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2, length(y)) 
                                             + numDeriv::grad(function(x) sta.mu$mu.eta(sta.mu$linkfun(pmin(pmax(x, sqrt(.Machine$double.eps)), 1-sqrt(.Machine$double.eps)))), mu[r])*dv1Y$mu[,r] ) } )
      dvy$sg = sapply(1:qp, function(r){( rep(-(4/(sg[r]^6)) * ((mu[r]^2) * trigamma(alp[r]) + ((1 - mu[r])^2) * trigamma(bet[r]) - trigamma(alp[r] + bet[r])) * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2, length(y))
                                             + numDeriv::grad(function(x) sta.sg$mu.eta(sta.sg$linkfun(pmin(pmax(x, sqrt(.Machine$double.eps)), 1-sqrt(.Machine$double.eps)))), sg[r])*dv1Y$sg[,r] ) } ) }
    return(dvy)
  }
  
  dvCy <- function(i,y,b,ghQ,info){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    alp = mu*(1-sg^2)/sg^2
    bet = (1-mu)*(1-sg^2)/sg^2
    y <- y.(y)
    dvy <- list(mu = list(sg = NULL))
    dvy$mu$sg = sapply(1:qp, function(r){( rep((2 * (1 - sg[r]^2)/(sg[r]^5)) * (mu[r] * trigamma(alp[r]) - (1 - mu[r]) * trigamma(bet[r])) * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])), length(y)) ) } )
    return(dvy)
  }
  
  sfun <- function(i,n,b,Z){
    mu = sta.mu$linkinv(c(as.matrix(Z$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(Z$sig)%*%b$sig[i,]))
    alp = mu*(1-sg^2)/sg^2
    bet = (1-mu)*(1-sg^2)/sg^2
    rbeta(n, shape1 = alp, shape2 = bet)
  }
  
  structure(list(family = "Beta", npar = 2, pars = c("mu","sigma"),
                 iuse = quote(BE()), dY = dy, sf = sfun,
                 dv1Y = dv1y, dv2Y = dv2y, dvCY = dvCy,
                 link.mu = sta.mu$name, linkf.mu = sta.mu$linkfun,
                 link.sg = sta.sg$name, linkf.sg = sta.sg$linkfun),
            class = "dist_glvmlss")
}

Poisson <- function(mu.link = "log"){
  
  sta.mu <- make.link(mu.link)
  
  dy <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sapply(1:nrow(ghQ$points), function(r) dpois(y,mu[r],log = T))
  }
  
  dv1y <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    qp = length(mu)
    dvy <- list(mu = NULL)
    dvy$mu = sapply(1:qp, function(r){ (y - mu[r])/mu[r] * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) } )
    return(dvy)
  }
  
  dv2y <- function(i,y,b,ghQ,info, dv1Y = NULL){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    qp = length(mu)
    dvy <- list(mu = NULL)
    if(info == "Fisher"){
      dvy$mu = sapply(1:qp, function(r){ rep(-1/mu[r] * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2, length(y)) } )
    } else {
      dvy$mu = sapply(1:qp, function(r){ -y/mu[r] * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2 + dv1Y$mu[,r] } )
    }
    return(dvy)
  }
  
  dvCy <- function(i,y,b,ghQ,info){ NULL }
  
  sfun <- function(i,n,b,Z){
    mu = sta.mu$linkinv(c(as.matrix(Z$mu)%*%b$mu[i,]))
    rpois(n,mu)
  }
  
  structure(list(family = "Poisson", npar = 1, pars = c("mu"),
                 iuse = quote(PO()), dY = dy, sf = sfun,
                 dv1Y = dv1y, dv2Y = dv2y, dvCY = dvCy,
                 link.mu = sta.mu$name, linkf.mu = sta.mu$linkfun),
            class = "dist_glvmlss")
}