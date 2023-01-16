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
      dvy$sg = sapply(1:qp, function(r){ rep(1/(sg[r]*(sg[r]-1)) * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2 , length(y)) } )
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

Lognormal <- function(mu.link = "identity", sg.link = "log"){
  
  sta.mu <- make.link(mu.link)
  sta.sg <- make.link(sg.link)
  
  dy <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    sapply(1:nrow(ghQ$points), function(r) dlnorm(y,mu[r],sg[r],log = T))
  }
  
  dv1y <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    y = log(y)
    qp = length(mu)
    dvy <- list(mu = NULL, sg = NULL)
    dvy$mu = sapply(1:qp, function(r) (y-mu[r])/(sg[r]^2) ) #  * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))
    dvy$sg = sapply(1:qp, function(r) (((y - mu[r])^2 - sg[r]^2)/(sg[r]^3)) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])))
    return(dvy)
  }
  
  dv2y <- function(i,y,b,ghQ,info,dv1Y = NULL){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    y = log(y)
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
    y = log(y)
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
    rlnorm(n,mu,sg)
  }
  
  structure(list(family = "Lognormal", npar = 2, pars = c("mu","sigma"),
                 iuse = quote(LOGNO()), dY = dy, sf = sfun,
                 dv1Y = dv1y, dv2Y = dv2y, dvCY = dvCy,
                 link.mu = sta.mu$name, linkf.mu = sta.mu$linkfun,
                 link.sg = sta.sg$name, linkf.sg = sta.sg$linkfun),
            class = "dist_glvmlss")
}

Beta0 <- function(mu.link = "logit", sg.link = "log"){

  sta.mu <- make.link(mu.link)
  sta.sg <- make.link(sg.link)

  dy <- function(i,y,b,ghQ){
    y <- y.(y)
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    alp = mu*sg
    bet = (1-mu)*sg
    sapply(1:nrow(ghQ$points), function(r) dbeta(y, shape1 = alp[r], shape2 = bet[r],log = T))
  }

  dv1y <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    alp = mu*sg
    bet = (1-mu)*sg
    y <- y.(y)
    dvy <- list(mu = NULL, sg = NULL)
    dvy$mu = sapply(1:qp, function(r){ ( -sg[r]*(2*atanh(1-2*y) + digamma(alp[r]) + digamma(bet[r])) ) * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) })
    dvy$sg = sapply(1:qp, function(r){ ( -2*mu[r]*atanh(1-2*y) + digamma(sg[r]) - mu[r]*digamma(alp[r]) - (1-mu[r])*digamma(bet[r]) ) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])) } )
    return(dvy)
  }

  dv2y <- function(i,y,b,ghQ,info, dv1Y = NULL){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    alp = mu*sg
    bet = (1-mu)*sg
    y <- y.(y)
    dvy <- list(mu = NULL, sg = NULL)
    if(info == "Fisher"){
      dvy$mu = sapply(1:qp, function(r){ rep((-sg[r]^2 * (-trigamma(alp[r]) + trigamma(bet[r]))) * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2, length(y)) } )
      dvy$sg = sapply(1:qp, function(r){ rep(( trigamma(sg[r]) - mu[r]^2*trigamma(alp[r]) - (1-mu[r])^2*trigamma(bet[r])) * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2, length(y)) } )
    } else {
      dvy$mu = sapply(1:qp, function(r){( rep(-sg[r]^2 * (-trigamma(alp[r]) + trigamma(bet[r])) * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2, length(y))
                                          + numDeriv::grad(function(x) sta.mu$mu.eta(sta.mu$linkfun(pmin(pmax(x, sqrt(.Machine$double.eps)), 1-sqrt(.Machine$double.eps)))), mu[r])*dv1Y$mu[,r] ) } )
      dvy$sg = sapply(1:qp, function(r){( rep(( trigamma(sg[r]) - mu[r]^2*trigamma(alp[r]) - (1-mu[r])^2*trigamma(bet[r])) * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2, length(y)) + dv1Y$sg[,r] ) } ) }
    return(dvy)
  }

  dvCy <- function(i,y,b,ghQ,info){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    alp = mu*sg
    bet = (1-mu)*sg
    y <- y.(y)
    dvy <- list(mu = list(sg = NULL))
    dvy$mu$sg = sapply(1:qp, function(r){ ( -2*atanh(1-2*y) - digamma(alp[r]) - digamma(bet[r]) - alp[r]*trigamma(alp[r])
                                            + bet[r]*trigamma(bet[r]) ) * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])) } )
    return(dvy)
  }

  sfun <- function(i,n,b,Z){
    mu = sta.mu$linkinv(c(as.matrix(Z$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(Z$sig)%*%b$sig[i,]))
    alp = mu*sg
    bet = (1-mu)*sg
    rbeta(n, shape1 = alp, shape2 = bet)
  }

  structure(list(family = "Beta0", npar = 2, pars = c("mu","sigma"),
                 iuse = quote(BEo()), dY = dy, sf = sfun,
                 dv1Y = dv1y, dv2Y = dv2y, dvCY = dvCy,
                 link.mu = sta.mu$name, linkf.mu = sta.mu$linkfun,
                 link.sg = sta.sg$name, linkf.sg = sta.sg$linkfun),
            class = "dist_glvmlss")
}

SkewNormal <- function(mu.link = "identity", sg.link = "log", nu.link = "identity"){
  
  sta.mu <- make.link(mu.link)
  sta.sg <- make.link(sg.link)
  sta.nu <- make.link(nu.link)
  
  dy <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    nu = sta.nu$linkinv(c(as.matrix(ghQ$out$nu)%*%b$nu[i,]))
    dSKN <- Vectorize( function(Y,mu,sigma,nu,log = T){
      tryCatch({z <- (Y - mu)/sigma
      w <- nu*z
      z2 <- (z^2)/2
      w2 <- (w^2)/2
      lpdf <- -0.5*log(2*pi) - z2
      lcdf_fun <- function(idx,w2.,w.){
        w2. <- w2.[idx]
        w. <- w.[idx]
        if(w2. == 0){
          lcdf <- lpnorm(0, log.p = T)
        } else lcdf <- log(0.5 * (1 + pgamma(w2., shape = 1/2, scale = 1) * sign(w.))) #  + .Machine$double.eps^10 
        return(lcdf) }
      lcdf <- sapply(1:length(w2), lcdf_fun, w2. = w2, w. = w)
      lf <- lpdf + lcdf + log(2) - log(sigma)
      if(log == T) return(lf) else return(exp(lf))}, error = function(e) NA) } )
    sapply(1:nrow(ghQ$points), function(r) dSKN(y,mu[r],sg[r],nu[r],log = T))
  }
  
  eta1 = Vectorize( function(x) exp( dnorm(x,log=T) - pnorm(x,log.p = T) ) )
  eta2 = Vectorize( function(x) -x*eta1(x) - eta1(x)^2 )
  d1mu = function(y, mu, sigma, nu){
    z <- (y - mu)/sigma
    w <- nu * z
    dldm <- z/sigma - (nu/sigma)*eta1(w)
    dldm
  }
  d1sg = function(y, mu, sigma, nu){
    z <- (y - mu)/sigma
    w <- nu * z
    dldd <- -eta1(w)*w/sigma + (z^2 - 1)/sigma
    dldd
  }
  d1nu = function(y, mu, sigma, nu){
    z <- (y - mu)/sigma
    w <- nu * z
    dldv <- eta1(w) * z
    dldv
  }
  d2mu = function(y, mu, sigma, nu){
    z <- (y - mu)/sigma
    w <- nu * z
    d2ldm2 <- -1/sigma^2 - (nu/sigma)^2*eta2(w)
    d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
    d2ldm2
  }
  d2sg = function(y, mu, sigma, nu){
    z <- (y - mu)/sigma
    w <- nu * z
    d2ldd2 <- sigma^(-2)*(-1 -3*z^2 -2*w*eta1(w) - w^2*eta2(w))
    d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
    d2ldd2
  }
  d2nu = function(y, mu, sigma, nu){
    z <- (y - mu)/sigma
    w <- nu * z
    d2ldv2 <- eta2(w)*z^2
    d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2, -1e-15)
    d2ldv2
  }
  dcms = function(y, mu, sigma, nu){
    z <- (y - mu)/sigma
    w <- nu * z
    d2ldmdd <- 1/sigma^2*(-2*z + nu*eta1(w) + nu^2*eta2(w)*z)
    d2ldmdd
  }
  dcmn = function(y, mu, sigma, nu){
    z <- (y - mu)/sigma
    w <- nu * z
    d2ldmdv <- -1/sigma*eta1(w) - nu/sigma*eta2(w)*z
    d2ldmdv
  }
  dcsn = function(y, mu, sigma, nu){
    z <- (y - mu)/sigma
    w <- nu * z
    d2ldddv <- -eta1(w)*z/sigma -eta2(w)*w*z/sigma
    d2ldddv
  }
  
  dv1y <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    nu = sta.nu$linkinv(c(as.matrix(ghQ$out$nu)%*%b$nu[i,]))
    qp = length(mu)
    dvy <- list(mu = NULL, sg = NULL, nu = NULL)
    dvy$mu = sapply(1:qp, function(r) d1mu(y,mu[r],sg[r],nu[r])) #  * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) # Identity
    dvy$sg = sapply(1:qp, function(r) d1sg(y,mu[r],sg[r],nu[r]) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])))
    dvy$nu = sapply(1:qp, function(r) d1nu(y,mu[r],sg[r],nu[r])) # * sta.nu$nu.eta(sta.nu$linkfun(nu[r])) # Identity
    return(dvy)
  }
  
  dv2y <- function(i,y,b,ghQ,info,dv1Y = NULL){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    nu = sta.nu$linkinv(c(as.matrix(ghQ$out$nu)%*%b$nu[i,]))
    qp = length(mu)
    dvy <- list(mu = NULL, sg = NULL, nu = NULL)
    if(info == "Fisher"){
      dvy$mu = sapply(1:qp, function(r) d2mu(y,mu[r],sg[r],nu[r]) ) # * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2
      dvy$sg = sapply(1:qp, function(r) d2sg(y,mu[r],sg[r],nu[r]) * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2)
      dvy$nu = sapply(1:qp, function(r) d2nu(y,mu[r],sg[r],nu[r]) ) # * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2
    } else {
      dvy$mu = sapply(1:qp, function(r) d2mu(y,mu[r],sg[r],nu[r]) ) # * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2
      dvy$sg = sapply(1:qp, function(r) d2sg(y,mu[r],sg[r],nu[r]) * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2 + dv1Y$sg[,r])
      dvy$nu = sapply(1:qp, function(r) d2nu(y,mu[r],sg[r],nu[r]) ) # * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2
    }
    return(dvy)
  }
  
  dvCy <- function(i,y,b,ghQ,info){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    nu = sta.nu$linkinv(c(as.matrix(ghQ$out$nu)%*%b$nu[i,]))
    qp = length(mu)
    dvy <- list(mu = list(sg = NULL, nu = NULL), sg = list(nu = NULL))
    if(info == "Fisher"){
      dvy$mu$sg = sapply(1:qp, function(r) dcms(y,mu[r],sg[r],nu[r]) )
      dvy$mu$nu = sapply(1:qp, function(r) dcmn(y,mu[r],sg[r],nu[r]) )
      dvy$sg$nu = sapply(1:qp, function(r) dcsn(y,mu[r],sg[r],nu[r]) )
    } else {
      dvy$mu$sg = sapply(1:qp, function(r){ dcms(y,mu[r],sg[r],nu[r]) * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])) })
      dvy$mu$nu = sapply(1:qp, function(r){ dcmn(y,mu[r],sg[r],nu[r]) * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) * sta.nu$mu.eta(sta.nu$linkfun(nu[r])) })
      dvy$sg$nu = sapply(1:qp, function(r){ dcsn(y,mu[r],sg[r],nu[r]) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])) * sta.nu$mu.eta(sta.nu$linkfun(nu[r])) })
    }
    return(dvy)
  }
  
  sfun <- function(i,n,b,Z){
    mu = sta.mu$linkinv(c(as.matrix(Z$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(Z$sig)%*%b$sig[i,]))
    nu = sta.sg$linkinv(c(as.matrix(Z$nu)%*%b$nu[i,]))
    gamlss.dist::rSN1(n,mu,sg,nu)
  }
  
  structure(list(family = "SkewNormal", npar = 3, pars = c("mu","sigma","nu"),
                 iuse = quote(SN1()), dY = dy, sf = sfun,
                 dv1Y = dv1y, dv2Y = dv2y, dvCY = dvCy,
                 link.mu = sta.mu$name, linkf.mu = sta.mu$linkfun,
                 link.sg = sta.sg$name, linkf.sg = sta.sg$linkfun,
                 link.nu = sta.nu$name, linkf.nu = sta.nu$linkfun),
            class = "dist_glvmlss")
}