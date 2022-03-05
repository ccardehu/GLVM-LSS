Normal <- function(mu.link = "identity", sg.link = "log"){
  
  sta.mu <- make.link(mu.link)
  sta.sg <- make.link(sg.link)
  
  dy <- function(i,y,b,ghQ){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    sapply(1:nrow(ghQ$points), function(r) dnorm(y,mu[r],sg[r],T))
  }
  
  dvy <- function(i,y,b,ghQ,info){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(ghQ$out$sig)%*%b$sig[i,]))
    qp = length(mu)
    dvy <- list(d1 = NULL, d2 = NULL, dc = list(mu = NULL))
    dvy$d1$mu = sapply(1:qp, function(r) (y-mu[r])/(sg[r]^2) ) #  * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))
    dvy$d2$mu = sapply(1:qp, function(r) rep(-1/sg[r]^2 , length(y)) ) # * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2
    dvy$d1$sg = sapply(1:qp, function(r) (((y - mu[r])^2 - sg[r]^2)/(sg[r]^3)) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])))
    if(info == "Fisher"){
      dvy$d2$sg = sapply(1:qp, function(r) rep(-2/(sg[r]^2) * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2 , length(y)) )
      dvy$dc$mu$sg = sapply(1:qp, function(r) rep(0, length(y)) )
    } else {
      dvy$d2$sg = sapply(1:qp, function(r){ (-3*(y-mu[r])^2/sg[r]^4 + 1/sg[r]^2) * sta.sg$mu.eta(sta.sg$linkfun(sg[r]))^2 + dvy$d1$sg[,r] })
      dvy$dc$mu$sg = sapply(1:qp, function(r){ -2*(y-mu[r])/sg[r]^3 * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) * sta.sg$mu.eta(sta.sg$linkfun(sg[r])) })
    }
    return(dvy)
  }
  
  sfun <- function(i,n,b,Z){
    mu = sta.mu$linkinv(c(as.matrix(Z$mu)%*%b$mu[i,]))
    sg = sta.sg$linkinv(c(as.matrix(Z$sig)%*%b$sig[i,]))
    rnorm(n,mu,sg)
  }
  
  structure(list(family = "Normal", npar = 2, pars = c("mu","sigma"),
                 iuse = quote(NO()),
                 dY = dy, dvY = dvy, sf = sfun,
                 link.mu = sta.mu$name, linkf.mu = sta.mu$linkfun,
                 link.sg = sta.sg$name, linkf.sg = sta.sg$linkfun),
            class = "dist_glvmlss")
  
}

Binomial <- function(n = 1, mu.link = "logit"){
  
  sta.mu <- make.link(mu.link)
  
  dy <- function(i,y,b,ghQ){
    mu = c(sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,])))
    sapply(1:nrow(ghQ$points), function(r) dbinom(y,n,mu[r],T))
  }
  
  dvy <- function(i,y,b,ghQ,info){
    mu = sta.mu$linkinv(c(as.matrix(ghQ$out$mu)%*%b$mu[i,]))
    qp = length(mu)
    dvy <- list(d1 = NULL, d2 = NULL)
    dvy$d1$mu = sapply(1:qp, function(r){ (y - mu[r])/(mu[r] * (1 - mu[r])) * sta.mu$mu.eta(sta.mu$linkfun(mu[r])) } )
    if(info == "Fisher"){
      dvy$d2$mu = sapply(1:qp, function(r) rep(-(1/(mu[r] * (1 - mu[r]))) * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2, length(y)))
    } else {
      dvy$d2$mu = sapply(1:qp, function(r){ ( -(mu[r]^2 - 2*y*mu[r] + y)/((mu[r]-1)*mu[r])^2 * sta.mu$mu.eta(sta.mu$linkfun(mu[r]))^2 
      + numDeriv::grad(function(x) sta.mu$mu.eta(sta.mu$linkfun(pmin(pmax(x,0),1))), mu[r])*dvy$d1$mu[,r] ) } ) 
    }
    return(dvy)
  }
  
  sfun <- function(i,n.,b,Z){
    mu = sta.mu$linkinv(c(as.matrix(Z$mu)%*%b$mu[i,]))
    rbinom(n.,n,mu)
  }
  
  structure(list(family = "Bernoulli", npar = 1, pars = c("mu"),
                 iuse = quote(BI()),
                 dY = dy, dvY = dvy, sf = sfun,
                 link.mu = sta.mu$name, linkf.mu = sta.mu$linkfun),
                 class = "dist_glvmlss")
}

