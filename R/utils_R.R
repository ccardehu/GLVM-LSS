fam4bstart <- function(family){
  out <- vector(mode = "list",length = length(family))
  for(i in 1:length(out)){
    if(family[i] == "normal") out[[i]] <- quote(NO())
    if(family[i] == "binomial") out[[i]] <- quote(BI())
    if(family[i] == "beta") out[[i]] <- quote(BE())
    if(family[i] == "sn") out[[i]] <- quote(SN1())
  }
  return(out)
}

b_ini <- function(Y,family,form){

  Y <- Y[complete.cases(Y),]
  lvar <- unique(unlist(lapply(1:length(form), function(i) all.vars(form[[i]]))))
  lvar <- grep("Z", lvar, fixed = T, value = T)
  q <- length(lvar)
  p <- ncol(Y)
  Z <- scale(princomp(Y,cor = T)$scores)[,1:q, drop = F]
  colnames(Z) <- paste0("Z", 1:q)

  bstart <- array(0,dim = c(p,q+1,length(form)))
  fstart <- fam4bstart(family)
  colnames(bstart) <- c("(Intercept)", lvar)
  rownames(bstart) <- colnames(Y)
  dimnames(bstart)[[3]] <- names(form)
  eq <- form

  for(i in 1:ncol(Y)){
    tmpY <- Y[,i]
    if(family[i] == "beta") tmpY <- fixy(tmpY)
    for(ii in 1:length(eq)){ eq[[ii]] <- update(eq[[ii]], tmpY ~ .)  }
    tmp <- try(gamlss::gamlss(eq$mu, sigma.formula = eq$sg, nu.formula = eq$nu, tau.formula = eq$ta,
                              family = fstart[[i]], data = as.data.frame(cbind(tmpY,Z)),
                              control = gamlss::gamlss.control(trace = F)), silent = T)

    if(!inherits(tmp, "try-error")){
      tmpB <- gamlss::coefAll(tmp)
      K <- length(tmpB)
      for(k in 1:K){
        bstart[i,names(tmpB[[k]]),k] <- tmpB[[k]]
      }
    } else {
      bstart[i,,] <- 0.01
    }
  }
  return(bstart)
}

fixy <- function(x){
  if (any(x == 1, na.rm = T)) x[x == 1] <- 1 - sqrt(.Machine$double.eps)
  if (any(x == 0, na.rm = T)) x[x == 0] <- sqrt(.Machine$double.eps)
  return(x)
}

genrb <- function(b){
  rb <- ifelse(b == 0 | b == 1, b, NA)
  return(rb)
}

genrR <- function(R){
  rR <- ifelse(R == 0 | R == 1, R, NA)
  return(rR)
}

rb2b <- function(b,rb){
  b[!is.na(rb)] <- rb[!is.na(rb)]
  return(list(b = b, rb = rb, irb = sum(!is.na(rb[,,1]))))
}

rR2R <- function(R,rR){
  if(!isSymmetric(rR)){
    rR[upper.tri(rR)] <- t(rR)[upper.tri(rR)]
  }
  R[!is.na(rR)] <- rR[!is.na(rR)]
  return(list(R = R, rR = rR, irR = sum(!is.na(rR[lower.tri(rR,diag = T)]))))
}

is.diag <- function(m){
  all(m[lower.tri(m)] == 0, m[upper.tri(m)] == 0)
}
