y. <- function(x){
 if (any(x == 1)) x[x == 1] <- 1 - .Machine$double.eps^(1/3)
 if (any(x == 0)) x[x == 0] <- .Machine$double.eps^1/3
 return(x)
}

fyz <- function(Y,ghQ,b,famL){
  Y <- Y$Y
  A1 <- lapply(1:length(famL), function(i) famL[[i]]$dY(i,Y[,i],b,ghQ))
  A1 <- array(unlist(A1),dim = c(nrow(Y),nrow(ghQ$points),length(famL)))
  A2 <- c(exp(rowSums(A1,dim = 2,na.rm = T))%*%ghQ$weights)
  return(list(pD = exp(rowSums(A1,dim = 2, na.rm = T))/A2, ll = sum(log(A2))))
}

ll <- function(cb,Y,ghQ,bg,famL,info,rb){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  f1 <- -d1ll(Y,ghQ,b,famL,info,f0$pD,rb)
  f2 <- -d2ll(Y,ghQ,b,famL,info,f0$pD,rb)
  return(list(value = -f0$ll, gradient = f1, hessian = f2))
}

lazyll <- function(cb,Y,ghQ,bg,famL,info,rb){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  f1 <- -d1ll(Y,ghQ,b,famL,info,f0$pD,rb)
  f2 <- ad2ll(Y,ghQ,b,famL,info,f0$pD,rb)
  return(list(value = -f0$ll, gradient = f1, hessian = f2))
}

pll <- function(cb,Y,ghQ,bg,famL,info,rb,pen.control){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  pMat <- pM(b,rb, penalty = pen.control$penalty,
             lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
  f1 <- -d1ll(Y,ghQ,b,famL,info,f0$pD,rb) + c(nrow(Y$Y)*crossprod(cb,pMat$full))
  f2 <- -d2ll(Y,ghQ,b,famL,info,f0$pD,rb) + nrow(Y$Y)*pMat$full 
  return(list(value = -f0$ll + 0.5*nrow(Y$Y)*crossprod(cb,pMat$full)%*%cb,
              gradient = f1, hessian = f2))
}

lazypll <- function(cb,Y,ghQ,bg,famL,info,rb,pen.control){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  pMat <- pM(b,rb, penalty = pen.control$penalty,
             lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
  f1 <- -d1ll(Y,ghQ,b,famL,info,f0$pD,rb) + c(nrow(Y$Y)*crossprod(cb,pMat$full))
  f2 <- ad2ll(Y,ghQ,b,famL,info,f0$pD,rb) + nrow(Y$Y)*pMat$full 
  return(list(value = -f0$ll + 0.5*nrow(Y$Y)*crossprod(cb,pMat$full)%*%cb,
              gradient = f1, hessian = f2))
}

d1ll <- function(Y,ghQ,b,famL,info,pd,rb){
  na.idx <- Y$na.idx
  Y <- Y$Y
  sc <- NULL
  w2 <- pd
  w3 <- c(ghQ$weights)
  for(i in 1:ncol(Y)){
    for(k in 1:famL[[i]]$npar){
      w1 <- famL[[i]]$dv1Y(i,Y[na.idx[,i],i],b,ghQ)[[k]]
      sc <- c(sc,colSums(ghQ$out[[k]]*colSums(w1*w2[na.idx[,i],])*w3)[rb[[k]][i,] == T])
    }
  }
  return(unname(sc))
}

d2ll <- function(Y,ghQ,b,famL,info,pd,rb){
  na.idx <- Y$na.idx
  Y <- Y$Y
  w2 <- pd
  w3 <- c(ghQ$weights)
  H <- matrix(0,sum(lb2cb(rb) == T),sum(lb2cb(rb) == T))
  idH <- 0
  for(i in 1:ncol(Y)){
    pY <- seq.int(from = max(idH+1,i), len = sum(lb2mb(rb)[i,] == T))
    h <- matrix(NA,sum(lb2mb(rb)[i,] == T),sum(lb2mb(rb)[i,] == T))
    dvO <- famL[[i]]$dv1Y(i,Y[na.idx[,i],i],b,ghQ)
    dvP <- famL[[i]]$dv2Y(i,Y[na.idx[,i],i],b,ghQ,info,dvO)
    dvC <- famL[[i]]$dvCY(i,Y[na.idx[,i],i],b,ghQ,info)
    H3 <- NULL
    id0 <- 0
    for(k in 1:famL[[i]]$npar){
      Zk <- as.matrix(ghQ$out[[k]])[,rb[[k]][i,] == T]
      dim(Zk) <- c(nrow(ghQ$points), sum(rb[[k]][i,] == T))
      pk <- seq.int(from = max(id0+1,k), len = ncol(Zk))
      w1 <- dvP[[k]]
      wa <- colSums(w1*w2[na.idx[,i],])*w3
      h1 <- lapply(1:nrow(ghQ$points), function(r) tcrossprod(Zk[r,]) * wa[r])
      h1 <- Reduce("+", h1)
      # ~~~~~~~~~~~
      wo <- dvO[[k]]
      wb <- colSums(wo^2*w2[na.idx[,i],])*w3
      h2 <- lapply(1:nrow(ghQ$points), function(r) tcrossprod(Zk[r,]) * wb[r])
      h2 <- Reduce("+", h2)
      # ~~~~~~~~~~~
      wc <- sweep(wo*w2[na.idx[,i],],2,w3,"*")
      h3 <- array(sapply(1:nrow(w1), function(r) Zk * wc[r,]),
                  dim = c(nrow(ghQ$points),ncol(Zk),nrow(w1)))
      h3 <- abind::abind(H3,h3,along = 2)
      h[pk,pk] <- h1 + h2
      id0 <- max(pk)
      if(tryCatch({!is.null(dvC[[k]])}, error = function(e) F)){
        id1 <- id0
        for(o in 1:length(dvC[[k]])){
          Zo <- as.matrix(ghQ$out[[k + o]])[,rb[[k + o]][i,] == T]
          dim(Zo) <- c(nrow(ghQ$points), sum(rb[[k + o]][i,] == T))
          po <- seq.int(from = max(id1+1,o), len = ncol(Zo))
          w1i <- dvC[[k]][[o]]
          w0i <- sweep(w1i*w2[na.idx[,i],],2,w3,"*")
          wai <- colSums(w0i)
          h1i <- lapply(1:nrow(ghQ$points), function(r) tcrossprod(Zk[r,],Zo[r,]) * wai[r])
          h1i <- Reduce("+", h1i)
          # ~~~~~~~~~~~~~~~~~~
          woi <- dvO[[k + o]]
          wbi <- colSums(wo*woi*w2[na.idx[,i],])*w3
          h2i <- lapply(1:nrow(ghQ$points), function(r) tcrossprod(Zk[r,],Zo[r,]) * wbi[r])
          h2i <- Reduce("+",h2i)
          # ~~~~~~~~~~~~~~~~~~
          h[pk,po] <- h1i + h2i
          h[po,pk] <- t(h1i + h2i) 
          id1 <- max(po)
        }
      }
      H3 <- h3
    }
    h3 <- t(sapply(1:nrow(w1), function(r) colSums(H3[,,r,drop = T])))
    h3 <- Reduce("+", lapply(1:nrow(w1), function(r) tcrossprod(h3[r,])))
    H[pY,pY] <- h - h3
    idH <- max(pY)
  }
  return(H)
}

ad2ll <- function(Y,ghQ,b,famL,info,pd,rb){
  na.idx <- Y$na.idx
  Y <- Y$Y
  w2 <- pd
  w3 <- c(ghQ$weights)
  H <- matrix(0,sum(lb2cb(rb) == T),sum(lb2cb(rb) == T))
  idH <- 0
  for(i in 1:ncol(Y)){
    pY <- seq.int(from = max(idH+1,i), len = sum(lb2mb(rb)[i,] == T))
    dvO <- famL[[i]]$dv1Y(i,Y[na.idx[,i],i],b,ghQ)
    H3 <- NULL
    id0 <- 0
    for(k in 1:famL[[i]]$npar){
      Zk <- as.matrix(ghQ$out[[k]])[,rb[[k]][i,] == T]
      dim(Zk) <- c(nrow(ghQ$points), sum(rb[[k]][i,] == T))
      pk <- seq.int(from = max(id0+1,k), len = ncol(Zk))
      wo <- dvO[[k]]
      wc <- sweep(wo*w2[na.idx[,i],],2,w3,"*")
      h3 <- array(sapply(1:nrow(wo), function(r) Zk * wc[r,]),
                  dim = c(nrow(ghQ$points),ncol(Zk),nrow(wo)))
      h3 <- abind::abind(H3,h3,along = 2)
      id0 <- max(pk)
      H3 <- h3
    }
    h3 <- t(sapply(1:nrow(wo), function(r) colSums(H3[,,r,drop = T])))
    h3 <- Reduce("+", lapply(1:nrow(wo), function(r) tcrossprod(h3[r,])))
    H[pY,pY] <- h3
    idH <- max(pY)
  }
  return(H)
}

d2llEM <- function(Y,ghQ,b,famL,info,pd,rb){
  na.idx <- Y$na.idx
  Y <- Y$Y
  w2 <- pd
  w3 <- c(ghQ$weights)
  H <- matrix(0,sum(lb2cb(rb) == T),sum(lb2cb(rb) == T))
  idH <- 0
  for(i in 1:ncol(Y)){
    pY <- seq.int(from = max(idH+1,i), len = sum(lb2mb(rb)[i,] == T))
    h <- matrix(NA,sum(lb2mb(rb)[i,] == T),sum(lb2mb(rb)[i,] == T))
    if(info == "Hessian") dvO <- famL[[i]]$dv1Y(i,Y[na.idx[,i],i],b,ghQ) else dvO <- NULL
    dvP <- famL[[i]]$dv2Y(i,Y[na.idx[,i],i],b,ghQ,info,dvO)
    dvC <- famL[[i]]$dvCY(i,Y[na.idx[,i],i],b,ghQ,info)
    id0 <- 0
    for(k in 1:famL[[i]]$npar){
      Zk <- as.matrix(ghQ$out[[k]])[,rb[[k]][i,] == T]
      dim(Zk) <- c(nrow(ghQ$points), sum(rb[[k]][i,] == T))
      pk <- seq.int(from = max(id0+1,k), len = ncol(Zk))
      w1 <- dvP[[k]]
      wa <- colSums(w1*w2[na.idx[,i],])*w3
      h1 <- lapply(1:nrow(ghQ$points), function(r) tcrossprod(Zk[r,]) * wa[r])
      h1 <- Reduce("+", h1)
      h[pk,pk] <- h1
      id0 <- max(pk)
      if(tryCatch({!is.null(dvC[[k]])}, error = function(e) F)){
        id1 <- id0
        for(o in 1:length(dvC[[k]])){
          Zo <- as.matrix(ghQ$out[[k + o]])[,rb[[k + o]][i,] == T]
          dim(Zo) <- c(nrow(ghQ$points), sum(rb[[k + o]][i,] == T))
          po <- seq.int(from = max(id1+1,o), len = ncol(Zo))
          w1i <- dvC[[k]][[o]]
          w0i <- sweep(w1i*w2[na.idx[,i],],2,w3,"*")
          wai <- colSums(w0i)
          h1i <- lapply(1:nrow(ghQ$points), function(r) tcrossprod(Zk[r,],Zo[r,]) * wai[r])
          h1i <- Reduce("+", h1i)
          h[pk,po] <- h1i
          h[po,pk] <- t(h1i) 
          id1 <- max(po)
        }
      }
    }
    H[pY,pY] <- h
    idH <- max(pY)
  }
  return(H)
}

lla <- function(cb,Y,ghQ,bg,famL,info,rb){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  return(-f0$ll)
}

plla <- function(cb,Y,ghQ,bg,famL,info,rb,pen.control){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  pMat <- pM(b,rb, penalty = pen.control$penalty,
             lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
  return(-f0$ll + 0.5*nrow(Y$Y)*crossprod(cb,pMat$full)%*%cb)
}

d1lla <- function(cb,Y,ghQ,bg,famL,info,rb){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  f1 <- -d1ll(Y,ghQ,b,famL,info,f0$pD,rb)
  return(f1)
}

d1plla <- function(cb,Y,ghQ,bg,famL,info,rb,pen.control){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  pMat <- pM(b,rb, penalty = pen.control$penalty,
             lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
  f1 <- -d1ll(Y,ghQ,b,famL,info,f0$pD,rb) + c(nrow(Y$Y)*crossprod(cb,pMat$full))
  return(f1)
}

ghq <- function(n){ # Function from GLMMADAPTIVE (2021)
  m <- trunc((n + 1) / 2)
  x <- w <- rep(-1, n)
  for (i in seq_len(m)) {
    z <- if (i == 1) {
      sqrt(2 * n + 1) - 1.85575 * (2 * n + 1)^(-0.16667)
    } else if (i == 2) {
      z - 1.14 * n^0.426/z
    } else if (i == 3) {
      1.86 * z - 0.86 * x[1]
    } else if (i == 4) {
      1.91 * z - 0.91 * x[2]
    } else {
      2 * z - x[i - 2]
    }
    for (its in seq_len(10)) {
      p1 <- 0.751125544464943
      p2 <- 0
      for (j in seq_len(n)) {
        p3 <- p2
        p2 <- p1
        p1 <- z * sqrt(2 / j) * p2 - sqrt((j - 1) / j) * p3
      }
      pp <- sqrt(2 * n) * p2
      z1 <- z
      z <- z1 - p1/pp
      if (abs(z - z1) <= 3e-14)
        break
    }
    x[i] <- z
    x[n + 1 - i] <- -z
    w[i] <- 2 / (pp * pp)
    w[n + 1 - i] <- w[i]
  }
  list(x = x, w = w)
}

ibeta <- function(Y,famL,form){
  
  Y <- Y$Y[complete.cases(Y$Y),]
  lvar <- unique(unlist(lapply(1:length(form), function(i) all.vars(form[[i]]))))
  lvar <- grep("Z", lvar, fixed = T, value = T)
  q <- length(lvar)
  Z <- scale(princomp(Y,cor = T)$scores)[,1:q, drop = F]
  colnames(Z) <- paste0("Z", 1:q)
  
  sZ <- NULL
  bstart <- NULL
  for(r in names(form)){
    sZ[[r]] <- as.data.frame(model.matrix(form[[r]],as.data.frame(Z)))
    bstart[[r]] <- matrix(NA, nrow = ncol(Y), ncol = ncol(sZ[[r]]))
    dimnames(bstart[[r]]) <- list(colnames(Y),colnames(sZ[[r]]))
  }
  
  for(i in 1:ncol(Y)){
    eq <- form
    if(famL[[i]]$family == "Beta" | famL[[i]]$family == "Beta0") tmpY <- y.(Y[,i]) else tmpY <- Y[,i]
    for(p in 1:length(eq)){ eq[[p]] <- update(eq[[p]], tmpY ~ .)  }
    if(famL[[i]]$family != "Beta0"){
      tmp <- try(gamlss::gamlss(eq$mu, sigma.formula = eq$sigma, tau.formula = eq$tau, nu.formula = eq$nu,
                                family = famL[[i]]$iuse, data = as.data.frame(cbind(tmpY,Z)),
                                control = gamlss::gamlss.control(trace = F)), silent = T)
    } else {
      eqtmp <- as.formula(paste0(deparse(eq$mu), " | ", deparse(eq$sigma[[3]])))
      tmp <- try(betareg::betareg(formula = eqtmp, data = as.data.frame(cbind(tmpY,Z)),
                                  link = famL[[i]]$link.mu, link.phi = famL[[i]]$link.sg, type = "BC"))
    }
    for(p in 1:famL[[i]]$npar){ 
      if(!inherits(tmp, "try-error")){ 
        if(famL[[i]]$family == "Beta0"){ 
          bstart[[p]][i,] <- tmp$coefficients[[p]]*0
          } else bstart[[p]][i,] <- coef(tmp,famL[[i]]$pars[p])
      } else {
        bstart[[p]][i,] <- 0.01
      }
    }
  }
  
  for(r in names(bstart)){
    for(j in which(grepl("Z",colnames(bstart[[r]])))){
      c <- as.integer(substr(colnames(bstart[[r]])[j],nchar(colnames(bstart[[r]])[j]),nchar(colnames(bstart[[r]])[j])))
      if(bstart[[r]][c,j] < 0) bstart[[r]][,j] <- -bstart[[r]][,j]
    } }
  return(bstart)
}

rmat <- function(res,b){
  rest <- l2rm(res)
  frparm <- rb <- b
  for(k in names(b)){ rb[[k]] <- b[[k]]*NA  }
  for(i in 1:nrow(rest)){
    frparm[[as.character(rest[i,1])]][as.character(rest[i,2]),as.character(rest[i,3])] <- as.numeric(rest[i,4])
    rb[[as.character(rest[i,1])]][as.character(rest[i,2]),as.character(rest[i,3])] <- as.numeric(rest[i,4])
  }
  for(k in names(rb)){ rb[[k]] <- is.na(rb[[k]])  }
  return(list(b = frparm, rb = rb))
}

l2rm <- function(rlist){
  if(!is.list(rlist)) stop("Restrictions should be a list with element(s), each of the type c('parameter',item,'restricted variable',value)")
  r <- as.data.frame(matrix(unlist(rlist), ncol=4, byrow = T))
  return(r)
}

mvghQ <- function(n, mu, psi, formula = ~ Z1 + Z2) {
  nl <- unique(unlist(lapply(1:length(formula), function(i) all.vars(as.formula(formula[[i]])))))
  nl <- nl[grep("Z", nl, fixed = T)]
  if(missing(mu)) mu <- rep(0,length(nl)); 
  if(missing(psi)) psi <- diag(length(nl))
  if(!all(dim(psi) == length(mu))) stop("muZ and PsiZ have nonconformable dimensions (in ghQ)")
  if(length(mu) > 0) dm  <- length(mu) else dm <- 1
  if(missing(n)) n <- round(100^(1/dm))
  gh  <- list("points" = ghq(n)$x, "weights" = ghq(n)$w)
  gh$weights <-  gh$weights * (2*pi)^(-1/2)*exp(gh$points^2/2)
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh$points[idx],nrow(idx),dm)
  wts <- as.matrix(apply(matrix(gh$weights[idx],nrow(idx),dm), 1, prod))
  wts <- wts + (1-sum(wts))/(length(wts))
  # Rotating if mu & psi != NULL
  eig <- eigen(psi) 
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  colnames(pts) <- nl
  out <- NULL
  for(i in names(formula)){
    out[[i]] <- as.data.frame(model.matrix(as.formula(formula[[i]]), as.data.frame(pts)))
  }
  return(list(points = pts, weights = c(wts), out = out))
}

lb2mb <- function(b){
  l2m <- NULL
  for(i in 1:length(b)){ l2m <- cbind(l2m,b[[i]]) }
  return(unname(l2m))
}

lb2cb <- function(b){
  return(c(t(lb2mb(b))))
}

cb2lb <- function(cb,b){
  bm <- matrix(cb,byrow=T,nrow=nrow(b$mu))
  c2l <- vector(mode = "list", length = length(b))
  names(c2l) <- names(b)
  for(i in seq_along(names(b))){
    if(i == 1) c2l[[i]] <- bm[,1:ncol(b[[i]]), drop = F]
    else c2l[[i]] <- bm[,seq(from = (ncol(b[[i-1]])+1), length = ncol(b[[i]])), drop = F]
    colnames(c2l[[i]]) <- colnames(b[[i]]); rownames(c2l[[i]]) <- rownames(b[[i]])
  }
  return(c2l)
}

mb2lb <- function(mb,b){
  m2l <- vector(mode = "list", length = length(b))
  names(m2l) <- names(b)
  for(i in seq_along(names(b))){
    if(i == 1) m2l[[i]] <- mb[,1:ncol(b[[i]]), drop = F]
    else m2l[[i]] <- mb[,seq(from = (ncol(b[[i-1]])+1), length = ncol(b[[i]])), drop = F]
    colnames(m2l[[i]]) <- colnames(b[[i]]); rownames(m2l[[i]]) <- rownames(b[[i]])
  }
  return(m2l)
}

m2pdm <- function(mat){
  eS <- eigen(mat, symmetric = TRUE)
  e.val <- eS$values
  e.vec <- eS$vectors
  check.eigen <- any(e.val <= 0)
  if(check.eigen == T){
    n.e.val <- e.val[e.val <= 0]
    s <- sum(e.val[n.e.val])*2  
    t <- s^2*100 + 1
    p <- min(e.val[(e.val <= 0) == F])
    e.val[e.val <= 0] <- p*(s - n.e.val)^2/t
    D <- diag(e.val)
    D.inv <- diag(1/e.val)
    res <- e.vec %*% D %*% t(e.vec) 
    res.inv <- e.vec %*% D.inv %*% t(e.vec)
    res.sqr <- e.vec %*% sqrt(D)
  } else { res <- mat; res.inv <- e.vec %*% diag(1/e.val) %*% t(e.vec) ; res.sqr <- e.vec %*% diag(sqrt(e.val))} 
  res.inv <- (res.inv + t(res.inv) ) / 2 
  return(list(mat = res, inv.mat = res.inv, sqr.mat = t(res.sqr), is.PDM = !check.eigen))
}

pM <- function(b, rb, penalty = "lasso", lambda = 0.1, w.alasso = NULL, a = NULL){
  penidx <- function(rb){
    for(i in names(rb)){
      rb[[i]][,colnames(rb[[i]]) == "(Intercept)"] <- F
      nQ <- sum(grepl("Z",colnames(rb[[i]]), fixed = T))
      rb[[i]][seq_len(nQ), ] <- F
    }
    return(rb)
  }
  lambda2lb <- function(lamlist,b){
    for(i in 1:length(b)) b[[i]] <- matrix(lamlist[i], nrow = nrow(b[[i]]), ncol = ncol(b[[i]]),
                                           dimnames = dimnames(b[[i]]))
    return(b)
  }
  id <- lb2cb(penidx(rb))
  param <- lb2cb(b)[id]
  if(!is.null(w.alasso) && is.list(w.alasso)){ w.alasso <- lb2cb(w.alasso)[id] }
  if(length(lambda) != length(b)) lambda <- rep(lambda, length(b))
  lambda <- lb2cb(lambda2lb(lambda,b))[id]
  S. <- Sraw <- diag(0,length(id))
  eps = 1e-7 # sqrt(.Machine$double.eps) # protective tolerance level
  if(penalty == "ridge"){
    A1 <-  c(lambda)*rep(1, length(param))
  }
  if(penalty == "lasso"){
    A1 <- c(lambda)/sqrt(param^2 + eps) 
  }
  if(penalty == "alasso"){
    if(is.null(a)) a = 2
    if( is.null(w.alasso) ) w.alasso <- 1
    w.al <- abs(w.alasso)^a
    A1 <- c(lambda)/(w.al*sqrt(param^2 + eps))
  }
  if(penalty == "scad"){
    if(is.null(a)) a = 3.7
    theta <- abs(param)
    f1 <- sapply(1:length(theta), function(i) { max(a*c(lambda[i]) - theta[i], 0)/((a-1)*c(lambda[i]) + eps) })
    f.d <- ((theta <= c(lambda)) + f1 * (theta > c(lambda)))
    A1 <- c(lambda)* f.d / ( sqrt(param^2 + eps) )
  }
  if(penalty == "mcp"){
    if(is.null(a)) a = 2.5
    theta <- abs(param) 
    f.d <- (c(lambda)-theta/a)*(theta < c(lambda)*a)
    A1 <- f.d / ( sqrt(param^2 + eps) )
  }
  if(length(A1) == 1) S <- matrix(A1) else S <- diag(A1)
  diag(S.)[id] <- A1
  diag(Sraw)[id] <- A1/lambda
  S. <- S.[lb2cb(rb),lb2cb(rb)]
  Sraw <- Sraw[lb2cb(rb),lb2cb(rb)]
  return(list(full = S., red = S, raw = Sraw))
}

SSE <- function(loglambda,b,gra,hes,rb,pml,Y){
  Y <- Y$Y
  lambda = exp(loglambda)
  S <- pM(b = b, rb = rb, penalty = pml$penalty, lambda = lambda, w.alasso = pml$w.alasso, a = pml$a)
  th <- lb2cb(b)[lb2cb(rb)]
  sJ <- m2pdm(hes)$sq
  K = sJ%*%th + solve(t(sJ))%*%gra
  A <- sJ%*%m2pdm(hes + nrow(Y)*S$full)$inv.mat%*%t(sJ)
  sse <- c(crossprod(K-A%*%K)) + 2*pml$gamma*sum(diag(A)) - sum(lb2cb(rb))
  return(sse)
}

op.lambda <- function(Y,ghQ,b,famL,info,rb,pen.control){
  fyz_c <- fyz(Y,ghQ,b,famL)
  gra <- d1ll(Y,ghQ,b,famL,info,fyz_c$pD,rb)
  hes <- -d2ll(Y,ghQ,b,famL,info,pd = fyz_c$pD,rb)
  lambda_ <- pen.control$lambda
  if(length(lambda_) != length(b)) lambda_ <- rep(lambda_, length(b))
  sse <- numeric(pen.control$iter.lim)
  lambda <- matrix(NA,ncol = length(lambda_),nrow = pen.control$iter.lim)
  i <- 1
  sse[i] <- SSE(log(lambda_),b,gra,hes,rb,pen.control,Y)
  lambda[i,] <- lambda_
  opt <- T
  
  while(opt){
    i <- i + 1
    d1sse <- numDeriv::grad(SSE,log(lambda_),b = b, gra = gra, hes = hes,
                            rb = rb, pml = pen.control, Y = Y)
    d2sse <- numDeriv::hessian(SSE,log(lambda_),b = b, gra = gra, hes = hes,
                            rb = rb, pml = pen.control, Y = Y)
    lambda[i,] <- lambda_ <- exp(log(lambda_) - solve(d2sse, d1sse))
    sse[i] <- SSE(log(lambda_),b,gra,hes,rb,pen.control,Y)
    if(sse[i] - sse[i-1] > 0){ grdes <- T; opt <- F ; i <- i - 1}
    if(i == pen.control$iter.lim){ grdes <- opt <- F }
  }
  
  lambda_ <- lambda[i,]; 
  
  while(grdes){
    i <- i + 1
    d1sse <- numDeriv::grad(SSE,log(lambda_),b = b, gra = gra, hes = hes,
                            rb = rb, pml = pen.control, Y = Y)
    lambda[i,] <- lambda_ <- exp(log(lambda_) - (0.01/i)*c(d1sse))
    sse[i] <- SSE(log(lambda_),b,gra,hes,rb,pen.control,Y)
    if(sse[i] - sse[i-1] > 0){ stdes <- T; grdes <- F ; i <- i - 1}
    if(i == pen.control$iter.lim){ stdes <- grdes <- F }
  }
  
  lambda_ <- lambda[i,]; 
  
  while(stdes){
    i <- i + 1
    d1sse <- numDeriv::grad(SSE,log(lambda_),b = b, gra = gra, hes = hes,
                            rb = rb, pml = pen.control, Y = Y)
    lambda[i,] <- lambda_ <- exp(log(lambda_) - c(d1sse))
    sse[i] <- SSE(log(lambda_),b,gra,hes,rb,pen.control,Y)
    if(sse[i] - sse[i-1] > 0){ stdes <- F; i <- i - 1}
    if(i == pen.control$iter.lim){ stdes <- F }
  }
  
  return(list(lambda = lambda[i,], miditer = i, sse = sse[i]))
}

op.lambda1 <- function(Y,ghQ,b,famL,info,rb,pen.control){
  fyz_c <- fyz(Y,ghQ,b,famL)
  gra <- d1ll(Y,ghQ,b,famL,info,fyz_c$pD,rb)
  hes <- -d2ll(Y,ghQ,b,famL,info,pd = fyz_c$pD,rb)
  lambda_ <- pen.control$lambda
  if(length(lambda_) != length(b)) lambda_ <- rep(lambda_, length(b))
  D1SSE <- function(loglambda_, b = b, gra = gra, hes = hes,
                    rb = rb, pml = pen.control, Y = Y){
      d1sse <- numDeriv::grad(SSE,loglambda_,b = b, gra = gra, hes = hes,
                              rb = rb, pml = pen.control, Y = Y)
      return(d1sse) }
  
  # out_ <- optim(par = log(lambda_), fn = SSE, gr = D1SSE, method = "L-BFGS-B",
  #               b = b, gra = gra, hes = hes, rb = rb, pml = pen.control, Y = Y)
  # 
  # return(list(lambda = exp(out_$par), miditer = out_$counts[1], sse = out_$value))
  
  D2SSE <- function(loglambda_, b = b, gra = gra, hes = hes,
                    rb = rb, pml = pen.control, Y = Y){
      d2sse <- numDeriv::hessian(SSE,loglambda_,b = b, gra = gra, hes = hes,
                              rb = rb, pml = pen.control, Y = Y)
      return(d2sse) }
  SSEtrust <- function(loglambda_t, bt, grat, hest,
                       rbt, pmlt, Yt){
      valuer <- SSE(loglambda_t, bt, grat, hest, rbt, pmlt,Yt)
      gradientr <- D1SSE(loglambda_t, bt, grat, hest, rbt, pmlt,Yt)
      hessianr <- D2SSE(loglambda_t, bt, grat, hest, rbt, pmlt,Yt)
      return(list(value = valuer, gradient = gradientr, hessian = hessianr)) }

  out__ <- trust::trust(objfun = SSEtrust, parinit = log(lambda_), rinit = 1, rmax = 5,
                       fterm = sqrt(.Machine$double.eps), iterlim = 100,
                       bt = b, grat = gra, hest = hes, rbt = rb, pmlt = pen.control, Yt = Y)

  return(list(lambda = exp(out__$arg), miditer = out__$iter, sse = out__$value))
}

GAIC <- function(mod){
  if(class(mod) != "glvmlss") stop("Model ('mod') object should be of class 'glvmlss'")
  ll <- mod$unploglik
  H <- mod$hes$H
  iH <- tryCatch({solve(H)}, error = function(e){m2pdm(H)$inv})
  
  if(!is.null(mod$hes$Hp)){
    iHp <- tryCatch({solve(mod$hes$Hp)}, error = function(e){m2pdm(mod$hes$Hp)$inv})
    GIC <- -2*ll + 2*sum(diag(iHp%*%H)) } else {
      GIC <- -2*ll + 2*sum(diag(iH%*%H))
    }
  return(GIC)
}

GBIC <- function(mod){
  if(class(mod) != "glvmlss") stop("Model ('mod') object should be of class 'glvmlss'")
  ll <- mod$unploglik
  H <- mod$hes$H
  iH <- tryCatch({solve(H)}, error = function(e){m2pdm(H)$inv})
  
  if(!is.null(mod$hes$Hp)){
    iHp <- tryCatch({solve(mod$hes$Hp)}, error = function(e){m2pdm(mod$hes$Hp)$inv})
    GIC <- -2*ll + log(mod$n)*sum(diag(iHp%*%H)) } else {
      GIC <- -2*ll + log(mod$n)*sum(diag(iH%*%H))
    }
  return(GIC)
}