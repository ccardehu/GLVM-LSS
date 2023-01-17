y. <- function(x){
 if (any(x == 1, na.rm = T)) x[x == 1] <- 1 - 1e-3 #.Machine$double.eps^(1/3)
 if (any(x == 0, na.rm = T)) x[x == 0] <- 1e-3 #.Machine$double.eps^(1/3)
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
  # Solve for Rz here? (i.e., update ghQ)
  f0 <- fyz(Y,ghQ,b,famL)
  f1 <- -d1ll(Y,ghQ,b,famL,info,f0$pD,rb)
  f2 <- -d2ll(Y,ghQ,b,famL,info,f0$pD,rb)
  return(list(value = -f0$ll, gradient = f1, hessian = f2))
}

lazyll <- function(cb,Y,ghQ,bg,famL,info,rb){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  # Solve for Rz here? (i.e., update ghQ)
  f0 <- fyz(Y,ghQ,b,famL)
  f1 <- -d1ll(Y,ghQ,b,famL,info,f0$pD,rb)
  f2 <- ad2ll(Y,ghQ,b,famL,info,f0$pD,rb)
  return(list(value = -f0$ll, gradient = f1, hessian = f2))
}

pll <- function(cb,Y,ghQ,bg,famL,info,rb,bp,pen.control){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  pMat <- pM(b,rb,bp,penalty = pen.control$penalty,
             lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
  f1 <- -d1ll(Y,ghQ,b,famL,info,f0$pD,rb) + c(nrow(Y$Y)*crossprod(cb,pMat$full))
  f2 <- -d2ll(Y,ghQ,b,famL,info,f0$pD,rb) + nrow(Y$Y)*pMat$full 
  return(list(value = -f0$ll + 0.5*nrow(Y$Y)*crossprod(cb,pMat$full)%*%cb,
              gradient = f1, hessian = f2))
}

lazypll <- function(cb,Y,ghQ,bg,famL,info,rb,bp,pen.control){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  # Solve for Rz here? (i.e., update ghQ)
  f0 <- fyz(Y,ghQ,b,famL)
  pMat <- pM(b,rb,bp,penalty = pen.control$penalty,
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

plla <- function(cb,Y,ghQ,bg,famL,info,rb,bp,pen.control){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  pMat <- pM(b,rb,bp,penalty = pen.control$penalty,
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

d1plla <- function(cb,Y,ghQ,bg,famL,info,rb,bp,pen.control){
  b <- lb2cb(bg)
  b[lb2cb(rb) == T] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  pMat <- pM(b,rb,bp,penalty = pen.control$penalty,
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
      if (abs(z - z1) <= 3e-14){ break }
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
      tmp <- try(gamlss::gamlss(eq$mu, sigma.formula = eq$sigma, nu.formula = eq$nu, tau.formula = eq$tau,
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
          bstart[[p]][i,] <- rep(tmp$coefficients[[p]],nrow(bstart[[p]][i,,drop=F]))
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
    frparm[[as.character(rest[i,1])]][as.character(rest[i,2]),as.character(rest[i,3])] <- as.numeric(as.character(rest[i,4]))
    rb[[as.character(rest[i,1])]][as.character(rest[i,2]),as.character(rest[i,3])] <- as.numeric(as.character(rest[i,4]))
  }
  for(k in names(rb)){ rb[[k]] <- is.na(rb[[k]]) & !is.na(b[[k]])}
  penb <- rb
  for(k in names(penb)){ penb[[k]][,grepl("(Intercept)", colnames(penb[[k]]), fixed = T)] <- F  }
  return(list(b = frparm, rb = rb, penb = penb))
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
  if(!all(dim(psi) == length(mu))) stop("muZ and Rz have nonconformable dimensions (in ghQ)")
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
  ic = 1
  for(i in seq_along(names(b))){
    colsel = seq(from = ic, length = ncol(b[[i]]))
    c2l[[i]] <- bm[,colsel, drop = F]
    ic = max(colsel) + 1
    colnames(c2l[[i]]) <- colnames(b[[i]]); rownames(c2l[[i]]) <- rownames(b[[i]])
  }
  return(c2l)
}

mb2lb <- function(mb,b){
  m2l <- vector(mode = "list", length = length(b))
  names(m2l) <- names(b)
  ic = 1
  for(i in seq_along(names(b))){
    colsel = seq(from = ic, length = ncol(b[[i]]))
    m2l[[i]] <- mb[,colsel, drop = F]
    ic = max(colsel) + 1
    colnames(m2l[[i]]) <- colnames(b[[i]]); rownames(m2l[[i]]) <- rownames(b[[i]])
  }
  return(m2l)
}

m2pdm <- function(mat){
  eS <- eigen(mat, symmetric = TRUE)
  e.val <- eS$values
  e.vec <- eS$vectors
  check.eigen <- any(e.val <= 0)
  # if(check.eigen == T){
  #   n.e.val <- e.val[e.val <= 0]
  #   # s <- sum(e.val[!n.e.val])
  #   s <- sum(e.val[n.e.val])*2
  #   t <- s^2*100 + 1
  #   p <- min(e.val[(e.val <= 0) == F])
  #   # p <- min(e.val[(e.val <= 0) == F])
  #   e.val[e.val <= 0] <- p*(s - n.e.val)^2/t
  #   D <- diag(e.val)
  #   D.inv <- diag(1/e.val)
  #   res <- e.vec %*% D %*% t(e.vec)
  #   res.inv <- e.vec %*% D.inv %*% t(e.vec)
  #   res.sqr <- e.vec %*% sqrt(D)
  # } else { res <- mat; res.inv <- e.vec %*% diag(1/e.val) %*% t(e.vec) ; res.sqr <- e.vec %*% diag(sqrt(e.val))}
  # res <- (res + t(res) ) / 2
  # res.inv <- (res.inv + t(res.inv) ) / 2
  if(check.eigen == T){
    neg.e.val <- e.val[e.val <= 0]
    s <- sum(e.val[-1])
    t <- s^2*100 + 1
    p <- min(e.val[e.val > 0])
    e.val[e.val <= 0] <- p*(s - neg.e.val)^2/t
    D <- diag(e.val)
    D.inv <- diag(1/e.val)
    res <- e.vec %*% D %*% t(e.vec)
    res.inv <- e.vec %*% D.inv %*% t(e.vec)
    res.sqr <- e.vec %*% sqrt(D)
  } else { res <- mat; res.inv <- e.vec %*% diag(1/e.val) %*% t(e.vec) ; res.sqr <- e.vec %*% diag(sqrt(e.val))}
  res <- (res + t(res) ) / 2
  res.inv <- (res.inv + t(res.inv) ) / 2
  # tobj <- nearPD(mat)
  # res <- tobj$mat
  # e.vec <- tobj$vectors
  # e.val <- tobj$values
  # res.inv <- e.vec %*% diag(1/e.val) %*% t(e.vec) ; res.sqr <- e.vec %*% diag(sqrt(e.val)) 
  # } else {  res <- mat; res.inv <- e.vec %*% diag(1/e.val) %*% t(e.vec) ; res.sqr <- e.vec %*% diag(sqrt(e.val)) }
  return(list(mat = res, inv.mat = res.inv, sqr.mat = t(res.sqr), is.PDM = !check.eigen))
}

pM <- function(b, rb, bp, penalty = "lasso", lambda = 0.1, w.alasso = NULL, a = NULL){
  lambda2lb <- function(lamlist,b){
    for(i in 1:length(b)) b[[i]] <- matrix(lamlist[i], nrow = nrow(b[[i]]), ncol = ncol(b[[i]]),
                                           dimnames = dimnames(b[[i]]))
    return(b)
  }
  id <- lb2cb(bp)
  param <- lb2cb(b)[id]
  if(!is.null(w.alasso) && is.list(w.alasso)){ w.alasso <- lb2cb(w.alasso)[id] }
  if(length(lambda) != length(b)) lambda <- rep(lambda, length(b))
  lambda <- lb2cb(lambda2lb(lambda,b))[id]
  S. <- Sraw <- diag(0,length(id))
  eps = 1e-8
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
  lambdaraw <- lambda
  lambdaraw[lambdaraw == 0] <- 1
  diag(Sraw)[id] <- A1/lambdaraw
  S. <- S.[lb2cb(rb),lb2cb(rb)]
  Sraw <- Sraw[lb2cb(rb),lb2cb(rb)]
  return(list(full = S., red = S, raw = Sraw))
}

SSE <- function(loglambda,b,gra,hes,rb,bp,pml,Y){
  Y <- Y$Y
  lambda <- numeric(length(b))
  lambda[pml$lambda.auto] = exp(loglambda)
  S <- pM(b = b, rb = rb, bp = bp, penalty = pml$penalty, lambda = lambda, w.alasso = pml$w.alasso, a = pml$a)
  th <- lb2cb(b)[lb2cb(rb)]
  sJ <- m2pdm(hes)$sq
  K = sJ%*%th + solve(t(sJ))%*%gra
  A <- sJ%*%m2pdm(hes + nrow(Y)*S$full)$inv.mat%*%t(sJ)
  sse <- c(crossprod(K-A%*%K)) + 2*pml$gamma*sum(diag(A)) - sum(lb2cb(rb))
  return(sse)
}

op.lambda <- function(Y,ghQ,b,famL,info,rb,bp,pen.control){
  fyz_c <- fyz(Y,ghQ,b,famL)
  gra <- d1ll(Y,ghQ,b,famL,info,fyz_c$pD,rb)
  hes <- -d2ll(Y,ghQ,b,famL,info,pd = fyz_c$pD,rb)
  lambda_ <- numeric(length(b))
  lambda_[pen.control$lambda.auto] <- pen.control$lambda[pen.control$lambda.auto]
  
  d1SSE <- function(loglambda_, b = b, gra = gra, hes = hes,
                    rb = rb, bp = bp, pml = pen.control, Y = Y){
      d1sse <- numDeriv::grad(SSE,loglambda_,b = b, gra = gra, hes = hes,
                              rb = rb, bp = bp, pml = pen.control, Y = Y)
      return(d1sse) }
  
  # out_ <- optim(par = log(lambda_), fn = SSE, gr = D1SSE, method = "L-BFGS-B",
  #               b = b, gra = gra, hes = hes, rb = rb, pml = pen.control, Y = Y)
  # 
  # return(list(lambda = exp(out_$par), miditer = out_$counts[1], sse = out_$value))
  
  d2SSE <- function(loglambda_, b = b, gra = gra, hes = hes,
                    rb = rb, bp = bp, pml = pen.control, Y = Y){
      d2sse <- numDeriv::hessian(SSE,loglambda_,b = b, gra = gra, hes = hes,
                              rb = rb, bp = bp, pml = pen.control, Y = Y)
      return(d2sse) }
  
  SSEtrust <- function(loglambda_t, bt, grat, hest,
                       rbt, bpt, pmlt, Yt){
      valuer <- SSE(loglambda_t, bt, grat, hest, rbt, bpt, pmlt,Yt)
      gradientr <- d1SSE(loglambda_t, bt, grat, hest, rbt, bpt, pmlt,Yt)
      hessianr <- d2SSE(loglambda_t, bt, grat, hest, rbt, bpt, pmlt,Yt)
      return(list(value = valuer, gradient = gradientr, hessian = hessianr)) }

  out_ <- trust::trust(objfun = SSEtrust, parinit = log(lambda_[pen.control$lambda.auto]), rinit = 1, rmax = 5,
                       fterm = sqrt(.Machine$double.eps), iterlim = 100,
                       bt = b, grat = gra, hest = hes, rbt = rb, bpt = bp, pmlt = pen.control, Yt = Y)

  lambda_[pen.control$lambda.auto] <- exp(out_$arg)
  return(list(lambda = lambda_, miditer = out_$iter, sse = out_$value))
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

fixrb <- function(rb,b,tolb){
  for(i in names(b)){
    rb[[i]] <- (rb[[i]] & !(abs(b[[i]]) <= tolb))
  }
  return(rb)
}

nearPD <- function (M, eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08,
                    maxits = 100) {
  if (!(is.numeric(M) && is.matrix(M) && identical(M, t(M))))
    stop("Input matrix M must be square and symmetric.\n")
  inorm <- function(x) max(rowSums(abs(x)))
  n <- ncol(M)
  U <- matrix(0.0, n, n)
  X <- M
  iter <- 0
  converged <- FALSE
  while (iter < maxits && !converged) {
    Y <- X
    T. <- Y - U
    e <- eigen(Y, symmetric = TRUE)
    Q <- e$vectors
    d <- e$values
    D <- if (length(d) > 1) diag(d) else as.matrix(d)
    p <- (d > eig.tol * d[1])
    QQ <- Q[, p, drop = FALSE]
    X <- QQ %*% D[p, p, drop = FALSE] %*% t(QQ)
    U <- X - T.
    X <- (X + t(X)) / 2
    conv <- inorm(Y - X)/inorm(Y)
    iter <- iter + 1
    converged <- conv <= conv.tol
  }
  X <- (X + t(X)) / 2
  e <- eigen(X, symmetric = TRUE)
  d <- e$values
  Eps <- posd.tol * abs(d[1L])
  if (d[n] < Eps) {
    d[d < Eps] <- Eps
    Q <- e$vectors
    o.diag <- diag(X)
    X <- Q %*% (d * t(Q))
    D <- sqrt(pmax(Eps, o.diag) / diag(X))
    X[] <- D * X * rep(D, each = n)
  }
  X <- (X + t(X)) / 2
  e <- eigen(X, symmetric = TRUE)
  return(list(mat = X, values = e$values, vectors = e$vectors))
}

# # upSa <- function(ghQ,pD){ # pz,ghQ,pD,respz
# #  
# #  if(is.null(respz)){ corlst <- which(lower.tri(pz),arr.ind = T) } else corlst <- which(lower.tri(pz),arr.ind = T)[-respz,]
# #  V <- Reduce("+",lapply(1:nrow(pD), function(m){
# #   Reduce("+",lapply(1:nrow(ghQ$points), function(i) tcrossprod(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i])) }))
# #  sc <- -nrow(pD)/2*solve(pz) + 0.5*solve(pz)%*%V%*%solve(pz)
# #  sc <- 2*sc - diag(diag(sc))
# #  sc <- c(sc[lower.tri(sc)],sc[upper.tri(sc)])
# #  kR <- kronecker(solve(pz),solve(pz))
# #  for(y in 1:nrow(corlst)){
# #   idr <- (((y-1)*nrow(pz))+nrow(pz)+1):((y)*nrow(pz)+nrow(pz))
# #   idc <- (((y-1)*nrow(pz))+nrow(pz)+1):((y)*nrow(pz)+nrow(pz))
# #   he <- nrow(pD)/2*(kR[id,id]) - 0.5*(solve(pz)%*%solve(pz)%*%V%*%solve(pz) - solve(pz)%*%V%*%solve(pz)%*%solve(pz))
# #  }
# # 
# #  # he <- 2*he - diag(diag(he))
# #  # for(y in 1:nrow(corlst)){
# #  #  pz[corlst[y,1],corlst[y,2]] <- pz[corlst[y,1],corlst[y,2]] - solve(he)[corlst[y,1],corlst[y,2]]*sc[corlst[y,1],corlst[y,2]]
# #  # }
# #  # pz[upper.tri(pz)] <- pz[lower.tri(pz)]
# #   
# #  
# #  pz <- (1/nrow(pD))*Reduce("+",lapply(1:nrow(pD), function(m){
# #   Reduce("+",lapply(1:nrow(ghQ$points), function(i) tcrossprod(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i]))}))
# #  pz <- diag(1/sqrt(diag(pz)))%*%pz%*%diag(1/sqrt(diag(pz)))
# #  return(pz)
# # }
# 
# upS <- function(pz,ghQ,pD,respz){
#  if(is.null(respz)){ corlst <- which(lower.tri(pz),arr.ind = T) } else corlst <- which(lower.tri(pz),arr.ind = T)[-respz,]
#  sv <- pz[lower.tri(pz)]
#  gra <- matrix(0,nrow = length(sv))
#  hess <- matrix(0, nrow = length(sv), ncol = length(sv))
#  Iq <- diag(nrow(pz))
#  for(y in 1:nrow(corlst)){
#   Dy <- matrix(0,nrow = nrow(pz),ncol(pz)); Dy[corlst[y,1],corlst[y,2]] <- 1
#   Dy <- Dy + t(Dy) - Dy%*%Dy
#   Gy <- solve(pz)%*%Dy%*%solve(pz)
#   # Fy <- solve(pz)%*%Dy%*%Gy + Gy%*%Dy%*%solve(pz)
#   P1 <- -nrow(pD)/2*sum(diag((2*solve(pz) - solve(pz)*Iq)%*%Dy))
#   P2 <- Reduce("+",lapply(1:nrow(pD), function(m){sum(diag(
#         Gy%*%Reduce("+",lapply(1:nrow(ghQ$points), function(i) tcrossprod(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i])) )) }))
#   P3 <- Reduce("+",lapply(1:nrow(pD), function(m){
#         v <- Reduce("+",lapply(1:nrow(ghQ$points), function(i) matrix(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i]));
#         t(v)%*%Gy%*%v }))
#   gra[y] <- P1+P2+P3
#  }
#  pz[lower.tri(pz)] <- matrix(sv) - solve(hess)%*%gra
# 
#   #
#   # H1 <- nrow(pD)/2*sum(diag(Gy%*%Dy))
#   # H2 <- -0.5*Reduce("+",lapply(1:nrow(pD), function(m){sum(diag(
#   # Fy%*%Reduce("+",lapply(1:nrow(ghQ$points), function(i) tcrossprod(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i])) )) }))
#   # H3 <- -0.5*Reduce("+",lapply(1:nrow(pD), function(m){
#   #   v <- Reduce("+",lapply(1:nrow(ghQ$points), function(i) matrix(ghQ$points[i,])*pD[m,i]*c(ghQ$weights)[i]));
#   #   t(v)%*%Fy%*%v }))
#   # pz[corlst[y,1],corlst[y,2]] <- pz[corlst[y,1],corlst[y,2]] - c(P1+P2+P3)/c(H1+H2+H3)
#   #
#   # gradient <- c(gradient, P1+P2+P3)
#   # hessian <- c(hessian, H1+H2+H3)
#  # }
#  # if(length(hessian) == 1) hessian <- matrix(hessian) else hessian <- diag(hessian)
#  pz[upper.tri(pz)] <- pz[lower.tri(pz)]
#  return(list(pz = pz, gradient = c(gra), hessian = hess))
# }
 
# grL <- function(lvec,ghQf,pDf,qf){
#   
#   Lf <- diag(qf)
#   lvecorg <- lvec
#   lvec <- c(1,lvec)
#   Lf[upper.tri(Lf, diag = T)] <- lvec
#   Lf <- t(Lf)
#   gra <- vector(mode = "numeric", length(lvecorg))
#   
#   auxF <- function(m,ghQff,pDff,Gyf){
#     zm  <- Reduce("+",lapply(1:nrow(ghQff$points), function(i){ matrix(ghQff$points[i,])*pDff[m,i]*c(ghQff$weights)[i] }))
#     z2m <- Reduce("+",lapply(1:nrow(ghQff$points), function(i){ tcrossprod(matrix(ghQff$points[i,]) - zm)*pDff[m,i]*c(ghQff$weights)[i] }))
#     
#     P2 <- sum(diag(Gyf%*%z2m))
#     P3 <- t(zm)%*%Gyf%*%zm
#     return(P2 + P3)
#   }
#   
#   for(i in 1:length(lvecorg)){
#     Dy <- matrix(0, nrow = nrow(Lf), ncol(Lf))
#     Dy[which(Lf == lvecorg[i])] <- 1
#     Ay <- t(Lf)%*%solve(Lf%*%t(Lf))
#     Gy <- solve(Lf%*%t(Lf))%*%Dy%*%Ay
#     P1 <- -nrow(pDf)*sum(diag((Ay%*%Dy)))
#     P2P3 <- Reduce("+",lapply(1:nrow(pDf), function(m){auxF(m,ghQf,pDf,Gy) }))
#     gra[i] <- P1+P2P3
#   }
#   
#   return(gra)
# }

# rmA <- grL(lvecT, ghQf = ghQf ,pDf = pDf,qf = qf)
# rmB <- numDeriv::grad(fun = opRz, x = lvecT, ghQf = ghQf ,pDf = pDf,qf = qf)

# opRz <- function(lvec,ghQf,pDf,qf){
#   Lf <- diag(qf)
#   lvec <- c(1,lvec)
#   Lf[upper.tri(Lf, diag = T)] <- lvec
#   llkZ <- sum(sapply(1:nrow(pDf), function(m){ sum(mvnfast::dmvn(ghQf$points, mu = rep(0,qf), sigma = t(Lf)%*%Lf, log = T)*pDf[m,]*c(ghQf$weights)) } ) )
#   return(llkZ)
# }

# upRz <- function(cRz, ghQ, pD, q, max_it = 100){
#   tmp <- optim(par = cRz, fn = opRz, ghQf = ghQ, pDf = pD, qf = q,
#                method = "L-BFGS-B", lower = -0.95, upper = 0.95, control = list(fnscale = -1, maxit = max_it))
#   Rz <- diag(q)
#   Rz[lower.tri(Rz)] <- tmp$par
#   Rz[upper.tri(Rz)] <- c(t(Rz)[upper.tri(Rz)])
#   Rz <- diag(1/sqrt(diag(Rz)))%*%Rz%*%diag(1/sqrt(diag(Rz)))
#   return(Rz)
# }

newRz <- function(Rz,pD,ghQ,q){
  
  L <- chol(Rz) # Rz = t(L)%*%L

  objF <- function(lvec,colNf,Lf,ghQf,pDf,qf){
    tlvec <- vector(mode = "numeric", length = qf)
    tlvec[1:length(lvec)] <- lvec
    Lf[,colNf] <- tlvec
    llkZ <- sum(sapply(1:nrow(pDf), function(m){ sum(mvnfast::dmvn(ghQf$points[,], mu = rep(0,qf), sigma = t(Lf)%*%Lf, log = T)*pDf[m,]*c(ghQf$weights)) } ) )
    return(llkZ)
  }
  
  for(i in 2:q){
    lvec <- L[,i][upper.tri(L,diag = T)[,i]]
    tmp <- optim(lvec, fn = objF, colNf = i, Lf = L, ghQf = ghQ, pDf = pD, qf = q,
                 control = list(fnscale = -1, maxit = 100))
    lvec <- vector(mode = "numeric", length = q)
    lvec[1:length(tmp$par)] <- tmp$par/sqrt(c(crossprod(tmp$par)))
    L[,i] <- lvec
  }
  
  .Rz <- t(L)%*%L
  
  return(.Rz)
}

# numDeriv::grad(func = LDLRz, x = pvec, ghQ = ghQ, pD = pD, q = q)
# grRz(.crz = c(Rz[lower.tri(Rz)]), .ghQ = ghQ, .pD = dfyz_t$pD, .q = q) 
# crz. = c(Rz[lower.tri(Rz)]); ghQ. = ghQ; pD. = dfyz_t$pD; q. = q; maxit. = 100
