fyz <- function(Y,ghQ,b,famL){
  A1 <- lapply(1:length(famL), function(i) famL[[i]]$dY(i,Y[,i],b,ghQ))
  A1 <- array(unlist(A1),dim = c(nrow(Y),nrow(ghQ$points),length(famL)))
  A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights)
  return(list(pD = exp(rowSums(A1,dim = 2))/A2, ll = sum(log(A2))))
}

ll <- function(cb,Y,ghQ,bg,famL,info){
  b <- lb2cb(bg)
  b[b != 0] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  f1 <- -d1ll(Y,ghQ,b,famL,info,f0$pD)
  f2 <- -d2ll(Y,ghQ,b,famL,info,f0$pD)
  return(list(value = -f0$ll, gradient = f1, hessian = f2))
}

d1ll <- function(Y,ghQ,b,famL,info,pd){
  sc <- NULL
  w2 <- pd
  w3 <- c(ghQ$weights)
  for(i in 1:ncol(Y)){
    for(k in 1:famL[[i]]$npar){
      w1 <- famL[[i]]$dvY(i,Y[,i],b,ghQ,info)$d1[[k]]
      sc <- c(sc,colSums(ghQ$out[[k]]*colSums(w1*w2)*w3)[b[[k]][i,] != 0])
    }
  }
  return(unname(sc))
}

d2ll <- function(Y,ghQ,b,famL,info,pd){
  w2 <- pd
  w3 <- c(ghQ$weights)
  H <- matrix(0,sum(lb2cb(b) != 0),sum(lb2cb(b) != 0))
  idH <- 0
  for(i in 1:ncol(Y)){
    pY <- seq.int(from = max(idH+1,i), len = sum(lb2mb(b)[i,] != 0))
    h <- matrix(NA,sum(lb2mb(b)[i,] != 0),sum(lb2mb(b)[i,] != 0))
    dv <- famL[[i]]$dvY(i,Y[,i],b,ghQ,info)
    dvP <- dv$d2
    dvO <- dv$d1
    dvC <- dv$dc
    H3 <- NULL
    id0 <- 0
    for(k in 1:famL[[i]]$npar){
      Zk <- as.matrix(ghQ$out[[k]])[,b[[k]][i,] != 0]
      dim(Zk) <- c(nrow(ghQ$points), sum(b[[k]][i,] != 0))
      pk <- seq.int(from = max(id0+1,k), len = ncol(Zk))
      w1 <- dvP[[k]]
      wa <- colSums(w1*w2)*w3
      h1 <- lapply(1:nrow(ghQ$points), function(r) tcrossprod(Zk[r,]) * wa[r])
      h1 <- Reduce("+", h1)
      # ~~~~~~~~~~~
      wo <- dvO[[k]]
      wb <- colSums(wo^2*w2)*w3
      h2 <- lapply(1:nrow(ghQ$points), function(r) tcrossprod(Zk[r,]) * wb[r])
      h2 <- Reduce("+", h2)
      # ~~~~~~~~~~~
      wc <- sweep(wo*w2,2,w3,"*")
      h3 <- array(sapply(1:nrow(w1), function(r) Zk * wc[r,]),
                  dim = c(nrow(ghQ$points),ncol(Zk),nrow(w1)))
      h3 <- abind::abind(H3,h3,along = 2)
      h[pk,pk] <- h1 + h2
      id0 <- max(pk)
      if(tryCatch({!is.null(dvC[[k]])}, error = function(e) F)){
        id1 <- id0
        for(o in 1:length(dvC[[k]])){
          Zo <- as.matrix(ghQ$out[[k + o]])[,b[[k + o]][i,] != 0]
          dim(Zo) <- c(nrow(ghQ$points), sum(b[[k + o]][i,] != 0))
          po <- seq.int(from = max(id1+1,o), len = ncol(Zo))
          w1i <- dvC[[k]][[o]]
          w0i <- sweep(w1i*w2,2,w3,"*")
          wai <- colSums(w0i)
          h1i <- lapply(1:nrow(ghQ$points), function(r) tcrossprod(Zk[r,],Zo[r,]) * wai[r])
          h1i <- Reduce("+", h1i)
          # ~~~~~~~~~~~~~~~~~~
          woi <- dvO[[k + o]]
          wbi <- colSums(wo*woi*w2)*w3
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

ad2ll <- function(Y,ghQ,b,famL,info,pd){
  w2 <- pd
  w3 <- c(ghQ$weights)
  H <- matrix(0,sum(lb2cb(b) != 0),sum(lb2cb(b) != 0))
  idH <- 0
  for(i in 1:ncol(Y)){
    pY <- seq.int(from = max(idH+1,i), len = sum(lb2mb(b)[i,] != 0))
    dv <- famL[[i]]$dvY(i,Y[,i],b,ghQ,info)
    dvO <- dv$d1
    H3 <- NULL
    id0 <- 0
    for(k in 1:famL[[i]]$npar){
      Zk <- as.matrix(ghQ$out[[k]])[,b[[k]][i,] != 0]
      dim(Zk) <- c(nrow(ghQ$points), sum(b[[k]][i,] != 0))
      pk <- seq.int(from = max(id0+1,k), len = ncol(Zk))
      wo <- dvO[[k]]
      wc <- sweep(wo*w2,2,w3,"*")
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

lla <- function(cb,Y,ghQ,bg,famL,info){
  b <- lb2cb(bg)
  b[b != 0] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  return(-f0$ll)
}

d1lla <- function(cb,Y,ghQ,bg,famL,info){
  b <- lb2cb(bg)
  b[b != 0] <- cb
  b <- cb2lb(b,bg)
  f0 <- fyz(Y,ghQ,b,famL)
  f1 <- -d1ll(Y,ghQ,b,famL,info,f0$pD)
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
  
  lvar <- unique(unlist(lapply(1:length(form), function(i) all.vars(form[[i]]))))
  lvar <- grep("Z", lvar, fixed = T, value = T)
  q <- length(lvar)
  Z <- scale(princomp(Y,cor = T)$scores)[,1:q]
  colnames(Z) <- paste0("Z", 1:q)
  
  sZ <- NULL
  bstart <- NULL
  for(r in names(form)){
    sZ[[r]] <- as.data.frame(model.matrix(form[[r]],as.data.frame(Z)))
    bstart[[r]] <- matrix(0, nrow = ncol(Y), ncol = ncol(sZ[[r]]))
    dimnames(bstart[[r]]) <- list(colnames(Y),colnames(sZ[[r]]))
  }
  
  for(i in 1:ncol(Y)){
    eq <- form
    if(famL[[i]]$family == "Beta") tmpY <- y.(Y[,i]) else tmpY <- Y[,i]
    for(p in 1:length(eq)){ eq[[p]] <- update(eq[[p]], tmpY ~ .)  }
    tmp <- gamlss::gamlss(eq$mu, sigma.formula = eq$sigma, tau.formula = eq$tau, nu.formula = eq$nu,
                          family = famL[[i]]$iuse, data = as.data.frame(cbind(tmpY,Z)),
                          control = gamlss::gamlss.control(trace = F))
    for(p in 1:famL[[i]]$npar){ bstart[[p]][i,] <- coef(tmp,famL[[i]]$pars[p]) }
  }
  
  for(r in names(bstart)){
    for(j in 1:q){
      if(bstart[[r]][j, paste0("Z",j)] < 0) bstart[[r]][, paste0("Z",j)] <- -bstart[[r]][, paste0("Z",j)]
    } }
  return(bstart)
}

rmat <- function(res,b){
  rest <- l2rm(res)
  frparm <- b
  for(i in 1:nrow(rest)){
    frparm[[as.character(rest[i,1])]][as.integer(rest[i,2]),as.character(rest[i,3])] <- as.numeric(rest[i,4])
  }
  return(frparm)
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
  if(check.eigen == TRUE){
    n.e.val <- e.val[e.val <= 0]
    s <- sum(e.val[n.e.val])*2  
    t <- s^2*100 + 1
    p <- min(e.val[(e.val <= 0) == FALSE])
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

# Not necessary for package
# ~~~~~~~~~~~~~~~~~~~~~~~~~

glvmlss_parsimE1 <- function(nsim, saveRes = T){
  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  progress <- function(a) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), a)
  opts <- list(progress = progress)
  FCOL <- suppressWarnings(
    foreach(l = 1:nsim,
            .combine = rbind,
            .options.snow = opts,
            .export = ls(parent.frame()) ) %dopar% {
              tryCatch({
                A <- glvmlss_sim(n,famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
                c0 <- Sys.time()
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = A$b)
                c1 <- Sys.time()
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ex2 <- c(c1-c0, B$iter, length(lb2cb(B$b)))
                return(c(cof,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("E1n",n, "p", p,".Rds"))
  return(FCOL)
}

glvmlss_parsimpost <- function(X, p, out = c("MSE","AB"),
                               mu.eq = ~ Z1+Z2, sg.eq = NULL,
                               ta.eq = NULL, nu.eq = NULL){
  
  # str <- as.character(X)
  # p = as.integer(substr(str,(nchar(str)+1)-2,nchar(str))); rm(str)
  environment(prep_form) <- environment()
  form <- prep_form()
  
  nam <- NULL
  for(i in names(form)){ nam <- rbind(nam,expand.grid(i,1:p,stringsAsFactors = F)) }
  nam <- unlist(lapply(1:nrow(nam),function(i) paste0(nam[i,], collapse = "")))
  
  nM <- NULL
  for(i in 1:p){
    for(j in names(form)){
      nM <- append(nM,paste0(nam[grepl(j,as.character(nam))][i],".",c(0,seq_len(length(all.vars(as.formula(form[[j]])))))))
    }
  }; rm(nam)
  
  ip <- !grepl(".0",nM,fixed = T)
  li1 <- c("mu1.2","mu2.1","sigma1.2","sigma2.1","tau1.2","tau2.1","nu1.2","nu2.1")
  li2 <- lapply(li1, function(i) !grepl(i,nM))
  for(i in 1:length(li2)){ ip <- ip & li2[[i]] }; names(ip) <- nM; rm(list=c("li1","li2"))
  ip[grepl(".0",nM,fixed = T)] <- T; nM <- nM[ip]
  
  X  <- X[complete.cases(X),]; X <- X[apply(X,1,min)!=-999, ]
  ix <- mean(X[,ncol(X)])
  X <- X[,-ncol(X)]
  rR <- !(X[,ncol(X)] == 1000)
  Xb0 <- X[rR,seq_len(ix[1])[ip]]; colnames(Xb0) <- nM ;
  Xbe <- X[rR,((ix[1]+1):(2*ix[1]))[ip]]; colnames(Xbe) <- nM
  if("iter" %in% out) Xit <- X[, max((ix[1]+1):(2*ix[1]))+2]
  if("time" %in% out) Xti <- X[, max((ix[1]+1):(2*ix[1]))+1]
  
  imu0 <- grepl("mu",colnames(Xb0),fixed = T) & grepl(".0",colnames(Xb0),fixed = T)
  imul <- grepl("mu",colnames(Xb0),fixed = T) & !grepl(".0",colnames(Xb0),fixed = T)
  isg0 <- grepl("sigma",colnames(Xb0),fixed = T) & grepl(".0",colnames(Xb0),fixed = T)
  isgl <- grepl("sigma",colnames(Xb0),fixed = T) & !grepl(".0",colnames(Xb0),fixed = T)
  ili <- list("mu0" = imu0, "mu1" = imul,"sg0" = isg0, "sg1" = isgl)
  
  res <- nam <- NULL
  for(ii in names(ili)){
    if("MSE" %in% out) assign(paste0("AvMSE",ii), mean(colMeans((Xb0[,ili[[ii]]] - Xbe[,ili[[ii]]])^2)) )
    if("AB" %in% out) assign(paste0("AvAB",ii), mean(abs(colMeans(Xb0[,ili[[ii]]] - Xbe[,ili[[ii]]]))) )
    if("SB" %in% out) assign(paste0("AvSB",ii), mean((colMeans(Xb0[,ili[[ii]]] - Xbe[,ili[[ii]]]))^2) )
    if("RB" %in% out) assign(paste0("AvRB",ii), mean((colMeans(Xb0[,ili[[ii]]] - Xbe[,ili[[ii]]]))/abs(colMeans(Xb0[,ili[[ii]]]))) ) } 
  
  if("MSE" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvMSE",ii))); nam <- c(nam, paste0("AvMSE",ii)) }
  if("AB" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvAB",ii))); nam <- c(nam, paste0("AvAB",ii)) }
  if("SB" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvSB",ii))); nam <- c(nam, paste0("AvSB",ii)) }
  if("RB" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvRB",ii))); nam <- c(nam, paste0("AvRB",ii)) }
  names(res) <- nam
  
  if("iter" %in% out){ nam <- c(names(res),"AvIter"); res <- c(res,mean(Xit)); names(res) <- nam }
  if("time" %in% out){ nam <- c(names(res),"AvTime"); res <- c(res,mean(Xti)); names(res) <- nam }
  return(c(round(res,4)))
}

y. <- function(x){
 if (any(x == 1)) x[x == 1] <- 1 - .Machine$double.eps^(1/3)
 if (any(x == 0)) x[x == 0] <- .Machine$double.eps^1/3
 return(x)
}