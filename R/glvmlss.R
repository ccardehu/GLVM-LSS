glvmlss <- function(data, family = list(),
                    mu.eq = ~ Z1, sg.eq = NULL,
                    nu.eq = NULL, ta.eq = NULL,
                    control = list(), ...){
  
  Y <- prep_Y(data)
  if(!all(Y$na.idx)) cat("\n Argument 'data' has missing values (NA): Handling with full information MML.")
  p <- ncol(Y$Y)
  famL <- prep_fam(family)
  form <- prep_form(mu.eq, sg.eq, nu.eq, ta.eq)
  q <- prep_Z(form)
  control <- prep_cont(control,q,p,...)
  if(control$corr.lv){
    if(!is.null(control$Rz)) Rz <- control$Rz else { Rz <- diag(q); Rz[Rz == 0] <- 0.05 }
  } else { Rz <- diag(q) }
  ghQ <- prep_ghq(nQP = control$nQP, form = form, Rz = Rz)
  b <- prep_stva(control = control,form = form,ghQ = ghQ,Y = Y,p = p,famL = famL,q = q)
  rb <- b$rb   # 'rb' are free parameters! (i.e., T = to be estimated)
  bp <- b$penb # 'bp' are the penalised parameters! (i.e., T = to be penalised)
  b <- b$b

  if(control$penalty == "none" | is.null(control$penalty)){
    environment(glvmlss_fit) <- environment()
    fit <- glvmlss_fit() 
  } else {
    pen.control <- list(penalty = control$penalty,
                        w.alasso = control$w.alasso, a = control$a)
    if(!is.null(control$lambda)){
      if(any("auto" %in% control$lambda)){
        if(length(control$lambda) == 1 && control$lambda == "auto"){
          pen.control$lambda <- rep(sqrt(.Machine$double.eps),length(b)) #1/nrow(Y$Y)
          pen.control$lambda.auto <- rep(T, length(pen.control$lambda))
        }
        if(length(control$lambda) > 1){
          if(length(control$lambda) != length(b)) stop("Define as many values for lambda as location, shape, or scale parameters.")
          pen.control$lambda <- control$lambda
          pen.control$lambda.auto <- pen.control$lambda == "auto"
          pen.control$lambda[pen.control$lambda.auto] <- sqrt(.Machine$double.eps) #1/nrow(Y$Y)
          pen.control$lambda <- as.numeric(pen.control$lambda)
        }
        autoL <- T; cycle <- 0
      } else {
        if(length(control$lambda) != length(b)){
          pen.control$lambda <- rep(control$lambda, length(b))
        } else{ pen.control$lambda <- control$lambda }
        pen.control$lambda.auto <- rep(T, length(pen.control$lambda))
        autoL <- F; cycle <- 0
      }
    } else {
      pen.control$lambda <- rep(sqrt(.Machine$double.eps),length(b)) #1/nrow(Y$Y)
      pen.control$lambda.auto <- rep(T, length(pen.control$lambda))
      autoL <- F; cycle <- 0
    }
    if(!is.null(control$gamma)) pen.control$gamma <- control$gamma else pen.control$gamma <- log(nrow(Y$Y))/2
    pen.control$iter.lim <- control$iter.lim
    environment(glvmlss_penfit) <- environment()
    fit <- glvmlss_penfit()      
  }
  
  fit <- structure(fit, class = "glvmlss")
  if(!is.null(fit$hes$H) | !is.null(fit$hes$Hp)) fit$GAIC <- GAIC(fit) else fit$GAIC <- NULL
  if(!is.null(fit$hes$H) | !is.null(fit$hes$Hp)) fit$GBIC <- GBIC(fit) else fit$GBIC <- NULL
  if(!is.null(fit$GAIC) & !is.null(fit$GBIC)) fit$edf <- (fit$GBIC + 2*fit$unploglik)/log(fit$n)
  
  if(control$f.scores){
    ghQ <- prep_ghq(control$nQP, form, fit$Rz)
    fit$f.scores <- as.data.frame(fyz_c$pD%*%(c(ghQ$weights)*ghQ$points)[,,drop=F]) 
    names(fit$f.scores) <- colnames(ghQ$points) }
  
  return(fit)
}

glvmlss_fit <- function(){
  
  convg <- ifelse(control$EM_iter == 0, T, F)
  pb <- lb2cb(b)
  cb <- pb[lb2cb(rb)]
  info <- control$mat.info
  iter <- 0
  updSS <- vector(mode = "numeric",length = p)
  llk <- numeric(control$EM_iter + 1)
  dfyz_t <- fyz(Y,ghQ,b,famL)
  llk[1] <- -dfyz_t$ll
  EM_appHess <- control$EM_appHess
  posbc <- seq(1:length(lb2cb(b)))[lb2cb(rb)]
  posbM <- lb2mb(cb2lb(seq(1:length(lb2cb(b))), b = b))
  SS <- rep(control$EM_lrate, p)
  
  while(!convg){
    bold1 <- b
    iter <- iter + 1
    pb <- lb2cb(b)
    cb <- pb[lb2cb(rb)]
    d1ll_t <- -d1ll(Y,ghQ,b,famL,info,dfyz_t$pD,rb)
    if(control$EM_use2d){
      if(!EM_appHess){
        da2ll_t <- -d2llEM(Y,ghQ,b,famL,info,dfyz_t$pD,rb)
      } else {
        da2ll_t <- ad2ll(Y,ghQ,b,famL,info,dfyz_t$pD,rb)
      }
      tryCatch({cb <- cb - solve(da2ll_t, d1ll_t)},
               error = function(e){cb <- cb - c(m2pdm(da2ll_t)$inv%*%d1ll_t)})
      pb[lb2cb(rb)] <- cb
      b <- cb2lb(pb,b)
      if(!is.null(control$iden.res) && control$iden.res == "orthogonal") b <- fixOrthob(b)
      if(control$corr.lv){
        Rz <- newRz(Rz = Rz, ghQ = ghQ, pD = dfyz_t$pD, q = q, control)
        ghQ <- prep_ghq(control$nQP, form, Rz)
      }
      dfyz_t <- fyz(Y,ghQ,b,famL)
      llk[iter+1] <- -dfyz_t$ll
      eps1 <- llk[iter+1] - llk[iter]
      if(eps1 > 0){
        b <- bold1
        pb <- lb2cb(b)
        cb <- pb[lb2cb(rb)]
        iter <- iter - 1
        break
      }
      convg <- ifelse(abs(eps1) < control$tol | iter == control$EM_iter, T, F)
      if(control$verbose){
        if(iter == 1) cat(paste0("\n EM iter: ",iter,", marginal loglik.: ", round(-llk[iter+1],5))) else
          cat(paste0("\r EM iter: ",iter,", marginal loglik.: ", round(-llk[iter+1],5)))
      }
    } else {
      eps = rep(1,p)        
      for(ii in 1:p){
        cbold <- cb
        while(eps[ii] > 0){
          if(SS[ii] == sqrt(.Machine$double.eps)){eps[ii] = 0; break}
          idx <- posbc %in% posbM[ii,]
          cb[idx] <- cb[idx] - SS[ii]*d1ll_t[idx]
          pb[lb2cb(rb)] <- cb
          b <- cb2lb(pb,b)
          if(!is.null(control$iden.res) && control$iden.res == "orthogonal") b <- fixOrthob(b)
          dfyz_t <- fyz(Y,ghQ,b,famL)
          llktmp <- -dfyz_t$ll
          eps[ii] <- llktmp - llk[iter]
          if(eps[ii] > 0){
            cb[idx] <- cbold[idx]
            SS[ii] <- max(SS[ii]/2, sqrt(.Machine$double.eps))
            updSS[ii] <- updSS[ii] + 1
          }
        }
      }
      if(control$corr.lv){
        Rz <- newRz(Rz = Rz, ghQ = ghQ, pD = dfyz_t$pD, q = q, control)
        ghQ <- prep_ghq(control$nQP, form, Rz)
      }
      dfyz_t <- fyz(Y,ghQ,b,famL)
      llk[iter+1] <- -dfyz_t$ll
      eps1 <- llk[iter+1] - llk[iter]
      if(eps1 > 0){
        b <- bold1
        pb <- lb2cb(b)
        cb <- pb[lb2cb(rb)]
        iter <- iter - 1
        eps1 <- 0
      }
      convg <- ifelse(abs(eps1) < control$tol | iter == control$EM_iter, T, F)
      if(control$verbose & (iter == 1 | iter %% 10 == 0)){
        if(iter == 1) cat(paste0("\n EM iter: ",iter,", marginal loglik.: ", round(-llk[iter+1],5), " (# of step size updates: ",sum(updSS),")")) else
          cat(paste0("\r EM iter: ",iter,", marginal loglik.: ", round(-llk[iter+1],5), " (# of step size updates: ",sum(updSS),")"))
      }   
    }
  }

  if(control$solver == "trust"){
    if(control$verbose) cat("\n Direct MML estimation (trust-region) ...")
    if(control$lazytrust) trustll <- lazyll else trustll <- ll
    cb <- trust::trust(objfun = trustll, parinit = cb, rinit = 1, rmax = 5,
                       fterm = control$tol, iterlim = control$iter.lim,
                       Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info, rb = rb)
    pb[lb2cb(rb)] <- cb$argument
    b <- cb2lb(pb,b)
    if(!is.null(control$iden.res) && control$iden.res == "orthogonal") b <- fixOrthob(b)
    convg <- ifelse(cb$converged,1,0)
    if(control$corr.lv){
      dfyz_t <- fyz(Y,ghQ,b,famL)
      Rz <- newRz(Rz = Rz, ghQ = ghQ, pD = dfyz_t$pD, q = q, control)
      ghQ <- prep_ghq(control$nQP,form,Rz) }
    if(control$verbose) cat(paste0("\r Direct MML estimation (trust-region) ... Converged after ",
                                   iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " ", control$solver, ")",
                                   "\n (Approximate) Marginal loglikelihood: ", round(-cb$value,5))) }
  
  if(control$solver == "nlminb"){
    if(control$verbose) cat("\n Direct MML estimation (nlminb) ...")
    cb <- nlminb(start = cb, objective = lla, gradient = d1lla, control = list(iter.max = control$iter.lim),
                 Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info, rb = rb)
    pb[lb2cb(rb)] <- cb$par
    b <- cb2lb(pb,b)
    if(!is.null(control$iden.res) && control$iden.res == "orthogonal") b <- fixOrthob(b)
    convg <- ifelse(cb$convergence == 0,1,0)
    if(control$corr.lv){
      dfyz_t <- fyz(Y,ghQ,b,famL)
      Rz <- newRz(Rz = Rz, ghQ = ghQ, pD = dfyz_t$pD, q = q, control)
      ghQ <-prep_ghq(control$nQP,form,Rz) }
    if(control$verbose) cat(paste0("\r Direct MML estimation (nlminb) ... Converged after ",
                                   iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " nlminb)",
                                   "\n (Approximate) Marginal loglikelihood: ", round(-cb$objective,5))) } 
  
  if(control$solver != "nlminb" & control$solver != "trust"){
    if(control$verbose) cat(paste0("\n Direct MML estimation (", control$solver,") ..."))
    cb <- optim(cb, fn = lla, gr = d1lla, method = control$solver,
                control = list(maxit = control$iter.lim, fnscale = 1),
                Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info, rb = rb)
    pb[lb2cb(rb)] <- cb$par
    b <- cb2lb(pb,b)
    if(!is.null(control$iden.res) && control$iden.res == "orthogonal") b <- fixOrthob(b)
    convg <- ifelse(cb$convergence == 0,1,0)
    if(control$corr.lv){
      dfyz_t <- fyz(Y,ghQ,b,famL)
      Rz <- newRz(Rz = Rz, ghQ = ghQ, pD = dfyz_t$pD, q = q, control)
      ghQ <- prep_ghq(control$nQP,form,Rz) }
    if(control$verbose) cat(paste0("\r Direct MML estimation (", control$solver, ") ... Converged.",
                                   "\n (Approximate) Marginal loglikelihood: ", round(-cb$value,5))) }
  
  assign(x = "fyz_c",value = fyz(Y,ghQ,b,famL),envir = parent.frame())
  if(!is.null(control$est.ci)){
    if(control$est.ci == "Standard" | control$est.ci == "Bayesian"){
      if(control$verbose) cat("\n Computing standard errors ...")
      hes_u <- -d2ll(Y = Y,ghQ = ghQ,b = b, famL = famL, info = "Fisher", pd = fyz_c$pD, rb = rb); hes_c <- NULL
      if(control$corr.lv){
        hes_R <- hessRz(Rz,Y,b,famL,form,q,control)
        hes_u <- magic::adiag(hes_u,-hes_R$hess)
      }
      pdMhes <- m2pdm(hes_u)
      if(control$verbose){
        if(!pdMhes$is.PD){ 
          cat(paste0("\r Computing standard errors ... done!\n \r Warning: Fisher Information matrix is not positive definite at solution (fixed).\n "))
        } else { 
          cat(paste0("\r Computing standard errors ... done!\n ")) 
        } 
      }
      seb <- rep(NA,length(lb2cb(b)))
      seb[lb2cb(rb)] <- sqrt(diag(pdMhes$inv.mat))[1:sum(lb2cb(rb))]
      seb <- cb2lb(seb,b)
      if(control$corr.lv){
        iRz <- sort(ncol(hes_u):(ncol(hes_u)-sum(hes_R$rR)+1))
        seRz <- array(NA, dim = dim(hes_R$rR))
        seRz[hes_R$rR] <- sqrt(diag(pdMhes$inv.mat))[iRz]
      } else seRz <- NULL
    }
    if(control$est.ci == "Approximate"){
      if(control$verbose) cat("\n Computing (approximate) standard errors ...")
      hes_u <- ad2ll(Y = Y,ghQ = ghQ,b = b, famL = famL, info = "Fisher", pd = fyz_c$pD, rb = rb); hes_c <- NULL
      if(control$corr.lv){
        hes_R <- hessRzEM(Rz,ghQ,fyz_c$pD,q,control)
        hes_u <- magic::adiag(hes_u,-hes_R$hess)
      }
      pdMhes <- m2pdm(hes_u)
      if(control$verbose){
        if(!pdMhes$is.PD){
          cat(paste0("\r Computing (approximate) standard errors ... done!\n \r Warning: Fisher Information matrix is not positive definite at solution (fixed).\n "))
        } else {
          cat(paste0("\r Computing (approximate) standard errors ... done!\n "))
        }
      }
      seb <- rep(NA,length(lb2cb(b)))
      seb[lb2cb(rb)] <- sqrt(diag(pdMhes$inv.mat))[1:sum(lb2cb(rb))]
      seb <- cb2lb(seb,b)
      if(control$corr.lv){
        iRz <- sort(ncol(hes_u):(ncol(hes_u)-sum(hes_R$rR)+1))
        seRz <- array(NA, dim = dim(hes_R$rR))
        seRz[hes_R$rR] <- sqrt(diag(pdMhes$inv.mat))[iRz]
      } else seRz <- NULL
    }
  } else seb <- hes_u <- hes_c <- seRz <- NULL
  
  return(list(b = b, Rz = Rz, loglik = fyz_c$ll, unploglik = fyz_c$ll,
              convergence = convg, iter = iter + cb$iter,
              hes = list(H = hes_u, Hp = hes_c), SE = list(b = seb, Rz = seRz), n = nrow(Y$Y)))
}

glvmlss_penfit <- function(){
  
  penON <- T
  lhist <- matrix(NA,nrow = control$autoL_iter+1, ncol = length(pen.control$lambda))
  ssehist <- numeric(control$autoL_iter+1)
  rownames(lhist) <- names(ssehist) <- paste0("c",c(1:(control$autoL_iter+1)))
  colnames(lhist) <- paste0("lambda_",names(b))
  lhist[1,] <- pen.control$lambda
  
  while(penON){
    
    convg <- ifelse(control$EM_iter == 0, T, F)
    pb <- lb2cb(b)
    cb <- pb[lb2cb(rb)]
    info <- control$mat.info
    iter <- 0
    updSS <- vector(mode = "numeric",length = p)
    llk <- numeric(control$EM_iter + 1)
    dfyz_t <- fyz(Y,ghQ,b,famL)
    pMat <- pM(b,rb,bp,penalty = pen.control$penalty,
               lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
    llk[1] <- -dfyz_t$ll + 0.5*nrow(Y$Y)*crossprod(cb,pMat$full)%*%cb
    EM_appHess <- control$EM_appHess
    posbc <- seq(1:length(lb2cb(b)))[lb2cb(rb)]
    posbM <- lb2mb(cb2lb(seq(1:length(lb2cb(b))), b = b))
    SS <- rep(control$EM_lrate, p)
    
    while(!convg){
      bold1 <- b
      iter <- iter + 1
      pb <- lb2cb(b)
      cb <- pb[lb2cb(rb)]
      pMat <- pM(b,rb,bp,penalty = pen.control$penalty,
                 lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
      d1ll_t <- -d1ll(Y,ghQ,b,famL,info,dfyz_t$pD,rb) + c(nrow(Y$Y)*crossprod(cb,pMat$full))
      if(control$EM_use2d){
        if(!EM_appHess){
          da2ll_t <- -d2llEM(Y,ghQ,b,famL,info,dfyz_t$pD,rb) + nrow(Y$Y)*pMat$full
        } else {
          da2ll_t <- -ad2ll(Y,ghQ,b,famL,info,dfyz_t$pD,rb) + nrow(Y$Y)*pMat$full
        }
        tryCatch({cb <- cb - solve(da2ll_t, d1ll_t)},
                 error = function(e){cb <- cb - c(m2pdm(da2ll_t)$inv%*%d1ll_t)})
        if(cycle == 0){
          da2ll_t <- -d2llEM(Y,ghQ,b,famL,info,dfyz_t$pD,rb) + nrow(Y$Y)*pMat$full
          ssehist[1] <- SSE(loglambda = log(pen.control$lambda[pen.control$lambda.auto]),
                            b = b, gra = -d1ll_t, hes = da2ll_t, rb = rb, bp = bp, pml = pen.control, Y = Y)
        }
        pb[lb2cb(rb)] <- cb
        b <- cb2lb(pb,b)
        if(!is.null(control$iden.res) && control$iden.res == "orthogonal") b <- fixOrthob(b)
        if(control$corr.lv){
          Rz <- newRz(Rz = Rz, ghQ = ghQ, pD = dfyz_t$pD, q = q, control)
          ghQ <- prep_ghq(control$nQP,form,Rz)
        }
        dfyz_t <- fyz(Y,ghQ,b,famL)
        pMat <- pM(b,rb,bp,penalty = pen.control$penalty,
                   lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
        llk[iter+1] <- -dfyz_t$ll + 0.5*nrow(Y$Y)*crossprod(cb,pMat$full)%*%cb
        eps1 <- llk[iter+1] - llk[iter]
        if(eps1 > 0){
          b <- bold1
          pb <- lb2cb(b)
          cb <- pb[lb2cb(rb)]
          iter <- iter - 1
          break
        }
        convg <- ifelse(abs(eps1) < control$tol | iter == control$EM_iter, T, F)
        if(control$verbose){
          if(!autoL){
            if(iter == 1) cat(paste0("\n EM iter: ",iter,", penalised marginal loglik.: ", round(-llk[iter+1],5))) else
              cat(paste0("\r EM iter: ",iter,", penalised marginal loglik.: ", round(-llk[iter+1],5)))
          } else {
            if(iter == 1) cat(paste0("\n [Cycle: ",cycle + 1,"] EM iter: ",iter,", penalised marginal loglik.: ", round(-llk[iter+1],5))) else
              cat(paste0("\r [Cycle: ",cycle + 1,"] EM iter: ",iter,", penalised marginal loglik.: ", round(-llk[iter+1],5)))
          }
        }
        
      } else {
        if(cycle == 0){
          da2ll_t <- -d2llEM(Y,ghQ,b,famL,info,dfyz_t$pD,rb) + nrow(Y$Y)*pMat$full
          ssehist[1] <- SSE(loglambda = log(pen.control$lambda[pen.control$lambda.auto]),
                            b = b, gra = -d1ll_t, hes = da2ll_t, rb = rb, bp = bp, pml = pen.control, Y = Y)
        }
        eps = rep(1,p)        
        for(ii in 1:p){
          cbold <- cb
          while(eps[ii] > 0){
            if(SS[ii] == sqrt(.Machine$double.eps)){eps[ii] = 0; break}
            idx <- posbc %in% posbM[ii,]
            cb[idx] <- cb[idx] - SS[ii]*d1ll_t[idx]
            pb[lb2cb(rb)] <- cb
            b <- cb2lb(pb,b)
            if(!is.null(control$iden.res) && control$iden.res == "orthogonal") b <- fixOrthob(b)
            dfyz_t <- fyz(Y,ghQ,b,famL)
            pMat <- pM(b,rb,bp,penalty = pen.control$penalty,
                       lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
            llktmp <- -dfyz_t$ll + 0.5*nrow(Y$Y)*crossprod(cb,pMat$full)%*%cb
            eps[ii] <- llktmp - llk[iter]
            if(eps[ii] > 0){
              cb[idx] <- cbold[idx]
              SS[ii] <- max(SS[ii]/2, sqrt(.Machine$double.eps))
              updSS[ii] <- updSS[ii] + 1
            }
          }
        }
        if(control$corr.lv){
          Rz <- newRz(Rz = Rz, ghQ = ghQ, pD = dfyz_t$pD, q = q, control)
          ghQ <- prep_ghq(control$nQP, form, Rz)
        }
        dfyz_t <- fyz(Y,ghQ,b,famL)
        pMat <- pM(b,rb,bp,penalty = pen.control$penalty,
                   lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
        llk[iter+1] <- -dfyz_t$ll + 0.5*nrow(Y$Y)*crossprod(cb,pMat$full)%*%cb
        eps1 <- llk[iter+1] - llk[iter]
        if(eps1 > 0){
          b <- bold1
          pb <- lb2cb(b)
          cb <- pb[lb2cb(rb)]
          iter <- iter - 1
          eps1 <- 0
        }
        convg <- ifelse(abs(eps1) < control$tol || iter == control$EM_iter, T, F)
        if(control$verbose & (iter == 1 | iter %% 10 == 0)){
          if(!autoL){
            if(iter == 1) cat(paste0("\n EM iter: ",iter,", penalised marginal loglik.: ", round(-llk[iter+1],5), " (# of step size updates: ",sum(updSS),")")) else
              cat(paste0("\r EM iter: ",iter,", penalised marginal loglik.: ", round(-llk[iter+1],5), " (# of step size updates: ",sum(updSS),")"))
          } else {
            if(iter == 1) cat(paste0("\n [Cycle: ",cycle + 1,"] EM iter: ",iter,", penalised marginal loglik.: ", round(-llk[iter+1],5), " (# of step size updates: ",sum(updSS),")")) else
              cat(paste0("\r [Cycle: ",cycle + 1,"] EM iter: ",iter,", penalised marginal loglik.: ", round(-llk[iter+1],5), " (# of step size updates: ",sum(updSS),")"))
          }
        }
      }
    }
    
    if(control$solver == "trust"){
      if(control$verbose){
        if(!autoL) cat("\n Direct penalised MML estimation (trust-region) ...") else
          if(cycle == 0) cat(paste0("\n [Cycle: ",cycle + 1,"] Direct penalised MML estimation (trust-region) ...")) else
            cat(paste0("\n [Cycle: ",cycle + 1,"] Direct penalised MML estimation (trust-region) ...")) }
      if(control$lazytrust) trustll <- lazypll else trustll <- pll
      cb <- trust::trust(objfun = trustll, parinit = cb, rinit = 1, rmax = 5,
                         fterm = control$tol, iterlim = control$iter.lim,
                         Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info, rb = rb, bp = bp,
                         pen.control = pen.control)
      pb[lb2cb(rb)] <- cb$argument
      b <- cb2lb(pb,b)
      if(!is.null(control$iden.res) && control$iden.res == "orthogonal") b <- fixOrthob(b)
      convg <- ifelse(cb$converged,1,0)
      if(control$corr.lv){
        dfyz_t <- fyz(Y,ghQ,b,famL)
        Rz <- newRz(Rz = Rz, ghQ = ghQ, pD = dfyz_t$pD, q = q, control)
        ghQ <- prep_ghq(control$nQP,form,Rz) }
      if(control$verbose){
        if(autoL){ cat(paste0("\r [Cycle: ",cycle + 1,"] Direct penalised MML estimation (trust-region) ... Converged after ",
                              iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " ", control$solver, "),",
                              " pen. margllk: ", round(-cb$value,5)))
        } else {
          cat(paste0("\r Direct penalised MML estimation (trust-region) ... Converged after ",
                     iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " ", control$solver, ")",
                     "\n (Approximate) Penalised Marginal loglikelihood: ", round(-cb$value,5))) } } }
    
    if(control$solver == "nlminb"){
      if(control$verbose){
        if(!autoL) cat("\n Direct penalised MML estimation ('nlminb') ...") else 
          if(cycle == 0) cat(paste0("\n [Cycle: ",cycle + 1,"] Direct penalised MML estimation ('nlminb') ...")) else
            cat(paste0("\n [Cycle: ",cycle + 1,"] Direct penalised MML estimation ('nlminb') ...")) }
      cb <- nlminb(start = cb, objective = plla, gradient = d1plla, control = list(iter.max = control$iter.lim),
                   Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info, rb = rb, bp = bp,
                   pen.control = pen.control)
      pb[lb2cb(rb)] <- cb$par
      b <- cb2lb(pb,b)
      if(!is.null(control$iden.res) && control$iden.res == "orthogonal") b <- fixOrthob(b)
      convg <- ifelse(cb$convergence == 0,1,0)
      if(control$corr.lv){
        dfyz_t <- fyz(Y,ghQ,b,famL)
        Rz <- newRz(Rz = Rz, ghQ = ghQ, pD = dfyz_t$pD, q = q, control)
        ghQ <- prep_ghq(control$nQP,form,Rz) }
      if(control$verbose){
        if(!autoL){ cat(paste0("\r Direct penalised MML estimation ('nlminb') ... Converged after ",
                               iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " nlminb)",
                               "\n (Approximate) Penalised Marginal loglikelihood: ", round(-cb$objective,5)))
        } else {
          cat(paste0("\r [Cycle: ",cycle + 1,"] Direct penalised MML estimation ('nlminb') ... Converged after ",
                     iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " nlminb),",
                     " pen. margllk: ", round(-cb$objective,5))) } } }
    
    if(control$solver != "trust" & control$solver != "nlminb"){
      if(control$verbose){
        if(!autoL) cat(paste0("\n Direct penalised MML estimation (", control$solver,") ...")) else
          if(cycle == 0) cat(paste0("\n [Cycle: ",cycle + 1,"] Direct penalised MML estimation (", control$solver,") ...")) else
            cat(paste0("\n [Cycle: ",cycle + 1,"] Direct penalised MML estimation (", control$solver,") ...")) }
      cb <- optim(cb, fn = plla, gr = d1plla, method = control$solver,
                  control = list(maxit = control$iter.lim, fnscale = 1),
                  Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info, rb = rb, bp = bp,
                  pen.control = pen.control)
      pb[lb2cb(rb)] <- cb$par
      b <- cb2lb(pb,b)
      if(!is.null(control$iden.res) && control$iden.res == "orthogonal") b <- fixOrthob(b)
      convg <- ifelse(cb$convergence == 0,1,0)
      if(control$corr.lv){
        dfyz_t <- fyz(Y,ghQ,b,famL)
        Rz <- newRz(Rz = Rz, ghQ = ghQ, pD = dfyz_t$pD, q = q, control)
        ghQ <- prep_ghq(control$nQP,form,Rz) }
      if(control$verbose){
        if(!autoL){ cat(paste0("\r Direct penalised MML estimation (", control$solver, ") ... Converged.",
                               "\n (Approximate) Penalised Marginal loglikelihood: ", round(-cb$value,5)))
        } else {
          cat(paste0("\r [Cycle: ",cycle + 1,"] Direct penalised MML estimation (", control$solver,") ... Converged (",
                     "pen. margllk: ", round(-cb$value,5), ")")) } } }
    
    if(autoL){
      tObj <- op.lambda(Y,ghQ,b,famL,info,rb,bp,pen.control, (control$EM_appHess | control$est.ci == "Approximate") )
      eps2 <- max(abs(pen.control$lambda - tObj$lambda))
      if(eps2 < control$tol) autoL <- F
      cycle <- cycle + 1
      lhist[cycle+1,]  <- tObj$lambda
      ssehist[cycle+1] <- tObj$sse
      if(ssehist[cycle+1] - ssehist[cycle] > 0){ 
        if(control$verbose) cat(paste0("\n \n SSE did not improve from cycle ", cycle,
                                       " to ", cycle+1,". Stopping automatic selection of tuning parameters in cycle ",cycle," (marked by *)."))
        autoL <- F
        rownames(lhist)[cycle] <- names(ssehist)[cycle] <- paste0("c",cycle,"*")
        cycle = cycle - 1
      } else { pen.control$lambda <- tObj$lambda }
    }
    
    if(!autoL| cycle == control$autoL_iter) penON <- F
    
  }
  
  rb4se <- fixrb(rb,b,control$tolb)
  pMat <- pM(b,rb4se,bp,penalty = pen.control$penalty,
             lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
  b. <- lb2cb(b)[lb2cb(rb4se)]
  assign(x = "fyz_c",value = fyz(Y,ghQ,b,famL),envir = parent.frame())
  
  if(!is.null(control$est.ci)){
    if(control$est.ci == "Standard"){
      if(control$verbose) cat("\n Computing standard errors ...")
      hes_u <- -d2ll(Y,ghQ,b,famL,"Fisher", pd = fyz_c$pD, rb4se)
      hes_c <- hes_u + nrow(Y$Y)*pMat$full
      if(control$corr.lv){
        hes_R <- hessRz(Rz,Y,b,famL,form,q,control)
        hes_u <- magic::adiag(hes_u,-hes_R$hess)
        hes_c <- magic::adiag(hes_c,-hes_R$hess)
      }
      pdMhes <- m2pdm(hes_c)
      if(control$verbose){
        if(!pdMhes$is.PD){ cat(paste0("\r Computing standard errors ... done!\n \r Warning: Fisher Information matrix is not positive definite at solution (fixed).\n "))
        } else {
          cat(paste0("\r Computing standard errors ... done!\n "))
        }
      }
      seb <- rep(NA,length(lb2cb(b)))
      seb[lb2cb(rb4se)] <- sqrt(diag(pdMhes$inv.mat%*%hes_u%*%pdMhes$inv.mat))[1:sum(lb2cb(rb4se))]
      seb <- cb2lb(seb,b)
      if(control$corr.lv){
        iRz <- sort(ncol(hes_u):(ncol(hes_u)-sum(hes_R$rR)+1))
        seRz <- array(NA, dim = dim(hes_R$rR))
        seRz[hes_R$rR] <- sqrt(diag(pdMhes$inv.mat%*%hes_u%*%pdMhes$inv.mat))[iRz]
      } else seRz <- NULL
    }
    if(control$est.ci == "Bayesian"){
      if(control$verbose) cat("\n Computing (Bayesian) standard errors ...")
      hes_u <- -d2ll(Y,ghQ,b,famL,"Fisher", pd = fyz_c$pD, rb4se)
      hes_c <- hes_u + nrow(Y$Y)*pMat$full
      if(control$corr.lv){
        hes_R <- hessRz(Rz,Y,b,famL,form,q,control)
        hes_u <- magic::adiag(hes_u,-hes_R$hess)
        hes_c <- magic::adiag(hes_c,-hes_R$hess)
      }
      pdMhes <- m2pdm(hes_c)
      if(control$verbose){
        if(!pdMhes$is.PD){ cat(paste0("\r Computing (Bayesian) standard errors ... done!\n \r Warning: Fisher Information matrix is not positive definite at solution (fixed).\n "))
        } else {
          cat(paste0("\r Computing (Bayesian) standard errors ... done!\n "))
        }
      }
      seb <- rep(NA,length(lb2cb(b)))
      seb[lb2cb(rb4se)] <- sqrt(diag(pdMhes$inv.mat))[1:sum(lb2cb(rb4se))]
      seb <- cb2lb(seb,b)
      if(control$corr.lv){
        iRz <- sort(ncol(hes_u):(ncol(hes_u)-sum(hes_R$rR)+1))
        seRz <- array(NA, dim = dim(hes_R$rR))
        seRz[hes_R$rR] <- sqrt(diag(pdMhes$inv.mat))[iRz]
      } else seRz <- NULL
    }
    if(control$est.ci == "Approximate"){
      if(control$verbose) cat("\n Computing (approximate) standard errors ...")
      hes_u <- ad2ll(Y,ghQ,b,famL,"Fisher", pd = fyz_c$pD, rb4se)
      hes_c <- hes_u + nrow(Y$Y)*pMat$full
      if(control$corr.lv){
        hes_R <- hessRzEM(Rz,ghQ,fyz_c$pD,q,control)
        hes_u <- magic::adiag(hes_u,-hes_R$hess)
        hes_c <- magic::adiag(hes_c,-hes_R$hess)
      }
      pdMhes <- m2pdm(hes_c)
      if(control$verbose){
        if(!pdMhes$is.PD){
          cat(paste0("\r Computing (approximate) standard errors ... done!\n \r Warning: Fisher Information matrix is not positive definite at solution (fixed).\n "))
        } else {
          cat(paste0("\r Computing (approximate) standard errors ... done!\n "))
        }
      }
      seb <- rep(NA,length(lb2cb(b)))
      seb[lb2cb(rb4se)] <- sqrt(diag(pdMhes$inv.mat))[1:sum(lb2cb(rb4se))]
      seb <- cb2lb(seb,b)
      if(control$corr.lv){
        iRz <- sort(ncol(hes_u):(ncol(hes_u)-sum(hes_R$rR)+1))
        seRz <- array(NA, dim = dim(hes_R$rR))
        seRz[hes_R$rR] <- sqrt(diag(pdMhes$inv.mat))[iRz]
      } else seRz <- NULL
    }
  } else seb <- hes_u <- hes_c <- seRz <- NULL
  
  return(list(b = b, Rz = Rz, loglik = c(fyz_c$ll - 0.5*nrow(Y$Y)*crossprod(b.,pMat$full)%*%b.),
              unploglik = fyz_c$ll, convergence = convg,
              iter = iter + cb$iter, hes = list(H = hes_u, Hp = hes_c), SE = list(b = seb, Rz = seRz),
              cycle = cycle + 1, gamma = pen.control$gamma, n = nrow(Y$Y),
              sse = ssehist[ssehist != 0],
              lambda = lhist[complete.cases(lhist),],
              lambdahat = lhist[cycle + 1,,drop = F]))
}

glvmlss_sim <- function(n, family = list(),
                        mu.eq = ~ Z1 + Z2, sg.eq = NULL,
                        nu.eq = NULL, ta.eq = NULL,
                        control = list(), ...){
  
  famL <- prep_fam(family)
  form <- prep_form(mu.eq, sg.eq, nu.eq, ta.eq)
  p <- length(famL)
  q <- prep_Z(form)
  control <- sim_cont(control,q,...)
  Z <- sim_Z(control,q,n,form)
  b <- sim_stva(control,p,Z,form,q)$b
  if(!is.null(control$iden.res) && control$iden.res == "orthogonal") b <- fixOrthob(b)
  
  Y <- sapply(1:p, function(i) famL[[i]]$sf(i,n,b,Z$Zmod))
  Y <- as.data.frame(Y); colnames(Y) <- paste0("Y", 1:p)
  
  return(list(Y = Y, Z = Z$Z, b = b, Rz = control$Rz))
}

