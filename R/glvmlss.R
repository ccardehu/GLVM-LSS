glvmlss <- function(data, family = list(),
                    mu.eq = ~ Z1 + Z2, sg.eq = NULL,
                    ta.eq = NULL, nu.eq = NULL,
                    control = list(), pen.control = list(), ...){
  
  Y <- data; if(!is.matrix(Y)) Y <- as.matrix(Y)
  p <- ncol(Y)
  environment(prep_fam) <- environment(prep_form) <- environment()
  famL <- prep_fam()
  form <- prep_form()
  environment(prep_Z) <- environment()
  q <- prep_Z()
  environment(prep_cont) <- environment(prep_ghq) <- environment()
  environment(prep_stva) <- environment()
  control <- prep_cont()
  ghQ <- prep_ghq()
  b <- prep_stva()
  rb <- b$rb; b <- b$b # 'rb' are free parameters!
  
  if(control$penalty == "none"){
    environment(glvmlss_fit) <- environment()
    fit <- glvmlss_fit(...) 
  } else {
    pen.control <- list(penalty = control$penalty,
                        w.alasso = control$w.alasso, a = control$a)
    if(!is.null(control$lambda)){ 
      if(any("auto" %in% control$lambda)){ 
        pen.control$lambda <- rep(1/nrow(Y),length(b))
        autoL <- T; cycle <- 0
      } else {
        if(length(control$lambda) != length(b)){
          pen.control$lambda <- rep(control$lambda, length(b)) 
        } else{ pen.control$lambda <- control$lambda }
        autoL <- F; cycle <- 0
      }
    } else { pen.control$lambda <- rep(1/nrow(Y),length(b)); autoL <- F; cycle <- 0 }
    if(!is.null(control$gamma)) pen.control$gamma <- control$gamma else pen.control$gamma <- log(nrow(Y))/2
    pen.control$iter.lim <- control$iter.lim
    environment(glvmlss_penfit) <- environment()
    fit <- glvmlss_penfit(...)      
  }
  
  structure(fit, class = "glvmlss")
}

glvmlss_fit <- function(...){
  
  convg <- ifelse(control$EM_iter == 0, T, F)
  pb <- lb2cb(b)
  cb <- pb[lb2cb(rb) == T]
  info <- control$mat.info
  iter <- 0
  llk <- numeric(control$EM_iter + 1)
  dfyz_t <- fyz(Y,ghQ,b,famL)
  llk[1] <- -dfyz_t$ll
  EM_appHess <- control$EM_appHess
  
  while(!convg){
    iter <- iter + 1
    pb <- lb2cb(b)
    cb <- pb[lb2cb(rb) == T]
    d1ll_t <- -d1ll(Y,ghQ,b,famL,info,dfyz_t$pD,rb)
    if(control$EM_use2d){
      if(!EM_appHess){
        da2ll_t <- -d2llEM(Y,ghQ,b,famL,info,dfyz_t$pD,rb)
      } else {
        da2ll_t <- ad2ll(Y,ghQ,b,famL,info,dfyz_t$pD,rb)
      }
      cb <- cb - solve(da2ll_t, d1ll_t)
    } else {
      cb <- cb - control$EM_lrate*exp(-0.5*iter)*d1ll_t }
    pb[lb2cb(rb) == T] <- cb
    b <- cb2lb(pb,b)
    dfyz_t <- fyz(Y,ghQ,b,famL)
    llk[iter+1] <- -dfyz_t$ll
    eps <- llk[iter+1] - llk[iter]
    if(eps > 0) break
    convg <- ifelse(abs(eps) < control$tol*1e4 || iter == control$EM_iter, T, F)
    if(control$verbose){ if(iter == 1) cat(paste0("\n EM iter: ",iter,", marginal loglik.: ", round(-llk[iter+1],5))) else
      cat(paste0("\r EM iter: ",iter,", marginal loglik.: ", round(-llk[iter+1],5))) }
  }
  
  if(control$solver == "trust"){
    if(control$verbose) cat("\n Direct MML estimation using trust-region algorithm ...")
    cb <- trust::trust(objfun = ll, parinit = cb, rinit = 1, rmax = 5,
                       fterm = control$tol, iterlim = control$iter.lim,
                       Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info, rb = rb)
    pb[lb2cb(rb) == T] <- cb$argument
    b <- cb2lb(pb,b)
    if(control$verbose) cat(paste0("\r Direct MML estimation using trust-region algorithm ... Converged after ",
                                   iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " ", control$solver, ")",
                                   "\n (Approximate) Marginal loglikelihood: ", round(-cb$value,5))) }
  if(control$solver == "nlminb"){
    if(control$verbose) cat("\n Direct MML estimation using 'nlminb' function ...")
    cb <- nlminb(start = cb, objective = lla, gradient = d1lla, control = list(iter.max = control$iter.lim),
                 Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info, rb = rb)
    pb[lb2cb(rb) == T] <- cb$par
    b <- cb2lb(pb,b)
    if(control$verbose) cat(paste0("\r Direct MML estimation using 'nlminb' function ... Converged after ",
                                   iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " nlminb)",
                                   "\n (Approximate) Marginal loglikelihood: ", round(-cb$objective,5))) } 
  if(control$solver != "nlminb" & control$solver != "trust"){
    if(control$verbose) cat(paste0("\n Direct MML estimation using ", control$solver," algorithm ..."))
    cb <- optim(cb, fn = lla, gr = d1lla, method = control$solver,
                control = list(maxit = control$iter.lim, fnscale = 1),
                Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info, rb = rb)
    pb[lb2cb(rb) == T] <- cb$par
    b <- cb2lb(pb,b)
    if(control$verbose) cat(paste0("\r Direct MML estimation using ", control$solver, " algorithm ... Converged.",
                                   "\n (Approximate) Marginal loglikelihood: ", round(-cb$value,5))) }
  
  for(r in names(b)){
    for(j in which(grepl("Z",colnames(b[[r]])))){
      c <- as.integer(substr(colnames(b[[r]])[j],nchar(colnames(b[[r]])[j]),nchar(colnames(b[[r]])[j])))
      if(b[[r]][c,j] < 0) b[[r]][,j] <- -b[[r]][,j]
    } }
    
  fyz_c <- fyz(Y,ghQ,b,famL)
  if(control$est.ci){
    hes_u <- -d2ll(Y,ghQ,b,famL,info, pd = fyz_c$pD, rb); hes_c <- NULL
    if(!m2pdm(hes_u)$is.PD & control$verbose) cat(paste0("\n Unstable solution: Hessian matrix is not positive definite at solution."))
  } else if(control$solver == "trust"){ hes_u <- cb$hess; hes_c <- NULL }else hes_u <- hes_c <- NULL
  
  return(list(b = b, loglik = fyz_c$ll, unploglik = fyz_c$ll, iter = iter + cb$iter,
              hes = list(H = hes_u, Hp = hes_c), n = nrow(Y)))
}

glvmlss_penfit <- function(...){
  
  penON <- T
  lhist <- matrix(NA,nrow = control$autoL_iter+1, ncol = length(pen.control$lambda))
  ssehist <- numeric(control$autoL_iter+1)
  rownames(lhist) <- names(ssehist) <- paste0("c",c(1:(control$autoL_iter+1)))
  lhist[1,] <- pen.control$lambda
  
  while(penON){
  
  bold <- b
  convg <- ifelse(control$EM_iter == 0, T, F)
  pb <- lb2cb(b)
  cb <- pb[lb2cb(rb) == T]
  info <- control$mat.info
  iter <- 0
  llk <- numeric(control$EM_iter + 1)
  dfyz_t <- fyz(Y,ghQ,b,famL)
  pMat <- pM(b,rb, penalty = pen.control$penalty,
             lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
  llk[1] <- -dfyz_t$ll + 0.5*nrow(Y)*crossprod(cb,pMat$full)%*%cb
  EM_appHess <- control$EM_appHess
  
  while(!convg){
    iter <- iter + 1
    pb <- lb2cb(b)
    cb <- pb[lb2cb(rb) == T]
    pMat <- pM(b,rb, penalty = pen.control$penalty,
               lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
    d1ll_t <- -d1ll(Y,ghQ,b,famL,info,dfyz_t$pD,rb) + c(nrow(Y)*crossprod(cb,pMat$full))
    if(control$EM_use2d){
      if(!EM_appHess){
        da2ll_t <- -d2llEM(Y,ghQ,b,famL,info,dfyz_t$pD,rb) + nrow(Y)*pMat$full
      } else {
        da2ll_t <- -ad2ll(Y,ghQ,b,famL,info,dfyz_t$pD,rb) + nrow(Y)*pMat$full
      }
      cb <- cb - solve(da2ll_t, d1ll_t)
    } else {
      cb <- cb - control$EM_lrate*exp(-0.5*iter)*d1ll_t }
    if(cycle == 0) ssehist[1] <- SSE(loglambda = log(pen.control$lambda),b = b, gra = -d1ll_t, hes = da2ll_t, rb = rb, pml = pen.control, Y = Y)
    pb[lb2cb(rb) == T] <- cb
    b <- cb2lb(pb,b)
    dfyz_t <- fyz(Y,ghQ,b,famL)
    pMat <- pM(b,rb, penalty = pen.control$penalty,
               lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
    llk[iter+1] <- -dfyz_t$ll + 0.5*nrow(Y)*crossprod(cb,pMat$full)%*%cb
    eps1 <- llk[iter+1] - llk[iter]
    if(eps1 > 0) { iter <- iter - 1; break }
    convg <- ifelse(abs(eps1) < control$tol*1e4 || iter == control$EM_iter , T, F)
    if(control$verbose){ if(iter == 1) cat(paste0("\n EM iter: ",iter,", penalised marginal loglik.: ", round(-llk[iter+1],5))) else
      cat(paste0("\r EM iter: ",iter,", penalised marginal loglik.: ", round(-llk[iter+1],5))) }
  }
  
  if(control$solver == "trust"){
    if(control$verbose){ 
      if(!autoL) cat("\n Direct penalised MML estimation using trust-region algorithm ...") else
        if(cycle == 0) cat(paste0("\n Direct penalised MML estimation (trust-region, cycle ", cycle + 1,") ...")) else
          cat(paste0("\n Direct penalised MML estimation (trust-region, cycle ", cycle + 1,") ...")) }
    cb <- trust::trust(objfun = pll, parinit = cb, rinit = 1, rmax = 5,
                       fterm = control$tol, iterlim = control$iter.lim,
                       Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info, rb = rb,
                       pen.control = pen.control)
    pb[lb2cb(rb) == T] <- cb$argument
    b <- cb2lb(pb,b)
    if(control$verbose){ 
      if(autoL){ cat(paste0("\r Direct penalised MML estimation (trust-region, cycle ", cycle + 1,") ... Converged in ",
                            iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " ", control$solver, "),",
                            " pen. margllk: ", round(-cb$value,5)))
        } else {
          cat(paste0("\r Direct penalised MML estimation using trust-region algorithm ... Converged in ",
                     iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " ", control$solver, ")",
                     "\n (Approximate) Penalised Marginal loglikelihood: ", round(-cb$value,5))) } }
  } else {
    if(control$solver == "nlminb"){
      if(control$verbose){
        if(!autoL) cat("\n Direct penalised MML estimation using 'nlminb' function ...") else 
          if(cycle == 0) cat(paste0("\n Direct penalised MML estimation ('nlminb', cycle ", cycle + 1,") ...")) else
            cat(paste0("\n Direct penalised MML estimation ('nlminb', cycle ", cycle + 1,") ...")) }
      cb <- nlminb(start = cb, objective = plla, gradient = d1plla, control = list(iter.max = control$iter.lim),
                   Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info, rb = rb,
                   pen.control = pen.control)
      pb[lb2cb(rb) == T] <- cb$par
      b <- cb2lb(pb,b)
      if(control$verbose){
        if(!autoL){ cat(paste0("\r Direct penalised MML estimation 'nlminb' function ... Converged in ",
                               iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " nlminb)",
                               "\n (Approximate) Penalised Marginal loglikelihood: ", round(-cb$objective,5)))
        } else {
          cat(paste0("\r Direct penalised MML estimation ('nlminb', cycle ", cycle + 1,") ... Converged in ",
                     iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " nlminb),",
                     " pen. margllk: ", round(-cb$objective,5))) } }
    } else {
      if(control$verbose){ 
        if(!autoL) cat(paste0("\n Direct penalised MML estimation using ", control$solver," algorithm ...")) else
          if(cycle == 0) cat(paste0("\n Direct penalised MML estimation (", control$solver,", cycle ", cycle + 1,") ...")) else
            cat(paste0("\n Direct penalised MML estimation (", control$solver,", cycle ", cycle + 1,") ...")) }
      cb <- optim(cb, fn = plla, gr = d1plla, method = control$solver,
                  control = list(maxit = control$iter.lim, fnscale = 1),
                  Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info, rb = rb,
                  pen.control = pen.control)
      pb[lb2cb(rb) == T] <- cb$par
      b <- cb2lb(pb,b)
      if(control$verbose){
        if(!autoL){ cat(paste0("\r Direct penalised MML estimation using ", control$solver, " algorithm ... Converged.",
                               "\n (Approximate) Penalised Marginal loglikelihood: ", round(-cb$value,5)))
        } else {
          cat(paste0("\r Direct penalised MML estimation (", control$solver,", cycle ", cycle + 1,") ... Converged (",
                     "pen. margllk: ", round(-cb$value,5), ")")) } } } }
  
  for(r in names(b)){
    for(j in which(grepl("Z",colnames(b[[r]])))){
      c <- as.integer(substr(colnames(b[[r]])[j],nchar(colnames(b[[r]])[j]),nchar(colnames(b[[r]])[j])))
      if(b[[r]][c,j] < 0) b[[r]][,j] <- -b[[r]][,j]
    } }
  
  if(autoL){
    tObj <- op.lambda(Y,ghQ,b,famL,info,rb,pen.control)
    eps2 <- max(abs(pen.control$lambda - tObj$lambda))
    if(eps2 < control$tol) autoL <- F
    cycle <- cycle + 1
    lhist[cycle+1,]  <- tObj$lambda
    ssehist[cycle+1] <- tObj$sse
    if(ssehist[cycle+1] - ssehist[cycle] > 0){ 
      if(control$verbose) cat(paste0("\n \n SSE did not improve from cycle ", cycle, " to ", cycle+1,". Stopping automatic selection of tuning parameters in cycle ",cycle," (marked by *)."))
      autoL <- F
      b <- bold
      rownames(lhist)[cycle] <- names(ssehist)[cycle] <- paste0("c",cycle,"*")
    } else { pen.control$lambda <- tObj$lambda }
  }
  
  if(!autoL| cycle == control$autoL_iter) penON <- F
  
  }
  
  pMat <- pM(b,rb, penalty = pen.control$penalty,
             lambda = pen.control$lambda, w.alasso = pen.control$w.alasso, a = pen.control$a)
  b. <- lb2cb(b)[lb2cb(rb) == T]
  fyz_c <- fyz(Y,ghQ,b,famL)
  
  if(control$est.ci){
    hes_u <- -d2ll(Y,ghQ,b,famL,info, pd = fyz_c$pD, rb)
    hes_c <- hes_u + nrow(Y)*pMat$full
    if(!m2pdm(hes_c)$is.PD & control$verbose) cat(paste0("\n Unstable solution: Hessian matrix is not positive definite at solution."))
  } else if(control$solver == "trust"){ 
    hes_u <- cb$hess - nrow(Y)*pMat$full
    hes_c <- cb$hess } else hes_u <- hes_c <- NULL
  
  return(list(b = b, loglik = c(fyz_c$ll - 0.5*nrow(Y)*crossprod(b.,pMat$full)%*%b.),
              unploglik = fyz_c$ll,
              iter = iter + cb$iter, hes = list(H = hes_u, Hp = hes_c),
              # n = nrow(Y), sse = ssehist[cycle+1], lambda = lhist[cycle+1,]))
              n = nrow(Y), sse = ssehist[ssehist != 0], lambda = lhist[complete.cases(lhist),]))
}

glvmlss_sim <- function(n, family = list(),
                        mu.eq = ~ Z1 + Z2, sg.eq = NULL,
                        ta.eq = NULL, nu.eq = NULL,
                        control = list(), ...){
  
  environment(prep_fam) <- environment(prep_form) <- environment()
  famL <- prep_fam()
  form <- prep_form()
  p <- length(famL)
  environment(prep_Z) <- environment()
  q <- prep_Z()
  environment(sim_cont) <- environment(sim_Z) <- environment()
  control <- sim_cont()
  Z <- sim_Z()
  environment(sim_stva) <- environment()
  b <- sim_stva()$b
  
  Y <- sapply(1:p, function(i) famL[[i]]$sf(i,n,b,Z$Zmod))
  Y <- as.data.frame(Y); colnames(Y) <- paste0("Y", 1:p)
  
  return(list(Y = Y, Z = Z$Z, b = b))
}

