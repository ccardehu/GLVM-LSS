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
  
  environment(glvmlss_fit) <- environment()
  fit <- glvmlss_fit(...)
  
  structure(fit, class = "glvmlss")
}

glvmlss_fit <- function(...){
  
  convg <- ifelse(control$EM_iter == 0, T, F)
  pb <- lb2cb(b)
  cb <- pb[lb2cb(b) != 0]
  info <- control$mat.info
  iter <- 0
  llk <- numeric(control$EM_iter + 1)
  dfyz_t <- fyz(Y,ghQ,b,famL)
  llk[1] <- -dfyz_t$ll
  EM_appHess <- control$EM_appHess
  
  while(!convg){
    iter <- iter + 1
    pb <- lb2cb(b)
    cb <- pb[lb2cb(b) != 0]
    d1ll_t <- d1ll(Y,ghQ,b,famL,info,dfyz_t$pD)
    if(control$EM_use2d){
      if(!EM_appHess){
        da2ll_t <- -d2llEM(Y,ghQ,b,famL,info,dfyz_t$pD)
        # if(!m2pdm(da2ll_t)$is.P) da2ll_t <- da2ll_t - control$EM_lrate*diag(nrow(da2ll_t))
      } else {
        da2ll_t <- ad2ll(Y,ghQ,b,famL,info,dfyz_t$pD)
        # if(!m2pdm(da2ll_t)$is.P) da2ll_t <- da2ll_t + control$EM_lrate*diag(nrow(da2ll_t))
        }
      cb <- cb + solve(da2ll_t, d1ll_t)
    } else {
      cb <- cb + control$EM_lrate*exp(-0.5*iter)*d1ll_t }
    pb[lb2cb(b) != 0] <- cb
    b <- cb2lb(pb,b)
    llk[iter+1] <- -fyz(Y,ghQ,b,famL)$ll
    eps <- llk[iter+1] - llk[iter]
    # if(eps < 0 && control$EM_use2d) EM_appHess <- T
    convg <- ifelse(abs(eps) < control$tol*1e4 || iter == control$EM_iter , T, F)
    if(control$verbose){ if(iter == 1) cat(paste0("\n EM iter: ",iter,", marginal loglik.: ", round(-llk[iter+1],5))) else
      cat(paste0("\r EM iter: ",iter,", marginal loglik.: ", round(-llk[iter+1],5))) }
    # if(iter%%control$EM_upeach) dfyz_t <- fyz(Y,ghQ,b,famL)
    dfyz_t <- fyz(Y,ghQ,b,famL)
  }
  
  if(control$solver == "trust"){
    if(control$verbose) cat("\n Refining estimates using trust-region algorithm ...")
    cb <- trust::trust(objfun = ll, parinit = cb, rinit = 1, rmax = 5,
                       fterm = control$tol, iterlim = control$iter.lim,
                       Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info)
    pb[lb2cb(b) != 0] <- cb$argument
    b <- cb2lb(pb,b)
    if(control$verbose) cat(paste0("\r Refining estimates using trust-region algorithm ... Converged after ",
                        iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " ", control$solver, ")",
                        "\n (Approximate) Marginal loglikelihood: ", round(-cb$value,5))) } else {
    if(control$solver == "nlminb"){
    if(control$verbose) cat("\n Refining estimates using 'nlminb' function ...")
    cb <- nlminb(start = cb, objective = lla, gradient = d1lla, control = list(iter.max = control$iter.lim),
                 Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info)
    pb[lb2cb(b) != 0] <- cb$par
    b <- cb2lb(pb,b)
    if(control$verbose) cat(paste0("\r Refining estimates using 'nlminb' function ... Converged after ",
                        iter + cb$iter, " iterations (", iter," EM + ", cb$iter, " nlminb)",
                        "\n (Approximate) Marginal loglikelihood: ", round(-cb$objective,5))) } else {
    if(control$verbose) cat(paste0("\n Refining estimates using ", control$solver," algorithm ..."))
    cb <- optim(cb, fn = lla, gr = d1lla, method = control$solver,
                control = list(maxit = control$iter.lim, fnscale = 1),
                Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info)
    pb[lb2cb(b) != 0] <- cb$par
    b <- cb2lb(pb,b)
    if(control$verbose) cat(paste0("\r Refining estimates using ", control$solver, " algorithm ... Converged.",
                        "\n (Approximate) Marginal loglikelihood: ", round(-cb$value,5))) } }
  
  for(r in names(b)){
    for(j in which(grepl("Z",colnames(b[[r]])))){
      c <- as.integer(substr(colnames(b[[r]])[j],nchar(colnames(b[[r]])[j]),nchar(colnames(b[[r]])[j])))
      if(b[[r]][c,j] < 0) b[[r]][,j] <- -b[[r]][,j]
    } }
  
  fyz_c <- fyz(Y,ghQ,b,famL)
  hes_c <- -d2ll(Y,ghQ,b,famL,info, pd = fyz_c$pD)
  if(!m2pdm(hes_c)$is.PD & control$verbose) cat(paste0("\n Unstable solution: Hessian matrix is not positive definite at solution."))
  
  return(list(b = b, loglik = fyz_c$ll, iter = iter + cb$iter, hes = hes_c))
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
  b <- sim_stva()
  
  Y <- sapply(1:p, function(i) famL[[i]]$sf(i,n,b,Z$Zmod))
  Y <- as.data.frame(Y); colnames(Y) <- paste0("Y", 1:p)
  
  return(list(Y = Y, Z = Z$Z, b = b))
}

