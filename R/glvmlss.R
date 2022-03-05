glvmlss_fit <- function(...){
  
  convg <- F
  info <- control$mat.info
  iter <- 0
  llk <- numeric(control$EM_iter + 1)
  dfyz_t <- fyz(Y,ghQ,b,famL)
  llk[1] <- dfyz_t$ll
  use2d <- control$EM_use2d
  
  while(!convg){
    iter <- iter + 1
    pb <- lb2cb(b)
    cb <- pb[lb2cb(b) != 0]
    d1ll_t <- d1ll(Y,ghQ,b,famL,info,dfyz_t$pD)
    if(use2d){
      da2ll_t <- ad2ll(Y,ghQ,b,famL,info,dfyz_t$pD)
      cb <- cb + solve(da2ll_t, d1ll_t)
    } else {
      cb <- cb + 0.001*exp(-0.5*iter)*d1ll_t }
    pb[lb2cb(b) != 0] <- cb
    b <- cb2lb(pb,b)
    llk[iter+1] <- fyz(Y,ghQ,b,famL)$ll
    if(llk[iter+1] < llk[iter]) use2d <- F
    convg <- ifelse(abs(llk[iter] - llk[iter+1]) < sqrt(.Machine$double.eps) || iter == control$EM_iter ,
                    T, F)
    if(control$verbose){ if(iter == 1) cat(paste0("\n EM iter: ",iter,", loglk: ", round(llk[iter+1],5))) else
      cat(paste0("\r EM iter: ",iter,", loglk: ", round(llk[iter+1],5))) }
  }
  
  if(control$solver == "trust"){
    if(control$verbose) cat("\n Refining estimates using trust-region algorithm ...")
    cb <- trust::trust(objfun = ll, parinit = cb, rinit = 1, rmax = 5,
                       fterm = sqrt(.Machine$double.eps), iterlim = 1e3,
                       minimize = F, Y = Y, bg = b, ghQ = ghQ,
                       famL = famL, info = info)
    pb[lb2cb(b) != 0] <- cb$argument
    b <- cb2lb(pb,b); iter <- iter + cb$iter 
    if(control$verbose) cat(paste0("\r Refining estimates using trust-region algorithm ... Converged after ",
                        iter, " iterations (30 EM + ", cb$iter, " ", control$solver, ")",
                        "\n (Approximate) Marginal loglikelihood: ", round(cb$value,5))) } else {
    if(control$verbose) cat(paste0("\n Refining estimates using ", control$solver," algorithm ..."))
    cb <- optim(cb, fn = lla, gr = d1lla, method = control$solver,
                control = list(maxit = 1e3, fnscale = -1),
                Y = Y, bg = b, ghQ = ghQ, famL = famL, info = info)
                pb[lb2cb(b) != 0] <- cb$par
                b <- cb2lb(pb,b); iter <- iter + sum(cb$counts, na.rm = T)
    if(control$verbose) cat(paste0("\r Refining estimates using ", control$solver, " algorithm ... Converged after ",
                        iter, " iterations (30 EM + ", sum(cb$counts, na.rm = T)," ", control$solver, ")",
                        "\n (Approximate) Marginal loglikelihood: ", round(cb$value,5))) }
  
  for(r in names(b)){
    for(j in 1:q){
      if(b[[r]][j, paste0("Z",j)] < 0) b[[r]][, paste0("Z",j)] <- -b[[r]][, paste0("Z",j)]
    } }
  
  return(list(b = b, loglik = cb$value, iter = iter))
}

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

