

GLVM.fit <- function(Y, fam, form, loadmt, ghp = 10, iter.lim = 500,
                     tol = 1e-3, silent = F, icoefs = icoefs, useoptim = T, skipEM = T){
  
  # evaluate @ Y = simR$Y; fam = fam; form = form; silent = F; icoefs = lc; useoptim = T; ghp = 20; loadmt = l1; tol = 1e-7; iter.lim = 700; skipEM = T
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  parY <- unlist(unique(lapply(1:length(fam),function(i) pFun(i,fam[i]))))
  for(i in parY){form[[i]] <- as.formula(form[[i]])}
  lvar <- unique(unlist(lapply(1:length(form), function(i) all.vars(form[[i]]))))
  lvar <- lvar[grep("Z", lvar, fixed = T)]
  if(length(lvar) == 0) stop("Latent variables in formula should be represented by letter Z (capital)")
  q. <- length(lvar)
  p. <- ncol(Y)
  
  # Latent variables (Z. for each mu, sigma, etc.)
  # ______________________________________________
  if(missing(ghp)) ghp <- 6
  gr <- mvgH(n = ghp, formula = form)
  
  # Restriction matrices (for each mu, sigma, etc.)
  # ______________________________________________
  
  if(!missing(loadmt)){
    if(!is.list(loadmt)) stop("Provided `loadmt` should be a list with dim p*(q+intercept) restriction matrices for mu, sigma, tau, nu")
    if(length(loadmt) != length(parY)) stop("Number of matrices in `loadmt` should match number of parameters mu, sigma, tau, nu")
    names(loadmt) <- parY 
    for(i in parY) {
      colnames(loadmt[[i]]) <- colnames(gr$out[[i]]); rownames(loadmt[[i]]) <- colnames(Y)
      if(length(loadmt[[i]]) != p.*ncol(gr$out[[i]])) stop("Provided `loadmt` is not lenght p*(q+intercept), revise dimension")
      if(!is.matrix(loadmt[[i]])) loadmt[[i]] <- matrix(loadmt[[i]], nrow = p., ncol = ncol(gr$out[[i]]))
    }
  } else {
    warning("Restrictions matrix not supplied, simple structure assumed", call. = F)
    loadmt <- NULL
    for(i in parY){
      if(length(all.vars(form[[i]])) > 1){
        r. <- p.%/%(ncol(gr$out[[i]])-1)
        tmp1 <- matrix(0, nrow = p., ncol = ncol(gr$out[[i]]), dimnames = list(colnames(Y), colnames(gr$out[[i]])))
        for(j in seq_along(all.vars(form[[i]]))){
          tmp2 <- seq((j-1)*(r.)+1,(j)*(r.),length.out = r.)
          tmp1[tmp2,all.vars(form[[i]])[j]] <- 1
          if(j == length(all.vars(form[[i]])) && p.%%2 != 0) tmp1[(tail(tmp2,1)+1):nrow(tmp1),all.vars(form[[i]])[j]] <- 1
        }
        tmp1[,!(colnames(tmp1) %in% all.vars(form[[i]]))] <- 1
        loadmt[[i]] <- tmp1; rm(tmp1,tmp2)
      } else{
        loadmt[[i]] <- matrix(1,nrow = p., ncol = ncol(gr$out[[i]]))
        # loadmt[[i]][1,ncol(gr$out[[i]])] <- 0 # ! HERE
        if(ncol(gr$out[[i]])-1 > 1){
          warning(paste0("Nonlinear functions for ", i, ", identification restriction impossed in item 1"), call. = F)
          loadmt[[i]][1,ncol(gr$out[[i]])] <- 0
        }
      }
      colnames(loadmt[[i]]) <- colnames(gr$out[[i]])
      rownames(loadmt[[i]]) <- colnames(Y)  
    }
  }
  
  # Starting values
  # _______________
  iter <- 0
  if(missing(tol)) tol <- 1e-3
  eps1 <- tol + 1
  
  if(!missing(icoefs)){
    if(!is.list(icoefs)) stop("Provided initial `icoefs` should be a list with dim p*(q+intercept) loading matrices for mu, sigma, tau, nu")
    if(length(icoefs) != length(parY)) stop("Number of matrices in `icoefs` should match number of parameters mu, sigma, tau, nu")
    names(icoefs) <- parY 
    for(i in parY) {
      colnames(icoefs[[i]]) <- colnames(gr$out[[i]]); rownames(icoefs[[i]]) <- colnames(Y)
      if(length(icoefs[[i]]) != p.*ncol(gr$out[[i]])) stop("Provided `icoefs` is not lenght p*(q+intercept), revise dimension")
      if(!is.matrix(icoefs[[i]])) icoefs[[i]] <- matrix(icoefs[[i]], nrow = p., ncol = ncol(gr$out[[i]]))
      icoefs[[i]] <- icoefs[[i]] + runif(length(icoefs[[i]]), min = -tol, max = tol)
    }
    bold <- icoefs
  } else {
    warning("Initial Coefficients not supplied, Value = 1 assumed", call. = F)
    bold <- lapply(parY, function(i) matrix(1, nrow = ncol(Y), ncol = ncol(gr$out[[i]])))
    names(bold) <- parY
  }
  
  for(i in parY){
    bold[[i]] <- bold[[i]]*loadmt[[i]]
  }
  
  
  #if(missing(useoptim)) useoptim <- "BFGS"
  
  hist <- matrix(NA,nrow = iter.lim, ncol = length(unlist(bold)))
  
  # Estimation (EM algorithm)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~
  
  tryCatch({
    
    while(eps1 > tol & iter < iter.lim){
      
      llo <- sum(log(mfy(Y,bold,gr,fam)))
      bnew <- bold
      
      if(skipEM == F){
        
        if("mu" %in% parY){
          efy <- mfy(Y,bnew,gr,fam)
          EC <- sapply(1:nrow(gr$points), function(z) mfyz(z,Y,gr$out,bnew,fam))/efy
          Sm <- lapply(1:p, function(r) bmsc(r,Y,bnew,gr,fam,EC))
          Hm <- lapply(1:p, function(r) bmhe(r,Y,bnew,gr,fam,EC))
          bnew$mu <- t(sapply(1:p, function(i) as.matrix(bnew$mu[i,]) - solve(Hm[[i]])%*%Sm[[i]]))
          bnew$mu <- bnew$mu*loadmt$mu
        }
        
        if("sigma" %in% parY){
          efy <- mfy(Y,bnew,gr,fam)
          EC <- sapply(1:nrow(gr$points), function(z) mfyz(z,Y,gr$out,bnew,fam))/efy
          Ss <- lapply(1:p, function(r) bssc(r,Y,bnew,gr,fam,EC))
          Hs <- lapply(1:p, function(r) bshe(r,Y,bnew,gr,fam,EC))
          bnew$sigma <- t(sapply(1:p, function(i) as.matrix(bnew$sigma[i,]) - solve(Hs[[i]])%*%Ss[[i]]))
          if(nrow(bnew$sigma) == 1) bnew$sigma <- t(bnew$sigma)
          bnew$sigma <- bnew$sigma*loadmt$sigma
        }
        
        # Computing epsilon  
        
        lln <- sum(log(mfy(Y,bnew,gr,fam)))
        eps1 <- abs(lln - llo)
        
        iter = iter+1
        bold <- bnew
        hist[iter,] <- unlist(bnew)
        if(silent == F) cat("\r EM iteration = ", iter, ", LogLike. = ", round(lln,5), ", \U0394 LogLike. = ", lln-llo, sep = "") # "\n",
        catmsg <- paste0("(after ", iter," EM iterations)")
      } else{iter = 40} #skipEM
      
      if(useoptim == T & iter >= 40){
        if(skipEM == T) cat("\n Using `optim-BFGS` to find ML estimates") else cat(paste0("\n Using `optim-BFGS` to refine ML estimates", catmsg))
        optmres <- optim(unlist(bnew), loglik, beta = bnew, ghQ = gr, loadmt = loadmt, method = "BFGS", control = list(maxit = iter.lim, fnscale = -1))
        bnew <- coefmod(bet = optmres$par, beta = bnew, gr = gr,loadmt)
        lln <- sum(log(mfy(Y,bnew,gr,fam)))
        iter <- optmres$convergence
        if(optmres$convergence != 0) warning("Maximization did not converge properly")
        eps1 <- tol/2
      }
      # 
      # if(useoptim == "trustR" & iter >= 40){
      #   if(skipEM == T) cat("\n Using `trust` to find ML estimates") else cat(paste0("\n Using `trust` to refine ML estimates", catmsg))
      #   optmres <- trust::trust(loglik, beta = bnew, ghQ = gr, loadmt = loadmt, method = "BFGS", control = list(maxit = iter.lim, fnscale = -1))
      #   bnew <- coefmod(bet = optmres$par, beta = bnew, gr = gr,loadmt)
      #   lln <- sum(log(mfy(Y,bnew,gr,fam)))
      #   iter <- optmres$convergence
      #   if(optmres$convergence != 0) warning("Maximization did not converge properly")
      #   eps1 <- tol/2
      # }
      
    }
    
    Sco <- Hes <- NULL
    if(skipEM == F){
      if("mu" %in% parY) {Sco$mu <- Sm; Hes$mu <- Hm}
      if("sigma" %in% parY) {Sco$sigma <- Ss; Hes$sigma <- Hs}
    }
    
    hist <- hist[1:iter,]
    
    return(list(b = bnew, loglik = lln, loadmt = loadmt, iter = iter, gr = gr, Score = Sco, Hessian = Hes, cvgRes = hist))
    
  }, error = function(e){
    
    print(paste("\n Error in Estimation, proceeded with next simulation \n Error:",e))
    berr <- lapply(parY, function(i) matrix(-999, nrow = ncol(Y), ncol = ncol(gr$out[[i]])))
    names(berr) <- parY; for(i in parY){ berr[[i]] <- berr[[i]]*loadmt[[i]] ; loadmt[[i]] <- loadmt[[i]]*-999}
    
    return(list(b = berr, loglik = -999, loadmt = loadmt,iter = -999))
    
  })
  
}
