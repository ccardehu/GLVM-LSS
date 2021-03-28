

GLVM.fit <- function(Y, fam, form, loadmt, ghp = 10, iter.lim = 500,
                     tol = 1e-3, silent = F, icoefs = icoefs, useoptim = T, skipEM = T){
  
  # evaluate @ Y = simR$Y; fam = fam; form = form; silent = F; icoefs = lc; useoptim = T;
  # ghp = 10; loadmt = l1; tol = 1e-7; iter.lim = 700; skipEM = F
  
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  parY <- unique(unlist(lapply(1:length(fam),function(i) pFun(fam[i]))))
  for(i in parY){form[[i]] <- as.formula(form[[i]])}
  lvar <- unique(unlist(lapply(1:length(form), function(i) all.vars(form[[i]]))))
  lvar <- lvar[grep("Z", lvar, fixed = T)]
  if(length(lvar) == 0) stop("Latent variables in formula should be represented by letter Z (capital)")
  q. <- length(lvar)
  p. <- ncol(Y)
  pC <- vector(mode = "list", length = length(parY)); names(pC) <- parY
  for(i in parY){ for(j in 1:p.){ if(i %in% pFun(fam[j])) {pC[[i]] <- append(pC[[i]],j)} } }
  
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
   message("Restrictions matrix not supplied, simple structure assumed")
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
      message(paste0("Nonlinear functions for ", i, ", identification restriction impossed in item 1"))
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
  if(missing(tol)) tol <- sqrt(.Machine$double.eps)
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
   message("Initial Coefficients not supplied, Value = 1 assumed")
   bold <- lapply(parY, function(i) matrix(1, nrow = ncol(Y), ncol = ncol(gr$out[[i]])))
   names(bold) <- parY
  }
  for(i in parY){
   bold[[i]] <- bold[[i]]*loadmt[[i]]
   bold[[i]][-pC[[i]],] <- 0
  }
  
  hist <- matrix(NA,nrow = iter.lim, ncol = length(unlist(bold)))
  
  # Estimation (EM algorithm)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~
  
  tryCatch({
   while(eps1 > tol & iter < iter.lim){
    
    llo <- sum(log(mfy(Y,bold,gr,fam)))
    bnew <- bold
    
     if(skipEM == F){

      # bn <- NULL
      # HH <- NULL
      # bnn <- coefmod3(bnew)
      # bo1 <- decoefmod(bnew)
      # for(i in 1:ncol(Y)){
      # bo2 <- bnn[i,]
      # bNOT0 <- which(bo2 != 0)
      # tmp1 <- iScHesFun(i, bo1, Y, bnew, gr, fam)
      # Sm <- tmp1$score
      # Hm <- tmp1$hessian
      # bni <- c(as.matrix(bo2[bNOT0]) - solve(Hm)%*%as.matrix(Sm))
      # bn <- c(bn,bni)
      # }
      
      opttrst <- trust::trust(objfun = llkFun, parinit = decoefmod(bnew), rinit = 1, rmax = 10,
                              iterlim = 10, minimize = F, Y = Y, beta = bnew, ghQ = gr, fam = fam)
      
      bn <- opttrst$argument*decoefmod(loadmt)
      bnew <- coefmod2(bn,bold)
        
      # Computing epsilon  
      lln <- sum(log(mfy(Y,bnew,gr,fam)))
      # eps1 <- max((unlist(bnew) - unlist(bold))^2)
      eps1 <- abs(lln - llo) #/(0.1+abs(lln))
      iter = iter+1
      bold <- bnew
      hist[iter,] <- unlist(bnew)
      if(silent == F) cat("\n EM iteration: ", iter, ", LogLike. = ", round(lln,5), ", \U0394 LogLike. = ", lln-llo, sep = "")
      catmsg <- paste0("(after ", iter," EM iterations)")
     } else{iter = 100} #skipEM
   
    if(useoptim == T & iter >= 100){ #  
        if(skipEM == T) cat("\n Using `optim-BFGS` to find ML estimates") else cat(paste0("\n Using `optim-BFGS` to refine ML estimates ", catmsg))
        optmres <- optim(decoefmod(bnew), loglikFun, Y = Y, beta = bnew, ghQ = gr, fam = fam, method = "BFGS",
                         control = list(maxit = iter.lim, fnscale = -1))
        # opttrst <- trust::trust(objfun = llkFun, parinit = decoefmod(borg), rinit = 1, rmax = 10,
        #                         iterlim = 10, minimize = F, Y = Y, beta = bnew, ghQ = gr, fam = fam)
        bnew <- coefmod2(bet = optmres$par, beta = bnew)
        lln <- sum(log(mfy(Y,bnew,gr,fam)))
        iter <- optmres$convergence
        if(optmres$convergence != 0) warning("Maximization did not converge")
        eps1 <- tol/2
      }
    
   }

   Sco <- Hes <- NULL
   # if(skipEM == F){
   #  if("mu" %in% parY) {Sco$mu <- Sm; Hes$mu <- Hm}
   #  if("sigma" %in% parY) {Sco$sigma <- Ss; Hes$sigma <- Hs}
   # }
   hist <- hist[1:iter,]
    
   return(list(b = bnew, loglik = lln, loadmt = loadmt, iter = iter, gr = gr, Score = Sco, Hessian = Hes, cvgRes = hist,
               Y = as.data.frame(Y), fam = fam, formula = form, eps = eps1))
    
  }, error = function(e){
    
    print(paste("\n Error in Estimation, proceeded with next simulation \n Error:",e))
    berr <- lapply(parY, function(i) matrix(-999, nrow = ncol(Y), ncol = ncol(gr$out[[i]])))
    names(berr) <- parY; for(i in parY){ berr[[i]] <- berr[[i]]*loadmt[[i]] ; loadmt[[i]] <- loadmt[[i]]*-999}
    
    return(list(b = berr, loglik = -999, loadmt = loadmt,iter = -999))
    
  })
  
}
