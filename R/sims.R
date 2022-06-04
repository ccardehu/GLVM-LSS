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
                A <- glvmlss_sim(n, famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
                c0 <- Sys.time()
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = A$b,
                             nQP = 25, solver = "trust")
                c1 <- Sys.time()
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ex2 <- c(c1-c0, B$iter, length(lb2cb(B$b)))
                return(c(cof,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("E1n",n, "p", p,"_",Sys.Date(),".Rds"))
  return(FCOL)
}

glvmlss_parsimE2 <- function(nsim, saveRes = T){
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
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = A$b,
                             nQP = 25, solver = "trust")
                c1 <- Sys.time()
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ex2 <- c(c1-c0, B$iter, length(lb2cb(B$b)))
                return(c(cof,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("E2n",n, "p", p,"_",Sys.Date(),".Rds"))
  return(FCOL)
}

glvmlss_parsimE3 <- function(nsim, saveRes = T){
  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl,ls(parent.frame()))
  doSNOW::registerDoSNOW(cl)
  progress <- function(a) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), a)
  opts <- list(progress = progress)
  FCOL <- suppressWarnings(
    foreach(l = 1:nsim,
            .combine = rbind,
            .options.snow = opts) %dopar% { # .export = ls(parent.frame()),
              tryCatch({
                A <- glvmlss_sim(n,famt, mu.eq = ~ Z1, sg.eq = ~ Z1, start.val = lc)
                c0 <- Sys.time()
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, start.val = A$b,
                             nQP = 40, solver = "trust")
                c1 <- Sys.time()
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ex2 <- c(c1-c0, B$iter, length(lb2cb(B$b)))
                return(c(cof,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("E3n",n, "p", p,"_",Sys.Date(),".Rds"))
  return(FCOL)
}

glvmlss_parsimpost <- function(file, out = c("MSE","AB"), trim = 0,
                               mu.eq = ~ Z1+Z2, sg.eq = NULL,
                               ta.eq = NULL, nu.eq = NULL, plot = F, outs = T){
  
  colMeans_ <- function(X = X, trim = trim){
    return(c(apply(X,MARGIN = 2,mean, trim = trim)))
  }
  
  X <- readRDS(file)
  
  environment(prep_form) <- environment()
  form <- prep_form()
  environment(prep_Z) <- environment()
  q <- prep_Z()
  
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
  if(q != 1) li1 <- c("mu1.2","mu2.1","sigma1.2","sigma2.1","tau1.2","tau2.1","nu1.2","nu2.1") else li1 <- c("NA")
  li2 <- lapply(li1, function(i) !grepl(i,nM))
  for(i in 1:length(li2)){ ip <- ip & li2[[i]]}; names(ip) <- nM; rm(list=c("li1","li2"))
  ip[grepl(".0",nM,fixed = T)] <- T; nM <- nM[ip];
  
  for(i in 1:q){
    tmp <- grepl(paste0(".",i),nM,fixed = T)
    names(tmp) <- nM
    tmp <- tmp[nM]
    assign(paste0("ip",i), tmp)
  }
  
  Xor <- X
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
  ita0 <- grepl("tau",colnames(Xb0),fixed = T) & grepl(".0",colnames(Xb0),fixed = T)
  ital <- grepl("tau",colnames(Xb0),fixed = T) & !grepl(".0",colnames(Xb0),fixed = T)
  inu0 <- grepl("nu",colnames(Xb0),fixed = T) & grepl(".0",colnames(Xb0),fixed = T)
  inul <- grepl("nu",colnames(Xb0),fixed = T) & !grepl(".0",colnames(Xb0),fixed = T)
  ili <- list("mu0" = imu0, "mu1" = imul, "sg0" = isg0, "sg1" = isgl,
              "ta0" = ita0, "ta1" = ital, "nu0" = inu0, "nu1" = inul)
  for(i in names(ili)){ if(all(!ili[[i]])) ili[[i]] <- NULL }
  
  # Fixing sign
  # ~~~~~~~~~~~
  
  for(i in 1:nrow(Xbe)){
    for(j in 1:q){
      for(ii in names(ili)[grepl("1",names(ili),fixed = T)]){
        ij <- which.max(abs(Xbe[i,ili[[ii]] & get(paste0("ip",j))]))
        if(sign(Xbe[i,ili[[ii]] & get(paste0("ip",j))][ij]) != sign(Xb0[i,ili[[ii]] & get(paste0("ip",j))][ij])) Xbe[i,ili[[ii]] & get(paste0("ip",j))] <- -Xbe[i,ili[[ii]] & get(paste0("ip",j))]
      }
      # ij <- which.max(abs(Xbe[i,get(paste0("ip",j))]))
      # if(sign(Xbe[i,get(paste0("ip",j))][ij]) != sign(Xb0[i,get(paste0("ip",j))][ij])) Xbe[i,get(paste0("ip",j))] <- -Xbe[i,get(paste0("ip",j))]
    }
    # if(!all(sign(Xbe[i,ipp]) == sign(Xb0[i,ipp]))) Xbe[i,ipp] <- -Xbe[i,ipp]
    # if(sign(max(Xbe[i,ipp])) != sign(max(Xb0[i,ipp])) && sign(min(Xbe[i,ipp])) != sign(min(Xb0[i,ipp]))) Xbe[i,ipp] <- -Xbe[i,ipp]
  }
  
  res <- nam <- NULL
  for(ii in names(ili)){
    if("MSE" %in% out) assign(paste0("AvMSE",ii), mean(colMeans_((Xb0[,ili[[ii]]] - Xbe[,ili[[ii]]])^2, trim = trim)) )
    if("AB" %in% out) assign(paste0("AvAB",ii), mean(abs(colMeans_(Xb0[,ili[[ii]]] - Xbe[,ili[[ii]]], trim = trim))) )
    if("SB" %in% out) assign(paste0("AvSB",ii), mean((colMeans_(Xb0[,ili[[ii]]] - Xbe[,ili[[ii]]], trim = trim))^2) )
    if("RB" %in% out) assign(paste0("AvRB",ii), mean((colMeans_(Xb0[,ili[[ii]]] - Xbe[,ili[[ii]]], trim = trim))/abs(colMeans_(Xb0[,ili[[ii]]], trim = trim))) ) } 
  
  if("MSE" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvMSE",ii))); nam <- c(nam, paste0("AvMSE",ii)) }
  if("AB" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvAB",ii))); nam <- c(nam, paste0("AvAB",ii)) }
  if("SB" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvSB",ii))); nam <- c(nam, paste0("AvSB",ii)) }
  if("RB" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvRB",ii))); nam <- c(nam, paste0("AvRB",ii)) }
  names(res) <- nam
  
  if("iter" %in% out){ nam <- c(names(res),"AvIter"); res <- c(res,mean(Xit)); names(res) <- nam }
  if("time" %in% out){ nam <- c(names(res),"AvTime"); res <- c(res,mean(Xti)); names(res) <- nam }
  if("PVS" %in% out){ nam <- c("PVS",names(res)); res <- c(nrow(X)/nrow(Xor),res); names(res) <- nam }
  
  if(plot){
    nid <- NULL
    for(ii in names(ili)){
      nid <- c(nid,which(ili[[ii]]))
    }
    boxplot(Xb0[,nid] - Xbe[,nid],
            main = paste0("Bias (by type of parameter): Simulation ", EX, ", n = ",n,", p = ",p),
            xlab = "", ylab = "Bias", pch = 16,
            col = "gray90", xaxt = "n", cex = 0.5, outline = outs)
    mtext(side = 1, text = "Parameter index", line = 4)
    for(i in seq_along(nM)){ axis(1, at=i, labels = nM[nid][i] , las = 2, cex.axis = 0.7) }
    abline(h = 0, col = 2, lwd = 3, lty = 3)
    a = 0
    abline(v = a, lty = 2, col = "purple")
    for(ii in names(ili)){
      a = a + sum(ili[[ii]] == T)
      abline(v = a + 0.5, lty = 2, col = "purple")
    }
  }

  return(c(format(round(res,4),nsmall = 4,scientific = F)))
}

glvmlss_parsimE1_pml <- function(nsim, saveRes = T){
  cores <- ifelse(parallel::detectCores() > 32, 32, parallel::detectCores())
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
                
                A <- glvmlss_sim(n, famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
                t0 <- Sys.time()
                
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, solver = "trust")
                C1a <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2,
                               solver = "trust", penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 1)
                C1b <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2,
                               solver = "trust", penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 2)
                C1c <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2,
                               solver = "trust", penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 3)
                C1d <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2,
                               solver = "trust", penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 4)
                C1e <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2,
                               solver = "trust", penalty = "alasso", lambda = "auto", w.alasso = B$b)
                
                t1 <- Sys.time()
                cof <- c(lb2cb(A$b), lb2cb(B$b), lb2cb(C1a$b), lb2cb(C1b$b), lb2cb(C1c$b), lb2cb(C1d$b), lb2cb(C1e$b))
                ex2 <- c(B$GBIC, C1a$GBIC, C1b$GBIC, C1c$GBIC, C1d$GBIC, C1e$GBIC,
                         t1-t0, length(lb2cb(B$b)))
                
                return(c(cof,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("E1n",n, "p", p,"_",Sys.Date(),"_C2.Rds"))
  return(FCOL)
}

glvmlss_parsimE2_pml <- function(nsim, saveRes = T){
  cores <- ifelse(parallel::detectCores() > 32, 32, parallel::detectCores())
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
                
                A <- glvmlss_sim(n, famt, mu.eq = ~ Z1+Z2, start.val = lc)
                t0 <- Sys.time()
                
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, solver = "trust")
                C1a <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2,
                               solver = "trust", penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 1)
                C1b <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2,
                               solver = "trust", penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 2)
                C1c <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2,
                               solver = "trust", penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 3)
                C1d <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2,
                               solver = "trust", penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 4)
                C1e <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2,
                               solver = "trust", penalty = "alasso", lambda = "auto", w.alasso = B$b)
                
                t1 <- Sys.time()
                cof <- c(lb2cb(A$b), lb2cb(B$b), lb2cb(C1a$b), lb2cb(C1b$b), lb2cb(C1c$b), lb2cb(C1d$b), lb2cb(C1e$b))
                ex2 <- c(B$GBIC, C1a$GBIC, C1b$GBIC, C1c$GBIC, C1d$GBIC, C1e$GBIC,
                         t1-t0, length(lb2cb(B$b)))
                
                return(c(cof,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("E2n",n, "p", p,"_",Sys.Date(),"_C2.Rds"))
  return(FCOL)
}

glvmlss_parsimpost_pml <- function(file,
                                   out = c("MSE","AB"), trim = 0,
                               mu.eq = ~ Z1+Z2, sg.eq = NULL,
                               ta.eq = NULL, nu.eq = NULL){
  
  colMeans_ <- function(X = X, trim = trim){
    return(c(apply(X,MARGIN = 2,mean, trim = trim)))
  }
  
  X <- readRDS(file)
  
  environment(prep_form) <- environment()
  form <- prep_form()
  environment(prep_Z) <- environment()
  q <- prep_Z()
  
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
  if(q != 1) li1 <- c("mu1.2","mu2.1","sigma1.2","sigma2.1","tau1.2","tau2.1","nu1.2","nu2.1") else li1 <- c("NA")
  li2 <- lapply(li1, function(i) !grepl(i,nM))
  for(i in 1:length(li2)){ ip <- ip & li2[[i]]}; names(ip) <- nM; rm(list=c("li1","li2"))
  ip[grepl(".0",nM,fixed = T)] <- T; nM <- nM[ip];
  
  for(i in 1:q){
    tmp <- grepl(paste0(".",i),nM,fixed = T)
    names(tmp) <- nM
    tmp <- tmp[nM]
    assign(paste0("ip",i), tmp)
  }
  
  Xor <- X
  X  <- X[complete.cases(X),]; #X <- X[apply(X,1,min)!=-999, ]
  ix <- mean(X[,ncol(X)])
  X <- X[,-ncol(X)]
  rR <- !(X[,ncol(X)] == 1000)
  Xb0 <- X[rR,seq_len(ix[1])[ip]]; colnames(Xb0) <- nM ;
  Xbeg0 <- X[rR,((ix[1]+1):(2*ix[1]))[ip]]; colnames(Xbeg0) <- nM
  Xbeg1 <- X[rR,((2*ix[1]+1):(3*ix[1]))[ip]]; colnames(Xbeg1) <- nM
  Xbeg2 <- X[rR,((3*ix[1]+1):(4*ix[1]))[ip]]; colnames(Xbeg2) <- nM
  Xbeg3 <- X[rR,((4*ix[1]+1):(5*ix[1]))[ip]]; colnames(Xbeg3) <- nM
  Xbeg4 <- X[rR,((5*ix[1]+1):(6*ix[1]))[ip]]; colnames(Xbeg4) <- nM
  Xbegf <- X[rR,((6*ix[1]+1):(7*ix[1]))[ip]]; colnames(Xbegf) <- nM
  XGBIC <- X[, seq(from = (7*ix[1])+1, length.out = 6)]
  if("time" %in% out) Xti <- X[, ncol(X)]
  
  # Fixing sign
  
  for(i in 1:nrow(Xbeg0)){
    for(j in 1:q){
      ij <- which.max(abs(Xbeg0[i,get(paste0("ip",j))]))
      if(sign(Xbeg0[i,get(paste0("ip",j))][ij]) != sign(Xb0[i,get(paste0("ip",j))][ij])) Xbeg0[i,get(paste0("ip",j))] <- -Xbeg0[i,get(paste0("ip",j))]
      if(sign(Xbeg1[i,get(paste0("ip",j))][ij]) != sign(Xb0[i,get(paste0("ip",j))][ij])) Xbeg1[i,get(paste0("ip",j))] <- -Xbeg1[i,get(paste0("ip",j))]
      if(sign(Xbeg2[i,get(paste0("ip",j))][ij]) != sign(Xb0[i,get(paste0("ip",j))][ij])) Xbeg2[i,get(paste0("ip",j))] <- -Xbeg2[i,get(paste0("ip",j))]
      if(sign(Xbeg3[i,get(paste0("ip",j))][ij]) != sign(Xb0[i,get(paste0("ip",j))][ij])) Xbeg3[i,get(paste0("ip",j))] <- -Xbeg3[i,get(paste0("ip",j))]
      if(sign(Xbeg4[i,get(paste0("ip",j))][ij]) != sign(Xb0[i,get(paste0("ip",j))][ij])) Xbeg4[i,get(paste0("ip",j))] <- -Xbeg4[i,get(paste0("ip",j))]
      if(sign(Xbegf[i,get(paste0("ip",j))][ij]) != sign(Xb0[i,get(paste0("ip",j))][ij])) Xbegf[i,get(paste0("ip",j))] <- -Xbegf[i,get(paste0("ip",j))]
    }
  }
  
  imu0 <- grepl("mu",colnames(Xb0),fixed = T) & grepl(".0",colnames(Xb0),fixed = T)
  imul <- grepl("mu",colnames(Xb0),fixed = T) & !grepl(".0",colnames(Xb0),fixed = T)
  isg0 <- grepl("sigma",colnames(Xb0),fixed = T) & grepl(".0",colnames(Xb0),fixed = T)
  isgl <- grepl("sigma",colnames(Xb0),fixed = T) & !grepl(".0",colnames(Xb0),fixed = T)
  ili <- list("int" = imu0 | isg0, "loa" = imul | isgl)
  for(i in names(ili)){ if(all(!ili[[i]])) ili[[i]] <- NULL }
  
  # Check from here
  
  res <- nam <- NULL
  for(ii in names(ili)){
    if("MSE" %in% out){
      assign(paste0("AvMSEg0",ii), mean(colMeans_((Xb0[,ili[[ii]]] - Xbeg0[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEg1",ii), mean(colMeans_((Xb0[,ili[[ii]]] - Xbeg1[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEg2",ii), mean(colMeans_((Xb0[,ili[[ii]]] - Xbeg2[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEg3",ii), mean(colMeans_((Xb0[,ili[[ii]]] - Xbeg3[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEg4",ii), mean(colMeans_((Xb0[,ili[[ii]]] - Xbeg4[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEgf",ii), mean(colMeans_((Xb0[,ili[[ii]]] - Xbegf[,ili[[ii]]])^2, trim = trim)))
    }
    
    if("AB" %in% out){
      assign(paste0("AvABg0",ii), mean(abs(colMeans_(Xb0[,ili[[ii]]] - Xbeg0[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABg1",ii), mean(abs(colMeans_(Xb0[,ili[[ii]]] - Xbeg1[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABg2",ii), mean(abs(colMeans_(Xb0[,ili[[ii]]] - Xbeg2[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABg3",ii), mean(abs(colMeans_(Xb0[,ili[[ii]]] - Xbeg3[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABg4",ii), mean(abs(colMeans_(Xb0[,ili[[ii]]] - Xbeg4[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABgf",ii), mean(abs(colMeans_(Xb0[,ili[[ii]]] - Xbegf[,ili[[ii]]], trim = trim))))
    }
  }
  
  abc <- c(1,2,3,4,"f")
  for(jj in abc){
    for(ii in names(ili)){ 
      assign(paste0("gamma",ii),c(get(paste0("AvMSEg0",ii)),get(paste0("AvABg0",ii)),get(paste0("AvMSEg",jj,ii)),get(paste0("AvABg",jj,ii))))
    }
    res <- rbind(res,get(paste0("gamma",ii)))
  }
  
  # Where are here!
  
  if("MSE" %in% out) for(ii in names(ili)){res <- c(res, get(paste0("AvMSE",ii))); nam <- c(nam, paste0("AvMSE",ii)) }
  
  if("MSE" %in% out) for(ii in names(ili)){res <- c(res, get(paste0("AvMSE",ii))); nam <- c(nam, paste0("AvMSE",ii)) }
  if("AB" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvAB",ii))); nam <- c(nam, paste0("AvAB",ii)) }
  names(res) <- nam
  
  if("time" %in% out){ nam <- c(names(res),"AvTime"); res <- c(res,mean(Xti)); names(res) <- nam }
  if("PVS" %in% out){ nam <- c("PVS",names(res)); res <- c(nrow(X)/nrow(Xor),res); names(res) <- nam }
  
  return(c(format(round(res,4),nsmall = 4,scientific = F)))
}