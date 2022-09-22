glvmlss_parsimpost_pml <- function(file,
                                   out = c("MSE","AB"), trim = 0, dgs = 1,
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
  p <- as.numeric(sub(".Rds.*","", sub(".*0p","", file)))
  # p <- as.numeric(sub("C2.*","", sub(".*0p","", file))) # When using _trust/optim files
  
  nam <- NULL
  for(i in names(form)){ nam <- rbind(nam,expand.grid(i,1:p,stringsAsFactors = F)) }
  nam <- unlist(lapply(1:nrow(nam),function(i) paste0(nam[i,], collapse = "")))
  
  nM <- NULL
  for(i in 1:p){
    for(j in names(form)){
      nM <- append(nM,paste0(nam[grepl(j,as.character(nam))][i],".",c(0,seq_len(length(all.vars(as.formula(form[[j]])))))))
    }
  }; #rm(nam)
  
  ip <- !grepl(".0",nM,fixed = T)
  if(q != 1) li1 <- c("mu1.2","mu2.1","sigma1.2","sigma2.1","tau1.2","tau2.1","nu1.2","nu2.1") else li1 <- c("NA")
  li2 <- lapply(li1, function(i) !grepl(i,nM))
  for(i in 1:length(li2)){ ip <- ip & li2[[i]]}; names(ip) <- nM; rm(list=c("li1","li2"))
  ip[grepl(".0",nM,fixed = T)] <- T; nM <- nM[ip];
  
  for(i in 1:q){
    tmp <- grepl(paste0(".",i),nM,fixed = T)
    names(tmp) <- nM
    tmp <- tmp[nM]
    assign(paste0("ip",i), tmp); #rm(tmp)
  }
  
  abc <- c(1,2,3,4,"f") # Influence factor values
  Xor <- X
  X  <- X[complete.cases(X),];
  ix <- mean(X[,ncol(X)])
  X <- X[,-ncol(X)]
  
  Xl <- X[,(ncol(X)-(length(form)*5)+1):ncol(X)]
  colnames(Xl) <- c(outer(paste0("lambda_",names(form)),abc,paste0))
  X <- X[,-((ncol(X)-(length(form)*5)+1):ncol(X))]
  
  Xc <- X[,(ncol(X)-5):ncol(X)]
  colnames(Xc) <- paste0("conv",c(0,abc))
  X <- X[,-((ncol(X)-5):ncol(X))]
  
  Xb0   <- X[, seq_len(ix[1])[ip]]; colnames(Xb0) <- nM ;
  Xbeg0 <- X[Xc[,"conv0"] == 1, ((ix[1]+1):(2*ix[1]))[ip]]; colnames(Xbeg0) <- nM
  Xbeg1 <- X[Xc[,"conv1"] == 1, ((2*ix[1]+1):(3*ix[1]))[ip]]; colnames(Xbeg1) <- nM
  Xbeg2 <- X[Xc[,"conv2"] == 1, ((3*ix[1]+1):(4*ix[1]))[ip]]; colnames(Xbeg2) <- nM
  Xbeg3 <- X[Xc[,"conv3"] == 1, ((4*ix[1]+1):(5*ix[1]))[ip]]; colnames(Xbeg3) <- nM
  Xbeg4 <- X[Xc[,"conv4"] == 1, ((5*ix[1]+1):(6*ix[1]))[ip]]; colnames(Xbeg4) <- nM
  Xbegf <- X[Xc[,"convf"] == 1, ((6*ix[1]+1):(7*ix[1]))[ip]]; colnames(Xbegf) <- nM
  XGBIC <- X[, seq(from = (7*ix[1])+1, length.out = 6)]
  colnames(XGBIC) <- c("GBICg0","GBICg1","GBICg2",
                       "GBICg3","GBICg4","GBICgf")
  
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
  
  # Fix sign indeterminacy
  # ~~~~~~~~~~~~~~~~~~~~~~
  
  for(g in c(0,abc)){
    Xtmp <- get(paste0("Xbeg",g))
    for(i in 1:nrow(Xtmp)){
      for(j in 1:q){
        for(ii in names(ili)[grepl("1",names(ili),fixed = T)]){
          ij <- which.max(abs(Xtmp[i, ili[[ii]] & get(paste0("ip",j)) & (Xb0[1,] != 0) ]))
          if(sign(Xtmp[i,ili[[ii]] & get(paste0("ip",j)) & (Xb0[1,] != 0)][ij]) != 
             sign(Xb0 [1,ili[[ii]] & get(paste0("ip",j)) & (Xb0[1,] != 0)][ij])){
            Xtmp[i,ili[[ii]] & get(paste0("ip",j))] <- -Xtmp[i,ili[[ii]] & get(paste0("ip",j))] }
        }
      }
    }
    assign(paste0("Xbeg",g), Xtmp)
  }; rm(Xtmp)
  
  # Fix simulations with penalised anchor items
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # ia <- NULL
  # for(f in names(form)){ ia <- c(ia, paste0(nam[grepl(f,nam,fixed = T)][1:q],".",1:q)) }
  # 
  # for(g in abc){
  #   tmp <- get(paste0("Xbeg",g))[,ia] - Xb0[Xc[,paste0("conv",g)] == 1,ia]
  #   tmpid <- sapply(1:nrow(tmp), function(i) any(abs(tmp)[i,] > 1))
  #   namestmp <- rownames(tmp[tmpid,])
  #   Xc[namestmp, paste0("conv",g)] <- 0
  # }
  # 
  # Xbeg1 <- X[Xc[,"conv1"] == 1, ((2*ix[1]+1):(3*ix[1]))[ip]]; colnames(Xbeg1) <- nM
  # Xbeg2 <- X[Xc[,"conv2"] == 1, ((3*ix[1]+1):(4*ix[1]))[ip]]; colnames(Xbeg2) <- nM
  # Xbeg3 <- X[Xc[,"conv3"] == 1, ((4*ix[1]+1):(5*ix[1]))[ip]]; colnames(Xbeg3) <- nM
  # Xbeg4 <- X[Xc[,"conv4"] == 1, ((5*ix[1]+1):(6*ix[1]))[ip]]; colnames(Xbeg4) <- nM
  # Xbegf <- X[Xc[,"convf"] == 1, ((6*ix[1]+1):(7*ix[1]))[ip]]; colnames(Xbegf) <- nM
  
  # ~~~~~~~~~~~~~~~~~~
  # Compute MSE and AB
  # ~~~~~~~~~~~~~~~~~~
  
  for(ii in names(ili)){
    if("MSE" %in% out){
      assign(paste0("AvMSEg0",ii), mean(colMeans_((Xb0[Xc[,"conv0"] == 1, ili[[ii]]] - Xbeg0[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEg1",ii), mean(colMeans_((Xb0[Xc[,"conv1"] == 1, ili[[ii]]] - Xbeg1[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEg2",ii), mean(colMeans_((Xb0[Xc[,"conv2"] == 1, ili[[ii]]] - Xbeg2[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEg3",ii), mean(colMeans_((Xb0[Xc[,"conv3"] == 1, ili[[ii]]] - Xbeg3[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEg4",ii), mean(colMeans_((Xb0[Xc[,"conv4"] == 1, ili[[ii]]] - Xbeg4[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEgf",ii), mean(colMeans_((Xb0[Xc[,"convf"] == 1, ili[[ii]]] - Xbegf[,ili[[ii]]])^2, trim = trim)))
    }
    
    if("AB" %in% out){
      assign(paste0("AvABg0",ii), mean(abs(colMeans_(Xb0[Xc[,"conv0"] == 1, ili[[ii]]] - Xbeg0[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABg1",ii), mean(abs(colMeans_(Xb0[Xc[,"conv1"] == 1, ili[[ii]]] - Xbeg1[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABg2",ii), mean(abs(colMeans_(Xb0[Xc[,"conv2"] == 1, ili[[ii]]] - Xbeg2[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABg3",ii), mean(abs(colMeans_(Xb0[Xc[,"conv3"] == 1, ili[[ii]]] - Xbeg3[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABg4",ii), mean(abs(colMeans_(Xb0[Xc[,"conv4"] == 1, ili[[ii]]] - Xbeg4[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABgf",ii), mean(abs(colMeans_(Xb0[Xc[,"convf"] == 1, ili[[ii]]] - Xbegf[,ili[[ii]]], trim = trim))))
    }
  }

  res <- matrix(NA, nrow = length(c(0,abc)), ncol = sum(c("MSE" %in% out, "AB" %in% out))*length(ili)+1); rownames(res) <- paste0("g",c(0,abc))
  colnames(res) <- c(c(outer("AvMSE",names(ili),paste0)),c(outer("AvAB",names(ili),paste0)),"AvGBIC")
  
  for(ii in names(ili)){
    for(jj in c(0,abc)){
      res[paste0("g",jj), paste0("AvMSE",ii)] <- get(paste0("AvMSEg",jj,ii))
      res[paste0("g",jj), paste0("AvAB",ii) ] <- get(paste0("AvABg",jj,ii))
    }
  }
  
  for(jj in c(0,abc)){ res[paste0("g",jj), "AvGBIC"] <- mean(XGBIC[Xc[,paste0("conv",jj)] == 1, paste0("GBICg",jj)]) }
  
  # ~~~~~~~~~~~~~~~~~~~~~
  # Compute CER, TPR, FPR
  # ~~~~~~~~~~~~~~~~~~~~~
  
  jlj <- list("u" = grepl(".0",colnames(Xb0),fixed = T), "p" = !grepl(".0",colnames(Xb0),fixed = T)) # Which coefficients are penalised
  for(i in 1:nrow(Xbeg0)){
    for(ii in names(ili)[grepl("1",names(ili),fixed = T)]){
      # Xbeg0[i, ili[[ii]]][ abs(Xbeg0[i, ili[[ii]] ]) < min( abs( Xb0[1, ili[[ii]] & Xb0[1,] != 0] ) / 2 ) ] <- 0  # 10^(-dgs)
      Xbeg0[i, ili[[ii]]][ abs(Xbeg0[i, ili[[ii]] ]) < 0.1 ] <- 0  #
    }
    # tidx <- Xbeg0[i,] < min(abs(Xb0[i, Xb0[i,] != 0 & jlj$p ])) / 2
    # Xbeg0[i,tidx] <- 0
  }
  for(i in abc){ assign(paste0("Xbeg",i), round(get(paste0("Xbeg",i)), dgs)) } # Rounding
  
  for(ii in names(ili)[grepl("1",names(ili),fixed = T)]){ # We can change this for "L" to denote loadings
    
    t0 <- (Xb0[1, jlj$p & ili[[ii]] ] != 0) # Penalised estimates that are different from zero (in paper \mathbb{T}, Section 3.9)
    
    for(jj in c(0,abc)){
      
      # Correct Estimation Rate (CER)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      assign(paste0("CER",jj,ii), sapply(1:nrow(get(paste0("Xbeg",jj))),
                                         function(i) sum((get(paste0("Xbeg",jj))[i, jlj$p & ili[[ii]] ] != 0) == t0)/sum(jlj$p & ili[[ii]])))
      assign(paste0("AvCER",jj,ii), mean(get(paste0("CER",jj,ii))))
      
      # True Positive Rate (TPR)
      # ~~~~~~~~~~~~~~~~~~~~~~~~
      assign(paste0("TPR",jj,ii), sapply(1:nrow(get(paste0("Xbeg",jj))),
                                         function(i) sum(get(paste0("Xbeg",jj))[i, (Xb0[1,] != 0) & ili[[ii]] ] != 0)/sum(Xb0[1,ili[[ii]]] != 0)))
      assign(paste0("AvTPR",jj,ii), mean(get(paste0("TPR",jj,ii))))
      
      # False Positive Rate (FPR)
      # ~~~~~~~~~~~~~~~~~~~~~~~~
      assign(paste0("FPR",jj,ii), sapply(1:nrow(get(paste0("Xbeg",jj))),
                                         function(i) sum(get(paste0("Xbeg",jj))[i, (Xb0[1,] == 0) & ili[[ii]] ] != 0)/sum(Xb0[1,ili[[ii]]] == 0)))
      assign(paste0("AvFPR",jj,ii), mean(get(paste0("FPR",jj,ii))))
    }
  }
  
  names <- c(paste0("AvCER", names(ili)[grepl("1",names(ili),fixed = T)]),
             paste0("AvTPR", names(ili)[grepl("1",names(ili),fixed = T)]),
             paste0("AvFPR", names(ili)[grepl("1",names(ili),fixed = T)]))
  rtmp <- matrix(NA, nrow = length(c(0,abc)), ncol = length(names)); rownames(rtmp) <- paste0("g",c(0,abc))
  colnames(rtmp) <- names; rm(names)
  
  for(ii in names(ili)[grepl("1",names(ili),fixed = T)]){
    for(jj in c(0,abc)){
      rtmp[paste0("g",jj), paste0("AvCER",ii)] <- get(paste0("AvCER",jj,ii))
      rtmp[paste0("g",jj), paste0("AvTPR",ii)] <- get(paste0("AvTPR",jj,ii))
      rtmp[paste0("g",jj), paste0("AvFPR",ii)] <- get(paste0("AvFPR",jj,ii))
    }
  }
  
  res <- cbind(res,rtmp);
   
  if("PVS" %in% out){ 
    rtmp <- matrix(NA,nrow = nrow(res), ncol = 1); rownames(rtmp) <- rownames(res)
    colnames(rtmp) <- "PVS"
    for(i in c(0,abc)){
      rtmp[paste0("g",i),] <- nrow(get(paste0("Xbeg",i)))/300# nrow(Xor)
    }
    res <- cbind(rtmp, res); rm(rtmp)
  }
  
  if("lambda" %in% out){ 
    rtmp <- matrix(0,nrow = nrow(res), ncol = length(form)); rownames(rtmp) <- rownames(res)
    colnames(rtmp) <- paste0("AvL",names(form))
    for(i in abc){
      rtmp[paste0("g",i),] <- colMeans(Xl[Xc[,paste0("conv",i)] == 1, grepl(i,colnames(Xl),fixed = T), drop = F])
    }
    res <- cbind(res,rtmp); rm(rtmp)
  }
  
  return(res)
  # return(c(format(round(res,4),nsmall = 4,scientific = F)))
}

glvmlss_parsimE1_pml <- function(nsim, saveRes = T){
  cores <- ifelse(parallel::detectCores() > 32, 32, parallel::detectCores())
  cl <- parallel::makeCluster(cores)
  clusterExport(cl, ls(.GlobalEnv))
  doSNOW::registerDoSNOW(cl)
  progress <- function(a) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), a)
  opts <- list(progress = progress)
  FCOL <- suppressWarnings(
    foreach(l = 1:nsim,
            .combine = rbind,
            .options.snow = opts) %dopar% {
              tryCatch({
                
                set.seed(l)
                A <- glvmlss_sim(n, famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
                
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, solver = "trust", start.val = A$b)
                C1a <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 1, solver = "trust")
                C1b <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 2, solver = "trust")
                C1c <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 3, solver = "trust")
                C1d <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 4, solver = "trust")
                C1e <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, solver = "trust")
                
                cof <- c(lb2cb(A$b), lb2cb(B$b), lb2cb(C1a$b), lb2cb(C1b$b), lb2cb(C1c$b), lb2cb(C1d$b), lb2cb(C1e$b))
                sds <- c(lb2cb(B$S), lb2cb(C1a$S), lb2cb(C1b$S), lb2cb(C1c$S), lb2cb(C1d$S), lb2cb(C1e$S))
                ex2 <- c(B$GBIC, C1a$GBIC, C1b$GBIC, C1c$GBIC, C1d$GBIC, C1e$GBIC,
                         B$conv, C1a$conv, C1b$conv, C1c$conv, C1d$conv, C1e$conv,
                         c(C1a$lambdahat), c(C1b$lambdahat), c(C1c$lambdahat),
                         c(C1d$lambdahat), c(C1e$lambdahat), l, length(lb2cb(B$b)))
                
                return(c(cof,sds,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("C2E1n",n, "p", p,".Rds"))
  return(FCOL)
}

glvmlss_parsimE2_pml <- function(nsim, saveRes = T){
  cores <- ifelse(parallel::detectCores() > 32, 32, parallel::detectCores())
  cl <- parallel::makeCluster(cores)
  clusterExport(cl, ls(.GlobalEnv))
  doSNOW::registerDoSNOW(cl)
  progress <- function(a) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), a)
  opts <- list(progress = progress)
  FCOL <- suppressWarnings(
    foreach(l = 1:nsim,
            .combine = rbind,
            .options.snow = opts) %dopar% {
              tryCatch({
                
                set.seed(l)
                A <- glvmlss_sim(n, famt, mu.eq = ~ Z1+Z2, start.val = lc)
                
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, solver = "trust", start.val = A$b)
                C1a <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 1, solver = "trust")
                C1b <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 2, solver = "trust")
                C1c <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 3, solver = "trust")
                C1d <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 4, solver = "trust")
                C1e <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, solver = "trust")
                
                cof <- c(lb2cb(A$b), lb2cb(B$b), lb2cb(C1a$b), lb2cb(C1b$b), lb2cb(C1c$b), lb2cb(C1d$b), lb2cb(C1e$b))
                sds <- c(lb2cb(B$S), lb2cb(C1a$S), lb2cb(C1b$S), lb2cb(C1c$S), lb2cb(C1d$S), lb2cb(C1e$S))
                ex2 <- c(B$GBIC, C1a$GBIC, C1b$GBIC, C1c$GBIC, C1d$GBIC, C1e$GBIC,
                         B$conv, C1a$conv, C1b$conv, C1c$conv, C1d$conv, C1e$conv,
                         c(C1a$lambdahat), c(C1b$lambdahat), c(C1c$lambdahat),
                         c(C1d$lambdahat), c(C1e$lambdahat), l, length(lb2cb(B$b)))
                
                return(c(cof,sds,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("C2E2n",n, "p", p,".Rds"))
  return(FCOL)
}

glvmlss_parsimE3_pml <- function(nsim, saveRes = T){
  cores <- ifelse(parallel::detectCores() > 32, 32, parallel::detectCores())
  cl <- parallel::makeCluster(cores)
  clusterExport(cl, ls(.GlobalEnv))
  doSNOW::registerDoSNOW(cl)
  progress <- function(a) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), a)
  opts <- list(progress = progress)
  FCOL <- suppressWarnings(
    foreach(l = 1:nsim,
            .combine = rbind,
            .options.snow = opts ) %dopar% { # 
              tryCatch({
                
                set.seed(l)
                A <- glvmlss_sim(n, famt, mu.eq = ~ Z1, sg.eq = ~ Z1, start.val = lc)
                
                B <- glvmlss(data = A$Y, family = famt,mu.eq = ~ Z1, sg.eq = ~ Z1, solver = "trust", start.val = A$b)
                C1a <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 1, solver = "trust")
                C1b <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 2, solver = "trust")
                C1c <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 3, solver = "trust")
                C1d <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, gamma = 4, solver = "trust")
                C1e <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, start.val = B$b,
                               penalty = "alasso", lambda = "auto", w.alasso = B$b, solver = "trust")
                
                cof <- c(lb2cb(A$b), lb2cb(B$b), lb2cb(C1a$b), lb2cb(C1b$b), lb2cb(C1c$b), lb2cb(C1d$b), lb2cb(C1e$b))
                sds <- c(lb2cb(B$S), lb2cb(C1a$S), lb2cb(C1b$S), lb2cb(C1c$S), lb2cb(C1d$S), lb2cb(C1e$S))
                ex2 <- c(B$GBIC, C1a$GBIC, C1b$GBIC, C1c$GBIC, C1d$GBIC, C1e$GBIC,
                         B$conv, C1a$conv, C1b$conv, C1c$conv, C1d$conv, C1e$conv,
                         c(C1a$lambdahat), c(C1b$lambdahat), c(C1c$lambdahat),
                         c(C1d$lambdahat), c(C1e$lambdahat), l, length(lb2cb(B$b)))
                
                return(c(cof,sds,ex2)) },
                error = function(e) { return(NA) }  ) } ) #
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("C2E3n",n, "p", p,".Rds"))
  return(FCOL)
}


# 
glvmlss_parsimpost_pml2 <- function(file,
                                   out = c("MSE","AB"), trim = 0, dgs = 1,
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
  p <- as.numeric(sub(".Rds.*","", sub(".*0p","", file)))

  nam <- NULL
  for(i in names(form)){ nam <- rbind(nam,expand.grid(i,1:p,stringsAsFactors = F)) }
  nam <- unlist(lapply(1:nrow(nam),function(i) paste0(nam[i,], collapse = "")))
  
  nM <- NULL
  for(i in 1:p){
    for(j in names(form)){
      nM <- append(nM,paste0(nam[grepl(j,as.character(nam))][i],".",c(0,seq_len(length(all.vars(as.formula(form[[j]])))))))
    }
  }; #rm(nam)
  
  ip <- !grepl(".0",nM,fixed = T)
  if(q != 1) li1 <- c("mu1.2","mu2.1","sigma1.2","sigma2.1","tau1.2","tau2.1","nu1.2","nu2.1") else li1 <- c("NA")
  li2 <- lapply(li1, function(i) !grepl(i,nM))
  for(i in 1:length(li2)){ ip <- ip & li2[[i]]}; names(ip) <- nM; rm(list=c("li1","li2"))
  ip[grepl(".0",nM,fixed = T)] <- T; nM <- nM[ip];
  
  for(i in 1:q){
    tmp <- grepl(paste0(".",i),nM,fixed = T)
    names(tmp) <- nM
    tmp <- tmp[nM]
    assign(paste0("ip",i), tmp); #rm(tmp)
  }
  
  abc <- c(1,2,3,4,"f") # Influence factor values
  Xor <- X
  ix <- mean(X[,ncol(X)], na.rm = T)
  X <- X[,-ncol(X)]
  X  <- X[!is.na(X[,ncol(X)]),];
  # Xseed <- setdiff(1:nrow(Xor), X[,ncol(X)])
  X <- X[,-ncol(X)] # -Xseed
  
  Xl <- X[,(ncol(X)-(length(form)*5)+1):ncol(X)]
  colnames(Xl) <- c(outer(paste0("lambda_",names(form)),abc,paste0))
  X <- X[,-((ncol(X)-(length(form)*5)+1):ncol(X))]
  
  Xc <- X[,(ncol(X)-5):ncol(X)]
  colnames(Xc) <- paste0("conv",c(0,abc))
  X <- X[,-((ncol(X)-5):ncol(X))]
  
  Xb0   <- X[, seq_len(ix[1])[ip]]; colnames(Xb0) <- nM ;
  Xbeg0 <- X[Xc[,"conv0"] == 1, ((ix[1]+1):(2*ix[1]))[ip]]; colnames(Xbeg0) <- nM
  Xbeg1 <- X[Xc[,"conv1"] == 1, ((2*ix[1]+1):(3*ix[1]))[ip]]; colnames(Xbeg1) <- nM
  Xbeg2 <- X[Xc[,"conv2"] == 1, ((3*ix[1]+1):(4*ix[1]))[ip]]; colnames(Xbeg2) <- nM
  Xbeg3 <- X[Xc[,"conv3"] == 1, ((4*ix[1]+1):(5*ix[1]))[ip]]; colnames(Xbeg3) <- nM
  Xbeg4 <- X[Xc[,"conv4"] == 1, ((5*ix[1]+1):(6*ix[1]))[ip]]; colnames(Xbeg4) <- nM
  Xbegf <- X[Xc[,"convf"] == 1, ((6*ix[1]+1):(7*ix[1]))[ip]]; colnames(Xbegf) <- nM
  
  XSeg0 <- X[Xc[,"conv0"] == 1, ((7*ix[1]+1):(8*ix[1]))[ip]]; colnames(XSeg0) <- nM
  XSeg1 <- X[Xc[,"conv1"] == 1, ((8*ix[1]+1):(9*ix[1]))[ip]]; colnames(XSeg1) <- nM
  XSeg2 <- X[Xc[,"conv2"] == 1, ((9*ix[1]+1):(10*ix[1]))[ip]]; colnames(XSeg2) <- nM
  XSeg3 <- X[Xc[,"conv3"] == 1, ((10*ix[1]+1):(11*ix[1]))[ip]]; colnames(XSeg3) <- nM
  XSeg4 <- X[Xc[,"conv4"] == 1, ((11*ix[1]+1):(12*ix[1]))[ip]]; colnames(XSeg4) <- nM
  XSegf <- X[Xc[,"convf"] == 1, ((12*ix[1]+1):(13*ix[1]))[ip]]; colnames(XSegf) <- nM
  
  XGBIC <- X[, seq(from = (13*ix[1])+1, length.out = 6)]
  colnames(XGBIC) <- c("GBICg0","GBICg1","GBICg2",
                       "GBICg3","GBICg4","GBICgf")
  
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
  
  # Fix sign indeterminacy
  # ~~~~~~~~~~~~~~~~~~~~~~
  
  for(g in c(0,abc)){
    Xtmp <- get(paste0("Xbeg",g))
    for(i in 1:nrow(Xtmp)){
      for(j in 1:q){
        for(ii in names(ili)[grepl("1",names(ili),fixed = T)]){
          ij <- which.max(abs(Xtmp[i, ili[[ii]] & get(paste0("ip",j)) & (Xb0[1,] != 0) ]))
          if(sign(Xtmp[i,ili[[ii]] & get(paste0("ip",j)) & (Xb0[1,] != 0)][ij]) != 
             sign(Xb0 [1,ili[[ii]] & get(paste0("ip",j)) & (Xb0[1,] != 0)][ij])){
            Xtmp[i,ili[[ii]] & get(paste0("ip",j))] <- -Xtmp[i,ili[[ii]] & get(paste0("ip",j))] }
        }
      }
    }
    assign(paste0("Xbeg",g), Xtmp)
  }; rm(Xtmp)
  
  # Fix simulations with penalised anchor items
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # ~~~~~~~~~~~~~~~~~~
  # Compute MSE and AB
  # ~~~~~~~~~~~~~~~~~~
  
  for(ii in names(ili)){
    if("MSE" %in% out){
      assign(paste0("AvMSEg0",ii), mean(colMeans_((Xb0[Xc[,"conv0"] == 1, ili[[ii]]] - Xbeg0[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEg1",ii), mean(colMeans_((Xb0[Xc[,"conv1"] == 1, ili[[ii]]] - Xbeg1[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEg2",ii), mean(colMeans_((Xb0[Xc[,"conv2"] == 1, ili[[ii]]] - Xbeg2[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEg3",ii), mean(colMeans_((Xb0[Xc[,"conv3"] == 1, ili[[ii]]] - Xbeg3[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEg4",ii), mean(colMeans_((Xb0[Xc[,"conv4"] == 1, ili[[ii]]] - Xbeg4[,ili[[ii]]])^2, trim = trim)))
      assign(paste0("AvMSEgf",ii), mean(colMeans_((Xb0[Xc[,"convf"] == 1, ili[[ii]]] - Xbegf[,ili[[ii]]])^2, trim = trim)))
    }
    
    if("AB" %in% out){
      assign(paste0("AvABg0",ii), mean(abs(colMeans_(Xb0[Xc[,"conv0"] == 1, ili[[ii]]] - Xbeg0[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABg1",ii), mean(abs(colMeans_(Xb0[Xc[,"conv1"] == 1, ili[[ii]]] - Xbeg1[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABg2",ii), mean(abs(colMeans_(Xb0[Xc[,"conv2"] == 1, ili[[ii]]] - Xbeg2[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABg3",ii), mean(abs(colMeans_(Xb0[Xc[,"conv3"] == 1, ili[[ii]]] - Xbeg3[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABg4",ii), mean(abs(colMeans_(Xb0[Xc[,"conv4"] == 1, ili[[ii]]] - Xbeg4[,ili[[ii]]], trim = trim))))
      assign(paste0("AvABgf",ii), mean(abs(colMeans_(Xb0[Xc[,"convf"] == 1, ili[[ii]]] - Xbegf[,ili[[ii]]], trim = trim))))
    }
  }

  res <- matrix(NA, nrow = length(c(0,abc)), ncol = sum(c("MSE" %in% out, "AB" %in% out))*length(ili)+1); rownames(res) <- paste0("g",c(0,abc))
  colnames(res) <- c(c(outer("AvMSE",names(ili),paste0)),c(outer("AvAB",names(ili),paste0)),"AvGBIC")
  
  for(ii in names(ili)){
    for(jj in c(0,abc)){
      res[paste0("g",jj), paste0("AvMSE",ii)] <- get(paste0("AvMSEg",jj,ii))
      res[paste0("g",jj), paste0("AvAB",ii) ] <- get(paste0("AvABg",jj,ii))
    }
  }
  
  for(jj in c(0,abc)){ res[paste0("g",jj), "AvGBIC"] <- mean(XGBIC[Xc[,paste0("conv",jj)] == 1, paste0("GBICg",jj)]) }
  
  # ~~~~~~~~~~~~~~~~~~~~~
  # Compute CER, TPR, FPR
  # ~~~~~~~~~~~~~~~~~~~~~
  
  jlj <- list("u" = grepl(".0",colnames(Xb0),fixed = T), "p" = !grepl(".0",colnames(Xb0),fixed = T)) # Which coefficients are penalised
  for(i in 1:nrow(Xbeg0)){
    for(ii in names(ili)[grepl("1",names(ili),fixed = T)]){
      # Xbeg0[i, ili[[ii]]][ abs(Xbeg0[i, ili[[ii]] ]) < min( abs( Xb0[1, ili[[ii]] & Xb0[1,] != 0] ) / 2 ) ] <- 0  # 10^(-dgs)
      Xbeg0[i, ili[[ii]]][ abs(Xbeg0[i, ili[[ii]] ]) < 0.1 ] <- 0  #
    }
    # tidx <- Xbeg0[i,] < min(abs(Xb0[i, Xb0[i,] != 0 & jlj$p ])) / 2
    # Xbeg0[i,tidx] <- 0
  }
  for(i in abc){ assign(paste0("Xbeg",i), round(get(paste0("Xbeg",i)), dgs)) } # Rounding
  
  for(ii in names(ili)[grepl("1",names(ili),fixed = T)]){ # We can change this for "L" to denote loadings
    
    t0 <- (Xb0[1, jlj$p & ili[[ii]] ] != 0) # Penalised estimates that are different from zero (in paper \mathbb{T}, Section 3.9)
    
    for(jj in c(0,abc)){
      
      # Correct Estimation Rate (CER)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      assign(paste0("CER",jj,ii), sapply(1:nrow(get(paste0("Xbeg",jj))),
                                         function(i) sum((get(paste0("Xbeg",jj))[i, jlj$p & ili[[ii]] ] != 0) == t0)/sum(jlj$p & ili[[ii]])))
      assign(paste0("AvCER",jj,ii), mean(get(paste0("CER",jj,ii))))
      
      # True Positive Rate (TPR)
      # ~~~~~~~~~~~~~~~~~~~~~~~~
      assign(paste0("TPR",jj,ii), sapply(1:nrow(get(paste0("Xbeg",jj))),
                                         function(i) sum(get(paste0("Xbeg",jj))[i, (Xb0[1,] != 0) & ili[[ii]] ] != 0)/sum(Xb0[1,ili[[ii]]] != 0)))
      assign(paste0("AvTPR",jj,ii), mean(get(paste0("TPR",jj,ii))))
      
      # False Positive Rate (FPR)
      # ~~~~~~~~~~~~~~~~~~~~~~~~
      assign(paste0("FPR",jj,ii), sapply(1:nrow(get(paste0("Xbeg",jj))),
                                         function(i) sum(get(paste0("Xbeg",jj))[i, (Xb0[1,] == 0) & ili[[ii]] ] != 0)/sum(Xb0[1,ili[[ii]]] == 0)))
      assign(paste0("AvFPR",jj,ii), mean(get(paste0("FPR",jj,ii))))
    }
  }
  
  names <- c(paste0("AvCER", names(ili)[grepl("1",names(ili),fixed = T)]),
             paste0("AvTPR", names(ili)[grepl("1",names(ili),fixed = T)]),
             paste0("AvFPR", names(ili)[grepl("1",names(ili),fixed = T)]))
  rtmp <- matrix(NA, nrow = length(c(0,abc)), ncol = length(names)); rownames(rtmp) <- paste0("g",c(0,abc))
  colnames(rtmp) <- names; rm(names)
  
  for(ii in names(ili)[grepl("1",names(ili),fixed = T)]){
    for(jj in c(0,abc)){
      rtmp[paste0("g",jj), paste0("AvCER",ii)] <- get(paste0("AvCER",jj,ii))
      rtmp[paste0("g",jj), paste0("AvTPR",ii)] <- get(paste0("AvTPR",jj,ii))
      rtmp[paste0("g",jj), paste0("AvFPR",ii)] <- get(paste0("AvFPR",jj,ii))
    }
  }
  
  res <- cbind(res,rtmp);
   
  if("PVS" %in% out){ 
    rtmp <- matrix(NA,nrow = nrow(res), ncol = 1); rownames(rtmp) <- rownames(res)
    colnames(rtmp) <- "PVS"
    for(i in c(0,abc)){
      rtmp[paste0("g",i),] <- nrow(get(paste0("Xbeg",i)))/nrow(Xor)
    }
    res <- cbind(rtmp, res); rm(rtmp)
  }
  
  if("lambda" %in% out){ 
    rtmp <- matrix(0,nrow = nrow(res), ncol = length(form)); rownames(rtmp) <- rownames(res)
    colnames(rtmp) <- paste0("AvL",names(form))
    for(i in abc){
      rtmp[paste0("g",i),] <- colMeans(Xl[Xc[,paste0("conv",i)] == 1, grepl(i,colnames(Xl),fixed = T), drop = F])
    }
    res <- cbind(res,rtmp); rm(rtmp)
  }
  
  return(res)
}

