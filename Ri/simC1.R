glvmlss_parsimE1 <- function(nsim, saveRes = T){
  cores <- ifelse(parallel::detectCores() > 32, 32, parallel::detectCores())
  cl <- parallel::makeCluster(cores)
  clusterExport(cl, ls(.GlobalEnv))
  doSNOW::registerDoSNOW(cl)
  progress <- function(a) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), a)
  opts <- list(progress = progress)
  FCOL <- suppressWarnings(
    foreach(l = 1:nsim,
            .combine = rbind,
            .options.snow = opts) %dopar% { # , .export = ls(parent.frame()) 
              tryCatch({
                
                set.seed(l)
                A <- glvmlss_sim(n, famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc,iden.res = "eiv", Rz = matrix(c(1,.3,.3,1), nrow = 2))
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, nQP = 25, iden.res = "eiv", corr.lv = T,
                             start.val = A$b, Rz = matrix(c(1,.3,.3,1), nrow = 2))
                cof <- c(lb2cb(A$b), A$Rz[2,1], lb2cb(B$b), B$Rz[2,1])
                ste <- c(lb2cb(B$SE$b),B$SE$Rz[1,2])
                ex2 <- c(1,B$conv, B$iter, length(lb2cb(B$b)) + 1) # include the estimated correlation 
                return(c(cof,ste,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("C1E1n",n, "p", p,"_",Sys.Date(),".Rds"))
  return(FCOL)
}

glvmlss_parsimE2 <- function(nsim, saveRes = T){
  cores <- ifelse(parallel::detectCores() > 32, 32, parallel::detectCores())
  cl <- parallel::makeCluster(cores)
  clusterExport(cl, ls(.GlobalEnv))
  doSNOW::registerDoSNOW(cl)
  progress <- function(a) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), a)
  opts <- list(progress = progress)
  FCOL <- suppressWarnings(
    foreach(l = 1:nsim,
            .combine = rbind,
            .options.snow = opts) %dopar% { #, .export = ls(parent.frame()) 
              tryCatch({
                
                set.seed(l)
                A <- glvmlss_sim(n,famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc,iden.res = "eiv", Rz = matrix(c(1,.45,.45,1), nrow = 2))
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, nQP = 25, iden.res = "eiv", corr.lv = T, start.val = A$b)
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ste <- lb2cb(B$SE)
                ex2 <- c(1,B$conv, B$iter, length(lb2cb(B$b)))
                return(c(cof,ste,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("C1E2n",n, "p", p,"_",Sys.Date(),".Rds"))
  return(FCOL)
}

glvmlss_parsimE3 <- function(nsim, saveRes = T){
  cores <- ifelse(parallel::detectCores() > 32, 32, parallel::detectCores())
  cl <- parallel::makeCluster(cores)
  clusterExport(cl, ls(.GlobalEnv))
  doSNOW::registerDoSNOW(cl)
  progress <- function(a) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), a)
  opts <- list(progress = progress)
  FCOL <- suppressWarnings(
    foreach(l = 1:nsim,
            .combine = rbind,
            .options.snow = opts) %dopar% { # .export = ls(parent.frame()),
              tryCatch({
                
                set.seed(l)
                A <- glvmlss_sim(n,famt, mu.eq = ~ Z1, sg.eq = ~ Z1, start.val = lc)
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, nQP = 50, start.val = A$b)
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ste <- lb2cb(B$SE)
                ex2 <- c(0,B$conv, B$iter, length(lb2cb(B$b)))
                return(c(cof,ste,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("C1E3n",n, "p", p,"_",Sys.Date(),".Rds"))
  return(FCOL)
}

glvmlss_parsimE4 <- function(nsim, saveRes = T){
  cores <- ifelse(parallel::detectCores() > 32, 32, parallel::detectCores())
  cl <- parallel::makeCluster(cores)
  clusterExport(cl, ls(.GlobalEnv))
  doSNOW::registerDoSNOW(cl)
  progress <- function(a) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), a)
  opts <- list(progress = progress)
  FCOL <- suppressWarnings(
    foreach(l = 1:nsim,
            .combine = rbind,
            .options.snow = opts) %dopar% { # .export = ls(parent.frame()),
              tryCatch({
                
                set.seed(l)
                A <- glvmlss_sim(n,famt, mu.eq = ~ Z1, sg.eq = ~ Z1, nu.eq = ~ Z1, start.val = lc)
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, nu.eq = ~ Z1, nQP = 100,
                             est.ci = "Approximate", solver = "nlminb", iter.lim = 1000, start.val = A$b)
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ste <- lb2cb(B$SE)
                ex2 <- c(0, B$conv, B$iter, length(lb2cb(B$b)))
                return(c(cof,ste,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("C1E4(CP-AppSE)n",n, "p", p,"_",Sys.Date(),".Rds"))
  return(FCOL)
}

glvmlss_parsimEConf <- function(nsim, saveRes = T){
  cores <- ifelse(parallel::detectCores() > 32, 32, parallel::detectCores())
  cl <- parallel::makeCluster(cores)
  clusterExport(cl, ls(.GlobalEnv))
  doSNOW::registerDoSNOW(cl)
  progress <- function(a) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), a)
  opts <- list(progress = progress)
  FCOL <- suppressWarnings(
    foreach(l = 1:nsim,
            .combine = rbind,
            .options.snow = opts) %dopar% { # .export = ls(parent.frame()),
              tryCatch({
                
                set.seed(l)
                A <- glvmlss_sim(n,famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, nu.eq = ~ Z1+Z2, start.val = lc, Rz = Rz, iden.res = iRes)
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, nu.eq = ~ Z1+Z2, EM_use2d = F,
                             EM_iter = 500, iter.lim = 1000, start.val = A$b, iden.res = iRes, corr.lv = T, Rz = Rz)
                cof <- c(lb2cb(A$b)[lb2cb(B$rb)], A$Rz[lower.tri(A$Rz, diag = F)], lb2cb(B$b)[lb2cb(B$rb)], B$Rz[lower.tri(B$Rz, diag = F)])
                ste <- c(lb2cb(B$SE$b)[lb2cb(B$rb)], B$SE$Rz[lower.tri(B$Rz, diag = F)])
                ex2 <- c(n, B$conv, length(lb2cb(B$b)[lb2cb(B$rb)]) + length(B$Rz[lower.tri(B$Rz, diag = F)]))
                return(c(cof,ste,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("C1E5_Conf_n",n,"_",Sys.Date(),".Rds"))
  return(FCOL)
}


glvmlss_parsimpost <- function(file, out = c("MSE","AB"), trim = 0, iden.res = "eiv",
                               mu.eq = ~ Z1+Z2, sg.eq = NULL,
                               nu.eq = NULL, ta.eq = NULL, plot = F, outs = F){

  X <- readRDS(file)
  
  form <- prep_form(mu.eq = mu.eq, sg.eq = sg.eq, nu.eq = nu.eq, ta.eq = ta.eq)
  q <- prep_Z(form = form)
  
  nam <- NULL
  for(i in names(form)){ nam <- rbind(nam,expand.grid(i,1:p,stringsAsFactors = F)) }
  nam <- unlist(lapply(1:nrow(nam),function(i) paste0(nam[i,], collapse = "")))
  
  nM <- NULL
  for(i in 1:p){
    for(j in names(form)){
      nM <- append(nM,paste0(nam[grepl(j,as.character(nam))][i],".",c(0,seq_len(length(all.vars(as.formula(form[[j]])))))))
    }
  };
  
  ip <- !grepl(".0",nM,fixed = T)
  if(q != 1 & iden.res == "eiv") li1 <- c("mu1.2","mu2.1","sigma1.2","sigma2.1","tau1.2","tau2.1","nu1.2","nu2.1") else li1 <- c("NA")
  if(q != 1 & iden.res == "recursive") li1 <- c("mu1.2","sigma1.2","tau1.2","nu1.2") else li1 <- c("NA") # Depends on the type of identification!
  li2 <- lapply(li1, function(i) !grepl(i,nM))
  for(i in 1:length(li2)){ ip <- ip & li2[[i]]}; names(ip) <- nM; rm(list=c("li1","li2"))
  ip[grepl(".0",nM,fixed = T)] <- T; nM <- nM[ip];
  
  for(i in 1:q){
    tmp <- grepl(paste0(".",i),nM,fixed = T)
    names(tmp) <- nM
    tmp <- tmp[nM]
    assign(paste0("ip",i), tmp); #rm(tmp)
  }
  
  # Original simulations
  Xor <- X
  # Cleaning number of estimated parameters
  ix <- mean(X[,ncol(X)], na.rm = T)
  probseedX <- is.na(X[,ncol(X)])
  X  <- X[!probseedX,];
  X <- X[,-ncol(X)]
  # Matrix for number of iterations
  Xit <- X[, ncol(X), drop = F]
  colnames(Xit) <- "iter"
  X <- X[,-ncol(X)]
  # Matrix for convergence status
  Xc <- X[, ncol(X), drop = F]
  colnames(Xc) <- "conv"
  X <- X[,-ncol(X)]
  # Matrix for LV correlation flag 
  Xlvc <- X[, ncol(X), drop = F]
  colnames(Xlvc) <- "LVcor"
  X <- X[,-ncol(X)]
  if(mean(Xlvc) == 1){
    corFlag <- T
    names(corFlag) <- "cor(Z1,Z2)"
    ip <- c(ip,corFlag)
    nM <- c(nM, names(corFlag))
    for(i in 1:q){
      tmp <- grepl(paste0(".",i),nM,fixed = T)
      names(tmp) <- nM
      tmp <- tmp[nM]
      assign(paste0("ip",i), c(get(paste0("ip",i)), !corFlag)); #rm(tmp)
    }
  }
  
  Xb0 <- X[Xc == 1,seq_len(ix[1])[ip]]; colnames(Xb0) <- nM
  Xbe <- X[Xc == 1,((ix[1]+1):(2*ix[1]))[ip]]; colnames(Xbe) <- nM
  Xse <- X[Xc == 1,((2*ix[1]+1):(3*ix[1]))[ip]]; colnames(Xse) <- nM
  Xit <- Xit[Xc == 1,]

  imu0 <- grepl("mu",colnames(Xb0),fixed = T) & grepl(".0",colnames(Xb0),fixed = T)
  imul <- grepl("mu",colnames(Xb0),fixed = T) & !grepl(".0",colnames(Xb0),fixed = T)
  isg0 <- grepl("sigma",colnames(Xb0),fixed = T) & grepl(".0",colnames(Xb0),fixed = T)
  isgl <- grepl("sigma",colnames(Xb0),fixed = T) & !grepl(".0",colnames(Xb0),fixed = T)
  ita0 <- grepl("tau",colnames(Xb0),fixed = T) & grepl(".0",colnames(Xb0),fixed = T)
  ital <- grepl("tau",colnames(Xb0),fixed = T) & !grepl(".0",colnames(Xb0),fixed = T)
  inu0 <- grepl("nu",colnames(Xb0),fixed = T) & grepl(".0",colnames(Xb0),fixed = T)
  inul <- grepl("nu",colnames(Xb0),fixed = T) & !grepl(".0",colnames(Xb0),fixed = T)
  icor <- grepl("cor",colnames(Xb0),fixed = T)
  ili <- list("mu0" = imu0, "mu1" = imul, "sg0" = isg0, "sg1" = isgl,
              "ta0" = ita0, "ta1" = ital, "nu0" = inu0, "nu1" = inul, "Rz" = icor)
  for(i in names(ili)){ if(all(!ili[[i]])) ili[[i]] <- NULL }
  
  # Fixing sign
  # ~~~~~~~~~~~
  
  for(i in 1:nrow(Xbe)){
    for(j in 1:q){
      flagj = 0
      for(ii in names(ili)[grepl("1",names(ili),fixed = T)]){
        ij <- which.max(abs(Xbe[i,ili[[ii]] & get(paste0("ip",j))]))
        if(sign(Xbe[i,ili[[ii]] & get(paste0("ip",j))][ij]) != sign(Xb0[i,ili[[ii]] & get(paste0("ip",j))][ij])){
          Xbe[i,ili[[ii]] & get(paste0("ip",j))] <- -Xbe[i,ili[[ii]] & get(paste0("ip",j))]
          flagj = flagj + 1
        }
      }
      if(flagj == 2) Xbe[i,"cor(Z1,Z2)"] <- -Xbe[i,"cor(Z1,Z2)"]
    }
  }
  
  res <- nam <- NULL
  for(ii in names(ili)){
    if("MSE" %in% out) assign(paste0("AvMSE",ii), mean(colMeans((Xb0[,ili[[ii]], drop = F] - Xbe[,ili[[ii]], drop = F])^2)) )
    if("AB" %in% out) assign(paste0("AvAB",ii), mean(abs(colMeans(Xb0[,ili[[ii]], drop = F] - Xbe[,ili[[ii]], drop = F]))) )
    if("SB" %in% out) assign(paste0("AvSB",ii), mean((colMeans(Xb0[,ili[[ii]], drop = F] - Xbe[,ili[[ii]], drop = F]))^2) )
    if("RB" %in% out) assign(paste0("AvRB",ii), mean((colMeans(Xb0[,ili[[ii]], drop = F] - Xbe[,ili[[ii]], drop = F]))/abs(colMeans(Xb0[,ili[[ii]]]))) )
  } 

  if("MSE" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvMSE",ii))); nam <- c(nam, paste0("AvMSE",ii)) }
  if("AB" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvAB",ii))); nam <- c(nam, paste0("AvAB",ii)) }
  if("SB" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvSB",ii))); nam <- c(nam, paste0("AvSB",ii)) }
  if("RB" %in% out) for(ii in names(ili)){ res <- c(res, get(paste0("AvRB",ii))); nam <- c(nam, paste0("AvRB",ii)) }
  names(res) <- nam
  
  if("iter" %in% out){ nam <- c(names(res),"AvIter"); res <- c(res,mean(Xit)); names(res) <- nam }
  if("PVS" %in% out){ nam <- c("PVS",names(res)); res <- c(mean(Xc),res); names(res) <- nam }

  if(plot){
    nid <- NULL
    for(ii in names(ili)){
      nid <- c(nid,which(ili[[ii]]))
    }
    boxplot(Xb0[,nid] - Xbe[,nid],
            main = paste0("Bias (by type of parameter): Simulation, n = ",n,", p = ",p),
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
  
  return(res)
}
