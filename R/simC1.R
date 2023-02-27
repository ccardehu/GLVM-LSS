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
                A <- glvmlss_sim(n, famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, nQP = 25)
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ste <- lb2cb(B$SE)
                ex2 <- c(B$iter, length(lb2cb(B$b)))
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
                A <- glvmlss_sim(n,famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, nQP = 25)
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ste <- lb2cb(B$SE)
                ex2 <- c(B$iter, length(lb2cb(B$b)))
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
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, nQP = 40)
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ste <- lb2cb(B$SE)
                ex2 <- c(B$iter, length(lb2cb(B$b)))
                return(c(cof,ste,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("C1E3n",n, "p", p,"_",Sys.Date(),".Rds"))
  return(FCOL)
}

glvmlss_parsimE4_appSE <- function(nsim, saveRes = T){
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
                             est.ci = "Approximate", solver = "nlminb", iter.lim = 1000)
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ste <- lb2cb(B$SE)
                ex2 <- c(B$iter, length(lb2cb(B$b)))
                return(c(cof,ste,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("C1E4(CP-AppSE)n",n, "p", p,"_",Sys.Date(),".Rds"))
  return(FCOL)
}

glvmlss_parsimE4_stdSE <- function(nsim, saveRes = T){
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
                             solver = "nlminb", iter.lim = 1000, EM_appHess = T)
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ste <- lb2cb(B$SE)
                ex2 <- c(B$iter, length(lb2cb(B$b)))
                return(c(cof,ste,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("C1E4(CP-StdSE)n",n, "p", p,"_",Sys.Date(),".Rds"))
  return(FCOL)
}

glvmlss_parsimpost <- function(file, out = c("MSE","AB"), trim = 0,
                               mu.eq = ~ Z1+Z2, sg.eq = NULL,
                               nu.eq = NULL, ta.eq = NULL, plot = F, outs = T){
  
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
  # if(q != 1) li1 <- c("mu1.2","mu2.1","sigma1.2","sigma2.1","tau1.2","tau2.1","nu1.2","nu2.1") else li1 <- c("NA")
  if(q != 1) li1 <- c("mu1.2","sigma1.2","tau1.2","nu1.2") else li1 <- c("NA") # Depends on the type of identification!
  li2 <- lapply(li1, function(i) !grepl(i,nM))
  for(i in 1:length(li2)){ ip <- ip & li2[[i]]}; names(ip) <- nM; rm(list=c("li1","li2"))
  ip[grepl(".0",nM,fixed = T)] <- T; nM <- nM[ip];
  
  for(i in 1:q){
    tmp <- grepl(paste0(".",i),nM,fixed = T)
    names(tmp) <- nM
    tmp <- tmp[nM]
    assign(paste0("ip",i), tmp); #rm(tmp)
  }
  
  Xor <- X
  ix <- mean(X[,ncol(X)], na.rm = T)
  probseedX <- is.na(X[,ncol(X)])
  X  <- X[!probseedX,];
  X <- X[,-ncol(X)]
  
  rR <- !(X[,ncol(X)] == 330)
  Xb0 <- X[rR,seq_len(ix[1])[ip]]; colnames(Xb0) <- nM ;
  Xbe <- X[rR,((ix[1]+1):(2*ix[1]))[ip]]; colnames(Xbe) <- nM
  Xse <- X[rR,((2*ix[1]+1):(3*ix[1]))[ip]]; colnames(Xse) <- nM
  Xit <- X[rR, ncol(X)]
  
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
    }
  }
  
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
  if("PVS" %in% out){ nam <- c("PVS",names(res)); res <- c(nrow(Xbe)/nrow(Xor),res); names(res) <- nam }
  
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
