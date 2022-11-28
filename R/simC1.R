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
                c0 <- Sys.time()
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, nQP = 25)
                c1 <- Sys.time()
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ex2 <- c(c1-c0, B$iter, length(lb2cb(B$b)))
                return(c(cof,ex2)) },
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
                c0 <- Sys.time()
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, nQP = 25)
                c1 <- Sys.time()
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ex2 <- c(c1-c0, B$iter, length(lb2cb(B$b)))
                return(c(cof,ex2)) },
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
                c0 <- Sys.time()
                B <- glvmlss(data = A$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, nQP = 40)
                c1 <- Sys.time()
                cof <- c(lb2cb(A$b), lb2cb(B$b))
                ex2 <- c(c1-c0, B$iter, length(lb2cb(B$b)))
                return(c(cof,ex2)) },
                error = function(e) { return(NA) }  ) } )
  stopCluster(cl)
  if(saveRes) saveRDS(FCOL, file = paste0("C1E3n",n, "p", p,"_",Sys.Date(),".Rds"))
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
  # p <- as.numeric(sub(".Rds.*","",sub(".*p","", file)))
  EX <- as.numeric(sub(".*R/C1E","",sub("n.*","", file)))
  # n <- as.numeric(sub("p.*","",sub(".*n","", file)))
  
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
  if(q != 1) li1 <- c("mu1.2","sigma1.2","tau1.2","nu1.2") else li1 <- c("NA") # Depends on the type of identification!
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
  rR <- !(X[,ncol(X)] == 330)
  Xb0 <- X[rR,seq_len(ix[1])[ip]]; colnames(Xb0) <- nM ;
  Xbe <- X[rR,((ix[1]+1):(2*ix[1]))[ip]]; colnames(Xbe) <- nM
  Xit <- X[, max((ix[1]+1):(2*ix[1]))+2]
  Xti <- X[, max((ix[1]+1):(2*ix[1]))+1]
  
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
  if("PVS" %in% out){ nam <- c("PVS",names(res)); res <- c(nrow(Xbe)/nrow(Xor),res); names(res) <- nam }
  
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
  
  return(res)
  # return(c(format(round(res,4),nsmall = 4,scientific = F)))
}
