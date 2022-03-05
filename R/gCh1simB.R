# Code
# File ft_post.R: Post-processing for Simulation examples (Chapter 2)
# Code written by: Camilo Cardenas-Hurtado (c.a.cardenas-hurtado@lse.ac.uk)

Ch1pos <- function(FCOL,form,p){
  
  # TRIM <- F
  
  # Labels for parameters
  # ~~~~~~~~~~~~~~~~~~~~~
  nam <- NULL
  for(i in names(form)){ nam <- rbind(nam,expand.grid(i,1:p,stringsAsFactors = F)) }
  nam <- unlist(lapply(1:nrow(nam),function(i) paste0(nam[i,], collapse = "")))
  
  nM <- NULL
  for(i in 1:p){
   for(j in names(form)){
    nM <- append(nM,paste0(nam[grepl(j,as.character(nam))][i],".",c(0,seq_len(length(all.vars(as.formula(form[[j]])))))))
   }
  }; rm(nam)
  
  # Index for penalised parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ip <- !grepl(".0",nM,fixed = T)
  li1 <- c("mu1.2","mu2.1","sigma1.2","sigma2.1","tau1.2","tau2.1","nu1.2","nu2.1")
  li2 <- lapply(li1,function(i) !grepl(i,nM))
  for(i in 1:length(li2)){ ip <- ip & li2[[i]] }; names(ip) <- nM; # rm(list=c("li1","li2"))
  ip[grepl(".0",nM,fixed = T)] <- T; nM <- nM[ip]
  
  # Matrices of parameters (true, unpenalised and penalised)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # nX <- FCOL[!complete.cases(FCOL),]
  X  <- FCOL[complete.cases(FCOL),]
  X <- X[apply(X,1,min)!=-999,]
  ix <- mean(X[,ncol(X)]) # else ix <- c(colMeans(X[,(ncol(X)-1):ncol(X)]))
  X <- X[,-ncol(X)] # else X <- X[,-c((ncol(X)-1):ncol(X))]
  rR <- !(X[,ncol(X)] == 500 | X[,ncol(X)] == 300) #} else {
  #   rR <- !(X[,ncol(X)-1] == 500 | X[,ncol(X)-1] == 300);
  #   rA <- !(X[,ncol(X)] == 500 | X[,ncol(X)] == 300)
  # }
  # if(!CM){ 
  Xorg <- X[rR,seq_len(ix[1])[ip]]; colnames(Xorg) <- nM;
    #        Xorg <- Xorg[,!(colnames(Xorg) %in% li1)] } else {
    # Xorg <- X[rR,1:ix[1]];  XorA <- X[rA,1:ix[1]]; colnames(Xorg) <- colnames(XorA) <- nM;
    # Xorg <- Xorg[,!(colnames(Xorg) %in% li1)]; XorA <- XorA[,!(colnames(XorA) %in% li1)] }
  Xest <- X[rR,((ix[1]+1):(2*ix[1]))[ip]]; colnames(Xest) <- nM #; Xest <- Xest[,!(colnames(Xest) %in% li1)];
  # if(CM){ Xalt <- X[rA,(2*ix[1]+1):(2*ix[1]+ix[2])]; }
  # if(!CM){ 
  # Xinf <- X[rR,((2*ix[1])+1):ncol(X)]; colnames(Xinf) <- c("GBIC", "GIC","Iter") #} else {
    # Xinf1 <- X[rR,c(((2*ix[1]+ix[2])+1):(ncol(X)-4),ncol(X)-1)];
    # Xinf2 <- X[rA,c(((2*ix[1]+ix[2])+3):(ncol(X)-2),ncol(X))];
    # colnames(Xinf1) <- c("GBIC-1", "GIC-1","Iter1");
    # colnames(Xinf2) <- c("GBIC-2", "GIC-2","Iter2") }
  
  imu0 <- grepl("mu",colnames(Xorg),fixed = T) & grepl(".0",colnames(Xorg),fixed = T)
  imul <- grepl("mu",colnames(Xorg),fixed = T) & !grepl(".0",colnames(Xorg),fixed = T)
  isg0 <- grepl("sigma",colnames(Xorg),fixed = T) & grepl(".0",colnames(Xorg),fixed = T)
  isgl <- grepl("sigma",colnames(Xorg),fixed = T) & !grepl(".0",colnames(Xorg),fixed = T)
  
  # isg <- grepl("sigma",colnames(Xorg),fixed = T)# & ip
  # ita <- grepl("tau",colnames(Xorg),fixed = T)# & ip
  # inu <- grepl("nu",colnames(Xorg),fixed = T)# & ip
  
  ili <- list("mu0" = imu0, "mu1" = imul,
              "sg0" = isg0, "sg1" = isgl)
  
  # if(TRIM == T){
  #   AB <- MSE <- res <- NULL
  #   for(i in 1:ncol(Xorg)){
  #    AB <- c(AB,abs(colMeans(Xorg)[i] - mean(Xest[(Xest[,i] > quantile(Xest[,i],.05)) & (Xest[,i] < quantile(Xest[,i],.95)),i])))
  #    MSE <- c(MSE,mean((colMeans(Xorg)[i] - Xest[(Xest[,i] > quantile(Xest[,i],.05)) & (Xest[,i] < quantile(Xest[,i],.95)),i])^2))
  #   }
  #   for(ii in names(ili)){
  #    assign(paste0("AvAB",ii), mean(AB[ili[[ii]]]) )
  #    assign(paste0("AvMSE",ii), mean(colMeans((Xorg[,ili[[ii]]] - Xest[,ili[[ii]]])^2)) )
  #   }
  #   res <- c(res, get(paste0("AvMSE",ii)), get(paste0("AvAB",ii)))
  # } else {
  
  res <- NULL
  for(ii in names(ili)){
   assign(paste0("AvAB",ii), mean(abs(colMeans(Xorg[,ili[[ii]]] - Xest[,ili[[ii]]]))) )
   assign(paste0("AvMSE",ii), mean(colMeans((Xorg[,ili[[ii]]] - Xest[,ili[[ii]]])^2)) ) } #}
  
  for(ii in names(ili)){
    res <- c(res, get(paste0("AvMSE",ii)))
  }
  
  for(ii in names(ili)){
    res <- c(res, get(paste0("AvAB",ii)))
  }

  # if(print){
  #   paste0(round(res,4)," & ")
  # }
  
  # if(CM){
  #  for(ii in names(ili)){
  #  assign(paste0("AvAB2",ii), mean(abs(colMeans(XorA[,ili[[ii]]] - Xalt[,ili[[ii]]]))) )
  #  assign(paste0("AvMSE2",ii), mean(colMeans((XorA[,ili[[ii]]] - Xalt[,ili[[ii]]])^2)) )
  #  res <- c(res, get(paste0("AvMSE2",ii)), get(paste0("AvAB2",ii))) }
  # }
    
  # if(!CM){ GInC <- colMeans(Xinf[,-ncol(Xinf)]); conv <- sum(rR)/nrow(X) } else{
  #   GInC <- c(colMeans(Xinf1[,-ncol(Xinf1)]), colMeans(Xinf2[,-ncol(Xinf2)]));
  #   conv <- sum(rR)/nrow(X) }
  # 
  # if(!infC){
    return(c(round(res,4)))
  # } else return(c(round(res,4),round(c(GInC),2),round(conv*100,1)))
}
