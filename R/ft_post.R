
# Code
# File ft_post.R: Post-processing for Simulation examples (Chapter 2)
# Code written by: Camilo Cardenas-Hurtado (c.a.cardenas-hurtado@lse.ac.uk)

post <- function(FCOL,form,p, round = T){
  
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
  
  ip <- !grepl(".0",nM,fixed = T);
  li1 <- c("mu1.2","mu2.1","sigma1.2","sigma2.1")
  li2 <- lapply(li1,function(i) !grepl(i,nM))
  for(i in 1:length(li2)){ ip <- ip & li2[[i]] }; names(ip) <- nM; rm(list=c("li1","li2"))
  
  # Matrices of parameters (true, unpenalised and penalised)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # nX <- FCOL[!complete.cases(FCOL),]
  X  <- FCOL[complete.cases(FCOL),]
  
  ix <- mean(X[,ncol(X)])
  a = 1; i0  <- apply(X[,((a-1)*ix+1):(a*ix)],1,min) != -999; b0  <- X[i0,((a-1)*ix+1):(a*ix)]; colnames(b0)   <- nM
  a = 2; iu  <- apply(X[,((a-1)*ix+1):(a*ix)],1,min) != -999; bu  <- X[iu,((a-1)*ix+1):(a*ix)]; colnames(bu)   <- nM
  a = 3; ip1 <- apply(X[,((a-1)*ix+1):(a*ix)],1,min) != -999; bp1 <- X[ip1,((a-1)*ix+1):(a*ix)]; colnames(bp1) <- nM
  a = 4; ip2 <- apply(X[,((a-1)*ix+1):(a*ix)],1,min) != -999; bp2 <- X[ip2,((a-1)*ix+1):(a*ix)]; colnames(bp2) <- nM
  a = 5; ip3 <- apply(X[,((a-1)*ix+1):(a*ix)],1,min) != -999; bp3 <- X[ip3,((a-1)*ix+1):(a*ix)]; colnames(bp3) <- nM
  a = 6; ip4 <- apply(X[,((a-1)*ix+1):(a*ix)],1,min) != -999; bp4 <- X[ip4,((a-1)*ix+1):(a*ix)]; colnames(bp4) <- nM
  a = 7; ip5 <- apply(X[,((a-1)*ix+1):(a*ix)],1,min) != -999; bp5 <- X[ip5,((a-1)*ix+1):(a*ix)]; colnames(bp5) <- nM
  
  if(round == T){
   rmlist <- list(b0 = b0,bu = bu,bp1 = bp1,bp2 = bp2,bp3 = bp3,bp4 = bp4,bp5 = bp5)
   names(rmlist) <- c("b0","bu","bp1","bp2","bp3","bp4","bp5")
   for(i in names(rmlist)){ assign(paste0(i), round(rmlist[[i]],3)) }
   rm(rmlist) }
  
  # Absolute Bias (AB) & MSE
  # ~~~~~~~~~~~~~~~~~~~~~~~~
  
  abmle  <- colMeans(abs(b0[iu,]-bu))
  abp1mle <- colMeans(abs(b0[ip1,]-bp1))
  abp2mle <- colMeans(abs(b0[ip2,]-bp2))
  abp3mle <- colMeans(abs(b0[ip3,]-bp3))
  abp4mle <- colMeans(abs(b0[ip4,]-bp4))
  abp5mle <- colMeans(abs(b0[ip5,]-bp5))
  
  msemle  <- colMeans((b0[iu,]-bu)^2)
  msep1mle <- colMeans((b0[ip1,]-bp1)^2)
  msep2mle <- colMeans((b0[ip2,]-bp2)^2)
  msep3mle <- colMeans((b0[ip3,]-bp3)^2)
  msep4mle <- colMeans((b0[ip4,]-bp4)^2)
  msep5mle <- colMeans((b0[ip5,]-bp5)^2)
  
  AvABm0 <- mean(abmle[!ip]);    AvABm1 <- mean(abmle[ip])
  AvABp0a <- mean(abp1mle[!ip]); AvABp1a <- mean(abp1mle[ip])
  AvABp0b <- mean(abp2mle[!ip]); AvABp1b <- mean(abp2mle[ip])
  AvABp0c <- mean(abp3mle[!ip]); AvABp1c <- mean(abp3mle[ip])
  AvABp0d <- mean(abp4mle[!ip]); AvABp1d <- mean(abp4mle[ip])
  AvABp0e <- mean(abp5mle[!ip]); AvABp1e <- mean(abp5mle[ip])
  
  AvMSEm0 <- mean(msemle[!ip]);    AvMSEm1 <- mean(msemle[ip])
  AvMSEp0a <- mean(msep1mle[!ip]); AvMSEp1a <- mean(msep1mle[ip])
  AvMSEp0b <- mean(msep2mle[!ip]); AvMSEp1b <- mean(msep2mle[ip])
  AvMSEp0c <- mean(msep3mle[!ip]); AvMSEp1c <- mean(msep3mle[ip])
  AvMSEp0d <- mean(msep4mle[!ip]); AvMSEp1d <- mean(msep4mle[ip])
  AvMSEp0e <- mean(msep5mle[!ip]); AvMSEp1e <- mean(msep5mle[ip])
  
  t0 <- unname(colMeans(b0) != 0)
  tp1 <- bp1 != 0; tp2 <- bp2 != 0; tp3 <- bp3 != 0; tp4 <- bp4 != 0; tp5 <- bp5 != 0
  
  # Correct Estimation Rate (CER)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CER1 <- sapply(1:nrow(tp1),function(i) sum(tp1[i,ip] == t0[ip])/sum(ip)); AvCER1 <- mean(CER1)
  CER2 <- sapply(1:nrow(tp2),function(i) sum(tp2[i,ip] == t0[ip])/sum(ip)); AvCER2 <- mean(CER2)
  CER3 <- sapply(1:nrow(tp3),function(i) sum(tp3[i,ip] == t0[ip])/sum(ip)); AvCER3 <- mean(CER3)
  CER4 <- sapply(1:nrow(tp4),function(i) sum(tp4[i,ip] == t0[ip])/sum(ip)); AvCER4 <- mean(CER4)
  CER5 <- sapply(1:nrow(tp5),function(i) sum(tp5[i,ip] == t0[ip])/sum(ip)); AvCER5 <- mean(CER5)
  
  # True Positive Rate (TPR)
  # ~~~~~~~~~~~~~~~~~~~~~~~~
  tset <- t0 & unname(ip)
  TPR1 <- sapply(1:nrow(bp1), function(i) sum(bp1[i,tset] != 0)/sum(tset)); AvTPR1 <- mean(TPR1)
  TPR2 <- sapply(1:nrow(bp2), function(i) sum(bp2[i,tset] != 0)/sum(tset)); AvTPR2 <- mean(TPR2)
  TPR3 <- sapply(1:nrow(bp3), function(i) sum(bp3[i,tset] != 0)/sum(tset)); AvTPR3 <- mean(TPR3)
  TPR4 <- sapply(1:nrow(bp4), function(i) sum(bp4[i,tset] != 0)/sum(tset)); AvTPR4 <- mean(TPR4)
  TPR5 <- sapply(1:nrow(bp5), function(i) sum(bp5[i,tset] != 0)/sum(tset)); AvTPR5 <- mean(TPR5)
   
  # False Positive Rate (FPR)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~
  fset <- !t0 & unname(ip)
  FPR1 <- sapply(1:nrow(bp1), function(i) sum(bp1[i,fset] != 0)/sum(fset)); AvFPR1 <- mean(FPR1)
  FPR2 <- sapply(1:nrow(bp2), function(i) sum(bp2[i,fset] != 0)/sum(fset)); AvFPR2 <- mean(FPR2)
  FPR3 <- sapply(1:nrow(bp3), function(i) sum(bp3[i,fset] != 0)/sum(fset)); AvFPR3 <- mean(FPR3)
  FPR4 <- sapply(1:nrow(bp4), function(i) sum(bp4[i,fset] != 0)/sum(fset)); AvFPR4 <- mean(FPR4)
  FPR5 <- sapply(1:nrow(bp5), function(i) sum(bp5[i,fset] != 0)/sum(fset)); AvFPR5 <- mean(FPR5)
  
  # Choosing True Model (CTM)
  # ~~~~~~~~~~~~~~~~~~~~~~~~
  # CTM <- sapply(1:nrow(bp), function(i) (sum(bp[i,tset] != 0) + sum(bp[i,fset] == 0))/(sum(tset) + sum(fset)) )
  # PCTM <- mean(CTM)
  
  res <- matrix(NA,5,14)
  colnames(res) <- c("AvMSEm0","AvMSEm1","AvABm0","AvABm1","GBICu","AvMSEp0","AvMSEp1","AvABp0","AvABp1","GBICp","AvCER","AvTPR","AvFPR","lambda")
  rownames(res) <- c("ga1.0","ga1.4","ga2.0","ga3.0","ga4.0")
  
  res[,"AvMSEm0"] <- AvMSEm0; res[,"AvMSEm1"] <- AvMSEm1; res[,"AvABm0"] <- AvABm0; res[,"AvABm1"] <- AvABm1;
  res[,"AvMSEp0"] <- c(AvMSEp0a,AvMSEp0b,AvMSEp0c,AvMSEp0d,AvMSEp0e)
  res[,"AvMSEp1"] <- c(AvMSEp1a,AvMSEp1b,AvMSEp1c,AvMSEp1d,AvMSEp1e)
  res[,"AvABp0"] <- c(AvABp0a,AvABp0b,AvABp0c,AvABp0d,AvABp0e)
  res[,"AvABp1"] <- c(AvABp1a,AvABp1b,AvABp1c,AvABp1d,AvABp1e)
  res[,"AvCER"] <- c(AvCER1,AvCER2,AvCER3,AvCER4,AvCER5)
  res[,"AvTPR"] <- c(AvTPR1,AvTPR2,AvTPR3,AvTPR4,AvTPR5)
  res[,"AvFPR"] <- c(AvFPR1,AvFPR2,AvFPR3,AvFPR4,AvFPR5)
  
  R <- X[,-c(1:(7*ix))]
  colnames(R) <- c("GBICu", "GBICp1", "GBICp2", "GBICp3","GBICp4","GBICp5","lambda1","lambda2","lambda3","lambda4","lambda5","ix")
  
  res[,"GBICu"]  <- mean(R[iu,"GBICu"])
  res[,"GBICp"]  <- sapply(1:5, function(i){ mean(R[get(paste0("ip",i)),paste0("GBICp",i)]) })
  res[,"lambda"] <- sapply(1:5, function(i){ mean(R[get(paste0("ip",i)),paste0("lambda",i)]) })
  
 return(list(res = res))
  
}


# FCOLorg[,c(101:150,153)] <- rFCOL; FCOL <- FCOLorg
# FCOLorg <- FCOL # FCOL <- FCOLorg
# FCOL1 <- FCOL[complete.cases(FCOL),]
# 
# FCOLor <- FCOL1[,c(1:mean(FCOL1[,ncol(FCOL1)-1]-1))] # Oracle
# FCOLml <- FCOL1[,c(mean(FCOL1[,ncol(FCOL1)-1]):mean(FCOL1[,ncol(FCOL1)]-1))] # Unrestricted ML
# FCOLpl <- FCOL1[,c(mean(FCOL1[,ncol(FCOL1)]):(c(mean(FCOL1[,ncol(FCOL1)])+ncol(FCOLml)-1)))] # Penalised ML
# 
# round(cbind(cb2lb(colMeans(FCOLor),lc.)$mu, lc.$mu),5)
# round(cbind(cb2lb(colMeans(FCOLml),lc.)$mu, lc.$mu),5)
# round(cbind(cb2lb(colMeans(FCOLpl),lc.)$mu, lc.$mu),5)
# 
# round(cbind(cb2lb(colMeans(FCOLor),lc.)$sigma, lc.$sigma),5)
# round(cbind(cb2lb(colMeans(FCOLml),lc.)$sigma, lc.$sigma),5)
# round(cbind(cb2lb(colMeans(FCOLpl),lc.)$sigma, lc.$sigma),5)
# 
# lcmat <- matrix(lb2cb(lc.), nrow = nrow(FCOLor), ncol = ncol(FCOLor), byrow = T)
# 
# # FCOLt <- FCOLor
# # FCOLt <- FCOLml
# # FCOLt <- FCOLpl
# 
# bias  <- FCOLt - lcmat
# rbias <- abs(bias/lcmat)*100
# rmse  <- sqrt(colMeans(bias^2))
# Bias  <- colMeans(bias)
# RBias <- colMeans(rbias)
# boxplot(bias); abline(h= 0, col = 2, lwd = 2)
# vli <- which(lb2cb(lc.) == 0)
# abline(v = vli, lwd = 1, col = "blue")
# 
# abiasOR <- colMeans(abs(FCOLor - lcmat))
# abiasML <- colMeans(abs(FCOLml - lcmat))
# abiasPL <- colMeans(abs(FCOLpl - lcmat))
# abiasRL1 <- abiasPL/abiasML
# abiasRL2 <- abiasPL/abiasOR
# abiasRL3 <- abiasML/abiasOR
# plot(abiasRL1); abline(h= 1, col = 2, lwd = 2)
# plot(abiasRL2); abline(h= 1, col = 2, lwd = 2)
# plot(abiasRL3); abline(h= 1, col = 2, lwd = 2)
# 
# rmseOR <- colMeans((FCOLor - lcmat)^2)^0.5
# rmseML <- colMeans((FCOLml - lcmat)^2)^0.5
# rmsePL <- colMeans((FCOLpl - lcmat)^2)^0.5
# rmseRL1 <- rmsePL/rmseML
# rmseRL2 <- rmsePL/rmseOR
# rmseRL3 <- rmseML/rmseOR
# plot(rmseRL1); abline(h= 1, col = 2, lwd = 2)
# plot(rmseRL2); abline(h= 1, col = 2, lwd = 2)
# plot(rmseRL3); abline(h= 1, col = 2, lwd = 2)
# 
# # boxplot(bias[,01:10]); abline(h= 0, col = 2, lwd = 2)
# # boxplot(bias[,11:20]); abline(h= 0, col = 2, lwd = 2)
# # boxplot(bias[,21:30]); abline(h= 0, col = 2, lwd = 2)
# # boxplot(bias[,31:40]); abline(h= 0, col = 2, lwd = 2)
# # boxplot(bias[,41:50]); abline(h= 0, col = 2, lwd = 2)
# 
# # length(FCOL[FCOL[,ncol(FCOL)] != 700 & FCOL[,ncol(FCOL)] != -999, ncol(FCOL)])
# # summary(FCOL[FCOL[,ncol(FCOL)] != 700 & FCOL[,ncol(FCOL)] != -999, ncol(FCOL)])
# # 
# # For horizontal table in Appendix
# # ________________________________
# Mrmse <- matrix(rmse, nrow = p, byrow = T)
# Mbias <- matrix(Bias, nrow = p, byrow = T)
# Mrbias <- matrix(RBias, nrow = p, byrow = T)
# # MR <- matrix(NA, nrow = p, ncol = ncol(Mbias)*3)
# # for(i in 1:ncol(Mrmse)){MR[,((i-1)*3+1):((i)*3)] <- round(cbind(Mrmse[,i], Mbias[,i], Mrbias[,i]),3)}
# MR <- matrix(NA, nrow = p, ncol = ncol(Mbias)*2)
# for(i in 1:ncol(Mrmse)){MR[,((i-1)*2+1):((i)*2)] <- round(cbind(Mrmse[,i], Mbias[,i]),3)}
# rownames(MR) <- paste("Item",1:p)
# colnames(MR) <- rep(c("RMSE","Bias"),ncol(lb2mb(coefs))); MR
#xtable(MR,digits = 3)
# 
# # For vertical table in Appendix
# # ______________________________
# Morg <- matrix(unlist(lc)); rownames(Morg) <- names(unlist(lc)); colnames(Morg) <- "True"
# Mest <- matrix(colMeans(FCOL1)); colnames(Mest) <- "Avg. Est."
# Mbias1 <- matrix(Bias); colnames(Mbias1) <- "Bias"
# Mrbias1 <- matrix(RBias); colnames(Mrbias1) <- "RB (%)"
# Mrmse1 <- matrix(rmse); colnames(Mrmse1) <- "RMSE"
# Msd <- matrix(apply(FCOL1,2,sd)); colnames(Msd) <- "SD"
# MR <- round(cbind(Morg, Mest, Msd, Mrmse1),3); head(MR,15)
# MR <- round(cbind(Morg, Mbias1, Mrbias1, Mrmse1),3); head(MR,15)
# #xtable(unname(MR),digits = 3)
# 
# }