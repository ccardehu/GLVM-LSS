
library(foreach)
library(doSNOW)
library(parallel)
library(itertools)
library(xtable)

EM_GLVM <- function(Silent = T, coefs = lc, nsim = nsim, n = n, p = p, form = form, fam = fam, loadmt = l1, ghp = 50,
                    saveRes = T, cleanRes = F){

EMfun <- function(l, silent = Silent){ # CHANGE HERE!
  
  AA <- simGLVM(n = n, p = p, form = form, dist = fam, coefs = coefs, loadmt = loadmt)
  tmp1 <- GLVM.fit(Y = AA$Y, fam = fam, form = AA$form, silent = Silent, ghp = ghp, iter.lim = 700, tol = 1e-7, loadmt = loadmt, icoefs = coefs, useoptim = F, skipEM = F)
  if(silent == F) cat("Run = ", l, "Iter. to Convg.: ",tmp1$iter,"\n")
  return(c(unlist(tmp1$b),tmp1$iter))
  
}

cores <- detectCores()-2
cl<-makeCluster(cores)
registerDoSNOW(cl)
progress <- function(n) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), n)
opts <- list(progress = progress)

FCOL <- foreach(l = 1:nsim,
                .combine = rbind,
                .packages = c("mvtnorm", "fastGHQuad"),
                .options.snow = opts) %dopar% EMfun(l, silent = T)

stopCluster(cl)

if(cleanRes == T) FCOL <- FCOL[FCOL[,ncol(FCOL)] %!in% c(700,-999),-ncol(FCOL)]
if(saveRes == T) save(FCOL, file = paste0("nsim",nsim,"_n",n, "_p", p,"_ZIP(LinMu_QtSig).RData"))
return(FCOL)

}

# # FCOLorg <- FCOL # FCOL <- FCOLorg
# FCOL1 <- FCOL[FCOL[,ncol(FCOL)] %!in% c(700,-999),-ncol(FCOL)]
# #FCOLHT <- FCOL1[,1:40]
# #FCOLHM <- FCOL1[,51:90]
# lc1 <- NULL
# lc1$mu <- lc$mu*l1$mu
# lc1$sigma <- lc$sigma*l1$sigma
# 
# coefmod(colMeans(FCOL1),borg,gr = ex1$gr,loadmt = l1)$mu; borg$mu
# coefmod(colMeans(FCOL1),borg,gr = ex1$gr,loadmt = l1)$sigma; borg$sigma
# 
# lcmat <- matrix(unlist(lc), nrow = nrow(FCOL1), ncol = ncol(FCOL1), byrow = T)
# 
# bias  <- FCOL1 - lcmat
# rbias <- (bias/lcmat)*100
# rmse  <- sqrt(colMeans(bias^2))
# Bias  <- colMeans(bias)
# RBias <- colMeans(rbias)
# 
# boxplot(bias[,01:10]); abline(h= 0, col = 2, lwd = 2)
# boxplot(bias[,11:20]); abline(h= 0, col = 2, lwd = 2)
# boxplot(bias[,21:30]); abline(h= 0, col = 2, lwd = 2)
# boxplot(bias[,31:40]); abline(h= 0, col = 2, lwd = 2)
# boxplot(bias[,41:50]); abline(h= 0, col = 2, lwd = 2)
# boxplot(bias[,51:60]); abline(h= 0, col = 2, lwd = 2)
# boxplot(bias[,61:70]); abline(h= 0, col = 2, lwd = 2)
# boxplot(bias[,71:80]); abline(h= 0, col = 2, lwd = 2)
# 
# 
# length(FCOL[FCOL[,ncol(FCOL)] != 700 & FCOL[,ncol(FCOL)] != -999, ncol(FCOL)])
# summary(FCOL[FCOL[,ncol(FCOL)] != 700 & FCOL[,ncol(FCOL)] != -999, ncol(FCOL)])
# 
# # For horizontal table in Appendix
# # ________________________________
# Mrmse <- matrix(rmse, nrow = p, byrow = F)
# Mbias <- matrix(Bias, nrow = p, byrow = F)
# Mrbias <- matrix(RBias, nrow = p, byrow = F)
# MR <- matrix(NA, nrow = p, ncol = ncol(Mbias)*3)
# for(i in 1:ncol(Mrmse)){MR[,((i-1)*3+1):((i)*3)] <- round(cbind(Mrmse[,i], Mbias[,i], Mrbias[,i]),3)}
# rownames(MR) <- paste("Item",1:p); MR
# #xtable(MR,digits = 3)
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