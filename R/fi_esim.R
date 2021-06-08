
library(foreach)
library(doSNOW)
library(parallel)
library(itertools)
library(xtable)

sil = T; form = s.form; restr = l1.; coefs = lc.; saveRes = T; cleanRes = F
coefsf = lc; restrf = l1

splvm.simfit <- function(l,sil = sil){
 AA <- splvm.sim(n = n, form = form, fam = fam, constraints = restr, coefs = coefs)
 if(ncol(AA$Z) == 1) nqp <- 40 else nqp <- 20
 # tmp0 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, 
 #                   control = list(method = "EM", start.val = coefs, constraint = restr,
 #                   ghQqp = nqp, iter.lim = 1e3, full.hess = F, tol = sqrt(.Machine$double.eps),
 #                   silent = T, information = "Fisher"))
 # tmp1 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form,
 #                   control = list(method = "EM", start.val = coefsf, constraint = restrf,
 #                   ghQqp = nqp, iter.lim = 1e3, full.hess = F, tol = sqrt(.Machine$double.eps),
 #                   silent = T, information = "Fisher"))
 tmp2 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form,
                   control = list(method = "PEM", 
                   ghQqp = nqp, iter.lim = 1e3, full.hess = F, tol = sqrt(.Machine$double.eps),
                   silent = T, information = "Fisher",
                   pml.control = list(type = "mcp", lambda = 0.0605, a = 1.95, w.alasso = NULL, pen.load = F)))
 if(!sil) cat("Run = ", l, "Iter. to Convg.: ",tmp1$iter,"\n")
 # return(c(lb2cb(tmp0$b),lb2cb(tmp1$b),lb2cb(tmp2$b),tmp0$iter,tmp1$iter,tmp2$iter))
 return(c(lb2cb(tmp2$b),tmp2$iter))
}

cores <- detectCores()-2
cl<-makeCluster(cores)
registerDoSNOW(cl)
progress <- function(n) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), n)
opts <- list(progress = progress)

FCOL <- foreach(l = 1:nsim,
                .combine = rbind,
                .packages = c("mvtnorm", "fastGHQuad","gamlss"),
                .options.snow = opts) %dopar% splvm.simfit(l,sil = sil)

stopCluster(cl)
if(cleanRes) FCOL <- FCOL[FCOL[,ncol(FCOL)] %!in% c(1000,-999),-ncol(FCOL)]
if(saveRes) save(FCOL, file = paste0("nsim",nsim,"_n",n, "_p", p,"_Ex1.RData"))

# FCOLorg[,c(101:150,153)] <- rFCOL; FCOL <- FCOLorg
FCOLorg <- FCOL # FCOL <- FCOLorg
FCOL1 <- FCOL[!apply(FCOL,1,function(a) any(a %in% c(1000,-999))),]

FCOLor <- FCOL1[,c(1:40,121)] # Oracle
FCOLml <- FCOL1[,c(41:80,122)] # Unrestricted ML
FCOLpl <- FCOL1[,c(81:120,123)] # Penalised ML

round(cbind(cb2lb(colMeans(FCOLor[,-ncol(FCOLor)]),lc.)$mu, lc.$mu),5)
round(cbind(cb2lb(colMeans(FCOLml[,-ncol(FCOLml)]),lc.)$mu, lc.$mu),5)
round(cbind(cb2lb(colMeans(FCOLpl[,-ncol(FCOLpl)]),lc.)$mu, lc.$mu),5)

round(cbind(cb2lb(colMeans(FCOLor[,-ncol(FCOLor)]),lc.)$sigma, lc.$sigma),5)
round(cbind(cb2lb(colMeans(FCOLml[,-ncol(FCOLor)]),lc.)$sigma, lc.$sigma),5)
round(cbind(cb2lb(colMeans(FCOLpl[,-ncol(FCOLor)]),lc.)$sigma, lc.$sigma),5)
 
lcmat <- matrix(lb2cb(lc.), nrow = nrow(FCOLor), ncol = ncol(FCOLor)-1, byrow = T)

# FCOLt <- FCOLor[,-ncol(FCOLor)]
# FCOLt <- FCOLml[,-ncol(FCOLml)]
# FCOLt <- FCOLpl[,-ncol(FCOLpl)]

bias  <- FCOLt - lcmat
rbias <- abs(bias/lcmat)*100
rmse  <- sqrt(colMeans(bias^2))
Bias  <- colMeans(bias)
RBias <- colMeans(rbias)
boxplot(bias); abline(h= 0, col = 2, lwd = 2)

abiasOR <- colMeans(abs(FCOLor[,-ncol(FCOLor)] - lcmat))
abiasML <- colMeans(abs(FCOLml[,-ncol(FCOLml)] - lcmat))
abiasPL <- colMeans(abs(FCOLpl[,-ncol(FCOLpl)] - lcmat))
abiasRL1 <- abiasPL/abiasML
abiasRL2 <- abiasPL/abiasOR
abiasRL3 <- abiasML/abiasOR
plot(abiasRL1); abline(h= 1, col = 2, lwd = 2)
plot(abiasRL2); abline(h= 1, col = 2, lwd = 2)
plot(abiasRL3); abline(h= 1, col = 2, lwd = 2)

rmseOR <- colMeans((FCOLor[,-ncol(FCOLor)] - lcmat)^2)^0.5
rmseML <- colMeans((FCOLml[,-ncol(FCOLml)] - lcmat)^2)^0.5
rmsePL <- colMeans((FCOLpl[,-ncol(FCOLpl)] - lcmat)^2)^0.5
rmseRL1 <- rmsePL/rmseML
rmseRL2 <- rmsePL/rmseOR
rmseRL3 <- rmseML/rmseOR
plot(rmseRL1); abline(h= 1, col = 2, lwd = 2)
plot(rmseRL2); abline(h= 1, col = 2, lwd = 2)
plot(rmseRL3); abline(h= 1, col = 2, lwd = 2)

# boxplot(bias[,01:10]); abline(h= 0, col = 2, lwd = 2)
# boxplot(bias[,11:20]); abline(h= 0, col = 2, lwd = 2)
# boxplot(bias[,21:30]); abline(h= 0, col = 2, lwd = 2)
# boxplot(bias[,31:40]); abline(h= 0, col = 2, lwd = 2)
# boxplot(bias[,41:50]); abline(h= 0, col = 2, lwd = 2)

# length(FCOL[FCOL[,ncol(FCOL)] != 700 & FCOL[,ncol(FCOL)] != -999, ncol(FCOL)])
# summary(FCOL[FCOL[,ncol(FCOL)] != 700 & FCOL[,ncol(FCOL)] != -999, ncol(FCOL)])
# 
# For horizontal table in Appendix
# ________________________________
Mrmse <- matrix(rmse, nrow = p, byrow = T)
Mbias <- matrix(Bias, nrow = p, byrow = T)
Mrbias <- matrix(RBias, nrow = p, byrow = T)
# MR <- matrix(NA, nrow = p, ncol = ncol(Mbias)*3)
# for(i in 1:ncol(Mrmse)){MR[,((i-1)*3+1):((i)*3)] <- round(cbind(Mrmse[,i], Mbias[,i], Mrbias[,i]),3)}
MR <- matrix(NA, nrow = p, ncol = ncol(Mbias)*2)
for(i in 1:ncol(Mrmse)){MR[,((i-1)*2+1):((i)*2)] <- round(cbind(Mrmse[,i], Mbias[,i]),3)}
rownames(MR) <- paste("Item",1:p)
colnames(MR) <- rep(c("RMSE","Bias"),ncol(lb2mb(coefs))); MR
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