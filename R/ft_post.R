
# Code
# File ft_post.R: Post-processing for Simulation examples (Chapter 2)
# Code written by: Camilo Cardenas-Hurtado (c.a.cardenas-hurtado@lse.ac.uk)


# FCOLorg[,c(101:150,153)] <- rFCOL; FCOL <- FCOLorg
FCOLorg <- FCOL # FCOL <- FCOLorg
FCOL1 <- FCOL[complete.cases(FCOL),]

FCOLor <- FCOL1[,c(1:mean(FCOL1[,ncol(FCOL1)-1]-1))] # Oracle
FCOLml <- FCOL1[,c(mean(FCOL1[,ncol(FCOL1)-1]):mean(FCOL1[,ncol(FCOL1)]-1))] # Unrestricted ML
FCOLpl <- FCOL1[,c(mean(FCOL1[,ncol(FCOL1)]):(c(mean(FCOL1[,ncol(FCOL1)])+ncol(FCOLml)-1)))] # Penalised ML

round(cbind(cb2lb(colMeans(FCOLor),lc.)$mu, lc.$mu),5)
round(cbind(cb2lb(colMeans(FCOLml),lc.)$mu, lc.$mu),5)
round(cbind(cb2lb(colMeans(FCOLpl),lc.)$mu, lc.$mu),5)

round(cbind(cb2lb(colMeans(FCOLor),lc.)$sigma, lc.$sigma),5)
round(cbind(cb2lb(colMeans(FCOLml),lc.)$sigma, lc.$sigma),5)
round(cbind(cb2lb(colMeans(FCOLpl),lc.)$sigma, lc.$sigma),5)

lcmat <- matrix(lb2cb(lc.), nrow = nrow(FCOLor), ncol = ncol(FCOLor), byrow = T)

# FCOLt <- FCOLor
# FCOLt <- FCOLml
# FCOLt <- FCOLpl

bias  <- FCOLt - lcmat
rbias <- abs(bias/lcmat)*100
rmse  <- sqrt(colMeans(bias^2))
Bias  <- colMeans(bias)
RBias <- colMeans(rbias)
boxplot(bias); abline(h= 0, col = 2, lwd = 2)
vli <- which(lb2cb(lc.) == 0)
abline(v = vli, lwd = 1, col = "blue")

abiasOR <- colMeans(abs(FCOLor - lcmat))
abiasML <- colMeans(abs(FCOLml - lcmat))
abiasPL <- colMeans(abs(FCOLpl - lcmat))
abiasRL1 <- abiasPL/abiasML
abiasRL2 <- abiasPL/abiasOR
abiasRL3 <- abiasML/abiasOR
plot(abiasRL1); abline(h= 1, col = 2, lwd = 2)
plot(abiasRL2); abline(h= 1, col = 2, lwd = 2)
plot(abiasRL3); abline(h= 1, col = 2, lwd = 2)

rmseOR <- colMeans((FCOLor - lcmat)^2)^0.5
rmseML <- colMeans((FCOLml - lcmat)^2)^0.5
rmsePL <- colMeans((FCOLpl - lcmat)^2)^0.5
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