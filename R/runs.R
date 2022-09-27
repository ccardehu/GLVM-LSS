
rm(list = ls())
set.seed(1234)

library(doSNOW)
source("R/prep.R")
source("R/fams.R")
source("R/glvmlss.R")
source("R/misc.R")
source("R/simC1.R")
source("R/simC2.R")

n = 500     # Number of individuals
# nsim = 1000  # Number of simulations

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Simulation 1: Heteroscedastic Normal model:
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# p = 10
# famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- Normal()}
# lc <- NULL
# lc$mu <- matrix(1,ncol = 3, nrow = p)
# lc$sigma <- matrix(1, ncol = 3, nrow = p)
# lc$mu[,1] <- runif(p,1,2)
# lc$mu[,c(2,3)] <- runif(length(lc$mu[,c(2,3)]),0.5,1.5) * sample(c(-1,1),size = 2*p,replace = T)
# lc$sigma <- matrix(runif(length(lc$sigma),0.1,0.4), nrow = p)
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # For sparse factor loading matrices:
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lc$mu[sample(2:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(lc$mu)){ if(lc$mu[i,2] != 0) lc$mu[i,3] <- rbinom(1,1,0.3)*lc$mu[i,3] }
# lc$sigma[sample(2:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(lc$sigma)){ if(lc$sigma[i,2] != 0) lc$sigma[i,3] <- rbinom(1,1,0.3)*lc$sigma[i,3] }
# # # ~~~~~~~~~~~~
# # # Simulations:
# # # ~~~~~~~~~~~~
# AA <- glvmlss_sim(n = n, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b)

# round(cbind(AA$b$mu,AB$b$mu,AC$b$mu),4)
# round(cbind(AA$b$si,AB$b$si,AC$b$si),4)
# AB$GBIC; AC$GBIC
# AB$GAIC; AC$GAIC
# AC$lambda;# AD$lambda
# AC$sse;# AD$sse

# data = AA$Y; family = famt; mu.eq = ~ Z1+Z2; sg.eq = ~ Z1+Z2;ta.eq = NULL; nu.eq = NULL;
# control <- list(EM_iter = 30, EM_use2d = T, iter.lim = 300,
#               EM_appHess = F, EM_lrate = 0.001, est.ci = T,
#               solver = "trust", start.val = NULL, mat.info = "Hessian",
#               iden.res = NULL, tol = sqrt(.Machine$double.eps), corr.lv = FALSE, tolb = 1e-4,
#               nQP = if(q == 1) 40 else { if(q == 2) ifelse(p <= 10, 15, 25) else 10 },
#               verbose = TRUE, autoL_iter = 30,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b, gamma = NULL, a = NULL, lazytrust = F)

# r1 <- NULL
# for(p in c(5,10,20)){
#   for(n in c(200,500,1000,5000)){
#     r2 <- glvmlss_parsimpost(paste0("R/DefSim/C1E1n",n,"p",p,".Rds"),
#                              out = c("PVS","MSE","AB"),
#                              mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, plot = F, outs = F)
#     r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
#   };
# }; r1; rm(r1)
 
# r1 <- NULL
# for(p in c(10)){
#   for(n in c(200,500,1000)){
#     r2 <- glvmlss_parsimpost_pml2(paste0("R/C2E1n",n,"p",p,".Rds"),
#                                  out = c("MSE","AB","PVS","lambda"),
#                                  mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2)
#     r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
#   };
# }; round(r1,4);


# # ~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~~
# # Simulation 2: IRT model:
# # ~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~~
# p = 10
# famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- Binomial()}
# lc <- NULL
# lc$mu <- matrix(1,ncol = 3, nrow = p)
# lc$mu[,1] <- runif(length(lc$mu[,1]),-1,1)
# lc$mu[,c(2:3)] <- runif(length(lc$mu[,c(2:3)]),1.5,2.5) * sample(c(-1,1),size = 2*p,replace = T)
# if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3];
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # For sparse factor loading matrices:
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lc$mu[sample(2:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(lc$mu)){ if(lc$mu[i,2] != 0) lc$mu[i,3] <- rbinom(1,1,0.5)*lc$mu[i,3] }
# # ~~~~~~~~~~~~
# # Simulations:
# # ~~~~~~~~~~~~
# AA <- glvmlss_sim(n,famt, mu.eq = ~Z1+Z2, start.val = lc)
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, verbose = T)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, verbose = T,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b)
# 
# round(cbind(AA$b$mu,AB$b$mu,AC$b$mu),5)
# GBIC(AB); GBIC(AC)
# GAIC(AB); GAIC(AC)
# AC$lambda
# AC$sse

# r1 <- NULL
# for(p in c(10,20)){
#   for(n in c(200,500,1000)){
#     r2 <- glvmlss_parsimpost_pml(paste0("R/DefSim/C2E2n",n,"p",p,".Rds"),
#                                  out = c("MSE","AB", "PVS", "lambda"),
#                                  mu.eq = ~ Z1+Z2)
#     r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
#   };
# }; round(r1,4); #rm(r1)


# # ~~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~~~
# # Simulation 3: ZI-Poisson:
# # ~~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~~~
# p = 10
# famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- ZIpoisson()}
# lc <- NULL
# lc$mu <- matrix(1,ncol = 3, nrow = p)
# lc$sigma <- matrix(1, ncol = 3, nrow = p)
# lc$mu[,1] <- runif(p,2,3)
# lc$mu[,c(2,3)] <- matrix(runif(length(lc$mu[,c(2,3)]),0.2,0.6), nrow = p) * sample(c(-1,1),size = 2*p,replace = T)
# lc$sigma[,1] <- runif(p,-2,-1)
# lc$sigma[,c(2,3)] <- runif(length(lc$sigma[,c(2,3)]),1.5,2.5)
# if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # For sparse factor loading matrices:
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lc$mu[sample(2:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(lc$mu)){ if(lc$mu[i,2] != 0) lc$mu[i,3] <- rbinom(1,1,0.5)*lc$mu[i,3] }
# lc$sigma[sample(2:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(lc$sigma)){ if(lc$sigma[i,2] != 0) lc$sigma[i,3] <- rbinom(1,1,0.5)*lc$sigma[i,3] }
# # ~~~~~~~~~~~~
# # Simulations:
# # ~~~~~~~~~~~~
# AA <- glvmlss_sim(n,famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b)
# 
# round(cbind(AA$b$mu,AB$b$mu,AC$b$mu),5)
# round(cbind(AA$b$s,AB$b$s,AC$b$s),5)
# AB$GBIC; AC$GBIC
# AB$GAIC; AC$GAIC
# AC$lambda
# AC$sse

# r1 <- NULL
# for(p in c(5,10,20)){
#   for(n in c(200,500,1000,5000)){
#     r2 <- glvmlss_parsimpost(paste0("R/DefSim/C1E2n",n,"p",p,".Rds"),
#                              out = c("MSE","AB","PVS"),
#                              mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, plot = F, outs = F)
#     r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
#   };
# }; r1; rm(r1)

# # ~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~
# # Simulation 4: Beta:
# # ~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~
# p = 10
# famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- Beta()}
# lc <- NULL
# lc$mu <- matrix(1,ncol = 2, nrow = p)
# lc$sigma <- matrix(1, ncol = 2, nrow = p)
# lc$mu[,1] <- runif(length(lc$mu[,1]), -1.5, 1.5)
# lc$mu[,2] <- runif(length(lc$mu[,1]), 0.5, 1.5) * sample(c(-1,1), p ,replace = T)
# lc$sigma[,1] <- runif(length(lc$sigma[,1]), -2, -0.5)
# lc$sigma[,2] <- runif(length(lc$sigma[,1]), 0.3, 0.6)* sample(c(-1,1), p ,replace = T)
# if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]
# if(lc$sigma[1,2] < 0) lc$sigma[1,2] <- -lc$sigma[1,2]
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # For sparse factor loading matrices:
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lc$sigma[sample(2:p, floor((p-2)/2), replace = F),2] <- 0
# # # ~~~~~~~~~~~~
# # # Simulations:
# # # ~~~~~~~~~~~~
# AA <- glvmlss_sim(n, famt, mu.eq = ~ Z1, sg.eq = ~ Z1, start.val = lc)
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, verbose = T)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, verbose = T,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b)
# 
# round(cbind(AA$b$mu,AB$b$mu,AC$b$m),4)
# round(cbind(AA$b$s,AB$b$s,AC$b$s),4)
# AB$GBIC; AC$GBIC
# AB$GAIC; AC$GAIC
# AC$lambda
# AC$sse

# testm <- function(x){  return(plogis(lc$mu[item,][1] + lc$mu[item,][2]*x)) }
# tests <- function(x){  return(plogis(lc$sigma[item,][1] + lc$sigma[item,][2]*x)) }
# 
# item = 1
# # curve(testm, from = -4, to = 4, ylim = c(0,1))
# curve(tests, from = -4, to = 4, ylim = c(0,1))
# for(ii in 2:p){
#   item = ii
#   # lines(seq(from = -4, to = 4, length.out = 100), testm(seq(from = -4, to = 4, length.out = 100)), col = ii)
#   lines(seq(from = -4, to = 4, length.out = 100), tests(seq(from = -4, to = 4, length.out = 100)), col = ii)
# }
# legend("bottomright", legend = paste0("item",1:p), lwd = 3, col = 1:p, cex = 0.8)

# data = AA$Y; family = famt; mu.eq = ~ Z1; sg.eq = ~ Z1;ta.eq = NULL; nu.eq = NULL;
# control <- list(EM_iter = 30, EM_use2d = T, iter.lim = 300,
#               EM_appHess = F, EM_lrate = 0.001, est.ci = T,
#               solver = "trust", start.val = NULL, mat.info = "Hessian",
#               iden.res = NULL, tol = sqrt(.Machine$double.eps), corr.lv = FALSE,
#               nQP = if(q == 1) 40 else { if(q == 2) ifelse(p <= 10, 15, 25) else 10 },
#               verbose = TRUE, autoL_iter = 30,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b, gamma = NULL, a = NULL)

# r1 <- NULL
# for(p in c(5,10,20)){
#   for(n in c(200,500,1000,5000)){
#     r2 <- glvmlss_parsimpost(paste0("R/DefSim/C1E3n",n,"p",p,".Rds"),
#                              out = c("PVS","MSE","AB"), trim = 0.0,
#                              mu.eq = ~ Z1, sg.eq = ~ Z1, plot = T, outs = F)
#     r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
#   };
# }; r1; rm(r1)

r1 <- NULL
for(p in c(20)){
  for(n in c(1000)){
    r2 <- glvmlss_parsimpost_pml2(paste0("R/C2E3n",n,"p",p,".Rds"),
                             out = c("PVS","MSE","AB","lambda"), trim = 0.0,
                             mu.eq = ~ Z1, sg.eq = ~ Z1)
    r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
  };
}; round(r1,4);#rm(r1)

ulim <- 0.07

plot(x = colMeans(XSeg0,na.rm = T), y = apply(Xbeg0,2,sd,na.rm = T), col = as.factor(Xb0[1,] == 0), pch = 16, ylab = "Empirical SE", xlab = "Estimated (Asymptotic) SE", main = "MLE",
     ylim = c(0,ulim), xlim = c(0,ulim))
abline(a = 0,b = 1, col = 2)

plot(x = colMeans(XSeg1,na.rm = T), y = apply(Xbeg1,2,sd,na.rm = T), col = as.factor(Xb0[1,] == 0), pch = 16, ylab = "Empirical SE", xlab = "Estimated (Asymptotic) SE", main = "Gamma = 1",
     ylim = c(0,ulim), xlim = c(0,ulim))
abline(a = 0,b = 1, col = 2)

plot(x = colMeans(XSeg2,na.rm = T), y = apply(Xbeg2,2,sd,na.rm = T), col = as.factor(Xb0[1,] == 0), pch = 16, ylab = "Empirical SE", xlab = "Estimated (Asymptotic) SE", main = "Gamma = 2",
     ylim = c(0,ulim), xlim = c(0,ulim))
abline(a = 0,b = 1, col = 2)

plot(x = colMeans(XSegf,na.rm = T), y = apply(Xbegf,2,sd,na.rm = T), col = as.factor(Xb0[1,] == 0), pch = 16, ylab = "Empirical SE", xlab = "Estimated (Asymptotic) SE", main = "Gamma = log(n)/2",
     ylim = c(0,ulim), xlim = c(0,ulim))
abline(a = 0,b = 1, col = 2)

plot(x = colMeans(XSeg3,na.rm = T), y = apply(Xbeg3,2,sd,na.rm = T), col = as.factor(Xb0[1,] == 0), pch = 16, ylab = "Empirical SE", xlab = "Estimated (Asymptotic) SE", main = "Gamma = 3",
     ylim = c(0,ulim), xlim = c(0,ulim))
abline(a = 0,b = 1, col = 2)

plot(x = colMeans(XSeg4,na.rm = T), y = apply(Xbeg4,2,sd,na.rm = T), col = as.factor(Xb0[1,] == 0), pch = 16, ylab = "Empirical SE", xlab = "Estimated (Asymptotic) SE", main = "Gamma = 4",
     ylim = c(0,ulim), xlim = c(0,ulim))
abline(a = 0,b = 1, col = 2)
 
# hist(Xbegf[,"mu7.1"], breaks = 25)
# abline(v = mean(Xbegf[,"mu7.1"], na.rm = T) + sd(Xbegf[,"mu7.1"], na.rm = T), col = 2, lty = 2, lwd = 2)
# abline(v = mean(Xbegf[,"mu7.1"], na.rm = T) - sd(Xbegf[,"mu7.1"], na.rm = T), col = 2, lty = 2, lwd = 2)
# abline(v = mean(Xbegf[,"mu7.1"], na.rm = T) + mean(XSegf[,"mu7.1"], na.rm = T), col = 3, lwd = 2, lty = 2)
# abline(v = mean(Xbegf[,"mu7.1"], na.rm = T) - mean(XSegf[,"mu7.1"], na.rm = T), col = 3, lwd = 2, lty = 2)
# 
# hist(XSegf[,"mu7.1"], breaks = 25)
