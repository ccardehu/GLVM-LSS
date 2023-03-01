rm(list = ls())
set.seed(1234)

source("R/prep.R")
source("R/fams.R")
source("R/glvmlss.R")
source("R/misc.R")
source("R/simC1.R")
source("R/simC2.R")

# n = 500     # Number of individuals
# # nsim = 1000  # Number of simulations
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
# lc$mu[,-1] <- runif(length(lc$mu[,-1]),0.5,1.5) * sample(c(-1,1),size = 2*p,replace = T)
# lc$sigma <- matrix(runif(length(lc$sigma),0.1,0.4), nrow = p)
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # For sparse factor loading matrices:
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lc$mu[sample(2:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(lc$mu)){ if(lc$mu[i,2] != 0) lc$mu[i,3] <- rbinom(1,1,0.3)*lc$mu[i,3] }
# lc$sigma[sample(2:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(lc$sigma)){ if(lc$sigma[i,2] != 0) lc$sigma[i,3] <- rbinom(1,1,0.3)*lc$sigma[i,3] }
# ~~~~~~~~~~~~
# Simulations:
# ~~~~~~~~~~~~
# AA <- glvmlss_sim(n = n, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc, iden.res = "eiv", Rz = matrix(c(1,.3,.3,1), nrow = 2))
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T, iden.res = "eiv", solver = "nlminb", corr.lv = T)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T, iden.res = "eiv", solver = "nlminb",
#               corr.lv = T, penalty = "alasso", lambda = "auto", w.alasso = AB$b)
# 
# round(cbind(AA$b$mu,AB$b$mu,AC$b$mu),3)
# round(cbind(AA$b$si,AB$b$si,AC$b$si),3)
# round(cbind(AA$Rz,AB$Rz,AC$Rz),3)
# AB$GBIC; AC$GBIC;

# data = AA$Y; family = famt; mu.eq = ~ Z1+Z2; sg.eq = ~ Z1+Z2; ta.eq = NULL; nu.eq = NULL; q = 2;
# control <- list(EM_iter = 30, EM_use2d = T, iter.lim = 300,
#             EM_appHess = F, EM_lrate = 0.01, est.ci = "Standard",
#             solver = "nlminb", start.val = NULL, mat.info = "Hessian", lazytrust = F,
#             iden.res = "eiv", tol = sqrt(.Machine$double.eps), tolb = 1e-4,
#             corr.lv = T, Rz = NULL, var.lv = c(1,1),
#             nQP = if(q == 1) 40 else { if(q == 2) ifelse(p <= 10, 15, 25) else 10 },
#             verbose = T, autoL_iter = 30, f.scores = F,
#             penalty = "alasso", lambda = "auto", w.alasso = AB$b, gamma = NULL, a = NULL)

# r1 <- NULL
# for(p in c(5,10)){
#   for(n in c(200,500)){
#     r2 <- glvmlss_parsimpost(paste0("R/C1E1n",n,"p",p,"_2023-03-02.Rds"),iden.res = "eiv",
#                              out = c("PVS","MSE","AB"),
#                              mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, plot = T, outs = F)
#     r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
#   };
# }; r1; rm(r1)
 
# r1 <- NULL
# for(p in c(10,20)){
#   for(n in c(200,500,1000)){
#     r2 <- glvmlss_parsimpost_pml(paste0("R/C2E1n",n,"p",p,"_HetLFM.Rds"),
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
# AA <- glvmlss_sim(n,famt, mu.eq = ~Z1+Z2, start.val = lc, Rz = matrix(c(1,.45,.45,1),nrow = 2), iden.res = "eiv")
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, verbose = T)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, verbose = T,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b)
# AD <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, verbose = T, corr.lv = T, iden.res = "eiv")
# AE <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, verbose = T, iden.res = "eiv",
#               penalty = "alasso", lambda = "auto", w.alasso = AD$b, corr.lv = T)
# 
# round(cbind(AA$b$mu,AB$b$mu,AC$b$mu, AD$b$mu,AE$b$mu),3) #,AB$b$mu,AC$b$mu,
# round(cbind(AA$Rz,AB$Rz,AC$Rz,AD$Rz,AE$Rz),3) #AB$Rz,AC$Rz,
# AB$GBIC; AC$GBIC; AD$GBIC; AE$GBIC
# 
# round(cbind(AA$b$mu,AB$b$mu,AC$b$mu),5)
# GBIC(AB); GBIC(AC)
# GAIC(AB); GAIC(AC)
# AC$lambda
# AC$sse
# 
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
# AA <- glvmlss_sim(n,famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc, Rz = matrix(c(1,.45,.45,1),nrow = 2), iden.res = "eiv")
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b)
# AD <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T, corr.lv = T, iden.res = "eiv")
# AE <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T, iden.res = "eiv",
#               penalty = "alasso", lambda = "auto", w.alasso = AD$b, corr.lv = T)
# 
# round(cbind(AA$b$mu,AB$b$mu,AC$b$mu, AD$b$mu,AE$b$mu),3) #,AB$b$mu,AC$b$mu,
# round(cbind(AA$b$si,AB$b$si,AC$b$si,AD$b$si,AE$b$si),3) #,AB$b$si,AC$b$si
# round(cbind(AA$Rz,AB$Rz,AC$Rz,AD$Rz,AE$Rz),3) #AB$Rz,AC$Rz,
# AB$GBIC; AC$GBIC; AD$GBIC; AE$GBIC
# 
# round(cbind(AA$b$mu,AB$b$mu,AC$b$mu),5)
# round(cbind(AA$b$s,AB$b$s,AC$b$s),5)
# AB$GBIC; AC$GBIC
# AB$GAIC; AC$GAIC
# AC$lambda
# AC$sse
# 
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
# n = 500
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
# A0 <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ 1, verbose = T, f.scores = T)
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, verbose = T, f.scores = T)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, verbose = T, f.scores = T,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b)
#   
# round(cbind(AA$b$mu,AB$b$mu,AC$b$m),4)
# round(cbind(AA$b$s,AB$b$s,AC$b$s),4)
# 
# round(cbind(AB$SE$mu,AC$S$m),4)
# round(cbind(AB$SE$s,AC$S$s),4)
# 
# A0$GBIC; AB$GBIC; AC$GBIC
# A0$GAIC; AB$GAIC; AC$GAIC
# AC$lambda
# AC$sse
# 
# hist(A0$f.scores$Z1, breaks = 25, freq = F); lines(density(c(A0$f.scores$Z1)), col = 2)
# hist(AB$f.scores$Z1, breaks = 25, freq = F); lines(density(c(AB$f.scores$Z1)), col = 2)
# hist(AC$f.scores$Z1, breaks = 25, freq = F); lines(density(c(AC$f.scores$Z1)), col = 2)
# 
# confintr::ci_cor(AA$Z$Z1,A0$f.scores$Z1, method = "kendall", type = c("bootstrap"))
# confintr::ci_cor(AA$Z$Z1,AB$f.scores$Z1, method = "kendall", type = c("bootstrap"))
# confintr::ci_cor(AA$Z$Z1,AC$f.scores$Z1, method = "kendall", type = c("bootstrap"))
# 
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
# 
# data = AA$Y; family = famt; mu.eq = ~ Z1; sg.eq = ~ Z1;ta.eq = NULL; nu.eq = NULL;
# control <- list(EM_iter = 30, EM_use2d = T, iter.lim = 300,
#               EM_appHess = F, EM_lrate = 0.001, est.ci = T,
#               solver = "trust", start.val = NULL, mat.info = "Hessian",
#               iden.res = NULL, tol = sqrt(.Machine$double.eps), corr.lv = FALSE,
#               nQP = if(q == 1) 40 else { if(q == 2) ifelse(p <= 10, 15, 25) else 10 },
#               verbose = TRUE, autoL_iter = 30,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b, gamma = NULL, a = NULL)
# 
# r1 <- NULL
# for(p in c(10,20)){
#   for(n in c(200,500,1000,5000)){
#     r2 <- glvmlss_parsimpost(paste0("R/C1E3n",n,"p",p,"_HetBFM.Rds"),
#                              out = c("PVS","MSE","AB","iter"), trim = 0.0,
#                              mu.eq = ~ Z1, sg.eq = ~ Z1, plot = T, outs = F)
#     r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
#   };
# }; round(r1,4); #rm(r1)
#  
# r1 <- NULL
# for(p in c(10,20)){
#   for(n in c(200,500,1000)){
#     r2 <- glvmlss_parsimpost_pml(paste0("R/C2E3n",n,"p",p,"_HetBFM.Rds"),
#                              out = c("MSE","AB","lambda", "PVS"),
#                              mu.eq = ~ Z1, sg.eq = ~ Z1)
#     r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
#   };
# }; round(r1,4);#rm(r1)


# ulim <- 0.07
# plot(x = colMeans(XSeg0,na.rm = T), y = apply(Xbeg0,2,sd,na.rm = T), col = as.factor(Xb0[1,] == 0), pch = 16, ylab = "Empirical SE", xlab = "Estimated (Asymptotic) SE", main = "MLE",
#      ylim = c(0,ulim), xlim = c(0,ulim))
# abline(a = 0,b = 1, col = 2)
# plot(x = colMeans(XSeg1,na.rm = T), y = apply(Xbeg1,2,sd,na.rm = T), col = as.factor(Xb0[1,] == 0), pch = 16, ylab = "Empirical SE", xlab = "Estimated (Asymptotic) SE", main = "Gamma = 1",
#      ylim = c(0,ulim), xlim = c(0,ulim))
# abline(a = 0,b = 1, col = 2)
# plot(x = colMeans(XSeg2,na.rm = T), y = apply(Xbeg2,2,sd,na.rm = T), col = as.factor(Xb0[1,] == 0), pch = 16, ylab = "Empirical SE", xlab = "Estimated (Asymptotic) SE", main = "Gamma = 2",
#      ylim = c(0,ulim), xlim = c(0,ulim))
# abline(a = 0,b = 1, col = 2)
# plot(x = colMeans(XSeg3,na.rm = T), y = apply(Xbeg3,2,sd,na.rm = T), col = as.factor(Xb0[1,] == 0), pch = 16, ylab = "Empirical SE", xlab = "Estimated (Asymptotic) SE", main = "Gamma = 3",
#      ylim = c(0,ulim), xlim = c(0,ulim))
# abline(a = 0,b = 1, col = 2)
# plot(x = colMeans(XSeg4,na.rm = T), y = apply(Xbeg4,2,sd,na.rm = T), col = as.factor(Xb0[1,] == 0), pch = 16, ylab = "Empirical SE", xlab = "Estimated (Asymptotic) SE", main = "Gamma = 4",
#      ylim = c(0,ulim), xlim = c(0,ulim))
# abline(a = 0,b = 1, col = 2)
# plot(x = colMeans(XSeg5,na.rm = T), y = apply(Xbeg5,2,sd,na.rm = T), col = as.factor(Xb0[1,] == 0), pch = 16, ylab = "Empirical SE", xlab = "Estimated (Asymptotic) SE", main = "Gamma = 5",
#      ylim = c(0,ulim), xlim = c(0,ulim))
# abline(a = 0,b = 1, col = 2)

# # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # # Simulation 5: Skewed-Normal distribution:
# # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rm(list = ls())
# set.seed(1234)
# 
# source("R/prep.R")
# source("R/fams.R")
# source("R/glvmlss.R")
# source("R/misc.R")
# source("R/simC1.R")
# source("R/simC2.R")
# 
# n = 2000     # Number of individuals
# 
# p = 10
# famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- SkewNormal()}
# lc <- NULL
# lc$mu <- matrix(1,ncol = 2, nrow = p)
# lc$sigma <- matrix(1, ncol = 2, nrow = p)
# lc$nu <- matrix(1, ncol = 2, nrow = p)
# lc$mu[,1] <- runif(p,-1,1)
# lc$mu[,c(2)] <- runif(length(lc$mu[,c(2)]),0.5,1.5) * sample(c(-1,1),size = p,replace = T)
# lc$sigma <- matrix(runif(length(lc$sigma),0.2,0.4), nrow = p)
# lc$nu[,1] <- runif(p,-2,2)
# lc$nu[,2] <- matrix(runif(length(lc$nu[,2]),0.2,0.5), nrow = p) * sample(c(-1,1),size = p,replace = T)
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # For sparse factor loading matrices:
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lc$sigma[sample(1:p, floor(p/2), replace = F),2] <- 0
# lc$nu[sample(1:p, floor(p/2), replace = F),2] <- 0
# 
# AA <- glvmlss_sim(n = n, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, nu.eq = ~ Z1, start.val = lc)
# A1 <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, nu.eq = ~ Z1, verbose = T,
#               solver = "nlminb", iter.lim = 1000, EM_use2d = F, est.ci = "Approximate", EM_iter = 1000)
# A1a <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, nu.eq = ~ Z1, verbose = T,
#               solver = "nlminb", iter.lim = 1000, EM_use2d = F, est.ci = "Approximate", EM_iter = 1000,
#               penalty = "alasso", lambda = "auto", w.alasso = A1$b)
# 
# round(cbind(AA$b$mu, A1$b$mu, A1a$b$mu), 3)
# round(cbind(AA$b$si, A1$b$si, A1a$b$si), 3)
# round(cbind(AA$b$nu, A1$b$nu, A1a$b$nu), 3)
# round(A1a$lambda,5)
# round(cbind(A1$SE$b$mu, A1a$SE$b$mu),3)
# round(cbind(A1$SE$b$si, A1a$SE$b$si),3)
# round(cbind(A1$SE$b$nu, A1a$SE$b$nu),3)
# 
# r1 <- NULL
# for(p in c(5,10,20)){ # 
#   for(n in c(200,500,1000,5000)){ #,5000
#     r2 <- glvmlss_parsimpost(file = paste0("R/C1E4(CP-AppSE)n",n,"p",p,"_2023-02-16.Rds"),
#                              out = c("MSE","AB","PVS"),
#                              mu.eq = ~ Z1, sg.eq = ~ Z1, nu.eq = ~ Z1, plot = T, outs = T) # ta.eq = NULL
#     r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
#   };
# }; round(r1,4); #rm(r1)
# 
# r1[,"PVS"] <- r1[,"PVS"]*100
# xtable::xtable(r1, digits = 4)
# 
# r1[,1] <- r1[,1]*100
# 
# rownames(r1) <- c(paste0(paste0("n",c(500,1000,5000)),"p",10), paste0(paste0("n",c(500,1000,5000)),"p",20))
# 
# temp <- data.frame(x = apply(Xbe, MARGIN = 2, FUN = sd),
#                    y = colMeans(Xse),
#                    nu = as.factor(grepl("nu",colnames(Xse),fixed = T)))
# 
# plot(x = temp$x, y = temp$y, col = temp$nu, pch = 16,
#      xlab = "Empirical SE", ylab = "Estimated SE (outer gradient)",
#      main = "Empirical vs estimated SE (n = 200, p = 10)")
# abline(a = 0, b = 1, col = 4, lwd = 2)
#  
# 
# data = AA$Y; family = famt; mu.eq = ~ Z1; sg.eq = ~ Z1; nu.eq = ~Z1; ta.eq = NULL; q = 1;
# control <- list(EM_iter = 1000, EM_use2d = F, iter.lim = 300,
#             EM_appHess = T, EM_lrate = 0.001, est.ci = "Approximate",
#             solver = "nlminb", start.val = NULL, mat.info = "Hessian", lazytrust = F,
#             iden.res = NULL, tol = sqrt(.Machine$double.eps), tolb = 1e-4,
#             corr.lv = FALSE, Rz = NULL,
#             nQP = 40,
#             verbose = T, autoL_iter = 30, f.scores = F,
#             penalty = "alasso", lambda = "auto", w.alasso = A1$b, gamma = NULL, a = NULL)
# 
# testf <- function(arg, Y, ghQ, borg, famL){
#   b <- cb2lb(arg,borg)
#   return(-fyz(Y,ghQ,b,famL)$ll)
# }
# 
# numGrad <- numDeriv::grad(func = testf, x = lb2cb(lc),Y = Y, ghQ = ghQ, borg = lc, famL = famL)
# anaGrad <- -d1ll(Y,ghQ,lc,famL,info,fyz(Y,ghQ,lc,famL)$pD,rb)
# 
# plot.ts(numGrad, type = "o")
# lines(anaGrad, col = 2)
# 
# numHess <- numDeriv::hessian(func = testf, x = lb2cb(lc),Y = Y, ghQ = ghQ, borg = lc, famL = famL)
# anaHess <- -d2ll(Y,ghQ,lc,famL,"Hessian",fyz(Y,ghQ,lc,famL)$pD,rb)
# # appHess <- ad2ll(Y,ghQ,lc,famL,"Hessian",fyz(Y,ghQ,lc,famL)$pD, rb)
# 
# plot.ts(diag(numHess) - diag(anaHess), type = "o")
# # plot.ts(diag(numHess) - diag(appHess), type = "o")
# # plot.ts(diag(anaHess) - diag(appHess), type = "o")
# lines(diag(anaHess), col = 2)
