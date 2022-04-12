
rm(list = ls())
set.seed(1234)

library(doSNOW)
source("R/prep.R")
source("R/fams.R")
source("R/glvmlss.R")
source("R/misc.R")

n = 500     # Number of individuals
# nsim = 1000  # Number of simulations

# # Simulation 1: Heteroscedastic Normal model:
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p = 10
famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- Normal()}
lc <- NULL
lc$mu <- matrix(1,ncol = 3, nrow = p)
lc$sigma <- matrix(1, ncol = 3, nrow = p)
lc$mu[,1] <- runif(p,1,2)
lc$mu[,c(2,3)] <- runif(length(lc$mu[,c(2,3)]),0.1,1.5) * sample(c(-1,1),size = 2*p,replace = T)
lc$sigma <- matrix(runif(length(lc$sigma),0.1,0.4), nrow = p)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For sparse factor loading matrices:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lc$mu[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
for(i in 3:nrow(lc$mu)){ if(lc$mu[i,2] != 0) lc$mu[i,3] <- rbinom(1,1,0.5) }
lc$sigma[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
for(i in 3:nrow(lc$sigma)){ if(lc$sigma[i,2] != 0) lc$sigma[i,3] <- rbinom(1,1,0.5) }
# ~~~~~~~~~~~~
# Simulations:
# ~~~~~~~~~~~~
AA <- glvmlss_sim(n,famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T, solver = "trust")
AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T,
              penalty = "alasso", lambda = "auto", w.alasso = AB$b, solver = "trust")
round(cbind(AA$b$mu,AB$b$mu,AC$b$mu),5)
round(cbind(AA$b$si,AB$b$si,AC$b$si),5)
GBIC(AB); GBIC(AC);
GAIC(AB); GAIC(AC);
# ~~~~~~~~~~~~~~
# Application 1:
# ~~~~~~~~~~~~~~
# be <- NULL
# be$mu <- matrix(1,ncol = 2, nrow = p)
# be$sigma <- matrix(1, ncol = 2, nrow = p)
# be$mu[,1] <- runif(p,1,2)
# be$mu[,2] <- runif(length(be$mu[,2]),0.1,1.5) * sample(c(-1,1),size = p,replace = T)
# be$sigma <- matrix(runif(length(be$sigma),0.1,0.4), nrow = p)
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # For sparse factor loading matrices:
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# be$sigma[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
# BA <- glvmlss_sim(n,famt, mu.eq = ~ Z1, sg.eq = ~ Z1, start.val = be)
# B1 <- glvmlss(data = BA$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ 1, verbose = T)
# B2 <- glvmlss(data = BA$Y, family = famt, mu.eq = ~ Z1 + Z2, sg.eq = ~ 1, verbose = T)
# B3 <- glvmlss(data = BA$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, verbose = T)
# B4 <- glvmlss(data = BA$Y, family = famt, mu.eq = ~ Z1, sg.eq = ~ Z1, verbose = T,
#               penalty = "alasso", lambda = "auto", w.alasso = B3$b, solver = "trust")
# B5 <- glvmlss(data = BA$Y, family = famt, mu.eq = ~ Z1 + Z2, sg.eq = ~ Z1 + Z2, verbose = T)
# B6 <- glvmlss(data = BA$Y, family = famt, mu.eq = ~ Z1 + Z2, sg.eq = ~ Z1 + Z2, verbose = T,
#               penalty = "alasso", lambda = "auto", w.alasso = B5$b, solver = "trust")
# GBIC(B1); GBIC(B2); GBIC(B3); GBIC(B4); GBIC(B5); GBIC(B6)


# famA <- vector("list",ncol(ais[,-c(1,2)])); for(i in 1:ncol(ais[,-c(1,2)])){ famA[[i]] <- Normal()}
# testA1 <- glvmlss(data = ais[,-c(1,2)], family = famA, mu.eq = ~ Z1, sg.eq = ~ 1, verbose = T)
# testA2 <- glvmlss(data = ais[,-c(1,2)], family = famA, mu.eq = ~ Z1 + Z2, sg.eq = ~ 1, verbose = T)
# testB1 <- glvmlss(data = ais[,-c(1,2)], family = famA, mu.eq = ~ Z1, sg.eq = ~ Z1, verbose = T)
# testB2 <- glvmlss(data = ais[,-c(1,2)], family = famA, mu.eq = ~ Z1 + Z2, sg.eq = ~ Z1, verbose = T)
# testC1 <- glvmlss(data = ais[,-c(1,2)], family = famA, mu.eq = ~ Z1 + Z2, sg.eq = ~ Z1 + Z2, verbose = T)
# testC2 <- glvmlss(data = ais[,-c(1,2)], family = famA, mu.eq = ~ Z1 + Z2, sg.eq = ~ Z1 + Z2, verbose = T,
#                   penalty = "alasso", lambda = "auto", w.alasso = testC1$b)
# testD1 <- glvmlss(data = ais[,-c(1,2)], family = famA, mu.eq = ~ Z1 + Z2 + Z3, sg.eq = ~ Z1 + Z2 + Z3, verbose = T)
# testD2 <- glvmlss(data = ais[,-c(1,2)], family = famA, mu.eq = ~ Z1 + Z2 + Z3, sg.eq = ~ Z1 + Z2 + Z3, verbose = T,
#                   penalty = "alasso", lambda = "auto", w.alasso = testD1$b)
# GBIC(testA1); GBIC(testA2); GBIC(testB1); GBIC(testB2)
# GBIC(testC1); GBIC(testC2); GBIC(testD1); GBIC(testD2)
# round(cbind(testC1$b$mu, testD1$b$mu), 4)
# round(cbind(testC2$b$mu, testD2$b$mu), 4)


# data = AA$Y; family = famt; mu.eq = ~ Z1+Z2; sg.eq = ~ Z1+Z2;ta.eq = NULL; nu.eq = NULL;
# control <- list(EM_iter = 30, EM_use2d = T, iter.lim = 300,
#               EM_appHess = F, EM_lrate = 0.001, est.ci = T,
#               solver = "trust", start.val = NULL, mat.info = "Hessian",
#               iden.res = NULL, tol = sqrt(.Machine$double.eps), corr.lv = FALSE,
#               nQP = if(q == 1) 40 else { if(q == 2) ifelse(p <= 10, 15, 25) else 10 },
#               verbose = TRUE, autoL_iter = 30,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b, gamma = NULL, a = NULL)

# # FCOL <- glvmlss_parsimE1(nsim, T)
# r1 <- NULL
# for(p in c(10,20)){
# for(n in c(200,500,1000,5000)){
#   r2 <- glvmlss_parsimpost(readRDS(paste0("R/E1n",n,"p",p,".Rds")), p,
#                            out = c("MSE","AB","PVS"),
#                            mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2)
#   r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
# }; #r1; #r1 <- NULL
# }; r1; rm(r1)


# Simulation 2: IRT model:
# ~~~~~~~~~~~~~~~~~~~~~~~~
# p = 10
# famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- Binomial()}
# lc <- NULL
# lc$mu <- matrix(1,ncol = 3, nrow = p)
# lc$mu[,1] <- runif(length(lc$mu[,1]),-1,1)
# lc$mu[,c(2:3)] <- runif(length(lc$mu[,c(2:3)]),1.5,2.5) * sample(c(-1,1),size = 2*p,replace = T)
# # lc$mu[,4] <- runif(length(lc$mu[,c(4)]),-0.5,0.5) * sample(c(-1,1),size = p,replace = T)
# if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3];
# # irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0), c("mu",3,"Z1:Z2",0), c("mu",1,"Z1:Z2",0), c("mu",2,"Z1:Z2",0))
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # For sparse factor loading matrices:
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lc$mu[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(lc$mu)){ if(lc$mu[i,2] != 0) lc$mu[i,3] <- rbinom(1,1,0.5) }
# # ~~~~~~~~~~~~
# # Simulations:
# # ~~~~~~~~~~~~
# AA <- glvmlss_sim(n,famt, mu.eq = ~Z1+Z2, start.val = lc)
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, verbose = T, est.ci = T)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, verbose = T,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b, est.ci = T)
# round(cbind(AA$b$mu,AB$b$mu,AC$b$mu),5)
# GBIC(AB); GBIC(AC)
# GAIC(AB); GAIC(AC)
# AC$lambda
# AC$sse



# Simulation 3: ZI-Poisson:
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# p = 10
# famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- ZIpoisson()}
# famr <- vector("list",p); for(i in 1:p){ famr[[i]] <- Poisson()}
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
# lc$mu[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(lc$mu)){ if(lc$mu[i,2] != 0) lc$mu[i,3] <- rbinom(1,1,0.5) }
# lc$sigma[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(lc$sigma)){ if(lc$sigma[i,2] != 0) lc$sigma[i,3] <- rbinom(1,1,0.5) }
# # ~~~~~~~~~~~~
# # Simulations:
# # ~~~~~~~~~~~~
# AA <- glvmlss_sim(n,famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T, solver = "trust",
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b)
# round(cbind(AA$b$mu,AB$b$mu,AC$b$mu),5)
# round(cbind(AA$b$s,AB$b$s,AC$b$s),5)
# GBIC(AB); GBIC(AC)
# GAIC(AB); GAIC(AC)
# AC$lambda
# AC$sse

# sum(rowSums(AA$Y) == 0)/n
# Y <- AA$Y; b <- AA$b
# 10 -> i; hist(Y[,i],breaks = 25); rbind(b$mu[i,],b$sigma[i,]); (colSums(apply(AA$Y,2,function(x) x == 0))/n)[i]
# round(cbind(AA$b$mu,AC$b$mu,AD$b$mu),5)
# round(cbind(AA$b$si,AC$b$si,AD$b$si),5)
# 
# r1 <- NULL
# for(p in c(10,20)){
# for(n in c(200,500,1000,5000)){
#   r2 <- glvmlss_parsimpost(readRDS(paste0("R/E2n",n,"p",p,".Rds")), p,
#                            out = c("MSE","AB","PVS"),
#                            mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2)
#   r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
# }; #r1; #r1 <- NULL
# }; r1; #rm(r1)

# Simulation 4: Beta:
# ~~~~~~~~~~~~~~~~~~~
p = 10
famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- Beta()}
lc <- NULL
lc$mu <- matrix(1,ncol = 3, nrow = p)
lc$sigma <- matrix(1, ncol = 3, nrow = p)
lc$mu <- matrix(runif(length(lc$mu),0,1),nrow = p)
lc$sigma <- matrix(runif(length(lc$sigma),-1,1), nrow = p)
lc$mu[abs(lc$mu[,1]) < 0.5, 1] <- 0
lc$sigma[, 1] <- runif(p,-1,0)
if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2] ; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3]
if(lc$sigma[1,2] < 0) lc$sigma[1,2] <- -lc$sigma[1,2] ; if(lc$sigma[2,3] < 0) lc$sigma[2,3] <- -lc$sigma[2,3]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For sparse factor loading matrices:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lc$mu[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
for(i in 3:nrow(lc$mu)){ if(lc$mu[i,2] != 0) lc$mu[i,3] <- rbinom(1,1,0.5) }
lc$sigma[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
for(i in 3:nrow(lc$sigma)){ if(lc$sigma[i,2] != 0) lc$sigma[i,3] <- rbinom(1,1,0.5) }
# ~~~~~~~~~~~~
# Simulations:
# ~~~~~~~~~~~~
# AA <- glvmlss_sim(n, famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T,
#               penalty = "alasso", lambda = "auto", w.alasso = AB$b, solver = "trust")
# AD <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T,
#               penalty = "scad", lambda = 0.1, solver = "trust")
# round(cbind(AA$b$mu,AB$b$mu,AC$b$mu),5)
# round(cbind(AA$b$s,AB$b$s,AC$b$s),5)
# GBIC(AB); GBIC(AC); GBIC(AD)
# GAIC(AB); GAIC(AC); GAIC(AD)
# round(AC$lambda,4)
# AC$sse

# data = AA$Y; family = famt; mu.eq = ~ Z1+Z2; sg.eq = ~ Z1+Z2;ta.eq = NULL; nu.eq = NULL;
# control <- list(EM_iter = 30, EM_use2d = T, iter.lim = 300,
#               EM_appHess = F, EM_lrate = 0.001, est.ci = F,
#               solver = "trust", start.val = NULL, mat.info = "Hessian",
#               iden.res = NULL, tol = sqrt(.Machine$double.eps), corr.lv = FALSE,
#               nQP = if(q == 1) 40 else { if(q == 2) ifelse(p <= 10, 15, 25) else 10 },
#               verbose = TRUE, autoL_iter = 30,
#               penalty = "scad", lambda = 0.03, w.alasso = AB$b, gamma = NULL, a = NULL)


# Y <- AA$Y; b <- AA$b
# 15 -> i; hist(Y[,i],breaks = 25); rbind(b$mu[i,],b$sigma[i,]);# (colSums(apply(AA$Y,2,function(x) x == 0))/n)[i]
# cbind(AA$b$mu,AC$b$mu)
# cbind(AA$b$si,AC$b$si)
# r1 <- NULL
# for(p in c(10,20)){
# for(n in c(200,500,1000,5000)){
#   r2 <- glvmlss_parsimpost(readRDS(paste0("R/E3n",n,"p",p,".Rds")), p,
#                            out = c("MSE","AB","PVS"),
#                            mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2)
#   r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
# }; #r1; #r1 <- NULL
# }; r1; rm(r1)
