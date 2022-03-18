
rm(list = ls())
set.seed(1234)

library(doSNOW)
source("R/prep.R")
source("R/fams.R")
source("R/glvmlss.R")
source("R/misc.R")

n = 500     # Number of individuals
nsim = 1000  # Number of simulations

# Simulation 1: Heteroscedastic Normal model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# p = 10
# famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- Normal()}
# lc <- NULL
# lc$mu <- matrix(1,ncol = 3, nrow = p)
# lc$sigma <- matrix(1, ncol = 1, nrow = p)
# lc$mu[,1] <- runif(p,1,2)
# lc$mu[,c(2,3)] <- runif(length(lc$mu[,c(2,3)]),0.1,1.5) * sample(c(-1,1),size = 2*p,replace = T)
# lc$sigma <- matrix(runif(length(lc$sigma),0.1,0.4), nrow = p)
# # lc$sigma[,c(2,3)] <- lc$sigma[,c(2,3)] * sample(c(-1,1),size = 2*p,replace = T)
# lc$mu <- lc$mu[,-1];# lc$sigma <- lc$sigma[,-1]
# # irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0))
# AA <- glvmlss_sim(n,famt, mu.eq = ~ Z1+Z2-1, sg.eq = ~ 1, start.val = lc)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2-1, sg.eq = ~ 1, verbose = T)
# # rmA <- fa(scale(AA$Y), nfactors=2, fm="ml", rotate="varimax")
# cbind(AA$b$mu,AC$b$mu);
# cbind(AA$b$s,AC$b$s)
# # # FCOL <- glvmlss_parsimE1(nsim, T)
# # r1 <- NULL
# # for(p in c(10,20)){
# # for(n in c(200,500,1000,5000)){
# #   r2 <- glvmlss_parsimpost(readRDS(paste0("R/E1n",n,"p",p,".Rds")), p,
# #                            out = c("MSE","AB"),
# #                            mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2)
# #   r1 <- rbind(r1,r2,deparse.level = 0); rm(r2)
# # }; #r1; #r1 <- NULL
# # }; r1; rm(r1)


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
# AA <- glvmlss_sim(n,famt, mu.eq = ~Z1+Z2, start.val = lc, iden.res = irestr)
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, start.val = AA$b, verbose = T)
# cbind(AA$b$mu,AB$b$mu)


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
# AA <- glvmlss_sim(n,famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T, start.val = AA$b, nQP = 25, solver = "trust")
# 
# sum(rowSums(AA$Y) == 0)/n
# Y <- AA$Y; b <- AA$b
# 1 -> i; hist(Y[,i],breaks = 25); rbind(b$mu[i,],b$sigma[i,]); (colSums(apply(AA$Y,2,function(x) x == 0))/n)[i]
# cbind(AA$b$mu,AC$b$mu)
# cbind(AA$b$si,AC$b$si)

# Simulation 4: Beta:
# ~~~~~~~~~~~~~~~~~~~
p = 10
famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- Beta()}
lc <- NULL
lc$mu <- matrix(1,ncol = 3, nrow = p)
lc$sigma <- matrix(1, ncol = 3, nrow = p)
lc$mu <- matrix(runif(length(lc$mu),-1,1),nrow = p)
lc$sigma <- matrix(runif(length(lc$sigma),-1,1), nrow = p)
if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2] ; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3]
if(lc$sigma[1,2] < 0) lc$sigma[1,2] <- -lc$sigma[1,2] ; if(lc$sigma[2,3] < 0) lc$sigma[2,3] <- -lc$sigma[2,3]
AA <- glvmlss_sim(n, famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T, solver = "trust", nQP = 25)

Y <- AA$Y; b <- AA$b
1 -> i; hist(Y[,i],breaks = 25); rbind(b$mu[i,],b$sigma[i,]);# (colSums(apply(AA$Y,2,function(x) x == 0))/n)[i]
cbind(AA$b$mu,AC$b$mu)
cbind(AA$b$si,AC$b$si)
