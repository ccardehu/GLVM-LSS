
rm(list = ls())
set.seed(1234)

library(doSNOW)
source("prep.R")
source("fams.R")
source("glvmlss.R")
source("misc.R")

n = 500     # Number of individuals
nsim = 1000  # Number of simulations

# Simulation 1: Heteroscedastic Normal model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p = 10
famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- Normal()}
lc <- NULL
lc$mu <- matrix(1,ncol = 3, nrow = p)
lc$sigma <- matrix(1, ncol = 3, nrow = p)
lc$mu[,1] <- runif(p,1,2)
lc$mu[,c(2,3)] <- runif(length(lc$mu[,c(2,3)]),0.1,1.5)# * sample(c(-1,1),size = 2*p,replace = T)
lc$sigma <- matrix(runif(length(lc$sigma),0.1,0.4), nrow = p)
lc$sigma[,c(2,3)] <- lc$sigma[,c(2,3)]# * sample(c(-1,1),size = 2*p,replace = T)
# irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0))
# AA <- glvmlss_sim(n,famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, start.val = lc)
# AC <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, verbose = T, start.val = AA$b)

FCOL <- glvmlss_parsimE1(10, T)






# Simulation 2: IRT model:
# ~~~~~~~~~~~~~~~~~~~~~~~~
# p = 10
# famt <- vector("list",p); for(i in 1:p){ famt[[i]] <- Binomial()}
# lc <- NULL
# lc$mu <- matrix(1,ncol = 3, nrow = p)
# lc$mu[,1] <- runif(length(lc$mu[,1]),-1,1)
# lc$mu[,c(2:3)] <- runif(length(lc$mu[,c(2:3)]),1.5,2.5) * sample(c(-1,1),size = 2*p,replace = T)
# if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3];
# # irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0))
# AA <- glvmlss_sim(n,famt, mu.eq = ~Z1+Z2, start.val = lc)
# AB <- glvmlss(data = AA$Y, family = famt, mu.eq = ~ Z1+Z2, verbose = T)

# Simulation 3: ZI-Poisson:
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Copy code from fe_estf.R

 
# Simulation 4: Beta:
# ~~~~~~~~~~~~~~~~~~~
# 
# Copy code from fe_estf.R

