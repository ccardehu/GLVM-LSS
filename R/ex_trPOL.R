
rm(list = ls())
set.seed(1234)

library(haven)
library(mvtnorm)
library(numDeriv)

source("ddFun.R")
source("SimFA.R")
source("GHFun.R")
source("scFun.R")
source("GLVMfit.R")
source("graphFun.R")
source("EMFun.R")

trPOL <- read_dta("C:/Users/carde/Dropbox/Camilo and Irini/Research/Data/ESS/Immigration/trPOL.dta")
trPOL <- as.data.frame(trPOL[complete.cases(trPOL),-c(1,ncol(trPOL))])
summary(trPOL)
# pairs(trPOL)

# testa <- Jmisc::demean(trPOL[,-c(9:11)])
testa <- trPOL[,-c(9:11)]
psych::fa(testa,1, fm="mle", rotate = "none")$loadings
psych::fa(testa,1, fm="mle", rotate = "none")$uniquenesses
fam00 <- rep("normal", ncol(testa))
formd00 <- NULL; formd00$mu <- "~ Z1"; formd00$sigma <- "~ 1"
l001. <- l002. <- NULL
l001.$mu <- matrix(1,ncol = 2, nrow = ncol(testa))
l002.$mu <- matrix(0,ncol = 2, nrow = ncol(testa))
l001.$sigma <- matrix(1, ncol = 1, nrow = ncol(testa))
l002.$sigma <- matrix(0, ncol = 1, nrow = ncol(testa))
test00 <- GLVM.fit(Y = testa, fam = fam00, form = formd00, silent = F, ghp = 30, iter.lim = 1000, tol = 1e-7, loadmt = l001., icoefs = l002., useoptim = F, skipEM = F)
test00$b

testb <- trPOL[,c(9:11)]
fam10 <- rep("binomial", ncol(testb))
l101. <- l102. <- NULL
l101.$mu <- matrix(1,ncol = 2, nrow = ncol(testb))
l102.$mu <- matrix(1,ncol = 2, nrow = ncol(testb))
test10 <- GLVM.fit(Y = testb, fam = fam10, form = formd00 , silent = F, ghp = 30, iter.lim = 1000, tol = 1e-7, loadmt = l101., icoefs = l102., useoptim = F, skipEM = F)
round(test10$b$mu,4)
test10$loglik
ltm(testb ~ z1, IRT.param = F, control = list(GHk = 30))

p. <- ncol(trPOL)
famd <- c(rep("normal", 8), rep("binomial", 3))
l01. <- l02. <- l11. <- l12. <- l21. <- l22. <- l31. <- l32. <- NULL
formd0 <- NULL; formd0$mu <- "~ Z1"; formd0$sigma <- "~ 1"
formd1 <- NULL; formd1$mu <- "~ Z1"; formd1$sigma <- "~ Z1"
formd2 <- NULL; formd2$mu <- "~ Z1 + I(Z1^2)"; formd2$sigma <- "~ Z1"
formd3 <- NULL; formd3$mu <- "~ Z1 + Z2 "; formd3$sigma <- "~ Z1"
l01.$mu <- matrix(1,ncol = 2, nrow = p.)
l02.$mu <- matrix(0,ncol = 2, nrow = p.)
l01.$sigma <- matrix(1, ncol = 1, nrow = p.)
l02.$sigma <- matrix(0, ncol = 1, nrow = p.)
l11.$mu <- matrix(1,ncol = 2, nrow = p.)
l12.$mu <- matrix(0,ncol = 2, nrow = p.)
l11.$sigma <- matrix(1, ncol = 2, nrow = p.)
l12.$sigma <- matrix(0, ncol = 2, nrow = p.)
l21.$mu <- matrix(1,ncol = 3, nrow = p.)
l22.$mu <- matrix(0,ncol = 3, nrow = p.)
l21.$sigma <- matrix(1, ncol = 2, nrow = p.)
l22.$sigma <- matrix(0, ncol = 2, nrow = p.)
l31.$mu <- matrix(1,ncol = 3, nrow = p.)
l32.$mu <- matrix(0,ncol = 3, nrow = p.)
l31.$sigma <- matrix(1, ncol = 2, nrow = p.)
l32.$sigma <- matrix(0, ncol = 2, nrow = p.)

test0 <- GLVM.fit(Y = trPOL, fam = famd, form = formd0 , silent = F, ghp = 30, iter.lim = 1000, tol = 1e-7, loadmt = l01., icoefs = l02., useoptim = F, skipEM = F)
test1 <- GLVM.fit(Y = trPOL, fam = famd, form = formd1 , silent = F, ghp = 30, iter.lim = 1000, tol = 1e-7, loadmt = l11., icoefs = l12., useoptim = F, skipEM = F)
test2 <- GLVM.fit(Y = trPOL, fam = famd, form = formd2 , silent = F, ghp = 30, iter.lim = 1000, tol = 1e-7, loadmt = l21., icoefs = l22., useoptim = F, skipEM = F)
test3 <- GLVM.fit(Y = trPOL, fam = famd, form = formd3 , silent = F, ghp = 10, iter.lim = 1000, tol = 1e-7, loadmt = l31., icoefs = l32., useoptim = F, skipEM = F)

test0$b;
test1$b;
test2$b;
test3$b;
test0$loglik; test1$loglik; test2$loglik; test3$loglik

round(Matrix::tril(var(trPOL)),2)
round(Matrix::tril(cor(trPOL)),2)

