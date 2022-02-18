
rm(list = ls())
# set.seed(1234)
source("f0_prep.R")
# sourceCpp("C:/Users/carde/OneDrive/Desktop/LSE/Research/Research/GitHub/SPLVM/C/testFun.cpp")

n = 3000     # Number of individuals
nsim = 1000  # Number of simulations

# Simulation 1: Heteroscedastic Normal model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s.form <- list("mu" = "~ Z1 + Z2", "sigma" = "~ Z1 + Z2"); p = 30
fam <- rep("normal",p)
l1 <- lc <- NULL
l1$mu <- lc$mu <- matrix(1,ncol = 3, nrow = p)
l1$sigma <- lc$sigma <- matrix(1, ncol = 3, nrow = p)
lc$mu[,1] <- runif(p,1,2)
lc$mu[,c(2,3)] <- runif(length(l1$mu[,c(2,3)]),0.1,1.5)# * sample(c(-1,1),size = 2*p,replace = T)
if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3]
lc$sigma <- matrix(runif(length(l1$sigma),0.1,0.4), nrow = p)
lc$sigma[,c(2,3)] <- lc$sigma[,c(2,3)]# * sample(c(-1,1),size = 2*p,replace = T)
if(ncol(lc$sigma) >= 2 && lc$sigma[1,2] < 0) lc$sigma[1,2] <- -lc$sigma[1,2]
if(ncol(lc$sigma) >= 2 && lc$sigma[2,3] < 0) lc$sigma[2,3] <- -lc$sigma[2,3]


# Simulation 2: IRT model:
# ~~~~~~~~~~~~~~~~~~~~~~~~
# s.form <- list("mu" = "~ Z1+Z2"); p = 30
# fam <- rep("binomial",p)
# l1 <- lc <- NULL
# l1$mu <- lc$mu <- matrix(1,ncol = 3, nrow = p)
# lc$mu[,1] <- runif(length(l1$mu[,1]),-1,1)
# lc$mu[,c(2:3)] <- runif(length(l1$mu[,c(2:3)]),1.5,2.5) * sample(c(-1,1),size = 2*p,replace = T)
# if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3]; #if(lc$mu[3,4] < 0) lc$mu[3,4] <- -lc$mu[3,4]


# Simulation 3: ZI-Poisson:
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# s.form <- list("mu" = "~ Z1 + Z2", "sigma" = "~ Z1 + Z2"); p = 10
# fam <- rep("ZIpoisson", p)
# l1 <- lc <- NULL
# l1$mu <- matrix(1,ncol = 3, nrow = p)
# l1$sigma <- matrix(1, ncol = 3, nrow = p)
# lc$mu <- matrix(runif(length(l1$mu),-1,1),nrow = p)
# lc$sigma <- matrix(runif(length(l1$sigma),0.5,1.5), nrow = p)
# lc$mu[,1] <- runif(p,0.5,1)
# if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2] ; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3]
# lc$sigma[,1] <- runif(p,-3,-1)
 
# Simulation 4: Beta:
# ~~~~~~~~~~~~~~~~~~~
# s.form <- list("mu" = "~ Z1 + Z2", "sigma" = "~ Z1 + Z2"); p = 20
# fam <- rep("beta", p)
# l1 <- lc <- NULL
# l1$mu <- matrix(1,ncol = 3, nrow = p)
# l1$sigma <- matrix(1, ncol = 3, nrow = p)
# lc$mu <- matrix(runif(length(l1$mu),-1,1),nrow = p)
# # lc$mu <- matrix(rep(0.4,length(l1$mu)),nrow = p)
# lc$sigma <- matrix(runif(length(l1$sigma),-1,1), nrow = p)
# # lc$mu[,1] <- runif(p,2,3)
# # lc$sigma[,1] <- runif(p,0.5,1)
# if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2] ; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3]
# if(lc$sigma[1,2] < 0) lc$sigma[1,2] <- -lc$sigma[1,2] ; if(lc$sigma[2,3] < 0) lc$sigma[2,3] <- -lc$sigma[2,3]

# ~~~~~~~~~~~
# Starts here:
# ~~~~~~~~~~~

l1. <- l1

# Parameter restrictions:
# Simulation 1
# ~~~~~~~~~~~~
l1.$mu[1,c(3)] <- l1.$mu[2,c(2)] <- 0
# l1.$mu[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(l1.$mu)){ if(l1.$mu[i,2] != 0) l1.$mu[i,3] <- rbinom(1,1,0.5) }
l1.$sigma[1,3] <- l1.$sigma[2,2] <- 0
# l1.$sigma[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(l1.$sigma)){ if(l1.$sigma[i,2] == 0) l1.$sigma[i,3] <- rbinom(1,1,0.5) }

# Simulation 2
# ~~~~~~~~~~~~
# l1.$mu[1,c(3)] <- l1.$mu[2,c(2)] <- 0
# l1.$mu[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(l1.$mu)){ if(l1.$mu[i,2] != 0) l1.$mu[i,3] <- rbinom(1,1,0.5) }

# Simulation 3
# ~~~~~~~~~~~~
# l1.$mu[1,3] <- l1.$mu[2,2] <- 0
# l1.$mu[sample(3:p, floor((p-3)/3), replace = F),2] <- 0
# for(i in 3:nrow(l1.$mu)){ if(l1.$mu[i,2] != 0) l1.$mu[i,3] <- rbinom(1,1,0.5) }
# l1.$sigma[1,3] <- l1.$sigma[2,2] <- 0
# l1.$sigma[sample(3:p, floor((p-3)/3), replace = F),2] <- 0
# for(i in 3:nrow(l1.$sigma)){ if(l1.$sigma[i,2] != 0) l1.$sigma[i,3] <- rbinom(1,1,0.5) }

# Simulation 4
# ~~~~~~~~~~~~
# l1.$mu[1,c(3)] <- l1.$mu[2,c(2)] <- 0
# l1.$sigma[1,c(3)] <- l1.$sigma[2,c(2)] <- 0
# l1.$mu[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(l1.$mu)){ if(l1.$mu[i,2] != 0) l1.$mu[i,3] <- rbinom(1,1,0.5) }
# l1.$sigma[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
# for(i in 3:nrow(l1.$sigma)){ if(l1.$sigma[i,2] != 0) l1.$sigma[i,3] <- rbinom(1,1,0.5) }


# Simulation procedure:
# ~~~~~~~~~~~~~~~~~~~~~
lc. <- NULL; for(i in names(l1.)){ lc.[[i]] <- lc[[i]]*l1.[[i]] }; sZ <- diag(2) #matrix(c(1,.3,.3,1),2)# randcorr::randcorr(2)
simR <- splvm.sim(n = n, form = s.form, fam = fam, constraints = l1., coefs = lc., sigZ = sZ)
Y <- simR$Y
Z <- simR$Z
borg <- simR$b
e.form <- s.form

# library(xtable)
# xtable(unname(borg$sigma), digits = 3)

# Restrictions (l1.)
# ~~~~~~~~~~~~~~~~~~ 

# Simulation 1: Quadratic heteroscedastic model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0), c("sigma",1,"Z2",0), c("sigma",2,"Z1",0))

# Simulation 2: Bernoulli
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0))

# Simulation 3: ZI-Poisson
# ~~~~~~~~~~~~~~~~~~~~~~~~
# irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0), c("sigma",1,"Z2",0), c("sigma",2,"Z1",0))

# Simulation 4: Beta
# ~~~~~~~~~~~~~~~~~~~~~~~~
# irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0), c("sigma",1,"Z2",0), c("sigma",2,"Z1",0))

# profvis::profvis({
fu <- splvm.fit(Y,fam,e.form,control = list(method = "EM", silent = F, information = "Hessian", startl.val = lc.,
                                            constraint = irestr, corr.lv = F, iter.lim = 1e3))
# })

# a.form <- list("mu" = "~ Z1 + Z2")
# afam = rep("normal",p)
# arestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0))
# fa <- splvm.fit(Y,afam,e.form,control = list(method = "ML", silent = F, information = "Hessian",
#                                             constraint = irestr, corr.lv = F))

# profvis::profvis({
fp <- splvm.fit(Y,fam,e.form, control = list(method = "PML", information = "Hessian", constraint = irestr, start.val = fu$b,
           silent = F, corr.lv = F, pml.control = list(type = "alasso", lambda = "auto", a = 1, gamma = log(n)/2,
           w.alasso = fu$b)))
 # })

round(c(fu$log,fp$log),3)
for(i in names(s.form)){ print(round(cbind(fu$b[[i]],fp$b[[i]],borg[[i]]),3)) }

GBIC(fu); GBIC(fp);
GIC(fu); GIC(fp);


#######################
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~

zsco <- fscore(testc1) # check function

modt <- fu
b1 <- lb2cb(modt$b)[lb2cb(modt$b) != 0]
b2 <- lb2cb(lc.)[lb2cb(lc.) != 0]
num.score <- grad(func = llkf, x = b1, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$ghQ, fam = modt$fam, bg = modt$b)
num.hess <- hessian(func = llkf, x = b1, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$ghQ, fam = modt$fam, bg = modt$b)
ana.res <- sche.test(lb2cb(modt$b),modt)
ana.score <- ana.res$gradient
ana.hess <- ana.res$hessian

xlabn <- NULL
for(i in 1:ncol(modt$Y)){
  xlabn <- append(xlabn,paste0("mu[",i,",",which(modt$b$mu[i,] != 0),"]"))
  if("sigma" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("sigma[",i,",",which(modt$b$sigma[i,] != 0),"]"))
  if("tau" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("tau[",i,",",which(modt$b$tau[i,] != 0),"]"))
  if("nu" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("nu[",i,",",which(modt$b$nu[i,] != 0),"]"))
}

par("mar" = c(6, 4, 4, 2) + 0.1)
plot(ana.score, main = "Score (first deriv. of log-likelihood) @ MLE", xlab = "", ylab = "Value", pch = 16,
     col = "gray50", xaxt = "n")
mtext(side = 1, text = "Parameter index", line = 4)
for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
abline(h = 0, col = "forestgreen", lty = 2)
abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
               by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")
points(num.score[1:length(ana.score)], col = "red", pch = 10)
segments(x0 = seq_along(ana.score), y0 = unlist(ana.score), y1 = num.score, col = "blue", lty = 3, lwd = 1)
legend("topleft",legend=c("Analytical", "Numerical"), col=c("gray50", "red"),
       pch = c(16,10), cex=1, inset = 0.02, box.col = "white")

int <- 1:length(which(lb2cb(modt$b) != 0))^2
# int <- 1:sum(length(modt$b$mu), length(modt$b$sigma))^2
# int <- 1:((sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) * 2)
plot(c(ana.hess)[int], main = "Hessian (2nd deriv. of log-likelihood) @ MLE",xlab ="Hessian entry Index", ylab = "Value",
     cex = 0.5, pch = 16, col = "gray50", ylim = c(min(ana.hess,num.hess),max(ana.hess,num.hess)))#, xaxt = "n")
abline(h = 0, col = "forestgreen", lty = 2);
abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")
points(c(num.hess)[int], col = "red", pch = 16, cex = 0.5)
segments(x0 = seq_along(ana.hess[int]), y0 = unlist(ana.hess)[int], y1 = num.hess[int], col = "blue", lty = 3, lwd = 1)
legend(0,-40000,legend=c("Analytical", "Numerical (Hessian)"), col=c("gray50", "red",0.6),
       pch = c(16,16), cex=0.8, inset = 0.02, box.col = "white")
rm(modt)

#######################

testm <- fp

hist(Z$Z1, breaks = 100, freq = F, xlim = c(-5,5), border = "gray", ylim = c(0,0.7))
hist(Z$Z2, breaks = 100, freq = F, xlim = c(-5,5), border = "gray", ylim = c(0,0.7))
# lines(density(rmtestfac$score.dat$z1), col = 5, lwd = 2)
lines(density(Z$Z1), col = 1, lwd = 2)
lines(density(Z$Z2), col = 1, lwd = 2)
lines(density(fscore(testm)$Z1), col = 2, lwd = 2)
lines(density(fscore(testm)$Z2), col = 2, lwd = 2)
plot(Z$Z1, fscore(testm)$Z1); abline(0,1, col = "red", lty = 3)
plot(Z$Z2, fscore(testm)$Z2); abline(0,1, col = "red", lty = 3)
# plot(Z$Z2, fscore(testa0)$Z2)
cor(Z$Z1, fscore(testm)$Z1,method = "kendall")
cor(Z$Z1, fscore(testm)$Z1,method = "pearson")
cor(Z$Z2, fscore(testm)$Z2,method = "pearson")
points(testm$ghQ$points, testm$ghQ$weights, col = "red", pch = 16)
lines(seq(-5,5,length.out = 100), dnorm(seq(-5,5,length.out = 100)), col = "blue", lwd = 2)

##################
#########################
##########################

Ex1 <- splvm.mcsim(3, sil = F, n = n, form = s.form, fam = fam, restr = l1, coefs = lc,
                   saveRes = F, cleanRes = F)
splvm.plot(1,testa1,
           control.plot = list(sim.mod = simR,
           plot.mean = T, plot.sd = F, qtl = c(0.025,0.975), plot.3D = T, bw = F))

splvm.plot(1,testa3,
           control.plot = list(sim.mod = simR,
           plot.mean = T, plot.sd = F, qtl = c(seq(0.1,0.9,by = 0.2)), plot.3D = T, bw = F))


##############
##############
##############

ftest <- function(par,xt){ dGU(xt,par[1],par[2],log = T)  } #ftest(c(-1,1),-1)
atest <- function(par,xt){ 
 atea <- (exp((xt - par[1])/par[2]) - 1)/par[2]
 ateb <- (((xt - par[1])/par[2]^2) * (exp((xt - par[1])/par[2]) - 1) - (1/par[2]))
return(list(a = atea, b = ateb)) 
}

num.score <- grad(func = ftest, x = c(1,1), method = "Richardson", method.args=list(r = 6, v = 2),
                  xt = 1)
atest(c(-1,1),-1)
optim(par = c(1,1),fn = ftest, gr = atest, xt = -1,method = "L-BFGS-B", lower = 0)
