
rm(list = ls())
set.seed(1234)
# rm(.Random.seed, envir=globalenv())

library(MASS)
library(abind)
library(splines)
library(mvtnorm)
library(numDeriv)
library(gamlss)
library(foreach)
library(doSNOW)
library(parallel)
library(itertools)
library(xtable)

# `%!in%` <- Negate(`%in%`)

source("fa_main.R")
source("fb_auxf.R")
source("fc_mfit.R")
source("fd_mghq.R")
source("ff_msim.R")
source("fg_grph.R")
# source("fi_esim.R")

n = 1000     # Number of individuals
p = 10       # Number of items
nsim = 1000  # Number of simulations
s.form <- list("mu" = "~ Z1*Z2")
fam <- rep("binomial",p)

l1 <- NULL
l1$mu <- matrix(1,ncol = 4, nrow = p)
# l1$sigma <- matrix(1, ncol = 2, nrow = p)

l1. <- l1
# l1.$mu[1,3] <- 0
# l1.$mu[6,2] <- 0
# l1.$mu[sample(1:p, p/2, replace = F),3] <- 0
l1.$mu[sample(1:p, p/2, replace = F),4] <- 0
# l1.$sigma[sample(1:p, p/2, replace = F),2] <- 0
# l1.$sigma[,3] <- 0
# l1.$sigma[sample(which(l1.$sigma[,2] == 1), sum(l1.$sigma[,2] == 1)/2, replace = F),3] <- 1

lc <- NULL
# lc$mu <- matrix(runif(length(l1$mu),-1.5,1.5),nrow = p) # for Normal
lc$mu <- matrix(runif(length(l1$mu),-1.5,1.5),nrow = p) # for Bernoulli
# lc$mu[,1] <- runif(p,1,2) # for Normal
# lc$mu[,2] <- runif(p,-2,2) # for Normal
if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2] # for Normal
lc$mu[abs(lc$mu) < 0.3] <- runif(sum(abs(lc$mu) < 0.3), 0.5,1.0) # for Normal
# lc$mu[,3] <- runif(p,-1,0.5)
# lc$mu[,4] <- runif(p,-0.3,-0.1)
# lc$sigma <- matrix(runif(length(l1$sigma), min = 0.3, max = 0.7), nrow = p) # for Normal
# if(ncol(lc$sigma) >= 2 && lc$sigma[1,2] < 0) lc$sigma[1,2] <- -lc$sigma[1,2] # for Normal
# lc$sigma[abs(lc$sigma) < 0.2] <- runif(sum(abs(lc$sigma) < 0.2), 0.2,0.5) # for Normal
# lc$sigma[,1] <- runif(p,-0.5,0.5)
# lc$sigma[,2] <- runif(p,-0.5,1)
# # lc$sigma[,3] <- runif(p,-0.3,-0.1)

lc. <- NULL
lc.$mu <- lc$mu * l1.$mu
# lc.$sigma <- lc$sigma * l1.$sigma

simR <- splvm.sim(n = n, form = s.form, fam = fam, constraints = l1., coefs = lc.)
Y <- simR$Y
Z <- simR$Z
borg <- simR$b
e.form <- s.form

# Restrictions (l1.)
restr <- list(c("mu",1,"Z1:Z2",0),c("mu",4,"Z1:Z2",0),c("mu",5,"Z1:Z2",0),
              c("mu",6,"Z1:Z2",0),c("mu",10,"Z1:Z2",0))

testa0 <- splvm.fit(Y,fam,e.form,
          control = list(method = "EM", constraint = restr, start.val = lc.,
          tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher"))

testa1 <- splvm.fit(Y,fam,e.form,control = list(method = "EM",
          tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher"))

testa2 <- splvm.fit(Y,fam,e.form, control = list(method = "PEM",
          tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher",
          pml.control = list(type = "alasso", lambda = "auto", w.alasso = testa1$b, pen.load = F)))

testa3 <- splvm.fit(Y,fam,e.form, control = list(method = "PEM",
          tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher",
          pml.control = list(type = "mcp", lambda = 0.5, pen.load = F)))

round(cbind(testa0$b$mu,testa1$b$mu,testa2$b$mu,testa3$b$mu,borg$mu),4)
round(cbind(testa0$b$sigma,testa1$b$sigma,testa2$b$sigma,testa3$b$sigma,borg$sigma),4)

GBIC(testa0);GBIC(testa1); GBIC(testa2); GBIC(testa3)
GIC(testa0);GIC(testa1); GIC(testa2); GIC(testa3)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test for different values of lambda in PML
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lambda <- seq(0.01,0.1,by = 0.01)
pas <- seq(1.5,2.5,by = 0.5)

test <- matrix(NA,ncol=length(pas),nrow = length(lambda))
colnames(test) <- pas; rownames(test) <- lambda
for(i in seq_along(lambda)){
 for(j in seq_along(pas)){
 test[i,j] <- GBIC(splvm.fit(Y = simR$Y, fam = fam, form = simR$form,
                   control = list(method = "PEM", constraint = restrf,
                   ghQqp = 20, iter.lim = 1e3, full.hess = F, tol = sqrt(.Machine$double.eps),
                   silent = T, information = "Fisher",
                   pml.control = list(type = "mcp", lambda = lambda[i], a = pas[j], w.alasso = NULL, pen.load = F))))
 cat(paste0("\n a = ",pas[j]," [",j,"] ","; \U03BB = ",lambda[i]," [",i,"] "))
 }
}

testor <- test

cores <- detectCores()-2
cl<-makeCluster(cores)
registerDoSNOW(cl)
progress <- function(n) setTxtProgressBar(txtProgressBar(max = length(lambda), style = 3), n)
opts <- list(progress = progress)
for(j in seq_along(pas)){
test[,j] <- foreach(i = 1:length(lambda),
                .combine = rbind,
                .packages = c("mvtnorm", "fastGHQuad","gamlss"),
                .options.snow = opts) %dopar% {GBIC(splvm.fit(Y = simR$Y, fam = fam, form = simR$form,
                   control = list(method = "PEM", constraint = restrf,
                   ghQqp = 20, iter.lim = 1e3, full.hess = F, tol = sqrt(.Machine$double.eps),
                   silent = T, information = "Fisher",
                   pml.control = list(type = "mcp", lambda = lambda[i], a = pas[j], w.alasso = NULL, pen.load = F))))}
}
stopCluster(cl)


which(test == min(test), arr.ind = T)

round(cbind(test[[which(lambda == 0.0150)]]$b$mu,
            test[[which(lambda == 0.0125)]]$b$mu,
            test[[which(lambda == 0.0100)]]$b$mu,
            test[[which(lambda == 0.0075)]]$b$mu,
            test[[which(lambda == 0.0050)]]$b$mu),3)

c(GBIC(test[[which(lambda == 0.0150)]]),
GBIC(test[[which(lambda == 0.0125)]]),
GBIC(test[[which(lambda == 0.0100)]]),
GBIC(test[[which(lambda == 0.0075)]]), 
GBIC(test[[which(lambda == 0.0050)]]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test: Comparison vs penfa & lavaan
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
library(penfa)
data(ccdata)

ccdata.subset <- ccdata[ccdata$country == "LEB",2:8]
### Single-group analysis (no mean-structure, unit factor variances)
syntax = 'Z1 =~ h1 + h2 + h3 + h4 + h5 + h6 + h7 
h1 + h2 + h3 + h4 + h5 + h6 + h7 ~ 1'
alasso_fit <- penfa(model = syntax, data = ccdata.subset, std.lv = TRUE,
pen.shrink = "alasso", information = "fisher",
eta = list(shrink = c("lambda" = 0.1), diff = c("none" = 0)),
strategy = "auto", meanstructure = T)
summary(alasso_fit)
penfaParEstim(alasso_fit)
penmat(alasso_fit)

library(lavaan)
CFA.model <- ' Z1 =~ h1 + h2 + h3 + h4 + h5 + h6 + h7 '
rmtest <- cfa(CFA.model, data = ccdata.subset, orthogonal = T, meanstructure = T, std.lv = TRUE, estimator="ML",
              information = "observed"); # coef(rmtest)

fam.ex <- rep("normal",ncol(ccdata.subset))
# ex.form <- list("mu" = "~ Z1-1", "sigma" = "~ 1")
# test <- t(t(ccdata.subset)-colMeans(ccdata.subset)) 
ex.form <- list("mu" = "~ Z1", "sigma" = "~ 1")
test <- ccdata.subset

ex.test.unp <- splvm.fit(test,fam.ex,ex.form,control = list(method = "EM",ghQqp = 50,
               tol = sqrt(.Machine$double.eps), silent = F,information = "Fisher"))

# rmtes2 <- ERP::emfa(as.matrix(ccdata.subset), 1, min.err = .Machine$double.eps, verbose = FALSE, svd.method = c("irlba")) #"irlba" "fast.svd"
# ex.EMfa <- NULL
# ex.EMfa$mu <- matrix(rmtes2$B, dimnames = dimnames(ex.test.unp$b$mu))
# ex.EMfa$sigma <- matrix(log(sqrt(rmtes2$Psi)),dimnames = dimnames(ex.test.unp$b$sigma))

# # ~~~~~~~~~~~~~
ex.ulb <- mb2lb(matrix(coef(rmtest),nrow = ncol(test),byrow = F)[,c(3,1,2)],ex.test.unp$b) #[,c(3,1,2)]
ex.ulb$sigma <- log(sqrt(ex.ulb$sigma))
c(logLik(rmtest)); ex.test.unp$loglik
round(cbind(ex.ulb$mu,ex.test.unp$b$mu),3)
round(cbind(ex.ulb$sigma,ex.test.unp$b$sigma),3)
 
ex.lb <- mb2lb(matrix(coef(alasso_fit),nrow = ncol(test),byrow = F)[,c(2,1,3)],ex.test.unp$b) # [,c(2,1,3)]
ex.lb$sigma <- log(sqrt(ex.lb$sigma))

ex.test.pen <- splvm.fit(test,fam.ex,ex.form,control = list(method = "PEM",
               tol = .Machine$double.eps^(1/2), silent = F,information = "Fisher",
               pml.control = list(type = "alasso", w.alasso = ex.lb, lambda = "auto", pen.load = T)))

# alasso_fit@Options$eta$shrink
# Y = test; fam = fam.ex; form = ex.form
# ex.test.pen$pml.control$lambda

# Penalised comparison
round(cbind(ex.lb$mu,ex.test.pen$b$mu),3)
round(cbind(ex.lb$sigma,ex.test.pen$b$sigma),3)

# Unpenalised comparison
round(cbind(ex.ulb$mu,ex.test.unp$b$mu),3)
round(cbind(ex.ulb$sigma,ex.test.unp$b$sigma),3)

# Penalised vs. Unpenalised
round(cbind(lb2mb(ex.ulb),lb2mb(ex.test.unp$b)),3)
round(cbind(lb2mb(ex.lb),lb2mb(ex.test.pen$b)),3)

GBIC(ex.test.pen); GIC(ex.test.pen)
GBIC(ex.test.unp); GIC(ex.test.unp)
ex.test.unp$loglik; ex.test.pen$loglik

c(alasso_fit@Optim$logl.unpen);ex.test.unp$loglik
c(alasso_fit@Optim$logl.pen); ex.test.pen$loglik;

# ~~~~~~~~~~~~~~~~~~~~~~~
# Test: Comparison vs ltm
# ~~~~~~~~~~~~~~~~~~~~~~~

library(ltm)
test <- WIRS
# test <- Mobility
# data(test)
rmltm <- ltm::ltm(test ~ z1*z2, IRT.param = F,control = list(GHk = 15, iter.em = 350))
# rmltm <- ltm(test ~ z1 + I(z1^2), IRT.param = F,control = list(GHk = 50, iter.em = 350))
# coef(rmltm)
# rmltm$log.Lik

ex.lb <- mb2lb(matrix(coef(rmltm),nrow = ncol(test),byrow = F),ex.test.unp$b)

fam.ex <- rep("binomial",ncol(test))
ex.form <- list("mu" = "~ Z1*Z2")
# ex.form <- list("mu" = "~ Z1 + I(Z1^2)")
ex.test.unp <- splvm.fit(test,fam.ex,ex.form, control = list(method = "EM",
               tol = .Machine$double.eps^0.5, silent = F,information = "Fisher"))

round(cbind(ex.test.unp$b$mu,ex.lb$mu),3)
rmltm$log.Lik; ex.test.unp$loglik

ex.test.pen <- splvm.fit(test,fam.ex,ex.form,
           control = list(method = "PEM", full.hess = F,start.val = ex.test.unp$b,
           ghQqp = 15, iter.lim = 1e3, tol = .Machine$double.eps^0.5, silent = F,information = "Fisher",
           pml.control = list(type = "scad", w.alasso = ex.test.unp$b, lambda = 0.1, pen.load = T)))

round(cbind(ex.test.pen$b$mu,ex.test.unp$b$mu,ex.lb$mu),3)
rmltm$log.Lik; ex.test.unp$loglik; ex.test.pen$loglik

GBIC(ex.test.pen)
GBIC(ex.test.unp)

mod <- NULL
mod$b <- ex.lb
mod$Y <- ex.test.unp$Y
mod$ghQ <- ex.test.unp$ghQ
mod$fam <- ex.test.unp$fam
mod$loglik <- llkf(lb2cb(mod$b),mod$Y,mod$ghQ,mod$b,mod$fam) # rmltm$log.Lik


# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~

zsco <- fscore(testc1) # check function

modt <- ex.test.pen
b1 <- lb2cb(modt$b)
b2 <- lb2cb(lc)
b3 <- rep(1,length(b1))*lb2cb(modt$loadmt)
b4 <- lb2cb(lc.)
num.score <- grad(func = llkf, x = b1, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$ghQ, fam = modt$fam, bg = modt$b)
num.hess <- hessian(func = llkf, x = b1, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$ghQ, fam = modt$fam, bg = modt$b)
ana.res <- sche.test(b1,modt)
ana.score <- ana.res$gradient
ana.hess <- ana.res$hessian

xlabn <- NULL
for(i in 1:ncol(modt$Y)){
  xlabn <- append(xlabn,paste0("mu[",i,",",1:ncol(modt$b$mu),"]"))
  if("sigma" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("sigma[",i,",",1:ncol(modt$b$sigma),"]"))
  if("tau" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("tau[",i,",",1:ncol(modt$b$tau),"]"))
  if("nu" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("nu[",i,",",1:ncol(modt$b$nu),"]"))
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

int <- 1:sum(length(modt$b$mu), length(modt$b$sigma))^2
# int <- 1:((sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) * 2)
plot(c(ana.hess)[int], main = "Hessian (2nd deriv. of log-likelihood) @ MLE",xlab ="Hessian entry Index", ylab = "Value",
     cex = 0.5, pch = 16, col = "gray50", ylim = c(min(ana.hess,num.hess),max(ana.hess,num.hess)))#, xaxt = "n")
abline(h = 0, col = "forestgreen", lty = 2);
abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")
points(c(num.hess)[int], col = "red", pch = 16, cex = 0.5)
segments(x0 = seq_along(ana.hess[int]), y0 = unlist(ana.hess)[int], y1 = num.hess[int], col = "blue", lty = 3, lwd = 1)
legend(0,-2500,legend=c("Analytical", "Numerical (Hessian)"), col=c("gray50", "red",0.6),
       pch = c(16,16), cex=0.8, inset = 0.02, box.col = "white")
rm(modt)


rbenchmark::benchmark(
  "test20" = {testa <- splvm.fit(Y,fam,e.form,
          control = list(method = "EM", constraint = l1, full.hess = F, start.val = lc,
          ghQqp = 20, iter.lim = 200, tol = sqrt(.Machine$double.eps), silent = F))},
  "test40" = {testa <- splvm.fit(Y,fam,e.form,
          control = list(method = "EM", constraint = l1, full.hess = F, start.val = lc,
          ghQqp = 40, iter.lim = 200, tol = sqrt(.Machine$double.eps), silent = F))},
  replications = 3, columns = c("test", "replications", "elapsed", "relative", "user.self", "sys.self")
)


hist(Z$Z1, breaks = 100, freq = F, xlim = c(-5,5), border = "gray", ylim = c(0,0.7))
hist(Z$Z2, breaks = 100, freq = F, xlim = c(-5,5), border = "gray", ylim = c(0,0.7))
# lines(density(rmtestfac$score.dat$z1), col = 5, lwd = 2)
lines(density(Z$Z1), col = 1, lwd = 2)
lines(density(Z$Z2), col = 1, lwd = 2)
lines(density(fscore(testa0)$Z1), col = 2, lwd = 2)
lines(density(fscore(testa0)$Z2), col = 2, lwd = 2)
plot(Z$Z1, fscore(testa0)$Z1); abline(0,1, col = "red", lty = 3)
plot(Z$Z2, fscore(testa0)$Z2); abline(0,1, col = "red", lty = 3)
# plot(Z$Z2, fscore(testa0)$Z2)
cor(Z$Z1, fscore(testa0)$Z1,method = "kendall")
cor(Z$Z1, fscore(testa0)$Z1,method = "pearson")
cor(Z$Z2, fscore(testa0)$Z2,method = "pearson")
points(testa0$ghQ$points, testa0$ghQ$weights, col = "red", pch = 16)
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
