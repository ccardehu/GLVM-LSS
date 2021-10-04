# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Example No.0 : Holzinger & Swineford data #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

rm(list = ls())
set.seed(1234)
source("f0_prep.R")

data.0 <- lavaan::HolzingerSwineford1939
data.0 <- data.0[data.0$school == "Grant-White",]
data.0 <- data.0[,-c(1:6)]
colnames(data.0) <- paste0("Y",1:ncol(data.0)); rownames(data.0) <- NULL; data.0 <- as.matrix(data.0)

fam.0 <- rep("normal",ncol(data.0))
form.0a <- list(mu = "~ Z1 + Z2 + Z3", sigma = "~ 1")
irst.0a0 <- list(c("mu",1,"Z2",0), c("mu",1,"Z3",0),c("mu",2,"Z2",0), c("mu",2,"Z3",0),c("mu",3,"Z2",0), c("mu",3,"Z3",0),
                 c("mu",4,"Z1",0), c("mu",4,"Z3",0),c("mu",5,"Z1",0), c("mu",5,"Z3",0),c("mu",6,"Z1",0), c("mu",6,"Z3",0),
                 c("mu",7,"Z1",0), c("mu",7,"Z2",0),c("mu",8,"Z1",0), c("mu",8,"Z2",0),c("mu",9,"Z1",0), c("mu",9,"Z2",0))
irst.0a1 <- list(c("mu",1,"Z2",0), c("mu",1,"Z3",0),
                 c("mu",4,"Z1",0), c("mu",4,"Z3",0),
                 c("mu",7,"Z1",0), c("mu",7,"Z2",0))
fit.0au0 <- splvm.fit(data.0,fam.0,form.0a,control = list(method = "ML", silent = F, information = "Fisher", constraints = irst.0a0, ghQp = 5))
fit.0au1 <- splvm.fit(data.0,fam.0,form.0a,control = list(method = "ML", silent = F, information = "Fisher", constraints = irst.0a1, ghQp = 5))
fit.0ap1 <- splvm.fit(data.0,fam.0,form.0a,control = list(method = "PML", silent = F, information = "Fisher", constraint = irst.0a1, ghQp = 5,
                      pml.control = list(type = "alasso", lambda = "auto", w.alasso = fit.0au1$b, gamma = 4)))

fit.0au0$loglik; fit.0au1$loglik; fit.0ap1$loglik
GIC(fit.0au0); GIC(fit.0au1); GIC(fit.0ap1);
GBIC(fit.0au0); GBIC(fit.0au1); GBIC(fit.0ap1);
round(cbind(fit.0au0$b$mu%*%diag(c(1,1,-1,-1)),fit.0au1$b$mu%*%diag(c(1,1,-1,-1)),fit.0ap1$b$mu%*%diag(c(1,1,-1,-1))),4)
round(cbind(fit.0au0$b$si,fit.0au1$b$si,fit.0ap1$b$si),4)

# Model with heteroscedasticity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

form.0b <- list(mu = "~ Z1 + Z2 + Z3", sigma = "~ Z1 + Z2 + Z3")
irst.0b0 <- list(c("mu",1,"Z2",0), c("mu",1,"Z3",0),c("mu",2,"Z2",0), c("mu",2,"Z3",0),c("mu",3,"Z2",0), c("mu",3,"Z3",0),
                 c("mu",4,"Z1",0), c("mu",4,"Z3",0),c("mu",5,"Z1",0), c("mu",5,"Z3",0),c("mu",6,"Z1",0), c("mu",6,"Z3",0),
                 c("mu",7,"Z1",0), c("mu",7,"Z2",0),c("mu",8,"Z1",0), c("mu",8,"Z2",0),c("mu",9,"Z1",0), c("mu",9,"Z2",0),
                 c("sigma",1,"Z2",0), c("sigma",1,"Z3",0),c("sigma",2,"Z2",0), c("sigma",2,"Z3",0),c("sigma",3,"Z2",0), c("sigma",3,"Z3",0),
                 c("sigma",4,"Z1",0), c("sigma",4,"Z3",0),c("sigma",5,"Z1",0), c("sigma",5,"Z3",0),c("sigma",6,"Z1",0), c("sigma",6,"Z3",0),
                 c("sigma",7,"Z1",0), c("sigma",7,"Z2",0),c("sigma",8,"Z1",0), c("sigma",8,"Z2",0),c("sigma",9,"Z1",0), c("sigma",9,"Z2",0))
irst.0b1 <- list(c("mu",1,"Z2",0), c("mu",1,"Z3",0),
                 c("mu",4,"Z1",0), c("mu",4,"Z3",0),
                 c("mu",7,"Z1",0), c("mu",7,"Z2",0),
                 c("sigma",1,"Z2",0), c("sigma",1,"Z3",0),
                 c("sigma",4,"Z1",0), c("sigma",4,"Z3",0),
                 c("sigma",7,"Z1",0), c("sigma",7,"Z2",0))
fit.0bu0 <- splvm.fit(data.0,fam.0,form.0b,control = list(method = "ML", silent = F, information = "Fisher", constraints = irst.0b0, ghQp = 5))
fit.0bp0 <- splvm.fit(data.0,fam.0,form.0b,control = list(method = "PML", silent = F, information = "Fisher", constraint = irst.0b0, ghQp = 5,
                      pml.control = list(type = "alasso", lambda = "auto", w.alasso = fit.0bu0$b, gamma = 3)))
fit.0bu1 <- splvm.fit(data.0,fam.0,form.0b,control = list(method = "ML", silent = F, information = "Fisher", constraints = irst.0b1, ghQp = 5))
fit.0bp1 <- splvm.fit(data.0,fam.0,form.0b,control = list(method = "PML", silent = F, information = "Fisher", constraint = irst.0b1, ghQp = 5,
                      pml.control = list(type = "alasso", lambda = "auto", w.alasso = fit.0bu1$b, gamma = 3), iter.lim = 200))

fit.0bu0$loglik; fit.0bp0$loglik
fit.0bu1$loglik; fit.0bp1$loglik
GBIC(fit.0bu0); GBIC(fit.0bp0);
GBIC(fit.0bu1); GBIC(fit.0bp1);
round(cbind(fit.0bu0$b$mu%*%diag(c(1,1,-1,-1)),fit.0bp0$b$mu%*%diag(c(1,1,-1,-1))),4)
round(cbind(fit.0bu0$b$si%*%diag(c(1,1,-1,-1)),fit.0bp0$b$si%*%diag(c(1,1,-1,-1))),4)
round(cbind(fit.0bu1$b$mu%*%diag(c(1,1,-1,-1)),fit.0bp1$b$mu%*%diag(c(1,1,-1,-1))),4)
round(cbind(fit.0bu1$b$si%*%diag(c(1,1,-1,-1)),fit.0bp1$b$si%*%diag(c(1,1,-1,-1))),4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Example No.1 : European Social Survey - ESS #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rm(list = ls())
set.seed(1234)
source("f0_prep.R")
library(haven)

data.1 <- read_dta("C:/Users/carde/Dropbox/Camilo and Irini/Research/Data/ESS/Immigration/trPOL.dta")
data.1 <- as.data.frame(data.1[complete.cases(data.1),-c(1,ncol(data.1))])
ex.d <- c(1:4,8:13)
data.1 <- data.1[,ex.d]; #colnames(data.1) <- paste0("Y",1:ncol(data.1))

plotV <- T
if(plotV) { X11(); pairs(data.1,
      upper.panel = panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex =  cex.cor * (1 + r) / 2)
},
      diag.panel  = panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)
} ) }

fam.1a <- c(rep("ZIpoisson", ncol(data.1)))
form.1a <- list(mu = "~ Z1", sigma = "~ Z1")
form.1b <- list(mu = "~ Z1 + Z2", sigma = "~ Z1 + Z2")
irst.1b0 <- list(c("mu",1,"Z2",0), c("mu",5,"Z1",0), c("sigma",1,"Z2",0), c("sigma",5,"Z1",0))

# ZIP Model
# ~~~~~~~~~

fit.1au0 <- splvm.fit(data.1,fam.1a,form.1a,control = list(method = "EM", silent = F))
fit.1au1 <- splvm.fit(data.1,fam.1a,form.1a,control = list(method = "ML", silent = F, information = "Hessian"))

fit.1bu0 <- splvm.fit(data.1,fam.1a,form.1b,control = list(method = "EM", silent = F, constraints = irst.1b0))
fit.1bu1 <- splvm.fit(data.1,fam.1a,form.1b,control = list(method = "ML", silent = F, information = "Fisher", constraints = irst.1b0))
fit.1bp1 <- splvm.fit(data.1,fam.1a,form.1b,control = list(method = "PML", silent = F, information = "Fisher", constraints = irst.1b0,
                      pml.control = list(type = "alasso", lambda = "auto", w.alasso = fit.1bu1$b, gamma = 2)))

fit.1au1$loglik; fit.1bu1$loglik; fit.1bp1$loglik
GBIC(fit.1au1); GBIC(fit.1bu1); GBIC(fit.1bp1)
round(cbind(fit.1au1$b$mu,fit.1bu1$b$mu,fit.1bp1$b$mu),4)
round(cbind(fit.1au1$b$si,fit.1bu1$b$si,fit.1bp1$b$si),4)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Example No.2 : Rizopoulos & Moustaki (2008) : WIRS #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rm(list = ls())
set.seed(1234)

source("f0_prep.R")

library(xtable)
library(ltm)
data.2 <- WIRS
rmltm <- ltm::ltm(data.2 ~ z1*z2, IRT.param = F,control = list(GHk = 40, iter.em = 350),constraint = cbind(1:6,4,-1.64))
# constraint = cbind(1:6,4,-1.64)
# rmltm0 <- ltm::ltm(test ~ z1+z2, IRT.param = F,control = list(GHk = 31, iter.em = 350))
# coef(rmltm)
# coef(rmltm0)
# rmltm$log.Lik
# rmltm0$log.Lik
sd(rep(factor.scores(rmltm, method = "EAP")$score.dat[,"z1"],factor.scores(rmltm, method = "EAP")$score.dat[,"Obs"]))
mod <- NULL
mod$b <- ex.lb
mod$Y <- ex.test.pen[[1]][[1]]$Y
mod$ghQ <- ourmod2$ghQ#ex.test.pen[[1]][[1]]$ghQ
mod$fam <- ex.test.pen[[1]][[1]]$fam
mod$pml.control <- ex.test.pen[[1]][[1]]$pml.control
mod$uploglik <- llkf(lb2cb(mod$b),mod$Y,mod$ghQ,mod$b,mod$fam) # Difference with penfa; c(alasso_fit@Optim$logl.unpen)
mod$info = "Hessian"
mod$hessian <- sche.test(lb2cb(mod$b),mod)$hessian
GBIC(mod)

ex.lb <- mb2lb(matrix(coef(rmltm),nrow = ncol(data.2),byrow = F),ex.test.unp$b)
# ex.lb0 <- NULL; ex.lb0$mu <- ex.lb$mu[,-4]

fam.ex <- rep("binomial",ncol(data.2))

# ex.form0 <- list("mu" = "~ Z1+Z2")
# ex.test.unp0 <- splvm.fit(test,fam.ex,ex.form0,
#                control = list(method = "EM", full.hess = F, start.val = ex.lb0,
#                ghQqp = 15, iter.lim = 2e3, tol = .Machine$double.eps^0.5, silent = F,information = "Fisher"))

ex.form <- list("mu" = "~ Z1*Z2")
ex.test.unp <- splvm.fit(data.2,fam.ex,ex.form, control = list(method = "ML", ghQp = 15,
               iter.lim = 2e3, silent = F))
ex.test.pen <- splvm.fit(data.2,fam.ex,ex.form, control = list(method = "PML",
               iter.lim = 2e3, silent = F,
               pml.control = list(type = "alasso", w.alasso = ex.test.unp$b, lambda = "auto", gamma = 3.5)))

round(cbind(ex.lb$mu,ex.test.unp$b$mu,ex.test.pen$b$mu),3)
c(rmltm$log.Lik,ex.test.unp$loglik,ex.test.pen$loglik)
GBIC(ex.test.unp);GBIC(ex.test.pen)
GIC(ex.test.unp);GIC(ex.test.pen)


ourmod <- ex.test.pen[[2]][[7]]
round(cbind(ourmod2$b$mu,mod$b$mu),2)
splvm.plot(2, ourmod2, control.plot = list(plot.3D = T, bw = F))
splvm.plot(2, mod, control.plot = list(plot.3D = T, bw = F))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Example No.3 : Niku et al. (2019) : AntTraits #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rm(list = ls())
set.seed(1234)

source("f0_prep.R")

library(mvabund)
data(antTraits)
data.3 <- unname(as.matrix(antTraits$abund)); rm(antTraits)
data.3 <- data.3[,-9]
fitp <- gllvm::gllvm(data.3, family = "poisson", num.lv = 2); # summary(fitp)

fam.3 <- rep("poisson", ncol(data.3))
form.3 <- list(mu = "~ Z1 + Z2")
sv.3 <- list(mu = matrix(c(unname(fitp$params$beta0), unname(fitp$params$theta)), nrow = ncol(data.3)))#,
irestr.3 <- list(c("mu",1,"Z2",0))

fit.3u <- splvm.fit(Y = data.3, fam = fam.3, form = form.3,
                    control = list(method = "ML", silent = F, information = "Fisher", constraint = irestr.3, start.val = sv.3))
fit.3p <- splvm.fit(data.3,fam.3,form.3,control = list(method = "PML", silent = F, information = "Fisher", constraint = irestr.3, start.val = sv.3,
                    pml.control = list(type = "alasso", lambda = "auto", w.alasso = fit.3u$b, gamma = 2.5)))

logLik(fitp)[1]; fit.3u$loglik; fit.3p$loglik
BIC(fitp); GBIC(fit.3u); GBIC(fit.3p);
AIC(fitp); GIC(fit.3u); GIC(fit.3p)
round(cbind(sv.3$mu,fit.3u$b$mu,fit.3p$b$mu),4)

pairs(data.3,
      upper.panel = panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex =  cex.cor * (1 + r) / 2)
},
      diag.panel  = panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)
} )
