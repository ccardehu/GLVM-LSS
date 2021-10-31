
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# Simulation 1: Heteroscedastic Normal model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
source("f0_prep.R")
s.form <- list("mu" = "~ Z1 + Z2", "sigma" = "~ Z1 + Z2")

# n = 200; p = c(10,30)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("OldSim/Fa/Ex1_nsim500_n200_p10.RData"); p = 10
r200p10 <- post(FCOL,s.form,p,T)
load("OldSim/Fa/Ex1_nsim500_n200_p30.RData"); p = 30
r200p30 <- post(FCOL,s.form,p,T)
# n = 500; p = 10
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("OldSim/Fa/Ex1_nsim500_n500_p10.RData"); p = 10
r500p10 <- post(FCOL,s.form,p,T)
load("OldSim/Fa/Ex1_nsim500_n500_p30.RData"); p = 30
r500p30 <- post(FCOL,s.form,p,T)
# n = 1000; p = 10
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("OldSim/Fa/Ex1_nsim500_n1000_p10.RData"); p = 10
r1000p10 <- post(FCOL,s.form,p,T)
load("OldSim/Fa/Ex1_nsim500_n1000_p30.RData"); p = 30
r1000p30 <- post(FCOL,s.form,p,T)

rbind(r200p10$res,r500p10$res, r1000p10$res)
rbind(r200p30$res,r500p30$res, r1000p30$res)

# write.table(rbind(r200p10$res,r500p10$res, r1000p10$res), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
# write.table(rbind(r200p30$res,r500p30$res, r1000p30$res), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 2: IRT model
# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
source("f0_prep.R")
s.form <- list("mu" = "~ Z1+Z2")

# n = 200; p = c(10,30)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("Ex2_nsim500_n200_p10.RData"); p = 10
r200p10 <- post(FCOL,s.form,p,T)
load("Ex2_nsim500_n200_p30.RData"); p = 30
r200p30 <- post(FCOL,s.form,p,T)
# n = 500; p = 10
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("Ex2_nsim500_n500_p10.RData"); p = 10
r500p10 <- post(FCOL,s.form,p,T)
load("Ex2_nsim500_n500_p30.RData"); p = 30
r500p30 <- post(FCOL,s.form,p,T)
# n = 1000; p = 10
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("Ex2_nsim500_n1000_p10.RData"); p = 10
r1000p10 <- post(FCOL,s.form,p,T)
load("Ex2_nsim500_n1000_p30.RData"); p = 30
r1000p30 <- post(FCOL,s.form,p,T)

rbind(r200p10$res,r500p10$res, r1000p10$res)
rbind(r200p30$res,r500p30$res, r1000p30$res)

# write.table(rbind(r200p10$res,r500p10$res, r1000p10$res), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
# write.table(rbind(r200p30$res,r500p30$res, r1000p30$res), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)