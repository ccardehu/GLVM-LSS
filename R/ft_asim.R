
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# Simulation 1: Heteroscedastic Normal model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
source("f0_prep.R")
form <- list("mu" = "~ Z1 + Z2", "sigma" = "~ Z1 + Z2")

# n = 200; p = c(10,30)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("Ex1_nsim500_n200_p10_a1.RData"); p = 10
r200p10a1 <- post(FCOL,form,p,T,1)
load("Ex1_nsim500_n200_p10_a2.RData"); p = 10
r200p10a2 <- post(FCOL,form,p,T,1)
load("Ex1_nsim500_n200_p10_a3.RData"); p = 10
r200p10a3 <- post(FCOL,form,p,T,1)

load("Ex1_nsim500_n200_p30_a1.RData"); p = 30
r200p30a1 <- post(FCOL,form,p,T,1)
load("Ex1_nsim500_n200_p30_a2.RData"); p = 30
r200p30a2 <- post(FCOL,form,p,T,1)
load("Ex1_nsim500_n200_p30_a3.RData"); p = 30
r200p30a3 <- post(FCOL,form,p,T,1)

# n = 500; p = 10
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("Ex1_nsim500_n500_p10_a1.RData"); p = 10
r500p10a1 <- post(FCOL,form,p,T,1)
load("Ex1_nsim500_n500_p10_a2.RData"); p = 10
r500p10a2 <- post(FCOL,form,p,T,1)
load("Ex1_nsim500_n500_p10_a3.RData"); p = 10
r500p10a3 <- post(FCOL,form,p,T,1)

load("Ex1_nsim500_n500_p30_a1.RData"); p = 30
r500p30a1 <- post(FCOL,form,p,T,1)
load("Ex1_nsim500_n500_p30_a2.RData"); p = 30
r500p30a2 <- post(FCOL,form,p,T,1)
load("Ex1_nsim500_n500_p30_a3.RData"); p = 30
r500p30a3 <- post(FCOL,form,p,T,1)

# n = 1000; p = 10
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("Ex1_nsim500_n1000_p10_a1.RData"); p = 10
r1000p10a1 <- post(FCOL,form,p,T,1)
load("Ex1_nsim500_n1000_p10_a2.RData"); p = 10
r1000p10a2 <- post(FCOL,form,p,T,1)
load("Ex1_nsim500_n1000_p10_a3.RData"); p = 10
r1000p10a3 <- post(FCOL,form,p,T,1)

load("Ex1_nsim500_n1000_p30_a1.RData"); p = 30
r1000p30a1 <- post(FCOL,form,p,T,1)
load("Ex1_nsim500_n1000_p30_a2.RData"); p = 30
r1000p30a2 <- post(FCOL,form,p,T,1)
load("Ex1_nsim500_n1000_p30_a3.RData"); p = 30
r1000p30a3 <- post(FCOL,form,p,T,1)

round(rbind(cbind(r200p10a1$res,r200p10a1$cr),cbind(r500p10a1$res,r500p10a1$cr),cbind(r1000p10a1$res,r1000p10a1$cr)),4)
round(rbind(cbind(r200p10a2$res,r200p10a2$cr),cbind(r500p10a2$res,r500p10a2$cr),cbind(r1000p10a2$res,r1000p10a2$cr)),4)
round(rbind(cbind(r200p10a3$res,r200p10a3$cr),cbind(r500p10a3$res,r500p10a3$cr),cbind(r1000p10a3$res,r1000p10a3$cr)),4)

round(rbind(cbind(r200p30a1$res,r200p30a1$cr),cbind(r500p30a1$res,r500p30a1$cr),cbind(r1000p30a1$res,r1000p30a1$cr)),4)
round(rbind(cbind(r200p30a2$res,r200p30a2$cr),cbind(r500p30a2$res,r500p30a2$cr),cbind(r1000p30a2$res,r1000p30a2$cr)),4)
round(rbind(cbind(r200p30a3$res,r200p30a3$cr),cbind(r500p30a3$res,r500p30a3$cr),cbind(r1000p30a3$res,r1000p30a3$cr)),4)


rbind(r200p30a2$res,r500p30a2$res, r1000p30a2$res)
rbind(r200p30a2$res,r500p30a2$res, r1000p30a2$res)
rbind(r200p30a3$res,r500p30a3$res, r1000p30a3$res)


# write.table(rbind(r200p10a1$res,r500p10a1$res, r1000p10a1$res,r200p30a1$res,r500p30a1$res, r1000p30a1$res), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
# write.table(rbind(r200p10a2$res,r500p10a2$res, r1000p10a2$res,r200p30a2$res,r500p30a2$res, r1000p30a2$res), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
# write.table(rbind(r200p10a3$res,r500p10a3$res, r1000p10a3$res,r200p30a3$res,r500p30a3$res, r1000p30a3$res), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 2: IRT model
# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
source("f0_prep.R")
form <- list("mu" = "~ Z1+Z2")

# n = 200; p = c(10,30)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("Ex2_nsim500_n200_p10_a1.RData"); p = 10
r200p10a1 <- post(FCOL,form,p,T,1)
load("Ex2_nsim500_n200_p10_a2.RData"); p = 10
r200p10a2 <- post(FCOL,form,p,T,1)
load("Ex2_nsim500_n200_p10_a3.RData"); p = 10
r200p10a3 <- post(FCOL,form,p,T,1)

load("Ex2_nsim500_n200_p30_a1.RData"); p = 30
r200p30a1 <- post(FCOL,form,p,T,1)
load("Ex2_nsim500_n200_p30_a2.RData"); p = 30
r200p30a2 <- post(FCOL,form,p,T,1)
load("Ex2_nsim500_n200_p30_a3.RData"); p = 30
r200p30a3 <- post(FCOL,form,p,T,1)

# n = 500; p = 10
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("Ex2_nsim500_n500_p10_a1.RData"); p = 10
r500p10a1 <- post(FCOL,form,p,T,1)
load("Ex2_nsim500_n500_p10_a2.RData"); p = 10
r500p10a2 <- post(FCOL,form,p,T,1)
load("Ex2_nsim500_n500_p10_a3.RData"); p = 10
r500p10a3 <- post(FCOL,form,p,T,1)

load("Ex2_nsim500_n500_p30_a1.RData"); p = 30
r500p30a1 <- post(FCOL,form,p,T,1)
load("Ex2_nsim500_n500_p30_a2.RData"); p = 30
r500p30a2 <- post(FCOL,form,p,T,1)
load("Ex2_nsim500_n500_p30_a3.RData"); p = 30
r500p30a3 <- post(FCOL,form,p,T,1)

# n = 1000; p = 10
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("Ex2_nsim500_n1000_p10_a1.RData"); p = 10
r1000p10a1 <- post(FCOL,form,p,T,1)
load("Ex2_nsim500_n1000_p10_a2.RData"); p = 10
r1000p10a2 <- post(FCOL,form,p,T,1)
load("Ex2_nsim500_n1000_p10_a3.RData"); p = 10
r1000p10a3 <- post(FCOL,form,p,T,1)

load("Ex2_nsim500_n1000_p30_a1.RData"); p = 30
r1000p30a1 <- post(FCOL,form,p,T,1)
load("Ex2_nsim500_n1000_p30_a2.RData"); p = 30
r1000p30a2 <- post(FCOL,form,p,T,1)
load("Ex2_nsim500_n1000_p30_a3.RData"); p = 30
r1000p30a3 <- post(FCOL,form,p,T,1)

round(rbind(cbind(r200p10a1$res,r200p10a1$cr),cbind(r500p10a1$res,r500p10a1$cr), cbind(r1000p10a1$res,r1000p10a1$cr)),4)
round(rbind(cbind(r200p10a2$res,r200p10a2$cr),cbind(r500p10a2$res,r500p10a2$cr), cbind(r1000p10a2$res,r1000p10a2$cr)),4)
round(rbind(cbind(r200p10a3$res,r200p10a3$cr),cbind(r500p10a3$res,r500p10a3$cr), cbind(r1000p10a3$res,r1000p10a3$cr)),4)

round(rbind(cbind(r200p30a1$res,r200p30a1$cr),cbind(r500p30a1$res,r500p30a1$cr), cbind(r1000p30a1$res,r1000p30a1$cr)),4)
round(rbind(cbind(r200p30a2$res,r200p30a2$cr),cbind(r500p30a2$res,r500p30a2$cr), cbind(r1000p30a2$res,r1000p30a2$cr)),4)
round(rbind(cbind(r200p30a3$res,r200p30a3$cr),cbind(r500p30a3$res,r500p30a3$cr), cbind(r1000p30a3$res,r1000p30a3$cr)),4)


rbind(r200p30a2$res,r500p30a2$res, r1000p30a2$res)
rbind(r200p30a2$res,r500p30a2$res, r1000p30a2$res)
rbind(r200p30a3$res,r500p30a3$res, r1000p30a3$res)


# write.table(rbind(r200p10a1$res,r500p10a1$res, r1000p10a1$res,r200p30a1$res,r500p30a1$res, r1000p30a1$res), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
# write.table(rbind(r200p10a2$res,r500p10a2$res, r1000p10a2$res,r200p30a2$res,r500p30a2$res, r1000p30a2$res), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
# write.table(rbind(r200p10a3$res,r500p10a3$res, r1000p10a3$res,r200p30a3$res,r500p30a3$res, r1000p30a3$res), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
