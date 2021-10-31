
# Code
# File ft_sim1.R: Simulation 1 (Chapter 2)
# Code written by: Camilo Cardenas-Hurtado (c.a.cardenas-hurtado@lse.ac.uk)

rm(list = ls())
source("f0_prep.R")

n = 1000     # Number of individuals
p = 30      # Number of items
nsim = 500  # Number of simulations

# Simulation 2: IRT model:
# ~~~~~~~~~~~~~~~~~~~~~~~~
s.form <- list("mu" = "~ Z1+Z2")
fam <- rep("binomial",p)
l1 <- lc <- NULL
l1$mu <- lc$mu <- matrix(1,ncol = 3, nrow = p)
lc$mu[,1] <- runif(length(l1$mu[,1]),-1,1)
lc$mu[,c(2,3)] <- runif(length(l1$mu[,c(2,3)]),1.5,2.5) * sample(c(-1,1),size = 2*p,replace = T)
if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3]

# ~~~~~~~~~~~
# Starts here:
# ~~~~~~~~~~~

l1. <- l1
# Parameter restrictions:
# Simulation 2
# ~~~~~~~~~~~~
l1.$mu[1,c(3)] <- l1.$mu[2,c(2)] <- 0
l1.$mu[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
for(i in 3:nrow(l1.$mu)){ if(l1.$mu[i,2] != 0) l1.$mu[i,3] <- rbinom(1,1,0.5) }

lc. <- NULL; for(i in names(l1.)){ lc.[[i]] <- lc[[i]]*l1.[[i]] }

# Restrictions (l1.)
# ~~~~~~~~~~~~~~~~~~ 

# Simulation 1: Heteroscedastic NLFM:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# restr <- list(c("mu",2,"Z1",0), c("mu",3,"Z1",0), c("mu",4,"Z1",0), c("mu",5,"Z1",0), c("mu",10,"Z1",0),
#               c("mu",1,"Z2",0), c("mu",8,"Z2",0), c("mu",9,"Z2",0),
#               c("sigma",2,"Z1",0), c("sigma",4,"Z1",0), c("sigma",6,"Z1",0), c("sigma",8,"Z1",0), c("sigma",10,"Z1",0),
#               c("sigma",1,"Z2",0), c("sigma",10,"Z2",0))

irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0))

cores <- detectCores()-2
cl<-makeCluster(cores)
registerDoSNOW(cl)
progress <- function(n) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), n)
opts <- list(progress = progress)
saveRes = T

FCOL <- foreach(l = 1:nsim,
                .combine = rbind,
                .packages = c("mvtnorm","fastGHQuad","gamlss","trust"),#,
                .options.snow = opts
                ) %dopar% splvm.simfit(l)

stopCluster(cl)
# if(cleanRes) FCOL <- FCOL[FCOL[,ncol(FCOL)] %!in% c(1000,-999),-ncol(FCOL)]
if(saveRes) save(FCOL, file = paste0("Ex2_","nsim",nsim,"_n",n, "_p", p,".RData"))
