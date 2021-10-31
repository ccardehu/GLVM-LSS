
# Code
# File ft_sim2.R: Simulation 2 (Chapter 2)
# Code written by: Camilo Cardenas-Hurtado (c.a.cardenas-hurtado@lse.ac.uk)

source("f0_prep.R")

n = 1000     # Number of individuals
nsim = 1000  # Number of simulations

# Simulation 2: Binomial model with interactions:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s.form <- list("mu" = "~ Z1+Z2"); p = 10
fam <- rep("binomial",p)
l1 <- lc <- NULL
l1$mu <- lc$mu <- matrix(1,ncol = 3, nrow = p)
lc$mu[,1] <- matrix(runif(length(l1$mu[,1]),-1,1),nrow = p)
lc$mu[,c(2,3)] <- matrix(runif(length(l1$mu[,c(2,3)]),1.5,2.5),nrow = p) * sample(c(-1,1),size = 2*p,replace = T)
if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]

# ~~~~~~~~~~~
# Starts here:
# ~~~~~~~~~~~

l1. <- l1
# Parameter restrictions:
# Simulation 2: Binomial
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
l1.$mu[sample(1:p, floor(p/2), replace = F),2] <- 0
for(i in 1:nrow(l1.$mu)){ if(l1.$mu[i,2] != 0) l1.$mu[i,3] <- rbinom(1,1,0.5) } #rbinom(1,1,0.5)
lc. <- NULL; for(i in names(l1.)){ lc.[[i]] <- lc[[i]]*l1.[[i]] }

# Restrictions (l1.)
# ~~~~~~~~~~~~~~~~~~ 

# Simulation 2: Binomial
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
restr <- list(c("mu",3,"Z1",0), c("mu",5,"Z1",0), c("mu",8,"Z1",0), c("mu",9,"Z1",0), c("mu",10,"Z1",0),
              c("mu",2,"Z2",0), c("mu",4,"Z2",0), c("mu",6,"Z2",0))

cores <- detectCores()
cl<-makeCluster(cores)
registerDoSNOW(cl)
progress <- function(n) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), n)
opts <- list(progress = progress)
saveRes = T

FCOL <- foreach(l = 1:nsim,
                .combine = rbind,
                .packages = c("mvtnorm", "fastGHQuad","gamlss","trust"),
                .options.snow = opts
                ) %dopar% splvm.simfit(l,3.5)

stopCluster(cl)
if(saveRes) save(FCOL, file = paste0("nsim",nsim,"_n",n, "_p", p,"_Ex2.RData"))
