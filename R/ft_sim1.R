
# Code
# File ft_sim1.R: Simulation 1 (Chapter 2)
# Code written by: Camilo Cardenas-Hurtado (c.a.cardenas-hurtado@lse.ac.uk)

source("f0_prep.R")

n = 1000     # Number of individuals
nsim = 1000  # Number of simulations

# Simulation 1: Heteroscedastic Normal model with interactions:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s.form <- list("mu" = "~ Z1*Z2", "sigma" = "~ Z1"); p = 15
fam <- rep("normal",p)
l1 <- lc <- NULL
l1$mu <- matrix(1,ncol = 4, nrow = p)
l1$sigma <- matrix(1, ncol = 2, nrow = p)
lc$mu <- matrix(runif(length(l1$mu),-1.5,1.5),nrow = p)
lc$sigma <- matrix(runif(length(l1$sigma),0.1,0.3), nrow = p)
lc$mu[,1] <- runif(p,1,2); if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]
lc$mu[abs(lc$mu) < 0.5] <- runif(sum(abs(lc$mu) < 0.5),0.5,1.5)
if(ncol(lc$sigma) >= 2 && lc$sigma[1,2] < 0) lc$sigma[1,2] <- -lc$sigma[1,2]

# ~~~~~~~~~~~
# Starts here:
# ~~~~~~~~~~~

l1. <- l1
# Parameter restrictions:
# Simulation 1: Heteroscedastic Normal model with interactions:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
l1.$mu[1:5,3:4] <- 0; l1.$mu[6:10,4] <- 0
l1.$sigma[1:10,2] <- 0
lc. <- NULL; for(i in names(l1.)){ lc.[[i]] <- lc[[i]]*l1.[[i]] }

# Restrictions (l1.)
# ~~~~~~~~~~~~~~~~~~ 

# Simulation 1: Heteroscedastic Normal model with interactions:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
restr <- list(c("mu",1,"Z2",0),c("mu",2,"Z2",0),c("mu",3,"Z2",0),
              c("mu",4,"Z2",0),c("mu",5,"Z2",0),
              c("mu",1,"Z1:Z2",0),c("mu",2,"Z1:Z2",0),c("mu",3,"Z1:Z2",0),
              c("mu",4,"Z1:Z2",0),c("mu",5,"Z1:Z2",0),
              c("mu",6,"Z1:Z2",0),c("mu",7,"Z1:Z2",0),c("mu",8,"Z1:Z2",0),
              c("mu",9,"Z1:Z2",0),c("mu",10,"Z1:Z2",0),
              c("sigma",1,"Z1",0),c("sigma",2,"Z1",0),c("sigma",3,"Z1",0),
              c("sigma",4,"Z1",0),c("sigma",5,"Z1",0),c("sigma",6,"Z1",0),
              c("sigma",7,"Z1",0),c("sigma",8,"Z1",0),c("sigma",9,"Z1",0),c("sigma",10,"Z1",0))

cores <- detectCores()-2
cl<-makeCluster(cores)
registerDoSNOW(cl)
progress <- function(n) setTxtProgressBar(txtProgressBar(max = nsim, style = 3), n)
opts <- list(progress = progress)
saveRes = T

FCOL <- foreach(l = 1:5,
                .combine = rbind,
                .packages = c("mvtnorm","fastGHQuad","gamlss","trust"),#,
                .options.snow = opts
                ) %dopar% splvm.simfit(l,4)

stopCluster(cl)
# if(cleanRes) FCOL <- FCOL[FCOL[,ncol(FCOL)] %!in% c(1000,-999),-ncol(FCOL)]
if(saveRes) save(FCOL, file = paste0("nsim",nsim,"_n",n, "_p", p,"_Ex1.RData"))
