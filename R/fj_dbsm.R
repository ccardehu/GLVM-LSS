
rm(list = ls())
source("f0_prep.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 1: Heteroscedastic Normal model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n = c(200,500,1000)     # Number of individuals
p = c(10,30)      # Number of items
nsim = 500  # Number of simulations

s.form <- list("mu" = "~ Z1+Z2", "sigma" = "~ Z1+Z2")
irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0), c("sigma",1,"Z2",0), c("sigma",2,"Z1",0))

for(p. in p){
 fam <- rep("normal",p.)
 l1 <- lc <- NULL
 l1$mu <- lc$mu <- matrix(1,ncol = 3, nrow = p.)
 l1$sigma <- lc$sigma <- matrix(1, ncol = 3, nrow = p.)
 lc$mu[,1] <- runif(p.,1,2)
 lc$mu[,c(2,3)] <- runif(length(l1$mu[,c(2,3)]),0.5,1.5)# * sample(c(-1,1),size = 2*p,replace = T)
 if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3]
 lc$sigma <- matrix(runif(length(l1$sigma),0.1,0.4), nrow = p.)
 lc$sigma[,c(2,3)] <- lc$sigma[,c(2,3)]# * sample(c(-1,1),size = 2*p,replace = T)
 if(ncol(lc$sigma) >= 2 && lc$sigma[1,2] < 0) lc$sigma[1,2] <- -lc$sigma[1,2]
 if(ncol(lc$sigma) >= 2 && lc$sigma[2,3] < 0) lc$sigma[2,3] <- -lc$sigma[2,3]
 
 # Parameter restrictions:
 # ~~~~~~~~~~~~~~~~~~~~~~~
 l1. <- l1
 l1.$mu[1,c(3)] <- l1.$mu[2,c(2)] <- 0
 l1.$mu[sample(3:p., floor((p.-2)/2), replace = F),2] <- 0
 for(i in 3:nrow(l1.$mu)){ if(l1.$mu[i,2] != 0) l1.$mu[i,3] <- rbinom(1,1,0.5) }
 l1.$sigma[1,3] <- l1.$sigma[2,2] <- 0
 l1.$sigma[sample(3:p., floor((p.-2)/2), replace = F),2] <- 0
 for(i in 3:nrow(l1.$sigma)){ if(l1.$sigma[i,2] == 0) l1.$sigma[i,3] <- rbinom(1,1,0.5) }
  
 lc. <- NULL; for(i in names(l1.)){ lc.[[i]] <- lc[[i]]*l1.[[i]] }
   
 for(n. in n){
  sdb <- array(NA,dim = c(n.,p.,nsim))
  for(i in 1:nsim){
   sdb[,,i] <- as.matrix(splvm.sim(n = n., form = s.form, fam = fam, constraints = l1., coefs = lc.)$Y)
  }
  save.image(paste0("DB1_n",n.,"_p",p.,".RData"))
 }
}

# ~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 2: IRT model
# ~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())
source("f0_prep.R")

n = c(200,500,1000)     # Number of individuals
p = c(10,30)      # Number of items
nsim = 500  # Number of simulations

s.form <- list("mu" = "~ Z1+Z2")
irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0)) 

for(p. in p){
 fam <- rep("binomial",p.)
 l1 <- lc <- NULL
 l1 <- lc <- NULL
 l1$mu <- lc$mu <- matrix(1,ncol = 3, nrow = p.)
 lc$mu[,1] <- runif(length(l1$mu[,1]),-1,1)
 lc$mu[,c(2,3)] <- runif(length(l1$mu[,c(2,3)]),1.5,2.5) * sample(c(-1,1),size = 2*p.,replace = T)
 if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3]
 
 # Parameter restrictions:
 # ~~~~~~~~~~~~~~~~~~~~~~~
 l1. <- l1
 l1.$mu[1,c(3)] <- l1.$mu[2,c(2)] <- 0
 l1.$mu[sample(3:p., floor((p.-2)/2), replace = F),2] <- 0
 for(i in 3:nrow(l1.$mu)){ if(l1.$mu[i,2] != 0) l1.$mu[i,3] <- rbinom(1,1,0.5) }

 lc. <- NULL; for(i in names(l1.)){ lc.[[i]] <- lc[[i]]*l1.[[i]] }
   
 for(n. in n){
  sdb <- array(NA,dim = c(n.,p.,nsim))
  for(i in 1:nsim){
   sdb[,,i] <- as.matrix(splvm.sim(n = n., form = s.form, fam = fam, constraints = l1., coefs = lc.)$Y)
  }
  save.image(paste0("DB2_n",n.,"_p",p.,".RData"))
 }
}

