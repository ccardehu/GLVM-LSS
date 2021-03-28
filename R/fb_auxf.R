
###############################################
## SPLVM: Auxiliary functions for estimation ##
###############################################

probs <- function(x){
 pr <- plogis(x)
 if (any(ind <- pr == 1)) pr[ind] <- 1 - sqrt(.Machine$double.eps)
 if (any(ind <- pr == 0)) pr[ind] <- sqrt(.Machine$double.eps)
 return(pr)
}

expit <- function(x){
   exp(x)/(1+exp(x)) }

logit <- function(x){
   log(x/(1-x)) }

lb2mb <- function(b){

# Goal: List betas to Matrix (cbind) betas
# Input : b (loadings matrix)
# Output: Matrix of pxQ (Q = # total parameters for item i)
# Testing: b = borg

l2m <- NULL
for(i in 1:length(b)){ l2m <- cbind(l2m,b[[i]]) }
return(unname(l2m))
}

lb2cb <- function(b){
   
# Goal: List betas to vector betas (by item)
# Input : b (loadings matrix)
# Output: vector of lenght R (R = # total parameters in the model)
# Testing: b = borg
   
return(c(t(lb2mb(b))))
}

cb2lb <- function(cb,b){

# Goal: vector betas to list betas
# Input : cb (loadings), b (list of loadings)
# Output: List of loadings by parameter
# Testing: b = borg; cb = lb2cb(b)

bm <- matrix(cb,byrow=T,nrow=nrow(b$mu))
c2l <- vector(mode = "list", length = length(b))
names(c2l) <- names(b)
for(i in seq_along(names(b))){
 if(i == 1) c2l[[i]] <- bm[,1:ncol(b[[i]]), drop = F]
 else c2l[[i]] <- bm[,seq(from = (ncol(b[[i-1]])+1), length = ncol(b[[i]])), drop = F]
 colnames(c2l[[i]]) <- colnames(b[[i]]); rownames(c2l[[i]]) <- rownames(b[[i]])
}
return(c2l)
}

mb2lb <- function(mb,b){

# Goal: matrix betas to list betas
# Input : mb (matrix loadings cbind by parameter), b (list of loadings)
# Output: List of loadings by parameter
# Testing: b = borg; mb = lb2mb(b)

m2l <- vector(mode = "list", length = length(b))
names(m2l) <- names(b)
for(i in seq_along(names(b))){
 if(i == 1) m2l[[i]] <- mb[,1:ncol(b[[i]]), drop = F]
 else m2l[[i]] <- mb[,seq(from = (ncol(b[[i-1]])+1), length = ncol(b[[i]])), drop = F]
 colnames(m2l[[i]]) <- colnames(b[[i]]); rownames(m2l[[i]]) <- rownames(b[[i]])
}
return(m2l)
}

pFun <- function(f){ # eval @ fam[1]
 if(f == "normal") pars <- c("mu","sigma")
 if(f == "lognormal") pars <- c("mu","sigma")
 if(f == "poisson") pars <- c("mu")
 if(f == "gamma") pars <- c("mu", "sigma")
 if(f == "binomial") pars <- c("mu")
 if(f == "ZIpoisson") pars <- c("mu", "sigma")
 # Add other parameters for other distributions
 return(pars)
}

loglkf <- function(cb,Y,ghQ,bg,fam){

# Goal: To compute log-likelihood function (to be used in trust function)
# Input : cb (non-zero loadings), Y (matrix of items), ghQ (GHQ object),
#         bg (guide matrix), fam (distributions)
# Output: List with log-likelihood, gradient & full Hessian
# Testing: cb = lb2cb(borg); Y = simR$Y; ghQ = gr; bg = borg; fam = fam; 

b1 <- lb2cb(bg)
b1[lb2cb(bg) != 0] <- cb # complete with cb including zeros
b1 <- cb2lb(b1,bg)
A1 <- dY(Y,ghQ,b1,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
pD1 <- exp(rowSums(A1,dim = 2))/A2
ll <- sum(log(A2)) # log-likelihood
dvL1 <- dvY(Y,ghQ,b1,fam)
res <- sche(ghQ,b1,fam,dvL1,pD1)
return(list(value = ll, gradient = res$gradient, hessian = res$hessian))
}

rFun <- function(n,i,fam,Z,b){ #
if(fam == "normal"){
mu = as.matrix(Z$mu)%*%matrix(b$mu[i,])
sigma = exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))
fyz <- rnorm(n, c(mu), c(sigma))
}
if(fam == "lognormal"){
mu = as.matrix(Z$mu)%*%matrix(b$mu[i,])
sigma = exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))
fyz <- rlnorm(n, drop(mu), drop(sigma))
}
if(fam == "poisson"){
mu = exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
fyz <- rpois(n, c(mu)) 
}
if(fam == "gamma"){
mu = exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
sigma = exp(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))
fyz <- rgamma(n,shape = mu, scale = sigma)
}
if(fam == "binomial"){
mu = probs(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
fyz <- rbinom(n,1,prob = mu)
}  
if(fam == "ZIpoisson"){
mu = exp(as.matrix(Z$mu)%*%matrix(b$mu[i,]))
sigma = probs(as.matrix(Z$sigma)%*%matrix(b$sigma[i,]))
fyz <- ifelse(rbinom(n,1,1-sigma) == 0, 0, rpois(n,mu))
}
# Add other distributions in SimFA::fod
return(fyz)
}
