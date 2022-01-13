
###############################################
## SPLVM: Auxiliary functions for estimation ##
###############################################

probs <- function(x){
 pr <- plogis(x)
 if (any(pr == 1)) pr[pr == 1] <- 1 - sqrt(.Machine$double.eps)
 if (any(pr == 0)) pr[pr == 0] <- sqrt(.Machine$double.eps)
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

loglkf <- function(cb,Y,ghQ,bg,fam,info,res,full){

# Goal: To compute log-likelihood function (to be used in trust function)
# Input : cb (non-zero loadings), Y (matrix of items), ghQ (GHQ object),
#         bg (guide matrix), fam (distributions)
# Output: List with log-likelihood, gradient & full Hessian
# Testing: cb = lb2cb(bold); Y = simR$Y; ghQ = gr; bg = bold; fam = fam; res = loadmt2

b1 <- lb2cb(bg)
b1[t(lb2mb(res))] <- cb # complete with cb including zeros
b1 <- cb2lb(b1,bg)
A1 <- dY(Y,ghQ,b1,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
ll <- sum(log(A2)) # log-likelihood old
pD <- exp(rowSums(A1,dim = 2))/A2 # this is EC (posterior density)
dvL <- dvY(Y,ghQ,b1,fam,info) # List with all the necessary derivatives
A3 <- sche(ghQ,b1,res,fam,dvL,pD,info,full) # score & Hessian object
return(list(value = ll, gradient = A3$gradient, hessian = A3$hessian))
}

penloglkf <- function(cb,Y,ghQ,bg,fam,pen.idx,info,res,full,
                      pml.control = list(type = "lasso", lambda = 1, w.alasso = NULL,a = NULL)) {

# Goal: To compute penalised log-likelihood function (to be used in trust function)
# Input : cb (non-zero loadings), Y (matrix of items), ghQ (GHQ object),
#         bg (guide matrix), fam (distributions), pen.idx (penalty index / T if penalty)
# Output: List with log-likelihood, gradient & full Hessian
# Testing: cb = lb2cb(bold); Y = Y; ghQ = ghQ; bg = bold; fam = fam; 
#          pml.control = list(type = "lasso", lambda = 1, w.alasso = NULL, gamma = 1, a = 3.7)

res <- uplm(bg,res)
b1 <- lb2cb(bg)
b1[t(lb2mb(res))] <- cb # complete with cb including zeros
b1 <- cb2lb(b1,bg)
A1 <- dY(Y,ghQ,b1,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
A2a <- lb2pM(b1,Y,pen.idx,res,pml.control)
ll <- c(sum(log(A2)) - A2a$lp) # penalized log-likelihood old
pD <- exp(rowSums(A1,dim = 2))/A2 # this is EC (posterior density)
dvL <- dvY(Y,ghQ,b1,fam,info) # List with all the necessary derivatives
A3 <- sche(ghQ,b1,res,fam,dvL,pD,info,full) # score & Hessian object
return(list(value = c(ll), gradient = A3$gradient - c(tcrossprod(cb,A2a$pM)), hessian = A3$hessian - A2a$pM))
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

dFun <- function(i,Y,Z.,b,fam){
  
# Goal: To compute (simulated) log-likelihood
# Input i (item), Y (original data), Z. (Z output in ff_msim), b (betas), fam (family)
# Output: Log-likelihood
# Testing: i = 1; Y = simR$Y; Z = simR$Zout; b = borg; fam = fam

for(a in names(Z.)){ if(!is.matrix(Z.[[a]])) Z.[[a]] <- as.matrix(Z.[[a]]) }
if(fam[i] == "normal"){ llg <- sum(dnorm(Y[,i], c(Z.$mu%*%matrix(b$mu[i,])), c(exp(Z.$sigma%*%matrix(b$sigma[i,]))),log = T)) }
if(fam[i] == "lognormal"){ llg <- sum(dlnorm(Y[,i], c(Z.$mu%*%matrix(b$mu[i,])), c(exp(Z.$sigma%*%matrix(b$sigma[,i]))),log = T)) }
if(fam[i] == "poisson"){ llg <- sum(dpois(Y[,i], c(exp(Z.$mu%*%matrix(b$mu[i,]))), log = T)) }
if(fam[i] == "gamma"){ llg <- sum(dgamma(Y[,i], shape = c(exp(Z.$mu%*%matrix(b$mu[i,]))), scale = c(exp(Z.$sigma%*%matrix(b$sigma[i,]))),log = T)) }
if(fam[i] == "binomial"){ llg <- sum(dbinom(Y[,i],1, c(probs(Z.$mu%*%matrix(b$mu[i,]))), log = T)) }
if(fam[i] == "ZIpoisson"){
 dZIpoisson <- function(Y,mu,sigma,log = T){
 u <- as.numeric(Y == 0)
 lf <- u*c(log(sigma + (1-sigma)*exp(-mu))) + (1-u)*(c(log(1-sigma)) - c(mu) + Y*c(log(mu)) - lfactorial(Y))
 if(log == T) return(lf) else return(exp(lf)) }
 llg <- sum(dZIpoisson(Y[,i], c(exp(Z.$mu%*%matrix(b$mu[i,]))), c(probs(Z.$sigma%*%matrix(b$sigma[i,]))), log = T)) }
return(llg)
}

llkf <- function(cb,Y,ghQ,bg,fam){

# Goal: To compute log-likelihood function
# Input : cb (non-zero loadings), Y (matrix of items), ghQ (GHQ object),
#         bg (guide matrix), fam (distributions)
# Output: List with log-likelihood, gradient & full Hessian
# Testing: cb = lb2cb(borg); Y = simR$Y; ghQ = ghQ; bg = borg; fam = fam; 

b1 <- lb2cb(bg)
b1[lb2cb(bg) != 0] <- cb # complete with cb including zeros
b1 <- cb2lb(b1,bg)
A1 <- dY(Y,ghQ,b1,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
ll <- sum(log(A2)) # log-likelihood
return(ll)
}

sche.test <- function(cb,mod){

# Goal: To compute scores & hessian from estimated model
# Input : vector of betas (cb), mod (splvm.fit object)
# Output: List with gradient (score) & full Hessian computed @ cb2lb(cb)
# Testing: cb = lb2cb(borg); mod = testUP

b <- cb2lb(cb,mod$b)
A1 <- dY(mod$Y,mod$ghQ,b,mod$fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%mod$ghQ$weights) # this is efy
mod.pD <- exp(rowSums(A1,dim = 2))/A2 # this is EC (posterior density)
mod.dvL <- dvY(mod$Y,mod$ghQ,b,mod$fam,mod$info) # List with all the necessary derivatives
res <- sche(ghQ = mod$ghQ, b = b, rs = mod$loadmt, fam = mod$fam, dvL = mod.dvL, pD = mod.pD, info = mod$info, full.hess = T)    
return(res)
}

m2pdm <- function(mat){
   
# Goal: checks & corrects if matrix object (mat) is Positive Definite
# Input : matrix (mat)
# Output: corrected (if necessary) positive definite matrix
# Testing: i = 1; ghQ = ghQ; b = l1; fam = fam; dvL = dvY(Y,ghQ,b,fam);
#          pD = exp(rowSums(dY(Y,ghQ,b,fam),dim = 2)) /
#          c(exp(rowSums(dY(Y,ghQ,b,fam),dim = 2))%*%ghQ$weights)
#          mat = hess(i,i,ghQ,b,fam,dvL,pD,"Fisher")$hessian

eS    <- eigen(mat, symmetric = TRUE)
e.val <- eS$values
e.vec <- eS$vectors
check.eigen <- any(e.val <= 0)
if(check.eigen == TRUE){
n.e.val <- e.val[e.val <= 0]
s <- sum(e.val[n.e.val])*2  
t <- s^2*100 + 1
p <- min(e.val[(e.val <= 0) == FALSE])
e.val[e.val <= 0] <- p*(s - n.e.val)^2/t
D <- diag(e.val)
D.inv <- diag(1/e.val)
res <- e.vec %*% D %*% t(e.vec) 
res.inv <- e.vec %*% D.inv %*% t(e.vec)
res.sqr <- e.vec %*% sqrt(D)
} else { res <- mat; res.inv <- e.vec %*% diag(1/e.val) %*% t(e.vec) ; res.sqr <- e.vec %*% diag(sqrt(e.val))} 
res.inv <- (res.inv + t(res.inv) ) / 2 
return(list(mat = res, inv.mat = res.inv, sqr.mat = t(res.sqr), is.PDM = !check.eigen))
}

lb2pM <- function(lb,Y,pen.idx,loadmt,
                  pml.control = list(type = "lasso", lambda = 1, w.alasso = NULL,a = NULL)){

# Goal: To compute the penalty part for the log-likelihood
# Input : lb (list betas), Y (matrix of items), pml.control (options for Penalty)
# Output: penalty for log-likelihood
# Testing: lb = bold; Y = simR$Y; pml.control = list(type = "lasso", lambda = 1, w.alasso = NULL, a = NULL)
   
b <- lb2cb(lb)[t(lb2mb(loadmt))]
P <- nrow(Y)*penM(lb2cb(lb),type = pml.control$type, id = pen.idx, rs = loadmt, lambda = pml.control$lambda,
                  w.alasso = pml.control$w.alasso,a = pml.control$a)$full
return(list(lp = 0.5*crossprod(b,P)%*%b, pM = P))
}

pidx <- function(lb,loadmt,pen.load){

# Goal: Index of penalised columns of betas 
# Input : lb (loadings list), loadmt (free parameters == T) pen.load (penalise loadings?)
# Output: Vector of indices (for columns to be used in lb2mb)
# Testing: lb = borg; loadmt; pen.load = F

idx <- NULL
for(i in names(lb)){
 tmp <- array(T,dim = dim(lb[[i]]), dimnames = list(NULL,colnames(lb[[i]])))
 if("(Intercept)" %in% colnames(tmp)) { tmp[,"(Intercept)"] <- F }
 tmp <- tmp & loadmt[[i]] & (lb[[i]] != 0); 
 if(i == "mu" && !pen.load){ tmp[,c(colnames(tmp) %in% c("(Intercept)","Z1","Z2","Z3"))] <- F }
 idx <- cbind(idx, tmp)
}
return(unname(idx))
}

ini.par <- function(Y,fam,form,pC,q.){

# Goal: To produce initial parameters for splvm
# Input : Y (matrix of items), fam (distributions),  form (formulas for mu, sigma, etc.)
#         pC (list of whether item i has mu or sigma, etc.), q. (# latent variables)
# Output: List of parameters
# Testing: Y = simR$Y; fam = fam; form = form; pC = pC; q. = q. # (get form fc_mfit first)

# Z. <- princomp(Y,cor = T)$scores[,q.:1,drop=F]%*%diag(c(rep(-1,q.-1),1))

Z. <- princomp(Y,cor = T)$scores[,1:q.,drop=F]%*%diag(c(-1,rep(ifelse(all(fam == "normal"),1,-1),q.-1)),nrow = length(c(-1,rep(ifelse(all(fam == "normal"),1,-1),q.-1))))
colnames(Z.) <- paste0("Z", 1:q.)
sZ <- NULL
bstart <- NULL
for(r in names(form)){
 sZ[[r]] <- as.data.frame(model.matrix(form[[r]],as.data.frame(Z.)))
 bstart[[r]] <- matrix(0, nrow = ncol(Y), ncol = ncol(sZ[[r]]))
 dimnames(bstart[[r]]) <- list(colnames(Y),colnames(sZ[[r]]))
}

if(!is.matrix(Y)) Y <- as.matrix(Y)
for(i in 1:ncol(Y)){
 eq <- NULL
 tmpY <- Y[,i]
 eq$mu <- update(form$mu, tmpY ~ .)
 if(i %in% pC$sigma) eq$sigma <- update(form$sigma, tmpY ~ .)
 if(i %in% pC$tau) eq$tau <- update(form$tau, tmpY ~ .)
 if(i %in% pC$nu) eq$nu <- update(form$nu, tmpY ~ .)
 if(fam[i] == "normal"){
  tmp <- gamlss(eq$mu,sigma.fo = eq$sigma,data = as.data.frame(cbind(tmpY,Z.)), control = gamlss.control(trace = F))
  bstart$mu[i,] <- coef(tmp,"mu")
  bstart$sigma[i,] <- coef(tmp,"sigma")
 }
 if(fam[i] == "lognormal"){ 
  tmp <- gamlss(eq$mu,sigma.fo = eq$sigma,family = LOGNO(),data = as.data.frame(cbind(tmpY,Z.)), control = gamlss.control(trace = F))
  bstart$mu[i,] <- coef(tmp,"mu")
  bstart$sigma[i,] <- coef(tmp,"sigma")
 }
 if(fam[i] == "poisson"){ 
  tmp <- gamlss(eq$mu,family = PO(),data = as.data.frame(cbind(tmpY,Z.)), control = gamlss.control(trace = F))
  bstart$mu[i,] <- coef(tmp,"mu")
 }
 if(fam[i] == "gamma"){ 
  tmp <- gamlss(eq$mu,sigma.fo = eq$sigma,family = GA(),data = as.data.frame(cbind(tmpY,Z.)), control = gamlss.control(trace = F))
  bstart$mu[i,] <- coef(tmp,"mu")
  bstart$sigma[i,] <- coef(tmp,"sigma")
 }
 if(fam[i] == "binomial"){ 
  tmp <- gamlss(eq$mu,family = BI(),data = as.data.frame(cbind(tmpY,Z.)), control = gamlss.control(trace = F))
  bstart$mu[i,] <- coef(tmp,"mu")
 }
 if(fam[i] == "ZIpoisson"){ 
  tmp <- gamlss(eq$mu,sigma.fo = eq$sigma,family = ZIP(),data = as.data.frame(cbind(tmpY,Z.)), control = gamlss.control(trace = F))
  bstart$mu[i,] <- coef(tmp,"mu")
  bstart$sigma[i,] <- coef(tmp,"sigma")
 }
}
for(r in names(form)){ if("Z1" %in% colnames(bstart[[r]]) && bstart[[r]][1,"Z1"] < 0){ bstart[[r]][,"Z1"] <- -bstart[[r]][,"Z1"] } }
return(bstart)
}

GBIC <- function(mod){

if(mod$loglik == -999){ return(-999) } else {
ll <- mod$uploglik
Hm <- mod$hessian$unp
iHm <- tryCatch({solve(mod$hessian$unp)}, error = function(e){-m2pdm(-mod$hessian$unp)$inv})

if(!is.null(mod$pml.control)){
ipS <- tryCatch({solve(mod$hessian$pen)}, error = function(e){-m2pdm(-mod$hessian$pen)$inv})
GBIC <- -2*ll + log(nrow(mod$Y))*sum(diag(ipS%*%Hm)) } else {
  GBIC <- -2*ll + log(nrow(mod$Y))*sum(diag(iHm%*%Hm))
}
return(GBIC) }
}

GIC <- function(mod){

if(mod$loglik == -999){ return(-999) } else {
  
ll <- mod$uploglik
Hm <- mod$hessian$unp
iHm <- tryCatch({solve(mod$hessian$unp)}, error = function(e){-m2pdm(-mod$hessian$unp)$inv})

if(!is.null(mod$pml.control)){
ipS <- tryCatch({solve(mod$hessian$pen)}, error = function(e){-m2pdm(-mod$hessian$pen)$inv})
GIC <- -2*ll + 2*sum(diag(ipS%*%Hm)) } else {
  GIC <- -2*ll + 2*sum(diag(iHm%*%Hm))
}
return(GIC) }
}

fscore <- function(mod){
   
# Goal: To compute Expected a-posterior Scores for model 'mod'
# Input : mod (splvm model)
# Output: Matrix of factor scores
# Testing: mod = testa
 
A1 <- dY(mod$Y,mod$ghQ,mod$b,mod$fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%mod$ghQ$weights) 
pD <- exp(rowSums(A1,dim = 2))/A2 # posterior density
zscore <- pD%*%(c(mod$ghQ$weights)*mod$ghQ$points)[,,drop=F]
zscore <- as.data.frame(zscore); names(zscore) <- colnames(mod$ghQ$points)  
return(zscore)
}

rmat <- function(res,b){

# Goal: To build a restriction matrix 
# Input: res (dataframe), b (betas matrix lb2mb form)
# Output: matrix of restrictions (l2bmb form)
# Testing: res = l2rm(list(c("sigma",2,"Z1",0),c("mu", 1,"(Intercept)",1))); b = bold
 
rest <- l2rm(res)
frparm <- vector(mode = "list",length = length(b))
names(frparm) <- names(b)
for(i in names(b)){
 frparm[[i]] <- array(NA, dim = dim(b[[i]]))
 dimnames(frparm[[i]]) <- dimnames(b[[i]])
}
for(i in 1:nrow(rest)){
  frparm[[as.character(rest[i,1])]][as.integer(rest[i,2]),as.character(rest[i,3])] <- as.numeric(rest[i,4])
}
return(frparm)
}

l2rm <- function(rlist){

# Goal: List to Restrictions matrix
# Input: List of restrictions
# Output: Dataframe of restrictions
# Testing: rlist = list(c("sigma",2,"Z1",0),c("mu", 1,"(Intercept)",1))

if(!is.list(rlist)) stop("Restrictions should be a list with element(s), each of the type c('parameter',item,'restricted variable',value)")
r <- as.data.frame(matrix(unlist(rlist), ncol=4, byrow = T))
return(r)
}
