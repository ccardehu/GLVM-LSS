

splvm.sim <- function(n,fam,form,constraints,coefs){

# Goal: Simulate semi-parametric LVM (using constraints, formulas and distributions)
# Input : n (sample size), fam (distributions),#         
#         form (list of formulas for measurement eqs. for each parameter),
#         constraints (matrix of coefs. constraints), coefs (loadings)
# Output: List of Y (simulated items), Z (simulated latent variables), b (original parameters),
#         constraints (restrictions), coefs (coefficients for simulation)
# Testing: n = 100; fam = fam; form = s.form; constraints = l1; coefs = lc

p <- length(fam)  
parY <- unique(unlist(lapply(1:length(fam),function(i) pFun(fam[i]))))
if(!is.list(form)) stop("Argument `form` should be a list with elements mu, sigma, tau and/or nu")
if(!all(names(form) == parY)) stop("Error in formula, check for mu, sigma, tau or nu")
for(i in parY){ form[[i]] <- as.formula(form[[i]]) }

lvar <- unique(unlist(lapply(1:length(form), function(i) all.vars(form[[i]]))))
lvar <- lvar[grep("Z", lvar, fixed = T)]
if(length(lvar) == 0) stop("Latent variables in formula should be represented by letter Z (capital)")
q. <- length(lvar)

# Latent variables
# ~~~~~~~~~~~~~~~~
muZ <- c(rep(0,q.))
sigZ <- diag(q.)
Z <- as.data.frame(rmvnorm(n, muZ, sigZ)); colnames(Z) <- paste0("Z", 1:q.)
Z. <- t. <- NULL
for(i in parY){
 Z.[[i]] <- as.data.frame(model.matrix(form[[i]],Z))
 t.[[i]] <- ncol(Z.[[i]]) # q. + 1 if intercept = T, or q. if intercept = F
}

# Loadings 
# ~~~~~~~~
if(!missing(coefs)){
 if(!is.list(coefs)) stop("Provided `coefs` should be a list with dim p*(q+intercept) loading matrices for mu, sigma, tau, nu")
 names(coefs) <- parY
 for(i in parY){
  if(length(coefs[[i]]) != p*t.[[i]]) stop("Provided `coefs` is not lenght p*(q+intercept), revise dimension")
  if(!is.matrix(coefs[[i]])) coefs[[i]] <- matrix(coefs[[i]], nrow = p, ncol = t.[[i]])
  colnames(coefs[[i]]) <- colnames(Z.[[i]]); rownames(coefs[[i]]) <- paste0("Y", 1:p)
 }
} else {
 warning("Coefficients for simulations not supplied, check code for assumed", call. = F)
 coefs <- NULL
 for(i in parY){
  tmp <- c(rep(0.4,p), rep(0.8,p*(t.[[i]]-1))) #tmp <- runif(n =  p*t.[[i]], 0.2,1.5)
  coefs[[i]] <- matrix(tmp, nrow = p, ncol = t.[[i]])
  colnames(coefs[[i]]) <- colnames(Z.[[i]]); rownames(coefs[[i]]) <- paste0("Y", 1:p)
  rm(tmp)
 }
}

# Restrictions on factor loadings
# _______________________________
if(!missing(constraints)){
 if(!is.list(constraints)) stop("Provided `constraints` should be a list with dim p*(q+intercept) restriction matrices for mu, sigma, tau, nu")
 names(constraints) <- parY
 if(length(constraints) != length(coefs)) stop("Number of matrices in `constraints` should match number of matrices in `coefs`")
 for(i in parY){
  if(length(constraints[[i]]) != p*t.[[i]]) stop("Provided `constraints` is not lenght p*(q+intercept), revise dimension")
  if(!is.matrix(constraints[[i]])) constraints[[i]] <- matrix(constraints[[i]], nrow = p, ncol = t.[[i]])
 }
} else {
 warning("Constraints matrix not supplied, simple structure assumed \n", call. = F)
 constraints <- NULL
 for(i in parY){
  if(length(all.vars(form[[i]])) > 1){
   r. <- p%/%(t.[[i]]-1)
   tmp1 <- matrix(0, nrow = nrow(coefs[[i]]), ncol = ncol(coefs[[i]]), dimnames = dimnames(coefs[[i]]))
   for(j in seq_along(all.vars(form[[i]]))){
    tmp2 <- seq((j-1)*(r.)+1,(j)*(r.),length.out = r.)
    tmp1[tmp2,all.vars(form[[i]])[j]] <- 1
    if(j == length(all.vars(form[[i]])) && p%%2 != 0) tmp1[(tail(tmp2,1)+1):nrow(tmp1),all.vars(form[[i]])[j]] <- 1
   }
   tmp1[,!(colnames(tmp1) %in% all.vars(form[[i]]))] <- 1
   constraints[[i]] <- tmp1; rm(tmp1,tmp2,r.)
  } else{
   constraints[[i]] <- matrix(1,nrow = p, ncol = t.[[i]])
   if(t.[[i]]-1 > 1){
    warning(paste0("Nonlinear functions for ", i, ", identification restriction impossed in item 1"), call. = F)
    constraints[[i]][1,t.[[i]]] <- 0
   }
  }
 }
}

borg <- NULL
for(i in parY){borg[[i]] <- coefs[[i]]*constraints[[i]]}
names(borg) <- parY

Y <- sapply(1:p, function(i) rFun(n,i,fam[i], Z.,borg))
Y <- as.data.frame(Y); colnames(Y) <- paste0("Y", 1:p)

return(list("Y" = Y, "Z" = Z, "Zout" = Z., "b" = borg, "formula" = form, "constraints" = constraints))
}
