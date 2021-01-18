
# set.seed(1234)

simGLVM <- function(n,p,form,dist,loadmt,coefs){

# To evaluate 
# n = 500; p = 10; form = form; dist = fam; 
# loadmt = l1; coefs = lc;
# _________________________________________________________
parY <- unlist(unique(lapply(1:length(dist),function(i) pFun(i,dist[i]))))

if(!is.list(form)) stop("Argument `form` should be a list with elements mu, sigma, tau or nu")
if(!all(names(form) == parY)) stop("Missing formula, check for mu, sigma, tau or nu")
for(i in parY){ form[[i]] <- as.formula(form[[i]]) }

lvar <- unique(unlist(lapply(1:length(form), function(i) all.vars(form[[i]]))))
lvar <- lvar[grep("Z", lvar, fixed = T)]
if(length(lvar) == 0) stop("Latent variables in formula should be represented by letter Z (capital)")
q. <- length(lvar)

# Latent variables (Z. for each mu, sigma, etc.)
# ______________________________________________
muZ <- c(rep(0,q.))
sigZ <- diag(q.)
Z <- as.data.frame(rmvnorm(n, muZ, sigZ)); colnames(Z) <- paste0("Z", 1:q.)
Z. <- t. <- f. <- NULL
for(i in parY){
 Z.[[i]] <- as.data.frame(model.matrix(form[[i]],Z))
 t.[[i]] <- ncol(Z.[[i]]) # q. + 1 if intercept = T, or q. if intercept = F
 f.[[i]] <- nrow(attributes(terms(form[[i]]))$factors) # Number of LVs in form[[i]]
 if(t.[[i]] != q. || sum(Z.[[i]][,1]) == n) colnames(Z.[[i]])[1] <- "Int"
}


# Loadings for mu and sigma (measurement equations)
# _________________________________________________
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
if(!missing(loadmt)){
 if(!is.list(loadmt)) stop("Provided `loadmt` should be a list with dim p*(q+intercept) restriction matrices for mu, sigma, tau, nu")
 names(loadmt) <- parY
 if(length(loadmt) != length(coefs)) stop("Number of matrices in `loadmt` should match number of matrices in `coefs`")
 for(i in parY){
  if(length(loadmt[[i]]) != p*t.[[i]]) stop("Provided `loadmt` is not lenght p*(q+intercept), revise dimension")
  if(!is.matrix(loadmt[[i]])) loadmt[[i]] <- matrix(loadmt[[i]], nrow = p, ncol = t.[[i]])
 }
} else {
 warning("Restrictions matrix not supplied, simple structure assumed", call. = F)
 loadmt <- NULL
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
   loadmt[[i]] <- tmp1; rm(tmp1,tmp2,r.)
  } else{
   loadmt[[i]] <- matrix(1,nrow = p, ncol = t.[[i]])
   if(t.[[i]]-1 > 1){
    warning(paste0("Nonlinear functions for ", i, ", identification restriction impossed in item 1"), call. = F)
    loadmt[[i]][1,t.[[i]]] <- 0
   }
  }
 }
}

# borg <- coefmod(bmu = coefs$mu, bsg = 0.5*coefs$sigma, loadmt = loadmt) # bsg = coefs$sigma

borg <- NULL
for(i in parY){borg[[i]] <- coefs[[i]]*loadmt[[i]]}
names(borg) <- parY

Y <- sapply(1:p, function(i) rFun(n,i,dist[i], Z.,borg))
Y <- as.data.frame(Y); colnames(Y) <- paste0("Y", 1:p)

return(list("Y" = Y, "Z" = Z., "FRes" = loadmt, "borg" = borg, "Z.mu" = muZ, "Z.Sigma" = sigZ,
            "formula" = form))
}
