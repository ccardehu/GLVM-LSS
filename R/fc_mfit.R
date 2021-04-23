
splvm.fit <- function(Y, fam, form,
                      control = list(method = c("ML","EM","hybrid","PML","PEM"), 
                                     start.val = NULL, constraint = NULL,
                                     ghQqp = 15, iter.lim = 150, full.hess = F, EM.iter.lim = 20,
                                     tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher",
                                     pml.control = list(type = "lasso", lambda = 1, w.alasso = NULL,
                                         gamma = 1, a = 3.7, pen.load = F)) )
                      {

# Goal: Fits semi-parametric LVM
# Input : Y (item matrix), fam (distributions),
#         form (list of formulas for measurement eqs. for each parameter),
#         control (list of controls)
# Output: Estimated semi-parametric LVM
# Testing: Y = simR$Y; fam = fam; form = e.form;
#          control = list(method = "PEM", start.val = lc, constraint = l1,
#          ghQqp = 15, iter.lim = 150, tol = sqrt(.Machine$double.eps), silent = F, full.hess = F, information = "Fisher")
#          control$pml.control = list(type = "scad", lambda = 1, w.alasso = 1, gamma = 1, a = 3.7, pen.load = T)

if(!is.matrix(Y)) Y <- as.matrix(Y)
parY <- unique(unlist(lapply(1:length(fam),function(i) pFun(fam[i]))))
for(i in parY){form[[i]] <- as.formula(form[[i]])}
lvar <- unique(unlist(lapply(1:length(form), function(i) all.vars(form[[i]]))))
lvar <- lvar[grep("Z", lvar, fixed = T)]
if(length(lvar) == 0) stop("\n Latent variables in formula should be represented by letter Z (capital)")
q. <- length(lvar)
p. <- ncol(Y)
pC <- vector(mode = "list", length = length(parY)); names(pC) <- parY
for(i in parY){ for(j in 1:p.){ if(i %in% pFun(fam[j])) {pC[[i]] <- append(pC[[i]],j)} } }

# Gaussian Hermite Quadrature
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(is.null(control$ghQqp)) control$ghQqp <- 15
ghQ <- mvghQ(n = control$ghQqp, formula = form)

# Constraint matrices (for each parameter)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

loadmt <- control$constraint
if(!is.null(loadmt)){
  if(!is.list(loadmt)) stop("Constraint should be a list with dim p*(q+intercept) restriction matrices for mu, sigma, tau, nu")
  if(length(loadmt) != length(parY)) stop("Number of matrices in `loadmt` should match number of parameters mu, sigma, tau, nu")
  names(loadmt) <- parY 
  for(i in parY) {
   colnames(loadmt[[i]]) <- colnames(ghQ$out[[i]]); rownames(loadmt[[i]]) <- colnames(Y)
   if(length(loadmt[[i]]) != p.*ncol(ghQ$out[[i]])) stop("Provided `loadmt` is not lenght p*(q+intercept), revise dimension")
   if(!is.matrix(loadmt[[i]])) loadmt[[i]] <- matrix(loadmt[[i]], nrow = p., ncol = ncol(ghQ$out[[i]]))
  }
} else {
 cat("\n Constraint matrices not supplied. No identification restrictions assumed.\n")
 loadmt <- NULL
 for(i in parY){
  if(length(all.vars(form[[i]])) > 1){
   r. <- p.%/%(ncol(ghQ$out[[i]])-1)
   tmp1 <- matrix(1, nrow = p., ncol = ncol(ghQ$out[[i]]), dimnames = list(colnames(Y), colnames(ghQ$out[[i]]))) # 0
   for(j in seq_along(all.vars(form[[i]]))){
    tmp2 <- seq((j-1)*(r.)+1,(j)*(r.),length.out = r.)
    tmp1[tmp2,all.vars(form[[i]])[j]] <- 1
    if(j == length(all.vars(form[[i]])) && p.%%2 != 0) tmp1[(tail(tmp2,1)+1):nrow(tmp1),all.vars(form[[i]])[j]] <- 1
   }
   tmp1[,!(colnames(tmp1) %in% all.vars(form[[i]]))] <- 1
   loadmt[[i]] <- tmp1; rm(tmp1,tmp2)
  } else{
   loadmt[[i]] <- matrix(1,nrow = p., ncol = ncol(ghQ$out[[i]]))
   # loadmt[[i]][1,ncol(gr$out[[i]])] <- 0 # ! HERE
   # if(ncol(ghQ$out[[i]])-1 > 1){
    # message(paste0("Nonlinear functions for ", i, ", identification restriction imposed in item 1"))
    # loadmt[[i]][1,ncol(ghQ$out[[i]])] <- 0
   # }
  }
  colnames(loadmt[[i]]) <- colnames(ghQ$out[[i]])
  rownames(loadmt[[i]]) <- colnames(Y)  
 }
}

# Starting values
# ~~~~~~~~~~~~~~~
iter <- 0
if(!is.null(control$tol)) tol <- control$tol else tol <- sqrt(.Machine$double.eps)
eps <- tol + 1
icoefs <- control$start.val
if(!is.null(icoefs)){
  if(!is.list(icoefs)) stop("\n Provided starting values should be a list with dim p*(q+intercept) loading matrices for mu, sigma, tau, nu")
  if(length(icoefs) != length(parY)) stop("\n Number of matrices in starting values should match number of parameters mu, sigma, tau, nu")
  names(icoefs) <- parY 
  for(i in parY) {
    colnames(icoefs[[i]]) <- colnames(ghQ$out[[i]]); rownames(icoefs[[i]]) <- colnames(Y)
    if(length(icoefs[[i]]) != p.*ncol(ghQ$out[[i]])) stop("\n Provided starting values is not lenght p*(q+intercept), revise dimension")
    if(!is.matrix(icoefs[[i]])) icoefs[[i]] <- matrix(icoefs[[i]], nrow = p., ncol = ncol(ghQ$out[[i]]))
    icoefs[[i]] <- icoefs[[i]] + runif(length(icoefs[[i]]), min = -2*tol, max = 2*tol)
  }
  bold <- icoefs
} else {
 cat("\n Starting values not supplied. Using 'gamlss' & PCA for initial values .... (Done)\n")
 # bold <- lapply(parY, function(i) matrix(1, nrow = ncol(Y), ncol = ncol(ghQ$out[[i]])))
 # names(bold) <- parY
 bold <- ini.par(Y,fam,form,pC,q.)
 }
for(i in parY){
 bold[[i]] <- bold[[i]]*loadmt[[i]]
 bold[[i]][-pC[[i]],] <- 0
}
  
# Estimation
# ~~~~~~~~~~
if(sum(lb2cb(loadmt)) > p.*(p.-1)/2) cat("\n Warning: Number of parameters is < p(p-1)/2. Model might be under-indentified.\n")
method <- control$method
pml.control <- control$pml.control
if(is.null(pml.control$pen.load)) pen.load <- F else pen.load <- pml.control$pen.load
if(is.null(control$information)) control$information <- "Fisher"

tryCatch({
while(eps > tol && iter < control$iter.lim){

pen.idx <- pidx(bold,pen.load)

if(method == "EM"){
A1 <- dY(Y,ghQ,bold,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
llo <- sum(log(A2)) # log-likelihood old
pD <- exp(rowSums(A1,dim = 2))/A2 # this is EC (posterior density)
dvL <- dvY(Y,ghQ,bold,fam,control$information) # List with all the necessary derivatives
A3 <- upB(bold,ghQ,fam,dvL,pD,full.hess = control$full.hess,information = control$information) # updated betas
# Compute log-likelihood for comparison, etc.
bnew <- A3$b
A1 <- dY(Y,ghQ,bnew,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
lln <- sum(log(A2)) ; dlln <- round(lln,3) # log-likelihood new
eps <- abs(lln-llo)
iter <- iter + 1
bold <- bnew
if(control$silent == F) cat("\r EM iter: ", iter, ", loglk: ", dlln, ", \U0394 loglk: ", lln-llo, sep = "")
mod.grad <- A3$gradient
mod.hess <- A3$hessian
}
  
if(method == "PEM"){
A1 <- dY(Y,ghQ,bold,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
llo <- c(sum(log(A2)) - lb2pM(bold,Y,pen.idx,pml.control)) # penalized log-likelihood old
pD <- exp(rowSums(A1,dim = 2))/A2 # this is EC (posterior density)
dvL <- dvY(Y,ghQ,bold,fam,control$information) # List with all the necessary derivatives
A3 <- upB.pen(bold,ghQ,fam,dvL,pD,full.hess = control$full.hess,pen.idx,
                pml.control = pml.control,information = control$information) # updated betas
# Compute log-likelihood for comparison, etc.
bnew <- A3$b
pen.idx <- pidx(bnew,pen.load)
A1 <- dY(Y,ghQ,bnew,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
lln <- c(sum(log(A2)) - lb2pM(bnew,Y,pen.idx,pml.control)) # penalized log-likelihood new
dlln <- round(lln,3) # penalized log-likelihood new (printing)
eps <- abs(lln-llo)
iter <- iter + 1
bold <- bnew
if(control$silent == F) cat("\r (Penalised) EM iter: ", iter, ", loglk: ", dlln, ", \U0394 loglk: ", lln-llo, sep = "")
mod.grad <- A3$gradient
mod.hess <- A3$hessian
}

if(method == "ML"){
b2r <- lb2cb(bold)
btr <- b2r[lb2cb(bold) != 0]
if(control$silent == F){
  if(exists("catmsg")) cat(paste0("\n Using trust-region algorithm to refine ML estimates ", catmsg))
  else cat("\n Using trust-region algorithm to find ML estimates") }
r1 <- trust::trust(objfun = loglkf, parinit = btr, rinit = 1, rmax = 10,
                   iterlim = control$iter.lim, minimize = F, Y = Y, bg = bold, ghQ = ghQ, fam = fam, info = control$information)
b2r[lb2cb(bold) != 0] <- r1$argument
bnew <- cb2lb(b2r,bold)
lln <- r1$value
iter <- control$iter.lim + 1; eps <- 0
if(control$silent == F) cat("\n Converged after ", r1$iter, " iterations (loglk: ", round(r1$value,3),")", sep = "")
mod.grad <- r1$gradient
mod.hess <- r1$hessian
}
  
if(method == "PML"){
b2r <- lb2cb(bold)
btr <- b2r[lb2cb(bold) != 0]
if(control$silent == F){
  if(exists("catmsg")) cat(paste0("\n Using trust-region algorithm to refine Penalised ML estimates ", catmsg))
  else cat("\n Using trust-region algorithm to find Penalised ML estimates") }
r1 <- trust::trust(objfun = penloglkf, parinit = btr, rinit = 1, rmax = 10,
                   iterlim = control$iter.lim, minimize = F, Y = Y, bg = bold, ghQ = ghQ, fam = fam,
                   pml.control = pml.control, pen.idx = pen.idx, info = control$information)
rep <- r1$argument; rep[abs(rep) < 1*sqrt(.Machine$double.eps)] <- 0
b2r[lb2cb(bold) != 0] <- rep; rm(rep)
bnew <- cb2lb(b2r,bold)
lln <- r1$value
iter <- control$iter.lim + 1; eps <- 0
if(control$silent == F) cat("\n Converged after ", r1$iter, " iterations (loglk: ", round(r1$value,3),")", sep = "")
mod.grad <- r1$gradient
mod.hess <- r1$hessian
}

if(method == "hybrid"){
A1 <- dY(Y,ghQ,bold,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
llo <- sum(log(A2)) # log-likelihood old
pD <- exp(rowSums(A1,dim = 2))/A2 # this is EC (posterior density)
dvL <- dvY(Y,ghQ,bold,fam,control$information) # List with all the necessary derivatives
A3 <- upB(bold,ghQ,fam,dvL,pD,full.hess = control$full.hess,information = control$information) #updated betas
# Compute log-likelihood for comparison, etc.
bnew <- A3$b
A1 <- dY(Y,ghQ,bnew,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
lln <- sum(log(A2)) ; dlln <- round(lln,3) # log-likelihood new
eps <- abs(lln-llo)
iter <- iter + 1
bold <- bnew
if(iter >= control$EM.iter.lim){ method <- "ML"; catmsg <- paste0("(after ", iter," Hybrid-EM iterations)") }
if(control$silent == F) cat("\r Hybrid-EM iter: ", iter, ", loglk: ", dlln, ", \U0394 loglk: ", lln-llo, sep = "")
mod.grad <- A3$gradient
mod.hess <- A3$hessian
}

if(method == "P-hybrid"){
A1 <- dY(Y,ghQ,bold,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
llo <- c(sum(log(A2)) - lb2pM(bold,Y,pen.idx,pml.control)) # log-likelihood old
pD <- exp(rowSums(A1,dim = 2))/A2 # this is EC (posterior density)
dvL <- dvY(Y,ghQ,bold,fam,control$information) # List with all the necessary derivatives
A3 <- upB.pen(bold,ghQ,fam,dvL,pD,full.hess = control$full.hess,pen.idx,
                pml.control = pml.control,information = control$information) #updated betas
# Compute log-likelihood for comparison, etc.
bnew <- A3$b
pen.idx <- pidx(bnew,pen.load)
A1 <- dY(Y,ghQ,bnew,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
lln <- c(sum(log(A2)) - lb2pM(bnew,Y,pen.idx,pml.control)) ; dlln <- round(lln,3) # log-likelihood new
eps <- abs(lln-llo)
iter <- iter + 1
bold <- bnew
if(iter >= control$EM.iter.lim){ method <- "PML"; catmsg <- paste0("(after ", iter," (Penalised) Hybrid-EM iterations)") }
if(control$silent == F) cat("\r Penalised Hybrid-EM iter: ", iter, ", loglk: ", dlln, ", \U0394 loglk: ", lln-llo, sep = "")
mod.grad <- A3$gradient
mod.hess <- A3$hessian
}
    
}
return(list(b = bnew, loglik = lln, loadmt = loadmt, iter = iter, ghQ = ghQ,
            Y = as.data.frame(Y), fam = fam, formula = form, eps = eps, method = method,
            pml.control = pml.control, gradient = mod.grad, hessian = mod.hess, info = control$information))
},
error = function(e){
cat(paste("\n Error in Estimation, proceeded with next simulation \n Error:",e))
berr <- lapply(parY, function(i) matrix(-999, nrow = ncol(Y), ncol = ncol(ghQ$out[[i]])))
names(berr) <- parY; for(i in parY){ berr[[i]] <- berr[[i]]*loadmt[[i]] ; loadmt[[i]] <- loadmt[[i]]*-999}
return(list(b = berr, loglik = -999, loadmt = loadmt,iter = -999)) } )
}

