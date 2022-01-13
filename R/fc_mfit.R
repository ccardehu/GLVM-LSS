
splvm.fit <- function(Y, fam, form,
                      control = list(method = c("ML","EM","hybrid","PML","PEM"), 
                                     start.val = NULL, constraint = NULL,
                                     ghQqp = 15, iter.lim = 150, full.hess = F, EM.iter.lim = 20,
                                     tol = sqrt(.Machine$double.eps), silent = F, information = "Fisher",
                                     pml.control = list(type = "lasso", lambda = 1, w.alasso = NULL,
                                         a = NULL, pen.load = F)) )
                      {

# Goal: Fits semi-parametric LVM
# Input : Y (item matrix), fam (distributions),
#         form (list of formulas for measurement eqs. for each parameter),
#         control (list of controls)
# Output: Estimated semi-parametric LVM
# Testing: Y = simR$Y; fam = fam; form = e.form;
#          control = list(method = "PEM", #start.val = lc, constraint = l1,
#          ghQqp = 15, iter.lim = 150, tol = sqrt(.Machine$double.eps), silent = F, full.hess = F, information = "Fisher")
#          control$pml.control = list(type = "alasso", lambda = 0.01, w.alasso = testa1$b, pen.load = F)

if(!is.matrix(Y)) Y <- as.matrix(Y)
parY <- unique(unlist(lapply(1:length(fam),function(i) pFun(fam[i]))))
for(i in parY){form[[i]] <- as.formula(form[[i]])}
lvar <- unique(unlist(lapply(1:length(form), function(i) all.vars(form[[i]]))))
lvar <- grep("Z", lvar, fixed = T, value = T)
if(length(lvar) == 0) stop("\n Latent variables in formula should be represented by letter Z (capital)")
q. <- length(lvar)
p. <- ncol(Y)
pC <- vector(mode = "list", length = length(parY)); names(pC) <- parY
for(i in parY){ for(j in 1:p.){ if(i %in% pFun(fam[j])) {pC[[i]] <- append(pC[[i]],j)} } }

# Gaussian Hermite Quadrature
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(is.null(control$ghQqp)) { if(q. == 1) control$ghQqp <- 25 else control$ghQqp <- 10 }
ghQ <- mvghQ(n = control$ghQqp, formula = form)

# Starting values
# ~~~~~~~~~~~~~~~
cycl <- iter <- 0
if(!is.null(control$tol)) tol <- control$tol else tol <- sqrt(.Machine$double.eps)
eps <- eps2 <- tol + 1; tol2 <- min(tol,1e-5) # max(tol,1e-5)
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
 if(!control$silent) cat("\n Argument 'control$start.val' not supplied: Initial guess defined using 'gamlss' & PCA (for factor scores).\n")
 bold <- suppressWarnings(ini.par(Y,fam,form,pC,q.))
}

# Constraints matrices (for each parameter)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
restr <- control$constraint
if(!is.null(restr)){
  if(!is.list(restr)) stop("Argument 'control$constraint' should be a list with element(s), each of the type c('parameter',item,'restricted variable',value)")
  restr <- rmat(restr,bold)
  loadmt <- vector(mode = "list",length = length(restr)); names(loadmt) <- names(restr)
  for(i in parY){loadmt[[i]] <- is.na(restr[[i]]); bold[[i]][!loadmt[[i]]] <- restr[[i]][!loadmt[[i]]] }
} else {
 if(!control$silent) cat("\n Argument 'control$constraint' not supplied: No (identification) restrictions assumed.\n")
 loadmt <- vector(mode = "list", length = length(parY)); names(loadmt) <- parY
 restr <- vector(mode = "list", length = length(loadmt)); names(restr) <- parY
 for(i in parY){
  if(length(all.vars(form[[i]])) > 1){
   r. <- p.%/%(ncol(ghQ$out[[i]])-1)
   tmp1 <- matrix(T, nrow = p., ncol = ncol(ghQ$out[[i]]), dimnames = list(colnames(Y), colnames(ghQ$out[[i]]))) # 0
   for(j in seq_along(all.vars(form[[i]]))){
    tmp2 <- seq((j-1)*(r.)+1,(j)*(r.),length.out = r.)
    tmp1[tmp2,all.vars(form[[i]])[j]] <- T
    if(j == length(all.vars(form[[i]])) && p.%%2 != 0) tmp1[(tail(tmp2,1)+1):nrow(tmp1),all.vars(form[[i]])[j]] <- T
   }
   tmp1[,!(colnames(tmp1) %in% all.vars(form[[i]]))] <- T
   loadmt[[i]] <- tmp1; rm(tmp1,tmp2)
  } else{
   loadmt[[i]] <- matrix(T,nrow = p., ncol = ncol(ghQ$out[[i]]))
  }
  colnames(loadmt[[i]]) <- colnames(ghQ$out[[i]])
  rownames(loadmt[[i]]) <- colnames(Y)
  restr[[i]] <- array(NA,dim = dim(loadmt[[i]]), dimnames = dimnames(loadmt[[i]]))
 }
}
loadmt2 <- loadmt

# Fixing values for restrictions parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in parY){ bold[[i]][-pC[[i]],] <- 0 }
  
# Estimation
# ~~~~~~~~~~
# if(sum(lb2cb(loadmt)) > p.*(p.-1)/2 && !control$silent) cat("\n Warning: Number of free parameters is less than p(p-1)/2. Model might be under-indentified.\n")
method <- control$method; pml.control <- control$pml.control
if(is.null(pml.control$pen.load)) {if(q. > 1){ pen.load <- T } else { pen.load <- F } } else pen.load <- pml.control$pen.load
if(is.null(control$information)) control$information <- "Fisher"
if(is.null(control$full.hess)) control$full.hess <- F
if(is.null(control$iter.lim)) control$iter.lim <- 3e2

stop.crit <- autoL <- F;
SSE <- olObj <- NULL; iiter <- "NA"
if(!is.null(pml.control$lambda)){
 if(pml.control$lambda == "auto" & !pml.control$type %in% c("lasso","alasso")) stop("\n Penalty type should be 'Alasso' or 'Lasso' if lambda = 'auto'")
 if(pml.control$lambda == "auto"){ pml.control$lambda <- 1/nrow(Y); autoL <- T } # starting lambda
 if(is.null(pml.control$gamma)){ pml.control$gamma <- 1.4 }
 if(is.null(pml.control$a)){pml.control$a <- 2 } 
}

hist.lambda <- ifelse(autoL, pml.control$lambda, NA) # starting lambda

tryCatch({
  
while(!stop.crit){

pen.idx <- pidx(bold,loadmt2,pen.load)

if(method == "EM"){
A1 <- dY(Y,ghQ,bold,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
llo <- sum(log(A2)) # log-likelihood old
pD <- exp(rowSums(A1,dim = 2))/A2 # this is EC (posterior density)
dvL <- dvY(Y,ghQ,bold,fam,control$information) # List with all the necessary derivatives
A3 <- sche(ghQ,bold,loadmt2,fam,dvL,pD,control$information,control$full.hess) # score & Hessian object
A4 <- upB(bold,A3,loadmt2) # updated betas
# Compute log-likelihood for comparison, etc.
bnew <- A4$b
A1 <- dY(Y,ghQ,bnew,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
upll <- lln <- sum(log(A2)) ; dlln <- round(lln,3) # log-likelihood new
eps <- abs(lln-llo) #/(0.1+abs(lln))
iter <- iter + 1
bold <- bnew
if(control$silent == F) cat("\r EM iter: ", iter, ", loglk: ", format(round(dlln, digits = 5), nsmall = 3), ", \U0394 loglk: ", format(round(lln-llo, digits = 3), scientific = T, nsmall = 3), sep = "")
mod.grad <- list(unp = A4$gradient)
mod.hess <- list(unp = A4$hessian)
adsol <- A4$adsol
for(r in names(bnew)){
  if("Z1" %in% colnames(bnew[[r]]) && bnew[[r]][1,"Z1"] < 0) bnew[[r]][,"Z1"] <- -bnew[[r]][,"Z1"]
  }
}
  
if(method == "PEM"){
A1 <- dY(Y,ghQ,bold,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
A2a <- lb2pM(bold,Y,pen.idx,loadmt2,pml.control)
llo <- c(sum(log(A2)) - A2a$lp) # penalized log-likelihood old
pD <- exp(rowSums(A1,dim = 2))/A2 # this is EC (posterior density)
dvL <- dvY(Y,ghQ,bold,fam,control$information) # List with all the necessary derivatives
A3 <- sche(ghQ,bold,loadmt2,fam,dvL,pD,control$information,control$full.hess) # score & Hessian object
A4 <- upB.pen(bold,A3,A2a,loadmt2) # updated betas
# Compute log-likelihood for comparison, etc.
bnew <- A4$b
loadmt2 <- uplm(bnew,loadmt2)
pen.idx <- pidx(bnew,loadmt2,pen.load)
A1 <- dY(Y,ghQ,bnew,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
A2a <- lb2pM(bnew,Y,pen.idx,loadmt2,pml.control)
lln <- c(sum(log(A2)) - A2a$lp) # penalized log-likelihood new
pD <- exp(rowSums(A1,dim = 2))/A2 # this is EC (posterior density)
dvL <- dvY(Y,ghQ,bnew,fam,control$information) # List with all the necessary derivatives
A3 <- sche(ghQ,bnew,loadmt2,fam,dvL,pD,control$information,control$full.hess) # score & Hessian object
dlln <- round(lln,3) # penalized log-likelihood new (printing)
eps <- abs(lln-llo) #/(0.1+abs(lln))
iter <- iter + 1
bold <- bnew
mod.grad <- list(pen = A4$gradient, unp = A3$gradient)
mod.hess <- list(pen = A4$hessian, unp = A3$hessian)
adsol <- A4$adsol
upll <- sum(log(A2))
for(r in names(bnew)){
  if("Z1" %in% colnames(bnew[[r]]) && bnew[[r]][1,"Z1"] < 0) bnew[[r]][,"Z1"] <- -bnew[[r]][,"Z1"]
}

if(!control$silent){
 if(!autoL){
  cat("\r Penalised EM iter: ", iter, ", loglk: ", format(round(dlln, digits = 5), nsmall = 3), ", \U0394 loglk: ", format(round(lln-llo, digits = 3), scientific = T, nsmall = 3), sep = "") } else {
    if(iter == 1){
  cat("\n (Automatic) Penalised EM iter: ", iter, ", loglk: ", format(round(dlln, digits = 5), nsmall = 3), ", \U0394 loglk: ", format(round(lln-llo, digits = 3), scientific = T, nsmall = 3), ", (\U03bb: ", format(round(pml.control$lambda, digits = 5), nsmall = 5), ", mid-cycle: ", cycl+1, ", iters: ", iiter,")",sep = "") } else {
  cat("\r (Automatic) Penalised EM iter: ", iter, ", loglk: ", format(round(dlln, digits = 5), nsmall = 3), ", \U0394 loglk: ", format(round(lln-llo, digits = 3), scientific = T, nsmall = 3), ", (\U03bb: ", format(round(pml.control$lambda, digits = 5), nsmall = 5), ", mid-cycle: ", cycl+1, ", iters: ", iiter,")",sep = "") }
 } }
}

if(method == "ML"){
b2r <- lb2cb(bold)
btr <- b2r[t(lb2mb(loadmt2))]
# control$full.hess <- T;  control$information = "Hessian"
if(!control$silent){
  if(exists("catmsg")) cat(paste0("\n Using trust-region algorithm to refine ML estimates ", catmsg))
  else cat("\n Using trust-region algorithm to find ML estimates ...") }
r1 <- trust::trust(objfun = loglkf, parinit = btr, rinit = 1, rmax = 5,fterm = tol,
                   iterlim = control$iter.lim, minimize = F, Y = Y, bg = bold, ghQ = ghQ, fam = fam,
                   info = control$information, res = loadmt2, full = control$full.hess)
b2r[t(lb2mb(loadmt2))] <- r1$argument
bold <- bnew <- cb2lb(b2r,bold)
lln <- r1$value
iter <- r1$iter; eps <- 0
if(!control$silent) cat("\n Converged after ", r1$iter, " iterations (loglk: ", round(r1$value,3),")", sep = "")
mod.grad <- list(unp = r1$gradient)
mod.hess <- list(unp = r1$hessian)
upll <- lln
adsol <- m2pdm(-mod.hess$unp)$is.PDM
for(r in names(bnew)){
  if("Z1" %in% colnames(bnew[[r]]) && bnew[[r]][1,"Z1"] < 0) bnew[[r]][,"Z1"] <- -bnew[[r]][,"Z1"]
  }
}
  
if(method == "PML"){
b2r <- lb2cb(bold)
btr <- b2r[t(lb2mb(loadmt2))]
if(!control$silent){
  if(exists("catmsg")) { cat(paste0("\n Using trust-region algorithm to refine Penalised ML estimates ", catmsg)) } else {
    if(cycl == 0) cat("\n Using trust-region algorithm to find Penalised ML estimates ...") } }
r1 <- trust::trust(objfun = penloglkf, parinit = btr, rinit = 1, rmax = 5,fterm = tol,
                   iterlim = control$iter.lim, minimize = F, Y = Y, bg = bold, ghQ = ghQ, fam = fam,
                   pml.control = pml.control, pen.idx = pen.idx, info = control$information,
                   res = loadmt2, full = control$full.hess)
rep <- r1$argument; # rep[abs(rep) < 1e5*sqrt(.Machine$double.eps)] <- 0
b2r[t(lb2mb(loadmt2))] <- rep; rm(rep)
bold <- bnew <- cb2lb(b2r,bold)
loadmt2 <- uplm(bnew,loadmt2)
pen.idx <- pidx(bnew,loadmt2,pen.load)
lln <- r1$value
iter <- r1$iter; eps <- 0
if(!control$silent){
 if(!autoL){
  cat("\n PMLE converged after ", r1$iter, " iterations (loglk: ", format(round(r1$value,3), nsmall = 3),")", sep = "") } else {
  cat("\n (Automatic) PMLE converged after ", r1$iter, " iterations (loglk: ", format(round(r1$value,3), nsmall = 3), ", \U03bb: ", format(round(pml.control$lambda, digits = 5), nsmall = 5), ", mid-cycle: ", cycl+1, ", iters: ", iiter,")",sep = "")
  } }
A1 <- dY(Y,ghQ,bnew,fam); A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights)
upll <- sum(log(A2))
pD <- exp(rowSums(A1,dim = 2))/A2 # this is EC (posterior density)
dvL <- dvY(Y,ghQ,bnew,fam,control$information) # List with all the necessary derivatives
A3 <- sche(ghQ,bnew,loadmt2,fam,dvL,pD,control$information,control$full.hess) # score & Hessian object
mod.grad <- list(pen = r1$gradient, unp = A3$gradient)
mod.hess <- list(pen = r1$hessian, unp = A3$hessian)
adsol <- m2pdm(-mod.hess$unp)$is.PDM
for(r in names(bnew)){
  if("Z1" %in% colnames(bnew[[r]]) && bnew[[r]][1,"Z1"] < 0) bnew[[r]][,"Z1"] <- -bnew[[r]][,"Z1"]
  }
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
if(!control$silent) cat("\r Hybrid-EM iter: ", iter, ", loglk: ", dlln, ", \U0394 loglk: ", lln-llo, sep = "")
mod.grad <- A3$gradient
mod.hess <- A3$hessian
upll <- lln
for(r in names(bnew)){ if("Z1" %in% colnames(bnew[[r]]) && bnew[[r]][1,"Z1"] < 0) bnew[[r]][,"Z1"] <- -bnew[[r]][,"Z1"] }
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
pen.idx <- pidx(bnew,loadmt,pen.load)
A1 <- dY(Y,ghQ,bnew,fam)
A2 <- c(exp(rowSums(A1,dim = 2))%*%ghQ$weights) # this is efy
lln <- c(sum(log(A2)) - lb2pM(bnew,Y,pen.idx,pml.control)) ; dlln <- round(lln,3) # log-likelihood new
eps <- abs(lln-llo)
iter <- iter + 1
bold <- bnew
if(iter >= control$EM.iter.lim){ method <- "PML"; catmsg <- paste0("(after ", iter," (Penalised) Hybrid-EM iterations)") }
if(!control$silent) cat("\r Penalised Hybrid-EM iter: ", iter, ", loglk: ", dlln, ", \U0394 loglk: ", lln-llo, sep = "")
mod.grad <- A3$gradient
mod.hess <- A3$hessian
upll <- sum(log(A2))
for(r in names(bnew)){ if("Z1" %in% colnames(bnew[[r]]) && bnew[[r]][1,"Z1"] < 0) bnew[[r]][,"Z1"] <- -bnew[[r]][,"Z1"]  }
}

if(!control$silent & cycl+1 == 20) { cat("\n Note: Automatic selection of \U03bb reached maximum number of mid-cycles (20)", sep = "") }

if(autoL && (eps < tol | iter >= control$iter.lim) && (eps2 > tol2 && cycl < 19)){
 cycl <- cycl + 1
 # pD <- exp(rowSums(A1,dim = 2))/A2 # this is EC (posterior density)
 # dvL <- dvY(Y,ghQ,bnew,fam,control$information) # List with all the necessary derivatives
 # A3 <- sche(ghQ,bnew,loadmt2,fam,dvL,pD,control$information,control$full.hess) # score & Hessian object
 olObj <- op.lambda(bnew,Y,pen.idx,loadmt2,pml.control,A3,control$iter.lim,tol2)
 laold <- pml.control$lambda; SSE <- olObj$sse
 lanew <- pml.control$lambda <- olObj$lambda ; iiter <- olObj$iter
 hist.lambda <- c(hist.lambda,lanew)
 eps2 <- abs(lanew - laold) # /(0.1+abs(lanew))
 eps <- 1
 iter <- 0 }

stop.crit <- (eps < tol || iter >= control$iter.lim)

}
  
conv <- ifelse(iter == control$iter.lim | cycl == 19 | !adsol, F, T)
  
for(r in names(bnew)){ for(j in 1:q.){
  if(paste0("Z",j) %in% colnames(bnew[[r]]) && bnew[[r]][j,paste0("Z",j)] < 0 && sum(bnew[[r]][j,-1] != 0) == 1) bnew[[r]][,paste0("Z",j)] <- -bnew[[r]][,paste0("Z",j)]
} }

return(list(b = bnew, loglik = lln, uploglik = upll, loadmt = loadmt2, iter = iter, iiter= cycl, ghQ = ghQ,
            Y = as.data.frame(Y), fam = fam, formula = form, eps = eps, method = method,
            pml.control = pml.control, gradient = mod.grad, hessian = mod.hess, info = control$information,
            pen.idx = pen.idx, hist = hist.lambda, conv = conv, sse = SSE))
},
error = function(e){
if(!control$silent) cat(paste("\n Error in Estimation, proceeded with next simulation \n Error: ",e))
berr <- lapply(parY, function(i) matrix(-999, nrow = ncol(Y), ncol = ncol(ghQ$out[[i]])))
names(berr) <- parY; for(i in parY){ berr[[i]] <- berr[[i]]*loadmt[[i]] ; loadmt[[i]] <- loadmt[[i]]*-999}
return(list(b = berr, loglik = -999, loadmt = loadmt,iter = -999, pml.control = pml.control, conv = 0)) } )
}

