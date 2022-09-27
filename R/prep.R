prep_cont <- function(){
  # Control settings
  # ~~~~~~~~~~~~~~~~
  con <- list(EM_iter = 30, EM_use2d = T, iter.lim = 300,
              EM_appHess = F, EM_lrate = 0.001, est.ci = T,
              solver = "trust", start.val = NULL, mat.info = "Hessian", lazytrust = F,
              iden.res = NULL, tol = sqrt(.Machine$double.eps), tolb = 1e-4,
              corr.lv = FALSE,
              nQP = if(q == 1) 40 else { if(q == 2) ifelse(p <= 10, 15, 25) else 10 },
              verbose = FALSE, autoL_iter = 30, f.scores = F,
              penalty = "none", lambda = NULL, w.alasso = NULL, gamma = NULL, a = NULL)
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(namc[!namc %in% namC]) > 0)
    warning("Unknown names in control: ", paste(namc[!namc %in% namC], collapse = ", "))
  return(con)
}

prep_form <- function(){
  # Formula settings
  # ~~~~~~~~~~~~~~~~
  is.formula <- function(x){ inherits(x,"formula") }
  if(!is.formula(mu.eq) & !is.null(mu.eq)) mu.eq <- as.formula(mu.eq)
  if(!is.formula(sg.eq) & !is.null(sg.eq)) sg.eq <- as.formula(sg.eq)
  if(!is.formula(ta.eq) & !is.null(ta.eq)) ta.eq <- as.formula(ta.eq)
  if(!is.formula(nu.eq) & !is.null(nu.eq)) nu.eq <- as.formula(nu.eq)
  form <- list("mu" = mu.eq, "sigma" = sg.eq, "tau" = ta.eq, "nu" = nu.eq)
  for(i in names(form)){ if(is.null(form[[i]])) form[[i]] <- NULL } 
  return(form)
}

prep_fam <- function(){
  # Family settings
  # ~~~~~~~~~~~~~~~
  if(!is.list(family)) stop("Argument `family' must be a list of valid distributions.")
  if(!all(sapply(1:length(family), function(i) class(family[[i]]) == "dist_glvmlss"))){
    stop("Arguments in `family' must be of class `dist_glvmlss'.")
  }
  return(family)
}

prep_Z <- function(){
  # Factor settings
  # ~~~~~~~~~~~~~~~
  lvar <- unique(unlist(lapply(1:length(form), function(i) all.vars(form[[i]]))))
  lvar <- grep("Z", lvar, fixed = T, value = T)
  if(length(lvar) == 0) stop("\n Latent variables in formula should be represented by capital letter Z")
  return(length(lvar))
}

prep_Y <- function(data){
  # Item settings
  # ~~~~~~~~~~~~~
  Y <-list(Y = NULL, na.idx = NULL)
  Y$Y <- data; if(!is.matrix(Y$Y)) Y$Y <- as.matrix(Y$Y)
  if(is.null(colnames(Y$Y))) colnames(Y$Y) <- paste0("Y",1:ncol(Y$Y))
  Y$na.idx <- !is.na(Y$Y)
  return(Y)
}

prep_ghq <- function(){
  # GHQ setting
  # ~~~~~~~~~~~
  if(control$corr.lv == T){ PsiZ <- diag(q); PsiZ[lower.tri(sigZ)] <- PsiZ[upper.tri(PsiZ)] <- 0.2 } else PsiZ <- diag(q)
  return(mvghQ(n = control$nQP, formula = form, psi = PsiZ))
}

prep_stva <- function(){
  # Starting betas
  # ~~~~~~~~~~~~~~
  if(!is.null(control$start.val)){
    if(!is.list(control$start.val)) stop("\n Provided starting values should be a list with dim p*(q+intercept) loading matrices for mu, sigma, tau, nu")
    if(length(control$start.val) != length(form)) stop("\n Number of matrices in starting values should match number of parameters mu, sigma, tau, nu")
    names(control$start.val) <- names(form)
    for(i in names(control$start.val)) {
      colnames(control$start.val[[i]]) <- colnames(ghQ$out[[i]]); rownames(control$start.val[[i]]) <- colnames(Y$Y)
      if(length(control$start.val[[i]]) != p*ncol(ghQ$out[[i]])) stop("\n Provided starting values is not lenght p*(q+intercept), revise dimension")
      if(!is.matrix(control$start.val[[i]])) control$start.val[[i]] <- matrix(control$start.val[[i]], nrow = p, ncol = ncol(ghQ$out[[i]]))
    }
    b <- control$start.val
  } else {
    if(control$verbose) cat("\n Argument 'control$start.val' not supplied: Starting values via 'gamlss' and PCA.")
    b <- suppressWarnings(ibeta(Y,famL,form))
  }
  
  # Constraints matrices (for each parameter)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(control$iden.res)){
    if(q == 1){
      rb <- penb <- b
      for(i in 1:length(b)){
        rb[[i]] <- penb[[i]] <- !is.na(rb[[i]])
        penb[[i]][,grepl("(Intercept)", colnames(penb[[i]]), fixed = T)] <- F 
        penb[[i]][1,2] <- F }
      b <-  list(b = b, rb = rb, penb = penb)
      return(b) }
    if(control$verbose & q != 1) cat("\n Argument 'control$iden.res' not supplied: Using recursive identification restriction.")
    control$iden.res <- "recursive" }
  if(is.character(control$iden.res)){
    if(!control$iden.res %in% c("orthogonal", "eiv", "recursive")) stop("Argument 'control$iden.res' should be one of c(orthogonal, eiv, recursive).")
    if(control$iden.res == "orthogonal"){
      rb <- penb <- b
      for(i in 1:length(b)){
        b[[i]][,!grepl("(Intercept)", colnames(b[[i]]), fixed = T)] <- svd(b[[i]][,!grepl("(Intercept)", colnames(b[[i]]), fixed = T)])$u
        rb[[i]] <- penb[[i]] <- !is.na(rb[[i]])
        penb[[i]][,grepl("(Intercept)", colnames(penb[[i]]), fixed = T)] <- F } }
    if(control$iden.res == "eiv"){
      rb <- penb <- b
      for(i in 1:length(b)){
        for(r in 1:q){ b[[i]][r, !colnames(b[[i]]) %in% c("(Intercept)", paste0("Z",r))] <- 0 }
        rb[[i]] <- penb[[i]] <- b[[i]] != 0 
        penb[[i]][,grepl("(Intercept)", colnames(penb[[i]]), fixed = T)] <- F 
        diag(penb[[i]][,!grepl("(Intercept)", colnames(penb[[i]]), fixed = T)]) <- F } }
    if(control$iden.res == "recursive"){
      rb <- penb <- b
      for(i in 1:length(b)){
        b[[i]][,!colnames(b[[i]]) %in% c("(Intercept)")][upper.tri(b[[i]][,!colnames(b[[i]]) %in% c("(Intercept)")])] <- 0
        rb[[i]] <- penb[[i]] <- b[[i]] != 0 
        penb[[i]][,grepl("(Intercept)", colnames(penb[[i]]), fixed = T)] <- F
        diag(penb[[i]][,!grepl("(Intercept)", colnames(penb[[i]]), fixed = T)]) <- F} }
    b <-  list(b = b, rb = rb, penb = penb)
  } else {
    if(q == 1){
      if(control$verbose) cat("\n Argument 'control$iden.res' not needed for one-factor models (q == 1).")
      rb <- penb <- b
      for(i in 1:length(b)){
        rb[[i]] <- penb[[i]] <- !is.na(rb[[i]])
        penb[[i]][,grepl("(Intercept)", colnames(penb[[i]]), fixed = T)] <- F 
        penb[[i]][1,2] <- F }
      b <-  list(b = b, rb = rb, penb = penb)
      return(b) }
    if(!is.list(control$iden.res)) stop("Argument 'control$iden.res' should be a list with element(s), each of the type c('parameter',item,'restricted variable',value)")
    b <- rmat(control$iden.res,b) }
  # Sign
  # ~~~~
  for(r in names(b$b)){
    for(j in which(grepl("Z",colnames(b$b[[r]])))){
      if(b$b[[r]][,j][b$b[[r]][,j] != 0][1] < 0) b$b[[r]][,j] <- -b$b[[r]][,j]
    }
  }
  return(b)
}

sim_cont <- function(){
  # Control settings
  # ~~~~~~~~~~~~~~~~
  con <- list(muZ = rep(0,q), PsiZ = diag(q), start.val = NULL, iden.res = NULL,
              verbose = FALSE)
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(namc[!namc %in% namC]) > 0)
    warning("Unknown names in control: ", paste(namc[!namc %in% namC], collapse = ", "))
  return(con)
}

sim_Z <- function(){
  if(sum(diag(control$PsiZ)) != q){ warning("Argument `control$PsiZ' standardised to a correlation matrix")
    control$PsiZ <- diag(1/sqrt(diag(control$PsiZ)))%*%control$PsiZ%*%diag(1/sqrt(diag(control$PsiZ)))
  }
  Z <- as.data.frame(mvtnorm::rmvnorm(n, control$muZ, control$PsiZ))
  colnames(Z) <- paste0("Z", 1:q)
  Zmod <- vector("list", length = length(form)); names(Zmod) <- names(form)
  for(i in names(form)){ Zmod[[i]] <- as.data.frame(model.matrix(form[[i]],Z)) }
  return(list(Z = Z, Zmod = Zmod))
}

sim_stva <- function(){
  # Starting betas
  # ~~~~~~~~~~~~~~
  if(!is.null(control$start.val)){
    if(!is.list(control$start.val)) stop("\n Provided starting values should be a list with dim p*(q+intercept) loading matrices for mu, sigma, tau, nu")
    if(length(control$start.val) != length(form)) stop("\n Number of matrices in starting values should match number of parameters mu, sigma, tau, nu")
    names(control$start.val) <- names(form)
    for(i in names(control$start.val)){
      if(length(control$start.val[[i]]) != p*ncol(Z$Zmod[[i]])) stop("Provided `control$start.val` is not lenght p*(q+intercept), revise dimension")
      if(!is.matrix(control$start.val[[i]])) control$start.val[[i]] <- matrix(control$start.val[[i]], nrow = p, ncol = ncol(Z$Zmod[[i]]))
      colnames(control$start.val[[i]]) <- colnames(Z$Zmod[[i]]); rownames(control$start.val[[i]]) <- paste0("Y", 1:p)
    }
    b <- control$start.val
  } else {
    if(control$verbose) cat("\n Argument `control$start.val` not supplied: check code for assumed")
    control$start.val <- vector("list",length = length(form))
    names(control$start.val) <- names(form)
    for(i in names(form)){
      tmp <- c(rep(0.4,p), rep(0.8,p*(ncol(Z$Zmod[[i]])-1)))
      control$start.val[[i]] <- matrix(tmp, nrow = p, ncol = ncol(Z$Zmod[[i]]))
      colnames(control$start.val[[i]]) <- colnames(Z$Zmod[[i]]); rownames(control$start.val[[i]]) <- paste0("Y", 1:p)
      rm(tmp) }
    b <- control$start.val }
  # Constraints matrices (for each parameter)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(control$iden.res)){
    if(control$verbose) cat("\n Argument 'control$iden.res' not supplied: Using recursive identification restriction.")
    control$iden.res <- "recursive" }
  if(is.character(control$iden.res)){
    if(!control$iden.res %in% c("orthogonal", "eiv", "recursive")) stop("Argument 'control$iden.res' should be one of c(orthogonal, eiv, recursive).")
    if(control$iden.res == "orthogonal"){
      rb <- b
      for(i in 1:length(b)){
        b[[i]][,!grepl("(Intercept)", colnames(b[[i]]), fixed = T)] <- svd(b[[i]][,!grepl("(Intercept)", colnames(b[[i]]), fixed = T)])$u
        rb[[i]] <- !is.na(rb[[i]]) } }
    if(control$iden.res == "eiv"){
      rb <- b
      for(i in 1:length(b)){
        for(r in 1:q){ b[[i]][r, !colnames(b[[i]]) %in% c("(Intercept)", paste0("Z",r))] <- 0 }
        rb[[i]] <- b[[i]] != 0 } }
    if(control$iden.res == "recursive"){
      rb <- b
      for(i in 1:length(b)){
        b[[i]][,!colnames(b[[i]]) %in% c("(Intercept)")][upper.tri(b[[i]][,!colnames(b[[i]]) %in% c("(Intercept)")])] <- 0
        rb[[i]] <- b[[i]] != 0 } }
    b <-  list(b = b, rb = rb)
  } else {
    if(!is.list(control$iden.res)) stop("Argument 'control$iden.res' should be a list with element(s), each of the type c('parameter',item,'restricted variable',value)")
    b <- rmat(control$iden.res,b) }
  # Sign
  # ~~~~
  for(r in names(b$b)){
   for(j in which(grepl("Z",colnames(b$b[[r]])))){
     if(b$b[[r]][,j][b$b[[r]][,j] != 0][1] < 0) b$b[[r]][,j] <- -b$b[[r]][,j]
   }
  }
  return(b)
}