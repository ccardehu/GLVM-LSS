prep_cont <- function(control,q,p,...){
  # Control settings
  # ~~~~~~~~~~~~~~~~
  con <- list(EM_iter = 30, EM_use2d = T, iter.lim = 300, DirectMaxFlag = T,
              EM_appHess = F, EM_lrate = 0.001, est.ci = "Approximate",
              solver = "nlminb", start.val = NULL, mat.info = "Hessian", lazytrust = F,
              iden.res = NULL, tol = sqrt(.Machine$double.eps), tolb = 1e-4,
              corr.lv = FALSE, Rz = NULL, var.lv = rep(1,q),
              nQP = if(q == 1) 40 else { if(q == 2) ifelse(p <= 10, 18, 25) else 10 },
              verbose = FALSE, autoL_iter = 30, f.scores = F, seed = 1234,
              penalty = "none", lambda = NULL, w.alasso = NULL, gamma = NULL, a = NULL)
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(namc[!namc %in% namC]) > 0)
    warning("Unknown names in control: ", paste(namc[!namc %in% namC], collapse = ", "))
  return(con)
}

prep_form <- function(mu.eq, sg.eq, nu.eq, ta.eq){
  # Formula settings
  # ~~~~~~~~~~~~~~~~
  is.formula <- function(x){ inherits(x,"formula") }
  if(!is.formula(mu.eq) & !is.null(mu.eq)) mu.eq <- as.formula(mu.eq)
  if(!is.formula(sg.eq) & !is.null(sg.eq)) sg.eq <- as.formula(sg.eq)
  if(!is.formula(nu.eq) & !is.null(nu.eq)) nu.eq <- as.formula(nu.eq)
  if(!is.formula(ta.eq) & !is.null(ta.eq)) ta.eq <- as.formula(ta.eq)
  form <- list("mu" = mu.eq, "sigma" = sg.eq, "nu" = nu.eq, "tau" = ta.eq)
  for(i in names(form)){ if(is.null(form[[i]])) form[[i]] <- NULL } 
  return(form)
}

prep_fam <- function(family){
  # Family settings
  # ~~~~~~~~~~~~~~~
  if(!is.list(family)) stop("Argument `family' must be a list of valid distributions.")
  if(!all(sapply(1:length(family), function(i) class(family[[i]]) == "dist_glvmlss"))){
    stop("Arguments in `family' must be of class `dist_glvmlss'.")
  }
  return(family)
}

prep_Z <- function(form){
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

prep_ghq <- function(nQP,form,Rz){
  # GHQ setting
  # ~~~~~~~~~~~
  # Rz <- control$Rz
  # if(control$corr.lv == T){ 
  #   if(is.null(Rz)){ Rz <- diag(q); Rz[lower.tri(Rz)] <- Rz[upper.tri(Rz)] <- 0.5 }
  #   } else Rz <- diag(q)
  return(mvghQ(n = nQP, formula = form, psi = Rz))
}

prep_stva <- function(control,form,ghQ,Y,p,famL,q){
  # Starting betas
  # ~~~~~~~~~~~~~~
  if(!is.null(control$start.val)){
    if(length(control$start.val) == 1 & is.character(control$start.val) & any(control$start.val == "random")){
      b <- vector(mode = "list", length = length(form))
      names(b) <- names(form)
      tmpnames <- c("(Intercept)", paste0("Z", 1:q))
      set.seed(control$seed)
      for(i in names(b)){
        size_ <- p*(sum(attr(terms(form[[i]]),"order")) + 1)
        b[[i]] <- matrix(runif(size_) * sample(c(1,-1), size = size_, replace = T), nrow = p)
        rownames(b[[i]]) <- colnames(Y$Y)
        colnames(b[[i]]) <- tmpnames[1:ncol(b[[i]])]
      }
    } else {
      if(!is.list(control$start.val)) stop("\n Provided starting values should be a list with dim p*(q+intercept) loading matrices for mu, sigma, nu, tau")
      if(length(control$start.val) != length(form)) stop("\n Number of matrices in starting values should match number of parameters mu, sigma, nu, tau")
      names(control$start.val) <- names(form)
      for(i in names(control$start.val)) {
        colnames(control$start.val[[i]]) <- colnames(ghQ$out[[i]]); rownames(control$start.val[[i]]) <- colnames(Y$Y)
        if(length(control$start.val[[i]]) != p*ncol(ghQ$out[[i]])) stop("\n Provided starting values is not lenght p*(q+intercept), revise dimension")
        if(!is.matrix(control$start.val[[i]])) control$start.val[[i]] <- matrix(control$start.val[[i]], nrow = p, ncol = ncol(ghQ$out[[i]]))
      }
      b <- control$start.val
    }
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
        penb[[i]][1,] <- F }
      b <-  list(b = b, rb = rb, penb = penb)
      return(b) }
    if(q != 1 & !control$corr.lv){
      if(control$verbose) cat("\n Argument 'control$iden.res' not supplied: Using recursive identification restriction (uncorrelated factors).")
      control$iden.res <- "recursive" }
    if(q != 1 & control$corr.lv){
      if(control$verbose) cat("\n Argument 'control$iden.res' not supplied: Using errors-in-variable identification restriction (correlated factors).")
      control$iden.res <- "eiv" } }
  if(is.character(control$iden.res)){
    if(!control$iden.res %in% c("orthogonal", "eiv", "recursive", "none")) stop("Argument 'control$iden.res' should be one of c(orthogonal, eiv, recursive, none).")
    if(control$iden.res == "orthogonal"){
      rb <- penb <- b
      b <- fixOrthob(b)
      for(i in 1:length(b)){
        rb[[i]] <- penb[[i]] <- !is.na(rb[[i]])
        penb[[i]][,grepl("(Intercept)", colnames(penb[[i]]), fixed = T)] <- F
      }
    }
    if(control$iden.res == "none"){
      rb <- penb <- b
      for(i in 1:length(b)){
        rb[[i]] <- penb[[i]] <- !is.na(rb[[i]])
        penb[[i]][,grepl("(Intercept)", colnames(penb[[i]]), fixed = T)] <- F
      }
    }
    if(control$iden.res == "eiv"){
      rb <- penb <- b
      for(i in names(b)){
        if(i == "mu") for(r in 1:q){ b[[i]][r, !colnames(b[[i]]) %in% c("(Intercept)", paste0("Z",r))] <- 0 }
        rb[[i]] <- penb[[i]] <- b[[i]] != 0 & !is.na(b[[i]])
        penb[[i]][,grepl("(Intercept)", colnames(penb[[i]]), fixed = T)] <- F 
        if(i == "mu") diag(penb[[i]][,!grepl("(Intercept)", colnames(penb[[i]]), fixed = T)]) <- F # Not penalise loadings in diagonal (to give direction)
      }
    }
    if(control$iden.res == "recursive"){
      rb <- penb <- b
      for(i in names(b)){
        if(i == "mu") b[[i]][,!colnames(b[[i]]) %in% c("(Intercept)")][upper.tri(b[[i]][,!colnames(b[[i]]) %in% c("(Intercept)")])] <- 0
        rb[[i]] <- penb[[i]] <- b[[i]] != 0 & !is.na(b[[i]])
        penb[[i]][,grepl("(Intercept)", colnames(penb[[i]]), fixed = T)] <- F
        if(i == "mu") diag(penb[[i]][,!grepl("(Intercept)", colnames(penb[[i]]), fixed = T)]) <- F # Not penalise loadings in diagonal (to give direction)
      }
    }
    b <-  list(b = b, rb = rb, penb = penb)
  } else {
    if(q == 1){
      if(control$verbose) cat("\n Argument 'control$iden.res' not needed for one-factor models (q == 1).")
      rb <- penb <- b
      for(i in names(b)){
        rb[[i]] <- penb[[i]] <- !is.na(rb[[i]])
        penb[[i]][,grepl("(Intercept)", colnames(penb[[i]]), fixed = T)] <- F 
        penb[[i]][1,2] <- F }
      b <-  list(b = b, rb = rb, penb = penb)
      return(b) }
    if(!is.list(control$iden.res)) stop("Argument 'control$iden.res' should be a list with element(s), each of the type c('parameter',item,'restricted variable',value)")
    b <- rmat(control$iden.res,b) }
  # Sign
  # ~~~~
  # for(j in paste0("Z",1:q)){
  #   if(j %in% colnames(b$b[[1]])){
  #     if(b$b[[1]][,j][b$b[[1]][,j] != 0][1] < 0){
  #       for(r in names(b$b)){
  #         if(j %in% colnames(b$b[[r]])) b$b[[r]][,j] <- -b$b[[r]][,j] }
  #     }
  #   }
  # }
  return(b)
}

sim_cont <- function(control,q,...){
  # Control settings
  # ~~~~~~~~~~~~~~~~
  con <- list(muZ = rep(0,q), Rz = diag(q), start.val = NULL, iden.res = NULL,
              verbose = FALSE)
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(namc[!namc %in% namC]) > 0)
    warning("Unknown names in control: ", paste(namc[!namc %in% namC], collapse = ", "))
  return(con)
}

sim_Z <- function(control,q,n,form){
  # if(sum(diag(control$Rz)) != q){ warning("Argument `control$Rz' standardised to a correlation matrix")
  #   control$Rz <- diag(1/sqrt(diag(control$Rz)))%*%control$Rz%*%diag(1/sqrt(diag(control$Rz)))
  # }
  Z <- as.data.frame(mvtnorm::rmvnorm(n, control$muZ, control$Rz))
  colnames(Z) <- paste0("Z", 1:q)
  Zmod <- vector("list", length = length(form)); names(Zmod) <- names(form)
  for(i in names(form)){ Zmod[[i]] <- as.data.frame(model.matrix(form[[i]],Z)) }
  return(list(Z = Z, Zmod = Zmod))
}

sim_stva <- function(control,p,Z,form,q){
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
    if(q == 1){
      rb <- b
      for(i in 1:length(b)){
        rb[[i]] <- !is.na(rb[[i]])
      }
      b <-  list(b = b, rb = rb)
    }
    if(q != 1 & (sum(control$Rz) == 2)){
      if(control$verbose) cat("\n Argument 'control$iden.res' not supplied: Using recursive identification restriction (uncorrelated factors).")
      control$iden.res <- "recursive" }
    if(q != 1 & !(sum(control$Rz) == 2)){
      if(control$verbose) cat("\n Argument 'control$iden.res' not supplied: Using errors-in-variable identification restriction (correlated factors).")
      control$iden.res <- "eiv" } 
  }
  if(is.character(control$iden.res)){
    if(!control$iden.res %in% c("orthogonal", "eiv", "recursive","none")) stop("Argument 'control$iden.res' should be one of c(orthogonal, eiv, recursive, none).")
    if(control$iden.res == "none"){
      rb <- b
      for(i in names(b)){
        rb[[i]] <- !is.na(rb[[i]]) } }
    if(control$iden.res == "orthogonal"){
      rb <- b
      b <- fixOrthob(b)
      for(i in names(b)){
        rb[[i]] <- !is.na(rb[[i]]) } }
    if(control$iden.res == "eiv"){
      rb <- b
      for(i in names(b)){
        if(i == "mu") for(r in 1:q){ b[[i]][r, !colnames(b[[i]]) %in% c("(Intercept)", paste0("Z",r))] <- 0 }
        rb[[i]] <- b[[i]] != 0 } }
    if(control$iden.res == "recursive"){
      rb <- b
      for(i in names(b)){
        if(i == "mu") b[[i]][,!colnames(b[[i]]) %in% c("(Intercept)")][upper.tri(b[[i]][,!colnames(b[[i]]) %in% c("(Intercept)")])] <- 0
        rb[[i]] <- b[[i]] != 0 } }
    b <-  list(b = b, rb = rb)
    flagSign <- T
  } else {
    if(q == 1){
      if(control$verbose) cat("\n Argument 'control$iden.res' not needed for one-factor models (q == 1). Z1 is Standard Normal.")
      # rb <- b
      # for(i in names(b)){
      #   rb[[i]] <- !is.na(rb[[i]])
      # }
      # b <-  list(b = b, rb = rb)
    } else {
      if(!is.list(control$iden.res)) stop("Argument 'control$iden.res' should be a list with element(s), each of the type c('parameter',item,'restricted variable',value)")
      b <- rmat(control$iden.res,b)
    }
    flagSign <- F}
  # Sign
  # ~~~~
  if(flagSign){
    for(r in names(b$b)){
      for(j in which(grepl("Z",colnames(b$b[[r]])))){
        if(b$b[[r]][,j][b$b[[r]][,j] != 0][1] < 0) b$b[[r]][,j] <- -b$b[[r]][,j]
      }
    }
  }
  return(b)
}
