pr_form <- function(mu.eq, sg.eq, nu.eq, ta.eq){

  is.formula <- function(x){ inherits(x,"formula") }
  i2formula <- function(x){
    if(!is.formula(x) & !is.null(x)) x <- as.formula(x)
    return(x)
  }
  form <- lapply(list(mu.eq, sg.eq, nu.eq, ta.eq), i2formula)
  names(form) <- c("mu","sg","nu","ta")
  for(i in names(form)){ if(is.null(form[[i]])) form[[i]] <- NULL }
  return(form)
}

pr_control <- function(control, ...){

  con <- list("EM.maxit" = 0, "DM.maxit" = 500, "R.maxit" = 50,
              "nQP" = 25, "Fisher" = T, "verbose" = T,
              "EM.useGD" = F, "LR" = 0.001,
              "start.b" = NULL, "start.R" = NULL, "restr" = NULL,
              "corLV" = F, "tolerance" = sqrt(.Machine$double.eps))
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(namc[!namc %in% namC]) > 0)
    warning("Unknown names in control: ", paste(namc[!namc %in% namC], collapse = ", "))
  return(con)
}

pr_q <- function(form){
  lvar <- unique(unlist(lapply(1:length(form), function(i) all.vars(form[[i]]))))
  lvar <- grep("Z", lvar, fixed = T, value = T)
  q <- length(lvar)
  return(q)
}

pr_b <- function(Y,family,form,control,q){
  if(!is.null(control$start.b)){
    b <- control$start.b
    if(!is.null(control$restr)) b <- b*control$restr
    return(b)
  } else {
    if(control$verbose) cat("\n Defining starting values (gamlss + PCA; missing control$start.b) ...")
    b <- suppressWarnings(b_ini(Y,family,form))
    if(control$verbose) cat("\r Defining starting values (gamlss + PCA; missing control$start.b) ... Done!")
    if(!is.null(control$restr)){
      b <- b*control$restr
      return(b)
    }
    if(control$corLV){
      if(control$verbose & q != 1) cat("\n Imposing errors-in-variables restrictions (missing control$restr; control$corLV == T).")
      for(r in 1:q){
        b[r, !colnames(b) %in% c("(Intercept)", paste0("Z",r)), 1] <- 0
      }
      return(b)
    } else {
      if(control$verbose & q != 1) cat("\n Imposing recursive restrictions (missing control$restr; control$corLV == F).")
      b[,!colnames(b) %in% c("(Intercept)"),1][upper.tri(b[, !colnames(b) %in% c("(Intercept)"), 1])] <- 0
      return(b)
    }
  }
}

pr_bsim <- function(control,q,p,form){
  if(is.null(control$start.b)) stop("\n Initial values not reported (missing argument: control$start.b)")
  b <- control$start.b
  colnames(b) <- c("(Intercept)", rep(paste0("Z",1:q)))
  rownames(b) <- paste0("Y",1:p)
  dimnames(b)[[3]] <- names(form)
  if(q == 1){
    return(b)
  } else {
    if(!is.null(control$restr)){
      b <- b*control$restr
      return(b)
    } else {
      if(control$corLV){
        if(control$verbose) cat("\n Imposing errors-in-variables restrictions (missing control$restr; control$corLV == T).")
        for(r in 1:q){
          b[r, !colnames(b) %in% c("(Intercept)", paste0("Z",r)), 1] <- 0
        }
        return(b)
      } else {
        if(control$verbose) cat("\n Imposing recursive restrictions (missing control$restr; control$corLV == F).")
        b[,!colnames(b) %in% c("(Intercept)"),1][upper.tri(b[, !colnames(b) %in% c("(Intercept)"), 1])] <- 0
        return(b)
      }
    }
  }
}
