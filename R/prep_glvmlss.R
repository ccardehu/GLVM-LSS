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

pr_control <- function(control,q,...){

  con <- list("EM.maxit" = 20, "DM.maxit" = 500, "R.maxit" = 50,
              "nQP" = 25, "Fisher" = T, "verbose" = T,
              "EM.useGD" = F, "LR" = 1e-3,
              "start.b" = NULL, "start.R" = NULL, "restr.b" = NULL, "restr.R" = NULL,
              "estim.R" = T, "sign.restr.b" = NULL, "updt.R.every" = 1,
              "corLV" = F, "tolerance" = 1e-4) # sqrt(.Machine$double.eps)
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(namc[!namc %in% namC]) > 0)
    warning("Unknown names in control: ", paste(namc[!namc %in% namC], collapse = ", "))
  # Sanity checks:
  if(!is.null(con$start.R)){
    if(is.diag(con$start.R)) con$corLV <- F else con$corLV <- T
  }
  if(q == 1){
    con$estim.R <- F
    con$corLV <- F
  } else {
    con$nQP <- 15
  }
  if(!is.integer(con$updt.R.every)) con$updt.R.every <- as.integer(round(con$updt.R.every))
  return(con)
}

pr_q <- function(form){
  lvar <- unique(unlist(lapply(1:length(form), function(i) all.vars(form[[i]]))))
  lvar <- grep("Z", lvar, fixed = T, value = T)
  q <- length(lvar)
  return(q)
}

pr_bR <- function(Y,family,form,control,q){

  nostartb <- is.null(control$start.b)
  norestrb <- is.null(control$restr.b)
  nostartR <- is.null(control$start.R)
  norestrR <- is.null(control$restr.R)

  if(all(!nostartb,!norestrb,!nostartR,!norestrR)){
    bL <- rb2b(control$start.b, control$restr.b)
    RL <- rR2R(control$start.R, control$restr.R)
    if(bL$irb + RL$irR < q^2) stop(paste0("\n Not enough identification restrictions (",bL$irb," in control$restr.b + ",RL$irR," in control$restr.R)"))
    return(list(b = bL$b, rb = bL$rb, R = RL$R, rR = RL$rR))
  }

  if(all(nostartb,norestrb,nostartR,norestrR)){
    if(!control$corLV){
      control$start.R <- control$restr.R <- diag(q)
    } else {
      control$start.R <- matrix(0.05,q,q); diag(control$start.R) <- 1
      control$restr.R <- matrix(NA,q,q); diag(control$restr.R) <- 1
    }
    if(control$verbose) cat("\n Defining starting values (PCA + gamlss; missing control$start.b) ...")
    b <- suppressWarnings(b_ini(Y,family,form))
    if(control$verbose) cat("\r Defining starting values (PCA + gamlss; missing control$start.b) ... Done!")
    if(q > 1){
      if(control$corLV){
        if(control$verbose) cat("\n Imposing errors-in-variables restrictions (missing control$restr.b; control$corLV == T).")
        for(r in 1:q){
          b[r, !colnames(b) %in% c("(Intercept)", paste0("Z",r)), 1] <- 0
        }
      } else {
        if(control$verbose) cat("\n Imposing recursive restrictions (missing control$restr.b; control$corLV == F).")
        b[,!colnames(b) %in% c("(Intercept)"),1][upper.tri(b[, !colnames(b) %in% c("(Intercept)"), 1])] <- 0
      }
    }
    control$start.b <- b
    control$restr.b <- genrb(b)
    bL <- rb2b(control$start.b, control$restr.b)
    RL <- rR2R(control$start.R, control$restr.R)
    if(bL$irb + RL$irR < q^2) stop(paste0("\n Not enough identification restrictions (",bL$irb," in control$restr.b + ",RL$irR," in control$restr.R)"))
    return(list(b = bL$b, rb = bL$rb, R = RL$R, rR = RL$rR))
  }

  if(!nostartR & !nostartb){
    if(norestrR) control$restr.R <- genrR(control$start.R)
    if(norestrb) control$restr.b <- genrb(control$start.b)
    bL <- rb2b(control$start.b, control$restr.b)
    RL <- rR2R(control$start.R, control$restr.R)
    if(bL$irb + RL$irR < q^2) stop(paste0("\n Not enough identification restrictions (",bL$irb," in control$restr.b + ",RL$irR," in control$restr.R)"))
    return(list(b = bL$b, rb = bL$rb, R = RL$R, rR = RL$rR))
  }

  if(nostartR & nostartb){
    if(!control$corLV){
      control$start.R <- diag(1.1,q)
    } else {
      control$start.R <- matrix(0.05,q,q); diag(control$start.R) <- 1.1
    }
    if(norestrR){
      diag(control$start.R) <- 1
      control$restr.R <- genrR(control$start.R)
    }
    if(control$verbose) cat("\n Defining starting values (PCA + gamlss; missing control$start.b) ...")
    b <- suppressWarnings(b_ini(Y,family,form))
    if(control$verbose) cat("\r Defining starting values (PCA + gamlss; missing control$start.b) ... Done!")
    if(norestrb & q > 1){
      if(control$corLV){
        if(control$verbose) cat("\n Imposing errors-in-variables restrictions (missing control$restr.b; control$corLV == T).")
        for(r in 1:q){
          b[r, !colnames(b) %in% c("(Intercept)", paste0("Z",r)), 1] <- 0
        }
      } else {
        if(control$verbose) cat("\n Imposing recursive restrictions (missing control$restr.b; control$corLV == F).")
        b[,!colnames(b) %in% c("(Intercept)"),1][upper.tri(b[, !colnames(b) %in% c("(Intercept)"), 1])] <- 0
      }
      control$restr.b <- genrb(b)
    }
    control$start.b <- b
    bL <- rb2b(control$start.b, control$restr.b)
    RL <- rR2R(control$start.R, control$restr.R)
    if(bL$irb + RL$irR < q^2) stop(paste0("\n Not enough identification restrictions (",bL$irb," in control$restr.b + ",RL$irR," in control$restr.R)"))
    return(list(b = bL$b, rb = bL$rb, R = RL$R, rR = RL$rR))
  }

  if(nostartR & !nostartb){
    if(!control$corLV){
      control$start.R <- diag(1.1,q)
    } else {
      control$start.R <- matrix(0.05,q,q); diag(control$start.R) <- 1.1
    }
    if(norestrR){
      diag(control$start.R) <- 1
      control$restr.R <- genrR(control$start.R)
    }
    if(norestrb) control$restr.b <- genrb(control$start.b)
    bL <- rb2b(control$start.b, control$restr.b)
    RL <- rR2R(control$start.R, control$restr.R)
    if(bL$irb + RL$irR < q^2) stop(paste0("\n Not enough identification restrictions (",bL$irb," in control$restr.b + ",RL$irR," in control$restr.R)"))
    return(list(b = bL$b, rb = bL$rb, R = RL$R, rR = RL$rR))
  }

  if(!nostartR & nostartb){
    if(norestrR) control$restr.R <- genrR(control$start.R)
    if(control$verbose) cat("\n Defining starting values (PCA + gamlss; missing control$start.b) ...")
    b <- suppressWarnings(b_ini(Y,family,form))
    if(control$verbose) cat("\r Defining starting values (PCA + gamlss; missing control$start.b) ... Done!")
    if(norestrb & q > 1){
      if(control$corLV){
        if(control$verbose) cat("\n Imposing errors-in-variables restrictions (missing control$restr.b; control$corLV == T).")
        for(r in 1:q){
          b[r, !colnames(b) %in% c("(Intercept)", paste0("Z",r)), 1] <- 0
        }
      } else {
        if(control$verbose) cat("\n Imposing recursive restrictions (missing control$restr.b; control$corLV == F).")
        b[,!colnames(b) %in% c("(Intercept)"),1][upper.tri(b[, !colnames(b) %in% c("(Intercept)"), 1])] <- 0
      }
      control$restr.b <- genrb(b)
    }
    control$start.b <- b
    bL <- rb2b(control$start.b, control$restr.b)
    RL <- rR2R(control$start.R, control$restr.R)
    if(bL$irb + RL$irR < q^2) stop(paste0("\n Not enough identification restrictions (",bL$irb," in control$restr.b + ",RL$irR," in control$restr.R)"))
    return(list(b = bL$b, rb = bL$rb, R = RL$R, rR = RL$rR))
  }
}

pr_bRsim <- function(control,q,p,form){

  nostartb <- is.null(control$start.b)
  norestrb <- is.null(control$restr.b)
  nostartR <- is.null(control$start.R)
  norestrR <- is.null(control$restr.R)

  if(nostartb) stop("\n Missing argument: control$start.b")
  colnames(control$start.b) <- c("(Intercept)", paste0("Z",1:q))
  rownames(control$start.b) <- paste0("Y",1:p)
  dimnames(control$start.b)[[3]] <- names(form)
  if(norestrb) control$restr.b <- genrb(control$start.b)

  if(nostartR){
    if(!control$corLV){
      if(control$verbose) cat("\n Missing control$start.R: R set to diag(q) (control$corLV == F).")
      control$start.R <- diag(1.1,q)
    } else {
      if(control$verbose) cat("\n Missing control$start.R: Off-diagonal entries in R set to 0.3 (control$corLV == T).")
      control$start.R <- matrix(0.3,q,q); diag(control$start.R) <- 1.1
    }
    if(norestrR){
      diag(control$start.R) <- 1
      control$restr.R <- genrR(control$start.R)
    }
  } else {
    if(norestrR) control$restr.R <- genrR(control$start.R)
  }
  bL <- rb2b(control$start.b, control$restr.b)
  RL <- rR2R(control$start.R, control$restr.R)
  if(bL$irb + RL$irR < q^2) stop(paste0("\n Not enough identification restrictions (",bL$irb," in control$restr.b and ",RL$irR," in control$restr.R)"))
  colnames(RL$R) <- rownames(RL$R) <- paste0("Z",1:q)
  dimnames(RL$rR) <- dimnames(RL$R)
  return(list(b = bL$b, rb = bL$rb, R = RL$R, rR = RL$rR))
}
