
grph.fun1 <- function(item,mod,qtl = c(seq(0.1,0.9,by = 0.2)) ){ #fFun + gFun

# Goal: To compute Expected values, SD, fitted parameters
#       and quantiles of distribution for item i
# Input: item (item), model (mod), quantiles (qtl)
# Output: List with E(i), sd(i), theta(i) and Q_{qtl}(i)
# Testing: item = 1; mod = testmod; qtl = c(seq(0.1,0.9,by = 0.2))
 
if(is.null(qtl)) qtl.ret <- F else qtl.ret <- T # missing(qtl) || 

Z <- mod$ghQ$out

if(mod$fam[item] == "normal"){
 mu = c(unname(as.matrix(Z$mu)%*%matrix(mod$b$mu[item,])))
 sigma = c(unname(exp(as.matrix(Z$sigma)%*%matrix(mod$b$sigma[item,]))))
 EY = mu
 SY = sigma
 if(qtl.ret){
  qM = matrix(nrow = length(mu), ncol = length(qtl)); colnames(qM) <- paste0("q",qtl*100)
  for(j in seq_along(qtl)){ qM[,j] <- qnorm(qtl[j],mu,sigma) } }
 theta = list(mu = mu, sigma = sigma)
}
if(mod$fam[item] == "lognormal"){
 mu = c(unname(as.matrix(Z$mu)%*%matrix(mod$b$mu[item,])))
 sigma = c(unname(exp(as.matrix(Z$sigma)%*%matrix(mod$b$sigma[item,]))))
 EY = exp(mu + 0.5*sigma^2)
 SY = sqrt((exp(sigma^2)-1)*exp(2*mu+sigma^2))
 if(qtl.ret){
  qM = matrix(nrow = length(mu), ncol = length(qtl)); colnames(qM) <- paste0("q",qtl*100)
  for(j in seq_along(qtl)){ qM[,j] <- qlnorm(qtl[j],mu,sigma) } }
 theta = list(mu = mu, sigma = sigma)
}
if(mod$fam[item] == "poisson"){
 mu = c(unname(exp(as.matrix(Z$mu)%*%matrix(mod$b$mu[item,]))))
 EY = mu
 SY = sqrt(mu)
 if(qtl.ret){
  qM = matrix(nrow = length(mu), ncol = length(qtl)); colnames(qM) <- paste0("q",qtl*100)
  for(j in seq_along(qtl)){ qM[,j] <- qpois(qtl[j],mu) } }
 theta = list(mu = mu)
}
if(mod$fam[item] == "gamma"){
 mu = c(unname(exp(as.matrix(Z$mu)%*%matrix(mod$b$mu[item,]))))
 sigma = c(unname(exp(as.matrix(Z$sigma)%*%matrix(mod$b$sigma[item,]))))
 EY = mu*sigma
 SY = sqrt(mu*sigma^2)
 if(qtl.ret){
  qM = matrix(nrow = length(mu), ncol = length(qtl)); colnames(qM) <- paste0("q",qtl*100)
  for(j in seq_along(qtl)){ qM[,j] <- qgamma(qtl[j], shape = mu, scale = sigma) } }
 theta = list(mu = mu, sigma = sigma)
}
if(mod$fam[item] == "binomial"){
 mu = c(unname(probs(as.matrix(Z$mu)%*%matrix(mod$b$mu[item,]))))
 EY = mu
 SY = sqrt(mu*(1-mu))
 if(qtl.ret){
  qM = matrix(nrow = length(mu), ncol = length(qtl)); colnames(qM) <- paste0("q",qtl*100)
  for(j in seq_along(qtl)){ qM[,j] <- qbinom(qtl[j], 1, mu) } }
 theta = list(mu = mu)
}
if(mod$fam[item] == "ZIpoisson"){
 mu = c(unname(exp(as.matrix(Z$mu)%*%matrix(mod$b$mu[item,]))))
 sigma = c(unname(probs(as.matrix(Z$sigma)%*%matrix(mod$b$sigma[item,]))))
 qZIpoisson <- function (p,mu,sigma){ # https://stats.stackexchange.com/questions/463844/is-there-an-equation-for-the-median-and-percentile-of-a-zero-inflated-poisson
  ly <- max(length(p), length(mu), length(sigma))
  p <- rep(p, length = ly); sigma <- rep(sigma, length = ly); mu <- rep(mu, length = ly)
  pnew <- ((p - sigma)/(1 - sigma)) - (1e-07)
  pnew <- ifelse(pnew > 0, pnew, 0)
  q <- ifelse(p > sigma*(1-sigma)*dpois(0,mu), qpois(pnew, mu), 0)
  return(q)
 }
 EY = (1-sigma)*mu
 SY = sqrt(mu*(1-sigma)*(1+mu*sigma))
 if(qtl.ret){
  qM = matrix(nrow = length(mu), ncol = length(qtl)); colnames(qM) <- paste0("q",qtl*100)
  for(j in seq_along(qtl)){ qM[,j] <- qZIpoisson(qtl[j], mu, sigma) } }
 theta = list(mu = mu, sigma = sigma)
}
if(qtl.ret){
 return(list(mean = EY, sd = SY, theta = theta, quant = as.data.frame(qM))) } else {
 return(list(mean = EY, sd = SY, theta = theta)) }
}

grph.fun2 <- function(y,item,mod,sel){ #f2Fun # Only for 3D graphs
   
# Goal: To compute density (using fitted parameters) for item i @ certain selection
# Input: item (item), model (mod), selection
# Output: Density f(Y_i), for a selection of theta[sel]
# Testing: y = y; item = 1; mod = testmod; sel = sel1

g <- grph.fun1(item,mod)
if(mod$fam[item] == "normal"){ fyz <- dnorm(y, g$theta$mu[sel], g$theta$sigma[sel]) }
if(mod$fam[item] == "lognormal"){ fyz <- dlnorm(y, g$theta$mu[sel], g$theta$sigma[sel]) }
if(mod$fam[item] == "poisson"){ fyz <- dpois(ceiling(y), g$theta$mu[sel]) }
if(mod$fam[item] == "gamma"){ fyz <- dgamma(y,shape = g$theta$mu[sel], scale = g$theta$sigma[sel]) }
if(mod$fam[item] == "binomial"){ fyz <- dbinom(round(y), 1, prob = g$theta$mu[sel]) }  
if(mod$fam[item] == "ZIpoisson"){
dZIpoisson <- function(Y,mu,sigma,log=T){
 u <- as.numeric(Y == 0); lf <- u*c(log(sigma + (1-sigma)*exp(-mu))) + (1-u)*(c(log(1-sigma)) - c(mu) + Y*c(log(mu)) - lfactorial(Y))
 if(log == T) return(lf) else return(exp(lf))
}; fyz <- dZIpoisson(ceiling(y), g$theta$mu[sel], g$theta$sigma[sel], log = F) }
return(fyz)    
}

grph.fun3 <- function(item,mod,qtl = c(seq(0.1,0.9,by = 0.2)), mar = q){ #fFun + gFun

# Goal: To compute (MARGINAL) Expected values, SD, fitted parameters
#       and quantiles of distribution for item i
# Input: item (item), model (mod), quantiles (qtl), margin to leave (integrate out all other factors)
# Output: List with (MARGINAL) E(i), sd(i), theta(i) and Q_{qtl}(i)
# Testing: item = 1; mod = testmod; qtl = c(seq(0.1,0.9,by = 0.2)); mar = 1
 
if(is.null(qtl)) qtl.ret <- F else qtl.ret <- T # missing(qtl) || 
lvar <- unique(unlist(lapply(1:length(mod$form), function(i) all.vars(mod$form[[i]]))))
lvar <- grep("Z", lvar, fixed = T, value = T)
lvar <- grep(paste0("Z",mar),lvar, invert = T, fixed = T, value = TRUE)

Zn <- list("x" = fastGHQuad::gaussHermiteData(mod$ghQ$n)$x, "w" = fastGHQuad::gaussHermiteData(mod$ghQ$n)$w)
Zn$w <- Zn$w * (2*pi)^(-1/2)*exp(Zn$x^2/2)
EY <- SY <- array(NA,dim = c(length(Zn$x),1,length(Zn$x)))
qM <- array(NA,dim = c(length(Zn$x),length(qtl),length(Zn$x)))
colnames(qM) <- paste0("q",qtl*100)

# EY <- SY <- array(NA,dim = c(length(mod$ghQ$weights),1,length(mod$ghQ$weights)))
# qM <- array(NA,dim = c(length(mod$ghQ$weights),length(qtl),length(mod$ghQ$weights)))
# colnames(qM) <- paste0("q",qtl*100)
# rma <- vector(length = length(mod$ghQ$points[,mar]))
# rma[] <- fastGHQuad::gaussHermiteData(mod$ghQ$n)$w

for(r in 1:nrow(EY)){
 Zr <- NULL
 Zs <- mod$ghQ$points # 
 Z <- mod$ghQ$out
 for(a in names(mod$ghQ$out)){
  # Zr[[a]] <- Z[[a]][r,grepl(lvar,colnames(Z[[a]]))]
  # Z[[a]][,grepl(lvar,colnames(Z[[a]]))] <- c(Zr[[a]])
  Zr[[a]] <- Zn$x[r]# Zs[r,lvar]
  Zs[,lvar] <- Zr[[a]]
  Z[[a]] <- unique(as.data.frame(model.matrix(mod$formula[[a]],as.data.frame(Zs)))) }
  if(mod$fam[item] == "normal"){
   mu = c(unname(as.matrix(Z$mu)%*%matrix(mod$b$mu[item,])))
   sigma = c(unname(exp(as.matrix(Z$sigma)%*%matrix(mod$b$sigma[item,]))))
   EY[,,r] = mu
   SY[,,r] = sigma
   if(qtl.ret){
    for(j in seq_along(qtl)){ qM[,j,r] <- qnorm(qtl[j],mu,sigma) } }
  # theta = list(mu = mu, sigma = sigma)
  }
  if(mod$fam[item] == "lognormal"){
   mu = c(unname(as.matrix(Z$mu)%*%matrix(mod$b$mu[item,])))
   sigma = c(unname(exp(as.matrix(Z$sigma)%*%matrix(mod$b$sigma[item,]))))
   EY[,,r] = exp(mu + 0.5*sigma^2)
   SY[,,r] = sqrt((exp(sigma^2)-1)*exp(2*mu+sigma^2))
   if(qtl.ret){
    for(j in seq_along(qtl)){ qM[,j,r] <- qlnorm(qtl[j],mu,sigma) } }
   # theta = list(mu = mu, sigma = sigma)
  }
  if(mod$fam[item] == "poisson"){
   mu = c(unname(exp(as.matrix(Z$mu)%*%matrix(mod$b$mu[item,]))))
   EY[,,r] = mu
   SY[,,r] = sqrt(mu)
   if(qtl.ret){
    for(j in seq_along(qtl)){ qM[,j] <- qpois(qtl[j],mu) } }
   # theta = list(mu = mu)
  }
  if(mod$fam[item] == "gamma"){
   mu = c(unname(exp(as.matrix(Z$mu)%*%matrix(mod$b$mu[item,]))))
   sigma = c(unname(exp(as.matrix(Z$sigma)%*%matrix(mod$b$sigma[item,]))))
   EY[,,r] = mu*sigma
   SY[,,r] = sqrt(mu*sigma^2)
   if(qtl.ret){
    for(j in seq_along(qtl)){ qM[,j,r] <- qgamma(qtl[j], shape = mu, scale = sigma) } }
   # theta = list(mu = mu, sigma = sigma)
  }
  if(mod$fam[item] == "binomial"){
   mu = c(unname(probs(as.matrix(Z$mu)%*%matrix(mod$b$mu[item,]))))
   EY[,,r] = mu
   SY[,,r] = sqrt(mu*(1-mu))
   if(qtl.ret){
    for(j in seq_along(qtl)){ qM[,j,r] <- qbinom(qtl[j], 1, mu) } }
   # theta = list(mu = mu)
  }
  if(mod$fam[item] == "ZIpoisson"){
   mu = c(unname(exp(as.matrix(Z$mu)%*%matrix(mod$b$mu[item,]))))
   sigma = c(unname(probs(as.matrix(Z$sigma)%*%matrix(mod$b$sigma[item,]))))
   qZIpoisson <- function (p,mu,sigma){ # https://stats.stackexchange.com/questions/463844/is-there-an-equation-for-the-median-and-percentile-of-a-zero-inflated-poisson
     ly <- max(length(p), length(mu), length(sigma))
     p <- rep(p, length = ly); sigma <- rep(sigma, length = ly); mu <- rep(mu, length = ly)
     pnew <- ((p - sigma)/(1 - sigma)) - (1e-07)
     pnew <- ifelse(pnew > 0, pnew, 0)
     q <- ifelse(p > sigma*(1-sigma)*dpois(0,mu), qpois(pnew, mu), 0)
     return(q) }
   EY[,,r] = (1-sigma)*mu
   SY[,,r] = sqrt(mu*(1-sigma)*(1+mu*sigma))
   if(qtl.ret){
    for(j in seq_along(qtl)){ qM[,j,r] <- qZIpoisson(qtl[j], mu, sigma) } }
   # theta = list(mu = mu, sigma = sigma)
  }
  EY[,,r] <- (EY[,,r]*c(Zn$w)[r])
  SY[,,r] <- (SY[,,r]*c(Zn$w)[r])
  qM[,,r] <- (qM[,,r]*c(Zn$w)[r])
 }

EY <- rowSums(EY,dims = 2); SY <- rowSums(SY,dims = 2); qM <- rowSums(qM, dims = 2)

if(qtl.ret){
 return(list(mean = EY, sd = SY, quant = as.data.frame(qM))) } else {
 return(list(mean = EY, sd = SY)) }
}

grph.fun4 <- function(fam){ # should be evaluated at i = 1:p; fam = fam[i]
  if(fam == "normal") lims <- c(-Inf,Inf)
  if(fam == "lognormal") lims <- c(0  + .Machine$double.eps, Inf)
  if(fam == "poisson") lims <- c(0,Inf)
  if(fam == "gamma") lims <- c(0  + .Machine$double.eps, Inf)
  if(fam == "binomial") lims <- c(0,1)
  if(fam == "ZIpoisson") lims <- c(0,Inf)
  # Add other parameters for other distributions
  return(lims)
}

splvm.plot <- function(item, mod, 
                       control.plot = list(sim.mod = simR,
                       plot.mean = T, plot.sd = F,
                       qtl = c(seq(0.1,0.9,by = 0.2)),
                       plot.3D = F, bw = F)){

# Goal: Plot semi-parametric LVM
# Input : item (item), model (mod),
#         form (list of formulas for measurement eqs. for each parameter),
#         control (list of controls)
# Output: plot for estimated semi-parametric LVM (mod)
# Testing: item = 1; mod = testmod; control.plot = list(sim.mod = simR, plot.mean = T, plot.sd = F, qtl = c(seq(0.1,0.9,by = 0.2)), plot.3D = F)
  
if("scales" %in% rownames(installed.packages()) == F) install.packages("scales")

# item = 1; mod = ex1; morg = simR; plot.mean = T; plot.org = F # plot.sd = T; plot.addpoints = F; sep.plots = F, plot.3D = T # quant = c(0.05,0.2,0.4,0.6,0.8,0.95)

if(is.null(control.plot$plot.mean)) control.plot$plot.mean <- ifelse(is.null(control.plot$qtl),T,F) # missing(control.plot$qtl) || 
control.plot$plot.sd <- ifelse(is.null(control.plot$qtl),T,F) # missing(control.plot$qtl) || 
control.plot$plot.org <- ifelse(!is.null(control.plot$sim.mod),T,F) # !missing(control.plot$sim.mod) || 
control.plot$plot.quant <- ifelse(!is.null(control.plot$qtl),T,F) # !missing(control.plot$qtl) || 
qtl <- control.plot$qtl
fam <- mod$fam;
f.mod <- grph.fun1(item,mod)
dimp <- ncol(mod$ghQ$points); Zsc <- fscore(mod)

if(!control.plot$plot.3D){ par(mfrow=c((dimp+1)%/%2,ifelse(dimp == 1, 1, 2)))
 for(i in 1:dimp){
  if(control.plot$plot.org){
   zlab <- min(control.plot$sim.mod$Z[[paste0("Z",i)]])
   plot(control.plot$sim.mod$Z[[paste0("Z",i)]], mod$Y[[paste0("Y",item)]], pch = 16, col = scales::alpha("gray", 0.7),
        xlab = bquote(Z[.(i)]), ylab = bquote(Y[.(item)]))
   if(dimp != 1) f.modmg <- grph.fun3(item,mod,qtl,i) else f.modmg <- grph.fun1(item,mod,qtl)
   if(control.plot$plot.mean) lines(unique(mod$ghQ$points[,paste0("Z",i)]), f.modmg$mean, lwd = 2, col = 2)
   if(control.plot$plot.sd){
    lines(unique(mod$ghQ$points[,paste0("Z",i)]), f.modmg$mean + f.modmg$sd, lwd = 2, col = 2, lty = 3, type = "b", pch = ".")
    lines(unique(mod$ghQ$points[,paste0("Z",i)]), f.modmg$mean - f.modmg$sd, lwd = 2, col = 2, lty = 3, type = "b", pch = ".")
    mtext(side = 3, line = 1, at = zlab, adj = 0, cex = 0.8, text = "E(Y|z) +/- 2*SD")
   }
   if(control.plot$plot.quant){
    for(j in seq_along(control.plot$qtl)){ lines(unique(mod$ghQ$points[,paste0("Z",i)]), f.modmg$quant[,j], lwd = 2, col = 2, lty = 3, type = "b", pch = ".") }
      mtext(side = 3, line = 1, at = zlab, adj = 0, cex = 0.8, text = paste0("(Q: ", paste0(qtl, collapse = ", "),")"))
   }
  } else {
   zlab <- min(Zsc[[paste0("Z",i)]])
   plot(Zsc[[paste0("Z",i)]], mod$Y[[paste0("Y",item)]], pch = 16, col = scales::alpha("gray", 0.7),
        xlab = bquote(Z[.(i)]), ylab = bquote(Y[.(item)]))     
   if(dimp != 1) f.modmg <- grph.fun3(item,mod,qtl,i) else f.modmg <- grph.fun1(item,mod,qtl)
   if(control.plot$plot.mean) lines(mod$ghQ$points[,paste0("Z",i)], f.modmg$mean, lwd = 2, col = 2)
   if(control.plot$plot.sd){
    lines(mod$ghQ$points[,paste0("Z",i)], f.modmg$mean + f.modmg$sd, lwd = 2, col = 2, lty = 3, type = "b", pch = ".")
    lines(mod$ghQ$points[,paste0("Z",i)], f.modmg$mean - f.modmg$sd, lwd = 2, col = 2, lty = 3, type = "b", pch = ".")
    mtext(side = 3, line = 1, at = zlab, adj = 0, cex = 0.8, text = "E(Y|z) +/- 2*SD")
   }
   if(control.plot$plot.quant){
    for(j in seq_along(control.plot$qtl)){ lines(mod$ghQ$points[,paste0("Z",i)], f.modmg$quant[,j], lwd = 2, col = 2, lty = 3, type = "b", pch = ".") }
    mtext(side = 3, line = 1, at = zlab, adj = 0, cex = 0.8, text = paste0("(Q: ", paste0(qtl, collapse = ", "),")"))
   }
  }
  legend("topleft",legend=c("Fitted model"), col=c("red"), cex=0.8, lty = 1, inset = 0.02, bty = "n")
  if(dimp == 1){ mtext(side = 3, line = 2, at = zlab, adj = 0, cex = 1, paste0("Probability function P(Y",item," | Z",i,")")) } else {
     mtext(side = 3, line = 2, at = zlab, adj = 0, cex = 1, paste0("Marginal probability P(Y",item," | Z",i,")")) }
 }
} else { #else 3D
 if(dimp < 2){
  par(mfrow = c(1,1), mar = c(3, 3, 4, 2) + 0.1)
  # tmpx1 <- sort(unique(round(c(mod$ghQ$points),5)));
  tmpx1 <- sort(unique(mod$ghQ$points));
  if(control.plot$plot.org){
   sel1 <- (tmpx1 > min(control.plot$sim.mod$Z[[paste0("Z",dimp)]]) & tmpx1 < max(control.plot$sim.mod$Z[[paste0("Z",dimp)]]))
  } else {
   sel1 <- (tmpx1 > min(Zsc[[paste0("Z",dimp)]]) & tmpx1 < max(Zsc[[paste0("Z",dimp)]]))
  }
  tmpx <- tmpx1[sel1]      
  tmpy <- sort(unique(mod$Y[[paste0("Y",item)]]))
  vX <- seq(min(tmpx)-1,max(tmpx)+1,length = 2); vY <- seq(min(tmpy)-1,max(tmpy)+1,length = 2)
  zlim. <- min(max(grph.fun2(tmpy,item,mod))+0.2,1)
  
  res  <- persp(vX,vY,matrix(0,2,2),zlim=c(0,zlim.),theta= -57, phi = 30, expand = 0.9, ticktype = "detailed",
                box = T, xlab = paste0("\n \n Z",dimp), ylab = paste0("\n \n Y",item), zlab = paste0("\n \n f (Y",item," | Z1)"),
                nticks = 5, border = T, cex.axis = 0.8, las = 3)
  mtext(side = 3, line = 2, adj = 0, cex = 1, text = paste0("Fitted probability function f(Y",item," | Z1)"))
   
  if(control.plot$plot.org){
   C <- trans3d(control.plot$sim.mod$Z[[paste0("Z",dimp)]], mod$Y[[paste0("Y",item)]], rep(0,length(mod$Y[[paste0("Y",item)]])),res)
   points(C, pch = 20, col= scales::alpha("gray", 0.6), cex = 0.8)
  } else {
   C <- trans3d(Zsc[[paste0("Z",dimp)]], mod$Y[[paste0("Y",item)]],rep(0,length(mod$Y[[paste0("Y",item)]])),res)
   points(C, pch = 20, col= scales::alpha("gray", 0.6), cex = 0.8)   
  }
  f.modmg <- grph.fun1(item,mod,qtl)
  if(control.plot$plot.mean){ C <- trans3d(tmpx,unique(f.modmg$mean)[sel1],rep(0,length(tmpx)),res); lines(C,lwd=2, col = 2) }
  if(control.plot$plot.sd){
   C <- trans3d(tmpx,unique(f.modmg$mean + 2*f.modmg$sd)[sel1],rep(0,length(tmpx)),res); lines(C,lwd=2, lty = 3, col = 2)
   C <- trans3d(tmpx,unique(f.modmg$mean - 2*f.modmg$sd)[sel1],rep(0,length(tmpx)),res); lines(C,lwd=2, lty = 3, col = 2)
   C <- trans3d(c(tmpx,rev(tmpx)),c(unique(f.modmg$mean + 2*f.modmg$sd)[sel1],rev(unique(f.modmg$mean - 2*f.modmg$sd)[sel1])),rep(0,2*length(tmpx)),res)
   polygon(C,border=NA,col=scales::alpha("yellow", 0.5))
   mtext(side = 3, line = 1, adj = 0, cex = 0.8, text = "E(Y|z) +/- 2*SD")
  }
  # if(plot.dist == T){
  cuts <- 10
  vXa <- tmpx[round(seq(1,length(tmpx),length.out = cuts))]
   for(w in 1:cuts){
    if(length(tmpy) < 10) stp <- 10 else stp <- ifelse(round(length(tmpy)/10) < 10, 10, round(length(tmpy)/10))
    x <- rep(vXa[w],stp)
    liminf <- ifelse(f.modmg$quant[which(vXa[w]==tmpx1),paste0("q",100*min(control.plot$qtl))]-2*f.modmg$sd[which(vXa[w]==tmpx1)] < max(min(vY), grph.fun4(fam[[item]])[1]),
                     max(min(vY),grph.fun4(fam[[item]])[1]),f.modmg$quant[which(vXa[w]==tmpx1),paste0("q",100*min(control.plot$qtl))]-2*f.modmg$sd[which(vXa[w]==tmpx1)])
    limsup <- ifelse(f.modmg$quant[which(vXa[w]==tmpx1),paste0("q",100*max(control.plot$qtl))]+2*f.modmg$sd[which(vXa[w]==tmpx1)] > min(max(vY),grph.fun4(fam[[item]])[2]),
                     min(max(vY),grph.fun4(fam[[item]])[2]), f.modmg$quant[which(vXa[w]==tmpx1),paste0("q",100*max(control.plot$qtl))]+2*f.modmg$sd[which(vXa[w]==tmpx1)])
    z0 <- rep(0,stp)
     
    if(any(fam[[item]] == c("ZIpoisson","poisson","binom"))) {
     y <- round(seq(liminf,limsup, length = stp))
     z <- ifelse(grph.fun2(y,item,mod,w) > 1, 1, grph.fun2(y,item,mod,w))
     C <- trans3d(c(x,x),c(y,rev(y)),c(z,z0),res)
     for(m in 1:length(z0)){ lines(trans3d(x[m],y[m],c(z0[m],z[m]),res), col = scales::alpha("blue",0.4), lty = 3) }
     C <- trans3d(x,y,z0,res); lines(C,lty=3, col = scales::alpha("lightblue",0.6))
     C <- trans3d(x,y,z,res)
     points(C, pch = 16, col= scales::alpha("blue",0.4), cex = 0.5)
    } else {
     y <- seq(liminf,limsup, length = stp)
     z <- ifelse(grph.fun2(y,item,mod,which(vXa[w]==tmpx1)) > 1, 1, grph.fun2(y,item,mod,which(vXa[w]==tmpx1)))
     C <- trans3d(c(x,x),c(y,rev(y)),c(z,z0),res)
     polygon(C,border=NA,col=scales::alpha("lightblue",0.3))
     C <- trans3d(x,y,z0,res)
     lines(C,lty=3, col = scales::alpha("lightblue",0.6))
     C <- trans3d(x,y,z,res)
     lines(C,col=scales::alpha("blue",0.4)) }
     
   }; #legend("bottomleft",legend=c("Fitted distribution"), col= scales::alpha("blue",0.5), lty =1, inset = -0.5, bty = "n")
   # }
   
   if(control.plot$plot.quant){
    for(j in qtl){
       C <- trans3d(tmpx,f.modmg$quant[sel1,paste0("q",100*j)],rep(0,length(tmpx)),res)
       lines(C,lwd=2, lty = 3, col = 2)
       #C <- trans3d(c(tmpx,rev(tmpx)),c(c(fFun.mod$mean + fFun.mod$sM)[sel1],rev(c(fFun.mod$mean - fFun.mod$sM)[sel1])),rep(0,2*length(tmpx)),res)
       #polygon(Ck,border=NA,col=scales::alpha("yellow", 0.5))
    }
    mtext(side = 3, line = 1, adj = 0, cex = 0.8, text = paste0("(Q: ", paste0(qtl, collapse = ", "),")"))   
   }
}
 if(dimp == 2){
 if(control.plot$plot.org) par(mfrow = c(1,2),mar = c(3, 3, 4, 2) + 0.1) else par(mfrow = c(1,1),mar = c(3, 3, 4, 2) + 0.1)
 for(k in 1:dimp){
  tmpx <- tmpy <- sort(unique(round(c(mod$ghQ$points),5))); tmpz <- matrix(NA, nrow = length(tmpx), ncol = length(tmpy))
  f.mod <- grph.fun1(item,mod,control.plot$qtl)
  tmpB <- cbind(round(mod$ghQ$points,5), f.mod$mean)
  for(i in 1:length(tmpx)){
   for(j in 1:length(tmpy)){
       tmpz[i,j] <- tmpB[intersect(which(tmpB[,1] == tmpx[i]), which(tmpB[,2] == tmpy[j])),3]}}
  tmpzz <- tmpz[-1, -1] + tmpz[-1, -ncol(tmpz)] + tmpz[-nrow(tmpz), -1] + tmpz[-nrow(tmpz), -ncol(tmpz)]
  cutz <- cut(c(tmpzz), 10)
  if(control.plot$bw) color <- colorRampPalette(c("white"))(length(tmpx)) else color <- colorRampPalette(c("royalblue", "orange"))(length(tmpx))
    
  res <- persp(x = tmpx, y = tmpy, z = tmpz, theta = -57, phi = 30,  expand = 0.9,
               xlab = "\n \n Z1", ylab = "\n \n Z2", zlab = paste0("\n \n E(Y",item," | Z)"),#r = 5, 
               ticktype = "detailed", nticks = 5, border = T, cex.axis = 0.8, las = 3,
               ltheta = -45, lphi = 45, col = color[cutz], shade = 0.25)
  mtext(side = 3, line = 2, adj = 0, cex = 1.2, text = paste0("Fitted E(Item ",item," | Z)"))
  
  if(control.plot$plot.org){
   control.plot$sim.mod$ghQ <- mvghQ(mod$ghQ$n, formula = mod$formula)
   control.plot$sim.mod$fam <- fam
   f.morg <- grph.fun1(item, control.plot$sim.mod, qtl = control.plot$qtl)
   tmpB. <- cbind(round(mod$ghQ$points,5), f.morg$mean); tmpz. <- tmpz
   for(i in 1:length(tmpx)){
    for(j in 1:length(tmpy)){
        tmpz.[i,j] <- tmpB.[intersect(which(tmpB.[,1] == tmpx[i]), which(tmpB.[,2] == tmpy[j])),3]}}
   color. <- colorRampPalette(c("forestgreen","orange"))(length(tmpx))
   tmpzz. <- tmpz.[-1, -1] + tmpz.[-1, -ncol(tmpz.)] + tmpz.[-nrow(tmpz.), -1] + tmpz.[-nrow(tmpz.), -ncol(tmpz.)]
   cutz. <- cut(c(tmpzz.), 10)
    
   res <- persp(x = tmpx, y = tmpy, z = tmpz., theta = -57, phi = 30,  expand = 0.9,
                xlab = "\n \n Z1", ylab = "\n \n Z2", zlab = paste0("\n \n E(Y",item," | Z)"), # r = 5, 
                ticktype = "detailed", nticks = 5, border = T, cex.axis = 0.8, las = 3,
                ltheta = -45, lphi = 45, col = color.[cutz.], shade = 0.25)
   mtext(side = 3, line = 2, adj = 0, cex = 1.2, text = paste0("Simulated E(Item ",item," | Z)"))
   } } 
 }
 if(dimp > 2){
    
 stop("Code for graphics (q > 3) under development.")
 # cmb <- combn(dimp,2)
 # # if(control.plot$plot.org == T && sep.plots == F){X11(); par(mfrow=c(((ncol(cmb)*2+1)%/%2),2))} else{X11(); par(mfrow=c((ncol(cmb)+1)%/%2,2))}
 # for(k in 1:ncol(cmb)){
 # tmpx <- tmpy <- sort(unique(round(c(mod$ghQ$points),5))); tmpz <- matrix(NA, nrow = length(tmpx), ncol = length(tmpy))
 # f.mod <- fFun(item, fam[[item]], Z = mod$gr$out, b = mod$b, forms = mod$formula, lvp = c(cmb[,k]))
 # tmpB <- cbind(round(mod$gr$points[,paste0("Z",c(cmb[,k]))],5), fFun.mod$EY2M)
 # for(i in 1:length(tmpx)){
 #  for(j in 1:length(tmpy)){
 #      tmpz[i,j] <- unique(tmpB[intersect(which(tmpB[,1] == tmpx[i]), which(tmpB[,2] == tmpy[j])),3])}}
 # color <- colorRampPalette(c("lightblue", "purple"))(length(tmpx))
 # tmpzz <- tmpz[-1, -1] + tmpz[-1, -ncol(tmpz)] + tmpz[-nrow(tmpz), -1] + tmpz[-nrow(tmpz), -ncol(tmpz)]
 # cutz <- cut(c(tmpzz), 10)
 # 
 # res <- persp(x = tmpx, y = tmpy, z = tmpz, theta = -45, phi = 30,  expand = 0.8,
 #              xlab = paste0("\n Z",cmb[1,k]), ylab = paste0("\n Z",cmb[2,k]), zlab = paste0("\n E(Y",item," | Z)"),
 #              r = 5, ticktype = "detailed", nticks = 7, 
 #              ltheta = -45, lphi = 45, col = color[cutz], shade = 0.25,
 #              main = paste0("Fitted values for E(item ",item," | Z)"))
 # 
 # if(plot.org == T){
 # if(missing(morg)) stop("Argument `morg' missing.")
 # tmp <- mvgH(mod$gr$n, mu = morg$Z.mu, sigma = morg$Z.Sigma, formula = morg$formula)
 # fFun.morg <- fFun(item, fam[[item]], Z = tmp$out, b = morg$b, forms = morg$formula, lvp = c(cmb[,k]))
 # tmpB. <- cbind(round(mod$gr$points[,paste0("Z",c(cmb[,k]))],5), fFun.morg$EY2M); tmpz. <- tmpz
 # for(i in 1:length(tmpx)){
 #  for(j in 1:length(tmpy)){
 #      tmpz.[i,j] <- unique(tmpB.[intersect(which(tmpB.[,1] == tmpx[i]), which(tmpB.[,2] == tmpy[j])),3])}}
 # color. <- colorRampPalette(c("gold","royalblue1"))(length(tmpx))
 # tmpzz. <- tmpz.[-1, -1] + tmpz.[-1, -ncol(tmpz.)] + tmpz.[-nrow(tmpz.), -1] + tmpz.[-nrow(tmpz.), -ncol(tmpz.)]
 # cutz. <- cut(c(tmpzz.), 10)
 # 
 # res <- persp(x = tmpx, y = tmpy, z = tmpz., theta = -45, phi = 30,  expand = 0.8,
 #              xlab = paste0("\n Z",cmb[1,k]), ylab = paste0("\n Z",cmb[2,k]), zlab = paste0("\n E(Y",item," | Z)"),
 #              r = 5, ticktype = "detailed", nticks = 7, 
 #              ltheta = -45, lphi = 45, col = color.[cutz.], shade = 0.25,
 #              main = paste0("Original E(item ",item," | Z)")) } } }
 
#   if(plot.addpoints == T) points(trans3d(morg$Z$mu$Z1, morg$Z$mu$Z2, mod$Y[[paste0("Y",item)]], res), pch = 16, cex = 0.5, col = alpha("gray", 0.7))
}
} # else 3D
} # function

plot.score <- function(mod = ex1){
num.score <- grad(func = loglikFun, x = unlist(mod$b), method = "Richardson",
                  method.args=list(r = 8), Y = mod$Y, ghQ = mod$gr, fam = mod$fam, beta = mod$b)
tmp.min <- min(num.score,unlist(mod$Score)) - 1e-5
tmp.max <- max(num.score,unlist(mod$Score)) + 1e-5
t <- seq_along(num.score)
plot(x = t, y = unlist(mod$Score), pch = 2, col = "gray50", lwd = 2, main = "Analytical vs. Numerical",
        ylim = c(tmp.min,tmp.max), xlab = "Parameters (#)", ylab = "Score (@ ML estimates parameter vector)"); abline(h = 0, col = "red", lty = 2)
points(num.score, col = "red", lwd = 2, pch = 1)
segments(x0 = t, y0 = unlist(mod$Score), y1 = num.score, col = "blue", lty = 3, lwd = 1)
legend("bottomright",legend=c("Analytical score", "Numerical score"), col=c("gray50", "red"), pch = c(2,1), cex=1, inset = 0.02, bty = "n")
}



