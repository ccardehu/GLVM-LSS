
plotGLVM <- function(item = 1, mod = ex1, morg = simR,
                     plot.mean = T, plot.sd = F, plot.org = F,
                     quant = c(0.05,0.2,0.4,0.6,0.8,0.95),
                     plot.3D = F, sep.plots = T, plot.dist = T,
                     plot.addpoints = F
                     ){

if("scales" %in% rownames(installed.packages()) == F) install.packages("scales")
# item = 1; mod = ex1; morg = simR; plot.mean = T; plot.org = F
# plot.sd = T; plot.addpoints = F; sep.plots = F, plot.3D = T
# quant = c(0.05,0.2,0.4,0.6,0.8,0.95)

if(missing(sep.plots)) sep.plots <- F
if(missing(plot.dist)) plot.dist <- F
plot.quant <- ifelse(missing(quant),F,T)
fam <- mod$fam
fFun.mod <- fFun(item, fam[[item]], Z = mod$gr$out, b = mod$b, qnt = quant,forms = mod$formula)
dimp <- ncol(mod$gr$points); Zorg <- morg$Z$mu
 
if(plot.3D == F){
 if(sep.plots == F) par(mfrow=c((dimp+1)%/%2,ifelse(dimp == 1, 1, 2)))
 for(i in 1:dimp){
   plot(Zorg[[paste0("Z",i)]], mod$Y[[paste0("Y",item)]], pch = 16, col = scales::alpha("gray", 0.7),
        xlab = bquote(Z[.(i)]), ylab = bquote(Y[.(item)]))
   if(plot.mean == T) lines(mod$gr$points[,paste0("Z",i)], fFun.mod$eM[,paste0("Z",i)], lwd = 2, col = 2)
   plot.sd <- ifelse(plot.quant, F, plot.sd) # plot.quant overrides plot.sd
   if(plot.sd == T){
    lines(mod$gr$points[,paste0("Z",i)], fFun.mod$eM[,paste0("Z",i)] + fFun.mod$sM[,paste0("Z",i)], lwd = 2, col = 2, lty = 3, type = "b", pch = ".")
    lines(mod$gr$points[,paste0("Z",i)], fFun.mod$eM[,paste0("Z",i)] - fFun.mod$sM[,paste0("Z",i)], lwd = 2, col = 2, lty = 3, type = "b", pch = ".")
    #lines(mod$gr$points[,paste0("Z",i)], fFun.mod$eM[,paste0("Z",i)] + 1.96*fFun.mod$sM[,paste0("Z",i)], lwd = 2, col = 2, lty = 3)
    #lines(mod$gr$points[,paste0("Z",i)], fFun.mod$eM[,paste0("Z",i)] - 1.96*fFun.mod$sM[,paste0("Z",i)], lwd = 2, col = 2, lty = 3)
    mtext(side = 3, line = 1, at = -3.7, adj = 0, cex = 0.8, text = "E(Y|z) +1/-1 SD")
   }
   if(plot.quant == T){
    if(missing(quant)){warning("Argument `quant' is missing. Setting it to 0.05 and 0.95."); quant <- c(0.05,0.95)}
    for(j in seq_along(quant)){ lines(mod$gr$points[,paste0("Z",i)], fFun.mod$quantM[[paste0("Z",i)]][,j], lwd = 2, col = 2, lty = 3, type = "b", pch = ".") }
    mtext(side = 3, line = 1, at = -3.7, adj = 0, cex = 0.8, text = paste0("(Q: ", paste0(quant, collapse = ", "),")"))
   }
  
  if(plot.org == T){
   if(missing(morg)) stop("Argument `morg' missing.")
   tmp <- mvgH(mod$gr$n, mu = morg$Z.mu, sigma = morg$Z.Sigma, formula = morg$formula)
   fFun.org <- fFun(item, fam[[item]], Z = tmp$out, b = morg$borg, qnt = quant, forms = morg$formula)
   lines(mod$gr$points[,paste0("Z",i)], fFun.org$eM[,paste0("Z",i)], lwd = 2, col = "black")
   if(plot.sd == T){
    lines(mod$gr$points[,paste0("Z",i)], fFun.org$eM[,paste0("Z",i)] + fFun.org$sM[,paste0("Z",i)], lwd = 2, col = 1, lty = 3, type = "b", pch = ".")
    lines(mod$gr$points[,paste0("Z",i)], fFun.org$eM[,paste0("Z",i)] - fFun.org$sM[,paste0("Z",i)], lwd = 2, col = 1, lty = 3, type = "b", pch = ".")
   }
   if(plot.quant == T){
    if(missing(quant)){warning("Argument `quant' is missing. Setting it to 0.05 and 0.95."); quant <- c(0.05,0.95)}
    for(j in seq_along(quant)){ lines(mod$gr$points[,paste0("Z",i)], fFun.org$quantM[[paste0("Z",i)]][,j], lwd = 2, col = 1, lty = 3, type = "b", pch = ".") }
   }
   legend("topleft",legend=c("Fitted model", "Original model"), col=c("red", "black"), cex=0.8, lty = 1, inset = 0.02, bty = "n")
   mtext(side = 3, line = 2, at = -3.7, adj = 0, cex = 1,
         paste0("P(item ",item," | Z",i,"), with all other LV set to 0"))
  }
   else legend("topleft",legend=c("Fitted model"), col=c("red"), cex=0.8, lty = 1, inset = 0.02, lwd = 2, bty = "n")
   mtext(side = 3, line = 2, at = -3.7, adj = 0, cex = 1,
         paste0("P(item ",item," | Z",i,"), with all other LV set to 0"))
 }
} else { #else 3D
 if(dimp < 2){
   par(mfrow = c(1,1))
   tmpx1 <- sort(unique(round(c(mod$gr$points),5)));
   sel1 <- (tmpx1 > min(Zorg[[paste0("Z",dimp)]]) & tmpx1 < max(Zorg[[paste0("Z",dimp)]]))
   tmpx <- tmpx1[sel1]
   tmpy <- sort(unique(mod$Y[[paste0("Y",item)]]))
   vX <- seq(min(tmpx)-1,max(tmpx)+1,length = 2); vY <- seq(min(tmpy)-1,max(tmpy)+1,length = 2)
   
   res  <- persp(vX,vY,matrix(0,2,2),zlim=c(0,1),theta= -55, phi = 30, expand = 0.8, ticktype ="detailed",
                 box = T, xlab = "\n Z1", ylab = paste0("Y",item),zlab = paste0("\n P(Y",item," | Z)"),
                 nticks = 4, main = paste0("Fitted distribution for P(item ",item," | Z)"))
   
   if(plot.addpoints == T){ C <- trans3d(Zorg[[paste0("Z",dimp)]], mod$Y[[paste0("Y",item)]],rep(0,length(mod$Y[[paste0("Y",item)]])),res)
     points(C, pch = 20, col= scales::alpha("gray", 0.6), cex = 0.8) }
   fFun.mod <- fFun(i = item, fam = fam[[item]], Z = mod$gr$out, b = mod$b, forms = mod$formula, qnt = quant)
   plot.sd <- ifelse(plot.quant, F, plot.sd) # plot.quant overrides plot.sd
   if(plot.sd == T){
     C <- trans3d(tmpx,c(fFun.mod$mean + fFun.mod$sM)[sel1],rep(0,length(tmpx)),res); lines(C,lwd=2, lty = 3, col = 2)
     C <- trans3d(tmpx,c(fFun.mod$mean - fFun.mod$sM)[sel1],rep(0,length(tmpx)),res); lines(C,lwd=2, lty = 3, col = 2)
     C <- trans3d(c(tmpx,rev(tmpx)),c(c(fFun.mod$mean + fFun.mod$sM)[sel1],rev(c(fFun.mod$mean - fFun.mod$sM)[sel1])),rep(0,2*length(tmpx)),res)
     polygon(C,border=NA,col=scales::alpha("yellow", 0.5))
   }
   if(plot.mean == T){ C <- trans3d(tmpx,fFun.mod$mean[sel1],rep(0,length(tmpx)),res); lines(C,lwd=2, col = 2) }
   if(plot.dist == T){
   tmpd <- gFun(mod$gr$out,mod$b,fam[item],item)
   cuts <- 10
   vXa <- tmpx[round(seq(1,length(tmpx),length.out = cuts))]
     for(w in 1:cuts){
     if(length(tmpy) < 10) stp <- 10 else stp <- ifelse(round(length(tmpy)/10) < 10, 10, round(length(tmpy)/10))
     x <- rep(vXa[w],stp)
     # y <- seq(min(min(tmpy)-1, fFun.mod$quantM[[paste0("Z",dimp)]][,paste0("q",100*min(quant))]),
              # max(tmpy)+1, length = stp)
     liminf <- ifelse(fFun.mod$quant[which(vXa[w]==tmpx1),paste0("q",100*min(quant))]-2*fFun.mod$sd[which(vXa[w]==tmpx1)] < max(min(vY), p2Fun(fam[[item]])[1]),
                      max(min(vY),p2Fun(fam[[item]])[1]),fFun.mod$quant[which(vXa[w]==tmpx1),paste0("q",100*min(quant))]-2*fFun.mod$sd[which(vXa[w]==tmpx1)])
     limsup <- ifelse(fFun.mod$quant[which(vXa[w]==tmpx1),paste0("q",100*max(quant))]+2*fFun.mod$sd[which(vXa[w]==tmpx1)] > min(max(vY),p2Fun(fam[[item]])[2]),
                      min(max(vY),p2Fun(fam[[item]])[2]), fFun.mod$quant[which(vXa[w]==tmpx1),paste0("q",100*max(quant))]+2*fFun.mod$sd[which(vXa[w]==tmpx1)])
     z0 <- rep(0,stp)
     if(any(fam[[item]] == c("ZIpoisson","poisson","binom"))) {
     y <- round(seq(liminf,limsup, length = stp))
     z <- ifelse(f2Fun(y,fam[item],tmpd, which(vXa[w]==tmpx1)) > 1, 1, f2Fun(y,fam[item],tmpd, which(vXa[w]==tmpx1)))
     C=trans3d(c(x,x),c(y,rev(y)),c(z,z0),res)
     for(m in 1:length(z0)){ lines(trans3d(x[m],y[m],c(z0[m],z[m]),res), col = scales::alpha("blue",0.4), lty = 3) }
     C=trans3d(x,y,z0,res); lines(C,lty=3, col = scales::alpha("lightblue",0.6))
     C=trans3d(x,y,z,res)
     points(C, pch = 16, col= scales::alpha("blue",0.4), cex = 0.5)
     #lines(C,col="blue", lty = 3)
     } else {
     y <- seq(liminf,limsup, length = stp)
     z <- ifelse(f2Fun(y,fam[item],tmpd, which(vXa[w]==tmpx1)) > 1, 1, f2Fun(y,fam[item],tmpd, which(vXa[w]==tmpx1)))
     C=trans3d(c(x,x),c(y,rev(y)),c(z,z0),res)
     polygon(C,border=NA,col=scales::alpha("lightblue",0.3))
     C=trans3d(x,y,z0,res)
     lines(C,lty=3, col = scales::alpha("lightblue",0.6))
     C=trans3d(x,y,z,res)
     lines(C,col=scales::alpha("blue",0.4)) }
     }; legend("bottomleft",legend=c("Fitted distribution"), col= scales::alpha("blue",0.5), lty =1, inset = 0.05, bty = "n")
   }
   
   if(plot.quant == T){
     if(missing(quant)){warning("Argument `quant' is missing. Setting it to 0.05 and 0.95."); quant <- c(0.05,0.95)}
     for(j in quant){
       C <- trans3d(tmpx,fFun.mod$quant[sel1,paste0("q",100*j)],rep(0,length(tmpx)),res)
       lines(C,lwd=2, lty = 3, col = 2)
       #C <- trans3d(c(tmpx,rev(tmpx)),c(c(fFun.mod$mean + fFun.mod$sM)[sel1],rev(c(fFun.mod$mean - fFun.mod$sM)[sel1])),rep(0,2*length(tmpx)),res)
       #polygon(Ck,border=NA,col=scales::alpha("yellow", 0.5))
     }
   }
}
 if(dimp == 2){
 if(plot.org == T) par(mfrow = c(1,2)) else par(mfrow = c(1,1))
 for(k in 1:dimp){
  tmpx <- tmpy <- sort(unique(round(c(mod$gr$points),5))); tmpz <- matrix(NA, nrow = length(tmpx), ncol = length(tmpy))
  fFun.mod <- fFun(i = item, fam = fam[[item]], Z = mod$gr$out, b = mod$b, forms = mod$formula, lvp = c(1:dimp))
  tmpB <- cbind(round(mod$gr$points,5), fFun.mod$mean)
  for(i in 1:length(tmpx)){
   for(j in 1:length(tmpy)){
       tmpz[i,j] <- tmpB[intersect(which(tmpB[,1] == tmpx[i]), which(tmpB[,2] == tmpy[j])),3]}}
  color <- colorRampPalette(c("lightblue", "purple"))(length(tmpx))
  tmpzz <- tmpz[-1, -1] + tmpz[-1, -ncol(tmpz)] + tmpz[-nrow(tmpz), -1] + tmpz[-nrow(tmpz), -ncol(tmpz)]
  cutz <- cut(c(tmpzz), 10)
    
  res <- persp(x = tmpx, y = tmpy, z = tmpz, theta = -45, phi = 30,  expand = 0.8,
               xlab = "\n Z1", ylab = "\n Z2", zlab = paste0("\n E(Y",item," | Z)"),
               r = 5, ticktype = "detailed", nticks = 5, 
               ltheta = -45, lphi = 45, col = color[cutz], shade = 0.25,
               main = paste0("Fitted values for E(item ",item," | Z)"))
  
  # if(plot.addpoints == T){ C <- trans3d(Zorg[[paste0("Z",dimp)]], mod$Y[[paste0("Y",item)]],rep(0,length(mod$Y[[paste0("Y",item)]])),res)
  # points(C, pch = 20, col= scales::alpha("gray", 0.7)) }
  
  if(plot.org == T){
   if(missing(morg)) stop("Argument `morg' missing.")
   tmp <- mvgH(mod$gr$n, mu = morg$Z.mu, sigma = morg$Z.Sigma, formula = morg$formula)
   fFun.morg <- fFun(i = item, fam = fam[[item]], Z = tmp$out, b = morg$b, forms = morg$formula, lvp = c(1:dimp))
   tmpB. <- cbind(round(mod$gr$points,5), fFun.morg$mean); tmpz. <- tmpz
   for(i in 1:length(tmpx)){
    for(j in 1:length(tmpy)){
        tmpz.[i,j] <- tmpB.[intersect(which(tmpB.[,1] == tmpx[i]), which(tmpB.[,2] == tmpy[j])),3]}}
   color. <- colorRampPalette(c("gold","forestgreen"))(length(tmpx))
   tmpzz. <- tmpz.[-1, -1] + tmpz.[-1, -ncol(tmpz.)] + tmpz.[-nrow(tmpz.), -1] + tmpz.[-nrow(tmpz.), -ncol(tmpz.)]
   cutz. <- cut(c(tmpzz.), 10)
    
   res <- persp(x = tmpx, y = tmpy, z = tmpz., theta = -45, phi = 30,  expand = 0.8,
                xlab = "\n Z1", ylab = "\n Z2", zlab = paste0("\n E(Y",item," | Z)"),
                r = 5, ticktype = "detailed", nticks = 5, 
                ltheta = -45, lphi = 45, col = color.[cutz.], shade = 0.25,
                main = paste0("Original E(item ",item," | Z)"))
   } } 
 }
 if(dimp > 2){
 cmb <- combn(dimp,2)
 if(plot.org == T && sep.plots == F){X11(); par(mfrow=c(((ncol(cmb)*2+1)%/%2),2))} else{X11(); par(mfrow=c((ncol(cmb)+1)%/%2,2))}
 for(k in 1:ncol(cmb)){
 tmpx <- tmpy <- sort(unique(round(c(mod$gr$points),5))); tmpz <- matrix(NA, nrow = length(tmpx), ncol = length(tmpy))
 fFun.mod <- fFun(item, fam[[item]], Z = mod$gr$out, b = mod$b, forms = mod$formula, lvp = c(cmb[,k]))
 tmpB <- cbind(round(mod$gr$points[,paste0("Z",c(cmb[,k]))],5), fFun.mod$EY2M)
 for(i in 1:length(tmpx)){
  for(j in 1:length(tmpy)){
      tmpz[i,j] <- unique(tmpB[intersect(which(tmpB[,1] == tmpx[i]), which(tmpB[,2] == tmpy[j])),3])}}
 color <- colorRampPalette(c("lightblue", "purple"))(length(tmpx))
 tmpzz <- tmpz[-1, -1] + tmpz[-1, -ncol(tmpz)] + tmpz[-nrow(tmpz), -1] + tmpz[-nrow(tmpz), -ncol(tmpz)]
 cutz <- cut(c(tmpzz), 10)
 
 res <- persp(x = tmpx, y = tmpy, z = tmpz, theta = -45, phi = 30,  expand = 0.8,
              xlab = paste0("\n Z",cmb[1,k]), ylab = paste0("\n Z",cmb[2,k]), zlab = paste0("\n E(Y",item," | Z)"),
              r = 5, ticktype = "detailed", nticks = 7, 
              ltheta = -45, lphi = 45, col = color[cutz], shade = 0.25,
              main = paste0("Fitted values for E(item ",item," | Z)"))
 
 if(plot.org == T){
 if(missing(morg)) stop("Argument `morg' missing.")
 tmp <- mvgH(mod$gr$n, mu = morg$Z.mu, sigma = morg$Z.Sigma, formula = morg$formula)
 fFun.morg <- fFun(item, fam[[item]], Z = tmp$out, b = morg$b, forms = morg$formula, lvp = c(cmb[,k]))
 tmpB. <- cbind(round(mod$gr$points[,paste0("Z",c(cmb[,k]))],5), fFun.morg$EY2M); tmpz. <- tmpz
 for(i in 1:length(tmpx)){
  for(j in 1:length(tmpy)){
      tmpz.[i,j] <- unique(tmpB.[intersect(which(tmpB.[,1] == tmpx[i]), which(tmpB.[,2] == tmpy[j])),3])}}
 color. <- colorRampPalette(c("gold","royalblue1"))(length(tmpx))
 tmpzz. <- tmpz.[-1, -1] + tmpz.[-1, -ncol(tmpz.)] + tmpz.[-nrow(tmpz.), -1] + tmpz.[-nrow(tmpz.), -ncol(tmpz.)]
 cutz. <- cut(c(tmpzz.), 10)
 
 res <- persp(x = tmpx, y = tmpy, z = tmpz., theta = -45, phi = 30,  expand = 0.8,
              xlab = paste0("\n Z",cmb[1,k]), ylab = paste0("\n Z",cmb[2,k]), zlab = paste0("\n E(Y",item," | Z)"),
              r = 5, ticktype = "detailed", nticks = 7, 
              ltheta = -45, lphi = 45, col = color.[cutz.], shade = 0.25,
              main = paste0("Original E(item ",item," | Z)")) } } }
 
#   if(plot.addpoints == T) points(trans3d(morg$Z$mu$Z1, morg$Z$mu$Z2, mod$Y[[paste0("Y",item)]], res), pch = 16, cex = 0.5, col = alpha("gray", 0.7))
} # else 3D
} # function

plot.score <- function(mod = ex1){
num.score <- grad(func = loglik, x = unlist(mod$b), method = "Richardson",
                  method.args=list(r = 8), ghQ = mod$gr, loadmt = mod$loadmt, beta = mod$b)
tmp.min <- min(num.score,unlist(mod$Score)) - 1e-5
tmp.max <- max(num.score,unlist(mod$Score)) + 1e-5
t <- seq_along(num.score)
plot(x = t, y = unlist(mod$Score), pch = 2, col = "gray50", lwd = 2, main = "Analytical vs. Numerical",
        ylim = c(tmp.min,tmp.max), xlab = "Parameters (#)", ylab = "Score (@ ML estimates parameter vector)"); abline(h = 0, col = "red", lty = 2)
points(num.score, col = "red", lwd = 2, pch = 1)
segments(x0 = t, y0 = unlist(mod$Score), y1 = num.score, col = "blue", lty = 3, lwd = 1)
legend("bottomright",legend=c("Analytical score", "Numerical score"), col=c("gray50", "red"), pch = c(2,1), cex=1, inset = 0.02, bty = "n")
}



