glvmss_plot <- function(item, mod, 
                       control.plot = list(sim.mod = simR,
                       plot.mean = T, plot.sd = F,
                       qtl = c(seq(0.1,0.9,by = 0.2)),
                       plot.3D = F, bw = F)){
  
  family = list(),
  mu.eq = ~ Z1, sg.eq = NULL,
  nu.eq = NULL, ta.eq = NULL,
  
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

