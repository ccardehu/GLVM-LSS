
plotGLVM <- function(item = 1, mod = ex1, morg = simR,
                     plot.mean = T, plot.sd = F, plot.org = F,
                     quant = c(0.05,0.2,0.4,0.6,0.8,0.95),
                     plot.3D = F, sep.plots = T, plot.dist = T#,
                     #plot.addpoints = F
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
dimp <- ncol(mod$gr$points); 
 
if(plot.3D == F){
 if(sep.plots == F) par(mfrow=c((dimp+1)%/%2,ifelse(dimp == 1, 1, 2)))
 Zorg <- morg$Z$mu
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
 if(dimp < 2) stop("Change argument `plot.3D': Only one latent variable in the model.")
 if(dimp == 2){
 if(plot.org == T) par(mfrow = c(1,2)) else par(mfrow = c(1,1))
 for(k in 1:dimp){
  tmpx <- tmpy <- sort(unique(round(c(mod$gr$points),5))); tmpz <- matrix(NA, nrow = length(tmpx), ncol = length(tmpy))
  fFun.mod <- fFun(i = item, fam = fam[[item]], Z = mod$gr$out, b = mod$b, forms = mod$formula, lvp = c(1:dimp))
  tmpB <- cbind(round(mod$gr$points,5), fFun.mod$EY2M)
  for(i in 1:length(tmpx)){
   for(j in 1:length(tmpy)){
       tmpz[i,j] <- tmpB[intersect(which(tmpB[,1] == tmpx[i]), which(tmpB[,2] == tmpy[j])),3]}}
  color <- colorRampPalette(c("lightblue", "purple"))(length(tmpx))
  tmpzz <- tmpz[-1, -1] + tmpz[-1, -ncol(tmpz)] + tmpz[-nrow(tmpz), -1] + tmpz[-nrow(tmpz), -ncol(tmpz)]
  cutz <- cut(c(tmpzz), 10)
    
  res <- persp(x = tmpx, y = tmpy, z = tmpz, theta = -45, phi = 30,  expand = 0.8,
               xlab = "\n Z1", ylab = "\n Z2", zlab = paste0("\n E(Y",item," | Z)"),
               r = 5, ticktype = "detailed", nticks = 7, 
               ltheta = -45, lphi = 45, col = color[cutz], shade = 0.25,
               main = paste0("Fitted values for E(item ",item," | Z)"))
  if(plot.org == T){
   if(missing(morg)) stop("Argument `morg' missing.")
   tmp <- mvgH(mod$gr$n, mu = morg$Z.mu, sigma = morg$Z.Sigma, formula = morg$formula)
   fFun.morg <- fFun(i = item, fam = fam[[item]], Z = tmp$out, b = morg$b, forms = morg$formula, lvp = c(1:dimp))
   tmpB. <- cbind(round(mod$gr$points,5), fFun.morg$EY2M); tmpz. <- tmpz
   for(i in 1:length(tmpx)){
    for(j in 1:length(tmpy)){
        tmpz.[i,j] <- tmpB.[intersect(which(tmpB.[,1] == tmpx[i]), which(tmpB.[,2] == tmpy[j])),3]}}
   color. <- colorRampPalette(c("gold","forestgreen"))(length(tmpx))
   tmpzz. <- tmpz.[-1, -1] + tmpz.[-1, -ncol(tmpz.)] + tmpz.[-nrow(tmpz.), -1] + tmpz.[-nrow(tmpz.), -ncol(tmpz.)]
   cutz. <- cut(c(tmpzz.), 10)
    
   res <- persp(x = tmpx, y = tmpy, z = tmpz., theta = -45, phi = 30,  expand = 0.8,
                xlab = "\n Z1", ylab = "\n Z2", zlab = paste0("\n E(Y",item," | Z)"),
                r = 5, ticktype = "detailed", nticks = 7, 
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






# attach(cars)
# n=2
# X= cars$speed
# Y=cars$dist
# df=data.frame(X,Y)
# vX=seq(min(X)-2,max(X)+2,length=n)
# vY=seq(min(Y)-15,max(Y)+15,length=n)
# mat=persp(vX,vY,matrix(0,n,n),zlim=c(0,.1),theta=-30,ticktype ="detailed", box = FALSE)
# reggig=glm(Y~X,data=df,family=gaussian(link="identity"))
# x=seq(min(X),max(X),length=501)
# C=trans3d(x,predict(reggig,newdata=data.frame(X=x),type="response"),rep(0,length(x)),mat)
# lines(C,lwd=2)
# sdgig=sqrt(summary(reggig)$dispersion)
# x=seq(min(X),max(X),length=501)
# y1=qnorm(.95,predict(reggig,newdata=data.frame(X=x),type="response"), sdgig)
# C=trans3d(x,y1,rep(0,length(x)),mat)
# lines(C,lty=2)
# y2=qnorm(.05,predict(reggig,newdata=data.frame(X=x),type="response"), sdgig)
# C=trans3d(x,y2,rep(0,length(x)),mat)
# lines(C,lty=2)
# C=trans3d(c(x,rev(x)),c(y1,rev(y2)),rep(0,2*length(x)),mat)
# polygon(C,border=NA,col="yellow")
# C=trans3d(X,Y,rep(0,length(X)),mat)
# points(C,pch=19,col="red")
# n=8
# vX=seq(min(X),max(X),length=n)
# mgig=predict(reggig,newdata=data.frame(X=vX))
# sdgig=sqrt(summary(reggig)$dispersion)
# for(j in n:1){
#    stp=251
#    x=rep(vX[j],stp)
#    y=seq(min(min(Y)-15,qnorm(.05,predict(reggig,newdata=data.frame(X=vX[j]),type="response"), sdgig)),max(Y)+15,length=stp)
#    z0=rep(0,stp)
#    z=dnorm(y, mgig[j], sdgig)
#    C=trans3d(c(x,x),c(y,rev(y)),c(z,z0),mat)
#    polygon(C,border=NA,col="light blue",density=40)
#    C=trans3d(x,y,z0,mat)
#    lines(C,lty=2)
#    C=trans3d(x,y,z,mat)
#    lines(C,col="blue")}

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



