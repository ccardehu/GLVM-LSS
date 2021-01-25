
library(scales)

plotGLVM <- function(item = 1, mod = ex1, morg = simR,
                     plot.mean = T, plot.sd = F, plot.org = T,
                     quant = c(0.05,0.2,0.4,0.6,0.8,0.95), # plot.quant = T, 
                     plot.addpoints = F){

# item = 1; mod = ex1; morg = simR; plot.org = T; plot.quant = T; plot.addpoints = F; plot.mean = T; plot.sd = F
# quant = c(0.05,0.2,0.4,0.6,0.8,0.95)

 fam <- mod$fam;# Y <- mod$Y; 
 plot.quant <- ifelse(missing(quant),F,T)
 fFun.mod <- fFun(item, fam[[item]], Z = mod$gr$out, b = mod$b, qnt = quant)
 
 if(ncol(mod$gr$points) == 1){
  plot(morg$Z$mu$Z1, mod$Y[[paste0("Y",item)]], pch = 16, col = alpha("gray", 0.7),
       xlab = expression("Latent Variable Z"[1]), ylab = substitute("Y"[item]))
 
 if(plot.mean == T) lines(mod$gr$points, fFun.mod$mean, lwd = 2, col = 2)
 plot.sd <- ifelse(plot.quant, F, plot.sd)
 
 if(plot.sd == T){
  lines(mod$gr$points, fFun.mod$mean + fFun.mod$sd, lwd = 2, col = 2, lty = 3)
  lines(mod$gr$points, fFun.mod$mean - fFun.mod$sd, lwd = 2, col = 2, lty = 3)
 }
  
 if(plot.quant == T){
  if(missing(quant)){warning("Argument `quant' is missing. Setting it to 0.05 and 0.95."); quant <- c(0.05,0.95)}
  for(i in seq_along(quant)){ lines(mod$gr$points, fFun.mod$quant[,i], lwd = 2, col = 2, lty = 3) }
 }

 if(plot.org == T){
  if(missing(morg)) stop("Argument `morg' missing.")
  tmp <- mvgH(length(mod$gr$points), mu = morg$Z.mu, sigma = morg$Z.Sigma, formula = morg$formula)
  fFun.org <- fFun(item, fam[[item]], Z = tmp$out, b = morg$borg, qnt = quant)
   
  lines(mod$gr$points, fFun.org$mean, lwd = 2, col = "black")
  if(plot.sd == T){
   lines(mod$gr$points, fFun.org$mean + fFun.org$sd, lwd = 2, col = 1, lty = 3)
   lines(mod$gr$points, fFun.org$mean - fFun.org$sd, lwd = 2, col = 1, lty = 3)  
  }
   
  if(plot.quant == T){
   if(missing(quant)){warning("Argument `quant' is missing. Setting it to 0.05 and 0.95."); quant <- c(0.05,0.95)}
    for(i in seq_along(quant)){ lines(mod$gr$points, fFun.org$quant[,i], lwd = 2, col = 1, lty = 3) }
  }
   
 legend("topleft",legend=c("Fitted model", "Original model"), col=c("red", "black"), cex=0.8, lty = 1, inset = 0.02, bty = "n")
 }
else legend("topleft",legend=c("Fitted model"), col=c("red"), cex=0.8, lty = 1, inset = 0.02, lwd = 2, bty = "n")
}
  
 else{
  #if(ncol(mod$gr$points) >= 1){
    
  if(plot.org == F){
    
   par(mfrow=c(1,1))
    
   tmpx <- tmpy <- -round(unique(mod$gr$points[,1]), 5); tmpz <- matrix(NA, nrow = length(tmpx), ncol = length(tmpy))
   tmpB <- cbind(round(mod$gr$points,5), fFun(item, fam[[item]], Z = mod$gr$out, b = mod$b)$EY)
   for(i in 1:length(tmpx)){
    for(j in 1:length(tmpy)){
     tmpz[i,j] <- tmpB[intersect(which(tmpB[,1] == tmpx[i]), which(tmpB[,2] == tmpy[j])),3]
    }
   }
   
   figcol <- colorRampPalette(c("lightblue", "purple"))
   color <- figcol(length(tmpx))
   tmpzz <- tmpz[-1, -1] + tmpz[-1, -ncol(tmpz)] + tmpz[-nrow(tmpz), -1] + tmpz[-nrow(tmpz), -ncol(tmpz)]
   cutz <- cut(c(tmpzz), 10)
   
   #rownames(tmpz) <- tmpx; colnames(tmpz) <- tmpy 
   
   res <- persp(x = tmpx, y = tmpy, z = tmpz, theta = -45, phi = 30,  expand = 0.8,
                xlab = "\n Z1", ylab = "\n Z2", zlab = paste0("\n E(Y",item," | Z)"),
                r = 5, ticktype = "detailed", nticks = 7, 
                ltheta = -45, lphi = 45, col = color[cutz], shade = 0.25, main = paste0("Expected (fitted) value of item Y",item))
   
   if(plot.addpoints == T) points(trans3d(morg$Z$mu$Z1, morg$Z$mu$Z2, mod$Y[[paste0("Y",item)]], res), pch = 16, cex = 0.5, col = alpha("gray", 0.7))
  }
   
  if(plot.org == T){
   
   tmp <- mvgH(length(unique(mod$gr$points[,1])), mu = morg$Z.mu, sigma = morg$Z.Sigma, formula = morg$formula)
   
   par(mfrow=c(1,2))
    
   tmpx <- tmpy <- -round(unique(mod$gr$points[,1]), 5)
   tmpz1 <- matrix(NA, nrow = length(tmpx), ncol = length(tmpy))
   tmpz2 <- matrix(NA, nrow = length(tmpx), ncol = length(tmpy))
   tmpB <- round(cbind(mod$gr$points, fFun(item, fam[[item]], Z = mod$gr$out, b = mod$b)$EY), 5)
   tmpA <- round(cbind(mod$gr$points, fFun(item, fam[[item]], Z = tmp$out, b = morg$borg)$EY), 5)
   for(i in 1:nrow(tmpB)){tmpz1[which(tmpx == tmpB[i,1]), which(tmpy == tmpB[i,2])] <- tmpB[i,3]}
   for(i in 1:nrow(tmpB)){tmpz2[which(tmpx == tmpA[i,1]), which(tmpy == tmpA[i,2])] <- tmpA[i,3]}
   
   figcol <- colorRampPalette(c("lightblue", "purple"))
   color <- figcol(length(tmpx))
   tmpzz1 <- tmpz1[-1, -1] + tmpz1[-1, -ncol(tmpz1)] + tmpz1[-nrow(tmpz1), -1] + tmpz1[-nrow(tmpz1), -ncol(tmpz1)]
   tmpzz2 <- tmpz2[-1, -1] + tmpz2[-1, -ncol(tmpz2)] + tmpz2[-nrow(tmpz2), -1] + tmpz2[-nrow(tmpz2), -ncol(tmpz2)]
   cutz1 <- cut(c(tmpzz1), 10)
   cutz2 <- cut(c(tmpzz2), 10)
   
   res1 <- persp(x = tmpx, y = tmpy, z = tmpz1, theta = -45, phi = 30,  expand = 0.8,
                 xlab = "\n Z1", ylab = "\n Z2", zlab = paste0("\n E(Y",item," | Z)"),
                 r = 5, ticktype = "detailed", nticks = 5, 
                 ltheta = -45, lphi = 45, col = color[cutz1], shade = 0.25, main = paste0("Expected (fitted) value of item Y",item))
   if(plot.addpoints == T) points(trans3d(morg$Z$mu$Z1, morg$Z$mu$Z2, mod$Y[[paste0("Y",item)]], res1), pch = 16, cex = 0.5, col = alpha("gray", 0.7))
   
   res2 <- persp(x = tmpx, y = tmpy, z = tmpz2, theta = -45, phi = 30,  expand = 0.8,
                 xlab = "\n Z1", ylab = "\n Z2", zlab = paste0("\n E(Y",item," | Z)"),
                 r = 5, ticktype = "detailed", nticks = 5, 
                 ltheta = -45, lphi = 45, col = color[cutz2], shade = 0.25, main = paste0("Expected (original) value of item Y",item))
   if(plot.addpoints == T) points(trans3d(morg$Z$mu$Z1, morg$Z$mu$Z2, mod$Y[[paste0("Y",item)]], res2), pch = 16, cex = 0.5, col = alpha("gray", 0.7)) 
    
  }
    
  
   #lines (trans3d(x = tmpx, y = max(tmpy), z = 0.5, pmat = res), col = alpha("gray", 1), lty = 3)
   #lines (trans3d(x = max(tmpx), y = tmpy, z = 0.5, pmat = res), col = alpha("gray", 1), lty = 3)
  
 }
  
  
  
}






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
