
########################################################################
########################################################################
##### Update this file to include plots of zero coefficients, etc. #####
########################################################################
########################################################################


modt <- ex2

xlabn <- NULL
for(i in 1:ncol(Y)){
  xlabn <- append(xlabn,paste0("mu[",i,",",1:ncol(modt$b$mu),"]"))
  if("sigma" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("sigma[",i,",",1:ncol(modt$b$sigma),"]"))
  if("tau" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("tau[",i,",",1:ncol(modt$b$tau),"]"))
  if("nu" %in% pFun(modt$fam[i])) xlabn <- append(xlabn, paste0("nu[",i,",",1:ncol(modt$b$nu),"]"))
}

testH1 <- decoefmod(modt$b)
#if(sum(which(testH1 == 0)) == 0){ id0 <- NULL } else 
id0 <- which(testH1 == 0)
testH2 <- decoefmod(borg)
testH3 <- unname(rep(1,length(unlist(borg))))*decoefmod(l1)

num.score1 <- grad(func = loglikFun, x = testH1, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
num.score2 <- grad(func = loglikFun, x = testH2, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
num.score3 <- grad(func = loglikFun, x = testH3, method = "Richardson", method.args=list(r = 6, v = 2),
                   Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)

num.hess1 <- hessian(func = loglikFun, x = testH1, method = "Richardson", method.args=list(r = 6, v = 2),
                     Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b);
num.hess2 <- hessian(func = loglikFun, x = testH2, method = "Richardson", method.args=list(r = 6, v = 2),
                     Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b);
num.hess3 <- hessian(func = loglikFun, x = testH3, method = "Richardson", method.args=list(r = 6, v = 2),
                     Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b);

plot.jac = F
if(plot.jac == T){
  num.jac1 <- jacobian(func = ScoreFun, x = testH1, method = "Richardson", method.args=list(r = 6, v = 2),
                       Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
  num.jac2 <- jacobian(func = ScoreFun, x = testH2, method = "Richardson", method.args=list(r = 6, v = 2),
                       Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
  num.jac3 <- jacobian(func = ScoreFun, x = testH3, method = "Richardson", method.args=list(r = 6, v = 2),
                       Y = modt$Y, ghQ = modt$gr, fam = modt$fam, beta = modt$b)
}

ana1 <- ScHesFun(testH1, modt$Y, modt$b, modt$gr, modt$fam)
ana2 <- ScHesFun(testH2, modt$Y, modt$b, modt$gr, modt$fam)
ana3 <- ScHesFun(testH3, modt$Y, modt$b, modt$gr, modt$fam)

# profvis::profvis({t1 <- ScHesFun(testH1, modt$Y, modt$b, modt$gr, modt$fam)})

ana.score1 <- rep(0,length(testH1))
ana.score2 <- rep(0,length(testH2))
ana.score3 <- rep(0,length(testH3))
# if(id0 == integer(0))
ana.score1[-id0] <- ana1$score
ana.score2[-id0] <- ana2$score
ana.score3[-id0] <- ana3$score

ana.hess1 <- matrix(0,length(ana.score1),length(ana.score1))
ana.hess2 <- matrix(0,length(ana.score2),length(ana.score2))
ana.hess3 <- matrix(0,length(ana.score3),length(ana.score3))
ana.hess1[-id0,-id0] <- ana1$hess
ana.hess2[-id0,-id0] <- ana2$hess
ana.hess3[-id0,-id0] <- ana3$hess

num.score1[id0] <- 0
num.score2[id0] <- 0
num.score3[id0] <- 0
num.hess1[c(id0),] <- 0
num.hess1[,c(id0)] <- 0
num.hess2[c(id0),] <- 0
num.hess2[,c(id0)] <- 0
num.hess3[c(id0),] <- 0
num.hess3[,c(id0)] <- 0

par("mar" = c(6, 4, 4, 2) + 0.1)
plot(ana.score1, main = "Score (first deriv. of log-likelihood) @ MLE", xlab = "", ylab = "Value", pch = 16,
     col = "gray50", xaxt = "n")
mtext(side = 1, text = "Parameter index", line = 4)
for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
abline(h = 0, col = "forestgreen", lty = 2)
abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
               by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")
points(num.score1[1:length(ana.score1)], col = "red", pch = 10)
segments(x0 = seq_along(ana.score1), y0 = unlist(ana.score1), y1 = num.score1, col = "blue", lty = 3, lwd = 1)
legend("topleft",legend=c("Analytical", "Numerical"), col=c("gray50", "red"),
       pch = c(16,10), cex=1, inset = 0.02, box.col = "white")

# plot(ana.score1 - num.score1,main = "Differences in Score values: Analytical vs. Numerical (@ MLE)",xlab ="", ylab = "Value", pch = 16, col = "gray50", xaxt = "n")
# mtext(side = 1, text = "Parameter index", line = 4)
# for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
# abline(h = 0, col = "forestgreen", lty = 2)
# abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
#                by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")

plot(ana.score2, main = "Score (first deriv. of log-likelihood) @ Original Parameters", xlab = "", ylab = "value", pch = 16,
     col = "gray50", xaxt = "n")
mtext(side = 1, text = "Parameter index", line = 4)
for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
abline(h = 0, col = "forestgreen", lty = 2)
abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
               by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")
points(num.score2[1:length(ana.score2)], col = "red", pch = 10)
segments(x0 = seq_along(ana.score2), y0 = unlist(ana.score2), y1 = num.score2, col = "blue", lty = 3, lwd = 1)
legend("bottomleft",legend=c("Analytical", "Numerical"), col=c("gray50", "red"), pch = c(16,10), cex=1, inset = 0.02, box.col = "white")

# plot(ana.score2 - num.score2,main = "Differences in Score values: Analytical vs. Numerical (@ Original)",xlab ="", ylab = "Value", pch = 16, col = "gray50", xaxt = "n")
# mtext(side = 1, text = "Parameter index", line = 4)
# for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
# abline(h = 0, col = "forestgreen", lty = 2)
# abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
#                by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")

plot(ana.score3, main = "Score (first deriv. of log-likelihood) @ Vector of 1s", xlab = "", ylab = "value", pch = 16,
     col = "gray50", xaxt = "n")
mtext(side = 1, text = "Parameter index", line = 4)
for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
abline(h = 0, col = "forestgreen", lty = 2)
abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
               by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")
points(num.score3[1:length(ana.score3)], col = "red", pch = 10)
segments(x0 = seq_along(ana.score3), y0 = unlist(ana.score3), y1 = num.score3, col = "blue", lty = 3, lwd = 1)
legend("topleft",legend=c("Analytical", "Numerical"), col=c("gray50", "red"), pch = c(16,10), cex=1, inset = 0.02, box.col = "white")

# plot(ana.score3 - num.score3,main = "Differences in Score values: Analytical vs. Numerical (@ 1s)",xlab ="", ylab = "Value", pch = 16, col = "gray50", xaxt = "n")
# mtext(side = 1, text = "Parameter index", line = 4)
# for(i in seq_along(xlabn)){ axis(1, at=i, labels = xlabn[i] , las = 2, cex.axis = 0.7) }
# abline(h = 0, col = "forestgreen", lty = 2)
# abline(v = seq(from = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y))),length.out = ncol(Y)-1,
#                by = ((sum(length(modt$b$mu), length(modt$b$sigma)) %/% ncol(Y)))) + 0.5, lty = 3, col = "gray80")

int <- 1:sum(length(modt$b$mu), length(modt$b$sigma))^2
# int <- 1:((sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) * 2)
plot(c(ana.hess1)[int], main = "Hessian (2nd deriv. of log-likelihood) @ MLE",xlab ="Hessian entry Index", ylab = "Value",
     cex = 0.5, pch = 16, col = "gray50", ylim = c(min(ana.hess1,num.hess1),max(ana.hess1,num.hess1)))#, xaxt = "n")
abline(h = 0, col = "forestgreen", lty = 2);
abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")
points(c(num.hess1)[int], col = "red", pch = 16, cex = 0.5)
segments(x0 = seq_along(ana.hess1[int]), y0 = unlist(ana.hess1)[int], y1 = num.hess1[int], col = "blue", lty = 3, lwd = 1)
legend(0,-2500,legend=c("Analytical", "Numerical (Hessian)"), col=c("gray50", "red",0.6),
       pch = c(16,16), cex=0.8, inset = 0.02, box.col = "white")
if(plot.jac == T){ points(c(num.jac1)[int], col = scales::alpha("blue",0.6), pch = 16, cex = 0.5) 
  legend("bottomright",legend=c("Analytical", "Numerical (Hessian)", "Numerical (Jacobian)"),
         col=c("gray50",scales::alpha("red",0.6), scales::alpha("blue",0.6)), pch = c(16,16,16),
         cex=1, inset = 0.02, box.col = "white")}

# plot(c(num.hess1 - ana.hess1), main = "Differences in Hessian Analytical vs. Numerical (@ MLE)",xlab ="Hessian entry Index", ylab = "Value", pch = 16, col = "gray50")
# abline(h = 0, col = "forestgreen", lty = 2)
# abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
#               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")

plot(c(ana.hess2)[int], main = "Hessian (2nd deriv. of log-likelihood) @ Original Parameters",
     xlab = "Hessian entry Index", ylab = "value", cex = 0.5, pch = 16, col = "gray50")
abline(h = 0, col = "forestgreen", lty = 2);
abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")
points(c(num.hess2)[int], col = "red", pch = 16, cex = 0.5)
segments(x0 = seq_along(ana.hess2[int]), y0 = unlist(ana.hess2)[int], y1 = num.hess2[int], col = "blue", lty = 3, lwd = 1)
legend(0,-3000,legend=c("Analytical", "Numerical (Hessian)"), col=c("gray50", "red"),
       pch = c(16,16), cex=1, inset = 0.02, box.col = "white")
if(plot.jac == T){ points(c(num.jac2)[int], col = scales::alpha("blue",0.6), pch = 16, cex = 0.5) 
  legend("bottomright",legend=c("Analytical", "Numerical (Hessian)", "Numerical (Jacobian)"),
         col=c("gray50",scales::alpha("red",0.6), scales::alpha("blue",0.6)), pch = c(16,16,16),
         cex=1, inset = 0.02, box.col = "white")}

# plot(c(num.hess2 - ana.hess2),main = "Differences in Hessian Analytical vs. Numerical (@ Original)",xlab ="Hessian entry Index", ylab = "Value", pch = 16, col = "gray50")
# abline(h = 0, col = "forestgreen", lty = 2)
# abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
#               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")

plot(c(ana.hess3)[int], main = "Hessian (2nd deriv. of log-likelihood) @ Vector of 1s",
     xlab = "Hessian entry Index", ylab = "value", cex = 0.5, pch = 16, col = "gray50", ylim = c(min(ana.hess3,num.hess3),max(ana.hess3,num.hess3)))
abline(h = 0, col = "forestgreen", lty = 2);
abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")
points(c(num.hess3)[int], col = "red", pch = 16, cex = 0.5)
segments(x0 = seq_along(ana.hess3[int]), y0 = unlist(ana.hess3)[int], y1 = num.hess3[int], col = "blue", lty = 3, lwd = 1)
legend(2000,-850,legend=c("Analytical", "Numerical (Hessian)"), col=c("gray50", "red",0.6),
       pch = c(16,16), cex=0.8, inset = 0.02, box.col = "white")
if(plot.jac == T){ points(c(num.jac3)[int], col = scales::alpha("blue",0.6), pch = 16, cex = 0.5) 
  legend("bottomleft",legend=c("Analytical", "Numerical (Hessian)", "Numerical (Jacobian)"),
         col=c("gray50",scales::alpha("red",0.6), scales::alpha("blue",0.6)), pch = c(16,16,16),
         cex=1, inset = 0.02, box.col = "white")}

# plot(c(num.hess3 - ana.hess3),main = "Differences in Analytical Hessian vs. Numerical Hessian (@ 1s)",xlab ="Hessian entry Index", ylab = "Value", pch = 16, col = "gray50")
# abline(h = 0, col = "forestgreen", lty = 2)
# abline(v = seq(from = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y),length.out = ncol(Y)-1,
#               by = sum(length(modt$b$mu), length(modt$b$sigma))^2 %/% ncol(Y)) + 0.5, lty = 3, col = "gray80")