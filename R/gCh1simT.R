# source("gCh1simE0n500p10.R")
# source("gCh1simE0n500p30.R")

# source("gCh1simE2n200p30.R")
# source("gCh1simE2n500p30.R")
# source("gCh1simE2n1000p30.R")

source("gCh1simB.R")
form <- list("mu" = "~ Z1+Z2", "sigma" = "~ Z1+Z2");
load(file = "Ch1_Ex1_nsim500_n200_p10.RData"); p = 10; Ch1pos(FCOL,form,p)
load(file = "Ch1_Ex1_nsim500_n500_p10.RData"); p = 10; Ch1pos(FCOL,form,p)
load(file = "Ch1_Ex1_nsim500_n1000_p10.RData"); p = 10; Ch1pos(FCOL,form,p)
load(file = "Ch1_Ex1_nsim500_n5000_p10.RData"); p = 10; Ch1pos(FCOL,form,p)
load(file = "Ch1_Ex1_nsim500_n200_p30.RData"); p = 30; Ch1pos(FCOL,form,p)
load(file = "Ch1_Ex1_nsim500_n500_p30.RData"); p = 30; Ch1pos(FCOL,form,p)
load(file = "Ch1_Ex1_nsim500_n1000_p30.RData"); p = 30; Ch1pos(FCOL,form,p)
load(file = "Ch1_Ex1_nsim500_n5000_p30.RData"); p = 30; Ch1pos(FCOL,form,p)

form <- list("mu" = "~ Z1+Z2", "sigma" = "~ Z1+Z2")
load(file = "Ch1_Ex2_nsim500_n200_p10.RData"); p = 10; Ch1pos(FCOL,form,p)
load(file = "Ch1_Ex2_nsim500_n500_p10.RData"); p = 10; Ch1pos(FCOL,form,p)
load(file = "Ch1_Ex2_nsim500_n1000_p10.RData"); p = 10; Ch1pos(FCOL,form,p)
load(file = "Ch1_Ex2_nsim500_n5000_p10.RData"); p = 10; Ch1pos(FCOL,form,p)
# load(file = "Ch1_Ex2_nsim500_n200_p30.RData"); p = 30; round(Ch1pos(FCOL,form,p),4)
# load(file = "Ch1_Ex2_nsim500_n500_p30.RData"); p = 30; round(Ch1pos(FCOL,form,p),4)
# load(file = "Ch1_Ex2_nsim500_n1000_p30.RData"); p = 30; round(Ch1pos(FCOL,form,p),4)

form <- list("mu" = "~ Z1+Z2", "sigma" = "~ Z1+Z2")
load(file = "Ch1_Ex3_nsim500_n200_p10.RData"); p = 10; round(Ch1pos(FCOL,form,p),4)
load(file = "Ch1_Ex3_nsim500_n500_p10.RData"); p = 10; round(Ch1pos(FCOL,form,p),4)
load(file = "Ch1_Ex3_nsim500_n1000_p10.RData"); p = 10; round(Ch1pos(FCOL,form,p),4)
load(file = "Ch1_Ex3_nsim500_n200_p30.RData"); p = 30; round(Ch1pos(FCOL,form,p),4)
load(file = "Ch1_Ex3_nsim500_n500_p30.RData"); p = 30; round(Ch1pos(FCOL,form,p),4)
load(file = "Ch1_Ex3_nsim500_n1000_p30.RData"); p = 30; round(Ch1pos(FCOL,form,p),4)
 
# form <- list("mu" = "~ Z1+Z2")
# load(file = "Ch1_Ex0_nsim100_n500_p10.RData"); p = 10; round(Ch1pos(FCOL,form,p),4)
# load(file = "Ch1_Ex0_nsim100_n500_p30.RData"); p = 30; round(Ch1pos(FCOL,form,p),4)


plot.ts(abs(colMeans(Xorg - Xest))) #/colMeans(Xorg))
points(y = abs(colMeans(Xorg - Xest))[which(grepl("mu",colnames(Xorg),fixed = T))],
       x = which(grepl("mu",colnames(Xorg),fixed = T)), col = 2, pch = 19)
points(y = abs(colMeans(Xorg - Xest))[which(grepl("sigma",colnames(Xorg),fixed = T))],
       x = which(grepl("sigma",colnames(Xorg),fixed = T)), col = 3, pch = 19)
plot(colMeans(Xorg), colMeans(Xest)); abline(a = 0, b = 1, col = 2, lwd = 3)
boxplot(Xorg - Xest); abline(h = 0, col = 2, lwd = 3, lty = 3)


plot.ts(abs(colMeans(Xorg - Xest)))
points(y = abs(colMeans(Xorg - Xest))[which(grepl(".0",colnames(Xorg),fixed = T))],
       x = which(grepl(".0",colnames(Xorg),fixed = T)), col = 2, pch = 19)
points(y = abs(colMeans(Xorg - Xest))[which(grepl(".1",colnames(Xorg),fixed = T))],
       x = which(grepl(".1",colnames(Xorg),fixed = T)), col = 3, pch = 19)
points(y = abs(colMeans(Xorg - Xest))[which(grepl(".2",colnames(Xorg),fixed = T))],
       x = which(grepl(".2",colnames(Xorg),fixed = T)), col = 4, pch = 19)

plot.ts(colMeans((Xorg - Xest)^2))
points(y = colMeans((Xorg - Xest)^2)[which(grepl("mu",colnames(Xorg),fixed = T))],
       x = which(grepl("mu",colnames(Xorg),fixed = T)), col = 2, pch = 19)
points(y = colMeans((Xorg - Xest)^2)[which(grepl("sigma",colnames(Xorg),fixed = T))],
       x = which(grepl("sigma",colnames(Xorg),fixed = T)), col = 3, pch = 19)

plot.ts(colMeans((Xorg - Xest)^2))
points(y = colMeans((Xorg - Xest)^2)[which(grepl(".0",colnames(Xorg),fixed = T))],
       x = which(grepl(".0",colnames(Xorg),fixed = T)), col = 2, pch = 19)
points(y = colMeans((Xorg - Xest)^2)[which(grepl(".1",colnames(Xorg),fixed = T))],
       x = which(grepl(".1",colnames(Xorg),fixed = T)), col = 3, pch = 19)
points(y = colMeans((Xorg - Xest)^2)[which(grepl(".2",colnames(Xorg),fixed = T))],
       x = which(grepl(".2",colnames(Xorg),fixed = T)), col = 4, pch = 19)

