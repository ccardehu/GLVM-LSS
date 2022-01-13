
rm(list = ls())
# set.seed(1234)
source("f0_prep.R")

n = 500 # Number of individuals

# Simulation 1: Heteroscedastic Normal model:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s.form <- list("mu" = "~ Z1 + Z2"); p = 10
fam <- rep("binomial",p)
l1 <- lc <- NULL
l1$mu <- lc$mu <- matrix(1,ncol = 3, nrow = p)
# l1$sigma <- lc$sigma <- matrix(1, ncol = 1, nrow = p)
lc$mu[,1] <- runif(p,1,2)
lc$mu[,c(2,3)] <- runif(length(l1$mu[,c(2,3)]),0.5,2.5)# * sample(c(-1,1),size = 2*p,replace = T)
if(lc$mu[1,2] < 0) lc$mu[1,2] <- -lc$mu[1,2]; if(lc$mu[2,3] < 0) lc$mu[2,3] <- -lc$mu[2,3]

l1. <- l1

l1.$mu[1,c(3)] <- l1.$mu[2,c(2)] <- 0
l1.$mu[sample(3:p, floor((p-2)/2), replace = F),2] <- 0
for(i in 3:nrow(l1.$mu)){ if(l1.$mu[i,2] != 0) l1.$mu[i,3] <- rbinom(1,1,0.5) }

lc. <- NULL; for(i in names(l1.)){ lc.[[i]] <- lc[[i]]*l1.[[i]] }
simR <- splvm.sim(n = n, form = s.form, fam = fam, constraints = l1., coefs = lc.)
Y <- simR$Y
Z <- simR$Z
borg <- simR$b
e.form <- s.form

irestr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0))
restr <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0), c("mu",4,"Z1",0), c("mu",6,"Z1",0), c("mu",7,"Z1",0), c("mu",9,"Z1",0),
              c("mu",3,"Z2",0), c("mu",5,"Z2",0))

f0 <- splvm.fit(Y,fam,e.form,control = list(method = "ML", silent = F, information = "Hessian", constraint = restr))
fu <- splvm.fit(Y,fam,e.form,control = list(method = "ML", silent = F, information = "Hessian", constraint = irestr))

cbind(fu$b$mu,f0$b$mu,borg$mu)
GBIC(fu); GIC(fu);
GBIC(f0); GIC(f0);

lista <- c(1:2)
listg <- sort(c(1,1.4,seq(2,9,by = 1),log(n)/2))

GBICmat <- matrix(NA,nrow = length(listg), ncol = length(lista))
rownames(GBICmat) <- listg; colnames(GBICmat) <- lista
upLLmat <- pLLmat <- Lmat <- GICmat <- GBICmat

for(i in 1:length(lista)){
 for(j in 1:length(listg)){
   cat(paste0("\nComputing for a = ",lista[i],", gamma = ",listg[j]))
   tmp <- splvm.fit(Y,fam,e.form, control = list(method = "PML", information = "Hessian", constraint = irestr, silent = T,start.val = fu$b,
                    pml.control = list(type = "alasso", lambda = "auto", a = lista[i], gamma = listg[j], w.alasso = fu$b)))
   GICmat[j,i] <- GIC(tmp); 
   GBICmat[j,i] <- GBIC(tmp); 
   Lmat[j,i] <- tmp$pml.control$lambda
   upLLmat[j,i] <- tmp$uploglik; pLLmat[j,i] <- tmp$loglik
 }
}

# which(GBICmat == min(GBICmat), arr.ind = TRUE)
# which(upLLmat == max(upLLmat), arr.ind = TRUE)
# which(pLLmat == max(pLLmat), arr.ind = TRUE)

Lmat2 <- seq(from = max(0,min(Lmat)-1e-5), to=max(Lmat)+1e-2,length.out = 40) # unique(sort(union(seq(from = max(0,min(Lmat)-1e-5), to=max(Lmat)+1e-1,length.out = 14),Lmat)))
{upLLmat2 <- pLLmat2 <- GBICmat2 <- GICmat2 <-
 # upLLmat3 <- pLLmat3 <- GBICmat3 <- GICmat3 <- 
 # upLLmat4 <- pLLmat4 <- GBICmat4 <- GICmat4 <- 
    matrix(NA,nrow = length(Lmat2), ncol = ncol(GBICmat))}

# crsq <- function(a,Lmat){
#   tL <- seq(from = max(0,min(Lmat)-1e-5), to=max(Lmat)+1e-2,length.out = 39)
#   # tLu[1] <- max(tLu)
#   # tLu <- seq(from = a, length.out = 5, by = 1e-4)
#   # for(i in 2:length(tLu)){tLu[i] <- tLu[1] + 1.5*log(i)}
#   # for(i in 1:length(tLu)){tLu[i] <- tLu[1] + exp(0.75*i)*(max(tLl) - tLl[5-i+1])}
#   tL <- sort(union(tL,a))
#   return(unique(tL))
# }


for(i in seq_along(lista)){
 for(j in seq_along(Lmat2)){
  # if(k == 1) cat(paste0("\nComputing for a = ",lista[i],", gamma = ",format(round(listg[j],1), nsmall = 1),", \U03bb = ", tmpL[k])) else 
  # cat(paste0("\rComputing for a = ",lista[i],", gamma = ",format(round(listg[j],1), nsmall = 1),", \U03bb = ", tmpL[k]))
  cat(paste0("\nComputing for a = ",lista[i],", \U03bb = ",format(round(Lmat2[j],5), nsmall = 5)))
  tmp <- splvm.fit(Y,fam,e.form, control = list(method = "PML", information = "Hessian", constraint = irestr, silent = T,start.val = fu$b,
                   pml.control = list(type = "alasso", lambda = Lmat2[j], a = lista[i], w.alasso = fu$b)))
  GBICmat2[j,i] <- GBIC(tmp);
  GICmat2[j,i] <- GIC(tmp); 
  upLLmat2[j,i] <- tmp$uploglik
  pLLmat2[j,i] <- tmp$loglik
  
  # tmp <- splvm.fit(Y,fam,e.form, control = list(method = "PML", information = "Hessian", constraint = irestr, silent = T,start.val = fu$b,
  #                  pml.control = list(type = "mcp", lambda = Lmat2[j], a = lista[i])))
  # GBICmat3[j,i] <- GBIC(tmp);
  # GICmat3[j,i] <- GIC(tmp); 
  # upLLmat3[j,i] <- tmp$uploglik
  # pLLmat3[j,i] <- tmp$loglik
  # 
  # tmp <- splvm.fit(Y,fam,e.form, control = list(method = "PML", information = "Hessian", constraint = irestr, silent = T,start.val = fu$b,
  #                  pml.control = list(type = "scad", lambda = Lmat2[j], a = lista[i])))
  # GBICmat4[j,i] <- GBIC(tmp);
  # GICmat4[j,i] <- GIC(tmp); 
  # upLLmat4[j,i] <- tmp$uploglik
  # pLLmat4[j,i] <- tmp$loglik
 }
}


plot(0,0, xlim = c(min(Lmat2),max(Lmat2)), ylim = c(min(GBICmat)-5,max(GBICmat2)), type = "l", ylab = "GBIC", xlab = "Tuning parameter (\U03bb)",
     main = "GBIC values accross tuning parameters (\U03bb)")
for(i in 1:ncol(GBICmat2)){ lines(Lmat2,GBICmat2[,i], lty = i, lwd = 2) }
legend("topleft", legend=c("a = 1", "a = 2"), inset = 0.05, box.col="white",
       lty=1:2, cex=0.8, lwd = 2)

abline(h = min(GBICmat2[,1]), lty = 1, col = 1); abline(h = min(GBICmat2[,2]), col = 1, lty = 2)
abline(v = log(n)/(2*n), lty = 3)


for(i in 1:ncol(GBICmat)){
  for(j in 1:(nrow(GBICmat))){
    lines(Lmat[j,i], GBICmat[j,i], pch = i*5, col = i*5, type = "b",lwd = 1)
    text(Lmat[j,i], GBICmat[j,i], label = round(as.numeric(rownames(Lmat)[j]),1), pos = i, cex = 0.8, col = i*5)
  }
}

legend("bottomright", legend=c("a = 1", "a = 2"), inset = 0.05, box.col="white",
       pch=c(5,10), cex=0.8, col=c(5,10,15))

# plot(0,0, xlim = c(min(Lmat2),0.02), ylim = c(min(GICmat)-2, 5250), type = "l", ylab = "GIC", xlab = "lambda")
# for(i in 1:ncol(Lmat)){
#   lines(Lmat2,GICmat2[,i], col = i+1, lwd = 2)
# }
# abline(h = min(GICmat2[,1]), lty = 3, col = 2); abline(h = min(GICmat2[,2]), col = 3, lty = 3); abline(h = min(GICmat2[,3]), col = 4, lty = 3)


# plot(Lmat2[1,1,],upLLmat2[1,1,], xlim = c(min(Lmat2),max(Lmat2)), ylim = c(min(upLLmat2)-2, max(upLLmat2)+2), type = "l")
# # points(Lmat2[1,1,],GBICmat2[1,1,], xlim = c(min(Lmat2),max(Lmat2)), ylim = c(min(GBICmat2)-3, max(GBICmat2)+3))
# for(i in 1:ncol(Lmat)){
#  for(j in 1:nrow(Lmat)){
#   lines(Lmat2[j,i,],upLLmat2[j,i,], col = i)
#  }
# }
# points(Lmat[which.min(upLLmat[,1]),1], min(upLLmat[,1]), col = 1, pch = 3)
# points(Lmat[which.min(upLLmat[,2]),2], min(upLLmat[,2]), col = 2, pch = 3)
# points(Lmat[which.min(upLLmat[,3]),3], min(upLLmat[,3]), col = 3, pch = 3)
# 
# plot(Lmat2[1,1,],pLLmat2[1,1,], xlim = c(min(Lmat2),max(Lmat2)), ylim = c(min(pLLmat2)-2, max(pLLmat2)+2), type = "l")
# # points(Lmat2[1,1,],GBICmat2[1,1,], xlim = c(min(Lmat2),max(Lmat2)), ylim = c(min(GBICmat2)-3, max(GBICmat2)+3))
# for(i in 1:ncol(Lmat)){
#  for(j in 1:nrow(Lmat)){
#   lines(Lmat2[j,i,],pLLmat2[j,i,], col = i)
#  }
# }
# points(Lmat[which.min(pLLmat[,1]),1], min(pLLmat[,1]), col = 1, pch = 3)
# points(Lmat[which.min(pLLmat[,2]),2], min(pLLmat[,2]), col = 2, pch = 3)
# points(Lmat[which.min(pLLmat[,3]),3], min(pLLmat[,3]), col = 3, pch = 3)


plot(t(sapply(1:length(Lmat),function(i){crsq(Lmat[i],Lmat)}))[1,], type = "l", ylim = c(min(Lmat),max(Lmat)))
color <- rainbow(14)
for(i in 2:15){
  lines(1:9,t(sapply(2:(length(Lmat)),function(j){crsq(Lmat[j],Lmat)}))[i,], col = i-1)
}


plot(t(sapply(1:9,function(i){GBICmat2[,,i]}))[1,], type = "l", ylim = c(min(GBICmat2)-5,max(GBICmat2)+2))
color <- rainbow(8)
for(i in 1:8){
  lines(1:15,t(sapply(1:9,function(i){GBICmat2[,,i]}))[i+1,], col = i+1)
}
points(1:15, GBICmat, pch = 3, col = 2)


t(sapply(1:length(Lmat),function(i) crsq(Lmat[i])))
t(sapply(1:length(Lmat),function(i) crsq(Lmat[i])))[,5] == c(Lmat)

round(cbind(GBICmat,GBICmat2),3)
GBICmat < GBICmat2
round(cbind(Lmat,Lmat2),4)
Lmat > Lmat2
(GBICmat < GBICmat2) & (Lmat > Lmat2)

max(Lmat2 - Lmat)

testLmat <- Lmat[matrix(GBICmat %in% apply(GBICmat,2,min),nrow(GBICmat),ncol(GBICmat))]


or <- list(c("mu",1,"Z2",0), c("mu",2,"Z1",0), c("mu",3,"Z2",0), c("mu",4,"Z1",0), c("mu",5,"Z2",0), c("mu",6,"Z1",0), c("mu",7,"Z1",0),
           c("mu",9,"Z1",0))
fo <- splvm.fit(Y,fam,e.form,control = list(method = "ML", silent = F, information = "Hessian", constraint = or))
fn <- splvm.fit(Y,fam,e.form, control = list(method = "PML", information = "Hessian", constraint = irestr, silent = F,
                                            pml.control = list(type = "alasso", lambda = "auto", a = 2, gamma = 2, w.alasso = fu$b)))
# sval <- list("mu"= fn$b$mu + matrix(runif(length(fu$b$mu),-1,1),nrow = nrow(fu$b$mu), ncol = ncol(fu$b$mu)))
ff <- splvm.fit(Y,fam,e.form, control = list(method = "PML", information = "Hessian", constraint = irestr, silent = F, #start.val = sval,
                                            pml.control = list(type = "alasso", lambda = Lmat2[3,2], a = 2, gamma = 2, w.alasso = fu$b)))

fn$pml.control$lambda == Lmat[3,2] # test for automatic consistency
GBIC(fn); GBIC(ff); 
GBICmat[3,2]; GBICmat2[3,2]


GBIC(fo);GBIC(fu);GBIC(fn); GBIC(ff)
GIC(fo);GIC(fu);GIC(fn); GIC(ff)
round(c(fo$log,fu$log,fn$log,ff$log),3)

for(i in names(s.form)){ print(round(cbind(fo$b[[i]],fu$b[[i]],fn$b[[i]],ff$b[[i]],borg[[i]]),3)) }


library(ltm)

simTest <- function(l){

simR <- splvm.sim(n = n, form = s.form, fam = fam, constraints = l1., coefs = lc.)
fl <- ltm(simR$Y ~ z1 + z2, constr = rbind(c(1, 3, 0), c(2, 2, 0)), IRT.param = F)
ff <- splvm.fit(simR$Y,simR$fam,simR$form,control = list(method = "ML", silent = F, information = "Fisher", constraint = irestr))
fh <- splvm.fit(simR$Y,simR$fam,simR$form,control = list(method = "ML", silent = F, information = "Hessian", constraint = irestr))

r0 <- lb2cb(simR$b)
if(fl$coefficients[2,3] < 0){ r1 <- lb2cb(list(mu = cbind(fl$coefficients[,1:2],-fl$coefficients[,3]))) } else {
  r1 <- lb2cb(list(mu = fl$coefficients)) }
r2 <- lb2cb(ff$b)
r3 <- lb2cb(fh$b)

return(c(r0,r1,r2,r3))

}

cores <- detectCores()
cl<-makeCluster(cores)
registerDoSNOW(cl)
progress <- function(n) setTxtProgressBar(txtProgressBar(max = 50, style = 3), n)
opts <- list(progress = progress)

rmTest <- foreach(l = 1:50,
                .combine = rbind,
                .packages = c("mvtnorm","fastGHQuad","gamlss","trust","ltm"),#,
                .options.snow = opts
                ) %dopar% simTest(l)

stopCluster(cl)

nam <- NULL
for(i in names(s.form)){ nam <- rbind(nam,expand.grid(i,1:p,stringsAsFactors = F)) }
nam <- unlist(lapply(1:nrow(nam),function(i) paste0(nam[i,], collapse = "")))

nM <- NULL
for(i in 1:p){
 for(j in names(s.form)){
  nM <- append(nM,paste0(nam[grepl(j,as.character(nam))][i],".",c(0,seq_len(length(all.vars(as.formula(s.form[[j]])))))))
 }
}; 

colnames(rmTest) <- rep(nM,4)

rorg <- rmTest[,1:30]; rltm <- rmTest[,31:60];
rfis <- rmTest[,61:90]; rhes <- rmTest[,91:120]

AB0ltm <- mean(abs(colMeans(rltm[,grepl(".0",nM,fixed = T)]) - colMeans(rorg[,grepl(".0",nM,fixed = T)])))
AB1ltm <- mean(abs(colMeans(rltm[,!grepl(".0",nM,fixed = T)]) - colMeans(rorg[,!grepl(".0",nM,fixed = T)])))
AB0fis <- mean(abs(colMeans(rfis[,grepl(".0",nM,fixed = T)]) - colMeans(rorg[,grepl(".0",nM,fixed = T)])))
AB1fis <- mean(abs(colMeans(rfis[,!grepl(".0",nM,fixed = T)]) - colMeans(rorg[,!grepl(".0",nM,fixed = T)])))
AB0hes <- mean(abs(colMeans(rhes[,grepl(".0",nM,fixed = T)]) - colMeans(rorg[,grepl(".0",nM,fixed = T)])))
AB1hes <- mean(abs(colMeans(rhes[,!grepl(".0",nM,fixed = T)]) - colMeans(rorg[,!grepl(".0",nM,fixed = T)])))

AB0ltm; AB0fis; AB0hes
AB1ltm; AB1fis; AB1hes

MSE0ltm <- mean(colMeans((rorg[,grepl(".0",nM,fixed = T)] - rltm[,grepl(".0",nM,fixed = T)])^2))
MSE1ltm <- mean(colMeans((rorg[,!grepl(".0",nM,fixed = T)] - rltm[,!grepl(".0",nM,fixed = T)])^2))
MSE0fis <- mean(colMeans((rorg[,grepl(".0",nM,fixed = T)] - rfis[,grepl(".0",nM,fixed = T)])^2))
MSE1fis <- mean(colMeans((rorg[,!grepl(".0",nM,fixed = T)] - rfis[,!grepl(".0",nM,fixed = T)])^2))
MSE0hes <- mean(colMeans((rorg[,grepl(".0",nM,fixed = T)] - rhes[,grepl(".0",nM,fixed = T)])^2))
MSE1hes <- mean(colMeans((rorg[,!grepl(".0",nM,fixed = T)] - rhes[,!grepl(".0",nM,fixed = T)])^2))
 
MSE0ltm; MSE0fis; MSE0hes
MSE1ltm; MSE1fis; MSE1hes




# syntax = 'Z1 =~ Y1 + 0*Y2 + Y3 + Y4 + Y5 + Y6 + Y7 + Y8 + Y9 + Y10
#           Z2 =~ 0*Y1 + Y2 + Y3 + Y4 + Y5 + Y6 + Y7 + Y8 + Y9 + Y10'
# alasso_fit <- penfa::penfa(model  = syntax,
#                      data   = Y,
#                      std.lv = T,
#                      pen.shrink = "alasso",gamma = 2,a.alasso = 2,meanstructure = T,
#                      orthogonal = T,int.ov.free = T)
# penfa::coef(alasso_fit)
# penfa::summary(alasso_fit)
