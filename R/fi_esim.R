
# sil = T; form = s.form; restrl = restr; coefs = lc.; saveRes = T; cleanRes = F
# coefsf = lc; restrf = l1.

splvm.simfit <- function(l,ga){
 tryCatch({
 AA <- splvm.sim(n = n, form = s.form, fam = fam, constraints = l1., coefs = lc.)
 tmp0 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, control = list(method = "ML", constraint = restr, silent = T))
 tmp1 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, control = list(method = "ML", silent = T))
 tmp2 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, control = list(method = "PML", information = "Hessian", full.hess = T,
                   silent = T, pml.control = list(type = "alasso", lambda = "auto", w.alasso = tmp1$b, gamma = ga)))
 r0 <- c(lb2cb(tmp0$b)); r1 <- c(lb2cb(tmp1$b)); r2 <- c(lb2cb(tmp2$b))
 rlist <- c(tmp2$pml.control$lambda,tmp0$iter,tmp1$iter,tmp2$iter,tmp2$iiter,length(r0),length(r0)+length(r1),length(r0)+length(r1)+length(r2))
 return(c(r0,r1,r2,GBIC(tmp0),GBIC(tmp1),GIC(tmp2),GIC(tmp0),GIC(tmp1),GIC(tmp2),rlist))},
 error = function(e) {return(NA)})
}

# splvm.simfit <- function(l){ # ,sil = sil
# 
# tryCatch({
# # if(ncol(AA$Z) == 1) nqp <- 40 else nqp <- 20
# AA <- splvm.sim(n = n, form = form, fam = fam, constraints = restrf, coefs = coefs)
# tmp0 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, control = list(method = "EM", constraint = restrl, silent = T))
# tmp1 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, control = list(method = "EM", silent = T))
# tmp2 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, control = list(method = "PEM", 
#                                                               silent = T, pml.control = list(type = "alasso", lambda = "auto", w.alasso = tmp1$b, gamma = 3)))
# # if(!sil) cat("Run = ", l, "Iter. to Convg.: ",tmp1$iter,"\n")
# # return(c(lb2cb(tmp0$b),lb2cb(tmp1$b),lb2cb(tmp2$b),tmp0$iter,tmp1$iter,tmp2$iter))
# r0 <- c(lb2cb(tmp0$b)); r1 <- c(lb2cb(tmp1$b)); r2 <- c(lb2cb(tmp2$b))
# rlist <- c(tmp2$pml.control$lambda,tmp0$iter,tmp1$iter,tmp2$iter,tmp2$iiter,length(r0),length(r0)+length(r1),length(r0)+length(r1)+length(r2))
# return(c(r0,r1,r2,GBIC(tmp0),GBIC(tmp1),GIC(tmp2),GIC(tmp0),GIC(tmp1),GIC(tmp2),rlist))},
# error = function(e) {return(NA)})
# }
