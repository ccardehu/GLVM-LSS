
# sil = T; form = s.form; restrl = restr; coefs = lc.; saveRes = T; cleanRes = F
# coefsf = lc; restrf = l1.

splvm.simfit <- function(l){
 tryCatch({
 AA   <- splvm.sim(n = n, form = s.form, fam = fam, constraints = l1., coefs = lc.)
 tmu  <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, control = list(method = "ML", silent = T, constraint = irestr, information = "Hessian"))
 tmp1 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, control = list(method = "PML", information = "Hessian", constraint = irestr,
                   silent = T, pml.control = list(type = "alasso", lambda = "auto", w.alasso = tmu$b, gamma = 1.0, a = 2)))
 tmp2 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, control = list(method = "PML", information = "Hessian", constraint = irestr,
                   silent = T, pml.control = list(type = "alasso", lambda = "auto", w.alasso = tmu$b, gamma = 1.4, a = 2)))
 tmp3 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, control = list(method = "PML", information = "Hessian", constraint = irestr,
                   silent = T, pml.control = list(type = "alasso", lambda = "auto", w.alasso = tmu$b, gamma = 2.0, a = 2)))
 tmp4 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, control = list(method = "PML", information = "Hessian", constraint = irestr,
                   silent = T, pml.control = list(type = "alasso", lambda = "auto", w.alasso = tmu$b, gamma = 3.0, a = 2)))
 tmp5 <- splvm.fit(Y = AA$Y, fam = fam, form = AA$form, control = list(method = "PML", information = "Hessian", constraint = irestr,
                   silent = T, pml.control = list(type = "alasso", lambda = "auto", w.alasso = tmu$b, gamma = 4.0, a = 2)))
 r0 <- c(lb2cb(AA$b)); ru <- c(lb2cb(tmu$b)); rp1 <- c(lb2cb(tmp1$b)); rp2 <- c(lb2cb(tmp2$b)); rp3 <- c(lb2cb(tmp3$b))
 rp4 <- c(lb2cb(tmp4$b)); rp5 <- c(lb2cb(tmp5$b));
 
 rlist <- c(tmp1$pml.control$lambda, tmp2$pml.control$lambda, tmp3$pml.control$lambda,
            tmp4$pml.control$lambda, tmp5$pml.control$lambda, length(r0))
 return(c(r0,ru,rp1,rp2,rp3,rp4,rp5,
          GBIC(tmu),GBIC(tmp1),GBIC(tmp2),GBIC(tmp3),GBIC(tmp4),GBIC(tmp5),
          rlist))
 },
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
