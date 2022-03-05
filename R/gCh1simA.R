Ch1sim <- function(l,CM = F){
 tryCatch({
  AA <- splvm.sim(n = n, form = form, fam = fam, constraints = l1., coefs = lc.)
  infS <- ifelse(n == 200, "Fisher", "Hessian")
  iLim <- ifelse(n == 200, 1000, 1000)
  t1 <- splvm.fit(Y = AA$Y, fam = fam, form = form, control = list(method = "ML", start.val = lc.,
                  silent = T, constraint = irestr, information = infS, iter.lim = iLim))
  cof <- c(c(lb2cb(AA$b)),c(lb2cb(t1$b)))
  ex1 <- c(GBIC(t1),GIC(t1))
  ex2 <- c(t1$iter,length(c(lb2cb(t1$b))))
  if(CM){
    t2 <- splvm.fit(Y = AA$Y, fam = afam, form = aform, control = list(method = "ML",
                    silent = T, constraint = arestr, information = infS, iter.lim = iLim))
    cof <- c(cof,lb2cb(t2$b))
    ex1 <- c(ex1,GBIC(t2),GIC(t2))
    ex2 <- c(ex2,t2$iter,length(c(lb2cb(t2$b))))
    ex2 <- c(ex2[1],ex2[3],ex2[2],ex2[4])
  }
  return(c(cof,ex1,ex2))
  },
  error = function(e) {return(NA)}) # return(NA) stop(e)
}
