library(Rcpp)
rm(list = ls())
set.seed(1234)

sourceCpp("C:/Users/carde/OneDrive/Desktop/LSE/Research/Research/GitHub/SPLVM/C/testFun.cpp")


rbenchmark::benchmark(
  "Rt" = { dY(Y,ghQ,b,fam) },
  "Ct" = { DY(Y,ghQ,b,fam) },
  replications = 10,
  columns = c("test","elapsed","relative")
)


all.equal(dY(Y,ghQ,b,fam), DY(Y,ghQ,b,fam))
