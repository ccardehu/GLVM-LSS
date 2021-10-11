
# Code
# File f0_prep.R: Loading required packages & sourcing functions
# Code written by: Camilo Cardenas-Hurtado (c.a.cardenas-hurtado@lse.ac.uk)

# rm(list = ls())
set.seed(1234)

# package.list <- c("MASS","abind","splines","mvtnorm","numDeriv","gamlss",
                  # "foreach","doSNOW","parallel","itertools","xtable")
# package.check <- package.list %in% installed.packages()[,"Package"]
# if(any(!package.check)){ install.packages(package.list[package.check]) }
# rm(list = c(package.list, package.check))

# library(package.list)

library(MASS)
library(abind)
library(splines)
library(mvtnorm)
library(numDeriv)
library(gamlss)
library(foreach)
library(doSNOW)
library(parallel)
library(itertools)
library(xtable)
library(fastGHQuad)

source("fa_main.R")
source("fb_auxf.R")
source("fc_mfit.R")
source("fd_mghq.R")
source("ff_msim.R")
source("fi_esim.R")