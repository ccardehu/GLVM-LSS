library(lavaan)
library(numDeriv)
#1101 update
#weights for IRLS
#sgnmodeltext and Q matrix
#bic_ex
#sgnrefit_optim2 with mode
#fixed the time bug in irls
#1116 fixed 2 bugs for irls, 1 early stopping if check L1, 2 weights has A rather than L, 3 new irls_ep=1e-5
#1116 gpa add a step to check whether the line search is successful
#1121
#re-parametrise the function with B be a upper.tri with diagonals 1,
# since this parametrisation is different from the rotation problem, we first transform B back with unit columns and then use jacobian for the inference
# checked the results of obj value and CIs, which is unchanged
#added 4 auxilary functions for sgnrefit_optim2B2, rho,invrho,x2par,par2x
#have a new permutation and sign flip function pefl.LT that simultaneously have L and B(since we use Dinv, I added several line to ensure D is invertible
# i.e, can not have 1 on the same row or columns)
# 
#------------------------------------- irls algo function -------------------------------------- 
Q<-function(w,L){return(sum(w*L*L))};
dQ<-function(w,L){return(2*w*L)}
rho<-function(x){return(x%*%diag(diag(t(x)%*%x)^(-1/2)))};
Q1<- function(L,p){return(sum(abs(L)^p))}

gpa <-function(A,TR,W,iter,stop_ep,alpha0=1){
  # Q is the weighted least square function
  alpha<-alpha0;
  Ti<-solve(TR);
  L<-A%*%t(Ti);
  ft<-Q(W,L);
  G<- -t(t(L)%*%dQ(W,L)%*%Ti);
  flag=0
  for (j in 1:iter){   
    alpha<-2*alpha;
    # find a proper step size
    for( i in 1:20 ){
      X<-TR-alpha*G;
      T_new<-rho(X);
      L<-A%*%solve(t(T_new));
      if (Q(W,L)<ft){
        ft<-Q(W,L);
        flag=1
        break;
      }else{
        alpha<-alpha*0.5;
      }
    }
    if (flag==1){
      # update gradient and T
      TR<-  T_new;
      G<-  -t(t(L)%*%dQ(W,L)%*%solve(TR));}
    else{
      break;
    }
  }
  T_final<-TR;
  L_final<-A%*%solve(t(T_final));
  return(list(T_final,L_final,ft,j))
}

# iterative re-weighted least square algorithm
irls<-function(p,A,irls_ep0=1e-5,
               stop_ep=1e-6,
               irls_iter=10000,gpa_iter=5,
               T=diag(rep(1,ncol(A))),W=(A^2+irls_ep0)^(p/2-1),
               alpha0=1,power=1/2){
  
  j_sum<-0
  start_time<-proc.time()
  for (it in 1:irls_iter){
    # solve by gpa
    gpa_result <- gpa(A,T,W,gpa_iter,stop_ep,alpha0);
    T_new <-gpa_result[[1]]
    L <-gpa_result[[2]]
    ft <-gpa_result[[3]]
    j <-gpa_result[[4]] 
    j_sum<-j_sum+j
    #repeatedly update weights in the objective function
    W <- (L^2+irls_ep0)^(p/2-1) 
    if ((max(abs(L-A%*%solve(t(T))))<stop_ep)){
      break
    }
    T<-T_new;
  }
  return(list(T=T,L=L,
              ft=ft,
              it=it,gpa_it=j_sum/it,
              t=proc.time()[3]-start_time[3],
              obj=Q1(L,p)))
}
# # gradient projection algorithm
# gpa <-function(A,TR,W,iter,stop_ep,alpha0=1){
#   # Q is the weighted least square function
#   alpha<-alpha0;
#   Ti<-solve(TR);
#   L<-A%*%t(Ti);
#   ft<-Q(W,L);
#   G<- -t(t(L)%*%dQ(W,L)%*%Ti);
#   
#   for (j in 1:iter){   
#     alpha<-2*alpha;
#     # find a proper step size
#     for( i in 1:20 ){
#       X<-TR-alpha*G;
#       T_new<-rho(X);
#       L<-A%*%solve(t(T_new));
#       if (Q(W,L)<ft){
#         ft<-Q(W,L);
#         break;
#       }else{
#         alpha<-alpha*0.5;
#       }
#     }
#     
#     # update gradient and T
#     TR<-   T_new;
#     G<-  -t(t(L)%*%dQ(W,L)%*%solve(TR));
#   }
#   T_final<-TR;
#   L_final<-A%*%solve(t(T_final));
#   #Phi<-t(T)%*%T;
#   return(list(T_final,L_final,ft,j))
# }
# 
# 
# 
# 
# # iterative re-weighted least square algorithm
# irls<-function(p,A,irls_ep0=1e-3,
#                stop_ep=1e-6,
#                irls_iter=10000,gpa_iter=5,
#                T=diag(rep(1,ncol(A))),W=(A^2+irls_ep0)^(p/2-1), #(abs(A)+mean(abs(A))*irls_ep0)^(p-2), #W=pmax(abs(A),irls_ep0)^(p-2),#,#
#                alpha0=1,power=1/2){
#   
#   j_sum<-0
#   start_time<-proc.time()
#   for (it in 1:irls_iter){
#     #irls_ep0<-0.1/it
#     #irls_ep=0.001/it
#     #irls_ep=min(irls_ep0,1/it)
#     # solve by gpa
#     gpa_result <- gpa(A,T,W,gpa_iter,stop_ep,alpha0);
#     T_new <-gpa_result[[1]]
#     L <-gpa_result[[2]]
#     ft <-gpa_result[[3]]
#     j <-gpa_result[[4]] 
#     
#     #Phi <-gpa_result[[5]]
#     j_sum<-j_sum+j
#     #repeatedly update weights in the objective function
#     #W <- (abs(L)+irls_ep)^(p-2);
#     #W <- (abs(L)+(mean(abs(L))*irls_ep0)/(it)^(power))^(p-2);
#     W <- (A^2+irls_ep0)^(p/2-1) #(abs(L)+mean(abs(L))*irls_ep0)^(p-2)
#     if (Q1(A%*%solve(t(T)),p)-Q1(L,p)< 0){
#       L<-A%*%solve(t(T))
#       break
#     }
#     else if ((max(abs(L-A%*%solve(t(T))))<stop_ep)){
#       break
#     }
#     T<-T_new;
#   }
#   return(list(T=T,L=L,
#               ft=ft,
#               it=it,gpa_it=j_sum/it,
#               t=proc.time()[3]-start_time[3],
#               obj=Q1(L,p)))
# }

#------------------------------------- other functions -------------------------------------- 
#rho <- function(x){return(x%*%diag(diag(t(x)%*%x)^(-1/2)))};# normalization function

hard.each<-function(c,L_irls,L_old){
  L_new=abs(L_irls)>c
  TPR= sum(L_new[L_old==1])/(sum(L_old))
  TNR= 1-sum(L_new[L_old==0])/(length(L_old)-sum(L_old))
  TR=sum(L_new==L_old)/(nrow(L_old)*ncol(L_old))
  TRall=all(L_new==L_old)
  return(list(TPR=TPR,TNR=TNR,TR=TR,TRall=TRall))
}

hard<-function(c_list,L_irls,L_old){
  
  ans_hard.each=sapply(c_list,hard.each,L_irls=L_irls,L_old=L_old)
  TPR=unlist(ans_hard.each[1,]) #length of c_list
  TNR=unlist(ans_hard.each[2,])
  id <- 1:(length(TNR)+1)
  AUC <- sum(diff((1-c(1,TPR))[id])*rollmean(c(0,TNR)[id],2)) #length of 1
  return(list(TPR=TPR,TNR=TNR,AUC=AUC))
}

# # permutation and sign flip of the solution
# permfilp <- function(L_irls,L_list){
#   # sign filp
#   n_col=ncol(L_irls)
#   for( i in 1:n_col){
#     if(sum(L_irls[,i])<0){
#       L_irls[,i]<- -L_irls[,i]}}
#   # column swap
#   #library(gtools)
#   P<-permutations(n_col, n_col)
#   minf<-norm(L_list-L_irls,'f')
#   L_new0<-L_irls
#   for (i in 1:nrow(P)){
#     if (norm(L_list-L_irls[,P[i,]],'f')< minf){
#       minf<-norm(L_list-L_irls[,P[i,]],'f')
#       L_new0<-L_irls[,P[i,]]
#     }
#   }
#   return(L_new0)
# }


# parameterEstimates(fit)$ci.upper
modeltext<-function(L_tmp){
  n_row=nrow(L_tmp)
  n_col=ncol(L_tmp)
  main=''
  #std.lv=T
  for (i in 1:n_col){
    t=0
    for (j in 1:n_row)
    {if(L_tmp[j,i]!=0){
      if(t==0){
        tmp=paste0('y',i,'=~','NA*','x',j)
      }else{
        tmp=paste0('+','x',j)
      }
      t=1 
      main=paste0(main, tmp)
      
    } 
    } 
    main=paste0(main,'\n')
  }
  main=paste0(main, '# orthogonal factors\n')
  for ( i in 1:n_col){
    if (!all(!L_tmp[,i])){
      main = paste0(main,'y',i,'~~','1*y',i,'\n')}
  }
  main=paste0(main, '# constrains\n')
  for ( j in 1:n_row){
    main = paste0(main,'x',j,'~~','c',j,'*','x',j,'\n')
    main = paste0(main,'c',j,'>0','\n')}
  # cat(main)
  return(main)
}

Gen.SLP<-function(Sig,N,n_col){
  n_row<-nrow(Sig)
  X <-t(chol(Sig)) %*% matrix(rnorm(N*n_row),n_row,N)
  S <- cov(t(X))
  var_nam=c()
  for( j in 1:n_row){
    var_nam<-c(var_nam,paste0('x',j))}
  colnames(S)[1:n_row] <- var_nam
  rownames(S)[1:n_row] <- var_nam
  
  L_tmp<-matrix(1,n_row,n_col)
  L_tmp[upper.tri(L_tmp)]<-0
  main<-modeltext(L_tmp)
  #cat(main)
  fit <-cfa(main,
            sample.cov=S,
            sample.nobs=N,
            orthogonal=T,
            se="none")
  ncoef.L<-(n_row*n_col-(n_col-1)*n_col/2)
  L_cfa<-matrix(0,n_row,n_col)
  L_cfa[lower.tri(L_cfa,diag=TRUE)]<-coef(fit)[1:ncoef.L]
  Psi_cfa<-coef(fit)[ncoef.L+1:n_row]
  L_vari<-L_cfa %*% varimax(L_cfa)$rotmat
  return(list(S=S,
              L_vari=L_vari,
              Psi_cfa=Psi_cfa))
}

# refit.bic<-function(c,L_irls,S,N){
#   L_zeros<-(abs(L_irls)>c)*sign(L_irls)
#   main<-sgnmodeltext(L_zeros)
#   #cat(main)
#   fit <-cfa(main,
#             sample.cov=S,
#             sample.nobs=N,
#             se="none")
#   n_row=nrow(L_zeros)
#   n_col=ncol(L_zeros)
#   bic_ex= 0 # log(choose(n_row*n_col,sum(L_zeros)))
#   return(bic=BIC(fit)+bic_ex)
# }
# 
# refit.L<-function(c,L_irls,S,N){
#   L_zeros<-(abs(L_irls)>c)*sign(L_irls)
#   main<-sgnmodeltext(L_zeros)
#   #cat(main)
#   fit <-cfa(main,
#             sample.cov=S,
#             sample.nobs=N,
#             se="none")
#   L_cfa<-matrix(0,n_row,n_col)
#   L_cfa[which(L_zeros!=0)]<-coef(fit)[1:sum(L_zeros)]
#   return(L_cfa)
# }

# ## objective function (in vector)
# fx<-function(x,n_row,n_col,L_zeros,S){
#   #reconstruct Lambda,Psi,Phi from vector x
#   Lnco = sum(L_zeros!=0)
#   Lambda<-matrix(0,n_row,n_col)
#   Lambda[which(L_zeros!=0)] = x[1:Lnco]
#   Psi = x[Lnco+1:n_row]
#   Phi = diag(n_col)
#   Phi[lower.tri(Phi)] = x[Lnco+n_row+1:(n_col*(n_col-1)/2)]
#   Phi[upper.tri(Phi)] = t(Phi)[upper.tri(Phi)]
#   Sigma  = Lambda%*%Phi%*%t(Lambda)+diag(exp(Psi))
#   return(g(Sigma,S))}
# 
# ## gradient of objective function (in vector)
# grad_fx<-function(x,n_row,n_col,L_zeros,S){
#   Lnco = sum(L_zeros!=0)
#   Lambda = matrix(0,n_row,n_col)
#   Lambda[which(L_zeros!=0)] = x[1:Lnco]
#   Psi = x[Lnco+1:n_row]
#   Phi = diag(n_col)
#   Phi[lower.tri(Phi)] = x[Lnco+n_row+1:(n_col*(n_col-1)/2)]
#   Phi[upper.tri(Phi)] = t(Phi)[upper.tri(Phi)]
#   Sigma  = Lambda%*%Phi%*%t(Lambda)+diag(exp(Psi))
#   Sigma_inv  = solve(Sigma)
#   Q = Sigma_inv-Sigma_inv%*%S%*%Sigma_inv
#   grad_Lambda =  Subg_Lambda(Q,Lambda,Phi)
#   grad_Psi1 =  Subg_Psi(Q,Psi)
#   grad_Phi = Subg_Phi1(Q,Lambda)
#   return(c(grad_Lambda[L_zeros!=0],grad_Psi1,grad_Phi[lower.tri(grad_Phi)]))}
# 
# sgnrefit_optim2<-function(c0,L_irls,S,N,mod=0,k=-1,Psi0=rep(0,n_row),Phi0=rep(0,n_col*(n_col-1)/2)){
#   #mode=0 return bic, mode = 1 return L, mode=2 return confidence interval
#   L_zeros=(abs(L_irls)>c0)*sign(L_irls)
#   if_hessian=F
#   if (mod==2){
#     L_zeros[k,] = 2
#   if_hessian=T}
#   n_row = nrow(L_irls)
#   n_col = ncol(L_irls)
#   
#   # input a vector L,Psi,B
#   Lnco = sum(L_zeros!=0)
#   L0 = L_irls[L_zeros!=0]
#   x0 = c(L0,Psi0,Phi0)
#   
#   # set sign constraints, where L_zeros=2 means no constraints(used in CI calculation)
#   # L_zeros>0 means positive sign, L_zeros<0 means negative sign
#   lowerb = rep(-Inf,length(x0))
#   lowerb[which(((L_zeros[L_zeros!=0])>0)&(L_zeros[L_zeros!=0]!=2))] = 10^(-3)
#   lowerb[Lnco+1:n_row]=-3
#   lowerb[Lnco+n_row+1:(n_col*(n_col-1)/2)]=-1
#   upperb = rep(Inf,length(x0))
#   upperb[which((L_zeros[L_zeros!=0])<0)] = -10^(-3)
#   upperb[Lnco+n_row+1:(n_col*(n_col-1)/2)]=1
#   
#   # solve the MLE by BFGS with sign constraints
#   result<-optim(x0, fx, grad_fx, n_row=n_row, n_col=n_col,L_zeros=L_zeros,S=S*(N-1)/N,
#                 method = "L-BFGS-B",lower = lowerb,upper=upperb, hessian = if_hessian)
#   x = result$par
#   if (mod==0){
#     bic = N*result$value+log(N)*length(x0)
#     return(bic)}
#   else if (mod == 1){
#     Lambda = matrix(0,n_row,n_col)
#     Lambda[which(L_zeros!=0)] = x[1:Lnco]
#     Psi = x[Lnco+1:n_row]
#     Phi = diag(n_col)
#     Phi[lower.tri(Phi)] = x[Lnco+n_row+1:(n_col*(n_col-1)/2)]
#     Phi[upper.tri(Phi)] = t(Phi)[upper.tri(Phi)]
#     return(list(Lambda=Lambda,Psi=Psi,Phi=Phi))}
#   else if (mod ==2 ){
#     # construct confidence interval
#     L_up = L_low = matrix(0,n_row,n_col)
#     se_est = sqrt(diag(solve(N*result$hessian/2))) #standard error
#     L_up[which(L_zeros!=0)] = (x+1.96*se_est)[1:Lnco]
#     L_low[which(L_zeros!=0)] = (x-1.96*se_est)[1:Lnco]
#     return(list(u=L_up,l=L_low))
#   }
# }
# this is used to calculate Confidence interval
grad_fx<-function(x,n_row,n_col,L_zeros,S){
  Lnco = sum(L_zeros!=0)
  Lambda = matrix(0,n_row,n_col)
  Lambda[which(L_zeros!=0)] = x[1:Lnco]
  Psi = x[Lnco+1:n_row]
  Phi = diag(n_col)
  Phi[lower.tri(Phi)] = x[Lnco+n_row+1:(n_col*(n_col-1)/2)]
  Phi[upper.tri(Phi)] = t(Phi)[upper.tri(Phi)]
  Sigma  = Lambda%*%Phi%*%t(Lambda)+diag(exp(Psi))
  Sigma_inv  = solve(Sigma)
  Q = Sigma_inv-Sigma_inv%*%S%*%Sigma_inv
  grad_Lambda =  Subg_Lambda(Q,Lambda,Phi)
  grad_Psi1 =  Subg_Psi(Q,Psi)
  grad_Phi = Subg_Phi1(Q,Lambda)
  return(c(grad_Lambda[L_zeros!=0],grad_Psi1,grad_Phi[lower.tri(grad_Phi)]))}

Subg_Bvec <- function(Q,Lambda,B){
  grad_B <-2*B %*% t(Lambda) %*% Q %*% Lambda
  grad_B[lower.tri(grad_B,diag = TRUE)] <- 0
  return(grad_B[upper.tri(grad_B)])
}


fxB<-function(x,n_row,n_col,L_zeros,S){
  #reconstruct Lambda,Psi,Phi from vector x
  Lnco = sum(L_zeros!=0)
  Lambda<-matrix(0,n_row,n_col)
  Lambda[which(L_zeros!=0)] = x[1:Lnco]
  Psi = x[Lnco+1:n_row]
  B = diag(n_col)
  B[upper.tri(B)]=x[Lnco+n_row+1:(n_col*(n_col-1)/2)]
  Sigma  = Lambda%*%t(B)%*%B%*%t(Lambda)+diag(exp(Psi))
  return(g(Sigma,S))}

## gradient of objective function (in vector)
grad_fxB<-function(x,n_row,n_col,L_zeros,S){
  Lnco = sum(L_zeros!=0)
  Lambda = matrix(0,n_row,n_col)
  Lambda[which(L_zeros!=0)] = x[1:Lnco]
  Psi = x[Lnco+1:n_row]
  B = diag(n_col)
  B[upper.tri(B)]=x[Lnco+n_row+1:(n_col*(n_col-1)/2)]
  
  Sigma  = Lambda%*%t(B)%*%B%*%t(Lambda)+diag(exp(Psi))
  Sigma_inv  = solve(Sigma)
  Q = Sigma_inv-Sigma_inv%*%S%*%Sigma_inv
  grad_Lambda =  Subg_Lambda(Q,Lambda,t(B)%*%B)
  grad_Psi1 =  Subg_Psi(Q,Psi)
  grad_B = Subg_Bvec(Q,Lambda,B)
  return(c(grad_Lambda[L_zeros!=0],grad_Psi1,grad_B))}

sgnrefit_optim2B2<-function(c0,L_irls,S,N,mod=0,k=-1,Psi0=rep(0,n_row),B=diag(n_col)){
  library(numDeriv)
  #mode=0 return bic, mode = 1 return L, mode=2 return confidence interval
  L_zeros=(abs(L_irls)>c0)*sign(L_irls)
  if_hessian=F
  if (mod==2){
    L_zeros[k,] = 2
    # if_hessian=T
  }
  n_row = nrow(L_irls)
  n_col = ncol(L_irls)
  
  # input a vector of L,Psi,B, where we require B to be upper.tri with diagonal 1.
  # the 0.01 diag is added to ensure B has full rank
  Lnco = sum(L_zeros!=0)
  resLB0=invrho(L_irls,rho(B+0.01*diag(n_col)))
  x0 = par2x(resLB0$L0,Psi0,resLB0$B0,L_zeros)
  
  # set sign constraints, where L_zeros=2 means no constraints(used in CI calculation)
  # L_zeros>0 means positive sign, L_zeros<0 means negative sign
  lowerb = rep(-Inf,length(x0))
  lowerb[which(((L_zeros[L_zeros!=0])>0)&(L_zeros[L_zeros!=0]!=2))] = 10^(-3)
  lowerb[Lnco+1:n_row]=-3
  upperb = rep(Inf,length(x0))
  upperb[which((L_zeros[L_zeros!=0])<0)] = -10^(-3)
  
  # solve the MLE by BFGS with sign constraints
  result<-optim(x0, fxB, grad_fxB, n_row=n_row, n_col=n_col,L_zeros=L_zeros,S=S*(N-1)/N,
                method = "L-BFGS-B",lower = lowerb,upper=upperb, hessian = if_hessian)
  x = result$par
  if (mod==0){
    bic = N*result$value+log(N)*length(x0)
    return(bic)}
  else {
    respar=x2par(x,n_row,n_col,L_zeros)
    L1=respar$L
    B1=respar$B
    resrho=rhoLB(L1,B1)
    L2=resrho$L
    B2=resrho$B
    if (mod == 1){
      return(list(L=L2,Psi=respar$Psi,B=B2))}
    else if (mod ==2 ){
      # construct confidence interval
      L_up = L_low = matrix(0,n_row,n_col)
      Phi2=t(B2)%*%B2
      x1= c(L2[L_zeros!=0],respar$Psi,Phi2[lower.tri(Phi2)]) #par2x(L2,respar$Psi,B2,L_zeros)
      se_est = sqrt(diag(solve(N*jacobian(grad_fx,x1, method="Richardson",n_row=n_row,n_col=n_col,L_zeros=L_zeros,S=S)/2))) #standard error
      L_up[which(L_zeros!=0)] = (x1+1.96*se_est)[1:Lnco]
      L_low[which(L_zeros!=0)] = (x1-1.96*se_est)[1:Lnco]
      return(list(u=L_up,l=L_low))
    }
  }}


par2x<-function(L,Psi,B,L_zeros){
  x = c(L[L_zeros!=0],Psi,B[upper.tri(B)])
  return(x)
}

# fix the bug that the D might degenerate
pefl.LT<-function(L_irls,L,T_irls){
  K=ncol(L)
  D=matrix(0,K,K)
  Llong=cbind(L,-L)
  sel=c()
  for (j in 1:K)
  {dis=Llong-L_irls[,j]
  Dis=sapply(1:(2*K),function(x)norm(dis[,x],"2"))
  #print(Dis)
  left=setdiff(1:(2*K),sel)
  i=left[which.min(Dis[left])]
  sel=c(sel,i)
  if (i<=K){
    D[j,i]=1
    sel=c(sel,i+K)
  }else{
    D[j,i-K]=-1
    sel=c(sel,i-K)
  }
  }
  Lnew=L_irls%*%D
  Tnew=T_irls%*%t(solve(D))
  return(list(L=Lnew,B=Tnew,D=D))
}
rhoLB<-function(L1,B1){
  B2 = B1%*%diag(diag(t(B1)%*%B1)^(-1/2))
  L2 = L1%*%diag(diag(t(B1)%*%B1)^(1/2))
  return(list(L=L2,B=B2))
}
x2par<-function(x,n_row,n_col,L_zeros){
  Lnco = sum(L_zeros!=0)
  L = matrix(0,n_row,n_col)
  L[which(L_zeros!=0)] = x[1:Lnco]
  Psi = x[Lnco+1:n_row]
  B = diag(n_col)
  B[upper.tri(B)]=x[Lnco+n_row+1:(n_col*(n_col-1)/2)]
  return(list(L=L,Psi=Psi,B=B))
}
invrho<-function(L_irls,T_irls){
  B=chol(t(T_irls)%*%T_irls)
  L0=L_irls%*%diag(diag(B))
  B0=(B%*%diag(diag(B)^(-1)))
  return(list(L0=L0,B0=B0))}

CI3<-function(c0,L_irls,S,N,Psi0, B){
  t1 = proc.time()
  L_upper=L_lower=matrix(0,nrow(L_irls),ncol(L_irls))
  for ( k in 1:nrow(L_irls)){
    res = sgnrefit_optim2B2(c0,L_irls,S,N,mod=2,k=k,Psi0=Psi0,B=B)
    L_lower[k,] = res$l[k,]
    L_upper[k,] = res$u[k,]
  }
  L_class<-(L_lower<=L)&(L<=L_upper)
  accuracy<-sum(L_class)/length(L_lower)
  return(list(L_class=L_class,L_upper=L_upper, L_lower= L_lower,
              ci.accuracy=accuracy,na.flag=!all(!is.infinite(L_lower)),
              t=(proc.time()-t1)[3]))
}

# refit.Lbic<-function(c,L_irls,S,N){
#   L_zeros<-(abs(L_irls)>c)*sign(L_irls)
#   main<-sgnmodeltext(L_zeros)
#   #cat(main)
#   fit <-cfa(main,
#             sample.cov=S,
#             sample.nobs=N,
#             se="none")
#   L_cfa<-matrix(0,n_row,n_col)
#   L_cfa[which(L_zeros!=0)]<-coef(fit)[1:sum(L_zeros)]
#   n_row=nrow(L_zeros)
#   n_col=ncol(L_zeros)
#   bic_ex= 0 #log(choose(n_row*n_col,sum(L_zeros)))
#   return(list(bic=BIC(fit)+bic_ex,L_bic=L_cfa))
# }

# CI_line<-function(k,L_zeros,S,N){
#   L_tmp<-L_zeros
#   L_tmp[k,]<-2
#   main<-sgnmodeltext(L_tmp)
#   fit <-cfa(main,sample.cov=S,sample.nobs=N) 
#   L_uppertmp<-matrix(0,n_row,n_col)
#   L_uppertmp[which(L_tmp!=0)]<-parameterEstimates(fit)$ci.upper[1:sum(L_tmp)]
#   L_lowertmp<-matrix(0,n_row,n_col)
#   L_lowertmp[which(L_tmp!=0)]<-parameterEstimates(fit)$ci.lower[1:sum(L_tmp)]
#   return( c(L_lowertmp[k,],L_uppertmp[k,]))
# }
# 
# CI<-function(L_bic,S,N,L){
#   n_row=nrow(L)
#   n_col=ncol(L)
#   bounds=sapply(1:n_row,CI_line,L_zeros=sign(L_bic),S=S,N=N)
#   L_lower=t(bounds[1:n_col,])
#   L_upper=t(bounds[n_col+1:n_col,])
#   L_lower[which(is.na(L_lower))]<--Inf
#   L_upper[which(is.na(L_upper))]<-Inf
#   L_class<-(L_lower<=L)&(L<=L_upper)
#   accuracy<-sum(L_class)/length(L_lower)
#   
#   return(list(L_class=L_class,
#               L_upper=L_upper, L_lower= L_lower,
#               ci.accuracy=accuracy,na.flag=!all(!is.infinite(L_lower)) ))
# }

# lasso.path<-function(lambda_list,L_vari,B,Psi_cfa,S,N,L){
#   nlambda<-length(lambda_list)
#   TPR=TNR=rep(0,nlambda)
#   L0_list=list()
#   t=it=0
#   L0=L_vari
#   B0=B
#   Psi0=Psi_cfa
#   bic.min=Inf
#   for (j in 1:nlambda){
#     ans_est<-prox_grad(L0,B0,Psi0,S,lambda_list[j],1)
#     L0_list[[j]]=L0=permfilp(ans_est$L,L)
#     Psi0=ans_est$Psi
#     t=ans_est$t+t
#     it=ans_est$it+it
#     
#     #soft-thresholding
#     ans_hard=hard.each(0,L0,L_old)
#     TPR[j] =ans_hard$TPR
#     TNR[j] =ans_hard$TNR
#     
#     if(lambda_list[j]){
#       ans_bic=refit.Lbic(0,L0,S,N)
#       bic=ans_bic$bic
#       if(bic<bic.min){
#         L_bic=ans_bic$L_bic
#         #print(L_bic)
#         l=lambda_list[j]
#         bic.min=bic
#       }
#     }
#   }
#   L_bic.res=hard.each(0,L_bic,L_old)
#   id <- 1:(length(TNR)+1)
#   AUC <- sum(diff((1-c(1,TPR))[id])*rollmean(c(0,TNR)[id],2))
#   return(list(L=L0_list,t=t,it=it, 
#               TPR=TPR,TNR=TNR,AUC=AUC,
#               c=l,L_bic=L_bic,L_bic.res=L_bic.res))
# }
#------------------------------------- lasso functions -------------------------------------- 
## the main function f = g + h
f <- function(Sigma, S, Lambda,lambda){
  log(det(Sigma))+sum(diag(S %*% solve(Sigma))) + lambda*sum(abs(Lambda))
}

## smooth function g
g <- function(Sigma,S){
  log(det(2*pi*Sigma))+sum(diag(S %*% solve(Sigma)))
}

g2 <- function(Sigma,S){
  nrow(Sigma)*log(2*pi)+sum(log(eigen(Sigma)$values))+sum(diag(S %*% solve(Sigma)))
}

# lkhd1 <- function(Lambda,B,Psi,S){
#   
#   Sigma <-Lambda%*%t(B)%*%B%*%t(Lambda)+diag(exp(Psi))
#   log(det(2*pi*Sigma))+sum(diag(S %*% solve(Sigma)))
# }
# 
# lkhd2 <- function(Lambda,T_i,Psi_cfa,S){
#   Sigma <-Lambda%*%t(B)%*%B%*%t(Lambda)+diag(Psi_cfa)
#   log(det(2*pi*Sigma))+sum(diag(S %*% solve(Sigma)))
# }

## subgradient of smooth function
Subg_Lambda <- function(Q,Lambda,Phi){
  2*Q %*% Lambda %*% Phi
}

## we use Phi instead of T to parameterise the likelihood 
Subg_Phi1 <- function(Q,Lambda){
  grad_Phi <-2* t(Lambda) %*% Q %*% Lambda
  return(grad_Phi)
}

Subg_B <- function(Q,Lambda,B){
  grad_B <-2*B %*% t(Lambda) %*% Q %*% Lambda
  grad_B[lower.tri(grad_B)] <- 0
  return(grad_B)
}

Subg_Psi <- function(Q,Psi){
  diag(Q)*exp(Psi)
}

## proximal function of none smooth function h
prox_L1 <- function(Lambda, lambda){
  sign(Lambda) * pmax(abs(Lambda) - lambda, 0)
}



# proximal gradient descend 
prox_grad<-function(Lambda0,B0,Psi0,S,lambda,if_fixB=0,
                    maxiter=10000,stop_ep=10^(-6),t0=1){
  
  Lambda<-Lambda0
  B<-B0
  Psi<-Psi0
  Sigma <-Lambda%*%t(B)%*%B%*%t(Lambda)+diag(exp(Psi))
  
  # initilization
  t <- t0                # step size
  beta0 <- 0.5             # line search factor for t parameter
  #time_iter<-rep(0, maxiter)
  #objective_iter<-rep(0,maxiter)
  start_time<-proc.time()
  for(i in 1:maxiter){
    
    t<-2*t
    Sigma_inv <-solve(Sigma)
    
    Q=Sigma_inv-Sigma_inv%*%S%*%Sigma_inv
    Phi=t(B)%*%B
    
    # gradient
    grad_Lambda <- Subg_Lambda(Q,Lambda,Phi)
    grad_B <- Subg_B(Q,Lambda,B)
    grad_Psi <- Subg_Psi(Q,Psi)
    
    # proximal gradient descend stage and line search
    for( i_iter in 1:20 ){
      
      Lambda_new <- prox_L1(Lambda - t*grad_Lambda, t*lambda)
      if (if_fixB==1){
        B_new<-B0  
      }else{
        B_new <-B- t*grad_B
        B_new <-rho(B_new)
      }                       # GP algorithm
      Psi_new<- Psi- t*grad_Psi
      Psi_new[Psi_new > 3] <- 3
      Psi_new[Psi_new < -3]<- -3
      Sigma_new <-Lambda_new%*%t(B_new)%*%B_new%*%t(Lambda_new)+diag(exp(Psi_new)) 
      
      if(f(Sigma_new, S, Lambda_new,lambda) < f(Sigma, S, Lambda,lambda)) {
        break}
      t <- beta0*t
    }
    
    
    if(i > 1 && (max(abs(Lambda- Lambda_new)) < stop_ep)
       && (max(abs(B-B_new)) < stop_ep)&& (max(abs(Psi- Psi_new)) < stop_ep)) {
      #cat(i, fill = TRUE)
      break }
    #time_iter[i]<-proc.time()[3]-start_time[3]
    #objective_iter[i]<- f(Sigma_new, S, Lambda_new,lambda)
    # if(i > 1 && (abs(objective_iter[i] - objective_iter[i-1]) < 10^-6*5)) {
    #   break
    #   
    # }
    
    # store temp Lambda value
    Lambda <- Lambda_new
    B <- B_new
    Psi <- Psi_new
    Sigma <-Sigma_new
  }
  time=proc.time()[3]-start_time[3]
  #print(paste(paste('Lambda',lambda),paste('iter=',i)))
  return(list(L=Lambda,B=B,Psi=Psi,it=i,t=time))
}