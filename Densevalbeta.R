RC$null1=matrix(0,nrow=length(Wsim)+RC$n,ncol=2)
RC$null2=t(RC$null1)
#need median value from the former MCMC for mu
RC$mu=mu
# dist fylki fyrir joint vigur
Densevalmbeta <- function(param,RC,Wsim){
x=param[1:(RC$n+2)]
th=param[(RC$n+3):length(param)]
zeta=th[1]
phi_b=exp(th[3])
sig_b2=exp(th[2])
lambda = th[4:9]
varr = c(exp(RC$Bsim %*% lambda))

c=min(Wsim)-exp(zeta)
l=log(Wsim-c)
m=length(Wsim)
n=RC$n
#Wsim_rep=as.matrix(c(unique(RC$w),Wsim))%*%rep(1,m+n)
#dist=abs(t(Wsim_rep)-Wsim_rep)
W_all=c(RC$O,Wsim)
dist=abs(outer(W_all,W_all,FUN="-"))
sigma_all=sig_b2*(1 + sqrt(5)*dist/phi_b+(5*dist^2)/(3*phi_b^2))*exp(-sqrt(5)*dist/phi_b) + diag(length(W_all))*RC$nugget
sigma_11=sigma_all[1:n,1:n]
sigma_22=sigma_all[(n+1):(m+n),(n+1):(m+n)]
sigma_12=sigma_all[1:n,(n+1):(n+m)]
sigma_21=sigma_all[(n+1):(n+m),1:n]

mu_beta=sigma_21%*%solve(sigma_11,x[3:length(x)])
Sigma=(sigma_22-sigma_21%*%solve(sigma_11,sigma_12))
ncols <- ncol(Sigma)
beta_u=as.numeric(mu_beta) + rnorm(ncols) %*% chol(Sigma)
#rbeta=mvrnorm(n=1,mu=mu_beta,Sigma=Sigma)
#beta_u=rmvn(n=1,mu=mu_beta,Sigma=Sigma)
# sigma_eps = diag(varr)
X=cbind(rep(1,m),l,matrix(0,m,n),diag(l))
x=c(x,beta_u)
# #Adding linear constraint
# sigma_all=rbind(cbind(RC$Sig_ab,RC$null2),cbind(RC$null1,sigma_all))
# L = t(chol(X %*% sigma_all %*% t(X) + sigma_eps))
# W = solve(L, X %*% sigma_all)
# x_u=c(RC$mu,beta_u)
#sss = (X %*% x_u) - RC$y + rbind(sqrt(varr) * as.matrix(rnorm(RC$N)), 0)
#x = as.matrix(x_u - t(W) %*% solve(L, sss))
#yp = X %*% x
#yp = yp[1:RC$N, ]
ypo = X%*%x + as.matrix(rnorm(m)) * sqrt(varr)
return(ypo)
}







