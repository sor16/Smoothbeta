qvdata=read.table("V508.txt",skip=3,sep="|",dec=",")
qvdata=qvdata[,c(2:4,7)]
qvdata[,3:4]=qvdata[,4:3]
names(qvdata)=c("Date","Time","W","Q")
qvdata$Time=as.character(qvdata$Time)
qvdata$Date=as.Date(gsub("\\.","-",qvdata$Date),"%d-%m-%Y")
qvdata=qvdata[with(qvdata,order(W)),]
wq=as.matrix(qvdata[,3:4])
# dist fylki fyrir joint vigur
Densevalmbeta <- function(th,RC,Wsim){ 
zeta=th[1]
phi_b=th[2]
sig_b2=th[3]
lambda = th[4:9]
varr_sim = c(exp(RC$Bsim %*% lambda))

c=min(Wsim)-exp(zeta)
l=log(Wsim-c)
m=length(Wsim)
n=RC$n
#Wsim_rep=as.matrix(c(unique(RC$w),Wsim))%*%rep(1,m+n)
#dist=abs(t(Wsim_rep)-Wsim_rep)
W_all=c(RC$O,Wsim)
dist=abs(outer(W_all,W_all,FUN="-"))
sigma_all=sig_b2*(1+(sqrt(5)/exp(phi_b))*dist+(5/(3*exp(phi_b)))*dist^2)*exp(-(sqrt(5)/exp(phi_b))*dist) + diag(n+m) * RC$nugget
sigma_11=sigma_all[1:n,1:n]
#sigma_trick=sigma_11 +sqrt(varr_sim)*diag(rnorm(n))
sigma_22=sigma_all[(n+1):(m+n),(n+1):(m+n)]
sigma_12=sigma_all[1:n,(n+1):(n+m)]
sigma_21=sigma_all[(n+1):(n+m),1:n]
# Wsim_rep=as.matrix(Wsim)%*%rep(1,m)
# dist2=abs(t(Wsim_rep)-Wsim_rep)
# sigma_11a=sig_b2*(1+(sqrt(5)/exp(phi_b))*RC$dist+(5/(3*exp(phi_b)))*RC$dist^2)*exp(-(sqrt(5)/exp(phi_b))*RC$dist)
# sigma_22a=sig_b2*(1+(sqrt(5)/exp(phi_b))*dist2+(5/(3*exp(phi_b)))*dist2^2)*exp(-(sqrt(5)/exp(phi_b))*dist2)
# dist3=sapply(unique(RC$w),function(x) abs(x-Wsim))
# sigma_21a=sig_b2*(1+(sqrt(5)/exp(phi_b))*dist3+(5/(3*exp(phi_b)))*dist3^2)*exp(-(sqrt(5)/exp(phi_b))*dist3)
# sigma_12a=t(sigma_21a)



mu_beta=sigma_21%*%solve(sigma_11,mu[3:length(mu)])
Sigma=(sigma_22-sigma_21%*%solve(sigma_11,sigma_12))
ncols <- ncol(Sigma)
rbeta=as.numeric(mu_beta) + rnorm(ncols) %*% chol(Sigma)
return(rbeta)
}

rbeta=mvrnorm(n=1,mu=mu_beta,Sigma=Sigma)
rbeta=rmvn(n=1,mu=mu_beta,Sigma=Sigma)
# 
# #
# Sig_eps = diag(c(varr, 0))
# R_Beta = (1 + sqrt(5) * RC$dist/exp(phi_b) + 5 * RC$dist^2/(3 *exp(phi_b)^2)) * 
#     exp(-sqrt(5) * RC$dist/exp(phi_b)) +diag(RC$n) * RC$nugget
# 
# Sig_x = rbind(cbind(RC$Sig_ab, RC$m1), cbind(RC$m2, exp(sig_b2) * R_Beta))
# 
# L = t(chol(X %*% Sig_x %*% t(X) + Sig_eps))
# 
# 
# W = solve(L, X %*% Sig_x)
# x_u = RC$mu_x + t(chol(Sig_x)) %*% rnorm(RC$n + 2)
# sss = (X %*% x_u) - RC$y + rbind(sqrt(varr) * as.matrix(rnorm(RC$N)), 0)
# x = as.matrix(x_u - t(W) %*% solve(L, sss))
# yp = X %*% x
# yp = yp[1:RC$N, ]
# ypo = yp + as.matrix(rnorm(RC$N)) * sqrt(varr)
return(rbeta)
}







