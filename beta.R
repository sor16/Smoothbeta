RC=list()

RC$mu_a=3.20
RC$mu_b=2.29
RC$sig_a=sqrt(1.21)
RC$sig_b=sqrt(0.48)
RC$p_ab=-0.61
RC$mu_c=1.9000
RC$nugget=10^-8


RC$mu_sb=0.5
RC$mu_pb=0.5
RC$tau_pb2=0.25^2
RC$s=3
RC$v=5

#wq = as.matrix(read.table('15.txt'))
qvdata=read.table("V508.txt",skip=3,sep="|",dec=",")
qvdata=qvdata[,c(2:4,7)]
qvdata[,3:4]=qvdata[,4:3]
names(qvdata)=c("Date","Time","W","Q")
qvdata$Time=as.character(qvdata$Time)
qvdata$Date=as.Date(gsub("\\.","-",qvdata$Date),"%d-%m-%Y")
qvdata=qvdata[with(qvdata,order(W)),]
wq=as.matrix(qvdata[,3:4])



RC$y=rbind(as.matrix(log(wq[,2])),0)
RC$w=as.matrix(0.01*wq[,1])
RC$w_tild=RC$w-min(RC$w)
# 
H=RC$w
Q=wq[,2]
dat=data.frame(H,Q)

###

Adist1 <- Adist(RC$w)
RC$A=Adist1$A
RC$dist=Adist1$dist
RC$n=Adist1$n
RC$N=Adist1$N

RC$P=diag(nrow=5,ncol=5,6)-matrix(nrow=5,ncol=5,1)

RC$Sig_ab= rbind(c(RC$sig_a^2, RC$p_ab*RC$sig_a*RC$sig_b), c(RC$p_ab*RC$sig_a*RC$sig_b, RC$sig_b^2))

RC$mu_x=as.matrix(c(RC$mu_a,RC$mu_b, rep(0,RC$n))) #Setja i RC

RC$B=B_splines(t(RC$w_tild)/RC$w_tild[length(RC$w_tild)])
RC$Z=cbind(t(rep(0,2)),t(rep(1,RC$n)))

#t_median=t(apply(t, 1, quantile, probs = c(0.5),  na.rm = TRUE))
#t is a 9 * 50000 matrix. MCMC for V508.txt
t_median=c(-0.1391361,-2.275301, 0.4984314, -5.993821, -6.354939, -6.678067, -6.688161,-6.704671, -6.698719)
#t is a 9 * 50000 matrix. MCMC for 15.txt
#t_median=c(-1.630088, -2.313527, 0.4967125, -4.465228, -3.648009, -4.898924, -5.848946,-5.031323, -5.061868)
t_m=t_median
# t_m =Densmin$par
# H=Densmin$hessian



phi_b=t_m[3]
sig_b2=t_m[2]
zeta=t_m[1]
lambda=t_m[4:9]

l=log(RC$w_tild+exp(t_m[1])) #as.matrix

varr_m=exp(RC$B%*%lambda)
Sig_eps=diag(as.numeric(rbind(varr_m,0)))
R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(1,RC$n,RC$n)*RC$nugget
Sig_x=rbind(cbind(RC$Sig_ab,matrix(0,nrow=2,ncol=RC$n)),cbind(matrix(0,nrow=RC$n,ncol=2),exp(sig_b2)*R_Beta))

X=Matrix(rbind(cbind(matrix(1,dim(l)),l,Matrix(diag(as.numeric(l)),sparse=TRUE)%*%RC$A),RC$Z),sparse=TRUE)


L=t(chol(as.matrix(X%*%Sig_x%*%t(X)+Sig_eps)))

w=solve(L,(-RC$y+X%*%RC$mu_x))
mu=RC$mu_x-Sig_x%*%(t(X)%*%(solve(t(L),w)))

ymu=X%*%mu
ymu=ymu[1:RC$N]
plot(RC$w,exp(ymu))