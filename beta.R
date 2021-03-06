library(doParallel)
library(foreach)
library(RCmodels)
library(ggplot2)
qvdata=read.table("V316.txt",skip=3,sep="|",dec=",")
qvdata=qvdata[,c(2:4,7)]
qvdata[,3:4]=qvdata[,4:3]
names(qvdata)=c("Date","Time","W","Q")
qvdata$Time=as.character(qvdata$Time)
qvdata$Date=as.Date(gsub("\\.","-",qvdata$Date),"%d-%m-%Y")
qvdata=qvdata[with(qvdata,order(W)),]
wq=as.matrix(qvdata[,3:4])
qvdata$W=0.01*qvdata$W

Nit=20000

RC=priors("Iceland")

RC$nugget=10^-8
RC$mu_sb=0.5
RC$mu_pb=0.5
RC$tau_pb2=0.25^2
RC$s=3
RC$v=5

RC$y=rbind(as.matrix(log(wq[,2])),0)
RC$w=as.matrix(0.01*wq[,1])
RC$w_tild=RC$w-min(RC$w)

Adist1 <- Adist(RC$w)
RC$A=Adist1$A
RC$dist=Adist1$dist
RC$n=Adist1$n
RC$N=Adist1$N
RC$O=Adist1$O

RC$P=diag(nrow=5,ncol=5,6)-matrix(nrow=5,ncol=5,1)
RC$Sig_ab= rbind(c(RC$sig_a^2, RC$p_ab*RC$sig_a*RC$sig_b), c(RC$p_ab*RC$sig_a*RC$sig_b, RC$sig_b^2))
RC$mu_x=as.matrix(c(RC$mu_a,RC$mu_b, rep(0,RC$n))) #Setja i RC

RC$B=B_splines(t(RC$w_tild)/RC$w_tild[length(RC$w_tild)])
RC$Z=cbind(t(rep(0,2)),t(rep(1,RC$n)))

RC$m1=matrix(0,nrow=2,ncol=RC$n)
RC$m2=matrix(0,nrow=RC$n,ncol=2)
theta.init=rep(0,9)

Dens = function(th) {-Densevalm22(th,RC)$p}
Densmin=optim(par=theta.init,Dens,method="L-BFGS-B",hessian=TRUE)
t_m =Densmin$par
H=Densmin$hessian
phi_b=t_m[3]
sig_b2=t_m[2]
zeta=t_m[1]
lambda=t_m[4:9]
l=log(RC$w_tild+exp(t_m[1])) #as.matrix
varr_m=exp(RC$B%*%lambda)
Sig_eps=diag(as.numeric(rbind(varr_m,0)))
R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(1,RC$n,RC$n)*RC$nugget
Sig_x=rbind(cbind(RC$Sig_ab,matrix(0,nrow=2,ncol=RC$n)),cbind(matrix(0,nrow=RC$n,ncol=2),exp(sig_b2)*R_Beta))

X=rbind(cbind(matrix(1,dim(l)),l,diag(as.numeric(l))%*%RC$A),RC$Z)


L=t(chol(as.matrix(X%*%Sig_x%*%t(X)+Sig_eps)))

w=solve(L,(-RC$y+X%*%RC$mu_x))
mu=RC$mu_x-Sig_x%*%(t(X)%*%(solve(t(L),w)))
LH=t(chol(H))/0.8

cl <- makeCluster(4)
registerDoParallel(cl)
WFill=W_unobserved(RC$O,min=ceiling((min(RC$O)-exp(t_m[1]))*10)/10,max=ceiling(max(RC$O)*10)/10)
RC$W_u=WFill$W_u
RC$W_u_tild=WFill$W_u_tild
RC$Bsim=B_splines(t(RC$W_u_tild)/RC$W_u_tild[length(RC$W_u_tild)])

ptm <- proc.time()
MCMC <- foreach(i=1:4,.combine=cbind,.export=c("Densevalm22","Densevalm22_u")) %dopar% {
    ypo_obs=matrix(0,nrow=RC$N,ncol=Nit)
    param=matrix(0,nrow=9+RC$n+2,ncol=Nit)
    t_old=as.matrix(t_m)
    Dens<-Densevalm22(t_old,RC)
    p_old=Dens$p
    ypo_old=Dens$ypo
    x_old=Dens$x
    
    for(j in 1:Nit){
        t_new=t_old+solve(t(LH),rnorm(9,0,1))
        Densnew<-Densevalm22(t_new,RC)
        x_new=Densnew$x
        ypo_new=Densnew$ypo
        p_new=Densnew$p
        logR=p_new-p_old
        
        if (logR>log(runif(1))){
            t_old=t_new
            p_old=p_new
            ypo_old=ypo_new
            x_old=x_new
            
            
        }
        ypo_obs[,j]=ypo_old
        param[,j]=rbind(t_old,x_old)    
    }
    seq=seq(2000,Nit,5)
    ypo_obs=ypo_obs[,seq]
    param=param[,seq]
    unobserved=apply(param,2,FUN=function(x) Densevalm22_u(x,RC))
    #x_obs=param[10:nrow(param),]
    output=rbind(ypo_obs,unobserved)
    
    return(output)
}
proc.time() - ptm
stopCluster(cl)

ptm <- proc.time()
betasamples=apply(MCMC[(RC$N+length(RC$W_u)+1):nrow(MCMC),],2,FUN=function(x){x[2]+x[3:length(x)]})
yposamples=MCMC[1:(RC$N+length(RC$W_u)),]
ypodata=as.data.frame(t(apply(yposamples,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
names(ypodata)=c("lower","fit","upper")
ypodata$W=c(RC$w,RC$W_u)
ypodata=ypodata[with(ypodata,order(W)),]
betadata=as.data.frame(t(apply(betasamples,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
names(betadata)=c("lower","fit","upper")
betadata$W=c(RC$O,RC$W_u)
betadata=betadata[with(betadata,order(W)),]



smoothbeta=ggplot(data=betadata)+geom_line(aes(W,fit))+geom_line(aes(W,lower),linetype="dashed")+geom_line(aes(W,upper),linetype="dashed")

rcrealsmooth=ggplot(data=ypodata)+geom_line(aes(exp(fit),W))+geom_line(aes(exp(lower),W),linetype="dashed")+geom_line(aes(exp(upper),W),linetype="dashed")+geom_point(data=qvdata,aes(Q,W))
xout=seq(ceiling((min(RC$O)-exp(t_m[1]))*10)/10,ceiling(max(RC$O)*10)/10-0.01,by=0.01)
interpol=approx(ypodata$W,ypodata$fit,xout=xout)
if(length(interpol$x)%%10==0) {
    table=t(as.data.frame(split(x=interpol$y, f=ceiling(seq_along(interpol$y)/10))))
    rownames(table)=seq(min(interpol$x),max(interpol$x),by=0.1)*100
    colnames(table)=0:9
    table=exp(table)
}else  {
    stop("sequence has to be of length that adds up to 10")
}
proc.time()-ptm
