

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

X=rbind(cbind(matrix(1,dim(l)),l,diag(as.numeric(l))%*%RC$A,RC$Z))


L=t(chol(as.matrix(X%*%Sig_x%*%t(X)+Sig_eps)))

w=solve(L,(-RC$y+X%*%RC$mu_x))
mu=RC$mu_x-Sig_x%*%(t(X)%*%(solve(t(L),w)))

ymu=X%*%mu
ymu=ymu[1:RC$N]
plot(RC$w,exp(ymu))


v=seq(min(wq[,1]),max(wq[,1]),length.out=40)
seq=unique(c(wq[,1],v))
seq=sort(seq)
seq2=c(seq[2:length(seq)],1000)
dist=abs(seq-seq2)
mindist=0.1
Wsim=0.01*seq[which(dist>mindist)]
RC$w_tildsim=as.matrix(Wsim-min(Wsim))
RC$Bsim=B_splines(t(RC$w_tildsim)/RC$w_tildsim[length(RC$w_tildsim)])

cl <- makeCluster(4)
# Register cluster
registerDoParallel(cl)
#Find out how many
MCMC <- foreach(i=1:4, .combine=cbind,.export=c("Densevalm22")) %dopar% {
    beta_u=matrix(0,nrow=length(Wsim),ncol=Nit)
    
    t_old=as.matrix(t_m)
    Dens<-Densevalmbeta(t_old,RC,Wsim)
    p_old=Dens$p
    ypo_old=Dens$ypo
    
    for(j in 1:Nit){
        t_new=t_old+solve(t(LH),rnorm(9,0,1))
        Densnew<-Densevalm22(t_new,RC)
        ypo_new=Densnew$ypo
        p_new=Densnew$p
        logR=p_new-p_old
        
        if (logR>log(runif(1))){
            t_old=t_new
            p_old=p_new
            ypo_old=ypo_new
            
        }
        ypo[,j]=rbind(ypo_old,t_old)
    }
    
    seq=seq(2000,Nit,5)
    ypo=ypo[,seq]
    
    return(ypo)
}
quantmatrix=head(MCMC,nrow(MCMC)-9)
t=tail(MCMC,9)


