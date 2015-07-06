library(doParallel)
library(foreach)
library(RCmodels)
qvdata=read.table("V508.txt",skip=3,sep="|",dec=",")
qvdata=qvdata[,c(2:4,7)]
qvdata[,3:4]=qvdata[,4:3]
names(qvdata)=c("Date","Time","W","Q")
qvdata$Time=as.character(qvdata$Time)
qvdata$Date=as.Date(gsub("\\.","-",qvdata$Date),"%d-%m-%Y")
qvdata=qvdata[with(qvdata,order(W)),]
wq=as.matrix(qvdata[,3:4])

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

Adist1 <- Adist(sort(RC$w))
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
Densmin=optim(par=theta.init,Dens,method="BFGS",hessian=TRUE)
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
# Register cluster
registerDoParallel(cl)
#Find out how many
ptm <- proc.time()
MCMC <- foreach(i=1:4,.combine=cbind,.export=c("Densevalm22")) %dopar% {
    output=matrix(0,nrow=nrow(wq)+9,ncol=Nit)
    
    t_old=as.matrix(t_m)
    Dens<-Densevalm22(t_old,RC)
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
        output[,j]=rbind(ypo_old,t_old)
    }
    
    seq=seq(2000,Nit,5)
    output=output[,seq]
    
    return(output)
}
quantmatrix=head(MCMC,nrow(MCMC)-9)
t=tail(MCMC,9)
proc.time() - ptm
##########################
#BETA_U
##########################
size_grid=50
v=seq(min(wq[,1]),max(wq[,1]),length.out=size_grid)
#seq=unique(c(wq[,1],v))
Wsim=0.01*sort(v)
#seq2=c(seq[2:length(seq)],1000)
#filter=abs(seq-seq2)
#mindist=0.1
#Wsim=0.01*seq[which(filter>mindist)]
#RC$w_tildsim=as.matrix(Wsim-min(Wsim))
# v=seq(min(wq[,1]),max(wq[,1]),length.out=size_grid)
# seq=unique(c(wq[,1],v))
# seq=sort(seq)
# seq2=c(seq[2:length(seq)],1000)
# filter=abs(seq-seq2)
# mindist=0.1
# Wsim=0.01*seq[which(filter>mindist)]
RC$w_tildsim=as.matrix(Wsim-min(Wsim))
RC$Bsim=B_splines(t(RC$w_tildsim)/RC$w_tildsim[length(RC$w_tildsim)])

#MCMC <- foreach(i=1:4, .combine=cbind,.export=c("Densevalmbeta")) %dopar% {
    beta_u=apply(t,2,FUN=function(x) Densevalmbeta(x,RC,Wsim=Wsim))
    seq=seq(2000,Nit,5)
    beta_u=beta_u[,seq]
#}
beta_u_median=t(apply(beta_u,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T))
X=rbind(rep(1,m),l,matrix(0,m,n),diag(l))
x=rbind(mu[1],mu[2])

stopCluster(cl)

# all_values <- function(x) {
#     if(is.null(x)) return(NULL)
#     row <- qvdata[data$id == x$id, ]
#     paste0(names(row), ": ", format(row), collapse = "<br />")
# }




# base <- data %>% ggvis(x = ~exp(Q), y = ~W,key:= ~id) %>%
#     layer_points() %>% add_tooltip(all_values, "hover")%>%layer_lines(x= ~exp(fit),y= ~W) %>%
#     layer_lines(x= ~exp(lower),y= ~W,strokeDash:=6)%>%layer_lines(x= ~exp(upper),y= ~W,strokeDash:=6)
# base
# 
data=as.data.frame(t(apply(ypo,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
names(data)=c("lower","fit","upper")
CI=ggplot(data=data,aes(W,fit))+geom_line()+geom_line(aes(W,lower),linetype="dashed")+geom_line(aes(W,upper),linetype="dashed")
fit=ggplot(data=data,aes(W,fit))+geom_line()