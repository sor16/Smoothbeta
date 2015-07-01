
qvdata=read.table("V508.txt",skip=3,sep="|",dec=",")
qvdata=qvdata[,c(2:4,7)]
qvdata[,3:4]=qvdata[,4:3]
names(qvdata)=c("Date","Time","W","Q")
qvdata$Time=as.character(qvdata$Time)
qvdata$Date=as.Date(gsub("\\.","-",qvdata$Date),"%d-%m-%Y")
qvdata=qvdata[with(qvdata,order(W)),]
wq=as.matrix(qvdata[,3:4])


v=seq(min(wq[,1]),max(wq[,1]),length.out=100)
seq=unique(c(wq[,1],v))
seq=sort(seq)
seq2=c(seq[2:length(seq)],1000)
dist=abs(seq-seq2)
mindist=0.1
newseq=seq[which(dist>mindist)]