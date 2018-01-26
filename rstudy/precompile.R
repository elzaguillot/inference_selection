## script to expolore some outputs created by the vscompute programm
## Left here to give an idea of what is possible

setwd('vitalit')

dataMoran=read.table("s_values_250.txt")
dataWright=read.table("s_valuesWF_250.txt")

names(dataMoran)[1]="popsize1"
names(dataMoran)[2]="popsize2"
names(dataMoran)[3]="s"
names(dataMoran)[4]="nbgen"

names(dataWright)[1]="popsize1"
names(dataWright)[2]="popsize2"
names(dataWright)[3]="s"
names(dataWright)[4]="nbgen"

plot(unlist(dataMoran[1,5:dim(dataMoran)[2]]),ylim=c(0.01,1))
for(i in 1:dim(dataMoran)[1])
    lines(unlist(dataMoran[i,5:dim(dataMoran)[2]]),col=rainbow(10)[dataMoran$s[i]*10],lwd=3)
for(i in 1:dim(dataWright)[1])
    lines(unlist(dataWright[i,5:dim(dataMoran)[2]]),col=rainbow(10)[dataWright$s[i]*10],lwd=3,lty=2)


# Ancient allele frequency
plot(dataWright[,3],dataWright[,5],col="red",pch="+",ylim=c(0.5,1))
for(i in c(seq(100,1000,by=100),seq(1000,10000,by=1000)))
{
    name1=paste("s_values_",i,".txt",sep="",collapse="")
    name2=paste("s_valuesWF_",i,".txt",sep="",collapse="")
    dataMoran=read.table(name1)
    dataWright=read.table(name2)
    points(dataWright[,3],dataWright[,5],col="red",pch="+",ylim=c(0.6,1))
    points(dataMoran[,3],dataMoran[,5],col="blue", pch="o")
}

# singleton freq
plot(dataWright[,3],dataWright[,6],col="red",pch="+",ylim=c(0,0.0001))
for(i in seq(50,1050,by=50))
{
    name1=paste("s_values_",i,".txt",sep="",collapse="")
    name2=paste("s_valuesWF_",i,".txt",sep="",collapse="")
    dataMoran=read.table(name1)
    dataWright=read.table(name2)
    points(dataWright[,3],dataWright[,6],col="red",pch="+",ylim=c(0.6,1))
    points(dataMoran[,3],dataMoran[,6],col="blue", pch="o")
}
