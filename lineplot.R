Args    =       commandArgs(TRUE);
inputF1 =       Args[1];#input file; header, sample names;
output1 =       Args[2];#output pdf file name

a=read.table(inputF1,sep="\t",header=FALSE)
##class(a[1]) ##should not be data frame type
png(output1)
plot(a[,1],a[,2],type="l",col=2,xlab="read position",ylab="% of total(per read position)",lty=1,lwd=3,ylim=c(0,0.5),main="Nucletide Distribution")
for(i in 3:6){
    lines(a[,1],a[,i],type="l",col=i,lwd=3)
}
legend("topright",c("A","C","G","T","N"),col=2:6,lty=1,box.col="white",lwd=3)
dev.off()
