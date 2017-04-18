Args	=	commandArgs(TRUE);
inputF1	=	Args[1];#input file; header, sample names;
output1	=	Args[2];#output pdf file name
library("ggplot2");
library("scales");
library("reshape");
mx	=	read.table(inputF1,sep="\t",head=T,quote="",check.names=F);
mx[,1]	=	as.factor(mx[,1]);
out	=	melt(mx);
colnames(out)   =       c("X","feat","value");
out$X	=	factor(out$X);
out$feat=	factor(out$feat);

hah=ggplot(out,aes(x=as.factor(X),y=value,fill=feat))+
geom_bar(position=position_dodge(),stat="identity",colour="white")+
#coord_cartesian(ylim=c(0.0,1.0))+
#scale_fill_manual(breaks=Feat,values=my_col,labels=Feat,guide = guide_legend(nrow=1,keywidth=0.5))+
theme_bw() +
theme(axis.title.x = element_text(size=0),axis.text.x  = element_text(size=20),axis.title.y = element_text(size=0),axis.text.y  = element_text(size=20))+
theme(legend.title=element_blank(),legend.position="bottom")+
theme(legend.text = element_text(size =18));
pdf(output1,width=10,height=8)
print(hah);
dev.off()

