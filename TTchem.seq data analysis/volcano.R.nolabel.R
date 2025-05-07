argv <- commandArgs(TRUE)
library(ggrepel)
library(ggpubr)
library(ggplot2)
library(stringr)
library(gridExtra)
library(ggforce)
library(dplyr)
deg=read.table(argv[1],row.names=1,header=T)
xlab=colnames(deg)[length(colnames(deg))-4]

x=deg[length(colnames(deg))-4]
n=str_replace(rownames(x),"Mus_musculus_tRNA-","")
n=str_replace(rownames(x),"Homo_sapiens_tRNA-","")
y=-log10(deg$padj)
y[y>5]<-5
x[x< -3] <- -3
x[x> 3] <- 3
dat=data.frame(x,y)
colnames(dat)=c("x","y")
rownames(dat)=n

significant_dn=subset(dat, x> -3 &x <= -1  & ( y >= -log10(0.05) & y<5))
significant_up=subset(dat, x<3 &x >= 1  & ( y >= -log10(0.05) & y<5))

nochange=subset(dat,(x >-3 & x<3 & y < -log10(0.05)) | (x > -1 & x< 1))

high_sig_dn=subset(dat, x <= -1 & y >=5)
high_sig_up=subset(dat, x >= 1  & y >=5)
high_nosig=subset(dat, (x > -1 & x < 1 ) & y>=5)
x_extreme_up_sig=subset(dat, x==3 & y >=-log10(0.05))
x_extreme_dn_sig=subset(dat, x==-3 & y >=-log10(0.05))
x_extreme_up_nosig=subset(dat, x==3 & y < -log10(0.05))
x_extreme_dn_nosig=subset(dat, x==-3 & y < -log10(0.05))

significant_all=rbind(significant_dn,significant_up,high_sig_dn,high_sig_up,x_extreme_up_sig,x_extreme_dn_sig)

up=length(subset(dat,x>=1 & y>=-log10(0.05))$x)
up_2=length(subset(dat,x<1 & x>=0 & y>=-log10(0.05))$x)
dn=length(subset(dat,x<=-1 & y>=-log10(0.05))$x)
dn_2=length(subset(dat,x>-1 & x<0 & y>=-log10(0.05))$x)


volcano_plot = ggplot(nochange,aes(x,y))+geom_point(col=rgb(0.5,0.5,0.5,1))

if (length(significant_up$x)>0){
        volcano_plot=volcano_plot+geom_point(data=significant_up,mapping=aes(x=x,y=y),col="orange")
}
if (length(significant_dn$x)>0){
        volcano_plot=volcano_plot+geom_point(data=significant_dn,mapping=aes(x=x,y=y),col="#00BFC4")
}
if (length(x_extreme_dn_nosig$x)>0){
        volcano_plot=volcano_plot+geom_point(data=x_extreme_dn_nosig,mapping=aes(x=x,y=y),shape=25,colour=rgb(0.5,0.5,0.5,1),fill=rgb(0.5,0.5,0.5,1))
}
if (length(x_extreme_up_nosig$x)>0){
        volcano_plot=volcano_plot+geom_point(data=x_extreme_up_nosig,mapping=aes(x=x,y=y),shape=25,colour=rgb(0.5,0.5,0.5,1),fill=rgb(0.5,0.5,0.5,1))
}
if (length(high_nosig$x)>0){
        volcano_plot=volcano_plot+geom_point(data=high_nosig,mapping=aes(x=x,y=y),shape=25,colour=rgb(0.5,0.5,0.5,1),fill=rgb(0.5,0.5,0.5,1))
}
if (length(x_extreme_dn_sig$x)>0){
        volcano_plot=volcano_plot+geom_point(data=x_extreme_dn_sig,mapping=aes(x=x,y=y),shape=25,colour="#00BFC4",fill="#00BFC4")
}
if (length(x_extreme_up_sig$x)>0){
        volcano_plot=volcano_plot+geom_point(data=x_extreme_up_sig,mapping=aes(x=x,y=y),shape=25,colour="orange",fill="orange")
}
if (length(high_sig_up$x)>0){
        volcano_plot=volcano_plot+geom_point(data=high_sig_up,mapping=aes(x=x,y=y),shape=25,colour="orange",fill="orange")
}
if (length(high_sig_dn$x)>0){
        volcano_plot=volcano_plot+geom_point(data=high_sig_dn,mapping=aes(x=x,y=y),shape=25,colour="#00BFC4",fill="#00BFC4")
}

volcano_plot=volcano_plot+
        xlab(xlab)+ylab("-log10(padj)")+
        geom_segment(aes(x = 1, xend=1, y= -log10(0.05), yend=5.5), lty=4)+
        geom_segment(aes(x = 0, xend=0, y= -log10(0.05), yend=5.5), lty=4)+
        geom_segment(aes(x = -1, xend=-1, y= -log10(0.05), yend=5.5), lty=4)+
        geom_hline(yintercept=-log10(0.05),lty=4)+
#	geom_text_repel(data=significant_all,label=rownames(significant_all),vjust=-0.5,hjust=0.5,col="black")+
        annotate("text",label=paste("n=",dn,sep=""),x=-2,y=5.3)+
        annotate("text",label=paste("n=",dn_2,sep=""),x=-0.5,y=5.3)+
        annotate("text",label=paste("n=",up,sep=""),x=2,y=5.3)+
        annotate("text",label=paste("n=",up_2,sep=""),x=0.5,y=5.3)+
        scale_x_continuous(breaks=seq(-10,10,by=1))+
        ylim(0,5.5)+
        scale_y_continuous(breaks=seq(0,6,by=1))+
        theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank())

f=str_replace(argv[1],".xls","")
pdf(paste(f,"nolabel.pdf",sep="."))
volcano_plot
dev.off()
png(paste(f,"nolabel.png",sep="."),height=4,width=4,units="in",res=200)
volcano_plot
dev.off()

