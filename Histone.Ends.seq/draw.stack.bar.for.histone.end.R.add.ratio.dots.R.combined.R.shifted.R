library(stringr)
library(gridExtra)
library(ggplot2)
library(ggpubr)
argv=commandArgs(TRUE)


a_0h=read.table(argv[1],header=T)
a_3h=read.table(argv[2],header=T)
a_24h=read.table(argv[3],header=T)


a_0h$time="0h"
a_3h$time="3h"
a_24h$time="24h"

a_3h$new_pos=a_3h$pos-3
for (i in rownames(a_3h)){if(a_3h[i,]$nt != "3+"){ a_3h[i,]$new_pos =a_3h[i,]$pos-as.numeric(a_3h[i,]$nt)}}
a_24h$new_pos=a_24h$pos-3
for (i in rownames(a_24h)){if(a_24h[i,]$nt != "3+"){ a_24h[i,]$new_pos =a_24h[i,]$pos-as.numeric(a_24h[i,]$nt)}}
a_0h$new_pos=a_0h$pos-3
for (i in rownames(a_0h)){if(a_0h[i,]$nt != "3+"){ a_0h[i,]$new_pos =a_0h[i,]$pos-as.numeric(a_0h[i,]$nt)}}

draw_stack_barplot=function(i){
	df1=subset(a_0h,a_0h$gene == i)
	df2=subset(a_3h,a_3h$gene == i)
	df3=subset(a_24h,a_24h$gene == i)
	df=rbind(df1,df2,df3)
	df$time=factor(df$time,levels=c("0h","3h","24h"))
	df$nt=paste(df$nt," nt tails",sep="")
	df$nt=factor(df$nt,levels=c("3+ nt tails","2 nt tails","1 nt tails","0 nt tails"))
	df$new_pos=max(df$new_pos) - df$new_pos - 5
	p<-ggplot(df,mapping=aes(fill=nt,y=CPM,x=new_pos))+
		geom_bar(position="stack", stat="identity")+
		scale_x_reverse("Nucleotides from 3' end",limits=c(10,0),breaks=seq(10,0,by=-2))+
		ylab("CPM")+ggtitle(unique(df$gene))+
		facet_wrap(~time,scales = "free_x", strip.position = "bottom")+
		labs(fill="")+
		scale_fill_manual(values=c("#ee2025", "#f5a713", "#6dbb4d","#3953a7"))+
		theme(plot.title = element_text(hjust = 0.5),
			panel.grid.major = element_blank(),
#			legend.position="none",
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.text = element_text(colour = "black"),
			axis.ticks = element_line(colour = "black"),
			strip.background = element_rect(colour="black",fill="white"),
			strip.placement = "outside")
	return(p)
}
if (grepl("H3",argv[1])){
	glist=lapply(unique(a_0h$gene),draw_stack_barplot)
	ggsave(argv[4],width=20,height=10,units="in",marrangeGrob(glist,nrow=3,ncol=3,top=NULL))
}else{
	glist=lapply(unique(a_0h$gene),draw_stack_barplot)
	ggsave(argv[4],width=20,height=10,units="in",marrangeGrob(glist,nrow=4,ncol=4,top=NULL))
}
dev.off()


