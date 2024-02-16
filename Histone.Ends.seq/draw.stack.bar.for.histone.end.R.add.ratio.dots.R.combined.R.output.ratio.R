library(stringr)
library(gridExtra)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggpubr)
argv=commandArgs(TRUE)


a_0h=read.table(argv[1],header=T)
a_3h=read.table(argv[2],header=T)
a_24h=read.table(argv[3],header=T)
a_0h_r=read.table(argv[5],header=T)
a_3h_r=read.table(argv[6],header=T)
a_24h_r=read.table(argv[7],header=T)



a_0h$time="0h"
a_3h$time="3h"
a_24h$time="24h"
a_0h_r$time="0h"
a_3h_r$time="3h"
a_24h_r$time="24h"
glist=list()
j=0
output_df=data.frame(matrix(ncol = 3, nrow = 0))
colnames(output_df)=c("0h","3h","24h")
output_df_ratio=output_df
for (i in unique(a_0h$gene)){
	j =j+1
	df1=subset(a_0h,a_0h$gene == i)
	df2=subset(a_3h,a_3h$gene == i)
	df3=subset(a_24h,a_24h$gene == i)
	df4=subset(a_0h_r,a_0h_r$gene == i)
	df5=subset(a_3h_r,a_3h_r$gene == i)
	df6=subset(a_24h_r,a_24h_r$gene == i)

	df4=df4 %>% top_n (1,total_CPM)
	df5=df5 %>% top_n (1,total_CPM)
	df6=df6 %>% top_n (1,total_CPM)
	df=rbind(df1,df2,df3)
	df_r=rbind(df4,df5,df6)
	df$time=factor(df$time,levels=c("0h","3h","24h"))
	df_r$time=factor(df_r$time,levels=c("0h","3h","24h"))
	df_r$nt="0"
	start=min(df$pos)
	end=max(df$pos)
	coef=max(df_r$total_CPM)/max(df_r$non_tail_cpm_ratio)
	p<-ggplot(df,mapping=aes(fill=nt,y=CPM,x=pos))+
		geom_bar(position="stack", stat="identity")+
		scale_y_continuous(sec.axis = sec_axis(~./coef,name="Non-templated ratio"))+
		scale_x_continuous("position",limits=c(start,end),breaks=seq(start,end,by=20),labels=seq(start,end,by=20))+
		ylab("CPM")+ggtitle(unique(df$gene))+
		facet_wrap(~time,scales = "free_x", strip.position = "bottom")+
		theme(plot.title = element_text(hjust = 0.5),
			panel.grid.major = element_blank(),
#			legend.position="none",
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			strip.background = element_rect(colour="black",fill="white"),
			strip.placement = "outside")
	p=p+geom_point(df_r,mapping=aes(x=pos,y=non_tail_cpm_ratio*coef),col="black")
	output_df[i,]=c(df4$non_tail_cpm_ratio,df5$non_tail_cpm_ratio,df6$non_tail_cpm_ratio)
	output_df_ratio[i,]=c(1,df5$non_tail_cpm_ratio/df4$non_tail_cpm_ratio,df6$non_tail_cpm_ratio/df4$non_tail_cpm_ratio)
	glist[[j]]=p	
}
bk=seq(0.7,2,by=0.01)

if (grepl("H3",argv[1])){
	ggsave(argv[4],width=14,height=10,units="in",marrangeGrob(glist,nrow=3,ncol=3,top=NULL))
dev.off()
pdf(argv[8],height=2.5,width=2)
pheatmap(output_df,cluster_rows=F,cluster_cols=F,scale="none",col=colorRampPalette(c('green','white','red'))(100),angle_col=0,main="Non temp ratio",legend_breaks=c(0,0.2,0.4,0.6,0.8,1),legend_labels=c(0,0.2,0.4,0.6,0.8,1))
pheatmap(output_df,cluster_rows=F,cluster_cols=F,scale="none",col=colorRampPalette(c('white','red'))(100),angle_col=0,main="Non temp ratio",legend_breaks=c(0,0.2,0.4,0.6,0.8,1),legend_labels=c(0,0.2,0.4,0.6,0.8,1))
pheatmap(output_df_ratio,cluster_rows=F,cluster_cols=F,scale="none",col=colorRampPalette(c('green','white','red'))(length(bk)),angle_col=0,main="fold to 0h",breaks=bk)
pheatmap(output_df_ratio,cluster_rows=F,cluster_cols=F,scale="none",col=colorRampPalette(c('white','red'))(length(bk)),angle_col=0,main="fold to 0h",breaks=bk)
dev.off()
}else{
	ggsave(argv[4],width=14,height=10,units="in",marrangeGrob(glist,nrow=4,ncol=4,top=NULL))
dev.off()
pdf(argv[8],height=4,width=2)
pheatmap(output_df,cluster_rows=F,cluster_cols=F,scale="none",col=colorRampPalette(c('green','white','red'))(100),angle_col=0,main="Non temp ratio",legend_breaks=c(0,0.2,0.4,0.6,0.8,1),legend_labels=c(0,0.2,0.4,0.6,0.8,1))
pheatmap(output_df,cluster_rows=F,cluster_cols=F,scale="none",col=colorRampPalette(c('white','red'))(100),angle_col=0,main="Non temp ratio",legend_breaks=c(0,0.2,0.4,0.6,0.8,1),legend_labels=c(0,0.2,0.4,0.6,0.8,1))
pheatmap(output_df_ratio,cluster_rows=F,cluster_cols=F,scale="none",col=colorRampPalette(c('green','white','red'))(length(bk)),angle_col=0,main="fold to 0h",breaks=bk)
pheatmap(output_df_ratio,cluster_rows=F,cluster_cols=F,scale="none",col=colorRampPalette(c('white','red'))(length(bk)),angle_col=0,main="fold to 0h",breaks=bk)
dev.off()
}


