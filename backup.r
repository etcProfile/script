setwd("~/Downloads/达尔文兰形态计量/zfjy")
read.table("zfjy_grow_curve.txt",header = T,sep = "\t") -> Lenght
read.table("zfjy_grow_curve.txt",header = T,sep = "\t") -> Lenght
View(Lenght)
melt(Length,id='day') -> melt_data
library(tidy)
library(tidyverse)
melt(Length,id='day') -> melt_data
library(reshape)
read.table("zfjy_grow_curve.txt",header = T,sep = "\t") -> Lenght
melt(Length,id='Day') -> melt_data
read.table("zfjy_grow_curve.txt",header = T,sep = "\t") -> length
melt(Length,id='Day') -> melt_data
melt(length,id='Day') -> melt_data
View(length)
melt(length,id='X') -> melt_data
wirte.table(melt_data,file="melt.txt",sep="\t")
write.table(melt_data,file="melt.txt",sep="\t")
read.table("melt.txt",header = T) -> melt
View(melt)
logSS <- getInitial(value ~ SSlogis(day, alpha, xmid, scale), data = melt)
K_start <- logSS["alpha"]
R_start <- 1/logSS["scale"]
N0_start <- logSS["alpha"]/(exp(logSS["xmid"]/logSS["scale"])+1)
log_fit <- nls(value~logF(day,Linf,gninf,Ti),data=melt,start=list(Linf=K_start,gninf=R_start,Ti=N0_start))
f.starts <- vbStarts(value~day,data=melt)
vb_fit <- nls(value~vb(day,Linf,K,t0),data=melt,start=f.starts)
log_formula <- formula(n ~ K*N0*exp(R*m)/(K + N0*(exp(R*m) - 1)))
formu<-nls(log_formula, start = list(K = K_start, R = R_start, N0 = N0_start))
log_formula <- formula(value ~ K*N0*exp(R*day)/(K + N0*(exp(R*day) - 1)))
formu<-nls(log_formula, start = list(K = K_start, R = R_start, N0 = N0_start))
logF <- logisticFuns()
log_fit <- nls(value~logF(day,Linf,gninf,Ti),data=melt,start=list(Linf=K_start,gninf=R_start,Ti=N0_start))
log_formula <- formula(value ~ K*N0*exp(day*m)/(K + N0*(exp(day*m) - 1)))
formu<-nls(log_formula, start = list(K = K_start, R = R_start, N0 = N0_start))
formu<-nls(log_formula, data=melt,start = list(K = K_start, R = R_start, N0 = N0_start))
log_formula <- formula(value ~ K*N0*exp(R*day)/(K + N0*(exp(R*day) - 1)))
formu<-nls(log_formula, data=melt,start = list(K = K_start, R = R_start, N0 = N0_start))
f.boot1 <- Boot(formu)
confint(f.boot1)
days <- seq(0,20,by=0.2)
predict2 <- function(x) predict(x,data.frame(day=days))
predict2(formu)
vb_boot2 <- Boot(formu,f=predict2)
vb_preds1 <- data.frame(days,predict(formu,data.frame(day=days)),confint(vb_boot2))
names(vb_preds1) <- c("day","fit","LCI","UCI")
vbFitPlot <- ggplot() +
geom_ribbon(data=log_preds1,aes(x=day,ymin=LCI,ymax=UCI),fill="gray90") +
geom_point(data=melt,aes(y=value,x=day),size=2,alpha=0.1) +
geom_line(data=log_preds1,aes(y=fit,x=day),size=1) +
scale_y_continuous(name="Total Length (mm)",limits=c(0,700),expand=c(0,0)) +
scale_x_continuous(name="Age (years)",expand=c(0,0)) +
theme_bw() +
theme(panel.grid=element_blank())
vbFitPlot <- ggplot() +
geom_ribbon(data=vb_preds1,aes(x=day,ymin=LCI,ymax=UCI),fill="gray90") +
geom_point(data=melt,aes(y=value,x=day),size=2,alpha=0.1) +
geom_line(data=vb_preds1,aes(y=fit,x=day),size=1) +
scale_y_continuous(name="Total Length (mm)",limits=c(0,700),expand=c(0,0)) +
scale_x_continuous(name="Age (years)",expand=c(0,0)) +
theme_bw() +
theme(panel.grid=element_blank())
vbFitPlot
vbFitPlot <- ggplot() +
geom_ribbon(data=vb_preds1,aes(x=day,ymin=LCI,ymax=UCI),fill="gray90") +
geom_point(data=melt,aes(y=value,x=day),size=2,alpha=0.1) +
geom_line(data=vb_preds1,aes(y=fit,x=day),size=1) +
scale_y_continuous(name="Total Length (mm)",limits=c(0,200),expand=c(0,0)) +
scale_x_continuous(name="Age (years)",expand=c(0,0)) +
theme_bw() +
theme(panel.grid=element_blank())
vbFitPlot
vbFitPlot <- ggplot() +
geom_ribbon(data=vb_preds1,aes(x=day,ymin=LCI,ymax=UCI),fill="gray90") +
geom_point(data=melt,aes(y=value,x=day),size=2,alpha=0.1) +
geom_line(data=vb_preds1,aes(y=fit,x=day),size=1) +
scale_y_continuous(name="Total Length (mm)",limits=c(0,200),expand=c(0,0)) +
scale_x_continuous(name="Age (years)",expand=c(0,0)) +
theme_bw()
vbFitPlot
vbFitPlot <- ggplot() +
geom_ribbon(data=vb_preds1,aes(x=day,ymin=LCI,ymax=UCI),fill="gray90") +
geom_point(data=melt,aes(y=value,x=day),size=2,alpha=0.1) +
geom_line(data=vb_preds1,aes(y=fit,x=day),size=1) +
scale_y_continuous(name="Total Length (mm)",limits=c(0,200),expand=c(0,0)) +
scale_x_continuous(name="Days",expand=c(0,0)) +
theme_bw()
vbFitPlot
library(eoffice)
topptx(filename = "zfjy_growth_curve.pptx")
topptx(filename = "zfjy_growth_curve.pptx",width = 6,height = 5)


plot_cpoa <- function(p, d, centroids=NULL, col_var, pch_var, col_comp, pch_comp,
                      shapes, colors, file_name, svg=F, constraint, 
                      ci, var_tbl, cap_var, perm_anova, cex.main=1, cex=0.8) {
  
  col <- rep("black", dim(p)[1])
  pch <- rep(19, dim(p)[1])
  for (i in 1:length(col_var)){
    index <- grep(col_var[i], d[rownames(p), col_comp[i]]) 
    col[index] <- colors[i]
  }
  for (i in 1:length(pch_var)){
    index <- grep(pch_var[i], d[rownames(p), pch_comp[i]]) 
    pch[index] <- shapes[i]
  }
  xlab <- paste("Constrained PCoA 1 (", format(cap_var[1]*100, digits=4), " %)", sep="")
  ylab <- paste("Constrained PCoA 2 (", format(cap_var[2]*100, digits=4), " %)", sep="")
  main <- paste(constraint, ": [", format(var_tbl["constrained", "proportion"]*100, digits=2),
                "% of variance; P < ", format(perm_anova[1,4], digits=2),
                "; 95% CI = ", format(ci[1]*100, digits=2), 
                "%, ", format(ci[2]*100, digits=2), "%]", sep="")	
  if(svg) svg(file=file_name) else pdf(file=file_name)
  plot(rbind(centroids, p), col="white", pch=19, xlab=xlab, ylab=ylab, main=main,font.main=1,font=1, cex.main=cex.main)
  points(p, col=col, cex=cex, pch=pch)
  arrows(x0=c(0,0,0,0), y0=c(0,0,0,0), x1=centroids[,1], y1=centroids[,2], col=colors, pch=17, length=0.1, cex=1.5, lwd=2)
  abline(v=0, h=0, lty="dotted")
  legend("topright",bty="n",title="Group",inset=0.01,legend=col_var,col=colors,pch=c(16,18,),text.font = 3,cex = 0.5)
  dev.off()
  
}


read.table("area.txt",header = T,sep = "\t") -> area
area_mean <- area %>% group_by(group1, group2) %>% summarise(mean_area = mean(area), sd_area = sd(area))
area_mean$group2 <- factor(area_mean$group2,levels = c("1.5cm","4cm","6cm","9cm","12cm","mature"))
area_mean$group1 <- factor(area_mean$group1,levels = c("Base","BMM","Middle","TMM","Tip"))
ggplot(area_mean,aes(x=group1,y=mean_area,fill=group2))+geom_bar(stat="identity",position = 'dodge',width = 0.75)+scale_fill_brewer(name="Spur length",palette = "Blues")+geom_errorbar(aes(ymin = mean_area - sd_area, ymax = mean_area + sd_area), position = position_dodge(0.75),width=0.3,alpha=0.6)+theme_bw()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),panel.grid=element_blank())+ylab("Cell Area")+xlab("Spur of Area")+theme(legend.position = "none")
topptx(filename = "Spur_Area.pptx",width = 7,height = 7)
read.table("height.txt",header = T,sep = "\t") -> height
height_mean <- height %>% group_by(group1, group2) %>% summarise(mean_height = mean(height), sd_height = sd(height))
height_mean$group2 <- factor(height_mean$group2,levels = c("1.5cm","4cm","6cm","9cm","12cm","mature"))
height_mean$group1 <- factor(height_mean$group1,levels = c("Base","BMM","Middle","TMM","Tip"))
ggplot(height_mean,aes(x=group1,y=mean_height,fill=group2))+geom_bar(stat="identity",position = 'dodge',width = 0.75)+scale_fill_brewer(name="Spur length",palette = "Blues")+geom_errorbar(aes(ymin = mean_height - sd_height, ymax = mean_height + sd_height), position = position_dodge(0.75),width=0.3,alpha=0.6)+theme_bw()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),panel.grid=element_blank())+ylab("Cell Height")+xlab("Spur of Area")+theme(legend.position = c(0.98,0.98),legend.justification = c("right","top"))
topptx(filename = "Spur_Height.pptx",width = 7,height = 7)
read.table("width.txt",header = T,sep = "\t") -> width
width_mean <- width %>% group_by(group1, group2) %>% summarise(mean_width = mean(width), sd_width = sd(width))
width_mean$group2 <- factor(width_mean$group2,levels = c("1.5cm","4cm","6cm","9cm","12cm","mature"))
width_mean$group1 <- factor(width_mean$group1,levels = c("Base","BMM","Middle","TMM","Tip"))
ggplot(width_mean,aes(x=group1,y=mean_width,fill=group2))+geom_bar(stat="identity",position = 'dodge',width = 0.75)+scale_fill_brewer(name="Spur length",palette = "Blues")+geom_errorbar(aes(ymin = mean_width - sd_width, ymax = mean_width + sd_width), position = position_dodge(0.75),width=0.3,alpha=0.6)+theme_bw()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),panel.grid=element_blank())+ylab("Cell Width")+xlab("Spur of Area")+theme(legend.position = "none")
topptx(filename = "Spur_Width.pptx",width = 7,height = 7)
read.table("e.txt",header = T,sep = "\t") -> e
e_mean <- e %>% group_by(group1, group2) %>% summarise(mean_e = mean(e), sd_e = sd(e))
e_mean$group2 <- factor(e_mean$group2,levels = c("1.5cm","4cm","6cm","9cm","12cm","mature"))
e_mean$group1 <- factor(e_mean$group1,levels = c("Base","BMM","Middle","TMM","Tip"))
ggplot(e_mean,aes(x=group1,y=mean_e,fill=group2))+geom_bar(stat="identity",position = 'dodge',width = 0.75)+scale_fill_brewer(name="Spur length",palette = "Blues")+geom_errorbar(aes(ymin = mean_e - sd_e, ymax = mean_e + sd_e), position = position_dodge(0.75),width=0.3,alpha=0.6)+theme_bw()+theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),panel.grid=element_blank())+ylab("Cell Anisotropy")+xlab("Spur of Area")+theme(legend.position = c(0.98,0.98),legend.justification = c("right","top"))
topptx(filename = "Spur_E.pptx",width = 7,height = 7)
read.table("target_counts.txt",header = T,sep = "\t")
read.table("target_counts.txt",header = T,sep = "\t") -> datqa
read.table("target_counts.txt",header = T,sep = "\t") -> data
View(data)
tapply(data$counts,data[,c("gene_id","condition")],mean)
read.table("target_counts.txt",header = T,sep = "\t") -> data
tapply(data$counts,data[,c("gene_id","condition")],mean)
tapply(data$counts,data[,c("gene_id","condition")],mean)
read.table("target_counts.txt",header = T,sep = "\t") -> data
tapply(data$counts,data[,c("gene_id","condition")],mean)
tapply(data$counts,data[,c("gene_id","condition")],as)
tapply(data$counts,data[,c("gene_id","condition")],sd)
tapply(data$counts,data[,c("gene_id","condition")],sd) -> data_sd
View(data_sd)
tapply(data$counts,data[,c("gene_id","condition")],mean) -> data_mean
rownames(data_sd) <- ("Day_sd","Night_sd")
rownames(data_sd) <- c("Day_sd","Night_sd")
rownames(data_sd) <- c("Day_sd,Night_sd")
names(data_sd) <- c("Day_sd,Night_sd")
View(data_sd)
names(data_sd) <- c("Day_sd","Night_sd")
df_sd <- as.matrix(data_sd)
View(df_sd)
names(df_sd) <- c("Day_sd","Night_sd")
View(data_mean)
View(data_sd)
colnames(data_sd) <- c("Day_sd","Night_sd")
colnames(data_mean) <- c("Day_mean","Night_mean")
data.frame(data_mean,data_sd)
data.frame(data_mean,data_sd) -> df
View(df)
data.frame(data$gene_id,df)
data.frame(df,data$group)
read.table("target_counts.txt",sep = "\t",header = T) -> group
data.frame(df,group$group)
read.table("target_counts.txt",sep = "\t",header = T,row.names = 1) -> group
View(group)
ls -l
ls
View(df)
write.table(df,file="gene_count.txt",quote = F,sep = "\t")
read.table("gene_count.txt",header = T,sep = "\t") -> plot
View(plot)
library(ggplot)
library(ggplot2)
library(eoffice)
ggplot(plot,aes(x=gene_id,y=Day_mean))+geom_bar()
ggplot(plot,aes(x=gene_id,y=Day_mean))+geom_bar(stat='identity',position="dodge")
read.table("gene_count.txt",header = T,sep = "\t") -> plot
ggplot(plot,aes(x=gene_id,y=Day_mean,fill="group1"))+geom_bar(stat='identity',position="dodge")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge")+geom_errorbar(aes(x=group1,ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd))
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge")+geom_errorbar(aes(x=gene_id,ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd))
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge")+geom_errorbar(aes(x=gene_id,ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.1)
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge")+geom_errorbar(aes(x=gene_id,ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.1,position = position_dodge(0.5))
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge")+geom_errorbar(aes(x=gene_id,ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.1,position = position_dodge(0.75))
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge")+geom_errorbar(aes(x=gene_id,ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.3,position = position_dodge(0.75))
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge")+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.3,position = position_dodge(0.75))
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 0.75)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.3,position = position_dodge(0.75))
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.3,position = position_dodge(1))
View(data)
read.table("target_counts.txt",header = T,sep = "\t") -> data
read.table("target_counts.txt",header = T,sep = "\t") -> data
tapply(data$counts,data[,c("gene_id","condition")],mean)
tapply(data$counts,data[,c("gene_id","condition")],sd) -> data_sd
tapply(data$counts,data[,c("gene_id","condition")],mean) -> data_mean
colnames(data_sd) <- c("Day_sd","Night_sd")
colnames(data_mean) <- c("Day_mean","Night_mean")
data.frame(data_mean,data_sd) -> df
write.table(df,file="gene_count.txt",quote = F,sep = "\t")
read.table("gene_count.txt",header = T,sep = "\t") -> plot
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.3,position = position_dodge(1))
View(plot)
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.3,position = position_dodge(1))+facet_grid(~group)
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.3,position = position_dodge(1))+facet_warp(~group)
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.3,position = position_dodge(1))+facet_wrap(~group)
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.3,position = position_dodge(1))+facet_wrap(~group,sacles="free")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.3,position = position_dodge(1))+facet_wrap(~group,scales="free")
read.table("gene_count.txt",header = T,sep = "\t") -> plot
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.3,position = position_dodge(1))+facet_wrap(~group,scales="free")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(1))+facet_wrap(~group,scales="free")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(1))+facet_wrap(~group,scales="free")+theme_bw()
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(1))+facet_wrap(~group,scales="free")+theme_classic()
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(1))+facet_wrap(~group,scales="free")+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("raw read counts")+xlab("Gene ID")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(1))+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("raw read counts")+xlab("Gene ID")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(1))+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Raw read counts")+xlab("Gene ID")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(1))+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Raw read counts")+xlab("Gene ID")+scale_fill_brewer(name="Condition",palette = "Greys")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(1))+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Raw read counts")+xlab("Gene ID")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(1))+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Raw read counts")+xlab("Gene ID")+scale_fill_brewer(type='qual')
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 0.7)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(0.7))+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Raw read counts")+xlab("Gene ID")+scale_fill_brewer(type='qual')
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(1))+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Raw read counts")+xlab("Gene ID")+scale_fill_brewer(name="Condition",palette = "Dark2")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 0.7)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(0.7))+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Raw read counts")+xlab("Gene ID")+scale_fill_brewer(type='qual')
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 1)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(1))+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Raw read counts")+xlab("Gene ID")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 0.7)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.15,position = position_dodge(0.7))+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Raw read counts")+xlab("Gene ID")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 0.7)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.1,position = position_dodge(0.7))+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Raw read counts")+xlab("Gene ID")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 0.7)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.1,alpha=0.7,position = position_dodge(0.7)+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Raw read counts")+xlab("Gene ID")
ggplot(plot,aes(x=gene_id,y=Day_mean,fill=group1))+geom_bar(stat='identity',position="dodge",width = 0.7)+geom_errorbar(aes(ymin=Day_mean-Day_sd,ymax=Day_mean+Day_sd),width=0.1,alpha=0.7,position = position_dodge(0.7))+facet_wrap(~group,scales="free")+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Raw read counts")+xlab("Gene ID")
topptx(filename = "gene_bar_plot.pptx",width = 8,height7)
topptx(filename = "gene_bar_plot.pptx",width = 8,height=7)
topptx(filename = "gene_bar_plot.pptx",width = 13,height=7)
+geom_errorbar(aes(ymin = mean_e - sd_e, ymax = mean_e + sd_e), position = position_dodge(0.75),width=0.3,alpha=0.6)+theme_bw()+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),panel.grid=element_blank())+ylab("Cell Anisotropy")+xlab("Spur of Area")
topptx(filename = "Spur_E.pptx",width = 11,height = 6)
e_mean$group1 <- factor(e_mean$group1,levels = c("1.5cm","4cm","6cm","9cm","12cm","mature"))


