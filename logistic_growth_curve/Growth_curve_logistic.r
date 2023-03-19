library(FSAdata) # for data
library(FSA)     # for vbFuns(), vbStarts(), confint.bootCase()
library(car)     # for Boot()
library(dplyr)   # for filter(), mutate()
library(ggplot2)
library(reshape)
library(deSolve)
#use melt(reshape) for wide data to long data(point plot) 
read.csv("length.csv",header="TRUE") -> Length
melt(Length,id='day') -> melt_data
write.table(melt_data,file="melt.txt",sep="\t")
# here also could use a function named "getInitial" from the package "stats"  to get logistic nonlinear initial parameter
#logF <- logisticFuns()
#logSS <- getInitial(value ~ SSlogis(day, alpha, xmid, scale), data = melt)
#differenct between the logistic and get Initial.
K_start <- logSS["alpha"]
R_start <- 1/logSS["scale"]
N0_start <- logSS["alpha"]/(exp(logSS["xmid"]/logSS["scale"])+1)
log_formula <- formula(value ~ K*N0*exp(R*day)/(K + N0*(exp(R*day) - 1)))
formu<-nls(log_formula,data=melt, start = list(K = K_start, R = R_start, N0 = N0_start))
f.boot1 <- Boot(formu)
confint(f.boot1)
days <- seq(0,22,by=0.2)#22 is my time scale
predict2 <- function(x) predict(x,data.frame(day=days))
predict2(log_fit) 
log_boot <- Boot(log_fit,f=predict2)  
log_preds <- data.frame(days,predict(log_fit,data.frame(day=days)),confint(log_boot))
names(log_preds) <- c("day","fit","LCI","UCI")
loggrowthPlot <- ggplot() + 
  geom_ribbon(data=log_preds,aes(x=day,ymin=LCI,ymax=UCI),fill="gray90") +
  geom_point(data=melt,aes(y=value,x=day),size=2,alpha=0.1) +
  geom_line(data=log_preds,aes(y=fit,x=day),size=1) +
  scale_y_continuous(name="Total Length (mm)",limits=c(0,700),expand=c(0,0)) +
  scale_x_continuous(name="Day",expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid=element_blank())
loggrowthPlot
topptx(loggrowthPlot,filename="Log_grow_curve.pptx",width=8,height=7)