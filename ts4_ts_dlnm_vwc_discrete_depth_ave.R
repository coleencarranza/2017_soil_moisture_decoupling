####PLOTTING SCATTER PLOTS 40cm vs 5cm SOIL MOSITURE
##USE data plotting only
#CHOOSE DIRECTORY
setwd("/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/data/ITC_SM")
stations<-list.files(pattern = ".csv")
stations.list<-lapply(stations,function(x) read.csv(x, header=T,stringsAsFactors = F))
names<-c("SM05","SM09","SM13","SM20")
stations.list<-setNames(stations.list,names)
stations.list<-lapply(stations.list,function(x) {x[, "date"] <- as.Date(x[, "date"],
                                                                        format = "%Y-%m-%d");x})
#PACKAGE load
require("graphics");library(RColorBrewer);library(scales); library(ecp);
library(plyr);library(imputeTS)


######CALCULATE Depth-averaged values using Qiu, 2012
#data.ts<-data
stations.list<-lapply(stations.list, function(x)  {x$VWC5<-na.ma(x$VWC5, k=6, weighting = "linear");x})
stations.list<-lapply(stations.list, function(x)  {x$VWC10<-na.ma(x$VWC10, k=6, weighting = "linear");x})
stations.list<-lapply(stations.list, function(x)  {x$VWC20<-na.ma(x$VWC20, k=6, weighting = "linear");x})
stations.list<-lapply(stations.list, function(x)  {x$VWC40<-na.ma(x$VWC40, k=6, weighting = "linear");x})



###change formula-calculate using Qiu, 2012
stations.list<-lapply(stations.list, function(x) {x<- within(x,VWC_zone_Qiu <-(VWC5*5 +((VWC5+VWC10)/2)*(10-5)+((VWC10+VWC20)/2)*(20-10)
                          +((VWC40+VWC20)/2)*(40-20))/40);x})

#CALCULATE ZONE AVERAGE SOIL MOISTURE (Gretchen Miller et al., 2007)
#data<- within(data,VWC_zone <-(VWC5*7.5+VWC10*7.5 +VWC20*15 +VWC40*30)/60)
#data<- within(data,VWC_zonep <- VWC_zone*100)
#plot(data$date, data$VWC5, type="l")



#####Create ts for plotting Time series
#PLOTTING TIME SERIES
#XTS SERIES - ORIGINAL DATA
require(xts)
xts<-lapply(stations.list, function(x) xts(x[,-1], x[,1]))


pdf(paste("/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/figures/depth_ave_40cm_v2.pdf"),
height=9,width=13)
par(oma = c(3,4,0,1))
layout(matrix(1:8,ncol=2,byrow =FALSE), widths = c(7,3),height = c(1,1))
#layout.show(8)
#########PLOTTING TIME SERIES##########################------------------------------------------
for(i in 1:4){
  par(mar = c(1,1.5,0.5,0),mgp = c(1, 0.7, 0), tck = -0.015)
  myColors <- c("royalblue2","mediumblue")
  plot(as.zoo(xts[[i]][,c(7,11)]),yaxt="n",xlab="",ylab="", xaxt = "n",lty= c(6,1),
       xlim=c(as.Date('2014-01-01', format="%Y-%m-%d"),as.Date('2016-12-31', format="%Y-%m-%d")),
       col = alpha(myColors,0.75), lwd=c(2,1.5),screens=1, ylim = c(0,0.60)) 
  text(as.Date('2016-12-28', format="%Y-%m-%d"),0.02, labels = names(xts[i]),cex=2,font =2)
  axis(side=2,at=seq(0,0.6,by=0.1), labels = c("",0.1,"",0.3,"",0.5,""),cex.axis=1.8,mgp=c(0.75,0.35,0), tck=-0.01)
  if (i %in% 4)
    axis.Date(1, at = seq(as.Date('2013-01-01', format="%Y-%m-%d"), 
                          as.Date('2017-01-01', format="%Y-%m-%d"), by = "years"), cex.axis= 1.8)
}
mtext(expression("Soil Moisture" ~ (cm^3~cm^{-3})), side=2, line=1.35,outer = TRUE,cex = 1.3)
mtext(side = 1, line = 2, "Time", outer = FALSE,cex  = 1.35)
legend(x="bottomleft",legend = c("40cm","depth-average"), y.intersp = 0.75, x.intersp = 0.25,bty = "n",
       seg.len=1.5, lwd = 2,lty = c(6,1), col = myColors,cex=2)
############SCATTERPLOT + Correlation############################----------------------------------------------
for(i in 1:4){
  par(mar = c(1.75,6,1,4),mgp = c(2, 0.5, 0), tck = -0.015)
plot(stations.list[[i]]$VWC_zone_Qiu~stations.list[[i]]$VWC40,xlab="",ylab="", cex.axis = 1.3, cex = 1.3, col = alpha("grey38",0.5))
cor=round(cor(stations.list[[i]]$VWC40,stations.list[[i]]$VWC_zone_Qiu),3)
xx<-min(stations.list[[i]]$VWC40)
yy<-max(stations.list[[i]]$VWC_zone_Qiu)
text(x=xx,y=yy,labels=paste0("R = ",cor), cex = 1.2, offset=1)
if (i %in% 4)
  axis.Date(1, at = seq(as.Date('2013-01-01', format="%Y-%m-%d"), 
                        as.Date('2017-01-01', format="%Y-%m-%d"), by = "years"), cex.axis= 1.8)
}
mtext(expression("Depth-Average Soil Moisture" ~ (cm^3~cm^{-3})), side=4, line=-0.7,outer=TRUE,cex = 1.33)
mtext("Soil Moisture at 40 cm ", side = 1, outer = FALSE, cex=1.35, line=2.7)

dev.off()
