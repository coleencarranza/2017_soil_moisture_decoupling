#CHOOSE DIRECTORY
setwd("/media/coleen/DDrive/1_OWASIS/4_Papers/ts_dnlm/data")
stations<-list.files(pattern = "_imp.csv")
stations.list<-lapply(stations,function(x) read.csv(x, header=T,stringsAsFactors = F))
names<-c("SM05","SM09","SM13","SM20")
stations.list<-setNames(stations.list,names)
stations.list<-lapply(stations.list,function(x) {x[, "YYYYMMDD"] <- as.Date(x[, "YYYYMMDD"],
                       format = "%Y-%m-%d");x})

#PLOTTING TIME SERIES
#XTS SERIES - ORIGINAL DATA
require(xts)
library(zoo)
xts.list<-lapply(stations.list,function(x) xts(x[,-10],x[,10]))


##PLOTTING XTS TIME SERIES########
#pdf(paste("/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/figures/station_ts_rev_v2.pdf"),
#   height=11 ,width=13, onefile = TRUE) 
par(mfrow= c(4,1),mar = c(0,0,0,0) + 0.01,mgp=c(1.5,0.55,0),oma = c(3.75,3.75,1,3.75), tck=-0.01)
for(i in 1:4){
    myColors <- c("skyblue1", "dodgerblue4", "grey30", "darkblue", "darkviolet")
  #Plot rainfall first so Soil moisture TS is more visible
  plot(as.zoo(xts.list[[i]][,9]),type="h",lwd =3,  col="snow3", yaxt="n",xlab="",ylab="", xaxt = "n",
       xlim= c(as.Date('2014-01-01', format="%Y-%m-%d"),as.Date('2016-12-31', format="%Y-%m-%d")),
       ylim = rev(c(0,80)))
  axis(side=4,at = seq(0,80,by = 10),labels = c("",10,"",30,"",50,"",70,""),cex.axis=1.8,mgp=c(0.75,0.65,0), tck=-0.01)
  par(new = T)
  plot(as.zoo(xts.list[[i]][,c(1,3)]),yaxt="n",xlab="",ylab="", xaxt = "n",
       xlim=c(as.Date('2014-01-01', format="%Y-%m-%d"),
       as.Date('2016-12-31', format="%Y-%m-%d")),cex.lab= 1.7,
       col = myColors, lwd=2,screens=1, ylim = c(0,0.60)) 
  text(as.Date('2016-12-28', format="%Y-%m-%d"),0.02, labels = names(xts.list[i]),cex=2,font =2)
  axis(side=2,at=seq(0,0.6,by=0.1), labels = c("",0.1,"",0.3,"",0.5,""),cex.axis=1.8,mgp=c(0.75,0.35,0), tck=-0.01)
  if (i %in% 4)
  axis.Date(1, at = seq(as.Date('2013-01-01', format="%Y-%m-%d"), 
                as.Date('2017-01-01', format="%Y-%m-%d"), by = "years"), cex.axis= 1.8)
}
mtext(expression("Soil Moisture" ~ (cm^3~cm^{-3})), side=2, line=1.35,outer = TRUE,cex = 1.3)
mtext(side = 4, line=2.25, 'Rainfall (mm)', outer = TRUE,cex  = 1.35)
mtext(side = 1, line = 2, "Time", outer = TRUE,cex  = 1.35)
legend(x="bottomleft",legend = c("5cm","40cm", "rainfall"), y.intersp = 0.75, x.intersp = 0.25,bty = "n",
       seg.len=2, lwd = 3,lty = 1, col = myColors,cex=2)
#dev.off()



###PLOTTING CROSS CORRELATION###ORIGINAL
#pdf(paste("/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/figures/xcor_v2.pdf"),
 #   height=5,width=5)
par(mfrow = c(2,2),mar = c(0, 0, 1.2, 1), oma = c(3, 3, 0.25, 0.25),mgp = c(1, 0.4, 0), tck = -0.009)
for(i in seq_along(stations.list)){
  ccf<-ccf(as.ts(stations.list[[i]]$VWC5),as.ts(stations.list[[i]]$VWC40), lag.max = 10,xaxt="n",yaxt = "n",
           ylim = c(0.6,0.85), main="",col = "grey60",lwd=5.5,ylab = "",xlab = "")
  if (i %in% c(3, 4))
  axis(1, cex.axis = 1.2)
  if (i %in% c(1, 3))
   axis(2, cex.axis = 1.2)
  title( main =names(stations.list[i]), cex.main = 1.3)
  abline(v=0, col = "blue",lty= 5, lwd=1)  
}
mtext("Lag (in Days)", side = 1, outer = TRUE, cex=1.25, line=1.8)
mtext("Correlation", side = 2, outer = TRUE, cex =1.25, line =1.6)
#dev.off()



####PLOTTING SCATTER PLOTS 40cm vs 5cm SOIL MOSITURE
##USE data plotting only
stations<-list.files(pattern = "_plotting.csv")
stations.list<-lapply(stations,function(x) read.csv(x, header=T,stringsAsFactors = F))
names<-c("SM05","SM09","SM13","SM20")
stations.list<-setNames(stations.list,names)
stations.list<-lapply(stations.list,function(x) {x[, "YYYYMMDD"] <- as.Date(x[, "YYYYMMDD"],
              format = "%Y-%m-%d");x})

#PACKAGE load
require("graphics");library(RColorBrewer);library(scales); library(ecp);
library(plyr)

#pdf(paste("/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/figures/scatter_v2.pdf"),
#   height=5.5,width=18)
#LINEAR AND LOESS REGRESSION
colB<-palette(colorRampPalette(brewer.pal(11,"Spectral"))(12))
par(mfrow = c(1,4),mar = c(0,0,2,0.5), oma = c(5,6,0,0),mgp = c(1, 1, 0), tck = -0.009)

for(i in seq_along(stations.list)){
#reg fit
lo<-loess(stations.list[[i]]$VWC40 ~ stations.list[[i]]$VWC5, span=0.75,degree=2,se=TRUE)### span can be changed for smoothness
lm<-lm(stations.list[[i]]$VWC40~stations.list[[i]]$VWC5)
plot(c(0,0.5), c(0,0.6), type="n", cex.main= 1.5,xaxt="n",yaxt = "n",ylab = "",xlab = "")
points(stations.list[[i]]$VWC40~stations.list[[i]]$VWC5,pch=21, bg=alpha(colB[stations.list[[i]]$mon],0.7),
       col=alpha("black",0.35),cex=ceiling(stations.list[[i]]$RDmm/5),lwd =0.5)
#fit loess line
j <- order(lo$x)
lines(lo$x[j],lo$fitted[j], lwd=2.5)

#cond mean +var+ bars
mean.agg<-aggregate(stations.list[[i]]$VWC40,by=list(round(stations.list[[i]]$VWC5,2)),mean,na.rm=TRUE)
var.agg<-aggregate(stations.list[[i]]$VWC40,by=list(round(stations.list[[i]]$VWC5,2)),var,na.rm=TRUE)
sd.agg<-aggregate(stations.list[[i]]$VWC40,by=list(round(stations.list[[i]]$VWC5,2)),sd,na.rm=TRUE)
sd.plus<- mean.agg$x + sd.agg$x
sd.minus<-mean.agg$x - sd.agg$x
#lm fit line
lines(lm$model$`stations.list[[i]]$VWC5`,lm$fitted.values,col = "red3", lwd = 2, lty=3)
#correlation grouped VWC values
x<-data.frame(cbind(stations.list[[i]]$VWC5,stations.list[[i]]$VWC40))
x$x<-round(x$X1/.05)*.05 #round(x$X1,2)
y<-ddply(x, "x", summarize, corr=cor(X1, X2,method = "spearman",use="pairwise.complete.obs"))#use="pairwise.complete.obs"
frequ<-data.frame(count(round(stations.list[[i]]$VWC5/0.05)*0.05))
corr.list<-merge(frequ,y, by="x")
arrows(mean.agg$Group.1, sd.minus, mean.agg$Group.1, sd.plus, lwd=1.25,col=alpha("black",0.7),
       length=0.05, angle=90, code=3)
points(mean.agg, col=alpha("black",0.7), pch=15, cex=1.5)

#add text
cor<-round(cor(stations.list[[i]]$VWC5,stations.list[[i]]$VWC40, method = "spearman",use="pairwise.complete.obs"),digits=3)
sigma.lm<-round(sqrt(sum(lm$residuals^2)/lm$df.residual),3)
sigma.lo<-round(sqrt(sum(lo$residuals^2)/lo$n),3)
points(0.7,0.4, pch="")
text(0.08,0.605,labels = bquote('R'[s]==.(cor)), cex = 2, col = "grey38")
text(0.08,0.57,labels = bquote('RSE'[italic(lm)]==.(sigma.lm)), cex = 2, col = "grey38")
text(0.08,0.535,labels = bquote('RSE'[italic(lo)]==.(sigma.lo)), cex = 2, col = "grey38")
lm_coef <- round(coef(lm), 3) # extract coefficients lm

##AXES scales
if (i %in% c(1:4))
  axis(1, cex.axis = 2)
if (i %in% (1))
  axis(2, cex.axis = 2)
if (i %in% (4)){
  legend(x=0.37,y=0.12,legend = month.abb[1:12], y.intersp = 0.67, x.intersp = 0.5,bty = "n",
         pt.cex=1.75, cex=1.5, col = alpha("black",0.55), pt.bg =alpha(colB,0.75), pch=21, ncol=2)
  legend(x=0,y=0.055,legend = c("linear fit", "loess fit","conditional mean +/- sd"),lty =c(6,1,1),cex=1.5,
         lwd=c(2,3,1), seg.len=2.5,x.intersp = 0.75,y.intersp = 0.7,bty = "n", col=c("red3","black",alpha("black",0.7)))
  c<-sort(unique(ceiling(stations.list[[2]]$RDmm/10)),decreasing = FALSE)[-1]
  ##spacing for rainfall should be adjusted every time! - remove  some numbers from sequence to avoid duplicate points)
  ## for SM09 - x.intersp=c(1.85,1.85,1.85,1.85,1.85,2.7),y.intersp = c(1,1,1.2,1.4,1.6,1.9)
  ##for other 3 stations -- x.intersp=c(1.85,1.85,1.85,1.85),y.intersp = c(1,1,1.2,1.4)
  legend(x=0.42,y=0.56,legend=c*10,pch=21,pt.bg=alpha("grey",0.6),col = "black", pt.cex=c*2, bty="n", cex=1.5,
         x.intersp=2,y.intersp = c(1,1,1.17,1.37,1.59,1.92))
  text(0.45,0.56, labels = "Rainfall (mm)", cex=1.5)
}

title( main =names(stations.list[i]), cex.main=2)

}
mtext(expression("Soil Moisture at 5cm " ~(cm^3~cm^{-3})), side = 1, outer = TRUE, cex=1.75, line=4.5)
mtext(expression("Soil Moisture at 40cm "~(cm^3~cm^{-3})), side = 2, outer = TRUE, cex=1.75, line=2.8)
#dev.off()

