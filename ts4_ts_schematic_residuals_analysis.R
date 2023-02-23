###Generate correlated points fron mvrnorm
library(MASS)
##SET1
mean1<-c(15,18.8)
cor1 = 0.6
#covar = cor* (sqrt(var1),sqrt(var2))
var1 = 6.5*6.5
var2 = 6*6
cov1<-cor1*(6*6.5)
vcov1<-matrix(c(var1,cov1,cov1,var2),nrow =2,byrow = T)

trek1<-mvrnorm(n=400,mu=mean1,Sigma =vcov1)
plot(trek1[,1],trek1[,2])


#SET2
mean2<-c(28,24.5)
cor2 = 0.85
#covar = cor* (sqrt(var1),sqrt(var2))
var3 = 6.5*6.5
var4 = 3*3
cov2<-cor2*(6.5*3)
vcov2<-matrix(c(var3,cov2,cov2,var4),nrow =2,byrow = T)

trek2<-mvrnorm(n=400,mu=mean2,Sigma =vcov2)
plot(trek2[,1],trek2[,2])


all<-rbind(trek1,trek2)
all<-all[all[,1]>0,]
all<-all/100 ######## make values not in percent

pdf(paste("/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/figures/sche_resid_analy", ".pdf", sep = ""),
    height=5,width=10)
##SCATTER PLOT
layout(matrix(1:2,ncol=2), widths = c(5,7),height = c(1,1))
par(mar = c(3,3,0.5,0.5), mgp = c(2, 0.35, 0), tck = -0.01)
plot(all[,1],all[,2], pch = 19, col = "grey60", cex = 0.7,yaxt = "n", xlab = "Surface soil moisture",ylab = "Subsurface soil moisture",cex.lab= 1.2)
axis(side=2, cex.axis = 1, las =2)  

###Add conditional mean + st. dev bars
#cond mean +var+ bars
mean.agg<-aggregate(all[,2],by=list(round(all[,1],2)),mean,na.rm=TRUE)
var.agg<-aggregate(all[,2],by=list(round(all[,1],2)),var,na.rm=TRUE)
sd.agg<-aggregate(all[,2],by=list(round(all[,1],2)),sd,na.rm=TRUE)
sd.plus<- mean.agg$x + sd.agg$x
sd.minus<-mean.agg$x - sd.agg$x

arrows(mean.agg$Group.1, sd.minus, mean.agg$Group.1, sd.plus, lwd=1,col="grey20",
       length=0.05, angle=90, code=3)
points(mean.agg, col="grey20", pch =15)

####Fit a loess function
lo<-loess(all[,2]~all[,1], span=0.4,degree=2,se=TRUE)### span can be changed for smoothness
j <- order(lo$x)
lines(lo$x[j],lo$fitted[j], lwd=2.5, col = "red4")


##Plot residuals from LOESS
library(scales); library(ecp);
library(plyr)

lo.res<-resid(lo)
###Function for rounding to different sequences##
mround <- function(x,base){
  base*round(x/base)}
#cond mean +var+ bars
mean.lores<-aggregate(lo.res,by=list(mround(lo$x,0.02)),mean,na.rm=TRUE)
#var
var.lores<-aggregate(lo.res,by=list(mround(lo$x,0.02)),var,na.rm=TRUE)
var.lores[is.na(var.lores)]<-0
var.lores$cum<-cumsum(var.lores$x)
x<-as.matrix(var.lores$cum)
infl<-e.divisive(diff(x), sig.lvl=.01,k=NULL,min.size=2,alpha = 1, R=199)


##Plot Residuals
par(mar = c(3,3,0.5,4), mgp = c(2, 0.35, 0), tck = -0.01)
plot(lo$x[j],lo.res[j],yaxt = "n",ylab = "Residuals", xlab="Surface soil moisture",pch =19, col = "grey60", cex = 0.75, cex.lab =1.2)
axis(side=2, cex.axis = 1, las =2)    
abline(0,0, lwd=1)
text(x = 0.235, y = 0.18,labels = expression(paste(italic(theta[c]))), cex = 1.25)

par(new = T)
plot(var.lores$Group.1,var.lores$cum,type="l",lwd=2.2,xlab="",ylab="",xlim = c(0, max(lo$x)),
     xaxt="n",yaxt = "n",ylim = c(0, max(var.lores$cum)))
points(var.lores$Group.1,var.lores$cum,pch = 15, cex = 0.65)
abline(0,0, lwd = 0.5)
axis(side=4, cex.axis = 1, las =2)
lines(var.lores$Group.1,var.lores$x, type = "h", lwd=3.5,col="grey38")
#add ecp estimates
inf<-infl$estimates
abline(v =var.lores$Group.1[c(inf)][-c(1,length(inf))], col = "grey38", lty = 3, lwd = 2.5)
inf.list<-var.lores$Group.1[c(inf)]
mtext('Cumulative Residual Variance', cex= 1.2, side =4, line=3)

dev.off()

