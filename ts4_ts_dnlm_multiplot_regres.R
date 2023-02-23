####PLOTTING SCATTER PLOTS 40cm vs 5cm SOIL MOSITURE
##USE data plotting only
#CHOOSE DIRECTORY
setwd("/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/data")
stations<-list.files(pattern = "_plotting.csv")
stations.list<-lapply(stations,function(x) read.csv(x, header=T,stringsAsFactors = F))
names<-c("SM05","SM09","SM13","SM20")
stations.list<-setNames(stations.list,names)
stations.list<-lapply(stations.list,function(x) {x[, "YYYYMMDD"] <- as.Date(x[, "YYYYMMDD"],
                       format = "%Y-%m-%d");x})
#PACKAGE load
require("graphics");library(RColorBrewer);library(scales); library(ecp);
library(plyr)

#pdf(paste("/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/figures/scatter_regres.pdf"),
 #height=12,width=15)
#LINEAR AND LOESS REGRESSION
colB<-palette(colorRampPalette(brewer.pal(11,"Spectral"))(12))
par(oma = c(5,5,0,5))
layout(matrix(1:8,ncol=2,byrow =FALSE), widths = c(3,7),height = c(1,1))
#layout.show(8)

for(i in seq_along(stations.list)){
par(mar = c(0,1.5,0.5,0),mgp = c(1, 0.7, 0), tck = -0.015)
  #reg fit
  lo<-loess(stations.list[[i]]$VWC40 ~ stations.list[[i]]$VWC5, span=0.75,degree=2,se=TRUE)### span can be changed for smoothness
  lm<-lm(stations.list[[i]]$VWC40~stations.list[[i]]$VWC5)
  plot(c(0,0.5), c(0,0.6), type="n", cex.main= 1.5,xaxt="n",yaxt = "n",ylab = "",xlab = "")
  points(stations.list[[i]]$VWC40~stations.list[[i]]$VWC5,pch=21, bg=alpha(colB[stations.list[[i]]$mon],0.7),
         col=alpha("black",0.25),cex=ceiling(stations.list[[i]]$RDmm/5),lwd =0.25)
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
  points(mean.agg, col=alpha("black",0.85), pch=15, cex=1)
  
  #add text
  cor<-round(cor(stations.list[[i]]$VWC5,stations.list[[i]]$VWC40, method = "spearman",use="pairwise.complete.obs"),digits=3)
  sigma.lm<-round(sqrt(sum(lm$residuals^2)/lm$df.residual),3)
  sigma.lo<-round(sqrt(sum(lo$residuals^2)/lo$n),3)
  points(0.7,0.4, pch="")
  text(0.08,0.59,labels = bquote('R'[s]==.(cor)), cex = 2, col = "grey38")
  text(0.08,0.54,labels = bquote('RSE'[italic(lm)]==.(sigma.lm)), cex = 2, col = "grey38")
  text(0.08,0.49,labels = bquote('RSE'[italic(lo)]==.(sigma.lo)), cex = 2, col = "grey38")
  lm_coef <- round(coef(lm), 3) # extract coefficients lm
  text(0.46,0.02, labels = names(stations.list[i]),cex=2,font =2)
  
  ##AXES scales
  if (i %in% c(1:4))
    axis(2,at=seq(0,0.6,by=0.1),labels=c(0.0,"",0.2,"",0.4,"",0.6), cex.axis = 1.75, las =2)
  if (i %in% (4)){
  axis(1, cex.axis = 1.75)
  legend(x=0,y=0.075,legend = c("linear fit", "loess fit","conditional mean +/- sd"),lty =c(6,1,1),cex=1.2,
         lwd=c(2,3,1), seg.len=2.5,x.intersp = 0.75,y.intersp = 0.68,bty = "n", col=c("red3","black",alpha("black",0.7)))
}
  #title( main =names(stations.list[i]), cex.main=2)
}
mtext("Soil Moisture at 5cm ", side = 1, outer = TRUE, cex=1.75, line=3.2)
mtext("Soil Moisture at 40cm ", side = 2, outer = TRUE, cex=1.75, line=2.35)


############LOESS function Residuals#########
for(i in seq_along(stations.list)){
par(mar = c(0,8,0.5,2),mgp = c(1, 0.7, 0), tck = -0.015)
#resid lo fit
  lo<-loess(stations.list[[i]]$VWC40 ~ stations.list[[i]]$VWC5, span=0.75,degree=2,se=TRUE)### span can be changed for smoothness\
  lo.res<-resid(lo)
  ###Function for rounding to different sequences##
  mround <- function(x,base){
    base*round(x/base)}
  #cond mean +var+ bars
  mean.lores<-aggregate(lo.res,by=list(mround(lo$x,0.01)),mean,na.rm=TRUE)
  #var
  var.lores<-aggregate(lo.res,by=list(mround(lo$x,0.01)),var,na.rm=TRUE)
  var.lores[is.na(var.lores)]<-0
  var.lores$cum<-cumsum(var.lores$x)
  x<-as.matrix(var.lores$cum)
  infl<-e.divisive(diff(x), sig.lvl=.05,k=NULL,min.size=2,alpha = 1, R=199)
  
  plot(lo$x,lo.res, main = "",pch=21, bg=alpha(colB[stations.list[[i]]$mon],0.7),
       col=alpha("grey38",0.25), ylab="", xlab="", yaxt="n",
       xaxt="n", cex=ceiling(stations.list[[i]]$RDmm/5), lwd=0.5, ylim=c(-0.2, 0.2), 
       xlim = c(0,0.5))
  abline(0,0, lwd=0.5)
  if (i %in% c(1:4))
  axis(2, cex.axis = 1.75, las =2)
  if (i %in% (4)){
    axis(1, cex.axis = 1.75,at=seq(0,0.5,by=0.05))
  }
  
  par(new = T)
  plot(var.lores$Group.1,var.lores$cum,type="l",lwd=1.25,xlab="",ylab="",
       xaxt="n",yaxt = "n",ylim = c(0, max(var.lores$cum)),xlim = c(0,0.5))
  points(var.lores$Group.1,var.lores$cum,pch = 15, cex = 0.5)
  abline(0,0, lwd = 0.5)
  axis(side=4, cex.axis = 1.35, las =2)
  lines(var.lores$Group.1,var.lores$x, type = "h", lwd=3.5,col="grey38")
  #add ecp estimates
  inf<-infl$estimates
  inf<-inf[-c(1,length(inf))]
  abline(v =var.lores$Group.1[c(inf)], col = "grey38", lty = 3, lwd = 2)
  inf.list<-var.lores$Group.1[c(inf)]
  if (i %in% (4)){
    legend(x=0.36,y=0.0075,legend = month.abb[1:12], y.intersp = 0.67, x.intersp = 0.5,box.lwd = 1,box.col = "white",bg = "white",
           pt.cex=2, cex=1.8, col = alpha("black",0.55), pt.bg =alpha(colB,0.8), pch=21, ncol=2)
    c<-sort(unique(ceiling(stations.list[[2]]$RDmm/10)),decreasing = FALSE)[-1]
    legend(x=0.45,y=0.0148,legend=c*10,pch=21,pt.bg="whitesmoke",col = "grey40", pt.cex=c*2, cex=1.5,
           x.intersp=3,y.intersp = c(1,1,1.17,1.37,1.59,1.92),box.lwd = 0,box.col = "white",bg = "white")
    text(0.47,0.0148, labels = "Rainfall (mm)", cex=1.5)
  }

}
mtext('Cumulative Residual Variance',outer = TRUE, cex= 1.75, side =4, line=3)
#dev.off()








###PLotting only residual plots--------------------------------------
#######LOESS function Residuals#########
pdf(paste("/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/figures/resid_v3.pdf"),
    height=7.5,width=13)
par(oma = c(2,3,0,4))
par(mfrow = c(2,2),mar = c(2.5,3.5,1,2),mgp = c(1, 0.7, 0), tck = -0.015)
for(i in seq_along(stations.list)){
  
  #resid lo fit
  lo<-loess(stations.list[[i]]$VWC40 ~ stations.list[[i]]$VWC5, span=0.75,degree=2,se=TRUE)### span can be changed for smoothness
  lo.res<-resid(lo)
  ###Function for rounding to different sequences##
  mround <- function(x,base){
    base*round(x/base)}
  #cond mean +var+ bars
  mean.lores<-aggregate(lo.res,by=list(mround(lo$x,0.01)),mean,na.rm=TRUE)
  #var
  var.lores<-aggregate(lo.res,by=list(mround(lo$x,0.01)),var,na.rm=TRUE)
  var.lores[is.na(var.lores)]<-0
  var.lores$cum<-cumsum(var.lores$x)
  x<-as.matrix(var.lores$cum)
  infl<-e.divisive(diff(x), sig.lvl=.05,k=NULL,min.size=2,alpha = 1, R=199)
  
  plot(lo$x,lo.res, main = "",pch=21, bg=alpha(colB[stations.list[[i]]$mon],0.7),
       col=alpha("grey38",0.25), ylab="", xlab="", yaxt="n",
       xaxt="n", cex=ceiling(stations.list[[i]]$RDmm/7.5), lwd=0.5, ylim=c(-0.2, 0.2),  xlim = c(0,0.5))
  abline(0,0, lwd=0.5)
  
  if (i %in% c(1,3)){
    axis(2, cex.axis = 1.5, las =2)
  }
  if (i %in% c(2,4)){
    axis(side=2,labels=F) 
  }
  
  if (i %in% c(3,4)){
    axis(1, cex.axis = 1.5,at=seq(0,0.5,by=0.1))
  }
  if (i %in% c(1,2)){
    axis(side=1,labels=F) 
  }
  
  par(new = T)
  plot(var.lores$Group.1,var.lores$cum,type="l",lwd=1.25,xlab="",ylab="",
       xaxt="n",yaxt = "n",ylim = c(0, max(var.lores$cum)),xlim = c(0,0.5))
  points(var.lores$Group.1,var.lores$cum,pch = 15, cex = 0.5)
  abline(0,0, lwd = 0.5)
  axis(side=4, cex.axis = 1.5, las =2)
  lines(var.lores$Group.1,var.lores$x, type = "h", lwd=3.5,col="grey38")
  #add ecp estimates
  inf<-infl$estimates
  inf<-inf[-c(1, length(inf))]
  abline(v =var.lores$Group.1[c(inf)], col = "grey38", lty = 3, lwd = 2)
  inf.list<-var.lores$Group.1[c(inf)]
  text(x=var.lores$Group.1[c(inf)]+0.015, y=max(var.lores$cum)*0.975,expression(italic(theta[c])), cex = 1.5)
}
mtext('Cumulative Residual Variance',outer = TRUE, cex= 1.75, side =4, line=2.75)
mtext('Residuals',outer = TRUE, cex= 1.75, side =2, line=0)
mtext(expression("Soil Moisture at 5cm " ~(cm^3~cm^{-3})), side = 1, outer = TRUE, cex=1.75, line=1.2)
dev.off()
