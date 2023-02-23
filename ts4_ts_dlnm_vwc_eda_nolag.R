#CHOOSE DIRECTORY
setwd("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/2014-2016")
data<-read.csv("itc20_win_day_dec16_sub.csv",   header =T, stringsAsFactors = FALSE)


#PACKAGE load
require(dlnm); require(splines);require(zoo); require(xts);require(lubridate);library(mgcv) ; library(splines) ; library(tsModel);
require("graphics");library(RColorBrewer);library(scales); library(ecp);
library(plyr)


#CHECK VERSION DLNM
if(packageVersion("dlnm")<"2.2.0")
  stop("update dlnm package to version >= 2.2.0")

###1.0 SETTING UP DATASETS
#ADD DATE COLUMNS, FIX DATE FORMATE
data$date<-as.POSIXct(data$date,format = "%Y-%m-%d %H:%M:%S")
#data <-data[data$date<="2016-10-18 23:45:00", ]
data$year<-year(data$date)
data$mon <-month(data$date)
data$doy<-yday(data$date)

##2.0 CLEANING DATA FOR TIME SERIES
#IMPUTE NAs  -multiple imputation
library(imputeTS)
library(Amelia)
library(missForest)

#bds <- matrix(c(2, 0.2, 0.3), nrow = 1, ncol = 3)
#amelia_fit <- amelia(data[,-c(10:11)], m=5, parallel = "multicore",ts = "date",polytime = 1,
empri = .07*nrow(data), autopri = 0.5,boot.type= "ordinary",
pplinetime = 0,noms = "year", bounds = bds, max.resample = 20)


#x<-amelia_fit$imputations[[2]]
#plot(x$date,x$VWC40, type = "l")

###convert to % values###
data$VWC5p <- data$VWC5*100
data$VWC40p <- data$VWC40*100
data$VWC80p <- data$VWC80*100
summary(data)
plot(data$date, data$VWC40p, type = "l")
options(na.action="na.exclude")


#CALCULATE ZONE AVERAGE SOIL MOISTURE (Gretchen Miller et al., 2007)
data<- within(data,VWC_zone <-(VWC5*7.5+VWC10*7.5 +VWC20*15 +VWC40*30)/60)
data<- within(data,VWC_zonep <- VWC_zone*100)
plot(data$date, data$VWC5, type="l")

###RAINFAL KNMI##
#TWENTHE_670
#REKKEN_674
#LOCHEM_663
#HOLTEN_687
#KNMI RAINFALL DOWNLOAD- change to closest KNMI station from SM station
setwd("/media/coleen/DDrive1/1_OWASIS/0_ExternalDATA/KNMI")
temp <- tempfile()
download.file("http://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/monv_reeksen/neerslaggeg_HOLTEN_687.zip",temp,mode="wb")
precip<- read.table(unz(temp,"neerslaggeg_HOLTEN_687.txt"), sep=",",skip=23, header=T)[ ,1:3]
head(precip)
precip$RDmm <-precip$RD/10
precip$YYYYMMDD<- strptime(precip$YYYYMMDD, "%Y%m%d")
##choice of 06:00:00 is abitrary - only so I can convert it to POSIXct
precip$date <- paste(precip$YYYYMMDD,"23:45:00")
precip$date <- as.POSIXct(precip$date,format = "%Y-%m-%d %H:%M:%S")

#count rainfall days and dry days
hist(precip$RDmm)
precip$rainy <-ifelse(precip$RDmm >=5,1,0)
precip$dry<-ifelse(precip$RDmm ==0,1,0)
precip$rainy.days<-(precip$rainy) * unlist(lapply(rle(precip$rainy)$lengths, seq_len))
precip$dry.days<-(precip$dry) * unlist(lapply(rle(precip$dry)$lengths, seq_len))


#MERGE SOIL MOISTURE AND RAINFALL TO ONE DATAFRAME
datetime1 <- strptime(data$date,format = "%Y-%m-%d %H:%M:%S")
datetime2 <- strptime(precip$date, format= "%Y-%m-%d %H:%M:%S")

data<-merge(  cbind(round(datetime1, "days"), data),
              cbind(round(datetime2, "days"), precip),
              by=1, suffixes=c("_dts1", "_dts2"))[-1]

head(data)
data.ts<-data[,c(2,3,8,10:14,20)]
#dec16 
data.ts<-data[,c(2,3,8,12:16,23)]
ts <-xts(data.ts,data$YYYYMMDD)
plot(as.zoo(ts[,c(2,1,3,9)]), col='black')

###IMPUTE MISSING VALUES##
#data.ts<-data
data.ts$VWC5<-na.ma(data.ts$VWC5, k=6, weighting = "linear")
data.ts$VWC10 <-na.seadec(data.ts$VWC10,algorithm = "ma")
data.ts$VWC20 <-na.seadec(data.ts$VWC20,algorithm = "ma")
data.ts$VWC40 <-na.seadec(data.ts$VWC40,algorithm = "ma")
data.ts$Temp5 <-na.seadec(data.ts$Temp5,algorithm = "ma")
ts <-xts(data.ts,data$YYYYMMDD)
plot(as.zoo(ts[,c(2,1,3,9)]), col='blue')


#ADDRESS SEASONALITY + TREND
library(forecast)
temp.ts<-ts(data.ts$Temp5,frequency = 365)
Temp.stl <- stl(temp.ts, s.window = "periodic", robust = FALSE)
plot(Temp.stl)
Temp.sa <- seasadj(Temp.stl)  # de-seasonalize
plot(Temp.sa)

VWC5.ts<-ts(data.ts$VWC5,frequency = 365)
VWC5.stl <- stl(VWC5.ts, s.window = "periodic", robust = FALSE)
plot(VWC5.stl)
VWC5.sa <- seasadj(VWC5.stl)  # de-seasonalize
plot(VWC5.sa)

VWC40.ts<-ts(data.ts$VWC40,frequency = 365)
VWC40.stl <- stl(VWC40.ts, s.window = "periodic", robust = FALSE)
plot(VWC40.stl)
VWC40.sa <- seasadj(VWC40.stl)  # de-seasonalize
plot(VWC40.sa)

#MERGE DECOMPOSED DATA INTO A DATAFRAME
temp5<- as.numeric(Temp.sa)
VWC5<- as.numeric(VWC5.sa)
VWC40<- as.numeric(VWC40.sa)
data.stl <-data.frame(cbind(temp5,VWC5,VWC40, data.ts[,c(4:6,9)]))
data.stl$VWC5p<-data.stl$VWC5*100
data.stl$VWC40p<-data.stl$VWC40*100
data.stl$YYYYMMDD<-data$YYYYMMDD

#########GETTING first >DIFFERENCE
dVWC5<-diff(data.ts$VWC5,1)
dVWC40<-diff(data.ts$VWC40,1)
dTemp5<-diff(data.ts$Temp5,1)
data.diff <-data.frame(cbind(dTemp5,dVWC5,dVWC40, data.ts[-1,c(4:6,9)]))

#PLOTTING TIME SERIES
#XTS SERIES - ORIGINAL DATA
require(xts)
ots.xts<-xts(data.ts, data$YYYYMMDD)
plot(ots.xts)

#pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/ITC05_TS_rev.pdf"),
#    height=4,width =12)
par(mfrow = c(1,1), mar = c(4,4,2,4) + 0.01,mgp=c(1.5,0.35,0), tck=-0.015)
myColors <- c("skyblue1", "dodgerblue4", "grey30", "darkblue", "darkviolet")
#Plot rainfall first so Soil moisture TS is more visible
plot(as.zoo(ots.xts[,9]),type="h",lwd =2.5,  col="snow3",yaxt="n",xlab="",ylab="", xaxt = "n",
     xlim= c(as.POSIXct('2014-01-01 23:45:00', format="%Y-%m-%d %H:%M:%S"),as.POSIXct('2016-12-31 23:45:00', format="%Y-%m-%d %H:%M:%S")),
     ylim = rev(c(0,80)))
axis(side=4,cex.axis = 1.4,mgp=c(0.5,0.35,0))
mtext(side = 4, line = 2, 'Rainfall (mm)', cex  = 1.3)
legend(x="bottomleft",legend = c("5cm","40cm", "rainfall"), y.intersp = 0.75, x.intersp = 0.25,bty = "n",
       seg.len =1.25, lwd = c(4,4,4.5),lty = 1, col = myColors)
par(new = T)
plot(as.zoo(ots.xts[,c(1,3)]),cex.lab = 1.3, xlab = "Time", ylab =expression("Soil Moisture" ~ (cm^3~cm^{-3})), 
     xlim= c(as.POSIXct('2014-01-01 23:45:00', format="%Y-%m-%d %H:%M:%S"),as.POSIXct('2016-12-31 23:45:00', format="%Y-%m-%d %H:%M:%S")),
     main = "SM09",cex.axis = 1.4, col = myColors, lwd=2,screens=1, ylim = c(0,0.60))
#dev.off()


#XTS SERIES - DECOMPOSED DATA
ts.xts<-xts(data.stl, data$YYYYMMDD)
#Plot rainfall first so Soil moisture TS is more visible
plot(as.zoo(ts.xts[,7]),type="h",lwd =1,  col="lightgrey",yaxt="n",xlab="",ylab="")
axis(side=4)
mtext(side = 4, line = 2, 'Rainfall (mm)',cex = 1.2)
legend(x="topleft",legend = c("5cm","40cm", "rainfall"), y.intersp = 0.5, x.intersp = 0.25,bty = "n",
       seg.len = 0.5, lwd = 2,lty = 1, col = myColors)
par(new = T)
plot(as.zoo(ts.xts[,2:3]),cex.lab = 1.2, xlab = "Time", ylab = expression("Soil Moisture " (cm^3/cm^3)), 
     main = "SM13 - decomposed", col = myColors, lwd=1.5,screens=1, ylim = c(0,0.60))


##cross correlation data - ORIGINAL
#pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/ITC20_ccf.pdf"),
#    height=3,width=3)
par(mfrow= c(1,1),mar = c(2,2,1,1),mgp = c(1,0.15,0), tck = -0.009)
ccf<-ccf(as.ts(data.ts$VWC5),as.ts(data.ts$VWC40), lag.max = 10,ylim = c(0.6,0.85),cex.lab= 0.8, cex.axis=0.7,
         main="", xlab = "Lag (Day)", ylab = "Correlation",
         col = "grey60",lwd=5.5)
title( main = "SM20", cex.main = 1)
abline(v=0, col = "blue",lty= 5, lwd=1)
#dev.off()

##cross correlation data - DeSeason
#pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/ITC20_ccf_deseas.pdf"),
#height=7,width = 7)
par(mfrow= c(1,1),mar = c(3,3,1.5,1),mgp = c(1.5,0.3,0), tck = -0.005)
ccf<-ccf(as.ts(data.stl$VWC5),as.ts(data.stl$VWC40), lag.max = 30,ylim = c(0,0.85),cex.lab= 1.3, cex.axis=1.3, cex = 1.5,
         main="", xlab = "Lag in Days", ylab = "Correlation",
         col = "grey40",lwd=6)
title( main = "SM05", cex.main = 1.5)
abline(v=0, col = "blue",lty= 4, lwd = 1.5)
#dev.off()




##WARNING-check if export functions will be executed - this will change you files in DDrive!!!
##3.0 EXPLORATORY DATA ANALYSIS
#Original Series
#HISTORGRAMS
hist(data$VWC5, xlim = c(0,0.5), main = "ITC Sm05 VWC@5cm",col = "grey", xlab="VWC 5cm", cex.main =2.25, cex.lab = 1.5, cex.axis=1.5)
hist(data$VWC40, xlim = c(0,0.5),main = "ITC Sm05 VWC@40cm",col = "grey", xlab="VWC 40cm", cex.main =2.25, cex.lab = 1.5, cex.axis=1.5)

##BOXPLOTS
box<-data.stl[,c(2,3)]
box$group<-round(box$VWC5,2)
boxplot(VWC40~group,data=box)
y<-ddply(box, "group", summarise, corr=cor(VWC5, VWC40,use="pairwise.complete.obs"))

mround <- function(x,base){
  base*round(x/base)
}


#######DECOMPOSED SERIES##########
#HISTORGRAMS
hist(data.stl$VWC5, xlim = c(0,0.5), main = "ITC SM20 VWC@5cm",col = "grey", xlab="VWC 5cm", cex.main =2.5, cex.lab = 2, cex.axis=2)
hist(data.stl$VWC40, xlim = c(0,0.5),main = "ITC SM20 VWC@40cm",col = "grey", xlab="VWC 40cm", cex.main =2.5, cex.lab = 2, cex.axis=2)
#CREATE LAG DATA

lag.df<-cbind(lag_VWC5,data.stl[,c(2:7)])

lag.df<-data.stl[,c(2:7)]

#############ORIGINAL DATA############
#CREATE LAG DATA
lag.df<-data[,c(2,8,12:14,23)]

#display.brewer.all() 
#LINEAR AND LOESS REGRESSION
rdm<-max(log(lag.df$RDmm))
colB<-palette(colorRampPalette(brewer.pal(11,"Spectral"))(12))
cm<-
#reg fit
  lo<- loess(lag.df$VWC40 ~ lag.df$VWC5, span=0.75,degree=2,se=TRUE)### span can be changed for smoothness
  lm<-lm(lag.df$VWC40~lag.df$VWC5)
  #exporting to file
  #pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/ITC05_vwc5-40_reg1.pdf", sep = ""),
 # height=7,width=6)
  #type= "cairo",res=300, width = 8*300, height=9*300, units = "px", pointsize = 12) ## height =1014 for SM20, other height =865
  # blank plot- ylim for SM20
  par(mfrow = c(1,1), mar = c(4.5,4.5,2,2) + 0.0,mgp=c(2.25,0.6,0), tck = -0.01)
  plot(c(0,0.5), c(0,0.6), type="n", xlab=expression("Soil Moisture at 5cm " ~ (cm^3~cm^{-3})), cex.main= 1.5,cex.lab = 1.4,
       cex.axis=1.4,ylab=expression("Soil Moisture at 40cm "~(cm^3~cm^{-3})), main= bquote("SM05 Lag"== 0))
  
  #SM points
  points(lag.df$VWC40~lag.df$VWC5,pch=21, bg=alpha(colB[lag.df$mon],0.75),
         col=alpha("black",0.55),cex=ceiling(lag.df$RDmm/10), lwd=0.5)
  
  #fit loess line
  j <- order(lo$x)
  lines(lo$x[j],lo$fitted[j], lwd=2.5)
  
  #cond mean +var+ bars
  mean.agg<-aggregate(lag.df$VWC40,by=list(round(lag.df$VWC5,2)),mean,na.rm=TRUE)
  var.agg<-aggregate(lag.df$VWC40,by=list(round(lag.df$VWC5,2)),var,na.rm=TRUE)
  sd.agg<-aggregate(lag.df$VWC40,by=list(round(lag.df$VWC5,2)),sd,na.rm=TRUE)
  sd.plus<- mean.agg$x + sd.agg$x
  sd.minus<-mean.agg$x - sd.agg$x
  #lm fit line
  lines(lm$model$`lag.df$VWC5`,lm$fitted.values,col = "red3", lwd = 2, lty=3)
  
  #correlation grouped VWC values
  x<-data.frame(cbind(lag.df$VWC5,lag.df$VWC40))
  x$x<-round(x$X1/.05)*.05 #round(x$X1,2)
  y<-ddply(x, "x", summarize, corr=cor(X1, X2,method = "spearman",use="pairwise.complete.obs"))#use="pairwise.complete.obs"
  frequ<-data.frame(count(round(lag.df$VWC5/0.05)*0.05))
  corr.list<-merge(frequ,y, by="x")
  
  arrows(mean.agg$Group.1, sd.minus, mean.agg$Group.1, sd.plus, lwd=1,col=alpha("black",0.7),
         length=0.05, angle=90, code=3)
  points(mean.agg, col=alpha("black",0.7), pch =15)
  
  
  #legend
  legend(x=0.4,y=0.095,legend = month.abb[1:12], y.intersp = 0.75, x.intersp = 0.5,bty = "n",
         pt.cex = 1, cex=0.95, col = alpha("black",0.55), pt.bg =alpha(colB,0.75), pch=21, ncol=2)
  legend(x=0,y=0.05,legend = c("linear fit", "loess fit","conditional mean +/- sd"),lty =c(6,1,1),cex=0.95,
         lwd=c(2,3,1), seg.len=1.5,x.intersp = 0.75,y.intersp = 0.75,bty = "n", col=c("red3","black",alpha("black",0.7)))
  f<-unique(floor(lag.df$RDmm/10))
  c<-sort(unique(ceiling(lag.df$RDmm/10)),decreasing = FALSE)[-1]
  ##spacing for rainfall should be adjusted every time! - remove  some numbers from sequence to avoid duplicate points)
  ## for SM09 - x.intersp=c(1.85,1.85,1.85,1.85,1.85,2.7),y.intersp = c(1,1,1.2,1.4,1.6,1.9)
  ##for other 3 stations -- x.intersp=c(1.85,1.85,1.85,1.85),y.intersp = c(1,1,1.2,1.4)
  legend(x=0.43,y=0.26,legend=c*10,pch=21,pt.bg=alpha("grey",0.6),col = "black",x.intersp=c(1.85,1.85,1.85,1.85),y.intersp = c(1,1,1.2,1.4),
         pt.cex=c, bty="n", cex=0.95)
  
  #add text
  cor<-round(cor(lag.df$VWC5,lag.df$VWC40, method = "spearman",use="pairwise.complete.obs"),digits=3)
  sigma.lm<-round(sqrt(sum(lm$residuals^2)/lm$df.residual),3)
  sigma.lo<-round(sqrt(sum(lo$residuals^2)/lo$n),3)
  points(0.7,0.4, pch="")
  text(0.05,0.60,labels = bquote('R'[s]==.(cor)), cex = 1, col = "grey38")
  text(0.05,0.58,labels = bquote('RSE'[italic(lm)]==.(sigma.lm)), cex = 1, col = "grey38")
  text(0.05,0.56,labels = bquote('RSE'[italic(lo)]==.(sigma.lo)), cex = 1, col = "grey38")
  lm_coef <- round(coef(lm), 3) # extract coefficients lm
  #text(0.07,0.375,labels=bquote(y(lm) == .(lm_coef[2])*x + .(lm_coef[1])), cex=1.25)
  text(0.45,0.26, labels = "Rainfall (mm)", cex = 1)
  #dev.off()
  
  #resid lm and lo fit
  lm.res<-resid(lm)
  lo.res<-resid(lo)
  count<-count(round(lag.df$VWC5,2))
  
  ###Function for rounding to different sequences##
  mround <- function(x,base){
    base*round(x/base)
  
  
  #LMRES
  #cond mean +var+ bars
  mean.lmres<-aggregate(lm.res,by=list(mround(lag.df$VWC5,0.01)),mean,na.rm=TRUE)
  #var
  var.lmres<-aggregate(lm.res,by=list(mround(lag.df$VWC5,0.01)),var,na.rm=TRUE)
  var.lmres[is.na(var.lmres)]<-0
  var.lmres$cum<-cumsum(var.lmres$x)
  x<-as.matrix(var.lmres$cum)
  infl<-e.divisive(diff(x),sig.lvl=.01,k=NULL,min.size=2,alpha = 1, R=199)
  #sd
  sd.lmres<-aggregate(lm.res,by=list(mround(lag.df$VWC5,0.01)),sd,na.rm=TRUE)
  sd.lmres.plus<- mean.agg$x + sd.agg$x
  sd.lmres.minus<-mean.agg$x - sd.agg$x
  
  #exporting
  #pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/ITC20_vwc5-40reslm.pdf", sep = ""),
  #height=6,width=9)
  #type= "cairo",res=300, width =8*300, height=6*300, units = "px", pointsize = 12)
  #plot resid lm
  par(mfrow = c(1,1), mar = c(4.5,3.8,3.5,4) + 0.0,mgp=c(2.5,0.7,0), tck = -0.01)
  plot(lag.df$VWC5,lm.res, main = bquote("SM20 Linear fit residuals Lag"  == 0),pch=21, bg=alpha(colB[lag.df$mon],0.7),
       col=alpha("black",0.55), ylab="Residuals", xlab=expression("Soil Moisture at 5cm"~(cm^3~cm^{-3})),cex.main= 2,cex.lab = 1.6,
       cex.axis=1.4, cex=ceiling(lag.df$RDmm/10), lwd=0.5, ylim=c(-0.2, 0.2), xlim = mround(range(lag.df$VWC5, na.rm=TRUE),0.01))
  abline(0,0)
  par(new = T)
  plot(var.lmres$Group.1,var.lmres$cum,type="l",lwd =3,cex.axis=1.4,  col="darkgrey",yaxt="n",xlab="",ylab="", ylim = c(0, max(var.lmres$cum)))
  lines(var.lmres$Group.1,var.lmres$x, type = "h", lwd=2.5)
  #add ecp estimates
  inf<-infl$estimates
  abline(v =var.lmres$Group.1[c(inf)], col = "grey38", lty = 3, lwd = 2)
  axis(side=4, cex.axis = 1.4)
  mtext(side = 4, line = 2, 'Cumulative Residual Variance',cex= 1.4)
  inf.list<-var.lmres$Group.1[c(inf)]
  # dev.off()
  
  
  #LORES
  #cond mean +var+ bars
  mean.lores<-aggregate(lo.res,by=list(mround(lag.df$VWC5,0.01)),mean,na.rm=TRUE)
  #var
  var.lores<-aggregate(lo.res,by=list(mround(lag.df$VWC5,0.01)),var,na.rm=TRUE)
  var.lores[is.na(var.lores)]<-0
  var.lores$cum<-cumsum(var.lores$x)
  x<-as.matrix(var.lores$cum)
  infl<-e.divisive(diff(x), sig.lvl=.01,k=NULL,min.size=2,alpha = 1, R=199)
  #sd
  sd.lores<-aggregate(lo.res,by=list(mround(lag.df$VWC5,0.01)),sd,na.rm=TRUE)
  sd.lores.plus<- mean.agg$x + sd.agg$x
  sd.lores.minus<-mean.agg$x - sd.agg$x
  
  #exporting
  #pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/ITC05_vwc5-40reslo_v2.pdf", sep = ""),
  # height=7,width=10)
  #type= "cairo",res=300, width = 8*300, height=6*300, units = "px", pointsize = 12)
  #plot resid lm
  par(mfrow = c(1,1), mar = c(3.7,3.5,3,3.3) + 0.0,mgp=c(2.2,0.5,0), tck = -0.01)
  plot(lag.df$VWC5,lo.res, main = bquote("SM05 Loess fit residuals Lag"  == 0),pch=21, bg=alpha(colB[lag.df$mon],0.7),
       col=alpha("grey38",0.5), ylab="Residuals", xlab=expression("Soil Moisture at 5cm"~(cm^3~cm^{-3})),cex.main= 1.5,cex.lab = 1.4,
       cex.axis=1.4, cex=ceiling(lag.df$RDmm/10), lwd=0.5, ylim=c(-0.2, 0.2), xlim = mround(range(lag.df$VWC5, na.rm=TRUE),0.01))
  abline(0,0)
  par(new = T)
  plot(var.lores$Group.1,var.lores$cum,type="l",lwd=1.25,cex.axis=1.4,yaxt="n",xlab="",ylab="", ylim = c(0, max(var.lores$cum)))
  points(var.lores$Group.1,var.lores$cum,pch = 15, cex = 0.5)
  abline(0,0, lwd = 0.5)
  axis(side=4, cex.axis = 1.4)
  lines(var.lores$Group.1,var.lores$x, type = "h", lwd=3.5,col="grey38")
  #add ecp estimates
  inf<-infl$estimates
  abline(v =var.lores$Group.1[c(inf)], col = "grey38", lty = 3, lwd = 2)
  
  mtext(side = 4, line = 2, 'Cumulative Residual Variance',cex= 1.4)
  inf.list<-var.lores$Group.1[c(inf)]
  #dev.off()
  
  }
  
######PLOT counts vs Variance#######
colnames(count)<-c("Interval","Samples")
var.lores<-var.lores[,1:2]
colnames(var.lores)<-c("Interval","Var")
var.lores$Var<-var.lores$Var
var.count05<-merge(count,var.lores, by = "Interval")### put number with station
m <-matrix(var.count05$Var)
var.count05$Var<-apply(m,MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

var.count05$Station = 05
min.var<- min(var.count05$Interval)
max.var<-max(var.count05$Interval)
var.count05$Interval.norm<- (var.count05$Interval - min.var)/(max.var-min.var)

min.var<- min(var.count05$Var)
max.var<-max(var.count05$Var)
var.count05$Var.norm<- (var.count05$Var - min.var)/(max.var-min.var)
#######merge all four data frames#####
var.count.all <-rbind(var.count05,var.count09,var.count13,var.count20)

##save var.count.all
write.csv(var.count.all, "/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/figures/var_count_all.csv",row.names = FALSE)

##read saved var.count.all if applicable
var.count.all<-read.csv("/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/figures/var_count_all.csv",
                        header =T, stringsAsFactors = FALSE)
###PLOT 2 layout column
pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/var.samples.all", ".pdf", sep = ""),
    height=5,width=7)


layout(matrix(1:2,ncol=2), widths = c(11,1),height = c(1,1))

##First plot
par(mar=c(3,3,2,0), mgp = c(1.6,0.4,0),tck = -0.01)  
colC<-palette(colorRampPalette(brewer.pal(9,"Spectral"))(100))
plot(Samples~Var.norm,data = var.count.all, pch = c(21,22,23,24)[as.factor(Station)], ylab = "Number of samples",xlab = "Residuals Variance ", cex.lab =1, cex.axis=1,
      cex.main=1.5,main = "Samples vs. Residuals Variance",col = alpha("black",0.75),bg=alpha(colC[(Interval.norm*100)], 0.75), cex =2.7)
legend(x=0.85,y=118,legend=paste("SM",sort(unique(var.count.all$Station),decreasing = FALSE)),pch=c(21,22,23,24),
       y.intersp=1.05,x.intersp = 0.75,bty = "n", pt.cex=1.8,cex=1, pt.bg =  alpha("black",0.55))

##Second plot
par(mar=c(11,0.5,5,0.8)) 
legend_image <- as.raster(alpha(colC,0.65))
plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=0.7, y = seq(0.01,0.99,l=2), labels = c("WET","DRY"), cex=0.56)
rasterImage(legend_image, 0.0,0.0,0.45,1)

dev.off()

#######################Version 2 variance calculation#################NOT USED
######PLOT counts vs Variance#######
colnames(count)<-c("Interval","Samples")
var.lores<-var.lores[,1:2]
colnames(var.lores)<-c("Interval","Var")
var.count20<-merge(count,var.lores, by.y  = "Interval",all=TRUE)### put number with station
var.count20<-var.count20[complete.cases(var.count20$Interval),]
var.count20$Station = 20
var.count.all <-rbind(var.count05,var.count09,var.count13,var.count20)

##save var.count.all
write.csv(var.count.all, "/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/figures/var_count_all_v2.csv",row.names = FALSE)
##read saved var.count.all if applicable
var.count.all<-read.csv("/media/coleen/DDrive1/1_OWASIS/4_Papers/ts_dnlm/figures/var_count_all.csv",
                        header =T, stringsAsFactors = FALSE)
max.int<-max(var.count.all$Interval)
min.int<-min(var.count.all$Interval)
var.count.all$Interval.norm<-round(((var.count.all$Interval-min.int)/(max.int-min.int))*max.int,2)
plot(var.count.all$Interval,var.count.all$Interval.norm)  

max.var<-max(var.count.all$Var)
min.var<-min(var.count.all$Var)
var.count.all$Var.norm<-round(((var.count.all$Var-min.var)/(max.var-min.var))*max.var,5)
var.count.all$Var.norm<-((var.count.all$Var-min.var)/(max.var-min.var))

###PLOT 2 layout column
#pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/var.samples.all", ".pdf", sep = ""),
#    height=7,width=11)

layout(matrix(1:2,ncol=2), widths = c(11,1),height = c(1,1))
Interval.norm<-var.count.all$Interval.norm
##First plot
par(mar=c(5,5,4,1), mgp = c(2.5,0.6,0))  
colC<-palette(colorRampPalette(brewer.pal(11,"Spectral"))(50))
plot(Samples~Var.norm,data = var.count.all, pch = c(21,22,23,24)[as.factor(Station)], ylab = "Number of samples",xlab = "Variance of Residuals", cex.lab = 1.5, cex.axis=1.5,
     cex.main=1.5,main = "Samples vs. Variance of Residuals",col = alpha("black",0.75),
     bg=alpha(colC[(var.count.all$Interval.norm*100)], 0.75), cex = 3.5)
legend(x=0.04,y=118,legend=paste("SM",sort(unique(var.count.all$Station),decreasing = FALSE)),pch=c(21,22,23,24),
       y.intersp=1,x.intersp = 0.75,bty = "n", pt.cex=2,cex=1.2, pt.bg =  alpha("black",0.55))

##Second plot
par(mar=c(10,0.2,11,1)) 
legend_image <- as.raster(alpha(colC,0.75))
plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=0.65, y = seq(0.01,0.99,l=2), labels = c("DRY","WET"), cex=0.7)
rasterImage(legend_image, 0,0,0.4,1)

#dev.off()