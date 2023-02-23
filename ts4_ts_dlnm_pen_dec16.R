#CHOOSE DIRECTORY
setwd("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/2014-2016")
data<-read.csv("itc13_win_day_dec16_sub.csv",   header =T, stringsAsFactors = FALSE)## add _sub for SM13 and SM20
gw<-read.csv("/media/coleen/DDrive1/1_OWASIS/0_ExternalDATA/GW_Dino/B34A0111001_1_SM13_gwsub.csv",   header =T, stringsAsFactors = FALSE)
#B28C0429001_1_SM20_gwsub
#B28G0409001_1_SM05_gwsub
#B34A0111001_1_SM13_gwsub
#B34H0057001_1_SM09_gwsub

#PACKAGE load
require(dlnm); require(splines);require(zoo); require(xts);require(lubridate);
library(mgcv) ; library(splines) ; library(tsModel);
require("graphics");library(RColorBrewer);library(scales); library(ecp);
library(plyr);library(imputeTS)
source("findmin.R")

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
download.file("http://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/monv_reeksen/neerslaggeg_LOCHEM_663.zip",temp,mode="wb")
precip<- read.table(unz(temp,"neerslaggeg_LOCHEM_663.txt"), sep=",",skip=23, header=T)[ ,1:3]
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

data.ts<-data[,c(2,3,8,12:16,23)]
data.ts$gw <-gw[,5] #-c(1:114) for rows of GW13 only
ts <-xts(data.ts,data$YYYYMMDD)
plot(as.zoo(ts[,c(2,1,3,9)]), col='black')


##2.0 CLEANING DATA FOR TIME SERIES
#IMPUTE LOCHEM_663NAs 
library(imputeTS)
data.ts$VWC5<-na.seadec(data.ts$VWC5, algorithm = "ma")
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
data.stl <-data.frame(cbind(temp5,VWC5,VWC40, data.ts[,c(4:6,9,10)]))
data.stl$VWC5p<-data.stl$VWC5*100
data.stl$VWC40p<-data.stl$VWC40*100


#PLOTTING TIME SERIES
#XTS SERIES - ORIGINAL DATA
require(xts)
ots.xts<-xts(data.ts, data$YYYYMMDD)
plot(ots.xts)

#pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/ITC05_TS.pdf"),
#   height=4,width =12)
par(mfrow = c(1,1), mar = c(4,4,2,4) + 0.01,mgp=c(1.5,0.35,0), tck=-0.015)
myColors <- c("lightblue2", "royalblue2", "grey", "darkblue", "darkviolet")
#Plot rainfall first so Soil moisture TS is more visible
plot(as.zoo(ots.xts[,9]),type="h",lwd =1.5,  col="grey",yaxt="n",xlab="",ylab="", xaxt = "n",
     xlim= c(as.POSIXct('2014-01-01 23:45:00', format="%Y-%m-%d %H:%M:%S"),as.POSIXct('2016-12-31 23:45:00', format="%Y-%m-%d %H:%M:%S")),
     ylim = rev(c(0,80)))
axis(side=4,cex.axis = 1.4,mgp=c(0.5,0.35,0))
mtext(side = 4, line = 2, 'Rainfall (mm)', cex  = 1.3)
legend(x="tbottomleft",legend = c("5cm","40cm", "rainfall"), y.intersp = 0.75, x.intersp = 0.25,bty = "n",
       seg.len = 1.25, lwd = c(2,2,1),lty = 1, col = myColors)
par(new = T)
plot(as.zoo(ots.xts[,c(1,3)]),cex.lab = 1.3, xlab = "Time", ylab =expression("Soil Moisture " (cm^3/cm^3)), 
     xlim= c(as.POSIXct('2014-01-01 23:45:00', format="%Y-%m-%d %H:%M:%S"),as.POSIXct('2016-12-31 23:45:00', format="%Y-%m-%d %H:%M:%S")),
     main = "SM05",cex.axis = 1.4, col = myColors, lwd=1.5,screens=1, ylim = c(0,0.60))
#dev.off()

#XTS SERIES - Seasonally adjusted DATA stl
ts.xts<-xts(data.stl, data$YYYYMMDD)

#Plot rainfall first so Soil moisture TS is more visible
plot(as.zoo(ts.xts[,7]),type="h",lwd =1,  col="lightgrey",yaxt="n",xlab="",ylab="",
     xlim= c(as.POSIXct('2014-01-01 23:45:00', format="%Y-%m-%d %H:%M:%S"),as.POSIXct('2016-12-31 23:45:00', format="%Y-%m-%d %H:%M:%S")))
axis(side=4)
mtext(side = 4, line = 2, 'Rainfall (mm)',cex = 1.4)
legend("topleft",legend = c("5cm","40cm", "rainfall"), y.intersp = 0.25, x.intersp = 0.25,bty = "n",
       seg.len = 0.35, lwd = 3,lty = 1, col = myColors, cex = 1.25)
par(new = T)
plot(as.zoo(ts.xts[,2:3]),cex.lab = 1.4, xlab = "Time", ylab = "Soil moisture [cm3/cm3]", 
     main =  "Station #05", col = myColors, lwd=1.5,screens=1, ylim = c(0,0.50),
     xlim= c(as.POSIXct('2014-01-01 23:45:00', format="%Y-%m-%d %H:%M:%S"),as.POSIXct('2016-12-31 23:45:00', format="%Y-%m-%d %H:%M:%S")))


##cross correlation seasonally adjusted data
par(mfrow= c(1,1),mgp = c(2.1,0.5,0), tck = -0.01)
ccf<-ccf(as.ts(data.ts$VWC5),as.ts(data.ts$VWC40), lag.max = 120)
abline(v=0, col = "blue",lty= 2)

##SPearman Rank cor then ccf###
VWC5.rank<-order(data.stl$VWC5)
VWC40.rank<-order(data.stl$VWC40)

ccf<-ccf(as.ts(data.ts$VWC5),as.ts(data.ts$VWC40), lag.max = 30,ylim = c(0.6,0.85),cex.lab= 1.3, cex.axis=1.3, cex = 1.5,
         main="", xlab = "Lag in Days", ylab = "Correlation",
         col = "grey40",lwd =6)
title( main = "SM20", cex.main = 1.5)
abline(v=0, col = "blue",lty= 2, lwd = 2)

###CCF###
ccf<-ccf(VWC5.rank,VWC40.rank, lag.max = 30,ylim = c(0,0.85),cex.lab= 1.3, cex.axis=1.3, cex = 1.5,
         main="", xlab = "Lag in Days", ylab = "Correlation",
         col = "grey40",lwd=6)
title( main = "SM05", cex.main = 1.5)
abline(v=0, col = "blue",lty= 4, lwd = 1.5)


#######5.2 SEASONAL ADJUSTMENT DATA########
##GAM WITH DOUBLY VARYING PENALTY ON THE LAG -internal GAM
# DEFINE MATRICES TO BE INCLUDED AS TERMS IN THE SMOOTHER
Q <- Lag (ts.xts$VWC5p,0:30)
R <- Lag (ts.xts$RDmm,0:30)
S <- Lag (ts.xts$gw,0:30)
L <- matrix(0:30,nrow(Q),ncol(Q),byrow=TRUE)

#DEFINE THE DOUBLY VARYING PENALTY MATRICES
# VARYING DIFFERENCE PENALTY APPLIED TO LAGS (EQ. 8b)
C <- do.call('onebasis',c(list(x=0:30,fun="ps",df=10,intercept=T)))
D <- diff(diag(30+1),diff=2)
P <- diag((seq(0,30-2))^2)
Slag1 <- t(C) %*% t(D) %*% P %*% D %*% C
#VARYING RIDGE PENALTY APPLIED TO COEFFs OF LAGS(Eq. 7a)
Slag2 <- diag(rep(0:1,c(6,4))) ### not sure why 10. needs to be same as df in C

# RUN THE GAM MODEL AND PREDICT (TAKES ~15sec IN A 2.4 GHz PC)
# PS: EXCLUDE THE DEFAULT PENALTY FROM THE LAG-RESPONSE FUNCTION WITH fx=T
# ADD ADDITIONAL PENALTY MATRICES FOR THE LAG SPACE IN addSlag OBJECT IN xt
xt <- list(addSlag=list(Slag1,Slag2))
#+s(R,L,bs="cb",k=10, fx=c(F,T),xt=xt)
#+ns(gw,knots=10 )
gam1 <-gam(VWC40p~s(Q,L,bs="cb",k=10, fx=c(F,T),xt=xt)+s(R,L,bs="cb",k=10, fx=c(F,T),xt=xt)
           ,family=quasipoisson(),data.stl,method='REML')
#+s(S,L,bs="cb",k=10, fx=c(F,T),xt=xt)
#+ns(gw,knots=10 )
pred3dgam1 <- crosspred("Q",gam1,at=0:50, cen=FALSE,cumul=F)## change range for SM20 - 0:30
predslgam1 <- crosspred("Q",gam1,at=0:50, cen=FALSE,cumul=F)

# CHECK CONVERGENCE, SMOOTHING PARAMETERS AND EDF
gam1$converged
gam1$sp
summary(gam1)$edf

#PLOTS
##Color ramp for 3D plots
z1 <- exp(pred3dgam1$matfit)
nrz <- nrow(z1)
ncz <- ncol(z1)
jet.colors <- colorRampPalette( c("lightpink1", "lightblue1") )
nbcol <- 100
color <- jet.colors(nbcol)
zfacet <- z1[-1, -1] + z1[-1, -ncz] + z1[-nrz, -1] + z1[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)


#pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/ITC13_dlnm_3D_col.pdf"),
 #   height=25/2.54,width=25/2.54)
par(mar=c(2,2,2,0.5),mgp = c(3,2,0), tck = -0.005)
plot(pred3dgam1,xlab="Surface Soil Moisture (in %)",zlab= "Relative Influence",ylab="Lag (Days)",main="SM13",phi=30, theta=35,cex.main = 1.5,r=sqrt(10),d=2,
     expand =1.5,col=color[facetcol], cex.lab = 1.25, cex.axis=1.25, zlim = c(0.4,1.2),## change c(0.4,1.2 for Sm20  to %50)
     ltheta=75,shade=0.25)## change to c(0.9,1.1)
#dev.off()

percentiles <- round(quantile(data.stl$VWC5p,c(0.001,0.05,0.95,0.99), na.rm = TRUE),0)
plot(predslgam1,var=percentiles,lag=c(2,5,7,10),ylab = "Relative Influence")
#pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/ITC13_dlnm_all1.pdf"),
#    height=10/2.54,width=16/2.54)
par(mar=c(4,5,3,2),mgp = c(2.2,0.75,0), tck = -0.01)
plot(predslgam1,"overall",ylab=expression(paste("Relative influence (", italic(hat(beta)),")")), cex.main=1.75,
     xlab="5cm Soil Moisture (in %)",lwd=1.5,main="SM13",
     ylim = c(0,3), cex.axis = 1.5,cex.lab = 1.5)
#dev.off()


###Get the RRfit all values for overall beta####
beta_overall<-predslgam1$allRRfit
beta_decoup<-round(beta_overall[beta_overall<1],2)
beta_decoup


####Select range of decoupled values from time series data#######
data09<-data.frame(data.ts)
data09$YYYYMMDD<-data$YYYYMMDD
data09$ratio5_40<-data09$VWC5/data09$VWC40
#data09<-data09[complete.cases(data09$VWC5),]
data09_decoup<-data09


data09_decoup$decoup_VWC5<-ifelse(data09$VWC5>=0.25,NA,data09$VWC5) ##for time series plotting
data09_decoup$group<-ifelse(data09$VWC5>=0.25 ,"coupled", "decoupled")
##SM13---- data13$VWC5>=0.07 & data13$VWC5<=0.34,NA,data09$VWC5)
##SM20---- data20$VWC5>=0.15 & data09$VWC5<=0.24,NA,data09$VWC5)
##SM05---- data05$VWC5<=0.24,NA,data05$VWC5)
##SM09----data09$VWC5>=0.28,NA,data09$VWC5)
#data09<-data09[complete.cases(data09$VWC5),]
decoup09<-data09[data09$VWC5<0.25  ,]#which(data05$VWC5<0.17|data05$VWC5>0.24),
couple09<-data09[data09$VWC5>=0.25,]
data09_decoup$group<-as.factor(data09_decoup$group)


#PLOTTING TIME SERIES - fix!
require(xts)
plot(data09_decoup$RDmm, type = "h")
decoup09.ts<-zoo(data09_decoup, data09_decoup$YYYYMMDD)

#pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/ITC20_TSdecoup.pdf"),
#    height=4,width =12)
par(mfrow = c(1,1), mar = c(4,4,2,4) + 0.01,mgp=c(1.5,0.55,0), tck=-0.01)
myColors <- c("grey38" ,"royalblue2","lightblue2","grey" , "darkviolet")

#Plot rainfall first so Soil moisture TS is more visible
plot(decoup09.ts$RDmm,type="h",lwd=1.5,  col="lightgrey",yaxt="n",xlab="",ylab="", ylim = rev(c(0,80)),xaxt = "n", 
     xlim= c(as.POSIXct('2014-01-01 23:45:00', format="%Y-%m-%d %H:%M:%S"),as.POSIXct('2016-12-31 23:45:00', format="%Y-%m-%d %H:%M:%S")))
axis(side=4,cex.axis = 1.4,mgp=c(0.5,0.35,0))
mtext(side = 4, line = 2, 'Rainfall (mm)',cex = 1.4)
legend(x="bottomleft",legend = c("5cm","40cm", "rainfall"), y.intersp = 0.75, x.intersp = 0.25,bty = "n",
       seg.len = 1.25, lwd = 2,lty = 1, col = c("lightblue2","royalblue2","grey"))
par(new = T)
plot(decoup09.ts[,c(13,3,1)],cex.lab = 1.3, xlab = "Time", ylab =expression("Soil Moisture " (cm^3/cm^3)),cex.axis=1.4,
     xlim= c(as.POSIXct('2014-01-01 23:45:00', format="%Y-%m-%d %H:%M:%S"),as.POSIXct('2016-12-31 23:45:00', format="%Y-%m-%d %H:%M:%S")),
     main = "SM20 ", col = myColors, lwd=c(7,1.5,1.5),screens=1, ylim = c(0,0.6))
#dev.off()


##BEAN PLOT
require(beanplot)
VWC5<-data09_decoup$VWC5p
VWC40<-data09_decoup$VWC40p
ratio5_40<-data09_decoup$ratio5_40
group<-data09_decoup$group


##BOXPLOT
#pdf(paste("/media/coleen/DDrive1/1_OWASIS/2_Statistical_analysis/ts_twente/dlnm/eda_dlnm_figs/figs_met_dec2016/ITC20_bean.pdf"),
 #   height=17/2.54,width=17/2.54)
par(mfrow = c(1,2),mar = c(1,3,2,1) + 0.01,mgp=c(1.8,0.5,0), tck=-0.02)
boxplot(ratio5_40, ylim = c(0, 2.0),boxwex = 0.5, ylab = "Ratio 5cm/40cm soil moisture",cex.lab =1.25,
        col ="grey50")
###Beanplot ratio
beanplot(ratio5_40~group, ll = 0.02,bw = 0.15, kernel = "gaussian", beanlines = "mean",  axes = F, log="",
         overallline = "median", method = "stack",
         main = "Beanplot SM20", side = "both", ylab="Ratio 5cm/40cm soil moisture", frame.plot = F, cex.lab =1.25, 
         col = list("grey25", c("grey90", "black")))
#axis(1)
axis(2 , seq(0,2, by = 0.5))
legend(x="topright", fill = c("grey25", "grey90"), bty = "n",x.intersp = 0.2,
       legend = c("Coupled", "Decoupled"), box.lty=0)
#dev.off()

#####T-test#########
### P-value less than α = 0.05, we reject the null hypothesis H0###
#test variance first
var.test(data09_decoup$ratio5_40~data09_decoup$group)
t.test(data09_decoup$ratio5_40~data09_decoup$group, paired = FALSE, var.equal = FALSE)

## normality assumption not followed fro SM20


###Wilcoxon-Mann-Whitney non parametric####
wilcox.test(data09_decoup$ratio5_40~data09_decoup$group, alternative = "two.sided" )


#######ZTEST########## 
## population variance
pop.var <- function(x) var(x) * (length(x)-1) / length(x)
z.test2sam = function(a, b, var.a, var.b){
    n.a = length(a)
    n.b = length(b)
    zeta = (mean(a) - mean(b)) / (sqrt(var.a/n.a + var.b/n.b))
    return(zeta)
  }

decoup<-decoup09$ratio5_40
var.de<-pop.var(decoup)
coupled<-couple09$ratio5_40
var.cou<-pop.var(coupled)


z.test2sam(decoup, coupled, var.de, var.cou)

############RAINFALL plots############
#1DAY rain
data09_decoup$D_ratio<-c(NA,diff(data09_decoup$ratio5_40,lag = 1))
data09_decoup<-data09_decoup[data09_decoup$rainy.days==1,]

2DAY rain
data09_decoup$D_ratio<-c(NA,NA,diff(data09_decoup$ratio5_40,lag = 2))
data09_decoup$RDmm2<-append(NA,rollapply(data09_decoup$RDmm,2,sum))
data09_decoup<-data09_decoup[data09_decoup$rainy.days==2,]

rain20<-data09_decoup[data09_decoup$VWC5<=0.2 ,]
rain30<-data09_decoup[data09_decoup$VWC5>0.2 & data09_decoup$VWC5<=0.3 ,]
rain40<-data09_decoup[data09_decoup$VWC5>0.3 & data09_decoup$VWC5<=0.4 ,]
rain50<-data09_decoup[data09_decoup$VWC5>0.4 & data09_decoup$VWC5<=0.5 ,]

#plot ratio 5/40 and rainfall day 1
plot(rain20$D_ratio,rain20$RDmm,col = "red", pch = 19, xlim = c(-0.1, 0.3),ylim = c(0,30))
points(rain30$D_ratio,rain30$RDmm,col = "orange", pch = 19)
points(rain40$D_ratio,rain40$RDmm, col = "green", pch = 19)
points(rain50$D_ratio,rain50$RDmm, col = "blue", pch = 19)


