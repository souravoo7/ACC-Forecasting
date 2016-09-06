#===================ADD LIBRARIES===================================
library('TSA');
library('forecast');
library('tseries');
library('moments');
library('lmtest');
library('sandwich');
library('qcc');
library('rugarch');
library('FinTS');
library('xts');
#=====================Data File & ADF test =========================
x=read.csv(file.choose(), sep=",", h=T)#read _ohlc files
y=read.csv(file.choose(), sep=",", h=T)# read the file containing 5 minute price data from 920 t0 1530 Hrs. for the month of June, 2014#
#z=read.csv(file.choose(), sep=",", h=T)# read the file containing close-to-close returns data for the month of June, 2014#
#calculate the returns 
returns<-diff(log(x$Adj.Close), lag=1);#price log returns
returns_hf<-diff(log(y$Close),lag=1);#high frequency price log returns
adf<-adf.test(returns);
par(mfrow=c(1,2));
acf(returns);
pacf(returns);
#================Determine the mean equation========================
m1=auto.arima(returns);
m2=arima(returns, order=c(1,0,0));
m3=arima(returns, order=c(0,0,1));
m4=arima(returns, order=c(1,0,1));
tsdiag(m1);
tsdiag(m2);
tsdiag(m3);
tsdiag(m4);
#futher analysis of means
m8=arima(returns, order=c(6,0,0), fixed=c(0,0,0,0,0,NA,NA));
m9=arima(returns, order=c(0,0,6), fixed=c(0,0,0,0,0,NA,NA));
m10=arima(returns, order=c(6,0,6), fixed=c(0,0,0,0,0,NA,0,0,0,0,0,NA,NA));
tsdiag(m8);
tsdiag(m9);
tsdiag(m10);
#============Returns Forecast======================================== 
results=matrix("na",nrow = 22, ncol = 1);
#forecasting ARMA(0,0)
for(i in 1:22){m=arima(returns[(i+1):(1216+i)], c(0,0,0)); results[i]=predict(m)$pred[1]};
write.csv(results, "m1.csv");
#forecasting ARMA(6,0)
for(i in 1:22){m=arima(returns[(i+1):(1216+i)], c(6,0,0)); results[i]=predict(m)$pred[1]};
write.csv(results, "m2.csv");
#forecasting ARMA(0,6)
for(i in 1:22){m=arima(returns[(i+1):(1216+i)], c(0,0,6)); results[i]=predict(m)$pred[1]};
write.csv(results, "m3.csv");
#forecasting ARMA(6,6)
for(i in 1:22){m=arima(returns[(i+1):(1216+i)], c(6,0,6)); results[i]=predict(m)$pred[1]};
write.csv(results, "m4.csv");

#=============Volatility Modelling===================================
#condunt the ARCH test
arch<-ArchTest(returns);
#======================================================================
  #sGARCH MODEL FITTING+normal distribution
s1<- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1),submodel="GARCH"), mean.model=list(armaOrder=c(1,0), include.mean=T), distribution="norm");
mv1<-ugarchfit(s1,returns);
#sGARCH MODEL FITTING+sGED distribution
s2<- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1),submodel="GARCH"), mean.model=list(armaOrder=c(1,0), include.mean=T), distribution="sged");
mv2<-ugarchfit(s2,returns);
#======================================================================
#eGARCH MODEL FITTING+norm distribution
s3<- ugarchspec(variance.model=list(model="eGARCH", garchOrder=c(1,1),submodel="GARCH"), mean.model=list(armaOrder=c(1,0), include.mean=T), distribution="norm");
mv3<-ugarchfit(s3,returns);
#eGARCH MODEL FITTING+sGED distribution
s4<- ugarchspec(variance.model=list(model="eGARCH", garchOrder=c(1,1),submodel="GARCH"), mean.model=list(armaOrder=c(1,0), include.mean=T), distribution="sged");
mv4<-ugarchfit(s4,returns);
#======================================================================
#gjrGARCH MODEL FITTING+norm distribution
s5<- ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1),submodel="GARCH"), mean.model=list(armaOrder=c(1,0), include.mean=T), distribution="norm");
mv5<-ugarchfit(s5,returns);
#gjrGARCH MODEL FITTING+sGED distribution
s6<- ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1),submodel="GARCH"), mean.model=list(armaOrder=c(1,0), include.mean=T), distribution="sged");
mv6<-ugarchfit(s6,returns);
#======================================================================
#EWMA estimation
s7<-ugarchspec(mean.model=list(armaOrder=c(1,0), include.mean=TRUE), variance.model=list(model="iGARCH",garchOrder=c(1,1),submodel="GARCH"), fixed.pars=list(alpha1=1-0.94, omega=0));
mv7<-ugarchfit(s7,returns);
#======================================================================
#SSRV Volatility Benchmark
rv=matrix("na", 21, 1);
for(i in 0:21){rv[i+1]=sum((returns_hf[((75*i)+1):((i*75)+75)])^2)};
#developing scalining factor.#
realised_variance=sum(as.numeric(rv[1:21]));
close_returns_variation=sum((returns_hf[1:22]-mean(returns_hf[1:22]))^2);
scale=close_returns_variation/realised_variance;
ssrv=as.numeric(rv)*scale;
annualised_ssrv=sqrt(ssrv*1);
write.csv(annualised_ssrv, "annualised ssrv.csv");
#======================================================================
#Forcasting volatility using models
for(i in 1:22 ){mvol=ugarchfit(s1, data=returns[(1+i):(1204+i)]); results[i]=sigma(ugarchforecast(mvol, n.ahead=1))[1]};
write.csv(results, 'mv1.csv');

for(i in 1:22 ){mvol=ugarchfit(s2, data=returns[(1+i):(1204+i)]); results[i]=sigma(ugarchforecast(mvol, n.ahead=1))[1]};
write.csv(results, 'mv2.csv');

for(i in 1:22 ){mvol=ugarchfit(s3, data=returns[(1+i):(1204+i)]); results[i]=sigma(ugarchforecast(mvol, n.ahead=1))[1]};
write.csv(results, 'mv3.csv');

for(i in 1:22 ){mvol=ugarchfit(s4, data=returns[(1+i):(1204+i)]); results[i]=sigma(ugarchforecast(mvol, n.ahead=1))[1]};
write.csv(results, 'mv4.csv');

for(i in 1:22 ){mvol=ugarchfit(s5, data=returns[(1+i):(1204+i)]); results[i]=sigma(ugarchforecast(mvol, n.ahead=1))[1]};
write.csv(results, 'mv5.csv');

for(i in 1:22 ){mvol=ugarchfit(s6, data=returns[(1+i):(1204+i)]); results[i]=sigma(ugarchforecast(mvol, n.ahead=1))[1]};
write.csv(results, 'mv6.csv');

for(i in 1:22 ){mvol=ugarchfit(s7, data=returns[(1+i):(1204+i)]); results[i]=sigma(ugarchforecast(mvol, n.ahead=1))[1]};
write.csv(results, 'mv7.csv');