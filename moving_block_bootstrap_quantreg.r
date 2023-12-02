########################################################################
#moving block bootstrap for quantile regression
########################################################################
library(quantreg)

mlo=read.csv("mlo.csv")
head(mlo)

#Deseasonalization based on the 2000-2020 climatology
base=mlo[mlo$year>=2000 & mlo$year<=2020,]
seasonality=predict(lm(y~sin(2*pi*month/12)+cos(2*pi*month/12)+sin(2*pi*month/6)+cos(2*pi*month/6), data=base), newdata=data.frame(month=1:12))
mlo=merge(mlo, data.frame(month=1:12, yd=seasonality), by="month")
mlo$yd=mlo$y-mlo$yd
mlo=mlo[order(mlo$x),]

##moving block bootstrap function
mbfun=function(formula,data,tau){
 #data=rbind(data,data)
 n=nrow(data)
 b=ceiling(n^0.25) 
 nblocks=ceiling(n/b)
 blocks=lapply(seq_len(n-b+1), function(i) seq(i, i+b-1))
 bn=sample(1:length(blocks),nblocks,replace=T)
 samp_data=data[unlist(blocks[bn]), ]  
 mod=rq(formula, data=samp_data, tau=tau)
 coef(mod)
}

set.seed(2013)
formula=yd~x
##note: one can directly use formula=y~x+sin(2*pi*month/12)+cos(2*pi*month/12)+sin(2*pi*month/6)+cos(2*pi*month/6)
##if data are not deseasonalized in the previous step (but this means that data are deseasonalized using the entire period)
fit=coef(rq(formula=formula, data=mlo, tau=0.5)) #intercept and slope
op=t(replicate(1000, mbfun(formula=formula,data=mlo,tau=0.5)))
fit_se=apply(op, 2, sd, na.rm=TRUE) #MBB standard error for intercept and slope
fit_pv=2*pt(q=abs(fit/fit_se), df=nrow(mlo)-2, lower.tail=FALSE) #MBB p value for intercept and slop

#print summary
data.frame("fit_coef"=fit, "standard_error"=fit_se, "SNR"=fit/fit_se, "p_value"=fit_pv)
##The unit of trend value and SE is ppbv per month, 
##one can convert it to ppbv per year by multiplying a factor of 12.

#multiple quantiles 
#op=t(replicate(1000, mbfun(formula=yd~x,data=mlo,tau=seq(0.1,0.9,by=0.1))[2,]))
