library(quantmod)
library(rugarch)


# DATA
amazon<-getSymbols("^DJI",src="yahoo", from="2012-01-01",to="2019-01-01", auto.assign = F)
amazonlog<-as.data.frame(log10(amazon[,6]))
N=length(amazonlog[,1])
amazon_ret<-100*(amazonlog[2:N,]-amazonlog[1:N-1,])

# RUGARCH

modelx<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                   mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), 
                   distribution.model = "norm")
setstart(modelx) <- list(omega=0.2, alpha1 = 0.8, beta1 = 0.2)
fitx = ugarchfit(data =amazon_ret , spec = modelx, method="BFGS")
coef(fitx)

#MANUAL
garch_loglik<-function(para,x){
  mu = mean(x)
  # Parameters
  omega0=para[1]
  alpha=para[2]
  beta=para[3]
  # Volatility and loglik initialisation
  loglik=0
  h=var(x)
  # Start of the loop
  vol=c()
  for (i in 2:length(x)){
    h=omega0+alpha*(x[i-1]-mu)^2+beta*h
    loglik=loglik+dnorm(x[i],mu,sqrt(h),log=TRUE)
  }
  print(para)
  return(-loglik)
}

para=c(0.2,0.8,0.2)

mlef<-optim(para, garch_loglik, gr = NULL,method = c("BFGS"),amazon_ret)

mlef$par

typeof(amazon_ret)
class(amazon_ret)
