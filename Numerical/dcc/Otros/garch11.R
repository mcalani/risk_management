library(zoo)
library(xts)
library(rugarch)
library(tseries)
library(rmgarch)

# ---------------------------------------------------------------------------- #
#                                   LOAD DATA                                  #
# ---------------------------------------------------------------------------- #

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data_csv <- read.csv("tc24febrero.csv")
data <- na.omit(data_csv)
day <- as.Date(data[,1], format = "%d/%m/%y")
TC  <- as.xts(zoo(data$Dolar.observado,order.by=day))
ret_tc = diff(log(as.numeric(TC)))
#plot.ts(ret_tc)

# ---------------------------------------------------------------------------- #
#                               RUGARCH ESTIMATE                               #
# ---------------------------------------------------------------------------- #

# Specify a standard GARCH model with constant mean
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                   mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), 
                   distribution.model = "norm")

# Estimate the model

fit <- ugarchfit(data = ret_tc, spec = spec, method="BFGS")

fit
coef(fit)
# Use the method sigma to retrieve the estimated volatilities 
garchvol <- c(mean(sigma(fit)),sigma(fit))

# Plot the volatility for 2017
vol <- zoo(garchvol, order.by = day)

# plot(vol, xaxt = "n")
# axis(1, day, format(day, "%m/%Y"), cex.axis = .7)

# ---------------------------------------------------------------------------- #
#                               Manual estimation                              #
# ---------------------------------------------------------------------------- #

garch_loglik<-function(para,x,mu){
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
  -loglik
}

# Optim Routine
para=c(50*var(ret_tc),0.2,0.9)
mlef<-optim(para, garch_loglik, gr = NULL,method = c("BFGS"),ret_tc,0)
mlef$par

# NLMINB Routine
garch_loglik2<-function(para,x){
  # Parameters
  #mu = mean(x)
  mu    = para[1]
  omega = para[2]
  alpha = para[3]
  beta  = para[4]
  # Volatility and loglik initialisation
  loglik=0
  h=var(x)
  # Start of the loop
  #vol=c()
  for (i in 2:length(x)){
    h      = omega + alpha*(x[i-1]-mu)^2 + beta*h
    loglik = loglik + dnorm(x[i],mu,sqrt(h),log=TRUE)
  }
  print(para)
  -loglik
}

LLike <- function(para) { garch_loglik2(para,ret_tc) }

S=1e-06
params = c(mu= mean(ret_tc), omega0 = 0.1 * var(ret_tc), alpha = 0.1, beta = 0.8)
lowerBounds = c(mu= S**2 * abs(mean(ret_tc)),omega0 = S^2, alpha = S, beta = S)
upperBounds = c(mu=10 * abs(mean(ret_tc)),omega0 = 100 * var(ret_tc), alpha = 1 - S, beta = 1 - S) 
mlef2 <-nlminb(start = params, objective = LLike,lower = lowerBounds, upper = upperBounds,
                control = list(abs.tol = 1e-30))

mlef2$par
coef(fit)
mean(ret_tc)

LLike(mlef2$par)
LLike(coef(fit))



class(ret_tc)

## WITH EXTERNAL REGRESSORS
#fb4 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1),
#                                      submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
#                mean.model = list(armaOrder = c(0, 0), external.regressors = vix),
#                distribution.model = "std", start.pars = list(), fixed.pars = list())



