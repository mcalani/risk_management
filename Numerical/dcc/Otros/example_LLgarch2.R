
# ---------------------------------------------------------------------------- #
#                              ESTIMATE GARCH(1,1)                             #
# ---------------------------------------------------------------------------- #

### FUNCTION

garch11Fit = function(x) {
  require(numDeriv)
  # Step 1: Initialize Model Parameters and Bounds:
  Mean = mean(x)
  Var = var(x)
  S = 1e-06
  params = c(mu = Mean, omega = 0.1 * Var, alpha = 0.1, beta = 0.8)
  lowerBounds = c(mu = -10 * abs(Mean), omega = S^2, alpha = S, beta = S)
  upperBounds = c(mu = 10 * abs(Mean), omega = 100 * Var, alpha = 1 - S, beta = 1 - S)
  # Step 2: Set Conditional Distribution Function:
  garchDist = function(z, hh) {
    dnorm(x = z/hh)/hh
  }
  # Step 3: Compose log-Likelihood Function:
  garchLLH = function(parm) {
    mu = parm[1]
    omega = parm[2]
    alpha = parm[3]
    beta = parm[4]
    z = (x - mu)
    Mean = mean(z^2)
    # Use Filter Representation:
    e = omega + alpha * c(Mean, z[-length(x)]^2)
    h = filter(e, beta, "r", init = Mean)
    hh = sqrt(abs(h)) -sum(log(garchDist(z, hh)))  #llh
    
  }
  # print(garchLLH(params)) Step 4: Estimate Parameters and Compute
  # Numerically Hessian:
  fit = nlminb(start = params, objective = garchLLH, lower = lowerBounds, 
               upper = upperBounds)
  Hessian <- numDeriv::hessian(func = garchLLH, x = fit$par)
  # Step 5: Create and Print Summary Report:
  se.coef = sqrt(diag(solve(Hessian)))
  tval = fit$par/se.coef
  matcoef = cbind(fit$par, se.coef, tval, 2 * (1 - pnorm(abs(tval))))
  dimnames(matcoef) = list(names(tval), c(" Estimate", " Std. Error", " t value", 
                                          "Pr(>|t|)"))
  cat("\nCoefficient(s):\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
}


## EXAMPLE: AMAZON

library(tseries) #to extract quotes from internet
#Obtain Microsoft series from Internet
msft.prices = get.hist.quote(
  instrument = "MSFT", 
  quote = "Close",  #adjusted does not exist anymore
  provider = c("yahoo"), origin = "1999-12-30", start="2000-01-01", end="2010-01-01",
  retclass = c("zoo"), quiet = FALSE, drop = FALSE)

#Transform to log-returns
msft <- as.data.frame(msft.prices)
N <- length(msft[, 1])
msft.returns <- 100*(log(msft[2:N, ])-log(msft[1:(N-1), ]))

# Fit with function akin to that found in library(fGarch)
garch11_model <- garch11Fit(msft.returns)

## EXAMPLE: Exchange rate Chile
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data_csv <- read.csv("tc24febrero.csv")
data <- na.omit(data_csv)
day <- as.Date(data[,1], format = "%d/%m/%y")
TC  <- as.xts(zoo(data$Dolar.observado,order.by=day))
ret_tc = diff(log(as.numeric(TC)))

garch11_TC <- garch11Fit(ret_tc)

# ---------------------------------------------------------------------------- #
#                                     OTROS                                    #
# ---------------------------------------------------------------------------- #

kurtosis_garch11 <- function(alpha1, beta1) {
    3 * (1 + alpha1 + beta1) * (1 - alpha1 - beta1)/(1 - beta1^2 - 2 * alpha1 * 
        beta1 - 3 * alpha1^2)
}
kurtosis_garch11(garch11_model[3, 1], garch11_model[4, 1])

library(rugarch)
# Specification of GARCH(1, 1) model using rugarch workhorse with ARMA(1,0)
# + mean, normal errors
model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 
    1)), mean.model = list(armaOrder = c(1, 0), include.mean = TRUE), distribution.model = "norm")
# Model fitting
model_fit <- ugarchfit(spec = model, data = msft.returns, solver = "nloptr", 
    solver.control = list(solver = 9))
# Did the optimization routine converge?
convergence(model_fit)  #0 == TRUE, indicating convergence

plot(model_fit, which = 1)  # series with 95% conf. int (+/- 2 conditional std. dev.)
plot(model_fit, which = 2)  #VaR 99%