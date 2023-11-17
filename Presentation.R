library(reshape2)
library(ggplot2)
library(quantmod)
library(timeSeries)
library(car)
library(rugarch)
library(ggExtra)

rm(list=ls())

data <- get(getSymbols("BAC", from = "2003-11-03", to = "2023-11-01", frequency = "daily", src = "yahoo", auto.assign=TRUE))
#Check if there is NA value
sum(is.na(data))

#calculate returns and delete N/A value
returns <- diff(log(data$BAC.Adjusted))
names(returns) <- c("returns")
returns <- returns[complete.cases(returns)]

#plot returns
barChart(returns,theme = 'white')

r_BAC <- as.vector(returns)
plot(r_BAC, main = "BCA", type = "l")

acf(returns, main = "Autocorrelation of of returns")
acf(returns^2, main = "Autocorrelation of of returns squared")


# Create the SPECIFICATIONS. ----------------------------------------------------------

#GARCH in Mean
spec_garch_m <- ugarchspec(
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE, archm = TRUE, archpow = 1),
  variance.model = list(garchOrder = c(1, 1)),
  distribution.model = "norm"
)

# Fit the model
GARCH_Mean <- ugarchfit(spec = spec_garch_m, data = returns)

#GARCH in Mean with ARMA
spec_garch_ARMA <- ugarchspec(
  variance.model = list(garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(1,0), include.mean = TRUE,archm = TRUE, archpow = 1),
  distribution.model = "norm"
)

#Fit the model
GARCH_ARMA <- ugarchfit(spec = spec_garch_ARMA, data = returns)


# Plot --------------------------------------------------------------------
plot(GARCH_Mean, which = "all")

plot(GARCH_ARMA, which = "all")

#Plotting the estimated conditional volatility
plot(GARCH_Mean@fit$sigma, type = "l", main = "Esimated conditional volatility for BAC(Mean)", ylab="Conditional volatility")

plot(GARCH_ARMA@fit$sigma, type = "l", main = "Esimated conditional volatility for BAC(ARMA)", ylab="Conditional volatility")

#-----------------------------------------------------------------

#Log-likelihood ratio test-------------------

# Number of restrictions GARCH_Mean vs GARCH_ARMA
coef(GARCH_Mean)
coef(GARCH_ARMA)

# Creating the LR statistic - using cat() for output display
cat(" Likelihood of GARCH_Mean: ", round(likelihood(GARCH_Mean),2), "\n", 
    "Likelihood of GARCH_ARMA: ", round(likelihood(GARCH_ARMA),2), "\n",
    "2 * (log Lu - log Lr): ", round(2*(likelihood(GARCH_Mean)-likelihood(GARCH_ARMA)),2))

# 95% quantile of chi-square
qchisq(p = 0.95, df = 1)
#3.841459>3.8, hence the unrestricted model adds value!


# Adding the p-value
p.value  <- 1 - pchisq(2*(likelihood(GARCH_Mean)-likelihood(GARCH_ARMA)), df = 1)
cat(" Likelihood of GARCH in Mean: ", round(likelihood(GARCH_Mean),2), "\n", 
    "Likelihood of GARCH in Mean with ARMA: ", round(likelihood(GARCH_ARMA),2), "\n",
    "2 * (Lu - Lr): ", round(2*(likelihood(GARCH_Mean)-likelihood(GARCH_ARMA)),2), "\n",
    "p-value:", p.value)

#Residuals of Garch in mean--------------------------------------
# Performing residual analysis
residuals_mean <- GARCH_Mean@fit$residual

# Autocorrelation of residuals
acf(residuals_mean, main = "Autocorrelation of residuals of Garch in mean")

# Autocorrelation of residuals squared
acf(residuals_mean^2, main = "Autocorrelation of residuals squared of Garch in mean")

# Jarque-Bera test for Normality
jarque.bera.test(residuals_mean) #pvalue small, hence reject normality

# Ljung-Box test for residuals
Box.test(residuals_mean, type = "Ljung-Box")
# Ljung-Box test for residuals squared
Box.test(residuals_mean^2, type = "Ljung-Box")

# qqPlot for residuals under a normal distribution
qqPlot(residuals_mean)

# qqPlot for residuals under a Student-t with 3 degrees of freedom
qqPlot(residuals_mean, distribution = "t", df = 3, envelope = FALSE)

#Residuals of Garch in mean with ARMA--------------------------------------
# Performing residual analysis
residuals_ARMA <- GARCH_ARMA@fit$residual

# Autocorrelation of residuals
acf(residuals_ARMA, main = "Autocorrelation of of residuals of Garch in mean with ARMA")

# Autocorrelation of residuals squared
acf(residuals_ARMA^2, main = "Autocorrelation of residuals squared of of Garch in mean with ARMA")

# Jarque-Bera test for Normality
jarque.bera.test(residuals_ARMA) #pvalue small, hence reject normality

# Ljung-Box test for residuals
Box.test(residuals_ARMA, type = "Ljung-Box")
# Ljung-Box test for residuals squared
Box.test(residuals_ARMA^2, type = "Ljung-Box")

# qqPlot for residuals under a normal distribution
qqPlot(residuals_ARMA)

# qqPlot for residuals under a Student-t with 3 degrees of freedom
qqPlot(residuals_ARMA, distribution = "t", df = 3, envelope = FALSE)
