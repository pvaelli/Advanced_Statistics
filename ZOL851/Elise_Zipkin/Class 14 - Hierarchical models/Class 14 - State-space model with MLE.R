#Structural time series models (Commandeur et al. 2011, Section 2.3) can be easily implemented
#in R through the function StructTS by B.D. Ripley, included in the base package stats.
#StructTS has the advantage of being of simple usage and quite reliable. It gives the main tools 
#for fitting a structural model for a time series by maximum likelihood; the options "level", 
#"trend", "BSM" are used to t a local level model, a local linear trend, or a local trend with 
#an additional seasonal component. The Nile river data are included in any standard distribution 
#of R as a time series object (i.e., a vector containing the data together with information about 
#start/end time and sampling frequency); a detailed description of the data is given in the help file ?Nile.

################################################################################

#Load the packages that we will need
library("dlm")
library("KFAS")
library("forecast")

fitNile <- StructTS(Nile, "level")

#The maximum likelihood estimates (MLEs) of the level and observation error variances, 1469
#and 15099, respectively, are included in the output, fitNile as fitNile$coef. Asymptotic standard errors are not provided.

plot(Nile, type = "o")
lines(fitted(fitNile), lty = "dashed", lwd = 2)
lines(tsSmooth(fitNile), lty = "dotted", lwd = 2)

plot(forecast(fitNile, level = c(50, 90), h = 10), xlim = c(1950, 1980))

#A polynomial DLM (a local level model is a polynomial DLM of order 1, a local linear trend is a polynomial 
#DLM of order 2), is easily defined in dlm through the function dlmModPoly.

#MLEs are computed in package dlm by the function dlmMLE, which is a wrapper to optim,
#the general-purpose optimizer included in package stats. The main arguments are the data,
#a starting value for the unknown (vector) parameter, and a function that sets up a dlm
#object using the unknown parameter. This is internally combined with a call to dlmLL,
#which evaluates the negative log likelihood, and passed to optim, which performs the actual
#optimization.

buildNile <- function(theta) {
   dlmModPoly(order = 1, dV = theta[1], dW = theta[2])
}

#The MLEs are then obtained by calling dlmMLE as follows:

fit <- dlmMLE(Nile, parm = c(100, 2), buildNile, lower = rep(1e-4, 2))

modNile <- buildNile(fit$par)
smoothNile <- dlmSmooth(Nile, modNile)

filterNile <- dlmFilter(Nile, modNile)

#Plot the residuals
plot(residuals(filterNile, sd = FALSE), type = "o",
         ylab = "Standardized prediction error")
abline(h = 0)

#Forcast the results
foreNile <- dlmForecast(filterNile, nAhead = 10)
 attach(foreNile)
hwidth <- qnorm(0.25, lower = FALSE) * sqrt(unlist(Q))
fore <- cbind(f, as.vector(f) + hwidth %o% c(-1, 1))
rg <- range(c(fore, window(Nile, start = c(1951, 1))))
plot(fore, type = "o", pch = 16, plot.type = "s", lty = c(1, 3, 3),
         ylab = "Nile level", xlab = "", xlim = c(1951, 1980), ylim = rg)
lines(window(Nile, start = c(1951, 1)), type = 'o')
lines(window(smoothNile$s, start = c(1951,1)), lty = 5)
abline(v = mean(c(time(f)[1], tail(time(Nile), 1))),
           lty = "dashed", col = "darkgrey")
legend("topleft", lty = c(1, 5, 1, 3), pch = c(1, NA, 16, 16), bty = "n",
           legend = c("observed level", "smoothed level", "forecasted level",
                       "50% probability limits"))
detach(foreNile)

####################################################################
#Using the KFAS package
#The code below only works with  Version 0.6.1 of the KFAS function

#Write out the log likelihood
logLik <- function(theta) {
  lik <- kf(yt = Nile, Zt = 1, Tt = 1, Rt = 1, Ht = theta[1],
        Qt = theta[2], a1 = 0, P1 = 1e7)
   return(-lik$lik)
  }
 
fit <- optim(par = c(100, 2), fn = logLik, lower = rep(1e-4, 2))

filterNile <- kf(yt = Nile, Zt = 1, Tt = 1, Rt = 1, Ht = fit$par[1],
                  Qt = fit$par[2], a1 = 0, P1 = 1e7)
smoothNile <- ks(filterNile)
attach(smoothNile)
hwidth <- qnorm(0.05, lower = FALSE) * sqrt(drop(Vt))
sm <- cbind(drop(ahat), as.vector(ahat) + hwidth %o% c(-1, 1))
sm <- ts(sm, start = start(Nile))
plot(sm, plot.type = "s", type = "l", lty = c(1, 5, 5),
        + ylab = "Level", xlab = "", ylim = range(Nile))
lines(Nile, type = "o", col = "darkgrey")
legend("bottomleft", col = c("darkgrey", rep("black", 2)),
          + lty = c(1, 1, 5), pch = c(1, NA, NA), bty = "n", legend =
            + c("data", "smoothed level", "90% probability limits"))
detach(smoothNile)

residNile <- drop(filterNile$vtuni / sqrt(filterNile$Ftuni))

foreNile <- forecast(filterNile, fc = 9)
attach(foreNile)
hwidth <- qnorm(0.25, lower = FALSE) * sqrt(drop(Pt.fc))
fore <- ts(cbind(drop(at.fc), drop(at.fc) + hwidth %o% c(-1, 1)),
               start = 1 + end(Nile)[1])
rg <- range(c(fore, window(Nile, start = c(1951, 1))))
detach(foreNile)




