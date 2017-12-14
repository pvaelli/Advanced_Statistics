#Pretend we take a Swiss survey of the Wallcreeper, a little cliff-inhabiting bird 
#that seems to have declined greatly in Switzerland in recent years. Assume that we 
#had data on the proportion of sample quadrats in which the species was observed 
#in Switzerland for the years 1990-2005 and that we were willing to assume that the 
#random deviations about a linear time trend were normally distributed. 

#This is for illustration only, usually, we would use logistic regression or a 
#site-occupancy model to make inference about such data that have to do with the 
#distribution of a species.

#Simulate the data
n <- 16  				# Number of years
a = 40					# Intercept
b = -1.5					# Slope
sigma2 = 25					# Residual variance

x <- 1:16 					# Values of covariate: here, year
eps <- rnorm(n, mean = 0, sd = sqrt(sigma2))	# Draw residuals from normal
y <- a + b*x + eps
plot((x+1989), y, xlab = "Year", las = 1, ylab = "Prop. occupied by Wallcreeper (%)", 
     cex = 1.2)


########################################################

#Analysis using JAGS (Bayesian)

#Load the correct library
library("R2jags")

# Write BUGS model
sink("linreg.txt")
cat("
    model {
    # Priors
    alpha ~ dnorm(0,0.001) # Note: small precision = huge variance (prec=1/variance)
    beta ~ dnorm(0,0.001)
    sigma ~ dunif(0, 100)   	# Note: Large variance = Small precision
    
    # Likelihood
    tau <- 1/ (sigma * sigma)
    for (i in 1:n) {
    y[i] ~ dnorm(mu[i], tau) 
    mu[i] <- alpha + beta*x[i]
    }
    
    # Derived quantities
    p.decline <- 1-step(beta)	# Probability of decline: step(beta) = if beta >=0; 0 otherwise
    
    }
    ",fill=TRUE)
sink()

# Package data
data <- list("x","y", "n")

# Inits function
inits <- function(){ list(alpha=rnorm(1), beta=rnorm(1), sigma = rlnorm(1))}

# Parameters to estimate
params <- c("alpha","beta", "p.decline", "sigma")

# MCMC settings
nc = 3  ;  ni=3500  ;  nb=500  ;  nt=3

# Start Gibbs sampler

out <- jags(data, inits, params, "linreg.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

#Look at the results
jagsout <- as.mcmc.list(out$BUGSoutput)
plot(jagsout)
summary(jagsout)


#Determining the trend
hist(out$BUGSoutput$sims.list$beta, 
     main = "Posterior distribution \nof Swiss Wallcreeper population trend", 
     col = "grey", xlab = "Occupancy-based population trend", xlim=c(-3.0,1.0))
abline(v = 0, col = "red")


#Plotting the results with credible intervals
plot((x+1989), y, xlab = "Year", las = 1, ylab = "Prop. occupied by Wallcreeper (%)", cex = 1.2)
abline(lm(y~ I(x+1989)), col = "blue", lwd = 2)
pred.y <- out$BUGSoutput$mean$alpha + out$BUGSoutput$mean$beta * x
points(1990:2005, pred.y, type = "l", col = "red", lwd = 2)
text(1994, 12, labels = "blue - ML; red - MCMC", cex = 1.2)

predictions <-array(dim=c(length(x),length(out$BUGSoutput$sims.list$alpha)))
for(i in 1:length(x)){
  predictions[i,]<-out$BUGSoutput$sims.list$alpha+out$BUGSoutput$sims.list$beta*i
}

LPB <-apply(predictions,1,quantile,probs=0.025)#Lowerbound
UPB <-apply(predictions,1,quantile,probs=0.975)#Upperbound

plot((x+1989),y,xlab="Year",las=1,ylab="Prop.occupied(%)",cex=1.2)
points(1990:2005,out$BUGSoutput$mean$alpha+out$BUGSoutput$mean$beta*x,type="l",col="black",
       lwd=2)
points(1990:2005,LPB,type="l",col="grey",lwd=2)
points(1990:2005,UPB,type="l",col="grey",lwd=2)









