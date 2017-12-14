##title: "Homework 5 - Class 14"

###Name: Patric Vaelli

###Due: October 27, 11:59pm

## 5 points total

#A population of ibex with initial size of 30 individuals was monitored during 25 years and grew 
#annually by about 2%. Assume that that population growth rate was normally distributed with a 
#mean of 2% and a variance (not standard deviation) of 2%.

#Shortly after the start of the annual snowmelt,ibexes are counted by the aid of a telescope from
#the slope opposite to the mountain ridge which constitutes the core area of the population. The 
#population survey is not perfect, and we assume that the variance of the observation error is 25. 

#Simulate the initial data parameters
n.years <- 25           # Number of years
N1 <- 30                # Initial population size
mean.lambda <- 1.02       # Mean annual population growth rate
sigma2.lambda <- 0.02     # Process (temporal) variation of the growth rate
sigma2.y <-  25         # Variance of the observation error

#Simulate the population sizes under this assumption of exponential growth
#Create a vector called N that will hold the true (unobserved) population sizes in years 1 through 25
N <- 1:25

#Define the population size in the first year of sampling
N[1] <- 30

#Define the year specific growth rates and put them in a vector called lambda
lambda <- rnorm(n.years, mean.lambda, sigma2.lambda)

#Write a loop that calculates the population sizes for N[2] through N[25] using the 
#population size from N[t-1] and the annual growth rates
for(i in N[2]:N[25]){
  pop.size <- N[i]*lambda[i]
}
pop.size  
  
  

#Simulate the observed data conditional on the true population sizes
#Call the vector of observed counts y
  
y <-   

#Create a loop to fill in the values for the y.    
#Note that the way we specified the basic state-space model, y can contain non integer values.  
#This is fine.  You can interpret these values as density rather than absolute counts.


  



######################################################################################

#Fit the model using JAGS.  All code is here for this but it's a good idea to review it so that
#you can start to understand the components that go into a Bayesian analysis in JAGS.

#Load the correct library
library("R2jags")

# Specify model in BUGS language
sink("ssm.txt")
cat("
    model {
    # Priors and constraints
    N.est[1] ~ dunif(0, 500)            #Prior for initial population size
    mean.lambda ~ dunif(0, 10)          #Prior for mean growth rate
    
    sigma.proc ~ dunif(0, 10)           #Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    
    sigma.obs ~ dunif(0, 100)           #Prior for sd of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # Stateprocess
    for (t in 1:(TT-1)) {
    lambda[t] ~ dnorm(mean.lambda, tau.proc)
    N.est[t+1]<-N.est[t]*lambda[t]
    }
    
    # Observationprocess
    for (t in 1:TT) {
    y[t] ~dnorm(N.est[t],tau.obs)
    }
    }
    ",fill =TRUE)
sink()

# Package the data
bugs.data <- list(y=y, TT=n.years)

# Generate Initial values
inits <- function() { list(sigma.proc=runif(1,0,5), mean.lambda=runif(1, 0.1,2), 
                           sigma.obs=runif(1,0,10), N.est=c(runif(1,20,40), 
                                                            rep(NA,(n.years-1))))  }

# Define the parameters to monitor
params <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc",
            "N.est")

# MCMC settings
nc <- 3    # Number of chains
ni <- 30000		# Number of draws from posterior
nb <- 20000		# Number of draws to discard as burn-in
nt <- 10		# Thinning rate

ssm <- jags(bugs.data, inits, params, "ssm.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, 
            working.directory = getwd())

#Look at the results
jagsout <- as.mcmc.list(ssm$BUGSoutput)
plot(jagsout)
summary(jagsout)

######################################################################################################

#Plot a figure of the data and the results.  Again, all code to do this is already included.

# Define function to draw a graph to summarize results
graph.ssm <- function(ssm, N, y){
  fitted <- lower <- upper <- numeric()
  n.years <- length(y)
  for (i in 1:n.years){
    fitted[i] <- mean(ssm$BUGSoutput$sims.list$N.est[,i])
    lower[i] <- quantile(ssm$BUGSoutput$sims.list$N.est[,i], 0.025)
    upper[i] <- quantile(ssm$BUGSoutput$sims.list$N.est[,i], 0.975)}
  m1 <- min(c(y, fitted, N, lower))
  m2 <- max(c(y, fitted, N, upper))
  par(mar = c(4.5, 4, 1, 1), cex = 1.2)
  plot(0, 0, ylim = c(m1, m2), xlim = c(0.5, n.years), ylab = "Population size", xlab = "Year", las = 1, 
       col = "black", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
  axis(2, las = 1)
  axis(1, at = seq(0, n.years, 5), labels = seq(0, n.years, 5))
  axis(1, at = 0:n.years, labels = rep("", n.years + 1), tcl = -0.25)
  polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
  points(N, type = "l", col = "red", lwd = 2)
  points(y, type = "l", col = "black", lwd = 2)
  points(fitted, type = "l", col = "blue", lwd = 2)
  legend(x = 1, y = m2, legend = c("True", "Observed", "Estimated"), lty = c(1, 1, 1), lwd = c(2, 2, 2), 
         col = c("red", "black", "blue"), bty = "n", cex = 1)
}

# Execute function: Produce figure 
graph.ssm(ssm, N, y)

