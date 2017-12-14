#Assume five groups with 10 snakes measured in each with SVL averages of 
#50,40,45,55, and 60. This corresponds to a base line population mean of 
#50 and effects of populations 2-5 of -10, -5, 5, and 10. We choose a 
#residual standard deviation of SVL of 3 and assemble everything

#Simulate the data
ngroups <-5 # Number of populations
nsample <-10 # Number of snakes in each
pop.means <-c(50,40,45,55,60) # Population mean SVL
sigma <-3 # Residualsd
n <-ngroups*nsample # Total number of datapoints
eps <-rnorm(n,0,sigma) # Residuals
x <-rep(1:5,rep(nsample,ngroups)) # Indicator for population
means <-rep(pop.means,rep(nsample,ngroups))
X <-as.matrix(model.matrix(~as.factor(x)-1)) #Create design matrix

X # Inspect the data
y <- as.numeric(X%*%as.matrix(pop.means)+eps)  
                # assemble--NOTE:as.numericESSENTIALforWinBUGS
#Plot the data
boxplot(y~x,col="grey",xlab="Population",ylab="SVL",main="",las=1)

#######################################################################

#Fit fixed effect ANoVA with JAGS
#Load the correct library
library("R2jags")

#Here we will fit a means parameterization of the model

# Write model
sink("anova.txt")
cat("
  model {
  # Priors
   for (i in 1:5) {  #Implicitly define pop.mean as a vector
    pop.mean[i] ~ dnorm(0,0.001)
     }
  sigma ~ dunif(0,100)
  
  # Likelihood
   for (i in 1:50) {
    y[i] ~ dnorm(mean[i],tau.residual)
      mean[i] <- pop.mean[ x[i] ]
    }

  # Derived quantities
  # Done to determine whethere there are differences in SVL between pairwise comparisons of the populations  
    tau.residual <-1/(sigma*sigma)
    effe2 <-pop.mean[2]-pop.mean[1]
    effe3 <-pop.mean[3]-pop.mean[1]
    effe4 <-pop.mean[4]-pop.mean[1]
    effe5 <-pop.mean[5]-pop.mean[1]

    }

    ",fill=TRUE)
sink()

# Bundledata
win.data <-list("y","x")

# Initsfunction
inits <-function(){list(pop.mean=rnorm(5,mean=mean(y)),sigma=rlnorm(1))}

# Parameters to estimate
params<-c("pop.mean","sigma","effe2","effe3","effe4","effe5")

# MCMC settings
ni <-1200
nb <-200
nt <-2
nc <-3


out <- jags(win.data, inits, params, "anova.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb, 
             working.directory = getwd())

# Inspect estimates
jagsout <- as.mcmc.list(out$BUGSoutput)
plot(jagsout)
summary(jagsout)


#######################################################################
#######################################################################

#Modify the model assuming that the groups are random samples from a larger
#population

#Simulate the data
npop <-10 # Number of populations: now choose 10 rather than 5
nsample <-12 # Number of snakes in each
n <-npop*nsample # Total number of data points

pop.grand.mean <- 50 #Grand mean SVL
pop.sd <-5 # sd of population effects about mean
pop.means <- rnorm(n=npop,mean=pop.grand.mean,sd=pop.sd)
sigma <-3 # Residual sd
eps <-rnorm(n,0,sigma) # Draw residuals

x <- rep(1:npop,rep(nsample,npop))
X <- as.matrix(model.matrix(~as.factor(x)-1))
y <- as.numeric(X%*%as.matrix(pop.means)+eps) # as.numeric is ESSENTIAL
boxplot(y ~x,col="grey",xlab="Population",ylab="SVL",main="",las=1)
# Plot of generated data
abline(h =pop.grand.mean)


#######################################################################

#Fit random effect ANoVA with JAGS

#Write the model
sink("re.anova.txt")
cat("

model {
  # Priors and some derived things
    for (i in 1:npop){
      pop.mean[i] ~ dnorm(mu,tau.group) # Prior for population means
      effe[i] <- pop.mean[i] - mu  #Population effects as derived quant's
    }

    mu ~ dnorm(0,0.001) # Hyperprior for grand mean svl
    sigma.group ~ dunif(0,10) # Hyperprior for sd of population effects
    sigma.residual ~ dunif(0,10) # Prior for residual sd
    
  # Likelihood
    for (i in 1:n) {
      y[i] ~ dnorm(mean[i],tau.residual)
      mean[i] <- pop.mean[x[i]]
    }

  # Derived quantities
    tau.group <-1/(sigma.group*sigma.group)
    tau.residual <-1/(sigma.residual*sigma.residual)
    }
    ",fill=TRUE)
sink()

# Bundle data
win.data <-list(y=y, x=x, npop=npop, n=n)

# Inits function
inits<-function () {
  list(mu = runif(1,0,100), sigma.group = rlnorm(1), sigma.residual = rlnorm(1))  }

# Params to estimate
parameters <-c("mu","pop.mean","effe","sigma.group","sigma.residual")

# MCMCsettings
ni <-1200
nb <-200
nt <-2
nc <-3

# Start JAGS

out.RE <- jags(win.data, inits, parameters, "re.anova.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb, 
             working.directory = getwd())


#Look at the results
jagsout2 <- as.mcmc.list(out.RE$BUGSoutput)
plot(jagsout2)
summary(jagsout2)

save.image()









