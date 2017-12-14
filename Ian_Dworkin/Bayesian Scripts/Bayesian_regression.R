# A working example of a linear regression from a Bayesian perspective using several "canned" functions

####Chunk 1: Inputting and checking the data
setwd("/Users/ian/R/R scripts/Dll data/") 
dll.data = read.csv("dll.txt", header=TRUE)   #data frame input

# Quality Control of the data frame - making sure that it is inputted correctly
summary(dll.data)  # basic summary information for variables in the object
dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
dll.data$genotype <- relevel(dll.data$genotype, "wt")  # Setting the wild-type (wt) as reference
str(dll.data)  # check that this fixed the problem
dll.data <- na.omit(dll.data)  # I am only doing this to look at some diagnostic plots, which get all bothered by the missing data.

dll.data$tarsusCentered <- dll.data$tarsus - mean(dll.data$tarsus) # Center the tarsus length

require(MCMCpack)
require(arm)


lm.1 <- lm(SCT ~ 1+ tarsusCentered + genotype + temp + genotype:temp, data=dll.data)
summary(lm.1)
confint(lm.1)

# using MCMCregress in MCMCpack

lm.MCMC <- MCMCregress(SCT~ 1 + tarsusCentered + genotype + temp + genotype:temp, 
            data=dll.data, burnin=1000, mcmc=10000, thin=1, 
            #b0 are the hyper-parameters for the prior means 
            # (i.e. the location parameter for the priors for intercept and slopes)
            b0=c(0,0, 0, 0, 0 ), 
            #B0 are the hyper-paramereters for the prior precision (1/variance). 
            #i.e. The scale of the prior distributions for the intercept and slope 
            #(how flat or peaked are the priors)
            B0=c(0,0,0,0,0), c0 = 0.001, d0 = 0.001)


#c0 & d0 are the shape hyper-parameters for the prior distribution for the residual variation of the model (A gamma distribution is used I believe) 

par(mfrow=c(1,1))
plot(lm.MCMC)
HPDinterval(lm.MCMC)
summary(lm.MCMC)

# check some diagnostics
raftery.diag(lm.MCMC)
effectiveSize(lm.MCMC)
acf(lm.MCMC)




# When it comes to using pre-built functions, it is always good to double check your results against another function with comparable capabilities, or wrte out the code yourself.


# THe arm library (for the Gelman and Hill book) has a related function


lm.MCMC.arm <- bayesglm(SCT~ 1 + tarsusCentered + genotype + temp + genotype:temp, 
                 data=dll.data, family=gaussian, 
                 prior.mean=0, prior.mean.for.intercept=11.13, 
                 prior.scale.for.intercept=10^-5,
                 n.iter=50000, scaled=T)
                 
summary(lm.MCMC.arm)

