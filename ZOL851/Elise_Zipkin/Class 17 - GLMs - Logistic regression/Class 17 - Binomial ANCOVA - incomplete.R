#The adder has an all-black and a zigzag morph, where females are brown and males are gray. It has been 
#hypothesized that the black color confers a thermal advantage, and therefore, the proportion of black
#individuals should be greater in cooler or wetter habitats. We simulate data that bear on this question
#and "study," by simulation, 10 adder populations each in the Jura mountains, the BlackForest, and the Alps.
#We capture a number of snakes in these populations and record the proportion of black adders. Then, 
#we relate these proportions to the mountain range and to a combined index of low temperature, wetness, 
#and northerliness of the site. Our expectation will of course be that there are relatively more black 
#adders at cool and wet sites.

#1. Distribution: Ci ~ Binomial(pi, Ni)
#2. Link function: logit (pi) = linear predictor
#3. Linear predictor: alpha_Jura + Beta1 * xBlackF + Beta2 * xAlps + Beta3 * xwet + Beta4 *xwet * xBlackF + Beta5 * xwet * xAlps


#Simulate the data
n.groups = 3
n.sample = 10
n = n.groups *n.sample
x = rep(1:n.groups,rep(n.sample,n.groups))
pop = factor(x, labels=c("Jura","BlackForest","Alps"))

#We construct a continuous wetness index: 0 denotes wet sites lacking sun and 1 is the converse. For
#ease of presentation, we sort this covariate; this has no effect on the analysis.
wetness.Jura = sort(runif(n.sample,0,1))
wetness.BlackF <- sort(runif(n.sample,0,1))
wetness.Alps <- sort(runif(n.sample,0,1))
wetness <- c(wetness.Jura,wetness.BlackF,wetness.Alps)

#We also need the number of adders examined in each population (Ni), i.e., the binomial totals, also called
#sample or trial size of the binomial distribution. 

N = round(runif(n,10,50)) #Get discrete Uniform values

#Build the design matrix by combining population and wetness
Xmat = 
print(Xmat,dig=2)

#Select the parameter values, i.e., choose values for alpha_Jura, beta1 - beta5.
beta.vec = c(-4, 1,2,6,2, -5)

#We assemble the number of black adders captured in each population in three steps:
#  1. We add up all components of the linear model to get the value of the linear predictor,
#  2. We apply the inverse logit transformation to get the expected proportion (pi) of black 
#adders in each population
#  3. We add binomial noise, i.e., use pi and Ni to draw binomial random numbers representing the 
#count of black adders in each sample of Ni

lin.pred <-                    #Value of lin. predictor
exp.p <-                       #Expected proportion
C <-                           # Add binomial noise

hist(C)                   # Inspect what we've created

par(mfrow =c(2,1))
matplot(cbind(wetness[1:10],wetness[11:20],wetness[21:30]), cbind(exp.p[1:10],
          exp.p[11:20], exp.p[21:30]), ylab="Expected proportion of black morphs", xlab="", col=c("red","green","blue"),
        pch=c("J","B","A"), lty="solid", type="b", las=1, cex =1.2,main="", lwd=1)

matplot(cbind(wetness[1:10],wetness[11:20],wetness[21:30]),
        cbind(C[1:10]/N[1:10],C[11:20]/N[11:20],C[21:30]/N[21:30]), ylab="Observed proportion of black morphs", 
        xlab="Wetness index",col=c("red","green","blue"),pch= c("J","B","A"), las=1,cex=1.2,main="")
par(mfrow =c(1,1))

####################################################################################

#Analyze the data in R with MLE

#Note the glm function require the number of failures, not the total number of observed N.
#Thus we create a vector indicating the number of zigzag morphs that we counted at each site.

Z=N-C

bin.ancova = 

#Examine the results and compare the effect sizes with the data generating effect sizes


#Pull out the parameter estimates and their confidence intervals.

  
  
#################################################################################

# Define model
sink("bin_glm.txt")
cat("
    model {
    
    # Priors - use a means parameterization
    for (i in 1:n.groups){
    alpha[i] ~ 
    beta[i] ~ 
    }
    
    # Likelihood - use a means parameterization
    #Is this a means or an effect parameterization?
    for (i in 1:n){
    C[i] ~ 
    logit(p[i]) <- 
    }
    
    # Derived quantities
    # Recover the effects relative to baseline level (no.1)
    a.effe2 <- alpha[2] - alpha[1] #Intercept Black Forest vs. Jura
    a.effe3 <- alpha[3] - alpha[1] #Intercept Alps vs. Jura
    b.effe2 <- beta[2] - beta[1] # Slope Black Forest vs. Jura
    b.effe3 <- beta[3] - beta[1] # Slope Alps vs. Jura
    
    }
    ",fill=TRUE)
sink()

# Bundle data
data = list(    )

# Inits function
inits = function(){ list(alpha=rnorm(n.groups), beta=rnorm(n.groups))}

# Parameters to estimate
params = c(   )

# MCMC settings
ni = 7000
nb = 2000
nt = 5
nc = 3

#Load the r2jags library
library("R2jags")

# Start JAGS
out = jags(data =data, inits=inits, parameters.to.save=params, model.file = "bin_glm.txt", n.thin=nt, 
           n.chains=nc, n.burnin=nb,n.iter=ni)

plot(out)
jagsout <- as.mcmc.list(out$BUGSoutput)
plot(jagsout)

# Examine the results from the Bayesian analysis
print(out, dig = 3) 

#Which parameters correspond to those that we used to generate the data?
beta.vec        # Truth in data generation

#Can you create a figure plotting the expected proportion of black morphs relative to the wetness indicator in 
#each of the three regions including the 95% credible intervals?


