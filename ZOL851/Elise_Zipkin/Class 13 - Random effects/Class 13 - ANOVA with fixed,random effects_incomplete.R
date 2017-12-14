##Recall: One-way ANOVA with fixed effects

#Assume 8 groups with 15 snakes measured in each with SVL averages of 
#50,40,45,55,60,35,75,80. 
#We choose a residual standard deviation of SVL of 3 and assemble everything

#Simulate the data
ngroups <-8 # Number of populations
nsample <-15 # Number of snakes in each
pop.means <-c(50,40,45,55,60,35,65,50) # Population mean SVL
sigma <-3 # Residualsd
n <-ngroups*nsample # Total number of datapoints
eps <-rnorm(n,0,sigma) # Residuals
x <-rep(1:ngroups,rep(nsample,ngroups)) # Indicator for population
means <-rep(pop.means,rep(nsample,ngroups))
X <-as.matrix(model.matrix(~as.factor(x)-1)) #Create design matrix

X # Inspect the data
y <- as.numeric(X%*%as.matrix(pop.means)+eps)  
                # assemble--NOTE:as.numericESSENTIALforWinBUGS
#Plot the data
boxplot(y~x,col="grey",xlab="Population",ylab="SVL",main="",las=1)

#What is the mean and SD of the data?
mean(y); sd(y)

#######################################################################

# Likelihood model for an AVOVA
#
#  y[i] ~ dnorm(pop.means[x[i]],sigma)
#
# i is each sample 1:n
# x is groups 1:ngroups

#######################################################################

#Fit a fixed effect ANOVA in R
#Remember R fits an effects parameterization

fit=lm(y~as.factor(x))

#Look at the results
anova(fit)
summary(fit)

#Can also use the aov function, which is not as useful because it does not give the population means,
#just whether or not at least one mean is different from the others
fit2=aov(y~as.factor(x))
summary(fit2)

#######################################################################
#######################################################################

#Modify the model assuming that the groups are random samples from a larger
#population

#Simulate the data
npop <-8 # Number of populations: 8
nsample2 <-12 # Number of snakes in each
n2 <-npop*nsample2 # Total number of data points

pop.grand.mean <- 50 #Grand mean SVL
pop.sd <-5 # sd of population effects about mean

pop.means <- rnorm(npop, pop.grand.mean, pop.sd)

sigma2 <-3 # Residual sd
eps2 <-   rnorm(n2, 0, sigma2) # Draw residuals

x2 <- rep(1:npop, rep(nsample2,npop))
X2 <- as.matrix(model.matrix(~as.factor(x2)-1))
y2 <- as.numeric(X2 %*% as.matrix(pop.means) + eps2) # as.numeric is ESSENTIAL
boxplot(y2 ~x2,col="grey",xlab="Population",ylab="SVL",main="",las=1)
# Plot of generated data
abline(h =pop.grand.mean)

#What is the overall mean and sd of the data?

mean(y2); sd(y2)

#######################################################################

# Likelihood model for an random effects AVOVA 
#
#  y2[i] ~ dnorm(pop.means[x2[i]],sigma)
#  pop.means[x] ~ dnorm(pop.grand.mean,pop.sd) 
#
# i is each sample 1:n2
# x2 is groups 1:ngroups

#######################################################################

#Fit a random effect ANOVA in R
pop <- as.factor(x2) #Define x as a factor and call it pop

#Need the lme4 library for this analysis
library('lme4')
lme.fit <- lmer(y2 ~ 1 + (1 | pop))

# Inspect results

lme.fit
summary(lme.fit)
#Print random effects and remember they are relative to the grand mean

  
#Compare this to population level means


#Can run this model with the aov function if the data are balanced (which they are) 
#Again results may not be as useful
fit3=aov(y2 ~ Error(pop) )
fit3

#What are we interested in when using a random-effects ANOVA compared to using a fixed-effects ANOVA












