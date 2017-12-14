# Last modified October 10th 2011

# Chunk 1: The basics of a power analysis using R pre-built functions.
# Chunk 2: Power analysis across multiple levels of "power" using sapply.
# Chunk 3: Using Monte Carlo methods for power analysis: basics
# Chunk 4: The basic function for a linear model to generate monte carlo samples under a model, and extract a p-value.
# Chunk 5: Some additional thoughts on Monte Carlo methods for power analysis
# Chunk 6: Extracting p-values for the whole model.


###Chunk 1: The basics of a power analysis using R pre-built functions.

# Let's say we wanted to know what value of "t" from a two-tailed t-test we would need to exceed to "reject" the null hypothesis at alpha=0.05

# Let's start assuming that we have sample size n=25 (df=n-1)
qt(p=0.975,df=24)

# In general this single point will not be particularly useful to us. 
# We can ask how this value of t changes as we increase sample sizes.
curve(qt(p=0.975,df=x),2,25, xlab = "Sample Size", ylab= "t value for two sided test, alpha=0.05")


# Of course this is also of pretty limited utility since we usually want to ask these questions with respect to parameters we have estimated (like  difference between group1 and group2, and their standard deviations), and start with (like sample size and alpha).


# Remember all that a t-test is essentially t = (mean(A) - mean(B))/(SD/sqrt(n))
# This t is also called the non-centrality parameter (ncp)
# This test statistic is the non-central t distribution

# We can now ask for a given set of parameters what the power will be for the t-test.
# Usually we will be interested in something like what sample size we will need to achieve a certain power level

pow1 <- power.t.test(delta=0.5, sd=2, sig.level=0.05, power=0.8)
pow1$n # So we need ~ 252 observations to achieve a power of 0.8



### Chunk 2: Power analysis across multiple levels of "power" using sapply.
# We can also do this to see how sample size requirements change with changing power.
power = seq(from=0.1, to=0.9, by=0.01) # creates a vector of values from 0.1 to 0.9

# This function is just so that we can compute the sample size we need for varying levels of power
pow.test <- function(x){
	pow2 <- power.t.test(delta=0.5, sd=2, sig.level=0.05, power=x) # We only allow power to vary.
	return(pow2$n) # This pulls out the sample size we need
	}

power.n <- sapply(power, pow.test) # This just uses one of the apply functions to repeat the function (pow.test) for each element of the vector "power". Thus for each value in the vector "power" (from 0.1 to 0.9), it inputs this value into pow.test() and then returns the estimated n (# of observations needed to achieve this power).


plot(power.n ~ power, ylab="# of Observations", main=" sample size vs power for a t-test with delta=0.5, sd=2 and alpha=0.95")




##### Chunk 3: Using Monte Carlo methods for power analysis: basics
# As Described in Ch.8 of Dalgaard there are many functions that can figure out the power of the test, or sample size etc...   However these are generally only useful for very simple designs (t-tests, 1 or 2 way ANOVA's simple linear regression). Our problems do not fall into these categories. 

 # What can we do? WE CAN SIMULATE OUR EXPERIMENTS!!!! (See chapters 7 & 8 of Gelman and Hill, and Chapter 5 of Bolker)
 # This allows us to not only ask how often we can reject the null hypothesis at a given alpha for #that model (and set of parameters), but we can easily change it to test a variety of parameters, #sample sizes and experimental designs. 
 
# See the Monte Carlo methods for inference script and lecture for a basic reminder, before proceeding through the rest of this script. 

############# Quick reminder for monte carlo regression.
# let's have a regression  Y~ N(a+b*x, sd)
a=5
b=0.7
x <- seq(2,20)
y_fixed <- a + b*x
plot(y_fixed~x)

# We need to add some variation, let's say sd=2
y.sim.1 <- rnorm(length(x),mean=y_fixed, sd=2) # Explain length(x)
plot(y.sim.1~x)
abline(a=5, b=0.7) # Expected relationship based on the parameters we used.

# But what is the actual parameter estimates for the regression?
y.sim.1.lm <- lm(y.sim.1 ~ x)
summary(y.sim.1.lm)  # notice parameter estimates and RSE!
confint(y.sim.1.lm) # does it include our expected values
abline(reg=y.sim.1.lm,lty=2) # estimated values based on simulated data.
# The point is, we have now done a little simulation of our regression model. 
###############


#### Chunk 4: The basic function for a linear model to generate monte carlo samples under a model, and extract a p-value.
##### Simple power analysis
# Let's write a little function for the simulation to extract the p_values for the slope
y.sim <- function(a=2,b=0.5, sd=4, x=2:20){
	  y_fixed <- a + b*x
	  y.sim.1 <- rnorm(length(x), mean=y_fixed, sd=sd)
	  y.sim.1.lm <- lm(y.sim.1 ~ x)
	  #return(summary(y.sim.1.lm))
	  return(summary(y.sim.1.lm)$coef[2,4]) # This just pulls out the p-value we want for the slope.
	  # Note 
	  }
	  
# let's talk about what we have done.
#################



# Let's look at some examples from this simulation
y.sim()
y.sim() # And again

# Now we have the tools to re-run the simulation 
n.sims = 100 # How many times do we want to simulate from this distribution.
p.lm.1 <- replicate(n.sims,y.sim()) 
p.lm.1  # all 100 (n.sims) p-values we generated
length(p.lm.1[p.lm.1 < 0.05])/n.sims  # Power -- how many sims have p values less than alpha as a proportion of number of simulations. Remember that a zero really means "less than 1/n.sims" for a p-value.

#  n.sims should generally be at least 1000, but for now it would take too much time
####################









########## 
# let's look at a null model now. This will again help us understand what we mean by alpha.
n.sims=1000
p.lm.null <- replicate(n.sims, y.sim(b=0))  
length(p.lm.null[p.lm.null < 0.05])/n.sims
# Why is this Approximately 0.05? 
# it looks like our level of alpha is what we would expect (or close to it). That is we still expect, by chance 5 out of 100 simulations to have p <0.05 just by chance, BY DEFINITION of how we specify alpha.

hist(p.lm.null) # THis distribution is approximately uniform between 0 and 1. This is what we expect when simulating under a null model.
#################



####Chunk 5: Some additional thoughts on Monte Carlo methods for power analysis

# We can use this framework to systematically alter parameters for our power analysis.

# More importantly we can do it for ANY MODEL WE WANT!!!!

#### Written after class, add to Thursday's lecture.
# This type of simulation is not only important for power analyses, but also for checking how well your model (and script) works with "known data", especially useful for complex data.





# A couple of thoughts

#1 - First, make sure that your simulated data makes sense. i.e. the values are reasonable (for instance you may not want negative values for your variables). So take a look.


# 2 There are a couple of ways to set up your explanatory variable, one is to use seq(from, to, by) as above
#  You may want to add noise to your explanatory variable  x <- rnorm(length(x),x, sd) will add some scatter
# You can directly use the observed value of your covariates, which is probably best!!!
#  The problem with using the seq() approach to modeling the explanatory variable is that your variable is spread out equally and regularly, unlike real data in many cases  (for instances if you were modeling fitness as a function of size).
# Then you can just specify the explanatory variable as a random variable itself with mean and sd (for a normally distributed data set)
### i.e.  x <- rnorm(sample_size,mean, sd)  .. and specify the y <- a + b*x as always

##3   ###### A BRIEF ASIDE ON WHAT P value you should be grabbing in your function !!!!!!! #############
# We should also think about the p-value (or whatever test statistic we want to get)
# For the above example we extracted the p-value directly associated with a particular coefficient.
# Sometimes we want something else, like the p value for the overall model
# We can extract the f statistic like
model.fstats <- summary(y.sim.1.lm)[["fstatistic"]] # For a univariate model this will be a vector of length 3

#then compute the p-value directly from the F distribution
pf(model.fstats[1], model.fstats[2], model.fstats[3], lower.tail=F) # Comparing to upper tail..


# Depending on your question, you could also consider using variance explained, RSE, etc.... So THINK HARD BEFORE STARTING!!!!

######################



###########ANOVA ANOVA ANOVA ANOVA ANOVA ANOVA ANOVA ANOVA ANOVA ANOVA ANOVA 
# Using this framework for a classic "ANOVA" power analyses
#4  this method can be used for ANOVA's if you use dummy variables for the categorial variables (you may need to change how you extract p values)

### x <- rep(c(0,1),sample_size/2) # This will create a dummy vector of 0's and 1's (representing two treatment levels)
# sample size is divided by 2 since for each increment of one on sample size you get two observations (0 & 1)
# This specifies things for balanced designs only.

# For unbalanced designs you could use
#  x <- c(rep(0, sample_size1), rep(1, sample_size2))

# The estimated b represents a co-efficient which shows the difference of group "1" from the overall mean of group "0"
# This is also how the lm() function spits it out, so it is useful for consistency.

# For a treatment with m groups, you need m-1 dummy variables
# so for three treatment levels you need 2 dummy variables (x1, and x2) 
# x1 accouts for the second level compared to the first (0's for level 1 &3, 1's for level 2)
# x2 accouts for the third level compared to the first (0's for level 1 &2, 1's for level 3)
# x1 <- c(rep(0, sample_size1), rep(1, sample_size2), rep(0, sample_size3))
# x2 <- c(rep(0, sample_size1), rep(0, sample_size2), rep(1, sample_size3))


# You can also use model.matrix(your real data) if you want to use the design matrix from your own study. Be careful of missing data though.
# One convenient way to generate simulated factors is to use gl()  # generate levels function

# Here is an example
factor.sim <- gl(n=4,k=5)  # Single factor with n=4 levels, each level has k=5 observations
factor.sim  # To show you what it looks like
design.sim <- model.matrix( ~ factor.sim)
design.sim  # To show you the design matrix, this is how linear models are constructed.
# You still need to generate a vector of the co-efficients for each of these parameters...
# you can do this with a for loop or writing a function and using replicate

sim.coefficients <- c(2, 2.5, 5, 3) # The simulated means for our 4 groups

simulated.values <- matrix(NA, nrow=5, ncol=4)
for (i in 1:nlevels(factor.sim)){	
	simulated.values[,i] <- rnorm(5, mean = sim.coefficients[i], sd = 2 )
}

simulated.values <- as.numeric(simulated.values)  # back into a vector

sim.data <- data.frame(treatment = factor.sim, simulated.values = simulated.values)
#########



#### Chunk 6: Extracting p-values for the whole model.
###### Extracting p-values from the whole model, instead of for a particular parameter.

# What would you need to do if you wanted to extract the p-value from the overall model, and not just a particular parameter as we have done here? This turns out to be a bit trickier, only because it is not a stored directly in an object. However, the information we need is.

summary(y.sim.1.lm)$fstatistic # This outputs the F value, the numerator df and denominator df for the whole model. 
# This all we need to compute a p-value (if you do not remember, go back and look at the notes on probability distributions, and particular on the family of "p" functions.)

# We will use pf(), to get the p value from the f distribution, using lower.tail=F (since for an F we want the upper tail)
p.value <- pf(summary(y.sim.1.lm)$fstatistic[1], summary(y.sim.1.lm)$fstatistic[2],summary(y.sim.1.lm)$fstatistic[3], lower.tail=F)
p.value # So you can use this instead for the model p-value.