
# ZOL851 - The absolute basics of performing Monte Carlo simulations in R
# Last modified October 4th 2012 ID
 
#  What do I mean by simulate?

 
norm.simulated <- rnorm(n = 100, mean  = 5, sd = 2)
par(mfrow=c(3,1))
plot(norm.simulated)
hist(norm.simulated)
hist(norm.simulated, freq=F)
curve(dnorm(x, mean=5, sd=2), from = 0, to = 10, add=T, col="red")

#Let's look at some summary statistics from this simulation
mean(norm.simulated)
sd(norm.simulated)
# While the mean and sd are close to the values we inputted, they are not exactly 5 & 2.

# let's repeat this experiment a few times
norm.simulated1 <- rnorm(n = 100, mean  = 5, sd = 2)
norm.simulated2 <- rnorm(n = 100, mean  = 5, sd = 2)
norm.simulated3 <- rnorm(n = 100, mean  = 5, sd = 2)
norm.simulated4 <- rnorm(n = 100, mean  = 5, sd = 2)

par(mfrow=c(2,2))
hist(norm.simulated1)
hist(norm.simulated2)
hist(norm.simulated3)
hist(norm.simulated4)
mean(norm.simulated1)
mean(norm.simulated2)
mean(norm.simulated3)
mean(norm.simulated4)

# All close to, but not exactly equal to 5.

# So for each round (iteration) of the simulation, we have effectively sampled 100 observations from a population with a normal distribution with mean =5 and sd =2. If we wanted to change the number of observations per iteration we simply change n. Similarly we could change the values for mean (or make it some function like a line), the standard deviation and the type of probability distribution we are sampling from.

#### more efficient performing the simulation 
# Of course performing this one iteration at a time is not particularly useful, so instead we find a way to repeat it efficiently.

# One way to do this in R is to use the replicate function, which just repeats the function call n times. Notice that there is the n=4 in the replicate which refers to how many times to repeat the call to the function (in this case the rnorm). the n=100 in rnorm is how many samples are drawn from the population (distribution).

norm.sim.all.2 <- replicate(n=4, rnorm(n=100, mean = 5, sd =2)) 
# And if we want more than 4 replications n=whatever we need for replicate
par(mfrow=c(2,2))
apply(X = norm.sim.all.2, MARGIN = 2, FUN=hist)
apply(X = norm.sim.all.2, MARGIN = 2, FUN=mean)
apply(X = norm.sim.all.2, MARGIN = 2, FUN=sd)


#let's do this again, but repeat the sampling process (number of iterations) 2000 times.
norm.sim.all.3 <- replicate(n=2000, rnorm(n=100, mean = 5, sd =2)) 

# We can ask about the spread of the means from our simulated trials
sd(apply(X = norm.sim.all.3, MARGIN = 2, FUN = mean)) # Spread of our means

# compare this to the expected standard error
2/sqrt(100) # Because this sd/sqrt(sample size) (the standard error of the mean) approximates the standard deviation of the sampling distribution

# The simulated sampling distribution for the means
par(mfrow=c(1,1))
hist(apply(X = norm.sim.all.3, MARGIN = 2, FUN = mean),
    main = "simulated sampling distribution for means",
    xlab = "estimated mean") 

# for loops  
### We can also run such simulations using for loops, which can be more efficient when the number of simulations (or parameters being estimated) is large

# How many iterations of the simulation do we want to do? We will call this "N"
N <- 2000

# First we initialize a variable to store the means
simulated_means <- rep(NA, N) # Create an empty vector of length N

head(simulated_means)  # Filled with "missing data" until we fill it

# Now we use the for loop to repeat this from i=1 until i=N (which is 2000 in this case )
for (i in 1:N){
	sim_data <- rnorm(n = 100, mean  = 5, sd = 2)
	simulated_means[i] <- mean(sim_data)
	rm(sim_data) # While this gets over-written each iteration of the loop, 
	   #we want to remove it after the last iteration.
}

# the vector "simulated_means" now contains all of the means from the N iterations of the simulation.
hist(simulated_means)
sd(simulated_means)


####
# Let's repeat this experiment with lower sample size (say 25)
norm.sim.all.4 <- replicate(n=2000, rnorm(n=25, mean=5,sd=2))
#apply(norm.sim.all.4,2,sd)
sd(apply(norm.sim.all.4,2,mean))
hist(apply(norm.sim.all.4,2,mean))

# How about if we have much larger sample size (say 1000)
norm.sim.all.5 <- replicate(n=2000, rnorm(1000,5,2))
#apply(norm.sim.all.5,2,sd)
sd(apply(norm.sim.all.5,2,mean))
hist(apply(norm.sim.all.5,2,mean))

# Do you see a pattern for the precision of the mean?
sd(apply(norm.sim.all.2,2,mean)) # n=100
sd(apply(norm.sim.all.4,2,mean)) # n=25
sd(apply(norm.sim.all.5,2,mean)) # n=1000

# This hopefully reinforces why the sample size influences the estimated the standard error so clearly!


# We have just learned about monte carlo methods. When we perform simulations from that distributions, as our n (sample size) increase we will get closer to the expected value (which in this case we provided). That is the law of large numbers.


# Instead of simulating just a mean and sd, let's specify something more interesting.. 
# How about a line/regression!!!!!

# let's have a regression  Y~ N(a+b*x, sd)
a=5
b=0.7
x <- seq(2,20) # our predictor variable
y_fixed <- a + b*x  # the deterministic part of our model

par(mfrow=c(2,1))
plot(y_fixed~x)  # The deterministic relationship.

# But our response is randomly sampled (conditioned on x), so we need to be able to incorporate this random variation. So let's say sd=2.5
y.sim.1 <- rnorm(length(x), mean=y_fixed, sd=2.5) 
    #length(x) just looks at the number of elements in the vector
# So we have now simulated values of our response variable y given a deterministic and stochastic component to the model!

plot(y.sim.1~x)
abline(a=5, b=0.7, lwd=2) # Expected relationship based on the "known" parameters we used.

# But what are the actual parameter estimates for the regression based on the simulated data (for y)?
# run a linear regression
y.sim.1.lm <- lm(y.sim.1 ~ x)
summary(y.sim.1.lm)  # notice parameter estimates and RSE!
confint(y.sim.1.lm) # does it include our expected values
abline(reg=y.sim.1.lm, lty=2, col="red", lwd=2) # estimated values based on simulated data.
# The point is, we have now done a little simulation of our regression model. 

