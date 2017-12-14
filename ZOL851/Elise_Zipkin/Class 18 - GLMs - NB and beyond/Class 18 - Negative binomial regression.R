#Want to estimate the number of water bird species at each location i in three different water depths: "Shallow", "Medium",
#"Deep".  We also think that turbidity may negatively influence the number of species (e.g., more species at lower values).
#Water birds tend to aggregate so we expect that they may not look Poisson distributed

#Load the libraries
require(foreign)
require(ggplot2)
require(MASS)

#Read in the data and look at a summary
birds = read.csv("waterbirddata.csv")

summary(birds)

#Plot a histogram of the data by depth (and overall)
ggplot(birds, aes(species, fill = depth)) +
  geom_histogram(binwidth=1) +
  facet_grid(depth ~ ., margins=TRUE, scales="free")

#Look at a summary - mean and standard deviation by depth
with(birds, tapply(species, depth, function(x) {
  sprintf("Mean (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))


#What would we expect from a poisson distribution?

#patric# We would expect mean = variance in poisson, which is why we need negative binomial

#Create a dummy variable called "pois" with species values that are Poissonly distributed
#with the observed mean

lambda = tapply(birds$species, birds$depth, mean)

for (i in 1:dim(birds)[1]) {
  a = pmatch(birds$depth[i],names(lambda))
birds$pois[i] = rpois(1,lambda[a])
}

#Plot a histogram of the data by depth (and overall)
ggplot(birds, aes(pois, fill = depth)) +
  geom_histogram(binwidth=1) +
  facet_grid(depth ~ ., margins=TRUE, scales="free")

#Run the above code several times to get an idea of what Poissonly distributed data with means of 
#lambda would look like.  Do the real data look like they have a Poisson distribution?

###########################################################################################

#Do an anlaysis of the data using a negative binomial regression that accounts for both depth and turbidity.
#The glm.nb function is in the MASS package
NB.regression <- glm.nb(species ~ turbidity + depth, data = birds)
summary(NB.regression)

#What do we find?
#Pull out the regression coefficients and their confidence intervals
est <- cbind(Estimate = coef(NB.regression), confint(NB.regression))

#Look at them and also on the inverse log scale
est
exp(est)

#The variable turbitidity has a coefficient of -0.006, which is statistically significant. This means 
#that for each one-unit increase in turbitidy, the expected log count of the number of species decreases by 0.006.
#The indicator variable shown as depthMedium is the expected difference in log count between group 2 and 
#the reference group (depth=1 - deep). The expected log count for level 2 of depth is 0.44 lower than the
#expected log count for level 1. 

####################################################################################

#Let's examine the predicted values
#Create new datasets with values of depth and turbitidy and then use the predict command to calculate 
#the predicted number of species.  First, we can look at predicted counts at each depth while holding 
#turbitidy at its mean. To do this, we create a new dataset with the combinations of depth and turbitidy 
#for which we would like to find predicted values, then use the predict command.

newdata1 <- data.frame(turbidity = mean(birds$turbidity),
                       depth = factor(1:3, levels = 1:3, labels = levels(birds$depth)))
newdata1$phat <- predict(NB.regression, newdata1, type = "response")
newdata1

#What does newdata1 tell us?  What if we had standardized turbitidy in our analysis?  Would our parameter estimates
#change?  In what ways?

#############################################

#Now we will obtain the mean predicted number of species for values of turbitidy across the 
#entire range for each depth and graph these.

newdata2 <- data.frame(
  turbidity = rep(seq(from = min(birds$turbidity), to = max(birds$turbidity), length.out = 100), 3),
  depth = factor(rep(1:3, each = 100), levels = 1:3, labels =
                  levels(birds$depth)))

newdata2 <- cbind(newdata2, predict(NB.regression, newdata2, type = "link", se.fit=TRUE))

newdata2 <- within(newdata2, {
  species <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

#The graph shows the expected number of species across the range of observed turbitidy values, 
#for each depth along with 95 percent confidence intervals. Note that the lines are not straight 
#because this is a log linear model, and what is plotted are the expected values, not the log of the expected values.

ggplot(newdata2, aes(turbidity, species)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = depth), alpha = .25) +
  geom_line(aes(colour = depth), size = 2) +
  labs(x = "Turbidity", y = "Predicted Species Richness")

