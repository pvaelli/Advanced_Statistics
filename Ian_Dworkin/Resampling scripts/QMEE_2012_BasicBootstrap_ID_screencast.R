# Written by Ian Dworkin - last modified on October 26th 2012



# Chunk 1: the basic non-parametric bootstrap, using explicit function calls
# Chunk 2: A bias in the percentile confidence intervals produced by bootstrapping...
# Chunk 3: An example of a basic bootstrapping approach for a statistical test (for clarity)

# Read in, and set up data as usual.
setwd("/Users/ian/R/R scripts/Dll data/") 
dll.data = read.csv("dll.csv", header=TRUE) 

#dll.data = read.csv("http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv",
#   header=TRUE)

dll.data <- na.omit(dll.data)  
dll.data$temp <- as.factor(dll.data$temp)
dll.data$replicate <- as.factor(dll.data$replicate)

####Chunk 1: the basic non parametric bootstrap to generate CIs for the mean
# This function is in fact longer than it needs to be, but this is for clarity.
# It is very simple, but allows us to bootstrap the mean

BootstrapMean <- function(X=dll.data$SCT){
	x.boot <- sample(X, size=length(X), replace=T)
	mean(x.boot) }

N <- 1000 # How many iterations we are doing.

# Use replicate to do N iterations
boot.replicate <- replicate(1000, BootstrapMean() )


# Or we can do it as a for loop
stat = rep(NA, N) # initialize vector to store bootstrap means

#Set up a loop to generate a series of bootstrap samples
for (i in 1:N){  
    stat[i] <- BootstrapMean()}

# quick plot from both replicate and for loop
par(mfrow=c(1,2))    
hist(stat)
abline(v=mean(dll.data$SCT), lwd=2, col="red")
hist(boot.replicate)
abline(v=mean(dll.data$SCT), lwd=2, col="red")

# SE of estimate (parametric)
with(dll.data, sd(SCT)/sqrt(length(SCT)) )


#SE of mean based on bootstrapping
sd(stat) # remember the sd of the means of the bootstrapped values = SE
sd(boot.replicate)

 
# parametric confidence intervals based on assumptions of normality
CI_L <- with(dll.data, mean(SCT) - 1.96*(sd(SCT)/sqrt(length(SCT))))
CI_U <- with(dll.data, mean(SCT) + 1.96*(sd(SCT)/sqrt(length(SCT))))
CI_L; CI_U   # 95% CI

# bootstrapped percentile confidence intervals
quantile(stat, probs = c(0.025, 0.975)) # empirically determined confidence intervals
quantile(boot.replicate, probs = c(0.025, 0.975))




#### Chunk 10: A bias in the percentile confidence intervals produced by bootstrapping...
# BUT THIS APPROACH CAN BE BIASED. The Bias corrected and accelerated (BCa) intervals are preferred. I will make the script available demonstrating the calculations for BCa, but they are a bit fiddly. Instead we will use a powerful library, called boot.
#

library(boot)
# First we write a function that generates a vector of means of SCT for bootstrapping
# Note the "id" variable in the index. This resamples from the index of the object (vector or matrix), not the elements of the object itself!
f1 <- function(y=dll.data$SCT, id)  {mean(y[id])} 

boot.out <- boot(dll.data$SCT, f1, 3000) # calling the boot() function, for 3000 replicates
boot.out
boot.ci(boot.out, conf = 0.95, type = c("basic", "bca", "perc")) # basic and adjusted CI

# Sometimes if you do not do enough bootstrap replicates (at least as many as the number of observations) you can get an error with
# "estimated adjustment 'a' is NA". Try doing more.


# NOT BIASED in this case!!!
##############################


## Chunk 11: An example of a basic bootstrapping approach for a statistical test.
# A non-parametric bootstrap analysis to test for differences in mean # of sex comb teeth between genotypes. This is not really a great way to do it, but in Thursday's lecture I will demonstrate a more general way.
names(dll.data)
levels(dll.data$genotype)
dll.SCT = dll.data[dll.data$genotype=="Dll", "SCT"]   # notice " == "
# selecting observations for the Dll genotype only
dll.SCT = na.omit(dll.SCT)  # removing NA's
l1 = length(dll.SCT)  # length of vector of values

 # and for wild type (wt)
wt.SCT = dll.data[dll.data$genotype=="wt", "SCT"]
wt.SCT = na.omit(wt.SCT)
l2 = length(wt.SCT)

diff.SCT = mean(dll.SCT) - mean(wt.SCT)   # mean differences

      #bootstrap resampling for the means, to compare difference
N= 1000 # number of replicates
stat = numeric(N)
for (i in 1:N){  
    SCT_Boot.dll = BootstrapMean(X=dll.SCT)  
    SCT_Boot.wt = BootstrapMean(X=wt.SCT)
    stat[i] = mean(SCT_Boot.dll) - mean(SCT_Boot.wt)
    }
    
hist(stat) # look at the distribution of resampled differences in means. 
# nothing overlapping with zero

# the pseudo p-value approach is as follows
length(stat[stat<=0])/N   # i.e. how many of the resampled differences are less.. 
 # ..than or equal to zero.  empirical p < 0.0002 (for 5000 resampling events)
 mean(stat) # the mean difference between Dll and wt in #SCT is
  
# make same plot as above clearer by adding some other information
hist(stat, xlim = c(0,1), xlab = " Dll - wt") 
arrows(0.1,200,0,0, lwd =2) # arrow to point to zero, the null hypothesis
text(0.1,240,"Ho = Diff of 0",offset=0, cex = 1.3)
arrows(0.875,880,0.71,600, lwd = 2) # arrow to point to zero, the null hypothesis
text(0.85,920,"Distribution of bootstrapped differences",offset=0, cex = 1.3, col= "red")

# Remember to get an actual p-value we need to do the sampling under the null (which is why the above are just pseudo p-values). You can get this from a permutation test or from sampling with replacement from the differences from random permutations of the data.

##############

# little function to compute the sd of bootstraps using the variance ( a double check)
#datx is the matrix of bootstrapped data, index is which column to use from the matrix
sd.boot.test <- function(datx,index) {
  var.boot <-sum((datx[index,] - mean(datx[index,]))^2)/(N-1)
  print("bootstrap variance & SD", quote=F)
  return(c(var.boot, sqrt(var.boot)))}   