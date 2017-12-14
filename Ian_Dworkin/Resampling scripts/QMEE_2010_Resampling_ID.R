# Written by Ian Dworkin - last modified on October 3rd 2011


# Chunk 1: read in data as usual
# Chunk 2: Overview of computation for resampling.
# Chunk 3: Fitting the basic model with lm() before attempting the resampling
# Chunk 4: A basic permutation test using a for() loop
# Chunk 5: A basic permutation test using an explicit functional call, and replicate() 
# Chunk 6: Distribution of p values under permutation
# Chunk 7: Some notes on the computation of permuted data sets.
# Chunk 8: the basic non-parametric bootstrap, using a for loop
# Chunk 9: the basic non-parametric bootstrap, using explicit function calls
# Chunk 10: A bias in the percentile confidence intervals produced by bootstrapping...
# Chunk 11: An example of a basic bootstrapping approach for a statistical test (for clarity)


# Chunk 1: Read in, and set up data as usual.
setwd("/Users/ian/R/R scripts/Dll data/") 
dll.data = read.csv("dll.csv", header=TRUE) 
dll.data <- na.omit(dll.data)  
dll.data$temp <- as.factor(dll.data$temp)
dll.data$replicate <- as.factor(dll.data$replicate)


### Chunk 2: Overview of computation for resampling  ############################

 # For performing the randomization/permutation test (or bootstrap) there are two important pieces of code to know. 1) the sample() function and 2) how to perform an iterative loop (with a for(){} loop.
 
 #  the function that performs the resampling is sample()
 #  sample(x, size, replace=False)
 #  where x is the vector to be resampled from, size is the length of the vector to be generated, and
 #  replace is whether or not sampling with replacement is to take place
 
 # the default value of size is length(x) which is what we want for both permutations and bootstrapping.


   # iterative loops using for (){}
   # for (i in 1:N) { do something here!!!}  # where i is the index variable, and N is the number of iterations.
   
   
   # If all of this code gives you a headache, just follow below, and remember that all you need to do is change N for number of resampling events, and change your model and test statistic.
########




####### Chunk 3: Examining the basic model we are fitting using lm (prior to getting to the resampling)
dll.anova = lm(SCT ~ genotype, data=dll.data) 
summary(dll.anova) # look at model summary statistics
summary.dll.anova = summary(dll.anova) 
summary.dll.anova   # prints this object to the screen
dll.anova.fstat  = summary.dll.anova$fstat # just pulls out the F statistic
dll.anova.fstat   # Fstats = 60.22
dll.anova.fstat1 = dll.anova.fstat[1]



#### Chunk 4:  A basic permutation test
# and here is the actual code for the permutation test
N = 500 # # of replicate resampling events
fstats = numeric(N) # initializes a numeric vector to input the resampled F stats
for (i in 1:N){  
    dll.resample = lm(sample(SCT, replace=F) ~ genotype, data=dll.data) #permutes SCT, and then
    # runs the model  
    fstats[i] = summary(dll.resample)$fstat[1] # place all of the fstats in this vector 
    }

hist(fstats) #let's remind ourselves again about what F really is.

# Pretty plot to look at the distribution of the F stastistics from permutation.
par(mfrow=c(2,1))
hist(fstats, xlim = c(0, 65),ylim = c(0,1), xlab = "F Statistic", main = " F Statistics: Empirically derived via randomization", freq=F )
arrows(dll.anova.fstat1, 0.6, dll.anova.fstat1, 0, lwd = 3)


# plot of F values based on F distribution

f.theor <- rf(1000,1,1946)
hist(f.theor, xlim = c(0,65),ylim = c(0,1), xlab = "F Statistic", freq=F,
main = " F values from F-distribution (1,1946)")
curve(df(x, 1,1946), 0, 60, add=T)


length(fstats[fstats>= 60.22])/N   # gives the empirical P value, 
 # input the observed F for your model where the 60.22 is.
 # in this case none of the permuted values is greater than the observed
 # so the empirical p-value <0.001 (for 1000 permutations)

length(fstats[fstats>= dll.anova.fstat1])/N # alternatively you do not need to directly
 # input the f value into this.
# Plot of  the distribution of F statistics for the permuted data sets.



##### Chunk 5: The permutation test using an explicit functional call.
##### Explicit functional call for the permutation test, and repeated using replicate()
# We can (and should) use a function like replicate() and write the permutation test as an explicit function
# Below is a function specifically for this data set. It is actually pretty straightforward to write this more generally.
permute.function <- function(){
	dll.resample = lm(sample(SCT, replace=F) ~ genotype, data=dll.data) #permutes SCT, and then
    # runs the model  
    fstats = summary(dll.resample)$fstat[1] # place all of the fstats in this vector 
    }

permute.N <- replicate(N, permute.function())

hist(permute.N)



#### Chunk 6:  Distribution of p values under permutation... Make sure you understand this.
# We can also use this to look at the distribution of p values under this model
permute.function.p <- function(){
	dll.resample = lm(sample(SCT, replace=F) ~ genotype, data=dll.data) #permutes SCT, and then
    # runs the model  
    anova(dll.resample)[1,5]
     # place all of the fstats in this vector 
    }
permute.N.p <- replicate(N, permute.function.p())

hist(permute.N.p)




#############
# Chunk 7: Some notes on the computation of permuted data sets.


 # For univariate tests use sample on the dependent variable only. For models with a single explanatory variable it will not matter, but for models with a number of explanatory variables if you use sample() on the explanatory variables it will break up the intra-observation correlation structure, which is not generally something you want.
 
 #  For multivariate tests, you can not use sample() quite as easily. For a simple model where you have multiple dependent variables, and a single explanatory variable use sample() on the explanatory variable. For complex models:
 
 # the basic idea is for matrix "a" with m rows
 # a[sample(m,m),]  # allows indexing

#####################




### There is a package (coin), which apparently implements a variety of permutation tests. I have not tried it.




####Chunk 8: the basic non parametric bootstrap
#This program selects 1000 bootstrap samples from your data
# This code provides the basic framework for boostrapping which we can use later.
SCT = dll.data[,"SCT"]
n = length(SCT) # length (# observations) of vector
N = 100 # # of replicates -creates an empty vector to store results
#the elements of the vector will be numbered from 1 to N

stat = numeric(N)  # Generates an empty vector where we will put the bootstrapped values


#Set up a loop to generate a series of bootstrap samples
for (i in 1:N){  
    #bootstrap sample counterparts to observed samples are denoted with "B"
    SCT_Boot= sample (SCT, n, replace=T) # notice that replace=T!!! 
    #SCT_Boot contains the vectors of bootstrapped  samples
    stat[i] = mean(SCT_Boot) # the mean for each set of resampled vectors
    }
hist(stat)
mean(SCT);   mean(stat) # compare bootrapped and observed means
sd(SCT) # standard deviation
sd(SCT)/sqrt(length(SCT)) # SE of estimate (parametric) 


#SE of mean based on bootstrapping
sd(stat) # remember the sd of the means of the bootstrapped values = SE
  
# parametric confidence intervals based on assumptions of normality
CI_L <- mean(SCT) - 1.96*(sd(SCT)/sqrt(length(SCT)))
CI_U <- mean(SCT) + 1.96*(sd(SCT)/sqrt(length(SCT)))
CI_U; CI_L   # 95% CI

# bootstrapped percentile confidence intervals
quantile(stat, probs = c(0.025, 0.975)) # empirically determined confidence intervals


####### Chunk 9: bootstrap with explicit functional call
# This function is in fact longer than it needs to be, but this is for clarity.
bootstrap.function <- function(X=SCT){
	x.boot <- sample(X, size=length(X), replace=T)
	mean(x.boot) }

boot.replicate <- replicate(1000, bootstrap.function() )


#### Really quick one liner approach...
bootstrap.call <- replicate(1000, mean(sample(SCT,size=length(SCT), replace=T))) # Samples from a distribution of data...




#### Chunk 10: A bias in the percentile confidence intervals produced by bootstrapping...
# BUT THIS APPROACH CAN BE BIASED.
#
# USE THE boot library

library(boot)
f1 <- function(SCT, id)  {mean(SCT[id])} # A function that generates a vector of means of SCT for bootstrapping
boot.out <- boot(SCT, f1, 3000) # calling the boot() function, for 3000 replicates
boot.out
boot.ci(boot.out, conf = 0.95, type = c("basic", "bca", "perc")) # basic and adjusted CI

# Sometimes if you do not do enough bootstrap replicated you can get an error with
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
    #bootstrap sample counterparts to observed samples are denoted with "B"
    SCT_Boot.dll= sample (dll.SCT, l1, replace=T) # notice that replace=T!!! 
    #SCT_Boot contains the vectors of bootstrapped  samples
    SCT_Boot.wt= sample (wt.SCT, l2, replace=T)
    stat[i] = mean(SCT_Boot.dll) - mean(SCT_Boot.wt)
    }
hist(stat) # look at the distribution of resampled differences in means. 
# nothing near zero, suggesting a highly significant difference

# the formal approach is as follows
length(stat[stat<=0])/N   # i.e. how many of the resampled differences are less.. 
 # ..than or equal to zero.  empirical p < 0.0002 (for 5000 resampling events)
 mean(stat) # the mean difference between Dll and wt in #SCT is
  
# make same plot as above clearer by adding some other information
hist(stat, xlim = c(0,1), xlab = " Dll - wt") 
arrows(0.1,400,0,0, lwd =2) # arrow to point to zero, the null hypothesis
text(0.1,440,"Ho = Diff of 0",offset=0, cex = 1.3)
arrows(0.875,880,0.71,600, lwd = 2) # arrow to point to zero, the null hypothesis
text(0.85,920,"Distribution of bootstrapped differences",offset=0, cex = 1.3, col= "red")



##############

# little function to compute the sd of bootstraps using the variance ( a double check)
#datx is the matrix of bootstrapped data, index is which column to use from the matrix
sd.boot.test <- function(datx,index) {
  var.boot <-sum((datx[index,] - mean(datx[index,]))^2)/(N-1)
  print("bootstrap variance & SD", quote=F)
  return(c(var.boot, sqrt(var.boot)))}   