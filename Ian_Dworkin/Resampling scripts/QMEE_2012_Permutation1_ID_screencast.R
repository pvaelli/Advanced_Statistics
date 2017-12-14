# Written by Ian Dworkin - last modified on October 25th 2012


# Chunk 1: read in data as usual
# Chunk 2: Overview of computation for resampling.
# Chunk 3: Fitting the basic model with lm() before attempting the resampling
# Chunk 4: A basic permutation test using a for() loop
# Chunk 5: A basic permutation test using an explicit functional call, and replicate() 
# Chunk 6: Distribution of p values under permutation
# Chunk 7: Some notes on the computation of permuted data sets.

# Chunk 1: Read in, and set up data as usual.
setwd("/Users/ian/R/R scripts/Dll data/") 
dll.data = read.csv("dll.csv", header=TRUE) 

#dll.data = read.csv("http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv",
#   header=TRUE)
dll.data <- na.omit(dll.data)  
dll.data$temp <- as.factor(dll.data$temp)
dll.data$replicate <- as.factor(dll.data$replicate)


### Chunk 2: Overview of computation for resampling  ############################

 # For performing the randomization/permutation test (or bootstrap) there are two important pieces of code to know. 1) the sample() function and 2) how to perform an iterative loop (with a for loop or replicate).
 
 #  the function that performs the resampling is sample()
 #  sample(x, size, replace=False)
 #  where x is the vector to be resampled from, size is the length of the vector to be generated, and
 #  replace is whether or not sampling with replacement is to take place
 
 # the default value of size is length(x) which is what we want for both permutations and bootstrapping.

   # If all of this code gives you a headache, just follow below, and remember that all you need to do is change N for number of resampling events, and change your model and test statistic.
########




####### Chunk 3: Examining the basic model we are fitting using lm (prior to getting to the resampling)
dll.anova = lm(SCT ~ genotype, data=dll.data) 
summary(dll.anova) # look at model summary statistics
 
dll.anova.fstat  = summary(dll.anova)$fstat # just pulls out the F statistic
dll.anova.fstat   
dll.anova.fstat1 = dll.anova.fstat[1]

anova(dll.anova)

#### Chunk 4:  A basic permutation test
# and here is the actual code for the permutation test
N = 500 # # of replicate resampling events


PermuteFunction <- function(y=dll.data$SCT, x=dll.data$genotype){
	model.resample = lm(sample(y, replace=F) ~ x) 
	#permutes response, then runs model
   fstats = summary(model.resample)$fstat[1] # place all of the fstats in this vector 
   return(fstats)}


# Try calling it a few times.
PermuteFunction()

# To perform the iterations of the permutations we can use either replicate()
permute.N <- replicate(N, PermuteFunction())
hist(permute.N)

# or a for loop (useful for big data or lots of iterations)
fstats = numeric(N) # initializes a numeric vector to input the resampled F stats

for (i in 1:N) {
	fstats[i] <- PermuteFunction()}

hist(fstats) 

#let's remind ourselves again about what F really is.

# Pretty plot to look at the distribution of the F stastistics from permutation.
par(mfrow=c(2,1))
hist(fstats, xlim = c(0, 65),ylim = c(0,1), xlab = "F Statistic", main = " F Statistics: Empirically derived via randomization", freq=F )
arrows(dll.anova.fstat1, 0.6, dll.anova.fstat1, 0, lwd = 3)


# plot of F values based on F distribution

f.theor <- rf(1000,1,1946)
hist(f.theor, xlim = c(0,65),ylim = c(0,1), xlab = "F Statistic", freq=F,
main = " F values from F-distribution (1,1946)")
curve(df(x, 1,1946), 0, 60, add=T)


length(fstats[fstats>= dll.anova.fstat1])/N   # gives the empirical P value, 
 # in this case none of the permuted values is greater than the observed
 # so the empirical p-value <0.001 (for 1000 permutations)


#### Chunk 6:  Distribution of p values under permutation... Make sure you understand this.
# We can also use this to look at the distribution of p values under this model
PermutePvals <- function(){
	dll.resample = lm(sample(SCT, replace=F) ~ genotype, data=dll.data) 
	anova(dll.resample)[1,5]}
permute.N.p <- replicate(N, PermutePvals())

hist(permute.N.p)




#############
# Chunk 7: Some notes on the computation of permuted data sets.


 # For univariate tests use sample on the dependent variable only. For models with a single explanatory variable it will not matter, but for models with a number of explanatory variables if you use sample() on the explanatory variables it will break up the intra-observation correlation structure, which is not generally something you want.
 
 #  For multivariate tests, you can not use sample() quite as easily. For a simple model where you have multiple dependent variables, and a single explanatory variable use sample() on the explanatory variable. For complex models:
 
 # the basic idea is for matrix "a" with m rows
 # a[sample(m,m),]  # allows indexing

#####################




### There is a package (coin), which apparently implements a variety of permutation tests. I have not tried it.

