# October 23rd 2014
# ZOL851  In Class exercise for resampling.

# It is often the case that a genetic or environemntal (i.e. temperature or nutritional) stressor/perturbation will not only influence the mean value of a trait, but the variance (both environmental and genetic) for the trait as well. There is a considerable literature that examines how and why this occurs. However, approropraite measures to assess variation are somewhat trickier than for estimating "mean" effect sizes. 

# A commonly used approach is to use the coefficient of variation
CV <- function(x) {sd(x)/mean(x)}

# However while the intepretation of CV is pretty straight forward, statistical inference for it is not. (There are in fact better approaches based on Levene's statistic among others).

# For the dll.data dataset I want you to assess how much if any additional variance (as measured using CV) the mutants have compared with the wild types. Pick one line (only one) AND one temperture (25 or 30).

# I suggest calculating CV for each combination of line and genotype (we will only use data at one temperature and ignore tarsus length for this)

 

dll.data = read.csv("http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv", header=TRUE)
 
dll.data = na.omit(dll.data)
dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
dll.data$genotype <-relevel(dll.data$genotype, ref="wt")
dll.data$tarsus.scaled <- scale(dll.data$tarsus, center=TRUE, scale=FALSE)

with(dll.data, table(genotype, line, temp))

#1 - pick one line and temperature combination and generate a data subset to use. You should try to use something with reasonable sample sizes for each group.

dll.sub <- subset(dll.data, line=="line-21" & temp=="25")

#2 Compute the CV for each group.

wt <- dll.sub[which(dll.sub$genotype=="wt"),]
Dll <- dll.sub[which(dll.sub$genotype=="Dll"),]
with(wt, CV(SCT))
with(Dll, CV(SCT))

with(dll.sub[dll.sub$genotype=="wt",], CV(SCT))
with(dll.sub[dll.sub$genotype=="Dll",], CV(SCT))

#3 Use a np bootstrap to compare CV between each group.

N <- 1000

BootCV <- function(x){
	var1 <- sample(x, length(x), replace=TRUE)
	CV(var1)
}

CVwt <- with(dll.sub[dll.sub$genotype=="wt",], replicate(N, BootCV(SCT)))
hist(CVwt, col="blue")

CVDll <- with(dll.sub[dll.sub$genotype=="Dll",], replicate(N, BootCV(SCT)))
hist(CVDll, add=T, col="green")

#wt.npboot <- replicate(N, with(dll.sub[dll.sub$genotype=="wt",], sample(SCT, length(SCT), replace=T,)))
#Dll.npboot <- replicate(N, with(dll.sub[dll.sub$genotype=="Dll",], sample(SCT, length(SCT), replace=T,)))

#with(dll.sub[dll.sub$genotype=="wt",], CV(SCT))
#CV(wt.npboot)

#with(dll.sub[dll.sub$genotype=="Dll",], CV(SCT))
#CV(Dll.npboot)



#4 Can you think of a way to implement a permutation test for this?

N <- 1000
wt.permutation <- function {
	replicate(N, with(dll.sub[dll.sub$genotype=="wt",], sample(SCT, length(SCT), replace=F,)))
	
Dll.permutation <- replicate(N, with(dll.sub[dll.sub$genotype=="Dll",], sample(SCT, length(SCT), replace=F,)))

with(dll.sub[dll.sub$genotype=="wt",], CV(SCT))
CV(wt.npboot)
CV(wt.permutation)

CV(Dll.permutation)

# 5 clearly we want to test this for the whole dataset. How might you approach this same question to address it for all lines? Write down steps.



#### Note a more common approach is to use Levene's test as you can deviates among individuals (see Dworkin 2005)
### Below I have written a general purpose function to generate Levene's deviates
LeveneDeviates <- function(y, group, med=TRUE, log_trans=TRUE) {
    
    #log transform data?
    if (log_trans==TRUE)
        y = log(y)
    
    # allows for use of mean or median as measure of central tendency
    if (med==TRUE)
        meds <- tapply(y, group, median, na.rm=TRUE)
    else 
        meds <- tapply(y, group, mean, na.rm=TRUE) 
    
    # calculates deviates for each observation from a measure of central tendency for a "group"
    abs(y - meds[group])}
    
lev_stat <- with(dll.data, 
    LeveneDeviates(y = SCT, group = genotype:temp:line, med=TRUE, log_trans=FALSE))

ls.lm <- lm(lev_stat ~ genotype*temp*line, data=dll.data)
anova(ls.lm)  # of course this ANOVA is not particularly valid given the lack of normality, use bootstrapping or permutation.


