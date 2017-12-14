### libraries
require(bbmle)

# Input and check data as normal
#dll.data = read.csv("http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv",
#   header=TRUE)  
setwd("/Users/ian/R/R scripts/Dll data/")  
dll.data = read.csv("dll.txt", header=TRUE)   #data frame input
##dll.data = read.csv("http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv", header=TRUE)

dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
dll.data$genotype <- relevel(dll.data$genotype, "wt") 
dll.data <- na.omit(dll.data)
dll.data$tarsus.centered <- scale(dll.data$tarsus, scale=F)



# We know from our previous examinations of this data set that there might be quite different variances for the two different genotype (the wild-type and the mutant Dll.)

CV_by_line <- with(dll.data, 
    tapply(SCT, INDEX=list(genotype,line), function(x) sd(x)/mean(x)))

plot(CV_by_line[1,], CV_by_line[2,], 
    pch=20, xlim=c(0.05, 0.25), ylim=c(0.05, 0.25),
    xlab="CV for wt", ylab="CV for Dll mutant")
abline(a=0, b=1, col="grey")    

# We can also focus on the residuals of the models.
linear_model_2 <- lm(SCT ~ tarsus.centered + genotype + temp + genotype:temp, data=dll.data)

# We can look at the residual variances broken down by genotype and temp 
tapply(resid(linear_model_2), INDEX=list(dll.data$genotype, dll.data$temp), sd)


# This is definitely suggestive of more "unexplained" variation for samples from the mutant (Dll genotype). Not only does this violate the assumption of identically distributed distributions of residual variation across groups, but it is potentially really interesting biologically. 


# That is, the mutant not only influences the trait mean, but may also increase the phenotypic variance of a trait, which relates directly to evolutionary theories of canalization and robustness.

# So how might we evaluate this?

# A traditional approach (see Dworkin 2005) might use a variant of Levene's test in combination with either a permutation test or a non-parametric (random/pairs) bootstrap to evaluate the differences. However that requires us to seperate out the model of changes in the "means" (location effects) with changes in the variances (scale effects). What if we want to model it all together?

# We want to model a situation where we have heterogeneity of variances for the two different genotypes modeled. For the moment I will focus just on the genotypic effects, but it is clear that we could also consider temperature (and maybe the interaction between genotype and temperature)

design_matrix <- model.matrix( ~ tarsus.centered + genotype + temp + genotype:temp, data=dll.data)
head(design_matrix)


# 2 - write support function/ NLL calculator function

# our parameters (b0-b4) are for the location effects (i.e. treatment contrasts, and slope of the relationship with tarsus).

# Our s0 and s1 are instead going to be analogous to our location effects but for the variances instead.

# s0 will be the residual variance for the first level of genotype (which is wt)
# s1 will be the "treatment contrast" for the residual variance for the mutant (dll) genotypic effect. THat is how much greater (or potentially less than) is the unexplained variance for the dll mutant over and above that for the wild-type.
# Since the genotypic effects are the third column of the design mantrix we will multiply s1*X[,3] for this

SE <- function(x){ 
	sd(x)/sqrt(length(x))}
	
glm_MLE_support <- function(b0 = mean(dll.data$SCT), 
     b1=0, b2=0, b3=0, b4=0, 
     s0= 1.5, s1=0.25, X=design_matrix, y = dll.data$SCT, sigma = SE(dll.data$SCT)) {
 
    deterministic_part <- b0*X[,1] + b1*X[, 2] + b2*X[, 3] + b3*X[,4] + b4*X[,5]
    random_part <- s0*X[,1] + s1*X[,3]
    -sum(dnorm( y, mean = deterministic_part, sd =  random_part, log=T))}




#4 try a fit	
glm_MLE_fit_1 <- mle2(glm_MLE_support, 
            start=list(b0 = 0, b1 = 0, b2 = 0, 
            b3 =0, b4=0, s0 =SE(dll.data$SCT), s1=0 ), method="BFGS")

summary(glm_MLE_fit_1)
-logLik(glm_MLE_fit_1)

prof_model_1 <- profile(glm_MLE_fit_1)
plot(prof_model_1)
confint(prof_model_1)

# What does this tell us? The "treatment contrast" for the unexplained variance for the mutant genotype is positive (with 95% CIs that do not include 0), which suggests that this may be an important effect.

# I am not sure what a sensible effect size might be right off the bat, but we can at least look at how much proportionally larger the unexplained variance is
 
sum(coef(glm_MLE_fit_1)[6:7])/coef(glm_MLE_fit_1)[6] # ~ 1.44 times as large

# I think it is reasonable to assess this model VS a model where we do not include different variance terms.

glm_MLE_support_2 <- function(b0 = mean(dll.data$SCT), 
     b1=0, b2=0, b3=0, b4=0, 
     X=design_matrix, y = dll.data$SCT, sigma = sd(dll.data$SCT)) {
 
    deterministic_part <- b0*X[,1] + b1*X[, 2] + b2*X[, 3] + b3*X[,4] + b4*X[,5]
    random_part <- sigma
    -sum(dnorm( y, mean = deterministic_part, sd =  random_part, log=T))}

glm_MLE_fit_2 <- mle2(glm_MLE_support_2, 
            start=list(b0 = 0, b1 = 0, b2 = 0, 
            b3 =0, b4=0, sigma =SE(dll.data$SCT) ), method="BFGS")

# As a reminder, this fit
summary(glm_MLE_fit_2)

# is the same as 
summary(linear_model_2)

# But what we are interested in is how the fit changes (as well as changes in coefficents and CIs)
-logLik(glm_MLE_fit_2)
-logLik(glm_MLE_fit_1)

# We can compute the LRT, but we need to consider if it is valid! 
Likelihood_ratio <- -2*(logLik(glm_MLE_fit_2) - logLik(glm_MLE_fit_1))
pchisq(Likelihood_ratio[1], df = 1,lower.tail=F) # highly significant.



# We may want to use a parametric bootstrap to test how robust this conclusion is to the assumptions. However, this is pretty small p-value, and in combination with the CIs and effect size I would at first blush argue this is probably an important biological effect and consistent with expectations under most models of canalization and robustness.

 # really ugly code I wrote last night to compute parametric bootstrap!
ParBootMLE <- function(Model1 = glm_MLE_fit_1, RestrictedModel = glm_MLE_fit_2, X = design_matrix){
	
	# to get simulated response under the restricted model (all in the fact that there is one sigma term)
	coef_restricted <- matrix(coef(RestrictedModel)[1:5], ncol=1)
	simulated_response <- rnorm(nrow(X), mean = X %*% coef_restricted, sd=coef(RestrictedModel)[6])
	
	# support function for restricted model
	glm_MLE_support_restr <- function(b0 , b1, b2, b3, b4, 
	y = simulated_response, sigma ) {
 
    deterministic_part <- b0*X[,1] + b1*X[, 2] + b2*X[, 3] + b3*X[,4] + b4*X[,5]
    random_part <- sigma
    -sum(dnorm( y, mean = deterministic_part, sd =  random_part, log=T))}

    # support function for full model
   glm_MLE_support_full <- function(b0, b1, b2, b3, b4, 
     s0, s1, y = simulated_response) {
    deterministic_part <- b0*X[,1] + b1*X[, 2] + b2*X[, 3] + b3*X[,4] + b4*X[,5]
    random_part <- s0*X[,1] + s1*X[,3]
    -sum(dnorm( y, mean = deterministic_part, sd =  random_part, log=T))}

    # MLE fit under restricted model
	glm_MLE_fit_restr <- mle2(glm_MLE_support_restr, 
            start=list(b0 = coef(RestrictedModel)[1], b1 = coef(RestrictedModel)[2], 
            b2 = coef(RestrictedModel)[3], b3 = coef(RestrictedModel)[4], 
            b4=coef(RestrictedModel)[5], sigma =coef(RestrictedModel)[6]), 
            method="BFGS")
            
    # MLE fit under full model        
    glm_MLE_fit_full <- mle2(glm_MLE_support_full, 
            start=list(b0 = coef(RestrictedModel)[1], b1 = coef(RestrictedModel)[2], 
            b2 = coef(RestrictedModel)[3], b3 = coef(RestrictedModel)[4], 
            b4=coef(RestrictedModel)[5], s0 =coef(RestrictedModel)[6], s1=0), method="BFGS")
    
    # compute Likelihood ratio and return numeric        
    LR <- -2*(logLik(glm_MLE_fit_restr) - logLik(glm_MLE_fit_full))        
    return(as.numeric(LR[1]))    
}

LR_Par_Sim <- replicate(200, ParBootMLE())
hist(LR_Par_Sim)

# our approximate p-value
(length(LR_Par_Sim[LR_Par_Sim > Likelihood_ratio]) + 1)/length(LR_Par_Sim)


# Using the canned gls() in nlme. WHile it is somewhat simpler to write out the function, I often have no idea how to extract the elements I want (so I stick with writing my own support function).
# I actually find the methods and the ability to extract features really painful in nlme, and lme4 (both primarly written by the same author). However there are many examples that use them, so some familiarity is probably reasonable.

require(nlme)
vf2 <- varIdent(form = ~ 1 | genotype)
GLS_model_1 <- gls(SCT ~ tarsus.centered + genotype + temp + genotype:temp,
    weights = vf2, data=dll.data, method="ML")

summary(GLS_model_1)

# I am not sure why it scales one of the groups to variance of one (and despite our relevel, it is treating dll as the base level), but to give you a sense that this is actually fitting the same model, let's do the same scaling.

# So it provides the residual standard error (for the first level, which for gls is the mutant)
GLS_model_1$sigma

# Which should be the same as 
sum(coef(glm_MLE_fit_1)[6:7]) # close enough computationally, although I am not sure why they are not identical.

# gls then provides the values for the variance function, scaled by the sigma above (in the variance function)
coef(glm_MLE_fit_1)[6]/sum(coef(glm_MLE_fit_1)[6:7])    

# Or examine the log likelihoods
logLik(glm_MLE_fit_1)
logLik(GLS_model_1)
