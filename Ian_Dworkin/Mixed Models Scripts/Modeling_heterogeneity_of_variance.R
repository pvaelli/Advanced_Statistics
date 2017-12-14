# Modeling heterogeneity in residual variance.
# updated Nov 2012
# This is a brief tutorial demonstrating how to model heterogeneity in variance, first by writing out the support function, and then using the pre-built gls() in the nlme library

library(nlme)
require(bbmle)
# library(AED) # This has the example dataset. THIS IS CURRENTLY NOT WORKING FOR  R 2.10.X
Squid <- read.table("/Users/ian/ZOL851_FS2010/Lectures/Module_15 Mixed Effects_AM/Squid.txt", h=T)
# http://highstat.com/Book2/AED_1.0.zip
# Install from local binary package

### Let us motivate this with an example
# This data set was to investigate factors that influenced squid sexual maturity
#data(Squid) # data is in the AED package
str(Squid) # DML is dorsal mantle length
Squid$fMONTH <- factor(Squid$MONTH)


M.1 <- lm (Testisweight ~ DML*fMONTH, data=Squid) # full model, we are not fitting it
op <- par(mfrow=c(2,2), mar=c(4,4,2,2))
plot(M.1, which = c(1))
plot(Squid$fMONTH, resid(M.1), xlab="Month", ylab="Residuals")
plot(Squid$DML, resid(M.1), xlab="Dorsal Mantle Length", ylab="Residuals")
par(op)

summary(M.1)
# Clearly variance increases with DML, and the residuals show a great deal of heterogeneity.

###  Explicit MLE approach (writing out our own Support function)

# Just as we have done for the "fixed" or deterministic part of a model using MLE, we can model the heterogeneity in the variance.

# In this case, since it looks like the heterogeneity in residual variation increases with increasing  DML, we can start there.
design.M.1 <- model.matrix(~ DML, data=Squid)  # Design matrix for this

NLL.M.1.HET.VAR <- function(b0, b1, sig){
	het.sig <- sig*sqrt(design.M.1[,2])  # Here is what's new. Modeling the variance as the sqrt of DML. See Chapter 5.2 in Pinheiro and Bates (page 209) for why we use this "link" function for the variance.
	det <- b0 + b1*design.M.1[,2]  # We write out our deterministic part like before.
	-sum(dnorm(Squid$Testisweight, mean=det, sd=het.sig, log=T))
	}

M.1.Het <- mle2(NLL.M.1.HET.VAR , start=list(b0=0, b1=0, sig=0.01)) 
summary(M.1.Het)
logLik(M.1.Het)
# Remember that the variance is no longer a constant, so that "sig" represents the slope of the fit


# There is  a pre-built function in the nlme library, gls() : GENERALIZED LEAST SQUARES
# Compare this to gls
vf1Fixed <- varFixed(~DML)
M.1.gls <- gls(Testisweight ~ DML, weights = vf1Fixed, data=Squid, method="ML") # The default method is REML, which is generally better, but for comparison with our ML approach above, we need to use method="ML"

summary(M.1.gls)
logLik(M.1.gls)
M.1.gls$sigma # The parameter for the variance.
plot(M.1.gls)

M.1.lm <- lm(Testisweight ~ DML, data=Squid)
summary(M.1.lm)

# Now that you get the general idea, let's stick with the gls() so we can keep the script short!

vf1Fixed <- varFixed(~DML)
M.2.gls <- gls(Testisweight ~ DML*fMONTH, weights = vf1Fixed, data=Squid) # using default "REML" method/
summary(M.2.gls)
M.2.gls$sigma