
# How to deal with temporal or spatial correlation (time series data)
#EXAMPLE FROM ZUUR ET AL. 2009

# library(AED) # CURRENTLY NOT WORKING IN  R 2.10.X.
library(nlme)

Hawaii <- read.table("/Users/ian/ZOL851_FS2010/Lectures/Module_15 Mixed Effects_AM/Hawaii.txt", h=T)
#data(Hawaii)
#Abundace of three species of birds collected from 1956-2003.
Hawaii <- na.omit(Hawaii)

Hawaii$Birds <- sqrt(Hawaii$Moorhen.Kauai) # for variance stabilization. Of course we could just use a different variance structure, but let's keep this simple

plot(Hawaii$Year, Hawaii$Birds, xlab="Year", ylab="abundance",type="b")

# What if we wanted to understand the relationship between annual rainfall and abundance, how do we take into account the increasing abundance over time?

M.0 <- lm(Birds ~ Year + Rainfall,  data=Hawaii)
summary(M.0)
# Rainfall seems to have no effect but is that only because year so large an effect? We have not really accounted for the fact that the abundance in one year could be correlated to the abundance in the next year.

par(mfrow=c(2,1))
plot(resid(M.0)~ Hawaii$Year, type="b")

# We can see this
acf(resid(M.0))

# Clearly we have violated the assumption of independence. More importantly we are not accounting for a biologically important source of information!!!!




# Compound symmetry structure in R using gls()
M.1 <- gls(Birds ~ Rainfall + Year, data= Hawaii, correlation = corCompSymm( form= ~Year))

M.2 <- gls(Birds ~ Rainfall + Year, data= Hawaii, correlation = corAR1( form = ~ Year))

summary(M.1)
summary(M.2)
AIC(M.1)
AIC(M.2)
require(bbmle)

AICctab(M.1,M.2, nobs=nrow(Hawaii))

acf(resid(M.2))