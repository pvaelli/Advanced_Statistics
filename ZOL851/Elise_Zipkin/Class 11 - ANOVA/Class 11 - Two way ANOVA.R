#TWO WAY ANOVA
#We will use the mourning cloak butterfly as an illustration for two way ANOVAs.
#Assume that we measured wing length of butterflies in each of three elevation
#classes in five different populations and that the effects of these factors
#interact. 

# Generate the data
# Choose sample size
n.pop <-5
n.elev <-3
nsample <-12
n <-n.pop*nsample

# Create factor levels
pop <-gl(n=n.pop,k=nsample,length=n)
elev <-gl(n=n.elev,k=nsample/n.elev,length=n)

# Choose effects
baseline <-40 #Intercept
pop.effects = c(-10, -5, 5,10) #Population effects
elev.effects  = c(5,10) #Elev effects
interaction.effects <- c(-2, 3,0,4,4,0,3, -2) #Interaction effects
all.effects <- c(baseline,pop.effects,elev.effects,interaction.effects)
sigma <- 3
eps <-rnorm(n,0,sigma) #Residuals
X <-as.matrix(model.matrix(~pop*elev)) #Create design matrix
X #Have a look at the design matrix

#Use matrix multiplication to assemble all components of the data
wing <-as.numeric(as.matrix(X)%*%as.matrix(all.effects)+eps)

#Plot the generated data
boxplot(wing ~elev*pop, col="grey",xlab="Elevation by Population", ylab=
          "Wing length",main="Simulated data set", las=1,ylim=c(20,70))
abline(h =40)

#Explore the data further by looking at a conditioning plot to examine
#the relationship between population and elevation on the data

library("lattice") # Load the lattice library
xyplot(wing ~elev|pop,ylab="Winglength",xlab="Elevation",main=
         "Population-specific relationship between wing and elevation class")

xyplot(wing ~pop|elev,ylab="Winglength",xlab="Population",main=
         "Elevation-specific relationship between wing and population")

#########
#Compare to the effects we used to generate the data
all.effects
#To the coefficient estimates  of the two-way interaction ANOVA model 
lm(wing ~ pop*elev)

#Are the parameter estimates the same?  Why or why not?

#As a check, let's create 1000 data sets and run the analysis 
#on all of those.

n.iter <-1000 # Desired number of iterations
estimates <-array(dim=c(n.iter,length(all.effects))) #Data structure to hold results

for(i in 1:n.iter){  #Run simulation n.iter times
  eps <- rnorm(n,0,sigma) #Residuals
  y <- as.numeric(as.matrix(X) %*% as.matrix(all.effects)+eps) #Assemble data
  fit.model <- lm(y ~ pop*elev)  #Break down data
  estimates[i,] <- fit.model$coefficients  #Save estimates of coefs.
}

print(apply(estimates,2,mean),dig=2)  #What does the apply function do?
all.effects

#What does this tell us?  Why did we do this?


####################
#Try fitting the model with only main effects first
mainfit = lm(wing ~elev+pop)
mainfit

#Let's now fit the interaction model with a means parameterization
intfit <-lm(wing ~ elev*pop -1 -pop -elev)
intfit

