#The code below is designed to estimate species-specific occupancy and detection,
#constant across sampling locations using the community occupancy model.
#The data are found in the file "occ data.csv".

#Read in the occurence data
data1 <- read.table("occ data.csv", header=TRUE,sep=",",na.strings=TRUE)
data1$Occ <- rep(1, dim(data1)[1])
#See the first ten lines of data
data1[1:10,]
#How many citings for each species
total.count = tapply(data1$Occ, data1$Species, sum)

#Find the number of unique species
uspecies = as.character(unique(data1$Species))
#n is the number of observed species
n=length(uspecies)

#Find the number of unique sampling locations
upoints = as.character(unique(data1$Point))
#J is the number of sampled points
J=length(upoints)

#Reshape the data using the R package "reshape"
library(reshape)

#The detection/non-detection data is reshaped into a three dimensional
#array X where the first dimension, j, is the point; the second
#dimension, k, is the rep; and the last dimension, i, is the species.
junk.melt=melt(data1,id.var=c("Species", "Point", "Rep"), measure.var="Occ")
X=cast(junk.melt, Point ~ Rep ~ Species)

#Add in the missing lines with NAs
for (i in 1: dim(X)[3]) {
   b = which(X[,,i] > 0)
   X[,,i][b] = 1
   X[,,i][-b] = 0
   X[,,i][1:36,4] = NA;  X[,,i][38:56,4] = NA;
   X[,,i][59:61,4] = NA;  X[,,i][66:70,4] = NA;
}

#Create all zero encounter histories to add to the detection array X
#as part of the data augmentation to account for additional
#species (beyond the n observed species).

#X.zero is a matrix of zeroes, including the NAs for when a point has not been sampled
  X.zero = matrix(0, nrow=70, ncol=4)
  X.zero[1:36,4] = NA;  X.zero[38:56,4] = NA;
  X.zero[59:61,4] = NA;  X.zero[66:70,4] = NA;
  
#K is a vector of length J indicating the number of reps at each point j
KK <- X.zero
a=which(KK==0); KK[a] <- 1
K=apply(KK,1,sum, na.rm=TRUE)
K=as.vector(K)

################

#Write the model code to a text file (used to run WinBUGS)
cat("
	model{

#Define prior distributions for community-level model parameters

u.mean ~ dunif(0,1)
mu.u <- log(u.mean) - log(1-u.mean)

v.mean ~ dunif(0,1)
mu.v <- log(v.mean) - log(1-v.mean)

tau.u ~ dgamma(0.1,0.1)
tau.v ~ dgamma(0.1,0.1)

for (i in 1:n) {

#Create priors for species i from the community level prior distributions
    u[i] ~ dnorm(mu.u, tau.u)
    v[i] ~ dnorm(mu.v, tau.v)

#Create a loop to estimate the Z matrix (true occurrence for species i
#at point j.
   for (j in 1:J) {
       logit(psi[j,i]) <- u[i]
        Z[j,i] ~ dbern(psi[j,i])

#Create a loop to estimate detection for species i at point k during #sampling period k.
     for (k in 1:K[j]) {
    	logit(p[j,k,i]) <-  v[i]
       mu.p[j,k,i] <- p[j,k,i]*Z[j,i]
       X[j,k,i] ~ dbern(mu.p[j,k,i])
}   	}		}


#Finish writing the text file into a document called basicmodel.txt
}
",file="basicmodel.txt")


#Create the necessary arguments to run the bugs() command
#Load all the data
sp.data = list(n=n, J=J, K=K, X=X)

#Specify the parameters to be monitored
sp.params = c("u", "v", "mu.u", "mu.v", "tau.u", "tau.v")

#Create a matrix of latent Z values to intialize the model
Z.guess = matrix(NA, nrow=70, ncol=58)
for (j in 1:J) {
  Z.guess[j,] = apply(X[j,,],2,max, na.rm=T)
}

#Specify the initial values
sp.inits = function() {
    list(u=rnorm(n), v=rnorm(n), Z = Z.guess)
           }

#Load the R2jags library
library(R2jags)

#Run the model and call the results "fit"
fit = jags(sp.data, sp.inits, sp.params, "basicmodel.txt", 
           n.chains=3, n.iter=3000, n.burnin=1000, n.thin=4)

#If model is not converged, run it longer.

#See a summary of the parameter estimates
jagsout <- as.mcmc.list(fit$BUGSoutput)
plot(jagsout)
summary(jagsout)

#See baseline estimates of species-specific occupancy and detection in one of
#the habitat types (CATO)
species.occ = fit$BUGSoutput$sims.list$u
species.det = fit$BUGSoutput$sims.list$v

#Show occupancy and detection estimates for all observed species 1:n
psi = plogis(species.occ)
p   = plogis(species.det)

occ.matrix <- cbind(apply(psi,2,mean),apply(psi,2,sd))
colnames(occ.matrix) = c("mean occupancy", "sd occupancy")
rownames(occ.matrix) = uspecies
det.matrix <- cbind(apply(p,2,mean),apply(p,2,sd))
colnames(det.matrix) = c("mean detection", "sd detection")

round(occ.matrix, digits=2)
