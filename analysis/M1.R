
########################################################################
#        M1: SURVIVAL  ~ Lifetime provisioning performance             #
#                                                                      #
# Date last modified: 10 August 2022                                   #
########################################################################


library(coda)
library(rjags)
library(shinystan)
library(tidyverse)

rootdir<- "~/Desktop/survival analysis"
datadir<- file.path(rootdir, "data")
date<-Sys.Date()
outdir<-file.path(rootdir, "output", date)
dir.create(outdir)

load(file.path(datadir, "jags_dataPP.RData"))        # A list of data for JAGS survival model.
#                                         nind = number of individuals, 
#                                         nyears = time series length,  
#                                         ageb = age, discretized into age groups 1-3, for each individual and time step
#                                         ap = age at first reproduction for each individual 
#                                         ch = capture histories for each individual over time series. Coded as 1 for sighted and 0 for not sighted

# List of data for embedded regression model on pup weaning mass. 
# Variables with a "v" after them are in the form of vectors for the regression 
# on pup weaning mass because we did not have data for every individual every year 
# (due to skip breeding or failed to collect mass).So,     
#                                         n.obs = number of observations, 
#                                         pwm_v = pup weaning masses, 
#                                         age_v = standardized age  
#                                         age_v2 = standardized age^2
#                                         par_v = parity
#                                         sex_v =  pup sex
#                                         afr_v = age at first reproduction 
#                                         amo_v = Atlantic Multidecadal Oscillation, averaged over preceeding 3 years
#                                         ind_v = individual identifier 
#                                         year_v = year
#                                         age_v2 = standardized age^2
#                                         age_v2 = standardized age^2



known.state.cjs <- function(ch){          # creates a matrix in which only the known states are coded in. So, as we only see 
  state <- ch                             #      individuals when ch = 1, whenever ch = 0 this new "state" matrix is NA.
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1)) 
    n2 <- max(which(ch[i,]==1)) 
    state[i,n1:n2] <- 1 
    state[i,n1] <- NA
  }
  state[state==0] <- NA 
  return(state)
}
jags.data[["z"]]=known.state.cjs(jags.data[["ch"]])            # adding this state matrix to the data list to be fed into JAGS.


get.first <- function(x) min(which(x!=0))                         # extracts the time step each individual was first seen and 
jags.data[["f"]] <- apply(jags.data[["ch"]], 1, get.first)        # adds to data list.  



#Model file 
jags.model.txt <- "model{

## Stage 1: Pup weaning mass mixed effects model 
### Parameters: pi, regression coefficients
#               alpha_mass, individual intercepts
#               epsilon_mass, year effect 

for(j in 1:n.obs){
      pwm_v[j] ~ dnorm(mu[j], tau)
      
      mu[j] <- intercept + pi[1]*age_v[j] + pi[2]*age_v2[j] + pi[3]*par_v[j] + pi[4]*sex_v[j] + pi[5]*afr_v[j] 
                + pi[6]*amo_v[j] + alpha_mass[ind_v[j]] + epsilon_mass[year_v[j]]
}


intercept ~ dnorm(50, 0.01)               # Mean of prior reflects prior knowledge on population (see Bowen et al. 2007, along with many others)

for(k in 1:6){
  pi[k] ~ dnorm(0,0.001)
}

for(i in 1:nind){
alpha_mass[i] ~ dnorm(0, tau.i.mass)
}

for(t in 1:nyears){
epsilon_mass[t] ~ dnorm(0,tau.t.mass)
}


sigma ~ dunif(0,10)
tau <- pow(sigma, -2)

sigma.t.mass ~ dunif(0,10)
tau.t.mass<- pow(sigma.t.mass,-2)
sigma2.t.mass<-pow(sigma.t.mass,2)

sigma.i.mass ~ dunif(0,10)
tau.i.mass <- pow(sigma.i.mass, -2)
sigma2.i.mass<-pow(sigma.i.mass,2)


## Stage 2: Survival analysis
### Parameters: beta, regression coefficients
#               alpha, individual intercepts for survival 
#               epsilon, year effect on survival 
#               p, detection probability 
#               



for(i in 1:nind){
for(t in f[i]:(nyears-1)){
  logit(phi[i,t]) <- beta[ageb[i,t]] + beta[4]*amo[t] + beta[5]*afr[i] + beta[6]*alpha_mass[i] + epsilon[t] + alpha[i]
  p[i,t] <-mean.p
}
}

mean.p ~ dbeta(1,1)


for(i in 1:nind){
alpha[i] ~ dnorm(0, tau.i)
}

for(t in 1:(nyears-1)){
epsilon[t] ~ dnorm(0,tau.t)
}

for(i in 1:6){
beta[i] ~ dnorm(0,0.001)T(-10,10)
}

sigma.t ~ dunif(0,10)
tau.t<- pow(sigma.t,-2)
sigma2.t<-pow(sigma.t,2)

sigma.i ~ dunif(0,10)
tau.i <- pow(sigma.i, -2)
sigma2.i<-pow(sigma.i,2)

for(i in 1:nind){
z[i,f[i]]<-1

for(t in (f[i]+1):nyears){
z[i,t] ~ dbern(mu1[i,t])
mu1[i,t] <- phi[i,t-1]*z[i,t-1]

ch[i,t] ~ dbern(mu2[i,t])
mu2[i,t] <- p[i,t-1]*z[i,t]

log_lik[i,t]<-logdensity.bern(ch[i,t],mu2[i,t])
}
lli[i]<- sum(log_lik[i,((f[i]+1):nyears)])

}
}" 

jags.file <- "JAGS_survival_lifePP.JAG"
# save the JAGS model syntax to a local file
sink(file=file.path(outdir, jags.file)) # open connection
cat(jags.model.txt,fill=TRUE) # send model syntax to the file
sink() # close connection


cjs.init.z <- function(ch,f){                    # initialize the state matrix.
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1)) 
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){ 
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

z.init<-cjs.init.z(jags.data[["ch"]],jags.data[["f"]])

inits<-function(){list(z=z.init, 
                       sigma.i=runif(1,0,5), 
                       sigma.t=runif(1,0,5), 
                       beta=runif(6,-5,5), 
                       mean.p=runif(1,0,1))}

n.chains<-2
n.adapt<-1000
n.burnin<-5000
n.iter<-20000
n.thin<-10


# COMPILE THE JAGS MODEL
m <- jags.model(file=file.path(outdir,jags.file), data=jags.data, n.chains=n.chains,
                n.adapt=n.adapt)

# BURN-IN PHASE
update(m, n.burnin)
# SAMPLE FROM POSTERIORS
post <- coda.samples(m,variable.names=c("intercept",
                                        "pi",
                                        "sigma.t.mass",
                                        "sigma.i.mass",
                                        "beta",
                                        "sigma.i",
                                        "sigma.t",
                                        "mean.p",
                                        "lli"
),
n.iter=n.iter,thin=n.thin)
saveRDS(post, file.path(outdir,"post_survival_lifePP.RData"))




k<-as.shinystan(post)
launch_shinystan(k)

#WAIC approximation function
waic1.jags <- function(post,n.obs,loglike.grep.txt = "lli"){
  logsumexp <- function(x){ max.x <- max(x); max.x - log(length(x)) +log(sum(exp(x-max.x)))} # this is equivalient to log( 1/S * sum_s[exp{ x }]))    
  # get the columns indices that have 'lli' in their name (for loglike-inidividual)
  ll.col <-grep(loglike.grep.txt,colnames(post[[1]]))
  # get the mcmc values for the log values
  if(length(post)>1){ lli.samp <- as.matrix(do.call("rbind",post)[,ll.col]) } else { lli.samp <- as.matrix(post[[1]][,ll.col])}
  # get the complexity penality (WAIC1; from Gelman et al 2014)
  logElike <- apply(lli.samp,2,logsumexp)
  Eloglike <- apply(lli.samp,2,mean)
  p_waic1 <- 2*sum(logElike-Eloglike)
  # get the lppd log pointwaise predictive density
  # use the log-sum-exp trick to handle underflow issues with exp(loglike) ~= 0
  logsumf.i <- apply(lli.samp,2,logsumexp)
  lppd <- sum(logsumf.i)
  waic1 <- -2*(lppd-p_waic1)
  return(waic1)
}
waic1.jags(post, nind)
