
library(ggplot2)
library(nimble)
library(MCMCvis)

#* Marked individuals
m.counts <- read.csv(file.path("data","encounter-history_counts.csv"),check.names=F)
# Marking occasions
mark.occ <-  read.csv(file.path("data","mark_occasions.csv"),check.names=F)

#* Population counts
pop.counts <- read.table(file.path("data","PopulationCounts.txt"),h=T)

pop.counts$data <- factor(pop.counts$data, levels=pop.counts$data)

n.visits <- as.numeric(nrow(pop.counts))

# Organize data
detected <- marked <- n.dets <- double <- numeric()

for(j in 1:n.visits){
  # detection
  ids <- which(!is.na(mark.occ[,j+1]))
  detected[j] <- sum(!is.na(m.counts[ids,j+1]))
  marked[j] <- length(ids)
  # double counts
  m.tot <- sum(!is.na(m.counts[,j+1]))
  n.dets[j] <- sum(m.counts[,j+1], na.rm=T)
  double[j] <- n.dets[j] - m.tot
}
marked[8] <- detected[8] <- n.dets[8] <- double[8] <- NA

counts <- pop.counts$nesting + pop.counts$walking


# Model 1: only detection -------------------------------------------------

#* Define model
mod1 <-  nimbleCode({
  #* Priors
  p ~ dunif(0,1) # detection prob.
  
  for(t in 1:n.visits){
    #* Marked inds.
    # Detection process
    detected[t] ~ dbin(p, marked[t])
    
    
    #* Counts
    N[t] <- counts[t] / p
  } #t
  Ntot <- sum(N[1:n.visits])
})

# Bundle data for nimble
dat1 <- list(
  detected=detected,
  marked=marked,
  counts = counts,
  n.visits = n.visits 
)

# Initial values
inits1 <-  function() list(
  p=runif(1)
)

# Parameters monitored
params1 <- c("p","Ntot","N")

# MCMC settings
ni <- 20000; nt <- 1; nb <- 10000; nc <- 3

# Run Nimble!
out1 <- nimbleMCMC(
  code=mod1,
  constants=dat1,
  inits=inits1,
  monitors=params1,
  niter = ni,
  nburnin = nb,
  nchains = nc
)

# Summarize results
resu1 <- MCMCsummary(out1)
MCMCtrace(out1,pdf=F,Rhat=T,n.eff=T)


# Model 2: detection and doubles ------------------------------------------
#* Define model
mod2 <-  nimbleCode({
  #* Priors
  p ~ dunif(0,1) # detection prob.
  omega ~ dunif(0,1) # proportion of double counts
  
  for(t in 1:n.visits){
    #* Marked inds.
    # Detection process
    detected[t] ~ dbin(p, marked[t])
    # Double counts
    doubles[t] ~ dbin(omega, n.dets[t])
    
    #* Counts
    N[t] <- (counts[t]*(1-omega)) / p
  } #t
  Ntot <- sum(N[1:n.visits])
})

# Bundle data for nimble
dat2 <- list(
  detected=detected,
  marked=marked,
  doubles=double,
  n.dets=n.dets,
  counts = counts,
  n.visits = n.visits 
)

# Initial values
inits2 <-  function() list(
  p=runif(1),
  omega=runif(1)
)

# Parameters monitored
params2 <- c("p","omega","Ntot","N")

# MCMC settings
ni <- 20000; nt <- 1; nb <- 10000; nc <- 3

# Run Nimble!
out2 <- nimbleMCMC(
  code=mod2,
  constants=dat2,
  inits=inits2,
  monitors=params2,
  niter = ni,
  nburnin = nb,
  nchains = nc
)

# Summarize results
resu2 <- MCMCsummary(out2)
MCMCtrace(out2,pdf=F,Rhat=T,n.eff=T)
