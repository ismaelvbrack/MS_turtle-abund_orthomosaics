

library(nimble)
library(MCMCvis)

# Import and Arrange data -------------------------------------------------

#* Population counts
pop.counts <- read.table(file.path("data","PopulationCounts.txt"),h=T)

pop.counts$data <- factor(pop.counts$data, levels=pop.counts$data)

# Number of occasions
J <- nrow(pop.counts)

#* Encounter history matrix with 3 states
states <- read.csv(file.path("data","encounter-history_states.csv"),check.names=F)
#* Counts of marked individuals
m.counts <- read.csv(file.path("data","encounter-history_counts.csv"),check.names=F)

m.walk <- m.double <- numeric()
for(j in 2:(J+1)){
  countj <- m.counts[which(states[,j]==1),j]
  m.walk[j-1] <- sum(countj)
  m.double[j-1] <- sum(countj) - length(countj)
}
#m.walk[8] <- m.double[8] <- NA

cbind(m.walk,m.double)

#* Identified and Unidentified marked individuals
m.ids <- read.table(file.path("data","mark_identification.txt"),h=T)

# Marking occasions
mark.occ <-  read.csv(file.path("data","mark_occasions.csv"),check.names=F)
fo <- apply(mark.occ[,2:(J+1)], 1, function(x) which(x==1))

# Replace with NA before the marking occasion
for(i in 1:nrow(states)){
  if(fo[i]==1){next}
  else{
    states[i,2:fo[i]] <- NA
  }
}

# Exclude occasion 8 (no mark resight data)
states <- states[,-9]
J <- J-1
fo1 <- fo
fo[which(fo %in% (9:12))] <- fo[which(fo %in% (9:12))]-1


# NIMBLE model ------------------------------------------------------------

modelMR5 <- source(file.path("R",
                             "nimble_MR5 theta(.)phi(occ)delta(.)omega(.).R"))

# Organize data for NIMBLE ------------------------------------------------

dat <- list(
  y=as.matrix(states[,2:(J+1)]), 
  m.unids=m.ids[-8,"unidentified"],
  m.detect=(m.ids[-8,"identified"] + m.ids[-8,"unidentified"]),
  m.walk=m.walk[-8] , m.double=m.double[-8],
  J=J,M=length(fo),
  fo=fo
)

# Initial values ----------------------------------------------------------

# Initial values for z can be very tricky... 
z.in <- matrix(NA,nrow=nrow(states),ncol=J)

for(i in 1:nrow(z.in)){
  foi <- fo[i]
  noi <- length(fo[i]:J)
  
  # If all detections are 3, put a 2 in the first
  if(sum(dat$y[i,]==3,na.rm=T)==noi){
    z.in[i,foi:J] <- c(2,rep(3,noi-1))
  }
  # If there is a 2 detection
  if(any(dat$y[i,]==2,na.rm=T)){
    t2 <- which(dat$y[i,]==2) # position of 2
    if(t2>foi){z.in[i,foi:(t2-1)] <- 1} # fill with 1 until the 2
    z.in[i,t2] <- 2 # 2
    if(t2<J){z.in[i,(t2+1):J] <- 3} # fill with 3 after the 2
  }
  # If there is a 1 but and a 2
  if(any(dat$y[i,]==1,na.rm=T) & !any(dat$y[i,]==2,na.rm=T)){
    t1 <- max(which(dat$y[i,]==1)) # position of the last 1
    z.in[i,foi:t1] <- 1
    if(t1<J){z.in[i,t1+1] <- 2} # fill with 2 if there is a 3 after the last 1
    if(t1<(J-1)){z.in[i,(t1+2):J] <- 3} # fill with 3 if we filled with a 2 
  }
  
}
inits <- function() list(
  phi=runif(2),
  theta=runif(1),
  pepa=runif(1,.9,1),
  delta=runif(1),
  omega=runif(1),
  
  z=z.in
)


# Run NIMBLE --------------------------------------------------------------

# Parameters monitored
params <- c("phi","theta","pepa","omega","delta")

# MCMC settings
ni <- 120000; nt <- 1; nb <- 20000; nc <- 3

##* RUN!
out5 <- nimbleMCMC(
  code=modelMR5,
  constants=dat,
  inits=inits,
  monitors=params,
  niter = ni,
  nburnin = nb,
  nchains = nc,
  summary=T
)

beepr::beep(8)

MCMCsummary(out5$samples)

MCMCtrace(out5$samples,ind=T,pdf=F,
          Rhat=T,n.eff=T,
          params=c("phi","theta","pepa","delta","omega"))

samps5 <- rbind(out5$samples$chain1,
                out5$samples$chain2,
                out5$samples$chain3)

effectiveSize(samps5)

# Export samples
write.table(samps5,file.path("outputs","MR5_postMCMC.txt"),row.names=F)

