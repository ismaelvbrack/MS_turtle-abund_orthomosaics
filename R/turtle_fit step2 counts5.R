
library(nimble)
library(MCMCvis)
library(foreach)
library(doParallel)

# Import data -------------------------------------------------------------
# Samples from the 1st step
samps1 <- read.table(file.path("outputs","MR5_postMCMC.txt"),h=T)

#* Population counts
pop.counts <- read.table(file.path("data","PopulationCounts.txt"),h=T)

pop.counts$data <- factor(pop.counts$data, levels=pop.counts$data)

# Number of occasions
J <- nrow(pop.counts)

#Initial values
# Counts
Cd.in <- rbinom(J,pop.counts$walking, 0.2)
Ctw.in <- pop.counts$walking - Cd.in

Nn.in <- pop.counts$nesting+ rev(round(seq(5000,20000,length.out=12)))
Nw.in <- Ctw.in+ rev(round(seq(5000,20000,length.out=12)))
N.in  <- Nn.in + Nw.in

B.in <- c(N.in[1],N.in[2:J] - Nw.in[1:(J-1)])


inits2 <- function() list(
  
  C.truewalk=Ctw.in,
  C.double=Cd.in,
  N.nest=Nn.in,
  
  #alpha0=rnorm(1),
  #alpha1=rnorm(1),
  b=rdirch(1,rep(1,J)),
  
  B=c(B.in[1:(J-1)], NA),
  Ntot2=sum(B.in)/10000
  
)

modelC5 <- source(file.path("R",
                            "nimble_counts b(T)theta(.)phi(.)delta(.)omega(.).R"))

# Mean (using only mean estimates from step1) --------------------------------------------------------------------
# Parameters monitored
params <- c("Ntot","b","N","N.nest","N.walk","C.truewalk","C.double","B")

mean.est <- apply(samps1,2,mean)

# Data
dat2 <- list(
  J=J, # number of occasions
  
  phi=mean.est["phi.1."],
  theta=mean.est["theta"],
  omega=mean.est["omega"],
  
  b.pri=rep(1,J), # constant prior for the Dirichlet
  # Count data
  C.tot=rowSums(pop.counts[,c("walking","nesting")]),
  C.nest=pop.counts$nesting, C.walk=pop.counts$walking
)

# MCMC settings
ni <- 40000; nt <- 10; nb <- 10000; nc <- 3

out5.2 <- nimbleMCMC(
  code=modelC5,
  constants=dat2,
  inits=inits2,
  monitors=params,
  niter = ni,
  nburnin = nb,
  nchains = nc
)

resu <- MCMCsummary(out5.2)

MCMCtrace(out5.2,ind=T,pdf=F,
          Rhat=T,n.eff=T,
          params=c("Ntot","B","b"))

# Multiple samples --------------------------------------------------------

# Parameters monitored
params <- c("Ntot","b","N","N.nest","N.walk","B") #,"alpha0","alpha1"

# MCMC settings
ni <- 40000; nt <- 100; nb <- 20000; nc <- 3


n.samps=881 # number of random samples from the 1st step

samps1 <- samps1[sample(1:nrow(samps1), n.samps),]

##** Running in parallel
registerDoParallel(makeCluster(5))

#out2 <- 
foreach(s=1:n.samps, .errorhandling = 'remove') %dopar% {
  
  require(nimble)
  # Data
  dat2 <- list(
    J=J, # number of occasions
    time=as.vector(scale(1:J)),
    
    phi=samps1[s,"phi.1."],
    theta=samps1[s,"theta"],
    omega=samps1[s,"omega"],
    b.pri=rep(1,J),
    # Count data
    C.tot=rowSums(pop.counts[,c("walking","nesting")]),
    C.nest=pop.counts$nesting, C.walk=pop.counts$walking
  )
  
  # Fit
  out <- nimbleMCMC(
    code=modelC5,
    constants=dat2,
    inits=inits2,
    monitors=params,
    niter = ni,
    nburnin = nb,
    thin=nt,
    nchains = nc
  )
  samps <- rbind(out$chain1,
                 out$chain2,
                 out$chain3)
  
  # Save samples
  write.table(samps,file.path("outputs","Counts5_phi1_postMCMC",paste0("s",s,".txt")),row.names=F)
  # Note that we create a folder to receive all the outputs, from each analysis fitted for each posterior sample of step1
  
  # Clean cache
  gc()
  nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)
  nimble:::clearCompiled(model) 
} # parallel s

