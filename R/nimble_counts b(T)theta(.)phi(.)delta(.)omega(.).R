

nimbleCode({
  # Priors ------------------------------------------------------------------

  # Entry probs.
  b[1:J] ~ ddirch(b.pri[1:J])

  # Total pop. size
  Ntot2 ~ dunif(0,20)
  Ntot <- round(Ntot2*10000)
  
  # Population counts -------------------------------------------------------
  # Stick-breaking binomials to represent multinomial entries
  # First occasion
  B[1] ~ dbin(b[1],Ntot)
  
  for(t in 2:(J-1)){
    # remaining Ntot to entry
    r.Ntot[t] <- Ntot - sum(B[1:(t-1)])
    
    # conditional entry prob.
    b.cond[t] <- b[t] / (1 - sum(b[1:(t-1)])) 
    
    B[t] ~ dbin(b.cond[t], r.Ntot[t])
    
  }
  # Last occasion
  B[J] <- Ntot - sum(B[1:(J-1)])    
  
  # Latent population dynamics
  N[1] <- B[1]
  for(t in 2:J){
    N[t] <- N.walk[t-1] + B[t]
  } #t 
  
  for(t in 1:J){
    # Nesting and Walking latent variables
    N.nest[t] ~ dbin(theta,N[t])
    N.walk[t] <- N[t] - N.nest[t]
    
    # Observation process for counts
    # Availability
    C.nest[t] ~ dbin(phi,N.nest[t])
    C.truewalk[t] ~ dbin(phi, N.walk[t])
    
    # Count errors
    C.double[t] ~ dbin(omega, C.walk[t])
    
    # Total counts
    C.tot[t] = C.nest[t] + C.truewalk[t] + C.double[t]
    
  } # t
  
}) # model

