


nimbleCode({
  # Priors ------------------------------------------------------------------
  
  # availability prob.
  phi[1] ~ dunif(0,1) # first occasion
  phi[2] ~ dunif(0,1) # following occasions
  
  omega ~ dunif(0,1) # prop. of repeated counts
  delta ~ dunif(0,1) # prop. of unidentified marked individuals
  
  # Nesting prob.
  theta ~ dunif(0,1)
  
  pepa ~ dbeta(20,1)
  
  # z transition matrix ------------------------------------------------------
  # z[t-1]=1
  psi[1,1] <- 1-theta
  psi[1,2] <- theta
  psi[1,3] <- 0
  # z[t-1]=2
  psi[2,1] <- 0
  psi[2,2] <- 1-pepa
  psi[2,3] <- pepa
  # z[t-1]=3
  psi[3,1] <- 0
  psi[3,2] <- 0
  psi[3,3] <- 1
  
  # y observation matrix ------------------------------------------------------
  for(i in 1:2){
    
    # z=1
    p[1,1,i] <- phi[i]*delta
    p[1,2,i] <- 0
    p[1,3,i] <- (1-phi[i]) + phi[i]*(1-delta)
    # z=2
    p[2,1,i] <- 0
    p[2,2,i] <- phi[i]*delta
    p[2,3,i] <- (1-phi[i]) + phi[i]*(1-delta)
    # z=3
    p[3,1,i] <- 0
    p[3,2,i] <- 0
    p[3,3,i] <- 1
  }
  
  # Mark-resight likelihood -------------------------------------------------
  
  # Multi-state CJS for theta and phi
  for(i in 1:M){
    
    #* First occasion after marking
    # latent state
    z[i,fo[i]] ~ dcat(psi[1,1:3])
    # observation
    y[i,fo[i]] ~ dcat(p [z[i,fo[i]] ,1:3, 1])
    
    #* Following occasions
    for(t in (fo[i]+1):J){
      # latent state transition
      z[i,t] ~ dcat(psi[z[i,t-1], 1:3])
      # observation
      y[i,t] ~ dcat(p[z[i,t], 1:3, 2])
      
    } # t  
  } # i
  
  # 
  for(t in 1:J){
    # Number of repeated detections
    m.double[t] ~ dbin(omega, m.walk[t])
    # Number of unidentified individuals
    m.unids[t] ~ dbin((1-delta), m.detect[t])
  }
}) # model
