# TODO
# edit model

# TODO
# what to do about the unknown age thing

YNAL.IPM <- nimbleCode({ 
  
  # OBSERVATION MATRIX ######
  
  # observed 
  #            |--------------------------------------- Observed event ---------------------------------------|
  #true state     Juv   Loaf-In   Breed-In-S  Breed-In-F  Breed-Out  Unobs;
  #Juvenile      
  #Pre-breed
  #Breed-S
  #Breed-F
  #Skip
  #Dead
  
  # TODO
  # it might be faster to make the array 1 year at a time
  # whichever involves updating fewer nodes at once
  # ie block wrt to year?
  # TODO
  # individuals are never observed - this is a change from a previous version so check correct
  pi[1:nstates, 1:nstates, 1:(nyear-1)] <- GetDetection(p[1:4, 1:3, 1:(nyear-1)], 
                                                        nyear, 
                                                        nstates)
  
  #OBSERVATION MODELS AND PRIORS 
  for(t in 1:(nyear-1)){
    # just setting up which transitions are possible, ever
    p[1,2,t] <- 0  
    p[1,3,t] <- 0 
    p[2,3,t] <- 0
    p[4,2,t] <- 0
    p[4,3,t] <- 0
    
    # pre breeder to loafer
    p[1,1,t] <- 1/(1+exp(-p.link[1,t])) 
    p.link[1,t] ~ dnorm(int.p[1],sd = sigma.p[1])
    
    # skip breeder to loafer
    p[4,1,t] <- 1/(1+exp(-p.link[2,t])) 
    p.link[2,t] ~ dnorm(int.p[2],sd = sigma.p[2])
  }
  
  # for years when effort was zero out of colony
  for(t in 1:nyear.in){
    
    # successful breeder to successful breeder
    p[2,1,yrs.in[t]] <- 1/(1+exp(-p.link[3,t])) 
    p.link[3,t] ~ dnorm(int.p[3],sd = sigma.p[3])
    
    # successful breeder to breeder out - impossible in years when no out obs were made
    p[2,2,yrs.in[t]] <- 0
    
    # failed breeder to loafer
    p[3,1,yrs.in[t]] <- exp(p.link[4,t])/(1+exp(p.link[4,t])+exp(p.link[5,t]))
    p.link[4,t] ~ dnorm(int.p[4],sd = sigma.p[4])
    
    # failed breeder to failed breeder
    p[3,2,yrs.in[t]] <- exp(p.link[5,t])/(1+exp(p.link[4,t])+exp(p.link[5,t]))
    p.link[5,t] ~ dnorm(int.p[5],sd = sigma.p[5])
    
    # failed breeder to breeder out - impossible in years when no out obs were made
    p[3,3,yrs.in[t]] <- 0
    
  }
  
  # for years when effort was nonzero out of colony
  for(t in 1:nyear.out){
    
    # successful breeder to successful breeder
    p[2,1,yrs.out[t]] <- exp(p.link[6,t])/(1+exp(p.link[6,t])+exp(p.link[7,t]))
    p.link[6,t] ~ dnorm(int.p[6],sd = sigma.p[6])
    
    # successful breeder to breeder out
    p[2,2,yrs.out[t]] <- exp(p.link[7,t])/(1+exp(p.link[6,t])+exp(p.link[7,t]))
    p.link[7,t] ~ dnorm(int.p[7],sd = sigma.p[7])
    
    # failed breeder to loafer
    p[3,1,yrs.out[t]] <- exp(p.link[8,t])/(1+exp(p.link[8,t])+exp(p.link[9,t])+exp(p.link[10,t]))
    p.link[8,t] ~ dnorm(int.p[8],sd = sigma.p[8])
    
    # failed breeder to failed breeder
    p[3,2,yrs.out[t]] <- exp(p.link[9,t])/(1+exp(p.link[8,t])+exp(p.link[9,t])+exp(p.link[10,t]))
    p.link[9,t] ~ dnorm(int.p[9],sd = sigma.p[9])
    
    # failed breeder to breeder out
    p[3,3,yrs.out[t]] <- exp(p.link[10,t])/(1+exp(p.link[8,t])+exp(p.link[9,t])+exp(p.link[10,t]))
    p.link[10,t] ~ dnorm(int.p[10],sd = sigma.p[10])
    
  }
  
  # priors
  for(g in 1:10){
    int.p[g] ~ dunif(-15,15)
    sigma.p[g] ~ dunif(0,10)
  }
  
  # END OBSERVATION MATRIX ###########
  
  # STATE TRANSITION MATRIX #########
  
  #S = survive
  #R = recruit (breed for first time)
  #B = breed (breed again)
  #F = fledge
  
  # TODO
  # it might be faster to make the array 1 year or 1 individual at a time
  # whichever involves updating fewer nodes at once
  # should it actually be indexed by age?? rather than nind
  S.t[1:nstates, 1:nstates, 1:(nyear-1), 1:nages] <- GetState(
    S[1:2, 1:(nyear-1)],
    R.age[1:nages, 1:(nyear-1)],
    F[1:(nyear-1)],
    B[1:2, 1:(nyear-1)],
    nages,
    nyear,
    nstates
  )

  # TODO
  # start checking here
  #PROCESS MODELS AND PRIORS
  for(t in 1:(nyear-1)){
    for(m in 1:2){
      S[m,t] <- 1/(1+exp(-S.rand[m,t]))
      S.rand[m,t] ~ dnorm(int.S[m],sd = sigma.S[m])
    }
  }
  for(m in 1:2){
    int.S[m] ~ dunif(-15,15)
    sigma.S[m] ~ dunif(0,10)
  }
  # TODO
  # does this work
  R.age[1, 1:(nyear-1)] <- 0
  R.age[2, 1:(nyear-1)] <- 0
  R.age[3, 1:(nyear-1)] <- 0
  R.rand[1:3] <- 0
  for(g in 4:14){
    R.age[g, 1:(nyear-1)] <- 1/(1+exp(-R.rand[g]))
    R.rand[g] ~ dnorm(int.R,sd = sigma.R)
  }
  for(g in 15:nages){ # TODO - why make this age assumption. Oh - because this is max age over data period
    # TODO remove hard coding here
    # TODO - why not have everything above this age be 1 age class
    R.age[g, 1:(nyear-1)] <- 0
  }
  int.R ~ dunif(-15,15)
  sigma.R ~ dunif(0,10)

  for(s in 1:2){
    for(t in 1:(nyear-1)){
      B[s,t] <- 1/(1+exp(-B.rand[s,t]))
      B.rand[s,t] ~ dnorm(int.B[s],sd = sigma.B[s])
    }
    int.B[s] ~ dunif(-15,15)
    sigma.B[s]  ~ dunif(0,10)
  }

  for(t in 1:(nyear-1)){
    F[t] <- 1/(1+exp(-F.rand[t]))
    F.rand[t] ~ dnorm(int.F,sd = sigma.F)
  }
  int.F ~ dunif(-15,15)
  sigma.F ~ dunif(0,10)

  # END STATE TRANSITION MATRIX ###########

  # LIKELIHOOD FOR MS MODEL ######

  for(i in 1:nind){
    #z[i,first[i]] <- z.first[i]
    
    zero[i] ~ dMultiStateHMM(y = Y[i,1:nyear], 
                             zfirst = z.first[i], 
                             first = first[i], 
                             ps = S.t[1:nstates, 1:nstates, 1:(nyear-1), 1:nages],
                             po = pi[1:nstates, 1:nstates, 1:(nyear-1)],
                             age = age[i, 1:nyear],
                             n.occasions = nyear,
                             nStates = nstates
                             )
    
    # TODO this needs editing if you want to do it this way
    # for(t in (first[i]+1):nyear){
    #   #state equation
    #   z[i,t] ~ dcat(S.t[z[i,(t-1)],t-1,i,1:nstates]) #### THIS NEEDS TO BE FIXED TO TRANSLATE TO NIMBLE
    #   #observation equation
    #   Y[i,t] ~ dcat(pi[z[i,t],t-1,1:nstates])        #### THIS NEEDS TO BE FIXED TO TRANSLATE TO NIMBLE
    #   
    # }
  }
  
})  # End nimble Model code
