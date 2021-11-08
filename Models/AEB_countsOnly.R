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
  
  # for(i in 1:nind){
  #   zero[i] ~ dMultiStateHMM(y = Y[i,1:nyear], 
  #                            zfirst = z.first[i], 
  #                            first = first[i], 
  #                            ps = S.t[1:nstates, 1:nstates, 1:(nyear-1), 1:nages],
  #                            po = pi[1:nstates, 1:nstates, 1:(nyear-1)],
  #                            age = age[i, 1:nyear],
  #                            n.occasions = nyear,
  #                            nStates = nstates
  #   )
  # }
  
  ###################################################
  # Likelihood for IPM
  ###################################################
  
  # TODO
  # 
  
  for(r in 1:nyear){
    nF[r] ~ dnorm(n[1,r],sd = sigma.Fl) 
    nP[r] ~ dnorm(n.breeders[r],sd = sigma.Pr)
  }
  sigma.Fl ~ dunif(0,1) # AEB changed this - was really really really precise
  sigma.Pr ~ dunif(0,3)
  
  for(i in 1:18){
    # TODO - fix some of these dimension issues
    n[i,1] <- stable[i]      ## changed to matrix with 1 column because that's how the data appeared in nimble
  }
  # TODO
  # is this correct??
  n.breeders[1] <- n[16,1] + n[17,1]
  
  for(r in 1:(nyear-1)){
    emm[r] ~ dunif(0,5) # TODO - very restrictive, is this reasonable?
  }
  
  for(r in 1:(nyear-1)){
    
    ### all age and breeding groups 
    n[1,r+1]  <- n[16,r+1]                         #total fledglings                 
    n[2,r+1]  <- n[1,r]*S[1,r]*0.5                 #1yr NB - females only  
    n[3,r+1]  <- n[2,r]*S[1,r]                     #2yr NB
    n[4,r+1]  <- n[3,r]*S[1,r]                     #3yr NB
    n[5,r+1]  <- n[4,r]*S[1,r]*(1-R.age[4,r])        #4yr NB
    n[6,r+1]  <- n[5,r]*S[1,r]*(1-R.age[5,r])        #5yr NB
    n[7,r+1]  <- n[6,r]*S[1,r]*(1-R.age[6,r])        #6yr NB
    n[8,r+1]  <- n[7,r]*S[1,r]*(1-R.age[7,r])        #7yr NB
    n[9,r+1]  <- n[8,r]*S[1,r]*(1-R.age[8,r])        #8yr NB
    n[10,r+1] <- n[9,r]*S[1,r]*(1-R.age[9,r])        #9yr NB
    n[11,r+1] <- n[10,r]*S[1,r]*(1-R.age[10,r])      #10yr NB
    n[12,r+1] <- n[11,r]*S[1,r]*(1-R.age[11,r])      #11yr NB
    n[13,r+1] <- n[12,r]*S[1,r]*(1-R.age[12,r])      #12yr NB
    n[14,r+1] <- n[13,r]*S[1,r]*(1-R.age[13,r])      #13yr NB
    n[15,r+1] <- n[14,r]*S[1,r]*(1-R.age[14,r])      #14yr NB
    
    n.breeders[r+1] <- n[4,r]*S[1,r]*R.age[4,r] +    #Breeders 
      n[5,r]*S[1,r]*R.age[5,r] +
      n[6,r]*S[1,r]*R.age[6,r] +
      n[7,r]*S[1,r]*R.age[7,r] +
      n[8,r]*S[1,r]*R.age[8,r] +
      n[9,r]*S[1,r]*R.age[9,r] +
      n[10,r]*S[1,r]*R.age[10,r] +
      n[11,r]*S[1,r]*R.age[11,r] +
      n[12,r]*S[1,r]*R.age[12,r] +
      n[13,r]*S[1,r]*R.age[13,r] +
      n[15,r]*S[1,r]*R.age[15,r] +
      n[16,r]*S[2,r]*B[1,r] +
      n[17,r]*S[2,r]*B[2,r] +
      n[18,r]*S[2,r] + 
      emm[r]
    
    n[16,r+1] <- n.breeders[r+1]*F[r] 
    n[17,r+1] <- n.breeders[r+1]*(1-F[r])
    
    n[18,r+1] <- n[16,r]*S[2,r]*(1-B[1,r]) +                #Skipped breeders
      n[17,r]*S[2,r]*(1-B[2,r]) 
  }
  
  
  ###################################################
  # Derived parameter for population growth rate
  ###################################################
  for(r in 1:(nyear-1)){
    lambda.t[r]<- n.breeders[r+1]/n.breeders[r]   # CHANGED TO GROWTH RATE rather than difference in number of breeding pairs between years
    loglam.t[r]<-log(lambda.t[r]) ## for calculating geometric mean of overall population growth rate
  }
  
  #### DERIVED PARAMETER: OVERALL POPULATION GROWTH RATE  #########
  mean.lambda<-exp((1/(nyear-1))*sum(loglam.t[1:(nyear-1)]))   # Geometric mean 
  
})  # End nimble Model code
