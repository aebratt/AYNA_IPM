GetState <- nimbleFunction(
  run = function(s = double(2), 
                 r = double(2), 
                 f = double(1),
                 b = double(2),
                 nages = double(0),
                 nyear = double(0), 
                 nstates = double(0)
  ) 
  {
    returnType(double(4))
    
    state.mat <- nimArray(0, dim = c(nstates, nstates, nyear-1, nages))
    
    #            |--------------------------------------- true state --------------|
    #true state     Juvenile   Pre-breed    Breed-S    Breed-F      Skip      Dead
    #Juvenile      
    #Pre-breed
    #Breed-S
    #Breed-F
    #Skip
    #Dead
    
    for(i in 1:nages){
      for(t in 1:(nyear-1)){
        state.mat[1, 1:nstates, t, i] <- c(0, s[1, t],            0,                  0,                      0,                 1 - s[1, t])
        state.mat[2, 1:nstates, t, i] <- c(0, s[1, t]*(1-r[i,t]), s[1,t]*r[i,t]*f[t], s[1,t]*r[i,t]*(1-f[t]), 0,                 1 - s[1, t])
        state.mat[3, 1:nstates, t, i] <- c(0, 0,                  s[2,t]*b[1,t]*f[t], s[2,t]*b[1,t]*(1-f[t]), s[2,t]*(1-b[1,t]), 1 - s[2, t])
        state.mat[4, 1:nstates, t, i] <- c(0, 0,                  s[2,t]*b[2,t]*f[t], s[2,t]*b[2,t]*(1-f[t]), s[2,t]*(1-b[2,t]), 1 - s[2, t])
        state.mat[5, 1:nstates, t, i] <- c(0, 0,                  s[2,t]*f[t],        s[2,t]*(1-f[t]),        0,                 1 - s[2, t])
        state.mat[6, 1:nstates, t, i] <- c(0, 0,                  0,                  0,                      0,                 1)
      }
    }
    
    return(state.mat)
  }
)

CGetState <- compileNimble(GetState)


GetDetection <- nimbleFunction(
  run = function(p = double(3), 
                 nyear = double(0), 
                 nstates = double(0)
  ) 
  {
    returnType(double(3))
    
    det.mat <- nimArray(0, dim = c(nstates, nstates, nyear-1))
    
    #            |--------------------------------------- Observed event ---------------------------------------|
    #true state     Juv   Loaf-In   Breed-In-S  Breed-In-F  Breed-Out  Unobs;
    #Juvenile      
    #Pre-breed
    #Breed-S
    #Breed-F
    #Skip
    #Dead
    
    for(t in 1:(nyear-1)){
      det.mat[1, 1:nstates, t] <- c(0, 0,        0,        0,        0,        1) # AEB change - individuals never observed. This was a vector of all 0 before
      det.mat[2, 1:nstates, t] <- c(0, p[1,1,t], 0,        0,        0,        (1-p[1,1,t]))
      det.mat[3, 1:nstates, t] <- c(0, 0,        p[2,1,t], 0,        p[2,2,t], (1-p[2,1,t]-p[2,2,t]))
      det.mat[4, 1:nstates, t] <- c(0, p[3,1,t], 0,        p[3,2,t], p[3,3,t], (1-p[3,1,t]-p[3,2,t]-p[3,3,t]))
      det.mat[5, 1:nstates, t] <- c(0, p[4,1,t], 0,        0,        0,        (1-p[4,1,t]))
      det.mat[6, 1:nstates, t] <- c(0, 0,        0,        0,        0,        1)
    }
    
    return(det.mat)
  }
)

CGetDetection <- compileNimble(GetDetection)

##### DHMM ######

# TODO - these are the functions from the murre model - need to be modified a bit
# this is basically a modification of dDHMMo from the nimbleEcology package
# but the transition matrix has an additional dimension - age
# therefore this is VERY similar to the murre model
# except we *know* age - and we do not reduce to age classes

## define custom distribution where we calcualte log-likelihood given ps and po matrices

# TODO
# check indexing!!! pretty sure it is correct but could use another set of eyes
deregisterDistributions("dMultiStateHMM") # just in case changes are made
dMultiStateHMM <- nimbleFunction(
  run = function(x = double(0), 
                 y = double(1), # vector of observations for individual i
                 zfirst = double(0), # first state at first observation, ind i
                 first = double(0), # first observation, ind i
                 ps = double(4), # state transition matrix
                 po = double(3), # observation matrix
                 age = double(1), # vector of age for individual i
                 n.occasions = double(0), # number of **observation occasions**
                 nStates = double(0), # number of possible states
                 log = double(0)) {
    
    returnType(double(0)) # returns log likelihood

    logLike <- 0 # placeholder
    
    # prob of state at k
    tmp <- matrix(0, n.occasions, nStates) # placeholder
    
    # we assume state is known at capture (and everything before first cap is ignored). 
    tmp[first, zfirst] <- 1
    
    # start in state zfirst then transition
    for(k in (first+1):n.occasions) {
      # tmp[k,s] is prob in state s at k, which is 
      # prob that any state at k-1 (tmp[(k-first), 1:nStates]) transtions to s during interval k-1 at age[k-1]
      # and detection at k give age at k and observation at k
      tmp[k,1:nStates] <-          
        (tmp[k-1, 1:nStates] %*% # prob of state at k-1 
         ps[1:nStates, 1:nStates, k-1, age[k-1]]) * # prob of transtiion from states at k-1 to state s duing interval k-1 if age[k-1]
        po[1:nStates, y[k], k-1]    # prob of detection across all states at occasion k if age[k] given what was oberved at k
      # po[,,k-1,] b/c first resight starts at occasion 2. Thus if k=2, we grab dimension 1 (k-1)
    }
    
    # total prob given y 
    prob <- sum(tmp[n.occasions, 1:nStates]) 
    
    # log likelihood
    logLike <- log(prob)
    
    return(logLike)
  }
)

##### RHMM ######
## nimble requires a random value generator for custom distributions
rMultiStateHMM <- nimbleFunction(
  run = function(n = integer(0), 
                 y = double(1), # vector of observations for individual i
                 zfirst = double(0), # first state at first observation, ind i
                 first = double(0), # first observation, ind i
                 ps = double(4), # state transition matrix
                 po = double(3), # observation matrix
                 age = double(1), # vector of age for individual i
                 n.occasions = double(0), # number of **observation occasions**
                 nStates = double(0) # number of possible states
                 ) { 
    
    returnType(double(0))
    return(1)
  }
)

