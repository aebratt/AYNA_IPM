##########################################################################
#
# ATLANTIC YELLOW-NOSED ALBATROSS INTEGRATED POPULATION MODEL 2000-2018
#
##########################################################################
# based on Kery and Schaub 2012, Chapter 11
# modified by Steffen oppel, December 2018
# implemented in NIMBLE by Abby Bratt, September 2021
# code for survival and trend SSM in separate files
# this IPM is much simpler than previous attempts by Horswill, Converse etc.

# TODO
# add complexity
# lots of posterior predictive samplers
# lots of problems with the rounded things it seems
# and is therefore not converging


#### LOAD LIBRARIES #####
library(tidyverse)
library(data.table)
library(nimble)
library(here)
library(coda)

#### SOURCE FILES ####
# CMR data
AYNA <- read_csv(here("Data", "AYNA_simple_encounter_history_1982_2018.csv"))
names(AYNA)
CH<-as.matrix(AYNA[,3:39], dimnames=F)
# TODO - is this a reasonable assumption???
AYNA$AGE[is.na(AYNA$AGE)]<-1  ## set all NA as 'adult'
apply(CH,2,sum) ### check that there are contacts in every season

# COUNT DATA FOR POPULATION TREND 
AYNA.pop <- read_csv(here("Data", "AYNA_pop_counts_1982_2018.csv"))
AYNA.pop<-subset(AYNA.pop,Year>1999)	## reduce data set to remove NA in 4 years
n.years<-dim(AYNA.pop)[1]		## defines the number of years
n.sites<-dim(AYNA.pop)[2]-1 ## defines the number of study areas
str(AYNA.pop)
names(AYNA.pop)

# BREEDING SUCCESS DATA FOR FECUNDITY 
AYNA.bs <- read_csv(here("Data", "AYNA_breed_success_1982_2017.csv"))
AYNA.bs<-subset(AYNA.bs,Year>1999)	## reduce data set to specified time period
J<-as.integer(AYNA.bs$n_nests*AYNA.bs$BREED_SUCC)
R<-AYNA.bs$n_nests

# DATA MANIPULATION #####

# CREATE MATRIX OF AGE FOR EACH OCCASION AND INDIVIDUAL
## this matrix will relate to the survival parameter estimates chosen in the model
## simple model only has 2 survival parameters:
## 1 - juvenile and immature survival (years 1-5)
## 2 - adult survival (birds >5 years old)

# COMPUTE VECTO WITH OCCASION OF FIRST CAPTURE
get.first <- function(x) min(which(x==1))
f <- apply(CH, 1, get.first)

## REMOVE BIRDS THAT ARE TOO YOUNG TO HAVE HAD A CHANCE TO RETURN

tooyoung<-ifelse(f>(dim(CH)[2]-5), # if first is greater than nyears-5
                 ifelse(AYNA$AGE==0,1,0), # and age was 0 - bird is too young
                 0) # otherwise, bird was not too young
CH<-CH[tooyoung==0,]  ## removes individuals that were ringed as chicks <5 years before end of time series
f <- apply(CH, 1, get.first)
ages<-AYNA$AGE[tooyoung==0]# AEB changed

toolate<-ifelse(f==dim(CH)[2],1,0)
CH<-CH[toolate==0,]  ## removes individuals ringed in last occasion end of time series
f <- apply(CH, 1, get.first)
ages<-ages[toolate==0]

## CREATE BLANK AGE MATRIX
AGEMAT<-matrix(2,nrow=nrow(CH),ncol=ncol(CH))
n.occ<-ncol(CH)

## LOOP OVER EACH BIRD RINGED AND SET PRE-CAPTURE DATA TO NA AND ADJUST AGE
for (l in 1:nrow(AGEMAT)){
  firstocc<-get.first(CH[l,])
  lastjuv<-firstocc+4
  lastjuv<-ifelse(lastjuv>n.occ,n.occ,lastjuv)
  young<-ages[l]
  if(firstocc>1){AGEMAT[l,1:(firstocc-1)]<-NA}  ## sets everything before first contact to NA
  if(young==0){AGEMAT[l,firstocc:lastjuv]<-1}  ## sets all juvenile years to 1
}

### CHECK WHETHER IT LOOKS OK ###
head(AGEMAT)
head(CH)

# RE-ARRANGE DATA TO REMOVE CONTACTS BEFORE 2000 (and individuals with no contacts after 2000)
## remove encounter occasions before 2000
rCH<-CH[,c(19:37)]	      ## reduce data set to exclude years before 2000
AGEMAT<-AGEMAT[,c(19:37)]	## reduce data set to exclude years before 2000
dim(rCH)
exclude<- apply(rCH,1,sum)
rCH<-rCH[exclude>0,]        ## removes individuals that were not observed in the last 19 years - could also set >1 to remove transients
AGEMAT<-AGEMAT[exclude>0,]  ## removes individuals that were not observed in the last 19 years
ages<-ages[exclude>0]
dim(rCH)
dim(AGEMAT)
head(rCH)
head(AGEMAT)

## PREPARE CONSTANTS
n.ind<-dim(rCH)[1]		## defines the number of individuals
n.years<-dim(rCH)[2]  ## defines the number of years
f <- apply(rCH, 1, get.first)

## CREATE MATRIX for INITIAL STATE Z FOR SURVIVAL MODEL
zinit<-rCH
for (l in 1:nrow(zinit)){
  firstocc<-get.first(zinit[l,])
  zinit[l,1:firstocc]<-NA  ## sets everything up to first contact to NA
  zinit[l,(firstocc+1):n.years]<-1  ## alive after first contact
}
dim(zinit)

## CREATE MATRIX for COUNTS
N.init=matrix(NA, nrow=n.years,ncol=n.sites)
N.init[1,]<-as.matrix(AYNA.pop[1,2:12])

#### MODEL CODE ####
code <- nimbleCode({
  #-------------------------------------------------
  # integrated population model for the Gough AYNA population
  # - age structured model with 6 age classes 
  # - adult survival based on CMR ringing data
  # - pre breeding census, female-based assuming equal sex ratio & survival
  # - productivity based on Area 1 nest monitoring data
  # - simplified population process with informed prior for adults skipping breeding and uninformed immatures recruiting
  # -------------------------------------------------
  
  #-------------------------------------------------  
  # 1. PRIORS FOR ALL DATA SETS
  #-------------------------------------------------
  
  # -------------------------------------------------        
  # 1.1. Priors and constraints FOR FECUNDITY
  # -------------------------------------------------
  for (t in 1:T){  
    ann.fec[t] ~ dunif(0.2,0.8)           # Priors on fecundity can range from 0-1 chicks per pair (constrained based on our data)
    imm.rec[t]~dunif(0.05,0.95)                ## RECRUITMENT PROBABILITY COULD SET MORE INFORMATIVE PRIOR HERE
    skip.prob[t]~dunif(0.15,0.45)              ## PRIOR FOR ADULT BREEDER SKIPPING PROBABILITY from Cuthbert paper that reported breeding propensity of 0.66
  } #t
  
  # -------------------------------------------------        
  # 1.2. Priors and constraints FOR POPULATION COUNTS
  # -------------------------------------------------
  for (s in 1:n.sites){			### start loop over every study area
    N.est[1,s] ~ dunif(0,200)   ## draw random value from a uniform distribution between 0 and 200 for initial population size
    mean.lambda[s] ~ dunif(0,10)	#Prior for mean growth rate
    sigma.proc[s] ~ dunif(0,10)	#Prior for SD of state process (annual variation in pop size)
    sigma.obs[s] ~ dunif(0,100)	#Prior for SD of observation process (variation in detectability)
  }
  
  # -------------------------------------------------        
  # 1.3. Priors and constraints FOR SURVIVAL
  # -------------------------------------------------
  
  ### RECAPTURE PROBABILITY
  mean.p ~ dunif(0, 1)                          # Prior for mean recapture
  logit.p <- log(mean.p / (1-mean.p))           # Logit transformation
  
  for (t in 1:T){
    logit(p[t]) <- logit.p  + capt.raneff[t]
    capt.raneff[t] ~ dnorm(0, sd = sigma.capt)
  }
  
  ### SURVIVAL PROBABILITY
  for (i in 1:nind){
    for (t in f[i]:(T-1)){
      logit(phi[i,t]) <- mu[AGEMAT[i,t]] + surv.raneff[t]
    } #t
  } #i
  
  ## AGE-SPECIFIC SURVIVAL 
  for (age in 1:2){
    beta[age] ~ dunif(0, 1)                         # Priors for age-specific survival
    mu[age] <- log(beta[age] / (1-beta[age]))       # Logit transformation
  }
  
  ## RANDOM TIME EFFECT ON SURVIVAL 
  for (t in 1:(T-1)){
    surv.raneff[t] ~ dnorm(0, sd = sigma.surv)
  }
  
  ### PRIORS FOR RANDOM EFFECTS
  sigma.surv ~ dunif(0, 10)                     # Prior for standard deviation of survival
  sigma.capt ~ dunif(0, 10)                     # Prior for standard deviation of capture
  
  #-------------------------------------------------  
  # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
  #-------------------------------------------------
  
  # -------------------------------------------------        
  # 2.1. System process: female based matrix model
  # -------------------------------------------------
  
  for (tt in 2:T){
    
    ## THE PRE-BREEDING YEARS ##
    nestlings[tt] <- max(1, ann.fec[tt] * 0.5 * Ntot.breed[tt])   ### number of locally produced FEMALE chicks
    JUV[tt] ~ dpois(nestlings[tt])                                ### need a discrete number otherwise dbin will fail, dpois must be >0
    N1[tt]  ~ dbin(ann.surv[1,tt-1], JUV[tt-1])                   ### number of 1-year old survivors 
    N2[tt] ~ dbin(ann.surv[1,tt-1], N1[tt-1])                     ### number of 2-year old survivors
    N3[tt] ~ dbin(ann.surv[1,tt-1], N2[tt-1])                     ### number of 3-year old survivors
    N4[tt] ~ dbin(ann.surv[1,tt-1], N3[tt-1])                     ### number of 4-year old survivors
    N5[tt] ~ dbin(ann.surv[1,tt-1], N4[tt-1])                     ### number of 5-year old survivors
    
    ## THE POTENTIAL RECRUITING YEARS ##
    N6[tt] ~ dbin(ann.surv[1,tt-1], round(N5[tt-1]))                                     ### number of 6-year old survivors that are ready for recruitment
    N.notrecruited[tt] ~ dbin(ann.surv[2,tt-1], max(1,non.recruits[tt-1]))      ### number of not-yet-recruited birds surviving from previous year
    non.recruits[tt]<-N6[tt]+N.notrecruited[tt]-ann.recruits[tt]                      ## number of birds that do not recruit is the sum of all available minus the ones that do recruit
    
    ## THE BREEDING YEARS ##
    Ntot.breed[tt] ~ dpois(pop.size[tt])                                           ### the annual number of breeding birds is the estimate from the count SSM
    ann.recruits[tt] ~ dbin(imm.rec[tt],Ntot.breed[tt]-Nold.breed[tt])          ### this total number comprises a bunch of new recruits, which is the number of total breeders that are not old breeders
    Nold.breed[tt]<- N.pot.breed[tt]-N.non.breed[tt]                              ### number of old breeders is survivors from previous year minus those that skip a year of breeding
    N.pot.breed[tt] ~ dbin(ann.surv[2,tt-1], Ntot.breed[tt-1]+N.non.breed[tt-1])   ### number of potential old breeders is the number of survivors from previous year breeders and nonbreeders
    N.non.breed[tt] ~ dbin(skip.prob[tt], N.pot.breed[tt])                             ### number of old nonbreeders (birds that have bred before and skip breeding) 
  } # tt
  
  ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on stable stage distribution from previous model
  JUV[1]<-round(Ntot.breed[1]*0.5*ann.fec[1])
  N1[1]<-round(Ntot.breed[1]*0.17574058)
  N2[1]<-round(Ntot.breed[1]*0.11926872)
  N3[1]<-round(Ntot.breed[1]*0.10201077)
  N4[1]<-round(Ntot.breed[1]*0.08725001)
  N5[1]<-round(Ntot.breed[1]*0.07462511)
  non.recruits[1]<-round(Ntot.breed[1]*0.3147774)
  Ntot.breed[1]<-sum(y.count[1,])
  N.non.breed[1]<- round(Ntot.breed[1]*0.12632740)
  
  # -------------------------------------------------        
  # 2.2. Observation process for population counts: state-space model of annual counts
  # -------------------------------------------------
  
  for (s in 1:n.sites){			### start loop over every study area
    
    ## State process for entire time series
    for (t in 1:(T-1)){
      lambda[t,s] ~ dnorm(mean.lambda[s], sd = sigma.proc[s])								# Distribution for random error of growth rate
      N.est[t+1,s]<-N.est[t,s]*lambda[t,s]										        # Linear predictor (population size based on past pop size and change rate)
    }														# run this loop over nyears
    
    ## Observation process
    for (t in 1:T){
      y.count[t,s] ~ dnorm(N.est[t,s], sd = sigma.obs[s])								# Distribution for random error in observed numbers (counts)
    }														# run this loop over t= nyears
  }		## end site loop
  
  # -------------------------------------------------        
  # 2.3. Likelihood for fecundity: Poisson regression from the number of surveyed broods
  # -------------------------------------------------
  for (t in 1:(T-1)){
    J[t] ~ dpois(rho.fec[t])
    rho.fec[t] <- R[t]*ann.fec[t]
  } #	close loop over every year in which we have fecundity data
  
  # -------------------------------------------------        
  # 2.4. Likelihood for adult and juvenile survival from CMR
  # -------------------------------------------------
  
  # Likelihood 
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):T){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[t] * z[i,t]
    } #t
  } #i
  
  #-------------------------------------------------  
  # 3. DERIVED PARAMETERS FOR OUTPUT REPORTING
  #-------------------------------------------------
  
  ## DERIVED SURVIVAL PROBABILITIES PER YEAR 
  for (t in 1:(T-1)){
    for (age in 1:2){
      logit(ann.surv[age,t]) <- mu[age] + surv.raneff[t]
    }
  }
  
  ## DERIVED POPULATION SIZE PER YEAR 
  for (t in 1:T){
    pop.size[t]<-sum(N.est[t,1:n.sites])
  }
  
  ## DERIVED OVERALL POPULATION GROWTH RATE 
  pop.growth.rate <- mean(lambda[1:(T-1),1:n.sites])  				# Arithmetic mean for whole time series
  
})

#### DATA ####
dat <- list(y = rCH,
            y.count=as.matrix(AYNA.pop[,2:12]),
            J=J,
            R=R) 

#### CONSTANTS ####
const <- list(f = f,
              T = n.years,
              nind = n.ind,
              AGEMAT=AGEMAT,
              n.sites=n.sites)



inits <- list(
  ann.fec = dunif(n.years, 0.2, 0.8),
  imm.rec = dunif(n.years, 0.05, 0.95),
  skip.prob = dunif(n.years, 0.15, 0.45),
  sigma.surv = runif(1, 0, 10),      
  sigma.capt = runif(1, 0, 10),
  beta = runif(2, 0, 1),
  z = zinit,
  mean.p = runif(1, 0, 1),
  sigma.proc=runif(n.sites,0,5),
  mean.lambda=runif(n.sites,0.1,2),
  sigma.obs=runif(n.sites,0,10),
  N.est=N.init
  )

#### PARAMETERS TO MONITOR ####
params <- c("pop.size",
            "ann.fec",
            "skip.prob",
            "imm.rec",
            "ann.surv",
            "beta",
            "lambda",
            "pop.growth.rate")

#### MCMC SETTINGS ####
nb <- 10000 #burn-in
ni <- nb + nb #total iterations
nt <- 1  #thin
nc <- 3  #chains


#### COMPILE CONFIGURE AND BUILD ####
Rmodel <- nimbleModel(code = code, constants = const, data = dat, 
                       check = FALSE, calculate = FALSE, inits = inits)
conf <- configureMCMC(Rmodel, monitors = params, thin = nt, 
                       control = list(maxContractions = 1000)) 
Rmcmc <- buildMCMC(conf)  
Cmodel <- compileNimble(Rmodel, showCompilerOutput = FALSE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

#### RUN MCMC ####
t.start <- Sys.time()
#sink("whyisitnotworking.txt")
out <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc, inits = inits,
                setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  
#sink()
t.end <- Sys.time()
(runTime <- t.end - t.start)

saveRDS(out, "AYNAresults.RDS")
out <- readRDS("AYNAresults.RDS")
gelman.diag(out, multivariate = FALSE)

out2 <- rbind(out$chain1, out$chain2, out$chain3) %>% as.data.frame()

juv.surv <- out2 %>% select(contains("ann.surv[1,")) %>% pivot_longer(everything()) 
round(quantile(juv.surv$value, c(0.025, 0.5, 0.975)), 3)

adult.surv <- out2 %>% select(contains("ann.surv[2,")) %>% pivot_longer(everything()) 
round(quantile(adult.surv$value, c(0.025, 0.5, 0.975)), 3)

growth <- out2 %>% select(contains("lambda")) %>% pivot_longer(everything()) %>% 
  mutate(Year = parse_number(str_remove(name, "lambda"))) %>% 
  group_by(Year) %>% 
  summarise(value = median(value))
round(quantile(growth$value, c(0.025, 0.5, 0.975)), 3)

fec <- out2 %>% select(contains("ann.fec")) %>% pivot_longer(everything()) 
round(quantile(fec$value, c(0.025, 0.5, 0.975)), 3)

skip <- out2 %>% select(contains("skip")) %>% pivot_longer(everything()) 
round(quantile(skip$value, c(0.025, 0.5, 0.975)), 3)

rec <- out2 %>% select(contains("rec")) %>% pivot_longer(everything()) 
round(quantile(rec$value, c(0.025, 0.5, 0.975)), 3)
