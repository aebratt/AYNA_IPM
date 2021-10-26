# TODO
# reorganize repo
# figure out data issue
# adapt model
# compile todolist throughout these files

library(tidyverse)
library(data.table)
library(nimble)
library(devtools)
library(mcmcplots)
library(coda)
library(here)
library(nimbleEcology)
library(foreach)
library(doParallel)
source(here("Models", "AEB_prepData.R"))
source(here("Models", "AEB_msOnly.R"))
source(here("Models", "AEB_nimbleFunctions.R"))

#### MCMC SETTINGS ####
nb <- 10000 #burn-in
ni <- 100000 + nb #total iterations
nt <- 1  #thin
nc <- 3  #chains
adaptInterval = 100

YNAL.Consts <- list(z.first=z.first,
                    first=first,
                    nind=nind,
                    nyear=nyear,
                    nyear.in=nyear.in,
                    yrs.in=yrs.in,
                    nyear.out=nyear.out,
                    yrs.out=yrs.out, 
                    #stable=stable
                    nstates = 6,
                    nages = 34 # TODO - remove hardcoding
)


zero <- numeric(dim(Y)[1])

# TODO - should these be data or constants
YNAL.Data <- list(Y=Y,
                  age=age,
                  zero = zero
                  #z=z.data#, 
                  #nF=nF, 
                  #nP=nP
)

# TODO - add more here
parameters <- c("S",
                "B",
                "R.age",
                "F",
                "p"
                #"n.breeders",
                #"lambda.t",
                #"mean.lambda"
) 

# 3. Specify initialisation values
S.rand.st <- matrix(runif(2*(nyear-1),2,3),nrow=2,ncol=nyear-1)
B.rand.st <- matrix(runif(2*(nyear-1),2,3),nrow=2,ncol=nyear-1)

inits <- list (#z=z.start,
               S.rand=S.rand.st,
               int.S=runif(2,2,4),
               sigma.S=runif(2),
               int.B=runif(2,0,2),
               sigma.B=runif(2),
               B.rand=B.rand.st,
               F.rand=runif(nyear-1),
               int.F=runif(1,1,3),
               sigma.F=runif(1)
               )

# cores=detectCores()
# cl <- makeCluster(nc, setup_strategy = "sequential") #not to overload your computer
# registerDoParallel(cl)

# foreach(i = 1:nc) %dopar% { #scenarios picked
#   library(nimble)
#   library(here)
#   source(here("Models", "AEB_nimbleFunctions.R"))

#### COMPILE CONFIGURE AND BUILD ####
Rmodel <- nimbleModel(code = YNAL.IPM, constants = YNAL.Consts, data = YNAL.Data, 
                      check = FALSE, calculate = FALSE, inits = inits)
conf <- configureMCMC(Rmodel, monitors = parameters, thin = nt, 
                      control = list(maxContractions = 1000, adaptInterval = adaptInterval)) 
# lots of initial model checking you can do by exploring conf1
# if you wanted to change samplers this is where you would do that

# MODEL EXPLORATION

# TODO
# what is going on here
# p.link 6, 7, 8, 9, t = 9
conf$printSamplers(type = "posterior") # check sampler defaults


##### CONSIDER BLOCK SAMPLING OF PARAMETERS THAT WILL ALWAYS BE EVALUATED TOGETHER ####

# blocking is smart for parameters that are correlated

# TODO
# could block wrt to time
# using RW block samplers
# TODO
# consider using AF slice samplers??
# TODO
# consider blocking on another dimension???
  # p link[1:10, 1-(nyears-1)]
  # S.rand[1:2,  1-(nyears-1)]
  # F rand[1-(nyears-1)]
  # B rand[1:2,  1-(nyears-1)
  # R rand[4:14]

  ## block all the capture prob variables
  conf$removeSamplers('int.p')
  conf$addSampler(target = "int.p[1:10]", type="AF_slice")
  conf$printSamplers("int.p")
# 
  conf$removeSamplers('sigma.p')
  conf$addSampler(target = "sigma.p[1:10]", type="AF_slice")
  conf$printSamplers("sigma.p")

  # TODO
  # for some reason this creates an error?
#   conf$removeSamplers('p.link')
# for(link in 1:10){
#     conf$addSampler(target = paste0("p.link[",link,", 1:",(nyear-1),"]"), type="AF_slice")
# }#s
#   conf$printSamplers("p.link")
#   
#   ## block all the survival variables

  conf$removeSamplers('int.S')
  conf$addSampler(target = "int.S[1:2]", type="AF_slice")
  conf$printSamplers("int.S")

  conf$removeSamplers('sigma.S')
  conf$addSampler(target = "sigma.S[1:2]", type="AF_slice")
  conf$printSamplers("sigma.S")
#   
  conf$removeSamplers('S.rand')
  for(s in 1:2){
    conf$addSampler(target = paste0("S.rand[",s,", 1:",(nyear-1),"]"), type="AF_slice")
  }#s
  conf$printSamplers("S.rand")
#   
#   ## block all the fecundity variables  
  conf$removeSamplers('F.rand')
  conf$addSampler(target = paste0("F.rand[1:",(nyear-1),"]"), type="AF_slice")
  conf$printSamplers("F.rand")
# 
#   ## block all the breeding variables  
#   
  conf$removeSamplers('int.B')
  conf$addSampler(target = "int.B[1:2]", type="AF_slice")
  conf$printSamplers("int.B")

  conf$removeSamplers('sigma.B')
  conf$addSampler(target = "sigma.B[1:2]", type="AF_slice")
  conf$printSamplers("sigma.B")
#   
  conf$removeSamplers('B.rand')
  for(b in 1:2){
    conf$addSampler(target = paste0("B.rand[",b,", 1:",(nyear-1),"]"), type="AF_slice")
  }#s
  conf$printSamplers("B.rand")
#   
#   ## block all the recruiting variables  

  conf$removeSamplers('R.rand')
  conf$addSampler(target = "R.rand[4:14]", type="AF_slice")
  conf$printSamplers("R.rand")

  ## block all the fledgling count variables  
  #conf$removeSamplers('nF')
  #conf$removeSamplers('tau.Fl')
  #conf$addSampler(target = c('nF','tau.Fl'), type = "AF_slice")
  ## block all the pair count variables  
  #conf$removeSamplers('nP')
  #conf$removeSamplers('tau.Pr')
  #conf$removeSamplers('sigma.Pr')
  #conf$addSampler(target = c('nP','tau.Pr','sigma.Pr'), type = "AF_slice")
  conf

############

Rmcmc <- buildMCMC(conf)  
Cmodel <- compileNimble(Rmodel, showCompilerOutput = FALSE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

#### RUN MCMC ####
t.start <- Sys.time()
out <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc, inits = inits,
               setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  
t.end <- Sys.time()
(runTime <- t.end - t.start)

# we ran the murres for 110000
# takes 1 minute to run 1 chain of 10 iterations
# est 110000/10 * 1 = 11000 minutes to run to completion (1 chain)
# est 7.5 days to run to completion (1 chain)
# obviously going to be more when adding the count data, 
# and any kind of structure to the model
# so definitely need to pursue blocking

# saveRDS(out, here("Models", paste("out-",i,".RDS", sep = "")))
#t.end <- Sys.time()
#(runTime <- t.end - t.start)

# } # foreach - scenarios picked (i)
# stopCluster(cl)
