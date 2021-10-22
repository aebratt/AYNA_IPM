##########################################################################
#
# YELLOW-NOSED ALBATROSS MULTI-STATE INTEGRATED POPULATION MODEL
#
##########################################################################
## written August 2016 by Sarah Converse 
## On top of ME model from October 2015 
## Variation on model by Cat Horswill

## last version YNAL_multievent_IPM_Sarah_v4_allSkippersBreed.R from Sarah on 1 Jan 2017

## modified to run in nimble by Steffen Oppel, April 2018

# # NEED TO DO: (1) streamline model code to reduce number of latent states
# (2) block samplers for parameters to make it more efficient
# (3) experiment with customised distributions for multistate from workshop

## AEB started working on this October 2021

# TODO
# try to get model working in NIMBLE with the data that sarah used
# figure out the data issue with Steffen


# Load necessary libraries
library(tidyverse)
library(data.table)
library(nimble)
library(here)
library(coda)
library(devtools)
library(nimble)
library(mcmcplots)
library(coda)


#############################################################
#
# LOAD AND MANIPULATE DATA
#
##############################################################

# TODO 
# need to get to the bottom of the issues with the data with Steffen
# why does he not think that multievent model is reasonable
# are these data correct???

# simple EH
AYNA <- read_csv(here("Data", "AYNA_simple_encounter_history_1982_2018.csv"))
names(AYNA)
CH<-as.matrix(AYNA[,3:39], dimnames=F)
AYNA$AGE[is.na(AYNA$AGE)]<-1  ## set all NA as 'adult'
apply(CH,2,sum) ### check that there are contacts in every season

# COUNT DATA FOR POPULATION TREND 
AYNA.pop <- read_csv(here("Data", "AYNA_pop_counts_1982_2018.csv")) %>% select(-1)
n.years<-dim(AYNA.pop)[1]		## defines the number of years, starting 1982
n.sites<-dim(AYNA.pop)[2]-1 ## defines the number of study areas
site.year.visited <- (!is.na(AYNA.pop)) + 0

#Number of breeding pairs and chicks fledged from census data
############################################################################# missing counts in 2015 
counts <- read.csv(here("Data", "YNAL_counts.csv"))
nP<-counts[,2]  #breeding pairs
nF<-counts[,3]   #young females (assuming 1:1 sex ratio) 

#Stable age distribution
stable.rate <- c(0.075,0.070,0.064,0.060,0.045,0.034,0.026,0.019,0.015,0.011,0.009,0.007,0.005,0.004,0.003,0.150,0.107,0.296)
stable <- rmultinom(1,450,stable.rate)

######### READ IN ENCOUNTER HISTORY ############################
#read in MARK EH
eh <- read.csv(here("Data", "YNdata_thru15_working.csv"))
eh <- as.matrix(eh)
eh <- eh[,-(1:6)]

YNstate <- eh
colnames(YNstate) = colnames(eh)

# Events are
#0 = unobserved                        -> #6  Unobs
#1 = loaf in colony,                   -> #2  Loaf In
#2 = loaf out of colony,               -> #6  Unobs
#3 = successful breed in colony,       -> #3  Breed-in-S
#4 = breed outside colony,             -> #5  Breed-Out
#5 = failed breed in colony,           -> #4  Breed-in-F
#7 = confirmed dead,                   -> #6  Unobs
#8 = juvenile hatched in colony,       -> #1  Juv
#9 = juvenile hatched out of colony    -> #6  Unobs

# TODO
# i dont think this is doing what it is supposed to be doing
#Condition on being captured first in the study area as a breeder or a juvenile
condcap <- rep(NA,nrow(YNstate))
for(i in 1:nrow(YNstate)){
  # returns the first occasion in which individual was captured in state 3, 5, or 8 ELSE ncol(Ynstate+1)
  condcap[i] <- min(c(which(YNstate[i,] == 3),which(YNstate[i,] == 5),which(YNstate[i,] == 8),(ncol(YNstate)+1)))
  if(condcap[i] == '1'){ # fenceposting
    YNstate[i,] <- YNstate[i,] 
  } else {
    # before it was observed in state 3, 5, or 8 consider it unobserved
    # if it was never oberved in these states it's gonna get omitted from analysis (effectively unobserved throughout)
    YNstate[i,(1:(condcap[i]-1))] <- 0 
  }
}  

#Reassign events as
#1 = Juvenile  
#2 = Loaf In
#3 = Breed In S 
#4 = Breed In F
#5 = Breed Out
#6 = no observation

#reassign events according to matrices in model        
for(i in 1:nrow(YNstate)){
  for(j in 1:ncol(YNstate)){
    if(YNstate[i,j] == '0'){
      YNstate[i,j] <- '6'
    }else if(YNstate[i,j] == '1'){
      YNstate[i,j] <- '2'
    }else if(YNstate[i,j] == '2'){
      YNstate[i,j] <- '6'
    }else if(YNstate[i,j] == '3'){
      YNstate[i,j] <- '3'
    }else if(YNstate[i,j] == '4'){
      YNstate[i,j] <- '5'
    }else if(YNstate[i,j] == '5'){
      YNstate[i,j] <- '4'
    }else if(YNstate[i,j] == '7'){
      YNstate[i,j] <- '6'
    }else if(YNstate[i,j] == '8'){
      YNstate[i,j] <- '1'
    }else if(YNstate[i,j] == '9'){
      YNstate[i,j] <- '6'
    }
  }
}

class(YNstate) <- 'numeric'

#pull out birds never observed in the admissable events (i.e. never NOT 6)
admit <- function(x) length(which(x!=6))
not.admit <- apply(YNstate,1,admit)
Y <- YNstate[-c(which(not.admit==0)),]

#number of individuals and number of years
nind <- nrow(Y)
nyear <- ncol(Y)

#determine first capture occassion for bounding likelihood
get.first <- function(x) min(which(x<6))
first <- apply(Y,1,get.first)

#Determine event at first release for all birds (should be 1, 2 or 3) 
first.event <- rep(NA,nrow(Y))
for(i in 1:nrow(Y)){
  first.event[i] <- Y[i,min(which(Y[i,]<6))]
}
table(first.event)
#get birds that were captured as juveniles
juv <- Y[which(first.event==1),]
#get when first bred
first.breed <- rep(NA,nrow(juv))
for(i in 1:nrow(juv)){
  first.breed[i] <- min(which(juv[i,]==3 | juv[i,]==4 | juv[i,] ==5),(nyear+1))
}
#get those indviduals that recruited
juv.rec <- juv[which(first.breed <(nyear+1)),]
first.breed <- first.cap <- rep(NA,nrow(juv.rec))
#get the age at observed recruitment for those individuals that recruited
for(i in 1:nrow(juv.rec)){
  first.cap[i] <- which(juv.rec[i,]==1)
  first.breed[i] <- min(which(juv.rec[i,]==3 | juv.rec[i,]==4 | juv.rec[i,] ==5),(nyear+1))
}
recruit.obs <- first.breed-first.cap

#age for known age birds 
age <- matrix(data=NA,nrow=nrow(Y),ncol=ncol(Y))
condcap3 <- rep(NA,nrow(Y))
for(i in 1:nrow(Y)){
  condcap3[i] <- min(c(which(Y[i,] == 1),(ncol(Y)+1)))
  for(j in 1:ncol(age)){
    if(j < condcap3[i]){
      age[i,j] <- 'NA'
    }else if(j == condcap3[i]){
      age[i,j] <- 'NA'
    } else {
      age[i,j] <- j-condcap3[i]
    }
  }
}

class(age) <- 'numeric'

# TODO
# why do this
age[is.na(age)]<- nyear			### changed from 33


#Deal with variable years for captures 

# TODO
# where does the effort data come from
#ps which are 0 out of colony
yrs.out <- c(18,22,23,24,25,26,27,32,34)
yrs.in <- c(1:(nyear-1))
yrs.in <- yrs.in[-yrs.out]
nyear.out <- length(yrs.out)
nyear.in <- length(yrs.in)

for(i in 1:nrow(Y)){
  for(j in c(1,yrs.in+1)){
    if(Y[i,j] == 5){
      Y[i,j] <- 6
    }
  }
}

#give this as data to initialize the first time step
# i.e. z.first == first.event
z.first <- rep(NA,nind)
for(i in 1:nind){
  if(first.event[i]==1){
    z.first[i] <- 1
  }else if(first.event[i]==3){
    z.first[i] <- 3
  }else if(first.event[i]==4){
    z.first[i] <- 4
  }  
}

# True states
# 1 - Juvenile      
# 2 - Pre-breed
# 3- Breed-S
# 4- Breed-F
# 5- Skip
# 6- Dead

#Initialize process matrix
z.start <- Y
for(i in 1:nind){
  if(first[i]>1){
    z.start[i,1:(first[i])] <- NA # everything before first is NA
  }
}
first.breed <- rep(NA,nind)
for(i in 1:nind){
  # when is first breeding
  first.breed[i] <- min(c(which(Y[i,] == 3),which(Y[i,] == 4),which(Y[i,] == 5)),(nyear+1)) 
}
for(i in 1:nind){
  for(t in first[i]:nyear){
    if(Y[i,t] == 1){ # when Y is known - z should be NA
      z.start[i,t] <- NA
    }else if (Y[i,t] == 2){ # if we observe you as a loafer
      if(first.breed[i]<t){  # and you are a known breeder
        z.start[i,t] <- 4  # assume z is failed breeder
      }else z.start[i,t] <- 2 # if we've never observed you breeding, assume pre-breeder
    }else if (Y[i,t] == 3){ 
      z.start[i,t] <- NA
    }else if (Y[i,t] == 4){
      z.start[i,t] <- NA
    }else if (Y[i,t] == 5){ # if we observe you as breeding out
      z.start[i,t] <- 4 # assume z is failed breeder
    }else if (Y[i,t] == 6){ # if we don't observe you
      if(first.breed[i]<t){ # and you are a known breeder
        z.start[i,t] <- 4 # assume you are a failed breeder
      }else z.start[i,t] <- 2 # if we've never observed you breeding, assume pre-breeder
    }
  }
}
for(i in 1:nind){
  for(t in first[i]:nyear){
    if(t < first.breed[i]){ 
      if(age[i,t] > 14 & age[i,t]<nyear){ # assume all individuals recruit into breeding pop by 15
        z.start[i,t] <- 4
      }
    }
  }
}
# just making sure
for(i in 1:nind){
  z.start[i,first[i]] <- NA
}

# when Y known - Z data is known too
# skip and dead are unobserved
z.data <- z.start
for(i in 1:nind){
  for(t in first[i]:nyear){
    if(Y[i,t] == 1){
      z.data[i,t] <- 1
    }else if (Y[i,t] == 2){
      z.data[i,t] <- NA
    }else if (Y[i,t] == 3){
      z.data[i,t] <- 3
    }else if (Y[i,t] == 4){
      z.data[i,t] <- 4
    }else z.data[i,t] <- NA 
  }
}
for(i in 1:nind){
  z.data[i,first[i]] <- NA
}