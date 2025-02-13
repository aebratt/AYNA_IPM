############################################################################
######## DATA PREPARATION FOR INTEGRATED POPULATION MODEL     ##############
############################################################################

### written by Steffen Oppel in November 2021 - adapted from TRAL code https://github.com/steffenoppel/TRAL_IPM
### steffen.oppel@rspb.org.uk
### uses existing CMR and breeding database to extract data

library(tidyverse)
library(lubridate)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select



#############################################################################
##   1. SPECIFY THE SPECIES AND START YEAR FOR WHICH YOU WANT A SUMMARY ####
#############################################################################
## SPECIFY THE SPECIES AND START YEAR FOR SURVIVAL MODEL
SP<-"AYNA"
start<-1978  ## for CMR data
IPMstart<-2007 ## for count and breeding success data - only start reliably in 2008


###################################################################################
##   2. READ IN DATA FROM DATABASE AND FILTER DATA FOR SPECIES OF INTEREST ####
###################################################################################

## run the RODBC import of nest and count data in a 32-bit version of R
#system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database\\RODBC_count_import.r")), wait = TRUE, invisible = FALSE, intern = T)
#system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database\\RODBC_count_import.r")), wait = TRUE, invisible = FALSE, intern = T)
#try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database"), silent=T)
load("GOUGH_seabird_data.RData")


## filter data for the selected species
nests<-nests %>% filter(Species==SP)
counts<-counts %>% filter(Species==SP)


## look at data
head(nests)  ## nest monitoring data
head(counts)  ## seabird count data



### CHECK DATA
unique(nests$Species)
unique(nests$Colony)
nests %>% filter (is.na(Colony))
nests<-nests %>% mutate(Colony=if_else(is.na(Colony),"Area 1",as.character(Colony)))
unique(nests$Colony)
unique(nests$StageFound)
unique(nests$LastStage)

unique(counts$Species)
unique(counts$Colony)
unique(counts$Breed_Stage)
counts %>% filter(Breed_Stage=="FLED")
unique(counts$Cohort)



#############################################################################
##   3. PREPARE THE BREEDING SUCCESS DATA FROM NEST RECORDS #################
#############################################################################
#try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM"), silent=T)


### remove only partially monitored nests
exclude <- nests %>%
  filter(LastStage=="INCU") %>%
  filter(SUCCESS==1)


### summary of breeding success per year from nests
FECUND<-nests %>% filter(Year<2021) %>% filter(Species==SP) %>% mutate(count=1) %>%
  mutate(Colony=if_else(is.na(Colony),"Area 1",as.character(Colony))) %>%
  filter(!NestID %in% exclude$NestID) %>%
  group_by(Year,Colony) %>%
  summarise(n_nests=sum(count),BREED_SUCC=mean(SUCCESS, na.rm=T))

FECUND<-FECUND %>% filter(Colony=="Area 1")


### PLOT TO SPOT ANY OUTLIERS OF BREEDING SUCCESS
ggplot(FECUND, aes(x=Year,y=BREED_SUCC)) +geom_point(size=2, color='darkred')+geom_smooth(method='lm') 



#############################################################################
##   4. PREPARE THE POPULATION COUNT DATA FROM COUNT RECORDS ################
#############################################################################
#try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdBreedingSuccess"), silent=T)

head(counts)

### summary of population counts of breeding pairs per year and colony
POPSIZE<-counts %>% filter(Species==SP) %>%
  mutate(Colony= as.character(Colony)) %>%
  mutate(Year=year(Date)) %>%
  filter(Year>IPMstart) %>%
  filter(Breed_Stage=="INCU") %>%
  filter(Cohort %in% c("INCU","TERR","AON")) %>%
  group_by(Year,Colony) %>%
  summarise(N=sum(Number, na.rm=T)) %>%
  spread(key=Colony, value=N)
POPSIZE

### PLOT TO SPOT ANY OUTLIERS OF COUNT DATA
POPSIZE %>% gather(key="Colony",value="N",-Year) %>% group_by(Year) %>% summarise(TOT=sum(N, na.rm=T)) %>% #summarise(avg=median(TOT))
ggplot(aes(x=Year,y=TOT)) +geom_point(size=2, color='darkred')+
  geom_smooth(method="lm",formula = y ~ x + I(x^2), colour="red", size=1.5, fill="indianred1", alpha=0.2) +
  theme(panel.background=element_rect(fill="white", colour="black"),  
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),  
        strip.text.x=element_text(size=18, color="black"),  
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_blank()) 


### summary of population counts of fledglings per year and colony
CHICKCOUNT<-counts %>% filter(Species==SP) %>%
  mutate(Colony= as.character(Colony)) %>%
  mutate(Year=year(Date)) %>%
  filter(Year>IPMstart) %>%
  filter(month(Date)<6) %>% ## need to exclude small chick counts in Dec
  filter(Breed_Stage %in% c("CHIC","FLED")) %>%
  filter(Cohort %in% c("CHIC","FLED")) %>%
  group_by(Year,Colony) %>%
  summarise(N=sum(Number, na.rm=T)) %>%
  ungroup() %>%
  spread(key=Colony, value=N)
CHICKCOUNT


### MAKE SURE BOTH MATRICES HAVE IDENTICAL DIMENSIONS
## adult counts in 11 areas, chick counts only in 4 areas
## note offset in year - adult count in year X is related to chick count in year X+1
dim(POPSIZE)
dim(CHICKCOUNT)

ADCOUNT <- POPSIZE %>% select(Year,`Area 1`,`Area 2`,`Area 3`,`Area 10`) %>%
  filter(Year<2021) %>% ### this cohort is not done yet
	arrange(Year)

CHICKCOUNT <- CHICKCOUNT %>% select(Year,`Area 1`,`Area 2`,`Area 3`,`Area 10`) %>%
  filter(Year>min(POPSIZE$Year)) %>%
  mutate(Year=Year-1) %>%
	arrange(Year)

### fill in missing counts from 2012
CHICKCOUNT<-counts %>% filter(Species==SP) %>%
  mutate(Colony= as.character(Colony)) %>%
  mutate(Year=year(Date)) %>%
  filter(Year==2011) %>%
  filter(month(Date)==12) %>% ## no large chick counts in March 2012, therefore take chick counts in Dec 2011
  filter(Breed_Stage %in% c("CHIC","FLED")) %>%
  filter(Cohort %in% c("CHIC","FLED")) %>%
  group_by(Year,Colony) %>%
  summarise(N=sum(Number, na.rm=T)) %>%
  spread(key=Colony, value=N) %>%
  select(Year,`Area 1`,`Area 2`,`Area 3`,`Area 10`) %>%
  bind_rows(CHICKCOUNT) %>%
  arrange(Year)

BS<-CHICKCOUNT/ADCOUNT

dim(ADCOUNT)
dim(CHICKCOUNT)


### PLOT TO SPOT ANY OUTLIERS OF COUNT DATA
BS %>% mutate(Year=ADCOUNT$Year) %>% gather(key="Colony",value="N",-Year) %>% group_by(Year) %>% summarise(SUCCESS=mean(N, na.rm=T)) %>% ##ungroup() %>% summarise(avg=median(SUCCESS, na.rm=T),low=min(SUCCESS, na.rm=T),hi=max(SUCCESS, na.rm=T))
ggplot(aes(x=Year,y=SUCCESS)) +geom_point(size=2, color='darkred')+geom_smooth(method='lm') + 
theme(panel.background=element_rect(fill="white", colour="black"),  
      axis.text=element_text(size=18, color="black"), 
      axis.title=element_text(size=20),  
      strip.text.x=element_text(size=18, color="black"),  
      strip.background=element_rect(fill="white", colour="black"), 
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.border = element_blank()) 





#############################################################################
##   5. PREPARE THE MARK-RECAPTURE DATA FOR SURVIVAL ANALYSIS ###############
#############################################################################

## run the RODBC import of CMR data in a 32-bit version of R
#system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\RODBC_CMR_import_AYNA.R")), wait = TRUE, invisible = FALSE, intern = T)
#system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\RODBC_CMR_import_AYNA.R")), wait = TRUE, invisible = FALSE, intern = T)
#try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM"), silent=T)
load("GOUGH_seabird_CMR_data.RData")


### COPIED FROM C:\STEFFEN\RSPB\UKOT\Gough\ANALYSIS\SeabirdSurvival\TRAL_survival_marray.r

## filter data for the selected species
contacts<-contacts %>% filter(SpeciesCode==SP) %>% ## %>% filter(Location %in% c("Area 1","Not Specified")) - filter by location optional
  mutate(Contact_Year=ifelse(!is.na(Date_Time) & (Contact_Year != year(Date_Time)),year(Date_Time),Contact_Year)) %>% ## fix years that were initially used as 'season'
  mutate(Contact_Year=ifelse(((Age %in% c("Chick","Fledgling") & is.na(Date_Time) & (Contact_Year == as.integer(substr(Contact_Season,1,4))))),  #
                             as.integer(Contact_Year)+1,as.integer(Contact_Year))) ## fix years that were initially used as 'season'


ages<-ages %>% filter(SpeciesCode==SP)
bands<-bands %>% filter(SpeciesCode==SP)

head(contacts)  ## CMR data
dim(contacts)
sort(unique(contacts$Contact_Year))
which(!((start:2021) %in% sort(unique(contacts$Contact_Year))))
(start:2021)[which(!(start:2021 %in% sort(unique(contacts$Contact_Year))))]
unique(contacts$Contact_Year)

contacts %>% filter(BirdID=="GO-16-18-556")


#############################################################################
##   6. AGE ASSIGNMENT OF BIRDS FOR SURVIVAL ANALYSIS ###############
#############################################################################

### EXTRACT AGE AT DEPLOYMENT FROM DATABASE
deploy_age<-contacts %>% arrange(BirdID, Date_Time,Contact_Year) %>%
  mutate(AGE=ifelse(Age %in% c("Chick","Fledgling"),0,1)) %>%
  arrange(BirdID,Contact_Year,Date_Time) %>%
  group_by(BirdID) %>%
  summarise(MIN_AGE=min(AGE,na.rm=T), MAX_AGE=max(AGE,na.rm=T), FIRST_AGE=first(Age,na.rm=T), FIRST_Date=first(Date_Time,na.rm=T), FIRST_YEAR=min(Contact_Year,na.rm=T)) %>% #filter(is.na(FIRST_YEAR))
  mutate(FIRST_AGE=ifelse(FIRST_AGE=="Unknown" & month(FIRST_Date)>6,"Adult", as.character(FIRST_AGE))) %>%  ### unknowns marked after June were not chicks
  mutate(FIRST_YEAR=ifelse(is.na(FIRST_YEAR),year(FIRST_Date),FIRST_YEAR))
head(deploy_age)
dim(deploy_age)
unique(deploy_age$FIRST_YEAR)


MISSAGE<-deploy_age %>% filter(is.na(FIRST_AGE)) %>%   left_join(bands, by="BirdID") %>%
  select(BirdID, Band_Number,MIN_AGE,FIRST_Date,FIRST_YEAR)
dim(MISSAGE)


### DEFINE A YEAR (=SEASON) FROM SEPT X to JUNE X+1 because birds breed during summer from Sept - April
## use Contact_Season as encounter occasion grouping variable

contacts %>% #filter(is.na(Date_Time)) %>%
  #filter(is.na(Contact_Year)) %>%
  filter(is.na(Contact_Season)) 


contacts<-contacts %>%
  mutate(Contact_Season=if_else(is.na(Contact_Season),if_else(Age=="Chick",paste(Contact_Year-1,"-",substr(Contact_Year,3,4), sep =""),
                                                              paste(Contact_Year,"-",as.integer(substr(Contact_Year,3,4))+1, sep = "")),
                                Contact_Season))
dim(contacts)
head(contacts)
sort(unique(contacts$Contact_Season))

### ASSIGN AGE TO BIRDS WHERE THIS IS NOT SPECIFIED
## include a column with continuous age 

contacts<-contacts %>%
  left_join(deploy_age, by="BirdID") %>%
  mutate(AGE=ifelse(Age=="Adult",1,ifelse(Age %in% c("Chick","Fledgling"),0,NA))) %>%    ### certain assignments based on provided age
  mutate(AGE=ifelse(is.na(AGE),ifelse(FIRST_AGE=="Adult",1,if_else(Contact_Year>FIRST_YEAR,1,0)),AGE)) %>%    ### certain assignments based on provided age
  mutate(AGE=ifelse(is.na(AGE), ifelse(Sex %in% c("Male","Female"),1,NA),AGE)) %>%       ### inferred assignment from sex info - only adults can be sexed
  mutate(ContAge=ifelse(FIRST_AGE %in% c("Chick","Fledgling"),Contact_Year-FIRST_YEAR,Contact_Year-FIRST_YEAR+5)) #%>%      ### continuous age since first deployment, at least 5 years for birds marked as 'adult'

contacts %>% filter(is.na(AGE))
contacts %>% filter(is.na(ContAge))
contacts %>% filter(ContAge==1)




################################################################################################
##   7. REMOVE BIRDS FROM OUTSIDE THE STUDY AREAS AND BEFORE 1978  ###############
##############################################################################################
unique(contacts$Location)
########## CREATE A LOOP OVER EVERY BIRD TO CHECK WHETHER THEY WERE EVER RECORDED IN STUDY AREAS
STUDY_AREAS<- c("Area 1","Not Specified","Between the base and seal beach","Area 8","Area 3","Area 3/10","Area 10", "Prion Cave","Tumbledown")
allbirds<-unique(contacts$BirdID)
fixed_contacts<-data.frame()

for (xid in allbirds){
  xcont<-contacts %>% filter(BirdID==xid) %>% arrange(Date_Time)
  xcont$INSIDE<-ifelse(xcont$Location %in% STUDY_AREAS,1,0)
  xcont$INSIDE<-ifelse(xcont$Location == "Not Specified" & xcont$Contact_Year>2014,1,xcont$INSIDE) 
  if(sum(xcont$INSIDE, na.rm=T)>0){fixed_contacts<-bind_rows(fixed_contacts ,xcont)}
}
dim(contacts)
dim(fixed_contacts)
length(unique(fixed_contacts$BirdID))
length(allbirds)


### CHECK WHAT BIRDS WERE RINGED AS CHICKS BEFORE 1978
oldchicks<-fixed_contacts %>% filter(Contact_Year<=start) %>% filter(ContAge<2)
fixed_contacts %>% filter(BirdID %in% oldchicks$BirdID)


### REMOVE RECORDS FROM BEFORE THE SET START YEAR AND BIRDS FIRST MARKED IN LAST YEAR
contacts<-fixed_contacts %>%
  filter(year(Date_Time)>=start)
dim(contacts)
unique(contacts$FIRST_AGE)

sort(unique(contacts$Contact_Season))

all.seasons <- paste(1978:2021, "-", 
                     c(79:99, 
                       "00", "01", "02", "03", "04", "05", "06", "07", "08", "09", 
                       10:22), 
                     sep = "")
all.seasons
which(!(all.seasons %in% sort(unique(contacts$Contact_Season))))
no.contact.seasons <- all.seasons[which(!(all.seasons %in% sort(unique(contacts$Contact_Season))))]
no.contact.seasons

################################################################################################
##   8. TRY TO QUANTIFY OBSERVATION EFFORT  ###############
##############################################################################################
## try to determine years with high and low detection probability

contacts %>% mutate(count=1) %>% group_by(Contact_Season) %>% summarise(n=sum(count)) %>%
  ggplot() + geom_bar(aes(x=Contact_Season,y=n), stat="identity")

## calculate number of individuals that had been marked by a given year
n_exist<-deploy_age %>% mutate(count=1) %>% rename(Contact_Year=FIRST_YEAR) %>%
  mutate(FIRST_AGE=if_else(FIRST_AGE=="Fledgling","Chick",FIRST_AGE)) %>%
  mutate(Contact_Year=ifelse(FIRST_AGE=="Chick",Contact_Year-1,Contact_Year)) %>%
  group_by(Contact_Year,FIRST_AGE) %>%
  summarise(N_marked=sum(count)) %>%
  arrange(Contact_Year) %>%
  spread(key=FIRST_AGE, value=N_marked, fill=0) %>%
  ungroup() %>%
  mutate(N_marked=Adult+Chick) %>%
  mutate(N_all = cumsum(N_marked)) %>%
  bind_rows(tibble(Contact_Year=2021,N_marked=0,N_all=0)) %>%
  mutate(N_all=if_else(Contact_Year==2021,dplyr::lag(N_all),N_all))
n_exist$N_all[1] = n_exist$N_marked[1]
n_exist$N_all[2] = (n_exist$N_all[1]*(0.92^4)) + n_exist$N_marked[2]
n_exist$N_all[3] = (n_exist$N_all[2]*(0.92^2)) + n_exist$N_marked[3]
n_exist$N_all[4] = (n_exist$N_all[3]*(0.92^3)) + n_exist$N_marked[4]
for (y in 5:dim(n_exist)[1]) {
  n_exist$N_all[y] = ((n_exist$Adult[y-1]+n_exist$N_all[y-2])*0.92) + (n_exist$Chick[y-1]*0.85) + n_exist$N_marked[y]
}
tail(n_exist)


goodyears<-contacts %>% group_by(Contact_Year) %>% summarise(n=length(unique(BirdID))) %>%
  left_join(n_exist, by='Contact_Year') %>%
  mutate(prop.seen=n/N_all) %>%
  mutate(p.sel=if_else(prop.seen>0.13,2,1))
tail(goodyears)

ggplot(goodyears) + geom_histogram(aes(x=prop.seen), binwidth=0.01)

goodyears %>% mutate(Effort=if_else(prop.seen<0.15,"low","high")) %>%
ggplot() + geom_bar(aes(x=Contact_Year,y=prop.seen, fill=Effort), stat="identity") + 
  
  labs(x = "Year",
       y = "Annual proportion of ringed birds recorded",
       fill = "Effort classification") +

  theme(panel.background=element_rect(fill="white", colour="black"),  
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),  
        strip.text.x=element_text(size=18, color="black"),  
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position=c(0.15, 0.90),
        panel.border = element_blank()) 

#ggsave("C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\AYNA_IPM\\FigS1.jpg", width=9, height=6)



#############################################################################
##   9. CREATE MATRIX OF ENCOUNTERS AND AGES ###############
#############################################################################
head(contacts)

sort(unique(contacts$Contact_Season))

### SIMPLE BINARY ENCOUNTER HISTORY FOR CHICKS AND ADULTS
AYNA_CHICK<- contacts %>% mutate(count=1) %>%
  filter(FIRST_AGE %in% c("Chick","Fledgling")) %>%
  group_by(BirdID,Contact_Season) %>%
  summarise(STATE=max(count)) %>%
  spread(key=Contact_Season, value=STATE, fill=0) %>%
  arrange(BirdID)
dim(AYNA_CHICK)  

## THIS MATRIX MUST BE PADDED BY YEARS WITH 0 CHICKS RINGED
colnames(AYNA_CHICK)[-1]
length(colnames(AYNA_CHICK)[-1])
which(!(all.seasons %in% colnames(AYNA_CHICK)[-1]))
all.seasons[which(!(all.seasons %in% colnames(AYNA_CHICK)[-1]))]
pad.chicks.vec <- all.seasons[which(!(all.seasons %in% colnames(AYNA_CHICK)[-1]))]
pad.chicks.vec

pad.chicks <- data.frame(matrix(0, nrow = dim(AYNA_CHICK)[1], ncol = length(pad.chicks.vec)))
colnames(pad.chicks) <- pad.chicks.vec
head(pad.chicks)

AYNA_CHICK.tmp <- cbind(AYNA_CHICK, pad.chicks)
colnames(AYNA_CHICK.tmp) <- c(colnames(AYNA_CHICK), pad.chicks.vec)
colnames(AYNA_CHICK.tmp)

AYNA_CHICK <- AYNA_CHICK.tmp %>% select("BirdID", sort(colnames(.)))
colnames(AYNA_CHICK)[-1]
dim(AYNA_CHICK)

### identify number of chicks ringed every year
phi.juv.possible<-AYNA_CHICK %>% gather(key='Year', value='count',-BirdID) %>% group_by(Year) %>% summarise(N=sum(count)) %>%
  mutate(JuvSurv=if_else(N<25,0,1)) %>%
  mutate(JuvSurv=if_else(Year>2018,0,JuvSurv)) ### not possible yet to estimate juvenile survival after 2018

phi.juv.possible$JuvSurv


AYNA_AD<- contacts %>% mutate(count=1) %>%
  group_by(BirdID,FIRST_AGE,Contact_Season) %>%
  summarise(STATE=max(count)) %>%
  spread(key=Contact_Season, value=STATE, fill=0) %>%
  filter(FIRST_AGE %in% c("Adult")) %>%    ### filter after spread to ensure that years without any adult contacts (1984, 2003, 2005) are included in matrix
  ungroup() %>%
  select(-FIRST_AGE) %>%
  arrange(BirdID)
dim(AYNA_AD)

colnames(AYNA_AD)[-1]
length(colnames(AYNA_AD)[-1])
which(!(all.seasons %in% colnames(AYNA_AD)[-1]))
all.seasons[which(!(all.seasons %in% colnames(AYNA_AD)[-1]))]
pad.adults.vec <- all.seasons[which(!(all.seasons %in% colnames(AYNA_AD)[-1]))]
pad.adults.vec

pad.adults <- data.frame(matrix(0, nrow = dim(AYNA_AD)[1], ncol = length(pad.adults.vec)))
colnames(pad.adults) <- pad.adults.vec
head(pad.adults)

AYNA_AD.tmp <- cbind(AYNA_AD, pad.adults)
colnames(AYNA_AD.tmp) <- c(colnames(AYNA_AD), pad.adults.vec)
colnames(AYNA_AD.tmp)

AYNA_AD <- AYNA_AD.tmp %>% select("BirdID", sort(colnames(.)))
colnames(AYNA_AD)[-1]
dim(AYNA_AD)

### CONVERT TO SIMPLE MATRICES WITHOUT BIRD ID COLUMN
CH.J<-as.matrix(AYNA_CHICK[,2:dim(AYNA_CHICK)[2]])
CH.A<-as.matrix(AYNA_AD[,2:dim(AYNA_CHICK)[2]])


### IDENTIFY WHICH CHICKS WERE EVER RECAPTURED AND PUT THOSE IN A SEPARATE ENCOUNTER HISTORY
cap <- apply(CH.J, 1, sum)
ind <- which(cap >= 2)
CH.J.R <- CH.J[ind,]    # Juvenile CH recaptured at least once
CH.J.N <- CH.J[-ind,]   # Juvenile CH never recaptured
dim(CH.J.R)

# FOR THOSE CHICKS THAT WERE RECAPTURED, Remove first capture and add the rest to the adult encounter history
first <- numeric()
for (i in 1:dim(CH.J.R)[1]){
  first[i] <- min(which(CH.J.R[i,]==1))
}
CH.J.R1 <- CH.J.R
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R1[i,first[i]] <- 0
}

# Add grown-up juveniles to adults FOR COMPLETE ENCOUNTER HISTORY TO BUILD ADULT MARRAY
CH.A.m <- rbind(CH.A, CH.J.R1)


# Create ENCOUNTER HISTORY matrix for juveniles, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.J.R1)[1]){
  second[i] <- min(which(CH.J.R1[i,]==1))
}
CH.J.R2 <- matrix(0, nrow = dim(CH.J.R)[1], ncol = dim(CH.J.R)[2])
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R2[i,first[i]] <- 1
  CH.J.R2[i,second[i]] <- 1
}



#############################################################################
##   10. CONVERT TO MARRAY ###############
#############################################################################

# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}


# Create TWO MARRAYS from the capture-histories: one for chicks and one for adults (including those chicks that were ever recaptured)

### ADULT MARRAY
adult.marray <- marray(CH.A.m)

# Create m-array for the chicks that were ultimately recaptured, but only up to the first recapture (the rest is included in the adult.marray)
CH.J.R.marray <- marray(CH.J.R2)
diag(CH.J.R.marray)  ## this should not contain any 1s because that would mean a chick is caught in the year after it was marked

# The last column ought to show the number of juveniles not recaptured again and should all be zeros, since all of them are released as adults
CH.J.R.marray[,dim(CH.J)[2]] <- 0

# Create the m-array for chicks never recaptured and add it to the recaptured chicks m-array
CH.J.N.marray <- marray(CH.J.N)
chick.marray <- CH.J.R.marray + CH.J.N.marray 



### CALCULATE THE PROPORTION OF CHICKS RECRUITING AT A CERTAIN AGE

recruit.age<-function(CH){
  n.occasions <- dim(CH)[2]-1
  total<-as.numeric()
  recruit.mat<- matrix(data = 0, ncol = n.occasions, nrow = n.occasions-1)
  # Calculate the age proportion of returned individuals at each time period
  for (t in 1:(n.occasions-1)){
    total[t] <- sum(CH[t,])
    for (col in (t+1):n.occasions){
      recruit.mat[t,col-t]<- CH[t,col]
    }
  }
  return(list(REC=recruit.mat,TOT=total))
}

RECRUIT.AGE.MAT<-recruit.age(CH.J.R.marray)
RECRUIT.AGE<-data.frame(age=seq(1,dim(RECRUIT.AGE.MAT$REC)[2],1),N=apply(RECRUIT.AGE.MAT$REC,2,sum)) %>%
  mutate(prop=N/sum(RECRUIT.AGE.MAT$TOT))

sum(RECRUIT.AGE.MAT$TOT)
sum(RECRUIT.AGE$N)

ggplot(RECRUIT.AGE) + geom_bar(aes(x=age,y=prop),stat='identity', fill='cornflowerblue') + 
  ylab("Proportion of AYNA first seen on Gough") + 
  xlab("Age in years") + 
  theme(panel.background=element_rect(fill="white", colour="black"),  
        axis.text=element_text(size=16, color="black"),  
        axis.title=element_text(size=18),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_blank()) 



### IDENTIFY THE CONTACTS OF AGE <4 TO DOUBLE-CHECK IN DATABASE
DOUBLE_CHECK<-contacts %>% filter(ContAge %in% c(1,2)) %>%
  arrange(BirdID,Date_Time)

contacts %>% filter(BirdID %in% unique(DOUBLE_CHECK$BirdID)) %>%
  arrange(BirdID,Date_Time)




#############################################################################
##   11. SAVE WORKSPACE ###############
#############################################################################
# TODO - change this
#setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM")
#save.image("AYNA_IPM_input.marray.RData")
