# contains some code to clean and prepare the data for analysis. Run before doing anything else

setwd("C:/Users/jamie/OneDrive/Documents/MMath Project/tommy_20200818_baseline_outcomes")
library("survival")
library("Hotelling")
library("mvtnorm")
library("tidyverse")


# ======= Load in Data ==========

demographics <- read.csv("demographicsnew.csv")
famhist <- read.csv("family_history.csv")
famhistmale <- read.csv("male_family_history_lifestyle.csv")
postpregnancies <- read.csv("PostClinicPregnancies.csv")
prepregnancies <- read.csv("PreClinicPregnancies.csv")
firstpregnancy <- read.csv("FirstPregnancyPostClinic.csv")
fertility <- read.csv("fertility.csv")
treatments <-read.csv("treatments.csv")
drugs <- read.csv("recreational_drugs.csv")


# ======= Change the weird IDMAT names =======

names(postpregnancies)[names(postpregnancies) == "ï..IDMAT"] <- "IDMAT"
names(demographics)[names(demographics) == "ï..IDMAT"] <- "IDMAT"
names(prepregnancies)[names(prepregnancies) == "ï..IDMAT"] <- "IDMAT"
names(famhist)[names(famhist) == "ï..IDMAT"] <- "IDMAT"
names(famhistmale)[names(famhistmale) == "ï..IDMAT"] <- "IDMAT"
names(firstpregnancy)[names(firstpregnancy) == "ï..IDMAT"] <- "IDMAT"
names(fertility)[names(fertility) == "ï..IDMAT"] <- "IDMAT"
names(treatments)[names(treatments) == "ï..IDMAT"] <- "IDMAT"
names(drugs)[names(drugs) == "ï..IDMAT"] <- "IDMAT"


# ======= Some Data Formatting =======

demographics$DateOfFirstConsultationMAT <- as.Date(demographics$DateOfFirstConsultationMAT, "%d/%m/%Y")
firstpregnancy$ConceptionDate <- as.Date(firstpregnancy$ConceptionDate, "%d/%m/%Y")
firstpregnancy$DateOfFirstConsultation <- as.Date(firstpregnancy$DateOfFirstConsultation, "%d/%m/%Y")
postpregnancies$ConceptionDate <- as.Date(postpregnancies$ConceptionDate, "%d/%m/%Y")
postpregnancies$expectedDueDate <- as.Date(postpregnancies$expectedDueDate, "%d/%m/%Y")

famhist$Height <- as.numeric(famhist$Height)
famhist$Weight <- as.numeric(famhist$Weight)
famhist$BMI <- as.numeric(famhist$BMI)
demographics$AgeMAT <- as.numeric(demographics$AgeMAT)
firstpregnancy$methodOfConception <- as.numeric(firstpregnancy$methodOfConception)

yesnotoindicator <- function(string){
  if(string=="Yes"){
    return(1)
  }
  if(string=="No"){
    return(0)
  }
  else(return(string))
}

demographics$SmokingMAT <- sapply(demographics$SmokingMAT, yesnotoindicator)
demographics$SmokingMale <- sapply(demographics$SmokingMale, yesnotoindicator)
demographics$AlcoholMAT <- sapply(demographics$AlcoholMAT, yesnotoindicator)
demographics$SmokingMale <- sapply(demographics$SmokingMale, yesnotoindicator)
fertility$FertilityTreatment <- sapply(fertility$FertilityTreatment, yesnotoindicator)
fertility$PeriodRegular <- sapply(fertility$PeriodRegular, yesnotoindicator)

nulltona <- function(string){
  if(!is.na(string) & string=="NULL"){
    return(NA)
  }
  else{
    return(string)
  }
}

demographics$AgeMAT <- sapply(demographics$AgeMAT, nulltona)
famhist$BMI <- sapply(famhist$BMI, nulltona)
demographics$SmokingMAT <- sapply(demographics$SmokingMAT, nulltona)
demographics$SmokingMAT <- sapply(demographics$SmokingMAT, nulltona)
demographics$AlcoholMAT <- sapply(demographics$AlcoholMAT, nulltona)
fertility$FertilityTreatment<- sapply(fertility$FertilityTreatment, nulltona)

demographics$SmokingMAT <- as.numeric(demographics$SmokingMAT)
demographics$AlcoholMAT <- as.numeric(demographics$AlcoholMAT)
demographics$SmokingMale <- as.numeric(demographics$SmokingMale)

demographics$BMIMale <- as.numeric(demographics$BMIMale)
demographics$AgePAT <- as.numeric(demographics$AgePAT)

fertility$FertilityTreatment <- as.numeric(fertility$FertilityTreatment)
fertility$PeriodRegular <- as.numeric(fertility$PeriodRegular)
fertility$PeriodIrregular <- 1 - fertility$PeriodRegular

nulltozero <- function(string){
  if(!is.na(string) & string=="NULL"){
    return(0)
  }
  else{
    return(string)
  }
}

demographics$VapingDayMAT <- sapply(demographics$VapingDayMAT, nulltozero)
demographics$VapingDayMAT <- as.numeric(demographics$VapingDayMAT)
demographics$doesvape <- as.numeric(demographics$VapingDayMAT>0)

for(i in 1:length(demographics$AlcoholUnitsMAT)){
  if(!is.na(demographics$AlcoholMAT[i]) & demographics$AlcoholMAT[i]==0){
    demographics$AlcoholUnitsMAT[i] <- 0 
  }
  if(is.na(demographics$AlcoholMAT[i])){
    demographics$AlcoholUnitsMAT[i] <- NA
  }
}

for(i in 1:length(demographics$CigarettesDayMAT)){
  if(!is.na(demographics$SmokingMAT[i]) & demographics$SmokingMAT[i]==0){
    demographics$CigarettesDayMAT[i] <- 0 
  }
  if(is.na(demographics$SmokingMAT[i])){
    demographics$CigarettesDayMAT[i] <- NA
  }
}

for(i in 1:length(demographics$CigarettesWeekMAT)){
  if(!is.na(demographics$SmokingMAT[i]) & demographics$SmokingMAT[i]==0){
    demographics$CigarettesWeekMAT[i] <- 0 
  }
  if(is.na(demographics$SmokingMAT[i])){
    demographics$CigarettesWeekMAT[i] <- NA
  }
}

fertility$FamilyMiscarriage <- as.numeric(famhist$FamilyMiscarriage)

demographics$AlcoholUnitsMAT <- sapply(demographics$AlcoholUnitsMAT, nulltona)
demographics$AlcoholUnitsMAT <- as.numeric(demographics$AlcoholUnitsMAT)
demographics$CigarettesDayMAT <- as.numeric(demographics$CigarettesDayMAT)
demographics$CigarettesWeekMAT <- as.numeric(demographics$CigarettesWeekMAT)


# ======= Create tiers =======

demographics$toomuchdrink <- rep(0, length(demographics$IDMAT))

for(i in 1:length(demographics$AlcoholMAT)){
  if(is.na(demographics$AlcoholUnitsMAT[i])){
    demographics$toomuchdrink[i] <-NA
  }
  else if(demographics$AlcoholUnitsMAT[i]>14){
    demographics$toomuchdrink[i] <- 1
  }
}

famhist$isoverweight <- rep(0, length(famhist$IDMAT))
famhist$isobese <- rep(0, length(famhist$IDMAT))

for(i in 1:length(famhist$BMI)){
  if(is.na(famhist$BMI[i])){
    famhist$isoverweight[i] <- NA
    famhist$isobese[i] <- NA
  }
  else if(famhist$BMI[i]>=25 & famhist$BMI[i]<=27.5){
    famhist$isoverweight[i] <- 1
  }
  else if(famhist$BMI[i]>27.5){
    famhist$isobese[i] <- 1
  }
}

demographics$heavysmoker <- rep(0, length(famhist$IDMAT))
demographics$lightsmoker <- rep(0, length(famhist$IDMAT))

for(i in 1:length(demographics$SmokingMAT)){
  if(is.na(demographics$CigarettesDayMAT[i])){
    demographics$heavysmoker[i] <- NA
    demographics$lightsmoker[i] <- NA
  }
  
  else if(demographics$CigarettesDayMAT[i]<7 & demographics$CigarettesDayMAT[i]>0){
    demographics$lightsmoker[i] <- 1
  }
  else if(demographics$CigarettesDayMAT[i]>=7){
    demographics$heavysmoker[i] <- 1
  }
}

demographics$middleage <- rep(0, length(demographics$IDMAT))
demographics$oldage <- rep(0, length(demographics$IDMAT))
demographics$agetier <- rep(0, length(demographics$IDMAT))

for(i in 1:length(demographics$AgeMAT)){
  if(is.na(demographics$AgeMAT[i])){
    demographics$middleage[i] <- NA
    demographics$oldage[i] <- NA
    demographics$agetier[i] <- NA
  }
  
  else if(demographics$AgeMAT[i]<=38 & demographics$AgeMAT[i]>34){
    demographics$middleage[i] <- 1
    demographics$agetier[i] <- 1
  }
  else if(demographics$AgeMAT[i]>40){
    demographics$oldage[i] <- 1
    demographics$agetier[i] <- 2
  }
}

demographics$agetier <- factor(demographics$agetier, levels=c(0,1,2), ordered=TRUE)


# ======= Create Natural Indicator ======= 

demographics$isnatural <- rep(0, length(demographics$IDMAT))
for(i in 1:length(firstpregnancy$methodOfConception)){
  id <- firstpregnancy$IDMAT[i]
  if(is.na(firstpregnancy$methodOfConception[i])){
    demographics$isnatural[demographics$IDMAT==id] <- NA
  }
  else if(firstpregnancy$methodOfConception[i]==1){
    demographics$isnatural[demographics$IDMAT==id] <- 1
  }
  else if(firstpregnancy$methodOfConception[i] %in% c(2,3,4,5,6)){
    demographics$isnatural[demographics$IDMAT==id] <- 0
  }
  
}


# ======= Create drugs lists =======

CannabisIDs <- c()
for(i in 1:length(drugs$IDMAT)){
  if(drugs$questionName[i]=="DrugsPlaceHolder" & drugs$type[i] %in% c("marijuana", "Marijuana", "MARIJUANA", "cannabis", "cannabis ", "Cannabis", "Cannabis ", "CANNABIS", "Canabis", "canabis", "CANABIS", "canibis", "CANIBIS", "cannibis", "CANNIBIS", "Cannibis", "marajuana", "Marajuana", "MARAJUANA", "Weed (Marijuna)", "canabis(vaporized)", "weed"))
    CannabisIDs <- c(CannabisIDs, drugs$IDMAT[i]) 
}

demographics$Cannabis <- as.numeric(demographics$IDMAT %in% CannabisIDs)


# ======= Create Adjusted BMI Column =======

famhist$BetterBMI <- famhist$Weight / (famhist$Height)^3


# ======= Create days since first consultation column =======

demographics$dayssincefirstconsultation <- as.numeric(as.Date("01/09/2020", "%d/%m/%Y") - demographics$DateOfFirstConsultationMAT)
firstpregnancy$dayssincefirstconsultation <- as.numeric(as.Date("01/09/2020", "%d/%m/%Y") - firstpregnancy$DateOfFirstConsultation)


# ======= Calculate the days till viable pregnancy and censoring information =======

firstpregnancy$timetofirstpregnancy <- as.numeric(firstpregnancy$ConceptionDate - firstpregnancy$DateOfFirstConsultation)

l <- length(demographics$IDMAT)

demographics$gotpreg <- rep(0, l)
demographics$daystillfirstpregnancy <- demographics$dayssincefirstconsultation
for(i in 1:l){
  if(demographics$IDMAT[i] %in% firstpregnancy$IDMAT){
    demographics$gotpreg[i] <- 1
    demographics$daystillfirstpregnancy[i] <- firstpregnancy$timetofirstpregnancy[firstpregnancy$IDMAT==demographics$IDMAT[i]]
  }
}

demographics$daystillviablepregnancy <- as.numeric(as.Date("01/09/2020", "%d/%m/%Y") - demographics$DateOfFirstConsultationMAT)
demographics$gotviablepreg <- rep(0, l) # 1 = have got viably pregnant
for(i in 1:l){
  # for those which have got pregnant since their first consultation, set daystillpregnancy as the time till 24 weeks into their first viable pregnancy
  if(demographics$IDMAT[i] %in% postpregnancies$IDMAT){
    subframe <- postpregnancies[postpregnancies$IDMAT==demographics$IDMAT[i],]
    k <- length(subframe$IDMAT)
    for(j in 1:k){
      # check if pregnancy was successful or has made it past 24 weeks
      if( (subframe$outcome[j]==1) | ( subframe$outcome[j]==10 & subframe$gestationWeeks[j]>=24 ) ){
        demographics$gotviablepreg[i] <- 1
        if(!is.na(subframe$ConceptionDate[j])){
          demographics$daystillviablepregnancy[i] <- as.numeric(subframe$ConceptionDate[j] + 24*7 - demographics$DateOfFirstConsultationMAT[i])
        }
        else{
          demographics$daystillviablepregnancy[i] <- as.numeric(subframe$expectedDueDate[j] - 12*7 - demographics$DateOfFirstConsultationMAT[i])
        }
        break
      }
    }
  }
}


# ======= Create Pregnancy History Totals =======

demographics$PostClinicTotalPregnancies <- rep(0, l)
demographics$PreClinicOutcome1 <- rep(0, l)
demographics$PreClinicOutcome2 <- rep(0, l)
demographics$PreClinicOutcome3 <- rep(0, l)
demographics$PreClinicOutcome4 <- rep(0, l)
demographics$PreClinicOutcome5 <- rep(0, l)
demographics$PreClinicOutcome6 <- rep(0, l)
demographics$PreClinicOutcome7 <- rep(0, l)
demographics$PreClinicOutcome8 <- rep(0, l)
demographics$PreClinicOutcome9 <- rep(0, l)
demographics$PreClinicOutcome10 <- rep(0, l)

for(i in 1:length(prepregnancies$IDMAT)){
  id <- prepregnancies$IDMAT[i]
  
  # add 1 to the total number of pre-clinic pregnancies for mother IDMAT
  demographics$PreClinicTotalPregnancies[demographics$IDMAT==id] <- demographics$PreClinicTotalPregnancies[demographics$IDMAT==id] + 1
  
  # add 1 to the relevant nmber of pre-clinic type pf pregnancies for mother IDMAT
  if(prepregnancies$outcome[i]==1){
    demographics$PreClinicOutcome1[demographics$IDMAT==id] <- demographics$PreClinicOutcome1[demographics$IDMAT==id] + 1
  }
  if(prepregnancies$outcome[i]==2){
    demographics$PreClinicOutcome2[demographics$IDMAT==id] <- demographics$PreClinicOutcome2[demographics$IDMAT==id] + 1
  }
  if(prepregnancies$outcome[i]==3){
    demographics$PreClinicOutcome3[demographics$IDMAT==id] <- demographics$PreClinicOutcome3[demographics$IDMAT==id] + 1
  }
  if(prepregnancies$outcome[i]==4){
    demographics$PreClinicOutcome4[demographics$IDMAT==id] <- demographics$PreClinicOutcome4[demographics$IDMAT==id] + 1
  }
  if(prepregnancies$outcome[i]==5){
    demographics$PreClinicOutcome5[demographics$IDMAT==id] <- demographics$PreClinicOutcome5[demographics$IDMAT==id] + 1
  }
  if(prepregnancies$outcome[i]==6){
    demographics$PreClinicOutcome6[demographics$IDMAT==id] <- demographics$PreClinicOutcome6[demographics$IDMAT==id] + 1
  }
  if(prepregnancies$outcome[i]==7){
    demographics$PreClinicOutcome7[demographics$IDMAT==id] <- demographics$PreClinicOutcome7[demographics$IDMAT==id] + 1
  }
  if(prepregnancies$outcome[i]==8){
    demographics$PreClinicOutcome8[demographics$IDMAT==id] <- demographics$PreClinicOutcome8[demographics$IDMAT==id] + 1
  }
  if(prepregnancies$outcome[i]==9){
    demographics$PreClinicOutcome9[demographics$IDMAT==id] <- demographics$PreClinicOutcome9[demographics$IDMAT==id] + 1
  }
  if(prepregnancies$outcome[i]==10){
    demographics$PreClinicOutcome10[demographics$IDMAT==id] <- demographics$PreClinicOutcome10[demographics$IDMAT==id] + 1
  }
}

demographics$PreClinicMiscarriages <- demographics$PreClinicOutcome2 + demographics$PreClinicOutcome3 + demographics$PreClinicOutcome4 + demographics$PreClinicOutcome5


# ======= Create Post Clinic Pregnancy Outcome Totals =======

demographics$PostClinicTotalPregnancies <- rep(0, l)
demographics$PostClinicOutcome1 <- rep(0, l)
demographics$PostClinicOutcome2 <- rep(0, l)
demographics$PostClinicOutcome3 <- rep(0, l)
demographics$PostClinicOutcome4 <- rep(0, l)
demographics$PostClinicOutcome5 <- rep(0, l)
demographics$PostClinicOutcome6 <- rep(0, l)
demographics$PostClinicOutcome7 <- rep(0, l)
demographics$PostClinicOutcome8 <- rep(0, l)
demographics$PostClinicOutcome9 <- rep(0, l)
demographics$PostClinicOutcome10 <- rep(0, l)

for(i in 1:length(postpregnancies$IDMAT)){
  id <- postpregnancies$IDMAT[i]
  
  # add 1 to the total number of pre-clinic pregnancies for mother IDMAT
  demographics$PostClinicTotalPregnancies[demographics$IDMAT==id] <- demographics$PostClinicTotalPregnancies[demographics$IDMAT==id] + 1
  
  # add 1 to the relevant nmber of pre-clinic type pf pregnancies for mother IDMAT
  if(postpregnancies$outcome[i]==1){
    demographics$PostClinicOutcome1[demographics$IDMAT==id] <- demographics$PostClinicOutcome1[demographics$IDMAT==id] + 1
  }
  if(postpregnancies$outcome[i]==2){
    demographics$PostClinicOutcome2[demographics$IDMAT==id] <- demographics$PostClinicOutcome2[demographics$IDMAT==id] + 1
  }
  if(postpregnancies$outcome[i]==3){
    demographics$PostClinicOutcome3[demographics$IDMAT==id] <- demographics$PostClinicOutcome3[demographics$IDMAT==id] + 1
  }
  if(postpregnancies$outcome[i]==4){
    demographics$PostClinicOutcome4[demographics$IDMAT==id] <- demographics$PostClinicOutcome4[demographics$IDMAT==id] + 1
  }
  if(postpregnancies$outcome[i]==5){
    demographics$PostClinicOutcome5[demographics$IDMAT==id] <- demographics$PostClinicOutcome5[demographics$IDMAT==id] + 1
  }
  if(postpregnancies$outcome[i]==6){
    demographics$PostClinicOutcome6[demographics$IDMAT==id] <- demographics$PostClinicOutcome6[demographics$IDMAT==id] + 1
  }
  if(postpregnancies$outcome[i]==7){
    demographics$PostClinicOutcome7[demographics$IDMAT==id] <- demographics$PostClinicOutcome7[demographics$IDMAT==id] + 1
  }
  if(postpregnancies$outcome[i]==8){
    demographics$PostClinicOutcome8[demographics$IDMAT==id] <- demographics$PostClinicOutcome8[demographics$IDMAT==id] + 1
  }
  if(postpregnancies$outcome[i]==9){
    demographics$PostClinicOutcome9[demographics$IDMAT==id] <- demographics$PostClinicOutcome9[demographics$IDMAT==id] + 1
  }
  if(postpregnancies$outcome[i]==10){
    demographics$PostClinicOutcome10[demographics$IDMAT==id] <- demographics$PostClinicOutcome10[demographics$IDMAT==id] + 1
  }
}

demographics$PostClinicMiscarriages <- demographics$PostClinicOutcome2 + demographics$PostClinicOutcome3 + demographics$PostClinicOutcome4 + demographics$PostClinicOutcome5











