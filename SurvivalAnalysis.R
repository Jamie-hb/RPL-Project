# -------------------
# Script to fit a Cox Regression to the data
# -------------------

library("survival")

# create viable ID set (i.e remove those without data)
# also remember to only include those with postive time till pregnancy! 
# i.e exclude those who were pregnant at the time of the first consultation
ViableIDs <- !is.na(data$AgeAtFirstConsultationMAT) & !is.na(data$BMIMAT) & !is.na(data$IsSmoker) & !is.na(data$Conception) & data$TimePositive & !is.na(data$CoffeeWeekMAT) & !is.na(data$PreClinicMiscarriages) 

# compute sample size 
SampleSize <- sum(ViableIDs)

# create explanatory, response and censoring variables
Age <- data$AgeAtFirstConsultationMAT[ViableIDs]
BMI <- data$BMIMAT[ViableIDs]
Smoking <- data$IsSmoker[ViableIDs]
Coffee <- data$CoffeeWeekMAT[ViableIDs]
Tea <- data$TeaWeekMAT[ViableIDs]
Caffeine <- Coffee 
PCM <- data$PreClinicMiscarriages[ViableIDs]
Time <- data$Time.to.conception.no.event[ViableIDs]
Censoring <- data$Conception[ViableIDs]

# run the Cox regression
PHFit <- coxph(Surv(time=Time, event=Censoring, type="right") ~ Age + BMI + Smoking + Caffeine)
print(PHFit)