# --------------------
# Script to fit a logistic regression to the data set
# --------------------

# create viable ID set (i.e remove those without data)
ViableIDs <- data$Conception==1 & !is.na(data$AgeAtFirstConsultationMAT) & !is.na(data$BMIMAT) & !is.na(data$IsSmoker) & !is.na(data$PolycysticOvaries) & !is.na(data$PreClinicMiscarriages)

# record sample size
SampleSize <- sum(ViableIDs)

# create explanatory and response variables
Age <- data$AgeAtFirstConsultationMAT[ViableIDs]
Age <- pmax(0, Age - 35) # note: 35 < 40 from before
BMI <- data$BMIMAT[ViableIDs]
Smoking <- data$IsSmoker[ViableIDs] 
PCOS <- data$PolycysticOvaries[ViableIDs]
PreClinicMiscarriages <- data$PreClinicMiscarriages[ViableIDs]
Outcome <- data$FirstPregSuccess[ViableIDs]

# do the LR
LRFit <- glm(Outcome ~ Age + BMI + Smoking + PCOS + PreClinicMiscarriages, family=binomial)
print(summary(LRFit)) # run this line to see the results of the LR