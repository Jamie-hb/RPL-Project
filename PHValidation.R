# ------------------------
# Validation for the PH model, done via LOOCV, using the fitted model to classify women into
# risk groups then plotting the KM curves for these groups
# ------------------------

# === Create viable ID set (i.e remove those without data) ===

ViableIDs <- !is.na(data$AgeAtFirstConsultationMAT) & !is.na(data$BMIMAT) & !is.na(data$IsSmoker) & !is.na(data$Conception) & data$TimePositive & !is.na(data$CoffeeWeekMAT) 
# record sample size
SampleSize <- sum(ViableIDs)

# create explanatory and response variables
Age <- data$AgeAtFirstConsultationMAT[ViableIDs]
BMI <- data$BMIMAT[ViableIDs]
Smoking <- data$IsSmoker[ViableIDs]
Caffeine <- data$CoffeeWeekMAT[ViableIDs]
Time <- data$Time.to.conception.no.event[ViableIDs]
Censoring <- data$Conception[ViableIDs]


# === Determine risk groups ===

PredictedRate <- rep(0, SampleSize)

for(i in 1:SampleSize){
  # train the LR on the data set excluding the ith individual
  PHFitNew <- coxph(Surv(time=Time[-i], event=Censoring[-i], type="right") ~ Age[-i] + BMI[-i] + Smoking[-i] + Caffeine[-i])
  
  # get the newly fitted coefficients 
  BetaNew <- as.vector(PHFitNew$coefficients)
  
  # make a prediction for the ith data point
  CovariatesNew <- c(Age[i] , BMI[i], Smoking[i], Caffeine[i])
  BTX <- BetaNew %*% CovariatesNew
  PredictedRate[i] <- exp(BetaNew %*% CovariatesNew)
}

LowRiskThreshold <- quantile(PredictedRate, 0.33)
MedRiskThreshold <- quantile(PredictedRate, 0.67)

RiskGroup <- 1*(PredictedRate<=LowRiskThreshold) + 2*(PredictedRate<=MedRiskThreshold & PredictedRate>LowRiskThreshold) + 3*(PredictedRate>MedRiskThreshold)


# === Plot KM curves for the different risk groups ===

plot(survfit.formula(Surv(time=Time[RiskGroup==1], event=Censoring[RiskGroup==1], type="right") ~ 1), col="blue", xlab="Time (Days)", ylab="Estimated S(t)", main="KM Curves for Risk Groups Predicted by the PH Model")
lines(survfit.formula(Surv(time=Time[RiskGroup==2], event=Censoring[RiskGroup==2], type="right") ~ 1), col="black")
lines(survfit.formula(Surv(time=Time[RiskGroup==3], event=Censoring[RiskGroup==3], type="right") ~ 1), col="red")
legend("topright", legend=c("Low Risk", "Medium Risk", "High Risk"), col=c("blue", "black", "red"), lty=c(1,1,1))

