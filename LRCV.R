# ----------------------
# Runs leave one out cross-validation on the training data for the LR model
# ----------------------

# === Create viable ID set (i.e remove those without data) ===

ViableIDs <- data$Conception==1 & !is.na(data$FirstPregSuccess) & !is.na(data$AgeAtFirstConsultationMAT) & !is.na(data$BMIMAT) & !is.na(data$PolycysticOvaries) & !is.na(data$PreClinicMiscarriages)
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

# set the threshold of predicted probability for predicting a successful pregnancy
Threshold <- 0.7


# === Build the confusion matrix ===

ConfusionMatrix <- matrix(0, 2, 2)
# The format used for the confusion matrix here is
#               Actual 0    Actual 1
# Predicted 0     TN          FN  
# Predicted 1     FP          TP

PredictedProbs <- rep(0, SampleSize)
Predictions <- rep(0, SampleSize)

for(i in 1:SampleSize){
  # train the LR on the data set excluding the ith individual
  LRFitNew <- glm(Outcome[-i] ~ Age[-i] + BMI[-i] + PCOS[-i] + PreClinicMiscarriages[-i], family=binomial)
  
  # get the newly fitted coefficients 
  BetaNew <- as.vector(LRFitNew$coefficients)
  
  # make a prediction for the ith data point
  CovariatesNew <- c(1, Age[i], BMI[i],  PCOS[i], PreClinicMiscarriages[i])
  BTX <- BetaNew %*% CovariatesNew
  PredictedProbs[i] <- exp(BTX) / (1 + exp(BTX)) 
  Predictions[i] <- as.numeric(PredictedProbs[i] >= Threshold)
  
  # update the confusion matrix
  if(Predictions[i]==0 & Outcome[i]==0){
    ConfusionMatrix[1,1] <- ConfusionMatrix[1,1] + 1
  }
  else if(Predictions[i]==0 & Outcome[i]==1){
    ConfusionMatrix[1,2] <- ConfusionMatrix[1,2] + 1
  }
  else if(Predictions[i]==1 & Outcome[i]==0){
    ConfusionMatrix[2,1] <- ConfusionMatrix[2,1] + 1
  }
  else if(Predictions[i]==1 & Outcome[i]==1){
    ConfusionMatrix[2,2] <- ConfusionMatrix[2,2] + 1
  }
}


# === Compute some accuracy metrics ===

Accuracy <- sum(diag(ConfusionMatrix)) / SampleSize # TP + TN / n
Sensitivity <- ConfusionMatrix[2,2] / sum(ConfusionMatrix[,2]) # TP / TP + FN
Specificity <- ConfusionMatrix[1,1] / sum(ConfusionMatrix[,1]) # TN / TN + FP
Precision <- ConfusionMatrix[2,2] / sum(ConfusionMatrix[2,]) # TP / TP + FP
NPV <- ConfusionMatrix[1,1] / sum(ConfusionMatrix[1,]) # TN / TN + FN

# Precision and NPV are probably the ones we're interested in, we want patients to have confidence
# in our predictions!


# === Print confusion matrix and accuracy metrics ===

AccuracyMetrics <- c(Threshold, Accuracy, Sensitivity, Specificity, Precision, NPV)
names(AccuracyMetrics) <- c("Threshold", "Accuracy", "Sensitivity", "Specificity", "Precision", "NPV")

rownames(ConfusionMatrix) <- c("Predicted 0", "Predicted 1")
colnames(ConfusionMatrix) <- c("Actual 0", "Actual 1")

print(ConfusionMatrix)
print(AccuracyMetrics)


