# --------------------
# Plots ROC and Precision-Recall curves for the LR model
# --------------------

# === Create viable ID set (i.e remove those without data) ===

ViableIDs <- data$Conception==1 & !is.na(data$FirstPregSuccess) & !is.na(data$IsSmoker) & !is.na(data$AgeAtFirstConsultationMAT) & !is.na(data$BMIMAT) & !is.na(data$PolycysticOvaries) & !is.na(data$PreClinicMiscarriages)
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


# === Define a function to compute the accuracy metrics for a given threshold 
#     (uses same script from LRCV.R) ===

LOOCV <- function(Threshold){
  ConfusionMatrixTemp <- matrix(0, 2, 2)
  # The format used for the confusion matrix here is
  #               Actual 0    Actual 1
  # Predicted 0     TN          FN  
  # Predicted 1     FP          TP
  
  for(i in 1:SampleSize){
    # train the LR on the data set excluding the ith individual
    LRFitNew <- glm(Outcome[-i] ~ Age[-i] + BMI[-i] + Smoking[-i] + PCOS[-i] + PreClinicMiscarriages[-i], family=binomial)
    
    # get the newly fitted coefficients 
    BetaNew <- as.vector(LRFitNew$coefficients)
    
    # make a prediction for the ith data point
    CovariatesNew <- c(1, Age[i] , BMI[i], Smoking[i], PCOS[i], PreClinicMiscarriages[i])
    BTX <- BetaNew %*% CovariatesNew
    PredictedProbs <- exp(BTX) / (1 + exp(BTX)) 
    Predictions <- as.numeric(PredictedProbs >= Threshold)
    
    # update the confusion matrix
    if(Predictions==0 & Outcome[i]==0){
      ConfusionMatrixTemp[1,1] <- ConfusionMatrixTemp[1,1] + 1
    }
    else if(Predictions==0 & Outcome[i]==1){
      ConfusionMatrixTemp[1,2] <- ConfusionMatrixTemp[1,2] + 1
    }
    else if(Predictions==1 & Outcome[i]==0){
      ConfusionMatrixTemp[2,1] <- ConfusionMatrixTemp[2,1] + 1
    }
    else if(Predictions==1 & Outcome[i]==1){
      ConfusionMatrixTemp[2,2] <- ConfusionMatrixTemp[2,2] + 1
    }
  }
  
  # compute some accuracy metrics
  TPR <- ConfusionMatrixTemp[2,2] / sum(ConfusionMatrixTemp[,2]) # TP / TP + FN
  FPR <- ConfusionMatrixTemp[2,1] / sum(ConfusionMatrixTemp[,1]) # FP / TN + FP
  
  Precision <- ConfusionMatrixTemp[2,2] / sum(ConfusionMatrixTemp[2,]) # TP / TP + FP
  Recall <- ConfusionMatrixTemp[2,2] / sum(ConfusionMatrixTemp[,2]) # TP / TP + FN
  
  return(c(FPR, TPR, Recall, Precision))
  
}


# === Plot the ROC curve ===

grain <- 0.1
mesh <- seq(0, 1, grain)
points <- matrix(0, 4, length(mesh))

for(i in 1:length(mesh)){
  points[,i] <- LOOCV(mesh[i])
}

par(mfrow=c(1,2))

plot(points[1,], points[2,], type="l", col="blue", main="ROC Curves", xlab="FPR", ylab="TPR")
abline(c(0,0),c(1,1),col="red")

plot(points[3,], points[4,], type="l", col="blue", main="Precision-Recall Curve", xlab="Recall", ylab="Precision")


