""" 10-Fold Cross-Validation for the cure rate model
        T ~ q(α^TX)Exp(β) + (1 - q(α^TX))∞ 
    with q = 1/1+exp(-α^TX) """
    
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import DataSetup
from DataSetup import *
#import csv
import random

# ======= Load in the data =======
# observed times
t = data['Time to conception/no event'].copy()
t = np.asarray(t, dtype='float64')  

# observed censoring 
delta = data['Conception'].copy()
delta = np.asarray(delta) 

# covariates 
Age = data['AgeAtFirstConsultationMAT'].copy()
Age = np.asarray(Age)
BMI = data['BMIMAT'].copy()
BMI = np.asarray(BMI)
Smoking  = data['IsSmoker'].copy()
Smoking = np.asarray(Smoking)
Coffee = data['CoffeeWeekMAT'].copy()
Coffee = np.asarray(Coffee)

# remove na entries:
# first find where the na entries are
workingset = np.logical_not( np.logical_or.reduce((np.isnan(t), np.isnan(delta), np.isnan(Age), np.isnan(BMI), np.isnan(Smoking), np.isnan(Coffee) )) )
workingset = np.logical_and(workingset, t>0)
# then remove them
t = t[workingset]
delta = delta[workingset]
Age = Age[workingset]
BMI = BMI[workingset]
Smoking = Smoking[workingset]
Coffee = Coffee[workingset]

# Compute sample size and problem dimension
SampleSize = t.size
NumPara = 6 # number of parameters to be estimated 

# Assemble data into a matrix, covariates[i,j] = covariate i for individual j
ones = np.repeat(1,SampleSize)
covariates = np.stack((ones, Age, BMI, Smoking, Coffee))


#  ======= Likelihood, Prior =======

def loglikelihood(alpha, beta, t, delta, Covariates):
    """ computes log-likelihood """
    
    ll = 0
    
    ll += np.sum(delta)*np.log(beta)
    ll -= np.sum(t*delta)*beta
    ll += np.sum((1 - delta)*( np.log(np.exp(-np.dot(alpha, Covariates)) + np.exp(-beta*t)) ) )
    ll -= np.sum(np.log(1 + np.exp(-np.dot(alpha, Covariates)) ))
       
    return ll

def logprior(alpha, beta):
    """ computes log prior likelihood
        α ~ N(0,Σ), 10000β ~ Gamma(6,1) """
    
    sigma_prior = np.array([10,0.3,0.3,3,1])
    a = 1
    b = 1
    
    lp = -np.sum((alpha**2) / (2*sigma_prior))
    lp += (a - 1)*np.log(1000*beta) - b*1000*beta
    
    return lp


# ======= MCMC =======

def MCMC(n, alphazero, betazero, sigma, t, delta, Covariates):
    """ runs MCMC sampler """
    
    accept = 0
    reject = 0
    
    theta = np.zeros((NumPara, n+1))
    theta[0:(NumPara - 1), 0] = alphazero
    theta[NumPara - 1, 0] = betazero
    
    for i in range(n):
        # propose a new value
        thetaprop = np.random.multivariate_normal(mean=theta[:,i], cov=sigma)
        
        # check for β* > 0
        if thetaprop[NumPara - 1]>0: 
            # set current and proposal alpha and beta values       
            alphacurrent = theta[0:(NumPara - 1), i]
            betacurrent = theta[NumPara - 1, i]
        
            alphaprop = thetaprop[0:(NumPara - 1)]
            betaprop = thetaprop[NumPara - 1]
            
            # compute log acceptance probability
            lap = loglikelihood(alphaprop, betaprop, t, delta, Covariates)
            lap -= loglikelihood(alphacurrent, betacurrent, t, delta, Covariates)
            lap += logprior(alphaprop, betaprop) - logprior(alphacurrent, betacurrent)
            
            # accept-reject step
            u = np.random.uniform(0,1)
            
            if np.log(u)<lap:
                theta[:,i+1] = thetaprop
                accept += 1
            else:
                theta[:,i+1] = theta[:,i]
                reject += 1
        else:
            theta[:,i+1] = theta[:,i]
            reject += 1
    
    return theta


def Fit(alphazero, betazero, sigma, t, delta, Covariates):
    """ Run the MCMC multiple times and return
        mean estimates of the last run """
        
    output = MCMC(10000, alphazero, betazero, sigma, t, delta, Covariates)
 
    sigmanew = 5.6644*np.cov(output) / NumPara
    thetazeronew = output[:,-1]
    alphazeronew = thetazeronew[0:(NumPara-1)]
    betazeronew = thetazeronew[NumPara-1]
    outputnew = MCMC(1000, alphazeronew, betazeronew, sigmanew, t, delta, Covariates)

    sigmanew2 = 5.6644*np.cov(outputnew) / NumPara
    thetazeronew2 = outputnew[:,-1]
    alphazeronew2 = thetazeronew2[0:(NumPara-1)]
    betazeronew2 = thetazeronew2[NumPara-1]
    outputnew2 = MCMC(1000, alphazeronew2, betazeronew2, sigmanew2, t, delta, Covariates)
        
    sigmanew3 = 5.6644*np.cov(outputnew2) / NumPara
    thetazeronew3 = outputnew2[:,-1]
    alphazeronew3 = thetazeronew3[0:(NumPara-1)]
    betazeronew3 = thetazeronew3[NumPara-1]
    outputnew3 = MCMC(10000, alphazeronew3, betazeronew3, sigmanew3, t, delta, Covariates)  
    
    # compute row-wise means of the final output
    res = np.mean(outputnew3, axis=1)
    
    return res
    

# alphazero = np.random.multivariate_normal(mean=np.array([0,0,0,0,0]), cov=np.diag([10,0.3,0.3,3,1]))
# betazero = np.random.gamma(6,1)/1000 

alphazero = np.array([0.82883734, -0.22993219,  0.0520644 , -1.15676514, 0.1])
betazero = 0.005 

sigma = np.array([[ 1.12660536e+00, -2.20321162e-02, -1.23968162e-02,
        -1.75385503e-01,  1.25278468e-02,  4.67710126e-05],
       [-2.20321162e-02,  5.57721444e-04,  9.51173429e-05,
         4.45922678e-03, -4.17533174e-04, -1.35718453e-06],
       [-1.23968162e-02,  9.51173429e-05,  3.24053474e-04,
         5.64609635e-04, -5.43844430e-05, -2.76171158e-07],
       [-1.75385503e-01,  4.45922678e-03,  5.64609635e-04,
         1.48867857e-01, -8.48026828e-03, -1.79307812e-05],
       [ 1.25278468e-02, -4.17533174e-04, -5.43844430e-05,
        -8.48026828e-03,  5.09209940e-03,  3.07777407e-06],
       [ 4.67710126e-05, -1.35718453e-06, -2.76171158e-07,
        -1.79307812e-05,  3.07777407e-06,  5.16316968e-08]])

# === CV ===
""" 10-fold CV """

# Create an array to store the predictions
# entry [i,j] = prediction for conception of individual i
# within Times[0], Times[1] and Times[2] days
Prediction = np.zeros((3,SampleSize))
Times = np.array([365,730,1000])

# Create validation sets
K = 10 # number of validation sets
L = np.floor(SampleSize / K).astype(int) # size of validation set

# create a matrix with each row being the indices of a validation set
Indexes = np.arange(SampleSize) # = [0,1,2,...,SampleSize]
random.shuffle(Indexes) # randomly permutes Indexes
IndexMatrix = np.zeros((K, L)).astype(int)
for i in range(10):
    IndexMatrix[i,:] = Indexes[(i*L):((i+1)*L)]
    
# the indices not in any validation set, these shall be added to the final validation set    
RemainderIndices = Indexes[(K*L):] 
FinalValidationSet = np.concatenate((IndexMatrix[9,:], RemainderIndices))
    

def ConceptionProb(alpha, beta, x, t):
    """ Computes probability that a woman with covariates x
        will achieve a pregnancy within t days.
        Uses fitted alpha and beta parameters """
        
    ECProb = 1 / (1 + np.exp(-np.dot(alpha, x))) # P(in Event Class) = 1/(1 + e^(-α.X))
    ExpCDF = 1 - np.exp(-beta*t) # P(T <= t | in Event Class) = 1 - e^(-βt)
    
    return ECProb * ExpCDF

# Make predictions
for j in range(9):
    # run the MCMC
    Coeffs = Fit(alphazero, betazero, sigma, np.delete(t,IndexMatrix[j,:]), np.delete(delta, IndexMatrix[j,:]), np.delete(covariates, IndexMatrix[j,:], axis=1))               
    print(Coeffs)
    
    # make predictions (note: the leading 1* is to turn True,False into 1,0)
    for i in IndexMatrix[j,:]:
        Prediction[0,i] = ConceptionProb(Coeffs[:5], Coeffs[5], covariates[:,i], Times[0]) 
        Prediction[1,i] = ConceptionProb(Coeffs[:5], Coeffs[5], covariates[:,i], Times[1]) 
        Prediction[2,i] = ConceptionProb(Coeffs[:5], Coeffs[5], covariates[:,i], Times[2]) 
    
# Some special code for the final validation set 
Coeffs = Fit(alphazero, betazero, sigma, np.delete(t,FinalValidationSet), np.delete(delta,FinalValidationSet), np.delete(covariates,FinalValidationSet,axis=1))               
print(Coeffs)    

# make predictions (note: the leading 1* is to turn True,False into 1,0)
for i in FinalValidationSet:
    Prediction[0,i] = ConceptionProb(Coeffs[:5], Coeffs[5], covariates[:,i], Times[0]) 
    Prediction[1,i] = ConceptionProb(Coeffs[:5], Coeffs[5], covariates[:,i], Times[1])
    Prediction[2,i] = ConceptionProb(Coeffs[:5], Coeffs[5], covariates[:,i], Times[2]) 
        
# ===== Test the Predictions =====

# First create the confusion matrices, in the format
#               Actual 0    Actual 1
# Predicted 0     TN          FN  
# Predicted 1     FP          TP
ConfusionMatrixOne = np.zeros((2,2))
ConfusionMatrixTwo = np.zeros((2,2))
ConfusionMatrixThree = np.zeros((2,2))

# GotPreg = np.zeros((3,SampleSize))
# GotPreg[0,:] = 1*(t <= Times[0])
# GotPreg[1,:] = 1*(t <= Times[1])
# GotPreg[2,:] = 1*(t <= Times[2])

# If δ = 1 we know T = t, if δ = 0 we know T >= t
# So if δ = 0 and T <= x, we do not know whether a 
# pregnancy was achieved within x days. I don't think
# this should lead to any bias in predictor evaluation, as this
# should affect both affirmative and negative predictions
# equally

for i in range(SampleSize):
    # = Test Times[0] days =
    # Predict P
    if Prediction[0,i]==1:
        # if t > Times[0], then FP (regrdless of δ)
        if t[i]>Times[0]:
            ConfusionMatrixOne[1,0] += 1
        # if t <= Times[0] and δ = 1, then TP
        if (t[i]<=Times[0]) & (delta[i]==1):
            ConfusionMatrixOne[1,1] += 1
        # if t <= Times[0] and δ = 0, correctness of prediction unknown
    
    # Predict N
    if Prediction[0,i]==0:
        # if t > Times[0], then TN (regardless of δ)
        if t[i]>Times[0]:
            ConfusionMatrixOne[0,0] += 1
        # if t <= Times[0] and δ = 1, then FN
        if (t[i]<=Times[0]) & (delta[i]==1):
            ConfusionMatrixOne[0,1] += 1
        # if t <= Times[0] and δ = 0, correctness of prediction unknown
        
            
    # = Test Times[1] days =
    # Predict P
    if Prediction[1,i]==1:
        # if t > Times[1], then FP (regrdless of δ)
        if t[i]>Times[1]:
            ConfusionMatrixTwo[1,0] += 1
        # if t <= Times[1] and δ = 1, then TP
        if (t[i]<=Times[1]) & (delta[i]==1):
            ConfusionMatrixTwo[1,1] += 1
        # if t <= Times[1] and δ = 0, correctness of prediction unknown
    
    # Predict N
    if Prediction[1,i]==0:
        # if t > Times[1], then TN (regardless of δ)
        if t[i]>Times[1]:
            ConfusionMatrixTwo[0,0] += 1
        # if t <= Times[1] and δ = 1, then FN
        if (t[i]<=Times[1]) & (delta[i]==1):
            ConfusionMatrixTwo[0,1] += 1
        # if t <= Times[1] and δ = 0, correctness of prediction unknown
        
            
    # = Test Times[2] days =
    # Predict P
    if Prediction[2,i]==1:
        # if t > Times[2], then FP (regrdless of δ)
        if t[i]>Times[2]:
            ConfusionMatrixThree[1,0] += 1
        # if t <= Times[2] and δ = 1, then TP
        if (t[i]<=Times[2]) & (delta[i]==1):
            ConfusionMatrixThree[1,1] += 1
        # if t <= Times[2] and δ = 0, correctness of prediction unknown
    
    # Predict N
    if Prediction[2,i]==0:
        # if t > Times[2], then TN (regardless of δ)
        if t[i]>Times[2]:
            ConfusionMatrixThree[0,0] += 1
        # if t <= Times[2] and δ = 1, then FN
        if (t[i]<=Times[2]) & (delta[i]==1):
            ConfusionMatrixThree[0,1] += 1
        # if t <= Times[2] and δ = 0, correctness of prediction unknown
            
            
# Compute accuracy metrics
AccuracyMetrics = np.zeros((3,5))

# Accuracy = TP + TN / TP + TN + FP + FN (what percentage of predictions were correct?)
# NOTE: cannot use sample size in the denominator, as we were not able to
# test all predictions due to censoring
AccuracyMetrics[0,0] = (ConfusionMatrixOne[0,0] + ConfusionMatrixOne[1,1]) / np.sum(ConfusionMatrixOne)
AccuracyMetrics[1,0] = (ConfusionMatrixTwo[0,0] + ConfusionMatrixTwo[1,1]) / np.sum(ConfusionMatrixTwo)
AccuracyMetrics[2,0] = (ConfusionMatrixThree[0,0] + ConfusionMatrixThree[1,1]) / np.sum(ConfusionMatrixThree)

# Sensitivty = TP / TP + FN (how good is the classifier at picking up positive cases?)
AccuracyMetrics[0,1] = ConfusionMatrixOne[1,1] / np.sum(ConfusionMatrixOne[:,1])
AccuracyMetrics[1,1] = ConfusionMatrixTwo[1,1] / np.sum(ConfusionMatrixTwo[:,1])
AccuracyMetrics[2,1] = ConfusionMatrixThree[1,1] / np.sum(ConfusionMatrixThree[:,1])

# Specificity = TN / TN + FP (how good is the classifier at picking up negative cases?)
AccuracyMetrics[0,2] = ConfusionMatrixOne[0,0] / np.sum(ConfusionMatrixOne[:,0])
AccuracyMetrics[1,2] = ConfusionMatrixTwo[0,0] / np.sum(ConfusionMatrixTwo[:,0])
AccuracyMetrics[2,2] = ConfusionMatrixThree[0,0] / np.sum(ConfusionMatrixThree[:,0])

# Precision = TP / TP + FP (how correct were the positive predictions?)
AccuracyMetrics[0,3] = ConfusionMatrixOne[1,1] / np.sum(ConfusionMatrixOne[1,:])
AccuracyMetrics[1,3] = ConfusionMatrixTwo[1,1] / np.sum(ConfusionMatrixTwo[1,:])
AccuracyMetrics[2,3] = ConfusionMatrixThree[1,1] / np.sum(ConfusionMatrixThree[1,:])

# NPV = TN / TN + FN (how correct were the negative predictions?)
AccuracyMetrics[0,4] = ConfusionMatrixOne[0,0] / np.sum(ConfusionMatrixOne[0,:])
AccuracyMetrics[1,4] = ConfusionMatrixTwo[0,0] / np.sum(ConfusionMatrixTwo[0,:])
AccuracyMetrics[2,4] = ConfusionMatrixThree[0,0] / np.sum(ConfusionMatrixThree[0,:])


# ======= For probabilistic predictions ===========

# BreakPoints = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
# ClassSize = np.zeros((3, BreakPoints.size))
# EmpiricalCount = np.zeros((3, BreakPoints.size))

# for i in range(SampleSize):
#     for k in range(3):
        
#         if t[i] > Times[k]:
#             Class = np.where(Prediction[k,i]<BreakPoints)[0][0]
#             ClassSize[k, Class] += 1
#             # EmpiricalCount[k, Class] += 0
            
#         if (t[i]<=Times[k]) & (delta[i]==1):
#             Class = np.where(Prediction[k,i]<BreakPoints)[0][0]
#             ClassSize[k, Class] += 1
#             EmpiricalCount[k, Class] += 1
            
# EmpiricalProb = EmpiricalCount / ClassSize
                
    

    
    














