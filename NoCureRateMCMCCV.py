""" Fits data from a simple cure rate model
        T ~ Exp( exp(beta^TX) ) 
    With right censoring using a MH sampler for 
    inference of beta"""

import pandas as pd
import DataSetup
from DataSetup import *
import numpy as np
import matplotlib.pyplot as plt
import csv
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
#BMI = np.maximum(0, BMI-35)
Smoking = data['IsSmoker'].copy()
Smoking = np.asarray(Smoking)
Caffeine = data['CoffeeWeekMAT'].copy()
Caffeine = np.asarray(Caffeine)

# remove na entries:
# first find where the na entries are
workingset = np.logical_not( np.logical_or.reduce((np.isnan(t), np.isnan(delta), np.isnan(Age), np.isnan(BMI), np.isnan(Smoking), np.isnan(Caffeine))) )
workingset = np.logical_and(workingset, t>0)
# then remove them
t = t[workingset]
delta = delta[workingset]
Age = Age[workingset]
BMI = BMI[workingset]
Smoking = Smoking[workingset]
Caffeine = Caffeine[workingset]

SampleSize = t.size

# Assemble data into a matrix
ones = np.repeat(1,t.size)
covariates = np.stack((ones, Age, BMI, Smoking, Caffeine))


# ======= Likelihood, Prior =======

def loglikelihood(beta, t, delta, Covariates):
    """ computes log-likelihood """
    
    ll = 0
    
    ll += np.sum(delta*np.dot(beta, Covariates))
    ll -= np.sum(np.exp(np.dot(beta, Covariates))*t)
    
    return ll

def logprior(beta):
    """ computes log prior likelihood """
    
    sigmas = np.array([10,10,10,10,10])
    
    lp = -np.sum((beta**2) / (2*sigmas))
    
    return lp

# ======= MCMC =======

def MCMC(n, betazero, sigma, t, delta, Covariates):
    """ runs MCMC sampler """
    
    accept = 0
    reject = 0
    
    theta = np.zeros((5, n+1))
    theta[:,0] = betazero
    
    for i in range(n):
        # propose a new value
        thetaprop = np.random.multivariate_normal(mean=theta[:,i], cov=sigma)
        
        # compute the log acceptance probability
        lap = loglikelihood(thetaprop, t, delta, Covariates)
        lap -= loglikelihood(theta[:,i], t, delta, Covariates)
        lap += logprior(thetaprop) - logprior(theta[:,i])
        
        # accept-reject step
        u = np.random.uniform(0,1)
        
        if np.log(u)<lap:
            theta[:,i+1] = thetaprop
            accept += 1
        else:
            theta[:,i+1] = theta[:,i]
            reject += 1
    
    return theta

# sigma = np.diag([0.001,0.001,0.001,0.001,0.001])

# betazero = np.array([-0.1,-0.1,-0.1,-0.1])
# qzero = 0.5       

sigma = np.diag([0.001, 0.001, 0.001, 0.001, 0.001])

betazero = np.array([-4.49523453, -0.03295658, -0.06429827, -1.06514478, -0.1])

# =============================================================================
# outputFile = open('TimeCovariatesData.csv','w',newline='')
# outputWriter = csv.writer(outputFile)
# =============================================================================


def Fit(betazero, sigma, t, delta, Covariates):
    """ Run the MCMC multiple times and return
        mean estimates of the last run """
        
    output = MCMC(10000, betazero, sigma, t, delta, Covariates)
 
    sigmanew = 5.6644*np.cov(output)/5
    betazeronew = output[:,-1]
    outputnew = MCMC(1000, betazeronew, sigmanew, t, delta, Covariates)

    sigmanew2 = 5.6644*np.cov(outputnew)/5
    betazeronew2 = outputnew[:,-1]
    outputnew2 = MCMC(1000, betazeronew2, sigmanew2, t, delta, Covariates)
        
    sigmanew3 = 5.6644*np.cov(outputnew2)/5
    betazeronew3 = outputnew2[:,-1]
    outputnew3 = MCMC(10000, betazeronew3, sigmanew3, t, delta, Covariates)   
    
    # compute row-wise means of the final output
    res = np.mean(outputnew3, axis=1)
    
    return res

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
    

def ConceptionProb(beta, x, t):
    """ Computes probability that a woman with covariates x
        will achieve a pregnancy within t days.
        Uses fitted beta parameters """
    
    return 1 - np.exp(-np.exp(np.dot(beta, x))*t)

# Make predictions
for j in range(9):
    # run the MCMC
    Coeffs = Fit(betazero, sigma, np.delete(t,IndexMatrix[j,:]), np.delete(delta, IndexMatrix[j,:]), np.delete(covariates, IndexMatrix[j,:], axis=1))               
    print(Coeffs)
    
    # make predictions (note: the leading 1* is to turn True,False into 1,0)
    for i in IndexMatrix[j,:]:
        Prediction[0,i] = 1*(ConceptionProb(Coeffs, covariates[:,i], Times[0]) > 0.5) 
        Prediction[1,i] = 1*(ConceptionProb(Coeffs, covariates[:,i], Times[1]) > 0.5) 
        Prediction[2,i] = 1*(ConceptionProb(Coeffs, covariates[:,i], Times[2]) > 0.5)
    
# Some special code for the final validation set 
Coeffs = Fit(betazero, sigma, np.delete(t,FinalValidationSet), np.delete(delta,FinalValidationSet), np.delete(covariates,FinalValidationSet,axis=1))               
print(Coeffs)    

# make predictions (note: the leading 1* is to turn True,False into 1,0)
for i in FinalValidationSet:
    Prediction[0,i] = 1*(ConceptionProb(Coeffs, covariates[:,i], Times[0]) > 0.5)
    Prediction[1,i] = 1*(ConceptionProb(Coeffs, covariates[:,i], Times[1]) > 0.5)
    Prediction[2,i] = 1*(ConceptionProb(Coeffs, covariates[:,i], Times[2]) > 0.5)
        
# Test predictions
# First create the confusion matrices, in the format
#               Actual 0    Actual 1
# Predicted 0     TN          FN  
# Predicted 1     FP          TP
ConfusionMatrixOne = np.zeros((2,2))
ConfusionMatrixTwo = np.zeros((2,2))
ConfusionMatrixThree = np.zeros((2,2))

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

