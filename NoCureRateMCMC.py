""" Fits data from a simple cure rate model
        T ~ Exp( exp(beta^TX) ) 
    With right censoring using a MH sampler for 
    inference of beta """

import pandas as pd
import DataSetup
from DataSetup import *
import numpy as np
import matplotlib.pyplot as plt
import csv

# ======= Simulate the data =======
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

# Assemble data into a matrix
ones = np.repeat(1,t.size)
Covariates = np.stack((ones, Age, BMI, Smoking, Caffeine))


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

def MCMC(n, betazero, sigma):
    """ runs MCMC sampler """
    
    accept = 0
    reject = 0
    
    theta = np.zeros((5, n+1))
    theta[:,0] = betazero
    
    for i in range(n):
        # propose a new value
        thetaprop = np.random.multivariate_normal(mean=theta[:,i], cov=sigma)
        
        # compute log acceptance probability
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
    
    acceptrate = accept/n    
    
    # plot trace plots and histograms for each parameter
    plt.figure()
    plt.plot(theta[0,], zorder=1)
    plt.title(r'$\beta^0$, acceptance rate = ' + str(acceptrate))
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\beta^0_t$')
        
    plt.figure()
    plt.hist(theta[0,], bins=50, zorder=1)
    plt.title(r'$\beta^0$')
    plt.xlabel(r'$\beta^0$')
    plt.ylabel('frequency')
        
    plt.figure()
    plt.subplot(2,2,1)
    plt.plot(theta[1,], zorder=1)
    plt.title(r'$\beta^0$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\beta^0_t$')
    
    plt.subplot(2,2,2)
    plt.plot(theta[2,], zorder=1)
    plt.title(r'$\beta^1$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\beta^1_t$')
    
    plt.subplot(2,2,3)
    plt.plot(theta[3,], zorder=1)
    plt.title(r'$\beta^2$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\beta^2_t$')
    
    plt.subplot(2,2,4)
    plt.plot(theta[4,], zorder=1)
    plt.title(r'$\beta^3$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\beta^3_t$')
    plt.tight_layout()   
 
    plt.figure()
    plt.subplot(2,2,1)
    plt.hist(theta[1,], bins=50, zorder=1)
    plt.title(r'$\beta^0$')
    plt.xlabel(r'$\beta^0$')
    plt.ylabel('frequency')
        
    plt.subplot(2,2,2)
    plt.hist(theta[2,], bins=50, zorder=1)
    plt.title(r'$\beta^1$')
    plt.xlabel(r'$\beta^1$')
    plt.ylabel('frequency')
            
    plt.subplot(2,2,3)
    plt.hist(theta[3,], bins=50, zorder=1)
    plt.title(r'$\beta^2$')
    plt.xlabel(r'$\beta^2$')
    plt.ylabel('frequency')
            
    plt.subplot(2,2,4)
    plt.hist(theta[4,], bins=50, zorder=1)
    plt.title(r'$\beta^3$')
    plt.xlabel(r'$\beta^3$')
    plt.ylabel('frequency')
    plt.tight_layout()  
    
    return theta

# sigma = np.diag([0.001,0.001,0.001,0.001,0.001])

# betazero = np.array([-0.1,-0.1,-0.1,-0.1])
# qzero = 0.5       

sigma = np.diag([0.001, 0.001, 0.001, 0.001, 0.001])

betazero = np.array([-4.49523453, -0.03295658, -0.06429827, -1.06514478, -0.1])

# ===== Script to store chains in Excel for convergence analysis =====

outputFile = open('TimeCovariatesData.csv','w',newline='')
outputWriter = csv.writer(outputFile)

output = MCMC(1000, betazero, sigma)

sigmanew = 5.6644*np.cov(output)/5
betazeronew = output[:,-1]
outputnew = MCMC(1000, betazeronew, sigmanew)

sigmanew2 = 5.6644*np.cov(outputnew)/5
betazeronew2 = outputnew[:,-1]
outputnew2 = MCMC(1000, betazeronew2, sigmanew2)
        
sigmanew3 = 5.6644*np.cov(outputnew2)/5
betazeronew3 = outputnew2[:,-1]
outputnew3 = MCMC(10000, betazeronew3, sigmanew3)  

for k in range(5):
        outputWriter.writerow(outputnew3[k,])
        
outputFile.close()

# =============================================================================
# outputFile = open('NoCureRateMCMC.csv','w',newline='')
# outputWriter = csv.writer(outputFile)
# 
# for i in range(10):
#     betazero = np.random.normal(0, [5,1,1,5], 4)
#     qzero = np.random.uniform(0,1)        
# 
#     output = MCMC(10000, qzero, betazero, sigma)
#     
#     for k in range(5):
#         outputWriter.writerow(output[k,])
#         
# outputFile.close()
# =============================================================================
 
 