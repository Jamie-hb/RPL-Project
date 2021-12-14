""" Simulates data from a simple cure rate model
        T ~ qExp( exp(β^TX) ) + (1 - q)∞ 
    With right censoring. Then tests a MH sampler for 
    inference of q and β"""

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
Age = np.maximum(0, Age-35)
BMI = data['BMIMAT'].copy()
BMI = np.asarray(BMI)
BMI = np.maximum(0, BMI-35)
Smoking = data['IsSmoker'].copy()
Smoking = np.asarray(Smoking)

# remove na entries:
# first find where the na entries are
workingset = np.logical_not( np.logical_or.reduce((np.isnan(t), np.isnan(delta), np.isnan(Age), np.isnan(BMI), np.isnan(Smoking))) )
workingset = np.logical_and(workingset, t>0)
# then remove them
t = t[workingset]
delta = delta[workingset]
Age = Age[workingset]
BMI = BMI[workingset]
Smoking = Smoking[workingset]

# ======= Likelihood, Prior =======

def loglikelihood(q, beta, t, delta, Age, BMI, Smoking):
    """ computes log-likelihood """
    
    ll = 0
    
    ll += np.sum(delta)*np.log(q)
    ll += np.sum(delta*(beta[0] + beta[1]*Age + beta[2]*BMI + beta[3]*Smoking))
    ll -= np.sum(delta*np.exp(beta[0] + beta[1]*Age + beta[2]*BMI + beta[3]*Smoking)*t)
    ll += np.sum((1 - delta)*np.log(1-q + q*np.exp(-np.exp(beta[0] + beta[1]*Age + beta[2]*BMI + beta[3]*Smoking)*t)))
    
    return ll

def logprior(q, beta):
    """ computes log prior likelihood """
    
    sigmas = np.array([10,10,10,10])
    
    lp = -np.sum((beta**2) / (2*sigmas))
    
    return lp

# ======= MCMC =======

def MCMC(n, qzero, betazero, sigma):
    """ runs MCMC sampler """
    
    accept = 0
    reject = 0
    
    theta = np.zeros((5, n+1))
    theta[0:4,0] = betazero
    theta[4,0] = qzero
    
    for i in range(n):
        # propose a new value
        thetaprop = np.random.multivariate_normal(mean=theta[:,i], cov=sigma)
        
        # set current and proposal values
        betacurrent = theta[0:4,i]
        qcurrent = theta[4,i]
        
        betaprop = thetaprop[0:4]
        qprop = thetaprop[4]
        
        # compute log acceptance probability
        lap = loglikelihood(qprop, betaprop, t, delta, Age, BMI, Smoking)
        lap -= loglikelihood(qcurrent, betacurrent, t, delta, Age, BMI, Smoking)
        lap += logprior(qprop, betaprop) - logprior(qcurrent, betacurrent)
        
        # accept-reject step
        u = np.random.uniform(0,1)
        
        if np.log(u)<lap:
            theta[:,i+1] = thetaprop
            accept += 1
        else:
            theta[:,i+1] = theta[:,i]
            reject += 1
    
    acceptrate = accept/n    
    
    betasample = theta[0:4,:]
    qsample = theta[4,:]
    
    plt.figure()
    plt.plot(qsample, zorder=1)
    plt.title(r'$q$, acceptance rate = ' + str(acceptrate))
    plt.xlabel(r'$t$')
    plt.ylabel(r'$q_t$')
        
    plt.figure()
    plt.hist(qsample, bins=50, zorder=1)
    plt.title(r'$q$')
    plt.xlabel(r'$q$')
    plt.ylabel('frequency')
        
    plt.figure()
    plt.subplot(2,2,1)
    plt.plot(betasample[0,], zorder=1)
    plt.title(r'$\beta^0$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\beta^0_t$')
    
    plt.subplot(2,2,2)
    plt.plot(betasample[1,], zorder=1)
    plt.title(r'$\beta^1$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\beta^1_t$')
    
    plt.subplot(2,2,3)
    plt.plot(betasample[2,], zorder=1)
    plt.title(r'$\beta^2$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\beta^2_t$')
    
    plt.subplot(2,2,4)
    plt.plot(betasample[3,], zorder=1)
    plt.title(r'$\beta^3$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\beta^3_t$')
    plt.tight_layout()   
 
    plt.figure()
    plt.subplot(2,2,1)
    plt.hist(betasample[0,], bins=50, zorder=1)
    plt.title(r'$\beta^0$')
    plt.xlabel(r'$\beta^0$')
    plt.ylabel('frequency')
        
    plt.subplot(2,2,2)
    plt.hist(betasample[1,], bins=50, zorder=1)
    plt.title(r'$\beta^1$')
    plt.xlabel(r'$\beta^1$')
    plt.ylabel('frequency')
            
    plt.subplot(2,2,3)
    plt.hist(betasample[2,], bins=50, zorder=1)
    plt.title(r'$\beta^2$')
    plt.xlabel(r'$\beta^2$')
    plt.ylabel('frequency')
            
    plt.subplot(2,2,4)
    plt.hist(betasample[3,], bins=50, zorder=1)
    plt.title(r'$\beta^3$')
    plt.xlabel(r'$\beta^3$')
    plt.ylabel('frequency')
    plt.tight_layout()  
    
    return theta

# sigma = np.diag([0.001,0.001,0.001,0.001,0.001])

betazero = np.array([-0.1,-0.1,-0.1,-0.1])
qzero = 0.5       

sigma = np.matrix([[ 2.11808832e-01, -6.01037185e-03, -7.28373028e-04,
        -8.48197294e-03,  2.56434160e-03],
       [-6.01037185e-03,  1.75957489e-04,  2.86458681e-05,
         2.24877669e-04, -1.13775757e-04],
       [-7.28373028e-04,  2.86458681e-05,  9.35238981e-04,
        -1.03977852e-03, -2.54804937e-04],
       [-8.48197294e-03,  2.24877669e-04, -1.03977852e-03,
         1.14239915e-01, -2.39218359e-03],
       [ 2.56434160e-03, -1.13775757e-04, -2.54804937e-04,
        -2.39218359e-03,  9.52041641e-04]])

betazero = np.array([-4.49523453, -0.03295658, -0.06429827, -1.06514478])
qzero = 0.7186742277757543      

# =============================================================================
# outputFile = open('TimeCovariatesData.csv','w',newline='')
# outputWriter = csv.writer(outputFile)
# =============================================================================

output = MCMC(10000, qzero, betazero, sigma)

# =============================================================================
# for k in range(5):
#         outputWriter.writerow(output[k,])
#         
# outputFile.close()
# 
# outputFile = open('TimeCovariatesDataDiagnostic.csv','w',newline='')
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
 
sigmanew = 5.6644*np.cov(output)/5
thetazeronew = output[:,-1]
betazeronew = thetazeronew[0:4]
qzeronew = thetazeronew[4]
outputnew = MCMC(1000, qzeronew, betazeronew, sigmanew)

sigmanew2 = 5.6644*np.cov(outputnew)/5
thetazeronew2 = outputnew[:,-1]
betazeronew2 = thetazeronew2[0:4]
qzeronew2 = thetazeronew2[4]
outputnew2 = MCMC(1000, qzeronew2, betazeronew2, sigmanew2)
        
sigmanew3 = 5.6644*np.cov(outputnew2)/5
thetazeronew3 = outputnew2[:,-1]
betazeronew3 = thetazeronew3[0:4]
qzeronew3 = thetazeronew3[4]
outputnew3 = MCMC(10000, qzeronew3, betazeronew3, sigmanew3)   