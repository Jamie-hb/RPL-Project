""" Runs a MH sampler for the model
        T ~ q(α^TX)Exp(β) + (1 - q(α^TX))∞ 
    with q = 1/1+exp(-α^TX)
    With right censoring. """
    
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import DataSetup
from DataSetup import *
import csv

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

# Assemble data into a matrix
ones = np.repeat(1,SampleSize)
Covariates = np.stack((ones, Age, BMI, Smoking, Coffee))


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

def MCMC(n, alphazero, betazero, sigma):
    """ runs MCMC sampler 
        n = chain length
        alphazero, betazero = initial parameters
        sigma = proposal covariance matrix """
    
    accept = 0
    reject = 0
    
    theta = np.zeros((NumPara, n+1))
    theta[0:(NumPara - 1), 0] = alphazero
    theta[NumPara - 1, 0] = betazero
    
    for i in range(n):
        # propose a new value
        thetaprop = np.random.multivariate_normal(mean=theta[:,i], cov=sigma)
        
        # check for β* > 0, o/w immediately reject
        if thetaprop[NumPara - 1]>0: 
            # set the current and proposal values
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
                            
    acceptrate = accept/n    
    print('acceptance rate = ' + str(acceptrate))
    
    alphasample = theta[0:(NumPara - 1),:]
    betasample = theta[NumPara-1,:]
    
    # plot trace plots and histograms for each parameter
    plt.figure()
    plt.plot(betasample, zorder=1)
    plt.title(r'$\beta$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\beta_t$, acceptance rate = ' + str(acceptrate))
        
    plt.figure()
    plt.hist(betasample, bins=100, zorder=1)
    plt.title(r'$\beta$')
    plt.xlabel(r'$\beta$')
    plt.ylabel('frequency')
        
    plt.figure()
    plt.subplot(2,2,1)
    plt.plot(alphasample[0,], zorder=1)
    plt.title('alpha0')
    plt.xlabel('t')
    plt.ylabel('alpha0_t')
    
    plt.subplot(2,2,2)
    plt.plot(alphasample[1,], zorder=1)
    plt.title('alpha1')
    plt.xlabel('t')
    plt.ylabel('alpha1_t')
    
    plt.subplot(2,2,3)
    plt.plot(alphasample[2,], zorder=1)
    plt.title('alpha2')
    plt.xlabel('t')
    plt.ylabel('alpha2_t')
    
    plt.subplot(2,2,4)
    plt.plot(alphasample[3,], zorder=1)
    plt.title('alpha3')
    plt.xlabel('t')
    plt.ylabel('alpha3_t')
    plt.tight_layout()
    
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(alphasample[4,], zorder=1)
    plt.title(r'$\alpha^4$')
    plt.xlabel('t')
    plt.ylabel(r'$\alpha^4_t$')
    
    plt.subplot(1,2,2)
    plt.hist(alphasample[4,], bins=100, zorder=1)
    plt.title('alpha4')
    plt.xlabel('alpah4')
    plt.ylabel('frequency')
    plt.tight_layout()
    
    plt.figure()
    plt.subplot(2,2,1)
    plt.hist(alphasample[0,], bins=100, zorder=1)
    plt.title('alpha0')
    plt.xlabel('alpha0')
    plt.ylabel('frequency')
    
    plt.subplot(2,2,2)
    plt.hist(alphasample[1,], bins=100, zorder=1)
    plt.title('alpha1, acceptance rate = ')
    plt.xlabel('t')
    plt.ylabel('alpha1_t')
        
    plt.subplot(2,2,3)
    plt.hist(alphasample[2,], bins=100, zorder=1)
    plt.title('alpha2')
    plt.xlabel('frequency')
        
    plt.subplot(2,2,4)
    plt.hist(alphasample[3,], bins=100, zorder=1)
    plt.title('alpha3')
    plt.xlabel('alpha3')
    plt.ylabel('frequency')
    plt.tight_layout()
    
    return theta

alphazero = np.array([0.82883734, -0.22993219,  0.0520644 , -1.15676514, 0.1])
betazero = 0.005

# alphazero = np.random.multivariate_normal(mean=np.array([0,0,0,0,0]), cov=np.diag([10,0.3,0.3,3,1]))
# betazero = np.random.gamma(6,1)/1000  

# sigma = np.matrix([[ 1.06951773e-02, -2.37195739e-03, -5.91441444e-04,
#         -6.93491844e-03, -5.19761156e-06],
#        [-2.37195739e-03,  3.58179203e-03, -3.66887098e-04,
#          2.89835546e-03, -9.09784796e-07],
#        [-5.91441444e-04, -3.66887098e-04,  9.49466377e-04,
#          5.30150724e-04, -1.66729300e-07],
#        [-6.93491844e-03,  2.89835546e-03,  5.30150724e-04,
#          9.06762376e-02, -2.09865151e-06],
#        [-5.19761156e-06, -9.09784796e-07, -1.66729300e-07,
#         -2.09865151e-06,  7.75189753e-08]])

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

output = MCMC(10000, alphazero, betazero, sigma)
 
sigmanew = 5.6644*np.cov(output) / NumPara
thetazeronew = output[:,-1]
alphazeronew = thetazeronew[0:(NumPara-1)]
betazeronew = thetazeronew[NumPara-1]
outputnew = MCMC(1000, alphazeronew, betazeronew, sigmanew)

sigmanew2 = 5.6644*np.cov(outputnew) / NumPara
thetazeronew2 = outputnew[:,-1]
alphazeronew2 = thetazeronew2[0:(NumPara-1)]
betazeronew2 = thetazeronew2[NumPara-1]
outputnew2 = MCMC(1000, alphazeronew2, betazeronew2, sigmanew2)
        
sigmanew3 = 5.6644*np.cov(outputnew2) / NumPara
thetazeronew3 = outputnew2[:,-1]
alphazeronew3 = thetazeronew3[0:(NumPara-1)]
betazeronew3 = thetazeronew3[NumPara-1]
outputnew3 = MCMC(10000, alphazeronew3, betazeronew3, sigmanew3) 

# ===== Some script to store the chains in Excel for convergence analysis ===== 

outputFile = open('EventClassMCMC.csv','w',newline='')
outputWriter = csv.writer(outputFile)

for k in range(6):
        outputWriter.writerow(outputnew3[k,])
        
outputFile.close()

# outputFile = open('FertileCovariatesDataDiagnostic.csv','w',newline='')
# outputWriter = csv.writer(outputFile)
# 
# for i in range(10):
#     alphazero = np.random.normal(0, [10,0.3,0.3,3], 4)
#     betazero = np.random.uniform(6,1)/1000      
# 
#     output = MCMC(10000, alphazero, betazero, sigma)
#     
#     for k in range(5):
#         outputWriter.writerow(output[k,])
#         
# outputFile.close()