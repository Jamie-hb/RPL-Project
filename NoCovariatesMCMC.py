""" Runs a Gibbs sampler for the model
        T ~ qExp(β) + (1 - q)∞ 
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

# remove na entries:
# first find where the na entries are
workingset = np.logical_not( np.logical_or.reduce((np.isnan(t), np.isnan(delta))) )
workingset = np.logical_and(workingset, t>0)
# then remove them
t = t[workingset]
delta = delta[workingset]

# ======= Define the priors =======

# β ~ Gamma(a,b)
a = 1
b = 0.005

# p ~ U[0,1]

# f_i ~ 1{gotpregnant}


# ======= Define some functions to help compute the posterior parameters =======

# beta
def betapostb(t, delta, f):
    """ compute Σ_i t_i(δ_i + f_i(1 - δ_i)) """
    
    s = t*(delta + f*(1 - delta))
    res = np.sum(s)
    return res

# q
def qpostalpha(f, delta, Np):
    """ compute Np + Σ_i f_i(1 - δ_i) """
    
    s = f*(1 - delta)
    res = Np + np.sum(s) + 1
    return res

def qpostbeta(f, delta):
    """ compute Σ_i (1 - f_i)(1 - δ_i) + 1 """
    
    s = (1 - f)*(1 - delta)
    res = np.sum(s) + 1
    return res

# f
def fprob(i, q, beta, t, delta):
    """ compute 1 if δ_i = 1, and
        qexp(-βt_i) / (qexp(-βt_i) + 1 - q) if δ_i = 0 """
        
    if delta[i]==1:
        return 1
    else:
        num = q*np.exp(-beta*t[i]) 
        denom = q*np.exp(-beta*t[i]) + 1 - q
        res = num / denom
        return res


    
# ======= Define the Gibbs sampler =======

def Gibbs_fit(t, delta, n):
    """ Implements the Gibbs sampler on given data 
        t = observed pregnancy times
        delta = observed censoring
        n = how many iterations to run the Gibbs sampler for """
    
    # compute the amount of data
    N = t.size
    
    # set up the result structures and set inital conditions
    beta = [0]*(n+1)
    beta[0] = 0.005
    q = [0]*(n+1)
    q[0] = np.random.uniform(0,1)
    
    # set up the fertility prior
    f = np.copy(delta)
    
    # compute how many women got pregnant Np = Σ_i δ_i
    Np = np.sum(delta)
    
    # Gibbs sampling
    for i in range(n):
        beta[i+1] = np.random.gamma(a + Np, 1 / (b + betapostb(t, delta, f)))
        q[i+1] = np.random.beta(qpostalpha(f, delta, Np), qpostbeta(f, delta))
        
        for j in range(N):
            if delta[j]==0:
                f[j] = np.random.binomial(1, fprob(j, q[i+1], beta[i+1], t, delta))
           
    # plot the results for q
    plt.figure(1)
    plt.plot(q, zorder=0)
    plt.title('q')
    plt.xlabel('t')
    plt.ylabel(r'$q_t$')
        
    plt.figure(2)
    plt.hist(q)
    plt.title('q')
    plt.xlabel('q')
    plt.ylabel('relative frequency')
        
    # plot the results for β
    plt.figure(3)
    plt.plot(beta, zorder=0)
    plt.title('beta')
    plt.xlabel('t')
    plt.ylabel(r'$\beta_t$')
        
    plt.figure(4)
    plt.hist(beta)
    plt.title('beta')
    plt.xlabel('beta')
    plt.ylabel('relative frequency')
        
    return np.stack((q, beta))


output = Gibbs_fit(t, delta, 100000)
coeff = np.mean(output, axis=1)

# ===== This bit is to put the chains into Excel for convergence analysis =======

# simpledata = open('simpleoutput.csv','w',newline='')
# simpledatawriter = csv.writer(simpledata)
# simpledatawriter.writerow(output[0,:])
# simpledatawriter.writerow(output[1,:])
# simpledata.close()

# ===== Script to test the predicitons made by the model =====

def ConceptionProb(q, beta, t):
    
    ECProb = q
    ExpCDF = 1 - np.exp(-beta*t)
    
    return ECProb * ExpCDF

Times = np.array([365, 730, 1000])

Prediction = np.array([0,0,0])
for i in range(3):
    Prediction[i] = 1*(ConceptionProb(coeff[0], coeff[1], Times[i]) > 0.5)

AccuracyMetrics = np.zeros((3,5))

for i in range(3):
    if Prediction[i]==0:
        temp = sum((t>Times[i])) / (sum(t>Times[i]) + sum((t<=Times[i]) & (delta==1)))
        AccuracyMetrics[i,0] = temp
        AccuracyMetrics[i,1] = 0
        AccuracyMetrics[i,2] = 1
        AccuracyMetrics[i,3] = -1
        AccuracyMetrics[i,4] = temp
    
    if Prediction[i]==1:
        temp = sum((t<=Times[i]) & (delta==1)) / (sum(t>Times[i]) + sum((t<=Times[i]) & (delta==1)))
        AccuracyMetrics[i,0] = temp
        AccuracyMetrics[i,1] = 1
        AccuracyMetrics[i,2] = 0
        AccuracyMetrics[i,3] = temp
        AccuracyMetrics[i,4] = -1