""" Simulates data from a simple cure rate model
        T ~ qExp(β) + (1 - q)∞ 
    With right censoring. Then tests a Gibbs sampler for 
    infernece of q and β"""

import numpy as np
import matplotlib.pyplot as plt

# ======= Simulate the data =======

def simulate(N, q, beta):
    """ simulates N data points from the distribution
        T ~ qExp(β) + (1 - q)∞, with right censoring
        C ~ Uniform[1,1200] """
    
    # simulate the actual data
    T = np.random.exponential(1/beta, N) # compute the times till pregnancy
    C = np.random.uniform(1,1200, N) # compute the truncation times
    f = np.random.binomial(1, q, N) # compute the fertilities
    
    # compute censoring information (for 100% fertility)
    event_happens = T < C
    event_happens = event_happens.astype(int)
    
    # compute the observed data
    # if f = 0 (<=> 1 - f = 1), we observe C
    # if f = 1, delta = 1, we observe T
    #         , delta = 0, we observe C
    t = (1 - f)*C + f*(event_happens*T + (1 - event_happens)*C)
    
    # computed observed censoring information
    delta = t < C
    delta = delta.astype(int)
    
    return np.stack((t, delta))


# ======= Define the priors =======

# β ~ Gamma(a,b)
a = 1
b = 0.005

# p ~ U[0,1]
# f_i ~ 1{gotpregnant}


# ======= Set up some functions to compute the posterior parameters =======

# beta
def betapostb(t, delta, f):
    """ compute Σi t_i(δ_i + f_i(1 - δ_i)) """
    
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

def Gibbs_fit(t, delta, n, betazero, qzero):
    """ Implements the Gibbs sampler on given data 
        t = observed pregnancy times
        delta = observed censoring
        n = how many iterations to run the Gibbs sampler for
        betazero = inital value for β
        qzero = initial value for q """
    
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
    
    for i in range(n):
        beta[i+1] = np.random.gamma(a + Np, 1 / (b + betapostb(t, delta, f)))
        q[i+1] = np.random.beta(qpostalpha(f, delta, Np), qpostbeta(f, delta))
        
        for j in range(N):
            if delta[j]==0:
                f[j] = np.random.binomial(1, fprob(j, q[i+1], beta[i+1], t, delta))
    
    return np.stack((q, beta))


# ======= Test the Gibbs sampler =======

def test(m, n, N):
    for i in range(m):
        # set the true values
        qtrue = 0.6
        betatrue = 0.005
        
        print('qtrue = ' + str(qtrue))
        print('betatrue = ' + str(betatrue))
        
        # simulate a data set
        data = simulate(N, qtrue, betatrue) 
        t = data[0,]
        delta = data[1,]
        #f = data[2,]
        
        # some plots of the data
        plt.figure(1)
        plt.hist(t)
        plt.title('data - T')
        plt.xlabel('time')
        plt.ylabel('frequency')
        
        plt.figure(2)       
        plt.hist(delta)
        plt.title('data - δ')
        plt.xlabel('delta')
        plt.ylabel('frequency')

        
        # run the Gibbs sampler
        sample = Gibbs_fit(t, delta, n, betatrue, qtrue)
        q = sample[0,]
        beta = sample[1,]
        
        # plot the results for q
        plt.figure(4)
        plt.plot(q, zorder=0)
        plt.hlines(qtrue, 0, n, colors='y', zorder=10)
        plt.title('q')
        plt.xlabel('t')
        plt.ylabel(r'$q_t$')
        
        plt.figure(5)
        plt.hist(q, zorder=0)
        plt.axvline(qtrue, color='y', zorder=1000)
        plt.title('q')
        plt.xlabel('q')
        plt.ylabel('relative frequency')
        
        # plot the results for β
        plt.figure(6)
        plt.plot(beta, zorder=0)
        plt.hlines(betatrue, 0, n, colors='y', zorder=1000)
        plt.title('beta')
        plt.xlabel('t')
        plt.ylabel(r'$\beta_t$')
        
        plt.figure(7)
        plt.hist(beta, zorder=0)
        plt.axvline(betatrue, color='y', zorder=1000)
        plt.title('beta')
        plt.xlabel('beta')
        plt.ylabel('relative frequency')


test(1, 10000, 1500)

