"""
Contains implementation of the Simple NoVaS and Exponential NoVaS transformations 
for volatility prediction.

Author: Nishant D. Gurnani
Date Created: 5/19/2018
Last Updated: 5/19/2018
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kurtosis

def simple_novas(returns, p):
    """
    Function performs NoVaS (normalizing and variance stabilizing transformation) on financial returns series
    as outlined in Politis (2007)
    
    input:
        returns - daily financial returns
        p - lag parameter (first p-1 entries of o)
    output:
        W_t - NoVaS transformed series
    """
    
    n = len(returns)
    X_t = np.zeros(2*p+n) # placeholder vector containing returns upto period t
    X_t[2*p:2*p+n] = returns
    X_t = pd.Series(X_t)
    
    # to calculate denominator term, sum lag p squared returns
    lagged_squared_returns = pd.Series(np.zeros(len(X_t))) 
    for i in range(p):
        lagged_squared_returns = lagged_squared_returns + (X_t**2).shift(-i)
    
    W_t = X_t/np.sqrt(lagged_squared_returns/p)
    W_t = W_t[2*p:2*p+n] # drop extra lag indices in array
    W_t = W_t.dropna() # drop any nan values
    return W_t

def simple_novas_prediction(returns, p):
    """
    Function performs NoVaS (normalizing and variance stabilizing transformation) on financial returns series
    as outlined in Politis (2007) and then compute the optimal one-step ahead L1 prediction
    
    input:
        returns - daily financial returns
        p - lag parameter (first p-1 entries of o)
    output:
        W_t - NoVaS transformed series
    """
    
    ## NoVaS transform
    
    n = len(returns)
    X_t = np.zeros(2*p+n) # placeholder vector containing returns upto period t
    X_t[2*p:2*p+n] = returns
    X_t = pd.Series(X_t)
    
    # to calculate denominator term, sum lag p squared returns
    lagged_squared_returns = pd.Series(np.zeros(len(X_t))) 
    for i in range(p):
        lagged_squared_returns = lagged_squared_returns + (X_t**2).shift(-i)
    
    W_t = X_t/np.sqrt(lagged_squared_returns/p)
    W_t = W_t[2*p:2*p+n] # drop extra lag indices in array
    W_t = W_t.dropna() # drop 
    
    # L1 prediction
    
    An2 = (lagged_squared_returns/p)
    An2 = An2[2*p:2*p+n]
    An2 = An2.dropna()
    An2 = An2.iloc[-1]

    a0 = 1/p
    W2 = W_t**2/(1-a0*(W_t**2))
    mu2 = np.nanmedian(W2)
    predictor = mu2*An2
    
    return predictor


def exponential_novas(x):
	pass

def main():
	pass

if __name__ == "__main__": main()



