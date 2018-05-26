"""
Contains implementation of the GARCH(1,1) model for volatility prediction

Author: Nishant D. Gurnani
Date Created: 5/26/2018
Last Updated: 5/26/2018
"""
import cvxopt
from functools import partial
import math
import numpy as np
import scipy
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.stattools import jarque_bera
import matplotlib.pyplot as plt

def compute_squared_sigmas(X, initial_sigma, theta):
	"""Function to compute the sigmas given the initial guess"""

    a0 = theta[0]
    a1 = theta[1]
    b1 = theta[2]
    
    T = len(X)
    sigma2 = np.ndarray(T)
    
    sigma2[0] = initial_sigma ** 2
    
    for t in range(1, T):
        # Here's where we apply the equation
        sigma2[t] = a0 + a1 * X[t-1]**2 + b1 * sigma2[t-1]
    
    return sigma2

def negative_log_likelihood(X, theta):
    """takes as input returns X and model parameters theta and returns the negative
       log likelihood"""

    T = len(X)
    
    # Estimate initial sigma squared
    initial_sigma = np.sqrt(np.mean(X ** 2))
    
    # Generate the squared sigma values
    sigma2 = compute_squared_sigmas(X, initial_sigma, theta)
    
    # Now actually compute
    return -sum(
        [-np.log(np.sqrt(2.0 * np.pi)) -
        (X[t] ** 2) / (2.0 * sigma2[t]) -
        0.5 * np.log(sigma2[t]) for
         t in range(T)]
    )

def garch11(X,theta):
	objective = partial(negative_log_likelihood, X)

	# Define the constraints for our minimizer
	def constraint1(theta):
	    return np.array([1 - (theta[1] + theta[2])])

	def constraint2(theta):
	    return np.array([theta[1]])

	def constraint3(theta):
	    return np.array([theta[2]])

	cons = ({'type': 'ineq', 'fun': constraint1},
	        {'type': 'ineq', 'fun': constraint2},
	        {'type': 'ineq', 'fun': constraint3})

	# Actually do the minimization
	result = scipy.optimize.minimize(objective, (1, 0.5, 0.5),
	                        method='SLSQP',
	                        constraints = cons)

	theta_mle = result.x
	return theta_mle

def prediction_interval():
	pass

