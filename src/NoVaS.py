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
	a = (np.ones(p+1))/(p+1)
	squared_returns = returns**2
	denominator = math.sqrt(np.dot(a, squared_returns.tail(p+1)[::-1]))	
	return np.divide(sp500_returns, denominator)

def exponential_novas(x):
	pass

def main():
	pass

if __name__ == "__main__": main()



