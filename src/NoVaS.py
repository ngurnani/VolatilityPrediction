"""
Contains implementation of the Simple NoVaS and Exponential NoVaS transformations 
for volatility prediction.

Author: Nishant D. Gurnani
Date Created: 5/19/2018
Last Updated: 5/19/2018
"""

import numpy as np
import pandas as pd
import math 


def empirical_kurtosis(x):
	""" Calculates empirical kurtosis of input array of data """
	n = len(x)
	numerator = (sum(np.subtract(x,np.mean(x)))**4)/n
	denominator = ((sum(np.subtract(x,np.mean(x)))**2)/n)**2
	return numerator/denominator


def simple_novas(x):
	pass

def exponential_novas(x):
	pass

def main():
	pass

if __name__ == "__main__": main()



