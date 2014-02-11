import numpy as np
from scipy.optimize import minimize
import scipy.stats as stats
import time
import matplotlib.pyplot as plt
import random

def binning(data, binN):
	return np.histogram(data, bins=binN)

def log_like(hist, a, W):
	"""N=number of points, a=smoothing factor, W=width"""
	Sum = np.nansum(hist + a - 1)
	return np.nansum(np.log((hist + a - 1) / (W* (Sum))) * hist)

def regress(params, data):
	h = binning(data, params[0])[0]
	W = (np.max(data) - np.min(data)) / len(h)
	return -log_like(h, params[1], W)

if __name__ == '__main__':
	np.random.seed(43)
	X = np.random.normal(size=(1000))
	P = minimize(regress, [100., 10], args=(X,), method='nelder-mead').x
	print P
	H = binning(X, P[0])
	plt.bar(np.arange(len(H[0])), H[0])
	plt.show()