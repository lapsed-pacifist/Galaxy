import lmfit as lm
from scipy.optimize import leastsq
import numpy as np
import storeData as SD
import fitting as F
from copy import deepcopy as copy
import time
import matplotlib.pyplot as plt

def unpack_params(P):
	d = {}
	for i, v in P.iteritems():
		d[i] = v.value
	return [d['MeB'], d['ReB'],d['nB'],d['MeD'],d['ReD']]

def sersic_exp(p, x, z):
	MeB, ReB, nB, MeD, ReD = p
	IeB = 10 ** ((MeB - z) / -2.5); IeD = 10 ** ((MeD - z) / -2.5)
	r_factorB = (x / ReB) ** (1 / nB); r_factorD = (x / ReD) ** (1 / 1.)
	exponentB = (r_factorB - 1) * F.get_b_n(nB); exponentD = (r_factorD - 1) * F.get_b_n(1.)
	modelB = IeB * np.exp(-1.*exponentB); modelD = IeD * np.exp(-1.*exponentD)
	return modelD + modelB

def residual(p, X, Y, W, Z):
	return (Y - sersic_exp(p, X, Z))/W

def bootstrap(fit_params, profile, sigma=0.68, fit_range=None, sample_size=500, full=False):
	N = len(profile.I)
	I0 = copy(profile.I)
	R0 = copy(profile.R)
	W0 = copy(profile.W)
	p = unpack_params(fit_params)
	out_array = np.zeros((len(p), sample_size))
	for i in range(sample_size):
		count = 0
		out = [5]
		while (out[-1] > 4) and (count < 5):
			sample = np.random.randint(0,N,(N,))
			out = leastsq(residual, p, args=(R0[sample], I0[sample], W0[sample], profile.zeropoint), full_output=True)
			out_array[:, i] = out[0]
			count += 1
	stds = np.std(out_array, 1)
	means = np.mean(out_array,1)
	if full:		
		return stds, means, out_array, p
	else:
		return out_array

def get_means(par_array, pars=None, sigma=0.68):
	if pars is None:
		pars = np.mean(par_array, 1)
	var_array = np.zeros((len(pars),2))
	for j in range(len(pars)):
		below = par_array[j,:][par_array[j,:] < pars[j]]
		above = par_array[j,:][par_array[j,:] >= pars[j]]
		b = np.floor((len(below) * sigma))
		a = np.floor((len(above) * sigma))
		var_array[j, 0] = below[-1*b]
		var_array[j, 1] = above[a]
	return pars, var_array

def store_confidence(profile, var_array, sigma):
	d = {}
	for i, key in enumerate(['MeB', 'ReB', 'nB', 'MeD', 'ReD']):
		d[key] = var_array[i,:]
	profile.confs.update({sigma: d})

if __name__ == '__main__':
	direct = r'C:\Users\User\Documents\Project work\GalaxyFitting\repository'
	gal_list, gal_names = SD.import_directory(direct)
	target = gal_list[34][0][0]
	fit = F.fit_sersic_exp(target, [1., target.R[-1]+1], True)

	lm.report_fit(fit.params, show_correl=False)
	st, mu, out, params = bootstrap(fit.params, target, sample_size=500, full=True)
	print get_means(out)
	print get_means(out, params)
