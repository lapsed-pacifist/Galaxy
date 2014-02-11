"""
Objective: find conditions for a galaxy to be included in the catalogue

create galaxy
add noise
fit galaxy 

test initial conditions at different noise levels
"""
import numpy as np
import storeData as SD
import lmfit as lm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

def fix_params(p, fixed, unfix=False):
	if type(fixed) is list:
		for n in fixed:
			p[n].vary = unfix
	else:
		for n, v in fixed.iteritems():
			p[n].vary = unfix
			if v is not None:
				p[n].value = float(v)

def copy_params(parameters, trans_vary=True):
	new = lm.Parameters()
	for p in parameters.values():
		new.add(p.name, value=p.value, min=p.min, max=p.max)
		if trans_vary:
			new[p.name].vary=p.vary
	return new

def convert_I(I, z):
	return z - (2.5 * np.log10(I))

def get_b_n(m):
	#   find sersic b_n coefficient'favoritecolor'
	#   this is more accurate than the usual b_n = 1.9992*n_sersic - 0.3271
	if m == 0.0: return -0.3271
	b_n = 2.0*m - 1.0/3.0 + 4.0/m/405.0 + 46.0/m/m/25515.0 \
				+ 131.0/m/m/m/1148175.0  \
				- 2194697.0/m/m/m/m/30690717750.0
	return b_n

def sersic(parameters, x, z):
	if type(parameters) is list:
		Me, Re, n = parameters
	else: 
		Me, Re, n = [parameters[i].value for i in ['MeD', 'ReD', 'nD']]
	Ie = 10 ** ((Me - z) / -2.5)
	r_factor = (x / Re) ** (1 / n)
	exponent = (r_factor - 1) * get_b_n(n)
	model = Ie * np.exp(-1.*exponent)
	return model

def sersic2(parameters, x, z, comp=False):
	MeB, ReB, nB, MeD, ReD, nD = [parameters[i].value for i in ['MeB', 'ReB', 'nB', 'MeD', 'ReD', 'nD']]
	IeB = 10 ** ((MeB - z) / -2.5); IeD = 10 ** ((MeD - z) / -2.5)
	r_factorB = (x / ReB) ** (1 / nB); r_factorD = (x / ReD) ** (1 / nD)
	exponentB = (r_factorB - 1) * get_b_n(nB); exponentD = (r_factorD - 1) * get_b_n(nD)
	modelB = IeB * np.exp(-1.*exponentB); modelD = IeD * np.exp(-1.*exponentD)
	if comp:
		return modelB + modelD, modelB, modelD
	else:
		return modelB + modelD
	

def residual(parameters, model_func, x, z, data, weights=None, fit_range=None):
	diff = (data - model_func(parameters, x, z))
	if weights is None:
		result = diff
	else: 
		result = diff / weights
	if fit_range is None:
		return result
	else:
		fit_i = [SD.translate_x(x, i) for i in fit_range]
		return result[fit_i[0]:fit_i[1]]

def fit(P, model, x, z, data, weights=None, fit_range=None, redchi_marker=None):
	pars = copy_params(P, False)
	fix_params(pars, {'nD':1., 'nB':None, 'ReB':None, 'MeB':None})
	exp_out = lm.minimize(residual, P, args=(sersic, x, z, data, weights, [10, x[-1]]))

	pars = copy_params(exp_out.params, False)
	fix_params(pars, {'nD':1.})
	out = lm.minimize(residual, P, args=(sersic2, x, z, data, weights, fit_range))
	if redchi_marker is not None:
		# total = sersic2(out.params, x, z, data, weights, fit_range, False)
		res_excl = residual(out.params, sersic2, x, z, data, weights, [0, redchi_marker])
		return out, res_excl
	else: return out


if __name__ == '__main__':
	# direct = 'C:\\Users\\User\\code\\Python\\Sersic_Plots'
	# gal_list, gal_names = SD.import_directory(direct)
	# temp_gal = SD.Galaxy('temp')
	# temp_gal.import_file('C:\\Users\\User\\code\\Python\\Sersic_Plots\\mega_1237665428090192022_cts.ascii', 'mega')
	# R = temp_gal[0][0].R
	# W = temp_gal[0][0].W
	zp = 30.
	R = np.linspace(0.1, 50, 500)
	P = lm.Parameters()
	P.add_many(('MeB', 26.), ('ReB', 0.5), ('nB', 4.), ('MeD', 21.), ('ReD', 9.), ('nD', 1.))
	total1, bulge1, disc1 = sersic2(P, R, zp, True)
	# noise = np.random.random(total1.shape) * 50.
	# prof = SD.Profile('test', total1+noise, R, W, zp, 0, 1)

	fig = plt.figure()
	# gs = gridspec.GridSpec(6,2)
	# ax = fig.add_subplot(gs[:4,0])
	# ax2 = fig.add_subplot(gs[:4,1])
	# ax2.xaxis.set_visible(False)
	# res1 = fig.add_subplot(gs[4:,0])
	# res2 = fig.add_subplot(gs[4:,1])
	# ax.plot(prof.R, prof.M, '.b')
	# ax.plot(prof.R, convert_I(bulge1, zp), 'g:')
	# ax.plot(prof.R, convert_I(disc1, zp), 'r--')
	# ax.set_ylim([35,15])

	# P['nD'].vary = False
	# P['ReD'].value = 10.; P['MeD'].value = 22.; P['nB'].value = 3.; P['ReB'].value = 4.; P['MeB'].value = 14.
	# out, res = fit(P, sersic2, R, zp, prof.I, redchi_marker=30.)
	# print np.sum(res) / out.nfree

	# total, bulge, disc = sersic2(out.params, R, zp, True)
	# ax2.plot(prof.R, prof.M, 'b.')
	# ax2.plot(prof.R,  convert_I(total, zp), 'k-')
	# ax2.plot(prof.R,  convert_I(disc, zp), 'r--')
	# ax2.plot(prof.R,  convert_I(bulge, zp), 'g:')
	# ax2.set_ylim(ax.get_ylim())
	# res1.plot(prof.R, prof.M - convert_I(total1, zp), 'b.')
	# res2.plot(prof.R, prof.M - convert_I(total, zp), 'b.')
	# res2.set_ylim([-2,2])
	# res1.set_ylim([-2,2])
	# lm.report_fit(out.params, show_correl=False)
	ax = fig.add_subplot(111)
	ax.plot(R, convert_I(total1, zp), 'k-')
	ax.plot(R, convert_I(disc1, zp), 'b--')
	ax.plot(R, convert_I(bulge1, zp), 'g:')

	ax.set_ylim([35, 15])

	plt.show()