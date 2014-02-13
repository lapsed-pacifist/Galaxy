"""
1) fit exp and then sersic-exp with fixed exp
"""
from scipy import stats
import storeData as SD
import numpy as np
import matplotlib.pyplot as plt
import lmfit as lm
import matplotlib.gridspec as gridspec

def add_bounds(p, bound):
	for n, v in bound.iteritems():
		if np.isinf(p[n].min):
			p[n].min = (1-v)*p[n].value
		else:
			p[n].min = max([(1-v)*p[n].value, p[n].min])
		if np.isinf(p[n].max):
			p[n].max = (1+v)*p[n].value
		else:
			p[n].max = min([(1+v)*p[n].value, p[n].max])

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

def fix_params(p, fixed, unfix=False):
	if type(fixed) is list:
		for n in fixed:
			p[n].vary = unfix
	else:
		for n, v in fixed.iteritems():
			p[n].vary = unfix
			if v is not None:
				p[n].value = v

def fit_sersic_exp(profile, fit_range=None, store=False):
	"""fits exp to middle then fixes and fits sersic-exp"""

	P = copy_params(profile.params)
	fix_params(P, {'nB':None, 'MeB':None, 'ReB':None, 'nD':1.}) # fit exp
	exp_fit = lm.minimize(residual, P, args=(sersic, profile.R, profile.zeropoint, profile.I, 1., [20, 40]))

	P = copy_params(exp_fit.params, False)
	fix_params(P, {'nD':1., 'MeD':None, 'ReD': None}) # fit temp bulge
	se_fit = lm.minimize(residual, P, args=(sersic2, profile.R, profile.zeropoint, profile.I, 1.))

	# P = copy_params(se_fit.params, False)
	# fix_params(P, {'nD':1., 'MeD':None, 'ReD': None, 'nB':4.}) # fit vaucouleurs
	# se_fit = lm.minimize(residual, P, args=(sersic2, profile.R, profile.zeropoint, profile.I, profile.W))

	# P = copy_params(se_fit.params, False)
	# fix_params(P, {'nD':1., 'MeD':None, 'ReD': None,}) # fit sersic bulge
	# se_fit = lm.minimize(residual, P, args=(sersic2, profile.R, profile.zeropoint, profile.I, profile.W))

	if store:
		profile.fits.update({'sersic+exp': [se_fit, fit_range]})
	return se_fit

def find_confs_se(profile, fit_name, conf_int=None, fit=None):
	if fit is not None:
		P = copy_params(fit.params)
	else:
		P = copy_params(profile.fits[fit_name][0].params, False)
	fix_params(P, {'ReB':None, 'MeB':None, 'nB':None, 'nD':1.}) # improve index
	disc_fit = lm.minimize(residual, P, args=(sersic2, profile.R, profile.zeropoint, profile.I, profile.W, profile.fits[fit_name][1]))

	if fit is not None:
		P = copy_params(fit.params)
	else:
		P = copy_params(profile.fits[fit_name][0].params, False)
	fix_params(P, {'ReD':None, 'MeD':None, 'nD':1.}) # improve index
	bulge_fit = lm.minimize(residual, P, args=(sersic2, profile.R, profile.zeropoint, profile.I, profile.W, profile.fits[fit_name][1]))

	if conf_int is None:
		ci_disc = lm.conf_interval(disc_fit, sigmas=[.68])
		ci_bulge = lm.conf_interval(bulge_fit, sigmas=[.68])
	else:
		ci_disc = lm.conf_interval(disc_fit, sigmas=[conf_int])
		ci_bulge = lm.conf_interval(bulge_fit, sigmas=[conf_int])

	# xnb, ynb, gridnb = lm.conf_interval2d(bulge_fit,'MeB','ReB',100,100)
	# fig = plt.figure()
	# fig.set_facecolor('white')
	# ax1 = fig.add_subplot(111)
	# ax1.contourf(xnb, ynb, gridnb, np.linspace(0,1,11))
	# ax1.set_xlabel('MeB/ mag arcsec$^{-1}$'); ax1.set_ylabel('ReB/ arcsec')
	# ax1.set_title('$\chi^2$ for varying parameters')
	# plt.show()

	confidence = {}
	for ci in [ci_bulge, ci_disc]:
		for key, item in ci.iteritems():
			if conf_int is None:
				confidence[key] = [item[1][1], item[1][1]]
			else:
				confidence[key] = [item[1][1]-item[0][1], item[2][1]-item[1][1]]
	if conf_int is None:
		confidence['nD'] = [1., 1.]
	else:
		confidence['nD'] = [0., 0.]
	profile.confs[fit_name] = confidence
	return confidence

def fit_vauc(profile, fit_range=None):
	"""fits exp to middle then fixes and fits sersic-exp"""
	P = copy_params(profile.params)
	fix_params(P, {'nB':None, 'MeB':None, 'ReB':None, 'nD':1.})
	exp_fit = lm.minimize(residual, P, args=(sersic, profile.R, profile.zeropoint, profile.I, profile.W, [10, 40]))

	P = copy_params(exp_fit.params, False)
	fix_params(P, {'nD':1., 'nB':4.})
	se_fit = lm.minimize(residual, P, args=(sersic2, profile.R, profile.zeropoint, profile.I, profile.W, fit_range))
	return se_fit

def multi_residuals(parameters, model_func, profile_list, fit_range):
	for i, prof in enumerate(profile_list):
		r = residual(parameters, model_func, prof.R, prof.zeropoint, prof.I, prof.W, fit_range)
		if i == 0:
			res_array = r
		else:
			res_array = np.concatenate((res_array, r))
	return res_array

def fit_2_sersic_exp(profile_list, fit_range=None, redchi=False):
	P = copy_params(profile_list[0].params)
	fix_params(P, {'nB':None, 'MeB':None, 'ReB':None, 'nD':1.})
	exp_fit = lm.minimize(multi_residuals, P, args=(sersic, profile_list, [10,40]))

	P = copy_params(exp_fit.params, False)
	fix_params(P, {'nD':1., 'MeD':None, 'ReD': None}) # fit temporary bulge
	se_fit = lm.minimize(multi_residuals, P, args=(sersic2, profile_list, fit_range))
	
	P = copy_params(se_fit.params, False)
	fix_params(P, {'nB':None, 'MeB': None, 'ReB':None}) # improve constants
	free_fit = lm.minimize(multi_residuals, P, args=(sersic2, profile_list, fit_range))

	P = copy_params(free_fit.params, False)
	# fix_params(P, {'MeB':None, 'MeD':None, 'ReB':None, 'ReD':None}) # improve index
	fix_params(P, {'nD':1., 'nB':None})
	co_fit = lm.minimize(multi_residuals, P, args=(sersic2, profile_list, fit_range))

	redchi = []
	for prof in profile_list:
		total = sersic2(free_fit.params, prof.R, prof.zeropoint)
		redchi.append( np.sum(((prof.I - total) / prof.W)**2.) / (len(prof.I)-free_fit.nvarys) )
	return free_fit, redchi

def sky_error_se(profile, sky_adj, fit_name):
	"""
	Adjusts the sky and generates fit for variables from the adjustment. 
	Can take the name of fit stored in the profile. Specify addition variables
	to fix in the dictionary fixed. Use overwrite to clear previous fixed.
	"""
	fit_data = profile.fits[fit_name]
	P = copy_params(fit_data[0].params, False)
	diff = {}
	profile.I += sky_adj*1.
	fit_sersic_exp(profile, fit_data[1], store=False)
	outup = find_confs_se(profile, fit_name)
	profile.I -= 1.*sky_adj
	outdown = find_confs_se(profile, fit_name)
	profile.I += 1.*sky_adj
	for key, val in outup.iteritems():
		up = outup[key][0] - P[key].value
		down = P[key].value - outdown[key][0]
		if up > 0 and down > 0:
			up = max((up, down))
			down = 0.
		elif up < 0 and down < 0:
			up = 0.
			down = min((up, down))
		elif up < 0 and down > 0:
			up, down = down, up
		diff.update({key: [up, -down]})
	profile.sky_fit[fit_name] = diff
	return diff

def find_redchi(profile, fit_data):
	total = sersic2(fit_data.params, profile.R, profile.zeropoint)
	# return stats.kstest((profile.I - total)/profile.W, 'norm')[1]
	return np.sum(((profile.I - total) / profile.W)**2.) / (fit_data.ndata)

def durbin_watson(residuals):
	i1 = np.array([residuals[i+1] - residuals[i] for i,v in enumerate(residuals[:-1])])
	return np.sum((i1) ** 2.) / np.sum(residuals ** 2.)

if __name__ == '__main__':
	direct = 'C:\\Users\\User\\code\\Python\\Sersic_Plots'
	gal_list, gal_names = SD.import_directory(direct)
	gal = gal_list[7][1][0]	
	# gal.params['nB'].max = 6.; gal.params['nD'].max = 6.

	# out, redchis = fit_2_sersic_exp([gal_list[2][0][0], gal_list[2][0][1]])
	out = fit_sersic_exp(gal, store=True)
	print find_redchi(gal, out)
	lm.report_fit(out.params, show_correl=False)
	find_confs_se(gal, 'sersic+exp', 0.68)
	# print sky_error_se(gal, 10., 'sersic+exp')


	fig = plt.figure()
	gs = gridspec.GridSpec(6,6)
	ax = fig.add_subplot(gs[0:4,0:]) # main plot
	ax.errorbar(gal.R, gal.M, yerr=gal.MW, fmt='b.')
	mod_R = np.linspace(0, gal.R[-1], 500)
	total, bulge, disc = sersic2(out.params, mod_R, gal.zeropoint, True)
	ax.plot(mod_R, convert_I(total, gal.zeropoint), 'k-')
	ax.plot(mod_R, convert_I(bulge, gal.zeropoint), 'g:', linewidth=2.)
	ax.plot(mod_R, convert_I(disc, gal.zeropoint), 'r--', linewidth=2.)
	ax.set_ylim([35,15])
	ax.set_title(gal.gal_name+gal.cam_name+gal.name)

	total = sersic2(out.params, gal.R, gal.zeropoint)
	res = fig.add_subplot(gs[4:,:], sharex=ax)
	ax.xaxis.set_visible(False)
	resi = gal.M - convert_I(total, gal.zeropoint)
	res.plot(gal.R, resi,  'b.')
	res.set_ylim([-2,2])

	# histo = fig.add_subplot(gs[4:,4:])
	# histo.yaxis.set_visible(False)
	# histo.hist(resi,10,range=(-2,2), orientation='horizontal')

	# figI = plt.figure()
	# axI = figI.add_subplot(111)
	# axI.errorbar(gal.R[-15:], gal.I[-15:], yerr=gal.W[-15:], fmt='b.')


	plt.show()