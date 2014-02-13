"""
fit disc component 
append bulge parameters 
fit both
"""
import storeData as SD
import numpy as np
import matplotlib.pyplot as plt
import lmfit as lm
import matplotlib.gridspec as gridspec

def copy_params(parameters, trans_vary=True):
	new = lm.Parameters()
	for p in parameters.values():
		new.add(p.name, value=p.value, min=p.min, max=p.max)
		if trans_vary:
			new[p.name].vary=p.vary
	return new

def fix_params(p, fixed, unfix=False):
	if type(fixed) is list:
		for n in fixed:
			p[n].vary = unfix
	else:
		for n, v in fixed.iteritems():
			p[n].vary = unfix
			if v is not None:
				p[n].value = v

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

def sersic2(parameters, x, z, comp='T'):
	MeB, ReB, nB, MeD, ReD, nD = [parameters[i].value for i in ['MeB', 'ReB', 'nB', 'MeD', 'ReD', 'nD']]
	IeB = 10 ** ((MeB - z) / -2.5); IeD = 10 ** ((MeD - z) / -2.5)
	r_factorB = (x / ReB) ** (1 / nB); r_factorD = (x / ReD) ** (1 / nD)
	exponentB = (r_factorB - 1) * get_b_n(nB); exponentD = (r_factorD - 1) * get_b_n(nD)
	modelB = IeB * np.exp(-1.*exponentB); modelD = IeD * np.exp(-1.*exponentD)
	comp_dict = {'D': modelD, 'B':modelB, 'T':modelB+modelD}
	ret = [comp_dict[i] for i in comp]
	if len(comp) == 1:
		return comp_dict[comp]
	return tuple(ret)

def res_func(parameters, model_func, x, z, data, weights=None, fit_range=None, comp='T'):
	"""takes model_func and generates residuals from data"""
	res = data - model_func(parameters, x, z, comp) 
	if weights is not None: 
		res /= weights
	if fit_range is not None: 
		fit_i = [SD.translate_x(x, i) for i in fit_range]
		return res[fit_i[0]:fit_i[1]]
	else:
		return res


def fit_sersic_exp(profile, fit_range=None, store=False):
	"""fits exp to middle then fixes and fits sersic-exp"""
	mid_section = [20., profile.R[-1] - 10.]

	P = copy_params(profile.params)
	fix_params(P, {'nB':None, 'MeB':None, 'ReB':None, 'nD':1.}) # fit exp
	exp_fit = lm.minimize(res_func, P, args=(sersic2, profile.R, profile.zeropoint, profile.I, profile.W, mid_section, 'D'))

	# P = copy_params(exp_fit.params, False)
	# fix_params(P, {'nD':1., 'MeD':None, 'ReD': None}) # fit temp bulge
	# se_fit = lm.minimize(res_func, P, args=(sersic2, profile.R, profile.zeropoint, profile.I, profile.W, None, 'T'))

	P = copy_params(exp_fit.params, False)
	fix_params(P, {'nD':1.})
	free_fit = lm.minimize(res_func, P, args=(sersic2, profile.R, profile.zeropoint, profile.I, profile.W, None, 'T'))

	if store:
		profile.fits.update({'sersic+exp': [free_fit, fit_range]})
	return free_fit

def find_redchi(profile, fit_data):
	total = sersic2(fit_data.params, profile.R, profile.zeropoint, 'T')
	return np.sum(((profile.I - total) / profile.W)**2.) / (fit_data.ndata)

if __name__ == '__main__':
	direct = 'repository'
	Gal_list, Gal_names = SD.import_directory(direct)
	G = Gal_list[0][0][1]
	plt.errorbar(G.R, G.I, G.W, 'b.')
	# P = G.params
	# X = np.linspace(0,60,1000)
	# out = fit_sersic_exp(G, store=True)
	# print find_redchi(G, out)
	# print out.redchi
	# lm.report_fit(out.params, show_correl=False)

	# fig = plt.figure()
	# gs = gridspec.GridSpec(6,6)
	# ax = fig.add_subplot(gs[0:4,0:]) # main plot
	# ax.errorbar(G.R, G.M, yerr=G.MW, fmt='b.')
	# mod_R = np.linspace(0, G.R[-1], 500)
	# total, bulge, disc = sersic2(out.params, mod_R, G.zeropoint, 'TBD')
	# ax.plot(mod_R, convert_I(total, G.zeropoint), 'k-')
	# ax.plot(mod_R, convert_I(bulge, G.zeropoint), 'g:', linewidth=2.)
	# ax.plot(mod_R, convert_I(disc, G.zeropoint), 'r--', linewidth=2.)
	# ax.set_ylim([35,15])
	# ax.set_title(G.gal_name+G.cam_name+G.name)

	# total = sersic2(out.params, G.R, G.zeropoint, 'T')
	# res = fig.add_subplot(gs[4:,:], sharex=ax)
	# ax.xaxis.set_visible(False)
	# resi = G.M - convert_I(total, G.zeropoint)
	# res.plot(G.R, resi,  'b.')
	# res.axhline(y=0, linestyle='--')
	# res.set_ylim([-2,2])
	plt.show()