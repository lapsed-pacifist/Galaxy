"""
fit disc component 
append bulge parameters 
fit both
"""
import storeData as SD
import numpy as np
import matplotlib.pyplot as plt

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
	return tuple(ret)

def res_func(parameters, model_func, data, x, z, weights=None, fit_range=None):
	"""takes model_func and generates residuals from data"""
	res = data - model_func(parameters, x, z) 
	if weights is not None: 
		res /= weights
	if fit_range is not None: 
		fit_i = [SD.translate_x(x, i) for i in fit_range]
		return res[fit_i[0]:fit_i[1]]
	else:
		return res


def BD_routine():
	pass

if __name__ == '__main__':
	direct = r'repository'
	gal_list, gal_names = SD.import_directory(direct)
	G = gal_list[0][0][0]
	P = G.params
	X = np.linspace(0,60,1000)
	tot = sersic2(P, X, G.zeropoint, 'T')
	res = res_func(P, sersic2, G.I, G.R, G.zeropoint, G.W, None)
	plt.plot(X, convert_I(tot, G.zeropoint))
	# plt.plot(X, convert_I(bulge, G.zeropoint))
	# plt.plot(X, convert_I(disc, G.zeropoint))
	plt.show()