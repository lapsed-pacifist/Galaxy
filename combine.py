"""
combine fit params using propagation and warn if values are not within a tolerance 


"""

import storeData as SD
import numpy as np
import bootstrap as b
import fitting as F
import matplotlib.pyplot as plt
import lmfit as lm
import plotData

def weighted_mean(*args):
	"""an arg is a list of [value, error]"""
	ws, xs = np.zeros((len(args),)), np.zeros((len(args),))
	for i, v in enumerate(args):
		ws[i] = 1./v[1] ** 2.
		xs[i] = v[0]
	xce = np.sum(ws * xs) / np.sum(ws)
	ace = 1. / np.sum(ws)
	return [xce, ace]

def combine_bootstrap(boot_list):
	"""combine each axis' camera and then all"""
	comb = None
	for i in boot_list:
		if comb is None:
			comb = i
		else:
			np.hstack((comb, i))

if __name__ == '__main__':
	# direct = r'C:\Users\User\Documents\Project work\GalaxyFitting\repository'
	# gal_list, gal_names = SD.import_directory(direct)
	# target = gal_list[34]
	# comb = None
	# for c in target:
	# 	for prof in c:
	# 		fit = F.fit_sersic_exp(prof, [1., prof.R[-1]+1], True)
	# 		lm.report_fit(fit.params, show_correl=False)
	# 		st, mu, out, params = b.bootstrap(fit.params, prof, sample_size=200, full=True)
	# 		if comb is None:
	# 			comb = out
	# 		else:
	# 			np.hstack((comb, out))
	# fig = plt.figure()
	# ax = fig.add_subplot(111)
	# hist, bins = np.histogram(comb[0,:], bins=80)
	# width = 0.95 * (bins[1] - bins[0])
	# center = (bins[:-1] + bins[1:]) / 2
	# ax.bar(center, hist, align='center', width=width)
	# print b.get_means(comb)

	# fig1 = plt.figure()
	# ax = fig.add_subplot(111)

	# plt.show()