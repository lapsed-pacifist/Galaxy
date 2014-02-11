import storeData as SD
from tests.basic_testing import *
import matplotlib.pyplot as plt

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
	print durbin_watson(resi)
	res.plot(gal.R, resi,  'b.')
	res.set_ylim([-2,2])

	# histo = fig.add_subplot(gs[4:,4:])
	# histo.yaxis.set_visible(False)
	# histo.hist(resi,10,range=(-2,2), orientation='horizontal')

	# figI = plt.figure()
	# axI = figI.add_subplot(111)
	# axI.errorbar(gal.R[-15:], gal.I[-15:], yerr=gal.W[-15:], fmt='b.')


	plt.show()