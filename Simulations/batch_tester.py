from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import lmfit as lm
import Simulate as S
import storeData as SD
import time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.misc import factorial
import csv
from scipy.special import gamma

def copy_params(parameters, trans_vary=True):
	new = lm.Parameters()
	for p in parameters.values():
		new.add(p.name, value=p.value, min=p.min, max=p.max)
		if trans_vary:
			new[p.name].vary=p.vary
	return new

def guess_params(true_params, percent):
	"""gets guess parameters differed from true ones by a percent, nD is fixed at 1"""
	P = copy_params(true_params, False)
	for i in ['Me', 'Re', 'n']:
		bulge = i+'B'
		disc = i+'D'
		P[bulge].value *= (1.- percent)
		P[disc].value *= (1.+ percent)
	P['nD'].value = 1.
	P['nD'].vary = False
	return P

def tester(param_template, bulge_params, R, zp, noise, weights, fit_range):
	"""returns fit_data and the initial guesses"""
	P = copy_params(param_template, False)
	P['MeB'].value = bulge_params['MeB'] # ready parameter set
	P['ReB'].value = bulge_params['ReB']
	P['nB'].value = bulge_params['nB']
	test_gal = S.sersic2(P, R, zp, False)
	test_gal += noise
	guesses = guess_params(P, 0.1)
	return S.fit(copy_params(guesses, True), S.sersic2, R, zp, test_gal, weights, fit_range, redchi_marker=None), guesses, test_gal

def preview(fit_data, fitted_profile, x, zp):
	"""shows a graph of fitting"""
	model_x = np.linspace(x[0], x[-1], 500)
	model, bulge, disc = S.sersic2(fit_data.params, model_x, zp, True)
	fig = plt.figure()
	gs = gridspec.GridSpec(6,1)
	ax = fig.add_subplot(gs[:4,0])
	res = fig.add_subplot(gs[4:,0])
	ax.plot(x, S.convert_I(fitted_profile, zp), 'b.')
	ax.plot(model_x, S.convert_I(model, zp), 'k-')
	ax.plot(model_x, S.convert_I(bulge, zp), 'g:')
	ax.plot(model_x, S.convert_I(disc, zp), 'r--')
	ax.set_ylim([35, 15])
	res_model = S.convert_I(S.sersic2(fit_data.params, x, zp, False), zp)
	res.plot(x, S.convert_I(fitted_profile, zp) - res_model, 'b.')
	plt.show()

def batch(filename, param_template, Me_range, Re_range, n_range, R, zp, noise, weights=None, fit_range=None):
	num = float(len(Me_range) * len(Re_range) * len(n_range))
	print num
	with open(filename, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		writer.writerow(['MeB_initial', 'ReB_initial', 'nB_initial', 'MeB_guess', 'ReB_guess', 'nB_guess', 'MeD_guess', 'ReD_guess'\
			, 'MeB_final', 'ReB_final', 'nB_final', 'MeD_final', 'ReD_final' \
			, 'MeB_stderr', 'ReB_stderr', 'nB_stderr', 'MeD_stderr', 'ReD_stderr', 'redchi'])
		count = 0
		for indme, Me in enumerate(Me_range):
			for indre, Re in enumerate(Re_range):
				for indn, n in enumerate(n_range):
					out, guess_p, test_gal = tester(param_template, {'MeB':Me, 'ReB':Re, 'nB':n}, R, zp, noise, weights, fit_range)
					# lm.report_fit(out.params)
					initials = [Me, Re, n]# + [param_template[i].value for i in ['MeD', 'ReD', 'nD']]
					guesses = [guess_p[i].value for i in ['MeB', 'ReB', 'nB', 'MeD', 'ReD']]
					finals = [out.params[i].value for i in ['MeB', 'ReB', 'nB', 'MeD', 'ReD']]
					std_errs = [out.params[i].stderr for i in ['MeB', 'ReB', 'nB', 'MeD', 'ReD']]
					redchi = out.redchi
					writer.writerow(initials+guesses+finals+std_errs+[redchi])
					perc = (100-(100*count/num))
					if (count%20) ==0:
						print '%.2f%% left' %  perc
					count += 1


					
def read_retrieve(filename):
	with open(filename, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		i = 0
		data = []
		for row in reader:
			if i == 0:
				headers = row
			else:
				data.append(row)
			i += 1
	return np.array(data).astype(float), headers

def get_originals(data_array_row, param_template):
	P_initial = copy_params(param_template, False)
	P_initial['MeB'].value = float(data_array_row[0]) # ready parameter set
	P_initial['ReB'].value = float(data_array_row[1])
	P_initial['nB'].value = float(data_array_row[2])

	P_fitted = copy_params(param_template, False)
	P_guess = copy_params(param_template, False)
	for i, name in enumerate(['MeB', 'ReB', 'nB', 'MeD', 'ReD']):
		P_fitted[name].value = data_array_row[8+i]
		P_guess[name].value = data_array_row[3+i]
	return P_initial, P_guess, P_fitted 

def redchi2(model, data, x, cutoff, weights=1, nvary=5):
	cutoff_R = SD.translate_x(x, cutoff)
	chi = (data[:cutoff_R] - model[:cutoff_R]) / weights
	return np.nansum(chi**2.) / (len(x[:cutoff_R]) - nvary)

def compare_residuals(truth, fitted, cutoff, x):
	"""compares residuals of truth+noise and the fitted model"""
	#truth_residual = truth + noise - truth = noise
	#fit_residual = truth + noise - fitted
	#diff = fit_residual - truth_residual
	cutoff_R = SD.translate_x(x, cutoff)
	return np.nansum((truth[:cutoff_R] - fitted[:cutoff_R])**2.)


def rate_fit(data_array_row, param_template, noise, cutoff, x):
	"""returns redchi2 for region before cutoff, truth-model, delta_params"""
	P_truth, P_guess, P_fit = get_originals(data_array_row, param_template)
	truth = S.sersic2(P_truth, R, 30., False)
	fitted = S.sersic2(P_fit, R, 30., False)
	noisy = truth + noise

	chi = redchi2(fitted, noisy, x, cutoff)
	diff = compare_residuals(truth, fitted, cutoff, x)

	delta_params = 0
	for p in P_truth.values():
		delta_params += (p.value - P_fit[p.name].value)
	return chi, diff, delta_params

def BD_ratio(params, zp):
	n = params['nB']
	bnb = S.get_b_n(n)
	bnd = S.get_b_n(params['nD'])
	pre = n * gamma(2*n) * np.exp(bnb) / (bnb ** (2 * n))
	h = params['ReD'] / ((bnd) ** params['nD'])
	Ie = 10 ** ((params['MeD'] - zp) / -2.5)
	I0 = Ie * np.exp(bnd)
	return pre * ((params['ReB'] / h)**2.) * (Ie / I0)


def extract_non_varying(data_array, array_headers, vary_dict):
	"""extracts only rows which have the values in varying"""
	h_index = [(array_headers.index(i), v) for i, v in vary_dict.iteritems()]
	conds = None 
	for i,v in h_index:
		if conds is None:
			conds = (data_array[:,i].astype(float) == float(v))
		else: 
			conds = (conds) & (data_array[:,i].astype(float) == float(v))
	return data_array[conds]

def add_data(array, headers, cutoff, zp, noise, x, param_template):
	chi_array = np.zeros((len(array), 3))
	BD_ratio_list = []
	for i, row in enumerate(array):
		chi_array[i,0], chi_array[i,1], chi_array[i,2] = rate_fit(row, param_template, noise, cutoff, x)
		initial, guess, final = get_originals(row, param_template)
		BD_ratio_list.append(BD_ratio(initial, zp))
	array = np.hstack((array, chi_array))
	array = np.hstack((array, np.array(BD_ratio_list).reshape(len(BD_ratio_list),1)))
	return array, headers+['cut_redchi', 'cut_diff', 'delta_p', 'BD_ratio']

def analysis_array(array, headers, cutoff, zp, noise, x, param_template):
	names = ['BD_ratio', 'BT_ratio', 'delta_params', 'params_ratio', 'cut_redchi', 'cut_diff']
	output = []
	for row in array:
		P_truth, P_guess, P_fit = get_originals(row, param_template)
		BD = BD_ratio(P_truth, zp)
		BT = 1. / (1 + ((1./BD)))
		deltas = {}
		for i in ['MeB', 'ReB', 'nB', 'MeD', 'ReD']:
			deltas['delta_'+i] = P_truth[i] - P_fit[i]
		Me_ratio = P_truth['MeB']/P_truth['MeD']
		Re_ratio = P_truth['ReB']/P_truth['ReD']
		total_t, bulge_t, disc_t = S.sersic2(P_truth, x, zp, True)
		total_f, bulge_f, disc_f = S.sersic2(P_fit, x, zp, True) 
		residual_excl =  total_t - total_f
		cutoff_ind = SD.translate_x(x, cutoff)
		redchi_excl = np.sum(residual_excl[:cutoff_ind]**2.) / (cutoff_ind - 5)
		cut_points_fit = len(np.where(np.diff(np.sign(bulge_f - disc_f)))[0]) - len(np.where((bulge_f - disc_f) == 0)[0])
		cut_points_initial = len(np.where(np.diff(np.sign(bulge_t - disc_t)))[0]) - len(np.where((bulge_t - disc_t) == 0)[0])

		d = {'BD_ratio':BD, 'BT_ratio':BT, 'Me_ratio':Me_ratio, 'Re_ratio':Re_ratio, 'redchi2_cut': redchi_excl, 'cross_initial': cut_points_initial,
		'cross_final': cut_points_fit}
		d.update(deltas)
		output.append(d)
	return output

def write_dict_to_csv(d_list, filename):
	w = csv.writer(open(filename, "wb"), delimiter=',')
	heads = [n for n,v in d_list[0].items()]
	w.writerow(heads)
	for i, d in enumerate(d_list):
		w.writerow([d[name] for name in heads])
	return
 
def classify_fit(fit_data_array, param_template, analysis_array, analysis_headers):
	classified = np.empty((len(fit_data_array),), dtype=int)
	
	stderr = fit_data_array[:,13:18]
	fitted = fit_data_array[:,8:13]
	delta_params = np.zeros((len(fit_data_array), 5))
	for i, name in enumerate(['MeB', 'ReB', 'nB', 'MeD', 'ReD']):
		delta_params[:, i] = analysis_array[:,analysis_headers.index('delta_'+name)]

	truth = np.zeros((len(fit_data_array), 5)) 
	truth[:,:3] = fit_data_array[:,:3]
	disc = [[param_template[i].value for i in ['MeD', 'ReD']]]*len(fit_data_array)
	truth[:,3:] = np.array(disc)
	bool_fit = (delta_params < fitted*0.01)# | (delta_params < 0.01) # 5 by N array

	in_range = np.zeros((len(fit_data_array)))
	for i, row in enumerate(bool_fit):
		in_range[i] = row.all()

	all_below = np.zeros((len(fit_data_array)))
	for i, row in enumerate(stderr):
		all_below[i] = (row <= (truth[i,:]*0.01)).all()

	for i in range(len(fit_data_array)): 
		if all_below[i] and in_range[i]:# good match 
			classified[i] = 1
		elif all_below[i] or in_range[i]: # disputed
			classified[i] = 0
		else:
			classified[i] = -1 # bad
	return classified

def histogram_fit(criterion, data_array, analysis_array, data_headers, analysis_headers, classified_array, bins, sep=1.0, 
	normalise=False, density=False, axis=None):

	

	if criterion in data_headers:
		plot_array = data_array[:, data_headers.index(criterion)]
	elif criterion in analysis_headers:
		plot_array = analysis_array[:, analysis_headers.index(criterion)]
	else:
		raise ValueError("Header not found")
	mask = data_array[:, data_headers.index('BT_ratio')] == 0.5
	plot_array = plot_array[mask]

	tot = len(plot_array)
	total, binning = np.histogram(plot_array, bins, density=density)
	good = plot_array[(classified_array[mask]==1)] 
	disp = plot_array[(classified_array[mask]==0)]
	bad = plot_array[(classified_array[mask]==-1)]

	if axis is None:
			fig = plt.figure()
			axis = fig.add_subplot(111)
	colour_list = ['g', 'y', 'r'][::-1]
	labs = ['Bad', 'Disputed', 'Good']
	cumul = np.zeros(total.shape)
	centres = (binning[1:] + binning[:-1])/2
	for ind, i in enumerate([bad, disp, good]):
		H, none = np.histogram(i, binning, density=density)
		H = H.astype(float)
		width = sep*(binning[1] - binning[0])
		center = (binning[:-1] + binning[1:]) / 2
		if normalise: H=100*H/total
		axis.bar(center, H, bottom=cumul, align='center', width=width, color=colour_list[ind], edgecolor=colour_list[ind], label=labs[ind])
		cumul += H
	if normalise:
		axis.set_ylabel('Percentage of galaxies')
		axis.set_ylim((0,100))
	else: 
		axis.set_ylabel('Number of galaxies')
	axis.set_xlim((binning[0], binning[-1]))
	return axis

if __name__ == '__main__':
	temp_gal = SD.Galaxy('temp')
	temp_gal.import_file('C:\\Users\\User\\code\\Python\\Sersic_Plots\\mega_1237665428090192022_cts.ascii', 'mega')
	R = temp_gal[0][0].R

	pars = lm.Parameters()
	pars.add_many(('MeB', None, True, 10.), ('ReB', None, True, 0.01), ('nB', None, True, 0.01),\
		('MeD', 21., True, 10.), ('ReD', 9., True, 0.01), ('nD', 1., False, 0.01))
	np.random.seed(42)
	noise_profile = np.random.random(R.shape)*30.

	# bulge_params = {'nB': 4., 'ReB':0.3, 'MeB':18.}

	# out, guess, gal = tester(pars, bulge_params, R, 30., noise_profile, None, None)
	# print out.params['ReD'].stderr
	# lm.report_fit(out.params, show_correl=False)
	# preview(out, gal, R, 30.)

	# Mes = [18.]
	# Res = np.arange(0.1, 12., 0.1)
	# ns = np.arange(0.5, 6, 0.1)

	# start = time.clock()
	# print 'running...'
	# batch('fit_data_staticMe.csv', pars, Mes, Res, ns, R, 30., noise_profile)
	# print 'took', time.clock() - start, 'seconds'

	array, headers = read_retrieve('combined.csv')
	print 'retrieved!'
	# save_dict = analysis_array(array, headers, 30., 30., noise_profile, R, pars)
	# print 'analysed!'
	# write_dict_to_csv(save_dict, 'combined_analysis_data.csv')
	
	ana_array, ana_headers = read_retrieve('combined_analysis_data.csv')
	print 'retrieved!'
	c = classify_fit(array, pars, ana_array, ana_headers)
	print 'classified!'
	print np.sum(c==1), np.sum(c==0), np.sum(c==-1)

	fig = plt.figure()
	fig.set_facecolor('white')
	ax1 = fig.add_subplot(221)
	ax2 = fig.add_subplot(222)
	ax3 = fig.add_subplot(223)
	ax4 = fig.add_subplot(224)
	
	histogram_fit('Re_ratio', array, ana_array, headers, ana_headers, c, 99, axis=ax1, normalise=True)
	ax1.set_xlabel('$Re_B/Re_D$')
	histogram_fit('Me_ratio', array, ana_array, headers, ana_headers, c, 50, axis=ax2, normalise=True)
	ax2.set_xlabel('$Me_B/Me_D$')
	histogram_fit('nB_initial', array, ana_array, headers, ana_headers, c, np.arange(0.5, 6, 0.1), axis=ax3, normalise=True)
	ax3.set_xlabel('$n_B/n_D$')
	histogram_fit('BT_ratio', array, ana_array, headers, ana_headers, c, 30, axis=ax4, normalise=True)
	ax4.set_xlabel('$B/T$')
	fig.suptitle('Goodness of fit versus changes in parameters (B/T frozen at 0.5)')
	for i in [ax1, ax2, ax3, ax4]:
		i.xaxis.labelpad = 0.1
	# plt.tight_layout()
	ax1.axvline(x=0.06, color='k', linestyle='--'); ax1.axvline(x=0.5, color='k', linestyle='--')
	ax2.axvline(x=0.9, color='k', linestyle='--')
	ax3.axvline(x=2.1, color='k', linestyle='--')
	ax4.axvline(x=0.008, color='k', linestyle='--'); ax4.axvline(x=0.17, color='k', linestyle='--'); ax4.axvline(x=0.61, color='k', linestyle='--')

	plt.show()

