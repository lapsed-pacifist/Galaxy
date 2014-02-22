import numpy as np
import lmfit as lm
import storeData as SD
import matplotlib.pyplot as plt
import fitting as F
import matplotlib.gridspec as gridspec
import sigfig as sf
from matplotlib import rc
from matplotlib.ticker import MaxNLocator, AutoMinorLocator 
from copy import deepcopy
from matplotlib.backends.backend_pdf import PdfPages
import sys
rc('font',**{'family':'serif','serif':['Times New Roman']})

params={'axes.labelsize':12,'xtick.labelsize':12,'ytick.labelsize':12,'figure.figsize':(8,6), 'figure.facecolor':'white'}
plt.rcParams.update(params)

def find_offaxis(Y, limits):
	nan_locator = [i for i, v in enumerate(Y) if np.isnan(v)]
	top = [i for i, v in enumerate(Y) if v > limits[1]]
	bottom = [i for i, v in enumerate(Y) if v < limits[0]]
	return top, bottom, nan_locator

def mark_offaxis(X, Y, limits, axis, nan_limit='up', color='k', len_scale=0.2):
	top, bottom, nans = find_offaxis(Y, limits)
	if nan_limit == 'up':
		top += nans
	elif nan_limit == 'down':
		bottom += nans
	length = (limits[1] - limits[0]) * len_scale	
	props = dict(arrowstyle="->", linewidth=2, color=color)
	for i in top:
		axis.annotate("", xy=(X[i], limits[1]), xytext=(X[i], limits[1]-length), arrowprops=props)
	for i in bottom:
		axis.annotate("", xy=(X[i], limits[0]), xytext=(X[i], limits[0]+length), arrowprops=props)

def extract_ci(ci, interval):
	zero = int((len(ci['nB'])/2.) + 0.5)
	conf = {}
	for key, v in ci.iteritems():
		up = v[zero + interval - 1][1]
		down = v[zero - interval - 1][1]
		conf[key] = [up, down]
	return conf

def plot_ci(profile, fit_data, ci, interval, axis):
	conf = extract_ci(ci, interval)
	Pup = deepcopy(fit_data.params)
	Pdown = deepcopy(fit_data.params)
	r = np.linspace(0, axis.get_xlim()[1], 1000)
	for key, v in conf.iteritems():
		if 'Me' not in key:
			Pdown[key].value = v[1]
			Pup[key].value = v[0]
		else:
			Pdown[key].value = v[0]
			Pup[key].value = v[1]
	up = F.convert_I(F.sersic2(Pup, r, profile.zeropoint, 'T'), profile.zeropoint)
	down = F.convert_I(F.sersic2(Pdown, r, profile.zeropoint, 'T'), profile.zeropoint)
	axis.fill_between(r, up, down, alpha=0.2)

def plot_ci_residual(profile, fit_data, ci, interval, axis):
	conf = extract_ci(ci, interval)
	Pup = deepcopy(fit_data.params)
	Pdown = deepcopy(fit_data.params)
	r = np.linspace(0, axis.get_xlim()[1], 1000)
	for key, v in conf.iteritems():
		if 'Me' not in key:
			Pdown[key].value = v[1]
			Pup[key].value = v[0]
		else:
			Pdown[key].value = v[0]
			Pup[key].value = v[1]
	total = F.convert_I(F.sersic2(fit_data.params, r, profile.zeropoint, 'T'), profile.zeropoint)
	up = total - F.convert_I(F.sersic2(Pup, r, profile.zeropoint, 'T'), profile.zeropoint)
	down = total - F.convert_I(F.sersic2(Pdown, r, profile.zeropoint, 'T'), profile.zeropoint) 
	axis.fill_between(r, up, down, alpha=0.1)

def add_sky_line(profile, axis):
	lims = axis.get_xlim()
	sky_error = F.convert_I(np.sqrt(profile.sky_average), profile.zeropoint)
	props = dict(arrowstyle="->", linewidth=2, color='k')
	axis.annotate("sky", xy=(0, sky_error), xytext=((lims[1]-lims[0])*0.05, sky_error), arrowprops=props, va='center')

def add_param_text(profile, fit_data, axis, sigs=1, perc=True, name_list=None, conf=None, align='right'):
	text = ''
	if name_list is None: 
		name_list = ['MeB', 'ReB', 'nB', 'MeD', 'ReD', 'nD']
	for key in name_list:
		name = list(key)
		if 'M' in name:
			name[name.index('M')] = '\\mu '
		name.insert(-1, '_')
		name = ''.join(name)
		if conf:
			error = conf[key]
		else:
			error = [fit_data.params[key].stderr]
		if not any(error): # if all 0
			text += '$%s = %.1f \\rm{ (fixed)}$\n' % (name, fit_data.params[key].value)
		else:
			val = fit_data.params[key].value
			errup = fit_data.params[key].stderr
			errdown = fit_data.params[key].stderr
			if perc:
				errup = errup * 100 / val
				errdown = errdown * 100 /val
				string = ['%.1f' % val, '%.1f' % errup, '%.1f' % errdown]
			else: string = sf.round_sig_error2(val, errup, errdown, sigs)
			if conf:
				if perc:
					text += '$%s = %s ^{+%s} _{-%s}\\%%$\n' % (name, string[0], string[1], string[2])
				else:
					text += '$%s = %s ^{+%s} _{-%s}$\n' % (name, string[0], string[1], string[2])
			else:
				if perc:
					text += '$%s = %s \pm %s\\%%$\n' % (name, string[0], string[1])
				else:
					text += '$%s = %s \pm %s$\n' % (name, string[0], string[1])
	text += '$\chi_{\\nu}^{2} = %.2f$' % fit_data.redchi
	if align == 'right':
		axis.text(0.99, 0.99, text, transform=axis.transAxes, ha='right', va='top')
	else:
		axis.text(0.01, 0.99, text, transform=axis.transAxes, ha='left', va='top')


def plot_fit(profile, axis, fit_data, conf=None, text=True, fmt='b.', text_align='right'):
	axis.errorbar(profile.R, profile.M, yerr=profile.MW, fmt=fmt)
 	mod_R = np.linspace(0, profile.R[-1], 500)
	total, bulge, disc = F.sersic2(fit_data.params, mod_R, profile.zeropoint, 'TBD')
	axis.plot(mod_R, F.convert_I(total, profile.zeropoint), 'k-')
	axis.plot(mod_R, F.convert_I(bulge, profile.zeropoint), 'g:', linewidth=2.)
	axis.plot(mod_R, F.convert_I(disc, profile.zeropoint), 'r--', linewidth=2.)
	axis.yaxis.set_minor_locator(AutoMinorLocator())
	# add_sky_line(profile, axis)
	if text:
		add_param_text(profile, fit_data, axis, 1, None, conf, align=text_align)
	mark_offaxis(profile.R, profile.M, [-100,100], axis, 'up', 'r')
	

def plot_residuals(profile, axis, fit_data, norm=False, fmt='b.', limits=None):
	model = F.convert_I(F.sersic2(fit_data.params, profile.R, profile.zeropoint, 'T'), profile.zeropoint)
	residual = profile.M - model
	if norm:
		residual /= profile.MW[0]
	axis.errorbar(profile.R, residual, yerr=profile.MW, fmt=fmt)
	if limits is None:
		limits = [-2, 2]
		axis.set_ylim(limits)
	else: axis.set_ylim(limits)
	mark_offaxis(profile.R, residual, limits, axis)
	axis.axhline(y=0, linestyle=':')

def prepare_single_figure():
	figure = plt.figure()
	gs = gridspec.GridSpec(6,6)
	axis = figure.add_subplot(gs[0:4,0:], label='main') # main plot
	axis.set_ylim([35,15])
	residuals = figure.add_subplot(gs[4:,0:], label='res') #residuals
	figure.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in figure.axes[:-1]], visible=False)
	return figure, axis, residuals

def prepare_double_figure():
	figure = plt.figure()
	gs = gridspec.GridSpec(6,2)
	axes, residuals = [], []
	for i in range(2):
		axis = figure.add_subplot(gs[0:4,i], label='main') # main plot
		axis.set_ylim([35,15])
		residual = figure.add_subplot(gs[4:,i], label='res') #residuals
		figure.subplots_adjust(hspace=0)
		plt.setp(axis.get_xticklabels(), visible=False)
		axes.append(axis)
		residuals.append(residual)
	figure.subplots_adjust(wspace=0)
	return figure, axes, residuals

def prepare_whole_figure():
	figure = plt.figure()
	gs = gridspec.GridSpec(12,2)
	axes, residuals = [], []
	for i in range(2):
		for j in range(2,7,4):
			axis = figure.add_subplot(gs[j:j+4,i], label='main') # main plot
			figure.subplots_adjust(hspace=0)
			plt.setp(axis.get_xticklabels(), visible=False)
			axes.append(axis)
	for i in range(2):
		for j in [0,10]:
			residual = figure.add_subplot(gs[j:j+2,i], label='res') #residuals
			residuals.append(residual)
	figure.subplots_adjust(wspace=0)
	figure.subplots_adjust(hspace=0)
	return figure, axes, residuals

def graph_profile(profile, fit, sky_error):
	ci = lm.conf_interval(fit)
	profile.confs.update(ci)
	fig, main, res = prepare_single_figure()
	plot_fit(profile, main, fit, sky_error)
	plot_residuals(profile, res, fit, False)
	# plot_residuals(profile, res, fit, True, 'g.')
	plot_ci(profile, fit, ci, 2, main)
	plot_ci_residual(profile, fit, ci, 2, res)
	main.yaxis.set_major_locator(MaxNLocator(prune='upper'))
	main.set_ylabel('$\mu$ [mag arcsec$^{-1}]$')
	res.set_xlabel('$R$ [arcsec]')
	res.set_ylabel('$\Delta\mu$ [mag arcsec$^{-1}$]')
	title = ''
	if profile.gal_name: title += profile.gal_name
	if profile.cam_name: title += profile.cam_name
	if profile.name: title += profile.name
	main.set_title(title)
	return fig

def graph_camera(camera, fit_name):
	fig, mains, resids = prepare_double_figure()
	for i, main in enumerate(mains):
		p = camera[i]
		fit = p.fits[fit_name][0]
		plot_fit(p, main, fit, p.sky_fit)
		plot_residuals(p, resids[i], fit, False)
		main.yaxis.set_major_locator(MaxNLocator(prune='upper'))
		main.set_ylabel('$\mu$ [mag arcsec$^{-1}]$')
		resids[i].set_xlabel('$R$ [arcsec]')
		resids[i].set_ylabel('$\Delta\mu$ [mag arcsec$^{-1}$]')
		title = ''
		if p.name: title += p.name
		main.set_title(title)
	fig.suptitle(camera.gal_name+' '+camera.name)
	mains[1].yaxis.tick_right()
	resids[1].yaxis.tick_right()
	mains[1].yaxis.set_label_position("right")
	resids[1].yaxis.set_label_position("right")
	return fig

def graph_galaxy(galaxy, fit_name):
	fig, mains, resids = prepare_double_figure()
	for i, main in enumerate(mains):
		p1, p2 = galaxy[0][i], galaxy[1][i]
		fit1, fit2 = p1.fits[fit_name][0], p2.fits[fit_name][0]
		if i ==1: al ='right'
		else: al = 'left'
		plot_fit(p1, main, fit1, p1.sky_fit, text=True, fmt='b.', text_align=al)
		plot_residuals(p1, resids[i], fit1, False, fmt='b.')
		plot_fit(p2, main, fit2, p2.sky_fit, text=False, fmt='g.')
		plot_residuals(p2, resids[i], fit2, False, fmt='g.')
		main.yaxis.set_major_locator(MaxNLocator(prune='upper'))
		main.set_ylabel('$\mu$ [mag arcsec$^{-1}]$')
		resids[i].set_xlabel('$R$ [arcsec]')
		resids[i].set_ylabel('$\Delta\mu$ [mag arcsec$^{-1}$]')
		title = ''
		if p1.name: title += p1.name
		main.set_title(title)
	fig.suptitle(galaxy.name)
	mains[1].yaxis.tick_right()
	resids[1].yaxis.tick_right()
	mains[1].yaxis.set_label_position("right")
	resids[1].yaxis.set_label_position("right")
	mains[0].set_xlim(mains[0].get_xlim()[::-1])
	resids[0].set_xlim(resids[0].get_xlim()[::-1])
	resids[0].xaxis.set_label_coords(1.0, -0.1)
	resids[1].set_xlabel('')
	return fig


if __name__ == '__main__':
	# n = 'sdss_1237665427552993436'
	# direct = r'repository/'+n+'_cts.ascii'
	# G = SD.Galaxy('first')
	# G.import_file(direct)
	# p = G[0][1]
	# p.add_ini_params(r'repository/'+n+'.ini')
	# fitted = F.fit_sersic_exp(p, store=True, fit_range=[3., p.R[-1]])
	# sky = F.sky_error(p, p.sky_var, 'sersic+exp')
	# graph_profile(p, fitted, sky)
	direct = 'repository'
	Gal_list, Gal_names = SD.import_directory(direct)
	gal = Gal_list[12]
	T = len(Gal_list) * 4
	cumul = 0
	gal = Gal_list[27]
	for i in range(2):
		for j in range(2):
			fitted = F.fit_sersic_exp(gal[i][j], store=True, fit_range=[3., gal[i][j].R[-1]])
			sky = F.sky_error(gal[i][j], gal[i][j].sky_var, 'sersic+exp')
			sys.stdout.write('\r')
			sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*cumul/T), (100*cumul/T)))
			sys.stdout.flush()
			cumul += 1

	# T = len(Gal_list)
	# cumul = 0
	# for gal in Gal_list:
	# 	figure = graph_galaxy(gal, 'sersic+exp')
	# 	plt.savefig('report/%s.png' % gal.ID)
	# 	plt.close()
	# 	sys.stdout.write('\r')
	# 	sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*cumul/T), (100*cumul/T)))
	# 	sys.stdout.flush()
	# 	cumul +=1

	plt.savefig('report/%s.png' % gal.ID)
	plt.close()