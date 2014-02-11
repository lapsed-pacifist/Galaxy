"""
plot components separately 
errorbars
residuals
residual errorbars
NaN/out of range marker arrows
confidence regions
fit region markers
sky region markers
sky level line
Re markers
text for fit data
"""
import storeData as SD
import lmfit as lm
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from copy import deepcopy
import fitGalaxy as F
import sigfig as sf
from matplotlib import rc
import basic_testing as bb
import convert_to_deluxe as CD

rc('font',**{'family':'serif','serif':['Times New Roman']})

params={'axes.labelsize':15,'xtick.labelsize':15,'ytick.labelsize':15,'figure.figsize':(8,6)}
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
	
def graph_models(profile, fit_data, R_model=None):
	if R_model is None:
		R_model = np.linspace(profile.R[0], profile.R[-1], 1000)
	P = fit_data.params
	total, bulge, disc = bb.sersic2(P, R_model, profile.zeropoint, True)
	return bb.convert_I(disc, profile.zeropoint), bb.convert_I(bulge, profile.zeropoint), bb.convert_I(total, profile.zeropoint), R_model

def plot_fit(profile, models, axis):
	axis.errorbar(profile.R, profile.M, yerr=profile.MW, fmt='b.')
	maxi = profile.M + profile.MW[1]
	mini = profile.M - profile.MW[0]
	lims = [min(mini)-0.5, max(maxi)+0.5]
	disc, bulge, total, R_model = models
	axis.plot(R_model, total, 'k-', linewidth=2, label='Total')
	axis.plot(R_model, disc, 'r--', linewidth=2, label='Disc')
	axis.plot(R_model, bulge, 'g:', linewidth=2, label='Bulge')
	mark_offaxis(profile.R, profile.M, lims, axis, 'up', 'r', 0.1)
	axis.set_title(profile.gal_name+profile.cam_name+profile.name)
	axis.set_ylim(lims[::-1])
	lg = axis.legend(loc=3, fancybox=True, borderaxespad=0)
	lg.draw_frame(False)
	add_fit_range(profile, axis)
	axis.set_xlim([0, profile.R[-1]+1])
	# ticks = np.arange(int(lims[1]), int(lims[0]), -5)
	# axis.set_yticks(ticks)

def plot_residuals(profile, fit_data, axis):
	disc, bulge, total, R_model = graph_models(profile, fit_data, profile.R)
	residual = total - profile.M
	axis.errorbar(profile.R, residual, yerr=profile.MW, fmt='b.')
	lims = [-1.5, 1.5]
	axis.set_ylim(lims)
	mark_offaxis(profile.R, residual, lims, axis, 'none')
	axis.set_yticks(np.arange(lims[0],lims[1]+1))
	axis.axhline(y=0, linestyle=':')

def plot_ci(profile, fit_data, ci, axis):
	Pup = deepcopy(fit_data.params)
	Pdown = deepcopy(fit_data.params)
	r = np.linspace(0, profile.R[-1], 1000)
	for key, v in ci.iteritems():
		if 'Me' not in key:
			Pdown[key].value += v[1]
			Pup[key].value -= v[0]
		else:
			Pdown[key].value -= v[0]
			Pup[key].value += v[1]
	up = bb.convert_I(bb.sersic2(Pup, r, profile.zeropoint, False), profile.zeropoint)
	down = bb.convert_I(bb.sersic2(Pdown, r, profile.zeropoint, False), profile.zeropoint)
	axis.fill_between(r, up, down, alpha=0.2)

def plot_ci_residual(profile, fit_data, ci, axis):
	Pup = deepcopy(fit_data.params)
	Pdown = deepcopy(fit_data.params)
	r = np.linspace(0, profile.R[-1], 1000)
	for key, v in ci.iteritems():
		if 'Me' not in key:
			Pdown[key].value += v[1]
			Pup[key].value -= v[0]
		else:
			Pdown[key].value -= v[0]
			Pup[key].value += v[1]
	total = bb.convert_I(bb.sersic2(fit_data.params, r, profile.zeropoint, False), profile.zeropoint)
	up = total - bb.convert_I(bb.sersic2(Pup, r, profile.zeropoint, False), profile.zeropoint)
	down = total - bb.convert_I(bb.sersic2(Pdown, r, profile.zeropoint, False), profile.zeropoint) 
	axis.fill_between(r, up, down, alpha=0.1)

def add_param_text(profile, fit_data, axis, sigs=1, perc=True, variables=None, conf=None):
	text = ''
	if variables is None: 
		variables = ['MeB', 'ReB', 'nB', 'MeD', 'ReD', 'nD']
	for key in variables:
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
			errup = conf[key][1]
			errdown = conf[key][0]
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
	text += '$\chi_{\\nu}^{2} = %.2f$' % bb.find_redchi(profile, fit_data)
	axis.text(0.99, 0.99, text, transform=axis.transAxes, ha='right', va='top')

def add_sky_line(profile, axis):
	lims = axis.get_xlim()
	sky_error = F.convert_I(np.sqrt(profile.sky_average), profile.zeropoint)
	props = dict(arrowstyle="->", linewidth=2, color='k')
	axis.annotate("sky", xy=(0, sky_error), xytext=((lims[1]-lims[0])*0.05, sky_error), arrowprops=props, va='center')

def add_Re_marker(profile, fit_data, axis):
	p = deepcopy(fit_data.params)
	rd = SD.translate_x(profile.R, p['ReD'].value)
	rb = SD.translate_x(profile.R, p['ReB'].value)
	axis.plot(p['ReD'].value, profile.M[rd], 'ko', markersize=10, alpha=0.5)
	axis.plot(p['ReB'].value, profile.M[rb], 'ko', markersize=10, alpha=0.5)

def add_fit_range(profile, axis):
	axis.axvline(x=profile.fit_range[0], linestyle='--')
	axis.axvline(x=profile.fit_range[1], linestyle='--')

def average_params(*params):
	P = lm.Parameters()
	for p in params[0].values():
		P.add(p.name, p.value, p.vary, min=p.min, max=p.max)
	for p in params[1:]:
		for i in p.values():
			P[i.name].value += i.value
			P[i.name].value /= 2
	return P


def plot_profile(profile, fit_name, figure, grid=None, side=None):
	# conf = profile.confs[fit_name]
	if grid is None:
		grid = gridspec.GridSpec(8,1)
	if side is None:
		ax = figure.add_subplot(grid[0:6])
		res = figure.add_subplot(grid[6:], sharex=ax)
	else:
		ax = figure.add_subplot(grid[0:6,side])
		res = figure.add_subplot(grid[6:,side], sharex=ax)
	fit = profile.fits[fit_name][0] #retrive fit_data
	mods = graph_models(profile, fit) # prepare models
	plot_fit(profile, mods, ax) #plot models and data
	plot_residuals(profile, fit, res) # plot residuals

	ci = profile.confs[fit_name]
	plot_ci(profile, fit, ci, ax) #plot ci area
	plot_ci_residual(profile, fit, ci, res) #plot ci residual area

	ax.xaxis.set_visible(False) # join axes
	figure.subplots_adjust(hspace=0)
	add_param_text(profile, fit, ax, conf=ci)
	add_sky_line(profile, ax)
	add_Re_marker(profile, fit, ax)
	return figure, ax, res

def plot_compare_profile(profile1, profile2, fit_name1, fit_name2, figure):
	gs = gridspec.GridSpec(8,2)
	figure, ax1, res1 = plot_profile(profile1, fit_name1, figure, gs, side=0)
	figure, ax2, res2 = plot_profile(profile2, fit_name2, figure, gs, side=1)
	ax2.yaxis.set_visible(False)
	res2.yaxis.set_visible(False)
	xticks = res2.xaxis.get_major_ticks()
	xticks[0].label1.set_visible(False)
	# ax2.set_xticklabels([1,2,3,4])
	figure.subplots_adjust(wspace=0)
	ax1.set_ylabel('$\mu / \\rm{mag\/arcsec}^{-1}$')
	res1.set_xlabel('R / arcsec')
	res2.set_xlabel('R / arcsec')
	ax2.set_ylim(ax1.get_ylim())
	return figure, ax1, res1, ax2, res2

def output_table(galaxy, fit_name, direct=None):
	T = CD.Table(galaxy.name, 'l|ccccccc')
	f = lambda x: str(x)
	var_names = ['MeB','ReB', 'nB', 'MeD', 'ReD','nD']
	h = map(f, [i for i in var_names]) #prepare headers
	col_h = ['Mega 1', 'Mega 2', 'SDSS 1', 'SDSS 2'] #first column
	h.insert(0, '')
	chi_col = []
	for c in ['mega', 'sdss']:
		camera = galaxy[c]
		for profile in camera:
			chi_col.append(bb.find_redchi(profile, profile.fits[fit_name][0]))
			values = []; errorup = []; errordown = []
			P = profile.fits[fit_name][0].params
			C = profile.confs[fit_name]
			S = profile.sky_fit[fit_name]
			for i in var_names:
				values.append(P[i].value)
				errorup.append(C[i][1])
				errordown.append(C[i][0])
			T.add_row(values, 1, errorup, errordown, True)
	T.add_header_column(col_h)
	T.add_header(h+ ['$\chi_{\\nu}^2$'])
	T.add_column(chi_col, 2)
	if direct:
		T.print_data(direct+T.name+'.tex')
	else:
		T.print_data('mydata.tex')

if __name__ == '__main__':
	direct = 'C:\\Users\\User\\code\\Python\\Sersic_Plots'
	gal_list, gal_names = SD.import_directory(direct)
	for i in range(len(gal_list)):
		for j in range(len(gal_list[i])):
			for g in gal_list[i][j]:
				print g.gal_name + g.cam_name + g.name
				out = bb.fit_sersic_exp(g, fit_range=None, store=True)
				bb.sky_error_se(g, 10., 'sersic+exp')
				bb.find_confs_se(g, 'sersic+exp', 0.95)

		prof1 = gal_list[i][0][0]
		prof2 = gal_list[i][0][1]
		prof3 = gal_list[i][1][0]
		prof4 = gal_list[i][1][1]
		
		d = 'C:\\Users\\User\\Documents\\Project work\\Latex Report\\Prelim\\Figs\\pre_plot\\'
		fig = plt.figure()
		fig.set_facecolor('white')
		plot_compare_profile(prof1, prof4, 'sersic+exp', 'sersic+exp', fig)
		# output_table(gal_list[i], 'sersic+exp', d)
		fig.set_size_inches(18.5,10.5)
		# plt.savefig(d+'%s.png' % (prof1.gal_name), dpi=300)
		plt.show()
	
