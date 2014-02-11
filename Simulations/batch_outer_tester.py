# import Queue
# import threading
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import gammainc
from scipy.special import gamma
import copy
import lmfit as lm
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.widgets import Slider
from storeData import translate_x
import Simulate as S

def get_b_n(m):
	#   find sersic b_n coefficient
	#   this is more accurate than the usual b_n = 1.9992*n_sersic - 0.3271
	if m == 0.0: return -0.3271
	b_n = 2.0*m - 1.0/3.0 + 4.0/m/405.0 + 46.0/m/m/25515.0 \
				+ 131.0/m/m/m/1148175.0  \
				- 2194697.0/m/m/m/m/30690717750.0
	return b_n

def find_outer_BD(params, cutoff, end=None):
	"""cutoff is in scale length multiples"""
	bnd = 1.67838865492
	cutoff_R = cutoff * params['ReD'].value / bnd 

	Ie_ratio = 10 ** ((params['MeD'].value - params['MeB'].value) / 2.5)
	Re_ratio = (params['ReB'].value / params['ReD'].value) ** 2.
	nB = params['nB'].value
	bnb = get_b_n(nB) 
	exp_ratio = np.exp(bnb - bnd)
	bn_ratio = (bnd*bnd / (bnb ** (2. * nB)))
	xb = bnb * ((cutoff_R / params['ReB'].value) ** (1./nB))
	xd = bnd * (cutoff_R / params['ReD'].value) 
	Bgamma = gamma(2.*nB)
	if end is None:
		gammas_B = (Bgamma - (gammainc(2. * nB, xb) * Bgamma))
		gammas_D = (1. - (gammainc(2., xd) * 1.))
	else:
		xb_end = bnb * ((end / params['ReB'].value) ** (1./nB))
		xd_end = bnd * (end / params['ReD'].value) 
		gammas_B = ((gammainc(2. * nB, xb_end) * Bgamma) - (gammainc(2. * nB, xb) * Bgamma))
		gammas_D = ((gammainc(2., xd_end) * 1.) - (gammainc(2., xd) * 1.))
	return Ie_ratio * Re_ratio * nB * exp_ratio * bn_ratio * gammas_B / gammas_D

def apply_func(M, R, N, pars, cutoff_R, end=None):
	new = copy.deepcopy(pars)
	new['MeB'].value = M; new['ReB'].value =R; new['nB'].value =N
	return find_outer_BD(new, cutoff_R, end)

def q_func(q, M, R, N, P, cutoff):
	q.put(apply_func(M, R, N, P, cutoff))

def apply_plot(M, R, N, pars, x):
	p = S.copy_params(pars)
	p['MeB'].value = M
	p['ReB'].value = R
	p['nB'].value = N
	return S.sersic2(p, x, 30., True)

if __name__ == '__main__':
	P = lm.Parameters()
	P.add_many(('MeB', None, True, 10.), ('ReB', None, True, 0.01), ('nB', None, True, 0.01),\
		('MeD', 21., True, 10.), ('ReD', 9., True, 0.01), ('nD', 1., False, 0.01))
	MeB = np.linspace(10,30, 49)
	ReB = np.linspace(0.1, P['ReD'].value*1.9, 50)
	nB = np.linspace(0.5, 5., 51)
	# nB = [4.]
	start_cut = 3.; end_cut = 50.
	BD_array = np.zeros((len(MeB),len(ReB),len(nB)))

	# print 'start!'
	# bar_string = ' '*100
	# T = len(MeB) * len(ReB) * len(nB)
	# cumul = 0
	# for i, m in enumerate(MeB):
	# 	for j, r in enumerate(ReB):
	# 		for k, n in enumerate(nB):
	# 			cumul += 1
	# 			sys.stdout.write('\r')
	# 			sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*cumul/T), (100*cumul/T)))
	# 			# sys.stdout.flush()
	# 			BD_array[i,j,k] = apply_func(m,r,n,P,start_cut,end_cut)
	# print '\nend!'

	# BT_array = 1. / (1. + ((1./BD_array)))

	# np.save('BT_file', BT_array)
	# print 'saved!'
	BT_array = np.load('BT_file.npy')
	# BT_array[BT_array < 0.] = np.nan
	# BD_array = 1. / (((1./BT_array)) - 1.)


	fig = plt.figure()
	
	ax = fig.add_subplot(111)
	plt.subplots_adjust(left=0.25, bottom=0.3)


	# frame = 5
	# X, Y = np.meshgrid(nB, ReB)
	# Z = BT_array[frame,:,:]/np.max(BT_array[frame,:,:])
	# l = ax.contourf(X, Y, Z)
	# cb = fig.colorbar(l,ax=ax)

	# axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
	# sframe = Slider(axframe, 'MeB', MeB[0], MeB[-1], valinit=MeB[0], valfmt='%.1f')

	# def update(val):
	# 	f = translate_x(MeB, sframe.val)
	# 	ax.cla()
	# 	l = ax.contourf(X, Y,  BT_array[f,:,:])
	# 	# cb.on_mappable_changed(l)
	# 	ax.set_xlabel('nB')
	# 	ax.set_ylabel('ReB')

		
	# sframe.on_changed(update)
	xr = np.linspace(0.1,100,600)
	fig_gal = plt.figure(); axg = fig_gal.add_subplot(111)
	total1, bulge1, disc1 = apply_plot(MeB[0], ReB[0], nB[0], P, xr)
	t, = axg.plot(xr, S.convert_I(total1, 30), 'k-')
	b, = axg.plot(xr, S.convert_I(bulge1, 30), 'g:')
	d, = axg.plot(xr, S.convert_I(disc1, 30), 'r--')
	l, = ax.plot(nB, BT_array[0,0,:])
	axg.set_ylim([35,15])
	ax.axvline(x=nB[0], linestyle=':', color='r')
	axg.axvline(x=start_cut* P['ReD'].value / 1.67838865492, linestyle='--')
	ax.set_ylim([0,1])
	lim = 0.05
	limit_line = MeB[BT_array[:,0,0] < lim]
	if len(limit_line) > 0:
		ax.axvline(x=limit_line[-1], linestyle='--')
	axframe1 = fig.add_axes([0.25, 0.1, 0.65, 0.03])
	axframe2 = fig.add_axes([0.25, 0.15, 0.65, 0.03])
	sframe1 = Slider(axframe1, 'MeB', MeB[0], MeB[-1], valinit=MeB[0], valfmt='%.1f')
	sframe2 = Slider(axframe2, 'ReB', ReB[0], ReB[-1], valinit=ReB[0], valfmt='%.1f')
	axframe3 = fig.add_axes([0.25, 0.20, 0.65, 0.03])
	sframe3 = Slider(axframe3, 'nB', nB[0], nB[-1], valinit=nB[0], valfmt='%.1f')

	def update(val):
		val1 = translate_x(MeB, sframe1.val)
		val2 = translate_x(ReB, sframe2.val)
		val3 = translate_x(nB, sframe3.val)
		ax.cla(); axg.cla()
		ax.plot(nB, BT_array[val1, val2, :])
		limit_line = nB[BT_array[val1, val2, :] < lim]
		if len(limit_line) > 0:
			ax.axvline(x=limit_line[-1], linestyle='--')
		total1, bulge1, disc1 = apply_plot(MeB[val1], ReB[val2], nB[val3], P, xr)
		ax.axvline(x=sframe3.val, linestyle=':', color='r')
		axg.plot(xr, S.convert_I(total1, 30), 'k-')
		axg.plot(xr, S.convert_I(bulge1, 30), 'g:')
		axg.plot(xr, S.convert_I(disc1, 30), 'r--')
		axg.axvline(x=start_cut* P['ReD'].value / 1.67838865492, linestyle='--')
		axg.set_ylim([35,15])
		ax.set_ylim([0,1])
		# ax.set_ylim([0., ax.get_ylim()[1]])
		# l.set_ydata(BT_array[val1, val2, :])

		fig_gal.canvas.draw_idle()

	sframe1.on_changed(update)
	sframe2.on_changed(update)
	sframe3.on_changed(update)


	# l = ax.contourf(nB, MeB, BT_array[:,0,:])
	# cb = fig.colorbar(l,ax=ax)

	plt.show()