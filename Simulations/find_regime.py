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
import matplotlib.patches as patches
from matplotlib.ticker import MultipleLocator

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

def mesh3d(*arrs):
	"""x,y,z, etc"""
	arrs = arrs[::-1]  #edit
	lens = map(len, arrs)
	dim = len(arrs)

	sz = 1
	for s in lens:
		sz*=s

	ans = []    
	for i, arr in enumerate(arrs):
		slc = [1]*dim
		slc[i] = lens[i]
		arr2 = np.asarray(arr).reshape(slc)
		for j, sz in enumerate(lens):
			if j!=i:
				arr2 = arr2.repeat(sz, axis=j) 
		ans.append(arr2)

	return tuple(ans)

def find_change(x):
	"""finds first change in value"""
	c = []
	for i, v in enumerate(x[1:]):
		if x[i] != x[i-1]:
			c.append(i)
	if len(c) == 0:
		return None
	else:
		return c[-1]

	
def give_limited_2d(vals, ind, BT, lim, choice='min'):
	"""displays as colour, the max value of vals[ind] which has BT < lim, on a 2d map """
	if ind == 0: # M
		a, b, rated = vals[1], vals[2], vals[0]
	elif ind == 1: #R
		a, b, rated = vals[0], vals[2], vals[1]
	else: #n
		a, b, rated = vals[0], vals[1], vals[2]
	arr = np.zeros((len(a), len(b)))
	for i in range(len(a)):
		for j in range(len(b)):
			if ind == 0:
				temp = BT[:,i,j]
			elif ind == 1:
				temp = BT[i,:,j]
			else:
				temp = BT[i,j,:]
			limits = (temp < lim)
			chg = find_change(limits)
			if chg is None:
				arr[i, j] = np.nan
			else:
				arr[i, j] = rated[chg]
	X, Y = np.meshgrid(a, b)
	return arr, X, Y

def give_2d_region(vals, ind, BT, lim, m='max'):
	if ind == 0: # M
		a, b, rated = vals[1], vals[2], vals[0]
	elif ind == 1: #R
		a, b, rated = vals[0], vals[2], vals[1]
	else: #n
		a, b, rated = vals[0], vals[1], vals[2]
	arr = np.zeros((len(a), len(b)))
	for i in range(len(a)):
		for j in range(len(b)):
			if ind == 0:
				temp = BT[:,i,j]
			elif ind == 1:
				temp = BT[i,:,j]
			else:
				temp = BT[i,j,:]
			arr[i,j] = find_such(rated, temp, lim, m)
	X, Y = np.meshgrid(a, b)
	return arr, X, Y

def find_such(values, BT_reduced, lim, m):
	"""in a 1d space of set BT[coord1,coord2] find max/min (in m) of values such that BT_reduced is less than lim"""
	a = values[BT_reduced <= lim]
	if len(a) == 0:
		return np.nan
	if m is 'max':
		return a[-1]
	else:
		return a[0]


if __name__ == '__main__':
	P = lm.Parameters()
	P.add_many(('MeB', None, True, 10.), ('ReB', None, True, 0.01), ('nB', None, True, 0.01),\
		('MeD', 21., True, 10.), ('ReD', 9., True, 0.01), ('nD', 1., False, 0.01))
	s = 50
	MeB = np.linspace(10,30, 10)
	ReB = np.linspace(0.1, P['ReD'].value*1.9, 5)
	nB = np.linspace(0.5, 5., 20)
	start_cut = 4.; end_cut = 50.
	# BD_array = np.zeros((len(MeB),len(ReB),len(nB)))
	# print 'start!'
	# bar_string = ' '*100
	# T = len(MeB) * len(ReB) * len(nB)
	# cumul = 0
	# for i, m in enumerate(MeB):
	# 	for j, r in enumerate(ReB):
	# 		for k, n in enumerate(nB):
	# 			cumul += 1
	# 			# sys.stdout.write('\r')
	# 			# sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*cumul/T), (100*cumul/T)))
	# 			# sys.stdout.flush()
	# 			BD_array[i,j,k] = apply_func(m,r,n,P,start_cut,end_cut)
	# print '\nend!'

	# BT_array = 1. / (1. + ((1./BD_array)))

	# np.save('BT_file', BT_array)

	BT_array = np.load('BT_file.npy')
	# sheet = (len(MeB)*len(ReB), len(nB))
	# BT_sheet = BT_array.reshape(sheet) # reshaped for nB in y-axis and 
	# mm, rr, nn = mesh3d(MeB,ReB,nB)
	# mr = rr*mm

	# mr_sheet = mr.reshape(sheet)
	# nn_sheet = nn.reshape(sheet)
	# fig = plt.figure(); ax = fig.add_subplot(111)
	# l = ax.contourf(mr_sheet, nn_sheet, BT_sheet)
	# fig.colorbar(l,ax=ax)
	# ax.set_ylabel('nB')
	# ax.set_xlabel('R/M')

	# find BT_array[M,R,n]
	# plt_array, X, Y = give_limited_2d([MeB, ReB, nB], 2, BT_array, 0.05) # find nB such that BT < 0.05
	plt_array, X, Y = give_2d_region([MeB, ReB, nB], 1, BT_array, 0.05, 'min')
	fig = plt.figure(); ax = fig.add_subplot(111)
	CS = ax.contourf(nB, MeB, plt_array,10,cmap=plt.cm.PuBu,origin='lower')
	# ax.contour(CS, levels=CS.levels[::2],colors = 'r', origin='lower', hold='on')
	
	xmin, xmax = ax.get_xlim()
	ymin, ymax = ax.get_ylim()
	xy = (xmin,ymin)
	width = xmax - xmin
	height = ymax - ymin
	p = patches.Rectangle(xy, width, height, hatch='\\\\\\\\\\', fill=None, zorder=-10)
	ax.add_patch(p)
	# ax.set_ylabel('ReB/ReD')
	# ax.set_xlabel('MeB/MeD')
	cbar = fig.colorbar(CS,ax=ax)
	cbar.ax.set_ylabel('nB')
	ax.set_title('Values of nB')
	plt.show()