from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import lmfit as lm
import Simulate as S
import storeData as SD
import time
import matplotlib.gridspec as gridspec
from scipy.misc import factorial
import csv
from scipy.special import gamma

def get_b_n(m):
	#   find sersic b_n coefficient'favoritecolor'
	#   this is more accurate than the usual b_n = 1.9992*n_sersic - 0.3271
	if m == 0.0: return -0.3271
	b_n = 2.0*m - 1.0/3.0 + 4.0/m/405.0 + 46.0/m/m/25515.0 \
				+ 131.0/m/m/m/1148175.0  \
				- 2194697.0/m/m/m/m/30690717750.0
	return b_n

def BD_ratio(params, zp):
	n = params['nB']
	bnb = S.get_b_n(n)
	bnd = S.get_b_n(params['nD'])
	pre = n * gamma(2*n) * np.exp(bnb) / (bnb ** (2 * n))
	h = params['ReD'] / ((bnd) ** params['nD'])
	Ie = 10 ** ((params['MeD'] - zp) / -2.5)
	I0 = Ie * np.exp(bnd)
	return pre * ((params['ReB'] / h)**2.) * (Ie / I0)

def BT_ratio(x):
	return 1. / (1 + ((1./x)))

if __name__ == '__main__':
	# P = lm.Parameters()
	# P.add_many(('MeB', None, True, 10.), ('ReB', None, True, 0.01), ('nB', None, True, 0.01),\
	# 	('MeD', 21., True, 10.), ('ReD', 9., True, 0.01), ('nD', 1., False, 0.01))

	# s = 30
	# Z = 30.
	# MeB = np.linspace(P['MeD']*0.1, P['MeD']*1.9, s)
	# ReB = np.linspace(P['ReD']*0.1, P['ReD']*1.9, s)
	# nB = np.linspace(P['nD']*0.1, P['nD']*1.9, s)

	# BD_array = np.zeros((s,s,s))

	# for i, m in enumerate(MeB):
	# 	for j, r in enumerate(ReB):
	# 		for k, n in enumerate(nB):
	# 			P['MeB'].value = m; P['ReB'].value =r; P['nB'].value =n
	# 			BD_array[i,j,k] = BD_ratio(P, Z)


	# fig = plt.figure()
	# ax1 = fig.add_subplot(131, projection='3d')
	# ax2 = fig.add_subplot(132, projection='3d')
	# ax3 = fig.add_subplot(133, projection='3d')

	# X, Y = np.meshgrid(MeB, ReB)
	# Z = BD_array[:,:,np.floor(len(BD_array)/2)]
	# ax1.plot_surface(X, Y, Z, alpha=0.3)
	# ax1.set_xlabel('MeB')
	# ax1.set_ylabel('ReB')
	# ax1.set_zlabel('B/T')

	# X, Y = np.meshgrid(MeB, nB)
	# Z = BD_array[:,np.floor(len(BD_array)/2),:]
	# ax2.plot_surface(X, Y, Z, alpha=0.3)
	# ax2.set_xlabel('MeB')
	# ax2.set_ylabel('nB')
	# ax2.set_zlabel('B/T')

	# X, Y = np.meshgrid(ReB, nB)
	# Z = BD_array[np.floor(len(BD_array)/2),:,:]
	# ax3.plot_surface(X, Y, Z, alpha=0.3)
	# ax3.set_xlabel('ReB')
	# ax3.set_ylabel('nB')
	# ax3.set_zlabel('B/T')

	# plt.show()

	print get_b_n(1)