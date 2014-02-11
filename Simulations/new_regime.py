import Simulate as S
import lmfit as lm
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
import time

def get_b_n(m):
	#   find sersic b_n coefficient
	#   this is more accurate than the usual b_n = 1.9992*n_sersic - 0.3271
	# if m == 0.0: return -0.3271
	b_n = 2.0*m - 1.0/3.0 + 4.0/m/405.0 + 46.0/m/m/25515.0 \
				+ 131.0/m/m/m/1148175.0  \
				- 2194697.0/m/m/m/m/30690717750.0
	return b_n

def mesh3d(*arrs):
	"""x,y,z, etc"""
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


def mag_diff(x, md, rd, mb, rb, nb, show='all'):
	disc = md + (2.5 * get_b_n(1.) * ((x /rd ) ** (1/1.) - 1)) / np.log(10)
	bulge =  mb + (2.5 * get_b_n(nb) * ((x / rb ) ** (1./nb) - 1.)) / np.log(10)
	if show is 'all':
		return bulge - disc, disc, bulge
	else:
		return bulge - disc

def mag_grad_diff(x, md, rd, mb, rb, nb, show='all'):
	disc =  (2.5 * get_b_n(1.) * (((1 / (rd * 1.)) * ((x / rd ) ** ((1-1.)/1.))))) / np.log(10)
	bulge =  (2.5 * get_b_n(nb) * (((1 / (rb * nb)) * ((x / rb ) ** ((1.-nb)/nb))))) / np.log(10)
	if show is 'all':
		return bulge - disc, disc, bulge
	else:
		return bulge - disc

def find_max_diff(md, rd, mb, rb, nb):
	"""finds where d delta_mu / dr = 0 and returns R and mu"""
	b_ratio = get_b_n(1.) / get_b_n(nb)
	rnb = rb ** (1. / nb)
	rnd = rd
	rn_ratio = rnb / rnd
	try:
		n_power = (nb) / (1. - nb)
	except ZeroDivisionError:
		return np.inf
	X = (b_ratio * rn_ratio * nb) ** n_power
	return X

def find_good_region(md, rd, mb, rb, nb, limit):
	xx = np.arange(0,1000)
	ans = mag_diff(md, rd, mb, rb, nb, xx)
	whr = np.where(ans[0] > limit)
	if whr[0].size == 0:
		return 0
	elif whr[0][-1] == 999:
		return np.nan
	else:
		return whr[0][-1]

def batch(md, rd, M_list, R_list, n_list, limit, m='max'):
	M, R, N = mesh3d(M_list, R_list, n_list)
	temp = np.ones(M.shape)
	MD, RD = temp*md, temp*rd
	max_locs = find_max_diff(MD, RD, M, R, N) + 1 # off-max locations
	grad_val = mag_grad_diff(MD, RD, M, R, N, max_locs)
	return grad_val[0]

	
def func(x, limit, mb, md, rd, rb, nb):
	disc = md + (2.5 * get_b_n(1.) * ((x /rd ) ** (1/1.) - 1)) / np.log(10)
	bulge =  mb + (2.5 * get_b_n(nb) * ((x / rb ) ** (1./nb) - 1.)) / np.log(10)
	return disc - bulge -limit

def batch_cross_finder(md, rd, M_list, R_list, n_list, limit, maxiter=50):
	M, R, N = mesh3d(M_list, R_list, n_list)
	temp = np.ones(M.shape)
	MD, RD = temp*md, temp*rd
	x = np.ones(MD.shape) * 1000.
	tol = 1
	new = 1
	for i in range(maxiter):
		diff = mag_diff(x, MD, RD, M, R, N)[0] - limit
		grad = mag_grad_diff(x, MD, RD, M, R, N)[0]
		new = x - (diff / grad)
		x = new
	return x

def single_newton(md,rd,mb,rb,nb,limit,maxiter=50):	
	x = 1000.
	tol = 0.1
	new = 1
	for i in range(maxiter):
		diff = mag_diff(x, md, rd, mb, rb, nb)[0] - limit
		grad = mag_grad_diff(x, md, rd, mb, rb, nb)[0]
		new = x - (diff / grad)
		x = new
	return new

if __name__ == '__main__':
	MeD, ReD = 21., 9.
	MeB = np.linspace(10,30, 49)
	ReB = np.linspace(0.1,18, 50)
	nB = np.linspace(0.5, 5., 51)
	BT_array = np.load('BT_file.npy')
	M,R,N=  mesh3d(MeB, ReB, nB) # out=nRM
	X = np.arange(0.01,1000,1)
	mb_test, rb_test, nb_test = 18., 3., 1.1

	diff, D, B = mag_diff(X, MeD,ReD, mb_test,rb_test,nb_test)
	grad_diff, grad_D, grad_B = mag_grad_diff(X, MeD,ReD,mb_test,rb_test,nb_test)
	plt.plot(X, D, 'r--')
	plt.plot(X, B, 'g:')
	plt.plot(diff)

	plt.plot(grad_diff)
	max_point = find_max_diff(MeD, ReD, mb_test, rb_test, nb_test)
	print max_point
	# print mag_grad_diff(MeD,ReD,mb_test,rb_test,nb_test, max_point+100)[0]
	print mag_diff(max_point+1, MeD,ReD, mb_test,rb_test,nb_test)[0] - \
	mag_diff(max_point, MeD,ReD, mb_test,rb_test,nb_test)[0]

	# plt.ylim([35,15])
	# plt.xlim([0,60])

	# print batch(MeD, ReD, MeB, ReB, nB, 1.)[25,25,25]
	# zeros = batch_cross_finder(MeD, ReD, MeB, ReB, nB, 1.)
	# print zeros.shape
	a = single_newton(MeD, ReD, mb_test, rb_test, nb_test, 1.)
	print a


	# b = newton(mag_diff, 1000., mag_grad_diff, args=(MeD,ReD,mb_test,rb_test,nb_test, 'not'))
	# print b
	plt.axvline(x=max_point)
	plt.show()
