"""
detect error in sky estimation and correct it
histogram
"""
import numpy as np
import matplotlib.pyplot as plt
import storeData as SD
import scipy.optimize as opt

def hist_area(x, Range):
	v = np.std(x[Range[0]:Range[1]:-1])
	mu = np.mean(x[Range[0]:Range[1]:-1])
	med = np.median(x[Range[0]:Range[1]:-1])
	synth = (3 * med) - (2 * mu)
	return v, mu, synth, med

def fit_sky(X, data):
	strt_res = lambda p, x, y: y - ( (p[0] * x) - p[1])
	return opt.leastsq(strt_res, [0., 0.], args=(X, data))

def iterate(I, R):
	ms, cs = [], []
	for i in range(len(R)-4,0,-1):
		sol = fit_sky(R[len(R):i:-1], I[len(R):i:-1])
		ms.append(sol[0][0])
		cs.append(sol[0][1])
	cs = np.array(cs)
	ms = np.array(ms)
	ans = np.min(np.abs(ms))
	return  R[0:-3][abs(ms[::-1]) == ans]


if __name__ == '__main__':
	direct = r'repository/mega_1237667322723631315_cts.ascii'
	G = SD.Galaxy('first')
	G.import_file(direct)
	p = G[0][0]
	loc = iterate(p.I, p.R)
	print loc
	plt.plot(p.R, p.I, 'b.')
	# plt.ylim([35,15])
	plt.axvline(x=loc)

	plt.show()