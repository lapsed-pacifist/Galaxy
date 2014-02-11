"""
make lists of variables to test M, R, N
make 3d array of meshgrids
find extrema for 3d array
find 2nd zero for 3d array

pipeline tests:
	1. Test for divergence at extremum
	2. B > D (at various points)
	3. Test for x > a point (bounds)

if divergence:
	if B > D-lim at extremum:
		reject
	else:
		accept
else:
	if B > D-lim beyond zero:
		if 

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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
		return bulge - disc # +ve if disc is greater in mag than bulge

def mag_grad_diff(x, md, rd, mb, rb, nb, show='all'):
	disc =  (2.5 * get_b_n(1.) * (((1 / (rd * 1.)) * ((x / rd ) ** ((1-1.)/1.))))) / np.log(10)
	bulge =  (2.5 * get_b_n(nb) * (((1 / (rb * nb)) * ((x / rb ) ** ((1.-nb)/nb))))) / np.log(10)
	if show is 'all':
		return bulge - disc, disc, bulge
	else:
		return bulge - disc

def find_extremum(md, rd, mb, rb, nb):
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

def find_zero(md, rd, mb, rb, nb, limit, maxiter=50):
	x = np.ones(md.shape) * 1000.
	tol = 1
	new = 1
	for i in range(maxiter):
		diff = mag_diff(x, md, rd, mb, rb, nb, 'none') - limit
		grad = mag_grad_diff(x, md, rd, mb, rb, nb, 'none')
		new = x - (diff / grad)
		x = new
	return x

def pipeline(md, rd, mb, rb, nb, limit, bounds):
	extreme_arr = find_extremum(md, rd, mb, rb, nb)
	diff_up = mag_diff(extreme_arr+1, md, rd, mb, rb, nb, 'none') - limit # finds B-D+lim up
	diff_down = mag_diff(extreme_arr-1, md, rd, mb, rb, nb, 'none') - limit # finds B-D+lim below
	diff_at = mag_diff(extreme_arr, md, rd, mb, rb, nb, 'none') - limit

	zeros = find_zero(md, rd, mb, rb, nb, limit) # finds furthest zeros of (B-D)
	diff_beyond_zero = mag_diff(zeros+1, md, rd, mb, rb, nb, 'none') - limit

	B_domin_at_extreme = diff_at >= 0 # is B - D > limit?
	divergence_at_extreme = (np.abs(diff_at - diff_down) >= 0) & (np.abs(diff_up - diff_at) >= 0) # does it diverge at extreme?
	disc_domin_after_zero = diff_beyond_zero >= 0 #(is disc dominated beyond 0?)
	below_lower_bound = zeros < bounds[0]
	above_upper_bound = zeros >= bounds[1]

	accepted = ((divergence_at_extreme) & (B_domin_at_extreme)) | ((disc_domin_after_zero) & (below_lower_bound)) | \
				((np.logical_not(disc_domin_after_zero)) & (above_upper_bound))
	return accepted

def view_2d(accepted_array, axis, mesh, m='max'):
	arr = accepted_array.astype(float) * mesh
	arr[arr==0] = np.nan
	if m is 'max':
		return np.nanmax(arr, axis=axis)
	elif m is 'min':
		return np.nanmin(arr, axis=axis)
	else:
		raise ValueError('Only max or min allowed')

def mask_minmax(map2d, mesh):
	Max, Min = np.max(mesh), np.min(mesh)
	map2d[(map2d == Max) | (map2d == Min)] = np.nan
	return map2d

if __name__ == '__main__':
	MeD, ReD = 21., 9.
	M_list = np.linspace(10,30, 49)
	R_list = np.linspace(0.1,18, 50)
	n_list = np.linspace(0.5, 5., 51)
	M, R, N = mesh3d(M_list, R_list, n_list)
	temp = np.ones(M.shape)
	MD, RD = temp*MeD, temp*ReD
	# ans = pipeline(MD, RD, M, R, N, 1., [10., 30.])
	# print 'percentage accepted: ', np.sum(ans) / float(ans.size)
	# np.save('accepted_file', ans)
	# print 'saved!'
	ans = np.load('accepted_file.npy')


	Map_N = view_2d(ans, 2, N, 'max') # axis 0, 1, 2 = M, R, N
	# Map_N[(Map_N > 5.) | (Map_N < 0.5)] = np.nan

	Map_R = view_2d(ans, 1, R, 'max') # axis 0, 1, 2 = M, R, N
	# Map_R[(Map_R > 18.) | (Map_R < 0.1)] = np.nan

	Map_M = view_2d(ans, 0, M, 'min') # axis 0, 1, 2 = M, R, N
	# Map_M[(Map_M >= 30.) | (Map_M <= 10.)] = np.nan

	fig1 = plt.figure(1); ax1 = fig1.add_subplot(111)
	fig2 = plt.figure(2); ax2 = fig2.add_subplot(111)
	fig3 = plt.figure(3); ax3 = fig3.add_subplot(111)

	CS1 = ax1.contourf(R_list, M_list, Map_N, 70)
	cbar1 = fig1.colorbar(CS1,ax=ax1)
	ax1.set_ylim(ax1.get_ylim()[::-1])

	CS2 = ax2.contourf(n_list, M_list, Map_R, 70)
	cbar2 = fig2.colorbar(CS2,ax=ax2)
	ax2.set_ylim(ax2.get_ylim()[::-1])

	CS3 = ax3.contourf(n_list, R_list/9., Map_M - 21., 70)
	cbar3 = fig3.colorbar(CS3,ax=ax3)
	# ax3.set_ylim(ax3.get_ylim()[::-1]

	# Ms = [15., 25.]; ns = [0.5, 2.5]

	# Ms_coord = [np.where(M_list == i)[0][0] for i in Ms]
	# ns_coord = [np.where(n_list == i)[0][0] for i in ns]

	# Ms_coord = [0, -1]
	# ns_coord = [0, -1]

	# print np.min(Map_R[Ms_coord[0]:Ms_coord[1],ns_coord[0]:ns_coord[1]])


	# ax2.axhline(y=M_list[Ms_coord[0]], color='w')
	# ax2.axvline(x=n_list[ns_coord[0]], color='w')

	# ax2.axhline(y=M_list[Ms_coord[1]], color='w')
	# ax2.axvline(x=n_list[ns_coord[1]], color='w')

	# fig = plt.figure()
	# ax = fig.add_subplot(111)
	
	m_2, r_2 = np.meshgrid(M_list, R_list)
	# ax.plot(((r_2/9.) / (21.-m_2)).ravel(), Map_N.ravel(), 'b.')

	R_choice = 0.2* 9
	whr = np.where(R_list < R_choice)[0][-1]
	print np.nanmax(Map_M[r_2 <=R_choice] - 21.) 
	for i in [fig1,fig2,fig3]:
		i.set_facecolor('white')
	ax3.set_ylabel('$R_{e,B}/R_{e,D}$', fontsize=16)
	ax3.set_xlabel('$n_B$', fontsize=16)
	cbar3.ax.set_ylabel('$\mu_{e,B} - \mu_{e,D} [\\rm{mag\/arcsecond^{-2}}$', fontsize=16)
	ax3.axhline(y = R_choice/9., color='k', linestyle='--')
	ax3.set_title('$Minimum\/\mu_{e,B} - \mu_{e,D}\/possible\/over\/a\/range\/of\/variables$')


	xmin, xmax = ax3.get_xlim()
	ymin, ymax = ax3.get_ylim()
	xy = (xmin,ymin)
	width = xmax - xmin
	height = ymax - ymin
	p = patches.Rectangle(xy, width, height, hatch='\\\\\\\\\\', fill=None, zorder=-10)
	ax3.add_patch(p)

	plt.show()
