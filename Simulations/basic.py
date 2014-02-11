import numpy as np
import matplotlib.pyplot as plt
import lmfit as lm
import matplotlib.patches as patches

def where_below():
	pass

	
def give_2d_region(vals, ind, BT, lim, m='max'):
	"""shows what value, on a 2d map, the chosen parameter must be below/above (min/max) to have a B/T below a limit.
	max for nB"""
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
			found = rated[temp <= 0.05]
			if len(found) == 0:
				arr[i,j] = np.nan
			elif m == 'max':
				arr[i,j] = max(found)
			elif m == 'min':
				arr[i,j] = min(found)
	return arr, a, b


if __name__ == '__main__':
	P = lm.Parameters()
	P.add_many(('MeB', None, True, 10.), ('ReB', None, True, 0.01), ('nB', None, True, 0.01),\
		('MeD', 21., True, 10.), ('ReD', 9., True, 0.01), ('nD', 1., False, 0.01))
	MeB = np.linspace(10,30, 49)
	ReB = np.linspace(0.1, P['ReD'].value*1.9, 50)
	nB = np.linspace(0.5, 5., 51)
	BT_array = np.load('BT_file.npy')
	print BT_array.shape
	fig = plt.figure()
	fig.set_facecolor('white')
	ax = fig.add_subplot(111)
	print len(BT_array[9,:,19])
	array, y, x = give_2d_region([MeB, ReB, nB], 2, BT_array, 0.05, 'max') # max for 1 and 2
	array[array >= 5] = np.nan
	CS = ax.contourf(x, y, array, 50)
	ax.contour(CS, levels=[CS.levels[0],CS.levels[-1]],colors = 'r', origin='lower', hold='on')
	cbar = fig.colorbar(CS,ax=ax)
	xmin, xmax = ax.get_xlim()
	ymin, ymax = ax.get_ylim()
	xy = (xmin,ymin)
	width = xmax - xmin
	height = ymax - ymin
	p = patches.Rectangle(xy, width, height, hatch='\\\\\\\\\\', fill=None, zorder=-10)
	ax.add_patch(p)
	ax.set_title('Maximum value of nB such that B/T < 0.05')
	ax.set_ylabel('MeB [mag arcsec$^{-1}$]')
	ax.set_xlabel('ReB [arcsec]')

	# fig2 = plt.figure()
	# fig2.set_facecolor('white')
	# ax2 = fig2.add_subplot(111)
	# for i in range(len(MeB))[::10]:
	# 	ax2.plot(ReB, array[i,:], label='%.1f' %MeB[i])#8


	# ax2.set_title('Maximum Value of nB at fixed MeB')
	# ax2.set_ylabel('Max nB')
	# ax2.set_xlabel('ReB')
	# leg = ax2.legend(loc=0)
	# leg.set_title('MeB Values')

	ax.set_ylim(ax.get_ylim()[::-1])
	plt.show()
