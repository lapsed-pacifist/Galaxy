import multiprocessing
import Queue
import Simulate as S
import lmfit as lm
import storeData as SD
from scipy import stats
import numpy as np
import csv
import time

def wrapper(queue, P, model, x, z, data, weights, fit_range, redchi_marker):
	out, res_excl = S.fit(P, model, x, z, data, weights, fit_range, redchi_marker)
	queue.put([out, res_excl])
	queue.close()

def Fitting(P, model, x, z, data, weights=None, fit_range=None, redchi_marker=None):
	queue = multiprocessing.Queue(1) # Maximum size is 1
	proc = multiprocessing.Process(target=wrapper, args=(queue, P, model, x, z, data, weights, fit_range, redchi_marker))
	proc.start()
	# Wait for 10 seconds
	try:
		result = queue.get(True, 10)
	except Queue.Empty:
		# Deal with lack of data somehow
		result = [None, None]
	finally:
		proc.terminate()
	return result[0], result[1]

				
def make_noise(array, multiple):
	noise = np.random.random(array.shape) * float(multiple)
	return noise

def start_routine(filename, P, Me_range, Re_range, n_range, noise_level, R, zp):
	noise = make_noise(R, noise_level)

	with open(filename, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		header = ['%s = %.2f' % (P[i].name, P[i].value) for i in ['MeD', 'ReD', 'nD']]
		writer.writerow(header + ['range= [%.1f, %.1f]' % (R[0], R[-1])] + ['noise level = %.1f' % noise_level])
		writer.writerow(['MeB_initial', 'ReB_initial', 'nB_initial', 'MeD_final', 'ReD_final', 'nD_final', 'MeB_final', 'ReB_final', 'nB_final',\
		 'redchi2_all', 'redchi2_excl', 'KS', 'KS_excl'])
		for nB in n_range:
			for ReB in Re_range:
				for MeB in Me_range:
					pars = S.copy_params(P, False)
					pars.add_many(('MeB', float(MeB), True, 1.), ('ReB', float(ReB), True, 0.01), ('nB', float(nB), True, 0.1))
					pars['nD'].vary = False
					test_gal = S.sersic2(pars, R, zp, False) + noise

					new_pars = S.copy_params(pars, False)
					
					fit_data, res_excl = S.fit(new_pars, S.sersic2, R, zp, test_gal, weights=None, fit_range=None, redchi_marker=30.)

					initials = [pars[i].value for i in ['MeB', 'ReB', 'nB']]
					if fit_data is None:
						writer.writerow(['N/A'] * 13) 
					else:
						finals = [new_pars[i].value for i in ['MeB', 'ReB', 'nB', 'MeD', 'ReD', 'nD']]
						redchi_excl = np.sum(res_excl) / fit_data.nfree
						KS, KS_excl = stats.kstest(fit_data.residual, 'norm')[1], stats.kstest(res_excl, 'norm')[1]
						writer.writerow(initials+finals+['%r' % (res_excl)])

if __name__ == '__main__':
	temp_gal = SD.Galaxy('temp')
	temp_gal.import_file('C:\\Users\\User\\code\\Python\\Sersic_Plots\\mega_1237665428090192022_cts.ascii', 'mega')
	R = temp_gal[0][0].R

	P = lm.Parameters()
	P.add('MeD', value=21, min=1.)
	P.add('ReD', value=9., min=0.01)
	P.add('nD', value=1., min=0.1)
	Mes = [18.]
	Res = np.arange(0.1, 10, 0.1)
	ns = [4.]

	timer = time.clock()
	start_routine('first_test.csv', P, Mes, Res, ns, 20., R, 30.)
	print 'time taken: %.1f' % (time.clock() - timer)
