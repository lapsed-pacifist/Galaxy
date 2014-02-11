from multiprocessing import Pool
from time import time
import lmfit as lm
import Simulate as S
import storeData as SD
import numpy as np

K = 50
def CostlyFunction((P, model, x, z, data, weights, fit_range)):
	return S.fit(P, model, x, z, data, weights, fit_range)
	

if __name__ == "__main__":
	temp_gal = SD.Galaxy('temp')
	temp_gal.import_file('C:\\Users\\User\\code\\Python\\Sersic_Plots\\mega_1237665428090192022_cts.ascii', 'mega')
	R = temp_gal[0][0].R
	W = temp_gal[0][0].W
	zp = 30.
	P = lm.Parameters()
	P.add_many(('MeB', 18.), ('ReB', 0.5), ('nB', 4.), ('MeD', 21.), ('ReD', 8.), ('nD', 1.))
	total1, bulge1, disc1 = S.sersic2(P, R, zp, True)
	noise = np.random.random(total1.shape) * 50.
	prof = SD.Profile('test', total1+noise, R, W, zp, 0, 1)

	P['nD'].vary = False
	P['ReD'].value = 10.; P['MeD'].value = 22.; P['nB'].value = 3.; P['ReB'].value = 4.; P['MeB'].value = 14.

	currtime = time()
	N = 10
	po = Pool()
	res = po.map_async(CostlyFunction,(P, S.sersic2, R, zp, prof.I, prof.W, None))
	print '2: parallel: time elapsed:', time() - currtime
	print res.get().params


