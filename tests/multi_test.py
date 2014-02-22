import multiprocessing as mp
import time
import pickle
import storeData as SD
import fitting as F

def process_gal(q, p):
	q.put([p.gal_name, p.cam_name, p.name, F.fit_sersic_exp(p, store=False, fit_range=[3., p.R[-1]])])

if __name__ == '__main__':
	G_list = pickle.load(open("saved_galaxies.p", 'r'))

	# pool = mp.Pool(processes = 8)
	# for G in G_list:
	# 	for C in G:
	# 		for P in C:
	# 			pool.apply_async(process_gal, args = (P, ))
	# pool.close()
	# start = time.clock()
	# pool.join()
	# result = pool.get()
	# print time.clock() - start

	result_queue = mp.Queue()
	fits = []
	for G in G_list:
		for C in G:
			for P in C:
				fits.append(process_gal(result_queue, P))
	jobs = [mp.Process(f) for f in fits]
	for job in jobs: job.start()
	for job in jobs: job.join()
	results = [result_queue.get() for f in fits]