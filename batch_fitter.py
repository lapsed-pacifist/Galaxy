import storeData as SD
import fitting as F
import sys
import numpy as np
import csv

direct = 'repository'
Gal_list, Gal_names = SD.import_directory(direct)
T = len(Gal_list) * 4
cumul = 0
pars = ['MeB', 'ReB', 'nB', 'MeD', 'ReD']
results = np.zeros([T/4, 4, len(pars)+1])
std_errs = np.zeros([T/4, 4, len(pars)])

for i, G in enumerate(Gal_list):
	gal_n = -1
	for j, C in enumerate(G):
		for k, P in enumerate(C):
			cumul += 1
			gal_n += 1
			ans = F.fit_sersic_exp(P, store=True, fit_range=[2., P.R[-1]])
			parameters = [ans.params[a].value for a in pars]
			parameters.append(ans.redchi)
			errs = [ans.params[a].stderr for a in pars]
			results[i, gal_n, :], std_errs[i, gal_n, :] = np.array(parameters), np.array(errs)
			sys.stdout.write('\r')
			sys.stdout.write("[%-20s] %.1f%% (evals=%-4i)" % ('='*int(20*cumul/T), (100*cumul/T), ans.nfev))
			sys.stdout.flush()

with open('fit_parameters.csv', 'wb') as f:
	writer = csv.writer(f)
	writer.writerow()