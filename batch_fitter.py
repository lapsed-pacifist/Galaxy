import storeData as SD
import plotData as plt
import fitting as F
import sys
from matplotlib.backends.backend_pdf import PdfPages
import shutil

direct = 'repository'
Gal_list, Gal_names = SD.import_directory(direct)
T = len(Gal_list) * 4
cumul = 0

figures = []
for G in Gal_list:
	for C in G:
		for P in C:
			sys.stdout.write('\r')
			sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*cumul/T), (100*cumul/T)))
			sys.stdout.write(" name: %s" % C.gal_name+C.name+P.name)
			sys.stdout.flush()
			F.fit_sersic_exp(P, store=True, fit_range=[3., P.R[-1]])
			F.sky_error(P, P.sky_var, 'sersic+exp')
			cumul += 1
		try:
			fig = plt.graph_camera(C, 'sersic+exp')
			fig.set_size_inches(18.5,10.5)
			fig.savefig('repository/graphs/'+C.gal_name+'_'+C.name+P.name+'.png', dpi=300, orientation='landscape', papertype='a4', pad_inches=0.)
		except ValueError:
			shutil.move('repository/'+P.cam_name+'_'+P.gal_name+'_cts.ascii', 'repository/bad')
			shutil.move('repository/'+P.cam_name+'_'+P.gal_name+'.ini', 'repository/bad')
			sys.stdout.write(' bad one moved :(')
		except IndexError:
			shutil.move('repository/'+P.cam_name+'_'+P.gal_name+'_cts.ascii', 'repository/bad')
			shutil.move('repository/'+P.cam_name+'_'+P.gal_name+'.ini', 'repository/bad')
			sys.stdout.write(' bad one moved :(')


# pdf_pages = PdfPages('my-fancy-document.pdf')
# for i in figures:
# 	pdf_pages.savefig(i)
# pdf_pages.close()

