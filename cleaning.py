"""
read gal_names
for each gal, show profile and image from SDSS
if rejected ask for reason, then mark in exclusions.txt
"""
import storeData as SD
import os
import matplotlib.pyplot as plt
import urllib
import cStringIO
from PIL import Image
import re
import sys
from matplotlib.backends.backend_pdf import PdfPages

def get_SDSS_image(ra, dec, scale):
	baseurl = 'http://skyservice.pha.jhu.edu/DR10/ImgCutout/getjpeg.aspx?ra=%s&dec=%s&scale=%s&width=600&height=600&opt=IG' % (ra, dec, scale)
	return Image.open(cStringIO.StringIO(urllib.urlopen(baseurl).read()))

def get_position(ID):
	url = 'http://skyserver.sdss3.org/dr10/en/tools/explore/summary.aspx?id=%s' % ID
	page = urllib.urlopen(url).read()
	regex = re.compile("(?<=javascript:gotochart\()(.*)(?=\);)")
	return regex.findall(page)[0].split(',')
	
def retrieve_image(ID, scale):
	RA, DEC = get_position(ID)
	return get_SDSS_image(RA, DEC, scale)

def display_galaxy(ID, directory):
	"""displays sky (bottom 25% of I), overall profiles from both cams, SDSS image"""
	fig = plt.figure()
	image1 = retrieve_image(ID, 0.6) #scale=0.6 (area around it)
	image2 = retrieve_image(ID, 0.2) #scale=0.2 (close up)
	G = SD.Galaxy('temp')
	G.import_file(directory+'/mega_'+ID+'_cts.ascii', 'mega')
	G.import_file(directory+'/sdss_'+ID+'_cts.ascii', 'sdss')
	ax = fig.add_subplot(1,2,0)
	ax2 = fig.add_subplot(1,2,1)
	ax.imshow(image1)
	ax2.imshow(image2)

	fig1 = plt.figure()

	count = 0
	for i in range(2):
		for j in range(2):
			mag = G[i][j].M
			rad = G[i][j].R
			L = len(rad) /2
			inten = G[i][j].I[L:-1]
			ax = fig1.add_subplot(4,2,count)
			ax.plot(rad[L:-1], inten, 'b.')
			ax.set_title('outer profile of %s' % G[i][j].cam_name+G[i][j].name)
			ax.set_ylim(ax.get_ylim()[0], 100)
			ax.axhline(y=0, linestyle='--')
			ax2  = fig1.add_subplot(4,2,count+1)
			ax2.plot(rad, mag, 'b.')
			ax2.set_title('magnitude of %s' % G[i][j].cam_name+G[i][j].name)
			ax2.set_ylim([35,15])
			count += 2
	fig.suptitle(ID)
	plt.tight_layout()
	return fig, fig1


if __name__ == '__main__':
	with open(r'repository/gal_names.txt', 'r') as f:
		name_list = f.read().split('\n')
	name_list = [line.split(' ') for line in name_list] #GMP is element 2 in sublist
	T = len(name_list)
	# image_list = []
	# plot_list = []
	# print 'retrieving information...'
	# error_list = []
	# for cumul, name in enumerate(name_list):
	# 	if os.path.isfile('repository/report/graphs_'+name[0]+'.png'):
	# 		print name[0], ' already exists'
	# 	else:
	# 		sys.stdout.write('\r')
	# 		sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*cumul/T), (100*cumul/T)))
	# 		sys.stdout.flush()
	# 		try:
	# 			a, b = display_galaxy(name[0], 'repository')
	# 			a.set_size_inches(18.5,10.5)
	# 			b.set_size_inches(18.5,10.5)
	# 			a.savefig('repository/report/images_'+name[0]+'.png', dpi=300, orientation='landscape', papertype='a4', pad_inches=0.)
	# 			b.savefig('repository/report/graphs_'+name[0]+'.png', dpi=300, orientation='landscape', papertype='a4', pad_inches=0.)
	# 			plt.close(a)
	# 			plt.close(b)
	# 		except IOError:
	# 			error_list.append(name[0])
	# print 'files not found: ', error_list
	im = None
	plt.ion()
	for cumul, name in enumerate(name_list):
		if im is None:
			im = plt.imread(r'C:\Users\User\Documents\Project work\GalaxyFitting\repository\report'+'\graphs_'+name[0]+'.png')
			fig = plt.figure()
			ax = fig.add_subplot(111)
			img = ax.imshow(im)
		else:
			img.set_data(im)
		plt.draw()
		accept = raw_input('OK?')