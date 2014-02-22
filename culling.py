"""
prepare csv
for every name in name_list:
	plot images in figure 1
	add button on screen to ask for acceptance 
	draw fig1
	if no then:
		add to exclusions list
		continue
	import SDSS+Mega using storeData
	draw double plot of SDSS intensities (<100) + mags in figure 2
	add sliders to figure 2 adjusting sky range and bulge-cutoff
	add button to confirm
	draw fig2
	store template guesses and recorded params in csv
"""

import storeData as SD
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import csv
from PIL import Image
import urllib
import cStringIO
import re
import sys

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

def iterate_config(galaxies, writer):
	fig1, fig2 = plt.figure(), plt.figure()
	zoom = fig1.add_subplot(111)
	ax1, ax2, ax3, ax4 = fig2.add_subplot(221), fig2.add_subplot(222),fig2.add_subplot(223),fig2.add_subplot(224)
	plt.ion()
	plt.show()
	T = len(galaxies)
	for i, g in enumerate(galaxies):
		sys.stdout.write('\r')
		sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*i/T), (100*i/T)))
		sys.stdout.flush()
		zoom.imshow('report/'+g.ID+'.png')
		acceptax = plt.axes([0.8, 0.025, 0.1, 0.04])
		denyax = plt.axes([0.9, 0.025, 0.1, 0.04])
		button = Button(acceptax, 'Deny', hovercolor='0.975')

		

def gather_data(galaxies):
	T = len(galaxies)
	for i, g in enumerate(galaxies):
		sys.stdout.write('\r')
		sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*i/T), (100*i/T)))
		sys.stdout.flush()
		fig = plt.figure()
		ax = fig.add_subplot(211)
		ax2 = fig.add_subplot(212)
		im1 = retrieve_image(g.ID, 0.1)
		im2 = retrieve_image(g.ID, 0.6)
		ax.imshow(im1)
		ax2.imshow(im2)
		fig.set_size_inches(18.5,10.5)
		plt.tight_layout()
		fig.savefig('report/'+g.ID+'.png', dpi=300, orientation='landscape', papertype='a4', pad_inches=0.)
		fig.clf()
		


if __name__ == '__main__':
	direct = r'C:\Users\User\Documents\Project work\GalaxyFitting\repository'
	gal_list, gal_names = SD.import_directory(direct)
	# with open('config.csv', 'wb') as f:
	# 	csvwriter = csv.writer(f)
	# 	iterate_config(gal_list, csvwriter)
	gather_data(gal_list)