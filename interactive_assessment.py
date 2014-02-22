import storeData as SD
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import csv
from PIL import Image
import urllib
import cStringIO
import re
import sys

class image_plot():
	"""contains updaters for an image plot assessment"""
	def __init__(self, gal_names, folder, fig):
		self.gal_names = gal_names
		self.folder = folder
		self.fig = fig
		self.current = 0
		self.L = len(gal_names)
		self.accepted =  [True] * self.L
		self.update()
		plt.show()


	def next(self, event):
		try:
			if self.current == self.L - 1:
				print 'Completed!'
				plt.close()
			else:
				self.current += 1
				self.update()
		except IOError:
			print 'Completed!'
			plt.close()

	def back(self, event):
		if self.current != 0:
			self.current -= 1
			self.update()

	def accept(self, label):
		if label == 'accept':
			self.accepted[self.current] = True
		else:
			self.accepted[self.current] = False

	def update(self):
		self.fig.clf()
		self.draw_objects()
		self.im1 = plt.imread(self.folder+'/'+self.gal_names[self.current]+'_1.png')
		self.im2 = plt.imread(self.folder+'/'+self.gal_names[self.current]+'_2.png')
		self.ax1.imshow(self.im1)
		self.ax2.imshow(self.im2)
		# plt.ion()
		self.backbutton = Button(self.backax, '<')
		self.nextbutton = Button(self.nextax, '>')
		self.backbutton.on_clicked(self.back)
		self.nextbutton.on_clicked(self.next)
		self.fig.suptitle(self.gal_names[self.current])
		self.rejcheck = RadioButtons(self.rejax,('accept', 'reject'))
		self.rejcheck.on_clicked(self.accept)
		self.ax1.text(0.9, 0.9, '%i/%i' % (self.current, self.L))

		for i in [self.ax1, self.ax2]:
			plt.setp(i.get_xticklabels(), visible=False)
			plt.setp(i.get_yticklabels(), visible=False)
		self.fig.canvas.draw()

	def draw_objects(self):
		self.ax1 = self.fig.add_axes([0.265, 0.05, 0.475, 0.475])
		self.ax2 = self.fig.add_axes([0.265, 0.48, 0.475, 0.475])
		self.rejax = self.fig.add_axes([0.1, 0.5, 0.1, 0.1])
		self.backax = self.fig.add_axes([0.0, 0.0, 0.5, 0.05])
		self.nextax = self.fig.add_axes([0.5, 0.0, 0.5, 0.05])
		self.count = self.fig.text(0.9, 0.9, '%i/%i' % (self.current+1, self.L-1))

class profile_plot():
	"""contains updaters for profile plots"""
	def __init__(self, gal_names, fig):
		self.gal_names = gal_names
		self.fig = fig
		self.current = 0
		self.L = len(gal_names)
		self.accepted =  [True] * self.L
		self.sky = {'SDSS':[0]* self.L, 'mega':[0]* self.L}
		self.bulge = {'SDSS':[0]* self.L, 'mega':[0]* self.L}
		self.update()
		plt.show()


	def next(self, event):
		if self.current == self.L - 1:
			print 'Completed!'
			plt.close()
		else:
			self.current += 1
			self.sky['SDSS'][self.current] = self.sdss_sky.val
			self.bulge['SDSS'][self.current] = self.sdss_bulge.val
			self.sky['mega'][self.current] = self.mega_sky.val
			self.bulge['mega'][self.current] = self.mega_bulge.val
			self.update()

	def back(self, event):
		if self.current != 0:
			self.current -= 1
			self.sky['SDSS'][self.current] = self.sdss_sky.val
			self.bulge['SDSS'][self.current] = self.sdss_bulge.val
			self.sky['mega'][self.current] = self.mega_sky.val
			self.bulge['mega'][self.current] = self.mega_bulge.val
			self.update()

	def accept(self, label):
		if label == 'accept':
			self.accepted[self.current] = True
		else:
			self.accepted[self.current] = False

	def update(self):
		self.fig.clf()
		self.draw_objects()
		sdss1 = self.gal_names[self.current]['sdss'][0]
		sdss2 = self.gal_names[self.current]['sdss'][1]
		mega1 = self.gal_names[self.current]['mega'][0]
		mega2 = self.gal_names[self.current]['mega'][1]
		self.SDSS1.plot(sdss1.R, sdss1.I, 'b.')
		self.SDSS2.plot(sdss2.R, sdss2.I, 'b.')
		self.mega1.plot(mega1.R, mega1.I, 'b.')
		self.mega2.plot(mega2.R, mega2.I, 'b.')
		for i in [self.SDSS1, self.SDSS2, self.mega1, self.mega2]:
			i.set_ylim([-20.,20.])

		self.SDSS1M.plot(sdss1.R, sdss1.M, 'r.')
		self.SDSS2M.plot(sdss2.R, sdss2.M, 'r.')
		self.mega1M.plot(mega1.R, mega1.M, 'r.')
		self.mega2M.plot(mega2.R, mega2.M, 'r.')
		for i in [self.SDSS1M, self.SDSS2M, self.mega1M, self.mega2M]:
			i.set_ylim([35,15])

		self.backbutton = Button(self.backax, '<')
		self.nextbutton = Button(self.nextax, '>')
		self.backbutton.on_clicked(self.back)
		self.nextbutton.on_clicked(self.next)
		self.fig.suptitle(self.gal_names[self.current].ID)
		self.rejcheck = RadioButtons(self.rejax,('accept', 'reject'))
		self.rejcheck.on_clicked(self.accept)
		self.fig.text(0.9, 0.9, '%i/%i' % (self.current, self.L))

		self.sdss_sky = Slider(self.sky_SDSS, "SDSS sky", sdss1.R[0], sdss1.R[-1], valinit=sdss1.R[-10], color='#AAAAAA')
		self.mega_sky = Slider(self.sky_M, "M sky", mega1.R[0], mega1.R[-1], valinit=mega1.R[-10], color='#AAAAAA')
		self.sdss_bulge = Slider(self.bulge_SDSS, "SDSS bulge", sdss1.R[0], sdss1.R[-1], valinit=sdss1.R[-10], color='#AAAAAA')
		self.mega_bulge = Slider(self.bulge_M, "M bulge", mega1.R[0], mega1.R[-1], valinit=mega1.R[-10], color='#AAAAAA')

		self.sdss_sky.on_changed(self.update_slider)

		for i in [self.SDSS1M, self.SDSS2M, self.mega1M, self.mega2M]:
			plt.setp(i.get_xticklabels(), visible=False)
			plt.setp(i.get_yticklabels(), visible=False)
		self.fig.canvas.draw()

	def update_slider(self, val):
		pass
		# self.SDSS1.axvline(x=self.sdss_sky.val, color='b', linestyle='-')
		# self.SDSS2.axvline(x=self.sdss_sky.val, color='b', linestyle='-')
		# self.mega1.axvline(x=self.mega_sky.val, color='b', linestyle='-')
		# self.mega2.axvline(x=self.mega_sky.val, color='b', linestyle='-')

		# self.SDSS1.axvline(x=self.sdss_bulge.val, color='g', linestyle='--')
		# self.SDSS2.axvline(x=self.sdss_bulge.val, color='g', linestyle='--')
		# self.mega1.axvline(x=self.mega_bulge.val, color='g', linestyle='--')
		# self.mega2.axvline(x=self.mega_bulge.val, color='g', linestyle='--')

	def draw_objects(self):
		self.SDSS1 = self.fig.add_axes([0.05, 0.5, 0.45, 0.4])
		self.SDSS2 = self.fig.add_axes([0.05, 0.1, 0.45, 0.4])
		self.mega1 = self.fig.add_axes([0.5, 0.5, 0.45, 0.4])
		self.mega2 = self.fig.add_axes([0.5, 0.1, 0.45, 0.4])

		self.SDSS1M = self.SDSS1.twinx()
		self.SDSS2M = self.SDSS2.twinx()
		self.mega1M = self.mega1.twinx()
		self.mega2M = self.mega2.twinx()

		self.rejax = self.fig.add_axes([0., 0.1, 0.05, 0.05])
		self.backax = self.fig.add_axes([0.0, 0.45, 0.05, 0.1])
		self.nextax = self.fig.add_axes([0.95, 0.45, 0.05, 0.1])

		self.sky_SDSS = self.fig.add_axes([0.05, 0., 0.45, 0.03])
		self.bulge_SDSS = self.fig.add_axes([0.05, 0.05, 0.45, 0.03])
		self.sky_M = self.fig.add_axes([0.5, 0., 0.45, 0.03])
		self.bulge_M = self.fig.add_axes([0.5, 0.05, 0.45, 0.03])
		

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

def gather_data(names):
	T = len(names)
	for i, g in enumerate(names): 
		sys.stdout.write('\r')
		sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*i/T), (100*i/T)))
		sys.stdout.flush()
		fig1 = plt.figure()
		fig2 = plt.figure()
		ax = fig1.add_subplot(111)
		ax2 = fig2.add_subplot(111)
		im1 = retrieve_image(g, 0.1)
		im2 = retrieve_image(g, 0.6)
		ax.imshow(im1)
		ax2.imshow(im2)
		plt.tight_layout()
		fig1.savefig('report/'+g+'_1.png', dpi=300, orientation='landscape', papertype='a4', pad_inches=0.)
		fig2.savefig('report/'+g+'_2.png', dpi=300, orientation='landscape', papertype='a4', pad_inches=0.)
		plt.close(fig1)
		plt.close(fig2)

if __name__ == '__main__':
	with open(r'repository/gal_names.txt', 'r') as f:
		name_list = f.read().split('\n')
	name_list = [line.split(' ') for line in name_list] #GMP is element 1 in sublist
	name_list = [name_list[i][0] for i,v in enumerate(name_list)]
	
	figure = plt.figure()
	I = image_plot(name_list, 'report', figure)
	with open(r'repository/exclusions.txt', 'wb') as f:
		for i,v in enumerate(name_list):
			if not I.accepted[i]:
				f.write(v+'\n')
	# direct = 'repository'
	# gals, names = SD.import_directory(direct)
	# I = profile_plot(gals, figure)
	# print I.sky['SDSS']
