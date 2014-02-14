import re
import asciitable
import numpy as np
import ConfigParser as cp
import matplotlib.pyplot as plt
from warnings import warn
import bisect
import lmfit as lm
import os
import fnmatch
import findSky 


def translate_x(x, point, dupl_pos=-1, round_method='right'):
	"""converts a point in ordered list x to the appropriate index of x. 
	For duplicates, it returns the last index in the set of duplicates by default."""
	if point < x[0]: return 0
	if point > x[-1]: return len(x)
	if round_method is 'left': #finds closest by rounding down 
		i = bisect.bisect_right(x, point) - 1
	elif round_method is 'right': #finds closest by rounding up 
		i = bisect.bisect_right(x, point)
	return i
	
def ConfigSectionMap(ParseObj, section):
	dict1 = {}
	options = ParseObj.options(section)
	for option in options:
		try:
			dict1[option] = ParseObj.get(section, option)
			if dict1[option] == -1:
				DebugPrint("skip: %s" % option)
		except:
			print("exception on %s!" % option)
			dict1[option] = None
	return dict1

def string_to_list(s):
	string = re.search(re.escape('[') + "(.*?)" + re.escape(']'), s).group(1).split()[0]
	return map(float, string.split(','))
	

def convert_ascii(table):
	"Converts ascii table to numpy array"
	x_len = len(table[0])
	y_len = len(table)
	ex = 0
	# check for method column:
	if type(table[0][1]) != 'numpy.float64': 
		ex = 1
	array = np.zeros([y_len, x_len - ex])
	for j in range(y_len):
		for i in range(ex, x_len):
			array[j, i-ex] = table[j][i]
	return array

def extract_headers(filename):
	with open(filename) as f:
		string = f.readline()
	headers = ['ID', 'scale', 'zp', 'sky']
	result = re.search(re.escape('###  ') + "(.*?)" + re.escape('\n'), string).group(1).split()
	d = {headers[i]: value for (i, value) in enumerate(result)}
	return d

class Galaxy(object):
	"""contains cameras and galaxy info"""
	def __init__(self, name):
		self.name = name
		self.camera_list = []

	def import_file(self, filename, camera_name=None):
		h = extract_headers(filename)
		T = asciitable.read(filename)
		table = convert_ascii(T)
		cam = Camera(camera_name, float(h['zp']), float(h['sky']), float(h['scale']), self.name) #create a camera from headers
		cam.add_data_table(table) #call add_table_data in camera
		self.camera_list.append(cam)
		self.ID = h['ID']
		# print re.search(re.escape('\\') + "(.*?)" + re.escape(self.ID), filename).group(1).split()
		# print re.match('.*?(\w+)(?<=.{%d})' % (0), filename).group(0)

	def __getitem__(self, index):
		if type(index) is str:
			names = [i.name for i in self.camera_list]
			return self.camera_list[names.index(index)]
		elif type(index) is int:
			return self.camera_list[index]
		else: raise ValueError('index must be a name (string) or an index (int)')

	def __len__(self):
		return len(self.camera_list)

class Camera(object):
	"""Contains Profiles and galaxy info"""
	def __init__(self, name, zeropoint, skylevel, scale, gal_name):
		self.profile_list = []
		self.name = name
		self.zeropoint = zeropoint
		self.skylevel = skylevel
		self.scale = scale
		self.gal_name = gal_name

	def __len__(self):
		return len(self.profile_list)

	def __getitem__(self, index):
		if type(index) is str:
			names = [i.name for i in self.profile_list]
			return self.profile_list[names.index(index)]
		elif type(index) is int:
			return self.profile_list[index]
		else: raise ValueError('index must be a name (string) or an index (int)')


	def add_data_table(self, data):
		one = Profile('1', data[:,1], data[:,0], data[:,2], self.zeropoint, self.skylevel, self.scale, self.gal_name, self.name) #create two profiles 1 and 2
		two = Profile('2', data[:,3], data[:,0], data[:,4], self.zeropoint, self.skylevel, self.scale, self.gal_name, self.name)
		self.profile_list.extend((one,two))

	

class Profile(object):
	"""Contains information for profile"""
	def __init__(self, name, I, R, W, zeropoint, skylevel, scale, gal_name=None, cam_name=None):
		"""stores data, converts I to counts per arcsec^2 and R to arcsec"""
		self.name = name
		self.skylevel = float(skylevel)
		self.zeropoint = float(zeropoint)
		self.scale = float(scale)
		self.I = (I - skylevel) / (scale ** 2)
		self.R = R * scale
		self.W = W #assigned when importing ini files
		self.M = self.zeropoint - (2.5*np.log10(self.I))
		self.MW = None
		self.params = lm.Parameters()
		self.fits = {}
		self.sky_fit = {}
		self.confs ={}
		self.gal_name = gal_name
		self.cam_name = cam_name

	def __len__(self):
		return len(self.fits)
		
	def add_ini_params(self, File, names=None):
		if os.path.isfile(File):
			Config = cp.ConfigParser()
			Config.read(File)
			d = {i:ConfigSectionMap(Config, i) for i in Config.sections()}
		else:
			#raise IOError("Ini File not found for %s" % File)
			print "Warning: ini file not found/blank for %s [creating a template]" % File
			f = open(File, 'w')
			Config = cp.ConfigParser()
			Config.add_section('Components')
			Config.set('Components','MeB',22)
			Config.set('Components','ReB',1)
			Config.set('Components','nB',4)
			Config.set('Components','MeD',20)
			Config.set('Components','ReD',8)
			Config.set('Components','nD',1)
			Config.add_section('Ranges')
			Config.set('Ranges', 'fitrange', '[0,60]')
			Config.set('Ranges', 'skyrange', '[35,60]')
			Config.write(f)
			f.close()
			Config = cp.ConfigParser()
			Config.read(File)
			d = {i:ConfigSectionMap(Config, i) for i in Config.sections()}
		comps = d['Components']
		self.fit_range = string_to_list(d['Ranges']['fitrange'])
		self.sky_range = string_to_list(d['Ranges']['skyrange'])
		self.W, self.sky_average = self.weight_method(self.sky_range)
		self.convert_errors()
		if names is None:
			names = ['MeB', 'ReB', 'nB', 'MeD', 'ReD', 'nD']
		for n in names:
			self.params.add(n, value=float(comps[n.lower()]), min=0.01)
		self.params['MeB'].max, self.params['MeD'].max = 40., 40.
		self.params['nB'].max = 8.

				
	def convert_errors(self):
		up = self.W + self.I
		down = -self.W + self.I
		down = [val if val > 0 else up[i] for i, val in enumerate(down)]
		mini = 0.1
		errup = self.zeropoint - (2.5 * np.log10(down))
		errdown = self.zeropoint - (2.5 * np.log10(up))
		self.MW = [self.M - errdown, abs(errup - self.M)]
		for i in range(len(up)):
			if self.MW[0][i] < mini: self.MW[0][i] = mini
			if self.MW[1][i] < mini: self.MW[1][i] = mini

	def preview(self):
		fig = plt.figure()
		fig.set_facecolor('white')
		ax = fig.add_subplot(211)
		# for i,v in enumerate(self.MW[1]):
		# 	if v < 0: print v, self.R[i]
		# 	self.MW[1][i] = abs(self.MW[1][i])
		if self.MW:
			ax.errorbar(self.R, self.M, yerr=self.MW, fmt='b.')
		else: 
			ax.plot(self.R, self.M,'b.')
		ax.set_ylim([31,17])
		ax.axhline(y=self.zeropoint - 2.5*np.log10(np.sqrt(self.sky_average)), linestyle='--')
		ax2 = fig.add_subplot(212)
		sky_ind = [translate_x(self.R, point) for point in self.sky_range]
		ax2.plot(self.R[sky_ind[0]:sky_ind[1]], self.I[sky_ind[0]:sky_ind[1]], 'b.')
		ax2.axhline(y=0, linestyle='--')
		ax2.axhline(y=np.sqrt(self.sky_average), linestyle=':')
		ax.set_title(self.gal_name+self.cam_name+self.name)
		# ax.set_xscale('log')
		plt.show()

	def weight_method(self, skyrange):
		"""calculates weighting by SNR based on sky in skyrange (in arcsec)
		There is a minimum weighting of sqrt(sky level)"""
		if self.cam_name is 'sdss':
			t, G = 54., 7.43
		elif self.cam_name is 'mega':
			t, G = 35., 1.7
		else: t, G = 1.,1.

		# loc = findSky.iterate(self.I, self.R)
		# skyrange = [int(loc[0]), self.fit_range[1]]
		sky_R = [translate_x(self.R, i) for i in skyrange]

		diff = sky_R[1] - sky_R[0]
		# if 0 > diff <= 2:
		# 	warn("Only %i points available for sky analysis in %s. Increase sky range for better results" % (diff, self.cam_name+self.gal_name))
		if diff == 0:
			# warn("No points available for sky analysis in %s between %.1f and %.1f. Increase sky range for better results\n\
			# 	...max_R = %.1f [automatically adjusting range]" % (self.cam_name+self.gal_name, skyrange[0], skyrange[1], self.R[-1]))
			sky_R = [-3, -1]
		sky_I = self.I[sky_R[0]:sky_R[1]] + self.skylevel
		sky_av = np.mean(sky_I)
		signal = (self.I + self.skylevel - sky_av) * t * G
		noise = np.sqrt(signal + (np.std(sky_I * t * G)**2.) + (self.skylevel * t * G))
		mini = np.sqrt(sky_av * t * G)
		if np.isnan(mini):
			warn("The error on averaged sky is below 0. Minimum weight is now tiny")
			mini = 1e-05
		Weighting = noise 
		for i, val in enumerate(Weighting):
			if (val < mini) or np.isnan(val):
				Weighting[i] = mini
		self.sky_var = np.std(sky_I)
		return Weighting / (t * G), sky_av / (t * G)

def import_directory(directory, names=None):
	files = [f.split('_') for f in os.listdir(directory) if fnmatch.fnmatch(f, '*.ascii')] #cam_number_cts.ascii
	galaxies = []
	name_list = []
	for f in files:
		if f[1] not in name_list:
			Gal = Galaxy(f[1])
			Gal.import_file(directory+'\\'+'_'.join(f), f[0])
			ini_name = '_'.join([f[0], f[1]]) + '.ini'
			Gal[f[0]][0].add_ini_params(directory+'\\'+ini_name, names)
			Gal[f[0]][1].add_ini_params(directory+'\\'+ini_name, names)
			galaxies.append(Gal)
			name_list.append(f[1])
		else:
			i = name_list.index(f[1])
			galaxies[i].import_file(directory+'\\'+'_'.join(f), f[0])
			ini_name = '_'.join([f[0], f[1]]) + '.ini'
			galaxies[i][f[0]][0].add_ini_params(directory+'\\'+ini_name, names)
			galaxies[i][f[0]][1].add_ini_params(directory+'\\'+ini_name, names)
	return galaxies, name_list
			
if __name__ == '__main__':
	direct = 'C:\\Users\\User\\code\\Python\\Sersic_Plots'
	gal_list, gal_names = import_directory(direct)
	# for i,v in enumerate(gal_list[0][0][0].MW[1]):
	# 	if v < 0: print v, gal_list[0][0][0].R[i]
	
	for gal in gal_list:
		for cam in gal:
			for prof in cam:
				prof.preview()
