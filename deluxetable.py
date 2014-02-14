"""
-takes data in form of array, headers, justs=None, multicolumn if duplicate, and errors for each column

write preamble(len(headers))
write tableheader
write data
write footer


create empty array of dimensions
have: columns in 1D array float format, up and down errors in 1D array float format
insert into function as add_col(column, sigfig, errorup=None, errordown=None)

if errorup=None:
	use round_sigfig(x,n)
	for each value in column:
		round and place in new array
		place this array into empty array

elif errordown=None and errorup !=None:
	use round_sigfig_error(x, err, n)
	for each value in column:
		round with one error
		place into empty array
elif errordown:
	use use round_sigfig_error2(x, err, err2 n)
	for each value in column:
		round with two errors
		place into empty array


"""
import sigfig as sf
import numpy as np
from copy import deepcopy

def cols_to_array(*columns):
	return np.vstack([i for i in columns]).T

class Table(object):
	"""Latex table stuff"""
	def __init__(self, name, justs):
		self.justs = justs
		self.array = None
		self.header = None
		self.header_columns = []
		self.name = name

	def add_header(self, headers):
		add = lambda x: '\\colhead{%s}' % x
		self.header = "\\tablehead{%s} \n" % ' & '.join(map(add, headers)) 

	def add_column(self, column, sigfig=None, errorup=None, errordown=None, percent=False):
		if sigfig is None:
			temp = np.array(column).reshape((len(temp), 1))
		else:
			temp = []
			for i, v in enumerate(column):
				if (errordown is not None) and (errorup is not None):
					rounded = sf.round_sig_error2(v, errorup[i], errordown[i], sigfig)
					string = "$%s ^{+%s}_{-%s}" % rounded
				elif (errordown is None) and (errorup is not None):
					rounded = sf.round_sig_error(v, errorup[i], sigfig)
					string = "$%s \\pm %s" % rounded
				elif (errorup is None):
					rounded = sf.round_sig(v, sigfig)
					string = "$%s" % rounded
				if percent:
					string += "\\%$"
				else: string += "$"
				temp.append(string)
			temp = np.array(temp).reshape((len(temp), 1))
		if self.array is None: 
			self.array = deepcopy(temp)
		else:
			self.array = np.hstack((self.array, temp))

	def add_row(self, row, sigfig, errorup=None, errordown=None, percent=False):
		if sigfig is None:
			temp = np.array(column).reshape((1, len(temp))) 
		else:
			temp = []
			for i, v in enumerate(row):
				if (errordown is not None) and (errorup is not None):
					rounded = sf.round_sig_error2(v, errorup[i], errordown[i], sigfig)
					string = "$%s ^{+%s}_{-%s}" % rounded
				elif (errordown is None) and (errorup is not None):
					rounded = sf.round_sig_error(v, errorup[i], sigfig)
					string = "$%s \\pm %s" % rounded
				elif (errorup is None):
					rounded = sf.round_sig(v, sigfig)
					string = "$%s" % rounded
				if percent:
					string += "\\%$"
				else: string += "$"
				temp.append(string)
			temp = np.array(temp).reshape((1, len(temp)))
		if self.array is None: 
			self.array = deepcopy(temp)
		else:
			self.array = np.vstack((self.array, temp))

	def add_header_column(self, column):
		self.header_columns.append(column)

	def print_data(self, filename, execute=False):
		with open(filename, 'w') as f:
			if execute:
				f.write("\\documentclass{aastex}\n")
				f.write("\\begin{document}\n")
			f.write("\\begin{deluxetable}{%s} \n\\tablecolumns{%i}\n" % (self.justs, len(self.justs)))
			f.write("\\tablewidth{0pt}\n")
			if self.name:
				f.write("\\tablecaption{%s\n\\label{tab:%s}}\n" % (self.name, self.name))
			f.write("\\centering\n")
			if self.header:
				f.write(self.header)
			f.write("\\startdata\n")
			
			for i in range(len(self.array)):
				preline = [str(j[i]) for j in self.header_columns]
				preline = ' & '.join(preline) + ' & '
				line = preline + ' & '.join(self.array[i]) + ' \\\\\n'
				f.write(line)
			f.write("\\enddata \n\\end{deluxetable}\n")
			if execute: f.write("\\end{document}")


if __name__ == '__main__':
	r = lambda x: np.random.random((x))
	col1 = r(10)*10
	col2 = r(10)
	col3 = r(10)
	col4 = [1,2,3,4,5,6,7,8,9,10]

	T = Table('cc')
	T.add_header(('a','b'))
	T.add_column(col1, 1, col2, col3)
	T.add_column(col4, 1)
	# T.columns.insert(0, ['s']*10)
	T.print_data('mydata.tex')
