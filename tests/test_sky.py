import storeData as SD
import matplotlib.pyplot as plt
import numpy as np

gal_list, gal_names = SD.import_directory(r'C:\Users\User\Documents\Project work\GalaxyFitting\repository')
G = gal_list[23][1][0]

# stds = np.array([np.std(G.I[::-1][:i+3]) for i, r in enumerate(G.R[-3::-1])])
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax2 = ax.twinx()

# ax.plot(G.R[-3::-1], stds)
# ax2.plot(G.R, G.I, 'b.')
# ax2.set_ylim([-10,20])


# point = G.R[stds[::-1]<10]
# a= np.arange(len(G.R))[stds[::-1]<10][3]
# print G.R[a]

# ax.axvline(x=point[3])
# G.preview()
# plt.show()

Irange = np.linspace(np.min(G.I)-1, G.I[len(G.I)/2], 100)
str_line = lambda x,c: np.ones(x.shape)*c
stds = []
for i in Irange:
	line = str_line(G.R, i)
	if i <= np.min(G.I):
		stds.append(np.nan)
	else:
		limit = range(len(G.I))[G.I <= i]
		error = np.std(G.I[limit[0]:]) / np.sqrt(len(G.I[limit[0]:]))
		stds.append(error)
plt.plot(Irange, stds, 'b.')
plt.show()

