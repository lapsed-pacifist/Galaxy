import shutil

with open('repository\\exclusions.txt', 'r') as f:
	exclusions = f.read().split('\n')
for i in exclusions:
	for j in ['mega', 'sdss']:
		shutil.move('repository\\'+ j +'_' + i + '_cts.ascii', 'repository\\bad\\')
		shutil.move('repository\\'+ j +'_' + i + '.ini', 'repository\\bad\\')