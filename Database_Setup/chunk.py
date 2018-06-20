import os
import sys

current = os.getcwd()
database = os.listdir(sys.argv[0])
os.chdir(sys.argv[0])

l = database
n = 100
lists = [l[i:i+n] for i in range(0, len(l), n)]

count = 1
for chunk in lists:
	thefile = open('PDB{}.list'.format(count), 'w')
	for item in chunk:
		line = '{}/{}\n'.format(sys.argv[2], item)
		thefile.write(line)
		print(line)
	thefile.close()
	count += 1
os.system('mv *.list ../')
