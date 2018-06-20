import os

current = os.getcwd()
database = os.listdir('Database')
os.chdir('Database')

l = database
n = 100
lists = [l[i:i+n] for i in range(0, len(l), n)]

count = 1
for chunk in lists:
	thefile = open('PDB{}.list'.format(count), 'w')
	for item in chunk:
		line = '/fefs1/generic/ssabban/Database/{}\n'.format(item)
		thefile.write(line)
		print(line)
	thefile.close()
	count += 1
os.system('mv *.list ../')
