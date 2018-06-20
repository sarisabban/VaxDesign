import os
import sys
import tqdm
from pyrosetta import *
from pyrosetta.toolbox import *
init()

def relax(filename):
	pose = pose_from_pdb(filename)
	scorefxn = get_fa_scorefxn()
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.apply(pose)
	pose.dump_pdb('_{}'.format(filename))

os.mkdir('Relaxed')
current = os.getcwd()
database = os.listdir(sys.argv[1])
os.chdir(sys.argv[1])
for i in tqdm.tqdm(database):
	relax(i)
	os.system('mv _* ../Relaxed')
