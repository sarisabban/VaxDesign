import os , numpy
from pyrosetta import *
from pyrosetta.toolbox import *
init()

TheProtein = sys.argv[1]
TheChain = sys.argv[2]
TheMotif = list(map(int , sys.argv[3:]))
#------------------------------------------------------------------------------
def Get(protein , chain , motif):
	pose = pose_from_rcsb(protein)
	FirstAA = int(pose.pdb_info().pose2pdb(1).split()[0]) - 1
	TheFile = open('{}.clean.pdb'.format(protein) , 'r')
	for line in TheFile:
		try:
			if line.split()[4] == chain:
				NewFile = open('temp.pdb' , 'a')
				NewFile.write(line)
				NewFile.close()
		except:
			pass
	os.remove('{}.clean.pdb'.format(protein))
	os.remove('{}.pdb'.format(protein))
	pdb = open('temp.pdb' , 'r')
	NewFile = open('Original.pdb' , 'w')
	count = 0
	num = 0
	AA2 = None
	for line in pdb:
		if not line.startswith('ATOM'):
			continue
		else:
			count += 1
			AA1 = line[23:27]
			if not AA1 == AA2:
				num += 1
			final_line = line[:7] + '{:4d}'.format(count) + line[11:17] + line[17:21] + 'A' + '{:4d}'.format(num) + line[26:]
			AA2 = AA1
			NewFile.write(final_line)
	NewFile.close()
	os.remove('temp.pdb')
	NewMotif = [numpy.absolute(AA - FirstAA) for  AA in motif]
	return(NewMotif)

def Design(filename , motif):
	pose = pose_from_pdb(filename)
	# A - Relax Original Structure
#	pyrosetta.rosetta.protocols.moves.AddPyMOLObserver(pose , False)
	scorefxn = get_fa_scorefxn()
	relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn)
	relax.apply(pose)
	# B - Preform Design
	for iteration in range(3):
		pack = standard_packer_task(pose)
		mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn , pack)
		for aa in motif:
			x = pose.residue(aa).name()
			if x == 'CYS:disulphide':
				continue
			else:
				pack.temporarily_set_pack_residue(aa , False)
		mover.apply(pose)
	# C - Relax Designed Structure
	relax.apply(pose)
	# D - Output Result
	pose.dump_pdb('Designed.pdb')
#------------------------------------------------------------------------------
NewMotif = Get(TheProtein , TheChain , TheMotif)
Design('Original.pdb' , NewMotif)