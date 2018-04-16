#!/usr/bin/python3

import os , numpy , Bio.PDB
from pyrosetta import *
from pyrosetta.toolbox import *
init()

TheProtein = sys.argv[1]
TheChain = sys.argv[2]
TheMotif = list(map(int , sys.argv[3:]))
#--------------------------------------------------------------------------
def SASA(filename):
	'''
	Calculates the different layers (Surface, Boundary, Core) of a 
	structure according its SASA (solvent-accessible surface area)
	Returns three lists Surface amino acids = [0] , Boundary amino 
	acids = [1] , Core amino acids = [2]
	'''
	parser = Bio.PDB.PDBParser()
	structure = parser.get_structure('X' , filename)
	dssp = Bio.PDB.DSSP(structure[0] , filename , acc_array = 'Wilke')
	lis = list()
	count = 0
	for x in dssp:
		if   x[1] == 'A' : sasa = 129 * (x[3])
		elif x[1] == 'V' : sasa = 174 * (x[3])
		elif x[1] == 'I' : sasa = 197 * (x[3])
		elif x[1] == 'L' : sasa = 201 * (x[3])
		elif x[1] == 'M' : sasa = 224 * (x[3])
		elif x[1] == 'P' : sasa = 159 * (x[3])
		elif x[1] == 'Y' : sasa = 263 * (x[3])
		elif x[1] == 'F' : sasa = 240 * (x[3])
		elif x[1] == 'W' : sasa = 285 * (x[3])
		elif x[1] == 'R' : sasa = 274 * (x[3])
		elif x[1] == 'C' : sasa = 167 * (x[3])
		elif x[1] == 'N' : sasa = 195 * (x[3])
		elif x[1] == 'Q' : sasa = 225 * (x[3])
		elif x[1] == 'E' : sasa = 223 * (x[3])
		elif x[1] == 'G' : sasa = 104 * (x[3])
		elif x[1] == 'H' : sasa = 224 * (x[3])
		elif x[1] == 'K' : sasa = 236 * (x[3])
		elif x[1] == 'S' : sasa = 155 * (x[3])
		elif x[1] == 'T' : sasa = 172 * (x[3])
		elif x[1] == 'D' : sasa = 193 * (x[3])
		lis.append((x[2] , sasa))
	surface = list()
	boundary = list()
	core = list()
	count = 0
	for x , y in lis:
		count = count + 1
		if y <= 25 and (x == '-' or x == 'T' or x == 'S'):
			core.append(count)
		elif 25 < y < 40 and (x == '-' or x == 'T' or x == 'S'):
			boundary.append(count)
		elif y >= 40 and (x == '-' or x == 'T' or x == 'S'):
			surface.append(count)
		elif y <= 15 and (x == 'G' or x == 'H' or x == 'I'):
			core.append(count)
		elif 15 < y < 60 and (x == 'G' or x == 'H' or x == 'I'):
			boundary.append(count)
		elif y >= 60 and (x == 'G' or x == 'H' or x == 'I'):
			surface.append(count)
		elif y <= 15 and (x == 'B' or x == 'E'):
			core.append(count)
		elif 15 < y < 60 and (x == 'B' or x == 'E'):
			boundary.append(count)
		elif y >= 60 and (x == 'B' or x == 'E'):
			surface.append(count)	
	return(surface , boundary , core)

def Get(protein , chain , motif):
	'''
	Get a protein and isolates the desired chain Generates the 
	Original.pdb file
	'''
	# A - Get Structure
	pose = pose_from_rcsb(protein)
	# B - Clean Structure
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
			final_line = line[:7] + '{:4d}'.format(count) \
			+ line[11:17] + line[17:21] + 'A' + \
			'{:4d}'.format(num) + line[26:]
			AA2 = AA1
			NewFile.write(final_line)
	NewFile.close()
	os.remove('temp.pdb')
	NewMotif = [numpy.absolute(AA - FirstAA) for  AA in motif]
	# C - Generate A Resfile of The Surface
	Surface = SASA('Original.pdb')[0]
	surface = [x for x in Surface if x not in NewMotif]
	Aminos = ['ALA' , 'CYS' , 'ASP' , 'GLU' , 'PHE' , 
		  'GLY' , 'HIS' , 'HIS_D' , 'ILE' , 'LYS' , 
		  'LEU' , 'MET' , 'ASN' , 'PRO' , 'GLN' , 
		  'ARG' , 'SER' , 'THR' , 'VAL' , 'TRP' , 'TYR']
	resfile = open('resfile.res' , 'a')
	resfile.write('ALLAA\nSTART\n')
	for aa  in surface:
		x = pose.residue(aa).name()
		if x == 'CYS:disulphide':
			continue
		else:
			newAminos = list()
			for res in Aminos:
				if res == x:
					continue
				else:
					newAminos.append(res)
		newAminos =['A' if x == 'ALA' else x for x in newAminos]
		newAminos =['C' if x == 'CYS' else x for x in newAminos]
		newAminos =['D' if x == 'ASP' else x for x in newAminos]
		newAminos =['E' if x == 'GLU' else x for x in newAminos]
		newAminos =['F' if x == 'PHE' else x for x in newAminos]
		newAminos =['G' if x == 'GLY' else x for x in newAminos]
		newAminos =['H' if x == 'HIS' else x for x in newAminos]
		newAminos =['' if x == 'HIS_D' else x for x in newAminos]
		newAminos =['I' if x == 'ILE' else x for x in newAminos]
		newAminos =['K' if x == 'LYS' else x for x in newAminos]
		newAminos =['L' if x == 'LEU' else x for x in newAminos]
		newAminos =['M' if x == 'MET' else x for x in newAminos]
		newAminos =['N' if x == 'ASN' else x for x in newAminos]
		newAminos =['P' if x == 'PRO' else x for x in newAminos]
		newAminos =['Q' if x == 'GLN' else x for x in newAminos]
		newAminos =['R' if x == 'ARG' else x for x in newAminos]
		newAminos =['S' if x == 'SER' else x for x in newAminos]
		newAminos =['T' if x == 'THR' else x for x in newAminos]
		newAminos =['V' if x == 'VAL' else x for x in newAminos]
		newAminos =['W' if x == 'TRP' else x for x in newAminos]
		newAminos =['Y' if x == 'TYR' else x for x in newAminos]
		Amins = ''.join(newAminos)
		TheLine = '{} A PIKAA {}\n'.format(str(aa) , str(Amins))
		resfile.write(TheLine)
	resfile.close()
	return(NewMotif , surface)

def Design(filename , motif , surface , resfile):
	'''
	Preforms the protein design Generates the Designed.pdb file
	'''
	pose = pose_from_pdb(filename)
	# A - Relax Original Structure
	#pyrosetta.rosetta.protocols.moves.AddPyMOLObserver(pose , False)
	scorefxn = get_fa_scorefxn()
	relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn)
	relax.apply(pose)
	# B - Preform Design
	for iteration in range(3):
		pack = standard_packer_task(pose)
		for aa in motif:
			x = pose.residue(aa).name()
			if x == 'CYS:disulphide':
				continue
			else:
				pack.temporarily_set_pack_residue(aa , \
				False)
		mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn , pack)
		mover.apply(pose)
	# C - Redesign Surface
	for x in range(3):
		pack = standard_packer_task(pose)
		pack.temporarily_fix_everything()
		pyrosetta.rosetta.core.pack.task.parse_resfile(pose , \
		pack , resfile)
		for aa  in surface:
			pack.temporarily_set_pack_residue(aa , True)
		mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn , pack)
		mover.apply(pose)
	# D - Relax Designed Structure
	relax.apply(pose)
	# E - Output Result
	pose.dump_pdb('Designed.pdb')
#--------------------------------------------------------------------------
NewMotif , Surface = Get(TheProtein , TheChain , TheMotif)
Design('Original.pdb' , NewMotif , Surface , 'resfile.res')
os.remove('resfile.res')
