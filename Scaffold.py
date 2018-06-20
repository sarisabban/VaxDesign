#!/usr/bin/python3

import os
import time
from pyrosetta import *
from pyrosetta.toolbox import *
init()

def Motif(PDB_ID, Chain, Motif_From, Motif_To):
	'''
	This function downloads a spesific protein from RCSB and isolates
	a specific user defined motif from it generating the motif.pdb file
	'''
	os.system('wget http://www.rcsb.org/pdb/files/{}.pdb'.format(PDB_ID))
	pdb = open('{}.pdb'.format(PDB_ID), 'r')
	Motif = open('motif.pdb', 'w')
	count = 0
	num = 0
	AA2 = None
	for line in pdb:
		if not line.startswith('ATOM'):
			continue
		if not line.split()[4] == Chain:
			continue
		try:
			if int(Motif_From) <= int(line.split()[5]) <= int(Motif_To):
				count += 1
				AA1 = line[23:27]
				if not AA1 == AA2:
					num += 1			
				final_line = line[:7]+'{:4d}'.format(count)+line[11:17]+line[17:21]+'A'+'{:4d}'.format(num)+line[26:]
				AA2 = AA1
				Motif.write(final_line)
		except:
			continue
	Motif.close()

def Receptor(PDB_ID, Chain):
	'''
	This function isolates a chain from a downloaded .pdb file and
	generates the receptor.pdb file
	'''
	pdb = open('{}.pdb'.format(PDB_ID), 'r')
	Receptor = open('receptor.pdb', 'w')
	for line in pdb:
		linesplit = line.split()
		if linesplit[0] == 'ATOM':
			if linesplit[4] == Chain:
				Receptor.write(line)
	Receptor.close()
	os.remove('{}.pdb'.format(PDB_ID))

def MotifGraft(receptor, motif, scaffold, RMSD):
	'''
	This script performs backbone motif grafting. It attempts to
	graft the motif onto all proteins in the database directory
	and returns all the successfully grafted structures the script
	returns a tuple [0] is the residue number where the motif starts
	and [1] where it ends
	'''
	scorefxn = get_fa_scorefxn()
	#Get motif hotspots
	motifpose = pose_from_pdb(motif)
	spots = list()
	for resi in range(motifpose.total_residue()):
		spots.append(str(resi+1))
	hotspots = ':'.join(spots)
	#Perform motif grafting
	mover = pyrosetta.rosetta.protocols.motif_grafting.movers.MotifGraftMover()
	mover.init_parameters(receptor, motif, RMSD, 2, 5, '0:0', '0:0', 'ALA', hotspots, True, False, True, False, False, True, False)
	'''
	context_structure				= 'receptor.pdb'
	motif_structure					= 'motif.pdb'
	RMSD_tolerance					= '2.0'
	NC_points_RMSD_tolerance			= '2'
	clash_score_cutoff				= '5'
	combinatory_fragment_size_delta			= '0:0'
	max_fragment_replacement_size_delta		= '0:0'
	clash_test_residue				= 'ALA'
	hotspots					= '1:2:3:4:5:6:7:8:9:10'
	full_motif_bb_alignment				= 'True'
	allow_independent_alignment_per_fragment	= 'False'
	graft_only_hotspots_by_replacement		= 'True'
	only_allow_if_N_point_match_aa_identity		= 'False'
	only_allow_if_C_point_match_aa_identity		= 'Flase'
	revert_graft_to_native_sequence			= 'True'
	allow_repeat_same_graft_output			= 'False'
	'''
	pose = pose_from_pdb(scaffold)
	mover.apply(pose)
	print(scorefxn(pose))

def main():
	#User inputs
	Protein		= sys.argv[1]
	RecChain	= sys.argv[2]
	Chain		= sys.argv[3]
	Motif_From	= sys.argv[4]
	Motif_To	= sys.argv[5]
	Directory	= sys.argv[6]

	hotspots = Motif(Protein, Chain, Motif_From, Motif_To)
	Receptor(Protein, RecChain)

	os.mkdir('Scaffolds')
	current = os.getcwd()
	database = os.listdir(Directory)
	os.chdir(Directory)
	print('\x1b[33m'+'Motif Grafting'+'\x1b[0m')
	for scaffold in database:
		try:
			MotifGraft('../receptor.pdb', '../motif.pdb', scaffold, 1.0)
			os.system('cp {} ../Scaffolds'.format(scaffold))
		except:
			continue
	os.remove('../receptor.pdb')
	os.remove('../motif.pdb')

if __name__ == '__main__':
	main()