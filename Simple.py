#!/usr/bin/python3

import os
import numpy
import Bio.PDB
from pyrosetta import *
from pyrosetta.toolbox import *
init()

def Get(protein , chain , motif):
	'''
	Get a protein and isolates the desired chain Generates the 
	original.pdb file
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
	NewFile = open('original.pdb' , 'w')
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
	return(NewMotif)

def BLAST(filename1 , filename2):
	'''
	Performs a BLAST alignment between two sequences and prints
	the sequences as well as the percentage of sequence
	similarity
	'''
	seq1 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET = True).get_structure('filename1' , filename1) , aa_only = True)[0].get_sequence()
	seq2 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET = True).get_structure('filename2' , filename2) , aa_only = True)[0].get_sequence()
	alignment = pairwise2.align.globalxx(seq1 , seq2)
	total = alignment[0][4]
	similarity = alignment[0][2]
	percentage = (similarity * 100) / total
	print(seq1)
	print(seq2)
	print('Sequence Similarity: {}%'.format(round(percentage , 3)))

def fixbb(filename , motif , relax_iters , design_iters):
	'''
	Performs the RosettaDesign protocol to change a structure's
	amino acid sequence while maintaining a fixed backbone.
	Generates the structure.pdb file
	'''
	#A - Relax original structure
	pose = pose_from_pdb(filename)
	chain = pose.pdb_info().chain(1)
	scorefxn = get_fa_scorefxn()
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	Rscore_before = scorefxn(pose)
	Rpose_work = Pose()
	Rpose_lowest = Pose()
	Rscores = []
	Rscores.append(Rscore_before)
	for nstruct in range(relax_iters):
		Rpose_work.assign(pose)
		relax.apply(Rpose_work)
		Rscore_after = scorefxn(Rpose_work)
		Rscores.append(Rscore_after)
		if Rscore_after < Rscore_before:
			Rscore_before = Rscore_after
			Rpose_lowest.assign(Rpose_work)
		else:
			continue
	pose.assign(Rpose_lowest)
	RFinalScore = scorefxn(pose)
	#B - Perform fixbb RosettaDesign without the motif residues
	packtask = standard_packer_task(pose)
	for aa in motif:
		packtask.temporarily_set_pack_residue(aa , False)
	packtask = standard_packer_task(pose)
	pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn , packtask)
	Dscore_before = 0
	Dpose_work = Pose()
	Dpose_lowest = Pose()
	Dscores = []
	Dscores.append(Dscore_before)
	for nstruct in range(design_iters):
		Dpose_work.assign(pose)
		pack.apply(Dpose_work)
		Dscore_after = scorefxn(Dpose_work)
		Dscores.append(Dscore_after)
		if Dscore_after < Dscore_before:
			Dscore_before = Dscore_after
			Dpose_lowest.assign(Dpose_work)
		else:
			continue
	pose.assign(Dpose_lowest)
	DFinalScore = scorefxn(pose)
	#C - Output Result
	pose.dump_pdb('structure.pdb')
	#D - Print report
	print('==================== Result Report ====================')
	print('Relax Scores:\n' , Rscores)
	print('Chosen Lowest Score:' , RFinalScore , '\n')
	print('Design Scores:\n' , Dscores)
	print('Chosen Lowest Score:' , DFinalScore , '\n')
	print('BLAST result, comparing the original structure to the designed structure:')
	BLAST(filename , 'structure.pdb')

def main():
	TheProtein = sys.argv[1]
	TheChain = sys.argv[2]
	TheMotif = list(map(int , sys.argv[3:]))
	Motif = Get(TheProtein , TheChain , TheMotif)
	fixbb('original.pdb' , Motif , 50 , 100)

if __name__ == '__main__':
	main()