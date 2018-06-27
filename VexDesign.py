#!/usr/bin/python3

#Import Modules
import sys
import os
import re
import time
import datetime
import subprocess
import random
import requests
import urllib.request
import bs4
import Bio.PDB
from Bio import pairwise2
from pyrosetta import *
from pyrosetta.toolbox import *
init()

#Terminal Text Colours
Black 	= '\x1b[30m'
Red	= '\x1b[31m'
Green	= '\x1b[32m'
Yellow	= '\x1b[33m'
Blue	= '\x1b[34m'
Purple	= '\x1b[35m'
Cyan	= '\x1b[36m'
White	= '\x1b[37m'
Cancel	= '\x1b[0m'

#Welcome Text
print(Green + '\n  ██╗   ██╗ █████╗ ██╗  ██╗\n  ██║   ██║██╔══██╗╚██╗██╔╝\n  ██║   ██║███████║ ╚███╔╝ \n  ╚██╗ ██╔╝██╔══██║ ██╔██╗ \n   ╚████╔╝ ██║  ██║██╔╝ ██╗\n    ╚═══╝  ╚═╝  ╚═╝╚═╝  ╚═╝\n                           \n  ██████╗ ███████╗███████╗██╗ ██████╗ ███╗   ██╗\n  ██╔══██╗██╔════╝██╔════╝██║██╔════╝ ████╗  ██║\n  ██║  ██║█████╗  ███████╗██║██║  ███╗██╔██╗ ██║\n  ██║  ██║██╔══╝  ╚════██║██║██║   ██║██║╚██╗██║\n  ██████╔╝███████╗███████║██║╚██████╔╝██║ ╚████║\n  ╚═════╝ ╚══════╝╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝\n                                                ' + Cancel)
print(Purple + '  ╔═╗┬ ┬┌┬┐┌─┐  ╔╦╗┌─┐┌─┐┬┌─┐┌┐┌  ╔═╗  ╦  ╦┌─┐┌─┐┌─┐┬┌┐┌┌─┐\n  ╠═╣│ │ │ │ │   ║║├┤ └─┐││ ┬│││  ╠═╣  ╚╗╔╝├─┤│  │  ││││├┤ \n  ╩ ╩└─┘ ┴ └─┘  ═╩╝└─┘└─┘┴└─┘┘└┘  ╩ ╩   ╚╝ ┴ ┴└─┘└─┘┴┘└┘└─┘' + Cancel)
print(Yellow + 'Authored by Sari Sabban on 31-May-2017 (sari.sabban@gmail.com)' + Cancel)
print(Cyan + '--------------------------------------------------------------' + Cancel)
#--------------------------------------------------------------------------------------------------------------------------------------
#The Functions

#1 - Extract Motif
def Motif(PDB_ID, Chain, Motif_From, Motif_To):
	'''
	This function downloads a spesific protein from RCSB
	and isolates a specific user defined motif from it
	Generates the motif.pdb file
	'''
	#Get the protein
	os.system('wget http://www.rcsb.org/pdb/files/' + PDB_ID + '.pdb')
	pdb = open(PDB_ID + '.pdb', 'r')
	#Isolate the motif
	Motif = open('motif.pdb', 'w')
	count = 0
	num = 0
	AA2 = None
	for line in pdb:
		if not line.startswith('ATOM'):														#Ignore all lines that do not start with ATOM
			continue
		if not line.split()[4] == Chain:													#Ignore all lines that do not have the specified chain (column 5)
			continue
		try:
			if int(Motif_From) <= int(line.split()[5]) <= int(Motif_To):									#Find residues between the user specified location
				count += 1														#Sequencially number atoms
				AA1 = line[23:27]													#Sequencially number residues
				if not AA1 == AA2:
					num += 1			
				final_line = line[:7] + '{:4d}'.format(count) + line[11:17] + line[17:21] + 'A' + '{:4d}'.format(num) + line[26:]	#Update each line of the motif to have its atoms and residues sequencially labeled, as well as being in chain A
				AA2 = AA1
				Motif.write(final_line)													#Write to new file called motif.pdb
		except:
			continue
	Motif.close()

#2 - Extract Receptor
def Receptor(PDB_ID, Chain):
	'''
	This function isolates a chain from a downloaded .pdb file 
	Generates the receptor.pdb file
	'''
	#Isolate the receptor
	pdb = open(PDB_ID + '.pdb', 'r')
	Receptor = open('receptor.pdb', 'w')
	for line in pdb:
		linesplit = line.split()
		if linesplit[0] == 'ATOM':
			if linesplit[4] == Chain:
				Receptor.write(line)
	Receptor.close()
	#Keep working directory clean, remove the protein's original file
	os.remove(PDB_ID + '.pdb')

#3 - RosettaRelax
def Relax(pose):
	'''
	Relaxes a structure
	Updates the original pose with the relaxed pose
	'''
	scorefxn = get_fa_scorefxn()
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.apply(pose)

#4 - Grafting
def Graft(receptor, motif, scaffold):
	'''
	Grafts a motif onto a protein scaffold structure
	Generates structure.pdb and returns a tuple [0] is the
	residue number where the motif starts and [1] where it ends
	'''
	scorefxn = get_fa_scorefxn()
	print(scorefxn(scaffold))
	mover = pyrosetta.rosetta.protocols.motif_grafting.movers.MotifGraftMover()
	#Setup motif hotspots
	motifpose = pose_from_pdb(motif)
	spots = list()
	for resi in range(motifpose.total_residue()):
		spots.append(str(resi + 1))
	hotspots = ':'.join(spots)
	#Setup grafting mover
	mover.init_parameters(receptor, motif, 1.0, 2, 5, '0:0', '0:0', 'ALA', hotspots, True, False, True, False, False, True, False)
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
	mover.apply(scaffold)
	print(scorefxn(scaffold))
	scaffold.dump_pdb('temp.pdb')
	#Extract grafted structure
	pdb = open('temp.pdb', 'r')
	Structure = open('temp2.pdb', 'w')
	for line in pdb:
		linesplit = line.split()
		if linesplit != []:
			if linesplit[0] == 'ATOM':
				if linesplit[4] == 'B':
					Structure.write(line)
	Structure.close()
	#Keep working directory clean
	os.remove('temp.pdb')
	#Renumber .pdb file starting at 1
	newpdb = open('temp2.pdb', 'r')
	thenewfile = open('grafted.pdb', 'w')
	count = 0
	num = 0
	AA2 = None
	for line in newpdb:
		count += 1														#Sequencially number atoms
		AA1 = line[23:27]													#Sequencially number residues
		if not AA1 == AA2:
			num += 1			
		final_line = line[:7] + '{:4d}'.format(count) + line[11:17] + line[17:21] + 'A' + '{:4d}'.format(num) + line[26:]	#Update each line of the motif to have its atoms and residues sequencially labeled, as well as being in chain A
		AA2 = AA1
		thenewfile.write(final_line)												#Write to new file called motif.pdb
	thenewfile.close()
	os.remove('temp2.pdb')
	#Identify start and finish residue number of inserted motif
	motifpose = pose_from_pdb('motif.pdb')												#Input motif structure as a pose
	graftpose = pose_from_pdb('grafted.pdb')											#Input graft structure as a pose
	MOTIF = motifpose.sequence()													#Get motif sequence
	GRAFT = graftpose.sequence()													#Get graft sequence
	start = GRAFT.index(MOTIF) + 1													#Identify start residue
	finish = start + len(MOTIF) - 1													#Identify end residue
	return((start, finish))														#Return values [0] = Motif_From [1] = Motif_To

#5 - RosettaDesign
class RosettaDesign():
	'''
	This class preforms RosettaDesign either fixed backbone 
	design (fixbb) or flexible backbone design (flxbb).
	It is preferred to perform the design many times and 
	select the best (lowest) scoring structure.
	'''
	def __init__(self):
		pass

	#5.1 - BLAST
	def BLAST(self, filename1, filename2):
		'''
		Performs a BLAST alignment between two sequences and prints
		the sequences as well as the percentage of sequence
		similarity
		'''
		seq1 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET=True).get_structure('filename1', filename1), aa_only=True)[0].get_sequence()
		seq2 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET=True).get_structure('filename2', filename2), aa_only=True)[0].get_sequence()
		alignment = pairwise2.align.globalxx(seq1, seq2)
		total = alignment[0][4]
		similarity = alignment[0][2]
		percentage = (similarity*100)/total
		print(seq1)
		print(seq2)
		print('Sequence Similarity: {}%'.format(round(percentage, 3)))

	#5.2 - Performs fixbb RosettaDesign
	def fixbb(self, filename, relax_iters, design_iters):
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
		#B - Perform fixbb RosettaDesign
		packtask = standard_packer_task(pose)
		pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, packtask)
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
		pose.dump_pdb('fixbb.pdb')
		#D - Print report
		print('==================== Result Report ====================')
		print('Relax Scores:\n', Rscores)
		print('Chosen Lowest Score:', RFinalScore, '\n')
		print('Design Scores:\n', Dscores)
		print('Chosen Lowest Score:', DFinalScore, '\n')
		print('BLAST result, comparing the original structure to the designed structure:')
		RosettaDesign.BLAST(self, filename, 'fixbb.pdb')

	#5.3 - Performs flxbb RosettaDesign
	def flxbb(self, filename, relax_iters, design_iters):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while allowing for a flexible backbone.
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
		#B - Perform flxbb RosettaDesign
		mover = pyrosetta.rosetta.protocols.flxbb.FlxbbDesign()
		Dscore_before = 0
		Dpose_work = Pose()
		Dpose_lowest = Pose()
		Dscores = []
		Dscores.append(Dscore_before)
		for nstruct in range(design_iters):
			Dpose_work.assign(pose)
			mover.apply(Dpose_work)
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
		pose.dump_pdb('flxbb.pdb')
		#D - Print report
		print('==================== Result Report ====================')
		print('Relax Scores:\n', Rscores)
		print('Chosen Lowest Score:', RFinalScore, '\n')
		print('Design Scores:\n', Dscores)
		print('Chosen Lowest Score:', DFinalScore, '\n')
		print('BLAST result, comparing the original structure to the designed structure:')
		RosettaDesign.BLAST(self, filename, 'flxbb.pdb')

	#5.4 - Rebuild loops
	def BDR(self, filename, refine_iters):
		#A - Generate constraints file
		structure = Bio.PDB.PDBParser(QUIET=True).get_structure('{}'.format(filename), filename)
		length = len(structure[0]['A'])
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		Type = ppb.build_peptides(structure, aa_only=False)
		model = Type
		chain = model[0]
		CST = []
		CST.append(0.0)
		for aa in range(1, length+1):
			try:
				residue1 = chain[0]
				residue2 = chain[aa]
				atom1 = residue1['CA']
				atom2 = residue2['CA']
				CST.append(atom1-atom2)
			except:
				pass
		atom = 1
		for cst in CST:
			line = 'AtomPair CA 1 CA '+str(atom)+' GAUSSIANFUNC '+str(cst)+' 1.0\n'
			thefile = open('structure.constraints', 'a')
			thefile.write(line)
			thefile.close()
			atom += 1
		#B - Generate blueprint file (remodeling only large loops)
		dssp = Bio.PDB.DSSP(structure[0], filename)
		SS = []
		SEQ = []
		for ss in dssp:
			if ss[2] == 'G' or ss[2] == 'H' or ss[2] == 'I':
				rename = 'HX'
			elif ss[2] == 'B' or ss[2] == 'E':
				rename = 'EX'
			else:
				rename = 'LX'
			SS.append(rename)
			SEQ.append(ss[1])
		buf = []
		items = []
		l_seen = 0
		for count, (ss, aa) in enumerate(zip(SS, SEQ), 1):
			buf.append((count, aa, ss))
			if 'LX' in {ss, aa}:
				l_seen += 1
				if l_seen >= 3:
					for count, aa, ss in buf:
						line = [str(count), aa, ss, '.' if ss in {'HX', 'EX'} else 'R']
						line = ' '.join(line)
						items.append(line)
					buf.clear()
			else:
				l_seen = 0
				for count, aa, ss in buf:
					line = [str(count), aa, ss, '.']
					line = ' '.join(line)
					items.append(line)
				buf.clear()
		if int(items[-1].split()[0]) != count:
			line = [str(count), aa, ss, '.']
			line = ' '.join(line)
			items.append(line)
		blueprint = open('structure.blueprint', 'a')
		for line in items:
			blueprint.write(line + '\n')
		blueprint.close()
		#C - Run BluePrint mover
		pose = pose_from_pdb(filename)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		secstr = pyrosetta.rosetta.protocols.fldsgn.potentials.SetSecStructEnergies(scorefxn,'structure.blueprint', True)
		secstr.apply(pose)
		BDR = pyrosetta.rosetta.protocols.fldsgn.BluePrintBDR()
		BDR.num_fragpick(200)
		BDR.use_fullmer(True)
		BDR.use_sequence_bias(False)
		BDR.max_linear_chainbreak(0.07)
		BDR.ss_from_blueprint(True)
		BDR.dump_pdb_when_fail('')
		BDR.set_constraints_NtoC(-1.0)
		BDR.use_abego_bias(True)
		#BDR.set_constraint_file('structure.constraints')
		BDR.set_blueprint('structure.blueprint')
		Dscore_before = 0
		Dpose_work = Pose()
		Dpose_lowest = Pose()
		Dscores = []
		Dscores.append(Dscore_before)
		for nstruct in range(refine_iters):
			Dpose_work.assign(pose)
			BDR.apply(Dpose_work)
			relax.apply(Dpose_work)
			Dscore_after = scorefxn(Dpose_work)
			Dscores.append(Dscore_after)
			if Dscore_after < Dscore_before:
				Dscore_before = Dscore_after
				Dpose_lowest.assign(Dpose_work)
			else:
				continue
		pose.assign(Dpose_lowest)
		DFinalScore = scorefxn(pose)
		#D - Output Result
		pose.dump_pdb('remodel.pdb')
		os.remove('structure.constraints')
		os.remove('structure.blueprint')
		#E - Print report
		print('==================== Result Report ====================')
		print('Design Scores:\n', Dscores)
		print('Chosen Lowest Score:', DFinalScore, '\n')

	#5.5 - Find amino acids in wrong layer
	def Layers(self, filename):
		'''
		This function will calculate the solvent-accessible surface area
		(SASA) and the secondary structure for each amino acid within a
		protein and points out the amino acids that are in the wrong layer.
		Returns a list of all positions to be mutated.
		'''
		#Identify SASA and secondary structure for each residue
		parser = Bio.PDB.PDBParser()
		structure = parser.get_structure('{}'.format(filename), filename)
		dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke')
		sasalist = []
		for x in dssp:
			if x[1] == 'A':
				sasa = 129*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'V':
				sasa = 174*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'I':
				sasa = 197*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'L':
				sasa = 201*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'M':
				sasa = 224*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'P':
				sasa = 159*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'Y':
				sasa = 263*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'F':
				sasa = 240*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'W':
				sasa = 285*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'R':
				sasa = 274*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'N':
				sasa = 195*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'C':
				sasa = 167*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'Q':
				sasa = 225*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'E':
				sasa = 223*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'G':
				sasa = 104*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'H':
				sasa = 224*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'K':
				sasa = 236*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'S':
				sasa = 155*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'T':
				sasa = 172*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'D':
				sasa = 193*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			if x[2] == 'G' or x[2] == 'H' or x[2] == 'I':
				ss = 'H'
			elif x[2] == 'B' or x[2] == 'E':
				ss = 'S'
			elif x[2] == 'S' or x[2] == 'T' or x[2] == '-':
				ss = 'L'
			sasalist.append((x[0], x[1], ss, sasa)) #(number, residue, secondary structure, SASA)
		#Identify residues in the wrong layer to mutate
		Resids = []
		SecStr = []
		SASAps = []
		MutPos = []
		Mutate = []
		for n, r, s, a in sasalist:
			if a == 'S' and s == 'L' and (	   r == 'P' or r == 'G' 
							or r == 'N' or r == 'Q'
							or r == 'S' or r == 'T'
							or r == 'D' or r == 'E'
							or r == 'R' or r == 'K'
							or r == 'H'):
				MutPos.append(' ')
			elif a=='B' and s=='L' and (	   r == 'A' or r == 'V'
							or r == 'I' or r == 'L'
							or r == 'F' or r == 'Y'
							or r == 'W' or r == 'G'
							or r == 'N' or r == 'Q'
							or r == 'S' or r == 'T'
							or r == 'P' or r == 'D'
							or r == 'E' or r == 'H'
							or r == 'R'):
				MutPos.append(' ')
			elif a=='C' and s=='L' and (	   r == 'A' or r == 'V'
							or r == 'I' or r == 'L'
							or r == 'P' or r == 'F'
							or r == 'W' or r == 'M'):
				MutPos.append(' ')
			elif a=='S' and s=='H' and (	   r == 'Q' or r == 'E'
							or r == 'K' or r == 'H'):
				MutPos.append(' ')
			elif a=='B' and s=='H' and (	   r == 'A' or r == 'V'
							or r == 'I' or r == 'L'
							or r == 'W' or r == 'Q'
							or r == 'E' or r == 'K'
							or r == 'F' or r == 'M'):
				MutPos.append(' ')
			elif a=='C' and s=='H' and (	   r == 'A' or r == 'V'
							or r == 'I' or r == 'L'
							or r == 'F' or r == 'W'):
				MutPos.append(' ')
			elif a=='S' and s=='S' and (	   r == 'Q' or r == 'T'
							or r == 'Y'):
				MutPos.append(' ')
			elif a=='B' and s=='S' and (	   r == 'A' or r == 'V'
							or r == 'I' or r == 'L'
							or r == 'F' or r == 'Y'
							or r == 'W' or r == 'Q'
							or r == 'T' or r == 'M'):
				MutPos.append(' ')
			elif a=='C' and s=='S' and (	   r == 'A' or r == 'V'
							or r == 'I' or r == 'L'
							or r == 'F' or r == 'W'
							or r == 'M'):
				MutPos.append(' ')
			else:
				MutPos.append('*')
				Mutate.append((n, r, s, a))		
			Resids.append(r)
			SASAps.append(a)
			SecStr.append(s)
		Resids=''.join(Resids)
		SASAps=''.join(SASAps)
		MutPos=''.join(MutPos)
		SecStr=''.join(SecStr)
		print('------------------------------')
		print('{}\n{}\n{}\n{}'.format(Resids, SecStr, SASAps, MutPos))
		return(Mutate)

	#5.6 - Refine structure
	def Refine(self, filename, mutations, refine_iters):
		'''
		This function takes the list of amino acids from the Layers()
		function that are in the wrong layer and mutates the structure by
		changing these position intothe preferred amino acids for the
		respective layer and secondary structure. Then refines the
		structure in an attempt to generate an ideal protein structure.
		Generates the refined.pdb file.
		'''
		pose = pose_from_pdb(filename)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		ideal = pyrosetta.rosetta.protocols.idealize.IdealizeMover()
		#A - Generate a resfile
		resfile = open('structure.res', 'a')
		resfile.write('NATRO\nSTART\n')
		for n, r, s, a in mutations:
			if s == 'L' and a == 'S':
				line = '{} A PIKAA PGNQSTDERKH\n'.format(n)
				resfile.write(line)
			elif s == 'H' and a == 'S':
				line = '{} A PIKAA QEKH\n'.format(n)
				resfile.write(line)
			elif s == 'S' and a == 'S':
				line = '{} A PIKAA QTY\n'.format(n)
				resfile.write(line)
			elif s == 'L' and a == 'B':
				line = '{} A PIKAA AVILFYWGNQSTPDEKR\n'.format(n)
				resfile.write(line)
			elif s == 'H' and a == 'B':
				line = '{} A PIKAA AVILWQEKFM\n'.format(n)
				resfile.write(line)
			elif s == 'S' and a == 'B':
				line = '{} A PIKAA AVILFYWQTM\n'.format(n)
				resfile.write(line)
			elif s == 'L' and a == 'C':
				line = '{} A PIKAA AVILPFWM\n'.format(n)
				resfile.write(line)
			elif s == 'H' and a == 'C':
				line = '{} A PIKAA AVILFWM\n'.format(n)
				resfile.write(line)
			elif s == 'S' and a == 'C':
				line = '{} A PIKAA AVILFWM\n'.format(n)
				resfile.write(line)
		resfile.close()
		#B - Refinement
		pack = standard_packer_task(pose)
		pack.temporarily_fix_everything()
		pyrosetta.rosetta.core.pack.task.parse_resfile(pose, pack, 'structure.res')
		for n, r, s, a in mutations:
			x = pose.residue(n).name()
			if x == 'CYS:disulphide':
				continue
			else:
				pack.temporarily_set_pack_residue(n, True) 
		print(pack)
		pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, pack)
		Dscore_before = 0
		Dpose_work = Pose()
		Dpose_lowest = Pose()
		Dscores = []
		Dscores.append(Dscore_before)
		for nstruct in range(refine_iters):
			Dpose_work.assign(pose)
			pack.apply(Dpose_work)
			ideal.apply(pose)
			relax.apply(Dpose_work)
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
		os.remove('structure.res')
		#D - Print report
		print('==================== Result Report ====================')
		print('Refine Scores:\n', Dscores)
		print('Chosen Lowest Score:', DFinalScore, '\n')
		RosettaDesign.BLAST(self, sys.argv[2], 'structure.pdb')
		RosettaDesign.Layers(self, 'structure.pdb')

	#5.7 - Preforms flxbb RosettaDesign for the whole protein except the motif
	def motif_fixbb(self, filename, Motif_From, Motif_To, relax_iters, design_iters):
		'''
		Applies RosettaDesign with a fixed back bone to 
		change the structure's amino acids (one layer at
		a time - like in the Design_Layer method) except
		for a desired continuous motif sequence while
		maintaining the same backbone.
		Just updates the pose with the new structure
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
		#B - Perform flxbb RosettaDesign without the motif residues
		packtask = standard_packer_task(pose)
		#Identify motif residues
		Motif = list(range(int(Motif_From), int(Motif_To) + 1))
		#Prevent motif residues from being designed
		for aa in Motif:
			packtask.temporarily_set_pack_residue(aa, False)
		pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, packtask)
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
		print('Relax Scores:\n', Rscores)
		print('Chosen Lowest Score:', RFinalScore, '\n')
		print('Design Scores:\n', Dscores)
		print('Chosen Lowest Score:', DFinalScore, '\n')
		print('BLAST result, comparing the original structure to the designed structure:')
		RosettaDesign.BLAST(self, filename, 'structure.pdb')

	#5.8 - Preforms fixbb RosettaDesign for the whole protein except the motif
	def motif_flxbb(self, filename, Motif_From, Motif_To, relax_iters, design_iters):
		'''
		Applies RosettaDesign with a flexible back bone to
		change the structure's amino acids (one layer at a
		time - like in the Design_Layer method) except for
		a desired continuous motif sequence while maintaining
		the same backbone.
		Just updates the pose with the new structure
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
		#B - Perform flxbb RosettaDesign without the motif residues
		packtask = standard_packer_task(pose)
		#Identify motif residues
		Motif = list(range(int(Motif_From), int(Motif_To) + 1))
		#Prevent motif residues from being designed
		for aa in Motif:
			packtask.temporarily_set_pack_residue(aa, False)
		mover = pyrosetta.rosetta.protocols.flxbb.FlxbbDesign()################# Must find a way to pass the packtask to this mover in order to prevent the redesigning of the motif sequence
		Dscore_before = 0
		Dpose_work = Pose()
		Dpose_lowest = Pose()
		Dscores = []
		Dscores.append(Dscore_before)
		for nstruct in range(design_iters):
			Dpose_work.assign(pose)
			mover.apply(Dpose_work)
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
		print('Relax Scores:\n', Rscores)
		print('Chosen Lowest Score:', RFinalScore, '\n')
		print('Design Scores:\n', Dscores)
		print('Chosen Lowest Score:', DFinalScore, '\n')
		print('BLAST result, comparing the original structure to the designed structure:')
		RosettaDesign.BLAST(self, filename, 'structure.pdb')

#6 - Fragment Generation and Identification
def Fragments(filename, username):
	'''
	Submits the pose to the Robetta server
	(http://www.robetta.org) for fragment generation that are
	used for the Abinitio folding simulation. Then measures the
	RMSD for each fragment at each position and chooses the
	lowest RMSD. Then averages out the lowest RMSDs. Then plots
	the lowest RMSD fragment for each positon.
	Generates the 3-mer file, the 9-mer file, the PsiPred file,
	the RMSD vs Position PDF plot with the averaged fragment
	RMSD printed in the plot
	'''
	#Make the 3-mer and 9-mer fragment files and the PSIPRED file using the Robetta server
	pose = pose_from_pdb(filename)
	sequence = pose.sequence()
	#Post
	web = requests.get('http://www.robetta.org/fragmentsubmit.jsp')
	payload = {
		'UserName':username,
		'Email':'',
		'Notes':'structure',
		'Sequence':sequence,
		'Fasta':'',
		'Code':'',
		'ChemicalShifts':'',
		'NoeConstraints':'',
		'DipolarConstraints':'',
		'type':'submit'}
	session = requests.session()
	response = session.post('http://www.robetta.org/fragmentsubmit.jsp', data=payload, files=dict(foo='bar'))		
	for line in response:
		line = line.decode()
		if re.search('<a href="(fragmentqueue.jsp\?id=[0-9].*)">', line):
			JobID = re.findall('<a href="(fragmentqueue.jsp\?id=[0-9].*)">', line)
	JobURL = 'http://www.robetta.org/' + JobID[0]
	#Check
	ID = JobID[0].split('=')
	print('Job ID: ' + str(ID[1]))
	while True:
		Job = urllib.request.urlopen(JobURL)
		jobdata = bs4.BeautifulSoup(Job, 'lxml')
		status = jobdata.find('td', string='Status: ').find_next().text
		if status == 'Complete':
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M'), 'Status:', status)
			break
		else:
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M'), 'Status:', status)
			time.sleep(1800)
			continue
	#Download
	sequence = pose.sequence()
	fasta = open('structure.fasta', 'w')
	fasta.write(sequence)
	fasta.close()
	time.sleep(1)
	os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/aat000_03_05.200_v1_3')
	os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/aat000_09_05.200_v1_3')
	os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/t000_.psipred_ss2')
	os.rename('aat000_03_05.200_v1_3', 'frags.200.3mers')
	os.rename('aat000_09_05.200_v1_3', 'frags.200.9mers')
	os.rename('t000_.psipred_ss2', 'pre.psipred.ss2')
	#Calculate the best fragment's RMSD at each position
	frag = open('frags.200.9mers', 'r')
	rmsd = open('temp.dat', 'w')
	for line in frag:
		if line.lstrip().startswith('position:'):
			line = line.split()
			size = line[1]
	frag.close()
	count = 0
	for x in range (int(size)):
		count +=1
		#Get the pose and make a copy of it to apply changes to
		pose_copy = pyrosetta.Pose()
		pose_copy.assign(pose)
		#Setup frame list
		frames = pyrosetta.rosetta.core.fragment.FrameList()
		#Setup the 9-mer fragment (9-mer is better than 3-mer for this analysis)
		fragset = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
		fragset.read_fragment_file('frags.200.9mers')
		fragset.frames(count, frames)
		#Setup the MoveMap
		movemap = MoveMap()
		movemap.set_bb(True)
		#Setup and apply the fragment inserting mover
		for frame in frames:
			for frag_num in range(1, frame.nr_frags() + 1 ):
				frame.apply(movemap, frag_num, pose_copy)
				#Measure the RMSD difference between the original pose and the new changed pose (the copy)
				RMSD = rosetta.core.scoring.CA_rmsd(pose, pose_copy)
				print(RMSD, '\t', count)
				rmsd.write(str(RMSD) + '\t' + str(count) + '\n')
				#Reset the copy pose to original pose
				pose_copy.assign(pose)
	rmsd.close()
	#Analyse the RMSD file to get the lowest RMSD for each position
	data = open('RMSDvsPosition.dat', 'w')
	lowest = {} 												#Mapping group number -> lowest value found
	for line in open('temp.dat'):
		parts = line.split()
		if len(parts) != 2:										#Only lines with two items on it
			continue
		first = float(parts[0])
		second = int(parts[1])
		if first == 0: 											#Skip line with 0.0 RMSD (this is an error from the 9-mer fragment file). I don't know why it happens
			continue
		if second not in lowest:
			lowest[second] = first
		else:
			if first < lowest[second]:
				lowest[second] = first
	for position, rmsd in lowest.items():
		#print(str(rmsd) + '\t' + str(position))
		data.write(str(position) + '\t' + str(rmsd) + '\n')
	data.close()
	#Calculate the average RMSD of the fragments
	data = open('RMSDvsPosition.dat', 'r')
	value = 0
	for line in data:
		line = line.split()
		RMSD = float(line[1])
		value = value + RMSD
		count = int(line[0])
	Average_RMSD = round(value / count, 2)
	#Plot the results
	gnuplot = open('gnuplot_sets', 'w')
	gnuplot.write("""
	reset\n
	set terminal postscript\n
	set output './plot_frag.pdf'\n
	set encoding iso_8859_1\n
	set term post eps enh color\n
	set xlabel 'Position'\n
	set ylabel 'RMSD (\\305)'\n
	set yrange [0:]\n
	set xrange [0:]\n
	set xtics auto\n
	set xtics rotate\n
	set grid front\n
	unset grid\n
	set title 'Fragment Quality'\n
	set key off\n
	set boxwidth 0.5\n
	set style fill solid\n
	set label 'Average RMSD = {}' at graph 0.01, graph 0.95 tc lt 7 font 'curior 12'\n
	plot 'RMSDvsPosition.dat' with boxes\n
	exit
	""".format(str(Average_RMSD)))
	gnuplot.close()
	os.system('gnuplot < gnuplot_sets')
	os.remove('gnuplot_sets')
	os.remove('temp.dat')
	return(Average_RMSD)
#--------------------------------------------------------------------------------------------------------------------------------------
#List of all functions and their arguments
'''
1   -	Motif(Protein, Chain, Motif_from, Motif_to)
2   -	Receptor(Protein, RecChain)
3   -	Relax(pose)
4   -	MotifPosition = Graft('receptor.pdb', 'motif.pdb', pose)
5   -	RD = RosettaDesign()
5.4 -	RD.motif_fixbb('grafted.pdb', MotifPosition[0], MotifPosition[1], 50, 100)
5.5 -	RD.Refine('fixbb.pdb', RD.Layers('fixbb.pdb'), 50)
6   -	Fragments(pose)
'''
#--------------------------------------------------------------------------------------------------------------------------------------
def protocol():
	#User inputs
	Protein		= sys.argv[1]
	RecChain	= sys.argv[2]
	Chain		= sys.argv[3]
	Motif_from	= sys.argv[4]
	Motif_to	= sys.argv[5]
	Scaffold	= sys.argv[6]
	UserName	= sys.argv[7]
	#1. Import scaffold
	pose = pose_from_pdb(Scaffold)
	#2. Isolate motif
	Motif(Protein, Chain, Motif_from, Motif_to)
	#3. Isolate receptor
	Receptor(Protein, RecChain)
	#4. Graft motif onto scaffold
	MotifPosition = Graft('receptor.pdb', 'motif.pdb', pose)
	#5. Sequence design the structure around the motif
	RD = RosettaDesign()
	RD.motif_fixbb('grafted.pdb', MotifPosition[0], MotifPosition[1], 50, 100)
	RD.Refine('fixbb.pdb', RD.Layers('fixbb.pdb'), 50)
	#6. Generate to test fragment quality and predict the Abinitio folding simulation success
	Fragments('structure.pdb', UserName)

if __name__ == '__main__':
	protocol()