#!/usr/bin/python3

'''
# VexDesign
A script that autonomously designs a vaccine. Authored by Sari Sabban on 31-May-2017 (sari.sabban@gmail.com).

## Requirements:
1. Make sure you install [PyRosetta](http://www.pyrosetta.org) as the website describes.
2. Use the following commands (in GNU/Linux) to install all nessesary programs and Python libraries for this script to run successfully:

`sudo apt update && sudo apt install python3-bs4 python3-biopython python3-lxml pymol DSSP gnuplot -y`

## How To Use:
1. Use the following command to run the script:

`python3 VaxDesign.py PDBID RCHAIN CHAIN FROM TO`

* PDBID = The protein's [Protein Data Bank](https://www.rcsb.org) identification name
* RCHAIN = The chain where your receptor resides within the protein .pdb file
* CHAIN = The chain where your target site resides (not part of the receptor) within the protein .pdb file
* FROM = The start of your target site
* TO = The end of your target site

Example:

`python3 VaxDesign.py 2y7q A B 420 429`

VaxDesign1.py works with PyRosetta4 python 3.6 release 176 and previous.

VaxDesign2.py works with PyRosetta4 python 3.6 release 177 onwards.

2. Calculation time is about 12 hours on a normal desktop computer for each structure.
3. Access to the internet is a requirement since the script will be sending and retrieving data from some servers.
4. Use this [Rosetta Abinitio](https://github.com/sarisabban/RosettaAbinitio) script to simulate the folding of the final designed vaccine's protein structure. An HPC (High Preformance Computer) and the original C++ [Rosetta](https://www.rosettacommons.org/) are required for this step.

## Description
This script autonomously designs a vaccine from a user specified target site. This is not artificial intellegance, you cannot just ask the the script to design "A" vaccine, you must understand what target site you want to develop antibodies against (make a liturature search and understand your disease and target site), then supply this target site to the script to build a protein structure around it so the final protein can be used as a vaccine. You must have prior understanding of Bioinformatics and Immunology in order to be able to understand what site to target and to supply it to the script. Once you identify a target site, the script will take it and run a long protocol, without the need for you to intervene, that will result in an ideal protein structure displaying your target site in its original 3D cofiguration. Thus the protien, theoretically, can be used as a vaccine against this site, and hopefully neutralise the disease you are researching. Everytime you run this script a different final protien structure will be generated, this is important to keep in mind, because if you want to generate many different structures to test or to use as boosts you can simply run the same target site again and you will end up with a different final structure.

This script has been last tested to work well with PyRosetta 4 Release 147 and using Python 3.5. If you use this script on a newer PyRosetta or Python version and it fails please notify me to get it updated.

Here is a [video](youtube.com/) that explains how to select a target site, how the script functions, and what results you sould get. If I did not make a video yet, bug me until I make one.

The script protocol is as follows:
1. Build Scaffold. --> STILL UNDER DEVELOPMENT --> I am having lots of trouble with De Novo Design (I have a very long temporary work around)
2. Isolate motif.
3. Isolate receptor.
4. Graft motif onto scaffold.
5. Sequence design the structure around the motif.
6. Generate fragments for Rosetta Abinitio folding simulation.
7. If average fragment RMSD is higher than 2Å repeat steps 5 and 6.

Output files are as follows:

|    | File Name               | Description                                                                                  |
|----|-------------------------|----------------------------------------------------------------------------------------------|
| 1  | DeNovo.pdb              | Scaffold structure                                                                           |
| 2  | motif.pdb	       | Original requested motif                                                                     |
| 3  | receptor.pdb            | Original receptor that binds motif                                                           |
| 4  | grafted.pdb             | Grafted motif to De Novo structure                                                           |
| 5  | structure.pdb           | Sequence designed structure                                                                  |
| 6  | structure.fasta         | Fasta of Rosetta Designed structure                                                          |
| 7  | frags.200.3mers         | 3-mer fragment of sequence designed structure from the Robetta server                        |
| 8  | frags.200.9mers         | 9-mer fragment of sequence designed structure from the Robetta server                        |
| 9  | pre.psipred.ss2         | PSIPRED secondary structure prediction of sequence designed structure from the Robetta server|
| 10 | plot_frag.pdf           | Plot of the fragment quality RMSD vs Position                                                |
| 11 | RMSDvsPosition.dat      | plot_frag.pdf's data                                                                         |
| 12 | FragmentAverageRMSD.dat | Average RMSD of the fragments                                                                |
'''
#--------------------------------------------------------------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------------------------------------------------------------
#The Functions

#1 - Extract Motif
def Motif(PDB_ID , Chain , Motif_From , Motif_To):
	'''
	This function downloads a spesific protein from RCSB
	and isolates a specific user defined motif from it
	Generates the motif.pdb file
	'''
	#Get the protein
	os.system('wget http://www.rcsb.org/pdb/files/' + PDB_ID + '.pdb')
	pdb = open(PDB_ID + '.pdb' , 'r')
	#Isolate the motif
	Motif = open('motif.pdb' , 'w')
	count = 0
	num = 0
	AA2 = None
	for line in pdb:
		if not line.startswith('ATOM'):					#Ignore all lines that do not start with ATOM
			continue
		if not line.split()[4] == Chain:				#Ignore all lines that do not have the specified chain (column 5)
			continue
		if int(Motif_From) <= int(line.split()[5]) <= int(Motif_To):	#Find residues between the user specified location
			count += 1						#Sequencially number atoms
			AA1 = line[23:27]					#Sequencially number residues
			if not AA1 == AA2:
				num += 1			
			final_line = line[:7] + '{:4d}'.format(count) + line[11:17] + line[17:21] + 'A' + '{:4d}'.format(num) + line[26:]	#Update each line of the motif to have its atoms and residues sequencially labeled, as well as being in chain A
			AA2 = AA1
			Motif.write(final_line)					#Write to new file called motif.pdb
	Motif.close()

#2 - Extract Receptor
def Receptor(PDB_ID , Chain):
	'''
	This function isolates a chain from a downloaded .pdb file 
	Generates the receptor.pdb file
	'''
	#Isolate the receptor
	pdb = open(PDB_ID + '.pdb' , 'r')
	Receptor = open('receptor.pdb' , 'w')
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
def Graft(receptor , motif , scaffold):
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
	mover.init_parameters(receptor , motif , 1.0 , 2 , 5 , '0:0' , '0:0' , 'ALA' , hotspots , True , False , True , False , False , True , False)
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
	pdb = open('temp.pdb' , 'r')
	Structure = open('temp2.pdb' , 'w')
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
	newpdb = open('temp2.pdb' , 'r')
	thenewfile = open('grafted.pdb' , 'w')
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
	return((start , finish))													#Return values [0] = Motif_From [1] = Motif_To

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
	def BLAST(self , filename1 , filename2):
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
		print('Sequence Similarity: {}%'.format(percentage))

	#5.2 - Preforms RosettaDesign
	def fixbb(self , filename , relax_iters , design_iters):
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
		pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn , packtask)
		backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
		backrub.pivot_residues(pose)
		GMC = pyrosetta.rosetta.protocols.monte_carlo.GenericMonteCarloMover()
		GMC.set_mover(backrub)
		GMC.set_scorefxn(scorefxn)
		GMC.set_maxtrials(500)
		GMC.set_temperature(1.0)
		GMC.set_preapply(False)
		GMC.set_recover_low(True)
		mover = pyrosetta.rosetta.protocols.moves.SequenceMover()
		mover.add_mover(pack)
		mover.add_mover(GMC)#####<--- problem here not accepting not rejecting moves
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
		print('Relax Scores:\n' , Rscores)
		print('Chosen Lowest Score:' , RFinalScore , '\n')
		print('Design Scores:\n' , Dscores)
		print('Chosen Lowest Score:' , DFinalScore , '\n')
		print('BLAST result, compairing the original structure to the designed structure:')
		RosettaDesign.BLAST(self , filename , 'structure.pdb')

	#5.3 - Preforms RosettaDesign for the whole protein except the motif
	def motif_fixbb(self , filename , Motif_From , Motif_To):
		'''
		Applies RosettaDesign to change the structure's
		amino acids (one layer at a time - like in the
		Design_Layer method) except for a desired
		continuous motif sequence while maintaining the
		same backbone.
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
		#B - Perform fixbb RosettaDesign without the motif residues
		packtask = standard_packer_task(pose)
		#Identify motif residues
		Motif = list(range(int(Motif_From) , int(Motif_To) + 1))
		#Prevent motif residues from being designed
		for aa in Motif:
			pack.temporarily_set_pack_residue(aa , False)
		pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn , packtask)
		backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
		backrub.pivot_residues(pose)
		GMC = pyrosetta.rosetta.protocols.monte_carlo.GenericMonteCarloMover()
		GMC.set_mover(backrub)
		GMC.set_scorefxn(scorefxn)
		GMC.set_maxtrials(500)
		GMC.set_temperature(1.0)
		GMC.set_preapply(False)
		GMC.set_recover_low(True)
		mover = pyrosetta.rosetta.protocols.moves.SequenceMover()
		mover.add_mover(pack)
		mover.add_mover(GMC)#####<--- problem here not accepting not rejecting moves
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
		print('Relax Scores:\n' , Rscores)
		print('Chosen Lowest Score:' , RFinalScore , '\n')
		print('Design Scores:\n' , Dscores)
		print('Chosen Lowest Score:' , DFinalScore , '\n')
		print('BLAST result, compairing the original structure to the designed structure:')
		RosettaDesign.BLAST(self , filename , 'structure.pdb')

#6 - Fragment Generation and Identification
def Fragments(filename):
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
		'UserName':'ac.research',
		'Email':'',
		'Notes':'structure',
		'Sequence':sequence,
		'Fasta':'',
		'Code':'',
		'ChemicalShifts':'',
		'NoeConstraints':'',
		'DipolarConstraints':'',
		'type':'submit'
	}
	session = requests.session()
	response = session.post('http://www.robetta.org/fragmentsubmit.jsp', data=payload , files=dict(foo='bar'))		
	for line in response:
		line = line.decode()
		if re.search('<a href="(fragmentqueue.jsp\?id=[0-9].*)">' , line):
			JobID = re.findall('<a href="(fragmentqueue.jsp\?id=[0-9].*)">' , line)
	JobURL = 'http://www.robetta.org/' + JobID[0]
	#Check
	ID = JobID[0].split('=')
	print('Job ID: ' + str(ID[1]))
	while True:
		Job = urllib.request.urlopen(JobURL)
		jobdata = bs4.BeautifulSoup(Job , 'lxml')
		status = jobdata.find('td', string='Status: ').find_next().text
		if status == 'Complete':
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M') , 'Status:' , status)
			break
		else:
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M') , 'Status:' , status)
			time.sleep(1800)
			continue
	#Download
	sequence = pose.sequence()
	fasta = open('structure.fasta' , 'w')
	fasta.write(sequence)
	fasta.close()
	time.sleep(1)
	os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/aat000_03_05.200_v1_3')
	os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/aat000_09_05.200_v1_3')
	os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/t000_.psipred_ss2')
	os.rename('aat000_03_05.200_v1_3' , 'frags.200.3mers')
	os.rename('aat000_09_05.200_v1_3' , 'frags.200.9mers')
	os.rename('t000_.psipred_ss2' , 'pre.psipred.ss2')
	#Calculate the best fragment's RMSD at each position
	frag = open('frags.200.9mers' , 'r')
	rmsd = open('temp.dat' , 'w')
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
		fragset.frames(count , frames)
		#Setup the MoveMap
		movemap = MoveMap()
		movemap.set_bb(True)
		#Setup and apply the fragment inserting mover
		for frame in frames:
			for frag_num in range( 1 , frame.nr_frags() + 1 ):
				frame.apply(movemap , frag_num , pose_copy)
				#Measure the RMSD difference between the original pose and the new changed pose (the copy)
				RMSD = rosetta.core.scoring.CA_rmsd(pose , pose_copy)
				print(RMSD , '\t' , count)
				rmsd.write(str(RMSD) + '\t' + str(count) + '\n')
				#Reset the copy pose to original pose
				pose_copy.assign(pose)
	rmsd.close()
	#Analyse the RMSD file to get the lowest RMSD for each position
	data = open('RMSDvsPosition.dat' , 'w')
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
	data = open('RMSDvsPosition.dat' , 'r')
	value = 0
	for line in data:
		line = line.split()
		RMSD = float(line[1])
		value = value + RMSD
		count = int(line[0])
	Average_RMSD = round(value / count , 2)
	#Plot the results
	gnuplot = open('gnuplot_sets' , 'w')
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
	set label 'Average RMSD = {}' at graph 0.01 , graph 0.95 tc lt 7 font 'curior 12'\n
	plot 'RMSDvsPosition.dat' with boxes\n
	exit
	""".format(str(Average_RMSD)))
	gnuplot.close()
	os.system('gnuplot < gnuplot_sets')
	os.remove('gnuplot_sets')
	os.remove('temp.dat')
	return(Average_RMSD)






#7 - DeNovo Design ###### MAYBE FORGET ABOUT GAN AND IMPLEMENT EPIGRAFTING   only of the required small dataset
def GAN():
	pass
#--------------------------------------------------------------------------------------------------------------------------------------
#List of All Functions And Their Arguments
'''
Remember: DeNovo() , Graft() , Design.Motif() do not export the pose, therefore you must call the pose from the exported .pdb file after each function so the pose can be used by the subsequent function.
Motif(Protein , Chain , Motif_from , Motif_to)
Receptor(Protein , RecChain)
Relax(pose)
SASA(pose)
MotifPosition = Graft('receptor.pdb' , 'motif.pdb' , pose)
Design.Pack(pose)
Design.Motif(pose , Motif_from , Motif_to)
Fragments(pose)
GAN()
'''
#--------------------------------------------------------------------------------------------------------------------------------------
def protocol():
	#User inputs
	Protein		= sys.argv[1]
	RecChain	= sys.argv[2]
	Chain		= sys.argv[3]
	Motif_from	= sys.argv[4]
	Motif_to	= sys.argv[5]

	#1. Build scaffold
	#GAN()										##### <-------- Requires Work
	pose = pose_from_pdb('generated.pdb')

	#2. Isolate motif
	Motif(Protein , Chain , Motif_from , Motif_to)

	#3. Isolate receptor
	Receptor(Protein , RecChain)

	#4. Graft motif onto scaffold
	MotifPosition = Graft('receptor.pdb' , 'motif.pdb' , pose)			##### <-------- Requires Work

	#5. Sequence design the structure around the motif
	RosettaDesign.motif_fixbb('grafted.pdb' , MotifPosition[0] , MotifPosition[1])	##### <-------- Requires Work

	#6. Generate to test fragment quality and predict the Abinitio folding simulation success
	Fragments('structure.pdb')

if __name__ == '__main__':
	protocol()
