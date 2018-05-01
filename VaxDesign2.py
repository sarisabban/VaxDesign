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

#Import Modules							   #Modules to download--->   #      #
import sys , os , re , time , datetime , subprocess , random , requests , urllib.request , Bio.PDB, bs4
from pyrosetta import *
from pyrosetta.toolbox import *
init()

#User Inputs
Protein		= sys.argv[1]
RecChain	= sys.argv[2]
Chain		= sys.argv[3]
Motif_from	= sys.argv[4]
Motif_to	= sys.argv[5]
#--------------------------------------------------------------------------------------------------------------------------------------
#The Functions

#1 - Extract Motif
def Motif(PDB_ID , Chain , Motif_From , Motif_To):
	''' This function downloads a spesific protein from RCSB and isolates a specific user defined motif from it '''
	''' Generates the motif.pdb file '''
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
	''' This function isolates a chain from a downloaded .pdb file '''
	''' Generates the receptor.pdb file '''
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
	''' Relaxes a structure '''
	''' Updates the original pose with the relaxed pose '''
	scorefxn = get_fa_scorefxn()
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.apply(pose)

#4 - SASA
def SASA(pose):
	''' Calculates the different layers (Surface, Boundary, Core) of a structure according its SASA (solvent-accessible surface area) '''
	''' Returns three lists Surface amino acids = [0] , Boundary amino acids = [1] , Core amino acids = [2] '''
	#Temporary generate a .pdb file of the pose to isolate the layers since it is not yet possible to do that using a Rosetta pose, this temporary .pdb file will be deleted after the layers are found
	pose.dump_pdb('ToDesign.pdb')
	#Standard script to setup biopython's DSSP to calculate SASA using Wilke constants
	parser = Bio.PDB.PDBParser()
	structure = parser.get_structure('X' , 'ToDesign.pdb')
	model = structure[0]
	dssp = Bio.PDB.DSSP(model , 'ToDesign.pdb' , acc_array = 'Wilke')
	#Loop to get SASA for each amino acid
	lis = list()
	count = 0
	for x in dssp:
		if x[1]   == 'A' : sasa = 129 * (x[3])
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
	#Label each amino acid depending on its SASA position according to the parameters highlighted in the paper by (Koga et.al., 2012 - PMID: 23135467). The parameters are as follows:
	#Surface:	Helix or Sheet: SASA => 60		Loop: SASA => 40
	#Boundry:	Helix or Sheet: 15 < SASA < 60		Loop: 25 < SASA < 40
	#Core:		Helix or Sheet: SASA =< 15		Loop: SASA =< 25	
	surface = list()
	boundary = list()
	core = list()
	count = 0
	for x , y in lis:
		count = count + 1
		if y <= 25 and (x == '-' or x == 'T' or x == 'S'):		#Loop (DSSP code is - or T or S)
			core.append(count)
		elif 25 < y < 40 and (x == '-' or x == 'T' or x == 'S'):	#Loop (DSSP code is - or T or S)
			boundary.append(count)
		elif y >= 40 and (x == '-' or x == 'T' or x == 'S'):		#Loop (DSSP code is - or T or S)
			surface.append(count)
		elif y <= 15 and (x == 'G' or x == 'H' or x == 'I'):		#Helix (DSSP code is G or H or I)
			core.append(count)
		elif 15 < y < 60 and (x == 'G' or x == 'H' or x == 'I'):	#Helix (DSSP code is G or H or I)
			boundary.append(count)
		elif y >= 60 and (x == 'G' or x == 'H' or x == 'I'):		#Helix (DSSP code is G or H or I)
			surface.append(count)
		elif y <= 15 and (x == 'B' or x == 'E'):			#Sheet (DSSP code is B or E)
			core.append(count)
		elif 15 < y < 60 and (x == 'B' or x == 'E'):			#Sheet (DSSP code is B or E)
			boundary.append(count)
		elif y >= 60 and (x == 'B' or x == 'E'):			#Sheet (DSSP code is B or E)
			surface.append(count)	
	os.remove('ToDesign.pdb')						#Keep working directory clean
	return(surface , boundary , core)

#5 - Grafting
def Graft(receptor , motif , scaffold):
	''' Grafts a motif onto a protein scaffold structure '''
	''' Generates structure.pdb and returns a tuple [0] is the residue number where the motif starts and [1] where it ends '''
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

#6 - RosettaDesign
class Design():
	#6.1 - Design the structure one layer at a time moving towards a tightly packed core with every loop
	def Pack(pose):
		''' Applies FastDesign to change the whole structure's amino acids (one layer at a time as well as designing towards an optimally packed core) while maintaining the same backbone. Should be faster than the Whole method and results in a better final structure than the Layer method '''
		''' Generates the Designed.pdb file '''
		#A - Relax original structure
		scorefxn = get_fa_scorefxn()												#Call the score function
		score1_original_before_relax = scorefxn(pose)										#Measure score before relaxing
		Relax(pose)														#Relax structure
		score2_original_after_relax = scorefxn(pose)										#Measure score after relaxing
		#B - FastDesign protocol												#Uses Generic Monte Carlo with PackStat as a filter to direct FastDesign towards an optimally packed structure core
		chain = pose.pdb_info().chain(1)											#Identify chain
		layers = [2 , 1 , 0]													#Layer Identity from SASA Surface = [0] , Boundary = [1] , Core = [2]
		for identity in layers:													#Loop through each layer
			#1 - Setup the PackStat filter
			filters = rosetta.protocols.simple_filters.PackStatFilter()
			#2 - Identify the layers
			sasa = SASA(pose)												#Re-calculate SASA every time because amino acid position can change from one layer to another during the design phase, therefore make sure to design the layer not the amino acid
			layer = sasa[identity]												#Changes every iteration to start with Core (sasa[2]) then Boundary (sasa[1]) then Surface (sasa[0])
			#3 - Generate the resfile											#Will generate a new Resfile for each layer (which is why it is deleted at the end of the loop)
			Resfile = open('Resfile.resfile' , 'w')
			Resfile.write('NATAA\n')
			Resfile.write('start\n')
			for line in layer:
				Resfile.write(str(line) + ' ' + chain + ' ALLAA\n')
			Resfile.close()
			#4 - Setup the FastDesign mover
			read = pyrosetta.rosetta.core.pack.task.operation.ReadResfile('Resfile.resfile')				#Call the generated Resfile
			task = pyrosetta.rosetta.core.pack.task.TaskFactory()								#Setup the TaskFactory
			task.push_back(read)												#Add the Resfile to the TaskFactory
			movemap = MoveMap()												#Setup the MoveMap
			movemap.set_bb(False)												#Do not change the phi and psi BackBone angles
			movemap.set_chi(True)												#Change the chi Side Chain angle
			mover = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()						#Call the FastDesign Mover
			mover.set_task_factory(task)											#Add the TaskFactory to it
			mover.set_movemap(movemap)											#Add the MoveMap to it
			mover.set_scorefxn(scorefxn)											#Add the Score Function to it
			#5 - Setup and apply the generic Monte Carlo mover
			MC = pyrosetta.rosetta.protocols.monte_carlo.GenericMonteCarloMover()						#Call Monter Carlo Class
			MC.set_mover(mover)												#Load The Mover
			MC.set_scorefxn(scorefxn)											#Set score function
			MC.set_maxtrials(10)												#Set number of monte carlo loops
			MC.set_temperature(1)												#Set temperature
			MC.set_preapply(False)												#To apply Boltzmann accept/reject to all applications of the mover (always use False)
			MC.set_drift(True)												#Make current pose = next iteration pose
			MC.set_sampletype('high')											#Move monte carlo to higher filter score
			#MC.recover_low(True)												#True - at the end of application, the pose is set to the lowest (or highest if sample_type="high") scoring pose
			#MC.stopping_condition()											#Stops before trials are done if a filter evaluates to true
			MC.add_filter(filters , False , 1.0 , 'high' , True)								#Add a filter (Filter Type , Adaptive , Temperature , Sample Type , Rank By)
			#MC.task_factory(task) #Causes an infinite loop									#Include a Task Factory
			#MC.boltzmann(pose) #For some reason hates a relaxed pose							#Evaulates a pose based on the scores/filters + temperatures
			MC.apply(pose)													#Apply Move
			os.remove('Resfile.resfile')											#To keep working directory clean, and to make sure each Resfile has the info for each layer only and they do not get mixed and appended together in one Resfile
		#C - Relax pose
		Relax(pose)														#Relax structure
		#D - Output result
		score3_of_design_after_relax = scorefxn(pose)										#Measure score of designed pose
		pose.dump_pdb('structure.pdb')												#Export final pose into a .pdb structure file
		print('---------------------------------------------------------')
		print('Original Structure Score:' , '\t' , score1_original_before_relax)
		print('Relaxed Original Score:' , '\t' , score2_original_after_relax)
		print('Relaxed Design Score:' , '\t\t' , score3_of_design_after_relax)
	#6.2 - Design the structure one layer at a time, except for a desired motif, moving towards a tightly packed core with every loop
	def Motif(pose , Motif_From , Motif_To):
		''' Applies RosettaDesign to change the structure's amino acids (one layer at a time - like in the Design_Layer method) except for a desired continuous motif sequence while maintaining the same backbone '''
		''' Just updates the pose with the new structure '''
		#A - Relax original structure
		Motif = list(range(int(Motif_From) , int(Motif_To) + 1))								#Identify motif residues
		scorefxn = get_fa_scorefxn()												#Call the score function
		score1_original_before_relax = scorefxn(pose)										#Measure score before relaxing
		Relax(pose)														#Relax structure
		score2_original_after_relax = scorefxn(pose)										#Measure score after relaxing
		#B - FastDesign protocol												#Uses Generic Monte Carlo with PackStat as a filter to direct FastDesign towards an optimally packed structure core
		chain = pose.pdb_info().chain(1)											#Identify chain
		layers = [2 , 1 , 0]													#Layer Identity from SASA Surface = [0] , Boundary = [1] , Core = [2]
		for identity in layers:													#Loop through each layer
			#1 - Setup the PackStat filter
			filters = rosetta.protocols.simple_filters.PackStatFilter()
			#2 - Identify the layers
			sasa = SASA(pose)												#Re-calculate SASA every time because amino acid position can change from one layer to another during the design phase, therefore make sure to design the layer not the amino acid
			layer = sasa[identity]												#Changes every iteration to start with Core (sasa[2]) then Boundary (sasa[1]) then Surface (sasa[0])
			#3 - Remove motif residues as to not get designed
			layer = list(set(layer) - set(Motif))
			#4 - Generate the resfile											#Will generate a new Resfile for each layer (which is why it is deleted at the end of the loop)
			Resfile = open('Resfile.resfile' , 'w')
			Resfile.write('NATAA\n')
			Resfile.write('start\n')
			for line in layer:
				Resfile.write(str(line) + ' ' + chain +' ALLAA\n')
			Resfile.close()
			#5 - Setup the FastDesign mover
			read = pyrosetta.rosetta.core.pack.task.operation.ReadResfile('Resfile.resfile')				#Call the generated Resfile
			task = pyrosetta.rosetta.core.pack.task.TaskFactory()								#Setup the TaskFactory
			task.push_back(read)												#Add the Resfile to the TaskFactory
			movemap = MoveMap()												#Setup the MoveMap
			movemap.set_bb(False)												#Do not change the phi and psi BackBone angles
			movemap.set_chi(True)												#Change the chi Side Chain angle
			mover = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()						#Call the FastDesign Mover
			mover.set_task_factory(task)											#Add the TaskFactory to it
			mover.set_movemap(movemap)											#Add the MoveMap to it
			mover.set_scorefxn(scorefxn)											#Add the Score Function to it
			#6 - Setup and apply the Generic Monte Carlo mover
			MC = pyrosetta.rosetta.protocols.monte_carlo.GenericMonteCarloMover()						#Call Monter Carlo Class
			MC.set_mover(mover)												#Load The Mover
			MC.set_scorefxn(scorefxn)											#Set score function
			MC.set_maxtrials(10)												#Set number of monte carlo loops
			MC.set_temperature(1)												#Set temperature
			MC.set_preapply(True)												#To apply Boltzmann accept/reject to all applications of the mover (always use False)
			MC.set_drift(True)												#Make current pose = next iteration pose
			MC.set_sampletype('high')											#Move monte carlo to higher filter score
			#MC.recover_low(True)												#True - at the end of application, the pose is set to the lowest (or highest if sample_type="high") scoring pose
			#MC.stopping_condition()											#Stops before trials are done if a filter evaluates to true
			MC.add_filter(filters , False , 1.0 , 'high' , True)								#Add a filter (Filter Type , Adaptive , Temperature , Sample Type , Rank By)
			#MC.task_factory(task) #Causes an infinite loop									#Include a Task Factory
			#MC.boltzmann(pose) #For some reason hates a relaxed pose							#Evaulates a pose based on the scores/filters + temperatures
			MC.apply(pose)													#Apply Move
			os.remove('Resfile.resfile')											#To keep working directory clean, and to make sure each Resfile has the info for each layer only and they do not get mixed and appended together in one Resfile
		#C - Relax pose
		Relax(pose)														#Relax structure
		#D - Output result
		score3_of_design_after_relax = scorefxn(pose)										#Measure score of designed pose
		pose.dump_pdb('structure.pdb')												#Export final pose into a .pdb structure file
		print('---------------------------------------------------------')
		print('Original Structure Score:' , '\t' , score1_original_before_relax)
		print('Relaxed Original Score:' , '\t' , score2_original_after_relax)
		print('Relaxed Design Score:' , '\t\t' , score3_of_design_after_relax)

#7 - Fragment Generation and Identification
def Fragments(pose):
	''' Submits the pose to the Robetta server (http://www.robetta.org) for fragment generation that are used for the Abinitio folding simulation. Then measures the RMSD for each fragment at each position and chooses the lowest RMSD. Then averages out the lowest RMSDs. Then plots the lowest RMSD fragment for each positon '''
	''' Generates the 3-mer file, the 9-mer file, the PsiPred file, the RMSD vs Position PDF plot with the averaged fragment RMSD printed in the plot '''
	#Make the 3-mer and 9-mer fragment files and the PSIPRED file using the Robetta server
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
	response = session.post('http://www.robetta.org/fragmentsubmit.jsp', data = payload , files = dict(foo = 'bar'))		
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
		status = jobdata.find('td', string = 'Status: ').find_next().text
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
	lowest = {} 									#Mapping group number -> lowest value found
	for line in open('temp.dat'):
		parts = line.split()
		if len(parts) != 2:							#Only lines with two items on it
			continue
		first = float(parts[0])
		second = int(parts[1])
		if first == 0: 								#Skip line with 0.0 RMSD (this is an error from the 9-mer fragment file). I don't know why it happens
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
	gnuplot.write("reset\nset terminal postscript\nset output './plot_frag.pdf'\nset encoding iso_8859_1\nset term post eps enh color\nset xlabel 'Position'\nset ylabel 'RMSD (\\305)'\nset yrange [0:]\nset xrange [0:]\nset xtics auto\nset xtics rotate\nset grid front\nunset grid\nset title 'Fragment Quality'\nset key off\nset boxwidth 0.5\nset style fill solid\nset label 'Average RMSD = " + str(Average_RMSD) + "' at graph 0.01 , graph 0.95 tc lt 7 font 'curior 12'\nplot 'RMSDvsPosition.dat' with boxes\nexit")
	gnuplot.close()
	os.system('gnuplot < gnuplot_sets')
	os.remove('gnuplot_sets')
	os.remove('temp.dat')
	return(Average_RMSD)

#8 - De Novo Design
def DeNovo():
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
DeNovo()
'''
#--------------------------------------------------------------------------------------------------------------------------------------
#The Protocol
#1. Build Scaffold
pose = pose_from_pdb('DeNovo.pdb')

#2. Isolate Motif
Motif(Protein , Chain , Motif_from , Motif_to)

#3. Isolate Receptor
Receptor(Protein , RecChain)

#4. Graft Motif onto Scaffold
MotifPosition = Graft('receptor.pdb' , 'motif.pdb' , pose)

#5. Sequence Design The Structure Around The Motif
home = os.getcwd()
for attempt in range(20):
	time.sleep(1)
	os.chdir(home)
	pose = pose_from_pdb('grafted.pdb')
	os.mkdir('Attempt_' + str(attempt + 1))
	os.chdir(home + '/Attempt_' + str(attempt + 1))
	Design.Motif(pose , MotifPosition[0] , MotifPosition[1])

	#6. Generate Fragments Locally To Test Fragment Quality And Predict Abinitio Fold Simulation Success
	pose = pose_from_pdb('structure.pdb')
	RMSD = Fragments(pose)

	#7. Average Fragment RMSD Should Be < 2Å - If Not Then Repeat
	if RMSD <= 2:
		break
	else:
		continue
