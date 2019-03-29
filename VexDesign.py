#!/usr/bin/python3

'''
MIT License

Copyright (c) 2017 Sari Sabban

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import os
import re
import bs4
import sys
import time
import glob
import numpy
import random
import Bio.PDB
import argparse
import requests
import datetime
import subprocess
import urllib.request
from pyrosetta import *
from pyrosetta.toolbox import *
print('\x1b[32m\n  ██╗   ██╗ █████╗ ██╗  ██╗\n  ██║   ██║██╔══██╗╚██╗██╔╝\n  ██║   ██║███████║ ╚███╔╝ \n  ╚██╗ ██╔╝██╔══██║ ██╔██╗ \n   ╚████╔╝ ██║  ██║██╔╝ ██╗\n    ╚═══╝  ╚═╝  ╚═╝╚═╝  ╚═╝\n                           \n  ██████╗ ███████╗███████╗██╗ ██████╗ ███╗   ██╗\n  ██╔══██╗██╔════╝██╔════╝██║██╔════╝ ████╗  ██║\n  ██║  ██║█████╗  ███████╗██║██║  ███╗██╔██╗ ██║\n  ██║  ██║██╔══╝  ╚════██║██║██║   ██║██║╚██╗██║\n  ██████╔╝███████╗███████║██║╚██████╔╝██║ ╚████║\n  ╚═════╝ ╚══════╝╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝\n                                                \x1b[0m\n\x1b[35m  ╔═╗┬ ┬┌┬┐┌─┐  ╔╦╗┌─┐┌─┐┬┌─┐┌┐┌  ╔═╗  ╦  ╦┌─┐┌─┐┌─┐┬┌┐┌┌─┐\n  ╠═╣│ │ │ │ │   ║║├┤ └─┐││ ┬│││  ╠═╣  ╚╗╔╝├─┤│  │  ││││├┤ \n  ╩ ╩└─┘ ┴ └─┘  ═╩╝└─┘└─┘┴└─┘┘└┘  ╩ ╩   ╚╝ ┴ ┴└─┘└─┘┴┘└┘└─┘\x1b[0m\n\x1b[33mAuthored by Sari Sabban on 31-May-2017 (sari.sabban@gmail.com)\x1b[0m\n\x1b[36m--------------------------------------------------------------\x1b[0m')
init()

parser = argparse.ArgumentParser(description='A script that autonomously designs a vaccine\nAuthored by Sari Sabban on 31-May-2017 (sari.sabban@gmail.com)\nhttps://sarisabban.github.io/VexDesign/')
parser.add_argument('-s', '--scaffold',		nargs='+', metavar='', help='search for scaffolds')
parser.add_argument('-p', '--protocol',		nargs='+', metavar='', help='Run full protocol')
parser.add_argument('-m', '--motif',		nargs='+', metavar='', help='Isolate motif')
parser.add_argument('-r', '--receptor',		nargs='+', metavar='', help='Isolate receptor')
parser.add_argument('-g', '--graft',		nargs='+', metavar='', help='Graft motif onto scaffold')
parser.add_argument('-f', '--ffl',			nargs='+', metavar='', help='Fold From Loop')
parser.add_argument('-d', '--design',		nargs='+', metavar='', help='Sequence design the structure around the motif')
parser.add_argument('-F', '--fragments',	nargs='+', metavar='', help='Generate and analyse fragments')
args = parser.parse_args()

def Motif(PDB_ID, Chain, Motif_From, Motif_To):
	'''
	This function downloads a spesific protein from RCSB and isolates a
	specific user defined motif from it Generates the motif.pdb file
	'''
	#Get the protein
	os.system('wget http://www.rcsb.org/pdb/files/{}.pdb'.format(PDB_ID))
	pdb = open('{}.pdb'.format(PDB_ID), 'r')
	#Isolate the motif
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
				final_line =	line[:7]+\
								'{:4d}'.format(count)+\
								line[11:17]+\
								line[17:21]+\
								'A'+\
								'{:4d}'.format(num)+\
								line[26:]
				AA2 = AA1
				Motif.write(final_line)
		except:
			continue
	Motif.close()

def Receptor(PDB_ID, Chain):
	'''
	This function isolates a chain from a downloaded .pdb file
	Generates the receptor.pdb file
	'''
	#Isolate the receptor
	pdb = open('{}.pdb'.format(PDB_ID), 'r')
	Receptor = open('receptor.pdb', 'w')
	for line in pdb:
		linesplit = line.split()
		if linesplit[0] == 'ATOM':
			if linesplit[4] == Chain:
				Receptor.write(line)
	Receptor.close()
	#Keep working directory clean, remove the protein's original file
	os.remove('{}.pdb'.format(PDB_ID))

def Graft(receptor, motif, scaffold):
	'''
	Grafts a motif onto a protein scaffold structure
	Generates grafted.pdb and returns a tuple [0] is the
	residue number where the motif starts and [1] where it ends
	'''
	scorefxn = get_fa_scorefxn()
	mover = pyrosetta.rosetta.protocols.motif_grafting.movers.MotifGraftMover()
	#Setup motif hotspots
	motifpose = pose_from_pdb(motif)
	spots = list()
	for resi in range(motifpose.total_residue()):
		spots.append(str(resi+1))
	hotspots = ':'.join(spots)
	#Setup grafting mover
	mover.init_parameters(
						receptor, 	# context_structure
						motif, 		# motif_structure
						1.0, 		# RMSD_tolerance
						2.0, 		# NC_points_RMSD_tolerance
						5, 			# clash_score_cutoff
						1, 			# min_fragment_size
						'0:0', 		# combinatory_fragment_size_delta
						'0:0', 		# max_fragment_replacement_size_delta
						'ALA', 		# clash_test_residue
						hotspots, 	# hotspots
						True, 		# full_motif_bb_alignment
						False, 		# allow_independent_alignment_per_fragment
						True, 		# graft_only_hotspots_by_replacement
						False, 		# only_allow_if_N_point_match_aa_identity
						False, 		# only_allow_if_C_point_match_aa_identity
						True, 		# revert_graft_to_native_sequence
						False,		# allow_repeat_same_graft_output
						1.0,		# output_cluster_tolerance
						GEO)		# output_filter
	mover.apply(scaffold)
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
		count += 1
		AA1 = line[23:27]
		if not AA1 == AA2:
			num += 1
		final_line =	line[:7]+\
						'{:4d}'.format(count)+\
						line[11:17]+\
						line[17:21]+\
						'A'+\
						'{:4d}'.format(num)+\
						line[26:]
		AA2 = AA1
		thenewfile.write(final_line)
	thenewfile.close()
	os.remove('temp2.pdb')
	#Identify start and finish residue number of inserted motif
	motifpose = pose_from_pdb('motif.pdb')
	graftpose = pose_from_pdb('grafted.pdb')
	MOTIF = motifpose.sequence()
	GRAFT = graftpose.sequence()
	start = GRAFT.index(MOTIF) + 1
	finish = start + len(MOTIF) - 1
	return((start, finish))

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
	#Make the 3-mer, 9-mer fragment files and PSIPRED files from Robetta server
	pose = pose_from_pdb(filename)
	sequence = pose.sequence()
	#Post
	web = requests.get('http://www.robetta.org/fragmentsubmit.jsp')
	payload = {	'UserName':				username,
				'Email':				'',
				'Notes':				filename.split('.')[0],
				'Sequence':				sequence,
				'Fasta':				'',
				'Code':					'',
				'ChemicalShifts':		'',
				'NoeConstraints':		'',
				'DipolarConstraints':	'',
				'type':					'submit'}
	session = requests.session()
	response = session.post('http://www.robetta.org/fragmentsubmit.jsp',\
							data=payload, files=dict(foo='bar'))
	for line in response:
		line = line.decode()
		if re.search('<a href="(fragmentqueue.jsp\?id=[0-9].*)">', line):
			JobID=re.findall('<a href="(fragmentqueue.jsp\?id=[0-9].*)">', line)
	JobURL = 'http://www.robetta.org/{}'.format(JobID[0])
	#Check
	ID = JobID[0].split('=')
	print('Job ID: {}'.format(str(ID[1])))
	while True:
		Job = urllib.request.urlopen(JobURL)
		jobdata = bs4.BeautifulSoup(Job, 'lxml')
		status = jobdata.find('td', string='Status: ').find_next().text
		if status == 'Complete':
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M'),\
															'Status:', status)
			break
		else:
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M'),\
															'Status:', status)
			time.sleep(900)
			continue
	#Download
	sequence = pose.sequence()
	fasta = open('structure.fasta', 'w')
	fasta.write(sequence)
	fasta.close()
	time.sleep(1)
	webserver = 'http://www.robetta.org/downloads/fragments'
	os.system('wget {}/{}/aat000_03_05.200_v1_3'.format(webserver, str(ID[1])))
	os.system('wget {}/{}/aat000_09_05.200_v1_3'.format(webserver, str(ID[1])))
	os.system('wget {}/{}/t000_.psipred_ss2'.format(webserver, str(ID[1])))
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
			for frag_num in range(1, frame.nr_frags()+1):
				frame.apply(movemap, frag_num, pose_copy)
				#Measure the RMSD between the original pose and the
				#new changed pose (the copy)
				RMSD = rosetta.core.scoring.CA_rmsd(pose, pose_copy)
				print(RMSD, '\t', count)
				rmsd.write('{}\t{}\n'.format(str(RMSD), str(count)))
				#Reset the copy pose to original pose
				pose_copy.assign(pose)
	rmsd.close()
	#Analyse the RMSD file to get the lowest RMSD for each position
	data = open('RMSDvsPosition.dat', 'w')
	lowest = {}
	for line in open('temp.dat'):
		parts = line.split()
		if len(parts) != 2:
			continue
		first = float(parts[0])
		second = int(parts[1])
		if first == 0:
			continue
		if second not in lowest:
			lowest[second] = first
		else:
			if first < lowest[second]:
				lowest[second] = first
	for position, rmsd in lowest.items():
		data.write(str(position)+'\t'+str(rmsd)+'\n')
	data.close()
	#Calculate the average RMSD of the fragments
	data = open('RMSDvsPosition.dat', 'r')
	value = 0
	for line in data:
		line = line.split()
		RMSD = float(line[1])
		value = value + RMSD
		count = int(line[0])
	Average_RMSD = round(value/count, 2)
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
	set label 'Average RMSD = {}' at graph 0.01,
	graph 0.95 tc lt 7 font 'curior 12'\n
	plot 'RMSDvsPosition.dat' with boxes\n
	exit
	""".format(str(Average_RMSD)))
	gnuplot.close()
	os.system('gnuplot < gnuplot_sets')
	os.remove('gnuplot_sets')
	os.remove('temp.dat')
	return(Average_RMSD)

def ScaffoldSearch(Protein, RecChain, Chain, Motif_From, Motif_To, Directory):
	'''
	This script searches a scaffold database for possible structures with
	successful grafting sites
	'''
	hotspots = Motif(Protein, Chain, Motif_From, Motif_To)
	Receptor(Protein, RecChain)
	os.mkdir('Scaffolds')
	current = os.getcwd()
	database = os.listdir(Directory)
	os.chdir(Directory)
	for scaffold in database:
		try:
			MotifGraft('../receptor.pdb', '../motif.pdb', scaffold, 1.0)
			os.system('cp {} ../Scaffolds'.format(scaffold))
		except:
			continue
	os.remove('../receptor.pdb')
	os.remove('../motif.pdb')

class RosettaDesign(object):
	def __init__(self, filename):
		''' Generate the resfile '''
		self.filename = filename
		parser = Bio.PDB.PDBParser()
		structure = parser.get_structure('{}'.format(filename), filename)
		dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke')
		sasalist = []
		for x in dssp:
			if x[1] == 'A':
				sasa = 129*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'V':
				sasa = 174*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'I':
				sasa = 197*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'L':
				sasa = 201*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'M':
				sasa = 224*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'P':
				sasa = 159*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'Y':
				sasa = 263*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'F':
				sasa = 240*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'W':
				sasa = 285*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'R':
				sasa = 274*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'N':
				sasa = 195*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'C':
				sasa = 167*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'Q':
				sasa = 225*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'E':
				sasa = 223*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'G':
				sasa = 104*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'H':
				sasa = 224*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'K':
				sasa = 236*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'S':
				sasa = 155*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'T':
				sasa = 172*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'D':
				sasa = 193*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			if x[2] == 'G' or x[2] == 'H' or x[2] == 'I':	ss = 'H'
			elif x[2] == 'B' or x[2] == 'E':				ss = 'S'
			elif x[2] == 'S' or x[2] == 'T' or x[2] == '-':	ss = 'L'
			sasalist.append((x[0], x[1], ss, sasa))
		resfile = open('resfile', 'a')
		resfile.write('NATRO\nSTART\n')
		for n, r, a, s in sasalist:
			if s == 'S' and a == 'L':
				line = '{} A PIKAA PGNQSTDERKH\n'.format(n)
				resfile.write(line)
			elif s == 'S' and a == 'H':
				line = '{} A PIKAA QEKH\n'.format(n)
				resfile.write(line)
			elif s == 'S' and a == 'S':
				line = '{} A PIKAA QTY\n'.format(n)
				resfile.write(line)
			elif s == 'B' and a == 'L':
				line = '{} A PIKAA AVILFYWGNQSTPDEKR\n'.format(n)
				resfile.write(line)
			elif s == 'B' and a == 'H':
				line = '{} A PIKAA AVILWQEKFM\n'.format(n)
				resfile.write(line)
			elif s == 'B' and a == 'S':
				line = '{} A PIKAA AVILFYWQTM\n'.format(n)
				resfile.write(line)
			elif s == 'C' and a == 'L':
				line = '{} A PIKAA AVILPFWM\n'.format(n)
				resfile.write(line)
			elif s == 'C' and a == 'H':
				line = '{} A PIKAA AVILFWM\n'.format(n)
				resfile.write(line)
			elif s == 'C' and a == 'S':
				line = '{} A PIKAA AVILFWM\n'.format(n)
				resfile.write(line)
		resfile.close()
	def __del__(self):
		''' Remove the resfile '''
		os.remove('resfile')
		for f in glob.glob('f[il]xbb.fasc'): os.remove(f)
	def choose(self):
		''' Choose the lowest scoring structure '''
		try:	scorefile = open('fixbb.fasc', 'r')
		except:	scorefile = open('flxbb.fasc', 'r')
		score = 0
		name = None
		for line in scorefile:
			line = json.loads(line)
			score2 = line.get('total_score')
			if score2 < score:
				score = score2
				name = line.get('decoy')
		os.system('mv {} structure.pdb'.format(name))
		for f in glob.glob('f[il]xbb_*'): os.remove(f)
	def fixbb(self):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while maintaining a fixed backbone.
		Generates the structure.pdb file
		'''
		pose = pose_from_pdb(self.filename)
		starting_pose = Pose()
		starting_pose.assign(pose)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		packtask = standard_packer_task(pose)
		pyrosetta.rosetta.core.pack.task.parse_resfile(pose,packtask,'resfile')
		fixbb = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, packtask, 10)
		job = PyJobDistributor('fixbb', 100, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			relax.apply(pose)
			fixbb.apply(pose)
			relax.apply(pose)
			job.output_decoy(pose)
		self.choose()
	def flxbb(self):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while allowing for a flexible backbone.
		Generates the structure.pdb file
		'''
		pose = pose_from_pdb(self.filename)
		starting_pose = Pose()
		starting_pose.assign(pose)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		resfile = rosetta.core.pack.task.operation.ReadResfile('resfile')
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()
		task.push_back(resfile)
		movemap = MoveMap()
		movemap.set_bb(True)
		movemap.set_chi(True)
		flxbb = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
		flxbb.set_task_factory(task)
		flxbb.set_movemap(movemap)
		flxbb.set_scorefxn(scorefxn)
		job = PyJobDistributor('flxbb', 100, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			relax.apply(pose)
			flxbb.apply(pose)
			relax.apply(pose)
			job.output_decoy(pose)
		self.choose()
	def fixbb_motif(self, Motif_From, Motif_To):
		'''
		Applies RosettaDesign with a fixed back bone to change the structure's
		amino acids (one layer at a time - like in the Design_Layer method)
		except for a desired continuous motif sequence while maintaining the
		same backbone. Just updates the pose with the new structure
		'''
		pose = pose_from_pdb(self.filename)
		starting_pose = Pose()
		starting_pose.assign(pose)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		packtask = standard_packer_task(pose)
		pyrosetta.rosetta.core.pack.task.parse_resfile(pose,packtask,'resfile')
		#Identify motif residues
		Motif = list(range(int(Motif_From), int(Motif_To) + 1))
		for aa in Motif:
			packtask.temporarily_set_pack_residue(aa, False)
		pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, packtask)
		print(packtask)
		fixbb = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, packtask, 10)
		job = PyJobDistributor('fixbb', 100, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			relax.apply(pose)
			fixbb.apply(pose)
			relax.apply(pose)
			job.output_decoy(pose)
		self.choose()
	def flxbb_motif(self, Motif_From, Motif_To):
		'''
		Applies RosettaDesign with a flexible back bone to
		change the structure's amino acids (one layer at a
		time - like in the Design_Layer method) except for
		a desired continuous motif sequence while maintaining
		the same backbone.
		'''
		# Remove motif residues from resfile
		resfile = open('resfile', 'r')
		resfile2 = open('resfile2', 'a')
		resfile2.write('NATRO\nSTART\n')
		next(resfile)
		next(resfile)
		for res in resfile:
			if not int(Motif_From) <= int(res.split()[0]) <= int(Motif_To):
				resfile2.write(res)
		resfile2.close()
		pose = pose_from_pdb(self.filename)
		starting_pose = Pose()
		starting_pose.assign(pose)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		resfile = rosetta.core.pack.task.operation.ReadResfile('resfile2')
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()
		task.push_back(resfile)
		movemap = MoveMap()
		movemap.set_bb(True)
		movemap.set_chi(True)
		flxbb = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
		flxbb.set_task_factory(task)
		flxbb.set_movemap(movemap)
		flxbb.set_scorefxn(scorefxn)
		job = PyJobDistributor('flxbb', 100, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			relax.apply(pose)
			flxbb.apply(pose)
			relax.apply(pose)
			job.output_decoy(pose)
		self.choose()
		os.remove('resfile2')
def FFL(Motif, Scaffold, Motif_From, Motif_To, username):
	'''
	Performs the Fold From Loops protocol
	'''
	print('FFL is not yet available in PyRosetta')

def protocol(Protein,RChain,Chain,Motif_from,Motif_to,Scaffold,Choice,UserName):
	#1. Import scaffold
	pose = pose_from_pdb(Scaffold)
	#2. Isolate motif
	Motif(Protein, Chain, Motif_from, Motif_to)
	#3. Isolate receptor
	Receptor(Protein, RChain)
	#4. Graft motif onto scaffold
	MotifPosition = Graft('receptor.pdb', 'motif.pdb', pose)
	#5. Fold From Loop
	#FFL('motif.pdb', 'grafted.pdb', MotifPosition, UserName)
	#6. RosettaDesign the structure around the motif
	if Choice == 'fixbb':
		RD = RosettaDesign('grafted.pdb')
		RD.fixbb_motif(MotifPosition[0], MotifPosition[1])
	elif Choice == 'flxbb':
		RD = RosettaDesign('grafted.pdb')
		RD.flxbb_motif(MotifPosition[0], MotifPosition[1])
	#7. Generate and analyse fragments
	Fragments('structure.pdb', UserName)

def main():
	if args.scaffold:		# Search for scaffolds
		ScaffoldSearch(	sys.argv[2],		# PDB ID
						sys.argv[3],		# Receptor chain
						sys.argv[4],		# Motif chain
						sys.argv[5],		# Motif from
						sys.argv[6],		# Motif to
						sys.argv[7])		# Directory of scaffolds
	elif args.protocol:		# Run full protocol
		protocol(		sys.argv[2],		# PDB ID
						sys.argv[3],		# Receptor chain
						sys.argv[4],		# Motif chain
						sys.argv[5],		# Motif from
						sys.argv[6],		# Motif to
						sys.argv[7],		# Scaffold PDB file name
						sys.argv[8],		# RosettaDesign choice
						sys.argv[9])		# Robetta server username
	elif args.motif:		# Isolate motif
		Motif(			sys.argv[2],		# PDB ID
						sys.argv[3],		# Motif chain
						sys.argv[4],		# Motif from
						sys.argv[5])		# Motif to
		os.remove('{}.pdb'.format(sys.argv[2]))
	elif args.receptor:		# Isolate receptor
		RCSB = 'http://www.rcsb.org/pdb/files'
		os.system('wget {}/{}.pdb'.format(RCSB, sys.argv[2]))
		Receptor(		sys.argv[2],		# PDB ID
						sys.argv[3])		# Receptor chain
	elif args.graft:		# Graft motif onto scaffold
		pose = pose_from_pdb(sys.argv[4])	# Scaffold PDB file name
		MotifPosition = Graft(sys.argv[2],	# Receptor PDB file name
							sys.argv[3],	# Motif PDB file name
							pose)
		print('Grafted between positions: {} and {}'.format(MotifPosition[0],
															MotifPosition[1]))
	elif args.ffl:			# Fold From Loop
		FFL(			sys.argv[2],		# Motif PDB file name
						sys.argv[3],		# Scaffold PDB file name
						sys.argv[4],		# Motif on scaffold from
						sys.argv[5],		# Motif on scaffold to
						sys.argv[6])		# Robetta server username
	elif args.design:		# Sequence design the structure around the motif
		if sys.argv[2] == 'fixbb':			# Choice
			RD = RosettaDesign(	sys.argv[3])# Scaffold PDB file name
			RD.fixbb_motif(		sys.argv[4],# Motif on scaffold from
								sys.argv[5])# Motif on scaffold to
		elif sys.argv[2] == 'flxbb':		# Choice
			RD = RosettaDesign(	sys.argv[3])# Scaffold PDB file name
			RD.flxbb_motif(		sys.argv[4],# Motif on scaffold from
								sys.argv[5])# Motif on scaffold to
	elif args.fragments:	# Generate fragments
		Fragments(		sys.argv[2],		# Filename
						sys.argv[3])		# Username

if __name__ == '__main__': main()

#2y7q A B 420 429
#Graft()
#ScaffoldSearch()
#protocol()
