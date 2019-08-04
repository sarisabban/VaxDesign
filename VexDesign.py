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
print('\x1b[32m\n  ██╗   ██╗ █████╗ ██╗  ██╗\n  ██║   ██║██╔══██╗╚██╗██╔╝\n  ██║   ██║███████║ ╚███╔╝ \n  ╚██╗ ██╔╝██╔══██║ ██╔██╗ \n   ╚████╔╝ ██║  ██║██╔╝ ██╗\n    ╚═══╝  ╚═╝  ╚═╝╚═╝  ╚═╝\n                           \n  ██████╗ ███████╗███████╗██╗ ██████╗ ███╗   ██╗\n  ██╔══██╗██╔════╝██╔════╝██║██╔════╝ ████╗  ██║\n  ██║  ██║█████╗  ███████╗██║██║  ███╗██╔██╗ ██║\n  ██║  ██║██╔══╝  ╚════██║██║██║   ██║██║╚██╗██║\n  ██████╔╝███████╗███████║██║╚██████╔╝██║ ╚████║\n  ╚═════╝ ╚══════╝╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝\n                                                \x1b[0m\n\x1b[35m  ╔═╗┬ ┬┌┬┐┌─┐  ╔╦╗┌─┐┌─┐┬┌─┐┌┐┌  ╔═╗  ╦  ╦┌─┐┌─┐┌─┐┬┌┐┌┌─┐\n  ╠═╣│ │ │ │ │   ║║├┤ └─┐││ ┬│││  ╠═╣  ╚╗╔╝├─┤│  │  ││││├┤ \n  ╩ ╩└─┘ ┴ └─┘  ═╩╝└─┘└─┘┴└─┘┘└┘  ╩ ╩   ╚╝ ┴ ┴└─┘└─┘┴┘└┘└─┘\x1b[0m\n\x1b[33mAuthored by Sari Sabban on 31-May-2017 (sari.sabban@gmail.com)\nRepository is found at https://sarisabban.github.io/VexDesign/\x1b[0m\n\x1b[36m--------------------------------------------------------------\x1b[0m')
init(' -out:level 0 -no_his_his_pairE -extrachi_cutoff 1 -multi_cool_annealer 10 -ex1 -ex2 -use_input_sc')
print('\x1b[36m--------------------------------------------------------------\x1b[0m')

parser = argparse.ArgumentParser(description='A script that autonomously designs a vaccine\nAuthored by Sari Sabban on 31-May-2017 (sari.sabban@gmail.com)\nhttps://sarisabban.github.io/VexDesign/')
parser.add_argument('-s', '--scaffold', nargs='+', metavar='', help='search for scaffolds')
parser.add_argument('-p', '--protocol', nargs='+', metavar='', help='Run full protocol')
parser.add_argument('-m', '--motif',    nargs='+', metavar='', help='Isolate motif')
parser.add_argument('-r', '--receptor', nargs='+', metavar='', help='Isolate receptor')
parser.add_argument('-g', '--graft',    nargs='+', metavar='', help='Graft motif onto scaffold')
parser.add_argument('-f', '--ffl',      nargs='+', metavar='', help='Fold From Loop')
parser.add_argument('-d', '--design',   nargs='+', metavar='', help='Sequence design the structure around the motif')
parser.add_argument('-F', '--fragments',nargs='+', metavar='', help='Generate and analyse fragments')
args = parser.parse_args()

def Motif(PDB_ID, Chain, Motif_From, Motif_To):
	'''
	This function downloads a spesific protein from RCSB and isolates a
	specific user defined motif from it Generates the motif.pdb file
	'''
	#Get the protein
	os.system('wget -q http://www.rcsb.org/pdb/files/{}.pdb'.format(PDB_ID))
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
	#Setup motif hotspots
	motifpose = pose_from_pdb(motif)
	spots = list()
	for resi in range(motifpose.total_residue()):
		spots.append(str(resi+1))
	hotspots = ':'.join(spots)
	#Setup score function
	scorefxn = get_fa_scorefxn()
	#Setup filters
	FLTR = rosetta.protocols.simple_filters.PackStatFilter()
	#Setup grafting mover
	graft = pyrosetta.rosetta.protocols.motif_grafting.movers.MotifGraftMover()
	graft.init_parameters(
						receptor,   # context_structure
						motif,      # motif_structure
						1.0,        # RMSD_tolerance
						2.0,        # NC_points_RMSD_tolerance
						5,          # clash_score_cutoff
						1,          # min_fragment_size
						'0:0',      # combinatory_fragment_size_delta
						'0:0',      # max_fragment_replacement_size_delta
						'ALA',      # clash_test_residue
						hotspots,   # hotspots
						True,       # full_motif_bb_alignment
						False,      # allow_independent_alignment_per_fragment
						True,       # graft_only_hotspots_by_replacement
						False,      # only_allow_if_N_point_match_aa_identity
						False,      # only_allow_if_C_point_match_aa_identity
						True,       # revert_graft_to_native_sequence
						False,      # allow_repeat_same_graft_output
						1.0,        # output_cluster_tolerance
						FLTR)       # output_filter
	graft.apply(scaffold)
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.apply(scaffold)
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
	Submits the pose to the Robetta server (http://www.robetta.org) for
	fragment generation that are used for the Abinitio folding simulation.
	Then measures the RMSD for each fragment at each position and chooses
	the lowest RMSD. Then averages out the lowest RMSDs. Then plots the
	lowest RMSD fragment for each positon. Generates the 3-mer file, the
	9-mer file, the PsiPred file, the RMSD vs Position PDF plot with the
	averaged fragment RMSD printed in the plot
	'''
	pose = pose_from_pdb(filename)
	sequence = pose.sequence()
	web = requests.get('http://www.robetta.org/fragmentsubmit.jsp')
	payload = {	'UserName':          username,
				'Email':             '',
				'Notes':             '{}'.format(filename.split('.')[0]),
				'Sequence':          sequence,
				'Fasta':             '',
				'Code':              '',
				'ChemicalShifts':    '',
				'NoeConstraints':    '',
				'DipolarConstraints':'',
				'type':              'submit'}
	session = requests.session()
	response = session.post('http://www.robetta.org/fragmentsubmit.jsp', data=payload, files=dict(foo='bar'))
	for line in response:
		line = line.decode()
		if re.search('<a href="(fragmentqueue.jsp\?id=[0-9].*)">', line):
			JobID = re.findall('<a href="(fragmentqueue.jsp\?id=[0-9].*)">', line)
	JobURL = 'http://www.robetta.org/' + JobID[0]
	ID = JobID[0].split('=')
	URL = 'http://robetta.org/fragmentqueue.jsp'
	page = urllib.request.urlopen(URL)
	data = bs4.BeautifulSoup(page, 'lxml')
	table = data.find('table', {'cellpadding':'3'})
	info = table.findAll('tr')
	print('''\u001b[32m{}
________      ______      __________
___  __ \________  /________  /__  /______ _
__  /_/ /  __ \_  __ \  _ \  __/  __/  __ `/
_  _, _// /_/ /  /_/ /  __/ /_ / /_ / /_/ /
/_/ |_| \____//_.___/\___/\__/ \__/ \__,_/

\u001b[34mRobetta.org Protein Structure Prediction Server\u001b[0m
'''.format('-'*57))
	print('\u001b[35m|JobID | Status | Length | User\u001b[0m')
	print('\u001b[33m-------------------------------------\u001b[0m')
	for item in info:
		try:
			job =      item.findAll('td')[0].findAll('a')[0].getText()
			status =   item.findAll('td')[1].getText()
			username = item.findAll('td')[3].getText()
			length =   item.findAll('td')[4].getText()
			target =   item.findAll('td')[5].getText()
		except: continue
		j = '\u001b[36m{}\u001b[33m'.format(job)
		s = '\u001b[36m{}\u001b[33m'.format(status)
		n = '\u001b[36m{}\u001b[33m'.format(username)
		l = '\u001b[36m{}\u001b[33m'.format(length)
		t = '\u001b[36m{}\u001b[33m'.format(target)
		if status == 'Queued':
			s =   '\u001b[31m{}\u001b[33m'.format(status)
			TBL = '\u001b[32m|{} |{}  | {}\t | {}\u001b[0m'.format(j, s, l, n)
		elif status == 'Active':
			s =   '\u001b[32m{}\u001b[33m'.format(status)
			TBL = '\u001b[33m|{} |{}  | {}\t | {}\u001b[0m'.format(j, s, l, n)
		else:
			TBL = '\u001b[33m|{} |{}| {}\t | {}\u001b[0m'.format(j, s, l, n)
		print(TBL)
	print('Fragments submitted to Robetta server [http://robetta.org/fragmentqueue.jsp]')
	print('Job ID: \u001b[32m{}\u001b[0m'.format(str(ID[1])))
	while True:
		Job = urllib.request.urlopen(JobURL)
		jobdata = bs4.BeautifulSoup(Job, 'lxml')
		status = jobdata.find('td', string='Status: ').find_next().text
		if status == 'Complete':
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M'), 'Status:', '\u001b[32m{}\u001b[0m'.format(status))
			break
		elif status == 'Active':
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M'), 'Status:', '\u001b[33m{}\u001b[0m'.format(status))
			time.sleep(180)
		else:
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M'), 'Status:', '\u001b[31m{}\u001b[0m'.format(status))
			time.sleep(300)
			continue
	sequence = pose.sequence()
	fasta = open('structure.fasta', 'w')
	fasta.write(sequence)
	fasta.close()
	time.sleep(1)
	os.system('wget -q http://www.robetta.org/downloads/fragments/'+str(ID[1])+'/aat000_03_05.200_v1_3')
	print('Downloading 3mer fragment file ...')
	os.system('wget -q http://www.robetta.org/downloads/fragments/'+str(ID[1])+'/aat000_09_05.200_v1_3')
	print('Downloading 9mer fragment file ...')
	os.system('wget -q http://www.robetta.org/downloads/fragments/'+str(ID[1])+'/t000_.psipred_ss2')
	print('Downloading PSIPRED file ...')
	os.rename('aat000_03_05.200_v1_3', 'frags.200.3mers')
	os.rename('aat000_09_05.200_v1_3', 'frags.200.9mers')
	os.rename('t000_.psipred_ss2', 'pre.psipred.ss2')
	frag = open('frags.200.9mers', 'r')
	data = open('RMSDvsPosition.dat', 'w')
	AVG = []
	for line in frag:
		if line.lstrip().startswith('position:'):
			line = line.split()
			size = line[1]
	frag.close()
	for i in range(1, int(size)):
		rmsd = []
		pose_copy = pyrosetta.Pose()
		pose_copy.assign(pose)
		frames = pyrosetta.rosetta.core.fragment.FrameList()
		fragset = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
		fragset.read_fragment_file('frags.200.9mers')
		fragset.frames(i, frames)
		movemap = MoveMap()
		movemap.set_bb(True)
		for frame in frames:
			for frag_num in range(1, frame.nr_frags()+1):
				frame.apply(movemap, frag_num, pose_copy)
				RMSD = rosetta.core.scoring.CA_rmsd(pose, pose_copy)
				rmsd.append(RMSD)
				lowest = min(rmsd)
				pose_copy.assign(pose)
		AVG.append(lowest)
		data.write(str(i)+'\t'+str(lowest)+'\n')
		print('\u001b[31mPosition:\u001b[0m {}\t\u001b[31mLowest RMSD:\u001b[0m {}\t|{}'.format(i, round(lowest, 3), '-'*int(lowest)))
	data.close()
	Average_RMSD = sum(AVG) / len(AVG)
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
	set label 'Average RMSD = {}' at graph 0.01 , graph 0.95 tc lt 7 font 'curior 12'\n
	plot 'RMSDvsPosition.dat' with boxes\n
	exit
	""".format(str(round(Average_RMSD, 3))))
	gnuplot.close()
	os.system('gnuplot < gnuplot_sets')
	os.remove('gnuplot_sets')
	print('\u001b[34mAverage RMSD:\u001b[0m {}'.format(round(Average_RMSD, 3)))
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
		''' Generate the resfile. '''
		AminoAcid = {	'A':129, 'P':159, 'N':195, 'H':224,
						'V':174, 'Y':263, 'C':167, 'K':236,
						'I':197, 'F':240, 'Q':225, 'S':155,
						'L':201, 'W':285, 'E':223, 'T':172,
						'M':224, 'R':274, 'G':104, 'D':193}
		self.filename = filename
		parser = Bio.PDB.PDBParser()
		structure = parser.get_structure('{}'.format(filename), filename)
		dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke')
		sasalist = []
		for aa in dssp:
			sasa = AminoAcid[aa[1]]*aa[3]
			if sasa <= 25:      sasa = 'C'
			elif 25 < sasa < 40:sasa = 'B'
			elif sasa >= 40:    sasa = 'S'
			if aa[2] == 'G' or aa[2] == 'H' or aa[2] == 'I':   ss = 'H'
			elif aa[2] == 'B' or aa[2] == 'E':                 ss = 'S'
			elif aa[2] == 'S' or aa[2] == 'T' or aa[2] == '-': ss = 'L'
			sasalist.append((aa[0], aa[1], ss, sasa))
		resfile = open('.resfile', 'a')
		#resfile.write('NATRO\nEX 1\nEX 2\nUSE_INPUT_SC\n')
		resfile.write('START\n')
		for n, r, a, s in sasalist:
			if s == 'S' and a == 'L':   line = '{} A PIKAA PGNQSTDERKH\n'.format(n)
			elif s == 'S' and a == 'H': line = '{} A PIKAA EHKQR\n'.format(n)
			elif s == 'S' and a == 'S': line = '{} A PIKAA DEGHKNPQRST\n'.format(n)
			elif s == 'B' and a == 'L': line = '{} A PIKAA ADEFGHIKLMNPQRSTVWY\n'.format(n)
			elif s == 'B' and a == 'H': line = '{} A PIKAA ADEHIKLMNQRSTVWY\n'.format(n)
			elif s == 'B' and a == 'S': line = '{} A PIKAA DEFHIKLMNQRSTVWY\n'.format(n)
			elif s == 'C' and a == 'L': line = '{} A PIKAA AFGILMPVWY\n'.format(n)
			elif s == 'C' and a == 'H': line = '{} A PIKAA AFILMVWY\n'.format(n)
			elif s == 'C' and a == 'S': line = '{} A PIKAA FILMVWY\n'.format(n)
			resfile.write(line)
		resfile.close()
		self.SASA = sasalist
		# aa_composition file
		with open('.comp', 'w')as comp:
			comp.write("""
				PENALTY_DEFINITION
				PROPERTIES AROMATIC
				NOT_PROPERTIES POLAR CHARGED
				FRACTION 0.1
				PENALTIES 100 0 100
				DELTA_START -1
				DELTA_END 1
				BEFORE_FUNCTION CONSTANT
				AFTER_FUNCTION CONSTANT
				END_PENALTY_DEFINITION
				""")
		# netcharge file
		with open('.charge', 'w')as comp:
			comp.write("""
				DESIRED_CHARGE 0
				PENALTIES_CHARGE_RANGE -1 1
				PENALTIES 10 0 10
				BEFORE_FUNCTION QUADRATIC
				AFTER_FUNCTION QUADRATIC
				""")
		self.pose = pose_from_pdb(self.filename)
		# pushback aa_composition
		comp = pyrosetta.rosetta.protocols.aa_composition.AddCompositionConstraintMover()
		comp.create_constraint_from_file('.comp')
		comp.apply(self.pose)
		# pushback netcharge
		charge = pyrosetta.rosetta.protocols.aa_composition.AddNetChargeConstraintMover()
		charge.create_constraint_from_file('.charge')
		charge.apply(self.pose)
		self.starting_pose = Pose()
		self.starting_pose.assign(self.pose)
		self.scorefxn = get_fa_scorefxn()
		self.scorefxn_G = get_fa_scorefxn()
		AAcomp      = pyrosetta.rosetta.core.scoring.ScoreType.aa_composition
		NETq        = pyrosetta.rosetta.core.scoring.ScoreType.netcharge
		AArep       = pyrosetta.rosetta.core.scoring.ScoreType.aa_repeat
		ASPpen      = pyrosetta.rosetta.core.scoring.ScoreType.aspartimide_penalty
		HBnet       = pyrosetta.rosetta.core.scoring.ScoreType.hbnet
		MHCep       = pyrosetta.rosetta.core.scoring.ScoreType.mhc_epitope
		VOIDpen     = pyrosetta.rosetta.core.scoring.ScoreType.voids_penalty
		ABurUnsatPen= pyrosetta.rosetta.core.scoring.ScoreType.approximate_buried_unsat_penalty
		BurUnsatPen = pyrosetta.rosetta.core.scoring.ScoreType.buried_unsatisfied_penalty
		#self.scorefxn_G.set_weight(AAcomp,      1.00)
		#self.scorefxn_G.set_weight(NETq,        1.00)
		#self.scorefxn_G.set_weight(HBnet,       1.00)
		#self.scorefxn_G.set_weight(VOIDpen,     0.10)
		self.scorefxn_G.set_weight(AArep,       1.00)
		self.scorefxn_G.set_weight(ASPpen,      1.00)
		self.scorefxn_G.set_weight(MHCep,       0.00)
		self.scorefxn_G.set_weight(BurUnsatPen, 1.00)
		self.scorefxn_G.set_weight(ABurUnsatPen,5.00)
		self.relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		self.relax.set_scorefxn(self.scorefxn)
	def __del__(self):
		''' Remove the resfile. '''
		os.remove('.resfile')
		os.remove('.comp')
		os.remove('.charge')
		for f in glob.glob('f[il]xbb.fasc'): os.remove(f)
	def choose(self):
		''' Choose the lowest scoring structure. '''
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
		Generates the structure.pdb file.
		'''
		resfile = pyrosetta.rosetta.core.pack.task.operation.ReadResfile('.resfile')
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()
		task.push_back(resfile)
		movemap = MoveMap()
		movemap.set_bb(False)
		movemap.set_chi(True)
		fixbb = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
		fixbb.set_task_factory(task)
		fixbb.set_movemap(movemap)
		fixbb.set_scorefxn(self.scorefxn_G)
		self.relax.apply(self.pose)
		job = PyJobDistributor('fixbb', 10, self.scorefxn)
		job.native_pose = self.starting_pose
		while not job.job_complete:
			self.pose.assign(self.starting_pose)
			fixbb.apply(self.pose)
			self.relax.apply(self.pose)
			job.output_decoy(self.pose)
	def flxbb(self):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while allowing for a flexible backbone.
		Generates the structure.pdb file.
		'''
		resfile = pyrosetta.rosetta.core.pack.task.operation.ReadResfile('.resfile')
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()
		task.push_back(resfile)
		movemap = MoveMap()
		movemap.set_bb(True)
		movemap.set_chi(True)
		flxbb = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
		flxbb.set_task_factory(task)
		flxbb.set_movemap(movemap)
		flxbb.set_scorefxn(self.scorefxn_G)
		self.relax.apply(self.pose)
		job = PyJobDistributor('flxbb', 10, self.scorefxn)
		job.native_pose = self.starting_pose
		while not job.job_complete:
			self.pose.assign(self.starting_pose)
			flxbb.apply(self.pose)
			self.relax.apply(self.pose)
			job.output_decoy(self.pose)
	def fixbb_motif(self, Motif_From, Motif_To):
		'''
		Applies RosettaDesign with a fixed back bone to
		change the structure's amino acids (one layer at a
		time - like in the Design_Layer method) except for
		a desired continuous motif sequence while maintaining
		the same backbone.
		'''
		# Remove motif residues from resfile
		resfile = open('.resfile', 'r')
		resfile2 = open('.resfile2', 'a')
		resfile2.write('START\n')
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
		resfile = rosetta.core.pack.task.operation.ReadResfile('.resfile2')
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()
		task.push_back(resfile)
		movemap = MoveMap()
		movemap.set_bb(False)
		movemap.set_chi(True)
		fixbb = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
		fixbb.set_task_factory(task)
		fixbb.set_movemap(movemap)
		fixbb.set_scorefxn(scorefxn)
		job = PyJobDistributor('fixbb', 10, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			fixbb.apply(pose)
			relax.apply(pose)
			job.output_decoy(pose)
		os.remove('.resfile2')
	def flxbb_motif(self, Motif_From, Motif_To):
		'''
		Applies RosettaDesign with a flexible back bone to
		change the structure's amino acids (one layer at a
		time - like in the Design_Layer method) except for
		a desired continuous motif sequence while maintaining
		the same backbone.
		'''
		# Remove motif residues from resfile
		resfile = open('.resfile', 'r')
		resfile2 = open('.resfile2', 'a')
		resfile2.write('START\n')
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
		resfile = rosetta.core.pack.task.operation.ReadResfile('.resfile2')
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()
		task.push_back(resfile)
		movemap = MoveMap()
		movemap.set_bb(True)
		movemap.set_chi(True)
		flxbb = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
		flxbb.set_task_factory(task)
		flxbb.set_movemap(movemap)
		flxbb.set_scorefxn(scorefxn)
		job = PyJobDistributor('flxbb', 10, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			flxbb.apply(pose)
			relax.apply(pose)
			job.output_decoy(pose)
		os.remove('.resfile2')
	def surf(self, motif_list):
		'''
		Applies RosettaDesign with a fixed backbone to change only the
		structure's surface amino acids except for the desired motif.
		'''
		pose = pose_from_pdb(self.filename)
		starting_pose = Pose()
		starting_pose.assign(pose)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		packtask = standard_packer_task(pose)
		pyrosetta.rosetta.core.pack.task.parse_resfile(pose, packtask, '.resfile')
		#Identify non-surface and motif residues
		Motif = motif_list
		for s in self.SASA:
			if s[3] != 'S':
				Motif.append(s[0])
		for aa in Motif:
			packtask.temporarily_set_pack_residue(int(aa), False)
		pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, packtask)
		fixbb = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, packtask, 10)
		job = PyJobDistributor('fixbb', 10, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			relax.apply(pose)
			fixbb.apply(pose)
			relax.apply(pose)
			job.output_decoy(pose)

def FFL(Motif, Scaffold, Motif_From, Motif_To, username):
	''' Performs the Fold From Loops protocol '''
	print('\x1b[31m[-] FFL is not yet fully available in PyRosetta\x1b[0m')

def protocol(Protein, RChain, Chain, Motif_from, Motif_to, Scaffold, Choice, UserName):
	#0. Make directory
	os.makedirs('Vaccine', exist_ok=True)
	os.chdir('Vaccine')
	print('\x1b[32m[+] Project directory created\x1b[0m')
	#1. Import scaffold
	os.system('mv ../{} .'.format(Scaffold))
	pose = pose_from_pdb(Scaffold)
	print('\x1b[32m[+] Imported scaffold\x1b[0m')
	#2. Isolate motif
	Motif(Protein, Chain, Motif_from, Motif_to)
	print('\x1b[32m[+] Isolated motif\x1b[0m')
	#3. Isolate receptor
	Receptor(Protein, RChain)
	print('\x1b[32m[+] Isolated receptor\x1b[0m')
	#4. Graft motif onto scaffold
	print('\x1b[32m[+] Grafting...\x1b[0m')
	MotifPosition = Graft('receptor.pdb', 'motif.pdb', pose)
	print('\x1b[32m[+] Grafted motif onto scaffold between positions: {} and {}\x1b[0m'.format(MotifPosition[0], MotifPosition[1]))
	#5. Fold From Loop
	#FFL('motif.pdb', 'grafted.pdb', MotifPosition, UserName)
	#print('\x1b[32m[+] Fold From Loop completed\x1b[0m')
	#6. RosettaDesign the structure around the motif
	if Choice == 'fixbb':
		print('\x1b[32m[+] Fixbb designing...\x1b[0m')
		RD = RosettaDesign('grafted.pdb')
		RD.fixbb_motif(MotifPosition[0], MotifPosition[1])
	elif Choice == 'flxbb':
		print('\x1b[32m[+] Flxbb designing...\x1b[0m')
		RD = RosettaDesign('grafted.pdb')
		RD.flxbb_motif(MotifPosition[0], MotifPosition[1])
	print('\x1b[32m[+] Design complete\x1b[0m')
	#7. Generate and analyse fragments
	X = [3, 3, 3, 3, 1, 3, 3, 2, 3, 3]
	for i in range(10):
		RMSD = Fragments('{}_{}.pdb'.format(f, str(i)), UserName)
		if RMSD < 2.0:
			for f in glob.glob('f[il]xbb_{}.pdb'.format(str(i))): os.system('mv {} structure.pdb'.format(f))
			for f in glob.glob('f[il]xbb_*'): os.remove(f)
			print('\x1b[32m[+++] Vaccine structure completed\x1b[0m')
			exit()
		else:
			for f in glob.glob('f[il]xbb_{}.pdb'.format(str(i))): os.remove(f)
	print('\x1b[31m[---] Vaccine structure failed\x1b[0m')

def main():
	if args.scaffold:   # Search for scaffolds
		print('\x1b[32m[+] Searching for scaffolds\x1b[0m')
		ScaffoldSearch(	sys.argv[2],         # PDB ID
						sys.argv[3],         # Receptor chain
						sys.argv[4],         # Motif chain
						sys.argv[5],         # Motif from
						sys.argv[6],         # Motif to
						sys.argv[7])         # Directory of scaffolds
	elif args.protocol: # Run full protocol
		protocol(		sys.argv[2],         # PDB ID
						sys.argv[3],         # Receptor chain
						sys.argv[4],         # Motif chain
						sys.argv[5],         # Motif from
						sys.argv[6],         # Motif to
						sys.argv[7],         # Scaffold PDB file name
						sys.argv[8],         # RosettaDesign choice
						sys.argv[9])         # Robetta server username
	elif args.motif:    # Isolate motif
		Motif(			sys.argv[2],         # PDB ID
						sys.argv[3],         # Motif chain
						sys.argv[4],         # Motif from
						sys.argv[5])         # Motif to
		os.remove('{}.pdb'.format(sys.argv[2]))
		print('\x1b[32m[+] Motif isolated\x1b[0m')
	elif args.receptor: # Isolate receptor
		RCSB = 'http://www.rcsb.org/pdb/files'
		os.system('wget -q {}/{}.pdb'.format(RCSB, sys.argv[2]))
		Receptor(		sys.argv[2],         # PDB ID
						sys.argv[3])         # Receptor chain
		print('\x1b[32m[+] Receptor isolated\x1b[0m')
	elif args.graft:    # Graft motif onto scaffold
		pose = pose_from_pdb(	sys.argv[4]) # Scaffold PDB file name
		MotifPosition = Graft(	sys.argv[2], # Receptor PDB file name
								sys.argv[3], # Motif PDB file name
								pose)
		print('\x1b[32m[+] Grafted motif onto scaffold between positions: {} and {}\x1b[0m'.format(MotifPosition[0], MotifPosition[1]))
	elif args.ffl:      # Fold From Loop
		FFL(			sys.argv[2],         # Motif PDB file name
						sys.argv[3],         # Scaffold PDB file name
						sys.argv[4],         # Motif on scaffold from
						sys.argv[5],         # Motif on scaffold to
						sys.argv[6])         # Robetta server username
		#print('\x1b[32m[+] Fold From Loop completed\x1b[0m')
	elif args.design:   # Sequence design the structure around the motif
		if sys.argv[2] == 'fixbb':           # Choice
			print('\x1b[32m[+] Fixbb designing\x1b[0m')
			RD = RosettaDesign(	sys.argv[3]) # Scaffold PDB file name
			RD.fixbb_motif(		sys.argv[4], # Motif on scaffold from
								sys.argv[5]) # Motif on scaffold to
			print('\x1b[32m[+] Design complete\x1b[0m')
		elif sys.argv[2] == 'flxbb':         # Choice
			print('\x1b[32m[+] Flxbb designing\x1b[0m')
			RD = RosettaDesign(	sys.argv[3]) # Scaffold PDB file name
			RD.flxbb_motif(		sys.argv[4], # Motif on scaffold from
								sys.argv[5]) # Motif on scaffold to
			print('\x1b[32m[+] Design complete\x1b[0m')
		elif sys.argv[2] == 'surface':       # Choice
			print('\x1b[32m[+] Surface designing\x1b[0m')
			RD = RosettaDesign(	sys.argv[3]) # Scaffold PDB file name
			RD.surf(			sys.argv[4:])# Motif amino acid list
			print('\x1b[32m[+] Design complete\x1b[0m')
	elif args.fragments:# Generate fragments
		Fragments(		sys.argv[2],         # Filename
						sys.argv[3])         # Username
		print('\x1b[32m[+] Fragments generated\x1b[0m')

if __name__ == '__main__': main()
