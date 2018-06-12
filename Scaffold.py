#!/usr/bin/python3

from pyrosetta import *
from pyrosetta.toolbox import *
init()

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
		try:
			if int(Motif_From) <= int(line.split()[5]) <= int(Motif_To):	#Find residues between the user specified location
				count += 1						#Sequencially number atoms
				AA1 = line[23:27]					#Sequencially number residues
				if not AA1 == AA2:
					num += 1			
				final_line = line[:7] + '{:4d}'.format(count) + line[11:17] + line[17:21] + 'A' + '{:4d}'.format(num) + line[26:]	#Update each line of the motif to have its atoms and residues sequencially labeled, as well as being in chain A
				AA2 = AA1
				Motif.write(final_line)					#Write to new file called motif.pdb
		except:
			continue
	Motif.close()

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



















def Graft_SC(receptor, motif, scaffold):
	'''
	Side chain grafting. Grafts a motif's side chains onto a protein
	scaffold structure.
	Generates structure.pdb and returns a tuple [0] is the
	residue number where the motif starts and [1] where it ends
	'''
	scorefxn = get_fa_scorefnx()
	#Task Operations
	pido = ProteinInterfaceDesign
	repack_chain1="1"
	repack_chain2="1"
	design_chain1="0"
	design_chain2="1"
	interface_distance_cutoff="8.0"



	hotspot_repack = OperateOnCertainResidues
		ResiduePDBInfoHasLabel
		property="HOTSPOT"		
		RestrictToRepackingRLT

	#Filters
	ddg = Ddg
	confidence="0"

	unsat = BuriedUnsatHbonds
	confidence="0"

	Sc = ShapeComplementarity
	confidence="0"

	#Movers

	motif_grafting = MotifGraft
	context_structure="context.pdb"
	motif_structure="motif.pdb"
	RMSD_tolerance="0.3"
	NC_points_RMSD_tolerance="0.5"
	clash_score_cutoff="5"
	clash_test_residue="GLY"
	hotspots="3:7"
	combinatory_fragment_size_delta="2:2"
	full_motif_bb_alignment="1"
	graft_only_hotspots_by_replacement="1"
	revert_graft_to_native_sequence="1"

	ala_pose = build_Ala_pose
	partner1="0"
	partner2="1"
	interface_cutoff_distance="8.0"
	task_operations="hotspot_repack"

	ppk = Prepack
	jump_number="0"

	design = PackRotamersMover
	task_operations="hotspot_repack,pido"

	rb_min = MinMover
	bb="0"
	chi="1"
	jump="1"

	#Protocol
	protocol = pyrosetta.rosetta.protocols.moves.SequenceMover()
	protocol.add_mover(motif_grafting)
	protocol.add_mover(ala_pose)
	protocol.add_mover(ppk)
	protocol.add_mover(design)
	protocol.add_mover(rb_min)
	protocol.add_mover(design)
	protocol.add_filter(unsat)
	protocol.add_filter(ddg)
	protocol.add_filter(Sc)
	protocol.apply(pose)








def Graft_BB(receptor, motif, scaffold):
	'''
	Back bone grafting. Grafts a motif's back bone onto a protein
	scaffold structure.
	Generates structure.pdb and returns a tuple [0] is the
	residue number where the motif starts and [1] where it ends
	'''



<ROSETTASCRIPTS>
<TASKOPERATIONS>
	<ProteinInterfaceDesign name="pido_far" interface_distance_cutoff="15.0"/>
	<ProteinInterfaceDesign name="pido_med" interface_distance_cutoff="12.0"/>
	<ProteinInterfaceDesign name="pido_near" interface_distance_cutoff="8.0"/>
	<OperateOnCertainResidues name="hotspot_repack">
		<ResiduePDBInfoHasLabel property="HOTSPOT"/>
		<RestrictToRepackingRLT/>
	</OperateOnCertainResidues>
	<SelectBySASA name="core" mode="sc" state="bound" probe_radius="2.2" core_asa="0" surface_asa="30" core="1" boundary="0" surface="0"/>
	<SelectBySASA name="core_and_boundary"  mode="sc" state="bound" probe_radius="2.2" core_asa="0" surface_asa="30" core="1" boundary="1" surface="0"/>
</TASKOPERATIONS>
<SCOREFXNS>
</SCOREFXNS>
<FILTERS>
	<Ddg name="ddg" confidence="0"/>
	<BuriedUnsatHbonds name="unsat" confidence="0"/>
	<ShapeComplementarity name="Sc" confidence="0"/>
</FILTERS>
<MOVERS>
	<MotifGraft name="motif_grafting" context_structure="context.pdb" motif_structure="motif.pdb" RMSD_tolerance="1.0" NC_points_RMSD_tolerance="1.0" 
	 clash_score_cutoff="5" clash_test_residue="GLY" hotspots="3:7" combinatory_fragment_size_delta="2:2" max_fragment_replacement_size_delta="-8:8" full_motif_bb_alignment="0" graft_only_hotspots_by_replacement="0"/>
	 <build_Ala_pose name="ala_pose" partner1="0" partner2="1" interface_cutoff_distance="8.0" task_operations="hotspot_repack"/>
	 <Prepack name="ppk" jump_number="0"/>
	 <PackRotamersMover name="design_core" task_operations="hotspot_repack,pido_far,core"/>
	 <PackRotamersMover name="design_boundary" task_operations="hotspot_repack,pido_med,core_and_boundary"/>
	 <PackRotamersMover name="design_interface" task_operations="hotspot_repack,pido_near"/>
	 <MinMover name="sc_min" bb="0" chi="1" jump="1"/>
</MOVERS>
<PROTOCOLS>
	<Add mover_name="motif_grafting"/>
	<Add mover_name="ala_pose"/>
	<Add mover_name="ppk"/>
	<Add mover_name="design_core"/>
	<Add mover_name="design_boundary"/>
	<Add mover_name="design_interface"/>
	<Add mover_name="sc_min"/>
	<Add filter_name="unsat"/>
	<Add filter_name="ddg"/>
	<Add filter_name="Sc"/>
</PROTOCOLS>
</ROSETTASCRIPTS>









	pose = pose_from_pdb(scaffold)
	scorefxn = get_fa_scorefxn()

	#Task Operations
	pido_far = ProteinInterfaceDesign()
	interface_distance_cutoff="15.0"

	pido_med = ProteinInterfaceDesign()
	interface_distance_cutoff="12.0"

	pido_near = ProteinInterfaceDesign()
	interface_distance_cutoff="8.0"

	hotspot_repack = OperateOnCertainResidues()
	ResiduePDBInfoHasLabel property="HOTSPOT
	RestrictToRepackingRLT

	<SelectBySASA name="core" mode="sc" state="bound" probe_radius="2.2" core_asa="0" surface_asa="30" core="1" boundary="0" surface="0"/>
	<SelectBySASA name="core_and_boundary"  mode="sc" state="bound" probe_radius="2.2" core_asa="0" surface_asa="30" core="1" boundary="1" surface="0"/>








	print(scorefxn(scaffold))
	mover = pyrosetta.rosetta.protocols.motif_grafting.movers.MotifGraftMover()
	#Setup motif hotspots
	motifpose = pose_from_pdb(motif)
	spots = list()
	for resi in range(motifpose.total_residue()):
		spots.append(str(resi+1))
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

#Motif('2y7q' , 'B' , 420 , 429)
#Receptor('2y7q' , 'A')
#Graft_SC('receptor.pdb', 'motif.pdb', 'scaffold.pdb')
#Graft_BB('receptor.pdb', 'motif.pdb', 'scaffold.pdb')