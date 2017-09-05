# VexDesign
A script that autonomously designs a vaccine. Authored by Sari Sabban on 31-May-2017 (sari.sabban@gmail.com).

## Requirements:
1. Make sure you install [PyRosetta](http://www.pyrosetta.org) as the website describes.
2. Use the following commands (in GNU/Linux) to install all necessary Python libraries for this script to run successfully:

`sudo apt install python3-pip pymol DSSP gnuplot && sudo python3 -m pip install zeep numpy biopython bs4`

3. Download the vall.jul19.2011.gz database (467 MB). This link is temporary until the database is included with PyRosetta:

`wget https://www.dropbox.com/s/4tcpq5vscqst5ww/vall.jul19.2011.gz`

## How To Use:
1. Use the following command to run the script:

`python3 VaxDesign.py PDBID RCHAIN CHAIN FROM TO VALL`

* PDBID = The protein's [Protein Data Bank](https://www.rcsb.org) identification name
* RCHAIN = The chain where your receptor resides within the protein .pdb file
* CHAIN = The chain where your target site resides (not part of the receptor) within the protein .pdb file
* FROM = The start of your target site
* TO = The end of your target site
* VALL = The path to the vall.jul19.2011.gz database

Example:

`python3 VaxDesign.py 2y7q A B 420 429 /home/acresearch/rosetta_src_2017.08.59291_bundle/tools/fragment_tools/vall.jul19.2011.gz`

2. Calculation time is about 720 hours on a normal desktop computer.
3. Access to the internet is a requirement since the script will be sending and retrieving data from some servers.
4. Use [this Rosetta Abinitio script](https://github.com/sarisabban/RosettaAbinitio) to simulate the folding of the final designed vaccine's protein structure. An HPC (High Preformance Computer) and the original C++ [Rosetta](https://www.rosettacommons.org/) are required for this step.

## Description
This script autonomously designs a vaccine from a user specified target site. This is not artificial intellegance, you cannot just ask the the script to design "A" vaccine, you must understand what target site you want to develop antibodies against (make a liturature search and understand your disease and target site), then supply this target site to the script to build a protein structure around it so the final protein can be used as a vaccine. You must have prior understanding of Bioinformatics and Immunology in order to be able to understand what site to target and to supply it to the script. Once you identify a target site, the script will take it and run a long protocol, without the need for you to intervene, that will result in an ideal protein structure displaying your target site in its original 3D cofiguration. Thus the protien, theoretically, can be used as a vaccine against this site, and hopefully neutralise the disease you are researching. Everytime you run this script a different final protien structure will be generated, this is important to keep in mind, because if you want to generate many different structures to test or to use as boosts you can simply run the same target site again and you will end up with a different final structure.

This script has been last tested to work well with PyRosetta 4 Release 147 and using Python 3.5. If you use this script on a newer PyRosetta or Python version and it fails please notify the author to get it updated.

Here is a [video](youtube.com/) that explains how to select a target site, how the script functions, and what results you sould get. If I did not make a video yet, bug me until I make one.

The script protocol is as follows:
1. Build Scaffold. --> STILL UNDER DEVELOPMENT --> I am having lots of trouble with De Novo Design (I have a very long temporary work around)
2. Isolate Motif.
3. Isolate Receptor.
4. Graft Motif onto Scaffold.
5. Sequence Design The Structure Around The Motif.
6. Generate Fragments for Rosetta Abinitio Folding Simulation.
7. If Average Fragment RMSD is Higher Than 2Å Repeat Steps 5 and 6.

Output files are as follows:

|    | File Name               | Description                                                                              |
|----|-------------------------|------------------------------------------------------------------------------------------|
| 1  | DeNovo.pdb              | Scaffold structure                                                                       |
| 2  | motif.pdb	             | Original requested motif                                                                 |
| 3  | receptor.pdb            | Original receptor that binds morif                                                       |
| 4  | grafted.pdb             | Grafted motif to De Novo structure                                                       |
| 5  | structure.pdb           | Sequence designed structure                                                              |
| 6  | structure.fasta         | Fasta of Rosetta Designed structure                                                      |
| 7  | frags.200.3mers         | 3-mer fragment of sequence designed structure from the Robetta server                    |
| 8  | frags.200.9mers         | 9-mer fragment of sequence designed structure from the Robetta server                    |
| 9  | pre.psipred.ss2     | PSIPRED secondary structure prediction of sequence designed structure from the Robetta server|
| 10 | plot_frag.pdf           | Plot of the fragment quality RMSD vs Position                                            |
| 11 | FragmentAverageRMSD.dat | Average RMSD of the fragments                                                            |
