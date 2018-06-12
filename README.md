# VexDesign
Scripts that autonomously designs a vaccine. Authored by Sari Sabban on 31-May-2017 (sari.sabban@gmail.com).

## Requirements:
1. Make sure you install [PyRosetta](http://www.pyrosetta.org) as the website describes.
2. Use the following commands (in GNU/Linux) to install all nessesary programs and Python libraries for this script to run successfully:

`sudo apt update && sudo apt install pymol gnuplot python3-bs4 python3-biopython python3-lxml -y`

## How To Use:
1. Go to [Robetta server](http://robetta.org/) and register a username, it is very easy and straightforward.
2. Generate a protein scaffold using the following command:

`python3 Scaffold.py PDBID RCHAIN CHAIN FROM TO`

* PDBID = The protein's [Protein Data Bank](https://www.rcsb.org) identification name
* RCHAIN = The chain where your receptor resides within the protein .pdb file
* CHAIN = The chain where your target site resides (not part of the receptor) within the protein .pdb file
* FROM = The start of your target site
* TO = The end of your target site

Example:

`python3 Scaffold.py 2y7q A B 420 429`

3. Generate a vaccine using the following command:

`python3 VexDesign.py PDBID RCHAIN CHAIN FROM TO SCAFFOLD USERNAME`

* PDBID = The protein's [Protein Data Bank](https://www.rcsb.org) identification name
* RCHAIN = The chain where your receptor resides within the protein .pdb file
* CHAIN = The chain where your target site resides (not part of the receptor) within the protein .pdb file
* FROM = The start of your target site
* TO = The end of your target site
* SCAFFOLD = The protein scaffold structure generated from the Scaffold.py script
* USERNAME = The [Robetta server](http://robetta.org/) username to generate and download the *Abinitio* fragment files and the PSIPRED secondary structure prediction file

Example:

`python3 VexDesign.py 2y7q A B 420 429 scaffold.pdb`

4. Calculation time is about 72 hours on a normal desktop computer.
5. Access to the internet is a requirement since the VexDesign.py script will be sending and retrieving data from some servers.
6. Use this [Rosetta Abinitio](https://github.com/sarisabban/RosettaAbinitio) script to simulate the folding of the final designed vaccine's protein structure. An HPC (High Preformance Computer) and the original C++ [Rosetta](https://www.rosettacommons.org/) are required for this step.

## Description
These scripts autonomously designs a vaccine from a user specified target site. This is not artificial intellegance, you cannot just ask these scripts to design "A" vaccine, you must understand what target site you want to develop antibodies against (make a liturature search and understand your disease and target site), then supply this target site to these scripts to build a protein structure around it so the final protein can be used as a vaccine. You must have prior understanding of Bioinformatics and Immunology in order to be able to understand what site to target and to supply to these scripts. Once you identify a target site, these scripts will take it and run a graft and desing protocol, without the need for you to intervene, that will result in an ideal protein structure displaying your target site in its original 3D cofiguration. Thus the protien, theoretically, can be used as a vaccine against this site, and hopefully neutralise the disease you are researching. Everytime you run the VexDesign.py script a different final protien sequence will be generated for the **same structure**, this is important to keep in mind, because if you want to generate many different structures with different sequences to test or to use as boosts you can simply run the same target site again and you will end up with the same structure but with a different sequence.

These scripts have been last tested to work well with PyRosetta 4 Release 177 using Python 3.6. If you use this script on a newer PyRosetta or Python version and it fails please notify me of the error and I will update it.

Here is a [video](youtube.com/) that explains how to select a target site, how the script functions, and what results you should get. If I did not make a video yet, bug me until I make one.

The protocol is as follows:

*Scaffold.py*

1. Generates a scaffold structure. **--> STILL UNDER DEVELOPMENT**

*VexDesign.py*

1. Isolates the motif.
2. Isolates the receptor.
3. Grafts the motif onto the scaffold.
4. Sequence designs the structure around the motif.
5. Generate fragments for Rosetta *Abinitio* folding simulation.

Output files are as follows:

|    | File Name               | Description                                                                                  |
|----|-------------------------|----------------------------------------------------------------------------------------------|
| 1  | scaffold.pdb            | Scaffold structure                                                                           |
| 2  | motif.pdb	       | Original requested motif                                                                     |
| 3  | receptor.pdb            | Original receptor that binds motif                                                           |
| 4  | grafted.pdb             | Grafted motif to scaffold structure                                                          |
| 5  | structure.pdb           | Sequence designed structure                                                                  |
| 6  | structure.fasta         | Fasta of Rosetta Designed structure                                                          |
| 7  | frags.200.3mers         | 3-mer fragment of sequence designed structure from the Robetta server                        |
| 8  | frags.200.9mers         | 9-mer fragment of sequence designed structure from the Robetta server                        |
| 9  | pre.psipred.ss2         | PSIPRED secondary structure prediction of sequence designed structure from the Robetta server|
| 10 | plot_frag.pdf           | Plot of the fragment quality RMSD vs Position                                                |
| 11 | RMSDvsPosition.dat      | plot_frag.pdf's data                                                                         |
| 12 | FragmentAverageRMSD.dat | Average RMSD of the fragments                                                                |

## Simple Design
If you have a strucutre that you want to target, but this structure is between 80 and 150 amino acids and is a ridgid structure (mainly helices and sheets and very little loops), you might get away with just re-designing the structure itself, while keeping the target site(s) conserved, without having to do anything complicated like isolate the motif and graft it onto a scaffold protein (as in the previous two scripts). The code for this protocol is very simple just use the following command:

`python3 Simple.py PDBID CHAIN NO_DESIGN_POSITIONS`

* PDBID = The protein's [Protein Data Bank](https://www.rcsb.org) identification name
* CHAIN = The chain of your target structure that you want to design
* NO_DESIGN_POSITIONS = The amino acid positions that you want to conserve (not re-design), i.e: your motif(s)

Example:

`python3 Simple.py 3WD5 A 19 20 21 22 23 24 65 66 67 68 140 141 143 144 145 146 147 148`

The script's protocol is as follows:
1. Gets the user defined protein structure and isolates the defined chain.
2. Sequence designs the structure (but not the user defined conserved amino acids). If the NO_DESIGN_POSITIONS is left empty the entire protein will be designed (good as a negative control).

The script will not access the Robetta server nor setup any of the *Abinitio* folding simulation files, it is meant to be used with more manual work from the user than automated work. Calculation time is about 24 hours.

Output files are as follows:

|    | File Name               | Description                                                                                  |
|----|-------------------------|----------------------------------------------------------------------------------------------|
| 1  | original.pdb            | The original non-designed structure (for comparisons)                                        |
| 2  | structure.pdb	       | The designed structure                                                                       |