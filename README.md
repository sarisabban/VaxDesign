# VaxDesign
Script that computationally designs a vaccine.

## Requirements
1. Download this script:

`wget https://raw.githubusercontent.com/sarisabban/VaxDesign/master/VaxDesign.py`

2. Use the following commands (in Ubuntu GNU/Linux) to install all nessesary programs and Python libraries for this script to run successfully:

`sudo apt update && sudo apt install gnuplot dssp python3-bs4 python3-biopython python3-lxml -y`

3. Make sure you install [PyRosetta](http://www.pyrosetta.org) as the website describes.

## Description
This script autonomously designs a vaccine from a user specified target site. It is best to use structures that contain a target site bound to a receptor or an antibody since these structures will be used to guide the vaccine development. Once you identify a target site, this script will graft it , without the need for you to intervene, onto a scaffold protein structure and design it resulting in an ideal protein structure displaying your target site in its original 3D cofiguration. Thus the protien, theoretically, can be used as a vaccine against this site, and hopefully neutralise the disease you are researching. Running the script willr esult in multiple protein structures with different sequences. It is best to use these structures and evaluate them (such as forward fold them using the *abinitio* protocol) to choose the best one to synthesise.

This script has been last tested to work well with PyRosetta 4 Release 223 using Python 3.7. If you use this script on a newer PyRosetta or Python version and it fails please notify me of the error and I will update it.

The script does the following:

*Scaffold Search*

1. Chooses ideal scaffold structures from a scaffold database.

*Vaccine Design*

1. Isolates the motif.

2. Isolates the receptor.

3. Grafts the motif onto a specified scaffold.

4. Performs the Fold From Loop protocol to optimise the grafted structure (NOT IMPLEMENTED YET).

5. RosettaDesign the structure around the motif to find the sequence that would fold the proteins to its designed structure.

6. Generate fragments for Rosetta AbinitioRelax folding simulation and analyses their quality.

Output files are as follows for each generated structure:

|    | File Name               | Description                                                                                  |
|----|-------------------------|----------------------------------------------------------------------------------------------|
| 1  | frags.200.3mers         | 3-mer fragment of sequence designed structure from the Robetta server                        |
| 2  | frags.200.9mers         | 9-mer fragment of sequence designed structure from the Robetta server                        |
| 3  | grafted.pdb             | Grafted motif to scaffold structure                                                          |
| 4  | motif.pdb               | Original requested motif                                                                     |
| 5  | pre.psipred.ss2         | PSIPRED secondary structure prediction of sequence designed structure from the Robetta server|
| 6  | plot_frag.pdf           | Plot of the fragment quality RMSD vs Position                                                |
| 7  | receptor.pdb            | Original receptor that binds motif                                                           |
| 8  | scaffold.pdb            | Scaffold structure                                                                           |
| 9  | structure.fasta         | Fasta of Rosetta Designed structure                                                          |
| 10 | **structure.pdb**       | Final RosettaDesigned vaccine structure                                                      |

## Manual
### Autonmous Structure Generation
1. Go to [Robetta server](http://robetta.org/) and register a username, it is very easy and straightforward.

2. Search potential protein scaffolds using the following command:

`python3 VaxDesign.py -s PDBID RCHAIN CHAIN FROM TO DATABASE` use `-s or --scaffold`

* PDBID = The protein's [Protein Data Bank](https://www.rcsb.org) identification name
* RCHAIN = The chain where the receptor resides within the protein .pdb file
* CHAIN = The chain where the target site resides (not part of the receptor) within the protein .pdb file
* FROM = The start of the target site
* TO = The end of the target site
* DATABASE = The directory that contains multiple .pdb scaffold structures to search through

3. The script will generate a directory called **Scaffolds** which will contain all the scaffold structures. Search through the directory and choose one scaffold (move it to the working directory next to the VaxDesign.py script).

4. Generate a vaccine using the following command:

`python3 VaxDesign.py -p PDBID RCHAIN CHAIN FROM TO SCAFFOLD PROTOCOL USERNAME` use `-p or --protocol`

* PDBID = The protein's [Protein Data Bank](https://www.rcsb.org) identification name
* RCHAIN = The chain where the receptor resides within the protein .pdb file
* CHAIN = The chain where the target site resides (not part of the receptor) within the protein .pdb file
* FROM = The start of the target site
* TO = The end of the target site
* SCAFFOLD = The protein scaffold structure generated from the Scaffold.py script
* PROTOCOL = Choose between fixbb ot flxbb RosettaDesign
* USERNAME = The [Robetta server](http://robetta.org/) username to generate and download the AbinitioRelax fragment files and the PSIPRED secondary structure prediction file

Calculation time is about 12 hours on a normal desktop computer. Access to the internet is a requirement since the VaxDesign.py script will be sending and retrieving data from some servers. The script will result in 20 vaccine structures.

5. Use this [Rosetta Abinitio](https://github.com/sarisabban/RosettaAbinitio) script to simulate the folding of all final designed vaccine structures. An HPC (High Preformance Computer) and the original C++ [Rosetta](https://www.rosettacommons.org/) are required for this step. This is just an evaluation step to check the design before you synthesise it, other folding simulation algorithms can be used, use what you find most reliable.

### Manual Structure Generation
You can run each step separatly, use the following commands to run each step:

* Isolate motif:

`python3 VaxDesign.py -m PDBID RCHAIN FROM TO` use `-m or --motif`

* Isolate receptor:

`python3 VaxDesign.py -r PDBID RCHAIN` use `-r or --receptor`

* Graft motif onto the scaffold structure

`python3 VaxDesign.py -g RECEPTOR.pdb MOTIF.pdb SCAFFOLD.pdb` use `-g or --graft`

* Run the Fold From Loop protocol (Not Implemented Yet)

`python3 VaxDesign.py -f motif.pdb grafted.pdb FROM TO USERNAME` use `-f or --ffd`

In this case the FROM and TO are the positions of the motif on the grafted structure and not the original protien structure

* RosettaDesign the structure

`python3 VaxDesign.py -d PROTOCOL grafted.pdb FROM TO` use `-d or --design`

In this case the FROM and TO are the positions of the motif on the grafted structure and not the original protien structure. PROTOCOL can be either `fixbb` or `flxbb`, with fixbb being the reccomended one.

You can also choose to RosettaDesign only the suface of the structure (without changing your motif), use `surface` as the protocol.

* Generate fragments

`python3 VaxDesign.py -F structure.pdb USERNAME` use `-F or --fragments`

* Help menu

`python3 VaxDesign.py -h` use `-h or --help`

### Tutorial
Here is a tutorial that walks you through how to use the script and the results that can be expected. Here is also a [video](https://youtu.be/rZl9yAI6jzc) that performs this tutorial

Search a database of structures for a scaffold:
`python3 VaxDesign.py -s 2y7q A B 332 337 Scaffold_Database`

Run the entire protocol:
`python3 VaxDesign.py -p 2y7q A B 332 337 scaffold.pdb fixbb siwa2` where scaffold.pdb is a cleaned structure of the PDB ID 3HZ7

*Breakdown of the protocol:*

Isolate the motif:
`python3 VaxDesign.py -m 2y7q B 332 337`

Isolate the receptor:
`python3 VaxDesign.py -r 2y7q A`

Graft the motif onto the scaffold using the receptor as reference:
`python3 VaxDesign.py -g receptor.pdb motif.pdb scaffold.pdb`

Fold From Loop protocol (Not Implemented Yet):
`python3 VaxDesign.py -f motif.pdb grafted.pdb 61 66 siwa2`

Sequence design the grafted structure:
`python3 VaxDesign.py -d fixbb grafted.pdb 61 66`

Generate fragments for the grafted structure (for an Abinitio simulation):
`python3 VaxDesign.py -F structure.pdb siwa2`

Sequence design of the surface only:
`python3 VaxDesign.py -d surface grafted.pdb 1 3 5 7 8 9 10 11 12 13 16 17 19 20 23 24 25 27 28 37 38 39 40 43 47 50 51 52 53 54 55 57 58 59 60 73 74`

## Reference
Please reference the following when using this script.
* Sari Sabban (2020) Computationally Grafting An IgE Epitope Onto A Scaffold. bioRxiv 2020.09.14.296392; doi: https://doi.org/10.1101/2020.09.14.296392 
