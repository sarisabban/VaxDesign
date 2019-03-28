# VexDesign
Scripts that computationally designs a vaccine. Authored by Sari Sabban on 31-May-2017 (sari.sabban@gmail.com).

## Requirements
1. Download this script:

`wget https://raw.githubusercontent.com/sarisabban/VexDesign/master/VexDesign.py`

2. Use the following commands (in Ubuntu GNU/Linux) to install all nessesary programs and Python libraries for this script to run successfully:

`sudo apt update && sudo apt install pymol gnuplot dssp python3-bs4 python3-biopython python3-lxml -y`

3. Make sure you install [PyRosetta](http://www.pyrosetta.org) as the website describes.

## Description
This script autonomously designs a vaccine from a user specified target site. This is not artificial intellegance, you cannot just ask these scripts to design "A" vaccine, you must understand what target site you want to develop antibodies against (make a liturature search and understand your disease and target site), then supply this target site to this script to build a protein structure around it so the final protein can be used as a vaccine. You must have prior understanding of Bioinformatics and Immunology in order to be able to understand what site to target and to supply to these scripts. Once you identify a target site, this script will take it and run a graft and desing protocol, without the need for you to intervene, that will result in an ideal protein structure displaying your target site in its original 3D cofiguration. Thus the protien, theoretically, can be used as a vaccine against this site, and hopefully neutralise the disease you are researching. Everytime you run the VexDesign.py script a different final protien sequence will be generated for the **same structure**, this is important to keep in mind, because if you want to generate many different structures with different sequences to target the same site you can simply run the script with the same parameters again and you will end up with the same structure but with a different sequence.

This script has been last tested to work well with PyRosetta 4 Release 182 using Python 3.6. If you use this script on a newer PyRosetta or Python version and it fails please notify me of the error and I will update it.

Here is a [video]() that explains how to select a target site, how the script functions, and what results you should get.

The script does the following:

*Scaffold Search*

1. Chooses ideal scaffold structures from a scaffold database (provided).

*Vaccine Design*

1. Isolates the motif.

2. Isolates the receptor.

3. Grafts the motif onto the scaffold.

4. Performs the Fold From Loop protocol to optimise the grafted structure.

5. Sequence designs the structure around the motif.

6. Generate fragments for Rosetta *Abinitio* folding simulation and analyses their quality.

Output files are as follows:

|    | File Name               | Description                                                                                  |
|----|-------------------------|----------------------------------------------------------------------------------------------|
| 1  | *PDBID_CHAIN*.pdb       | Scaffold structure                                                                           |
| 2  | FragmentAverageRMSD.dat | Average RMSD of the fragments                                                                |
| 3  | frags.200.3mers         | 3-mer fragment of sequence designed structure from the Robetta server                        |
| 4  | frags.200.9mers         | 9-mer fragment of sequence designed structure from the Robetta server                        |
| 5  | grafted.pdb             | Grafted motif to scaffold structure                                                          |
| 6  | motif.pdb	   		   | Original requested motif                                                                     |
| 7  | pre.psipred.ss2         | PSIPRED secondary structure prediction of sequence designed structure from the Robetta server|
| 8  | plot_frag.pdf           | Plot of the fragment quality RMSD vs Position                                                |
| 9  | receptor.pdb            | Original receptor that binds motif                                                           |
| 10 | RMSDvsPosition.dat      | plot_frag.pdf's data                                                                         |
| 11 | structure.fasta         | Fasta of Rosetta Designed structure                                                          |
| 12 | structure.pdb           | Sequence designed structure                                                                  |

## How To Use
### Automatic
1. Go to [Robetta server](http://robetta.org/) and register a username, it is very easy and straightforward.

2. Search potential protein scaffolds using the following command:

`python3 VexDesign.py -s PDBID RCHAIN CHAIN FROM TO DATABASE`

* PDBID = The protein's [Protein Data Bank](https://www.rcsb.org) identification name
* RCHAIN = The chain where the receptor resides within the protein .pdb file
* CHAIN = The chain where the target site resides (not part of the receptor) within the protein .pdb file
* FROM = The start of the target site
* TO = The end of the target site
* DATABASE = The directory that contains the database of scaffold structures to search through

Example:

`python3 VexDesign.py -s 2y7q A B 420 429 Database`

3. The script will generate a directory called **Scaffolds** which will contain all the scaffold structures. Search through the directory and choose one scaffold (move it to the working directory next to the VexDesign.py script).

4. Generate a vaccine using the following command:

`python3 VexDesign.py -p PDBID RCHAIN CHAIN FROM TO SCAFFOLD CHOICE USERNAME`

* PDBID = The protein's [Protein Data Bank](https://www.rcsb.org) identification name
* RCHAIN = The chain where the receptor resides within the protein .pdb file
* CHAIN = The chain where the target site resides (not part of the receptor) within the protein .pdb file
* FROM = The start of the target site
* TO = The end of the target site
* SCAFFOLD = The protein scaffold structure generated from the Scaffold.py script
* CHOICE = Choose between fixbb ot flxbb RosettaDesign
* USERNAME = The [Robetta server](http://robetta.org/) username to generate and download the *Abinitio* fragment files and the PSIPRED secondary structure prediction file

Example:

`python3 VexDesign.py -p 2y7q A B 420 429 scaffold.pdb fixbb siwa2`

Calculation time is about 72 hours on a normal desktop computer. Access to the internet is a requirement since the VexDesign.py script will be sending and retrieving data from some servers.

5. Use this [Rosetta Abinitio](https://github.com/sarisabban/RosettaAbinitio) script to simulate the folding of the final designed vaccine's protein structure. An HPC (High Preformance Computer) and the original C++ [Rosetta](https://www.rosettacommons.org/) are required for this step. This is just an evaluation step to check the design before you synthesise it, other folding simulation algorithms can be used, use what you find most reliable.

### Optimised
You can run each step separatly, use the following commands to run each step:

* Isolate motif:

`python3 VexDesign.py -m PDBID RCHAIN FROM TO` Example `python3 VexDesign.py -m 2y7q B 420 429`

* Isolate receptor:

`python3 VexDesign.py -r PDBID RCHAIN` Example `python3 VexDesign.py -r 2y7q A`

* Graft motif onto the scaffold structure

`python3 VexDesign.py -g receptor.pdb motif.pdb scaffold.pdb`

* Run the Fold From Loop protocol

`python3 VexDesign.py -f motif.pdb grafted.pdb FROM TO USERNAME` Example `python3 VexDesign.py -f motif.pdb grafted.pdb 8 17 acresearch`

In this case the FROM and TO are the positions of the motif on the grafted structure and not the original protien structure

* Sequence design the structure

`python3 VexDesign.py -d grafted.pdb FROM TO` Example `python3 VexDesign.py -d grafted.pdb 8 17`

In this case the FROM and TO are the positions of the motif on the grafted structure and not the original protien structure

* Generate fragments

`python3 VexDesign.py -F structure.pdb USERNAME` Example `python3 VexDesign.py -F structure.pdb acresearch`

* Help menu

`python3 VexDesign.py -h`

It is advised that each step is run multiple times (1000 times) and choose the lowest scoring structure to move on to the next step, thus these options are provided to run each step separatly (better on a high performace computer for speed).

### Simple Design
If you have a strucutre that you want to target, but this structure is between 80 and 150 amino acids and is a ridgid structure (mainly helices and sheets and very little loops), you might get away with just re-designing the structure itself, while keeping the target site(s) conserved, without having to do anything complicated like isolate the motif and graft it onto a scaffold protein (as discussed previously). The code for this protocol is very simple just use the following command:

`python3 VexDesign.py -S PDBID CHAIN NO_DESIGN_POSITIONS`

* PDBID = The protein's [Protein Data Bank](https://www.rcsb.org) identification name
* CHAIN = The chain of your target structure that you want to design
* NO_DESIGN_POSITIONS = The amino acid positions that you want to conserve (not re-design), i.e: your motif(s)

Example:

`python3 VexDesign.py -S 3WD5 A 19 20 21 22 23 24 65 66 67 68 140 141 143 144 145 146 147 148`

The script's protocol is as follows:
1. Gets the user defined protein structure and isolates the defined chain.
2. Sequence designs the structure (but not the user defined conserved amino acids). If the NO_DESIGN_POSITIONS is left empty the entire protein will be designed (good as a negative control).

The script will not access the Robetta server nor setup any of the *Abinitio* folding simulation files, it is meant to be used with more manual work from the user than automated work. Calculation time is about 24 hours.

Output files are as follows:

|    | File Name               | Description                                                                                  |
|----|-------------------------|----------------------------------------------------------------------------------------------|
| 1  | original.pdb            | The original non-designed structure (for comparisons)                                        |
| 2  | structure.pdb	       | The designed structure                                                                       |

## References
Please reference the following when using this script.
* 

## Thigns to do:
* Add a video
* Add scaffold database
* Add the references
