# VexDesign
A script that autonomously designs a vaccine. Authored by Sari Sabban on 31-May-2017 (sari.sabban@gmail.com).

## Requirements:
1. Make sure you install [PyRosetta](http://www.pyrosetta.org) as the website describes.
2. Use the following commands (in GNU/Linux) to install all nessesary Python libraries for this script to run successfully:

`sudo apt install python3-pip pymol DSSP gnuplot`

`sudo python3 -m pip install zeep numpy biopython bs4`

## How To Use:
1. Use the following command to run the script:

`python3 VaxDesign.py PDBID CHAIN FROM TO`

* PDBID = The protein's [Protein Data Bank](https://www.rcsb.org) identification name.
* CHAIN = The chain where your target site resides within the protein .pdb file.
* FROM = the start of your target site.
* TO = the end of your target site.

2. Calculation time is about 72 hours.
3. Access to the internet is a requirement since the script will be sending and retrieving data from several servers.

## Description
This script autonomously designs a vaccine from a user specified target site. This is not artificial intellegance, you cannot just ask the the script to design "A" vaccine, you must understand what target site you want to develop antibodies against (make a liturature search and understand your disease and target site), then supply this target site to the script to build a protein structure around it so the final protein can be used as a vaccine. You must have prior understanding of Bioinformatics and Immunology in order to be able to understand what site to target and to supply it to the script. Once you identify a target site, the script will take it and run a long protocol, without the need for you to intervene, that will result in an ideal protein structure displaying your target site in its original 3D cofiguration. Thus the protien, theoretically, can be used as a vaccine against this site, and hopefully neutralise the disease you are researching. Everytime you run this script a different final protien structure will be generated, this is important to keep in mind, because if you want to generate many different structures to test or to use as boosts you can simply run the same target site again and you will end up with a different final structure.

Here is a [video](youtube.com/) that explains how to select a target site, how the script functions, and what results you get. If I did not make a video yet, bug me until I make one.

The script protocol is as follows:
* STILL UNDERT DEVELOPMENT
