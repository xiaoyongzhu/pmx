.. _script_mutate:

Mutate
------

The ``pmx mutate`` script can be used to create a hybrid structure file, which
can be used for alchemical free energy calculations as shown in the :ref:`tutorials <tutorials>`. ::

    $ pmx mutate -h
    usage: mutate.py [-h] [-f infile] [-fB infileB] [-o outfile] [-ff ff]
                 [--script script] [--resinfo]

    This script applies mutations of residues in a structure file for subsequent
    free energy calculations. It supports mutations to protein, DNA, and RNA
    molecules.

    The mutation information and dummy placements are taken from the hybrid residue
    database "mutres.mtp". The best way to use this script is to take a pdb/gro file
    that has been written with pdb2gmx with all hydrogen atoms present.

    The program can either be executed interactively or via script. The script file
    simply has to consist of "resi_number target_residue" pairs.

    The script uses an extended one-letter code for amino acids to account for
    different protonation states. Use the --resinfo flag to print the dictionary.

    optional arguments:
        -h, --help       show this help message and exit
        -f infile        Input structure file in PDB or GRO format. Default is "protein.pdb"
        -fB infileB      Input structure file of the B state in PDB or GRO format (optional).
        -o outfile       Output structure file in PDB or GRO format. Default is "mutant.pdb"
        -ff ff           Force field to use. Available choices are:
                           amber99sb-star-ildn-mut
                           charmm36m-mut.ff
                           amber99sb-star-ildn-bsc1-mut
                           amber14sb-mut.
                         Default is "amber99sb-star-ildn-mut"
        --script script  Text file with list of mutations (optional).
        --resinfo        Show the list of 3-letter -> 1-letter residues
