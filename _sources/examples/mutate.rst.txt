.. _example_mutate:

Create Mutant Structure
-----------------------

Example on how to create a structure file containing hybrid residues for
an alchemical transformation. At the moment, **pmx** supports the mutation of
protein, DNA, and RNA residues.

This is how you can mutate a :download:`protein <peptide.pdb>` residue. ::

    >>> from pmx.model import Model
    >>> from pmx.alchemy import mutate

    >>> # load the PDB file
    >>> m = Model('peptide.pdb', for_gmx=True)
    >>> # perform mutation
    >>> m2 = mutate(m=m, mut_resid=9, mut_resname='R', ff='amber99sb-star-ildn-mut')
    >>> # save mutant PDB file
    >>> m2.write('mutant.pdb')


The ``Model`` can also be modified in place: ::

    >>> mutate(m=m, mut_resid=9, mut_resname='R', ff='amber99sb-star-ildn-mut', inplace=True)
    >>> m.write('mutant.pdb')

Similarly, you can also mutate :download:`DNA <dna.pdb>` with the same :func:`pmx.alchemy.mutate` function. ::

    >>> m = Model('peptide.pdb', for_gmx=True)
    >>> m2 = mutate(m=m, mut_resid=2, mut_resname='T', ff='amber99sb-star-ildn-bsc1-mut')
    >>> m2.write('mutant_dna.pdb')

With the ``verbose`` option you can get additional info printed to screen. ::

    >>> m2 = mutate(m=m, mut_resid=2, mut_resname='T', ff='amber99sb-star-ildn-bsc1-mut', verbose=True)
    log_> Residue to mutate: 2 | DA | A
    log_> Mutation to apply: A->T
    log_> Hybrid residue name: DAT
    log_> Inserted hybrid residue DAT at position 2 (chain A)
