
## Summary ##
This is a collection of the input files for 13 protein-ligand systems (482 ligand perturbations):
PDE2, CMET, Bace (Hunt et al), Bace p2, Galectin, JNK1, Tyk2, Bace, MCL1, CDK2, Thrombin, PTP1b, P38

## Citation: ##
  Large Scale Relative Protein Ligand Binding Affinities Using Non-Equilibrium Alchemy. 
  Gapsys, Pérez-Benito, Aldeghi, Seeliger, Van Vlijmen, Tresadern and de Groot. 
  2020. Chemical Science. 10.1039/C9SC03754C

## Folder structure: ##
	- ligands_gaff2: for every ligand two topology files are present. 
	ffMOL.itp contains the atomtypes, the rest of topology parameters are in MOL.itp. 
	Structure is in the mol_gmx.pdb file. Force field: Gaff v2.1
	
	- ligands_cgenff: topology and structure for cgenff. 
	CGenFF v3.0.1 with the atom typing based on MATCH was used for most of the ligands, 
	except for Bace inhibitors, where paramchem.org web-server was used in combination with CGenFF v4.1.
	
	- protein_amber: structure and topology for protein and, 
	if available, co-crystallized waters and ions. Force field: amber99sb*ILDN
	
	- protein_charmm: structure and topology for Charmm36m force field
	
	- transformations_gaff2: edge information and hybrid structures/topologies for gaff.
	
	- transformations_cgenff: edge information for cgenff. 
	

