use strict;
use warnings;
use Cwd;
use helper_functions;

my $base = getcwd;


############
# mcl1 #
############
my $case = "mcl1";
system("mkdir $base/$case");
# raw gaff data
my $ligsRef = get_ligands("/home/vgapsys/project/janssen/schroedinger_set/$case/top_gaff2_sigmahole");
my @ligs = @$ligsRef;
#print "@ligs\n";

system("mkdir $base/$case/ligands_gaff2");
copy_gaff(\@ligs,"$base/$case/ligands_gaff2","/home/vgapsys/project/janssen/schroedinger_set/$case/top_gaff2_sigmahole");

system("mkdir $base/$case/ligands_cgenff");
copy_cgenff(\@ligs,"$base/$case/ligands_cgenff","/home/vgapsys/project/janssen/schroedinger_set/$case/top_cgenff_sigmahole");

my $protPath = "/home/vgapsys/project/janssen/schroedinger_set/$case";
system("mkdir $base/$case/protein_amber");
system("cp $protPath/topol_amber.top $base/$case/protein_amber/.");
system("cp $protPath/protein/amber/*itp $base/$case/protein_amber/.");
system("cp $protPath/protein/amber/protein.pdb $base/$case/protein_amber/.");
system("cp $protPath/protein/water.pdb $base/$case/protein_amber/.");
#system("cp $protPath/protein/ions.pdb $base/$case/protein_amber/.");

system("mkdir $base/$case/protein_charmm");
system("cp $protPath/topol_charmm.top $base/$case/protein_charmm/.");
system("cp $protPath/protein/charmm/*itp $base/$case/protein_charmm/.");
system("cp $protPath/protein/charmm/protein.pdb $base/$case/protein_charmm/.");
system("cp $protPath/protein/water.pdb $base/$case/protein_charmm/.");
#system("cp $protPath/protein/ions.pdb $base/$case/protein_charmm/.");

system("mkdir $base/$case/transformations_gaff2");
system("cp /home/vgapsys/project/janssen/schroedinger_set/$case/ligands/edges.txt $base/$case/transformations_gaff2/.");
copy_transformations("$base/$case/transformations_gaff2","/home/vgapsys/project/janssen/schroedinger_set/$case/sim_gaff2_sigmahole/water");

system("mkdir $base/$case/transformations_cgenff");
system("cp /home/vgapsys/project/janssen/schroedinger_set/$case/ligands/edges.txt $base/$case/transformations_cgenff/.");
copy_transformations("$base/$case/transformations_cgenff","/home/vgapsys/project/janssen/schroedinger_set/$case/sim_cgenff_sigmahole/water");

