use strict;
use warnings;
use Cwd;
use helper_functions;

my $base = getcwd;


############
# pde2 #
############
my $case = "pde2";
system("mkdir $base/$case");
# raw gaff data
my $ligsRef = get_ligands("/home/vgapsys/project/janssen/top_gaff2_am1bcc_sigmahole_new");
my @ligs = @$ligsRef;
#print "@ligs\n";

system("mkdir $base/$case/ligands_gaff2");
copy_gaff(\@ligs,"$base/$case/ligands_gaff2","/home/vgapsys/project/janssen/top_gaff2_am1bcc_sigmahole_new");

system("mkdir $base/$case/ligands_cgenff");
copy_cgenff(\@ligs,"$base/$case/ligands_cgenff","/home/vgapsys/project/janssen/top_cgenff_match_sigmahole_new");

system("mkdir $base/$case/protein_amber");
system("cp /home/vgapsys/project/janssen/sim_gaff2_am1bcc_sigmahole_new/prot/topol.top $base/$case/protein_amber/topol_amber.top");
system("cp /home/vgapsys/project/janssen/For_vytas/top_amber/*itp $base/$case/protein_amber/.");
system("cp /home/vgapsys/project/janssen/For_vytas/top_amber/prot_ions.pdb $base/$case/protein_amber/.");
system("cp /home/vgapsys/project/janssen/For_vytas/water/water.pdb $base/$case/protein_amber/.");

system("mkdir $base/$case/protein_charmm");
system("cp /home/vgapsys/project/janssen/sim_cgenff_match_sigmahole_new/prot/topol.top $base/$case/protein_charmm/topol_charmm.top");
system("cp /home/vgapsys/project/janssen/For_vytas/top_charmm/*itp $base/$case/protein_charmm/.");
system("cp /home/vgapsys/project/janssen/For_vytas/top_charmm/prot_ions.pdb $base/$case/protein_charmm/.");
system("cp /home/vgapsys/project/janssen/For_vytas/water/water.pdb $base/$case/protein_charmm/.");

system("mkdir $base/$case/transformations_gaff2");
system("cp /home/vgapsys/project/janssen/For_vytas/edges.txt $base/$case/transformations_gaff2/.");
copy_transformations("$base/$case/transformations_gaff2","/netmount/energy/vgapsys/janssen/other_sets/pde2/sim_gaff2_sigmahole_scaleMass/water");

system("mkdir $base/$case/transformations_cgenff");
system("cp /home/vgapsys/project/janssen/For_vytas/edges.txt $base/$case/transformations_cgenff/.");
copy_transformations("$base/$case/transformations_cgenff","/netmount/energy/vgapsys/janssen/other_sets/pde2/sim_cgenff_sigmahole_scaleMass/water");

