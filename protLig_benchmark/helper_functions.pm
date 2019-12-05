sub get_ligands
{
    $path = $_[0];
    opendir(FOO,$path);
    my @cont = readdir(FOO);
    closedir(FOO);
    my @ligs = ();
    foreach my $dir(@cont)
    {
        if( $dir =~ /^lig_.*/ )
        { 
        	push(@ligs,$dir);
        }
    }
    return(\@ligs);
}

sub copy_gaff
{
	my $ligsRef = $_[0];
	my $dest = $_[1];
	my $source = $_[2];

	my @ligs = @$ligsRef;

	foreach my $lig(@ligs)
	{
		system("mkdir $dest/$lig");
		system("cp $source/$lig/MOL.itp $dest/$lig/.");
		system("cp $source/$lig/ffMOL.itp $dest/$lig/.");
		system("cp $source/$lig/mol_gmx.pdb $dest/$lig/.");
	}
}

sub copy_cgenff
{
	my $ligsRef = $_[0];
	my $dest = $_[1];
	my $source = $_[2];

	my @ligs = @$ligsRef;

	foreach my $lig(@ligs)
	{
		system("mkdir $dest/$lig");
		system("cp $source/$lig/MOL.itp $dest/$lig/.");
		system("cp $source/$lig/ffMOL.itp $dest/$lig/.");
		system("cp $source/$lig/mol_gmx.pdb $dest/$lig/.");
	}
}

sub copy_transformations
{
	my $dest = $_[0];
	my $source = $_[1];

	# read edges
        opendir(FOO,$source);
        my @cont = readdir(FOO);
        closedir(FOO);
        my @edges = ();
        foreach my $dir(@cont)
        {
            if( $dir =~ /^edge_.*/ )
            {
                push(@edges,$dir);
            }
        }

	# create folders and copy
	for my $edge(@edges)
	{
		system("mkdir $dest/$edge");
		system("cp $source/$edge/pairs.dat $dest/$edge/.");
		system("cp $source/$edge/mergedA.pdb $dest/$edge/.");
		system("cp $source/$edge/mergedB.pdb $dest/$edge/.");
		system("cp $source/$edge/merged.itp $dest/$edge/.");
		system("cp $source/$edge/ffMOL.itp $dest/$edge/.");
	}
}

1;
