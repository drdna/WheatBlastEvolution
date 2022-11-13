package MasterList;

# Reads in list of strains and population info and returns a hash keyed by strainID (values = population)

sub strains {

  my $infile = $_[0];

#  die "Check input file. Expecting name to start with Samples...\n" if $infile !~ /Samples/;

  open(LIST, $infile);

  while($L = <LIST>) {

    chomp($L);

    @List = split(/\t+|\s+/, $L);

    ($StrainID, $pop) = @List[0,1];

    $StrainHash{$StrainID} = $pop

  }

  close LIST;

  return(\%StrainHash)

}

1;
