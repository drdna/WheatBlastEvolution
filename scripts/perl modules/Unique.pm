package Unique;

my %alignHash = undef;

sub isit {

  my $Query = $_[0];

  my $inFile = "ALIGN_STRINGS_MASK/B71V3.$Query"."_alignments";

  open(ALIGNMENT, "$inFile"); #|| die "Can't open align file; $inFile\n";

  while(my $L = <ALIGNMENT>) {

    chomp($L);

    next if $L eq '';

    my($Chr, $alignString) = split(/\t/, $L);

    $alignHash{$Chr} = $alignString;

  }
  return \%alignHash;

  close ALIGNMENT

}

1; 
