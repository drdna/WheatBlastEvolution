package FetchGenome;

sub getAlign {

  my %GenomeHash;

  my $EOL = $/;

  $/ = undef;

  open(GENOME, "$_[0]") || die "Can't open genome file\n";

  my $wholeRecord = <GENOME>;

  my @GenomeList = split(/>/, $wholeRecord);

  shift @GenomeList;

  foreach my $fastaRecord (@GenomeList) {

    my ($ID, $Seq) = split(/\n/, $fastaRecord, 2);

    $AlignHash{$ID} = $Seq 

  }

  close GENOME;

  $/ = $EOL;

  return(\%AlignHash)

}

sub getSeqs {

  my %GenomeHash;

  my $EOL = $/;

  $/ = undef;

  open(GENOME, "$_[0]") || die "Can't open genome file\n";

  my $wholeRecord = <GENOME>;

  my @GenomeList = split(/>/, $wholeRecord);

  shift @GenomeList;

  foreach my $fastaRecord (@GenomeList) {

    my ($ID, $Seq) = split(/\n/, $fastaRecord, 2);

#    $ID =~ s/.+?(\d+)$/$1/;

#    print "$ID\n";

    $ID =~ s/ +.+//;

    $Seq =~ s/[^A-Za-z-]//g;

    $Seq =~ s/\r//g;

    $GenomeHash{$ID} = $Seq

  }

  close GENOME;

  $/ = $EOL;

  return(\%GenomeHash)

}


sub getLens {

  my %GenomeHash;

  my $EOL = $/;

  $/ = undef;

  open(GENOME, "$_[0]") || die "Can't open genome file\n";

  my $wholeRecord = <GENOME>;

  my @GenomeList = split(/>/, $wholeRecord);

  shift @GenomeList;

  foreach my $fastaRecord (@GenomeList) {

    my ($ID, $Seq) = split(/\n/, $fastaRecord, 2);

#    print "$ID\t";

#    $ID =~ s/.+?(\d+$)/$1/;

    # $ID =~ s/ +.+//;

    $Seq =~ s/[^A-Za-z-]//g;

    $RefSeqLen = length($Seq);

    $GenomeLenHash{$ID} = $RefSeqLen;

  }

  close GENOME;

  $/ = $EOL;

  return(\%GenomeLenHash)

}

1;

