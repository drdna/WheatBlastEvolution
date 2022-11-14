#!/usr/bin/perl

# Use continuous x scale and manually compute (center) position & width of each bar

# Open CP file

$Count = 0;

## open .cp file

open(FILE, $ARGV[0]);


## specify outfile names

$outfile1 = $ARGV[0];

$outfile2 = $ARGV[0];

$outfile1 =~ s/out/out\.forR/;

$outfile2 =~ s/out/out\.width/;


## read in entire file

$EOF= $/;

$/ = undef;

$F = <FILE>;

@HAPS = split(/\nHAP /, $F);

$Header = shift @HAPS;

$Header =~ s/ /\t/g;


## Loop through each haplotype block

foreach $haploType (@HAPS) {


  # Grab the haplotype ID

  @copyProbs = split(/\n/, $haploType);

  $hapInfo = shift @copyProbs;

  ($hapNum, $hapID) = split(/ /, $hapInfo);

  $outfile1 =~ s/.+\.chr/$hapID\.chr/;

  $outfile2 =~ s/.+\.chr/$hapID\.chr/;

  print "$outfile1, $outfile2\n";

  open(OUT1, '>', $outfile1);

  open(OUT2, '>', $outfile2);

  print OUT1 "$Header\n";

  @snpPositions = ();

  %snpPosHash = undef;

  foreach $copyProb (@copyProbs) {

    @haploType = split(/ /, $copyProb);

    $snpPosition = $haploType[0];

    $rndSnpPosition = sprintf("%.2f", $snpPosition/1000000);

#    print "$rndSnpPosition\n";

    next if exists($snpPosHash{$rndSnpPosition});

    push @snpPositions, $rndSnpPosition;

    $haploType[0] = $rndSnpPosition;

    $snpPosHash{$rndSnpPosition} = 1;
 
    print OUT1 join ("\t", @haploType), "\n";
  
  }

  close OUT1;

}

close FILE

