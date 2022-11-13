#!/usr/bin/perl

# Use continuous x scale and manually compute (center) position & width of each bar

# Open CP file

$Count = 0;

open(FILE, $ARGV[0]);

$outfile = $ARGV[0];

$outfile1 = $outfile;

$outfile2 = $outfile;

$outfile1 =~ s/out/out\.forR/;

$outfile2 =~ s/out/out\.width/;

$EOF= $/;

$/ = undef;

$F = <FILE>;

@HAPS = split(/\nHAP /, $F);

$Header = shift @HAPS;

$Header =~ s/ /\t/g;

# Loop through each haplotype block

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

  foreach $copyProb (@copyProbs) {

    @haploType = split(/ /, $copyProb);

    push @snpPositions, $haploType[0];
 
    print OUT1 join ("\t", @haploType), "\n";
  
  }

  close OUT1;

  & WIDTHS

}

sub WIDTHS {

  @snpPositions = reverse(@snpPositions);

  for($i = 0; $i <= @snpPositions - 1; $i++) {

    #determine width and center of last bar in series

    # Note: center is 3/4 distance to x from x-1; width = (x - x-1)/2

    if($i == @snpPositions-1) {

      $Width = ($snpPositions[$i] - $snpPositions[$i-1])/2;

      $newCenter = $snpPositions[$i-1] + 3*($snpPositions[$i] - $snpPositions[$i])/4;

      push @Widths, $Width;

      push @Centers, $newCenter

    }

    # determine width and center of first bar in series

    # note: start of window is at position of first x value; center is 1/4 distance to x+1

    elsif($i == 0) {

      $Width = ($snpPositions[$i+1] - $snpPositions[$i])/2;

      $newCenter = $snpPositions[$i] + (($snpPositions[$i+1] - $snpPositions[$i])/4);

      push @Widths, $Width;

      push @Centers, $newCenter

    }

    # calculate widths and new bar centers for "internal" datapoints

     else {

      $midpointPost = ($snpPositions[$i]+$snpPositions[$i+1])/2;

      $midpointPre = ($snpPositions[$i-1]+$snpPositions[$i])/2;

      $Width = $midpointPost - $midpointPre;

      $newCenter = $midpointPre + ($Width/2);

      push @Widths, $Width;

      push @Centers, $newCenter

    }

  }

  @Widths = reverse(@Widths);

  unshift @Widths, 'width';
 
  print OUT2 join(" ", @Widths), "\n";

  @Widths = ();

  @Centers = reverse(@Centers);

  unshift @Centers, 'center';
 
  print OUT2 join(" ", @Centers), "\n";

  @Centers = ();  

  close OUT2

}

close FILE
