#!/usr/bin/perl

die "Usage: perl ConstantSites.pl <ref-genome> <variant-sites-list> <VCF-dir>\n" if @ARGV < 3;

# load module

use FetchGenome;

#use warnings;

& MERGED_ALIGN_STRING;

& REF_GENOME;

& MASK_NON_ALIGNED_SITEs;

& MASK_VARIANT_SITEs;

& VCF_FILTER_LOG;

& MASK_REF_STRINGs;

& CONSTANT_SITEs;



# Read in merged alignment string

sub MERGED_ALIGN_STRING {

  open(MAS, "Merged_alignString.txt") || die "Can't find Merged_alignString.txt file\n";

  while($MAS = <MAS>) {

    chomp($MAS);

    if($MAS =~ /Merged_Alignment_String/) {

      $stringStarted = 'yes';

      next

    }

    if ($stringStarted eq 'yes') {

      my ($chr, $alignString) = split(/\t/, $MAS);

      $AlignStringHash{$chr} = $alignString

    }

  }

  close MAS

}


# grab ref genome sequences

sub REF_GENOME {

  $RefHashRef = FetchGenome::getSeqs($ARGV[0]);

  %RefHash = %$RefHashRef;

}



sub MASK_NON_ALIGNED_SITEs {

  @CHRs = qw(chr1 chr2 chr5);		# uncomment to include chromosome 5

  foreach my $chr (@CHRs) {

    while($AlignStringHash{$chr} =~ /(0+)/g) {

      ($startPos, $length) = ($-[0], length($1));

      substr($RefHash{$chr}, $startPos, $length, 'X' x $length)

    }

  }

}


sub MASK_VARIANT_SITEs {

  open(VS, $ARGV[1]);

  while($V = <VS>) {

    chomp($V);

    ($chr, $pos, $var, $numAlt, $altList) = split(/\t/, $V);

#    next if $chr =~ /chr5/;

    substr($RefHash{$chr}, $pos-1, 1, 'X')    

  }

  close VS

}


# read VCF filtering log files

sub VCF_FILTER_LOG {

  opendir(VCF, $ARGV[2]);

  @VCFs = readdir(DIR);

  foreach $VCF (@VCFs) {

    & MASK_REF_STRING if $VCF =~ /log$/;

  }


}
 

sub MASK_REF_STRINGs {

  open(V, $VCF);  

  while($L = <V>) {

    chomp($L);

    ($chr, $pos, $fail, $reason4fail) = split(/\t/, $L);

#    next if $chr =~ /chr5/;

    substr($RefHash{$chr}, $pos-1, 1, 'X')

  }

  close V;

}


sub CONSTANT_SITEs {

$SeqSegChr1 = substr($RefHash{chr1}, 2300000, 1100000);
$SeqSegChr2 = substr($RefHash{chr2}, 500000, 1000000);
$SeqSegChr5 = substr($RefHash{chr5}, 300000, 1700000);			# uncomment to include chromosome 5

(@SeqSegs) = ($SeqSegChr1, $SeqSegChr2, $SeqSegChr5);		# uncomment to include chromosome 5

foreach $Seg (@SeqSegs) {

  $InvariantPartition .= $Seg;

  while($Seg =~ /A/g) {

    $A++

  }

  while($Seg =~ /C/g) {

    $C++

  }

  while($Seg =~ /G/g) {

    $G++

  }

  while($Seg =~ /T/g) {

    $T++

  }

}

print "Contant Sites:- A: $A; C: $C; G: $G; T: $T\n";

print "$InvariantPartition\n"

}

