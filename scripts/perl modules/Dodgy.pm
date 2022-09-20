package Dodgy;

# Creates a hash of all sites where at least one sample failed the VCF filtering steps 

# subroutine


sub SNPs {

  $VCFDIR = $_[0];

  opendir(DIR, $VCFDIR);

  @VCFs = readdir(DIR);

  foreach $VCF (@VCFs) {

    next unless $VCF =~ /SSfilter/;

    open(VCF, "$VCFDIR/$VCF") || die "Can't open VCF file $VCF\n";

    while($V = <VCF>) {

      if($V =~ /FAIL/) {

        @Dodgy = split(/\t/, $V);

        $DodgyHash{$Dodgy[0]}{$Dodgy[1]} = 1

      }

    }

  }

  return (\%DodgyHash);

  closedir DIR;

}

1;

