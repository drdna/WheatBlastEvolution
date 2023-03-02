#!/usr/bin/perl

# Calculates whole genome proportions for each donor contribution to swarm

die "Usage: perl Whole_genome_proportions.pl <CP-dataframe>\n" if @ARGV != 1;

& GENOTYPES;

sub GENOTYPES {

  open(GT, "$ARGV[0]") || die "Can't find CP dataframe file\n";

  #write entries to list
  while($T = <GT>) {
    chomp($T);
    $T =~ s/"//g;
    $prevPos = $pos;
    ($pop, $pos, $prob, $width, $chr, $hapl, $yr, $haplyr, $haplyr2) = split(/,/, $T);
    $pos =~ s/\.\d+//;
    $chr =~ s/omosome //;
    $chr =~ s/C/c/;
    $genotypeHash{$pop}{$chr}{$pos} = 1;

    $AllHash{$chr}{$pos} = 1;

    $PoLHash{$chr}{$pos} = 1 if $hapl =~ /PoL/;

    $PoTHash{$chr}{$pos} = 1 if $hapl =~ /PoT/;

  }
  close GT
}


foreach $pop (keys %genotypeHash) {

  print "$pop\n";

  foreach $chr (keys %{$genotypeHash{$pop}}) {

    $numSites = keys %{$genotypeHash{$pop}{$chr}};

    print "$chr\t$numSites\n"

  }

}

foreach $chr (keys %posHash) {

  $numSites += keys %{$posHash{$chr}};

}

print "numSites = $numSites\n";


