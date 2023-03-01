#!/usr/bin/perl

# Code to generate a table of reciprocal crossover from the ChromoPainter dataframe used by plotting program 

die "Usage: perl Reciprocal.v3.pl <ChromoPainter-dataframe>\n" if @ARGV < 1;

# load List module

use List::MoreUtils qw(uniq);

# open chromopainter dataframe

open(F, $ARGV[0]);

while($L = <F>) {

  chomp($L);

  $L =~ s/"//g;

  $prevlineage = $lineage;

  $prevpos = $pos;

  $prevpos =~ s/\..+//;
  
  ($lineage, $pos, $prob, $pref, $chr, $hapl, $yr, $hapl2, $yr2)    = split(/,/, $L);

  $chr =~ s/[^0-9]+//g;

  $pos =~ s/\..+//;

  # look at predicted parentage of current SNP
  # hash the relevant details if there has been a crossover 
  # between it and the previous SNP

  if($prevlineage ne $lineage) {

    push @{$ReciprocalHash{$chr}{$prevpos}{$pos}{$prevlineage}{$lineage}}, $hapl;

  }

}

close F;


# if there has been a crossover between current SNP and previous one,
# cross-reference against the crossovers hash to see if any other haplotype has reciprocal parentage for the SNPs
# flanking the crossover


# print Table header lines

print "Chromosome\tPre-X-over\tPre-X-over\tPost-X-over\tPost-X-over\tHaplotype(s)\n";

print "\tSNP pos.\tallele\tSNP pos.\tallele\n";


# print reciprocal crossover details

foreach $chr (sort {$a cmp $b} keys %ReciprocalHash) {

  foreach $prevpos (sort {$a <=> $b} keys %{$ReciprocalHash{$chr}}) {

    foreach $pos (sort {$a <=> $b} keys %{$ReciprocalHash{$chr}{$prevpos}}) {

      foreach $prevlineage (sort {$a cmp $b} keys %{$ReciprocalHash{$chr}{$prevpos}{$pos}}) {

        foreach $lineage (sort {$a <=> $b} keys %{$ReciprocalHash{$chr}{$prevpos}{$pos}{$prevlineage}}) {

          next if exists($Analyzed{$chr}{$prevpos}{$pos}{$lineage}{$prevlineage});

          $Analyzed{$chr}{$prevpos}{$pos}{$prevlineage}{$lineage} = 1;

          if (exists($ReciprocalHash{$chr}{$prevpos}{$pos}{$lineage}{$prevlineage})) {

            @List1 = sort {$a cmp $b} uniq @{$ReciprocalHash{$chr}{$prevpos}{$pos}{$prevlineage}{$lineage}};
 
            @List2 = sort {$a cmp $b} uniq @{$ReciprocalHash{$chr}{$prevpos}{$pos}{$lineage}{$prevlineage}};

            print join ("\t", ($chr, $prevpos, $prevlineage, $pos, $lineage)), "\t", join(", ", @List1), "\n";

            print join ("\t", ('', $prevpos, $lineage, $pos, $prevlineage)), "\t", join(", ", @List2), "\n";

          }

        }

      }

    }

  }

}
