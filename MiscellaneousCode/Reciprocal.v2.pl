#!/usr/bin/perl

# Reads the final chromopainter dataframe that is used for plotting chromosomes in R.

die "Usage: perl Reciprocal.v2.pl <chromopainter-dataframe>\n" if @ARGV < 1;


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
  # hash the relavnt details if there has been a crossover 
  # between it and the previous SNP

  if($prevlineage ne $lineage) {

    $crossovers{$hapl}{$chr}{$pos}{$lineage} = $lineage;

    $crossovers{$hapl}{$chr}{$prevpos}{$prevlineage} = $prevlineage;

  }

}

close F;

# re-open chromopainter dataframe

open(F, $ARGV[0]);

while($L = <F>) {

  chomp($L);

  $L =~ s/"//g;

  $prevpos = $pos;

  $prevpos =~ s/\..+//;

  $prevlineage = $lineage;

  ($lineage, $pos, $prob, $pref, $chr, $hapl, $yr, $hapl2, $yr2)    = split(/,/, $L);

  $chr =~ s/[^0-9]+//g;

  $pos =~ s/\..+//;

  # if there has been a crossover between current SNP and previous one,
  # check if any other haplotype has reciprocal parentage for the SNPs
  # flanking the crossover

  if($prevlineage ne $lineage) {

    foreach $otherHapl (sort {$a cmp $b} keys %crossovers) { 

      if($otherHapl ne $hapl && exists($crossovers{$otherHapl}{$chr}{$pos}{$prevlineage}) && exists($crossovers{$otherHapl}{$chr}{$prevpos}{$lineage})) {

        @{$ReciprocalHash{$chr}{$prevpos}{$pos}{1}}[0] = $prevlineage unless @{$ReciprocalHash{$chr}{$prevpos}{$pos}{1}} > 0;

        push @{$ReciprocalHash{$chr}{$prevpos}{$pos}{1}}, $hapl;

        @{$ReciprocalHash{$chr}{$prevpos}{$pos}{2}}[0] = $lineage unless @{$ReciprocalHash{$chr}{$prevpos}{$pos}{2}} > 0;

        push @{$ReciprocalHash{$chr}{$prevpos}{$pos}{2}}, $otherHapl;

      }

    }

  }

}

close F;


print "Chromosome\tPre-X-over\tPre-X-over\tPost-X-over\tPost-X-over\tHaplotype(s)\n";

print "\tSNP pos.\tallele\tSNP pos.\tallele\n";

foreach $chr (sort {$a cmp $b} keys %ReciprocalHash) {

  foreach $pos1 (sort {$a <=> $b} keys %{$ReciprocalHash{$chr}}) {

    foreach $pos2 (sort {$a <=> $b} keys %{$ReciprocalHash{$chr}{$pos1}}) {

      $prevlineage = splice(@{$ReciprocalHash{$chr}{$pos1}{$pos2}{1}}, 0 , 1);

      $lineage = splice(@{$ReciprocalHash{$chr}{$pos1}{$pos2}{2}}, 0 , 1); 

      @List1 = sort {$a cmp $b} (uniq @{$ReciprocalHash{$chr}{$pos1}{$pos2}{1}});

      @List2 = sort {$a cmp $b} (uniq @{$ReciprocalHash{$chr}{$pos1}{$pos2}{2}});

      print join ("\t", ($chr, $pos1, $prevlineage, $pos2, $lineage)), "\t", join(", ", @List1), "\n";

      print join ("\t", ('', $pos1, $lineage, $pos2, $prevlineage)), "\t", join(", ", @List2), "\n";

    }

  }

}
