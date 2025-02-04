#!/usr/bin/perl

## written by M. Farman
## 09/29/2020

# reads pairwise distance data from iSNPcaller's "SNPcounts" output file 

print "Usage: perl Pairwise_distances_boxplot.pl boxplot.strain.idlist AllSNPcountsJan2021.txt\n" if @ARGV != 2;

open(STRAINIDs, $ARGV[0]);

while($S = <STRAINIDs>) {

  chomp($S);

  ($ID, $Pop, $Incl) = split(/\s+/, $S);

  $popHash{$ID} = $Pop

}

close(STRAINIDs);


open(SNP_SUMMARY, $ARGV[1]);

while($L = <SNP_SUMMARY>) {

    chomp($L);

    if($L =~ /masked/) {

        @LineList = split(/\t/, $L);

        ($Q, $S) = @LineList;

	$Q =~ s/Query: //;

        $Q =~ s/_.+//;

        $S =~ s/Subject: //;

        $S =~ s/_.+//;

      }
      if($L =~ /Weighted SNPs = (\d+)/) {

        $Value = sprintf("%.0f",$1);

        $Distances{$Q}{$S} = $Value;

    }

}

close SNP_summary;

foreach $Query (keys %Distances) {

# debugger

#  unless(exists($popHash{$Query})) {

#    print "Query: $Query is missing from popHash\n";

#  }

  foreach $Subject (keys %{$Distances{$Query}}) {

    print "$Query\t$Subject\t$Distances{$Query}{$Subject}\t$popHash{$Query}\n" if $popHash{$Query} eq $popHash{$Subject} && exists($popHash{$Query})

  }

}


