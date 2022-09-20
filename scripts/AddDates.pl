#!/usr/bin/perl

die "Usage: perl AddDates.pl <samples-list> <fasta-file>\n" if @ARGV != 2;
 
open(L, "$ARGV[0]");

while($L = <L>) {

  chomp($L);

  @List = split(/ +/, $L);

  foreach $strain (@List) {

    $Strains{$strain} = 1

  }

}

close L;


open(D, "WB_dates.txt");

while($D = <D>) {

  chomp($D);

  ($id, $date) = split(/\t/, $D);

  $Dates{$id} = $date

}

close D;

open(F, "$ARGV[1]");

while($F = <F>) {

  chomp($F);

  if($F =~ /^>(.+)/) {

    $print = 'no';

    my $id = $1;

    if(exists($Strains{$id})) {

      print ">$id"."_"."$Dates{$id}\n";

      $print = 'yes'

    }

  }

  elsif($print eq 'yes') {

    print "$F\n"

  }

}

close F;

    
