#!/usr/bin/perl

open(DATED, "Root-to-tip.txt");

while($D = <DATED>) {

  chomp($D);

  @List = split(/\t/, $D);

  $Hash{$List[1]} = $List[0];

}

close DATED;

open(UNDATED, "Root-to-tip2_data.txt");

while($U = <UNDATED>) {

  chomp($U);

  $count ++;

  if($count > 1) {

    @List = split(/\t/, $U);

    print "$Hash{$List[0]}\t$U\n";

  }

}

close UNDATED; 
