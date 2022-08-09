#!/usr/bin/perl


open(LOOKUP, "Root-to-tip_sample_lineages.txt") || die "Usage: perl Root-to-tip_lineages.pl\nPlace Root-to-tip_sample_lineages.txt file in working directory\n";

while($D = <LOOKUP>) {

  chomp($D);

  @List = split(/\t| +/, $D);

  $Hash{$List[0]} = $List[1];

}

close LOOKUP;


# open outfile

open(OUT, '>', "Root-to-tip_data_lineages.txt") || die "Can't create output file\n";


# open undated dataset and add a lineages column 

open(DATA, "Root-to-tip_data.txt") || die "Usage: perl Root-to-tip_lineages.pl\nPlace Root-to-tip_data.txt file in working directory\n";

while($U = <DATA>) {

  chomp($U);

  $count ++;

  if($count > 1) {

    @List = split(/\t| +/, $U);

    print OUT "$Hash{$List[0]}\t$U\n";

  }
  
  else {
  
    print OUT "$U\n"
    
  }

}

close DATA; 

close OUT;
