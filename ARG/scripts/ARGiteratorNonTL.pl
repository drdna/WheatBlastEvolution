#!/usr/bin/perl

## parse .arg files and determine distributions for
## time to most recent recombination events (TMRRE)
## for each PoT/PoL strain

##


# read file to find node IDs for each taxon/parent node

& READ_TAXON_NODES;


# iterate through arg file n times

for($i = 1; $i <= 10; $i++) {  
  & NODE_ITERATION;
  & READ_HASH;
}

& LINEAGES;

& READ_HASH;

& CHECK_TMRRE;

& TMRRE;

close F;


# read nodeHash after all iterations are complete

sub READ_HASH {

  # iterate through each non-PoT/PoL strain

  foreach $taxon1 (sort {$a cmp $b} keys %nodeHash) {
    next if $lineageHash{$taxon1} eq 'TL';

    # iterate through each non-PoT/PoL comparator

    foreach $taxon2 (sort {$a cmp $b} keys %nodeHash) {
      next if $taxon1 eq $taxon2 || $lineageHash{$taxon2} ne 'TL';

      # check nodes for each non-T/L member
      # if a given node is a recomb node and in path of non-T/L taxon
      # add node to array containing recombinant nodes common to the two taxa - ** keyed by ages of nodes **

      foreach $node (sort {$a <=> $b} keys %{$nodeHash{$taxon1}}) {
        if(exists($nodeHash{$taxon2}{$node}{recomb}) && $nodeHash{$taxon2}{$node}{recomb} == $nodeHash{$taxon1}{$node}{recomb}) {
          my $tmrre = $nodeHash{$taxon2}{$node}{recomb};

          # create a hash of taxa converging at recombination node and dates of convergence 
          push @{$tmrre{$taxon1}{$tmrre}{$taxon2}}, $node;
        }
      }
    }
  }
}

# code to check what's stored in TMRRE hash

sub CHECK_TMRRE {

  foreach $taxon (keys %tmrre) {
#    print "TAXON: $taxon\n"; 
    foreach $tmrre (sort {$a <=> $b} keys %{$tmrre{$taxon}}) {
      @compareTaxa = sort {$a cmp $b} keys %{$tmrre{$taxon}{$tmrre}};
#      print "@compareTaxa\n";
      foreach $taxon2 (@compareTaxa) {
#        print "$tmrre\t$taxon2: @{$tmrre{$taxon}{$tmrre}{$taxon2}}\n"
      }
    }
  }
}


sub TMRRE {

  # $tmrre{taxon1}{recomb node ages}{taxon2}
  # foreach taxon1, find youngest recombination node that intersects with a non-PoT/PoL taxon
  # report it

  foreach $taxon1 (sort {$a cmp $b} keys %tmrre) {
#    print "TAXON1: $taxon1\n";
    $success = 'no';
    @tmrres = (sort {$a <=> $b} keys %{$tmrre{$taxon1}});
#    print "@tmrres\n";

    foreach $tmrre (@tmrres) {
#      print "TMRRE: $tmrre\n";
      last if $success eq 'yes';
      if($tmrre eq '') {
        next
      }
      else {
        foreach $taxon2 (sort {$a <=> $b} keys %{$tmrre{$taxon1}{$tmrre}}) {
          if($lineageHash{$taxon2} eq 'TL') { 
#             print "$taxon2\n";
             push @taxons, $taxon2;
             push @nodes, @{$tmrre{$taxon1}{$tmrre}{$taxon2}}
          }
        }
      }        

#      print "Strain:$taxon1\nNearest PoT/PoL donor: @taxons\nTMRRE: $tmrre\nList of nodes with same age: @nodes\n\n";
      print join ("\t", ($taxon1, $tmrre, @taxons)), "\n";


      @taxons = ();
      @nodes = ();
      $success = 'yes'

    }
  }
}


# read file to find node IDs for each taxon/parent node

sub READ_TAXON_NODES {

  open(F, "$ARGV[0]");
  while($F = <F>) {
    chomp($F);
    ($node, $type, $age, $position, $parent, $child) = split(/\t/, $F);
    $age =~ s/^ //;
    next if $node =~ /name|start/;
    ($parent1, $parent2) = split(/,/, $parent) if $parent =~ /,/;
    if($node =~ /[A-Za-z]/ ) {
      if($parent =~ /,/) {
        foreach my $parent (@parents) {
          $nodeHash{$node}{$parent}{$type} = $age
        }
      }
      else {
        $nodeHash{$node}{$parent}{$type} = $age
      }
    }
  }
}


sub NODE_ITERATION {

  $newNode = 'yes';
  open(F, "$ARGV[0]");
  while($F = <F>) {
    chomp($F);
    next if $F =~ /name|start/;
    ($node, $type, $age, $position, $parent, $child) = split(/\t/, $F);
    $type = 'coal' if $type eq 'gene';
    next if $parent !~ /^\d/;
    $age =~ s/ //;
    ($parent1, $parent2) = split(/,/, $parent) if $parent =~ /,/;
    foreach $taxon (sort {$a cmp $b} keys %nodeHash) {
      if($parent =~ /,/) {
        if(exists($nodeHash{$taxon}{$node})) {
          $nodeHash{$taxon}{$parent1}{$type} = $age;
          $nodeHash{$taxon}{$parent2}{$type} = $age;
          $newNode = 'no'				# if match, change new node flag to "no"
        }
      }
      elsif(exists($nodeHash{$taxon}{$node})) {
        $nodeHash{$taxon}{$parent}{$type} = $age;
        $newNode = 'no'                                 # if match, change new node flag to "no"
      }
    }
  }
}


# read in lineage information for each strain

sub LINEAGES {
  open(L, "ARGstrain_lineages.txt");
  while($L = <L>) {
    chomp($L);
    ($strain, $index, $lineage) = split(/\t/, $L);
    $lineageHash{$strain} = $lineage;
  }
  @strainList = keys %lineageHash;
  $numStrains = @strainList;
  close L;
}
