#!/usr/bin/perl

die "Usage: Summarize_snps_no_dodgy.pl <VCF-dir>\n" if @ARGV < 1;

die "Program expects a DIRECTORY as argument\n" unless -d $ARGV[0];

use Dodgy;


& DODGY;

& CUT_SITES;

& ADD_BACK_STRAINS;

& READ_SITES_LIST;

& PRINT_OUTPUT;




## SUBROUTINES

# Run Dodgy module to identify sites that failed the VCF filter based on at least one criterion

sub DODGY {

  $DodgySNPsHashRef = Dodgy::SNPs($ARGV[0]);

  %DodgySNPsHash = %$DodgySNPsHashRef;

}



sub CUT_SITES {

  # open manually curated list of sites that escaped intial filters 

  open(CUT, "Chr1Chr2Chr5_disallowed_sites") || die "Can't find Chr1Chr2Chr5_disallowed_sites file\n";

  while($C = <CUT>) {

    chomp($C);

    my ($chr, $pos, $reason4cut) = split(/\t| +/, $C);

    $CutHash{$chr}{$pos} = $reason4cut;

  }

  close CUT;

}


sub ADD_BACK_STRAINS {


  # open list of strains to add to variant sites lists (usually filtered out due to poor coverage)

  open(ADD, "Chr1Chr2Chr5_add-backs") || die "Can't find Chr1Chr2Chr5_add-backs file\n";

  while($A = <ADD>) {

    chomp($A);

    my ($chr, $pos, $addons) = split(/\t| +/, $A, 3);

    @addons = split(/ /, $addons);

    @{$AddHash{$chr}{$pos}} = @addons;

  }

  close ADD;

}



sub READ_SITES_LIST {

  # read file list and create a list of strain IDs

  opendir(DIR, $ARGV[0]) || die "Can't open directory\n";

  @filesList = readdir(DIR);

  foreach $file (@filesList) {

    next if $file =~ /PY6025|PY6047|BdJes|T7-3|Py5020/;

    next unless $file =~ /SSfilter/;

    ($strainID = $file) =~ s/_.*//;

    $strainHash{$strainID} = 1;

  }

  # open SmartSNPsV2.pl-filtered VCF files and read SNP data

  foreach $file (@filesList) {

    next if $file =~ /PY6025|PY6047|BdJes|T7-3|Py5020/;

    next unless $file =~ /SSfilter/;

    open(F, "$ARGV[0]/$file")|| die "Can't open file\n";

    while($L=<F>) {

      next if $L =~ /^#/;

      next if $L =~ /FAIL/;	# skip over failed calls

      chomp($L);

      @L = split(/\t/, $L);

      if ($file =~ /^(.+)?_geno.+vcf/) {		# only process genotype files

        $strainID = $1;

        $var = $L[3].$L[4];

        push @{$Hash{$L[0]}{$L[1]}{$var}}, $strainID;

      }

    }

    close F

  }

}


sub PRINT_OUTPUT {

  @strains = (sort {$a cmp $b} keys %strainHash);

  print "STRAINS\n";

  print join (" ", @strains), "\nDATA\n";

  foreach $chr (sort {$a cmp $b} keys %Hash) {
 
    foreach $pos (sort {$a <=> $b} keys %{$Hash{$chr}}) {

      if(exists($DodgySNPsHash{$chr}{$pos}) || exists($CutHash{$chr}{$pos})) {			# skip over any sites listed in "dodgy" or "cut" hashes

      #  print "skipping $chr\$pos\n";

        next;

      }

      next if exists($AddHash{$chr}{$pos});		# skip over any additional sites listed in the "add" hash (super aggressive filtering to minimize false calls)

      foreach $var (sort {$a <=> $b} keys %{$Hash{$chr}{$pos}}) {

#        push @{$Hash{$chr}{$pos}{$var}}, @{$AddHash{$chr}{$pos}};	uncomment this line to add a strain that showed a SNP but didn't meet call threshold

        $arrayLen = @{$Hash{$chr}{$pos}{$var}};

        @array = sort {$a cmp $b} @{$Hash{$chr}{$pos}{$var}};

        print join("\t", ($chr, $pos, $var, $arrayLen)), "\t@array\n";

    }

  }

}

}
