#!/usr/bin/perl

##############################
#
# Generate_ARGsites.pl
#
# written by Mark L. Farman
#
# Purpose: Read SNPcaller outfiles, score all SNP positions, and record haplotypes for all individuals in a specified list
#
# Note: Must be re-run every time a new strain is added (unavoidable)
#
##############################

#use strict;

#use warnings;

use MasterList;

use FetchGenome;


die "Usage: perl Generate_ARGsites.pl <SAMPLES_LIST> <SNPFILE_DIR> <REF_GENOME_FASTA> <CHR#> <CHR-LENGTH> <skip>\n" if @ARGV != 5;

# SNP outfiles must be named according to the format: Ref_v_Subject_out 


### declare global variables

my ($Strains, $snpsDir, $refGenome, $chrNum, $chrLength, $nth) = @ARGV;

$StrainsHashRef = MasterList::strains($Strains);

%StrainsHash = %{$StrainsHashRef};

@StrainsList = keys %StrainsHash;

& SNP_ALLELES;

& REF_GENOME_SNPs;

& BASE_REF_HASH;

#my $InvariantRef = INVARIANT();

#@Invariant = @$InvariantRef;

& READ_WRITE_SNPs;

& PRINT_SITES;


## SUBROUTINES


sub SNP_ALLELES {

  # Loop through all SNP files and record positions of all SNPs relative to B71 reference genome

  #print "Identifying variant sites...\n";

  $snpsDir =~ s/\/$//;

  opendir(OUTFILES, $snpsDir);

  @snpFilesList = readdir(OUTFILES);

  #print "Identifying all relevant SNP loci\n";

  foreach $snpFile (@snpFilesList) {

    next unless $snpFile =~ /out$/;

    my($Q, $S) = $snpFile =~ /(.+)_v_(.+)_out/;

    unless(exists($StrainsHash{$S})) {

#      print "$S\n";

      next

    }

    open(FILE, "$snpsDir/$snpFile") || die "Problem\n";

    #print "$snpFile\n";

    while($L = <FILE>) {

## strip off next line if complete dataset needed

      $chromo = 'Chr'.$chrNum;

      next if $L !~ /$chromo/;

      next if $L =~ /Chr8|Chr9/;

      chomp($L);

      next if $L =~ /repeat/;

      @Data = split(/\t/, $L);

      if(@Data == 7) {

        ($B71Ref, $Other, $B71Pos, $OtherPos, $B71Nucl, $OtherNucl, $dir) = @Data;

        next if $B71Nucl !~ /^[AGTC]$/ || $OtherNucl !~ /^[AGTC]$/;

        $B71Ref =~ s/.+?(\d)$/$1/;	# strip off everything except contig identifier (at end)

        $B71Ref = "Chr"."$B71Ref";	# add back a prefix to signal a B71 chromosome ID

        ## record the presence of variant allele at each chromosome position

        $diffAllelesHash{$B71Ref}{$B71Pos} = 1

      }

    }

    close FILE  

  }

  close OUTFILES;

}


sub REF_GENOME_SNPs {

  #print "Retrieving fully masked reference genome\n";
  
  ## Hash masked reference to allow grabbing of nucleotides at variant positions

  # NB: genome masked for all repeats and positions not aligned in any one isolate


  $maskedGenome = "B71v2sh_fully_masked.fasta";

  $maskedGenomeHashRef = FetchGenome::getAlign($maskedGenome);

  %maskedGenomeHash = %$maskedGenomeHashRef;

}

sub BASE_REF_HASH {

  #print "Creating variant positions hash...\n";

  foreach my $Chr (keys %diffAllelesHash) {

    foreach my $pos (keys %{$diffAllelesHash{$Chr}} ) {

       my $RefNucl = substr($maskedGenomeHash{$Chr}, $pos-1, 1);

#       print "$RefNucl\n";

       next if $RefNucl =~ /[Nn]/;

       foreach my $Strain (@StrainsList) {

         $VariantHash{$Chr}{$pos}{$Strain} = $RefNucl

      }

    }

  }

}  


sub INVARIANT {

  #print "Retrieving inviariant nucleotides...\n";

  foreach my $chr (sort {$a cmp $b} keys %maskedGenomeHash) {

    $maskedSeq = $maskedGenomeHash{$chr};

    while($maskedSeq =~ /([AGTCagtc])/g) {

      $invNucl = $1;

      unless(exists($diffAllelesHash{$chr}{$+[0]})) {

        push @Invariant, $invNucl if $invNucl ne 'N'

      }

    }     

  }

  return(\@Invariant)

}


sub READ_WRITE_SNPs {

  #print "Reading genotypes of selected strains...\n";


  # read SNP reports again and examine genotypes at each possible SNP position:

  foreach my $File (@snpFilesList) {

    %subjSNPsHash = undef;

    next if($File !~ /out$/);

    ($Q, $S, $outsffx) = split(/_v_|_/, $File);			        # capture genome identifiers

    next unless exists($StrainsHash{$S});

    open(SNPs, "$snpsDir/$File") || die "Can't open SNPs file\n";

    ##print "Query: $Q\tSubject: $S\n";

    while(my $L = <SNPs>) {

      chomp($L);

      my @SNPs = split(/\t/, $L);

      ($qid, $sid, $qpos, $qend, $qnucl, $snucl, $dir) = @SNPs;

      next if @SNPs != 7;

      next if $qid =~ /scaf/;

      ($ChromoNum = $qid) =~ s/.+?(\d+)$/$1/;

       next if $ChromoNum !~ /^$chrNum$/;   

      next if $ChromoNum !~ /^[1-7]$/;					# inactivate for other projects

      $ChromoNum = 'Chr'.$ChromoNum;

      next if $qnucl !~ /^[AGTC]$/ || $snucl !~ /^[AGTC]$/;

      if(exists($VariantHash{$ChromoNum}{$qpos})) {

        $VariantHash{$ChromoNum}{$qpos}{$S} = $snucl

      }

    }

    close SNPs;
    
  }

}


sub PRINT_SITES {

  print "NAMES\t";

  print join ("\t", sort {$a cmp $b} @StrainsList), "\n";

# loop through VariantHash to record sites that need genotypes assigned

  foreach my $chr (sort {$a cmp $b} keys %VariantHash) {

    print "REGION\t$chr\t1\t$chrLength\n";

    foreach my $pos (sort {$a <=> $b} keys %{$VariantHash{$chr}}) {

      $varNum++;

      $everyNth = $nth+1;      

      if($varNum == 1 || $varNum % $everyNth == 0) {

        foreach my $strain (sort {$a cmp $b} @StrainsList) {

          $Data .= $VariantHash{$chr}{$pos}{$strain};

        }

        print "$pos\t$Data\n";

        $Data = '';

      }

    }

  }

}

