#!/usr/bin/perl

# Generates dataset listing alleles at SNP sites that are called in ALL strains

die "Usage: Create_fasta_alignment.pl <SNP-dataset> <align-string-directory\n" if @ARGV != 2;

use FetchGenome;

# Read align strings files and gather info on gaps and repeats in each B71 x test strain alignment


# get reference sequences

$refRef = FetchGenome::getSeqs("FASTAs/B71-VCF-corr2.fasta");


# create chr1 and chr2 hashes

$chrHashRef = REF_GENOME($refRef);


# create hash of gaps/repeats list 

$maskHashRef = LIST_GAPS_REPEATS();

$maskedSeqsRef = MASK_GENOMES($chrHashRef, $maskHashRef);

SNPs($maskedSeqsRef);


## subroutines

sub REF_GENOME {

  $refRef = $_[0];

  %refGenomeHash = %{$refRef};

  @keys = keys %refGenomeHash;

  foreach $value ('chr1', 'chr2', 'chr5') {

    $chrHash{$value} = $refGenomeHash{$value};

#    print "$chrHash{$value}\n";

  }

  $refGenomeHash = undef;

  return(\%chrHash)

} 

sub LIST_GAPS_REPEATS {

  opendir(ALIGNDIR, "$ARGV[1]") || die "Align files!\n";

  @alignFiles = readdir(ALIGNDIR);

  foreach $alignFile (@alignFiles) {

    next if $alignFile !~ /alignments$/;

    if($alignFile=~ /\.(.+)_alignments/) {

      $Strain = $1;

      $strainHash{$Strain}= 1;

      print "Working on $Strain\n";

    }
  
    open(ALIGN, "$ARGV[1]/$alignFile");

    while($L = <ALIGN>) {

      chomp($L);

      if($L =~ /chr1|chr2|chr5/) {

        ($chr, $alignString) = split(/\t/, $L);

        while($alignString =~ /(0+|2+)/g) {

          push @{$maskHash{$Strain}{$chr}}, ($-[0], length($1))

        }

      }

    }

    close ALIGN

  }

  return(\%maskHash)

}

sub MASK_GENOMES {

  ($chrHashRef, $maskHashRef) = @_;

  %chrHash = %{$chrHashRef};

  %maskHash = %{$maskHashRef};

  $maskedSeqsHash{B71}{chr1} = substr($chrHash{chr1}, 0, 3400000);

  $maskedSeqsHash{B71}{chr2} = substr($chrHash{chr2}, 0, 4500000);

  $maskedSeqsHash{B71}{chr5} = substr($chrHash{chr5}, 0, 200000);

  foreach $strain (sort {$a cmp $b} keys %maskHash) {

    foreach $chr (sort {$a cmp $b} keys %{$maskHash{$strain}}) {

      print "No align file for strain: $strain\n" unless exists($maskHash{$strain});

      @gapsRepeats = @{$maskHash{$strain}{$chr}};  

      # no need to make sequence any larger than it needs to be...

      $seq = substr($chrHash{$chr}, 0, 3400000) if $chr eq 'chr1';
      $seq = substr($chrHash{$chr}, 0, 4500000) if $chr eq 'chr2';
      $seq = substr($chrHash{$chr}, 0, 2000000) if $chr eq 'chr5';

      for($i = 0; $i <= @gapsRepeats-2; $i += 2) {
     
        ($index, $length) = @gapsRepeats[$i..$i+1];

        next if ($chr eq 'chr1' && $index >= 3400000);

        next if ($chr eq 'chr2' && $index >= 4500000);

        next if ($chr eq 'chr5' && $index >= 2000000);

        substr($seq, $index, $length) =~ tr/A|C|G|T/N/

      }

      $maskedSeqsHash{$strain}{$chr} = $seq;

    }

  }  

  return(\%maskedSeqsHash)

}


  
sub SNPs {

  $maskedSeqsRef = $_[0];

  %maskedSeqsHash = %{$maskedSeqsRef};

  open(SNPs, "$ARGV[0]");

  while($L = <SNPs>) {

    chomp($L);

    if($L =~ /STRAINS/) {

      $Start = 'strains';

      next

    }

    if($L =~ /DATA/) {

      $Start = 'data';

      next

    }

    if($Start eq 'strains') {

      @StrainsList = split(/ /, $L);

    }

    else {

      & PROCESS_VARIANT_STRAINS

    }

  }

}

sub PROCESS_VARIANT_STRAINS {

    $Lines ++;

    ($chr, $pos, $refAlt, $refStrain, $altStrains) = split(/\t/, $L);

    $chr =~ s/Chr/chr/;

    $ref = substr($refAlt, 0, 1);

    $RefSeq .= $ref;

    $Alt = substr($refAlt, 1, 1);

    @altStrains = split(/ /, $altStrains);

    foreach $strain (keys %strainHash) {

      $success = 'no';

      foreach $snpStrain (@altStrains) {

        if($strain eq $snpStrain) {
          $success = 'yes';
          next
        }
      }

      if($success eq 'yes') {
          $haplotypeHash{$strain}{$chr} .= $Alt;
      }

      else {
        $haplotypeHash{$strain}{$chr} .= $ref;
      }

    }   

}


open(L, "Dating_strain_list");

while($L = <L>) {

  chomp($L);

  $Strains{$L} = 1

}

close L;


open(D, "WB_dates.txt");

while($D = <D>) {

  chomp($D);

  ($id, $date) = split(/\t/, $D);

  $Dates{$id} = $date

}

close D;


($outfile = $ARGV[0]) =~ s/txt/fasta/;

open(OUT, '>', $outfile);

print OUT ">B71_2012\n$RefSeq\n";

$NumStrains = @StrainsList;

foreach $strain (@StrainsList) {

  foreach $chr (sort {$a cmp $b} keys %{$maskedSeqsHash{$strain}}) {

    if($chr eq 'chr1' || $chr eq 'chr2' || $chr eq 'chr5') {

      $Seq .= $haplotypeHash{$strain}{$chr};

      $SeqLen = length($Seq);

    }

  }

  print OUT ">$strain"."_"."$Dates{$strain}\n$Seq\n";

  $Seq = ''

}

