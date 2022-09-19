#!/usr/bin/perl

## SmartSNPs.pl: a program to filter illegal SNPs calls in VCF files

## written by Mark L. Farman


## Usage warnings

die "Usage: perl SmartSNPs.pl <VCF-file> <alignfile> <AltRefRatio> <ReadDepth> <outfile> [<cutoff>]\n" if @ARGV < 5;


## load module

use isUnique;


## Read arguments

($VCF, $alignFile, $ratio, $depth, $outfile, $cutoff) = @ARGV if @ARGV == 6;

($VCF, $alignFile, $ratio, $depth, $outfile) = @ARGV if @ARGV == 5;

$cutoff = 1 if $cutoff eq '';


## initialize variables

my $repeat = 0;


## read in alignment strings that record copy number at each site

$alignHashRef = isUnique::getAlignStrings($alignFile);

$alignHash = %$alignHashRef;


## create output filenames

($outfile = $VCF) =~ s/\.vcf/_SSfilter\.vcf/ unless $VCF =~ /SSfilter/ || $VCF =~ /log/;

($log = $VCF) =~ s/\.vcf/_failed\.log/ unless $VCF =~ /SSfilter/ || $VCF =~ /log/; 


## open output file for filtered data

open (OUT, '>', $outfile) || die "Can't open output file: $outfile\n";


## open VCF file and read sample info field

open(VCF, $VCF) || die "Can't open VCF file\n";

($VCFheader = $VCF) =~ s/.+\///;

print "$VCFheader\n";		# prints a VCF file ID header line

$dataStarted = 'no';

while($V = <VCF>) {

  chomp($V);

  if($V =~ /^\#CHROM/) {

    $dataStarted = 'yes';	# flag indicating data lines are starting

#    @VCF = split(/\t/, $V, 10);

#    @Strains = split(/\t/, $VCF[9]);

    print OUT "##source=SmartSNPsV2.pl @ARGV\n";

    print OUT "$V\n";    

    next

  }

  print OUT "$V\n" if $dataStarted eq 'no'; 		# uncomment this line if the complete VCF header is required in the output

  if($dataStarted eq 'yes') {

    $numRecords ++;			# count total variants in file

    @VCF = split(/\t/, $V, 10);

    my ($allowed, $repeat) = REPEATS(\@VCF, $alignHashRef);

    ($allowed, $heterozygote) = HETEROZYGOTES(\@VCF, $allowed) if $cutoff < 2 && $allowed eq 'yes';

    ($allowed, $lowCoverage) = LOW_COVERAGE(\@VCF, $allowed) if $allowed eq 'yes';

    $VCF[6] = 'PASS' if $allowed eq 'yes';

    $VCF[6] = 'FAIL' if $allowed eq 'no';

    #print join ("\t", @VCF), "\n";

    print OUT join ("\t", @VCF), "\n";

  }

}
   
& PRINT_SUMMARY;

& SUMMARY_STATISTICS;



## SUBROUTINES ##

# remove and count SNPs occurring in repeat regions

sub REPEATS {

  # look at SNP site in corresponding alignment string
  # reject site if # of alignments is greater than specified cutoff
  # increment the repeated sites counter

  ($Vref, $alignHashRef) = @_;

  @V = @{$Vref};

  %alignHash = %{$alignHashRef};

  ($ID, $pos, $varid, $refnucl, $altnucl, $qual, $filter, $info, $format, $genotypeInfo) = @V;

  ($IDnumeric = $ID) =~ s/.+?(\d+)/$1/;

  my $allowed = 'yes';

  $substr = substr($alignHash{$IDnumeric}, $pos-1, 1);

#  print "substr = $substr, $pos\n";

#  print "$alignHash{$IDnumeric}\n";

#  print "$V\n*$IDnumeric*, $pos, $substr\n" if $substr !~ /[12]/;

  if(substr($alignHash{$IDnumeric}, $pos-1, 1) > $cutoff || substr($alignHash{$IDnumeric}, $pos-1, 1) == 0) {

    @{$summaryHash{$ID}{$pos}} = ('FAIL', 'REPEAT');

    $allowed = 'no';

    $substr = substr($alignHash{$IDnumeric}, $pos-1, 1); 	# capture substring in case printing is needed for checking

    $repeat ++;

  }


  return $allowed, $repeat;
}


# remove and count SNPs that are in heterozygous (repeated) regions

sub HETEROZYGOTES {

  my ($Vref, $allowed) = @_;

  @V = @$Vref;

  # check for heterozygous genotype calls
  # reject site if number of alt allele calls is less than 50 times ref allele calls 
  # cutoff needs to be adjusted based on overall read depth

  ($ID, $pos, $varid, $refnucl, $altnucl, $qual, $filter, $info, $format, $genotypeInfo) = @V;

  @genotypeInfo = split(/:/, $genotypeInfo);

  ($ref, $alt) = split(/,/, $genotypeInfo[1]);

  if($genotypeInfo[0] == 1 && $ref > 0) {

    unless($alt/$ref >= $ratio) {

      @{$summaryHash{$ID}{$pos}} = ('FAIL', 'HETEROZYGOUS');

      $heterozygote ++;

      $allowed = 'no';

    }

  }

  return $allowed, $heterozygote;

}


# remove and count SNPs with low read coverage

sub LOW_COVERAGE {

  my ($Vref, $allowed) = @_;

  @V = @$Vref;

  # reject line if read depth used to infer genotype is < 30
  # adjust according to average read depth

  ($ID, $pos, $varid, $refnucl, $altnucl, $qual, $filter, $info, $format, $genotypeInfo) = @V;

  @genotypeInfo = split(/:/, $genotypeInfo);

  ($ref, $alt) = split(/,/, $genotypeInfo[1]);

  if($ref+$alt < $depth) {

    @{$summaryHash{$ID}{$pos}} = ('FAIL', 'LOW_COVERAGE');

    $lowCoverage ++;

    $allowed = 'no';

  }

  return $allowed, $lowCoverage;

}


## calculate and print summary statistics

sub SUMMARY_STATISTICS {

  $allowed = $numRecords - $repeat - $heterozygote - $lowCoverage;

  print "#NumRecords: $numRecords; Allowed: $allowed; Repeated: $repeat; Non-repeat heterozygotes: $heterozygote; Low coverage: $lowCoverage\n" if $cutoff == 1;

  print "#NumRecords: $numRecords; Allowed: $allowed; Repeated: $repeat; Low coverage: $lowCoverage\n" if $cutoff > 1;

  close VCF;

}

sub PRINT_SUMMARY {

  ## open summary output file

  open(LOG, '>', $log) || die "Can't open outfile\n";

  print LOG "$VCFheader\n";

  foreach my $SeqID (sort {$a cmp $b} keys %summaryHash) {

    foreach my $position (sort {$a <=> $b} keys %{$summaryHash{$SeqID}}) {

      print LOG "$SeqID\t$position\t", 

      join ("\t", @{$summaryHash{$SeqID}{$position}}), "\n";

    }

  }

  print LOG "#NumRecords: $numRecords; Allowed: $allowed; Repeated: $repeat; Non-repeat heterozygotes: $heterozygote; Low coverage: $lowCoverage\n" if $cutoff == 1;

  print LOG "#NumRecords: $numRecords; Allowed: $allowed; Repeated: $repeat; Low coverage: $lowCoverage\n" if $cutoff > 1;

  close LOG

}
