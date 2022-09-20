#!/usr/bin/perl

die "Usage: perl GATKviSNPcaller.pl <samples-list> <SNP-list> (from Summarize_snps script) <iSNPcaller-output-directory>\n" if @ARGV < 3;

use Unique;

## read arguments

($samplesList, $GATKsnpsList, $iSNPsOutdir) = map {$_ =~ s/\///} @ARGV;


## read sample IDs to a hash

open(SAMPLES, $samplesList);
while($S = <SAMPLES>) {
  chomp($S);
  @Samples = split(/ +/, $S);
  foreach $sample (@Samples) {  
    $SamplesHash{$sample} = 1
  }
}
close SAMPLES;


## Read in GATK SNP calls information

open(SNPList, $GATKsnpsList);
while($L = <SNPList>) {
  if($L =~ /DATA/) {
    $Data = 'yes';
    next
  }
  chomp($L);
#  print "$L\n";
  if($Data eq 'yes') {
    ($chr, $pos, $var, $num, $strains) = split(/\t/, $L);
    @Strains = split(/ /, $strains);
    @{$SNPHash{$chr}{$pos}{$var}} = @Strains
  }
}
close SNPList;
opendir(DIR, $iSNPsOutdir) || die "Can't open SNPs directory\n";
@SNPfiles = readdir(DIR);
foreach $file (@SNPfiles) {
  if($file =~ /_(.+)_out/) {
    $strain = $1;
    next unless(exists($SamplesHash{$strain}));
    next if $strain =~ /T3-1-[23]|PY6025|Py5020|T7-3|BdJes|PY6047/;	# ignore strains with poor coverage, or duplicate datasets
#    $outfilesHash{$strain} = 1;
  }

  $alignHashRef = Unique::isit($strain);
  %alignHash = %{$alignHashRef};

  open(FILE, "$iSNPsOutdir/$file") || die "Can't open SNPs file: $file\n";

#  print "processing: $strain\n";
  while(my $L = <FILE>) {
    $rep=undef;
    next if $L !~ /chr1|chr2|chr5/; # comment out to look at calls in other chromosomes
    chomp($L);
    @List = split(/\t/, $L);
    ($q, $s, $qpos, $spos, $qn, $sn, $dir, $rep) = @List if @List == 8;
    ($q, $s, $qpos, $spos, $qn, $sn, $dir) = @List if @List == 7;
    $var = $qn.$sn;
    if(exists($SNPHash{$q}{$qpos}{$var})) {
      push @{$indelHash{$q}{$qpos}{$var}}, $strain if $rep =~ /repeat/;
      push @{$VarHash{$q}{$qpos}{$var}}, $strain;
    }
    elsif(exists($SNPHash{$q}{$qpos})) {
      push @{$indelHash{$q}{$qpos}{$var}}, $strain if $rep =~ /repeat/;;
    }
    if(substr($alignHash{$q}, $qpos-1, 1) == 0) {
      push @{$nonAligned{$q}{$qpos}{$var}}, $strain;
    }
  }
}

& PRINT;

sub PRINT {

# Change loop hash to validate different sets

foreach my $chr (sort {$a cmp $b} keys %VarHash) {
  foreach my $pos (sort {$a <=> $b} keys %{$VarHash{$chr}}) {
    foreach my $var (sort {$a cmp $b} keys %{$VarHash{$chr}{$pos}}) {
      @sortedVCFstrains = sort {$a cmp $b} @{$SNPHash{$chr}{$pos}{$var}};
      @sortediSNPstrains = sort {$a cmp $b} @{$VarHash{$chr}{$pos}{$var}};
      $numVCF = @sortedVCFstrains;
      $numiSNP = @sortediSNPstrains;	
      print "$chr\t$pos\t$var\t$numVCF\t@sortedVCFstrains\n";
      print "$chr\t$pos\t$var\t$numiSNP\t@sortediSNPstrains\t | @{$nonAligned{$chr}{$pos}{$var}}|@{$indelHash{$chr}{$pos}{$var}}\n\n";

#      print "$chr\t$pos\t$var\t${$SNPHash{$chr}{$pos}{$var}}[0]\t@sortedStrains|@{$nonAligned{$chr}{$pos}{$var}}|@{$indelHash{$chr}{$pos}{$


    }
  }
}

}

## report if a strain is missing from the iSNP outfiles directory

sub MISSING {

foreach $strain (sort {$a cmp $b} keys %SamplesHash) {

  unless(exists($outfilesHash{$strain})) {

    print "Strain $strain is missing\n";

  }

}

}
