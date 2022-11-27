#!/usr/bin/perl

die "Usage: perl Create_align_stringsv2.pl <ref_genome> <blast-dir/blast-file> <out-dir>\n" if @ARGV < 3;

use FetchGenome;

# Get lengths of RefSeqs 


$GenomeLengthHashRef = FetchGenome::getLens($ARGV[0]);

%GenomeLengthHash = %$GenomeLengthHashRef;

# check for log file

my $outDir = $ARGV[2];

unless (-d "$outDir") {

  mkdir "$outDir" || die "can't make dir: $outDir\n";

}

# Read blast directory

($blast = $ARGV[1]) =~ s/\/$//;

if( -d $blast) {

  opendir(BLASTDIR, $blast) || die "Can't open blast directory for reading\n";

  & CHECKLOG if $ARGV[3] ne 'force';

  & READ_DIR;

}


if( -f $blast) {

  ALIGN_STRINGS($blast)

}


sub CHECKLOG {

  if(-f "$outDir/alignLog") {

    open(LOGFILE, "$outDir/alignLog") || die "Can't open alignLog file\n";

    while(my $L = <LOGFILE>) {

      chomp($L);

      $logHash{$L} = 1

    }

  }

}


sub READ_DIR {

  @Files = readdir(BLASTDIR);

  foreach my $File (@Files) {

    next if $File =~ /^\./;

    next if $File !~ /BLAST$/;

    print "Processing file: $File\n";

    ($Q, $S, $suff) = split(/\./, $File);

    next if(exists($logHash{$S}));

    ALIGN_STRINGS($File, $S);

  }

}

sub ALIGN_STRINGS {

  my ($blastFile, $sStrain) = @_;

  my %alignStringHash = undef;

  foreach my $Seq (keys %GenomeLengthHash) {

    $alignStringHash{$Seq} = "0" x $GenomeLengthHash{$Seq};

  }

  if (-f $blastFile) {

    open(BLAST, "$blastFile") || die "Can't open BLAST file\n";

  }

  elsif( -d $blast) {

   open(BLAST, "$blast/$blastFile") || die "Can't open BLAST files\n";

  }

  if($blastFile =~ /([A-Za-z0-9-_]+)\.([A-Za-z0-9-_]+)\.BLAST/) {

    @outfields = ($1, $2);

  }

  $outfile = join(".", @outfields), 

  $outfile = "$outfile"."_alignments";

  print "Writing outfile: $outfile\n";

  open(OUT, '>', "$ARGV[2]/$outfile") || die "Can't open outfile: $ARGV[2]/$outfile\n";

  while(my $blast = <BLAST>) {

    ($qid, $sid, $qs, $qe, $ss, $se, $btop) = split(/\t/, $blast);

    substr($alignStringHash{$qid}, $qs-1, ($qe-$qs)+1) =~ tr/01/12/;

  }

  foreach $record (sort {$a cmp $b} keys %alignStringHash) {

    next if $record eq '';

    print OUT "$record\t$alignStringHash{$record}\n"

  }

  close OUT;

  open (LOGFILE, '>>', "$outDir/alignLog") || die "Can't open log file for writing\n";

  print LOGFILE "$outfields[1]\n";

  close LOGFILE;

  close BLAST;

}     
