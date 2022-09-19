package StrictUnique4;

### This version does not assign a chromosome length hash based on B71 reference genome

### Reads BLAST reports and reports uniquely aligned nucleotide positions

### in each genome. Output format is a text file with each nt position represented as '0' (missing or repeated)

### or a '1' (uniquely aligned)

### Usage: StrictUnique.pm <inDirectory> <outDirectory> 

##############################
#
# SStrictUnique2.pm
#
# written by Mark L. Farman
#
# Purpose: Read BLAST reports, count and report SNPs in uniquely aligned genomic segments 
#
##############################

#use strict;
#use warnings;
use Cwd qw(cwd);
use File::Path;
use File::Copy;

# BLAST result files must be named according to the format: Query.Subject.BLAST
# The BLAST format used must be: -outfmt '6 qseqid sseqid qstart qend sstart send btop' 

###############################

  my @QArray;
  my @SArray;
  my $Q;
  my $S;
  my $BLAST;
  my %SNPsTempHash = ();
  my $Border = "@" x 40;
  my $Separator = "-" x 40;
  my $DateString;
  my $indir;
  my $outdir;
  my $logdir;
  my $cwd;

sub SNPs {

#  if(@ARGV != 2) { die "Usage: Run_SNPcaller.pl <indir> <outdir>\n"}
  ($indir, $outdir) = @{$_[0]};
  print "$indir, $outdir\n";
  $logdir = "$outdir/LOGS";

  if(-d "$outdir") {
    print "A directory named $outdir already exists. Are you sure you want to continue?\n";
   
    $answer = <STDIN>;

    chomp($answer);

    if($answer =~ /^Y|y/) {

      print   "\n$Border\n\n".
              "RUNNING THE SNP CALLER\n\n";

      & RUN_PROGRAM;
    }
    else {
      die "Run again using a different outfile name.\n"
    }
  }

  else {
    mkpath( "$outdir/LOGS");
    & RUN_PROGRAM
  }

  if(-f "SNPcall_tempfile") {
    close SNPTEMPFILE;
  }
}
1;

#################################
# SUBROUTINES

# 1. Run the main program

sub RUN_PROGRAM {

  ## Check for a tempfile from an interrupted previous run. Read file, hash information for back referencing

  my $BlastFile;
  my %GenomeHash;

  $cwd = cwd;

  chdir("$outdir") || die "Didn't cd\n";

  if(-f "SNPcall_tempfile") {
    open(SNPTEMPFILE, "$SNPcall_tempfile");
    my $LineCount = 0;

    while(my $L = <SNPTEMPFILE>) {
      chomp($L);
      $BlastFile = $L;
      $SNPsTempHash{$BlastFile} = 1;
    }

    close SNPTEMPFILE;
  }
  open(SNPTEMPFILE, '>>', "SNPcall_tempfile");

  unless(-d "SNP_COUNTS") {

        mkdir "SNP_COUNTS"
  }

  & OPEN_OUTFILES;					# create output directory and start writing main outfile

  my @FilesList = READ_DIR();				# read file list from input directory

  foreach my $BlastFile (@FilesList) {
    next if $BlastFile !~ /BLAST/;
    next if(exists($SNPsTempHash{$BlastFile}));

    if(($BlastFile !~ /^\./) && ($BlastFile =~ /BLAST/)) {	# capture genome identifiers 
      ($Q, $S, $BLAST) = split(/\./, $BlastFile);
      foreach $Genome ($Q, $S) {
        $Genome =~ s/_.*//;
        unless(exists($GenomeHash{$Genome})) {
          print LOGFILE "$Genome\n";
	  $GenomeHash{$Genome} = 1
        }
      }
    }


    print MAINREPORT "$Q\t$S\t";			# write info on compared genomes to outfile

    print CUMULATIVE "$Q\t$S\t";                        # write info on compared genomes to outfile

    PRE_SCREEN($BlastFile);					# prescreen BLAST report files and identify repeated DNA segments

    my ($SNPs, $MAligned, $MMAligned) = COUNT_SNPS($BlastFile);	# count SNPs in single-copy DNA segments
    
    & PRINT_SNPs;

    ($SNPs, $MAligned, $MMAligned, my $TotalAligned) = RESULTS($SNPs, $MAligned, $MMAligned, $BlastFile);	# write results summary 

  }

  & PRINT_ALIGN_STRINGS;

  close MAINREPORT;

  close CUMULATIVE;

  chdir($cwd)

}


# 2. Open and name the SNP summary and logfiles

sub OPEN_OUTFILES {

    my $DateTime = localtime();
    my @DateTime = split(/\s+/, $DateTime);
    shift(@DateTime);
    $DateString = join "_", @DateTime;
    $DateString =~ s/\:/-/g;
    my $fh1 = select MAINREPORT;
    $| = 1;
    select $fh1;
    open(MAINREPORT, '>>', "SNP_counts_out_"."$DateString") || die "Can't create snpsout file\n";
    my $fh2 = select LOGFILE;
    $| = 1;
    select $fh2;
    open(LOGFILE, '>>', "LOGS/logfile_"."$DateString") || die "Can't create logfile\n";
    print LOGFILE "SNPcaller run: $DateString\nGenomes analyzed:\n\n";	
    open(CUMULATIVE, '>>', "SNP_counts_cumulative.txt");
}

## 3. Read input directory

sub READ_DIR {

    opendir(INDIR, "$cwd/$indir") || die "Can't open input directory\n";
    my @FilesList = readdir(INDIR);
    return(@FilesList)
}


## 4. Pre-screen for repeated sequences in query and subject

sub PRE_SCREEN {

# Checks if query and subject nucleotides in current alignment are listed in the respective repeats hashes
# if yes, skips on to next nucleotide
# If no, checks alignment hashes to see if they are already in an alignment
# if yes, writes them to repeats hash and deletes alignment hash record
# if no, write to alignment hash

  my $k;
  my $j;
  my $qnum;
  my $snum;
  undef @QArray;
  undef @SArray;
  my ($File)  = @_;

  open(BLAST, "$cwd/$indir/$File") || die "Can't open BLAST file\n";
  print "$Separator\nQuery: $Q\tSubject: $S\n"."Creating hash table of single-copy regions\n";

  # the following examines each BLAST alignment and, for both query and subject, counts how often each base position occurs in an aligned segment

  while(my $L = <BLAST>) {
  chomp($L);
    my @BLAST = split(/\t/, $L);
    my($qid, $sid, $qpos, $qend, $spos, $send, $AlignSumm) = @BLAST;

#    print "$qid, $S\n";


# function to count times nucleotide alignments

    $countString{$qid}{$S} .= "0" x ($qend - length $countString{$qid}{$S}) if $qend > length($countString{$qid}{$S});

    substr($countString{$qid}{$S}, $qpos - 1, 1 + $qend - $qpos) =~ tr/01/12/;


    my @List = ($qid, $sid);
    my @List2 = map { $_ =~ s/.+?(\d+)$/$1/; $_} @List;
    ($qnum, $snum) = @List2;

    for($k = $qpos; $k <= $qend; $k++) {		# loop through each query base in current alignment                                               		
      $QArray[$qnum][$k] ++
    }

    if($spos < $send) {					# determine orientation of subject sequence
      for($j = $spos; $j <= $send; $j++) {		# loop through each subject base in current alignment
        $SArray[$snum][$j] ++
      }
    }

    else {						# repeat above actions, except count subject base positions backwards (inverted alignment)
      for($j = $send; $j <= $spos; $j++) {
        $SArray[$snum][$j] ++
      }
    }
  }
  close BLAST;
}


# 6. Read BLAST file again and count SNPs and number of uniquely aligned bases

	
sub COUNT_SNPS {

  ## assign local variables

  my $SNPs = 0;
  my $MAligned = 0;
  my $MMAligned = 0;
  my $TotalAligned = 0;
  my ($File) = @_;						# read arguments				
  open(BLAST, "$cwd/$indir/$File") || die "Can't open BLAST file\n";
  my $outfile = "$Q"."_v_"."$S"."_out";					# generate name for output file 
  open(OUT, '>', "$outfile");					# open output file for storing information on all SNPs
  print "Identifying SNPs between query and subject sequences\n";	# report information on run progress to screen


  ## start reading BLAST report  

  while(my $TheLine = <BLAST>) {
    chomp($TheLine);
    my($qid, $sid, $qpos, $qend, $spos, $send, $AlignSumm) = split(/\t/, $TheLine);
    ## run subroutine that parses the back trace operations portion of BLAST report
    ($SNPs, $MAligned, $MMAligned) = ANALYZE_ALIGNMENT($qid, $sid, $qpos, $qend, $spos, $send, $AlignSumm, $SNPs, $MAligned, $MMAligned);
  }

  close BLAST;
  return ($SNPs, $MAligned, $MMAligned)				# return results to main program

}


# 7. Parse back trace operations report


sub ANALYZE_ALIGNMENT {

  my $AlignedSegment = 0;
  my ($qid, $sid, $qpos, $qend, $spos, $send, $AlignSumm, $SNPs, $MAligned, $MMAligned) = @_;
  my @Align = split(/(\d+)/, $AlignSumm);			# convert BTOP output to list
  shift @Align;

  my @List = ($qid, $sid);
  my @List2 = map { $_ =~ s/.+?(\d+)$/$1/; $_} @List;
  my($qnum, $snum) = @List2;

  foreach my $Alignment (@Align) {				# loop through BTOP list;
    $AlignedSegment ++;

    if($Alignment =~ /\d+/) {					# if list member is an integer, analyze the alignment 
      ($qpos, $spos, $MAligned) = COUNT_ALIGNED_MATCHES($qid, $sid, $qpos, $spos, $send, $Alignment, $MAligned);

      if($AlignedSegment == 1) {
         # adjust nucleotide positions for the first SNP calls in the alignment
         $qpos --;

         if($spos < $send) {
           $spos --
         }

         else {
           $spos ++
         }
      }
    }

    if($Alignment =~ /\w+/) {					# if list members are alphabetic, analyze the mismatch
      ($qpos, $spos, $SNPs, $MMAligned) = ANALYZE_ALIGNED_MISMATCHES($qid, $sid, $qpos, $spos, $send, $Alignment, $SNPs, $MMAligned)
    }

  }
  return ($SNPs, $MAligned, $MMAligned)				# return results to COUNT_SNPS subroutine
}


# 8. Count aligned bases that do not occur in repeated DNA segments


sub COUNT_ALIGNED_MATCHES {

  my($qid, $sid, $qpos, $spos, $send, $Alignment, $MAligned) = @_;
  my @List = ($qid, $sid);
  my @List2 = map { $_ =~ s/.+?(\d+)$/$1/; $_} @List;
  my($qnum, $snum) = @List2;

  for(my $m = 0; $m <= $Alignment - 1; $m ++) {
                                                 # loop thro$
    if(($QArray[$qnum][$qpos] == 1) && ($SArray[$snum][$spos] == 1 )) {
        $MAligned ++;
    }

    $qpos ++;

    if($spos < $send) {                                                                     # keep trac$
        $spos ++
    }

    else {
        $spos --
    }

  }

  return($qpos, $spos, $MAligned)					# return results to ANALYZE_ALIGNMENT subroutine

}


# 9. Analyze, count and record mismatching bases 


sub ANALYZE_ALIGNED_MISMATCHES {

    my ($qid, $sid, $qpos, $spos, $send, $Alignment, $SNPs, $MMAligned) = @_;

    my @List = ($qid, $sid);

    my ($qnum, $snum) = map { $_ =~ s/.+?(\d+)$/$1/; $_ } @List;

    my @MutsList = split(//, $Alignment);				# Split mismatch string into a list of characters

    my $NumMutations = @MutsList;					# Count number of characters in list

    for(my $j = 0; $j <= $NumMutations-2; $j += 2) {			# Loop through mismatches, one pair of bases at a time

        my ($qnucl, $snucl) = ($MutsList[$j], $MutsList[$j+1]);		# Read query base and subject base in mismatch

        if(($qnucl =~ /[A-Za-z]/) && ($snucl =~ /[A-Za-z]/)) {  		# check that both are simple nucleotide substitutions
            $qpos ++;							# increment base position counters

            if($spos < $send) {
                $spos ++;
                $SNPsHash{$qid}{$qpos} = [($sid, $spos, $qnucl, $snucl, 'f')]
            }

            else {
                $spos --;
                $SNPsHash{$qid}{$qpos} = [($sid, $spos, $qnucl, $snucl, 'r')]
            }

            if(($QArray[$qnum][$qpos] < 2) && ($SArray[$snum][$spos] < 2) && ($qnucl =~ /[AGTCagtc]/) && ($snucl =~ /[AGTCagtc]/)) {	# make sure base positions are in uniquely aligned regions
                $SNPs ++;						# increment SNPs counter
                $MMAligned ++;						# increment aligned bases counter
            }

            elsif(($QArray[$qnum][$qpos] >= 2) && ($SArray[$snum][$spos] >= 2)) {
#                print "Both repeated\n";
                 push @{$SNPsHash{$qid}{$qpos}}, 'Both repeated'  
	    }

            elsif($QArray[$qnum][$qpos] < 2) {                  # add comment on repeat status
                push @{$SNPsHash{$qid}{$qpos}}, 'Subject repeated'
            }

            elsif($SArray[$snum][$spos] < 2) {                  # add comment on repeat status
                push @{$SNPsHash{$qid}{$qpos}}, 'Query repeated'
            }

        }

        elsif(($qnucl eq "-") && ($snucl =~ /[AGTCagtc]/))  {		# if subject has insertion	 					

            $snucl = ${$SNPsHash{$qid}{$qpos+1}}[3].$snucl;

            if($spos < $send) {						# in/de-crement subject base position tracker (mutation is NOT reported)
                $spos ++;
                $SNPsHash{$qid}{$qpos+1} = [($sid, $spos, $qnucl, $snucl, 'f')]
            }

            else {
                $spos --;
                $SNPsHash{$qid}{$qpos+1} = [($sid, $spos, $qnucl, $snucl, 'r')]
            }

            if(($QArray[$qnum][$qpos+1] >= 2) && ($SArray[$snum][$spos] >= 2)) {
                push @{$SNPsHash{$qid}{$qpos+1}}, 'Both repeated';
            }

            elsif($QArray[$qnum][$qpos+1] > 1) {                  # add comment on repeat status
                push @{$SNPsHash{$qid}{$qpos+1}}, 'Query repeated' 
            }

            elsif($SArray[$snum][$spos] > 1) {                  # add comment on repeat status
                push @{$SNPsHash{$qid}{$qpos+1}}, 'Subject repeated'
            }

        }  
	 	
        elsif(($qnucl =~ /[A-Za-z]/) && ($snucl eq "-")) {		# if subject has deletion

            $qpos ++;							# increment query base position tracker

            if($spos < $send) {
              $SNPsHash{$qid}{$qpos} = [($sid, $spos, $qnucl, $snucl, 'f')]
            }

            else {
              $SNPsHash{$qid}{$qpos} = [($sid, $spos, $qnucl, $snucl, 'r')]
            }

            if(($QArray[$qnum][$qpos] >= 2) && ($SArray[$snum][$spos] >= 2)) {
                push @{$SNPsHash{$qid}{$qpos}}, 'Both repeated';
            }

            elsif(($QArray[$qnum][$qpos] < 2) && ($SArray[$snum][$spos] >= 2)) {                  # add comment on repeat status
                push @{$SNPsHash{$qid}{$qpos}}, 'Subject repeated';
            }

            elsif(($SArray[$snum][$spos] < 2) && ($QArray[$qnum][$qpos] >= 2)) {                  # add comment on repeat status
                push @{$SNPsHash{$qid}{$qpos}}, 'Query repeated'
            }

	}


    }
    
    return ($qpos, $spos, $SNPs, $MMAligned)				# return results to ANALYZE_ALIGNMENT subroutine

}


# 10.  Print results to screen and file


sub RESULTS {

  my ($SNPs, $MAligned, $MMAligned, $BlastFile) = @_;

  print "Analysis complete.  Total SNPs = $SNPs\n";			# report total SNPs in uniquely aligned regions

  my $TotalAligned = $MAligned+$MMAligned;				# Sum up aligned based - in matched as well as mismatched alignments

  if($TotalAligned == 0) {

    print "Total SNPs = 0\nWeighted SNPs = 0\n";

    print MAINREPORT "0\t0\n";

    print CUMULATIVE "0\t0\n"

  }

  else {

    print  "Total aligned = $TotalAligned\n";				# print total uniquely aligned bases

    my $Weighted_SNPs = sprintf("%.0f", $SNPs/$TotalAligned * 1000000);	# generate a weighted SNPs count (SNPs per megabase uniquely aligned DNA)

    print  "Weighted SNPs (SNPs/Mb) = $Weighted_SNPs\n";				# report weighted SNPs

    print MAINREPORT "$SNPs\t";						# write above information to a file

    print CUMULATIVE "$SNPs\t";

    print MAINREPORT "$TotalAligned\t";

    print CUMULATIVE "$TotalAligned\t";

    print MAINREPORT "$Weighted_SNPs\n";

    print CUMULATIVE "$Weighted_SNPs\n";

    print SNPTEMPFILE "$BlastFile\n";

  }

  ($SNPs, $MAligned, $MMAligned, $TotalAligned) = (0, 0, 0, 0);		# reset counters ready for next BLAST file

  return($SNPs, $MAligned, $MMAligned, $TotalAligned)			# return value to main program, to ensure counters are reset

}

sub PRINT_ALIGN_STRINGS {


  foreach my $Gene (sort {$a cmp $b} keys %countString) {

#        print "*$Gene*\n";

    open(CHR, '>>', "$Gene"."_alignments") || die "Can't open for writing\n";

    foreach $Strain (sort {$a cmp $b} keys %{$countString{$Gene}}) {

      print CHR "$Strain\t$countString{$Gene}{$Strain}\n";

    }

  }

#  undef %countString;

  close CHR;

}

sub PRINT_SNPs {

  foreach $qid (sort {$a cmp $b} keys %SNPsHash) {

    foreach $qpos (sort {$a <=> $b} keys %{$SNPsHash{$qid}}) {

      @SNPsList = @{$SNPsHash{$qid}{$qpos}};
 
      $sid = shift @SNPsList;

      print OUT "$qid\t$sid\t$qpos\t";

      print OUT join ("\t", @SNPsList), "\n";

    }

  }

  close OUT;

  undef %SNPsHash

}
