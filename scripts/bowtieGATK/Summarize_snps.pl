#!/usr/bin/perl

die "Usage: Summarize_snps.pl <VCF-dir>\n" if @ARGV < 1;

die "Program expects a DIRECTORY as argument\n" unless -d $ARGV[0];


opendir(DIR, $ARGV[0]) || die "Can't open directory\n";

@filesList = readdir(DIR);

foreach $file (@filesList) {

  next unless $file =~ /SSfilter/;

  open(F, "$ARGV[0]/$file")|| die "Can't open file\n";

  while($L=<F>) {

    next if $L =~ /^#/;

    next if $L =~ /FAIL/;

    chomp($L);

    @L = split(/\t/, $L);

    if ($file =~ /^(.+)?_geno.+vcf/) {		# skip over data for non-wheat blast samples

      $ERR = $1;

      push @{$Hash{$L[0]}{$L[1]}}, $ERR;

    }

  }

}

  foreach $chr (sort {$a cmp $b} keys %Hash) {
 
    foreach $pos (sort {$a <=> $b} keys %{$Hash{$chr}}) {

      $arrayLen = @{$Hash{$chr}{$pos}};

      print "$chr\t$pos\t$arrayLen\t@{$Hash{$chr}{$pos}}\n";

  #      print "BGL SNPs: $BGLsnps\n" if $pop eq 'BGL'; 	# uncomment for total SNPs summaries

  #      print "ZMB SNPs: $ZMBsnps\n" if $pop eq 'ZMB';         # uncomment for total SNPs summaries

    }

  }


