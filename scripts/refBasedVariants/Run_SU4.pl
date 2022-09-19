#!/usr/bin/perl

print "Usage: perl Run_SU4.pl <blast-dir> <snps-dir>\n" if @ARGV != 2;

use StrictUnique4;

StrictUnique4::SNPs(\@ARGV);

mkdir "$ARGV[1]/ALIGN";

#system("mv $ARGV[1]/*_alignments $ARGV[1]/ALIGN/");
