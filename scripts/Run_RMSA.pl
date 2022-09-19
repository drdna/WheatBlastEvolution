#!/usr/bin/perl

use RMSA;

die "Usage: perl Run_RMSA.pl <genome_to_mask>\n" if @ARGV != 1;

RMSA::RUN($ARGV[0]);


