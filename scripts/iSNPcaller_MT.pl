#!/usr/bin/perl

##############################
#
# written by Mark L. Farman
#
# Purpose: Mask genome sequences, perform reciprocal pairwise blast searches against each other 
#
# Note: works in incremental fasion: only newly added genomes are processed when main program is run
#
# Usage: perl iSNPcaller_MT.pl <Project_name> [hit run when prompted]
#
##############################

# use strict;

use warnings;

#use File::Find::Rule;

use File::Copy;

use File::Path;

use constant MAX_PROCESSES => 48;

use Parallel::ForkManager;

use Repeatmaskerfast;

use Snpcallerfast;

##############################

# Declare global variables;

my @To_be_BLASTed;

my %NewGenomesHash;

##############################

# Read command line argument

if(@ARGV) {

    my $Run = $ARGV[0];

    if($Run =~ /[RUN|run]/) {

      chdir("IMPB");

      & RUN_PROGRAM

    }

}

else {

    & INITIAL_DIALOG;

}


##############################

# Initiate user dialog


sub INITIAL_DIALOG {

    print "Is this the first time you have run the Incremental_Multipairwise_Blast script?\n\n";

    my $Answer = <STDIN>;

    chomp($Answer);

   if($Answer =~ /^[Y|y]/) {

        & MAKE_DIRECTORIES

    }

   elsif($Answer =~ /^[N|n]/) {

        chdir("IMPB");

        & RUN_PROGRAM

    }

    else {

        while($Answer !~ /^[YyNn]/) {

            print "Unrecognized entry.  Please try again\n\n";

            & INITIAL_DIALOG;

        }

    }

}


# Create directories for inputs and results

sub MAKE_DIRECTORIES {

        mkpath( 'IMPB/GENOMES/PROCESSED');

        mkpath( 'IMPB/NEW_GENOMES/FASTA');
	 
        mkpath( 'IMPB/NEW_GENOMES/BLAST_RESULTS');
	
        mkpath( 'IMPB/PROCESSED/FASTA');
	
        mkpath( 'IMPB/PROCESSED/BLAST_RESULTS');

        print 	"Directories used by the program have been created under the IMPB top directory. ".
		"Copy your genome .fasta files into the IMPB/GENOMES/ directory and then ".
		"type: run\n\n";

        chdir("IMPB"); 

        & GET_RUN_INPUT;

}


sub GET_RUN_INPUT {

    my $Answer = <STDIN>;

    chomp($Answer);

    if($Answer =~ /^[R|r]/) {

	print "\nNOTE: Next time you run this program, you can bypass the initial dialog by typing: ".

              "perl Multipairwise_BLAST.pl run\n\n";

        & RUN_PROGRAM

    }

    else {

        while($Answer !~ /^[R|r]/) {

	    print "Unrecognized entry.  Please try again\n\n";

	    & GET_RUN_INPUT;

	}

    }

}

sub RUN_PROGRAM {

    # check to make sure required directories are in place

    # set base directory as current working directory
    my $base_dir = shift // 'IMPB';

#    my $find_rule = File::Find::Rule->new;

    # Do not descend past the first level
    #$find_rule->maxdepth(1);

    # Only return directories
    #$find_rule->directory;

    # Apply the rule and retrieve the subdirectories
    #my @sub_dirs = $find_rule->in(".");

    #unless(grep(/NEW_GENOMES/, @sub_dirs)) {

                die "Required NEW_GENOMES and/or PROCESSED directories apear to be missing. "
                ."Please create them in the working directory and populate them with "
                ."the appropriate genome sequences."

    }

    & MASK_GENOMES;

    & ADD_NEW_GENOMES;

    & BLAST_NEW_GENOMES;                    # Blast the new genomes against one another

    & BLAST_NEW_V_OLD;		               # Blast the new genomes against the ones that have already been processed

    & MOVE_NEW_FASTAs_TO_PROCESSED_DIR;     # Move new .fasta files to FASTA folder in PROCESSED directory

    print   "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n".
            "Running the SNP caller\n";

    Snpcallerfast::SNPs;

    & MOVE_OLD_BLAST_RESULTS;               # move old BLAST results from NEW_GENOME
                                                # directory into PROCESSED directory
}

sub MASK_GENOMES {

    Repeatmaskerfast::MASK();

    opendir(GENOMES, "GENOMES");

    my @FileList = readdir(GENOMES);

    foreach my $Masked_file (@FileList) {

        if($Masked_file =~ /masked.+fasta$/) {

            move("GENOMES/$Masked_file", "NEW_GENOMES/FASTA/$Masked_file")

        }

        elsif($Masked_file =~ /fasta$/) {

            move("GENOMES/$Masked_file", "GENOMES/PROCESSED")

        }

    }

}


sub ADD_NEW_GENOMES {

    my @New;

    my $New_genome;

    # open directory containing the new genome sequences

    opendir(NEW, "NEW_GENOMES/FASTA");

    @New = readdir(NEW);

    foreach $New_genome (@New) {

        if($New_genome =~ /\.fasta$|\.fsa$/) {
			
    	    $NewGenomesHash{$New_genome} = 1;

            push @To_be_BLASTed, $New_genome;               # add genome to the "To_be_BLASTed" list
			
	}

    }

}

sub BLAST_NEW_GENOMES {
	
	my %New_genome;

# Perform reciprocal pairwise BLASTs between all new genomes
	
# Takes each genome in To_be_BLASTed list and BLASTs against all others

 	my $Num_new_genomes = @To_be_BLASTed;

	if($Num_new_genomes > 1) {
	
		print "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n".
		      "BLASTing new genomes against one another\n\n";
	
		my $i = 0;

	        my $j = 0;

		### note: this algorithm needs to be rewritten to maximize multi-threading capability

                ### run the nested loops and write genome comparison list to a single array

                ### then loop through the single array                 

                for($i = 0; $i <= $Num_new_genomes-2; $i++) {

			$New_genome{$i} = 1;

                        my $pm = Parallel::ForkManager->new(MAX_PROCESSES);

			for ($j = $i + 1; $j <= $Num_new_genomes-1; $j++) {

				my $pid = $pm->start and next;

				my ($Q, $SUB) = ($To_be_BLASTed[$i], $To_be_BLASTed[$j]);

		                $Q =~ s/\.fasta//;

		                $SUB =~ s/\.fasta//;

		 		$New_genome{$j} = 1;
				
				print "BLASTing query $To_be_BLASTed[$i] against subject $To_be_BLASTed[$j]\n";
				
				exec("blastn -subject NEW_GENOMES/FASTA/$To_be_BLASTed[$j] -query NEW_GENOMES/FASTA/$To_be_BLASTed[$i]".
				 " -out NEW_GENOMES/BLAST_RESULTS/$Q.$SUB.BLAST".
				 " -max_target_seqs 20000 -evalue 1e-20 -outfmt '6 qseqid sseqid qstart qend sstart send btop'".
                 	        " 2>/dev/null >>exceptions.txt") || die $!
			
			}
			$pm->wait_all_children()
		}

	}

}


sub MOVE_NEW_FASTAs_TO_PROCESSED_DIR {
	
    # move new genome .fasta files to the PROCESSED directory
	
    foreach my $Genome (@To_be_BLASTed) {

	move("NEW_GENOMES/FASTA/$Genome", "PROCESSED/FASTA/$Genome")

    }

}

sub BLAST_NEW_V_OLD {

	print 	"\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n".
		"BLASTing new genome(s) against old ones\n\n";

# perform pairwise BLASTs between new and old genomes

# Reads PROCESSED/FASTA directory, reads To_be_BLASTed list
# performs blasts between "new" and "already processed" genomes

	opendir(PROCESSED, "PROCESSED/FASTA");
	
	my @Processed = readdir(PROCESSED);

	my $Fasta_files = "no";
		
	foreach my $Processed_genome (@Processed) {
			
		if(exists($NewGenomesHash{$Processed_genome})) {

			next

		}

		elsif($Processed_genome =~ /\.fasta$/) {

			$Fasta_files = "yes";

                        my $SUB = $Processed_genome;

	                $SUB =~ s/\.fasta$//;
			
			foreach my $New_genome (@To_be_BLASTed) {

	           		my $Query = $New_genome;

	                        $Query =~ s/\.fasta$//;

				if($New_genome eq $Processed_genome) {
					
					next
				
				}
				
				else  {
				
					print "BLASTing query $New_genome against subject $Processed_genome\n";

	                	   	system("blastn -subject PROCESSED/FASTA/$Processed_genome -query NEW_GENOMES/FASTA/$New_genome".
						" -out 'NEW_GENOMES/BLAST_RESULTS/$Query.$SUB.BLAST'".
						" -max_target_seqs 20000 -evalue 1e-20 -outfmt '6 qseqid sseqid qstart qend sstart send btop'".
                                 		" 2>/dev/null >>exceptions.txt");

                		}

            		}

	        }

    	}

	if($Fasta_files eq "no") {

		print "\nNo old genomes to analyze yet\n"

	}

}

sub MOVE_OLD_BLAST_RESULTS {

  system('mv NEW_GENOMES/BLAST_RESULTS/* PROCESSED/BLAST_RESULTS/')

}
