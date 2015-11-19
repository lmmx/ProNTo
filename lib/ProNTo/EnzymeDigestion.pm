#!/usr/bin/env perl

package ProNTo::enzymeDigestion;
######################################################################################
## enzyme_digestion.pl
#######################################################################################

use strict;
use diagnostics;
use warnings;
use Getopt::Long;

# Defining user options.
# If there is no user input, error message and usage are displayed and program closes.

my ($trypsin, $lysc, $argc, $gluc, $help_opt );

# Setting missed cleavages to a default value of 0.  
my $n_missed_cleavages = 0;

if (!GetOptions ("trypsin" => \$trypsin,
	    "lys-c" => \$lysc,
	    "arg-c" => \$argc,
	    "glu-c" => \$gluc,
	    "missed-cleavages=i" => \$n_missed_cleavages,
	    "help" => \$help_opt )) {
	    	error_out("No arguments recognized");
	    }


# Declaring arrays and variables used.
my (@sequences, @peptides, @mc_peptides, @sequences_ref);                       
my ($protein, $sequence, $peptide, $unusual_counter,$unknown_counter, $peptides_ref);
my ($n,$r,$s,$i,$j);

# Array of proteins for testing the program, should be silenced with a hash when using the proteins from task 1.
#my @proteins=qw(DAAAAATTLTTTAMTTTTTTCKMMFRPPPPPGGGGGGGGGGGG ALTAMCMNVWEITYHKGSDVNRRASFAQPPPQPPPPLLAIKPASDASD);

############################################################################
##############           Main subroutine                          ##########
############################################################################
sub main {  
	# If user selects option --help, goes to that subroutine and prints usage. 
	# Then exits the program 
	if ($help_opt) {
		help();
		exit;
	} 
	# Reading FASTA format 
	for $protein (@proteins) { 
		if ($protein !~ /^>/)   {     
			# Counting unusual and unknown aminoacids
#Only counting 1	if ($protein=~ m/[BOUJZ]/g) {$unusual_counter++;}
#per protein		if ($protein=~ m/X/g) {$unknown_counter++;}
			# Eliminating header and blank spaces    
			$protein =~ s/\s//g;
			# Saving protein sequences into a new array.
                        push @sequences, $protein;
        	}
#		print "$unusual_counter unusual aminoacids (BOUJZ) were found, $unknown_counter unknown aminoacids (X) were found\n";
	} 
	# Goes to the subroutine for digestion of the proteins and returns and array of peptides
	# obtained from the digestion.
	@peptides=digestion(\@sequences);
	# Goes to subroutine for printing the peptides obtained in the digestion.			
	print_peptides(\@sequences,\@peptides);  
	# If user selects a number of missed cleavages, goes to subroutine and prints an
	# array with new peptides formed with $n misscleavages, where $n goes from 1 to the 
	# max allowed cleavages. I.e: if $n_missed_cleavages=3, should return array of peptides
	# with $n=1,2 and 3 misscleavages.
	if ($n_missed_cleavages > 0) {@mc_peptides=missed_cleavages(\@peptides,\@sequences,$n_missed_cleavages);}     
	# Exits the program.
	exit;
}
##### END of main sub #####

# Subroutine for protein digestion, receives the protein array from the main sub and processes it with the selected enzyme
sub digestion { 
	# First parameter is the reference of the sequences array
	my $sequences_ref=shift;
	# Declaring a new array for the peptides
	my @peptides= ();
	# Dereferencing the sequences array
	my @sequences_ref= @{$sequences_ref};
	for (@sequences_ref) {
		
		# Checks user input and runs the matching pattern of the specified enzyme for digestion
		# If pattern matches, ads an = sign and splits there to divide the protein in peptides.
		# Then pushes peptide to peptide array.
		if ($trypsin) {
			# Pattern for Trypsine, it matches if the sequence has K or R, without P after.
			if ($_ =~ s{ (?<= [KR]) (?! P) }{=}xmsg) {                       
				push @peptides, split ('=', $_);
			}
		}
		elsif ($lysc) {
			# Pattern for Endoproteinase Lys-C, it matches if the sequence has K , without P after.
			if ($_ =~ s{ (?<= [K]) (?! P) }{=}xmsg) 	{
				push @peptides, split ('=', $_);
			}
		}
		elsif ($argc) {
			# Pattern for Endoproteinase Arg-C, it matches if the sequence has R, without P after.
			if ($_ =~ s{ (?<= [R]) (?! P) }{=}xmsg) {
				push @peptides, split ('=', $_);
			}
		}
		elsif ($gluc) {
			# Pattern for V8 proteinase, it matches if the sequence has E, without P after.
			if ($_ =~ s{ (?<= [E]) (?! P) }{=}xmsg) {
				push @peptides, split ('=', $_);
			}
		}
		else {
			# If no enzymes are selected, error message and usage are printed on the screen
			error_out("No enzyme selected");                           
		}
	}
	# Returns array of peptides to main subroutine.
	return @peptides;
}

# This subroutine is called if the user selects at least one missed cleavage for the digestion.
# Joins concatenated peptides from the same protein to simulate the enzyme missing a splitting spot.
sub missed_cleavages {
	# Fist parameter is the normal peptides array.
	# Second parameter is the proteins sequence array.
	# Third parameter is the number of max missed cleavaes. I.e: if $n_missed_cleavages is set to 2,
	# Will calculate new peptides for 1 and 2 missed cleavages.
	my $peptide_ref= shift;
	my $sequence_ref=shift;
	my $missed_cleavages_ref=shift;
	# Dereferencing both arrays
	my @sequences_ref=@{$sequence_ref};
	my @peptides_ref=@{$peptide_ref};
	# Creating an array containing allowed missed cleavages.
	my @n_cleav = (1...$missed_cleavages_ref);

	# This loop creates and prints the missed cleavage peptides
	print "\tPeptides generated by misscleavage:\n";
	for my $k(0...$#sequences_ref){
		# Defining variable to start printing proteins from 1
		my $printed_k=$k+1;
		my $sequence_ref = $sequences_ref[$k];
		for my $n (@n_cleav) {
			# Pep is the starting peptide for the cleavage.
			for my $pep (0...@peptides_ref-1) {
				# Defining variable to start printing peptides from 1
				my $printed_pep=$pep+1;
				# This is the ending position of the missed cleavage, the sum of the first peptide used
				# And the number of missed_cleavages
				my $end_pep = $pep + $n;
				# Creating new array to store new peptides
				my $mc_pep;
				# If the sum of the starting peptide and number of mc is bigger than the number of peptides
				# of the protein, this becomes the last step...
				if ($end_pep>=@peptides_ref) {
					last;
				}
				# ... If not, this variable takes the value of the starting peptide, the last, and
				# The ones between them, and concatenates them one by one, in sequence order
				for my $pos ($pep...$end_pep) {
					$mc_pep = $mc_pep.$peptides[$pos];
				}
				push @mc_peptides, $mc_pep;
				# Printing the new peptides
				print "> Protein: $printed_k|Peptide:$printed_pep| $n cleavages\n";
        			print "$mc_peptides[$pep]\n\n";
			}
		}
	}
	# Return to main sub
	return;
}

# This sub prints the command line options if the user calls it from the terminal
# It is also used if the user writes a non defined argument or doesn't select any
sub help {
	print "\n\n Protein Digestion Options\n";
	print "========================================\n\n";
	print "Options:\n";
	print "[--trypsin]: digests protein with trypsin\n";
	print "[--lys-c]: digests protein with endoproteinase Lys-C\n";
	print "[--arg-c]: digests protein with endoproteinase Arg-C\n";
	print "[--glu-c]: digests protein with V8 proteinase\n";
	print "[--missed-cleavages]: number of allowed missed cleavages. Default value is 0.\n";
	print "[--help]: prints help\n\n";
	return;
}

# If there is no command line option selected or the ones written are not recognised
# Prints error and help.
sub error_out {
	my $error_msg = shift;
	help;
	die "$error_msg";
}

# This subroutine does the printing for the normal peptides
sub print_peptides {
        # Parameter 1 is array of proteins (i.e. the sequences)
        # Param 2 is the array of peptides generated from digestion of the proteins
	my $sequence_ref = shift;
        my $peptide_ref = shift;
	# Dereferencing
	my @sequences_ref = @{$sequence_ref};
	my @peptides_ref = @{$peptide_ref};

	# Number of total proteins
	my $n_proteins = scalar @sequences_ref;
	print "\tProcessing $n_proteins proteins\n\n";
	
	# Loop for printing the peptides of each protein
	# Follows same order than the printing of mc_peptides
	for my $i (0..$#sequences_ref) {
		my $printed_i = $i+1;
    		my $sequence_ref = $sequences_ref[$i];
    		my $j;
    		print "\t". scalar @peptides_ref . " peptides for protein $printed_i\n";
    			for my $j (0..$#peptides_ref) {
        			my $printed_j = $j+1;
        			print "> Protein: $printed_i|Peptide:$printed_j| No cleavages\n";
        			print "$peptides_ref[$j]\n\n";
    			}
	}
	return;
    
}

# Calls main subroutine to run the program
main();

# End of module evaluates to true

1;

__END__

# End of file evaluates to false

0;
