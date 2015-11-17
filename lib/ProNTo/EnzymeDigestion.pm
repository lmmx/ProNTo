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
	    "missed-cleavages=i" => $n_missed_cleavages,
	    "help" => \$help_opt )) {
	    	error_out("No arguments recognized");
	    }


# Declaring arrays and variables used.
my (@sequences, @peptides, @mc_peptides);                       
my ($protein, $sequence, $peptide, $unusual_counter,$unknown_counter, $protein_size, $peptide_size, $mc_peptide_size);
my ($n,$r,$s,$i,$j);

# Array of proteins for testing the program, should be silenced with a hash when using the proteins from task 1.
my @proteins=qw(DAAAAATTLTTTAMTTTTTTCKMMFRPPPPPGGGGGGGGGGGG ALTAMCMNVWEITYHKGSDVNRRASFAQPPPQPPPPLLAIKPASDASD);

############################################################################
##############           Main subroutine                          ##########
############################################################################
sub main {  
	# If user selects option --help, goes to that subroutine and prints usage.  
	if ($help_opt) {
		help();
		exit;
	} 
	# Reading FASTA format 
	for $protein (@proteins) { 
		if ($protein !~ /^>/)   {     
			# Counting unusual and unknown aminoacids
#Only counting 1	if ($protein=~ m/(BOUJZ)/g) {$unusual_counter++;}
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
	# If user selects a number of missed cleavages, goes to subroutine and returns an
	# array with new peptides formed with $n misscleavages, where $n goes from 1 to the 
	# max allowed cleavages. I.e: if $n_missed_cleavages=3, should return array of peptides
	# with $n=1,2 and 3 misscleavages.
	if ($n_missed_cleavages > 0) {@mc_peptides=missed_cleavages(\@peptides);}  
	# Getting the lenght of the arrays containing the proteins, the peptides and the 
	# misscleavage peptides in order to use them in the subroutine to print results.	
	$protein_size=scalar(@sequences);
	$peptide_size=scalar(@peptides);
	$mc_peptide_size=scalar(@mc_peptides);
	# Goes to subroutine for printing the peptides obtained in the digestion.			
	print_peptides(\@sequences,\@peptides,\@mc_peptides,\$n);     
	# Exits the program.
	exit;
}
##### END of main sub #####

sub digestion { 
	my $sequences_ref=shift;
	my @peptides= ();
	my @sequences_ref= @{$sequences_ref};
	for (@sequences_ref) {
		
		#Checks user input and runs the matching pattern of the specified enzyme for digestion
		#If pattern matches, ads an = sign and splits there to divide the protein in peptides.
		if ($trypsin) {
			#Patern for Trypsine, it matches if the sequence has K or R, without P after.
			if ($_ =~ s{ (?<= [KR]) (?! P) }{=}xmsg) {                       
				push @peptides, split ('=', $_);
			}
		}
		elsif ($lysc) {
			#Patern for Endoproteinase Lys-C, it matches if the sequence has K , without P after.
			if ($_ =~ s{ (?<= [K]) (?! P) }{=}xmsg) 	{
				push @peptides, split ('=', $_);
			}
		}
		elsif ($argc) {
			#Patern for Endoproteinase Arg-C, it matches if the sequence has R, without P after.
			if ($_ =~ s{ (?<= [R]) (?! P) }{=}xmsg) {
				push @peptides, split ('=', $_);
			}
		}
		elsif ($gluc) {
			#Patern for V8 proteinase, it matches if the sequence has E, without P after.
			if ($_ =~ s{ (?<= [E]) (?! P) }{=}xmsg) {
				push @peptides, split ('=', $_);
			}
		}
		else {
			# If no enzymes are selected, error message and usage are printed on the screen
			error_out("No enzyme selected");                           
		}
	}
	# Returns array of digested peptides to main subroutine.
	return @peptides;
}

sub missed_cleavages {
	my $peptides_ref= shift;
	my @peptides_ref={@$peptides_ref};
	for $n(1...$n_missed_cleavages){
		my $counter = 0;
		for $counter($counter<scalar(@peptides)-$n,$counter++){
			my @mc_peptides = join ('',@peptides[$counter...$counter+$n]);
		}
	}
return @mc_peptides;
return $n;
}

sub help {
	print "\n\n Protein Digestion Options\n";
	print "========================================\n\n";
	print "Options:\n";
	print "[-t]: digests protein with trypsin\n";
	print "[-l]: digests protein with endoproteinase Lys-C\n";
	print "[-a]: digests protein with endoproteinase Arg-C\n";
	print "[-v]: digests protein with V8 proteinase\n";
	print "[-c]: number of allowed missed cleavages. Default value is 0.\n";
	print "[-h]: prints help\n\n";
	return;
}

sub error_out {
	my $error_msg = shift;
	help;
	die "$error_msg";
}

sub print_peptides {
        # Usage:
        #       parameter 1 is array of proteins (i.e. the sequences)
        #       param 2 is the array of peptides generated from digestion of the proteins
	#	param 3 is the array of additional peptides generated from missing the cleavages
        #       param 4 is the number of allowed misscleavages

	my $sequence_ref = shift;
        my $peptide_ref = shift;
        my $mc_peptide_ref = shift;
        my $n_ref= shift;
	my @sequences_ref = @{$sequence_ref};
	my @peptides_ref = @{$peptide_ref};
	my @mc_peptides_ref = @{$mc_peptide_ref};

	my $n_proteins = scalar @sequences_ref;
	print "Processing $n_proteins proteins\n";

	for my $i (0..$#sequences_ref) {
		my $printed_i = $i+1;
		print "Protein: $printed_i\n";
    		my $sequence_ref = $sequences_ref[$i];
    		my $j;
    		print scalar @peptides_ref . " peptides for protein $printed_i\n";
    			for my $j (0..$#peptides_ref) {
        			my $printed_j = $j+1;
        			print "> Protein: $printed_i|Peptide:$printed_j| $n_ref cleavages\n";
        			print "$peptides_ref[$j]\n";
    			}
	}
	return;
    
}

main();

# End of module evaluates to true

1;

__END__

# End of file evaluates to false

0;
