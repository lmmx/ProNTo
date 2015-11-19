#!/usr/bin/env perl

package ProNTo::enzymeDigestion;

use strict;
#use diagnostics;
use warnings;
use Getopt::Long;

use File::Spec; # for class methods related to file path parsing
use File::Basename; # for dirname function

# Initialising variables in scope, then assign some of them in BEGIN block

my ($abslibfile, $libdirname, $n_missed_cleavages, $trypsin, $lysc, $argc, $gluc, $help_opt);

BEGIN {
    # __FILE__ provides a relative path to the bin file
    # According to http://stackoverflow.com/a/90721/2668831 this is the most robust method
    $abslibfile = File::Spec->rel2abs(__FILE__);
    # use dirname x2 to go "up" 2 directories from bin file, then descend into the lib path:
    $libdirname = dirname($abslibfile);
    # set some default variable values
	$n_missed_cleavages = 0;
};

#load modules in the [path]/lib/ProNTo directory
use lib $libdirname;
use Methods;

__PACKAGE__->main(@ARGV) unless caller;

# Defining user options.
# If there is no user input, error message and usage are displayed and program closes.

sub handleopts {
	if (!GetOptions ("trypsin" => \$trypsin,
			"lys-c" => \$lysc,
			"arg-c" => \$argc,
			"glu-c" => \$gluc,
			"missed-cleavages=i" => $n_missed_cleavages,
			"help" => \$help_opt )
		){
			die ProNTo::Methods::help('NO_ENZ');
		}
	if (defined $help_opt) {
		print STDERR ProNTo::Methods::help();
	}
	
}

# Declaring arrays and variables used.
my (@sequences, @peptides, @mc_peptides);                       
my ($protein, $sequence, $peptide, $unusual_counter,$unknown_counter, $protein_size, $peptide_size, $mc_peptide_size);
my ($n,$r,$s,$i,$j);

# Array of proteins for testing the program, should be silenced with a hash when using the proteins from task 1.
my @proteins=qw(
	DAAAAATTLTTTAMTTTTTTCKMMFRPPPPPGGGGGGGGGGGG
	ALTAMCMNVWEITYHKGSDVNRRASFAQPPPQPPPPLLAIKPASDASD
	);

sub main {
	handleopts();
#	reader();
	print "Testing...";
	exit;
}

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
			ProNTo::Methods::help('noenz');
		}
	}
	# Returns array of digested peptides to main subroutine.
	return @peptides;
}

sub reader {
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
	print_peptides();
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
}

sub print_peptides {
	print "$unusual_counter unusual aminoacids (BOUJZ) were found, $unknown_counter unknown aminoacids (X) were found";
	for ( $i=1, $i<$protein_size, $i++){
			for ( $j=1, $j<$peptide_size, $j++){
					print "Protein[$i] Peptide[$j] No cleavages\n";
					print " $peptides[$j]\n";
			}
	}
	for ( $r=1, $r<$protein_size, $r++){
			for ( $s=1, $s<$mc_peptide_size, $s++){
					print "Protein[$r] Peptide[$s] $n cleavages\n";
					print " $mc_peptides[$s]\n";
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