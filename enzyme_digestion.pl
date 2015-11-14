#!/usr/bin/env perl
######################################################################################
## enzyme_digestion.pl
#######################################################################################

use warnings;
use strict;
use diagnostics;

#Defining user options
use vars qw($opt_t $opt_l $opt_a $opt_v $opt_c $opt_h);       
use Getopt::Std;

#Setting missed cleavages to a default value of 0   
my $n_missed_cleavages = (defined $opt_c) ? $opt_c : 0;    

#Declaring arrays and variables
my (@proteins, @sequences, @peptides, @mc_peptides);                       
my ($protein, $sequence, $peptide, $unusual_counter,$unknown_counter, $protein_size, $peptide_size, $mc_peptide_size);
my ($n,$r,$s,$i,$j);
############################################################################
##############           Main subroutine                          ##########
############################################################################
@proteins=("DAAAAATTLTTTAMTTTTTTCKMMFRPPPPPGGGGGGGGGGGG","ALTAMCMNVWEITYHKGSDVNRRASFAQPPPQPPPPLLAIKPASDASD");
main();
sub main {    
	# If there is no user input, error message and usage are displayed and program closes
	if (!getopts ('tlavc:h')) {
		error_out("No arguments recognized");
	}
	if ($opt_h) {
		help();
		exit;
	} 
	# Reading FASTA file and eliminating header and blank spaces    
	foreach $protein (@proteins) { 
		if ($protein !~ /^>/)   {     
			# Counting unusual and unknown aminoacids
			if ($protein=~ m/(BOUJZ)/g) {$unusual_counter++;}
			if ($protein=~ m/X/g) {$unknown_counter++;}
			$protein =~ s/\s//g;                
                        @sequences = chomp ($protein);
         }
	} 
	foreach $sequence (@sequences) {                            
	    #Checks user input and goes to the specified enzyme subroutine for every sequence (Protein)	                                                         
		if ($opt_t) {@peptides = trypsine(\$sequence);}
		elsif ($opt_l) {@peptides = endo_lysc(\$sequence);}
		elsif ($opt_a) {@peptides = endo_argc(\$sequence);}
		elsif ($opt_v) {@peptides = _v8(\$sequence);}
		else {
			# If no enzymes are selected, error message and usage are printed on the screen
			error_out("No enzyme selected");                           
		}
	}
	#If user selects a number of missed cleavages, goes to subroutine and returns and array with new peptides 
	if ($n_missed_cleavages > 0) {@mc_peptides=missed_cleavages(\@peptides);}  
	$protein_size=scalar(@sequences);
	$peptide_size=scalar(@peptides);
	$mc_peptide_size=scalar(@mc_peptides);
	print_peptides();     
	exit;
}
##### END of main sub #####

# Subroutine for trypsine, looks for regex, and splits sequence into peptides
sub trypsine {                                                
    my $sequence = shift; 
	if ($sequence =~ m/(KR)!P/g) {                       
		$sequence =~ s/(KR)/$1=/;
		@peptides = split ('=',$sequence);
	}
	return @peptides;                           
}

sub endo_lysc {                                                 
    my $sequence = shift; 
	if ($sequence =~ m/K!P/g) 	{
		$sequence =~ s/K/$1=/;
		@peptides = split ('=', $sequence);
	}
	return @peptides;
}

sub endo_argc {
    my $sequence = shift; 
	if ($sequence =~ m/R!P/g) {
		$sequence =~ s/R/$1=/;
		@peptides = split ('=', $sequence);
	}
	return @peptides;
}

sub _v8 {
	my $sequence = shift; 
    if ($sequence =~ m/E!P/g) {
		$sequence =~ s/E/$1=/;
		@peptides = split ('=', $sequence);
	}
	return @peptides;
}

sub missed_cleavages {
my $peptides= shift;
my @peptides={@$peptides};
foreach $n(1...$n_missed_cleavages){
	my $counter = 0;
	foreach $counter($counter<scalar(@peptides)-$n,$counter++){
	my @mc_peptides = join ('',@peptides[$counter...$counter+$n]);
	}
}
return @mc_peptides;
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
	print "USAGE:$0\n";
	print "\t[-t]: digests protein with trypsin\n";
	print "\t[-l]: digests protein with endoproteinase Lys-C\n";
	print "\t[-a]: digests protein with endoproteinase Arg-C\n";
	print "\t[-v]: digests protein with V8 proteinase\n";
	print "\t[-c]: number of allowed missed cleavages. Default value is 0.\n";
	print "\t[-h]: prints help\n\n";
	die "$error_msg";
}

sub print_peptides {
print "$unusual_counter unusual aminoacids (BOUJZ) were found, $unknown_counter unknown aminoacids (X) were found";
foreach ( $i=1, $i<$protein_size, $i++){
        foreach ( $j=1, $j<$peptide_size, $j++){
                print "Protein[$i] Peptide[$j] No cleavages\n";
                print " $peptides[$j]\n";
        }
foreach ( $r=1, $r<$protein_size, $r++){
        foreach ( $s=1, $s<$mc_peptide_size, $s++){
                print "Protein[$r] Peptide[$s] $n cleavages\n";
                print " $mc_peptides[$s]\n";
        }
}
return;
}
}
