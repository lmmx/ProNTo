#!/usr/bin/env perl
######################################################################################
## enzyme_digestion.pl
#######################################################################################

use warnings;
use strict;
use diagnostics;

#Defining user options
# If there is no user input, error message and usage are displayed and program closes
use Getopt::Long;
if (!GetOptions ("trypsin=s" => \$trypsin,
	    "lys-c=s" => \$lysc,
	    "arg-c=s" => \$argc,
	    "glu-c=s" => \$gluc,
	    "missed-cleavages=i" => $n_missed_cleavages,
	    "help=s" => \$help_opt )) {
	    	error_out("No arguments recognized");
	    }
#Setting missed cleavages to a default value of 0   
my $n_missed_cleavages = (defined $n_missed_cleavages) ? $n_missed_cleavages : 0;    

#Declaring arrays and variables
my (@sequences, @peptides, @mc_peptides);                       
my ($protein, $sequence, $peptide, $unusual_counter,$unknown_counter, $protein_size, $peptide_size, $mc_peptide_size);
my ($n,$r,$s,$i,$j);
############################################################################
##############           Main subroutine                          ##########
############################################################################
my @proteins=qw(DAAAAATTLTTTAMTTTTTTCKMMFRPPPPPGGGGGGGGGGGG ALTAMCMNVWEITYHKGSDVNRRASFAQPPPQPPPPLLAIKPASDASD);
main();
sub main {    
	if ($help_opt) {
		help();
		exit;
	} 
	# Reading FASTA file and eliminating header and blank spaces    
	for $protein (@proteins) { 
		if ($protein !~ /^>/)   {     
			# Counting unusual and unknown aminoacids
			if ($protein=~ m/(BOUJZ)/g) {$unusual_counter++;}
			if ($protein=~ m/X/g) {$unknown_counter++;}
			$protein =~ s/\s//g;
                        push @sequences, $protein;
        	}
	} 
	for $sequence (@sequences) {                            
	    #Checks user input and goes to the specified enzyme subroutine for every sequence (Protein)	                                                         
		if ($trypsin) {@peptides = trypsin(\$sequence);}
		elsif ($lysc) {@peptides = endo_lysc(\$sequence);}
		elsif ($argc) {@peptides = endo_argc(\$sequence);}
		elsif ($gluc) {@peptides = _v8(\$sequence);}
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

# Subroutine for trypsin, looks for regex, and splits sequence into peptides
sub trypsin {                                                
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
for $n(1...$n_missed_cleavages){
	my $counter = 0;
	for $counter($counter<scalar(@peptides)-$n,$counter++){
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
	help;
	die "$error_msg";
}

sub print_peptides {
print "$unusual_counter unusual aminoacids (BOUJZ) were found, $unknown_counter unknown aminoacids (X) were found";
for ( $i=1, $i<$protein_size, $i++){
        for ( $j=1, $j<$peptide_size, $j++){
                print "Protein[$i] Peptide[$j] No cleavages\n";
                print " $peptides[$j]\n";
        }
for ( $r=1, $r<$protein_size, $r++){
        for ( $s=1, $s<$mc_peptide_size, $s++){
                print "Protein[$r] Peptide[$s] $n cleavages\n";
                print " $mc_peptides[$s]\n";
        }
}
return;
}
}
