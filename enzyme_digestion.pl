#!/usr/bin/perl -w
use strict;

use vars qw($opt_t $opt_l $opt_a $opt_v $opt_c $opt_h);       #Defining user options
use Getopt::Std;
#getopts ('tlavc:h');

my $missed_cleavages = (defined $opt_c) ? $opt_c : 0;         #Setting missed cleavages to a default value of 0     
my (@proteins, @sequences, @peptides);                        #Declaring arrays and variables
my ($protein, $sequence, $peptide,$unusual_counter,$unknown_counter);

############################################################################
##############           Main subroutine                          ##########
############################################################################
@proteins=("DAAAAATTLTTTAMTTTTTTCKMMFRPPPPPGGGGGGGGGGGG","ALTAMCMNVWEITYHKGSDVNRRASFAQPPPQPPPPLLAIKPASDASD");
main();
sub main
{                                                             
	if (!getopts ('tlavc:h'))              # If there is no user input, error message and usage are displayed and program closes
	{
		print "\nNo arguments recognized\n";
		error_message();
		exit;
	}
	if ($opt_h) {help();}                  
	foreach $protein (@proteins)          # Reading FASTA file and eliminating header and blank spaces       
	{ 
                while $protein
			{
			if (=~ m/(BOUJZ)/g) {$unusual_counter++;}
			if (=~ m/X/g) {$unknown_counter++;}
 			}
		if ($protein !~ /^>/)   
		{         
			$protein =~ s/\s//g;                
                        @sequences = chomp ($protein);

                }
#        print "holacaracola";
	} 
	foreach $sequence (@sequences)                             #Checks user input and goes to the specified enzyme subroutine for 
	{                                                          # every sequence(Protein)
	if ($opt_t) {@peptides = trypsine(\$sequence);}
	elsif ($opt_l) {@peptides = endo_lysc(\$sequence);}
	elsif ($opt_a) {@peptides = endo_argc(\$sequence);}
	elsif ($opt_v) {@peptides = _v8(\$sequence);}
	else                                                      # If no enzymes are selected, error message and usage are printed on the screen
	{
		print "\nNo enzyme selected\n\n";
		error_message();                           
	}
	}	                                                                       
	if ($missed_cleavages > 0) {@mc_peptides=missed_cleavages(\@peptides);}        #If user selects a number of missed cleavages, goes to 
	print_peptides();                                                              #subroutine and returns and array with new peptides
}

##### END of main sub #####
sub trypsine 
{                                                
	my ($sequence)= @_;        
	if ($sequence =~ m/(KR)!P/g)             # Subroutine for trypsine, looks for regex, and splits sequence into peptides
	{                       
		$sequence =~ s/(KR)/$1=/;
		@peptides = split ('=',$sequence);
	}
	return @peptides;                           
}

sub endo_lysc 
{                                                 
        my ($sequence) = @_; 
	if ($sequence =~ m/K!P/g) 
	{
		$sequence =~ s/K/$1=/;
		@peptides = split ('=', $sequence);
	}
	return @peptides;
}

sub endo_argc 
{
        my ($sequence) = @_; 
	if ($sequence =~ m/R!P/g) 
	{
		$sequence =~ s/R/$1=/;
		@peptides = split ('=', $sequence);
	}
	return @peptides;
}

sub _v8 
{
	my ($sequence) = @_; 
        if ($sequence =~ m/E!P/g) 
	{
		$sequence =~ s/E/$1=/;
		@peptides = split ('=', $sequence);
	}
	return @peptides;
}

sub missed_cleavages 
{
my ($peptides)= @_;
my $peparray=@$peptides;
foreach my $n(1...$missed_cleavages)
{
	my $counter = 1;
	for $counter($counter<length($peparray)-$n,$counter++)
	{
	my @mc_peptides = join @peptides ('',$counter,$counter(1...$n));
	}
}
return @mc_peptides;
}
sub help 
{
print "\n\n Protein Digestion Options\n";
print "========================================\n\n";
print "Options:\n";
print "[-t]: digests protein with trypsin\n";
print "[-l]: digests protein with endoproteinase Lys-C\n";
print "[-a]: digests protein with endoproteinase Arg-C\n";
print "[-v]: digests protein with V8 proteinase\n";
print "[-c]: number of allowed missed cleavages. Default value is 0.\n";
print "[-h]: prints help\n\n";
exit;
}

sub error_message
{
print "USAGE:$0\n";
print "\t[-t]: digests protein with trypsin\n";
print "\t[-l]: digests protein with endoproteinase Lys-C\n";
print "\t[-a]: digests protein with endoproteinase Arg-C\n";
print "\t[-v]: digests protein with V8 proteinase\n";
print "\t[-c]: number of allowed missed cleavages. Default value is 0.\n";
print "\t[-h]: prints help\n\n";
exit;
}

sub print_peptides
{
print "@peptides @mc_peptides @unusual_counter @unknown_counter";
}
