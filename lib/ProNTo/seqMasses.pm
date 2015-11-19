#!/usr/bin/env perl

package ProNTo::seqMasses;

use diagnostics;
use strict;
use warnings;
use Getopt::Long;

# must be 2 arguments:
#  --mass-type : average or monoisotopic (accepted values: "a" or "m" - without quotes)
#  --pep-data : path to peptide FASTA file

my ($z,
	$pep_data,
	$mass_type
    );
BEGIN {
        $z = 1;
        $pep_data = $ARGV[0];
        $mass_type = "m";
};

__PACKAGE__->main(@ARGV) unless caller;

sub main {
        handleOpts();
        readPepFASTA();
}

sub handleOpts {
        # Defining bin size and parameters
        if (!GetOptions (
                "z=i" => \$z,
                "pep-data=s" => \$pep_data,
				"mass-type=s" => \$mass_type)
			){
                print "Help_page\n"; # to be replaced by generic help subroutine
                die "Error, no arguments passed to ProNTo::seqMasses::handleOpts\n";
        }
		
		if (!defined $mass_type || ($mass_type ne 'a' and $mass_type ne 'm') ) {
			die ("Please select either m for monoisotopic, or a for average masses");
		}
}

sub readPepFASTA{

	my @header = ("Protein", "Peptide Number", "Missed Cleavages", "Sequence", "m/z", "z");
	printf ("%-11s \t %-11s \t %-s \t %-s \t %-8s \t %s \n", @header);
	#The above two lines print the headers of each column for the table in correct format
	
	open(SEQFILE,$pep_data) or die "Unable to open file $pep_data\n";
	
	while ( <SEQFILE> ) {
		processSeqLine($_);
	}
	close SEQFILE or die "Unable to close file $pep_data\n";
}
	
# The subroutine is not specifically required for this program but it makes the while loop
# much easier to read and follow

sub processSeqLine {
	my $currentline = shift;
	chomp $currentline;
	if ( /^>/ ) {
		# add the header line in FASTA file to the table for each peptide
		my @first_column = substr((split(/|/,$currentline)),1);
		printf ("%-11s \t %-11s %-11s \t", @first_column);	  # and knocks off the >
	} elsif ( /^A-Z/ ) {
		# Obtain m/z values
		my $mz = MZconverter($currentline);
		my @seqdata = ($currentline, $mz, $z);
		printf ("%-s \t %-8s \t %s \n", @seqdata);
	} else {
		die "$pep_data in wrong format";
	} 
}

sub massWater {
	if ($mass_type eq 'a') { #average
		return 18.0106;
	} else { # monoisotopic
		return 18.0153;
	}
}

sub massTable {
	my $mass_type = shift;
	# This is for the user definition in the command line to which masses are used
	if ($mass_type eq 'a') {
		my %mass_table = ( # average masses
			   A =>  71.08, C => 103.14, D => 115.09, E => 129.12,
			   F => 147.18, G =>  57.05, H => 137.14, I => 113.16,
			   K => 128.17, L => 113.16, M => 131.19, N => 114.10,
			   P =>  97.12, Q => 128.13, R => 156.19, S =>  87.08,
			   T => 101.10, V =>  99.13, W => 186.21, Y => 163.18,
			   '\s' => 0.0, "*" => 0.0
		   );
		return %mass_table;
	} else { # monoisotopic
		my %mass_table = (
			   A =>  71.0371, C => 103.0092, D => 115.0269, E => 129.0426,
			   F => 147.0684, G =>  57.0215, H => 137.0589, I => 113.0841,
			   K => 128.0950, L => 113.0841, M => 131.0405, N => 114.0429,
			   P =>  97.0528, Q => 128.0586, R => 156.1011, S =>  87.0320,
			   T => 101.0477, V =>  99.0684, W => 186.0793, Y => 163.0633,
			   '\s' => 0.0, "*" => 0.0
		   );
		return %mass_table;
	}
}

sub MZconverter {
	# pass in current line as a parameter
	my $currentline = shift;
	my %mass_ref = massTable($mass_type);
	my $water_mass = massWater($mass_type);
	my @residues = split(//,$currentline); 			   #splits sequence data into parts
	my @residue_masses = map { $mass_ref{$_} } @residues; 	   #uses hashes to find masses
	my $total_residue_mass = join ("+", @residue_masses);
	return(($total_residue_mass + $water_mass)/$z);
}

# End of module evaluates to true

1;

__END__

# End of file evaluates to false

0;
