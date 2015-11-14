#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

# Defining bin size and parameters

my $data =  "table.tsv";
my $AAfilterContains = "None";
my $AAfilterTerm = "None";
my ($binwidth, $binvalue);

handleopts();

sub handleopts {
        if (!GetOptions (
                "bin-width=i" => \$binwidth,
                "file=s" => \$data,
                "bin-value=i" => \$binvalue,
                "aa-filter-contains=s" => \$AAfilterContains,
                "aa-filter-term=s" => \$AAfilterTerm)){
                print STDERR "Error, failed to obtain information from user\n";
                print "Help_page\n";
                exit(0);
        }
}

open(IFILE, $data) or die "Error, could not open input file\n";

# READING THE DATA

my @lines = <IFILE>;
my @MZvalues;
foreach my $line (@lines) {
        my $MZvalue = processMzValue($line, $AAfilterContains, $AAfilterTerm);
        push @MZvalues, $MZvalue;
}
sub processMzValue {
        # take the first parameter to the subroutine as a scalar variable, called line
        my $line = shift;
        my $AAfilterContains = shift;
        my $AAfilterTerm = shift;
        if ($AAfilterContains ne "None") {processAA($AAfilterContains);}
        if ($AAfilterTerm ne "None") {

                # Specify the terminus and residues in the format:
                #       --aa-filter-term N:RKH
                # to indicate N-terminal R, K, or H residues
                # or
                #       --aa-filter-term C:WQD
                # to indicate C-terminal W, Q, or D residues
                my $terminus =~ /(^.):(.+)/;
                # The above regular expression has 2 match groups:
                #       $1 is the terminus (C or N)
                #       $2 is one or more residues

                processAA($2, $1);
        }

        next if ($line =~ /^$/);
        # this error-proofs the read so that at every loop it will ignore the spaces

        my @column = split( /\t|\n/, $line );
        # splits the column

        my $mzcol = $column[2];
        if ( $mzcol =~ /(\d+)/ ) {
                my $mzcolmatch = $1;
                # capture m/z value for the current row in table

                print "$mzcolmatch\n" or die "Nope";
        }
}

sub processAA {
        # Filter for amino acid composition/terminal residue to generate stats
}

print "\n\nAll done\n";
exit;
