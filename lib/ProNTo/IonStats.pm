#!/usr/bin/env perl

package ProNTo::IonStats;

use diagnostics;
use strict;
use warnings;
use Getopt::Long;

use Exporter;
our @ISA = qw(Exporter);

# #################################
# TBC : which variables to export ?
# #################################
our @EXPORT_OK = qw();

# Initialising variables in scope, then assign some of them

my ($binwidth,
    $massmin,
    $massmax,
    $data,
    $AAfilterContains,
    $AAfilterTerm,
    @sortedmzvals
    );

BEGIN {
        $data =  "mztable.tsv";
        $AAfilterContains = "None";
        $AAfilterTerm = "None";
        $binwidth = 10;
        $massmin = 1000;
        $massmax = 1500;
};

__PACKAGE__->main(@ARGV) unless caller;

sub main {
        handleopts();
        print "Bin width is $binwidth m/z; range is $massmin to $massmax m/z units.\n";
        reader();
        print "\nAll done\n";
}

sub handleopts {
        # Defining bin size and parameters
        if (!GetOptions (
                "bin-width=f" => \$binwidth,
                "file=s" => \$data,
                "bin-min=f" => \$massmin,
                "bin-max=f" => \$massmax,
                "aa-filter-contains=s" => \$AAfilterContains,
                "aa-filter-term=s" => \$AAfilterTerm)){
                print "Help_page\n"; # to be replaced by generic help subroutine
                die "Error, failed to obtain information from user\n";
        }
}

sub reader {
        open(IFILE, "../../data/$data") or die "Error, could not open input file\n";

        # READING THE DATA
        
        my @lines = <IFILE>;
        my @MZvalues;
        foreach my $line (@lines) {
                my $MZvalue = processMzValue($line, $AAfilterContains, $AAfilterTerm);
                push @MZvalues, $MZvalue;
        }
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

        my @columns = split( /\t|\n/, $line );
        # splits the row into columns

        my $mzcol = $columns[2];
        if ( $mzcol =~ /(\d+)/ ) {
                my $mzcolmatch = $1;
                # capture m/z value for the current row in table
                print "$mzcolmatch\n";
        }
}

sub processAA {
        # Filter for amino acid composition/terminal residue to generate stats
}

sub generateBins {
        my @MZvalues = shift;
        my @sortedMZvalues = sort @MZvalues;
        print "\n@sortedMZvalues\n";
        return;
}

# End of module evaluates to true

1;

__END__

# End of file evaluates to false

0;