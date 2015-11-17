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
    @mzdata,
    @sortedmzvals,
    $usage
    );

BEGIN {
        $data =  "../../data/mztable.tsv";
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
        @mzdata = reader();
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
                $usage = "Error, failed to obtain sufficient arguments to calculate ion statistics.";
                die $usage;
        }
}

sub reader {
        open(IFILE, "$data") or die "Error, could not open input file\n";

        # READING THE DATA
        
        my @lines = <IFILE>;
        my @MZvalues;
        foreach my $line (@lines) {
                my $MZvalue = LineToMZ($line, $AAfilterContains, $AAfilterTerm);
                push @MZvalues, $MZvalue;
        }
}

sub LineToMZ {
        # take the first parameter to the subroutine as a scalar variable, called line
        my $line = shift;
        my $AAfilterContains = shift;
        my $AAfilterTerm = shift;
        if ($AAfilterContains ne "None") {
                # check if contains, or doesn't contain, a given set of residues
                processAA($line, $AAfilterContains);
        }
        if ($AAfilterTerm ne "None") {
                # Specify the terminus and residues in the format:
                #       --aa-filter-term N:M
                # to indicate an N-terminal Met residue
                # or
                #       --aa-filter-term C:WQD
                # to indicate C-terminal W, Q, or D residues
                $AAfilterTerm =~ /(^.):(.+)/;
                # The above regular expression has 2 match groups:
                #       $1 is the terminus (C or N)
                #       $2 is one or more residues
                my $terminus = $1;
                my $term_residues = $2;
                processAA($line, $term_residues, $terminus);
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
        my $line = shift;
        my $num_params = scalar @_;
        if ($num_params == 2) {
                # using $line and $AAfilterContains
                $AAfilterContains = shift;
        } elsif ($num_params == 3){
                
        } elsif ($num_params == 5) {
                
        } else {
                die "Incorrect number of parameters passed to processAA";
        }
        next if ($line =~ /^$/);
        my @AAcolumn = split( /\t|\n/, $line );
        my $AAfilcolumn = $AAcolumn[2];
        if ( $AAfilcolumn =~ /$AAfilterContains/) {
                my $mz = $AAcolumn[3];
                my $AAmatch = $1;
                print "Here's one\n";
        }
}

sub generateBins {
        my @MZvalues = shift;
        my @sortedMZvalues = sort @MZvalues;
        print "\n@sortedMZvalues\n";
        return;

        my $bin_number = ($massmax - $massmin)/$binwidth;
        my $bin = int ($MZ - $massmin)/$binwidth;
        foreach $MZ (@MZvalues) {
                if ($MZ > $massmin and $MZ <= $massmax and $MZ <= $bin)
        }
}

# End of module evaluates to true

1;

__END__

# End of file evaluates to false

0;
