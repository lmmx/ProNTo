#!/usr/bin/env perl

package ProNTo::IonStats;

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
    $binvalue,
    $data,
    $AAfilterContains,
    $AAfilterTerm);
BEGIN {
        $data =  "mztable.tsv";
        $AAfilterContains = "None";
        $AAfilterTerm = "None";
};

__PACKAGE__->main(@ARGV) unless caller;

sub main {
        handleopts();
        reader();
        print "\n\nAll done\n";
}

sub handleopts {
        # Defining bin size and parameters
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

# End of module evaluates to true

1;

__END__

# End of file evaluates to false

0;