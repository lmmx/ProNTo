#!/usr/bin/perl

use diagnostics;
use strict;
use warnings;
use Getopt::Long;

# Defining bin size and parameters

print "Enter file name:\n";
my $data = <>;
#my $AAfilterContains = <>;
print "Enter m/z minimum range limit:\n";
my $massmin = <>;
print "Enter m/z maximum range limit:\n";
my $massmax = <>;
print "Enter binwidth:\n";
my $binwidth = <>;
my @sorted;
my $mzcolmatch;

handleopts();

sub handleopts {
        if (!GetOptions (
                "bin-width=i" => \$binwidth,
                "file=s" => \$data,
                "mass-min=i" => \$massmin,
                "mass-max=i" => \$massmax,
#                "aa-filter-contains=s" => \$AAfilterContains
                )) {
                print STDERR "Error, failed to obtain information from user\n";
                print "Help_page\n";
                exit(0);
        }
}
open(IFILE, $data) or die "Error, could not open input file\n";

my @lines = <IFILE>;
my @MZvalues;
my $AA;
my $line;
foreach my $line (@lines) {
#        if ($AAfilterContains =~ /[a-zA-z]/) {
#                my $MZvalue = processAA($line);
#        } else {
                my $MZvalue = processMZValue($line);
#        }
}
sub processMZValue {
        my $line = shift;
        next if ($line =~ /^$/);
        # this error-proofs the read so that at every loop it will ig$
        my @column = split( /\s+/, $line );
        # splits the column
        my $mzcol = $column[2];
        if ( $mzcol =~ /(\d+)/ ) {
                my $mzcolmatch = $1;
                # capture m/z value for the current row in table
                push @MZvalues, $mzcolmatch;
        }
}

#Failed attempt at screening for peptides with certain AAs;
#sub processAA {
#        my $line = shift;
#       chomp($AAfilterContains);
#                next if ($line =~ /^$/);
#                my @AAcolumn = split( /\t|\n/, $line);
#                my $AAfilcolumn = $AAcolumn[2];
#                if ( $AAfilcolumn =~ /$AAfilterContains/) {
#                       next if ($line =~ /^$/);
#                       my @column = split( /\t|\n/, $line );
#                       my $mzcol = $column[3];
#                       my $mzcolmatch = $1;
#                       print "$mzcolmatch\n" or die "Nope";
#                       push @MZvalues, $mzcolmatch;

#                }
#}
print "$data was the file used.\n The binwidth is $binwidth.\n";
generateBins();
sub generateBins {
        my $i;
        my @finalCount;
        my $MZvalue;
        my $binline;
        my $binmin = $massmin;
        my $binmax = $massmin + $binwidth;
        my $num_bins = ($massmax - $massmin)/$binwidth;
        my $bin_id = 0; # start at bin 0
        # start loop with $binmin = $massmin and $binmax = $massmin + $binwidth
        for $MZvalue (@MZvalues) {
#               my @binned_mz_values;

#               SORT:
#                if ($MZvalue >= $binmin and $MZvalue <= $binmax) {

#                        addToCurrentBin($mzValue);
#                } else {
#                       while ($MZvalue > $binmax) {
#                               incrementBinBoundaries(); # increase $binmin and $binmax and return, exit if $binm$
#                               goto SORT;
#                       }
                if ($MZvalue >= $binmin and $MZvalue <= $massmax) {
                push @sorted, (int(($MZvalue-$massmin)/$binwidth));

                }

        }
        my @sortSorted = sort @sorted;
        for ($i = 0; $i < $num_bins; $i++) {
                $finalCount[$i]=0;
                foreach $binline (@sortSorted) {
                        if ($binline == $i) {
                                $finalCount[$i]++;
                        }
                }
                my $printA = $massmin + $binwidth * $i;
                my $printB = $massmin + $binwidth * ($i+1);
                print "BIN $printA\-$printB: $finalCount[$i]\n";
        }
}

# Failed attempt at sliding window approach

#sub addToCurrentBin {
#       push @bin_id ($MZvalue);
#       push @{sorted{}{}{}}, $MZvalue;
#}

#sub incrementBinBoundaries {
#       if ($binmax > $massmax) { exit; }
#       else {
#       $binmin = $binmin + $binwidth;
#       $binmax = $binmin + $binwidth;
#       $bin_id++;
#       $sorted{}->{"bin_id}++;
#       }
#}

#foreach my $sortedBins (sort keys %sorted) {
#       print scalar($sorted{}{});

exit;