#!/usr/bin/env perl

package ProNTo::ORFprediction;

use strict;
use warnings;
use Getopt::Long;

use Exporter;
our @ISA = qw(Exporter);

# #################################
# TBC : which variables to export ?
# #################################
our @EXPORT_OK = qw();

__PACKAGE__->main(@ARGV) unless caller;

# Read in command line args

my $fh = shift;

# if file not provided then show help message

# read FASTA DNA nucleotide sequence into array of text strings
# stripping off the ">" at the start of the line
my @FASTAseqs;

# Find all possible ORFs beginning with one of @ORFstartStrings, set as just "ATG"
my @ORFstartStrings = qw(ATG);
# Ending in @ORFendStrings set as "TGA", "TAA", or "TAG"
my @ORFendStrings = qw(TGA TAA TAG);

# find all bounding coordinates
my @ORFstartCoords;
my @ORFendCoords;
my $prevEnd = 0;

sub RetrieveORF {
    my $startCoord = shift; # this is where to begin from
    my $skipPast = shift; # the previous position seen: use next end coord in @ORFendCoords
    # ...
    # when finished bump $prevEnd:
    #  - get index of end coordinate in @ORFendCoords,
    #  - increment (++)
    #  - get coordinate at this index, set as new $prevEnd
}

my @ORFpreds;

for my $putativeStart (@ORFstartCoords) {
    my $retrievedORF = RetrieveORF($putativeStart, $prevEnd);
    push @ORFpreds, $retrievedORF;
}

# finish with a list of predicted ORFs, @ORFpreds

sub TranslateFrame {
    my $framestring = shift;
    # use loaded genetic code
    my $transl;
    my @splitframe = ( $framestring =~ m/.../g );
    for my $codon (@splitframe) {
        
    }
    # concatenate each in @splitframe to $transl
    return $transl;
}

for my $ORFpred (@ORFpreds) {
	TranslateFrame($ORFpred)
}

# End of module evaluates to true

1;

__END__

# End of file evaluates to false

0;