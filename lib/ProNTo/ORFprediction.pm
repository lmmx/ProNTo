#!/usr/bin/env perl

package ProNTo::ORFprediction;

use diagnostics;
use strict;
use warnings;
use Getopt::Long;

use File::Spec; # for class methods related to file path parsing
use File::Basename; # for dirname function

# declare these to be assigned in the BEGIN block
my ($absbin, $libdirname);

use Exporter;

our @ISA = qw(Exporter);

# #################################
# TBC : which variables to export ?
# #################################
our @EXPORT_OK = qw();

# Initialising variables in scope, then assign some of them

my ($sequence_input, $usage, $abslibfile, $libdirname);

BEGIN {
    # __FILE__ provides a relative path to the bin file
    # According to http://stackoverflow.com/a/90721/2668831 this is the most robust method
    $abslibfile = File::Spec->rel2abs(__FILE__);
    # use dirname x2 to go "up" 2 directories from bin file, then descend into the lib path:
    $libdirname = dirname($abslibfile);
    
    # set some default variable values
    $sequence_input = '../../data/genome.fasta';
};
use lib $libdirname;

__PACKAGE__->main(@ARGV) unless caller;

sub main {
    $usage = ProNTo::Methods::help();
    &handleopts();
    &reader();
}

sub handleopts {
        if (!GetOptions (
                "sequence-file=s" => \$sequence_input)
            ){
                print "Help_page\n"; # to be replaced by generic help subroutine
                die $usage;
        }
}

sub reader {
    my $inputusage = "Error, could not open input file\n";
    open(IFILE, "$sequence_input") or die $inputusage;
    my (@outputarray, $outputvalue);
    my @lines = <IFILE>;
    foreach my $line (@lines) {
            my $newvariable = processLine($line);
            push @outputarray, $outputvalue;
    }
}

sub processLine {
        # take the first parameter to the subroutine as a scalar variable, called line
        my $line = shift;
        next if ($line =~ /^$/);
        my @columns = split( /\t|\n/, $line );
        # splits the row into columns

        my $mzcol = $columns[2];
        if ( $mzcol =~ /(\d+)/ ) {
                my $mzcolmatch = $1;
                # capture m/z value for the current row in table
                print "$mzcolmatch\n";
        }
}


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