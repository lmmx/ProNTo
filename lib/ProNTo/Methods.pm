#!/usr/bin/env perl

package ProNTo::Methods;

#use diagnostics;
use strict;
use warnings;
use Getopt::Long;

use File::Spec; # for class methods related to file path parsing
use File::Basename; # for dirname function

use Exporter;

our @ISA = qw(Exporter);

# #################################
# TBC : which variables to export ?
# #################################
our @EXPORT_OK = qw(help);

# Initialising variables in scope, then assign some of them

my ($abslibfile, $libdirname);

BEGIN {
    # __FILE__ provides a relative path to the bin file
    # According to http://stackoverflow.com/a/90721/2668831 this is the most robust method
    $abslibfile = File::Spec->rel2abs(__FILE__);
    # use dirname to go "up" a directory from lib file to the folder containing other lib files, then descend into the lib path:
    $libdirname = dirname($abslibfile);
    
    # set some default variable values here?
};

use lib $libdirname;

# This module exists only to provide/expose common methods to other modules
# It has no main subroutine, nor a PACKAGE line specifying any default function call.

sub helpMessages {
    my $helptype = shift;
    my %helpMessageHash = (
        'BAD_SEQ_FILE' => 'Error, could not open input file.',
        'NO_SEQ_FILE' => 'Error, could not open input file: no file specified.',
        'NO_ENZ' => 'No enzyme selected',
        'BAD_MIN_ORF' => 'Error, minimum ORF size must be a non-negative integer.',
        'ORF_SELECT' => 'Error, please supply either a comma-separated list of reading frames to use (1 to 6) or set the use-longest-orf flag'
        );
    if (!defined $helptype) {
        # give back undefined if called without an argument
        return undef;
    } elsif (exists $helpMessageHash{$helptype}) {
        # give back help message if possible
        return $helpMessageHash{$helptype};
    } else {
        # give back undefined if no such help message
        return undef;
    }
}

sub usageStatement {
	my @helpstrings = (
	"Protein Digestion Options",
	"=====================================",
	"Options:",
    "    [-t]: digests protein with trypsin",
    "    [-l]: digests protein with endoproteinase Lys-C",
    "    [-a]: digests protein with endoproteinase Arg-C",
    "    [-v]: digests protein with V8 proteinase",
    "    [-c]: number of allowed missed cleavages. Default value is 0.",
    "    [-h]: prints help",
    ""
	);
	my $helpstr = join("\n", @helpstrings);
	return $helpstr;
}

sub help {
    my $helptype = shift;
    my $usagehelp = usageStatement; # provide generic usage guide from subroutine in this module
    my $helpmessage;
    if (!defined $helptype) {
        # Provide generic help message for unparametrised call to help()
        return $usagehelp;
    } else {
        # generic handling for a hash of potential help types in another subroutine
        $helpmessage = helpMessages($helptype);
    }
    
    # if the call to helpMessages produced a message, prepend to usage output
    if (defined $helpmessage) {
        $helpmessage = join("\n", ($helpmessage, $usagehelp));
    } else {
        # Unknown option to the help function, ignore it (shouldn't happen - use for debugging!)
        $helpmessage = $usagehelp . "\nWARNING! Unknown internal parameter to helpMessages(): $helptype";
    }
    return $helpmessage;
}

# End of module evaluates to true

1;

__END__

# End of file evaluates to false

0;