#!/usr/bin/env perl

package ProNTo::Methods;

use diagnostics;
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

# declare these to be assigned in the BEGIN block
my ($absbin, $libdirname);

# Initialising variables in scope, then assign some of them

my ($usage, $abslibfile, $libdirname);

BEGIN {
    # __FILE__ provides a relative path to the bin file
    # According to http://stackoverflow.com/a/90721/2668831 this is the most robust method
    $abslibfile = File::Spec->rel2abs(__FILE__);
    # use dirname x2 to go "up" 2 directories from bin file, then descend into the lib path:
    $libdirname = dirname($abslibfile);
    
    # set some default variable values here?
};

use lib $libdirname;

# This module exists only to provide/expose common methods to other modules

sub help {
    return "hello world i'm helping\n";
}

# End of module evaluates to true

1;

__END__

# End of file evaluates to false

0;