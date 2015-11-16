#!/usr/bin/env perl

package ProNTo::Template;

use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter);

# #################################
# TBC : which variables to export ?
# #################################
our @EXPORT_OK = qw();

=head1 NAME
ProNTo::Template - Testing out module usage
=head1 VERSION
Version 1.0
=cut

our $VERSION = '1.0';

=head1 SYNOPSIS
C<ProNTo::Template>, short for Perl Proteomics Pipeline,
is a template Perl module , within a pipeline used for
processing genomic sequence to produce a statistical
report for proteomics analyses.
Features include:
=over 4
=item * List the main features 
=item * Then
=item * List the "bonus" features
=back

PROteomics eNzymology TOolkit (ProNTo)

ProNTo: interpretation of genome sequence to asssess suitability of
digestion enzymes for proteomic mass spectral analysis.
    use ProNTo;
    # ...
    
ProNTo is well suited for use in the laboratory, either with a known enzyme,
or to assess among different enzymes.
    # ...

=head1 IMPORTANT LINKS
=over 4
=item * L<https://github.com/lmmx/p-p-pipeline/issues>
The queue for bugs & enhancements in ProNTo.
=back
=cut

print "Hello world";

# End of module evaluates to true

1;

__END__

# End of file evaluates to false

0;