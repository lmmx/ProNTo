#!/usr/bin perl

use strict;
use warnings;

# for is the same as foreach

  for my $protein (@proteins) {

  # regular expression to match either Lys or Arg unless next is Proline:
  # something like: [(K|R)^P]
  # i.e. either K or R, and then "not P"

  # $trypsin_regex = [(K|R)^P]

  # array of digested peptides is generated from the protein, by splitting on this regular expression
  # When there are no missed cleavages, the peptides are known as "limit peptides" (create in an array called @limit_peps)

  # this will currently only allow trypsin: make a subroutine to handle a generic enzyme and its regular expression

  if ($enzyme eq 'trypsin') {
    my @limit_peps = split($protein_sequence, $trypsinregex);
  } else {
    exit; # this is to remind you it only works with trypsin! Make a subroutine to handle others
  }
  # you now have a list of peptides
  
  # Keep a list of all the peptides, @all_peptides: made by joining "limit peptides" (i.e. with 0 missed cleavages), and the ones with up to "n missed cleavages"
  
  my @all_digested_peptides;
  my @missed_cleavage_peptides;

  for my $counter ($counter < $length_of_limit_peptides) {
    # loop through the A, B, C, D, E, and join them to simulate "1 missed cleavage" (or more, given by the user as $n)
    # use better variable names than n though!
    my $missed_cleavage_peptide = join @limit_peps[$counter,$counter+$n]
    # e.g. join @limit_peps[0,1] would give join ("A", "B"), which would give the first "missed cleavage" peptide, "AB"
    my $missed_cleavage_header = $description; # header line for sequence separated by pipes
    my $missed_cleavage_entry = $missed_cleavage_peptide . "\n" . $missed_cleavage_header; #concatenate the header line and the sequence
    push @missed_cleavage_peptides, $missed_cleavage_entry;
  }

  push @all_digested_peptides, (@limit_peptides, @missed_cleavage_peptides);


for $peptide (@all_digested_peptides) {
  print ">" . $peptide_header . "\n" . $peptide_sequence;
  # make a unique identifier for each peptide, which will allow accountability for each sequence
  # separate sections with pipe '|' symbols
  # this will generate FASTA format output to STDOUT for all peptides that have been stored in the @all_digested_peptides array
}

# NOTE - the FASTA format must have unique identifiers
