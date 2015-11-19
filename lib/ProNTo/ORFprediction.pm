#!/usr/bin/env perl

package ProNTo::ORFprediction;

#use diagnostics;
use strict;
use warnings;
use Getopt::Long;

use File::Spec; # for class methods related to file path parsing
use File::Basename; # for dirname function

# Initialising variables in scope, then assign some of them in BEGIN block

my ($sequence_input, $abslibfile, $libdirname, $input_fh, $help_opt, $debug_option, $summary_option, $verbose,
    $read_sequence, $seq_direction, $seq_is_reversed, $fasta_width, %codon_table, $translate_off, $FASTAheader,
    $overlapping_ORFs, $biggest_orf_only, @start_positions_list, @end_positions_list, @uncalled_positions_list, @orf_cli_list,
    @ORFendStrings, %ORF_lengths, %maxORF, %ORF_start_end_pairs, $match_string_length, $min_orf_size, @selected_rf);

BEGIN {
    # __FILE__ provides a relative path to the bin file
    # According to http://stackoverflow.com/a/90721/2668831 this is the most robust method
    $abslibfile = File::Spec->rel2abs(__FILE__);
    # use dirname x2 to go "up" 2 directories from bin file, then descend into the lib path:
    $libdirname = dirname($abslibfile);
    
    # set some default variable values
    $read_sequence = '';
    $FASTAheader = 'Sequence';
    $fasta_width = 60;
    $translate_off = 0; # 1 is translate off: i.e. leave ORFs as DNA in output
    $debug_option = 0; # 0 is off, i.e. just ORF summary to STDERR (and FASTA to STDOUT)
    $summary_option = 0; # default 0, switch on via --summary: provide a summary of ORFs to STDERR...
                           # (frames listed in decreasing order of total length with no. of ORFs & specifies longest)
    $verbose = 0; # gives user more feedback on progress (to STDERR), set with -V flag
    $seq_direction = '';
    @selected_rf = (1...6); # default to use all reading frames
    $biggest_orf_only = 0; # default 0 (off) - user can select to use only the reading frame with greatest total length
    $overlapping_ORFs = 1; # default 1: overlapping ORFs are OK (can start new ORF within boundaries of previous)
    $seq_is_reversed = 0; # 0 means not reversed, 1 means reversed
    %ORF_lengths = (1 => 0, 2 => 0, 3 => 0, 4 => 0, 5 => 0, 6 => 0); # keep track of total ORF length in each R.F.

    # if this goes in a subroutine it'll get regenerated again and again, so set up once at startup
    # Stop codons have been removed from table and supplied in an array @ORFendStrings below
    %codon_table = (
	"AAA" => "K", "AAC" => "N", "AAG" => "K", "AAT" => "N",
	"ACA" => "T", "ACC" => "T", "ACG" => "T", "ACT" => "T",
	"AGA" => "R", "AGC" => "S", "AGG" => "R", "AGT" => "S",
	"ATA" => "I", "ATC" => "I", "ATG" => "M", "ATT" => "I",
	"CAA" => "Q", "CAC" => "H", "CAG" => "Q", "CAT" => "H",
	"CCA" => "P", "CCC" => "P", "CCG" => "P", "CCT" => "P",
	"CGA" => "R", "CGC" => "R", "CGG" => "R", "CGT" => "R",	
	"CTA" => "L", "CTC" => "L", "CTG" => "L", "CTT" => "L",
	"GAA" => "E", "GAC" => "D", "GAG" => "E", "GAT" => "D",
	"GCA" => "A", "GCC" => "A", "GCG" => "A", "GCT" => "A",
	"GGA" => "G", "GGC" => "G", "GGG" => "G", "GGT" => "G",
	"GTA" => "V", "GTC" => "V", "GTG" => "V", "GTT" => "V",
	"TAA" => "*", "TAC" => "Y", "TAG" => "*", "TAT" => "Y",
	"TCA" => "S", "TCC" => "S", "TCG" => "S", "TCT" => "S",
	"TGA" => "*", "TGC" => "C", "TGG" => "W", "TGT" => "C",
	"TTA" => "L", "TTC" => "F", "TTG" => "L", "TTT" => "F" );
    
    # Codon match sequences are always 3 nt, regex match index is for start of a match,
    # so increment +2 to get the end nucleotide's index if getting end positions
    $match_string_length = 3;
    $min_orf_size = 150; # default 150 nucleotide minimum ORF size
};
use lib $libdirname;
use Methods;

__PACKAGE__->main(@ARGV) unless caller;

sub main {
    handleopts();
    reader();
    mapSequence();
    outputFASTA();
}

sub handleopts {
    GetOptions("help" => \$help_opt,
               "debug" => \$debug_option,
               "V" => \$verbose,
               "orf-list=i{1,6}" => \@orf_cli_list, # specify between 1 and 6 reading frames
               "use-longest-orf" => \$biggest_orf_only,
               "t" => \$translate_off,
               "summary" => \$summary_option,
               "V" => \$verbose,
               "min-orf=i" => \$min_orf_size,
               "sequence-file=s" => \$sequence_input);
    if (defined $help_opt and $help_opt eq 1) {
        print STDERR ProNTo::Methods::help();
        exit(1); # printed to STDERR but this was the behaviour desired by the user, so success exit code
    }
    if (defined $min_orf_size and $min_orf_size lt 0) {
        die ProNTo::Methods::help('BAD_MIN_ORF');
    }
    
    if ($biggest_orf_only == 1 and scalar @orf_cli_list > 0) {
        # can't use both
        die ProNTo::Methods::help('ORF_SELECT');
    }
    
    if (@orf_cli_list) {
        @selected_rf = @orf_cli_list;
    }
    
    
    if (!defined $sequence_input) {
            # if file not provided then show help message
            die ProNTo::Methods::help('NO_SEQ_FILE');
    }
    
}

sub reader {
    my $input_fh;
    open $input_fh, "<", $sequence_input
         or die ProNTo::Methods::help('BAD_SEQ_FILE');
    while (my $line = <$input_fh>) {
        # dot assignment operator would appear simpler, but Camel book says to prefer join
        # http://docstore.mik.ua/orelly/perl2/prog/ch24_02.htm
        # Convert sequence to upper case
        $read_sequence = join('', $read_sequence, processFASTA($line));
    }
    close $input_fh or die "$sequence_input: close failed: !";
    return;
}

sub mapSequence {
    mapORFs('forward');
    mapORFs('reverse');
    assessMap();
}

sub arrangeSequence {
    $seq_direction = shift;
    if ($seq_direction eq 'reverse' and $seq_is_reversed == 0) {
        # reverse the sequence - not already reversed
        $read_sequence = reverse $read_sequence;
        $seq_is_reversed = 1;
    } elsif ($seq_direction eq 'forward' and $seq_is_reversed == 1) {
        # reverse the sequence back to normal - it has been reversed (should not happen but just in case)
        $read_sequence = reverse $read_sequence;
        $seq_is_reversed = 0;
    }
}

sub mapORFs {
    $seq_direction = shift;
    arrangeSequence($seq_direction);
    
    # Could do all of this in a much simpler lookbehind regular expression, but that would be potentially memory intensive,
    # especially if the regex were to be expanded.
    #
    # Instead, this uses an algorithm to step through genomic coordinates of matches, finding the first 'stop' after each 'start'.
    
    #for my $rf_num (1...3) { processSomehow($read_sequence); }
    @start_positions_list = match_start_codons($read_sequence);
    @end_positions_list = match_end_codons($read_sequence);
    @uncalled_positions_list = match_non_acgt($read_sequence); # uncalled bases

    # rather than loop through the entire list of start/end positions, skip to a 'seen so far' position at this index:
    my $end_seen_index = my $start_seen_index = 0; # initialise at start of coordinate lists & no ORFs seen
    my %last_ORF_end_position = (1 => 0, 2 => 0, 3 => 0, 4 => 0, 5 => 0, 6 => 0);
    for my $start_position_index ($start_seen_index...@start_positions_list) {
        # print "Start at: $start_position\n";
        
        # The loop through @end_positions_list is about to start, keep track of how far into the list you pass until an ORF end
        my $end_index_advance_counter = 0;
                
        # loop from the index of the latest position seen to the end - don't waste loops or make new array subsets
        # for each index in the list of indices from 0 to N (where N is the number of positions in @end_positions_list, i.e. )
        for my $end_position_index ($end_seen_index...@end_positions_list - 1)  {
            my $start_position = $start_positions_list[$start_position_index];
            my $end_position = $end_positions_list[$end_position_index];
            my $ORFlength = $end_position - $start_position + 1;
            # use 4, 5, and 6 as the equivalent of reading frames 1, 2, and 3 in reverse direction
            my $reverse_frame_compensate = ($seq_is_reversed == 1) ? 3 : 0;
            my $reading_frame = ($start_position % 3) + 1 + $reverse_frame_compensate;
            if ($ORFlength % 3 == 0 # else OUT OF FRAME
                and
                $ORFlength >= $min_orf_size  # else TOO SHORT
                ) {
                if ($overlapping_ORFs == 1 # if default setting of "overlapping ORFs allowed" is being used...
                    or # ... or else if user has chosen to disallow ORFs within boundaries of other ORFs in same reading frame:
                        (   $overlapping_ORFs == 0
                        and
                        $start_position > $last_ORF_end_position{$reading_frame}
                        )
                    ) {
                    if (noUncalledBases($start_position, $end_position)
                        and
                        length (noInternalStopCodons($start_position, $end_position)) > 0) {
                        if ($debug_option == 1) {
                            printORFlogline($reading_frame, $ORFlength, $start_position_index, $start_position, $end_position_index, $end_position, ($last_ORF_end_position{$reading_frame}));
                        }
                        $ORF_lengths{$reading_frame} = $ORF_lengths{$reading_frame} + $ORFlength;
                        # print "New ORF total for $reading_frame:" . $ORF_lengths{$reading_frame} . "\n";
                        my @ORF_start_end_pair = ($start_position, $end_position);
                        push(@{ $ORF_start_end_pairs{$reading_frame} }, \@ORF_start_end_pair);
                        $last_ORF_end_position{$reading_frame} = $end_position;
                    }
                }
                # Start position was either assigned to an ORF pair, or was unsuitable for ORF-hood.
                # Now stop trying to match up this start position, advance to the next loop using it (within child loop so last not next)
                last;
            }
            # else just increment the $end_index_advance_counter and keep going
            
            # Not all end positions will be used, as they may occur in intergenic regions [i.e. after a STOP codon, before a START]
            # keep track of the index of the end position in the list of end positions:
            $end_index_advance_counter++;
        }
        # increment the index to start the coordinate comparison skipping past the ends already compared
        $end_seen_index = $end_seen_index + $end_index_advance_counter;
        # also record the actual position of this end, to skip past any start index less than the end just used
        
        # Don't reuse the one just matched up: increment the seen index to restart at the following end coordinate
        $end_seen_index++;
    }
}

sub printORFlogline {
    my ($log_rf, $orf_len, $startposi, $startpos, $endposi, $endpos, $lastorfendpos_logrf) = @_;
    print STDERR $log_rf . "|";
    print STDERR $orf_len . "nt.|Start: $startposi: $startpos\t";
    print STDERR "End: $endposi: $endpos ($lastorfendpos_logrf)\n";
    return;
}

sub noInternalStopCodons {
    # The check for minimum length has the side effect that internal stop codons may be allowed
    # Obviously these ORFs would not be viable if the substring to the first stop is below the threshold set
    my $startpos = shift;
    my $endpos = shift;
    my $ORF_string = substr($read_sequence, $startpos, $endpos - $startpos + 1) . "\n";
    return lengthCheckandTranslateFrame($ORF_string);
}

sub lengthCheckandTranslateFrame {
    my $checkstring = shift;
    my $peptide_string = TranslateFrame($checkstring);
    $peptide_string =~ /(^.+?)\*/;
    my $pep_to_stop = $1;
    my $min_codons_to_stop = ($min_orf_size / 3) - 1;
    if ((length($pep_to_stop)) > $min_codons_to_stop) {
        return $pep_to_stop;
    } else {
        return '';
    }
}

sub noUncalledBases {
    my $startpos = shift;
    my $endpos = shift;
    for my $uncalled_position (@uncalled_positions_list) {
        if ($uncalled_position > $startpos and $uncalled_position < $endpos) {
            if ($debug_option == 1) {
                print STDERR "Uh oh! $uncalled_position is between $startpos and $endpos - discard ORF candidate\n";
            }
            return 0; # False logical (failure), has uncalled bases - prevent use of ORF
        } else {
            return 1; # True logical (sucess), no uncalled bases - permit use of ORF
        }
    }
    # extension: could potentially look at position of base, given frame may be degenerate thus assignable to codon
}

sub assessMap {
    print STDERR "ORF summary:\n============\n" if $summary_option == 1;
    %maxORF = ('rf' => 'TBC',
               'val' => 'TBC');
    for my $orflength_rf (sort { $ORF_lengths{$b} <=> $ORF_lengths{$a} or $b cmp $a } keys %ORF_lengths) {
        if ($maxORF{'rf'} eq 'TBC') {
            # must be first iteration of loop
            $maxORF{'rf'} = $orflength_rf;
            $maxORF{'val'} = $ORF_lengths{$orflength_rf};
        }
        print STDERR "Frame " . $orflength_rf
            . ":" . $ORF_lengths{$orflength_rf}
            . " (" . scalar @{$ORF_start_end_pairs{$orflength_rf}}
            . " ORFs)\n"  if $summary_option == 1;
    }
    print STDERR "The largest is $maxORF{'val'} nucleotides, in reading frame $maxORF{'rf'} ("
        . scalar @{$ORF_start_end_pairs{$maxORF{'rf'}}}
        . " ORFs)\n"  if $summary_option == 1;
}

sub outputFASTA {
    if ($biggest_orf_only == 1) {
        @selected_rf = ($maxORF{'rf'}); # needed to wait until assessMap() ran to assign based on CLI flag setting
    }
    if ($verbose == 1){print scalar @selected_rf . " frame"; print "s" if scalar @selected_rf > 1; print " to use\n";};
    for my $output_rf (sort @selected_rf) {
        print STDERR " Generating FASTA for frame $output_rf\n" if ( $debug_option == 1 or $verbose == 1);
        print "\n" . generateFASTA($output_rf) . "\n";
    }
    return;
}

sub generateFASTA {
    my $fasta_rf = shift;
    my $fasta_counter = 0;
    if ($fasta_rf < 4) {
        arrangeSequence('forward');
    } else {
        arrangeSequence('reverse');
    }
    
    print STDERR "  There are " .
                  scalar @{$ORF_start_end_pairs{$fasta_rf}} .
                  " ORF start/end pairs for reading frame $fasta_rf\n" if ( $debug_option == 1 or $verbose == 1);
    for my $start_end_pair_index (0...scalar @{$ORF_start_end_pairs{$fasta_rf}} - 1) {
        my $start_end_pair_ref = @{$ORF_start_end_pairs{$fasta_rf}}[$start_end_pair_index];
        my @start_end_pair = @$start_end_pair_ref;
        my $start_coord = $start_end_pair[0];
        my $end_coord = $start_end_pair[1];
        my $ORF_string = substr($read_sequence, $start_coord, $end_coord - $start_coord + 1) . "\n";
        if ($translate_off == 1) { # for debugging
            print "\n$ORF_string\n";
            next;
        }
        my $peptide_string = lengthCheckandTranslateFrame($ORF_string);
        if (length($peptide_string) > 0) {
            my $peptide_FASTA = $FASTAheader . "|reading_frame:". $fasta_rf . "|pep_id:" . $fasta_counter . "\n"
            . join("\n", unpack("(A$fasta_width)*", $peptide_string))
            . "\n\n";
            print "$peptide_FASTA";
            $fasta_counter++;
        }
    }
}

sub match_start_codons {
    my $string = shift;
    my $regex = 'ATG';
    my @ret;
    # see comment at BEGIN block for $match_string_length
    while ($string =~ /$regex/g) {
        push @ret, $-[0];
   #     print $-[0] . "|" . substr($string , $-[0], 3) . "\n";
    }
    return @ret;
}

sub match_end_codons {
    # Pass in 2 parameters: regular expression and string to match on
    # generic function via http://stackoverflow.com/a/87504/2668831
    # ... plus 3rd parameter, whether providing index of start or end of matching string
    my $string = shift;
    # Ending in @ORFendStrings set as "TGA", "TAA", or "TAG"
    @ORFendStrings = qw(TGA TAA TAG);
    my $regex = join('|', @ORFendStrings);
    my @ret;
    # see comment at BEGIN block for $match_string_length
    while ($string =~ /$regex/g) {
        push @ret, $-[0] + 2;
    #    print $-[0] . "|" . substr($string , $-[0], 3) . "|" . substr($string, $-[0] + 2, 1) . "\n";
    }
    return @ret;
}

sub match_non_acgt {
    my $dna_string = shift;
    my @ret;
    while ($dna_string =~ /[^ACGT]/g) {
        push @ret, $-[0];
        # print $-[0] . "|" . substr($dna_string , $-[0], 1) . "\n";
    }
    return @ret;
}

sub processFASTA {
    my $line = shift;
    if ($line =~ /^>\s*(.+)/) {
        $FASTAheader = "> " . $1;
        print STDERR "Processing FASTA file ==$FASTAheader\n";
        return;
    } else {
        next if ($line =~ /^$/); # skip blank lines
        chomp($line);
        return uc "$line";
    }
}

sub TranslateFrame {
    my $framestring = shift;
    chomp $framestring;
    # use loaded genetic code
    my $transl = my $translcodon = '';
    # split frame (already confirmed multiple of 3) into codons and translate each
    my @framecodons = ( $framestring =~ m/[ACGT]{3}/g );
    for my $framecodon (@framecodons) {
        my $translcodon = $codon_table{$framecodon};;
        $transl = join('', $transl, $translcodon);
    }
    # concatenate each in @framecodons into a string to return, $transl
    return $transl;
}

# End of module evaluates to true

1;

__END__

# End of file evaluates to false

0;