#!/usr/bin/env perl
use strict;
use warnings;

# usage: perl gff2RMout.pl genome.gff3 genome.out

# Open the GFF3 file for reading
open(my $gff3_file, "<", $ARGV[0]) or die "Cannot open input.gff3: $!";

# Open the RepeatMasker .out file for writing
open(my $rm_out_file, ">", $ARGV[1]) or die "Cannot open output.out: $!";

# Read and process the GFF3 file line by line
while (my $line = <$gff3_file>) {
    chomp $line;
    next if $line =~ /^#/;
    my @fields = split /\t/, $line;

    # Extract necessary information (adjust indices based on your GFF3 structure)
    my $seqName = $fields[0];
    my $source = $fields[1];
    my $featureType = $fields[2];
    my $start = $fields[3];
    my $end = $fields[4];
    my $div = $fields[5];
    my $strand = $fields[6];
    my $phase = $fields[7];
    my $attributes = $fields[8];
    my ($family, $classification) = ($1, $2) if $attributes =~ /Name=([^;]+);Classification=([^;]+)/i;

    # Format the output for RepeatMasker .out (modify according to the expected format)
    printf $rm_out_file "1000 0.001 0.001 0.001 $seqName $start $end NA $strand $family $classification\n";
}

# Close the filehandles
close($gff3_file);
close($rm_out_file);

