#!/bin/env perl

# This script takes annotations of transposable elements and converts it into a GTF file that is compatible with TEtranscripts.
# The input file requires the following information in a tab-separated file: genomic location of TE (chrom, start, end, strand) and TE name (element/subfamily, e.g. L1HS).
# For a clearer annotation in TEtranscripts results, we recommend also providing the TE family (e.g. L1) and class (e.g. LINE) information in the input file. Otherwise, the TE name would be used.

use strict;
use warnings;
use Getopt::Std;
use Carp;

# Parsing through the command line to ensure essential (and optional) parameters are set.

my ($chr,$start,$end,$strand,$name,$te_name,$te_fam,$te_class,$one_based,$score) = parse_cmd_line();

my $infile = shift @ARGV;

open my $ifh, "<", $infile or die $!;

# To account for multiple insertions of TE, this hash will keep track of the number of instances observed for each TE element/subfamily.
my %TE;

while(my $line = <$ifh>){
  next if $line =~ /^#/;   # Ignore all comment lines.
  chomp $line;
  # my @tmp = split "\t", $line;
  $_=$line;
  my @tmp = split;
  if($tmp[$strand] eq 'C'){$tmp[$strand]='-';}

# This step removes annotations that fall into repetitive "structural" RNA and low complexity genomic sequences. These annotations are commonly found in RepeatMasker output / UCSC rmsk tables.
  next if $tmp[$te_fam] =~ /(Low_complexity|Simple_repeat|rRNA|scRNA|snRNA|srpRNA|tRNA)/;

# This step ensures that the strand information is either "+" or "-". This is critical for TEtranscripts to assign read if RNA-seq library strand-specificity is invoked.
# If strand information is absent, the annotation is currently ignored.
  if($tmp[$strand] !~ /^(\+|-)/){
    print STDERR "Strand information on line $. ($tmp[$strand]) is not recognized. Only \"+\" & \"-\" are accepted\n";
    print STDERR "Line $. is skipped\n";
    next;
  }

# Adjust start position from 0-based indexing (e.g. BED or UCSC rmsk table) to 1-based indexing (GTF).
# Can be skipped if the input file uses 1-based indexing for genomic coordinates.
  $tmp[$start] ++ unless $one_based;

# This portion generates the name for the particular instance/insertion of the TE in the genome.
# It will check the TE hash to determine if the element/subfamily has been previously observed. If so, it will add the suffix "_dupX" to the particular instance, where X is the number of additional TE instances beyond the first. 
  my $te_instance;
  if(exists $TE{$tmp[$te_name]}){
    $te_instance = $tmp[$te_name] . "_dup$TE{$tmp[$te_name]}";
    $TE{$tmp[$te_name]} ++;
  }
  else{
    $te_instance = $tmp[$te_name];
    $TE{$tmp[$te_name]} = 1;
  }

# This step processes the TE family name to remove all "?" (common in UCSC rmsk annotations), and converts "Other" and "Unknown" to the TE subfamily name.
# However, it will retain the "Other" and "Unknown" TE class information in the output.  
  my $te_familyID = $tmp[$te_fam];
  my $te_classID = $tmp[$te_class];
  $te_familyID  =~ s/\?//;
  if($te_familyID =~ /(Other|Unknown)/){
    $te_classID = $te_familyID if $te_class == $te_fam;
    $te_familyID = $tmp[$te_name];
  }
  
# This generates the output following the defined GTF format (chrom, source, feature (exons in all cases), start position (1-based), end position, score ("." in all cases), strand, frame ("." in all cases) and attributes).
# TEtranscripts requires the following attributes in the TE GTF file: gene_id (TE subfamily/element name), transcript_id (TE instance name, e.g. L1HS_dup56), family_id (TE family name) and class_id (TE class name).
  print "$tmp[$chr]\t$name\texon\t$tmp[$start]\t$tmp[$end]";
  ($score < 0) ? (print "\t.") : (print "\t$tmp[$score]");
  print "\t$tmp[$strand]\t.\tgene_id \"$tmp[$te_name]\"; transcript_id \"$te_instance\"; family_id \"$te_familyID\"; class_id \"$te_classID\";\n";
}

# This subroutine parses the command line for arguments and checks that all parameters are set correctly.
sub parse_cmd_line{
  my %opts;
  getopts('c:s:e:o:n:t:f:C:S:1', \%opts);

# This portion checks for the presence of the essential parameters, and throws an error if they are absent.
  my $warnings;
  $warnings .= "No chromosome column provided\n" if !defined $opts{'c'};
  $warnings .= "No start position column provided\n" if !defined $opts{'s'};
  $warnings .= "No stop/end position column provided\n" if !defined $opts{'e'};
  $warnings .= "No strand column provided\n" if !defined $opts{'o'};
  $warnings .= "No column provided for TE name\n" if !defined $opts{'t'};
  $warnings .= "No input file provided\n" if scalar @ARGV < 1;
  if(defined $warnings){
    print STDERR "\n", $warnings;
    usage();
  }

# This portion checks that the columns provided are valid (i.e. numeric and greater than zero), and ensure that the columns for the start and end positions are not the same.  
  my $error = validate_column($opts{'c'}, "chromosome");
  $error .= validate_column($opts{'s'}, "start position");
  $error .= validate_column($opts{'e'}, "end position");
  $error .= validate_column($opts{'o'}, "strand");
  $error .= validate_column($opts{'t'}, "TE name");
  $error .= validate_column($opts{'f'}, "TE family") if defined $opts{'f'};
  $error .= validate_column($opts{'C'}, "TE class") if defined $opts{'C'};
  $error .= validate_column($opts{'S'}, "score") if defined $opts{'S'};
  if($error =~ /1/){
    usage();
  }
  my $chr = $opts{'c'} - 1;
  my $start = $opts{'s'} - 1;
  my $end = $opts{'e'} - 1;
  my $strand = $opts{'o'} - 1;
  if($start == $end){
    print STDERR "\nColumns for start position ($start) & end position ($end) cannot be the same\n";
    usage();
  }
  my $te_name = $opts{'t'} - 1;

# This portion looks for the presence of optional parameters.
# If columns for the TE family and/or class information is absent, it will switch to the default (using the subfamily/element name for the family, and family name for the class)
  my $name;
  (defined $opts{'n'}) ? ($name = $opts{'n'}) : ($name = "user_provided");
  my $te_fam;
  (defined $opts{'f'}) ? ($te_fam = $opts{'f'} - 1) : ($te_fam = $te_name);
  my $te_class;
  (defined $opts{'C'}) ? ($te_class = $opts{'C'} - 1) : ($te_class = $te_fam);
  my $one_based;
  (defined $opts{'1'}) ? ($one_based = 1) : ($one_based = 0);
  my $score;
  (defined $opts{'S'}) ? ($score = $opts{'S'} - 1) : ($score = -1);
  return ($chr,$start,$end,$strand,$name,$te_name,$te_fam,$te_class,$one_based, $score);
}

# This subroutine checks that the column id is valid (numeric & greater than 0).
sub validate_column{
  my ($num, $id) = shift @_;
  if($num < 1 || $num =~ /\D/){
    print STDERR "Column number for $id ($num) is invalid\n";
    return 1;
  }
  else{
    return 0;
  }
}  


sub usage{
  print <<EOF

 Usage: makeTEgtf.pl -c [chrom column] -s [start column] -e [stop/end column] 
                     -o [strand column] -n [source] -t [TE name column] 
                     (-f [TE family column] -C [TE class column] -1)
                     [INFILE]

 Output is printed to STDOUT

 Required parameters:
  -c [chrom column]     -    Column containing chromosome name
  -s [start column]     -    Column containing feature start position
  -e [stop/end column]  -    Column containing feature stop/end position
  -o [strand column]    -    Column containing strand information (+ or -)
  -t [TE name column]   -    Column containing TE name
  [INFILE]              -    File name to be processed into GTF

 Optional parameters:
  -n [source]           -    Source of the TE information 
                             (e.g. mm9_rmsk for RepeatMasker track from
                              mm9 mouse genome)
                             Defaults to "user-provided" if not specified
  -f [TE family column] -    Column containing TE family name. 
                             Defaults to TE name if not specified
  -C [TE class column]  -    Column containing TE class name. 
                             Defaults to TE family name if not specified
  -S [score column]     -    Column containing the score of the TE prediction
                             (e.g. score from RepeatMasker)
  -1                    -    Input coordinates uses 1-based indexing
                             This should be used if the input file uses
                             1-based coordinates. This should be invoked
                             if the genomic coordinates are obtained from
                             a GFF3, GTF, SAM or VCF file
                             Default: off if using BED, BAM or UCSC rmsk
                                      input files

EOF
    ;
  exit 1;
}
