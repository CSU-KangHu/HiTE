#!/usr/bin/perl -w

#######################################################
# Author  :  Aurelie K
# date    :  v1.0, Feb 2013
# email   :  4urelie.k@gmail.com
# Purpose :  Read a repeat masker output file => write TEs in a bed file
######################################################

use strict;
use warnings;

my $usage = "\nUsage:
    perl scriptname.pl <file.out> <base>
	
	This script reads a repeat masker output file => write TEs in a bed file
	<base> = base0 or base1; ie if base0, +1 to all starts +. put everything in base 1\n\n";

my $file = shift or die "$usage";
my $base = shift or die "$usage";

######################################################
# write TEs in bed format
######################################################

#0		1		2	3	4		5		6		7			8	9		10			11	12	13	14		15
#319	31.4	6.6	4.2	chr1	7349611	7349813	-189845619	+	MIRb	SINE/MIR	42	253	-15	6103	*

open(TE, "<$file");
open(BED, ">$file.bed");
while(<TE>) {
	chomp(my $line = $_);
	unless ($line =~ /Score|score|SW|^#/){
		if ($line =~ /\w/){
			$line =~ s/^\s+//;
			my @line = split('\s+',$line);
			$line[8] =~ s/C/-/;
			my ($chr,$start,$end,$strand) = ($line[4],$line[5],$line[6],$line[8]);
			my $ID = $line[0];
			for (my $i=1; $i<=$#line;$i++) {
				$ID = $ID.";".$line[$i];
			}
			if	($base eq "base0") {
				$start=$start+1;
			}
			print BED "$chr\t$start\t$end\t$ID\t.\t$strand\n"; #with ID being the whole line => easy to acces to RMoutput
		}
	}	
}
close TE;
close BED;
exit;
