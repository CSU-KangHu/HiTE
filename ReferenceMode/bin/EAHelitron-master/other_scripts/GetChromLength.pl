#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2018
# Get chromosome Length v1.0000 2018/07/19
# hukaining@gmail.com

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';

use re 'eval';
our $opfn="";
my $verbose;
#our $upstreml=3000;
#our $downstreml=500;
#our $seqfilename ='';
#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn,"verbose"=>\$verbose)
or die("Error in command line arguments\nUsage: perl GetChromeLegth.pl [-o outfileprefix] <input_fasta_file>\n");

if ($opfn eq ""){
$opfn="Out.ChromL";
print "Output file:$opfn.txt \n";
}else{
print "Output file:$opfn.txt \n";
}

open OUT, "> $opfn.txt" or die ("[-] Error: Can't open or creat $opfn.txt\n");


our $loadingstarttime=time();

print "Start loading genomeic sequence.\n";
#if(!open SEQFILENAME,"< $seqfilename"){
# die "File not Found\n";
#}
#if(!open SEQFILENAME,"< $_"){
# die "File not Found\n";
#}
our $Chri=0;
our @Chrname=();
our @Chrseq=();
#@ARGV = qw#''  Not_Find_a_File#;
#say @ARGV;
#say $0;
while(defined(our $seq = <>)){
#while(our $seq = <SEQFILENAME>){
	if ($seq =~ m/^.*>/) {
	$seq=~ m/^.*>([a-zA-Z0-9_.-]+) ?.*/;
	print "$1\n";
	 $Chrname[$Chri]= $1;
	$Chri++;
	}else{
		$seq =~ s/\s//;
		$seq =~ tr/MRWSYKVHDBmrwsykvhdb/CGTGCGGCGGCGTGCGGCGG/;#1.31add snp replace
	 $Chrseq[$Chri-1] .=$seq;
			}
}
for (our $ni=0;$ni<$Chri;$ni++){
	our $Chrnamelist =$Chrname[$ni];
  our $Chrl= length($Chrseq[$ni]);
  print OUT "$Chrnamelist\t$Chrl\n";
}

our $loadingendtime=time();
print "$Chri Sequences\n";
print "Finished loading!\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;