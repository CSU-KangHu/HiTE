#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2019
# Count Seq ATGCN (CountATGC)v1.2000 2019/05/27
# hukaining@gmail.com

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
#use Parallel::ForkManager;
#our $MAX_processes=2;
#my $pm=Parallel::ForkManager->new($MAX_processes);
use re 'eval';
our $opfn="";
my $verbose;
#our $upstreml=3000;
#our $downstreml=500;
#our $seqfilename ='';
#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn,"verbose"=>\$verbose)
or die("Error in command line arguments\nUsage: perl countATGC [-o outfileprefix] <inputfile>\n");
our $loadingstarttime=time();

print "Start loading genomeic sequence(s).\n";
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
#say $ARGV[0];
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
#close SEQFILENAME;
our $loadingendtime=time();
print "$Chri Sequences\n";
print "Finished loading!\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;

if ($opfn eq ""){
	$opfn="ATGCNout";
	print "Output file: ATout.txt \n";
}else{
	print "Output file: $opfn.txt \n";
}

open OUT, "> $opfn.txt" or die ("[-] Error: Can't open or creat $opfn.fa\n");
print OUT "Seq\tA\tT\tG\tC\tN\tAT\tGC\tAll\tATGC\tA\tT\tG\tC\tN\tAT\tGC\tATGC\n";
our $starttime=time();
#  our $countA=0;
#	our $countT=0;
#	our $countG=0;
#	our $countC=0;
#	our $oneseql=0;
#	our $oneseq="";
our $seqall="";
print"#################################\n#            Start              #\n#################################\n";
for (our $ni=0;$ni<$Chri;$ni++){
	#print "Deal with $Chrname[$ni].\n";
	our $countA=0;
	our $countT=0;
	our $countG=0;
	our $countC=0;
	our $countN=0;
	our $oneseql=0;
	our $oneseq = $Chrseq[$ni];
	our $countAT=0;
	our $countGC=0;
	our $countnN=0;
	$oneseq=~s/\s//;
   $oneseql=length($oneseq);
	 $countA = $oneseq=~ tr/A/A/;
	 $countT = $oneseq=~ tr/T/T/;
	 $countG = $oneseq=~ tr/G/G/;
	 $countC = $oneseq=~ tr/C/C/;
	 $countN = $oneseq=~ tr/N/N/;
	 $countnN = $countA + $countT + $countG +$countC;
	print "$Chrname[$ni]: A=$countA, T=$countT, G=$countG, C=$countC, N=$countN, All=$oneseql, AllnotN=$countnN\n";
	printf "$Chrname[$ni]: A=%g, T=%g, G=%g, C=%g, N=%g, AT=%g, GC=%g, ATGC=%g \n", $countA/$oneseql, $countT/$oneseql, $countG/$oneseql, $countC/$oneseql, $countN/$oneseql, ($countA+$countT)/$oneseql, ($countG+$countC)/$oneseql, $countnN/$oneseql;
  $countAT=$countA+$countT;
  $countGC=$countG+$countC;
  print OUT "$Chrname[$ni]\t$countA\t$countT\t$countG\t$countC\t$countN\t$countAT\t$countGC\t$oneseql\t$countnN";
  printf OUT "\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", $countA/$oneseql, $countT/$oneseql, $countG/$oneseql, $countC/$oneseql, $countN/$oneseql, ($countA+$countT)/$oneseql, ($countG+$countC)/$oneseql, $countnN/$oneseql;
   $seqall.=$oneseq;
}

our $countA=0;
	our $countT=0;
	our $countG=0;
	our $countC=0;
	our $countN=0;
	our $oneseql=0;
	our $oneseq = $seqall; #Whole genome.
	our $countAT=0;
	our $countGC=0;
	our $countnN=0;
	$oneseq=~s/\s//;
   $oneseql=length($oneseq);
	 $countA = $oneseq=~ tr/A/A/;
	 $countT = $oneseq=~ tr/T/T/;
	 $countG = $oneseq=~ tr/G/G/;
	 $countC = $oneseq=~ tr/C/C/;
	 $countN = $oneseq=~ tr/N/N/;
	 $countnN = $countA + $countT + $countG +$countC;
	print "All: A=$countA, T=$countT, G=$countG, C=$countC, N=$countN, All=$oneseql, AllnotN=$countnN\n";
	printf "All: A=%g, T=%g, G=%g, C=%g, N=%g, AT=%g, GC=%g, ATGC=%g\n", $countA/$oneseql, $countT/$oneseql, $countG/$oneseql, $countC/$oneseql, $countN/$oneseql, ($countA+$countT)/$oneseql, ($countG+$countC)/$oneseql, ($countnN)/$oneseql;
  $countAT=$countA+$countT;
  $countGC=$countG+$countC;
  print OUT "All\t$countA\t$countT\t$countG\t$countC\t$countN\t$countAT\t$countGC\t$oneseql\t$countnN";
  printf OUT "\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", $countA/$oneseql, $countT/$oneseql, $countG/$oneseql, $countC/$oneseql, $countN/$oneseql, ($countA+$countT)/$oneseql, ($countG+$countC)/$oneseql, $countnN/$oneseql;

	close OUT;
	our $endtime=time();

print"#################################\n#             End               #\n#################################\n";
  printf "Done! %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;
	#print OUT ""
