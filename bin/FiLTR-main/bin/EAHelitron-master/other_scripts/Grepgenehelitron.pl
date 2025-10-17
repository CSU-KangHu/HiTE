#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2018
# Grep gene/helitron in GFF/GTF for silding window analysis. v1.0000 2018/07/18
# hukaining@gmail.com
#
use strict;
use warnings;
#use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
#use Parallel::ForkManager;
#our $MAX_processes=2;
#my $pm=Parallel::ForkManager->new($MAX_processes);
use re 'eval';
our $opfn="";
our $ipfn="";
my $verbose;
our $refgtf="";
our $sortID="gene_id";
our $feature="gene";

#our $seqfilename ='';
#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"i=s"=>\$ipfn,"T"=>\$refgtf,"f=s"=>\$feature)
or die("Error in command line arguments\nUsage:\n	 perl Grepgenehelitron [options] -i <input GFF>\n
	 options:\n
	 [-o outprefix default: SWLinput]\n
	 [-T input file is ref genome GTF format]\n
	 [-f string|Specify feature type in GTF annotation.default: gene]\n
	 Note: Grep gene/helitron in GFF/GTF for silding window analysis. v1.0000 2018/07/18\n");

open IN, "$ipfn" or die ("[-] Error: Can't open the EAHelitron output GFF File or ref-genome GTF.");
# open REFGTF, "$refgtf" or die ("[-] Error: Can't open the RefGenome GTF File.");


if ($opfn eq ""){
	$opfn="SWLinput";
	print "Output file:$opfn.txt \n";
}else{
	print "Output file:$opfn.txt \n";
}

open OUT, "> $opfn.txt" or die ("[-] Error: Can't open or creat $opfn.txt\n");
print OUT "Chr\tFeature\tStart\tEnd\tMidPos\n";
#our $loadingstarttime=time();
#say "Start loading genomeic sequence.";
#print "Feature: $feature\n";
#print "SortID: $sortID\n";
our $Hcount=0;
#our @Chrname=();
#our @Chrseq=();
#our (%HChr,%Hpm,%Hstart,%Hend);
#our @HIDS=();

our $loadingstarttime=time();
if ($refgtf eq "") {
  print "Start loading input GFF.\n";
  print "Feature: H\\d\+\\.3 \n";
  while(defined(our $inrow = <IN>)){
#while(our $seq = <SEQFILENAME>){
	if ($inrow =~ m/^\#/) {next;}
	our @tmp = split (/\t/,$inrow);
	if ($tmp[8] =~ m/H\d+\.3/) {
		print OUT "$tmp[0]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t",($tmp[3]+$tmp[4])/2,"\n" ;
		$Hcount++;
		}else{
		  next;
		}
		
	}
}else {
  print "-T option on. Start loading input GTF.\n";
  print "Feature: $feature\n";
  while(defined(our $inrow = <IN>)){
#while(our $seq = <SEQFILENAME>){
	if ($inrow =~ m/^\#/) {next;}
	our @tmp = split (/\t/,$inrow);
	if ($tmp[2] =~ m/$feature/) {
		print OUT "$tmp[0]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t",($tmp[3]+$tmp[4])/2,"\n" ;
		$Hcount++;
		}else{
		  next;
		}
		
	}
}


close IN;
close OUT;
our $loadingendtime=time();

print "Done. $Hcount lines.\n";
printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;