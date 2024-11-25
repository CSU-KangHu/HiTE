#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2018
# Creat a random DNA seq with your ATGC rate v1.0000 2018/03/08
# hukaining@gmail.com
#
use strict;
#use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
#use Parallel::ForkManager;
#our $MAX_processes=2;
#my $pm=Parallel::ForkManager->new($MAX_processes);
use re 'eval';
our $opfn="";
our $fasheader="randomDNArate1.1.1.1";
my $verbose;
our $dnal=100;
our $atgcrate="1,1,1,1";
#our $strict="yes";


GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"l=i"=>\$dnal,"n=s"=>\$fasheader,"r=s" => \$atgcrate)
or die("Error in command line arguments\nUsage: 	perl RandomDNA.pl [options] \n
	 options:\n
	 [-o string|output prefix default: rnddna$dnal]\n
	 [-n string|Add fasta header default: \"randomDNA\" ]\n
	 [-l int|Length of the DNA. default: $dnal]\n
	 [-r string|A,T,G,C rate, sep by comma. default:1,1,1,1]\n
	 Note: Creat a random DNA seq with your ATGC rate v1.0000 2018/03/08\n");

chomp $dnal;
if ($dnal<=0) {
	die ("[-] Error: Can't creat $dnal bp.\n");
}
#print $ARGV[0],"\n";
#if (not $ARGV[0]) {
	#die ("[-] Error: Not find a inputfile.\n");
#}
#my @tmpfliename=@ARGV;
if ($fasheader eq ""){
  $fasheader="randomDNArate$atgcrate";
}
if ($opfn eq ""){
	$opfn="rnddna$dnal.r$atgcrate";
	print "Output file: $opfn.fa \nFasta header: $fasheader\n";
}else{
	print "Output file: $opfn.fa \nFasta header: $fasheader\n";
}

open OUT, "> $opfn.fa" or die ("[-] Error: Can't open or creat $opfn.fa\n");
print OUT ">$fasheader\n";
our @ATGC=split(",",$atgcrate)or die ("[-] Error: Can't read the -r ATGC-rate. 4 numbers sep by 3 commas. \n");
our $ratesum=$ATGC[0]+$ATGC[1]+$ATGC[2]+$ATGC[3];
if ($ratesum<0) {
  die ("[-] Error: rate sum less than zero.\n");
}
our $Arate=$ATGC[0];
our $Trate=$ATGC[1];
our $Grate=$ATGC[2];
our $Crate=$ATGC[3];
our $Arange=$Arate;
our $Trange=$Arange+$Trate;
our $Grange=$Trange+$Grate;
our $Crange=$ratesum;
print "ATGC rate: A:$Arate T:$Trate G:$Grate C:$Crate Sum:$ratesum\n";



our $starttime=time();

#our $tmpseqname="";
our $count1=0;
our $count2=0;
#our $cutseq="";
our $finaldna="";
######main

while($count1<$dnal){
#    my(@nts)=qw(A T G C);
#    my $newnt=$nts[rand @nts];
#    #$finaldna=$finaldna.$newnt;
#    print OUT "$newnt";
#    $count1++;
     my $newnt="";
     given(rand($ratesum)){
      when (0<=$_ and $_<=$Arange){$newnt="A"}
      when ($Arange<$_ and$_<=$Trange){$newnt="T"}
      when ($Trange<$_ and $_<=$Grange){$newnt="G"}
      default {$newnt="C"}
     }



     
     print OUT "$newnt";
     $count1++;


} #while End.

print OUT "\n";
print "Finished Creating a $dnal bp random DNA sequence.\n";
close OUT;
our $endtime=time();
printf "Done! %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;
