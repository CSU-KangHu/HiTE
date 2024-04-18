#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2018
# Sliding Window Count (SWcount)v2.1000 2018/07/18
# hukaining@gmail.com
#
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use re 'eval';

our $opfn="";
my $verbose;
our $step = 500000;
our $window = 1000000;
our $Chrlenlist = "Chrlength.txt";
#our $Fline= -1;

GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"s=i"=>\$step, "w=i"=>\$window, "l=s"=> \$Chrlenlist)
	or die("Error in command line arguments\nUsage: 
	\nperl WSAcount.pl [options] -l <Chromesome Length List(tab)> <inputFSIout>\n
	options:\n
	 [-o outprefix default: filtered.wsa.out]\n
	 [-s int|Step size default: 500000]\n
	 [-w int|window size default: 1000000]\n
	 Note: Sliding Window Count (SWcount)v2.1000 2018/07/18 \n");
	#open IN,"chrA09.snp.vcf";

print "Step Size: $step\nWindow Size: $window\n";
if ($step>$window) {
	die("[-] Error: -s <Step Size> must <= -w <Window size> .\n");
}

open LEN, "$Chrlenlist" or die ("[-] Error: Can't open the Chromesome Length List.");

if ($opfn eq ""){
$opfn="filtered.wsa.out";
print "Output file:$opfn.txt \n";
}else{
print "Output file:$opfn.txt \n";
}

open OUT, "> $opfn.count.txt" or die ("[-] Error: Can't open or creat $opfn.count.txt\n");
print OUT "#CHROM\tPOS_Start\tPOS_End\tPOS_Mid\tCounts\n";


our %chrl;
our @chrs=();
while (my $chrinfo=<LEN>) {
	 chomp ($chrinfo);
	 if ($chrinfo =~ m/^\#/) {next;}
	 my @chrlen = split (/\t/,$chrinfo);
	 push @chrs, $chrlen[0];
	 $chrl{$chrlen[0]}=$chrlen[1];
}

our @allsnp=();
while (defined(our $brow=<>)){
	chomp ($brow);
	push @allsnp, $brow;
	
}


############
##Main
############
foreach our $chrid (@chrs){
	our $len=$chrl{$chrid};
	print "Deal with $chrid \n";
	#my $count1=0;
	
	our @snps=();
	
	foreach my $row (@allsnp) {
	#while (my $row=<>) {
		chomp $row;
		my @FSI = split (/\t/,$row);
		if ($FSI[0] eq $chrid ){
			push @snps, $FSI[4];
		}else {
			next;
		}
		}
		
		if (scalar @snps ==0) {
			print "[-] Not load $chrid pos.\n";
			next;
		}
		
	
		for (my $posstart=0;$posstart <= $len;$posstart+=$step){
			my $posend = $posstart + $window;
			my $posmid = ($posstart+$posend)/2;
			my $count1=0;
			foreach my $snppos (@snps){
				if ($snppos>=$posstart and $snppos<=$posend){$count1++;}
				if ($snppos>$posend){
					print OUT "$chrid\t$posstart\t$posend\t$posmid\t$count1\n";
					last;
				}
			}
		}
		print "[+] $chrid done.\n";
		
}

close OUT;
print "All done.\n"
