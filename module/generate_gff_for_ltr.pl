#!/usr/bin/env perl

##Generate gff3 file from pass list of LTR_retriever
##Usage: perl generate_gff_for_ltr.pl LTR.pass.list


use warnings;
use strict;

my $usage = "\n\tperl generate_gff_for_ltr.pl LTR.pass.list\n\n";
open List, "<$ARGV[0]" or die "ERROR: $usage";
open GFF, ">$ARGV[0].gff3" or die "ERROR: $!";
open my $fh, '<', $ARGV[1] or die "ERROR: $usage";

# 创建一个空的哈希
my %dict;

# 读取文件内容并存储到哈希中
while (my $line = <$fh>) {
    chomp $line; # 去除换行符
    my ($key, $value) = split /\t/, $line; # 按照 \t 分割行
    $dict{$key} = $value; # 存储到哈希中
}
# 关闭文件句柄
close $fh;

my $date=`date -u`;
chomp ($date);
print GFF "##gff-version 3\n##date $date
##ltr_identity: Sequence identity (0-1) between the left and right LTR region.
##tsd: target site duplication\n";

my $annotator="HiTE";
my $i=1;
while (<List>){
	next if /^#/;
	next if /^\s+/;
	my ($chr, $element_start, $element_end, $element_length, $lLTR_start, $lLTR_end, $lLTR_length, $rLTR_start, $rLTR_end, $rLTR_length, $seq_ID, $loc, $motif, $TSD, $lTSD, $rTSD, $IN, $sim, $strand, $supfam);
	($loc, undef, $motif, $TSD, $lTSD, $rTSD, $IN, $sim, $strand, $supfam) = (split);
#chr03:19323647..19328532        pass    motif:TGCA      TSD:GTCGC       19323642..19323646      19328533..19328537      IN:19323854..19328325   99.03
	($chr, $lLTR_start, $rLTR_end)=($1, $2, $3) if $loc=~/^(.*):([0-9]+)..([0-9]+)$/;
	($lLTR_start, $rLTR_end)=($rLTR_end, $lLTR_start) if $rLTR_end<$lLTR_start;
	($lLTR_end, $rLTR_start)=($1-1, $2+1) if $IN=~/^IN:([0-9]+)..([0-9]+)$/;
	$motif=~s/motif://gi;
	$TSD=~s/TSD://gi;
	$lTSD=~s/\.\./\t/;
	$rTSD=~s/\.\./\t/;
	$element_start=(split /\s+/, $lTSD)[0];
	$element_end=(split /\s+/, $rTSD)[1];
	# 根据哈希中的内容替换 chr_name
    $chr = $dict{$chr};
	my $id = "$chr:$lLTR_start..$rLTR_end";
	my $so = "LTR";
#	$so = "Copia_LTR_retrotransposon" if $supfam eq "Copia";
#	$so = "Gypsy_LTR_retrotransposon" if $supfam eq "Gypsy";
	if ($TSD eq "NA"){
		$element_start=$lLTR_start;
		$element_end=$rLTR_end;
		}
	my $chr_ori=$chr;
#	my $info = "Name=$id;motif=$motif;tsd=$TSD;ltr_identity=$sim;Method=structural";
	my $info = "name=$id;classification=LTR/$supfam";
	my $info2 = "ltr_identity=$sim;motif=$motif;tsd=$TSD";
	# print GFF "$chr\t$annotator\trepeat_region\t$element_start\t$element_end\t.\t$strand\t.\tid=repeat_region_$i;$info;$info2\n";
	# print GFF "$chr\t$annotator\ttarget_site_duplication\t$lTSD\t.\t$strand\t.\tid=lTSD_$i;parent=repeat_region_$i;$info;$info2\n" unless $TSD eq "NA";
	print GFF "$chr\t$annotator\tlong_terminal_repeat\t$lLTR_start\t$lLTR_end\t.\t$strand\t.\tid=lLTR_$i;parent=repeat_region_$i;$info;$info2\n";
	print GFF "$chr\t$annotator\t$so\t$lLTR_start\t$rLTR_end\t.\t$strand\t.\tid=LTRRT_$i;parent=repeat_region_$i;$info;$info2\n";
	print GFF "$chr\t$annotator\tlong_terminal_repeat\t$rLTR_start\t$rLTR_end\t.\t$strand\t.\tid=rLTR_$i;parent=repeat_region_$i;$info;$info2\n";
	# print GFF "$chr\t$annotator\ttarget_site_duplication\t$rTSD\t.\t$strand\t.\tid=rTSD_$i;parent=repeat_region_$i;$info;$info2\n" unless $TSD eq "NA";
	print GFF "###\n";
	$i++;
	}