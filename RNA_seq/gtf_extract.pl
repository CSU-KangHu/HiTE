#!/usr/bin/env perl -w 
use strict;

die "perl $0 <gtf>  <id.txt(transcript_id\\tgene_id)>  > out.gtf\n" unless @ARGV == 2 ;

my $gtf = shift;
my $id = shift;

## read id table
my %gid;
my %tid;
open ID, $id  or die $!;
while(<ID>){
  chomp;
  my @t = split ;
  next if (scalar(@t) <2);
  $tid{$t[0]} = 1;
  $gid{$t[1]} = 1;
}
close ID;

## read gtf
open GTF, $gtf or die $!;
while( <GTF>){
  chomp ;
  my $flag = 0;
  my @t = split /\t/; 
  next if ( scalar(@t) != 9);
  ## keep gene line
  if($t[8] =~ /gene_id "(\S+?)"/){
    my $gid = $1;
    if($t[2] eq "gene" &&  exists $gid{$gid}){
      $flag = 1;
    }
  }
  ## keep lines with transcript_id
  if($t[8] =~ /transcript_id "(\S+?)"/){
    my $tid = $1;
    if( exists $tid{$tid}){
      $flag = 1;
    }
  }
  if($flag == 1){
    print "$_\n";
  }
}
close GTF;
