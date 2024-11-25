#!/bin/bash

QUERY=$1
OUTPUT=$1
protein_db=$2
getorf -sequence $QUERY --outseq $OUTPUT.TE.orfs -minsize 400 -reverse
grep '>' $OUTPUT.TE.orfs | awk '{print $1"\t"$2"\t"$4}' | sed 's/\[//g;s/\]//g;s/#/--/g;s/>//g' > $OUTPUT.TE.orfs.R
blastp -query $OUTPUT.TE.orfs -db $protein_db -outfmt 6 | sort -k1,1 -k12,12nr | sort -u -k1,1 | sed 's/#/--/g' > $OUTPUT.TE.blastp.out
#join orfs with their prot hit
  sort -k1,1 $OUTPUT.TE.orfs.R > $OUTPUT.file1
  sort -k1,1 $OUTPUT.TE.blastp.out > $OUTPUT.file2
  join -a1 -11 -21 $OUTPUT.file1 $OUTPUT.file2 | \
   sed -E 's/ /\t/g' | \
   sed 's/DNA\/Maverick/MAV\/Maverick/g;s/DNA\/Crypton/CRY\/Crypton/g;s/LINE\/Penelope/PLE\/Penelope/g;s/LTR\/DIRS/DIRS\/DIRS/g;s/DNA/TIR/g' | \
   awk '{if (NF == 3) {print $0"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"} else {print $0}}' | \
   awk '$NF != "NA" {if ($2 < $3) {print $2"\t"$3"\t"$4"\t"(($2 + (3 * ($9-1)) ))"\t"(($2 + (3 * ($10-1)) ))"\t+"} else {print $2"\t"$3"\t"$4"\t"(($3 + (3 * ($9-1)) ))"\t"(($3 + (3 * ($10-1)) ))"\t-"}}; $NF == "NA" {print $0}' | \
   sed -E 's/--/\t/g' | \
   awk -v LINEcol="3399ff" \
    -v SINEcol="800080" \
    -v TIRcol="ff6666" \
    -v LTRcol="00cc44" \
    -v RCcol="ff6600" \
    -v Low_complexitycol="d1d1e0" \
    -v Satellitecol="ff99ff" \
    -v Simple_repeatcol="#8686ac" \
    -v PLEcol="b2edba" \
    -v DIRScols="fce7bd" \
    -v CRYcols="8f1800" \
    -v MAVcols="669999" \
    -v Unknowncol="c9c9c9" \
  '/LINE\// {print $0"\t"LINEcol} \
  /SINE\// {print $0"\t"SINEcol} \
  /TIR\// {print $0"\t"TIRcol} \
  /LTR\// {print $0"\t"LTRcol} \
  /RC\// {print $0"\t"RCcol} \
  /Low_complexity/ {print $0"\t"Low_complexitycol} \
  /Satellite/ {print $0"\t"Satellitecol} \
  /Simple_repeat/ {print $0"\t"Simple_repeatcol} \
  /Penelope/ {print $0"\t"PLEcol} \
  /DIRS/ {print $0"\t"DIRScols} \
  /CRY/ {print $0"\t"CRYcols} \
  /MAV/ {print $0"\t"MAVcols} \
  !/LINE\// && !/SINE\// && !/TIR\// && !/LTR\// && !/RC\// && !/Low_complexity/ && !/Satellite/ && !/Simple_repeat/ && !/Penelope/ && !/DIRS/ && !/CRY/ && !/MAV/ {print $0"\t"Unknowncol}' | \
  cut -f 1-9 | sort -k2,2n | awk '$NF != "NA" {print $0}; $NF == "NA" {if ($3 < $4) {sub(/NA\tNA\tNA\tNA\tNA\tNA/, "NO\tHIT\t0\t0\t+\tFFFFFF", $0); print $0} else {sub(/NA\tNA\tNA\tNA\tNA\tNA/, "NO\tHIT\t0\t0\t-\twhite", $0); print $0}}' > $OUTPUT.orftetable

awk '{len = $6 - $5 + 1; print $0 "\t" len}' $OUTPUT.orftetable > $OUTPUT.orftetable.with_len
sort -k9,9nr $OUTPUT.orftetable.with_len > $OUTPUT.orftetable.with_len.sorted
awk '
{
    for (i=1; i<=NR-1; i++) {
        if ($5 >= starts[i] && $6 <= ends[i] && ($6 - $5 + 1) <= 0.95 * lens[i]) {
            next
        }
    }
    starts[NR] = $5
    ends[NR] = $6
    lens[NR] = $9
    print $0
}' $OUTPUT.orftetable.with_len.sorted > $OUTPUT.orftetable.clean

rm $OUTPUT.TE.orfs $OUTPUT.TE.orfs.R $OUTPUT.TE.blastp.out $OUTPUT.file1 $OUTPUT.file2 $OUTPUT.orftetable $OUTPUT.orftetable.with_len $OUTPUT.orftetable.with_len.sorted
