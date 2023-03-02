# EAHelitron    
<img src="https://github.com/dontkme/PersonalScripts/raw/master/helitron-mini-01.png"  align="right" />

[![releaseVersion](https://img.shields.io/badge/release%20version-1.5.4-green.svg?style=flat)](https://github.com/dontkme/EAHelitron/releases) [![Last-changedate](https://img.shields.io/badge/last%20change-2022--05--27-green.svg)](https://github.com/dontkme/EAHelitron/commit) ![perlVersion](https://img.shields.io/badge/perl-%3E%3D5.10-blue.svg?sytle=flat)


Easy-to-Annotate Helitrons Unix-like Command-Line.              


`EAHelitron` is written in Perl. Used the *Helitron* conservative structure traits: 5’ terminal with TC, 3’ terminal with CTAGt, and before CTAG 2-10 nt has a GC-rich hairpin loop. We used the Perl regular expression(RE) engine and its Embedded-Code Construct to find all the matches and then printed and converted them to a GFF3 format file. Using the above GFF3 file, we can visualize these *Helitrons* in genome visualization tools such as IGV, Gbrowse, and Jbrowse, and easily characterize the captured genes.

EAHelitron is a Unix-like program that you can run it on all Perl 5.10+ supported machines, and write commands in your shell scripts or through pipes. Linux, Mac OS, and Windows tests passed.





## Getting Started



### Prerequisites

Make sure you have Perl on your system.
Type these words into your system's terminal.
```
perl -v
```
If the terminal shows the Perl version information, then we can download the EAHelitron files. https://github.com/dontkme/EAHelitron/archive/master.zip



### Installing


Just unzip the zip file.


```
unzip EAHelitron-master.zip
```

And go to the extracted folder and run EAHelitron.

```
cd EAHelitron-master
perl EAHelitron -h
```

If the screen displays help and version information. It works.

### Running 

**Basic mode：**

Example: Predict Helitrons and search for their 5'TC terminals in the upstream 20,000bp range.

```
perl EAHelitron –o testEAHout –u 20000 teat.fas
```
Or 

```  
./EAHelitron –o testEAHout –u 20000 teat.fas
```   
#### Options:
        
         [-o string|outprefix Default: EAHeli_out]
         [-u int|upstream length Default: 3000]
         [-d int|downstream length Default: 500]
         Advanced options:
         [-T string|TC pattern. User's 5'TC pattern]
         [-H string|Hairpin pattern. User's Hairpin left pattern]
         [-r int[0-5]|CTRRt 3' terminal fuzzy level;
                 0: CTAGT
                 1: CT[AG]GT
                 2: CTA[AG]T
                 3: CT[AG]{2}T
                 4: CT[AG]{2}.{1}
                 5: CTAG.{1}
                 Default: 0]

We also provide **EAHelitron_P**, a multi-threaded version that can speed up running in a large genome.

(**Need Perl Parallel::ForkManager.** You could install it by command: cpan Parallel::ForkManager )

```
perl EAHelitron_P –p 8 –o testEAHout –u 20000 teat.fas
```
**-p**: How many threads to use. It is recommended not to exceed the number of sequences contained in the fasta file you input.

---

#### **Advanced options:** 
**-r**: CTRRt 3' terminal fuzzy level:
6 fuzzy levels of CTRRt terminal [0-5]

```
perl EAHelitron_P –p 8 -r 3 –o testEAHout_r3 teat.fas
```
Users can enter their own patterns (Perl RE) to predict Helitrons.
(**Warning**: Advanced options may significantly increase the false positive rate, only for exploring).


**-H**: Use Hairpin left sequence RE pattern:

Example:
1. Only use a GC as hairpin left sequence pattern:
```
perl EAHelitron_P –p 8 -H "GC" –o testEAHout_H_GC teat.fas
```
2. or Use 1 **GC**(G or C) 1 **AT**(A or T) and 5 **GC** or(**|**) 6 **GC** as haripin left sequence pattern:
```
perl EAHeliton_P –p 8 -H "[GC]{1}[AT]{1}[GC]{5}|[GC]{6}" –o testEAHout_H_GC teat.fas
```

**-T**: Use 5' TC sequence RE pattern:

Example:

1. Only use 'TC' as 5' TC sequence RE pattern:
```
perl EAHelitron_P –p 8 -T "TC" –o testEAHout_T_TC teat.fas
```
2. Use TC or(**|**) TCT.TACTA.T as 5' TC sequence RE pattern (The dot '.' to indicate any character, we can use [ATCGN] instead of '.', if you like):
```
perl EAHeliton_P –p 8 -T "TC|TCT.TACTA.T" –o testEAHout_T teat.fas
```

We can also use them in combination:
```
perl EAHelitron_P –p 8 -T "TC" -H "GC" –o testEAHout_T_H teat.fas
```
Or more complex:
```
perl EAHeliton_P –p 8 -H "[GC]{1}[AT]{1}[GC]{5}|[GC]{6}" -T "TC|TCT.TACTA.T" –o testEAHout_T_H teat.fas
```

---

### Outputs:
The outputs are named like **EAHout.3.txt EAHout.5.txt EAHout.5.fa EAHout.gff3 EAout.u20000.fas.** (The prefix 'EAHout' can be set with the -o option, 20000 is the value of your -u option)

***.3.txt**: All 3’ terminal sequences with a 10nt left flank and a 4nt right flank in fasta format. All sequences are named by their local chromosome name, an 'H' means Helitron, '.3' suffix to mark they are 3’ terminals. The minus strain terminals have a 'tr' prefix. (e.g. Chr1H10.3, trChr5H40.3).

***.5.txt**: All 5’terminal sequences were matched in the 3’ terminal’s upstream sequences, with a 5 nt left flank and a 20 nt right flank. The names of sequences have .5.1 suffix to mark they are 5’ terminal and the match orders numbers. (e.g. Chr1H10.5.1,trChr5H40.5.2) 

***.5.fa**: Possible full-length Helitron sequences which start with 5’ terminal and end with 3’ terminal. (Only for Helitrons)

***.u*.fas**: All 3’ terminal upstream sequences.  

***.d*.fas**: All 3’ terminal downstream sequences.

***.gff3**: Helitron location information in GFF3 format.

***.bed**: Helitron 3'-ends location in bed format.

***.len.txt**: Summary of genome sequences length, Helitron counts, and Helitron Densities

## History

(EAHelitron) v1.5400 2022/05/27 New default hairpin-left-sequence pattern (allow 2 [AT] in 5 [GC]). The default results for *A. thaliana* increased from 665 to 708, and the false positive rate increased slightly from 5.91% to 6.47%.

(EAHelitron) v1.5300 2021/06/25 Add a feature. User TC pattern and hairpin left sequence pattern options.

(EAHelitron) v1.5200 2020/09/25 Use a new regular expression, which is based on BioPerl, to get chromosome names, in order to adapt to more cases. Thanks to [Darcy Jones](https://github.com/darcyabjones)'s advice.

(EAHelitron) v1.5100 2019/06/10 Add a feature. Output a BED of 3’-ends.  (*.bed)

(EAHelitron) v1.5000 2019/03/08 Add a feature. Output a summary of genome sequences length, Helitron counts and Helitron Densities. (*.len.txt)

(EAHelitron) v1.4000 2018/08/31 Add CTAGT fuzzy level [0-5] option.

(EAHelitron) v1.3100 2017/09/19 Add snp switch.

(EAHelitron)v1.3000 2017/08/29 Add a feature for downstream sequences.

...

(EAHelitron) v1.0000 2016/09/22 first version upload to GitHub.

## Citation

If you've found EAHelitron useful in your work, please cite our [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2945-8):
> Hu, K. et al. Helitron distribution in Brassicaceae and whole genome Helitron density as a character for distinguishing plant species. BMC bioinformatics 20, 1-20, doi:ARTN 354
10.1186/s12859-019-2945-8 (2019).

## Authors

**Hu Kaining** - *Initial work* - [dontkme](https://github.com/dontkme)

## License
![licence](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details


