### Build Ground Truth

1. Establish a parent directory <directory>  Ex: zmays

2. Place all FASTA files (extension .fa) into <directory>/Fa   Ex: zmays/Fa

3. Place all RepeatMasker output files into <directory>/Out  Ex: zmays/Out
   All RepeatMasker files must start with the name of the chromosome (can exclude .fa)  Ex: chromosome1.fa -> chromosome1.out  or chrX.fa -> chrX.fa.out

4. Obtain the Repbase download of all LTRS for the given genome in FASTA format. Ex: Zea_mays_LTR_Retrotransposon.txt

5. ```python generateTruthPipeline.py <directory> <RepbaseFile> ```
Ex: ```python generateTruthPipeline.py zmays Zea_mays_LTR_Retrotransposon.txt```