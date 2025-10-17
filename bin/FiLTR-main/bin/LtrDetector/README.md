# LtrDetector

https://www.biorxiv.org/content/biorxiv/early/2018/10/22/448969.full.pdf



## Compiling source code

You must have the GNU compiler. Follow the steps to obtain a compatible version.

###  MacOS
Use HomeBrew:
	``` brew install gcc```

Run this command to see list of gcc files installed by HomeBrew:
```brew list gcc ```

Example:
joseph:LtrDetector jdv$ brew list gcc

/usr/local/Cellar/gcc/8.2.0/bin/c++-8

/usr/local/Cellar/gcc/8.2.0/bin/cpp-8

/usr/local/Cellar/gcc/8.2.0/bin/g++-8

On Makefile, change the CXX value to g++-Version as indicated by the listing. In the example above,
CXX = g++-8

### Linux (Ubuntu)
The GNU compiler comes standard. If you do not have it, run:

``` sudo apt-get install gcc ```

## Compile

To locate the source directory:
	``` cd src	```

The following command makes the required directories:
	``` make bin```

The following command makes the binary that is located under the ``bin`` directory:
	``` make tr -j```

To find the binary:
	```cd ../bin```

To find help information:
``` LtrDetector -help ```


### Example run with the defaults:

```LtrDetector -chromDir ~/Data/zmays/Fa -destDir ~/Data/zmays/detector```

### Example run on 4 cores:

```LtrDetector -chromDir ~/Data/thaliana/Fa -destDir ~/Data/thaliana/detector -nThreads 4```

### Output Explanation

LtrDetector output is formatted as a BED file with 18 columns, which correspond to:

1. Sequence identifier.
2. Retrotransposon start coordinate within sequence.
3. Retrotransposon end coordinate within sequence.
4. Left LTR start coordinate.
5. Left LTR end coordinate.
6. Right LTR start coordinate.
7. Right LTR end coordinate.
8. % Identity between left and right LTRs (0-100).
9. Left Target Site Duplication start coordinate.
10. Left Target Site Duplication end coordinate.
11. Right Target Site Duplication start coordinate.
12. Right Target Site Duplication end coordinate.
13. Polypurine Tract start coordinate.
14. Polupurine Tract end coordinate.
15. Strand on chromosome (+ or -).
16. Percentage of purines in Polypurine Tract (0-100).
17. TG motif start coordinate.
18. CA motif end coordinate.

Dashes on any column pertaining to the Target Site Duplication or Polypurine Tract indicate that no such element was found.

This output format provides up to three alternatives for the boundary of the element. Columns 2 and 3 are set to be immediately inside of the TSDs if they exist,
otherwise they are the boundaries denoted by the alignment adjustment step. Columns 4 and 7 report the alignment	
adjusted boundaries directly. Columns 17 and 18 report the boundaries based on the TG..CA box, if one exists.


### To run visualization tool:

Step 1: Run LtrDetector with these three flags -rawScores -cleanedScores -bedFormat

```LtrDetector -chromDir ~/Data/thaliana/Fa -destDir ~/Data/thaliana/detector -rawScores -cleanedScores -bedFormat```

Step 2: Create a virtual environment using venv module of Python3
``` python3 -m venv myEnvironment```
myEnvironment is the name of a folder that will house your virtual environment. You can name it whatever you want.

Step 3: Activate virtual environment
``` source myEnvironment/bin/activate ```
You must now run the following commands from the same terminal. If you want to use multiple terminals, the virtual environment must be activated in all of them.

Step 4: Install python dependencies using provided requirements.txt file.
``` pip install -r requirements.txt ```

Step 5: To run the visualization tool on chr1, pass LTR-RTs found in chr1 (.bed), output directory where the graphs will be stored, and the scores file (raw or cleaned).

 ```python visualize.py ~/Data/thaliana/detector/chr1Detector.bed ~/Data/thaliana/detector/test ~/Data/thaliana/detector/chr1RawScores.csv```

| -arg     | Description | Default |
| ---------------- | ----------- | ------- |
| -chromDir | Directory containing files to be scanned. Files must have .fa extension. | required |
| -destDir | Output directory where the results are stored. | required |
( IMPORTANT: Files under the output directory are deleted at the start of the program.)
| -minLen | Minimum length of complete LTR-RTs. Constrains scoring system and filters. | 400 |
| -maxLen |  Maximum length of complete LTR-RTs. Constrains scoring system and filters. | 22000 |
| -minLenLTR | Minimum length of LTR direct repeat. Constrains filters. | 100 |
| -maxLenLTR | Maximum length of LTR direct repeat. Constrains filters. | 6000 |
( Note run time is highly dependent on this parameter, as it provides an upper bound for alignment length in the boundary adjustment step.)               
| -id | Minimum identity [0-100] between 5' and 3' LTRs. | 85 |
| -k  | Length of k-mers to adjust scoring system. Tradeoff between noise and resistance to mutation. | 13 |
| -plateauSeed | Minimum length of plateaus to be initially considered 'Keep' in merging step. | 10 |
| -nThreads | Number of cores to be used. | 1 |
| -gapTol | Number of base pairs that two plateaus can differ by in height/distance. Affects both plateau merging and pairing steps. | 200 |
|-seqLevel| Forces parallel execution on sequences within multi-FASTA file. Loads all sequences into memory|disabled|
| -rawScores | prints the raw scores to a file called xxxxRawScores.txt under the output directory. | disabled |
| -cleanedScores | prints the scores after merging to a file called xxxxCleanedScores.txt under the output directory. | disabled |
| -nested | searches for nested elements. Results are stored in seperate files (marked as xxxxNestedDetector.bed) under the output directory | disabled |
| -bedFormat | prints BED format without additional annotations (PPT and TSD). | disabled |
| -help | prints this help message. | disabled |


### License

Academic/non-profit use: The software is provided as-is under the GNU GPLv3.

Any restrictions to use for-profit: License needed.
