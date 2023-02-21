# HiTE, an Ensemble Method for High-Precision Transposable Element Annotation
[![GitHub](https://img.shields.io/badge/python-3-blue)](https://www.python.org/)
[![GitHub](https://img.shields.io/badge/license-GPL--3.0-green)](https://github.com/CSU-KangHu/HiTE/blob/master/LICENSE)

<!-- ## <a name="introduction"></a>Introduction
We have developed an ensemble method for high-precision TE annotation, known as **HiTE**, 
which has undergone extensive benchmarking and has proven to be the best TE annotation tool available. 
HiTE achieved the highest precision and discovered the most gold standard TE models based on four model species: 
Oryza sativa, Caenorhabditis briggsae, Drosophila melanogaster, and Danio rerio. Furthermore, HiTE can discover 
novel TEs with low copy numbers that are not included in known libraries. -->

## Table of Contents

<!-- - [Introduction](#introduction)
  - [The workflow of HiTE](#pipeline)
  - [Performance](#performance)
-->
- [Installation](#install)
  - [One-step installation](#one-step)
  - [Step-by-step installation](#step-step)
- [Quick start](#start)
- [Usage](#cmd)
- [Input](#inputs)
- [Output](#outputs)
- [Replace the Dfam library in RepeatMasker](#classified)
- [Contact](#contact)

<!-- ### <a name="pipeline"></a>The workflow of HiTE
![输入图片说明](pic/Framework.png) 

### <a name="performance"></a>Performance comparison of general-purpose TE annotators based on benckmarking method of RepeatModeler2 and EDTA 
![输入图片说明](pic/RM2_results.png)

![输入图片说明](pic/EDTA_results.png) -->

## <a name="install"></a>Installation
### <a name="one-step"></a>Option 1. One-step installation (recommended)
#### Download project 
```
git clone https://github.com/CSU-KangHu/HiTE.git
```
Find the **yml** file in the project directory and run:
```
conda env create --name HiTE -f environment.yml

conda activate HiTE
```

#### Replacing the Dfam library in RepeatMasker (optional)
HiTE is ready to go, but if you require the TE library to be **comprehensively classified**, you need to configure RepeatMasker with the complete Dfam library.
Use `--classified 0` if you do not need classified TE models.

[The simplest way to replace the Dfam library](#classified)

### <a name="step-step"></a>Option 2. Step-by-step installation
<details>
<summary>Show more details</summary>

#### 1. Download project 
```
git clone https://github.com/CSU-KangHu/HiTE.git
```

#### 2. installing python3
HiTE requires python3, please ensure you run HiTE with python3.

#### 3. installing RMBlast
Please install  [RMBlast](https://www.repeatmasker.org/rmblast/).

#### <a name="genome_tools"></a>4. Installing genome tools
Download [Genome Tools](http://genometools.org/pub/binary_distributions/).

<!-- For example:
```
wget http://genometools.org/pub/binary_distributions/gt-1.6.2-Linux_x86_64-64bit-complete.tar.gz
tar zxvf gt-1.6.2-Linux_x86_64-64bit-complete.tar.gz
```
-->

#### <a name="ltrretriever"></a>5. Installing LTR_retriever
Install the latest release of LTR_retriever
from the [LTR_retriever Github page](https://github.com/oushujun/LTR_retriever).

#### <a name="repeatmodeler"></a>6. Installing RepeatMasker and RepeatModeler2 (Optional)
It is recommended to download the [RepeatMasker](http://www.repeatmasker.org/RepeatMasker/) with the complete Dfam 3.6 library.

HiTE uses RepeatClassifier from [RepeatModeler2](https://www.repeatmasker.org/RepeatModeler/) to classify the TE models. 
Please follow the documentation to install RepeatModeler2 and configure RepeatMasker.

If you do not need classified TE models, you can skip this step and run HiTE with `--classified 0`.

#### <a name="configure"></a>7. Configuring dependencies
```
cd /your_path_to/HiTE/ReferenceMode
vim ParamConfig.json
```
Change
- Genome_Tools_Home
- LTR_retriever_Home
- RMBlast_Home
- RepeatModeler_Home

to the actual installation directories of RepeatMasker, Genome_Tools, LTR_retriever, and RMBlast, respectively.

Then, run

```
cd /your_path_to/HiTE/ReferenceMode
python configure.py
```
to validate all configurations.
</details>


## <a name="start"></a>Quick start
```
cd /your_path_to/HiTE/ReferenceMode
python main.py -g ../demo/genome.fa -t 48 -o ../demo/test --plant 0
```

If the following files exist in the **demo/test** directory, it means the program runs successfully: 
* confident_helitron_0.fa (9.1 KB)
* confident_other_0.fa (0 KB)
* confident_TE.cons.fa (117 KB)
* confident_tir_0.fa (159 KB)
* genome_all.fa.harvest.scn (11 KB)
* genome.rename.fa.finder.combine.scn (1.3 KB)
* Genome.rename.fa.LTRlib.fa (48 KB)
* longest_repeats_0.fa (1.2 MB)
* longest_repeats_0.flanked.fa (1.5 MB)

Note:
1. Please make sure you execute the **main.py** script under the **/your_path_to/HiTE/ReferenceMode** directory.
2. To avoid automatic deletion of files, set the output path parameter ```-o``` to a new directory.
 
#### Possible issues
Type `RepeatMasker` in the conda HiTE virtual environment, and an error similar to the following occurred:

`undefined symbol: Perl_xs_apiversion_bootcheck`

This error is caused by the incompatibility between the local perl of the system and the perl installed by conda. The solution is:

`export PERL5LIB=/`

## <a name="inputs"></a>Inputs
HiTE works with genome assemblies in **fasta**, **fa**, and **fna** formats using `-g`.


For other optional parameters, please refer to [Commands and options](#cmd).

## <a name="outputs"></a>Outputs
HiTE outputs many temporary files, which allow you to quickly restore the previous 
running state in case of any interruption during the running process. If
the pipeline completes successfully, the following files are generated:

* longest_repeats_*.fa.
* confident_tir_*.fa
* confident_helitron_*.fa
* confident_other_*.fa
* genome_all.fa.harvest.scn
* ${ref_name}.finder.combine.scn
* ${ref_name}.LTRlib.fa
* confident_TE.cons.fa
* confident_TE.cons.fa.classified

Note that "*" represents the number of blocks that the genome is divided into.
For example, if the genome input is 400 MB and the chunk size input is set to 100,
then * is equal to 4 (400/100), and you can find 4 files: repeats_0.fa, repeats_1.fa,
repeats_2.fa, and repeats_3.fa in your output directory.

**confident_TE.cons.fa** and **confident_TE.cons.fa.classified** are the 
unclassified and classified TE libraries generated by HiTE, respectively.

It is worth noting that **confident_TE.cons.fa.classified** is generated by RepeatClassifier from 
RepeatModeler2, which depends on the Dfam library in RepeatMasker.

**confident_TE.cons.fa.classified** can be used directly as TE library in RepeatMasker by `-lib`.


## <a name="classified"></a>Replace the Dfam library in RepeatMasker
Since the Dfam library included in RepeatMasker by default is not complete, it will seriously affect the classification effect.
We recommend updating RepeatMasker with the complete Dfam 3.6 library as described at http://www.repeatmasker.org/RepeatMasker/, including download, unpack, and reconfiguration.
We also provide an optional method to avoid the big Dfam.h5.gz (15 GB) download and reconfiguration, as follows:
1. download **RepeatMasker_Lib.zip** from google drive: 
https://drive.google.com/file/d/1vQLamfINdJ5iDwggYigWKe7Gor4t6JMK/view?usp=sharing

2. upload **RepeatMasker_Lib.zip** to RepeatMasker/Libraries, where RepeatMasker is your installation directory of RepeatMasker.
   (e.g., ~/anaconda2/envs/HiTE/share/RepeatMasker)

3. `cd RepeatMasker/Libraries`

4. `unzip RepeatMasker_Lib.zip && mv RepeatMasker_Lib/* ./`

<!--
```
1. cd RepeatMasker/Libraries/
2. wget https://www.dfam.org/releases/Dfam_3.6/families/Dfam.h5.gz
3. gunzip Dfam.h5.gz
```

Run Configure Script
```
1. cd RepeatMasker 
2. perl ./configure
``` -->


## <a name="cmd"></a>Usage
Type `python main.py -h` for help.
```
The simplest command:
python main.py -g $genome_assembly -o $output_dir

Most frequently used commands:
python main.py -g $genome_assembly -o $output_dir -t 40 --chunk_size 400 --plant 0 --recover 1

uusage: main.py [-h] [-g Genome assembly] [-t thread num]
               [--chunk_size chunk_size] [--plant is_plant]
               [--remove_nested remove_nested] [--classified classified]
               [--recover recover] [--debug is_debug] [-o output dir]
               [--flanking_len flanking_len]
               [--fixed_extend_base_threshold fixed_extend_base_threshold]
               [--tandem_region_cutoff tandem_region_cutoff]
               [--max_repeat_len max_repeat_len]
               [--chrom_seg_length chrom_seg_length]
               [--global_flanking_filter global_flanking_filter]

########################## HiTE, version 2.0.3 ##########################

optional arguments:
  -h, --help            show this help message and exit
  -g Genome assembly    input genome assembly path
  -t thread num         input thread num, default = [ 40 ]
  --chunk_size chunk_size
                        the chunk size of large genome, default = [ 400 MB ]
  --plant is_plant      is it a plant genome, 1: true, 0: false. default = [ 1
                        ]
  --remove_nested remove_nested
                        Whether to clear the nested TE, 1: true, 0: false.
                        default = [ 1 ]
  --classified classified
                        Whether to classify TE models, HiTE uses
                        RepeatClassifier from RepeatModeler to classify TEs,
                        1: true, 0: false. default = [ 1 ]
  --recover recover     Whether to enable recovery mode to avoid repeated
                        calculations, 1: true, 0: false. default = [ 0 ]
  --debug is_debug      Open debug mode, 1: true, 0: false. default = [ 0 ]
  -o output dir         output dir
  --flanking_len flanking_len
                        the flanking length of repeat to find the true
                        boundary, default = [ 50 ]
  --fixed_extend_base_threshold fixed_extend_base_threshold
                        the base number of extend base, default = [ 1000 ]
  --tandem_region_cutoff tandem_region_cutoff
                        Cutoff of the raw masked repeat regarded as tandem
                        region, default = [ 0.5 ]
  --max_repeat_len max_repeat_len
                        the maximum length of repeat, default = [ 30000 ]
  --chrom_seg_length chrom_seg_length
                        the length of genome segments, default = [ 500000 ]
  --global_flanking_filter global_flanking_filter
                        Whether to filter false positives by global flanking
                        alignment, significantly reduce false positives but
                        require more memory, especially when inputting a large
                        genome. 1: true (require more memory), 0: false.
                        default = [ 1 ]
```

## <a name="contact"></a>Contact
Kang Hu kanghu@csu.edu.cn

Jianxin Wang jxwang@mail.csu.edu.cn
