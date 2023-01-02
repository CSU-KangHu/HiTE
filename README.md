# HiTE, a Ensemble Method for High-Precision Transposable Element Annotation
[![GitHub](https://img.shields.io/github/license/BioinformaticsCSU/HiTE)](https://gitee.com/kkanghu/HiTE/blob/master/LICENSE)


## Getting started
Clone the [latest release](https://github.com/yangao07/TideHunter/releases):
```
git clone https://github.com/BioinformaticsCSU/HiTE.git
```
## Table of Contents

- [Introduction](#introduction)
  - [Pipeline of HiTE](#pipeline)
  - [Genome coverage by each major subclass](#cover_genome)
- [Installation](#install)
  - [Installing RepeatMasker](#repeatmasker)
  - [Installing genome tools](#genome_tools)
  - [Installing LTR_retriever](#ltrretriever)
  - [Configuring dependencies](#configure)
- [Getting started with toy example in `demo`](#start)
- [Commands and options](#cmd)
- [Input](#input)
- [Output](#output)
  - [Genome annotation information](#repeatmasker_annotation_info)
- [Contact](#contact)

## <a name="introduction"></a>Introduction
HiTE is an efficient TE annotation tool for genome assemblies based on the masking of repeated kmers.

HiTE offers a more **comprehensive** ability to annotate TEs and achieves remarkable efficiency. e.g., more than **21** times faster than RepeatModeler2 in the rice genome. It can serve as a novel solution to the existing methods to promote TE annotation performance.

### <a name="pipeline"></a>Pipeline of HiTE
![输入图片说明](pic/Framework_1.png) 

### <a name="cover_genome"></a>Genome coverage by each major subclass 
![输入图片说明](pic/cover_genome_1.png)

## <a name="install"></a>Installation
### HiTE requires python3, please ensure you run HiTE with python3.

### <a name="repeatmasker"></a>Installing RepeatMasker
It is recommended to download the latest release of RepeatMasker 
from the [RepeatMasker Download page](http://www.repeatmasker.org/RepeatMasker/).

Please install  **RMBlast** as the Sequence Search Engine of RepeatMasker.

### <a name="genome_tools"></a>Installing genome tools
Download [Genome Tools](http://genometools.org/pub/binary_distributions/).

For example:
```
wget http://genometools.org/pub/binary_distributions/gt-1.6.2-Linux_x86_64-64bit-complete.tar.gz
tar zxvf gt-1.6.2-Linux_x86_64-64bit-complete.tar.gz
```

### <a name="ltrretriever"></a>Installing LTR_retriever
Install the latest release of LTR_retriever
from the [LTR_retriever Github page](https://github.com/oushujun/LTR_retriever).
Install [LTRharvest](http://genometools.org/pub/binary_distributions/) and [LTR_FINDER_parallel](https://github.com/oushujun/LTR_FINDER_parallel).

记得chmod +x tools_dir里的所有工具

运行成功，记得check每个阶段的文件是否存在。

### <a name="configure"></a>Configuring dependencies
```
cd /your_path_to/HiTE/ReferenceMode
vim ParamConfig.json
```
Change
- RepeatMasker_Home
- Genome_Tools_Home
- LTR_retriever_Home
- RMBlast_Home 

to the actual installation directories of RepeatMasker, Genome_Tools, LTR_retriever, and RMBlast, respectively.

Then, run

```
cd /your_path_to/HiTE/ReferenceMode
python configure.py
```
to validate all configurations.

## <a name="cmd"></a>Commands and options
```
python main.py $genome_assembly $alias_name

usage: main.py [-h] [-k kmer size] [-t thread num] [-s sensitive mode]
               [--fault_tolerant_bases fault_tolerant_bases] [-o output dir]
               [--min_ltr_complete_len min_ltr_complete_len]
               [--max_ltr_complete_len max_ltr_complete_len]
               [--min_ltr_direct_repeat_len min_ltr_direct_repeat_len]
               [--max_ltr_direct_repeat_len max_ltr_direct_repeat_len]
               [--min_tir_complete_len min_tir_complete_len]
               [--max_tir_complete_len max_tir_complete_len]
               [--min_tir_direct_repeat_len min_tir_direct_repeat_len]
               [--max_tir_direct_repeat_len max_tir_direct_repeat_len]
               [--long_repeat_threshold long_repeat_threshold]
               Genome assembly alias name

run HiTE...

positional arguments:
  Genome assembly       input genome assembly path
  alias name            input alias name

optional arguments:
  -h, --help            show this help message and exit
  -k kmer size          input kmer size, default = [ 31 ]
  -t thread num         input thread num
  -s sensitive mode     sensitive mode, default = [ 0 ]
  --fault_tolerant_bases fault_tolerant_bases
                        the base number of fault tolerant in repeated kmers
                        masking, default = [ 50 ]
  -o output dir         output dir
  --min_ltr_complete_len min_ltr_complete_len
                        Minimum complete LTR length, default = [ 400 ]
  --max_ltr_complete_len max_ltr_complete_len
                        Maximum complete LTR length, default = [ 22000 ]
  --min_ltr_direct_repeat_len min_ltr_direct_repeat_len
                        Minimum LTR direct repeat length, default = [ 100 ]
  --max_ltr_direct_repeat_len max_ltr_direct_repeat_len
                        Maximum LTR direct repeat length, default = [ 6000 ]
  --min_tir_complete_len min_tir_complete_len
                        Minimum complete TIR length, default = [ 1000 ]
  --max_tir_complete_len max_tir_complete_len
                        Maximum complete TIR length, default = [ 40000 ]
  --min_tir_direct_repeat_len min_tir_direct_repeat_len
                        Minimum TIR direct repeat length, default = [ 0 ]
  --max_tir_direct_repeat_len max_tir_direct_repeat_len
                        Maximum TIR direct repeat length, default = [ 1000 ]
  --long_repeat_threshold long_repeat_threshold
                        Threshold of long repeat, default = [ 2000 ]
```

## <a name="start"></a>Getting started with toy example in `demo`
```
cd /your_path_to/HiTE/ReferenceMode
python main.py ../demo/Ecoli_K12_Ref.fasta ecoli
```

## <a name="input"></a>Input
HiTE works with genome assemblies in FASTA, FA, and FNA formats.

## <a name="output"></a>Output
HiTE outputs an annotated consensus TE library in FASTA format.

### <a name="repeatmasker_annotation_info"></a>Genome annotation information
The annotated TE library is further input into RepeatMasker for annotating the whole genome, and the annotation information for the genome is also output.



## <a name="contact"></a>Contact
Kang Hu kanghu@csu.edu.cn

Jianxin Wang jxwang@mail.csu.edu.cn
