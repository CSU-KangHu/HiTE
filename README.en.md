# KmerRepFinder: a comprehensive kmer-based  _de novo_  annotation method of transposable elements
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ltr_retriever/badges/license.svg)](https://gitee.com/kkanghu/KmerRepFinder/blob/master/LICENSE)


## Getting started
Download the [latest release](https://github.com/yangao07/TideHunter/releases):
```
git clone https://gitee.com/kkanghu/KmerRepFinder.git
```
## Table of Contents

- [Introduction](#introduction)
- [Installation](#install)
  - [Installing RepeatMasker](#repeatmasker)
  - [Installing LTR_retriever](#ltrretriever)
  - [Configuring depedencies](#configure)
- [Getting started with toy example in `demo`](#start)
- [Commands and options](#cmd)
- [Input](#input)
- [Output](#output)
  - [Annotation of genome](#repeatmasker_annotation_info)
- [Contact](#contact)

## <a name="introduction"></a>Introduction
KmerRepFinderis an efficient TE annotation tool for genome assemblies based on the masking of repeated kmers.

KmerRepFinder offers a more comprehensive ability to annotate TEs and achieves remarkable efficiency, which is more than 21 times faster than RepeatModeler2 in the rice genome and is expected to serve as a novel solution to the existing methods to promote TE annotation performance.

## <a name="install"></a>Installation

### <a name="repeatmasker"></a>Installing RepeatMasker
It is recommended to download the latest release of RepeatMasker 
from the [RepeatMasker Download page](http://www.repeatmasker.org/RepeatMasker/).

Please install  **RMBlast** as Sequence Search Engine of RepeatMasker.

### <a name="ltrretriever"></a>Installing LTR_retriever
Install the latest release of LTR_retriever
from the [LTR_retriever Github page](https://github.com/oushujun/LTR_retriever).

### <a name="configure"></a>Configuring dependencies
modify ParamConfig.json,change the install path.
```
python configure.py
```

## <a name="start"></a>Getting started with toy example in `test_data`
```
python main.py -R ../demo/Ecoli_K12_Ref.fasta -a ecoli
```

## <a name="cmd"></a>Commands and options
```
usage: main.py [-h] [-R Reference Path] [-k kmer size] [-t thread num]
               [-a alias name] [-s sensitive mode]
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

run kmerRepFinder...

optional arguments:
  -h, --help            show this help message and exit
  -R Reference Path     input reference path
  -k kmer size          input kmer size, default = [ 31 ]
  -t thread num         input thread num
  -a alias name         input alias name
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


## <a name="input"></a>Input
KmerRepFinder works with genome assemblies in FASTA, FA, and FNA formats.

## <a name="output"></a>Output
KmerRepFinder can output an annotated consensus TE library in FASTA format.

### <a name="repeatmasker_annotation_info"></a>Annotation information using RepeatMasker
Annotation information using RepeatMasker



## <a name="contact"></a>Contact
Kang Hu kanghu@csu.edu.cn

Jianxin Wang jxwang@mail.csu.edu.cn