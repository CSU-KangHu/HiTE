# HiTE, a progressive Method for High-Precision Transposable Element Annotation
[![GitHub](https://img.shields.io/badge/python-3-blue)](https://www.python.org/)
[![GitHub](https://img.shields.io/badge/license-GPL--3.0-green)](https://github.com/CSU-KangHu/HiTE/blob/master/LICENSE)

## Table of Contents
- [Installation](#install)
  - [Run with Singularity](#install_singularity)
  - [Run with Docker](#install_docker)
  - [Run with Conda](#install_conda)
  <!-- - [Step-by-step installation](#step-step) -->
- [Quick start](#start)
- [Usage](#cmd)
- [Input](#inputs)
- [Output](#outputs)
- [Replace the Dfam library in RepeatMasker](#classified)
- [Question and Answer](#QA)

## <a name="install"></a>Installation
### <a name="install_singularity"></a>Option 1. Run with Singularity (recommended)
```sh
# pull singularity image (once for all). There will be a HiTE.sif file.
singularity pull HiTE.sif docker://kanghu/hite:2.0.3

# run HiTE
singularity run -B /home/user:/home/user --pwd /HiTE/ReferenceMode HiTE.sif python main.py \
 --genome genome.fa \
 --thread 40 \
 --outdir output_dir \
 --plant 1 [other parameters]
```

### <a name="install_docker"></a>Option 2. Run with Docker
```sh
# pull docker image (once for all).
docker pull docker://kanghu/hite:2.0.3

# run HiTE
docker run -v /home/user:/home/user kanghu/hite:2.0.3 python main.py \
 --genome genome.fa \
 --thread 40 \
 --outdir output_dir \
 --plant 1 [other parameters]
```

### <a name="install_conda"></a>Option 3. Run with conda
```sh
# Download project 
git clone https://github.com/CSU-KangHu/HiTE.git

# Find the **yml** file in the project directory and run
cd HiTE
conda env create --name HiTE -f environment.yml
conda activate HiTE

# run HiTE
cd HiTE/ReferenceMode
python main.py \
 --genome genome.fa \
 --thread 40 \
 --outdir output_dir \
 --plant 1 [other parameters]
```

#### Updating the Dfam library in RepeatMasker (optional)
* HiTE is ready to go!
* use `--classified 0` if you do not need classified TE models.
* If you require the TE library to be **comprehensively classified**, you need to configure RepeatMasker with the complete Dfam library.
[The simplest way to update the Dfam library](#classified)
* If you installed HiTE with Singularity or Docker, you can skip this step.


<!--
### <a name="step-step"></a>Option 4. Step-by-step installation
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
-->

## <a name="start"></a>Quick start
```
cd /your_path_to/HiTE/ReferenceMode
python main.py --genome ../demo/genome.fa --thread 48 --outdir ../demo/test --plant 0
```

If the following files exist in the **demo/test** directory, it means the program runs successfully: 
```text
demo/test/
├── confident_helitron_0.fa (9.1 KB)
├── confident_other_0.fa (0 KB)
├── confident_TE.cons.fa (117 KB)
├── confident_tir_0.fa (159 KB)
├── genome_all.fa.harvest.scn (11 KB)
├── genome.rename.fa.finder.combine.scn (1.3 KB)
├── Genome.rename.fa.LTRlib.fa (48 KB)
├── longest_repeats_0.fa (1.2 MB)
└── longest_repeats_0.flanked.fa (1.5 MB)
```

Note:
1. Please make sure you execute the **main.py** script under the **/your_path_to/HiTE/ReferenceMode** directory.
2. To avoid automatic deletion of files, set the output path parameter ```--outdir``` to a new directory.

## <a name="inputs"></a>Inputs
HiTE works with genome assemblies in **fasta**, **fa**, and **fna** formats using `--genome`.


For other optional parameters, please refer to [Commands and options](#cmd).

## <a name="outputs"></a>Outputs
HiTE outputs many temporary files, which allow you to quickly restore the previous 
running state in case of any interruption during the running process. If
the pipeline completes successfully, the output directory should look like the following:
```text
output_dir/
├── longest_repeats_*.fa
├── confident_tir_*.fa
├── confident_helitron_*.fa
├── confident_other_*.fa
├── genome_all.fa.harvest.scn
├── ${ref_name}.finder.combine.scn
├── ${ref_name}.LTRlib.fa
├── confident_TE.cons.fa
└── confident_TE.cons.fa.classified
```

1. **confident_TE.cons.fa** and **confident_TE.cons.fa.classified** are the 
unclassified and classified TE libraries generated by HiTE, respectively.
2. **confident_TE.cons.fa.classified** can be used directly as TE library in RepeatMasker by `-lib`.
3. It is worth noting that **confident_TE.cons.fa.classified** is generated by RepeatClassifier from 
RepeatModeler2, which depends on the Dfam library in RepeatMasker.
4. Note that "*" represents the number of blocks that the genome is divided into.
For example, if the genome input is 400 MB and the chunk size input is set to 100,
then * is equal to 4 (400/100), and you can find 4 files: repeats_0.fa, repeats_1.fa,
repeats_2.fa, and repeats_3.fa in your output directory.

## <a name="classified"></a>Replace the Dfam library in RepeatMasker
Since the Dfam library included in RepeatMasker by default is not complete, it will seriously affect the classification effect.
We recommend updating RepeatMasker with the complete Dfam 3.6 library as described at http://www.repeatmasker.org/RepeatMasker/, including download, unpack, and reconfiguration.
We also provide an optional method to avoid the big Dfam.h5.gz (15 GB) download and reconfiguration, as follows:
1. download **RepeatMasker_Lib.zip** from [google drive](https://drive.google.com/file/d/1vQLamfINdJ5iDwggYigWKe7Gor4t6JMK/view?usp=sharing) or [github](https://github.com/CSU-KangHu/TE_annotation):

2. upload **RepeatMasker_Lib.zip** to RepeatMasker/Libraries, where RepeatMasker is your installation directory of RepeatMasker.
   (e.g., ~/anaconda2/envs/HiTE/share/RepeatMasker)

3. `cd RepeatMasker/Libraries`

4. `unzip RepeatMasker_Lib.zip && mv RepeatMasker_Lib/* ./`


## <a name="cmd"></a>Usage
Type `python main.py -h` for help.
```
The simplest command:
python main.py --genome $genome_assembly --outdir $output_dir

Most frequently used commands:
python main.py --genome $genome_assembly --outdir $output_dir --thread 40 --chunk_size 400 --plant 0 --recover 1

usage: main.py [-h] [--genome genome] [--thread thread_num]
               [--chunk_size chunk_size] [--plant is_plant]
               [--remove_nested is_remove_nested] [--classified is_classified]
               [--recover is_recover] [--debug is_debug] [--outdir output_dir]
               [--flanking_len flanking_len]
               [--fixed_extend_base_threshold fixed_extend_base_threshold]
               [--tandem_region_cutoff tandem_region_cutoff]
               [--max_repeat_len max_repeat_len]
               [--chrom_seg_length chrom_seg_length]
               [--global_flanking_filter global_flanking_filter]

########################## HiTE, version 2.0.3 ##########################

optional arguments:
  -h, --help            show this help message and exit
  --genome genome       Input genome assembly path
  --thread thread_num   Input thread num, default = [ 40 ]
  --chunk_size chunk_size
                        The chunk size of large genome, default = [ 400 MB ]
  --plant is_plant      Is it a plant genome, 1: true, 0: false. default = [ 1
                        ]
  --remove_nested is_remove_nested
                        Whether to unwrap the nested TE, 1: true, 0: false.
                        default = [ 1 ]
  --classified is_classified
                        Whether to classify TE models, HiTE uses
                        RepeatClassifier from RepeatModeler to classify TEs,
                        1: true, 0: false. default = [ 1 ]
  --recover is_recover  Whether to enable recovery mode to avoid starting from
                        the beginning, 1: true, 0: false. default = [ 0 ]
  --debug is_debug      Open debug mode, and temporary files will be kept, 1:
                        true, 0: false. default = [ 0 ]
  --outdir output_dir   The path of output directory; It is recommended to use
                        a new directory to avoid automatic deletion of
                        important files.
  --flanking_len flanking_len
                        The flanking length of candidates to find the true
                        boundaries, default = [ 50 ]
  --fixed_extend_base_threshold fixed_extend_base_threshold
                        The length of variation can be tolerated during
                        pairwise alignment, default = [ 1000 ]
  --tandem_region_cutoff tandem_region_cutoff
                        Cutoff of the candidates regarded as tandem region,
                        default = [ 0.5 ]
  --max_repeat_len max_repeat_len
                        The maximum length of a single repeat, default = [
                        30000 ]
  --chrom_seg_length chrom_seg_length
                        The length of genome segments, default = [ 500000 ]
  --global_flanking_filter global_flanking_filter
                        Whether to filter false positives by global flanking
                        alignment, significantly reduce false positives but
                        require more memory, especially when inputting a large
                        genome. 1: true (require more memory), 0: false.
                        default = [ 1 ]
```

## <a name="QA"></a>Question and Answer
You may want to check out this [Q&A](https://github.com/CSU-KangHu/HiTE/wiki) page for best practices and get answered.
* [Issues with installation and usage](https://github.com/CSU-KangHu/HiTE/wiki/Issues-with-installation-and-usage)
* [Nextflow of HiTE](https://github.com/CSU-KangHu/HiTE/wiki/Nextflow-of-HiTE)
* [Case study](https://github.com/CSU-KangHu/HiTE/wiki/Case-study)