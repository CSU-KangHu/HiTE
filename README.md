# HiTE: a fast and accurate dynamic boundary adjustment approach for full-length Transposable Elements detection and annotation in Genome Assemblies

# HiTE

[![GitHub](https://img.shields.io/badge/python-3-blue)](https://www.python.org/)
[![GitHub](https://img.shields.io/badge/license-GPL--3.0-green)](https://github.com/CSU-KangHu/HiTE/blob/master/LICENSE)
[![DockerHub](https://img.shields.io/badge/Singularity-support-blue)](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)
[![DockerHub](https://img.shields.io/badge/Docker-support-orange)](https://hub.docker.com/repository/docker/kanghu/hite/general)
[![Conda](https://img.shields.io/badge/Conda-support-yellow)](https://docs.conda.io/en/latest/)
[![Nextflow](https://img.shields.io/badge/Nextflow-support-red)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10576547.svg)](https://doi.org/10.5281/zenodo.10576547)

`HiTE` is a Python software that uses a dynamic boundary adjustment approach to detect and annotate full-length Transposable Elements in Genome Assemblies. In comparison to other tools, HiTE demonstrates superior performance in detecting a greater number of full-length TEs.

Similar works include:
* [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler)
* [EDTA](https://github.com/oushujun/EDTA)
* [EarlGrey](https://github.com/TobyBaril/EarlGrey)


<!-- 
## Application

HiTE has been successfully applied to multiple practical applications, and you can refer to our most recent case for reference.

[TE annotation in einkorn assemblies using HiTE](https://github.com/CSU-KangHu/HiTE/wiki/TE-annotation-in-einkorn-assemblies-using-HiTE)
-->

## Table of Contents
- [Installation](#install)
  - [Dowload project](#download)
  - [Run with Conda](#install_conda)
  - [Run with Singularity](#install_singularity)
  - [Run with Docker](#install_docker)
  - [Run with nextflow](#install_nextflow)
  <!-- - [Step-by-step installation](#step-step) -->
- [Demo data](#demo)
- [Code Structure](#code)
- [Usage](#cmd)
- [Input](#inputs)
- [Output](#outputs)
- [Experiment reproduction](#ER)
- [More tutorials](#QA)

## <a name="install"></a>Installation
Recommended Hardware requirements: 40 CPU processors, 128 GB RAM. 

Recommended OS: (Ubuntu 16.04, CentOS 7, etc.)

### <a name="download"></a>Dowload project
```sh
git clone https://github.com/CSU-KangHu/HiTE.git
# Alternatively, you can download the zip file directly from the repository.

cd HiTE && chmod +x tools/* bin/NeuralTE/tools/* bin/LTR_FINDER_parallel-master/bin/LTR_FINDER.x86_64-1.0.7/ltr_finder
```

For common issues related to installation and usage, please visit: [https://github.com/CSU-KangHu/HiTE/wiki/Issues-with-installation-and-usage](https://github.com/CSU-KangHu/HiTE/wiki/Issues-with-installation-and-usage)

### <a name="install_conda"></a>Option 1. Run with conda
```sh
# Find the **yml** file in the project directory and run
cd HiTE
conda env create --name HiTE -f environment.yml
conda activate HiTE

# run HiTE
python main.py \
 --genome ${genome} \
 --thread ${thread} \
 --outdir ${output_dir} \
 [other parameters]
 
 # e.g., my command: python main.py 
 # --genome /home/hukang/HiTE/demo/genome.fa 
 # --thread 40 
 # --outdir /home/hukang/HiTE/demo/test/
```

### <a name="install_singularity"></a>Option 2. Run with Singularity
```sh
# pull singularity image (once for all). There will be a HiTE.sif file.
singularity pull HiTE.sif docker://kanghu/hite:3.2.0

# run HiTE
singularity run -B ${host_path}:${container_path} --pwd /HiTE ${pathTo/HiTE.sif} python main.py \
 --genome ${genome} \
 --thread ${thread} \
 --outdir ${output_dir} \
 [other parameters]
 
 # (1) The options "--genome" and "--outdir" need to be specified as absolute paths.
 # (2) The option "-B" is used to specify directories to be mounted.
 #     It is recommended to set ${host_path} and ${container_path} to your user directory, and ensure 
 #     that all input and output files are located within the user directory.
 # (3) "--pwd /HiTE" and "python main.py" do not need to be changed.
 
 # e.g., my command: singularity run -B /home/hukang:/home/hukang --pwd /HiTE /home/hukang/HiTE.sif python main.py 
 # --genome /home/hukang/HiTE/demo/genome.fa 
 # --thread 40 
 # --outdir /home/hukang/HiTE/demo/test/
```

### <a name="install_docker"></a>Option 3. Run with Docker
```sh
# pull docker image (once for all).
docker pull kanghu/hite:3.2.0

# run HiTE
docker run -v ${host_path}:${container_path} kanghu/hite:3.2.0 python main.py \
 --genome ${genome} \
 --thread ${thread} \
 --outdir ${output_dir} \
 [other parameters]
 
 # (1) The options "--genome" and "--outdir" need to be specified as absolute paths.
 # (2) The option "-v" is used to specify directories to be mounted.
 #     It is recommended to set ${host_path} and ${container_path} to your user directory, and ensure 
 #     that all input and output files are located within the user directory.
 
 # e.g., my command: docker run -v /home/hukang:/home/hukang kanghu/hite:3.2.0 python main.py 
 # --genome /home/hukang/HiTE/demo/genome.fa 
 # --thread 40 
 # --outdir /home/hukang/HiTE/demo/test/
```
For those unable to download images from Docker Hub, we have uploaded the Docker and Singularity images to Zenodo: [https://zenodo.org/records/14013132](https://zenodo.org/records/14013132).
```sh
# Load the Docker image
docker load -i hite_docker_3.2.0.tar
```

### <a name="install_nextflow"></a>Option 4. Run with nextflow
Nextflow is built on top of the popular programming language, Groovy, and supports the execution of workflows 
on a wide range of computing environments, including **local machines, clusters, cloud platforms, and HPC** systems.
It also provides advanced features such as **data provenance tracking, automatic parallelization, error handling**, 
and support for containerization technologies like **Docker** and **Singularity**.

We provide a [tutorial](https://github.com/CSU-KangHu/HiTE/wiki/Run-HiTE-with-Nextflow) on how to run HiTE with nextflow.


## <a name="demo"></a>Demo data

Check `HiTE/demo/genome.fa` for demo genome assembly, and run HiTE with demo data (e.g., Conda mode):
```sh
python ${pathTo/HiTE}/main.py \
 --genome ${pathTo/genome.fa} \
 --thread 40 \
 --outdir ${outdir}

 # e.g., my command: python /home/hukang/HiTE/main.py 
 # --genome /home/hukang/HiTE/demo/genome.fa 
 # --thread 40 
 # --outdir /home/hukang/HiTE/demo/test/
```

If the following files exist in the **demo/test** directory, it means the program runs successfully: 
```text
demo/test/
├── confident_helitron.fa
├── confident_other.fa
├── confident_non_ltr.fa
├── confident_tir.fa
├── confident_ltr_cut.fa.cons
└── confident_TE.cons.fa
```
Click on [Outputs](#outputs) for further details.

Note:
To avoid automatic deletion of files, set the output path parameter ```--outdir``` to an empty directory.

### Predicting conserved protein domains in TEs
To predict conserved protein domains in TEs, you need to add `--domain 1` parameter.

The output file is **confident_TE.cons.fa.domain**, which is shown as follows:
```sh
TE_name domain_name     TE_start        TE_end  domain_start    domain_end

N_111   Gypsy-50_SB_1p#LTR/Gypsy        164     4387    1       1410
...
```

### Obtaining full-length TE annotations
To obtain full-length TE annotations on the genome, you need to include the parameter `--intact_anno 1`.

The output file is **HiTE_intact.sorted.gff3**, which is shown as follows:
```sh
##gff-version 3
##date 2024-04-29 03:20:00 UTC
##ltr_identity: Sequence identity (0-1) between the left and right LTR region.
##tir_identity: Sequence identity (0-1) between the left and right TIR region.
##tsd: target site duplication
chr_0   HiTE    TIR     68394   68599   .       +       .       id=te_intact_352;name=TIR_89;classification=DNA/MULE;tir=1-68,129-206;tir_identity=0.882353;tsd=TA;tsd_len=2
chr_0   HiTE    Helitron        3534305 3534481 .       -       .       id=te_intact_517;name=Helitron_1;classification=RC/Helitron;hairpin_loop=GCGCCGAAGGCGC
chr_1   HiTE    Non_LTR 1036    1315    .       +       .       id=te_intact_472;name=Denovo_Non_LTR_0;classification=LINE/L1;polya_t=AAAAAA;tsd=AAAATTGA;tsd_len=8
chr_1   HiTE    repeat_region   139832  146291  .       -       .       id=repeat_region_1;name=chr_1:139837..146286;classification=LTR/Copia;ltr_identity=1.0000;motif=TGCA;tsd=GTATA
chr_1   HiTE    target_site_duplication 139832  139836  .       -       .       id=lTSD_1;parent=repeat_region_1;name=chr_1:139837..146286;classification=LTR/Copia;ltr_identity=1.0000;motif=TGCA;tsd=GTATA
chr_1   HiTE    long_terminal_repeat    139837  140815  .       -       .       id=lLTR_1;parent=repeat_region_1;name=chr_1:139837..146286;classification=LTR/Copia;ltr_identity=1.0000;motif=TGCA;tsd=GTATA
chr_1   HiTE    LTR     139837  146286  .       -       .       id=LTRRT_1;parent=repeat_region_1;name=chr_1:139837..146286;classification=LTR/Copia;ltr_identity=1.0000;motif=TGCA;tsd=GTATA
chr_1   HiTE    long_terminal_repeat    145317  146286  .       -       .       id=rLTR_1;parent=repeat_region_1;name=chr_1:139837..146286;classification=LTR/Copia;ltr_identity=1.0000;motif=TGCA;tsd=GTATA
chr_1   HiTE    target_site_duplication 146287  146291  .       -       .       id=rTSD_1;parent=repeat_region_1;name=chr_1:139837..146286;classification=LTR/Copia;ltr_identity=1.0000;motif=TGCA;tsd=GTATA
```

## <a name="inputs"></a>Inputs
**Required Parameters:**
* `--genome`. HiTE works with genome assemblies in **fasta**, **fa**, and **fna** formats using the `--genome` parameter.

**Useful Parameters:**
* `--curated_lib`. HiTE supports users providing a fully trusted curated library, which will be used to pre-mask highly homologous sequences in the genome, thereby reducing the computational load to some extent. We recommend using TE libraries from Repbase.
* `--annotate`. Use the TE library generated by HiTE to annotate the genome. This will produce annotation files such as `HiTE.out`, `HiTE.gff`, and `HiTE.tbl`. To generate more detailed information on genome annotation proportions, please refer to [https://github.com/CSU-KangHu/HiTE/issues/7](https://github.com/CSU-KangHu/HiTE/issues/7).

For other optional parameters, please refer to [Usage](#cmd).

## <a name="outputs"></a>Outputs
HiTE outputs many temporary files, which allow you to quickly restore the previous 
running state (use `--recover 1`) in case of any interruption during the running process. If
the pipeline completes successfully, the output directory should look like the following:
```shell
output_dir/
├── longest_repeats_*.fa
├── longest_repeats_*.flanked.fa
├── confident_tir_*.fa
├── confident_helitron_*.fa
├── confident_non_ltr_*.fa
├── confident_other_*.fa
├── confident_ltr_cut.fa.cons
├── confident_TE.cons.fa
├── HiTE.out (require `--annotate 1`)
├── HiTE.gff (require `--annotate 1`)
└── HiTE.tbl (require `--annotate 1`)
```

1. **confident_TE.cons.fa** are the classified TE libraries generated by HiTE, which can be used directly as TE library in RepeatMasker by `-lib`.
2. **longest_repeats_*.fa** represents the output of the FMEA algorithm, while **longest_repeats_*.flanked.fa** extends the sequences at both ends of **longest_repeats_*.fa**.
3. **confident_tir_*.fa**, **confident_helitron_*.fa**, **confident_non_ltr_*.fa** represent the identification results of the TIR, Helitron, and non-LTR modules in HiTE respectively, while **confident_other_*.fa** indicates the identification results of the homology-based non-LTR searching module. 
4. Note that "*" represents the number of blocks that the genome is divided into.
For example, if the genome input is 400 MB and the chunk size input is set to 100,
then * is equal to 4 (400/100), and you can find 4 files: repeats_0.fa, repeats_1.fa,
repeats_2.fa, and repeats_3.fa in your output directory. 
5. The **HiTE.out**, **HiTE.gff**, and **HiTE.tbl** files are generated using parameter `--annotate 1`. 
The **HiTE.out** and **HiTE.gff**, are genome annotation files, with **HiTE.gff** being visualizable 
in the IGV (Integrative Genomics Viewer). Additionally, **HiTE.tbl** offers statistical information 
on the proportion of each transposon type within the genome.

## <a name="code"></a>Code Structure
The code structure of HiTE is organized as follows:
```shell
Pipeline: main.py
    ├──LTR: judge_LTR_transposons.py
    ├──Homology-Non-LTR: judge_Other_transposons.py
    ├──split genome into chunks: split_genome_chunks.py
      ├──De novo TE searching: coarse_boundary.py
      ├──TIR: judge_TIR_transposons.py
      ├──Helitron: judge_Helitron_transposons.py
      └──De novo-Non-LTR: judge_Non_LTR_transposons.py
    ├──generate TE library: get_nonRedundant_lib.py
      └──unwrap nested TE: remove_nested_lib.py
    ├──genome annotation: annotate_genome.py
    ├──benchmarking reproduction: benchmarking.py
    └──clean temporary files: clean_lib.py
```

## <a name="cmd"></a>Usage
Type `python main.py -h` for help.
```
The simplest command:
python main.py --genome $genome_assembly --outdir $output_dir

Most frequently used commands:
python main.py --genome $genome_assembly --outdir $output_dir --thread 40 --plant 0 --recover 1 --annotate 1

usage: main.py [-h] --genome genome --outdir output_dir [--thread thread_num] [--chunk_size chunk_size] [--miu miu] [--plant is_plant] [--te_type te_type] [--curated_lib curated_lib]
               [--remove_nested is_remove_nested] [--domain is_domain] [--recover is_recover] [--annotate is_annotate] [--intact_anno intact_anno] [--search_struct search_struct] [--BM_RM2 BM_RM2]
               [--BM_EDTA BM_EDTA] [--BM_HiTE BM_HiTE] [--EDTA_home EDTA_home] [--coverage_threshold coverage_threshold] [--species species] [--skip_HiTE skip_HiTE] [--is_denovo_nonltr is_denovo_nonltr]
               [--debug is_debug] [--use_NeuralTE use_NeuralTE] [--is_wicker is_wicker] [--flanking_len flanking_len] [--fixed_extend_base_threshold fixed_extend_base_threshold]
               [--tandem_region_cutoff tandem_region_cutoff] [--max_repeat_len max_repeat_len] [--chrom_seg_length chrom_seg_length]

########################## HiTE, version 3.2 ##########################

optional arguments:
  -h, --help            show this help message and exit
  --genome genome       Input genome assembly path
  --outdir output_dir   The path of output directory; It is recommended to use a new directory to avoid automatic deletion of important files.
  --thread thread_num   Input thread num, default = [ 40 ]
  --chunk_size chunk_size
                        The chunk size of genome, default = [ 400 MB ]
  --miu miu             The neutral mutation rate (per bp per ya), default = [ 1.3e-08 ]
  --plant is_plant      Is it a plant genome, 1: true, 0: false. default = [ 1 ]
  --te_type te_type     Retrieve specific type of TE output [ltr|tir|helitron|non-ltr|all]. default = [ all ]
  --curated_lib curated_lib
                        Provide a fully trusted curated library, which will be used to pre-mask highly homologous sequences in the genome. We recommend using TE libraries from Repbase. default = [ None ]
  --remove_nested is_remove_nested
                        Whether to remove nested TE, 1: true, 0: false. default = [ 1 ]
  --domain is_domain    Whether to obtain TE domains, HiTE uses RepeatPeps.lib from RepeatMasker to obtain TE domains, 1: true, 0: false. default = [ 0 ]
  --recover is_recover  Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false. default = [ 0 ]
  --annotate is_annotate
                        Whether to annotate the genome using the TE library generated, 1: true, 0: false. default = [ 0 ]
  --intact_anno intact_anno
                        Whether to generate annotation of full-length TEs, 1: true, 0: false. default = [ 0 ]
  --search_struct search_struct
                        Is the structural information of full-length copies being searched, 1: true, 0: false. default = [ 1 ]
  --BM_RM2 BM_RM2       Whether to conduct benchmarking of RepeatModeler2, 1: true, 0: false. default = [ 0 ]
  --BM_EDTA BM_EDTA     Whether to conduct benchmarking of EDTA, 1: true, 0: false. default = [ 0 ]
  --BM_HiTE BM_HiTE     Whether to conduct benchmarking of HiTE, 1: true, 0: false. default = [ 0 ]
  --EDTA_home EDTA_home
                        When conducting benchmarking of EDTA, you will be asked to input EDTA home path.
  --coverage_threshold coverage_threshold
                        The coverage threshold of benchmarking methods.
  --species species     Which species you want to conduct benchmarking, six species support (dmel, rice, cb, zebrafish, maize, ath).
  --skip_HiTE skip_HiTE
                        Whether to skip_HiTE, 1: true, 0: false. default = [ 0 ]
  --is_denovo_nonltr is_denovo_nonltr
                        Whether to detect non-ltr de novo, 1: true, 0: false. default = [ 1 ]
  --debug is_debug      Open debug mode, and temporary files will be kept, 1: true, 0: false. default = [ 0 ]
  --use_NeuralTE use_NeuralTE
                        Whether to use NeuralTE to classify TEs, 1: true, 0: false. default = [1 ]
  --is_wicker is_wicker
                        Use Wicker or RepeatMasker classification labels, 1: Wicker, 0: RepeatMasker. default = [ 0 ]
  --flanking_len flanking_len
                        The flanking length of candidates to find the true boundaries, default = [ 50 ]
  --fixed_extend_base_threshold fixed_extend_base_threshold
                        The length of variation can be tolerated during pairwise alignment, default = [ 1000 ]
  --tandem_region_cutoff tandem_region_cutoff
                        Cutoff of the candidates regarded as tandem region, default = [ 0.5 ]
  --max_repeat_len max_repeat_len
                        The maximum length of a single repeat, default = [ 30000 ]
  --chrom_seg_length chrom_seg_length
                        The length of genome segments, default = [ 100000 ]
```

## <a name="ER"></a>Experiment reproduction
The quantitative experimental results from the HiTE paper can be reproduced following the [Experiment reproduction](https://github.com/CSU-KangHu/HiTE/wiki/Experiment-reproduction).
### Benchmarking method of HiTE (BM_HiTE)
```sh
# run BM_HiTE
cd HiTE && python module/lib_evaluation.py -g ${genome} \
 --standard_lib ${standard_lib} \
 --test_lib ${test_lib} \
 --work_dir ${out_dir} \
 --coverage_threshold [0.8/0.95/0.99] \
 --cat Total
```

## <a name="QA"></a>More tutorials
You may want to check out this [Wiki](https://github.com/CSU-KangHu/HiTE/wiki) page for more tutorials.
* [Issues with installation and usage](https://github.com/CSU-KangHu/HiTE/wiki/Issues-with-installation-and-usage)
* [How to make HiTE into a Docker image](https://github.com/CSU-KangHu/HiTE/wiki/How-to-make-HiTE-into-a-Docker-image)
* [Run HiTE with Nextflow](https://github.com/CSU-KangHu/HiTE/wiki/Run-HiTE-with-Nextflow)

## Citations
Please cite our paper if you find `HiTE` useful:

Hu, K., Ni, P., Xu, M. et al. HiTE: a fast and accurate dynamic boundary adjustment approach for full-length transposable element detection and annotation. Nat Commun 15, 5573 (2024). [https://doi.org/10.1038/s41467-024-49912-8](https://doi.org/10.1038/s41467-024-49912-8)