# HiTE Now Includes panHiTE Support!
![Version](https://img.shields.io/badge/version-1.0.0-orange)

## panHiTE
![image](https://github.com/user-attachments/assets/8773346a-9ed2-4602-ac29-3be35d675ca7)

We’re excited to announce a major update to HiTE, now featuring the powerful panHiTE functionality. This new workflow is specifically designed for population genomics, streamlining TE annotation and analysis across multiple genomes.

Key features in this release include:

1. **Enhanced LTR Detection**: We’ve replaced the original `LTR_retriever` with our new tool, `FiLTR`, which offers superior **accuracy and completeness** in LTR detection.
2. **Recovery of Low-Copy TEs**: Low-copy TEs identified in individual genomes are re-aligned to the pan-genomes, ensuring sufficient copy numbers and enabling recovery of the true TEs.
3. **TIDELs Detection**: TE-induced differential expression loci (TIDELs) contribute to individual-specific variations within populations. By integrating RNA-seq data, panHiTE can pinpoint gene expression changes associated with specific TE insertions, making it an invaluable tool for studying individual-specific traits.

## Table of Contents
- [panHiTE Tutorial](#tutorial)
  - [Tool Preparation](#tool_pre)
    - [Install Nextflow](#install_nf)
    - [Installation via Conda](#install_conda)
    - [Installation via Singularity](#install_singularity)
    - [Installation via Docker](#install_docker)
    - [Alternative Mirror Download](#bakup_images)
  - [Data Preparation](#data)
  - [Run pipeline](#run)
    - [Full Workflow](#full_workflow)
    - [Skipping TIDELs Detection Workflow](#skip_de)
    - [Only panTE detection Workflow](#only_panTE)
    - [Running on HPC Platform](#run_hpc)
    - [Checking the Output](#check_output)
    - [Checkpoint Recovery](#nextflow_restore)
  - [Usage](#cmd)
  - [Output Preview](#output_preview)

## <a name="tutorial"></a>panHiTE Tutorial

In this tutorial, we'll demonstrate how to use panHiTE with a demo data.

### <a name="tool_pre"></a>1. Tool Preparation

#### <a name="install_nf"></a>Install Nextflow
```bash
# 1. Install Nextflow
conda create -n nextflow -c conda-forge -c bioconda nextflow==24.10.3

# 2. Download HiTE
git clone https://github.com/CSU-KangHu/HiTE.git
```

#### <a name="install_conda"></a>Option 1. Installation via Conda
###### 1. Install Conda environment
```bash
# 1. Grant execution permissions for the required tools.
cd HiTE && python configure.py

# 2. Create the HiTE environment and record the environment path
source ~/.bashrc # or open a new terminal
conda env create --name HiTE -f environment.yml
conda activate HiTE
which python 
# For example, the output might be: /home/xxx/miniconda3/envs/HiTE/bin/python. 
# Record this path as it will be used later: /home/xxx/miniconda3/envs/HiTE
```

###### 2. Install Required R Packages

To install R packages, you may need to configure the CRAN mirror first. For example, to use the default CRAN mirror, add the following configuration to your `~/.Rprofile` file:

```R
options(repos = c(CRAN = "https://cloud.r-project.org/"))
```

We provide a script to install all necessary R dependencies in one step:  
```bash
cd HiTE
conda activate HiTE
Rscript RNA_seq/install_R_dependencies.R
```

If the installation completes successfully, the following R packages will be installed:  
- **argparser**  
- **tibble**  
- **dplyr**  
- **minpack.lm**  
- **readr**  
- **stringr**  
- **tidyr**  
- **Rsubread**
- **limma** 
- **edgeR**
- **ggplot2**

If the script fails to install the packages, you can also install them manually.
We strongly recommend that you follow the steps below to ensure the required R packages are correctly installed:
```sh
R
# R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"
# ...
> library(argparser)
> library(tibble)
> library(dplyr)
> library(minpack.lm)
> library(readr)
> library(stringr)
> library(tidyr)
> library(Rsubread)
> library(limma)
> library(edgeR)
> library(ggplot2)
```

#### <a name="install_singularity"></a>Option 2. Installation via Singularity
```sh
# pull singularity image (once for all). There will be a HiTE.sif file.
singularity pull HiTE.sif docker://kanghu/hite:3.3.2
```

#### <a name="install_docker"></a>Option 3. Installation via Docker
```sh
# pull docker image (once for all).
docker pull kanghu/hite:3.3.2
```

#### <a name="bakup_images"></a>Alternative Mirror Download
For those unable to download images from Docker Hub, we have uploaded the Docker and Singularity images to Zenodo: [https://zenodo.org/records/14881832](https://zenodo.org/records/14881832).
```sh
# Load the Docker image
docker load -i hite_docker_3.3.2.tar
```

### <a name="data"></a>2. Data Preparation

Download the [demo data](https://zenodo.org/records/14893504) from Zenodo.
A complete genome assembly, annotation, and RNA-seq reads data were downloaded from the publication: _Kang M, Wu H, Liu H, Liu W, Zhu M, Han Y, Liu W, Chen C, Song Y, Tan L, Yin K. *The pan-genome and local adaptation of Arabidopsis thaliana.* Nature Communications. 2023 Oct 6;14(1):6259_.

1. `pan_genomes_dir` (**Required**)  
   All genome assemblies should be stored in a single directory, which should be specified as the `pan_genomes_dir` parameter.

2. `genome_list` (**Required**)  
   A tab-delimited file with the following columns:  
   - Column 1: Genome assembly file name  
   - Column 2: Gene annotation file name (optional)  
   - Columns 3: is paired-end RNA-seq data flag
   - Columns 4 & 5: Paths to RNA-seq data (optional, single-end data in column `4, 5, ...`, paired-end data in columns `4 & 5`)  
    
    2.1. A Complete Example (_genome_list_):
    ```markdown
    # genome_name   gene_annotation_name    is_paired (1:True/0:False)      Illumina RNA-seq reads (tab-delimited)
    44.ket_10.fa    44.ket_10.gff   1       CRR624282_Chr1_f1.fq.gz CRR624282_Chr1_r2.fq.gz
    02.tibet.fa     02.tibet.gff    1       CRR624279_Chr1_f1.fq.gz CRR624279_Chr1_r2.fq.gz
    25.per_1.fa     25.per_1.gff    0       SRR748686_Chr1.fq.gz
    ```
    
    2.2 Example Without Gene Expression Analysis (_genome_list_no_RNA_):
    ```markdown
    # genome_name   gene_annotation_name
    44.ket_10.fa    44.ket_10.gff
    02.tibet.fa     02.tibet.gff
    25.per_1.fa     25.per_1.gff
    ```
    
    2.3 Example for panTE Detection and Annotation Only (_genome_list_no_RNA_no_gene_):
    ```markdown
    # genome_name
    44.ket_10.fa
    02.tibet.fa
    25.per_1.fa
    ```


3. `out_dir` (**Required**): Specify the output directory path.

4. `gene_dir` (Optional, required for analyses in example 2.1 and 2.2): Place all gene annotation files in a single directory and set this directory as the `gene_dir` parameter.  
   **Important**:
   1. Please ensure that the chromosome names in the genome assembly are consistent with those in the genome annotation GFF file to avoid errors. For example, using `CM072657.1` in the genome assembly while using `Chr1` in the annotation file is not allowed. It is recommended to use a unified naming convention, such as `Chr1`, across both files.
   2. Ensure that the gene_id in multiple gene annotation files has a **consistent name**, with the last element separated by an `underscore`.  
      For example, in the file `44.ket_10.gff`, a gene_id might be `ket_10_AT1G01010`, and in the file 02.tibet.gff, the gene_id should be `tibet_AT1G01010`.

5. `RNA_dir` (Optional, required for analysis in example 2.1): Set this as the parent directory for RNA-seq reads. The paths listed in columns 3 & 4 of the `genome_list` should be accessible via this directory.


### <a name="run"></a>3. Run pipeline
• Primary Installation Method: Conda (recommended)  
• Alternative Methods:  
  - Singularity  
  - Docker  
  Note: Container images may lag behind Conda in version updates

#### <a name="full_workflow"></a>3.1 Full Workflow
To run panHiTE from start to end, use the following command:
###### Installation via Conda
```bash
# 1. activate nextflow
conda activate nextflow

# 2. Write the execution script
# Make sure to replace the `xxx` placeholders with your actual absolute paths.
source_dir=xxx
pan_genomes_dir=xxx
gene_dir=xxx
RNA_dir=xxx
genome_list=xxx
out_dir=xxx
work_dir=xxx   # (default: /tmp on each compute node) When the default /tmp partition on compute nodes has insufficient space, specify a shared storage path that is universally accessible across all nodes in the cluster.
conda_name=/home/xxx/miniconda3/envs/HiTE  # You need to replace the previously recorded HiTE conda environment path here
cd $source_dir && /usr/bin/time -v nextflow run panHiTE.nf \
 -profile conda --conda_name ${conda_name} \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --genes_dir ${gene_dir} \
 --RNA_dir ${RNA_dir} \
 --out_dir ${out_dir} \
 --threads ${threads} \
 --work_dir ${work_dir}


# Example script:
#source_dir=/home/xxx/HiTE
#pan_genomes_dir=/home/xxx/demo/genomes
#gene_dir=/home/xxx/demo/gff_files
#RNA_dir=/home/xxx/demo/RNA_seq_files
#genome_list=/home/xxx/demo/genome_list
#out_dir=/home/xxx/demo/panHiTE_output
#cd $source_dir && /usr/bin/time -v nextflow run panHiTE.nf \
# -profile conda --conda_name /home/xxx/miniconda3/envs/HiTE \
# --pan_genomes_dir ${pan_genomes_dir} \
# --genome_list ${genome_list} \
# --genes_dir ${gene_dir} \
# --RNA_dir ${RNA_dir} \
# --out_dir ${out_dir} \
# --threads 40 \
# --miu 7e-9 
```

###### Installation via Singularity
```bash
# 1. activate nextflow
conda activate nextflow

# 2. Write the execution script
# Make sure to replace the `xxx` placeholders with your actual absolute paths.
source_dir=xxx
HiTE_sif_path=xxx
pan_genomes_dir=xxx
gene_dir=xxx
RNA_dir=xxx
genome_list=xxx
out_dir=xxx
/usr/bin/time -v nextflow run ${source_dir}/panHiTE.nf \
 -profile singularity --singularity_name ${HiTE_sif_path} \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --genes_dir ${gene_dir} \
 --RNA_dir ${RNA_dir} \
 --out_dir ${out_dir} \
 --threads ${threads}


# Example script:
#source_dir=/home/xxx/HiTE
#HiTE_sif_path=/home/xxx/hite_3.3.2.sif
#pan_genomes_dir=/home/xxx/demo/genomes
#gene_dir=/home/xxx/demo/gff_files
#RNA_dir=/home/xxx/demo/RNA_seq_files
#genome_list=/home/xxx/demo/genome_list
#out_dir=/home/xxx/demo/panHiTE_output
#threads=40
#/usr/bin/time -v nextflow run ${source_dir}/panHiTE.nf \
# -profile singularity --singularity_name ${HiTE_sif_path} \
# --pan_genomes_dir ${pan_genomes_dir} \
# --genome_list ${genome_list} \
# --genes_dir ${gene_dir} \
# --RNA_dir ${RNA_dir} \
# --out_dir ${out_dir} \
# --threads ${threads} \
# --miu 7e-9 
```

###### Installation via Docker
```bash
# 1. activate nextflow
conda activate nextflow

# 2. Write the execution script
# Make sure to replace the `xxx` placeholders with your actual absolute paths.
source_dir=xxx
pan_genomes_dir=xxx
gene_dir=xxx
RNA_dir=xxx
genome_list=xxx
out_dir=xxx
/usr/bin/time -v nextflow run ${source_dir}/panHiTE.nf \
 -profile docker --docker_name kanghu/hite:3.3.2 \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --genes_dir ${gene_dir} \
 --RNA_dir ${RNA_dir} \
 --out_dir ${out_dir} \
 --threads ${threads}


# Example script:
#source_dir=/home/xxx/HiTE
#pan_genomes_dir=/home/xxx/demo/genomes
#gene_dir=/home/xxx/demo/gff_files
#RNA_dir=/home/xxx/demo/RNA_seq_files
#genome_list=/home/xxx/demo/genome_list
#out_dir=/home/xxx/demo/panHiTE_output
#threads=40
#/usr/bin/time -v nextflow run ${source_dir}/panHiTE.nf \
# -profile docker --docker_name kanghu/hite:3.3.2 \
# --pan_genomes_dir ${pan_genomes_dir} \
# --genome_list ${genome_list} \
# --genes_dir ${gene_dir} \
# --RNA_dir ${RNA_dir} \
# --out_dir ${out_dir} \
# --threads ${threads} \
# --miu 7e-9 
```

#### <a name="skip_de"></a>3.2 Skipping TIDELs Detection Workflow

If you do not need to perform TE-induced differential expression loci (TIDELs) detection, you do not need to specify the `--RNA_dir` parameter. Additionally, the `--genome_list` input file only needs to include two columns: `genome_name` and `gene_annotation_name`. Refer to the example file `demo/genome_list_no_RNA`.  

###### Installation via Conda
```bash
source_dir=xxx
pan_genomes_dir=xxx
gene_dir=xxx
genome_list=xxx
out_dir=xxx
conda_name=/home/xxx/miniconda3/envs/HiTE  # Replace this with the path to your HiTE conda environment

cd $source_dir && /usr/bin/time -v nextflow run panHiTE.nf \
 -profile conda --conda_name ${conda_name} \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --genes_dir ${gene_dir} \
 --out_dir ${out_dir} \
 --threads ${threads}
```
###### Installation via Singularity
```bash
source_dir=xxx
HiTE_sif_path=xxx
pan_genomes_dir=xxx
gene_dir=xxx
genome_list=xxx
out_dir=xxx
/usr/bin/time -v nextflow run ${source_dir}/panHiTE.nf \
 -profile singularity --singularity_name ${HiTE_sif_path} \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --genes_dir ${gene_dir} \
 --out_dir ${out_dir} \
 --threads ${threads}
```

###### Installation via Docker
```bash
source_dir=xxx
pan_genomes_dir=xxx
gene_dir=xxx
genome_list=xxx
out_dir=xxx
/usr/bin/time -v nextflow run ${source_dir}/panHiTE.nf \
 -profile docker --docker_name kanghu/hite:3.3.2 \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --genes_dir ${gene_dir} \
 --out_dir ${out_dir} \
 --threads ${threads}
```


#### <a name="only_panTE"></a>3.3 Only panTE detection Workflow

If you only require the panTE library and TE annotation for each genome, you do not need to specify the `--RNA_dir` and `--genes_dir` parameters. Additionally, the `--genome_list` input file only needs to include one column: `genome_name`. Refer to the example file `demo/genome_list_no_RNA_no_gene`.  

###### Installation via Conda
```bash
source_dir=xxx
pan_genomes_dir=xxx
genome_list=xxx
out_dir=xxx
conda_name=/home/xxx/miniconda3/envs/HiTE  # Replace this with the path to your HiTE conda environment

cd $source_dir && /usr/bin/time -v nextflow run panHiTE.nf \
 -profile conda --conda_name ${conda_name} \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --out_dir ${out_dir} \
 --threads ${threads}
```

###### Installation via Singularity
```bash
source_dir=xxx
HiTE_sif_path=xxx
pan_genomes_dir=xxx
genome_list=xxx
out_dir=xxx
/usr/bin/time -v nextflow run ${source_dir}/panHiTE.nf \
 -profile singularity --singularity_name ${HiTE_sif_path} \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --out_dir ${out_dir} \
 --threads ${threads}
```

###### Installation via Docker
```bash
source_dir=xxx
pan_genomes_dir=xxx
genome_list=xxx
out_dir=xxx
/usr/bin/time -v nextflow run ${source_dir}/panHiTE.nf \
 -profile docker --docker_name kanghu/hite:3.3.2 \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --out_dir ${out_dir} \
 --threads ${threads}
```


#### <a name="run_hpc"></a>3.4 Running on HPC Platform
Run HiTE for each genome and annotate each genome using the panTE library. These two steps can be parallelized on the HPC platform to effectively reduce runtime. We tested this on an HPC platform managed by Slurm, and the key step is to provide the correct HPC configuration.

###### Modify the `HiTE/nextflow.config` file
```markdown
// HPC with singularity-loading, adapted from nanome
hpc { // general HPC configuration
    params {
        // HPC Slurm default parameters
        qos = 'resq'
        partition = 'resQ'
        queue = partition
        processors = 40
        memory = '100.GB'
        time = '50.h'
        ...
        account = 'xxx'
        queueSize = 5	// max number of job submit
        ...
    }
}

# Please consult your HPC platform to modify the above parameters
# If running a large and complex genome, you may want to set the `time` parameter a bit higher; otherwise, it may terminate unexpectedly due to exceeding the runtime limit.
```

###### Option 1. Run panHiTE using Conda and HPC
```bash
# Add the hpc option in the -profile section
source_dir=xxx
pan_genomes_dir=xxx
gene_dir=xxx
RNA_dir=xxx
genome_list=xxx
out_dir=xxx
conda_name=/home/xxx/miniconda3/envs/HiTE  # Here, replace with the previously recorded HiTE conda environment name
cd $source_dir && /usr/bin/time -v nextflow run panHiTE.nf \
 -profile conda,hpc --conda_name ${conda_name} \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --genes_dir ${gene_dir} \
 --RNA_dir ${RNA_dir} \
 --out_dir ${out_dir} \
 --threads ${threads}
```

###### Option 2. Run panHiTE using Singularity and HPC
```bash
# Add the hpc option in the -profile section
source_dir=xxx
HiTE_sif_path=xxx
pan_genomes_dir=xxx
gene_dir=xxx
RNA_dir=xxx
genome_list=xxx
out_dir=xxx
cd $source_dir && /usr/bin/time -v nextflow run panHiTE.nf \
 -profile singularity,hpc --singularity_name ${HiTE_sif_path} \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --genes_dir ${gene_dir} \
 --RNA_dir ${RNA_dir} \
 --out_dir ${out_dir} \
 --threads ${threads}
```

###### Option 3. Run panHiTE using Docker and HPC
```bash
# Add the hpc option in the -profile section
source_dir=xxx
pan_genomes_dir=xxx
gene_dir=xxx
RNA_dir=xxx
genome_list=xxx
out_dir=xxx
cd $source_dir && /usr/bin/time -v nextflow run panHiTE.nf \
 -profile docker,hpc --docker_name kanghu/hite:3.3.2 \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --genes_dir ${gene_dir} \
 --RNA_dir ${RNA_dir} \
 --out_dir ${out_dir} \
 --threads ${threads}
```

#### <a name="check_output"></a>3.5 Checking the Output

Please verify that the `${out_dir}/pan_run_hite_single/${genome}` directory contains all the necessary output files to ensure that HiTE has successfully run for each genome.

The expected structure of the output directory is as follows:  
```markdown
.
├── confident_helitron.fa     (1)
├── confident_ltr_cut.fa      (2)
├── confident_ltr.internal.fa (3)
├── confident_ltr.terminal.fa (4)
├── confident_non_ltr.fa      (5)
├── confident_other.fa        (6)
├── confident_TE.cons.fa      (7)
├── confident_tir.fa          (8)
├── intact_LTR.fa             (9)
├── intact_LTR.fa.classified  (10)
├── intact_LTR.list           (11)
├── helitron_low_copy.fa      (12)
├── non_ltr_low_copy.fa       (13)
└── tir_low_copy.fa           (14)
```

The complete output should include the following files:  
- Files (2), (3), (4), (9), (10), and (11) are results from the LTR module.  
- Files (8) and (14) is from the TIR module.  
- Files (1) and (12) is the output of the Helitron module.  
- Files (5), (6), and (13) are results from the non-LTR module.  
- File (7) contains the merged results of all TEs.  

If any file has a size of 0, such as `confident_tir.fa`, it indicates that HiTE did not detect any TIR elements in the genome. This could be due to one of two reasons:  
1. The genome genuinely lacks TIR elements.  
2. The TIR detection step in HiTE did not run correctly.  

To determine the cause, you can check if other similar genomes contain TIR elements. If you're unable to confirm, the safest approach is to delete the output directory for the affected genome and rerun the pipeline to ensure proper execution.

#### <a name="nextflow_restore"></a>3.6 Checkpoint Recovery

The panHiTE pipeline consists of 10 processes: 
* `pan_preprocess_genomes`
* `pan_run_hite_single`
* `pan_remove_redundancy`
* `pan_recover_low_copy_TEs`
* `pan_merge_TE_recover`
* `pan_generate_bam_for_RNA_seq`
* `pan_annotate_genomes`
* `pan_summarize_tes`
* `pan_gene_te_relation`
* `pan_detect_de_genes`

To ensure smooth execution, the pipeline stores the output files of all processes. If the program is interrupted for any reason, simply rerun the pipeline. Nextflow will automatically detect these intermediate files and skip the completed processes, making it particularly efficient for large-scale pangenome analyses.  

If you need to rerun a specific process and its downstream processes, you can delete the corresponding `process_name` directory within the output directory specified by the `--out_dir` parameter.

<div style="text-align: left;">
    <img src="https://github.com/user-attachments/assets/5b1c55b0-e2fc-4abc-8683-e930ce6b5376" alt="recovery" width="1200"/> 
</div>

### <a name="cmd"></a>4. Usage
Type `nextflow run panHiTE.nf --help` for help.
```
panHiTE - Nextflow PIPELINE (v1.0.0)
=================================
Usage:
The typical command is as follows:
nextflow run panHiTE.nf --pan_genomes_dir xxx --genome_list xxx --genes_dir xxx --RNA_dir xxx --out_dir xxx --threads 40 --skip_analyze 0 --miu 7e-9

Mandatory arguments:
  --pan_genomes_dir      A directory containing the pan-genomes
  --genome_list          A text file with genome and gene names. Each line represents a pair of genome and gene names, separated by a tab (	). The genome name is mandatory, while the gene name is optional. If a gene name is provided, the genes_dir parameter must also be specified.
  --out_dir              Output directory
General options:
  --softcore_threshold   occurrence of core_TE = num_of_genomes, softcore_threshold * num_of_genomes <= softcore_TE < num_of_genomes, 2 <= dispensable_TE < softcore_threshold * num_of_genomes, private_TE = 1. default = [ 0.8 ]
  --genes_dir            A directory containing the gene annotation files, gff format.
  --RNA_dir              A directory containing the RNA-seq files.
  --te_type              Retrieve specific type of TE output [ltr|tir|helitron|non-ltr|all]. default = [ all ]
  --threads              Input thread num. default = [ 10 ]
  --skip_analyze         Whether to skip analyze, only generate panTE library. default = [ 0 ]
  --miu                  The neutral mutation rate (per bp per ya). default = [ 1.3e-8 ]
  --work_dir             Temporary work dir. default = [ /tmp ]
```

### <a name="output_preview"></a>5. Output Preview
- **panHiTE.CorePan_fitmodel.pdf** and **panHiTE.CorePan_fitsmooth.pdf**
    - Saturation curves of Pan and Core TE family counts based on a fitted model and LOESS smoothing.
    <div style="text-align: left;">
        <img src="https://github.com/user-attachments/assets/92075ed8-d32a-41fb-8583-0e6e178ed1e4" alt="panHiTE CorePan_fitmodel" width="400"/>
    </div>

- **panTE.fa**  
A pan-TE library generated by panHiTE.

- **TE_summary.pdf**  
Statistical summaries of core, softcore, dispensable, and private TEs in the pan-genome, including:  

  - **(a) Proportions of TE families** based on full-length TE annotations:
  <div style="text-align: left; display: flex;">
    <img src="https://github.com/user-attachments/assets/7919a974-b699-4b84-a348-a5cff7d05034" alt="Full length TEs Ratio" width="400"/>
  </div>

  - **(b) Proportions of TE families** based on all TE annotations:
  <div style="text-align: left; display: flex;">
    <img src="https://github.com/user-attachments/assets/a5f4bcb5-61ba-4d17-a412-4b0212231568" alt="TEs Ratio" width="400"/>
  </div>

  - **(c) Genome coverage statistics** based on full-length TE annotations:
  <div style="text-align: left; display: flex;">
    <img src="https://github.com/user-attachments/assets/4b18cfdc-5568-47bd-9cbf-2c6ebcd3fdb2" alt="Full length Genome Coverage" width="800"/>
  </div>

  - **(d) Genome coverage statistics** based on all TE annotations:
  <div style="text-align: left; display: flex;">
    <img src="https://github.com/user-attachments/assets/ca7f5609-68ac-44b6-bf74-fd9f783a66d1" alt="Genome Coverage" width="800"/>
  </div>

  - **(e) Counts of different TE superfamilies** based on full-length TE annotations:
  <div style="text-align: left; display: flex;">
    <img src="https://github.com/user-attachments/assets/5f20d97b-b359-469b-bc3c-d02267729631" alt="Full length TE Classes Ratio" width="800"/>
  </div>

  - **(f) Counts of different TE superfamilies** based on all TE annotations:
  <div style="text-align: left; display: flex;">
    <img src="https://github.com/user-attachments/assets/566281fd-7350-40aa-a921-3e19bfd84464" alt="TE Classes Ratio" width="800"/>
  </div>

  - **(g) Genome coverage by TE superfamilies** based on full-length TE annotations:
  <div style="text-align: left; display: flex;">
    <img src="https://github.com/user-attachments/assets/f1608bd9-940f-43a0-ac86-d506c5c3898f" alt="Full length TE Classes Coverage" width="800"/>
  </div>

  - **(h) Genome coverage by TE superfamilies** based on all TE annotations:
  <div style="text-align: left; display: flex;">
    <img src="https://github.com/user-attachments/assets/2bf445d4-2db1-4932-bbe5-f4c7594313f8" alt="TE Classes Coverage" width="800"/>
  </div>

  - **(i) Insertion time distributions of Copia and Gypsy TEs** across samples (default mutation rate: 1.3e-8):
  <div style="text-align: left; display: flex;">
    <img src="https://github.com/user-attachments/assets/7b3f63b1-9968-4513-a771-d2a1ec95695c" alt="Insertion Time" width="800"/>
  </div>


- **DE_genes_from_TEs.tsv**  
This file contains TEs inserted upstream, downstream, or inside genes, which significantly alter gene expression levels across populations.  

```markdown
"Gene_name"	"Insert_type"	"fold_change"	"P_adjust_value"	"significant"	"direct"
"AT1G25155"	"Upstream"	-6.92725450753374	0.016414503776662	"Significant"	"down"
"AT2G42220"	"Upstream"	-4.45176810716045	2.32542377192524e-05	"Significant"	"down"
"AT3G23050"	"Upstream"	4.34397796164078	0.000293455145666126	"Significant"	"up"
"AT2G33850"	"Downstream"	3.94925482861562	0.00478988357806407	"Significant"	"up"
"AT2G14890"	"Downstream"	3.71530705794358	0.0152999573408698	"Significant"	"up"
"AT1G67810"	"Upstream"	-3.59258449763117	0.00455800117040596	"Significant"	"down"
"AT1G25210"	"Upstream"	-3.52314372296497	0.0495361618097377	"Significant"	"down"
"AT3G08030"	"Upstream"	-3.47714210543959	0.0067639204988128	"Significant"	"down"
"AT3G45010"	"Upstream"	3.3175163272979	0.00131016382787026	"Significant"	"up"
"AT1G10640"	"Upstream"	-3.14588823685148	0.0156025086537739	"Significant"	"down"
"AT5G65040"	"Downstream"	-3.09346141153513	0.0351020874076906	"Significant"	"down"
"AT1G68238"	"Downstream"	-3.06028547336873	0.0163103841689352	"Significant"	"down"
...
```

**_Column Descriptions_:**  
1. **`Gene_name`**: Name of the gene.  
2. **`Insert_type`**: Location of the TE insertion relative to the gene: Upstream, Inside, or Downstream.  
3. **`fold_change`**: Log2 fold change of gene expression in samples with TE insertion compared to those without insertion.  
4. **`P_adjust_value`**: Ajusted P-value. Significance was determined using the two-sample t-test with Benjamini-Hochberg p-value adjustment to control for false discovery rate.
5. **`significant`**: Indicates whether there is a significant difference in gene expression.  
6. **`direct`**: Indicates whether the TE insertion leads to upregulation or downregulation of gene expression.

The formula for fold change is:  
`logFoldChange = log2(Upstream + 1) - log2(No_Insertion + 1)`.


- **all_gene_TEs_details.tsv**  
This file details the positional relationships between TEs and genes across samples (e.g., upstream, downstream, inside, and distance), along with gene expression data.  

```markdown
"Gene_name"	"Genome_name"	"TE_name"	"Chromosome"	"TE_start"	"TE_end"	"Gene_start"	"Gene_end"	"Position"	"Species"	"Distance"	"expression"
"AT1G01020"	"12.li_of_095.fa"	"44-Helitron_0"	"Chr1"	19966	21685	9740	12084	"Upstream"	"12.li_of_095.fa"	7882	12.85
"AT1G01020"	"15.kelsterbach_2.fa"	"44-Helitron_0"	"Chr1"	20690	22409	10464	12810	"Upstream"	"15.kelsterbach_2.fa"	7880	8.29
"AT1G01020"	"01.col.fa"	"44-Helitron_0"	"Chr1"	20107	21826	9885	12227	"Upstream"	"01.col.fa"	7880	4.07
"AT1G01020"	"21.ms_0.fa"	"44-Helitron_0"	"Chr1"	20743	22462	10519	12861	"Upstream"	"21.ms_0.fa"	7882	24.16
"AT1G01020"	"41.sorbo.fa"	"44-Helitron_0"	"Chr1"	20516	22204	10259	12613	"Upstream"	"41.sorbo.fa"	7903	39.28
"AT1G01020"	"44.ket_10.fa"	"44-Helitron_0"	"Chr1"	19241	20919	8958	11367	"Upstream"	"44.ket_10.fa"	7874	12.65
"AT1G01020"	"38.dra_2.fa"	"44-Helitron_0"	"Chr1"	20412	21912	10182	12532	"Upstream"	"38.dra_2.fa"	7880	9.97
...
```

Column Descriptions:
  1. `Gene_name`: Gene name.
  2. `Genome_name`: Genome name.
  3. `TE_name`: TE name.
  4. `Chromosome`: Chromosome name.
  5. `TE_start` / `TE_end`: Start and end positions of the TE.
  6. `Gene_start` / `Gene_end`: Start and end positions of the gene.
  7. `Position`: TE position relative to the gene (e.g., upstream, downstream, inside).
  8. `Species`: Species name.
  9. `Distance`: Distance between the TE and the gene.
  10. `expression`: Gene expression level.
