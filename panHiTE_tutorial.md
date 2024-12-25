# HiTE Now Includes panHiTE Support!
![Version](https://img.shields.io/badge/version-0.0.1--beta-orange)

## panHiTE
![image](https://github.com/user-attachments/assets/63e511c2-aaef-43d8-b6e7-2ae221a3e25a)

We are excited to announce a significant update to HiTE, which now includes the powerful panHiTE functionality. This release introduces the following key features:

1. **Improved LTR Detection**: The original `LTR_retriever` has been replaced with our new tool, `FiLTR`, offering improved **accuracy and completeness** in LTR detection.
2. **Introducing panHiTE**: This new workflow is designed for population genomics, facilitating TE annotation and analysis across multiple genomes. panHiTE accepts `RNA-seq` data, enabling the detection of gene expression changes associated with **specific TE insertions**, making it ideal for exploring individual-specific traits.

---

## Table of Contents
- [panHiTE Tutorial](#tutorial)
  - [Installation via Conda](#install_conda)
  - [Install Required R Packages](#install_r)
  - [Data Preparation](#data)
  - [Running panHiTE through Nextflow](#run)
    - [Full Workflow](#full_workflow)
    - [Skipping Differential Gene Detection Workflow](#skip_de)
    - [Running on HPC Platform](#run_hpc)
    - [Checking the Output](#check_output)
    - [Checkpoint Recovery](#nextflow_restore)
  - [Usage](#cmd)
  - [Output Preview](#output_preview)
  - [Common issues](#common_issues)

## <a name="tutorial"></a>panHiTE Tutorial

In this tutorial, we'll demonstrate how to use panHiTE with a demo data.

### <a name="install_conda"></a>1. Installation via Conda
*(Docker and Singularity versions will be available once panHiTE stabilizes.)*

```bash
# 1. Install Nextflow
conda create -n nextflow -c conda-forge -c bioconda nextflow==24.10.3

# 2. Download HiTE
git clone https://github.com/CSU-KangHu/HiTE.git

# 3. Grant execution permissions for the required tools.
cd HiTE && python configure.py

# 4. Create the HiTE environment and record the environment path
source ~/.bashrc # or open a new terminal
conda env create --name HiTE -f environment.yml
conda activate HiTE
which python 
# For example, the output might be: /home/xxx/miniconda3/envs/HiTE/bin/python. 
# Record this path as it will be used later: /home/xxx/miniconda3/envs/HiTE
```

### <a name="install_r"></a>2. Install Required R Packages
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


### <a name="data"></a>3. Data Preparation

Download the [demo data](https://zenodo.org/records/14280542) from Zenodo.
A complete genome assembly, annotation, and RNA-seq reads data were downloaded from the publication: _Kang M, Wu H, Liu H, Liu W, Zhu M, Han Y, Liu W, Chen C, Song Y, Tan L, Yin K. *The pan-genome and local adaptation of Arabidopsis thaliana.* Nature Communications. 2023 Oct 6;14(1):6259_.

1. `pan_genomes_dir` (**Required**)  
   All genome assemblies should be stored in a single directory, which should be specified as the `pan_genomes_dir` parameter.

2. `genome_list` (**Required**)  
   A tab-delimited file with the following columns:  
   - Column 1: Genome assembly file name  
   - Column 2: Gene annotation file name (optional)  
   - Columns 3 & 4: Paths to RNA-seq data (optional, single-end data in column 3, paired-end data in columns 3 & 4)  
    
    2.1. A Complete Example (_genome_list_):
    ```markdown
    # genome_name   gene_annotation_name    is_paired (1:True/0:False)      Illumina RNA-seq reads (tab-delimited)
    44.ket_10.fa    44.ket_10.gff   1       CRR624282_Chr1_f1.fq.gz CRR624282_Chr1_r2.fq.gz
    02.tibet.fa     02.tibet.gff    1       CRR624279_Chr1_f1.fq.gz CRR624279_Chr1_r2.fq.gz
    25.per_1.fa     25.per_1.gff    0       SRR748686_Chr1.fq.gz    SRR748686_Chr1.copy.fq.gz
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


### <a name="run"></a>4. Running panHiTE through Nextflow

#### <a name="full_workflow"></a>4.1 Full Workflow
To run panHiTE from start to end, use the following command:  
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
conda_name=/home/xxx/miniconda3/envs/HiTE  # You need to replace the previously recorded HiTE conda environment path here
cd $source_dir && /usr/bin/time -v nextflow run panHiTE.nf \
 -profile conda --conda_name ${conda_name} \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --genes_dir ${gene_dir} \
 --RNA_dir ${RNA_dir} \
 --out_dir ${out_dir} \
 --threads ${threads}


# Example script:
#source_dir=/home/xxx/HiTE
#pan_genomes_dir=/home/xxx/ath_pan_genome/pan_genome/ath/ath_genome_and_annotation/genomes
#gene_dir=/home/xxx/ath_pan_genome/pan_genome/ath/ath_genome_and_annotation/gff_files
#RNA_dir=/home/xxx/ath_pan_genome/pan_genome/ath/ath_genome_and_annotation/RNA_seq_files
#genome_list=/home/xxx/ath_pan_genome/pan_genome/ath/ath_genome_and_annotation/genome_list
#out_dir=/home/xxx/ath_pan_genome/pan_genome/ath/panHiTE_output
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

#### <a name="skip_de"></a>4.2 Skipping Differential Gene Detection Workflow

If you only require the panTE library and TE annotation for each genome without performing differential gene detection caused by TEs, you do not need to specify the `--RNA_dir` parameter. Additionally, the `--genome_list` input file only needs to include two columns: `genome_name` and `gene_annotation_name`. Refer to the example file `demo/genome_list_no_RNA`.  

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

#### <a name="run_hpc"></a>4.3 Running on HPC Platform
Run HiTE for each genome and annotate each genome using the panTE library. These two steps can be parallelized on the HPC platform to effectively reduce runtime. We tested this on an HPC platform managed by Slurm, and the key step is to provide the correct HPC configuration.

1. Modify the `HiTE/nextflow.config` file
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
# It is recommended not to set the `queueSize` too high, as too much parallelism could lead to high disk I/O, potentially causing an unexpected termination.
```

2. Run panHiTE using conda and HPC
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
 --threads ${threads} \
 --skip_analyze 0
```

#### <a name="check_output"></a>4.4 Checking the Output

Please verify that the `${out_dir}/run_hite_single/${genome}` directory contains all the necessary output files to ensure that HiTE has successfully run for each genome.  

The complete output should include the following files:  
- Files (2), (3), (4), (9), (10), and (11) are results from the LTR module.  
- File (8) is from the TIR module.  
- File (1) is the output of the Helitron module.  
- Files (5) and (6) are results from the non-LTR module.  
- File (7) contains the merged results of all TEs.  

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
└── intact_LTR.list           (11)
```

If any file has a size of 0, such as `confident_tir.fa`, it indicates that HiTE did not detect any TIR elements in the genome. This could be due to one of two reasons:  
1. The genome genuinely lacks TIR elements.  
2. The TIR detection step in HiTE did not run correctly.  

To determine the cause, you can check if other similar genomes contain TIR elements. If you're unable to confirm, the safest approach is to delete the output directory for the affected genome and rerun the pipeline to ensure proper execution.

#### <a name="nextflow_restore"></a>4.5 Checkpoint Recovery

The panHiTE pipeline consists of 9 processes: 
* `preprocess_genomes`
* `run_hite_single`
* `merge_terminal_te`
* `pan_remove_redundancy`
* `annotate_genomes`
* `summarize_tes`
* `pan_gene_te_relation`
* `pan_generate_bam_for_RNA_seq`
* `pan_detect_de_genes`  

To ensure smooth execution, the pipeline stores the output files of all processes. If the program is interrupted for any reason, simply rerun the pipeline. Nextflow will automatically detect these intermediate files and skip the completed processes, making it particularly efficient for large-scale pangenome analyses.  

If you need to rerun a specific process and its downstream processes, you can delete the corresponding `process_name` directory within the output directory specified by the `--out_dir` parameter.

<div style="text-align: left;">
    <img src="https://github.com/user-attachments/assets/5b1c55b0-e2fc-4abc-8683-e930ce6b5376" alt="recovery" width="1200"/> 
</div>

### <a name="cmd"></a>5. Usage
Type `nextflow run panHiTE.nf --help` for help.
```
panHiTE - Nextflow PIPELINE (v0.0.1)
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
```

### <a name="output_preview"></a>6. Output Preview

- **DE_genes_from_TEs.tsv**  
This file contains TEs inserted upstream, downstream, or inside genes, which significantly alter gene expression levels across populations.  

```markdown
"Gene_name"	"NoInsertion (mean_TPM)"	"Upstream (mean_TPM)"	"Downstream (mean_TPM)"	"Inside (mean_TPM)"	"Upstream_vs_NoInsertion_logFoldChange"	"Inside_vs_NoInsertion_logFoldChange"	"Downstream_vs_NoInsertion_logFoldChange"
"AT4G34410"	562.145	NA	1.14	NA	NA	NA	-8.03975183205617
"AT3G48360"	150.99	0.16	16.46	NA	-7.033707790821	NA	-3.12185094231913
"AT4G01950"	389.61	2.33	17.3366666666667	NA	-6.87406289658701	NA	-4.4129255767341
"AT5G45890"	137.77	0.21	NA	NA	-6.84154485479315	NA	NA
"AT4G22880"	68.44	NA	0.07	NA	NA	NA	-6.02008424604333
...
```

_Column Descriptions_:
  1. `Gene_name`: Gene name.
  2. `NoInsertion (mean_TPM)`: Mean expression (TPM) of the gene with no TE insertion.
  3. `Upstream (mean_TPM)`: Mean expression (TPM) of the gene with a TE inserted upstream.
  4. `Downstream (mean_TPM)`: Mean expression (TPM) of the gene with a TE inserted downstream.
  5. `Inside (mean_TPM)`: Mean expression (TPM) of the gene with a TE inserted inside.
  6. `Upstream_vs_NoInsertion_logFoldChange`: Log2 fold change of upstream insertion compared to no insertion.
  7. `Inside_vs_NoInsertion_logFoldChange`: Log2 fold change of inside insertion compared to no insertion.
  8. `Downstream_vs_NoInsertion_logFoldChange`: Log2 fold change of downstream insertion compared to no insertion.  
Positive fold changes indicate upregulation, while negative values indicate downregulation. The formula for fold change is:  
`Upstream_vs_NoInsertion_logFoldChange = log2(Upstream + 1) - log2(NoInsertion + 1)`.

_Result Interpretation_:

"AT4G34410" 562.145 NA 1.14 NA NA NA -8.03975183205617  

The average expression level (TPM) of gene AT4G34410 without TE insertion is 562.145. In specific samples (one or more), a TE insertion is present downstream of this gene, which results in an average expression level of 1.14 after the insertion. The log fold change is -8.03, indicating that the gene is downregulated due to the TE insertion.  
You can further examine the `all_gene_TEs_details.tsv` file to identify which samples contain the TE insertion for this gene, the distance between the TE insertion and the gene, expression levels in specific samples, and other related information.


- **all_gene_TEs_details.tsv**  
This file details the positional relationships between TEs and genes across samples (e.g., upstream, downstream, inside, and distance), along with gene expression data.  

```markdown
"Gene_name"	"Genome_name"	"TE_name"	"Chromosome"	"TE_start"	"TE_end"	"Gene_start"	"Gene_end"	"Position"	"Species"	"Distance"	"expression"
"AT1G01020"	"44.ket_10.fa"	"44-Helitron_44"	"Chr1"	19241	20919	8958	11367	"Upstream"	"44.ket_10"	7874	12.65
"AT1G01030"	"44.ket_10.fa"	"44-Helitron_44"	"Chr1"	19241	20919	13881	15947	"Upstream"	"44.ket_10"	3294	0.98
"AT1G01040"	"44.ket_10.fa"	"44-Helitron_44"	"Chr1"	19241	20919	25358	33453	"Upstream"	"44.ket_10"	4439	31.61
"AT1G01050"	"44.ket_10.fa"	"25-TIR_69"	"Chr1"	44534	44768	33396	35397	"Upstream"	"44.ket_10"	9137	91.79
"AT1G01050"	"02.tibet.fa"	"25-TIR_69"	"Chr1"	47261	47494	36147	38160	"Upstream"	"02.tibet"	9101	102.58
"AT1G01050"	"25.per_1.fa"	"25-TIR_69"	"Chr1"	44736	44971	33603	35605	"Upstream"	"25.per_1"	9131	88.57
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

- **panHiTE.CorePan_fitmodel.pdf** and **panHiTE.CorePan_fitsmooth.pdf**  
Saturation curves of Pan and Core TE family counts based on a fitted model and LOESS smoothing.

<div style="text-align: left;">
    <img src="https://github.com/user-attachments/assets/92075ed8-d32a-41fb-8583-0e6e178ed1e4" alt="panHiTE CorePan_fitmodel" width="800"/>
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


### <a name="common_issues"></a>7. Common issues
...