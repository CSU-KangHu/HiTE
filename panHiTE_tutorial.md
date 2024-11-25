# HiTE Now Includes panHiTE Support!

We are excited to announce a significant update to HiTE, which now includes the powerful panHiTE functionality. This release introduces the following key features:

1. **Improved LTR Detection**: The original `LTR_retriever` has been replaced with our new tool, `FiLTR`, offering improved **accuracy and completeness** in LTR detection.
2. **Introducing panHiTE**: This new workflow is designed for population genomics, facilitating TE annotation and analysis across multiple genomes. panHiTE accepts `RNA-seq` data, enabling the detection of gene expression changes associated with **specific TE insertions**, making it ideal for exploring individual-specific traits.

---

## panHiTE Tutorial

In this tutorial, we'll demonstrate how to use panHiTE with three complete *Arabidopsis thaliana* assemblies, gene annotations, and sequencing data from the following study:
_Kang, Minghui, Haolin Wu, Huanhuan Liu, Wenyu Liu, Mingjia Zhu, Yu Han, Wei Liu et al. "The pan-genome and local adaptation of Arabidopsis thaliana." Nature Communications 14, no. 1 (2023): 6259._

### 1. Installation via Conda  
*(Docker and Singularity versions will be available once panHiTE stabilizes.)*

```bash
git clone https://github.com/CSU-KangHu/HiTE.git

cd HiTE && python configure.py

source ~/.bashrc # or open a new terminal

conda env create --name HiTE -f environment.yml
conda activate HiTE
```

---

### 2. Data Preparation  

Download the [demo data](https://zenodo.org/records/14189130) from Zenodo. Since the dataset includes RNA-seq data, it is relatively large, so please be patient.

1. **pan_genomes_dir (Required)**  
   All genome assemblies should be stored in a single directory, which should be specified as the `pan_genomes_dir` parameter.

2. **genome_list (Required)**  
   A tab-delimited file with the following columns:  
   - Column 1: Genome assembly file name  
   - Column 2: Gene annotation file name (optional)  
   - Columns 3 & 4: Paths to RNA-seq data (optional, single-end data in column 3, paired-end data in columns 3 & 4)  

   Example:  
   ![image](https://github.com/user-attachments/assets/e0ae6a3c-c7e5-44a5-ac11-8c8b47af331d)

3. **out_dir (Required)**  
   Specify the output directory path.

4. **gene_dir (Optional)**  
   Place all gene annotation files in a single directory and set this directory as the `gene_dir` parameter.  
   **Important**: Ensure that the gene_id in multiple gene annotation files has a consistent name, with the last element separated by an underscore.  
   For example, in the file 01.col.gtf, a gene_id might be `col_AT1G01010`, and in the file 02.tibet.gtf, the gene_id should be `tibet_AT1G01010`.

5. **RNA_dir (Optional)**  
   Set this as the parent directory for RNA-seq reads. The paths listed in columns 3 & 4 of the `genome_list` should be accessible via this directory.

---

### 3. Running panHiTE  

#### 3.1 Full Workflow  
To run panHiTE from start to end, use the following command:  
```bash
python panHiTE.py \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --genes_dir ${gene_dir} \
 --RNA_dir ${RNA_dir} \
 --out_dir ${out_dir} \
 --thread ${threads} \
 --recover 1 \
 --te_type all \
 --skip_analyze 0 \
 --miu ${miu}
 
# Example script:
#source_dir=/public/home/xxx/test/panHiTE
#pan_genomes_dir=/public/home/xxx/ath_pan_genome/pan_genome/ath/ath_genome_and_annotation/genomes
#gene_dir=/public/home/xxx/ath_pan_genome/pan_genome/ath/ath_genome_and_annotation/gtf_files
#RNA_dir=/public/home/xxx/ath_pan_genome/pan_genome/ath/ath_genome_and_annotation/RNA_seq_files
#genome_list=/public/home/xxx/ath_pan_genome/pan_genome/ath/ath_genome_and_annotation/genome_list
#out_dir=/public/home/xxx/ath_pan_genome/pan_genome/ath/panHiTE_serial_output
#cd $source_dir && /usr/bin/time -v python panHiTE.py \
# --pan_genomes_dir ${pan_genomes_dir} \
# --genome_list ${genome_list} \
# --genes_dir ${gene_dir} \
# --RNA_dir ${RNA_dir} \
# --out_dir ${out_dir} \
# --thread 40 \
# --recover 1 \
# --te_type all \
# --skip_analyze 0 \
# --miu 7e-9
```
- **`--recover 1`**: This option detects existing intermediate files and resumes from the last checkpoint, saving time.

#### 3.2 skip downstream analysis 
If you only need the panLTR library and the annotation for each genome, you can choose to skip the downstream analysis by setting `--skip_analyze 1`.
```bash
python panHiTE.py \
 --pan_genomes_dir ${pan_genomes_dir} \
 --genome_list ${genome_list} \
 --out_dir ${out_dir} \
 --thread ${threads} \
 --recover 1 \
 --te_type all \
 --skip_analyze 1 \
 --miu ${miu}

# Example script:
#source_dir=/public/home/xxx/test/panHiTE
#pan_genomes_dir=/public/home/xxx/ath_pan_genome/pan_genome/ath/ath_genome_and_annotation/genomes
#gene_dir=/public/home/xxx/ath_pan_genome/pan_genome/ath/ath_genome_and_annotation/gtf_files
#RNA_dir=/public/home/xxx/ath_pan_genome/pan_genome/ath/ath_genome_and_annotation/RNA_seq_files
#genome_list=/public/home/xxx/ath_pan_genome/pan_genome/ath/ath_genome_and_annotation/genome_list
#out_dir=/public/home/xxx/ath_pan_genome/pan_genome/ath/panHiTE_serial_output
#cd $source_dir && /usr/bin/time -v python panHiTE.py \
# --pan_genomes_dir ${pan_genomes_dir} \
# --genome_list ${genome_list} \
# --out_dir ${out_dir} \
# --thread 40 \
# --recover 1 \
# --te_type all \
# --skip_analyze 1 \
# --miu 7e-9
```
---

### 4. Features Under Development

The current panHiTE pipeline runs HiTE detection on multiple genomes sequentially, generates the panHiTE library, 
and then annotates each genome sequentially as well. In theory, these steps can be parallelized on an HPC platform. 
Therefore, we plan to develop an HPC-compatible panHiTE version based on Nextflow, enabling faster and more efficient panHiTE analysis.
