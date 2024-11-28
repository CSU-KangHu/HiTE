#!/usr/bin/env python
import os
import sys
import json
current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")
from Util import Logger, generate_bam_for_RNA_seq, quantitative_gene


if __name__ == "__main__":
    # 检查命令行参数数量
    if len(sys.argv) != 6:
        print("Usage: python detect_de_genes.py <genome_info_list_file> <threads> <recover (0 or 1)> <output_dir>")
        sys.exit(1)

    # 获取命令行参数
    genome_metadata = sys.argv[1]
    threads = int(sys.argv[2])
    recover = int(sys.argv[3])  # 0 为 False，1 为 True
    output_dir = sys.argv[4]
    gene_te_associations = sys.argv[5]

    # 设置日志
    if output_dir is None:
        output_dir = os.getcwd()
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    log = Logger(output_dir + '/detect_de_genes.log', level='debug')

    # Load the metadata
    with open(genome_metadata, 'r') as f:
        genome_data = json.load(f)
    genome_info_list = genome_data["genome_info"]
    gene_annotation_list = genome_data["gene_annotations"]
    is_RNA_analyze = genome_data["is_RNA_analyze"]

    if is_RNA_analyze:
        # Step 7.1: 生成 BAM 文件
        log.logger.info("Start generating BAM files for RNA-seq data...")
        RNA_seq_dir = os.path.join(project_dir, 'RNA_seq')
        new_batch_files = generate_bam_for_RNA_seq(genome_info_list, threads, recover, RNA_seq_dir, log)
        log.logger.info("BAM file generation completed.")

        # Step 7.2: 基因定量
        log.logger.info("Start gene quantification using featureCounts...")
        gene_express_dir = os.path.join(output_dir, 'gene_quantities')
        os.makedirs(gene_express_dir, exist_ok=True)
        gene_express_table = quantitative_gene(new_batch_files, gene_express_dir, threads, recover, log)
        log.logger.info(f"Gene quantification completed. Results saved to {gene_express_table}.")

        # Step 7.3: 差异表达基因检测
        log.logger.info("Start detecting DE genes associated with LTR insertions...")
        script_dir = os.path.join(project_dir, 'RNA_seq')
        detect_DE_genes_from_TEs_cmd = (
            f"cd {output_dir} && Rscript {script_dir}/detect_DE_genes_from_TEs.R {gene_express_table} {gene_te_associations}"
        )
        log.logger.debug(detect_DE_genes_from_TEs_cmd)
        exit_code = os.system(detect_DE_genes_from_TEs_cmd)
        if exit_code == 0:
            log.logger.info("DE gene detection completed successfully.")
        else:
            log.logger.error("Error occurred during DE gene detection.")
