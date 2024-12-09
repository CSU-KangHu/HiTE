#!/usr/bin/env python
import argparse
import os
import sys
import json
current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")
from Util import Logger, generate_bam_for_RNA_seq, quantitative_gene


if __name__ == "__main__":
    # 创建解析器
    parser = argparse.ArgumentParser(description="panHiTE detect de genes.")
    parser.add_argument("--genome_info", type=str, help="genome info json.")
    parser.add_argument("--RNA_dir", type=str, help="RNA sequence data directory.")
    parser.add_argument("--threads", type=int, help="Number of threads to use.")
    parser.add_argument("--output_dir", nargs="?", default=os.getcwd(),
                        help="Output directory (default: current working directory).")

    # 解析参数
    args = parser.parse_args()
    genome_info = args.genome_info
    RNA_dir = args.RNA_dir
    threads = args.threads

    # 处理输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    log = Logger(output_dir + '/detect_de_genes.log', level='debug')

    # Load the metadata
    with open(genome_info, 'r') as f:
        genome_info_list = json.load(f)


    # Step 7.1: 生成 BAM 文件
    log.logger.info("Start generating BAM files for RNA-seq data...")
    new_batch_files = generate_bam_for_RNA_seq(genome_info_list, threads, RNA_dir, log)
    log.logger.info("BAM file generation completed.")

    # Step 7.2: 基因定量
    log.logger.info("Start gene quantification using featureCounts...")
    gene_express_dir = os.path.join(output_dir, 'gene_quantities')
    os.makedirs(gene_express_dir, exist_ok=True)
    RNA_tool_dir = os.path.join(project_dir, 'RNA_seq')
    gene_express_table = quantitative_gene(genome_info_list, RNA_tool_dir, gene_express_dir, threads, log)
    log.logger.info(f"Gene quantification completed. Results saved to {gene_express_table}.")


