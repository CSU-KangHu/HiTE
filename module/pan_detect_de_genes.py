#!/usr/bin/env python
import argparse
import os
import shutil
import sys
import json
import uuid

current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")
from Util import Logger, quantitative_gene, file_exist, create_or_clear_directory, copy_files, quantitative_TE, \
    clean_old_tmp_files_by_dir

if __name__ == "__main__":
    # 创建解析器
    parser = argparse.ArgumentParser(description="panHiTE detect de genes.")
    parser.add_argument("--genome_info_for_bam_json", type=str, help="genome info json.")
    parser.add_argument("--gene_te_associations", type=str, help="gene te association file.")
    parser.add_argument("--RNA_dir", type=str, help="RNA sequence data directory.")
    parser.add_argument("--threads", type=int, help="Number of threads to use.")
    parser.add_argument("--output_dir", nargs="?", default=os.getcwd(),
                        help="Output directory (default: current working directory).")

    # 解析参数
    args = parser.parse_args()
    genome_info_for_bam_json = args.genome_info_for_bam_json
    gene_te_associations = args.gene_te_associations
    RNA_dir = args.RNA_dir
    threads = args.threads

    threads = max(1, threads - 4)

    # 处理输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    log = Logger(output_dir + '/detect_de_genes.log', level='debug')

    # Load the metadata
    with open(genome_info_for_bam_json, 'r') as f:
        genome_info_list = json.load(f)

    clean_old_tmp_files_by_dir('/tmp')

    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = '/tmp/pan_detect_de_genes_' + str(unique_id)
    try:
        create_or_clear_directory(temp_dir)

        # Step 7.1: TE定量
        log.logger.info("Start TE quantification using featureCounts...")
        TE_express_dir = os.path.join(temp_dir, 'TE_quantities')
        os.makedirs(TE_express_dir, exist_ok=True)
        RNA_tool_dir = os.path.join(project_dir, 'RNA_seq')
        TE_express_table = quantitative_TE(genome_info_list, RNA_tool_dir, TE_express_dir, temp_dir, threads, log)
        log.logger.info(f"TE quantification completed. Results saved to {TE_express_table}.")

        # Step 7.2: 基因定量
        log.logger.info("Start gene quantification using featureCounts...")
        gene_express_dir = os.path.join(temp_dir, 'gene_quantities')
        os.makedirs(gene_express_dir, exist_ok=True)
        RNA_tool_dir = os.path.join(project_dir, 'RNA_seq')
        gene_express_table = quantitative_gene(genome_info_list, RNA_tool_dir, gene_express_dir, temp_dir, threads, log)
        log.logger.info(f"Gene quantification completed. Results saved to {gene_express_table}.")

        # Step 7.3: 差异表达基因检测
        log.logger.info("Start detecting DE genes associated with LTR insertions...")
        if file_exist(gene_express_table) and file_exist(gene_te_associations):
            script_dir = os.path.join(project_dir, 'RNA_seq')
            detect_DE_genes_from_TEs_cmd = (
                f"Rscript {script_dir}/detect_DE_genes_from_TEs.R {gene_express_table} {gene_te_associations}"
            )
            log.logger.debug(detect_DE_genes_from_TEs_cmd)
            exit_code = os.system(detect_DE_genes_from_TEs_cmd)
            if exit_code == 0:
                log.logger.info("DE gene detection completed successfully.")
            else:
                log.logger.error("Error occurred during DE gene detection.")
        else:
            # 创建空的输出文件以跳过nextflow检查
            de_genes = os.path.join(temp_dir, "DE_genes_from_TEs.tsv")
            all_genes = os.path.join(temp_dir, "all_gene_TEs_details.tsv")
            os.system('touch ' + de_genes)
            os.system('touch ' + all_genes)

        # 计算完之后将结果拷贝回输出目录
        copy_files(temp_dir, output_dir)

    except Exception as e:
        # 如果出现异常，打印错误信息并删除临时目录
        print(f"An error occurred: {e}")
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        raise  # 重新抛出异常，以便上层代码可以处理

    else:
        # 如果没有异常，删除临时目录
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)