#!/usr/bin/env python
import os
import sys
import json
current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")

from Util import Logger, find_gene_relation_tes

if __name__ == "__main__":
    # 检查命令行参数数量
    if len(sys.argv) != 4:
        print("Usage: python find_gene_relation_tes.py <genome_info_list_file> <gene_annotation_list_file> <output_dir> <recover (0 or 1)>")
        sys.exit(1)

    # 获取命令行参数
    genome_metadata = sys.argv[1]
    output_dir = sys.argv[2]
    recover = int(sys.argv[3])  # 0 为 False，1 为 True

    # 设置输出目录
    if output_dir is None:
        output_dir = os.getcwd()
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    log = Logger(output_dir + '/find_gene_relation_tes.log', level='debug')

    # Load the metadata
    with open(genome_metadata, 'r') as f:
        genome_data = json.load(f)
    genome_info_list = genome_data["genome_info"]
    gene_annotation_list = genome_data["gene_annotations"]

    # 调用 find_gene_relation_tes 函数
    log.logger.info('Start finding gene-TE relations...')
    find_gene_relation_tes(genome_info_list, gene_annotation_list, output_dir, recover, log)
    log.logger.info('Gene-TE relation analysis completed.')
