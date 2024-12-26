#!/usr/bin/env python
import argparse
import os
import sys
import json
current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")

from Util import Logger, find_gene_relation_tes, create_or_clear_directory, copy_files

if __name__ == "__main__":
    # 创建解析器
    parser = argparse.ArgumentParser(description="panHiTE find gene te relations.")
    parser.add_argument("--genome_info_json", type=str, help="genome info json.")
    parser.add_argument("--softcore_threshold", type=float, default=0.8, help="The rate value of softcore.")
    parser.add_argument("--output_dir", nargs="?", default=os.getcwd(),
                        help="Output directory (default: current working directory).")

    # 解析参数
    args = parser.parse_args()
    genome_info_json = args.genome_info_json
    softcore_threshold = args.softcore_threshold

    # 处理输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    log = Logger(output_dir + '/find_gene_relation_tes.log', level='debug')

    # Load the metadata
    with open(genome_info_json, 'r') as f:
        genome_info_list = json.load(f)

    # 创建本地临时目录，存储计算结果
    temp_dir = '/tmp/pan_gene_te_relation'
    create_or_clear_directory(temp_dir)

    # 调用 find_gene_relation_tes 函数
    log.logger.info('Start finding gene-TE relations...')
    find_gene_relation_tes(genome_info_list, temp_dir, softcore_threshold, log)

    # 计算完之后将结果拷贝回输出目录
    copy_files(temp_dir, output_dir)