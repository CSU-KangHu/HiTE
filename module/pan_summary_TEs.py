#!/usr/bin/env python
import os
import sys
import json
current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")
from Util import Logger, summary_TEs

if __name__ == "__main__":
    # 使用 sys.argv 接收参数
    if len(sys.argv) != 7:
        print("Usage: python summary_TEs.py <genome_info_list_file> <pan_genomes_dir> <panTE_lib> <output_dir> <intact_ltr_paths_file> <recover>")
        sys.exit(1)

    # 获取命令行参数
    genome_metadata = sys.argv[1]
    pan_genomes_dir = sys.argv[2]
    panTE_lib = sys.argv[3]
    output_dir = sys.argv[4]
    intact_ltr_paths_file = sys.argv[5]
    recover = int(sys.argv[6])  # 0 为 False，1 为 True

    if output_dir is None:
        output_dir = os.getcwd()
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    log = Logger(output_dir + '/panHiTE.log', level='debug')

    # Load the metadata
    with open(genome_metadata, 'r') as f:
        genome_data = json.load(f)
    genome_info_list = genome_data["genome_info"]

    with open(intact_ltr_paths_file, 'r') as f:
        intact_ltr_paths = json.load(f)

    # 调用 summary_TEs 函数
    log.logger.info('Start analysing using TE annotation files...')
    summary_TEs(genome_info_list, pan_genomes_dir, panTE_lib, output_dir, intact_ltr_paths, recover, log)
