#!/usr/bin/env python
import argparse
import os
import sys
import json
current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")
from Util import Logger, summary_TEs


if __name__ == "__main__":
    # 创建解析器
    parser = argparse.ArgumentParser(description="panHiTE summary TEs.")
    parser.add_argument("--genome_info_json", type=str, help="genome info json.")
    parser.add_argument("--pan_genomes_dir", type=str, help="pan genomes directory.")
    parser.add_argument("--panTE_lib", type=str, help="panTE library.")
    parser.add_argument("--softcore_threshold", type=float, default=0.8, help="The rate value of softcore.")
    parser.add_argument("--recover", type=int, help="is recover.")
    parser.add_argument("--output_dir", nargs="?", default=os.getcwd(),
                        help="Output directory (default: current working directory).")

    # 解析参数
    args = parser.parse_args()
    genome_info_json = args.genome_info_json
    pan_genomes_dir = args.pan_genomes_dir
    panTE_lib = args.panTE_lib
    softcore_threshold = args.softcore_threshold
    recover = args.recover

    # 处理输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    log = Logger(output_dir + '/panHiTE.log', level='debug')

    # Load the metadata
    with open(genome_info_json, 'r') as f:
        genome_info_list = json.load(f)

    # 调用 summary_TEs 函数
    log.logger.info('Start analysing using TE annotation files...')
    summary_TEs(genome_info_list, panTE_lib, output_dir, softcore_threshold, recover, log)
