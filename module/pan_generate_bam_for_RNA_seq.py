#!/usr/bin/env python
import argparse
import os
import re
import shutil
import sys
import json
import uuid

current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")
from Util import Logger, generate_bam_for_RNA_seq, copy_files, create_or_clear_directory, clean_old_tmp_files_by_dir


def preprocess_RNA_seq(input_str):
    # 去掉方括号 []
    input_str = input_str.strip('[]')

    # 替换 'key:value' 格式为 '"key":"value"' 格式
    input_str = re.sub(r'(\w+):([^,]+)', r'"\1":"\2"', input_str)

    # 将 'true' 和 'false' 替换为 JSON 格式的 'true' 和 'false'
    input_str = input_str.replace('true', 'True').replace('false', 'False')

    # 如果有引号问题，手动调整双引号
    return f'{{{input_str}}}'

if __name__ == "__main__":
    # 创建解析器
    parser = argparse.ArgumentParser(description="panHiTE generate bam for RNA-seq reads.")
    parser.add_argument("--genome_name", type=str, help="genome name.")
    parser.add_argument("--reference", type=str, help="the reference path.")
    parser.add_argument("--RNA_seq", type=str, help="the RNA_seq json string.")
    parser.add_argument("--RNA_dir", type=str, help="RNA sequence data directory.")
    parser.add_argument("--threads", type=int, help="Number of threads to use.")
    parser.add_argument("--output_dir", nargs="?", default=os.getcwd(), help="Output directory (default: current working directory).")
    parser.add_argument('-w', '--work_dir', nargs="?", default='/tmp', help="The temporary work directory for HiTE.")

    # 解析参数
    args = parser.parse_args()
    genome_name = args.genome_name
    reference = args.reference
    RNA_dir = args.RNA_dir
    threads = args.threads
    work_dir = args.work_dir
    work_dir = os.path.abspath(work_dir)

    reference = os.path.realpath(reference)

    RNA_seq = {}  # 默认空字典

    if preprocess_RNA_seq and isinstance(preprocess_RNA_seq, str):  # 确保非空且是字符串
        try:
            RNA_seq = json.loads(preprocess_RNA_seq)
        except (json.JSONDecodeError, TypeError):  # 捕获 JSON 解析错误或类型错误
            RNA_seq = {}  # 解析失败时设为空字典

    # 处理输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    log = Logger(output_dir + '/generate_bam_for_RNA_seq.log', level='debug')

    # clean_old_tmp_files_by_dir('/tmp')

    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = os.path.join(work_dir, 'pan_generate_bam_for_RNA_seq_' + str(unique_id))
    try:
        create_or_clear_directory(temp_dir)
        # Step 7.1: 生成 BAM 文件
        log.logger.info("Start generating BAM files for RNA-seq data...")
        if 'Status' in RNA_seq:
            genome_name += '.' + RNA_seq['Status']
        generate_bam_for_RNA_seq(genome_name, reference, RNA_seq, threads, RNA_dir, temp_dir, log)
        log.logger.info("BAM file generation completed.")

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