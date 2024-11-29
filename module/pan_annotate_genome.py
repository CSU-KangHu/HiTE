#!/usr/bin/env python
import json
import os
import sys

current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")

from Util import Logger


def run_repeat_masker(result_file, output_dir, threads, panTE_lib, reference, genome_name, recover, log):
    # 如果文件不存在或需要恢复，运行 RepeatMasker
    if not recover or not os.path.exists(result_file):
        RepeatMasker_command = f'cd {output_dir} && RepeatMasker -e ncbi -no_is -norna -nolow -pa {threads} -gff -lib {panTE_lib} -cutoff 225 {reference}'
        log.logger.info(f"Running command: {RepeatMasker_command}")
        os.system(RepeatMasker_command)

        # 移动 RepeatMasker 生成的结果文件
        mv_file_command = f'mv {reference}.out {output_dir}/{genome_name}.out && mv {reference}.tbl {output_dir}/{genome_name}.tbl && mv {reference}.out.gff {output_dir}/{genome_name}.gff'
        log.logger.info(f"Running command: {mv_file_command}")
        os.system(mv_file_command)
    else:
        log.logger.info(f'{result_file} exists, skipping RepeatMasker.')

    return {
        "genome_name": f"{output_dir}/{genome_name}.gff"
    }

if __name__ == "__main__":
    # 使用 sys.argv 获取命令行参数
    if len(sys.argv) != 8:
        print("Usage: python pan_annotate_genome.py <output_dir> <threads> <panTE_lib> <reference> <genome_name> <recover>")
        sys.exit(1)

    result_file = sys.argv[1]
    output_dir = sys.argv[2]
    threads = int(sys.argv[3])
    panTE_lib = sys.argv[4]
    reference = sys.argv[5]
    genome_name = sys.argv[6]
    recover = int(sys.argv[7]) # 0 为 False，1 为 True

    if output_dir is None:
        output_dir = os.getcwd()
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    log = Logger(output_dir + '/panHiTE.log', level='debug')

    # 调用 RepeatMasker 函数
    result = run_repeat_masker(result_file, output_dir, threads, panTE_lib, reference, genome_name, recover, log)

    # 保存结果到 JSON
    result_path = os.path.join(output_dir, f"{genome_name}_annotation.json")
    with open(result_path, 'w') as f:
        json.dump(result, f, indent=4)
    log.logger.info(f"Result saved to: {result_path}")
