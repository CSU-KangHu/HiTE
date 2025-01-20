#!/usr/bin/env python
import argparse
import json
import os
import shutil
import subprocess
import sys
import uuid
from datetime import datetime
current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")

from Util import Logger, get_TE_class, get_full_length_copies_from_gff_v1, copy_files, create_or_clear_directory


def for_test(output_dir, threads, panTE_lib, reference, genome_name, log):
    full_length_copies_path = f'{output_dir}/{genome_name}.full_length.copies'
    full_length_gff_path = f'{output_dir}/{genome_name}.full_length.gff'
    gff_path = f'{output_dir}/{genome_name}.gff'
    os.system('touch ' + full_length_copies_path)
    os.system('touch ' + full_length_gff_path)
    os.system('touch ' + gff_path)

def run_repeat_masker(output_dir, threads, panTE_lib, reference, annotate, genome_name, log):
    if annotate == 1:
        RepeatMasker_command = f'RepeatMasker -e ncbi -no_is -norna -nolow -pa {threads} -gff -lib {panTE_lib} -cutoff 225 {reference}'
        log.logger.info(f"Running command: {RepeatMasker_command}")
        os.system(RepeatMasker_command)

        # 移动 RepeatMasker 生成的结果文件
        mv_file_command = f'cp -f {reference}.out.gff {output_dir}/{genome_name}.gff && cp -f {reference}.tbl {output_dir}/{genome_name}.tbl && cp -f {reference}.out {output_dir}/{genome_name}.out'
        log.logger.info(f"Running command: {mv_file_command}")
        os.system(mv_file_command)

        # 生成全长TE注释文件
        # 获取 TE_name 和 TE_class的对应关系
        te_classes, new_te_contigs = get_TE_class(panTE_lib)
        full_length_threshold = 0.95
        TE_gff = f'{output_dir}/{genome_name}.gff'
        full_length_copies, flank_full_length_copies, full_length_annotations = get_full_length_copies_from_gff_v1(panTE_lib, reference, TE_gff, full_length_threshold, te_classes)
        # 存储 full_length_copies
        full_length_copies_path = f'{output_dir}/{genome_name}.full_length.copies'
        with open(full_length_copies_path, 'w') as f:
            json.dump(full_length_copies, f)

        full_length_tmp_gff_path = f'{output_dir}/{genome_name}.full_length.tmp.gff'
        full_length_tmp_sort_gff_path = f'{output_dir}/{genome_name}.full_length.tmp.sort.gff'
        full_length_gff_path = f'{output_dir}/{genome_name}.full_length.gff'
        intact_count = 0
        with open(full_length_tmp_gff_path, 'w') as f_save:
            for query_name in full_length_annotations.keys():
                for copy_annotation in full_length_annotations[query_name]:
                    classification = copy_annotation[0]
                    chr_name = copy_annotation[1]
                    chr_start = str(copy_annotation[2] + 1)
                    chr_end = str(copy_annotation[3])
                    intact_count += 1
                    update_annotation = 'id=te_intact_' + str(
                        intact_count) + ';name=' + query_name + ';classification=' + classification
                    f_save.write(chr_name + '\t' + 'HiTE' + '\t' + classification + '\t' + chr_start + '\t' + chr_end + '\t' + '.\t' + str(copy_annotation[4]) + '\t' + '.\t' + update_annotation + '\n')
        os.system('sort -k1,1 -k4n ' + full_length_tmp_gff_path + ' > ' + full_length_tmp_sort_gff_path)
        gff_lines = []
        with open(full_length_tmp_sort_gff_path, 'r') as f_r:
            for line in f_r:
                if line.startswith('#'):
                    continue
                gff_lines.append(line)

        date = datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
        header = (
            "##gff-version 3\n"
            f"##date {date}\n"
        )
        with open(full_length_gff_path, "w") as gff_file:
            gff_file.write(header)
            for line in gff_lines:
                gff_file.write(line)
        os.remove(full_length_tmp_gff_path)
        os.remove(full_length_tmp_sort_gff_path)

if __name__ == "__main__":
    # 创建解析器
    parser = argparse.ArgumentParser(description="panHiTE annotate genome.")
    parser.add_argument("--panTE_lib", type=str, help="pan TE lib.")
    parser.add_argument("--reference", type=str, help="the reference path.")
    parser.add_argument("--genome_name", type=str, help="genome name.")
    parser.add_argument('--annotate', type=int, default=1,
                        help='Whether to annotate the genome using the TE library generated, 1: true, 0: false.')
    parser.add_argument("--threads", type=int, help="Number of threads to use.")
    parser.add_argument("--output_dir", nargs="?", default=os.getcwd(),
                        help="Output directory (default: current working directory).")

    # 解析参数
    args = parser.parse_args()
    panTE_lib = args.panTE_lib
    reference = args.reference
    genome_name = args.genome_name
    annotate = args.annotate
    threads = args.threads

    # 处理输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    log = Logger(output_dir + '/panHiTE.log', level='debug')

    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = '/tmp/pan_annotate_genome_' + str(unique_id)
    create_or_clear_directory(temp_dir)
    # 调用 RepeatMasker 函数
    run_repeat_masker(temp_dir, threads, panTE_lib, reference, annotate, genome_name, log)

    # 计算完之后将结果拷贝回输出目录
    copy_files(temp_dir, output_dir)

    # 删除临时目录
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
