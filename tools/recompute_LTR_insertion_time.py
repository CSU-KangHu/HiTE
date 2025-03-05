# -*- coding: utf-8 -*-
import argparse
import os
import shutil
import subprocess
import sys

import math

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)

from module.Util import store_fasta

def estimate_insert_time(identity, miu):
    """
    估计序列的插入时间
    :param identity: 比对相似性（0到1之间）
    :param miu: 突变率
    :return: 插入时间
    """
    d = 1 - identity
    K = -3 / 4 * math.log(1 - d * 4 / 3)
    T = K / (2 * miu)
    return T


def read_fasta(fasta_path):
    """
    读取FASTA文件
    :param fasta_path: FASTA文件路径
    :return: 序列名列表和序列字典
    """
    contignames = []
    contigs = {}
    if os.path.exists(fasta_path):
        with open(fasta_path, 'r') as rf:
            contigname = ''
            contigseq = ''
            for line in rf:
                if line.startswith('>'):
                    if contigname != '' and contigseq != '':
                        contigs[contigname] = contigseq
                        contignames.append(contigname)
                    contigname = line.strip()[1:].split(" ")[0].split('\t')[0]
                    contigseq = ''
                else:
                    contigseq += line.strip().upper()
            if contigname != '' and contigseq != '':
                contigs[contigname] = contigseq
                contignames.append(contigname)
    return contignames, contigs


def extract_ltr_terminal_candidates(blast_output, sequence_length, ltr_length_threshold=50):
    """
    从BLAST输出中提取最有可能是LTR终端的比对记录
    :param blast_output: BLAST输出内容
    :param sequence_length: 序列总长度
    :param ltr_length_threshold: LTR长度阈值
    :return: 候选比对记录列表
    """
    candidates = []
    lines = blast_output.split('\n')
    for line in lines:
        fields = line.split('\t')
        if len(fields) < 12:
            continue
        query_start = int(fields[6])
        query_end = int(fields[7])
        subject_start = int(fields[8])
        subject_end = int(fields[9])
        identity = float(fields[2]) / 100

        # 过滤掉自身的比对
        if query_start == subject_start or query_end == subject_end:
            continue

        # 检查比对是否位于序列的两端
        is_near_5prime = (query_start < ltr_length_threshold or subject_start < ltr_length_threshold)
        is_near_3prime = (query_end > sequence_length - ltr_length_threshold or
                          subject_end > sequence_length - ltr_length_threshold)

        # 比对区域到 5' 端的距离（取最小值）
        distance_to_5prime = min(query_start, subject_start)
        # 比对区域到 3' 端的距离（取最大值与序列长度的差值）
        distance_to_3prime = sequence_length - max(query_end, subject_end)
        # 总距离（取较小值）
        total_distance = distance_to_5prime + distance_to_3prime

        # 检查比对长度和相似性
        if is_near_5prime and is_near_3prime:
            candidates.append((identity, total_distance))

    return candidates


def main(intact_LTR_path, intact_LTR_list, miu, temp_dir):
    """
    主函数
    :param intact_LTR_path: intact_LTR.fa 文件路径
    :param intact_LTR_list: intact_LTR.list 文件路径
    :param miu: 突变率
    :param temp_dir: 临时文件目录
    """
    # 1. 读取 intact_LTR.fa
    intact_names, intact_contigs = read_fasta(intact_LTR_path)

    # 2. 解析 intact_LTR.list，获取LTR终端位置
    os.makedirs(temp_dir, exist_ok=True)
    new_lines = []

    with open(intact_LTR_list, 'r') as f_r:
        for cur_line in f_r:
            columns = cur_line.strip().split("\t")
            LTR_position = columns[0]  # 第1列
            if LTR_position not in intact_contigs:
                continue
            intact_LTR_seq = intact_contigs[LTR_position]
            cur_intact_path = os.path.join(temp_dir, f"{LTR_position}.fa")
            cur_intact_contigs = {LTR_position: intact_LTR_seq}
            store_fasta(cur_intact_contigs, cur_intact_path)

            # 3. 对完整的LTR进行blastn自比对
            blastn_command = f'blastn -query {cur_intact_path} -subject {cur_intact_path} -num_threads 1 -outfmt 6'
            result = subprocess.run(blastn_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, executable='/bin/bash')

            if result.returncode == 0:
                candidates = extract_ltr_terminal_candidates(result.stdout, len(intact_LTR_seq))
                if candidates:
                    # 按总距离排序，选择最佳记录
                    candidates.sort(key=lambda x: x[1])
                    identity = candidates[0][0]  # 选择最佳记录
                    insertion_time = int(estimate_insert_time(identity, float(miu)))

                    # 更新原文件的 identity 和插入时间
                    columns[7] = str(identity)
                    columns[12] = str(insertion_time)
                    new_line = '\t'.join(columns)
                    new_lines.append(new_line)

    # 4. 重新生成 intact_LTR.list
    adjust_intact_LTR_list = intact_LTR_list + ".adjust"
    with open(adjust_intact_LTR_list, 'w') as f_save:
        for line in new_lines:
            f_save.write(line + '\n')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)



if __name__ == '__main__':
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="计算LTR插入时间并调整intact_LTR.list文件")
    parser.add_argument('--intact_LTR_path', required=True, help="intact_LTR.fa 文件路径")
    parser.add_argument('--intact_LTR_list', required=True, help="intact_LTR.list 文件路径")
    parser.add_argument('--miu', required=True, type=float, help="突变率")
    args = parser.parse_args()

    # 设置临时文件目录
    temp_dir = os.path.join(os.path.dirname(args.intact_LTR_path), "temp")

    # 调用主函数
    main(args.intact_LTR_path, args.intact_LTR_list, args.miu, temp_dir)