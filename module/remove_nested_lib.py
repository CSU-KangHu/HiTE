#!/usr/bin/env python
import argparse
import os
import sys

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, multi_process_align


def remove_self_alignment(blastnResults_path, new_blastnResults_path):
    records = []
    with open(blastnResults_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            if query_name == subject_name:
                continue
            else:
                records.append(line)
    f_r.close()

    with open(new_blastnResults_path, 'w') as f_save:
        for line in records:
            f_save.write(line)
    f_save.close()

def remove_nest(blastnResults_path, query_path, subject_path, output, coverage = 0.95, identity_threshold = 95):
    minlen = 100
    new_subject_contigs = {}
    query_names, query_contigs = read_fasta(query_path)
    subject_names, subject_contigs = read_fasta(subject_path)
    with open(blastnResults_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            line = str(line).replace('\n', '')
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            if not query_contigs.__contains__(query_name) or not subject_contigs.__contains__(subject_name):
                continue
            if query_name == subject_name:
                continue

            query_len = len(query_contigs[query_name])
            subject_len = len(subject_contigs[subject_name])
            identity = float(parts[2])
            alignment_len = int(parts[3])
            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            if s_start > s_end:
                tmp = s_start
                s_start = s_end
                s_end = tmp
            if float(alignment_len)/query_len < coverage or identity < identity_threshold:
                continue
            subject_seq = subject_contigs[subject_name]
            sbj_seq_p1 = subject_seq[0: s_start]
            sbj_seq_p2 = subject_seq[s_end:]
            sbj_seq_new = sbj_seq_p1 + sbj_seq_p2
            sbj_len_new = len(sbj_seq_new)
            if sbj_len_new >= minlen:
                if not new_subject_contigs.__contains__(subject_name):
                    new_subject_contigs[subject_name] = sbj_seq_new
                else:
                    old_sbj = new_subject_contigs[subject_name]
                    if sbj_len_new < len(old_sbj):
                        new_subject_contigs[subject_name] = sbj_seq_new
    f_r.close()

    for subject_name in subject_contigs.keys():
        if new_subject_contigs.__contains__(subject_name):
            subject_contigs[subject_name] = new_subject_contigs[subject_name]
    store_fasta(subject_contigs, output)


if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE remove nested TEs')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number.')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--max_iter_num', metavar='max_iter_num',
                        help='Maximum number of iterations to resolve nested insertions.')
    parser.add_argument('--input1', metavar='input1',
                        help='The path of input1')
    parser.add_argument('--input2', metavar='input2',
                        help='The path of input2')
    parser.add_argument('--output', metavar='output',
                        help='The path of output')

    args = parser.parse_args()

    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    max_iter_num = args.max_iter_num
    input1 = args.input1
    input2 = args.input2
    output = args.output

    default_max_iter_num = 3
    if max_iter_num is None:
        max_iter_num = default_max_iter_num
    else:
        max_iter_num = int(max_iter_num)

    iter_num = 0
    while iter_num < max_iter_num:
        # remove nested TE
        blastnResults_path = tmp_output_dir + '/rm_nested.self.out'
        confident_TE_blast_dir = tmp_output_dir + '/rm_nested_blast'
        multi_process_align(input1, input2, blastnResults_path, confident_TE_blast_dir, threads)
        clean_output = output
        remove_nest(blastnResults_path, input1, input2, clean_output, coverage=0.95, identity_threshold=95)
        input = clean_output
        iter_num += 1



