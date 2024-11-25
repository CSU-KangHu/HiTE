import argparse
import copy
import os
import sys


current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

from configs import config
from Util import multi_process_align_blastx, read_fasta, get_overlap_len, merge_overlap_seq, store_fasta, get_domain_info

def merge_intervals(intervals):
    """
    合并重叠的区间。

    :param intervals: 包含多个起始和结束位置的列表，格式 [(start1, end1), (start2, end2), ...]
                      注意，位置是1基的，且区间包含start，但不包含end。
    :return: 合并后的区间列表
    """
    # 按起始位置排序区间
    intervals.sort(key=lambda x: x[0])

    merged = []
    for interval in intervals:
        # 如果 merged 为空或当前区间不与上一个区间重叠，直接添加
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
        else:
            # 否则，合并当前区间与上一个区间
            merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))

    return merged


def remove_segments(sequence, positions):
    positions = [(start - 1, end) for start, end in positions]

    # 排序 positions 以按起始位置进行处理
    positions = sorted(positions)

    # 初始化新的序列片段
    new_sequence = []
    current_pos = 0

    for start, end in positions:
        # 确保位置在序列范围内
        if start < len(sequence):
            # 添加当前段的前一部分
            if current_pos < start:
                new_sequence.append(sequence[current_pos:start])
            # 更新当前位置
            current_pos = min(end, len(sequence))

    # 添加最后一个片段后的部分
    if current_pos < len(sequence):
        new_sequence.append(sequence[current_pos:])

    return ''.join(new_sequence)


def purge_internal_seq(query_path, clean_query_path, blast_other):
    query_records = {}
    with open(blast_other, 'r') as f_r:
        for idx, line in enumerate(f_r):
            # print('current line idx: %d' % (idx))
            parts = line.split('\t')
            query_name = parts[0]
            identity = float(parts[2])
            alignment_len = int(parts[3])
            q_start = int(parts[6])
            q_end = int(parts[7])
            if q_start > q_end:
                tmp = q_start
                q_start = q_end
                q_end = tmp
            if not query_records.__contains__(query_name):
                query_records[query_name] = []
            subject_pos = query_records[query_name]
            subject_pos.append((q_start, q_end))
    f_r.close()

    query_names, query_contigs = read_fasta(query_path)
    clean_query_contigs = copy.deepcopy(query_contigs)
    for idx, query_name in enumerate(query_records.keys()):
        subject_pos = query_records[query_name]
        if query_name not in query_contigs:
            continue
        query_seq = query_contigs[query_name]
        connected_segs = merge_intervals(subject_pos)
        # print(query_name, connected_segs)
        clean_query_seq = remove_segments(query_seq, connected_segs)
        clean_query_contigs[query_name] = clean_query_seq
    store_fasta(clean_query_contigs, clean_query_path)

def purge_internal_seq_by_table(query_path, subject_path, clean_query_path, table):
    subject_names, subject_contigs = read_fasta(subject_path)
    query_records = {}
    with open(table, 'r') as f_r:
        for idx, line in enumerate(f_r):
            if idx <= 1:
                continue
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            subject_len = len(subject_contigs[subject_name])
            q_start = int(parts[2])
            q_end = int(parts[3])
            if q_start > q_end:
                tmp = q_start
                q_start = q_end
                q_end = tmp
            s_start = int(parts[4])
            s_end = int(parts[5])
            if float(abs(s_end-s_start))/subject_len >= 0.95:
                if not query_records.__contains__(query_name):
                    query_records[query_name] = []
                subject_pos = query_records[query_name]
                subject_pos.append((q_start, q_end))
    f_r.close()

    query_names, query_contigs = read_fasta(query_path)
    clean_query_contigs = copy.deepcopy(query_contigs)
    for idx, query_name in enumerate(query_records.keys()):
        subject_pos = query_records[query_name]
        if query_name not in query_contigs:
            continue
        query_seq = query_contigs[query_name]

        connected_segs = merge_intervals(subject_pos)
        clean_query_seq = remove_segments(query_seq, connected_segs)
        clean_query_contigs[query_name] = clean_query_seq
        # print(len(query_seq), len(clean_query_seq))
    store_fasta(clean_query_contigs, clean_query_path)


if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run clean LTR internal program')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number.')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--internal_seq', metavar='internal_seq',
                        help='The path of LTR internal sequence')

    args = parser.parse_args()

    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    internal_seq = args.internal_seq

    project_dir = config.project_dir
    src_dir = project_dir + '/src'
    tool_dir = project_dir + '/tools'

    internal_seq_filter = internal_seq + '.filter_tandem'
    # Step1. 清理全部由串联重复组成的LTR内部序列
    tandem_filter_command = 'python ' + src_dir + '/filter_tandem_repeats.py -f ' + internal_seq + ' > ' + internal_seq_filter
    os.system(tandem_filter_command)

    # Step2. 清理内部序列中包含的其余转座子蛋白质序列
    lib_dir = project_dir + '/databases'
    line_db = lib_dir + '/Tpases020812LINE'
    dna_db = lib_dir + '/Tpases020812DNA'
    other_protein_db = lib_dir + '/OtherPeps.lib'

    # line_blast_dir = tmp_output_dir + '/internal_line'
    # blastx_line = internal_seq_filter + '.line.out'
    # multi_process_align_blastx(internal_seq_filter, line_db, blastx_line, line_blast_dir, threads)

    # dna_blast_dir = tmp_output_dir + '/internal_dna'
    # blastx_dna = internal_seq_filter + '.dna.out'
    # multi_process_align_blastx(internal_seq_filter, dna_db, blastx_dna, dna_blast_dir, threads)
    #
    # plant_protein_blast_dir = tmp_output_dir + '/internal_plantP'
    # blastx_plantP = internal_seq_filter + '.plantP.out'
    # multi_process_align_blastx(internal_seq_filter, plant_protein_db, blastx_plantP, plant_protein_blast_dir, threads)

    # blast_other = internal_seq_filter + '.other.out'
    # # os.system('cat ' + blastx_line + ' ' + blastx_dna + ' ' + blastx_plantP + ' > ' + blast_other)
    # os.system('cat ' + blastx_line + ' > ' + blast_other)
    #
    # clean_internal_seq = internal_seq + '.clean'
    # purge_internal_seq(internal_seq_filter, clean_internal_seq, blast_other)

    # filter others
    clean_internal_seq = internal_seq_filter + '.noOther'
    output_table = clean_internal_seq + '.domain'
    temp_dir = tmp_output_dir + '/domain'
    get_domain_info(internal_seq_filter, other_protein_db, output_table, threads, temp_dir)
    purge_internal_seq_by_table(internal_seq_filter, other_protein_db, clean_internal_seq, output_table)
    internal_seq_filter = clean_internal_seq

    # filter LINEs
    clean_internal_seq = internal_seq_filter + '.noLINE'
    output_table = clean_internal_seq + '.domain'
    temp_dir = tmp_output_dir + '/domain'
    get_domain_info(internal_seq_filter, line_db, output_table, threads, temp_dir)
    purge_internal_seq_by_table(internal_seq_filter, line_db, clean_internal_seq, output_table)
    internal_seq_filter = clean_internal_seq

    # filter TIRs
    clean_internal_seq = internal_seq_filter + '.noTIR'
    output_table = clean_internal_seq + '.domain'
    temp_dir = tmp_output_dir + '/domain'
    get_domain_info(internal_seq_filter, dna_db, output_table, threads, temp_dir)
    purge_internal_seq_by_table(internal_seq_filter, dna_db, clean_internal_seq, output_table)
    internal_seq_filter = clean_internal_seq

    if os.path.exists(temp_dir):
        os.system('rm -rf ' + temp_dir)

    clean_internal_seq = internal_seq + '.clean'
    os.system('cp ' + internal_seq_filter + ' ' + clean_internal_seq)














