#-- coding: UTF-8 --
import argparse
import codecs
import os
import sys

import json

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, getReverseSequence, read_fasta_v2, get_copies, flanking_copies, \
    multi_process_align, multi_process_align_and_get_copies, flanking_copies_v2, rename_fasta


def getPolyASeq(output, longest_repeats_path):
    res = {}
    contigNames, contigs = read_fasta(longest_repeats_path)
    with open(output, 'r') as f_r:
        for i, line in enumerate(f_r):
            parts = line.split('\t')
            repeat_id = parts[0]
            query_name = parts[2]
            start_pos = int(parts[3])
            end_pos = int(parts[4])
            seq = contigs[query_name]
            if start_pos > end_pos:
                #反向互补
                polyAseq = getReverseSequence(seq[end_pos-1:])
            else:
                polyAseq = seq[0: end_pos]
            new_query_name = 'N_'+str(i)
            res[new_query_name] = polyAseq
    f_r.close()
    return res

def extract_sequence_from_db(db_path, features, store_path):
    store_contigs = {}
    db_contignames, db_contigs = read_fasta(db_path)
    for name in db_contignames:
        for feature in features:
            if name.__contains__(feature):
                db_seq = db_contigs[name]
                store_contigs[name] = db_seq
    store_fasta(store_contigs, store_path)

def preprocess():
    db_path = '/public/home/hpc194701009/repeat_detect_tools/RepeatMasker-4.1.2/RepeatMasker/Libraries/RepeatMasker.lib'
    features = ['#LINE', '#SINE', 'DIRS', '#DNA/Crypton']
    store_path = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa/non_LTR.lib'
    extract_sequence_from_db(db_path, features, store_path)


def extract_xsmall(ref_masked_path):
    contigs = {}
    ref_names, ref_contigs = read_fasta_v2(ref_masked_path)
    node_index = 0
    for name in ref_names:
        seq = ref_contigs[name]
        #提取小写字母
        cur_seq = ''
        for i in range(len(seq)):
            if seq[i] >= 'a' and seq[i] <= 'z':
                cur_seq += seq[i]
            else:
                if cur_seq != '' and len(cur_seq) >= 100:
                    query_name = 'Homology_'+str(node_index)
                    contigs[query_name] = cur_seq
                    cur_seq = ''
                    node_index += 1

        if cur_seq != '' and len(cur_seq) >= 100:
            query_name = 'H_' + str(node_index)
            contigs[query_name] = cur_seq
            node_index += 1
    return contigs


if __name__ == '__main__':
    #preprocess()
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('--seqs', metavar='seqs',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa/longest_repeats_0.flanked.fa')
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/dmel')
    parser.add_argument('--query_coverage', metavar='query_coverage',
                        help='e.g., 0.95')
    parser.add_argument('--subject_coverage', metavar='subject_coverage',
                        help='e.g., 0.95')
    parser.add_argument('--ref_index', metavar='ref_index',
                        help='e.g., 0')
    parser.add_argument('--library_dir', metavar='library_dir',
                        help='e.g., ')


    args = parser.parse_args()
    longest_repeats_flanked_path = args.seqs
    reference = args.g
    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    query_coverage = float(args.query_coverage)
    subject_coverage = float(args.subject_coverage)
    ref_index = args.ref_index
    library_dir = args.library_dir

    # 除了LTR、TIR、Helitron之外的其他转座子，包括LINE、SINE、DIRS、PLE(在Dfam中属于LINE)、Crypton它们缺少或具有复杂的终端结构特点，且没有稳定的TSD特征。
    # 想要根据结构特征去识别有困难，我们根据同源性搜索的方法去识别。

    # LINE （1000-7000bp），通常以poly(A)结尾和 SINE(100-600bp)，generate TSDs (5–15 bp)，通常以poly(T)结尾，发现也有polyA结尾。我们还需要考虑反向互补序列。
    non_LTR_lib = library_dir + '/non_LTR.lib'

    other_TE_dir = tmp_output_dir + '/other_TE_' + str(ref_index)
    os.system('rm -rf ' + other_TE_dir)
    if not os.path.exists(other_TE_dir):
        os.makedirs(other_TE_dir)

    # blastnResults_path = tmp_output_dir + '/non_LTR.lib.out'
    # multi_process_align(longest_repeats_path, other_TE_lib, blastnResults_path, blast_program_dir, other_TE_dir, threads)
    # all_copies = get_copies(blastnResults_path, longest_repeats_path, other_TE_lib,
    #                         query_coverage=0.9, subject_coverage=0.9, threads=threads)


    all_copies = multi_process_align_and_get_copies(non_LTR_lib, longest_repeats_flanked_path, other_TE_dir, 'other', threads, query_coverage=query_coverage, subject_coverage=subject_coverage)

    flanking_len = 0
    all_copies = flanking_copies_v2(all_copies, non_LTR_lib, longest_repeats_flanked_path, flanking_len, copy_num=1)

    confident_other_path = tmp_output_dir + '/confident_other_' + str(ref_index) + '.fa'
    confident_other_contigs = {}
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        if len(copies) >= 1:
            seq = copies[0][4]
            confident_other_contigs[query_name] = seq
    store_fasta(confident_other_contigs, confident_other_path)

    # for test
    confident_other_rename_path = tmp_output_dir + '/confident_other_' + str(ref_index) + '.rename.fa'
    rename_fasta(confident_other_path, confident_other_rename_path)

    confident_other_rename_consensus = tmp_output_dir + '/confident_other_' + str(ref_index) + '.rename.cons.fa'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_other_rename_path + ' -o ' + confident_other_rename_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)






    # RepeatMasker_command = 'cd ' + other_TE_masked_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(threads) \
    #                        + ' -lib ' + other_TE_lib + ' -q -no_is -nolow -div 40 -xsmall -dir ' + other_TE_masked_dir + ' ' + reference
    # os.system(RepeatMasker_command)
    #
    # (reference_dir, reference_filename) = os.path.split(reference)
    # (reference_name, reference_extension) = os.path.splitext(reference_filename)
    # ref_masked_path = other_TE_masked_dir + '/' + reference_filename + '.masked'
    # contigs = extract_xsmall(ref_masked_path)
    #
    # homology_TE_path = tmp_output_dir + '/homology_TE.fa'
    # store_fasta(contigs, homology_TE_path)





    # # 步骤：1.修改polyATail，使其能够识别我们序列中的PolyA，获得候选的polyA序列，如果可以的话再实现一个polyTTail。
    # # 2.获取候选polyA序列的拷贝，
    #
    # output = tmp_output_dir + '/longest_repeats.fa.polyA.set'
    # polyA_temp_dir = tmp_output_dir + '/polyA_temp'
    # multi_process_polyATail(longest_repeats_path, output, polyA_temp_dir, TRsearch_dir)
    # candidate_polyA_path = tmp_output_dir + '/candidate_polyA.fa'
    # polyA_contigs = getPolyASeq(output, longest_repeats_path)
    # store_fasta(polyA_contigs, candidate_polyA_path)
    #
    # blastnResults_path = tmp_output_dir + '/polyA.ref.out'
    # blast_program_dir = '/public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2'
    # polyA_blast_dir = tmp_output_dir + '/polyA_blast'
    # multi_process_align(candidate_polyA_path, reference, blastnResults_path, blast_program_dir, polyA_blast_dir, threads)
    # all_copies = get_copies(blastnResults_path, candidate_polyA_path)
    #
    # flanking_len = 50
    # # 过滤掉拷贝数小于2, flanking copies
    # ref_names, ref_contigs = read_fasta(reference)
    # new_all_copies = {}
    # for query_name in all_copies.keys():
    #     copies = all_copies[query_name]
    #     if len(copies) < 2:
    #         continue
    #     for copy in copies:
    #         ref_name = copy[0]
    #         copy_ref_start = int(copy[1])
    #         copy_ref_end = int(copy[2])
    #         direct = copy[4]
    #         copy_len = copy_ref_end - copy_ref_start + 1
    #         copy_seq = ref_contigs[ref_name][copy_ref_start - 1 - flanking_len: copy_ref_end + flanking_len]
    #         if direct == '-':
    #             copy_seq = getReverseSequence(copy_seq)
    #         if len(copy_seq) < 100:
    #             continue
    #         if not new_all_copies.__contains__(query_name):
    #             new_all_copies[query_name] = []
    #         copy_list = new_all_copies[query_name]
    #         copy_list.append((ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq))
    #
    # # 如果一条序列是转座子，那么它的拷贝两端序列是没有相似性的。
    # polyA_flank_blast_dir = tmp_output_dir + '/polyA_flank_blast'
    # filter_flank_similar(new_all_copies, polyA_flank_blast_dir, blast_program_dir)

    # 2.判断我们具有准确边界的LTR是否是真实的。
    # 条件：
    # ①.它要有多份拷贝（单拷贝的序列需要靠判断它是否出现在“连续性原件“的直接侧翼序列，如基因、CDS或另一个转座子，因此我们不考虑单拷贝）。
    # ②.我们在获得拷贝后，需要扩展50 bp范围。然后取所有以polyA结尾的拷贝。
    # ③.记录下有polyA结构的拷贝及数量（robust of the evidence）。
    # candidate_line_path = tmp_output_dir + '/candidate_line.fa'
    # flanking_len = 50
    # is_transposons(candidate_line_path, tmp_output_dir, flanking_len)




