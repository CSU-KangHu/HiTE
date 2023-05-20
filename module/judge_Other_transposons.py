#-- coding: UTF-8 --
import argparse
import codecs
import os
import sys

import json

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, getReverseSequence, read_fasta_v2, get_copies, flanking_copies, \
    multi_process_align, multi_process_align_and_get_copies, flanking_copies_v2, rename_fasta, Logger, file_exist, \
    get_seq_families


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
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--member_script_path', metavar='member_script_path',
                        help='e.g., /home/hukang/HiTE/tools/make_fasta_from_blast.sh')
    parser.add_argument('--subset_script_path', metavar='subset_script_path',
                        help='e.g., /home/hukang/HiTE/tools/ready_for_MSA.sh')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/dmel')
    parser.add_argument('--library_dir', metavar='library_dir',
                        help='e.g., ')
    parser.add_argument('--recover', metavar='recover',
                        help='e.g., 0')
    parser.add_argument('-r', metavar='Reference path',
                        help='input Reference path')


    args = parser.parse_args()
    threads = int(args.t)
    member_script_path = args.member_script_path
    subset_script_path = args.subset_script_path
    tmp_output_dir = args.tmp_output_dir
    library_dir = args.library_dir
    recover = args.recover
    reference = args.r

    reference = os.path.realpath(reference)

    is_recover = False
    recover = int(recover)
    if recover == 1:
        is_recover = True

    tmp_output_dir = os.path.abspath(tmp_output_dir)

    log = Logger(tmp_output_dir + '/HiTE_other.log', level='debug')

    # 一些思路：
    # 除了LTR、TIR、Helitron之外的其他转座子，包括LINE、SINE、DIRS、PLE(在Dfam中属于LINE)、Crypton它们缺少或具有复杂的终端结构特点，且没有稳定的TSD特征。
    # 想要根据结构特征去识别有困难，我们根据同源性搜索的方法去识别。
    # LINE （1000-7000bp），通常以poly(A)结尾和 SINE(100-600bp)，generate TSDs (5–15 bp)，通常以poly(T)结尾，发现也有polyA结尾。我们还需要考虑反向互补序列。

    confident_other_path = tmp_output_dir + '/confident_other.fa'
    resut_file = confident_other_path
    if not is_recover or not file_exist(resut_file):
        non_LTR_lib = library_dir + '/non_LTR.lib'
        other_TE_dir = tmp_output_dir + '/other_TE'
        os.system('rm -rf ' + other_TE_dir)
        if not os.path.exists(other_TE_dir):
            os.makedirs(other_TE_dir)


        confident_non_ltr_contigs = {}
        contignames, contigs = read_fasta(non_LTR_lib)

        # 1.将库比对到参考上，获取其拷贝
        flanking_len = 0
        copy_files = get_seq_families(non_LTR_lib, reference, member_script_path, subset_script_path, flanking_len,
                                      tmp_output_dir, threads)

        # 2.取最长拷贝当做发现的non-LTR元素
        for copy_file in copy_files:
            name = copy_file[0]
            member_names, member_contigs = read_fasta(copy_file[1])
            member_items = member_contigs.items()
            member_items = sorted(member_items, key=lambda x: len(x[1]), reverse=True)
            if len(member_items) > 0:
                seq = member_items[0][1]
                confident_non_ltr_contigs[name] = seq
        store_fasta(confident_non_ltr_contigs, confident_other_path)
        rename_fasta(confident_other_path, confident_other_path, 'Other')
    else:
        log.logger.info(resut_file + ' exists, skip...')
