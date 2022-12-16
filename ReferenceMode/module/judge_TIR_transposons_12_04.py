import argparse
import codecs
import os
import sys

import json
import time

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, read_fasta_v1, store_fasta, getReverseSequence, \
    Logger, calculate_max_min, get_copies, flanking_copies, \
    multi_process_tsd, multi_process_itr, filter_dup_itr, multi_process_align, flank_region_align_v1, store_copies_v1, \
    multi_process_tsd_v1, multi_process_TRF, multi_process_ltr


#from Util_11_14 import multi_process_tsd


def is_transposons(filter_dup_path, reference, threads, tmp_output_dir, flanking_len, blast_program_dir, TRsearch_dir):
    log.logger.info('determine true TIR')

    log.logger.info('------flank TIR copy and see if the flanking regions are repeated')
    starttime = time.time()
    # 我们将copies扩展50bp，一个orig_query_name对应一个文件，然后做自比对。
    # 解析每个自比对文件，判断C0与C1,C2...等拷贝的比对情况，如果有flanking区域包含在比对区域内，那么这条拷贝应该被抛弃，如果所有拷贝被抛弃，则该条序列应该是假阳性。
    flanking_len = 50
    similar_ratio = 0.2
    TE_type = 'tir'
    confident_copies = flank_region_align_v1(filter_dup_path, flanking_len, similar_ratio, reference, TE_type, tmp_output_dir, blast_program_dir, threads, log)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of flanking TIR copy and see if the flanking regions are repeated: %.8s s" % (dtime))

    log.logger.info('------store confident TIR sequences')
    filter_dup_names, filter_dup_contigs = read_fasta(filter_dup_path)
    confident_tir_path = tmp_output_dir + '/confident_tir_'+str(ref_index)+'.fa'
    confident_tir = {}
    for name in confident_copies.keys():
        copy_list = confident_copies[name]
        if len(copy_list) >= 2:
            confident_tir[name] = filter_dup_contigs[name]
    store_fasta(confident_tir, confident_tir_path)


def get_score(confident_TIR):
    # (tsd, te_len, cur_tir_info, contigs[name])
    TE_len_list = []
    tsd_len_list = []
    for info in confident_TIR:
        TE_len_list.append(info[1])
        tsd_len_list.append(len(info[0]))

    max_TE_len = max(TE_len_list)
    min_TE_len = min(TE_len_list)

    max_tsd_len = max(tsd_len_list)
    min_tsd_len = min(tsd_len_list)

    res_confident_TIR = []
    for info in confident_TIR:
        if max_TE_len == min_TE_len:
            TE_len_normal = 0
        else:
            TE_len_normal = calculate_max_min(info[1], max_TE_len, min_TE_len)
        if max_tsd_len == min_tsd_len:
            tsd_len_normal = 0
        else:
            tsd_len_normal = calculate_max_min(len(info[0]), max_tsd_len, min_tsd_len)
        score = 0.4*TE_len_normal+0.6*tsd_len_normal
        res_confident_TIR.append((score, info[2], info[3]))
    res_confident_TIR.sort(key=lambda x: -x[0])
    return res_confident_TIR

def quick_tir_search(tir_tsd_path, tir_out, filter_dup_path, quick_tir_dir):
    # ①首先我们将候选序列按照打分排序，形成 candidate_tir_dict = {query: [(c1, score1, seq), (c2, score2, seq)]}。
    # ②遍历candidate_tir_dict，对于每个query，选择打分最高的序列，形成round1.fa，并用一个 candidate_tir_index = {query: index}记录已选择query对应数组中的第index个数据，判断round1.fa中数量大小，如果小于10000，再选择新的序列加进来。
    # ③调用multi_process_itr。读取round1.fa.itr,对于每个query，删除 candidate_tir_dict中对应的记录，表示这个query已找到具有TSD+TIR结构的序列。
    # ④如果len(candidate_tir_dict)>0或者iter_num < max_iter_num，进入②
    contignames, contigs = read_fasta_v1(tir_tsd_path)
    candidate_tir_dict = {}
    for name in contignames:
        parts = name.split('-C_')
        orig_query_name = parts[0]
        tsd = parts[1].split(' ')[0].split('-tsd_')[1]
        te_len = len(contigs[name])
        cur_tir_info = parts[1]
        # te_len = int(parts[1].split('-')[1].split('_')[1])
        if not candidate_tir_dict.__contains__(orig_query_name):
            candidate_tir_dict[orig_query_name] = set()
        confident_TIR = candidate_tir_dict[orig_query_name]
        confident_TIR.add((tsd, te_len, cur_tir_info, contigs[name]))

    for name in candidate_tir_dict.keys():
        confident_TIR = candidate_tir_dict[name]
        confident_TIR = get_score(confident_TIR)
        candidate_tir_dict[name] = confident_TIR

    candidate_tir_index = {}
    max_iter_num = 10
    iter_num = 0

    while (iter_num < max_iter_num and len(candidate_tir_dict) > 0):
        print('quick tir search round:' + str(iter_num) + ' ...')
        quick_tir_path = quick_tir_dir + '/round_'+str(iter_num)+'.fa'
        quick_tirs = {}
        while len(quick_tirs) < 10000:
            for name in list(candidate_tir_dict):
                confident_TIR = candidate_tir_dict[name]
                if not candidate_tir_index.__contains__(name):
                    candidate_tir_index[name] = 0
                cur_index = candidate_tir_index[name]
                if cur_index >= len(confident_TIR):
                    del candidate_tir_dict[name]
                    continue
                cur_info = confident_TIR[cur_index]
                quick_tirs[name+'-C_'+cur_info[1]] = cur_info[2]
                candidate_tir_index[name] = cur_index+1
        store_fasta(quick_tirs, quick_tir_path)
        quick_tir_out = quick_tir_path + '.itr'
        multi_process_itr(quick_tir_path, quick_tir_out, tmp_output_dir, TRsearch_dir)

        quick_out_names, quick_out_contigs = read_fasta(quick_tir_out)
        for quick_name in quick_out_names:
            orig_name = quick_name.split('-C_')[0]
            del candidate_tir_dict[orig_name]
        iter_num += 1
        print('this round found sequence number:' + str(len(quick_out_names)) + ', remaining:'+str(len(candidate_tir_dict)))

    if os.path.exists(tir_out):
        os.remove(tir_out)

    for i in range(max_iter_num):
        quick_tir_path = quick_tir_dir + '/round_'+str(i)+'.fa'
        quick_tir_out = quick_tir_path + '.itr'
        if os.path.exists(quick_tir_out):
            os.system('cat '+quick_tir_out+' >> '+tir_out)

    filter_dup_path = tmp_output_dir + '/tir_tsd.filter_dup.fa'
    filter_dup_itr(tir_out, filter_dup_path)



if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('--copies', metavar='copies',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa/longest_repeats_0.copies')
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--TRsearch_dir', metavar='TRsearch_dir',
                        help='e.g., /public/home/hpc194701009/repeat_detect_tools/REPET_linux-x64-3.0/bin')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa')
    parser.add_argument('--blast_program_dir', metavar='blast_program_dir',
                        help='e.g., /public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2')
    parser.add_argument('--TRF_Path', metavar='TRF_Path',
                        help='e.g., /public/home/hpc194701009/repeat_detect_tools/trf409.linux64')
    parser.add_argument('--tandem_region_cutoff', metavar='tandem_region_cutoff',
                        help='e.g., 0.5')
    parser.add_argument('--flanking_len', metavar='flanking_len',
                        help='e.g., 20')
    parser.add_argument('--plant', metavar='plant',
                        help='e.g., 1')
    parser.add_argument('--ref_index', metavar='ref_index',
                        help='e.g., 0')

    args = parser.parse_args()
    reference = args.g
    longest_repeats_flanked_copies_file = args.copies
    threads = int(args.t)
    TRsearch_dir = args.TRsearch_dir
    tmp_output_dir = args.tmp_output_dir
    blast_program_dir = args.blast_program_dir
    flanking_len = int(args.flanking_len)
    plant = int(args.plant)
    TRF_Path = args.TRF_Path
    tandem_region_cutoff = float(args.tandem_region_cutoff)
    ref_index = args.ref_index

    log = Logger('HiTE.log', level='debug')

    confident_tir_path = tmp_output_dir + '/confident_tir_'+str(ref_index)+'.fa'
    if os.path.exists(confident_tir_path):
        os.remove(confident_tir_path)

    log.logger.info('loading ' + longest_repeats_flanked_copies_file)
    file = open(longest_repeats_flanked_copies_file, 'r')
    js = file.read()
    all_copies = json.loads(js)

    # 取20条全长拷贝两端flanking 50bp以包含TSD，因为我们需要靠拷贝中是否有相同长度的TSD数量来支持，所以拷贝数量不能太少

    # 对于每条序列而言，获取它具有TIR+TSD结构的拷贝，分析拷贝中出现最多次数的tsd_len，取第一条具有tsd_len的拷贝
    log.logger.info('------get TIR+TSD in copies of candidate TIR')
    starttime = time.time()
    tir_tsd_path = tmp_output_dir + '/tir_tsd_'+str(ref_index)+'.fa'
    tir_tsd_dir = tmp_output_dir + '/tir_tsd_temp'
    multi_process_tsd(all_copies, tir_tsd_path, tir_tsd_dir, flanking_len, threads, TRsearch_dir, plant)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of getting TSD in copies of candidate TIR: %.8s s" % (dtime))

    # 过滤掉串联重复
    (repeat_dir, repeat_filename) = os.path.split(tir_tsd_path)
    (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)
    repeats_path = tmp_output_dir + '/tir_tsd_'+str(ref_index)+'.filter_tandem.fa'
    trf_dir = tmp_output_dir + '/tir_trf_temp'
    # 去掉那些在终端20 bp、LTR、Internal中存在50%以上串联重复的序列
    multi_process_TRF(tir_tsd_path, repeats_path, TRF_Path, trf_dir, tandem_region_cutoff, threads=threads,
                      TE_type='tir')

    # # 2.发现我们识别的候选TIR元素里面有一部分是LTR，因此我们根据ltrsearch和tsd在4-6bp来过滤出这些LTR元素，然后把它补充到LTR模块中去。
    # ltr_from_itr_path = tmp_output_dir + '/confident_ltr_from_itr_'+str(ref_index)+'.fa'
    # ltr_from_itr_cut_path = tmp_output_dir + '/confident_ltr_from_itr_cut_' + str(ref_index) + '.fa'
    # temp_dir = tmp_output_dir + '/tir_ltr_temp'
    # multi_process_ltr(repeats_path, ltr_from_itr_path, temp_dir, TRsearch_dir, threads=48)
    #
    # # 2.1 在tir_tsd.filter_tandem.fa中删除ltr元素
    # ltr_contignames, ltr_contigs = read_fasta_v1(ltr_from_itr_path)
    # tir_contignames, tir_contigs = read_fasta(repeats_path)
    # for ltr_contigname in ltr_contignames:
    #     parts = ltr_contigname.split(' ')
    #     orig_name = parts[0]
    #     del tir_contigs[orig_name]
    # store_fasta(tir_contigs, repeats_path)
    #
    # # 2.2 只保留ltr中具有4-6 bp TSD的元素
    # node_index = 0
    # confident_ltr_contigs = {}
    # confident_ltr_cut_contigs = {}
    # for ltr_contigname in ltr_contignames:
    #     tsd = ltr_contigname.split('-tsd_')[1].split('-')[0]
    #     tsd_len = len(tsd)
    #     if tsd_len < 4 or tsd_len > 6:
    #         continue
    #     else:
    #         seq = ltr_contigs[ltr_contigname]
    #         parts = ltr_contigname.split(' ')
    #         query_name = 'N_' + str(node_index)
    #         # LTR(1,594)..(2154,2747)
    #         LTR_info = parts[1].replace('LTR', '').replace('(', '').replace(')', '')
    #         LTR_info_parts = LTR_info.split('..')
    #         LTR_left_pos_parts = LTR_info_parts[0].split(',')
    #         LTR_right_pos_parts = LTR_info_parts[1].split(',')
    #         lLTR_start = int(LTR_left_pos_parts[0])
    #         lLTR_end = int(LTR_left_pos_parts[1])
    #         rLTR_start = int(LTR_right_pos_parts[0])
    #         rLTR_end = int(LTR_right_pos_parts[1])
    #
    #         left_LTR = seq[lLTR_start - 1: lLTR_end]
    #         right_LTR = seq[rLTR_start - 1: rLTR_end]
    #
    #         LTR_internal = seq[lLTR_end: rLTR_start - 1]
    #
    #         left_LTR_name = query_name + '-lLTR'
    #         right_LTR_name = query_name + '-rLTR'
    #         internal_query_name = query_name + '-ILTR'
    #
    #         if len(left_LTR) > len(right_LTR):
    #             confident_ltr_cut_contigs[left_LTR_name] = left_LTR
    #         else:
    #             confident_ltr_cut_contigs[right_LTR_name] = right_LTR
    #         confident_ltr_cut_contigs[internal_query_name] = LTR_internal
    #
    #         confident_ltr_contigs[query_name] = ltr_contigs[ltr_contigname]
    #         node_index += 1
    # store_fasta(confident_ltr_contigs, ltr_from_itr_path)
    # store_fasta(confident_ltr_cut_contigs, ltr_from_itr_cut_path)


    # 3.判断我们具有准确边界的TIR是否是真实的。
    # 条件：
    # ①.它要有多份拷贝（单拷贝的序列需要靠判断它是否出现在“连续性原件“的直接侧翼序列，如基因、CDS或另一个转座子，因此我们不考虑单拷贝）。
    # ②.判断它的拷贝是否有相同长度的TSD。在通过比对获得拷贝边界时，经常由于不是整个序列的全比对，导致拷贝的准确边界无法识别。
    # 因此，我们在获得拷贝后，需要扩展50 bp范围，记录此时的边界s1, e1，并且在[0:s1, e1:]范围内搜索相同长度的TSD。
    # ③.判断以TSD为边界的TIR拷贝是否具有itr结构，记录下有TSD+TIR结构的拷贝及数量（robust of the evidence）。
    repeats_path = tmp_output_dir + '/tir_tsd_'+str(ref_index)+'.filter_tandem.fa'
    is_transposons(repeats_path, reference, threads, tmp_output_dir, flanking_len, blast_program_dir, TRsearch_dir)
