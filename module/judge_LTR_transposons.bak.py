import argparse
import codecs
import os
import sys

import json

# cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# sys.path.append(cur_dir)
import time

from Util import read_fasta, store_fasta, Logger, get_copies, multi_process_align, store_LTR_seq, getReverseSequence, \
    flank_region_align_v1, multi_process_TRF


def is_canonical_ltr(seq, tsd_len):
    # 判断序列是否具备经典的5 bp TSD immediately connected with the 5'-TG…CA-3
    first_2bp = seq[0:2]
    end_2bp = seq[-2:]
    if first_2bp == 'TG' and end_2bp == 'CA' and tsd_len == 5:
        return True
    return False


def is_transposons(candidate_ltr_path, reference, threads, tmp_output_dir, blast_program_dir, ref_index):
    log.logger.info('determine true LTR')

    log.logger.info('------flank LTR copy and see if the flanking regions are repeated')
    starttime = time.time()
    # 我们将copies扩展50bp，一个orig_query_name对应一个文件，然后做自比对。
    # 解析每个自比对文件，判断C0与C1,C2...等拷贝的比对情况，如果有flanking区域包含在比对区域内，那么这条拷贝应该被抛弃，如果所有拷贝被抛弃，则该条序列应该是假阳性。
    flanking_len = 50
    similar_ratio = 0.1
    TE_type = 'ltr'
    confident_copies = flank_region_align_v1(candidate_ltr_path, flanking_len, similar_ratio, reference, TE_type,
                                             tmp_output_dir, blast_program_dir, threads, ref_index, log)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of flanking LTR copy and see if the flanking regions are repeated: %.8s s" % (dtime))

    confident_ltr_path = tmp_output_dir + '/confident_ltr_'+str(ref_index)+'.fa'
    confident_ltr_cut_path = tmp_output_dir + '/confident_ltr_cut_'+str(ref_index)+'.fa'
    confident_ltr = {}
    confident_ltr_cut = {}
    ltr_names = set()
    ltr_names.update(set(confident_copies.keys()))

    candidate_names, candidate_contigs = read_fasta(candidate_ltr_path)

    for name in ltr_names:
        seq = candidate_contigs[name]
        # 如果有连续10个以上的N就过滤掉
        if seq.__contains__('NNNNNNNNNN'):
            continue

        #copy_list = confident_copies[name]
        #if len(copy_list) >= 2 or is_canonical_ltr(seq, tsd_len):
        #if len(copy_list) >= 2:

        confident_ltr[name] = seq

        # cut the intact LTR into left, right and internal sequences
        # name format: N_0-LTR-lLTRStart_100-lLTREnd_400-rLTRStart_800-rLTREnd_1100-tsd_atcg
        parts = name.split('-')
        lLTR_start = int(parts[2].split('_')[1])
        lLTR_end = int(parts[3].split('_')[1])

        rLTR_start = int(parts[4].split('_')[1])
        rLTR_end = int(parts[5].split('_')[1])

        left_LTR = seq[lLTR_start - 1: lLTR_end]
        right_LTR = seq[rLTR_start - 1: rLTR_end]

        LTR_internal = seq[lLTR_end: rLTR_start - 1]

        left_LTR_name = parts[0] + '-lLTR'
        right_LTR_name = parts[0] + '-rLTR'
        internal_query_name = parts[0] + '-ILTR'

        if len(left_LTR) > len(right_LTR):
            confident_ltr_cut[left_LTR_name] = left_LTR
        else:
            confident_ltr_cut[right_LTR_name] = right_LTR
        # confident_ltr_cut[left_LTR_name] = left_LTR
        # confident_ltr_cut[right_LTR_name] = right_LTR
        confident_ltr_cut[internal_query_name] = LTR_internal

    log.logger.info('------store confident LTR sequences')
    store_fasta(confident_ltr, confident_ltr_path)
    store_fasta(confident_ltr_cut, confident_ltr_cut_path)

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa')
    parser.add_argument('--blast_program_dir', metavar='blast_program_dir',
                        help='e.g., /public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2')
    # parser.add_argument('--ltrharvest_output', metavar='ltrharvest_output',
    #                     help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa/genome_0.fa.harvest.scn')
    parser.add_argument('--TRF_Path', metavar='TRF_Path',
                        help='e.g., /public/home/hpc194701009/repeat_detect_tools/trf409.linux64')
    parser.add_argument('--tandem_region_cutoff', metavar='tandem_region_cutoff',
                        help='e.g., 0.5')
    parser.add_argument('--flanking_len', metavar='flanking_len',
                        help='e.g., 20')
    parser.add_argument('--ref_index', metavar='ref_index',
                        help='e.g., 0')
    parser.add_argument('--debug', metavar='debug',
                        help='e.g., 0')

    args = parser.parse_args()

    reference = args.g
    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    blast_program_dir = args.blast_program_dir
    # ltrharvest_output = args.ltrharvest_output
    TRF_Path = args.TRF_Path
    tandem_region_cutoff = float(args.tandem_region_cutoff)
    flanking_len = int(args.flanking_len)
    ref_index = args.ref_index
    debug = int(args.debug)

    log = Logger('HiTE.log', level='debug')

    confident_ltr_path = tmp_output_dir + '/confident_ltr_'+str(ref_index)+'.fa'
    confident_ltr_cut_path = tmp_output_dir + '/confident_ltr_cut_'+str(ref_index)+'.fa'
    if os.path.exists(confident_ltr_path):
        os.remove(confident_ltr_path)
    if os.path.exists(confident_ltr_cut_path):
        os.remove(confident_ltr_cut_path)

    log.logger.info('storing candidate LTR')

    # 1.因为现在的LTR序列不能准确的确定边界，因此下面的步骤都是为了能够准确获得边界
    # 这里是LTR_harvest结果
    candidate_ltr_harvest_path = tmp_output_dir + '/candidate_ltr_harvest_'+str(ref_index)+'.fa'
    candidate_ltr_harvest_cut_path = tmp_output_dir + '/candidate_ltr_harvest_cut_'+str(ref_index)+'.fa'

    # 2.因为现在的LTR序列不能准确的确定边界，因此下面的步骤都是为了能够准确获得边界
    # 这里是LTR_Finder结果
    candidate_ltr_finder_path = tmp_output_dir + '/candidate_ltr_finder_'+str(ref_index)+'.fa'
    candidate_ltr_finder_cut_path = tmp_output_dir + '/candidate_ltr_finder_cut_'+str(ref_index)+'.fa'

    # # 3. 这里是来自TIR中过滤出来的LTR
    # ltr_from_itr_path = tmp_output_dir + '/confident_ltr_from_itr_' + str(ref_index) + '.fa'
    # ltr_from_itr_cut_path = tmp_output_dir + '/confident_ltr_from_itr_cut_' + str(ref_index) + '.fa'

    candidate_ltr_path = tmp_output_dir + '/candidate_ltr_'+str(ref_index)+'.fa'
    candidate_ltr_cut_path = tmp_output_dir + '/candidate_ltr_cut_'+str(ref_index)+'.fa'
    os.system('cat ' + candidate_ltr_harvest_path + ' > ' + candidate_ltr_path)
    os.system('cat ' + candidate_ltr_finder_path + ' >> ' + candidate_ltr_path)
    #os.system('cat ' + ltr_from_itr_path + ' >> ' + candidate_ltr_path)
    os.system('cat ' + candidate_ltr_harvest_cut_path + ' > ' + candidate_ltr_cut_path)
    os.system('cat ' + candidate_ltr_finder_cut_path + ' >> ' + candidate_ltr_cut_path)
    #os.system('cat ' + ltr_from_itr_cut_path + ' >> ' + candidate_ltr_cut_path)

    contignames, contigs = read_fasta(candidate_ltr_path)

    #过滤掉串联重复
    (repeat_dir, repeat_filename) = os.path.split(candidate_ltr_path)
    (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)
    repeats_path = tmp_output_dir + '/candidate_ltr_cut_'+str(ref_index)+'.filter_tandem.fa'
    trf_dir = tmp_output_dir + '/ltr_trf_temp'
    #去掉那些在终端50 bp、LTR、Internal中存在50%以上串联重复的序列
    multi_process_TRF(candidate_ltr_cut_path, repeats_path, TRF_Path, trf_dir, tandem_region_cutoff, threads=threads, TE_type='ltr')
    filter_tandem_contignames, filter_tandem_contigs = read_fasta(repeats_path)

    if debug == 0:
        #remove temp dir
        os.system('rm -rf ' + trf_dir)

    candidate_ltr_contigs = {}
    for name in contignames:
        header = name.split('-')[0]
        internal_name = header + '-ILTR'
        # 过滤之后的文件，既包含LTR又包含Internal序列，说明可以保留
        if filter_tandem_contigs.__contains__(name) and filter_tandem_contigs.__contains__(internal_name):
            candidate_ltr_contigs[name] = contigs[name]

    candidate_ltr_path = tmp_output_dir + '/candidate_ltr_'+str(ref_index)+'.filter_tandem.fa'
    store_fasta(candidate_ltr_contigs, candidate_ltr_path)


    # 2.判断我们具有准确边界的LTR是否是真实的。
    # 条件：
    # ①.它要有多份拷贝（单拷贝的序列需要靠判断它是否出现在“连续性原件“的直接侧翼序列，如基因、CDS或另一个转座子，因此我们不考虑单拷贝）。
    # ②.判断它的拷贝是否有相同长度的TSD。在通过比对获得拷贝边界时，经常由于不是整个序列的全比对，导致拷贝的准确边界无法识别。
    # 因此，我们在获得拷贝后，需要扩展50 bp范围。然后可以用LTRharvest继续识别边界。
    # ③.记录下有TSD+LTR结构的拷贝及数量（robust of the evidence）。
    is_transposons(candidate_ltr_path, reference, threads, tmp_output_dir, blast_program_dir, ref_index)




