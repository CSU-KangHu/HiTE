import argparse
import codecs
import os
import sys

import json
import time

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, multi_process_helitronscanner, get_copies, multi_process_align, \
    flank_region_align_v1, multi_process_EAHelitron, Logger, flanking_copies

if __name__ == '__main__':
    log = Logger('HiTE.log', level='debug')
    EAHelitron = '/public/home/hpc194701009/repeat_detect_tools/EAHelitron-master'
    blast_program_dir = '/public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2'
    reference = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/cb4/GCF_000004555.2_CB4_genomic.fna'
    threads = 48
    flanking_len= 50

    # 1.提取repbase中Helitron元素，形成一个repbase.helitron.fa文件
    repbase_file = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/repbase/cbrrep.ref'
    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/cb4/test'
    repbase_helitron = tmp_output_dir + '/repbase.helitron.fa'
    contignames, contigs = read_fasta(repbase_file)
    helitron_contigs = {}
    for name in contignames:
        if name.__contains__('HELITRON'):
            helitron_contigs[name] = contigs[name]
    store_fasta(helitron_contigs, repbase_helitron)

    # 2.获取repbase.helitron.fa的拷贝，存成文件。然后用EAHelitron看能识别多少
    blastnResults_path = tmp_output_dir + '/repbase_helitron.ref.out'
    repbase_helitron_blast_dir = tmp_output_dir + '/repbase_helitron_blast'
    multi_process_align(repbase_helitron, reference, blastnResults_path, blast_program_dir, repbase_helitron_blast_dir,
                        threads)
    all_copies = get_copies(blastnResults_path, repbase_helitron, reference, query_coverage=0.99, threads=threads)
    all_copies = flanking_copies(all_copies, repbase_helitron, reference, flanking_len, copy_num=20)

    repbase_helitron_copies = tmp_output_dir + '/repbase_helitron.copies.fa'
    with open(repbase_helitron_copies, 'w') as f_save:
        for query_name in all_copies.keys():
            f_save.write('>' + query_name + '\n')
            copy_list = all_copies[query_name]
            # (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq)
            for index, copy in enumerate(copy_list):
                f_save.write('\t' + copy[0] + ':' + str(copy[1]) + '_' + str(copy[2]) + '_' + str(copy[3]) + '\n' + copy[4] + '\n')
    f_save.close()

    #多线程运行EAHelitro
    candidate_helitron_path = tmp_output_dir + '/repbase_EAHelitron.fa'
    temp_dir = tmp_output_dir + '/repbase_helitron_tmp'
    multi_process_EAHelitron(all_copies, flanking_len, candidate_helitron_path, temp_dir, EAHelitron, threads)
    candidate_helitron_contignames, candidate_helitron_contigs = read_fasta(candidate_helitron_path)
    candidate_helitron_rename = tmp_output_dir + '/repbase_EAHelitron.rename.fa'
    node_index = 0
    with open(candidate_helitron_rename, 'w') as f_save:
        for name in candidate_helitron_contignames:
            seq = candidate_helitron_contigs[name]
            new_name = 'Helitron_' + str(node_index)
            f_save.write('>' + new_name + '\n' + seq + '\n')
            node_index += 1
    f_save.close()
    candidate_helitron_contignames, candidate_helitron_contigs = read_fasta(candidate_helitron_rename)

    # 2.flanking candidate helitron 50bp,看flanking region是否有高同源性
    confident_copies = {}
    if len(candidate_helitron_contigs) > 0:
        log.logger.info('------flank Helitron copy and see if the flanking regions are repeated')
        starttime = time.time()
        # 我们将copies扩展50bp，一个orig_query_name对应一个文件，然后做自比对。
        # 解析每个自比对文件，判断C0与C1,C2...等拷贝的比对情况，如果有flanking区域包含在比对区域内，那么这条拷贝应该被抛弃，如果所有拷贝被抛弃，则该条序列应该是假阳性。
        flanking_len = 50
        similar_ratio = 0.2
        TE_type = 'helitron'
        confident_copies = flank_region_align_v1(candidate_helitron_rename, flanking_len, similar_ratio, reference, TE_type,
                                                 tmp_output_dir, blast_program_dir, threads, log)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info(
            "Running time of flanking TIR copy and see if the flanking regions are repeated: %.8s s" % (dtime))

    confident_helitron_path = tmp_output_dir + '/repbase_confident_helitron.fa'
    confident_helitron = {}
    for name in confident_copies.keys():
        copy_list = confident_copies[name]
        if len(copy_list) >= 2:
            confident_helitron[name] = candidate_helitron_contigs[name]
    store_fasta(confident_helitron, confident_helitron_path)



