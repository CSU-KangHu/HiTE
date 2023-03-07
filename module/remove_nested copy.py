#-- coding: UTF-8 --
import argparse
import codecs
import multiprocessing
import os
import sys

import json
import time

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, multi_process_align_and_get_copies, remove_ltr_from_tir, rename_fasta, \
    flank_region_align_v1


def is_transposons(filter_dup_path, reference, threads, tmp_output_dir, ref_index, log):
    log.logger.info('determine true TIR')

    log.logger.info('------flank TIR copy and see if the flanking regions are repeated')
    starttime = time.time()
    # 我们将copies扩展50bp，一个orig_query_name对应一个文件，然后做自比对。
    # 解析每个自比对文件，判断C0与C1,C2...等拷贝的比对情况，如果有flanking区域包含在比对区域内，那么这条拷贝应该被抛弃，如果所有拷贝被抛弃，则该条序列应该是假阳性。
    flanking_len = 50
    similar_ratio = 0.1
    TE_type = 'tir'
    confident_copies = flank_region_align_v1(filter_dup_path, flanking_len, similar_ratio, reference, TE_type, tmp_output_dir, threads, ref_index, log)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of flanking TIR copy and see if the flanking regions are repeated: %.8s s" % (dtime))

    log.logger.info('------store confident TIR sequences')
    filter_dup_names, filter_dup_contigs = read_fasta(filter_dup_path)
    if ref_index == -1:
        confident_tir_path = tmp_output_dir + '/confident_tir.rename.cons.fa'
    else:
        confident_tir_path = tmp_output_dir + '/confident_tir_'+str(ref_index)+'.fa'
    confident_tir = {}
    for name in confident_copies.keys():
        copy_list = confident_copies[name]
        if len(copy_list) >= 2:
            confident_tir[name] = filter_dup_contigs[name]
    store_fasta(confident_tir, confident_tir_path)

if __name__ == '__main__':
    #preprocess()
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('--confident_ltr_cut', metavar='confident_ltr_cut',
                        help='e.g., ')
    parser.add_argument('--confident_tir', metavar='confident_tir',
                        help='e.g., ')
    parser.add_argument('--confident_helitron', metavar='confident_helitron',
                        help='e.g., ')
    parser.add_argument('--confident_other', metavar='confident_other',
                        help='e.g., ')
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/dmel')
    parser.add_argument('--global_flanking_filter', metavar='global_flanking_filter',
                        help='e.g., 1')
    parser.add_argument('--remove_nested', metavar='remove_nested',
                        help='e.g., 1')
    parser.add_argument('--test_home', metavar='test_home',
                        help='e.g., ')




    args = parser.parse_args()
    reference = args.g
    confident_ltr_cut_path = args.confident_ltr_cut
    confident_tir_path = args.confident_tir
    confident_helitron_path = args.confident_helitron
    confident_other_path = args.confident_other
    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    global_flanking_filter = int(args.global_flanking_filter)
    remove_nested = int(args.remove_nested)
    test_home = args.test_home

    log = Logger(tmp_output_dir + '/HiTE.log', level='debug')

    # 1.2 confident_ltr_cut_path比对到TIR候选序列上，并且过滤掉出现在LTR库中的TIR序列
    temp_dir = tmp_output_dir + '/tir_blast_ltr'
    all_copies = multi_process_align_and_get_copies(confident_ltr_cut_path, confident_tir_path, temp_dir, 'tir',
                                                    threads, query_coverage=0.8)
    remove_ltr_from_tir(confident_ltr_cut_path, confident_tir_path, all_copies)

    # 1.4 生成一致性tir序列
    confident_tir_rename_path = tmp_output_dir + '/confident_tir.rename.fa'
    rename_fasta(confident_tir_path, confident_tir_rename_path)

    confident_tir_rename_consensus = tmp_output_dir + '/confident_tir.rename.cons.fa'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_tir_rename_path + ' -o ' + confident_tir_rename_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    # # 如果切分成了多个块，TIR需要重新flank_region_align_v1到整个基因组，以过滤掉那些在分块中未能过滤掉的假阳性。
    # if global_flanking_filter == 1:
    #     ref_index = -1
    #     is_transposons(confident_tir_rename_consensus, reference, threads, tmp_output_dir, ref_index, log)

    # 1.5 解开TIR中包含的nested TE
    clean_tir_path = tmp_output_dir + '/confident_tir.clean.fa'
    remove_nested_command = 'cd ' + test_home + ' && python3 ' + test_home + '/remove_nested_lib.py ' \
                            + ' -t ' + str(threads) \
                            + ' --tmp_output_dir ' + tmp_output_dir + ' --max_iter_num ' + str(5) \
                            + ' --input1 ' + confident_tir_rename_consensus \
                            + ' --input2 ' + confident_tir_rename_consensus \
                            + ' --output ' + clean_tir_path
    os.system(remove_nested_command)

    # # cd-hit -aS 0.95 -c 0.8合并一些冗余序列
    # clean_tir_consensus = tmp_output_dir + '/confident_tir.clean.cons.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + clean_tir_path + ' -o ' + clean_tir_consensus + ' -T 0 -M 0'
    # os.system(cd_hit_command)

    # # 1.6 生成一致性other序列
    # confident_other_rename_path = tmp_output_dir + '/confident_other.rename.fa'
    # rename_fasta(confident_other_path, confident_other_rename_path)
    #
    # confident_other_rename_consensus = tmp_output_dir + '/confident_other.rename.cons.fa'
    # cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + confident_other_rename_path + ' -o ' + confident_other_rename_consensus + ' -T 0 -M 0'
    # os.system(cd_hit_command)

    # 合并所有的TE（TIR+Helitron+Other）
    confident_TE_path = tmp_output_dir + '/confident_TE.fa'
    os.system('cat ' + clean_tir_path + ' > ' + confident_TE_path)
    os.system('cat ' + confident_helitron_path + ' >> ' + confident_TE_path)
    os.system('cat ' + confident_other_path + ' >> ' + confident_TE_path)

    # 解开LTR内部包含的nested TE，然后把解开后的LTR合并到TE库中
    # 获取LTR的内部序列
    confident_ltr_terminal_path = tmp_output_dir + '/confident_ltr_cut.terminal.fa'
    confident_ltr_internal_path = tmp_output_dir + '/confident_ltr_cut.internal.fa'
    ltr_names, ltr_contigs = read_fasta(confident_ltr_cut_path)
    ltr_internal_contigs = {}
    ltr_terminal_contigs = {}
    for name in ltr_names:
        if name.__contains__('_INT#'):
            ltr_internal_contigs[name] = ltr_contigs[name]
        else:
            ltr_terminal_contigs[name] = ltr_contigs[name]
    store_fasta(ltr_internal_contigs, confident_ltr_internal_path)
    store_fasta(ltr_terminal_contigs, confident_ltr_terminal_path)

    clean_ltr_internal_path = tmp_output_dir + '/confident_ltr_cut.internal.clean.fa'
    if remove_nested == 1:
        starttime = time.time()
        log.logger.info('Start step2.4: remove nested TE in LTR internal')
        # 将所有的LTR序列暂时加入到TE序列中，用来解开nested TE
        temp_confident_TE_path = tmp_output_dir + '/confident_TE.temp.fa'
        os.system('cat ' + confident_ltr_cut_path + ' > ' + temp_confident_TE_path)
        os.system('cat ' + confident_TE_path + ' >> ' + temp_confident_TE_path)
        if os.path.getsize(temp_confident_TE_path) > 0 and os.path.getsize(confident_ltr_internal_path) > 0:
            remove_nested_command = 'cd ' + test_home + ' && python3 ' + test_home + '/remove_nested_lib.py ' \
                                    + ' -t ' + str(threads) \
                                    + ' --tmp_output_dir ' + tmp_output_dir + ' --max_iter_num ' + str(5) \
                                    + ' --input1 ' + temp_confident_TE_path \
                                    + ' --input2 ' + confident_ltr_internal_path \
                                    + ' --output ' + clean_ltr_internal_path
            os.system(remove_nested_command)

            os.system('cat ' + confident_ltr_terminal_path + ' > ' + confident_ltr_cut_path)
            os.system('cat ' + clean_ltr_internal_path + ' >> ' + confident_ltr_cut_path)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step2.4: %.8s s" % (dtime))
    os.system('cat ' + confident_ltr_cut_path + ' >> ' + confident_TE_path)

    # confident_ltr_cut_consensus = tmp_output_dir + '/confident_ltr_cut.cons.fa'
    # cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + confident_ltr_cut_path + ' -o ' + confident_ltr_cut_consensus + ' -T 0 -M 0'
    # os.system(cd_hit_command)




