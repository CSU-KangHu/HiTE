#-- coding: UTF-8 --
import argparse
import codecs
import os
import sys

import json
import time

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, determine_repeat_boundary_v3, multi_process_TRF, flanking_seq, \
    file_exist, determine_repeat_boundary_v4

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('--fixed_extend_base_threshold', metavar='fixed_extend_base_threshold',
                        help='The length of variation can be tolerated during pairwise alignment')
    parser.add_argument('--max_repeat_len', metavar='max_repeat_len', help='The maximum length of a single repeat')
    parser.add_argument('--thread', metavar='thread_num',
                        help='Input thread num')
    parser.add_argument('--flanking_len', metavar='flanking_len',
                        help='flanking_len')
    parser.add_argument('--tandem_region_cutoff', metavar='tandem_region_cutoff',
                        help='tandem_region_cutoff')
    parser.add_argument('--ref_index', metavar='ref_index',
                        help='ref_index')
    parser.add_argument('-r', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa')
    parser.add_argument('--recover', metavar='recover',
                        help='e.g., 0')
    parser.add_argument('--debug', metavar='recover',
                        help='e.g., 0')


    args = parser.parse_args()
    cut_reference = args.g
    fixed_extend_base_threshold = int(args.fixed_extend_base_threshold)
    max_repeat_len = int(args.max_repeat_len)
    thread = int(args.thread)
    ref_index = args.ref_index
    flanking_len = int(args.flanking_len)
    tandem_region_cutoff = float(args.tandem_region_cutoff)
    reference = args.r
    tmp_output_dir = args.tmp_output_dir
    recover = args.recover
    debug = args.debug

    if debug is None:
        debug = 0
    else:
        debug = int(debug)

    is_recover = False
    recover = int(recover)
    if recover == 1:
        is_recover = True
    
    tmp_output_dir = os.path.abspath(tmp_output_dir) 

    log = Logger(tmp_output_dir + '/HiTE_coarse.log', level='debug')


    repeats_path = cut_reference
    longest_repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.fa'
    resut_file = longest_repeats_path
    if not is_recover or not file_exist(resut_file):
        # -------------------------------Stage02: this stage is used to do pairwise comparision, determine the repeat boundary-------------------------------
        determine_repeat_boundary_v4(repeats_path, longest_repeats_path, fixed_extend_base_threshold, max_repeat_len,
                                     tmp_output_dir, thread, ref_index, debug)
    else:
        log.logger.info(resut_file + ' exists, skip...')

    # longest_repeats_cons = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.cons.fa'
    # resut_file = longest_repeats_cons
    # if not is_recover or not file_exist(resut_file):
    #     # 减少一些冗余序列
    #     cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                      + ' -G 0 -g 1 -A 80 -i ' + longest_repeats_path + ' -o ' + longest_repeats_cons + ' -T 0 -M 0'
    #     os.system(cd_hit_command)
    # else:
    #     log.logger.info(resut_file + ' exists, skip...')

    # repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.filter_tandem.fa'
    # resut_file = repeats_path
    # if not is_recover or not file_exist(resut_file):
    #     trf_dir = tmp_output_dir + '/trf_temp_' + str(ref_index)
    #     multi_process_TRF(longest_repeats_path, repeats_path, trf_dir, tandem_region_cutoff,
    #                       threads=thread)
    #     os.system('rm -rf ' + trf_dir)
    # else:
    #     log.logger.info(resut_file + ' exists, skip...')

    longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'
    resut_file = longest_repeats_flanked_path
    if not is_recover or not file_exist(resut_file):
        flanking_seq(longest_repeats_path, longest_repeats_flanked_path, reference, flanking_len)
    else:
        log.logger.info(resut_file + ' exists, skip...')


    # output
    # longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'






