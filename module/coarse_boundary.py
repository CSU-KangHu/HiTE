#-- coding: UTF-8 --
import argparse
import codecs
import os
import sys

import json
import time

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, determine_repeat_boundary_v3, multi_process_TRF, flanking_seq

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


    log = Logger(tmp_output_dir + '/HiTE.log', level='debug')


    repeats_path = cut_reference
    longest_repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.fa'
    # -------------------------------Stage02: this stage is used to do pairwise comparision, determine the repeat boundary-------------------------------
    determine_repeat_boundary_v3(repeats_path, longest_repeats_path, fixed_extend_base_threshold, max_repeat_len,
                                 tmp_output_dir, thread, ref_index, log)


    trf_dir = tmp_output_dir + '/trf_temp'
    (repeat_dir, repeat_filename) = os.path.split(longest_repeats_path)
    (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)
    repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.filter_tandem.fa'
    multi_process_TRF(longest_repeats_path, repeats_path, trf_dir, tandem_region_cutoff,
                      threads=thread)

    longest_repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.filter_tandem.fa'
    longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'
    flanking_seq(longest_repeats_path, longest_repeats_flanked_path, reference, flanking_len)

    # output
    # longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'






