#-- coding: UTF-8 --
import argparse
import codecs
import os
import sys

import json

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, convertToUpperCase_v1, multi_line

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/dmel')
    parser.add_argument('--chrom_seg_length', metavar='chrom_seg_length',
                        help='The length of genome segments')
    parser.add_argument('--chunk_size', metavar='chunk_size',
                        help='The chunk size of large genome')


    args = parser.parse_args()
    reference = args.g
    tmp_output_dir = args.tmp_output_dir
    chrom_seg_length = int(args.chrom_seg_length)
    chunk_size = float(args.chunk_size)

    tmp_output_dir = os.path.abspath(tmp_output_dir) 

    log = Logger(tmp_output_dir + '/HiTE.log', level='debug')

    log.logger.info('Start Splitting Reference into chunks')
    # using multiple threads to gain speed
    reference_pre = convertToUpperCase_v1(reference)

    # 将基因组切成更小的块，以提升后续的比对性能
    reference_tmp = multi_line(reference_pre, chrom_seg_length)
    cut_references = []
    cur_ref_contigs = {}
    cur_base_num = 0
    ref_index = 0
    with open(reference_tmp, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ref_name = parts[0].replace('>', '')
            start = parts[1]
            seq = parts[2]
            new_ref_name = ref_name + '$' + start
            cur_ref_contigs[new_ref_name] = seq
            cur_base_num += len(line)
            if cur_base_num >= chunk_size * 1024 * 1024:
                # store references
                cur_ref_path = tmp_output_dir + '/genome.cut' + str(ref_index) + '.fa'
                store_fasta(cur_ref_contigs, cur_ref_path)
                cut_references.append(cur_ref_path)
                cur_ref_contigs = {}
                cur_base_num = 0
                ref_index += 1
        if len(cur_ref_contigs) > 0:
            cur_ref_path = cur_ref_path = tmp_output_dir + '/genome.cut' + str(ref_index) + '.fa'
            store_fasta(cur_ref_contigs, cur_ref_path)
            cut_references.append(cur_ref_path)
    f_r.close()
    print(cut_references)






