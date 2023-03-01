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


if __name__ == '__main__':
    #preprocess()
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-g', metavar='Genome assembly cut',
                        help='input cut genome assembly path')
    parser.add_argument('--ref_name', metavar='ref_name',
                        help='e.g., ')
    parser.add_argument('--debug', metavar='debug',
                        help='e.g., 0')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/dmel')



    args = parser.parse_args()
    cut_references = args.g
    ref_name = args.ref_name
    debug = int(args.debug)
    tmp_output_dir = args.tmp_output_dir





