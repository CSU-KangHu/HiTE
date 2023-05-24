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
from Util import read_fasta, store_fasta, rename_reference, file_exist, Logger, run_LTR_harvest, run_LTR_retriever, \
    rename_fasta

if __name__ == '__main__':
    #preprocess()
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('--ltrfinder_home', metavar='ltrfinder_home',
                        help='e.g., ')
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/dmel')
    parser.add_argument('--recover', metavar='recover',
                        help='e.g., 0')
    parser.add_argument('--miu', metavar='miu',
                        help='e.g., ')


    args = parser.parse_args()
    reference = args.g
    LTR_finder_parallel_Home = args.ltrfinder_home
    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    recover = args.recover
    miu = args.miu

    tmp_output_dir = os.path.abspath(tmp_output_dir) 

    log = Logger(tmp_output_dir + '/HiTE_ltr.log', level='debug')

    is_recover = False
    recover = int(recover)
    if recover == 1:
        is_recover = True

    # 1.重命名reference文件
    ref_rename_path = tmp_output_dir + '/genome.rename.fa'
    rename_reference(reference, ref_rename_path)

    resut_file = tmp_output_dir + '/genome_all.fa.harvest.scn'
    if not is_recover or not file_exist(resut_file):
        log.logger.info('Start step2.1: Running LTR_harvest and LTR_finder_parallel')
        run_LTR_harvest(ref_rename_path, tmp_output_dir, threads, LTR_finder_parallel_Home, log)
    else:
        log.logger.info(resut_file + ' exists, skip...')

    # resut_file = ref_rename_path + '.finder.combine.scn'
    # if not is_recover or not file_exist(resut_file):
    #     starttime = time.time()
    #     log.logger.info('Start step2.2: Running LTR finder parallel to obtain candidate LTRs')

    #     # 运行LTR_finder_parallel来获取候选的LTR序列
    #     # 2.运行LTR_finder_parallel
    #     LTR_finder_parallel_command = 'perl ' + LTR_finder_parallel_Home + '/LTR_FINDER_parallel -harvest_out -seq ' + ref_rename_path + ' -threads ' + str(threads)
    #     log.logger.debug('cd ' + tmp_output_dir + ' && ' + LTR_finder_parallel_command + ' > /dev/null 2>&1')
    #     os.system('cd ' + tmp_output_dir + ' && ' + LTR_finder_parallel_command + ' > /dev/null 2>&1')

    #     endtime = time.time()
    #     dtime = endtime - starttime
    #     log.logger.info("Running time of LTR finder parallel: %.8s s" % (dtime))
    # else:
    #     log.logger.info(resut_file + ' exists, skip...')

    # 合并LTR_harvest+LTR_finder结果，输入到LTR_retriever
    ltrharvest_output = tmp_output_dir + '/genome_all.fa.harvest.scn'
    ltrfinder_output = ref_rename_path + '.finder.combine.scn'
    ltr_output = tmp_output_dir + '/genome_all.fa.rawLTR.scn'
    os.system('cat ' + ltrharvest_output + ' ' + ltrfinder_output + ' > ' + ltr_output)

    resut_file = ref_rename_path + '.LTRlib.fa'
    if not is_recover or not file_exist(resut_file):
        starttime = time.time()
        log.logger.info('Start step2.3: run LTR_retriever to get confident LTR')
        run_LTR_retriever(ref_rename_path, tmp_output_dir, threads, miu, log)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of step2.3: %.8s s" % (dtime))
    else:
        log.logger.info(resut_file + ' exists, skip...')

    confident_ltr_cut_path = tmp_output_dir + '/confident_ltr_cut.fa'
    rename_fasta(resut_file, confident_ltr_cut_path, 'LTR')
    #os.system('cp ' + resut_file + ' ' + confident_ltr_cut_path)




