#-- coding: UTF-8 --
import argparse
import codecs
import os
import sys

import json
import time

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, multi_process_helitronscanner, get_copies, multi_process_align, \
    flank_region_align_v1, multi_process_EAHelitron, Logger, flanking_copies, rename_fasta


def cut_reference(fasta_path, line_len):
    tmp_fasta_path = fasta_path + ".helitron.tmp"
    contigNames, contigs = read_fasta(fasta_path)
    with open(tmp_fasta_path, 'w') as f_w:
        for contigName in contigNames:
            contig = contigs[contigName]
            start = 0
            end = len(contig)
            while start < end:
                seg = contig[start:start+line_len]
                line = '>' + contigName + '-' + str(start) + '\n' + seg + '\n'
                f_w.write(line)
                start += line_len
    f_w.close()
    return tmp_fasta_path

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('--seqs', metavar='seqs',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa/longest_repeats_0.flanked.fa')
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    # parser.add_argument('--HSDIR', metavar='HSDIR',
    #                     help='e.g., /public/home/hpc194701009/repeat_detect_tools/TrainingSet')
    # parser.add_argument('--HSJAR', metavar='HSJAR',
    #                     help='e.g., /public/home/hpc194701009/repeat_detect_tools/HelitronScanner/HelitronScanner.jar')
    parser.add_argument('--EAHelitron', metavar='EAHelitron',
                        help='e.g., /public/home/hpc194701009/repeat_detect_tools/EAHelitron-master')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa')
    parser.add_argument('--flanking_len', metavar='flanking_len',
                        help='e.g., 50')
    parser.add_argument('--ref_index', metavar='ref_index',
                        help='e.g., 0')

    args = parser.parse_args()

    reference = args.g
    longest_repeats_flanked_path = args.seqs
    threads = int(args.t)
    # HSDIR = args.HSDIR
    # HSJAR = args.HSJAR
    tmp_output_dir = args.tmp_output_dir
    EAHelitron = args.EAHelitron
    ref_index = args.ref_index
    flanking_len = int(args.flanking_len)

    tmp_output_dir = os.path.abspath(tmp_output_dir) 

    log = Logger(tmp_output_dir+'/HiTE.log', level='debug')

    #sh_dir = os.getcwd() + '/'

    # log.logger.info('loading ' + longest_repeats_flanked_copies_file)
    # file = open(longest_repeats_flanked_copies_file, 'r')
    # js = file.read()
    # all_copies = json.loads(js)

    # 取10条全长拷贝两端flanking 50bp以包含Helitron边界

    #步骤： 1.多线程运行EAHelitron，获得候选Helitron序列
    candidate_helitron_path = tmp_output_dir + '/candidate_helitron_'+str(ref_index)+'.fa'
    temp_dir = tmp_output_dir + '/helitron_tmp_'+str(ref_index)
    multi_process_EAHelitron(longest_repeats_flanked_path, flanking_len, candidate_helitron_path, temp_dir, EAHelitron, threads)

    #运行helitronscanner
    # candidate_helitron_path = tmp_output_dir + '/candidate_helitronscanner_' + str(ref_index) + '.fa'
    # multi_process_helitronscanner(all_copies, candidate_helitron_path, sh_dir, temp_dir, HSDIR, HSJAR, threads)
    # candidate_helitron_contignames, candidate_helitron_contigs = read_fasta(candidate_helitron_path)

    # 生成一致性序列
    candidate_helitron_cons = tmp_output_dir + '/candidate_helitron_' + str(ref_index) + '.cons.fa'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + candidate_helitron_path + ' -o ' + candidate_helitron_cons + ' -T 0 -M 0'
    os.system(cd_hit_command)
    candidate_helitron_contignames, candidate_helitron_contigs = read_fasta(candidate_helitron_cons)


    # 2.flanking candidate helitron 50bp,看flanking region是否有高同源性
    confident_copies = {}
    if len(candidate_helitron_contigs) > 0:
        log.logger.info('------flank Helitron copy and see if the flanking regions are repeated')
        starttime = time.time()
        # 我们将copies扩展50bp，一个orig_query_name对应一个文件，然后做自比对。
        # 解析每个自比对文件，判断C0与C1,C2...等拷贝的比对情况，如果有flanking区域包含在比对区域内，那么这条拷贝应该被抛弃，如果所有拷贝被抛弃，则该条序列应该是假阳性。
        flanking_len = 50
        similar_ratio = 0.1
        TE_type = 'helitron'
        confident_copies = flank_region_align_v1(candidate_helitron_cons, flanking_len, similar_ratio, reference,
                                                 TE_type, tmp_output_dir, threads, ref_index, log)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of flanking Helitron copy and see if the flanking regions are repeated: %.8s s" % (dtime))

    confident_helitron_path = tmp_output_dir + '/confident_helitron_'+str(ref_index)+'.fa'
    confident_helitron = {}
    for name in confident_copies.keys():
        copy_list = confident_copies[name]
        if len(copy_list) >= 2:
            confident_helitron[name] = candidate_helitron_contigs[name]
    store_fasta(confident_helitron, confident_helitron_path)
    rename_fasta(confident_helitron_path, confident_helitron_path, 'Helitron')
    #
    # confident_helitron_path = tmp_output_dir + '/confident_helitron_'+str(ref_index)+'.fa'
    # confident_helitron_rename_path = tmp_output_dir + '/confident_helitron_'+str(ref_index)+'.rename.fa'
    # rename_fasta(confident_helitron_path, confident_helitron_rename_path)
    #
    # candidate_helitron_rename_consensus = tmp_output_dir + '/confident_helitron_'+str(ref_index)+'.rename.cons.fa'
    # cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + confident_helitron_rename_path + ' -o ' + candidate_helitron_rename_consensus + ' -T 0 -M 0'
    # os.system(cd_hit_command)


