#-- coding: UTF-8 --
import argparse
import codecs
import os
import sys

import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, multi_process_helitronscanner, get_copies, multi_process_align, \
    flank_region_align_v1, multi_process_EAHelitron, Logger, flanking_copies, rename_fasta, file_exist, \
    flank_region_align_v3, run_HelitronScanner_v1


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
    parser.add_argument('--seqs', metavar='seqs',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa/longest_repeats_0.flanked.fa')
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--HSDIR', metavar='HSDIR',
                        help='e.g., /home/hukang/repeat_detect_tools/TrainingSet')
    parser.add_argument('--HSJAR', metavar='HSJAR',
                        help='e.g., /home/hukang/repeat_detect_tools/HelitronScanner/HelitronScanner.jar')
    parser.add_argument('--sh_dir', metavar='sh_dir',
                        help='e.g., /home/hukang/HiTE/modcd ule')
    parser.add_argument('--member_script_path', metavar='member_script_path',
                        help='e.g., /home/hukang/HiTE/tools/make_fasta_from_blast.sh')
    parser.add_argument('--subset_script_path', metavar='subset_script_path',
                        help='e.g., /home/hukang/HiTE/tools/ready_for_MSA.sh')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa')
    parser.add_argument('--flanking_len', metavar='flanking_len',
                        help='e.g., 50')
    parser.add_argument('--ref_index', metavar='ref_index',
                        help='e.g., 0')
    parser.add_argument('--recover', metavar='recover',
                        help='e.g., 0')
    parser.add_argument('--debug', metavar='recover',
                        help='e.g., 1')
    parser.add_argument('-r', metavar='Reference path',
                        help='input Reference path')

    args = parser.parse_args()

    longest_repeats_flanked_path = args.seqs
    threads = int(args.t)
    HSDIR = args.HSDIR
    HSJAR = args.HSJAR
    sh_dir = args.sh_dir
    member_script_path = args.member_script_path
    subset_script_path = args.subset_script_path
    tmp_output_dir = args.tmp_output_dir
    ref_index = args.ref_index
    flanking_len = int(args.flanking_len)
    recover = args.recover
    debug = args.debug
    reference = args.r

    # 将软链接路径转换绝对路径
    longest_repeats_flanked_path = os.path.realpath(longest_repeats_flanked_path)
    reference = os.path.realpath(reference)

    if debug is None:
        debug = 0
    else:
        debug = int(debug)

    is_recover = False
    recover = int(recover)
    if recover == 1:
        is_recover = True

    tmp_output_dir = os.path.abspath(tmp_output_dir) 

    log = Logger(tmp_output_dir+'/HiTE_helitron.log', level='debug')

    # 取10条全长拷贝两端flanking 50bp以包含Helitron边界
    candidate_helitron_path = tmp_output_dir + '/candidate_helitron_' + str(ref_index) + '.fa'
    resut_file = candidate_helitron_path
    if not is_recover or not file_exist(resut_file):
        # 运行helitronscanner
        HS_temp_dir = tmp_output_dir + '/HS_temp'
        if not os.path.exists(HS_temp_dir):
            os.makedirs(HS_temp_dir)
        candidate_helitron_path = tmp_output_dir + '/candidate_helitron_' + str(ref_index) + '.fa'
        multi_process_helitronscanner(longest_repeats_flanked_path, candidate_helitron_path, sh_dir, HS_temp_dir, HSDIR, HSJAR, threads, debug)
        candidate_helitron_contignames, candidate_helitron_contigs = read_fasta(candidate_helitron_path)
        if not debug:
            os.system('rm -rf ' + HS_temp_dir)

    # 生成一致性序列
    candidate_helitron_cons = tmp_output_dir + '/candidate_helitron_' + str(ref_index) + '.cons.fa'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + candidate_helitron_path + ' -o ' + candidate_helitron_cons + ' -T 0 -M 0'
    os.system(cd_hit_command)
    #rename_fasta(candidate_helitron_cons, candidate_helitron_cons, 'Helitron')

    # 对HelitronScanner识别的结果使用同源性过滤方法过滤
    confident_helitron_path = tmp_output_dir + '/confident_helitron_' + str(ref_index) + '.fa'
    resut_file = confident_helitron_path
    if not is_recover or not file_exist(resut_file):
        flanking_len = 50
        similar_ratio = 0.2
        TE_type = 'helitron'
        # 多轮迭代是为了找到更加准确的边界
        iter_num = 3
        input_file = candidate_helitron_cons
        for i in range(iter_num):
            result_type = 'cons'
            output_file = tmp_output_dir + '/confident_helitron_' + str(ref_index) + '.r' + str(i) + '.fa'
            flank_region_align_v3(input_file, output_file, flanking_len, similar_ratio, reference, TE_type,
                                  tmp_output_dir, threads,
                                  ref_index, log, member_script_path, subset_script_path, 1, debug, result_type)
            input_file = output_file
        cur_confident_helitron_path = tmp_output_dir + '/confident_helitron_' + str(ref_index) + '.r' + str(iter_num - 1) + '.fa'
        rename_fasta(cur_confident_helitron_path, confident_helitron_path, 'Helitron')
    else:
        log.logger.info(resut_file + ' exists, skip...')


