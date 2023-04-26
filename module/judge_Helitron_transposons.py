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
    parser.add_argument('-g', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('--seqs', metavar='seqs',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa/longest_repeats_0.flanked.fa')
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--HSDIR', metavar='HSDIR',
                        help='e.g., /home/hukang/repeat_detect_tools/TrainingSet')
    parser.add_argument('--HSJAR', metavar='HSJAR',
                        help='e.g., /home/hukang/repeat_detect_tools/HelitronScanner/HelitronScanner.jar')
    parser.add_argument('--sh_dir', metavar='sh_dir',
                        help='e.g., /home/hukang/HiTE/module')
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

    args = parser.parse_args()

    reference = args.g
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

    # 将软链接路径转换绝对路径
    reference = os.path.realpath(reference)
    longest_repeats_flanked_path = os.path.realpath(longest_repeats_flanked_path)

    if debug is None:
        debug = 0
    else:
        debug = int(debug)

    is_recover = False
    recover = int(recover)
    if recover == 1:
        is_recover = True

    tmp_output_dir = os.path.abspath(tmp_output_dir) 

    log = Logger(tmp_output_dir+'/HiTE.log', level='debug')

    # 取10条全长拷贝两端flanking 50bp以包含Helitron边界
    candidate_helitron_path = tmp_output_dir + '/candidate_helitron_' + str(ref_index) + '.fa'
    resut_file = candidate_helitron_path
    if not is_recover or not file_exist(resut_file):
        # 运行helitronscanner
        HS_temp_dir = tmp_output_dir + '/HS_temp'
        if not os.path.exists(HS_temp_dir):
            os.makedirs(HS_temp_dir)
        candidate_helitron_path = tmp_output_dir + '/candidate_helitron_' + str(ref_index) + '.fa'
        multi_process_helitronscanner(longest_repeats_flanked_path, candidate_helitron_path, sh_dir, HS_temp_dir, HSDIR, HSJAR, threads)
        candidate_helitron_contignames, candidate_helitron_contigs = read_fasta(candidate_helitron_path)
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

        # # 用HelitronScanner再次过滤出至少有两个拷贝具有Helitron结构的序列
        # # 1. 先将序列获取拷贝，并扩展50bp
        # # 2. 用HelitronScanner运行拷贝文件，至少有两个以上的拷贝需要具有Helitron结构
        # temp_dir = tmp_output_dir + '/' + TE_type + '_copies_' + str(ref_index)
        # HS_temp_dir = tmp_output_dir + '/HS_temp'
        # os.system('rm -rf ' + HS_temp_dir)
        # if not os.path.exists(HS_temp_dir):
        #     os.makedirs(HS_temp_dir)
        #
        # jobs = {}
        # job_id = 0
        # candidate_count = 0
        # for name in os.listdir(temp_dir):
        #     if name.endswith('.blast.bed.fa'):
        #         copy_file = temp_dir + '/' + name
        #         prefix = name.split('.fa.blast.bed.fa')[0]
        #         jobs[job_id] = (copy_file, prefix)
        #         job_id += 1
        #
        # ex = ProcessPoolExecutor(threads)
        # objs = []
        # for job_id in jobs.keys():
        #     job = jobs[job_id]
        #     obj = ex.submit(run_HelitronScanner_v1, sh_dir, HS_temp_dir, job[0], HSDIR, HSJAR, job[1])
        #     objs.append(obj)
        # ex.shutdown(wait=True)
        #
        # confident_helitron_names = []
        # for obj in as_completed(objs):
        #     cur_candidate_Helitrons, prefix = obj.result()
        #     if len(cur_candidate_Helitrons) <= 1:
        #         continue
        #     # cur_candidate_file = temp_dir + '/' + prefix + '.HS.fa'
        #     # store_fasta(cur_candidate_Helitrons, cur_candidate_file)
        #     confident_helitron_names.append(prefix)
        #
        # raw_helitron_names, raw_helitron_contigs = read_fasta(candidate_helitron_cons)
        # confident_helitron_contigs = {}
        # for name in confident_helitron_names:
        #     seq = raw_helitron_contigs[name]
        #     confident_helitron_contigs[name] = seq
        # store_fasta(confident_helitron_contigs, confident_helitron_path)


    else:
        log.logger.info(resut_file + ' exists, skip...')


