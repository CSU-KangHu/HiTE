#!/usr/bin/env python
import argparse
import os
import shutil
import subprocess
import sys
import uuid

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, multi_process_helitronscanner, multi_process_EAHelitron, \
    Logger, rename_fasta, file_exist, flank_region_align_v5, create_or_clear_directory, copy_files, \
    clean_old_tmp_files_by_dir, update_prev_TE, lib_add_prefix, store_fasta


def run_Helitron_detection(tmp_output_dir, longest_repeats_flanked_path, prev_TE, ref_index, is_recover, threads, debug,
                           flanking_len, reference, split_ref_dir, all_low_copy_helitron, min_TE_len, log):
    HSDIR = cur_dir + '/bin/HelitronScanner/TrainingSet'
    HSJAR = cur_dir + '/bin/HelitronScanner/HelitronScanner.jar'
    sh_dir = cur_dir + '/bin'
    # EAHelitron = cur_dir + '/bin/EAHelitron-master'
    subset_script_path = cur_dir + '/tools/ready_for_MSA.sh'

    candidate_helitron_path = tmp_output_dir + '/candidate_helitron_' + str(ref_index) + '.fa'
    resut_file = candidate_helitron_path
    if not is_recover or not file_exist(resut_file):
        # run helitronscanner
        HS_temp_dir = tmp_output_dir + '/HS_temp'
        if not os.path.exists(HS_temp_dir):
            os.makedirs(HS_temp_dir)
        candidate_helitronscanner_path = tmp_output_dir + '/candidate_helitron_' + str(
            ref_index) + '.HelitronScanner.fa'
        multi_process_helitronscanner(longest_repeats_flanked_path, candidate_helitronscanner_path, sh_dir, HS_temp_dir,
                                      HSDIR, HSJAR, threads, debug)
        # candidate_helitron_contignames, candidate_helitron_contigs = read_fasta(candidate_helitron_path)
        if not debug:
            shutil.rmtree(HS_temp_dir)

        # # run EAHelitron
        # EA_temp_dir = tmp_output_dir + '/EA_temp'
        # if not os.path.exists(EA_temp_dir):
        #     os.makedirs(EA_temp_dir)
        # candidate_eahelitron_path = tmp_output_dir + '/candidate_helitron_' + str(ref_index) + '.EAHelitron.fa'
        # multi_process_EAHelitron(longest_repeats_flanked_path, flanking_len, candidate_eahelitron_path, EA_temp_dir,
        #                          EAHelitron, threads)
        # if not debug:
        #     os.system('rm -rf ' + EA_temp_dir)

        # Combine results from HelitronScanner and EAHelitron.
        with open(candidate_helitron_path, 'w') as outfile:
            subprocess.run(['cat', candidate_helitronscanner_path], stdout=outfile)

        # os.system('cat ' + candidate_helitronscanner_path + ' > ' + candidate_helitron_path)
        # os.system('cat ' + candidate_eahelitron_path + ' >> ' + candidate_helitron_path)

    candidate_helitron_cons = candidate_helitron_path + '.cons'
    resut_file = candidate_helitron_cons
    if not is_recover or not file_exist(resut_file):
        log.logger.info('------clustering candidate Helitron')
        # cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
        #                  + ' -G 0 -g 1 -A 80 -i ' + candidate_helitron_path + ' -o ' + candidate_helitron_cons + ' -T 0 -M 0' + ' > /dev/null 2>&1'
        # os.system(cd_hit_command)
        cd_hit_command = [
            'cd-hit-est',
            '-aS', '0.95',
            '-aL', '0.95',
            '-c', '0.8',
            '-G', '0',
            '-g', '1',
            '-A', '80',
            '-i', candidate_helitron_path,
            '-o', candidate_helitron_cons,
            '-T', '0',
            '-M', '0'
        ]
        subprocess.run(cd_hit_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Apply homology filtering to the results identified by HelitronScanner.
    confident_helitron_path = tmp_output_dir + '/confident_helitron_' + str(ref_index) + '.fa'
    resut_file = confident_helitron_path
    delete_files = []
    if not is_recover or not file_exist(resut_file):
        flanking_len = 50
        TE_type = 'helitron'
        # Multiple iterations are performed to find more accurate boundaries.
        iter_num = 3
        input_file = candidate_helitron_cons
        for i in range(iter_num):
            result_type = 'cons'
            output_file = tmp_output_dir + '/confident_helitron_' + str(ref_index) + '.r' + str(i) + '.fa'
            cur_resut_file = output_file
            delete_files.append(output_file)
            if not is_recover or not file_exist(cur_resut_file):
                flank_region_align_v5(input_file, output_file, flanking_len, reference, split_ref_dir, TE_type,
                                      tmp_output_dir, threads, ref_index, log, subset_script_path, 1, debug,
                                      i, all_low_copy_helitron, result_type)
            input_file = output_file
        cur_confident_helitron_path = tmp_output_dir + '/confident_helitron_' + str(ref_index) + '.r' + str(
            iter_num - 1) + '.fa'
        cur_confident_helitron_cons = tmp_output_dir + '/confident_helitron_' + str(ref_index) + '.r' + str(
            iter_num - 1) + '.cons.fa'
        delete_files.append(cur_confident_helitron_path)
        delete_files.append(cur_confident_helitron_cons)
        # clustering
        # cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
        #                  + ' -G 0 -g 1 -A 80 -i ' + cur_confident_helitron_path + ' -o ' + cur_confident_helitron_cons + ' -T 0 -M 0' + ' > /dev/null 2>&1'
        # os.system(cd_hit_command)
        cd_hit_command = [
            'cd-hit-est',
            '-aS', '0.95',
            '-aL', '0.95',
            '-c', '0.8',
            '-G', '0',
            '-g', '1',
            '-A', '80',
            '-i', cur_confident_helitron_path,
            '-o', cur_confident_helitron_cons,
            '-T', '0',
            '-M', '0'
        ]
        subprocess.run(cd_hit_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # 去除 Helitron header 中包含 #, and filter TE by length
        cur_names, cur_contigs = read_fasta(cur_confident_helitron_cons)
        new_cur_contigs = {}
        for name in cur_names:
            if len(cur_contigs[name]) >= min_TE_len:
                new_cur_contigs[name.split('#')[0]] = cur_contigs[name]
        store_fasta(new_cur_contigs, cur_confident_helitron_cons)
        rename_fasta(cur_confident_helitron_cons, confident_helitron_path, 'Helitron_' + str(ref_index))

        # remove temp files
        for delete_file in delete_files:
            if os.path.exists(delete_file):
                os.remove(delete_file)

    else:
        log.logger.info(resut_file + ' exists, skip...')

    raw_name = os.path.basename(reference).split('.')[0]
    if file_exist(resut_file):
        lib_add_prefix(resut_file, raw_name)

    update_prev_TE(prev_TE, resut_file)
    # os.system('cat ' + resut_file + ' >> ' + prev_TE)

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE Helitron module...')
    parser.add_argument('--seqs', metavar='seqs',
                        help='Please enter the result of de novo TE searching in HiTE, typically named longest_repeats_*.fa. Please provide the absolute path.')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number.')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--flanking_len', metavar='flanking_len',
                        help='The flanking length of candidates to find the true boundaries.')
    parser.add_argument('--ref_index', metavar='ref_index',
                        help='The current split genome index.')
    parser.add_argument('--recover', metavar='recover',
                        help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false.')
    parser.add_argument('--debug', metavar='recover',
                        help='Open debug mode, and temporary files will be kept, 1: true, 0: false.')
    parser.add_argument('-r', metavar='Reference path',
                        help='Input Reference path.')
    parser.add_argument('--split_ref_dir', metavar='Split Reference path',
                        help='Please enter the directory of the split genome.')
    parser.add_argument('--prev_TE', metavar='prev_TE',
                        help='TEs fasta file that has already been identified. Please use the absolute path.')
    parser.add_argument('--all_low_copy_helitron', metavar='all_low_copy_helitron',
                        help='all low copy helitron path, to recover helitron using pan-genome')
    parser.add_argument('--min_TE_len', metavar='min_TE_len',
                        help='The minimum TE length')
    parser.add_argument('-w', '--work_dir', nargs="?", default='/tmp', help="The temporary work directory for HiTE.")

    args = parser.parse_args()

    longest_repeats_flanked_path = args.seqs
    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    ref_index = args.ref_index
    flanking_len = int(args.flanking_len)
    recover = args.recover
    debug = args.debug
    reference = args.r
    split_ref_dir = args.split_ref_dir
    prev_TE = args.prev_TE
    all_low_copy_helitron = args.all_low_copy_helitron
    min_TE_len = int(args.min_TE_len)
    work_dir = args.work_dir
    work_dir = os.path.abspath(work_dir)

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

    if tmp_output_dir is None:
        tmp_output_dir = os.getcwd()

    tmp_output_dir = os.path.abspath(tmp_output_dir)
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)

    log = Logger(tmp_output_dir+'/HiTE_helitron.log', level='debug')
    all_low_copy_helitron = os.path.abspath(all_low_copy_helitron)
    if not os.path.exists(all_low_copy_helitron):
        os.system('touch ' + all_low_copy_helitron)

    # clean_old_tmp_files_by_dir('/tmp')

    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = os.path.join(work_dir, 'judge_Helitron_transposons_' + str(unique_id))
    try:
        create_or_clear_directory(temp_dir)

        run_Helitron_detection(temp_dir, longest_repeats_flanked_path, prev_TE, ref_index, is_recover, threads, debug,
                               flanking_len, reference, split_ref_dir, all_low_copy_helitron, min_TE_len, log)

        # 计算完之后将结果拷贝回输出目录
        copy_files(temp_dir, tmp_output_dir)

    except Exception as e:
        # 如果出现异常，打印错误信息并删除临时目录
        print(f"An error occurred: {e}")
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        raise  # 重新抛出异常，以便上层代码可以处理

    else:
        # 如果没有异常，删除临时目录
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)