#!/usr/bin/env python
import argparse
import os
import shutil
import sys
import time
import uuid

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, multi_process_tsd, rename_fasta, file_exist, \
    run_itrsearch, get_short_tir_contigs, flank_region_align_v5, multi_process_tsd_v1, remove_no_tirs, \
    create_or_clear_directory, copy_files, clean_old_tmp_files_by_dir, update_prev_TE, lib_add_prefix


def is_transposons(filter_dup_path, reference, threads, tmp_output_dir, ref_index, log, subset_script_path, plant,
                   debug, TRsearch_dir, split_ref_dir, all_low_copy_tir, min_TE_len, is_recover):
    log.logger.info('determine true TIR')
    log.logger.info('------flank TIR copy and see if the flanking regions are repeated')
    starttime = time.time()
    flanking_len = 50
    TE_type = 'tir'

    # Utilizing multi-alignment sequences with a sliding window mode to filter out false positive sequences.
    # Multiple iterations are performed to achieve more accurate boundary identification.
    delete_files = []
    iter_num = 3
    input_file = filter_dup_path
    for i in range(iter_num):
        result_type = 'cons'
        output_file = tmp_output_dir + '/confident_tir_' + str(ref_index) + '.r' + str(i) + '.fa'
        delete_files.append(output_file)
        resut_file = output_file
        if not is_recover or not file_exist(resut_file):
            flank_region_align_v5(input_file, output_file, flanking_len, reference, split_ref_dir, TE_type,
                                  tmp_output_dir, threads, ref_index, log, subset_script_path, plant, debug,
                                  i, all_low_copy_tir, result_type)
        input_file = output_file

    confident_tir_path = tmp_output_dir + '/confident_tir_' + str(ref_index) + '.r' + str(iter_num-1) + '.fa'
    delete_files.append(confident_tir_path)

    # filter TE by length
    confident_tir_names, confident_tir_contigs = read_fasta(confident_tir_path)
    filter_confident_tir_contigs = {}
    for name in confident_tir_contigs.keys():
        if len(confident_tir_contigs[name]) >= min_TE_len:
            filter_confident_tir_contigs[name] = confident_tir_contigs[name]
    store_fasta(filter_confident_tir_contigs, confident_tir_path)

    final_confident_tir_path = tmp_output_dir + '/confident_tir_' + str(ref_index) + '.fa'
    rename_fasta(confident_tir_path, final_confident_tir_path, 'TIR_' + str(ref_index))

    # remove temp files
    for delete_file in delete_files:
        if os.path.exists(delete_file):
            os.remove(delete_file)

    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of flanking TIR copy and see if the flanking regions are repeated: %.8s s" % (dtime))

def run_TIR_detection(tmp_output_dir, longest_repeats_flanked_path, reference, prev_TE, flanking_len, threads,
                      debug, split_ref_dir, all_low_copy_tir, plant, ref_index, min_TE_len, is_recover, log):
    TRsearch_dir = cur_dir + '/tools'
    subset_script_path = cur_dir + '/tools/ready_for_MSA.sh'

    tir_tsd_path = tmp_output_dir + '/tir_tsd_' + str(ref_index) + '.fa'
    resut_file = tir_tsd_path
    if not is_recover or not file_exist(resut_file):
        # Retrieve candidate sequences with TIR+TSD structure from De novo TE searching.
        log.logger.info('------get TIR+TSD in copies of candidate TIR')
        starttime = time.time()
        tir_tsd_dir = tmp_output_dir + '/tir_tsd_temp_' + str(ref_index)
        # multi_process_tsd(longest_repeats_flanked_path, tir_tsd_path, tir_tsd_dir, flanking_len, threads, TRsearch_dir, plant)
        multi_process_tsd_v1(longest_repeats_flanked_path, tir_tsd_path, tir_tsd_dir, flanking_len, threads,
                             TRsearch_dir, plant)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of getting TSD in copies of candidate TIR: %.8s s" % (dtime))

    tir_tsd_cons = tmp_output_dir + '/tir_tsd_' + str(ref_index) + '.cons.fa'
    resut_file = tir_tsd_cons
    if not is_recover or not file_exist(resut_file):
        log.logger.info('------clustering candidate TIR')
        starttime = time.time()
        cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                         + ' -G 0 -g 1 -A 80 -i ' + tir_tsd_path + ' -o ' + tir_tsd_cons + ' -T 0 -M 0' + ' > /dev/null 2>&1'
        os.system(cd_hit_command)
        endtime = time.time()
        dtime = endtime - starttime
        log.logger.info("Running time of clustering candidate TIR: %.8s s" % (dtime))

    confident_tir_path = tmp_output_dir + '/confident_tir_' + str(ref_index) + '.fa'
    resut_file = confident_tir_path
    if not is_recover or not file_exist(resut_file):
        # Utilize homologous boundary search method to determine the authenticity of TE sequences.
        is_transposons(tir_tsd_cons, reference, threads, tmp_output_dir, ref_index, log,
                       subset_script_path, plant, debug, TRsearch_dir, split_ref_dir, all_low_copy_tir, min_TE_len, is_recover)
    else:
        log.logger.info(resut_file + ' exists, skip...')


    raw_name = os.path.basename(reference).split('.')[0]
    if file_exist(resut_file):
        lib_add_prefix(resut_file, raw_name)

    update_prev_TE(prev_TE, resut_file)
    # os.system('cat ' + resut_file + ' >> ' + prev_TE)


if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE TIR module...')
    parser.add_argument('--seqs', metavar='seqs',
                        help='Please enter the result of de novo TE searching in HiTE, typically named longest_repeats_*.fa. Please provide the absolute path.')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number.')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--tandem_region_cutoff', metavar='tandem_region_cutoff',
                        help='Cutoff of the candidates regarded as tandem region.')
    parser.add_argument('--flanking_len', metavar='flanking_len',
                        help='The flanking length of candidates to find the true boundaries.')
    parser.add_argument('--plant', metavar='plant',
                        help='Is it a plant genome, 1: true, 0: false.')
    parser.add_argument('--ref_index', metavar='ref_index',
                        help='The current split genome index.')
    parser.add_argument('--recover', metavar='recover',
                        help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false.')
    parser.add_argument('--debug', metavar='debug',
                        help='Open debug mode, and temporary files will be kept, 1: true, 0: false.')
    parser.add_argument('-r', metavar='Reference path',
                        help='Input Reference path.')
    parser.add_argument('--split_ref_dir', metavar='Split Reference path',
                        help='Please enter the directory of the split genome.')
    parser.add_argument('--prev_TE', metavar='prev_TE',
                        help='TEs fasta file that has already been identified. Please use the absolute path.')
    parser.add_argument('--all_low_copy_tir', metavar='all_low_copy_tir',
                        help='all low copy tir path, to recover tir using pan-genome')
    parser.add_argument('--min_TE_len', metavar='min_TE_len',
                        help='The minimum TE length')
    parser.add_argument('-w', '--work_dir', nargs="?", default='/tmp', help="The temporary work directory for HiTE.")

    args = parser.parse_args()
    longest_repeats_flanked_path = args.seqs
    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    flanking_len = int(args.flanking_len)
    plant = int(args.plant)
    tandem_region_cutoff = float(args.tandem_region_cutoff)
    ref_index = args.ref_index
    recover = args.recover
    debug = args.debug
    reference = args.r
    split_ref_dir = args.split_ref_dir
    prev_TE = args.prev_TE
    all_low_copy_tir = args.all_low_copy_tir
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

    log = Logger(tmp_output_dir+'/HiTE_tir.log', level='debug')

    all_low_copy_tir = os.path.abspath(all_low_copy_tir)
    if not os.path.exists(all_low_copy_tir):
        os.system('touch ' + all_low_copy_tir)

    # clean_old_tmp_files_by_dir('/tmp')

    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = os.path.join(work_dir, 'judge_TIR_transposons_' + str(unique_id))
    try:
        create_or_clear_directory(temp_dir)

        run_TIR_detection(temp_dir, longest_repeats_flanked_path, reference, prev_TE, flanking_len, threads, debug,
                          split_ref_dir, all_low_copy_tir, plant, ref_index, min_TE_len, is_recover, log)

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
        if os.path.exists(temp_dir) and debug != 1:
            shutil.rmtree(temp_dir)