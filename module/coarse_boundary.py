#!/usr/bin/env python
import argparse
import os
import shutil
import sys
import uuid

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import Logger, flanking_seq, file_exist, determine_repeat_boundary_v5, create_or_clear_directory, copy_files, \
    clean_old_tmp_files_by_dir


def run_coarse_boundary(tmp_output_dir, cut_reference, reference, is_recover, ref_index,
                        prev_TE, fixed_extend_base_threshold, flanking_len, max_repeat_len,
                        thread, debug, log):
    repeats_path = cut_reference
    longest_repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.fa'
    resut_file = longest_repeats_path
    if not is_recover or not file_exist(resut_file):
        # -------------------------------This stage is used to do pairwise comparision, determine the repeat boundary-------------------------------
        determine_repeat_boundary_v5(repeats_path, longest_repeats_path, prev_TE, fixed_extend_base_threshold,
                                     max_repeat_len, tmp_output_dir, thread, ref_index, reference, debug)
    else:
        log.logger.info(resut_file + ' exists, skip...')

    longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'
    resut_file = longest_repeats_flanked_path
    if not is_recover or not file_exist(resut_file):
        flanking_seq(longest_repeats_path, longest_repeats_flanked_path, reference, flanking_len)
    else:
        log.logger.info(resut_file + ' exists, skip...')

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE De novo TE searching...')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='Input genome assembly path.')
    parser.add_argument('--prev_TE', metavar='prev_TE',
                        help='TEs fasta file that has already been identified. Please use the absolute path.')
    parser.add_argument('--fixed_extend_base_threshold', metavar='fixed_extend_base_threshold',
                        help='The length of variation can be tolerated during pairwise alignment.')
    parser.add_argument('--max_repeat_len', metavar='max_repeat_len', help='The maximum length of a single repeat.')
    parser.add_argument('--thread', metavar='thread_num',
                        help='Input thread num.')
    parser.add_argument('--flanking_len', metavar='flanking_len',
                        help='The flanking length of candidates to find the true boundaries.')
    parser.add_argument('--tandem_region_cutoff', metavar='tandem_region_cutoff',
                        help='Cutoff of the candidates regarded as tandem region.')
    parser.add_argument('--ref_index', metavar='ref_index',
                        help='The current split genome index.')
    parser.add_argument('-r', metavar='Genome assembly',
                        help='Input genome assembly path.')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--recover', metavar='recover',
                        help='Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false.')
    parser.add_argument('--debug', metavar='recover',
                        help='Open debug mode, and temporary files will be kept, 1: true, 0: false.')
    parser.add_argument('-w', '--work_dir', nargs="?", default='/tmp', help="The temporary work directory for HiTE.")

    args = parser.parse_args()
    cut_reference = args.g
    prev_TE = args.prev_TE
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
    work_dir = args.work_dir
    work_dir = os.path.abspath(work_dir)

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

    log = Logger(tmp_output_dir + '/HiTE_coarse.log', level='debug')

    # clean_old_tmp_files_by_dir('/tmp')

    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = os.path.join(work_dir, 'coarse_boundary_' + str(unique_id))
    try:
        create_or_clear_directory(temp_dir)
        # 运行计算任务
        run_coarse_boundary(temp_dir, cut_reference, reference, is_recover, ref_index,
                            prev_TE, fixed_extend_base_threshold, flanking_len, max_repeat_len,
                            thread, debug, log)

        # 计算完之后将结果拷贝回输出目录
        copy_files(temp_dir, tmp_output_dir)

    except Exception as e:
        # 如果出现异常，打印错误信息并删除临时目录
        print(f"An error occurred: {e}")
        # if os.path.exists(temp_dir):
        #     shutil.rmtree(temp_dir)
        raise  # 重新抛出异常，以便上层代码可以处理

    else:
        # 如果没有异常，删除临时目录
        if os.path.exists(temp_dir) and debug != 1:
            shutil.rmtree(temp_dir)



