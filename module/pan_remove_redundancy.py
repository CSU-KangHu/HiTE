#!/usr/bin/env python
import argparse
import os
import shutil
import time
import uuid

current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")

from Util import Logger, deredundant_for_LTR_v5, ReassignInconsistentLabels, split_internal_out, \
    create_or_clear_directory, copy_files


def remove_redundancy(pan_terminal_tmp_lib, pan_internal_tmp_lib, output_dir, threads, log):
    panTE_lib = os.path.join(output_dir, 'panTE.fa')
    """
    对所有基因组生成的 TE library 去冗余并合并最终库
    """
    # Step 3.1 Remove LTR terminal redundancy
    pan_terminal_tmp_lib_cons = pan_terminal_tmp_lib + '.cons'
    terminal_coverage_threshold = 0.95
    starttime = time.time()
    log.logger.info('Start Removing LTR terminal redundancy')
    # 去冗余处理 LTR terminal
    deredundant_for_LTR_v5(pan_terminal_tmp_lib, output_dir, threads, 'terminal', terminal_coverage_threshold, debug=0)
    endtime = time.time()
    log.logger.info(f"Running time of Removing LTR terminal redundancy: {endtime - starttime:.4f} s")

    # Step 3.2 Remove LTR internal redundancy
    pan_internal_tmp_lib_cons = pan_internal_tmp_lib + '.cons'
    internal_coverage_threshold = 0.8
    starttime = time.time()
    log.logger.info('Start Removing LTR internal redundancy')
    # 去冗余处理 LTR internal
    deredundant_for_LTR_v5(pan_internal_tmp_lib, output_dir, threads, 'internal', internal_coverage_threshold, debug=0)
    endtime = time.time()
    log.logger.info(f"Running time of Removing LTR internal redundancy: {endtime - starttime:.4f} s")

    # Step 3.3 合并最终 TE library
    log.logger.info('Merging final TE libraries')
    os.system(f'cat {pan_terminal_tmp_lib_cons} {pan_internal_tmp_lib_cons} > {panTE_lib}')

    # Reassign inconsistent classification labels
    ReassignInconsistentLabels(panTE_lib)
    log.logger.info('Redundancy removal and library merge completed')


if __name__ == "__main__":
    # 创建解析器
    parser = argparse.ArgumentParser(description="panHiTE remove redundancy.")
    parser.add_argument("--merge_te_file", type=str, help="merged pan te file.")
    parser.add_argument("--threads", type=int, help="Number of threads to use.")
    parser.add_argument("--output_dir", nargs="?", default=os.getcwd(),
                        help="Output directory (default: current working directory).")

    # 解析参数
    args = parser.parse_args()
    merge_te_file = args.merge_te_file
    threads = args.threads

    # 处理输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    log = Logger(output_dir + '/panHiTE.log', level='debug')

    # 创建本地临时目录，存储计算结果
    unique_id = uuid.uuid4()
    temp_dir = '/tmp/pan_remove_redundancy_' + str(unique_id)
    create_or_clear_directory(temp_dir)

    # 根据文件的header将LTR内部序列和其他元素区分开存储
    other_path, internal_path = split_internal_out(merge_te_file, temp_dir)
    # 调用冗余去除和库合并函数
    remove_redundancy(other_path, internal_path, temp_dir, threads, log)

    # 计算完之后将结果拷贝回输出目录
    copy_files(temp_dir, output_dir)

    # 删除临时目录
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)