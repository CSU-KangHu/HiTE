#!/usr/bin/env python
import argparse
import os
import time

current_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = os.path.join(current_folder, ".")

from Util import Logger, deredundant_for_LTR_v5, ReassignInconsistentLabels


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
    parser.add_argument("--pan_terminal_tmp_lib", type=str, help="pan terminal lib.")
    parser.add_argument("--pan_internal_tmp_lib", type=str, help="pan internal lib.")
    parser.add_argument("--threads", type=int, help="Number of threads to use.")
    parser.add_argument("--output_dir", nargs="?", default=os.getcwd(),
                        help="Output directory (default: current working directory).")

    # 解析参数
    args = parser.parse_args()
    pan_terminal_tmp_lib = args.pan_terminal_tmp_lib
    pan_internal_tmp_lib = args.pan_internal_tmp_lib
    threads = args.threads

    # 处理输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    log = Logger(output_dir + '/panHiTE.log', level='debug')

    # 调用冗余去除和库合并函数
    remove_redundancy(pan_terminal_tmp_lib, pan_internal_tmp_lib, output_dir, threads, log)
