import os
import time
import sys

from Util import Logger, deredundant_for_LTR_v5, ReassignInconsistentLabels


def remove_redundancy(pan_terminal_tmp_lib, pan_internal_tmp_lib, output_dir, threads, panTE_lib, log):
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
    # 使用 sys.argv 获取命令行参数
    if len(sys.argv) != 6:
        print(
            "Usage: python remove_redundancy.py <pan_terminal_tmp_lib> <pan_internal_tmp_lib> <output_dir> <threads> <panTE_lib>")
        sys.exit(1)

    pan_terminal_tmp_lib = sys.argv[1]
    pan_internal_tmp_lib = sys.argv[2]
    output_dir = sys.argv[3]
    threads = int(sys.argv[4])
    panTE_lib = sys.argv[5]

    if output_dir is None:
        output_dir = os.getcwd()
    output_dir = os.path.abspath(output_dir)
    log = Logger(output_dir + '/panHiTE.log', level='debug')

    # 调用冗余去除和库合并函数
    remove_redundancy(pan_terminal_tmp_lib, pan_internal_tmp_lib, output_dir, threads, panTE_lib, log)
