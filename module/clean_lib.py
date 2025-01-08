#!/usr/bin/env python
import argparse
import os
import re
import sys


cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import Logger



if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--debug', metavar='recover',
                        help='Open debug mode, and temporary files will be kept, 1: true, 0: false.')

    args = parser.parse_args()

    tmp_output_dir = args.tmp_output_dir
    debug = int(args.debug)

    log = Logger(tmp_output_dir+'/HiTE_clean.log', level='debug')


    # remove temp files and directories
    if debug == 0:
        keep_files_temp = []
        keep_files = [
            'chr_name\\.map',
            'genome\\.rename\\.fa',
            'genome\\.fa\\.masked',
            'genome\\.rename\\.fa\\.pass\\.list',
            '.*\\.scn',
            'all_TE\\.fa',
            'low_confident_TE\\.cons\\.fa',
            'genome\\.rename\\.fa\\.LTRlib\\.fa',
            'TE_merge_tmp\\.fa\\.classified',
            'confident_TE\\.cons\\.fa',
            'confident_TE\\.cons\\.fa\\.domain',
            'confident_ltr_cut\\.fa',
            'confident_ltr\\.internal\\.fa',
            'confident_ltr\\.terminal\\.fa',
            'intact_LTR\\.list',
            'intact_LTR\\.fa',
            'intact_LTR\\.fa\\.classified',
            'confident_TE\\.cons\\.fa\\.classified',
            'longest_repeats(_\\d+)?\\.flanked\\.fa',
            'longest_repeats(_\\d+)?\\.fa',
            'confident_tir(_\\d+)?\\.fa',
            'confident_helitron(_\\d+)?\\.fa',
            'confident_non_ltr(_\\d+)?\\.fa',
            'confident_other(_\\d+)?\\.fa',
            'repbase\\.out',
            'HiTE\\.out',
            'HiTE\\.tbl',
            'HiTE\\.gff',
            'HiTE\\.cat\\.gz',
            'HiTE_intact\\.sorted\\.gff3',
            'BM_RM2\\.log',
            'BM_EDTA\\.log',
            'BM_HiTE\\.log',
            # 以下规则匹配 Nextflow 的隐藏文件
            '\\.command\\.begin',
            '\\.command\\.run',
            '\\.command\\.trace',
            '\\.command\\.sh',
            '\\.command\\.out',
            '\\.command\\.err',
            '\\.exitcode',
            '\\.command\\.log',
            '.*_low_copy.fa'
        ]

        all_files = os.listdir(tmp_output_dir)
        for filename in all_files:
            is_del = True
            for keep_file in keep_files:
                is_match = re.match(keep_file+'$', filename)
                if is_match is not None:
                    is_del = False
                    break
            if is_del:
                os.system('rm -rf ' + tmp_output_dir + '/' + filename)


