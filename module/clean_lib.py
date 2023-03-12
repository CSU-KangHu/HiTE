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
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa')
    parser.add_argument('--debug', metavar='debug',
                        help='e.g., 1')
    parser.add_argument('--ref_name', metavar='ref_name',
                        help='e.g., ')

    args = parser.parse_args()

    tmp_output_dir = args.tmp_output_dir
    debug = int(args.debug)
    ref_name = args.ref_name

    log = Logger(tmp_output_dir+'/HiTE.log', level='debug')


    # remove temp files and directories
    if debug == 0:
        keep_files_temp = []
        keep_files = [ref_name + '\.rename\.fa', 
                    ref_name + '\.rename\.fa\.pass\.list', 
                    '.*\.scn',
                    ref_name + '\.rename\.fa\.LTRlib\.fa', 
                    'confident_TE\.cons\.fa', 
                    'confident_ltr_cut\.fa',
                    'confident_TE\.cons\.fa\.classified', 
                    'longest_repeats(_\d+)?\.flanked\.fa', 
                    'longest_repeats(_\d+)?\.fa',
                    'confident_tir(_\d+)?\.fa',
                    'confident_helitron(_\d+)?\.fa', 
                    'confident_other(_\d+)?\.fa']

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


