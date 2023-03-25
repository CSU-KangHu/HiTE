import argparse
import os
import re
import sys


cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, rename_fasta, multi_process_align_and_get_copies, remove_ltr_from_tir, \
    get_domain_info

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--confident_TE_consensus', metavar='confident_TE_consensus',
                        help='e.g., ')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa')
    parser.add_argument('--classified', metavar='classified',
                        help='e.g., 1')
    parser.add_argument('--domain', metavar='domain',
                        help='e.g., 1')
    parser.add_argument('--TEClass_home', metavar='TEClass_home',
                        help='e.g., ')
    parser.add_argument('--protein_path', metavar='protein_path',
                        help='e.g., ')
    parser.add_argument('--debug', metavar='debug',
                        help='e.g., 1')

    args = parser.parse_args()
    threads = int(args.t)
    confident_TE_consensus = args.confident_TE_consensus
    tmp_output_dir = args.tmp_output_dir
    classified = args.classified
    domain = args.domain
    TEClass_home = args.TEClass_home
    protein_path = args.protein_path
    debug = int(args.debug)

    tmp_output_dir = os.path.abspath(tmp_output_dir) 

    log = Logger(tmp_output_dir+'/HiTE.log', level='debug')

    sample_name = 'test'
    # 1.classify
    if classified is not None and int(classified) == 1:
        # # use RepeatClassifier to classify TE models
        # command = 'cd ' + tmp_output_dir + ' && ' + RepeatModeler_Home + '/RepeatClassifier -pa ' + str(threads) + ' -consensi ' + confident_TE_consensus
        # log.logger.debug(command)
        # os.system(command + ' > /dev/null 2>&1')

        TEClass_command = 'python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
                          + ' --consensus ' + confident_TE_consensus + ' --genome 1' \
                          + ' --thread_num ' + str(threads) + ' --split_num ' + str(threads) + ' -o ' + tmp_output_dir
        #log.logger.debug(TEClass_command)
        os.system(TEClass_command)

        # 把unknown的序列放在最后
        classified_TE_path = confident_TE_consensus + '.classified'
        names, contigs = read_fasta(classified_TE_path)
        names.sort(key=lambda x: x.split('#')[1])
        with open(classified_TE_path, 'w') as f_save:
            for name in names:
                f_save.write('>'+name+'\n'+contigs[name]+'\n')
        f_save.close()

    # 2.get domain
    if domain is not None and int(domain) == 1:
        output_table = confident_TE_consensus + '.domain'
        temp_dir = tmp_output_dir + '/domain'
        get_domain_info(confident_TE_consensus, protein_path, output_table, threads, temp_dir)

