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
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--classified_TE_consensus', metavar='classified_TE_consensus',
                        help='e.g., ')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa')
    parser.add_argument('--annotate', metavar='annotate',
                        help='e.g., 1')
    parser.add_argument('-r', metavar='reference',
                        help='e.g., Reference Path')

    args = parser.parse_args()
    threads = int(args.t)
    classified_TE_consensus = args.classified_TE_consensus
    tmp_output_dir = args.tmp_output_dir
    annotate = args.annotate
    reference = args.r

    reference = os.path.abspath(reference)

    log = Logger(tmp_output_dir+'/HiTE_annotate_genome.log', level='debug')

    # 3.annotate the genome
    if annotate is not None and int(annotate) == 1:
        RepeatMasker_command = 'cd ' + tmp_output_dir + ' && RepeatMasker -e ncbi -pa ' + str(threads) \
                               + ' -q -no_is -norna -nolow -div 40 -gff -lib ' + classified_TE_consensus + ' -cutoff 225 ' \
                               + reference
        log.logger.debug(RepeatMasker_command)
        os.system(RepeatMasker_command)

        mv_file_command = 'mv ' + reference + '.out ' + tmp_output_dir + '/HiTE.out && mv ' \
                          + reference + '.tbl ' + tmp_output_dir + '/HiTE.tbl && mv ' \
                          + reference + '.out.gff ' + tmp_output_dir + '/HiTE.gff'
        log.logger.debug(mv_file_command)
        os.system(mv_file_command)


