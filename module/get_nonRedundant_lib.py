import argparse
import os
import re
import sys


cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, rename_fasta, multi_process_align_and_get_copies, remove_ltr_from_tir



if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--confident_ltr_cut', metavar='confident_ltr_cut',
                        help='e.g., ')
    parser.add_argument('--confident_tir', metavar='confident_tir',
                        help='e.g., ')
    parser.add_argument('--confident_helitron', metavar='confident_helitron',
                        help='e.g., ')
    parser.add_argument('--confident_other', metavar='confident_other',
                        help='e.g., ')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa')

    args = parser.parse_args()

    threads = int(args.t)
    confident_ltr_cut_path = args.confident_ltr_cut
    confident_tir_path = args.confident_tir
    confident_helitron_path = args.confident_helitron
    confident_other_path = args.confident_other
    tmp_output_dir = args.tmp_output_dir

    tmp_output_dir = os.path.abspath(tmp_output_dir) 

    log = Logger(tmp_output_dir+'/HiTE.log', level='debug')

    # 1. confident_ltr_cut_path比对到TIR候选序列上，并且过滤掉出现在LTR库中的TIR序列
    temp_dir = tmp_output_dir + '/tir_blast_ltr'
    all_copies = multi_process_align_and_get_copies(confident_ltr_cut_path, confident_tir_path, temp_dir, 'tir',
                                                    threads, query_coverage=0.8)
    remove_ltr_from_tir(confident_ltr_cut_path, confident_tir_path, all_copies)

    # 2. 生成一致性tir序列
    confident_tir_rename_path = tmp_output_dir + '/confident_tir.rename.fa'
    rename_fasta(confident_tir_path, confident_tir_rename_path)

    confident_tir_rename_consensus = tmp_output_dir + '/confident_tir.rename.cons.fa'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_tir_rename_path + ' -o ' + confident_tir_rename_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    # 合并所有的TE（TIR+Helitron+Other）
    confident_TE_path = tmp_output_dir + '/confident_TE.fa'
    os.system('cat ' + confident_tir_rename_consensus + ' > ' + confident_TE_path)
    os.system('cat ' + confident_helitron_path + ' >> ' + confident_TE_path)
    os.system('cat ' + confident_other_path + ' >> ' + confident_TE_path)
    os.system('cat ' + confident_ltr_cut_path + ' >> ' + confident_TE_path)

    # 3.generate consensus
    sample_name = 'test'
    confident_TE_consensus = tmp_output_dir + '/confident_TE.cons.fa'

    rename_fasta(confident_TE_path, confident_TE_path)
    contignames, contigs = read_fasta(confident_TE_path)
    new_contigs = {}
    for name in contignames:
        seq = contigs[name]
        if len(seq) < 100:
            continue
        new_contigs[name] = seq
    store_fasta(new_contigs, confident_TE_path)

    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_TE_path + ' -o ' + confident_TE_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)



