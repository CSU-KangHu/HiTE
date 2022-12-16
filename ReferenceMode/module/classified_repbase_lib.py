import argparse
import os
import sys


cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, read_fasta_v1

if __name__ == '__main__':
    tools_dir = os.getcwd() + '/../tools'
    threads = 48

    dir_path = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/repbase'
    repbase_names = ['cbrrep.ref', 'drorep.ref', 'oryrep.ref', 'zebrep.ref']
    #repbase_names = ['cbrrep.ref']

    for name in repbase_names:
        repbase_path = dir_path + '/' + name

        repbase_contignames, repbase_contigs = read_fasta_v1(repbase_path)

        #记录name与后面的注释
        repbase_notes = {}
        for repbase_contigname in repbase_contignames:
            parts = repbase_contigname.split('\t')
            repbase_notes[parts[0]] = parts[1]

        TEClass_home = os.getcwd() + '/../classification'
        TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + name \
                          + ' --consensus ' + repbase_path + ' --genome 1' \
                          + ' --thread_num ' + str(threads) + ' --split_num ' + str(48) + ' -o ' + dir_path
        os.system(TEClass_command)

        classified_path = repbase_path + '.final.classified'
        classified_contignames, classified_contigs = read_fasta(classified_path)
        #追加注释
        noted_contigs = {}
        for classified_contigname in classified_contignames:
            orig_name = classified_contigname.split('#')[0]
            note = repbase_notes[orig_name]
            noted_contigs[classified_contigname + '\t' + note] = classified_contigs[classified_contigname]
        store_fasta(noted_contigs, classified_path)


