import argparse
import os
import random
import time

from Util import Logger, read_fasta, store_fasta

model_library = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/output/model_lib.fa'
output_library = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/output/krf_output/CRD.2022-05-15.21-16-3/family_model.fasta'
#output_library = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/output/rm2_output/family_model.fasta'
output_dir = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/output/krf_output'

similarity_cutoff = 0.9
length_difference_cutoff = 0.9

if __name__ == '__main__':
    log = Logger('GenomeSimulator.log', level='debug')
    model_contignames, model_contigs = read_fasta(model_library)
    output_contignames, output_contigs = read_fasta(output_library)
    threads = 48
    tool_dir = '/public/home/hpc194701009/CompleteRepeatDetection/ReferenceMode'
    blast_program_dir = tool_dir + '/tools/rmblast-2.9.0-p2'
    blastn2Results_path = output_dir + '/tmpBlastResults2.out'
    makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + model_library
    align_command = blast_program_dir + '/bin/blastn -db ' + model_library + ' -num_threads ' \
                    + str(threads) + ' -query ' + output_library + ' -outfmt 6 > ' + blastn2Results_path
    log.logger.debug(makedb_command)
    os.system(makedb_command)
    log.logger.debug(align_command)
    os.system(align_command)

    query_name_set = set()
    target_name_set = set()
    with open(blastn2Results_path, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            query_name = parts[0]
            target_name = parts[1]
            identity = float(parts[2])
            match_base = int(parts[3])
            q_start = int(parts[6])
            q_end = int(parts[7])
            t_start = int(parts[8])
            t_end = int(parts[9])

            query_len = len(output_contigs[query_name])
            target_len = len(model_contigs[target_name])
            long_len = query_len if query_len > target_len else target_len
            short_len = query_len if query_len < target_len else target_len

            similarity = float(match_base) / short_len

            len_diff = abs(query_len - target_len)
            length_difference = float(long_len-len_diff)/long_len
            if similarity >= similarity_cutoff and length_difference >= length_difference_cutoff:
                query_name_set.add(query_name)
                target_name_set.add(target_name)
            #print((similarity, length_difference))

    precision = float(len(query_name_set))/len(output_contigs)
    recall = float(len(target_name_set)) / len(model_contigs)
    f1_score = 2 * precision * recall / (precision + recall)

    print('true repeats: %d, total find repeats: %d, precision: %f' % (len(query_name_set), len(output_contigs), precision))
    print('recall: %f' % recall)
    print('f1_score: %f' % f1_score)
    #print(target_name_set)
