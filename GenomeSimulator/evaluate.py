import argparse
import os
import random
import time

from Util import Logger, read_fasta, store_fasta

model_library = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/no_simple_repeats/dmel_curated.fasta'
output_library = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-26.0-39-21/family_dmel.fasta'
output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-26.0-39-21'

# model_library = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/10M_low_freq_out/model_lib.fa'
# output_library = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/10M_low_freq_out/krf_output/CRD.2022-05-26.0-13-4/family_model.fasta'
# output_dir = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/10M_low_freq_out/krf_output/CRD.2022-05-26.0-13-4'

#output_library = '/public/home/hpc194701009/KmerRepFinder_test/library/rm2_lib/sort_lib/Dmel-families.fa'
#output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/rm2_lib/sort_lib'

similarity_cutoff = 0.95
length_difference_cutoff = 0.95

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
    query_cluster = {}
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
            key = query_name + '$' +target_name
            if not query_cluster.__contains__(key):
                query_cluster[key] = ([], -1, -1)
            tuple = query_cluster[key]
            cluster = tuple[0]
            if identity >= similarity_cutoff:
                cluster.append((q_start, q_end, t_start, t_end ))
            query_cluster[key] = (cluster, query_len, target_len)
            # long_len = query_len if query_len > target_len else target_len
            # short_len = query_len if query_len < target_len else target_len
            #
            # similarity = float(match_base) / short_len
            #
            # len_diff = abs(query_len - target_len)
            # length_difference = float(long_len-len_diff)/long_len
            # if similarity >= similarity_cutoff and length_difference >= length_difference_cutoff:
            #     query_name_set.add(query_name)
            #     target_name_set.add(target_name)
            # #print((similarity, length_difference))
    for key in query_cluster.keys():
        parts = key.split('$')
        query_name = parts[0]
        target_name = parts[1]

        tuple = query_cluster[key]
        query_len = tuple[1]
        target_len = tuple[2]
        query_array = ['' for _ in range(query_len)]
        target_array = ['' for _ in range(target_len)]
        query_masked_len = 0
        target_masked_len = 0
        for record in tuple[0]:
            qstart = record[0]
            qend = record[1]
            if qstart > qend:
                tmp = qend
                qend = qstart
                qstart = tmp
            for i in range(qstart, qend):
                query_array[i] = 'X'

            tstart = record[2]
            tend = record[3]
            if tstart > tend:
                tmp = tend
                tend = tstart
                tstart = tmp
            for i in range(tstart, tend):
                target_array[i] = 'X'
        for j in range(len(query_array)):
            if query_array[j] == 'X':
                query_masked_len += 1
        for j in range(len(target_array)):
            if target_array[j] == 'X':
                target_masked_len += 1
        # long_len = query_masked_len if query_masked_len > target_masked_len else target_masked_len
        # short_len = query_masked_len if query_masked_len < target_masked_len else target_masked_len
        # len_diff = abs(query_len - target_len)
        # length_difference = float(long_len-len_diff)/long_len
        # if key.__contains__('family-273#LTR/Gypsy'):
        #     print(query_masked_len, query_len, target_masked_len, target_len)
        #     print(float(query_masked_len)/query_len, float(target_masked_len)/target_len)
        if float(query_masked_len)/query_len >= length_difference_cutoff and float(target_masked_len)/target_len >= length_difference_cutoff:
            query_name_set.add(query_name)
            target_name_set.add(target_name)


    precision = float(len(query_name_set))/len(output_contigs)
    recall = float(len(target_name_set)) / len(model_contigs)
    f1_score = 2 * precision * recall / (precision + recall)

    print('true repeats: %d, total find repeats: %d, precision: %f' % (len(query_name_set), len(output_contigs), precision))
    print('recall: %f' % recall)
    print('f1_score: %f' % f1_score)
    #print(set(output_contignames)-query_name_set)
    #print(set(model_contignames) - target_name_set)

    #print(target_name_set)
