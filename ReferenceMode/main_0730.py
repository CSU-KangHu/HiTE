import argparse
import codecs
import random
import sys

import datetime
import json
import multiprocessing
import os
import threading
import time
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import math

#from command import run_bwa
from Util import convertToUpperCase, read_fasta, getReverseSequence, \
    Logger, store_fasta, get_candidate_repeats, divided_array, split2cluster_normal, convertToUpperCase_v1, multi_line

dict_encode = {'A': 0b001, 'C': 0b011, 'G': 0b010, 'T': 0b100, 'N': 0b101}
def three_bit_encode(str):
    bytes = 0b000
    for i in range(len(str)):
        bytes <<= 3
        base = str[i]
        if not dict_encode.__contains__(base):
            base = 'N'
        bytes += dict_encode[base]
    return bytes

def print_seqs(header, sequence, length, outfile):
    print('>' + header, file=outfile)
    while len(sequence) > 0:
        print(sequence[:length], file=outfile)
        sequence = sequence[length:]

def generate_candidate_repeats(cur_segments, k_num, unique_kmer_map, fault_tolerant_bases):
    cur_lines = []
    cur_masked_segments = []
    for line in cur_segments:
        cur_lines.append(line)
        masked_line = list(line)
        last_masked_pos = -1
        for i in range(len(line)-k_num+1):
            kmer = line[i: i+k_num]
            # get reverse complement kmer
            r_kmer = getReverseSequence(kmer)
            # filter invalid kmer, contains 'N'
            if "N" in r_kmer:
                continue
            unique_key = kmer if kmer < r_kmer else r_kmer

            if unique_kmer_map.__contains__(unique_key):
                # mask position
                if last_masked_pos == -1:
                    for j in range(i, i+k_num):
                        masked_line[j] = 'X'
                    last_masked_pos = i+k_num-1
                else:
                    # do not need to mask position which has been masked
                    start_mask_pos = i if i > last_masked_pos else last_masked_pos+1
                    end_mask_pos = i+k_num
                    for j in range(start_mask_pos, end_mask_pos):
                        masked_line[j] = 'X'
                    last_masked_pos = end_mask_pos - 1
        cur_masked_segments.append(masked_line)

    repeat_list = []
    cur_repeat_str = ''
    try_connect_str = ''
    last_end_pos = -1
    for seq_index, cur_masked_segment in enumerate(cur_masked_segments):
        for i in range(len(cur_masked_segment)):
            if cur_masked_segment[i] == 'X':
                if try_connect_str != '':
                    cur_repeat_str = try_connect_str
                    try_connect_str = ''
                cur_repeat_str = cur_repeat_str + cur_lines[seq_index][i]
                last_end_pos = i
            elif cur_repeat_str != '' and cur_masked_segment[i] != 'N':
                if (i - last_end_pos) <= fault_tolerant_bases:
                    if try_connect_str == '':
                        try_connect_str = cur_repeat_str
                    try_connect_str = try_connect_str + cur_lines[seq_index][i]
                else:
                    repeat_list.append(cur_repeat_str)
                    cur_repeat_str = ''
                    try_connect_str = ''
    return repeat_list

def generate_candidate_repeats_v1(contigs, k_num, unique_kmer_map, fault_tolerant_bases):
    cur_lines = []
    cur_masked_segments = {}
    for ref_name in contigs.keys():
        line = contigs[ref_name]
        cur_lines.append(line)
        masked_line = list(line)
        last_masked_pos = -1
        for i in range(len(line)-k_num+1):
            kmer = line[i: i+k_num]
            # get reverse complement kmer
            r_kmer = getReverseSequence(kmer)
            # filter invalid kmer, contains 'N'
            if "N" in r_kmer:
                continue
            unique_key = kmer if kmer < r_kmer else r_kmer

            if unique_kmer_map.__contains__(unique_key):
                # mask position
                if last_masked_pos == -1:
                    for j in range(i, i+k_num):
                        masked_line[j] = 'X'
                    last_masked_pos = i+k_num-1
                else:
                    # do not need to mask position which has been masked
                    start_mask_pos = i if i > last_masked_pos else last_masked_pos+1
                    end_mask_pos = i+k_num
                    for j in range(start_mask_pos, end_mask_pos):
                        masked_line[j] = 'X'
                    last_masked_pos = end_mask_pos - 1
        cur_masked_segments[ref_name] = masked_line

    repeat_dict = {}
    cur_repeat_str = ''
    try_connect_str = ''
    last_start_pos = -1
    last_end_pos = -1
    for seq_index, cur_masked_item in enumerate(cur_masked_segments.items()):
        ref_name = cur_masked_item[0]
        cur_masked_segment = cur_masked_item[1]
        if not repeat_dict.__contains__(ref_name):
            repeat_dict[ref_name] = []
        repeat_list = repeat_dict[ref_name]
        for i in range(len(cur_masked_segment)):
            if cur_masked_segment[i] == 'X':
                if last_start_pos == -1:
                    # record masked sequence start position
                    last_start_pos = i
                if try_connect_str != '':
                    # recover skip gap sequence
                    cur_repeat_str = try_connect_str
                    try_connect_str = ''
                cur_repeat_str = cur_repeat_str + cur_lines[seq_index][i]
                last_end_pos = i
            elif cur_repeat_str != '' and cur_masked_segment[i] != 'N':
                # meet unmasked base
                if (i - last_end_pos) <= fault_tolerant_bases:
                    # skip gap
                    if try_connect_str == '':
                        try_connect_str = cur_repeat_str
                    try_connect_str = try_connect_str + cur_lines[seq_index][i]
                else:
                    # can not skip gap
                    repeat_list.append((last_start_pos, last_end_pos, cur_repeat_str))
                    cur_repeat_str = ''
                    try_connect_str = ''
                    last_start_pos = -1
        # keep last masked sequence
        if cur_repeat_str != '':
            repeat_list.append((last_start_pos, last_end_pos, cur_repeat_str))
            cur_repeat_str = ''
            try_connect_str = ''
            last_start_pos = -1
        repeat_dict[ref_name] = repeat_list
    return repeat_dict, cur_masked_segments


def generate_kmer_map(cur_segments, k_num, unique_kmer_map, merged_repeats, partiton_index):
    log.logger.debug('partition %d process: %d segments' %(partiton_index, len(cur_segments)))
    kmer_info_map = {}
    cur_lines = []
    cur_masked_segments = []
    for line in cur_segments:
        parts = line.split("\t")
        contigName = parts[0].split(" ")[0]
        # remove ">"
        contigName = contigName[1:]
        start = int(parts[1])
        line = parts[2]
        cur_lines.append(line)
        masked_line = list(line)
        last_masked_pos = -1
        for i in range(len(line)-k_num+1):
            kmer = line[i: i+k_num]
            # get reverse complement kmer
            r_kmer = getReverseSequence(kmer)
            # filter invalid kmer, contains 'N'
            if "N" in r_kmer:
                continue
            # encode with binary
            kmer_num = three_bit_encode(kmer)
            r_kmer_num = three_bit_encode(r_kmer)
            unique_key = kmer_num if kmer_num < r_kmer_num else r_kmer_num

            pos = start + i
            if unique_kmer_map.__contains__(unique_key):
                # mask position
                if last_masked_pos == -1:
                    for j in range(i, i+k_num):
                        masked_line[j] = 'X'
                    last_masked_pos = i+k_num-1
                else:
                    # do not need to mask position which has been masked
                    start_mask_pos = i if i > last_masked_pos else last_masked_pos+1
                    end_mask_pos = i+k_num
                    for j in range(start_mask_pos, end_mask_pos):
                        masked_line[j] = 'X'
                    last_masked_pos = end_mask_pos - 1
        cur_masked_segments.append(masked_line)

    repeat_list = []
    cur_repeat_str = ''
    for seq_index, cur_masked_segment in enumerate(cur_masked_segments):
        for i in range(len(cur_masked_segment)):
            if cur_masked_segment[i] == 'X':
                cur_repeat_str = cur_repeat_str + cur_lines[seq_index][i]
            elif cur_repeat_str != '':
                repeat_list.append(cur_repeat_str)
                cur_repeat_str = ''
    merged_repeats[partiton_index] = repeat_list
    log.logger.debug('partition %d process finished' % (partiton_index))

def getUniqueKmer_v1(cur_segments, partiton_index):
    #print('partition %d process: %d segments' % (partiton_index, len(cur_segments)))
    unique_kmers = []
    for line in cur_segments:
        kmer = line.split(' ')[0]
        r_kmer = getReverseSequence(kmer)
        #unique_key = kmer if kmer < r_kmer else r_kmer
        unique_kmers.append(kmer)
        unique_kmers.append(r_kmer)
    return unique_kmers

def getUniqueKmer(unique_kmer_path):
    unique_kmer_map = {}
    with open(unique_kmer_path, 'r') as f_r:
        for line in f_r:
            kmer = line.split(' ')[0]
            r_kmer = getReverseSequence(kmer)
            kmer_num = three_bit_encode(kmer)
            r_kmer_num = three_bit_encode(r_kmer)
            unique_key = kmer_num if kmer_num < r_kmer_num else r_kmer_num
            unique_kmer_map[unique_key] = 1
    return unique_kmer_map


# A sequence may include multiple align position, e.g.,
# Node_0-len_5109 Node_0-len_5109 100.000 4651    0       0       459     5109    1       4651    0.0     8589
# Node_0-len_5109 Node_30444-len_20481    100.000 217     0       0       1       217     20265   20481   1.37e-110       401

def parse_self_blast_output(blastnResults_path, repeats_path, candidate_repeats_path):
    repeatContigNames, repeatContigs = read_fasta(repeats_path)
    # To facilite searching
    # strcuture: {'Node_1': {'Node_1': [(),(),], 'Node_2':[(),(),], ...}}
    # step1. construct blast records clustering by query name
    query_records = {}
    with open(blastnResults_path, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            query_name = parts[0]
            target_name = parts[1]
            identity = float(parts[2])
            match_base = int(parts[3])
            query_length = int(query_name.split('-')[1].split('_')[1])
            q_start = int(parts[6])
            q_end = int(parts[7])
            t_start = int(parts[8])
            t_end = int(parts[9])

            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            records = query_records[query_name]
            if not records.__contains__(target_name):
                records[target_name] = []
            same_target_records = records[target_name]
            same_target_records.append((identity, match_base, query_length, q_start, q_end, t_start, t_end))
            records[target_name] = same_target_records
            query_records[query_name] = records

    # strcuture: {'Node_1': {'Node_1': [(),(),], 'Node_2':[(),(),], ...}}
    # Node_0-len_5109 Node_0-len_5109 100.000 4651    0       0       459     5109    1       4651    0.0     8589
    # Node_0-len_5109 Node_30444-len_20481    100.000 217     0       0       1       217     20265   20481   1.37e-110       401

    # step2. splice sequence
    # a map is used to avoid adding redudant sequence
    candidate_family_repeat = []
    stored_repeat = {}
    for query_name in query_records.keys():
        query_length = int(query_name.split('-')[1].split('_')[1])
        records = query_records[query_name]
        for target_name in records.keys():
            for record in records[target_name]:
                # self alignment, identity < 80% and segment length < 80bp should be neglected
                if (query_name == target_name and record[1] == query_length) \
                        or record[0] < 80:
                    continue
                # jugde direction
                q_start = record[3]
                q_end = record[4]
                if q_start > q_end:
                    q_tmp = q_start
                    q_start = q_end
                    q_end = q_tmp
                seg_seq1 = repeatContigs[query_name][q_start-1: q_end]

                t_start = record[5]
                t_end = record[6]
                if t_start > t_end:
                    t_tmp = t_start
                    t_start = t_end
                    t_end = t_tmp
                seg_seq2 = repeatContigs[target_name][t_start-1: t_end]

                if not stored_repeat.__contains__(query_name):
                    stored_repeat[query_name] = []
                query_pos_records = stored_repeat[query_name]
                queryIsNewSequence = True
                for pos_record in query_pos_records:
                    if q_start >= pos_record[0] and q_end <= pos_record[1]:
                        queryIsNewSequence = False
                        break
                if queryIsNewSequence:
                    candidate_family_repeat.append(seg_seq1)
                    query_pos_records.append((q_start, q_end))
                stored_repeat[query_name] = query_pos_records

                if not stored_repeat.__contains__(target_name):
                    stored_repeat[target_name] = []
                target_pos_records = stored_repeat[target_name]
                targetIsNewSequence = True
                for pos_record in target_pos_records:
                    if t_start >= pos_record[0] and t_end <= pos_record[1]:
                        targetIsNewSequence = False
                        break
                if targetIsNewSequence:
                    # if spliced segments is not exactly equal, add to set
                    if not math.isclose(record[0], 100.000, rel_tol=1e-5):
                        candidate_family_repeat.append(seg_seq2)
                    target_pos_records.append((t_start, t_end))
                stored_repeat[target_name] = target_pos_records


    # step3. generate candidate repeats
    node_index = 0
    with open(candidate_repeats_path, 'w') as f_save:
        for sequence in candidate_family_repeat:
            f_save.write('>Node_'+str(node_index)+'-len_'+str(len(sequence))+'\n'+sequence+'\n')
            node_index += 1


def run_LTR_retriever_v1(Genome_Tools_Home, LTR_retriever_Home, reference, tmp_output_dir, threads):
    starttime = time.time()
    log.logger.debug('start LTR_harvest detection...')
    ltrharvest_command1 = Genome_Tools_Home + '/bin/gt suffixerator -db ' + reference + ' -indexname ' \
                          + reference + ' -tis -suf -lcp -des -ssp -sds -dna'
    ltrharvest_command2 = Genome_Tools_Home + '/bin/gt ltrharvest -index ' + reference \
                          + ' -seed 20 -minlenltr 100 -maxlenltr 7000 -similar 85 -mintsd 4 -maxtsd 6 ' \
                            '-motif TGCA -motifmis 1 -vic 10 -seqids yes > ' + tmp_output_dir + '/genome.fa.harvest.scn'

    os.system(ltrharvest_command1)
    os.system(ltrharvest_command2)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("LTR_harvest running time: %.8s s" % (dtime))

    starttime = time.time()
    log.logger.debug('start LTR_retriever detection...')
    LTR_retriever_command = 'cd ' + tmp_output_dir + ' && ' + LTR_retriever_Home + '/LTR_retriever -noanno -genome ' + reference \
                            + ' -inharvest ' + tmp_output_dir + '/genome.fa.harvest.scn -threads ' + str(threads)
    os.system(LTR_retriever_command)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("LTR_retriever running time: %.8s s" % (dtime))


def run_LTR_retriever(LTR_retriever_Home, reference, tmp_output_dir, threads):
    starttime = time.time()
    log.logger.debug('start LTR_retriever detection...')
    LTR_retriever_command = 'cd ' + tmp_output_dir + ' && ' + LTR_retriever_Home + '/LTR_retriever -genome ' + reference \
                            + ' -inharvest ' + tmp_output_dir + '/genome.fa.harvest.scn -threads ' + str(threads)
    os.system(LTR_retriever_command)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("LTR_retriever running time: %.8s s" % (dtime))

def run_LTR_harvest(Genome_Tools_Home, reference, tmp_output_dir):
    starttime = time.time()
    log.logger.debug('start LTR_harvest detection...')
    ltrharvest_command1 = Genome_Tools_Home + '/bin/gt suffixerator -db ' + reference + ' -indexname ' \
                          + reference + ' -tis -suf -lcp -des -ssp -sds -dna'
    ltrharvest_command2 = Genome_Tools_Home + '/bin/gt ltrharvest -index ' + reference \
                          + ' -seed 20 -minlenltr 100 -maxlenltr 7000 -similar 85 -mintsd 4 -maxtsd 6 ' \
                            '-motif TGCA -motifmis 1 -vic 10 -seqids yes > ' + tmp_output_dir + '/genome.fa.harvest.scn'

    os.system(ltrharvest_command1)
    os.system(ltrharvest_command2)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("LTR_harvest running time: %.8s s" % (dtime))



def run_GRF(GRF_Home, reference, tmp_output_dir, threads):
    grf_tir_command = GRF_Home + '/bin/grf-main -i ' + reference + ' -o ' + tmp_output_dir + ' -c 0 --min_tr 10 -t ' + str(threads)
    os.system(grf_tir_command)

    grf_mite_command = 'sh ' + GRF_Home + '/script/run_mite_detection.sh ' + reference + ' ' + tmp_output_dir + ' ' + str(threads)
    os.system(grf_mite_command)

def get_flanking_instances(repeats_path, blast_program_dir, tools_dir, extend_base_threshold, flanking_len):
    split_longest_repeats_path = repeats_path[0]
    reference_path = repeats_path[1]
    blastnResults_path = repeats_path[2]

    ref_names, ref_contigs = read_fasta(reference_path)

    makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + reference_path
    align_command = blast_program_dir + '/bin/blastn -db ' + reference_path + ' -num_threads ' \
                    + str(1) + ' -query ' + split_longest_repeats_path + ' -outfmt 6 > ' + blastnResults_path
    os.system(makedb_command)
    os.system(align_command)

    # query_records = {query_name: {subject_name: [(q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end)] }}
    query_records = {}
    with open(blastnResults_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            identity = float(parts[2])
            alignment_len = int(parts[3])
            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            if identity < 80:
                continue
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))

    instances = {}
    for query_name in query_records.keys():
        query_len = int(query_name.split('-len_')[1])
        subject_dict = query_records[query_name]
        repeat_list = []
        for subject_name in subject_dict.keys():
            ref_seq = ref_contigs[subject_name]
            #masked_ref_seq = list(ref_seq)
            subject_pos = subject_dict[subject_name]
            # alignment direction
            forward_subject_pos = []
            reverse_subject_pos = []
            for pos_item in subject_pos:
                s_start = pos_item[2]
                s_end = pos_item[3]
                if s_start <= s_end:
                    forward_subject_pos.append(pos_item)
                else:
                    reverse_subject_pos.append(pos_item)

            # connect repeats
            connected_repeat_list = []
            # start forward process
            forward_subject_pos.sort(key=lambda x: (x[2], x[3]))
            last_start_pos = -1
            last_end_pos = -1
            for repeat_item in forward_subject_pos:
                start_pos = repeat_item[2]
                end_pos = repeat_item[3]
                if last_start_pos != -1:
                    if (start_pos - last_end_pos) <= extend_base_threshold:
                        # connected repeat
                        last_end_pos = end_pos
                    else:
                        # not connected repeat
                        # keep last connected repeat
                        connected_repeat_list.append((last_start_pos, last_end_pos))
                        last_start_pos = -1
                        last_end_pos = -1
                if last_start_pos == -1:
                    # start a new connected repeat
                    last_start_pos = start_pos
                    last_end_pos = end_pos
            if last_start_pos != -1:
                connected_repeat_list.append((last_start_pos, last_end_pos))

            # start reverse process
            reverse_subject_pos.sort(key=lambda x: (-x[2], -x[3]))
            last_start_pos = -1
            last_end_pos = -1
            for repeat_item in reverse_subject_pos:
                start_pos = repeat_item[2]
                end_pos = repeat_item[3]
                if last_start_pos != -1:
                    if (last_end_pos - start_pos) <= extend_base_threshold:
                        # connected repeat
                        last_end_pos = end_pos
                    else:
                        # not connected repeat
                        # keep last connected repeat
                        connected_repeat_list.append((last_start_pos, last_end_pos))
                        last_start_pos = -1
                        last_end_pos = -1
                if last_start_pos == -1:
                    # start a new connected repeat
                    last_start_pos = start_pos
                    last_end_pos = end_pos
            if last_start_pos != -1:
                connected_repeat_list.append((last_start_pos, last_end_pos))

            for pos_item in connected_repeat_list:
                start_pos = pos_item[0]
                end_pos = pos_item[1]
                cur_instance_len = abs(end_pos - start_pos) + 1
                long_len = cur_instance_len if cur_instance_len > query_len else query_len
                short_len = cur_instance_len if cur_instance_len <= query_len else query_len
                if float(short_len) / long_len >= 0.9:
                    # extend flanking sequence
                    if start_pos <= end_pos:
                        extend_start = start_pos - flanking_len
                        extend_end = end_pos + flanking_len
                        # if cannot extend, ignore this instance
                        if extend_start < 0 or extend_end > len(ref_seq) - 1:
                            continue
                        extend_seq = ref_seq[extend_start-1: extend_end]
                        if not extend_seq.__contains__('N'):
                            repeat_list.append(extend_seq)
                    else:
                        extend_start = start_pos + flanking_len
                        extend_end = end_pos - flanking_len
                        # if cannot extend, ignore this instance
                        if extend_end < 0 or extend_start > len(ref_seq) - 1:
                            continue
                        extend_seq = ref_seq[extend_start-1: extend_end: -1]
                        if not extend_seq.__contains__('N'):
                            repeat_list.append(extend_seq)
        instances[query_name] = repeat_list
    return instances

# def filter_false_positive(repeats_path, blast_program_dir, tools_dir, extend_base_threshold, flanking_len):
#     instance_path = repeats_path[0]
#     reference_path = repeats_path[1]
#     blastnResults_path = repeats_path[2]
#
#     align_command = blast_program_dir + '/bin/blastn -db ' + reference_path + ' -num_threads ' \
#                     + str(1) + ' -query ' + instance_path + ' -outfmt 6 > ' + blastnResults_path
#     os.system(align_command)
#
#     tp_names = []
#     instance_names, instance_contigs = read_fasta(instance_path)
#     query_records = {}
#     with open(blastnResults_path, 'r') as f_r:
#         for idx, line in enumerate(f_r):
#             parts = line.split('\t')
#             query_name = parts[0]
#             subject_name = parts[1]
#             identity = float(parts[2])
#             alignment_len = int(parts[3])
#             q_start = int(parts[6])
#             q_end = int(parts[7])
#             s_start = int(parts[8])
#             s_end = int(parts[9])
#             if identity < 80:
#                 continue
#             if not query_records.__contains__(query_name):
#                 query_records[query_name] = {}
#             subject_dict = query_records[query_name]
#
#             if not subject_dict.__contains__(subject_name):
#                 subject_dict[subject_name] = []
#             subject_pos = subject_dict[subject_name]
#             subject_pos.append((q_start, q_end, s_start, s_end))
#
#     max_over_distance = 50
#     for query_name in query_records.keys():
#         subject_dict = query_records[query_name]
#         query_seq = instance_contigs[query_name]
#         query_real_start = flanking_len
#         query_real_end = len(query_seq)-flanking_len
#         over_count = 0
#         is_false_positive = False
#         for subject_name in subject_dict.keys():
#             if is_false_positive:
#                 break
#             subject_pos = subject_dict[subject_name]
#             for pos_item in subject_pos:
#                 if over_count >= 2:
#                     is_false_positive = True
#                     break
#                 q_start = pos_item[0]
#                 q_end = pos_item[1]
#                 if query_real_start - q_start > max_over_distance \
#                         or q_end - query_real_end > max_over_distance:
#                     over_count += 1
#         if not is_false_positive:
#             tp_names.append(query_name)
#
#
#         # masked_query_seq = [0 for i in range(len(query_seq))]
#         # for subject_name in subject_dict.keys():
#         #     subject_pos = subject_dict[subject_name]
#         #     for pos_item in subject_pos:
#         #         q_start = pos_item[0]
#         #         q_end = pos_item[1]
#         #         for j in range(q_start-1, q_end):
#         #             masked_query_seq[j] = masked_query_seq[j] + 1
#         #
#         # last_start_pos = -1
#         # last_end_pos = -1
#         # cur_repeat_str = ''
#         # try_connect_str = ''
#         # repeat_list = []
#         # for i in range(len(masked_query_seq)):
#         #     if masked_query_seq[i] >= 2:
#         #         if last_start_pos == -1:
#         #             # record masked sequence start position
#         #             last_start_pos = i
#         #         if try_connect_str != '':
#         #             # recover skip gap sequence
#         #             cur_repeat_str = try_connect_str
#         #             try_connect_str = ''
#         #         cur_repeat_str = cur_repeat_str + query_seq[i]
#         #         last_end_pos = i
#         #     elif cur_repeat_str != '':
#         #         # meet unmasked base
#         #         if (i - last_end_pos) <= extend_base_threshold:
#         #             # skip gap
#         #             if try_connect_str == '':
#         #                 try_connect_str = cur_repeat_str
#         #             try_connect_str = try_connect_str + query_seq[i]
#         #         else:
#         #             # can not skip gap
#         #             repeat_list.append((last_start_pos, last_end_pos, cur_repeat_str))
#         #             cur_repeat_str = ''
#         #             try_connect_str = ''
#         #             last_start_pos = -1
#         # # keep last masked sequence
#         # if cur_repeat_str != '':
#         #     repeat_list.append((last_start_pos, last_end_pos, cur_repeat_str))
#         #
#         # max_over_distance = 50
#         # if len(repeat_list) == 1:
#         #     repeat_instance = repeat_list[0]
#         #     if abs(repeat_instance[0]-query_real_start) > max_over_distance \
#         #             or abs(repeat_instance[1]-query_real_end) > max_over_distance:
#         #         continue
#         #     tp_names.append(query_name)
#
#     true_positive_path = single_tmp_dir + '/true_positive.fa'
#     with open(true_positive_path, 'w') as f_save:
#         for name in tp_names:
#             seq = instance_contigs[name]
#             f_save.write('>'+name+'\n'+seq+'\n')


def remove_chimerism(repeats_path, blast_program_dir, round_num, partition_index):
    pure_repeat_file = repeats_path[0]
    repeat_pool_file = repeats_path[1]
    blastn2Results_path = repeats_path[2]

    repeat_pool_names, repeat_pool_contigs = read_fasta(repeat_pool_file)

    if len(repeat_pool_contigs) > 0:
        makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + repeat_pool_file
        align_command = blast_program_dir + '/bin/blastn -db ' + repeat_pool_file + ' -num_threads ' \
                        + str(1) + ' -query ' + pure_repeat_file + ' -outfmt 6 > ' + blastn2Results_path
        os.system(makedb_command)
        os.system(align_command)

    # parse blastn output, determine the repeat boundary
    subject_records = {}
    redundant_repeat = {}
    with open(blastn2Results_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            # print('current line idx: %d' % (idx))
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            identity = float(parts[2])
            alignment_len = int(parts[3])
            query_len = int(query_name.split('-len_')[1])
            subject_len = int(subject_name.split('-len_')[1])
            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            if identity < 80 or float(alignment_len)/query_len < 0.95:
                continue
            if float(alignment_len) / subject_len >= 0.95:
                redundant_repeat[subject_name] = 1
                continue
            # alignment_len/query_len >= 0.95 and alignment_len / subject_len < 0.95
            # this is chimerism record
            if not subject_records.__contains__(subject_name):
                subject_records[subject_name] = {}
            query_dict = subject_records[subject_name]

            if not query_dict.__contains__(query_name):
                query_dict[query_name] = []
            subject_pos = query_dict[query_name]
            subject_pos.append((q_start, q_end, s_start, s_end))

    remove_chimerism_repeats = []
    for subject_name in subject_records.keys():
        query_dict = subject_records[subject_name]
        subject_seq = repeat_pool_contigs[subject_name]
        mask_subject_seq = list(subject_seq)
        redundant_repeat[subject_name] = 1
        for query_name in query_dict.keys():
            subject_pos = query_dict[query_name]
            # masking subject sequence
            for pos_item in subject_pos:
                s_start = pos_item[2]
                s_end = pos_item[3]
                for j in range(s_start-1, s_end):
                    mask_subject_seq[j] = 'X'
        remain_seq = ''
        for j in range(len(mask_subject_seq)):
            if mask_subject_seq[j] != 'X':
                remain_seq += mask_subject_seq[j]
        if len(remain_seq) >= 80:
            remove_chimerism_repeats.append(remain_seq)

    #print('round_num: %d, partition_index: %d, remove_chimerism_repeats size: %d' %(round_num, partition_index, len(remove_chimerism_repeats)))

    cur_repeat_seqs = []
    for name in repeat_pool_contigs.keys():
        if redundant_repeat.__contains__(name):
            continue
        cur_repeat_seqs.append(repeat_pool_contigs[name])
    for seq in remove_chimerism_repeats:
        cur_repeat_seqs.append(seq)

    return cur_repeat_seqs


def multiple_alignment(repeats_path, blast_program_dir, tools_dir):
    split_repeats_path = repeats_path[0]
    original_repeats_path = repeats_path[1]
    blastn2Results_path = repeats_path[2]

    makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + original_repeats_path
    align_command = blast_program_dir + '/bin/blastn -db ' + original_repeats_path + ' -num_threads ' \
                    + str(1) + ' -query ' + split_repeats_path + ' -outfmt 6 > ' + blastn2Results_path
    os.system(makedb_command)
    os.system(align_command)

    return blastn2Results_path

def multiple_alignment_blastx(repeats_path, blast_program_dir, tools_dir):
    split_repeats_path = repeats_path[0]
    protein_db_path = repeats_path[1]
    blastx2Results_path = repeats_path[2]

    blast_db_command = blast_program_dir + '/bin/makeblastdb -dbtype prot -in ' + protein_db_path
    align_command = blast_program_dir + '/bin/blastx -db ' + protein_db_path + ' -num_threads ' \
                    + str(1) + ' -query ' + split_repeats_path + ' -word_size 7 -outfmt 6 > ' + blastx2Results_path
    os.system(blast_db_command)
    os.system(align_command)

    return blastx2Results_path

def get_longest_repeats(repeats_path, blast_program_dir, tools_dir, extend_base_threshold, max_single_repeat_len):
    split_repeats_path = repeats_path[0]
    original_repeats_path = repeats_path[1]
    blastn2Results_path = repeats_path[2]
    split_repeats_names, split_repeats_contigs = read_fasta(split_repeats_path)
    (single_tmp_dir, split_repeats_filename) = os.path.split(split_repeats_path)

    makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + original_repeats_path
    align_command = blast_program_dir + '/bin/blastn -db ' + original_repeats_path + ' -num_threads ' \
                    + str(1) + ' -query ' + split_repeats_path + ' -outfmt 6 > ' + blastn2Results_path
    os.system(makedb_command)
    os.system(align_command)

    # parse blastn output, determine the repeat boundary
    # query_records = {query_name: {subject_name: [(q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end)] }}
    query_records = {}
    with open(blastn2Results_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            #print('current line idx: %d' % (idx))
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            identity = float(parts[2])
            alignment_len = int(parts[3])
            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            if query_name == subject_name or identity < 80:
                continue
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))

    longest_repeats = {}
    for idx, query_name in enumerate(query_records.keys()):
        print('total query size: %d, current query name: %s, idx: %d' % (len(query_records), query_name, idx))
        subject_dict = query_records[query_name]

        # if there are more than one longest query overlap with the final longest query over 90%,
        # then it probably the true TE
        longest_queries = []
        for subject_name in subject_dict.keys():
            subject_pos = subject_dict[subject_name]
            # subject_pos.sort(key=lambda x: (x[2], x[3]))

            # cluster all closed fragments, split forward and reverse records
            forward_pos = []
            reverse_pos = []
            for pos_item in subject_pos:
                if pos_item[2] > pos_item[3]:
                    reverse_pos.append(pos_item)
                else:
                    forward_pos.append(pos_item)
            forward_pos.sort(key=lambda x: (x[2], x[3]))
            reverse_pos.sort(key=lambda x: (-x[2], -x[3]))

            clusters = {}
            cluster_index = 0
            for k, frag in enumerate(forward_pos):
                if not clusters.__contains__(cluster_index):
                    clusters[cluster_index] = []
                cur_cluster = clusters[cluster_index]
                if k == 0:
                    cur_cluster.append(frag)
                else:
                    is_closed = False
                    for exist_frag in reversed(cur_cluster):
                        if (frag[2] - exist_frag[3] < extend_base_threshold):
                            is_closed = True
                            break
                    if is_closed:
                        cur_cluster.append(frag)
                    else:
                        cluster_index += 1
                        if not clusters.__contains__(cluster_index):
                            clusters[cluster_index] = []
                        cur_cluster = clusters[cluster_index]
                        cur_cluster.append(frag)

            cluster_index += 1
            for k, frag in enumerate(reverse_pos):
                if not clusters.__contains__(cluster_index):
                    clusters[cluster_index] = []
                cur_cluster = clusters[cluster_index]
                if k == 0:
                    cur_cluster.append(frag)
                else:
                    is_closed = False
                    for exist_frag in reversed(cur_cluster):
                        if (exist_frag[3] - frag[2] < extend_base_threshold):
                            is_closed = True
                            break
                    if is_closed:
                        cur_cluster.append(frag)
                    else:
                        cluster_index += 1
                        if not clusters.__contains__(cluster_index):
                            clusters[cluster_index] = []
                        cur_cluster = clusters[cluster_index]
                        cur_cluster.append(frag)

            for cluster_index in clusters.keys():
                cur_cluster = clusters[cluster_index]
                cur_cluster.sort(key=lambda x: (x[0], x[1]))

                cluster_longest_query_start = -1
                cluster_longest_query_end = -1
                cluster_longest_query_len = -1

                #print('subject pos size: %d' %(len(cur_cluster)))
                # record visited fragments
                visited_frag = {}
                for i in range(len(cur_cluster)):
                    # keep a longest query start from each fragment
                    origin_frag = cur_cluster[i]
                    if visited_frag.__contains__(origin_frag):
                        continue
                    cur_frag_len = origin_frag[1] - origin_frag[0]
                    cur_longest_query_len = cur_frag_len
                    longest_query_start = origin_frag[0]
                    longest_query_end = origin_frag[1]
                    longest_subject_start = origin_frag[2]
                    longest_subject_end = origin_frag[3]

                    visited_frag[origin_frag] = 1
                    # try to extend query
                    for j in range(i + 1, len(cur_cluster)):
                        ext_frag = cur_cluster[j]
                        if visited_frag.__contains__(ext_frag):
                            continue
                        # could extend
                        # extend right
                        if ext_frag[1] > longest_query_end:
                            # judge subject direction
                            if longest_subject_start < longest_subject_end and ext_frag[2] < ext_frag[3]:
                                # +
                                if ext_frag[3] > longest_subject_end:
                                    # forward extend
                                    if ext_frag[0] - longest_query_end < extend_base_threshold and ext_frag[
                                        2] - longest_subject_end < extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start < ext_frag[
                                            2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start

                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= extend_base_threshold:
                                        break
                            elif longest_subject_start > longest_subject_end and ext_frag[2] > ext_frag[3]:
                                # reverse
                                if ext_frag[3] < longest_subject_end:
                                    # reverse extend
                                    if ext_frag[0] - longest_query_end < extend_base_threshold and longest_subject_end - \
                                            ext_frag[2] < extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start > ext_frag[
                                            2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start

                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= extend_base_threshold:
                                        break
                    if cur_longest_query_len > cluster_longest_query_len:
                        cluster_longest_query_start = longest_query_start
                        cluster_longest_query_end = longest_query_end
                        cluster_longest_query_len = cur_longest_query_len
                # keep this longest query
                if cluster_longest_query_len != -1:
                    longest_queries.append((cluster_longest_query_start, cluster_longest_query_end, cluster_longest_query_len))

        # # strategy before 06/15, too much false positive sequences
        # # generate fasta, use cd-hit-est to cluster sequences
        # local_longest_query_records = {}
        # local_longest_query_file = single_tmp_dir + '/local_longest_query_' + str(idx) + '.fa'
        # l_idx = 0
        # with open(local_longest_query_file, 'w') as f_save:
        #     for item in longest_queries:
        #         local_longest_seq = split_repeats_contigs[query_name][item[0]-1: item[1]]
        #         local_query_name = 'L_' + str(l_idx)
        #         f_save.write('>' + local_query_name + '\n' + local_longest_seq + '\n')
        #         local_longest_query_records[local_query_name] = (item[0], item[1], len(local_longest_seq))
        #         l_idx += 1
        # local_longest_query_consensus = single_tmp_dir + '/local_longest_query_' + str(idx) + '.cons.fa'
        # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.9) + ' -aL ' + str(0.9) + ' -c ' + str(
        #     0.9) + ' -G 0 -g 1 -d 0 -A 80 -i ' + local_longest_query_file + ' -o ' + local_longest_query_consensus
        # os.system(cd_hit_command)
        #
        # cluster_file = local_longest_query_consensus + '.clstr'
        # cluster_idx = -1
        # clusters = {}
        # with open(cluster_file, 'r') as f_r:
        #     for line in f_r:
        #         line = line.replace('\n', '')
        #         if line.startswith('>'):
        #             cluster_idx = line.split(' ')[1]
        #         else:
        #             if not clusters.__contains__(cluster_idx):
        #                 clusters[cluster_idx] = []
        #             cur_cluster = clusters[cluster_idx]
        #             name = line.split(',')[1].split(' ')[1].strip()[1:]
        #             name = name[0: len(name) - 3]
        #             cur_cluster.append(name)
        #             if line.endswith('*'):
        #                 clusters['rep_' + str(cluster_idx)] = name
        #
        # # local_longest_names, local_longest_contigs = read_fasta(local_longest_query_file)
        # # new_longest_queries = []
        # # for cluster_idx in clusters.keys():
        # #     if cluster_idx.isdecimal():
        # #         cluster_size = len(clusters[cluster_idx])
        # #         if cluster_size <= 5:
        # #             continue
        # #         name = clusters['rep_' + str(cluster_idx)]
        # #         new_longest_queries.append((name, cluster_size))
        # # new_longest_queries.sort(key=lambda x: -x[1])
        #
        # # sorted_longest_queries = []
        # # for item in new_longest_queries:
        # #     name = item[0]
        # #     seq = local_longest_contigs[name]
        # #     sorted_longest_queries.append((item[1], seq))
        # #longest_repeats[query_name] = sorted_longest_queries
        #
        # new_longest_queries = []
        # for cluster_idx in clusters.keys():
        #     if cluster_idx.isdecimal():
        #         cluster_size = len(clusters[cluster_idx])
        #         if cluster_size <= 1:
        #             continue
        #         name = clusters['rep_' + str(cluster_idx)]
        #         new_longest_queries.append(local_longest_query_records[name])

        # we now consider, we should take some sequences from longest_queries to represent this query sequence.
        # we take the longest sequence by length, if the latter sequence overlap with the former sequence largely (50%),
        # continue find next sequence until the ratio of query sequence over 90% or no more sequences.
        longest_queries.sort(key=lambda x: -x[2])
        query_seq = split_repeats_contigs[query_name]
        masked_query_seq = list(query_seq)
        represent_queries = []
        for cur_longest_query in longest_queries:
            cur_start = cur_longest_query[0]
            cur_end = cur_longest_query[1]
            cur_len = cur_longest_query[2]

            overlap_count = 0
            for i in range(cur_start-1, cur_end):
                if masked_query_seq[i] == 'X':
                    overlap_count += 1
                else:
                    masked_query_seq[i] = 'X'
            if float(overlap_count)/cur_len >= 0.5:
                continue
            else:
                if cur_len <= max_single_repeat_len:
                    # cur_seq = ''
                    # for i in range(cur_start-1, cur_end):
                    #     cur_seq += query_seq[i]
                    # represent_queries.append(cur_seq)
                    represent_queries.append(query_seq[cur_start-1: cur_end])

            masked_count = 0
            for j in range(len(masked_query_seq)):
                if masked_query_seq[j] == 'X':
                    masked_count += 1

            if float(masked_count)/len(masked_query_seq) >= 0.9:
                break

        longest_repeats[query_name] = represent_queries

    # print(longest_repeats)
    return longest_repeats


def get_pure_repeat_from_pool(repeat_pool, pure_repeat_len):
    new_pure_repeat = {}
    new_repeat_pool = {}
    for name in repeat_pool.keys():
        seq = repeat_pool[name]
        if len(seq) <= pure_repeat_len:
            new_pure_repeat[name] = seq
        else:
            new_repeat_pool[name] = seq
    return new_pure_repeat, new_repeat_pool


def determine_repeat_boundary(repeats_path, blast_program_dir):
    repeatNames, repeatContigs = read_fasta(repeats_path)
    # parallel
    tmp_blast_dir = tmp_output_dir + '/tmp_blast'
    os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    (repeat_dir, repeat_filename) = os.path.split(repeats_path)
    (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)

    repeat_files = []
    segments_cluster = divided_array(list(repeatContigs.items()), threads)
    for partition_index, cur_segments in enumerate(segments_cluster):
        single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
        if not os.path.exists(single_tmp_dir):
            os.makedirs(single_tmp_dir)
        split_repeat_file = single_tmp_dir + '/repeats_split.fa'
        cur_contigs = {}
        for item in cur_segments:
            cur_contigs[item[0]] = item[1]
        store_fasta(cur_contigs, split_repeat_file)
        os.system('cp ' + repeats_path + ' ' + single_tmp_dir)
        repeat_files.append((split_repeat_file, single_tmp_dir + '/' + repeat_filename,
                             single_tmp_dir + '/repeat.pairwise.out'))

    ex = ProcessPoolExecutor(threads)
    longest_repeats = {}
    jobs = []
    for file in repeat_files:
        job = ex.submit(get_longest_repeats, file, blast_program_dir, tools_dir, extend_base_threshold,
                        max_single_repeat_len)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        cur_longest_repeats = job.result()
        for query_name in cur_longest_repeats.keys():
            longest_repeats[query_name] = cur_longest_repeats[query_name]

    # # store longest_repeats for testing
    # longest_repeats_file = tmp_output_dir + '/longest_repeats.csv'
    # with codecs.open(longest_repeats_file, 'w', encoding='utf-8') as f:
    #     json.dump(longest_repeats, f)

    longest_repeats_path = tmp_output_dir + '/longest_repeats.fa'
    node_index = 0
    with open(longest_repeats_path, 'w') as f_save:
        for query_name in longest_repeats.keys():
            seqs = longest_repeats[query_name]
            for longest_seq in seqs:
                if len(longest_seq) >= 80 and not longest_seq.__contains__('N'):
                    f_save.write('>N_' + str(node_index) + '-len_' + str(len(longest_seq)) + '\n' + longest_seq + '\n')
                    node_index += 1
    return longest_repeats_path


def filter_duplication(longest_repeats_path, reference, blast_program_dir):
    (ref_dir, ref_filename) = os.path.split(reference)
    # parallel alignment
    tmp_blast_dir = tmp_output_dir + '/tmp_blast'
    os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)
    longestRepeatNames, longestRepeatContigs = read_fasta(longest_repeats_path)
    longest_repeats_files = []
    segments_cluster = divided_array(list(longestRepeatContigs.items()), threads)
    for partition_index, cur_segments in enumerate(segments_cluster):
        single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
        if not os.path.exists(single_tmp_dir):
            os.makedirs(single_tmp_dir)
        split_longest_repeat_file = single_tmp_dir + '/longest_repeats_split.fa'
        cur_contigs = {}
        for item in cur_segments:
            cur_contigs[item[0]] = item[1]
        store_fasta(cur_contigs, split_longest_repeat_file)
        os.system('cp ' + reference + ' ' + single_tmp_dir)
        longest_repeats_files.append((split_longest_repeat_file, single_tmp_dir + '/' + ref_filename,
                                      single_tmp_dir + '/longest_repeat.ref.out'))

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeats_files:
        job = ex.submit(multiple_alignment, file, blast_program_dir, tools_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    blastn2Results_path = tmp_output_dir + '/longest_repeat.ref.out'
    if os.path.exists(blastn2Results_path):
        os.remove(blastn2Results_path)

    for job in as_completed(jobs):
        cur_blastn2Results_path = job.result()
        os.system('cat ' + cur_blastn2Results_path + ' >> ' + blastn2Results_path)

    # test_array = ['N_7703-len_377', 'N_11985-len_165', 'N_9567-len_372', 'N_16763-len_732', 'N_25963-len_353', 'N_7575-len_614', 'N_23910-len_677', 'N_23982-len_450', 'N_37880-len_7446', 'N_21427-len_4703', 'N_6510-len_223', 'N_5399-len_762', 'N_12660-len_397', 'N_16621-len_258', 'N_19696-len_5266', 'N_19999-len_334', 'N_18408-len_7424', 'N_13852-len_6730', 'N_20019-len_296', 'N_5588-len_814', 'N_15501-len_524', 'N_16585-len_250', 'N_26987-len_598', 'N_19407-len_1495', 'N_36194-len_716', 'N_30063-len_515', 'N_39778-len_203', 'N_29025-len_5170', 'N_2778-len_634', 'N_26332-len_607', 'N_15941-len_665', 'N_12025-len_216', 'N_4416-len_410', 'N_26031-len_369', 'N_11507-len_454', 'N_15470-len_626', 'N_40845-len_417', 'N_3673-len_344', 'N_37672-len_675', 'N_11132-len_392', 'N_30418-len_468', 'N_2178-len_556', 'N_38528-len_327', 'N_37812-len_299', 'N_22152-len_601', 'N_36161-len_609', 'N_16792-len_324', 'N_22378-len_561', 'N_34549-len_622', 'N_29746-len_523', 'N_18060-len_448', 'N_13023-len_245', 'N_27450-len_419', 'N_21551-len_395', 'N_33549-len_1106', 'N_39743-len_1437', 'N_7787-len_506', 'N_5293-len_296', 'N_40841-len_356', 'N_3254-len_673', 'N_16635-len_429', 'N_13937-len_1776', 'N_8389-len_508', 'N_25263-len_361', 'N_18695-len_400', 'N_251-len_505', 'N_15946-len_7409', 'N_10402-len_387', 'N_23462-len_573', 'N_29371-len_870', 'N_16694-len_314', 'N_1358-len_414', 'N_10465-len_432', 'N_34878-len_7372', 'N_33736-len_269', 'N_25293-len_336', 'N_11303-len_308', 'N_13128-len_1389', 'N_17380-len_378', 'N_13106-len_217', 'N_17486-len_431', 'N_21575-len_484', 'N_10459-len_425', 'N_13277-len_432', 'N_18057-len_478', 'N_38909-len_568', 'N_37844-len_378', 'N_26769-len_296', 'N_26058-len_469', 'N_18235-len_288', 'N_12426-len_983', 'N_14824-len_200', 'N_21882-len_272', 'N_16193-len_653', 'N_19020-len_452', 'N_13827-len_381', 'N_35401-len_371', 'N_22937-len_634', 'N_16280-len_7405', 'N_36502-len_4732', 'N_29377-len_539', 'N_33577-len_656', 'N_22745-len_483', 'N_33829-len_492', 'N_7729-len_545', 'N_8545-len_7028', 'N_827-len_202', 'N_29550-len_621', 'N_4200-len_563', 'N_40226-len_395', 'N_39485-len_656', 'N_18977-len_507', 'N_6485-len_582', 'N_37741-len_623', 'N_15504-len_406', 'N_13082-len_4670', 'N_12075-len_446', 'N_25001-len_337']
    # test_dict = {}

    full_len_records_dict = {}
    # filter false positive strategy 06/27
    # test_query_name = 'N_19696-len_5266'
    longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_tandem.fa'
    blastn2Results_path = tmp_output_dir + '/longest_repeat.ref.out'
    orig_names, orig_contigs = read_fasta(longest_repeats_path)
    query_records = {}
    false_positive_queries = {}
    last_query_name = ''
    full_len_records = 0
    total_records_num = 0
    total_records_num_dict = {}
    with open(blastn2Results_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            parts = line.split('\t')
            query_name = parts[0]
            query_len = int(query_name.split('-len_')[1])
            subject_name = parts[1]
            identity = float(parts[2])
            alignment_len = int(parts[3])
            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            total_records_num += 1
            # if not query_records.__contains__(query_name):
            #     query_records[query_name] = {}
            # subject_dict = query_records[query_name]
            #
            # if not subject_dict.__contains__(subject_name):
            #     subject_dict[subject_name] = []
            # subject_pos = subject_dict[subject_name]
            # subject_pos.append((q_start, q_end, s_start, s_end))

            if query_name != last_query_name:
                total_records_num_dict[last_query_name] = total_records_num
                total_records_num = 0
                # if last_query_name in test_array:
                #     if not test_dict.__contains__(last_query_name):
                #         test_dict[last_query_name] = []
                #     array = test_dict[last_query_name]
                #     array.append(full_len_records)
                #     #print('last_query_name: %s, full_len_records: %d' %(last_query_name, full_len_records))
                full_len_records_dict[last_query_name] = full_len_records
                # if full_len_records >= 10:
                #     false_positive_queries[last_query_name] = 1
                # if query_records.__contains__(last_query_name):
                #     del query_records[last_query_name]
                full_len_records = 0
            if abs(q_start - 1) <= 10 and abs(q_end - query_len) <= 10:
                full_len_records += 1

            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))
            last_query_name = query_name
        full_len_records_dict[last_query_name] = full_len_records
        # if full_len_records >= 10:
        #     false_positive_queries[last_query_name] = 1
        total_records_num_dict[last_query_name] = total_records_num

    longest_repeats = {}
    for idx, query_name in enumerate(query_records.keys()):
        subject_dict = query_records[query_name]
        query_len = int(query_name.split('-len_')[1])
        extend_base_threshold = int(extend_base_ratio * query_len)
        # if there are more than one longest query overlap with the final longest query over 90%,
        # then it probably the true TE
        longest_queries = []
        for subject_name in subject_dict.keys():
            subject_pos = subject_dict[subject_name]
            # subject_pos.sort(key=lambda x: (x[2], x[3]))

            # cluster all closed fragments, split forward and reverse records
            forward_pos = []
            reverse_pos = []
            for pos_item in subject_pos:
                if pos_item[2] > pos_item[3]:
                    reverse_pos.append(pos_item)
                else:
                    forward_pos.append(pos_item)
            forward_pos.sort(key=lambda x: (x[2], x[3]))
            reverse_pos.sort(key=lambda x: (-x[2], -x[3]))

            clusters = {}
            cluster_index = 0
            for k, frag in enumerate(forward_pos):
                if not clusters.__contains__(cluster_index):
                    clusters[cluster_index] = []
                cur_cluster = clusters[cluster_index]
                if k == 0:
                    cur_cluster.append(frag)
                else:
                    is_closed = False
                    for exist_frag in reversed(cur_cluster):
                        if (frag[2] - exist_frag[3] < extend_base_threshold):
                            is_closed = True
                            break
                    if is_closed:
                        cur_cluster.append(frag)
                    else:
                        cluster_index += 1
                        if not clusters.__contains__(cluster_index):
                            clusters[cluster_index] = []
                        cur_cluster = clusters[cluster_index]
                        cur_cluster.append(frag)

            cluster_index += 1
            for k, frag in enumerate(reverse_pos):
                if not clusters.__contains__(cluster_index):
                    clusters[cluster_index] = []
                cur_cluster = clusters[cluster_index]
                if k == 0:
                    cur_cluster.append(frag)
                else:
                    is_closed = False
                    for exist_frag in reversed(cur_cluster):
                        if (exist_frag[3] - frag[2] < extend_base_threshold):
                            is_closed = True
                            break
                    if is_closed:
                        cur_cluster.append(frag)
                    else:
                        cluster_index += 1
                        if not clusters.__contains__(cluster_index):
                            clusters[cluster_index] = []
                        cur_cluster = clusters[cluster_index]
                        cur_cluster.append(frag)

            for cluster_index in clusters.keys():
                cur_cluster = clusters[cluster_index]
                cur_cluster.sort(key=lambda x: (x[0], x[1]))

                # cluster_longest_query_start = -1
                # cluster_longest_query_end = -1
                # cluster_longest_query_len = -1

                # record visited fragments
                visited_frag = {}
                for i in range(len(cur_cluster)):
                    # keep a longest query start from each fragment
                    origin_frag = cur_cluster[i]
                    if visited_frag.__contains__(origin_frag):
                        continue
                    cur_frag_len = origin_frag[1] - origin_frag[0]
                    cur_longest_query_len = cur_frag_len
                    longest_query_start = origin_frag[0]
                    longest_query_end = origin_frag[1]
                    longest_subject_start = origin_frag[2]
                    longest_subject_end = origin_frag[3]
                    frag_num = 1

                    visited_frag[origin_frag] = 1
                    # try to extend query
                    for j in range(i + 1, len(cur_cluster)):
                        ext_frag = cur_cluster[j]
                        if visited_frag.__contains__(ext_frag):
                            continue
                        # could extend
                        # extend right
                        if ext_frag[1] > longest_query_end:
                            # judge subject direction
                            if longest_subject_start < longest_subject_end and ext_frag[2] < ext_frag[3]:
                                # +
                                if ext_frag[3] > longest_subject_end:
                                    # forward extend
                                    if ext_frag[0] - longest_query_end < extend_base_threshold and ext_frag[
                                        2] - longest_subject_end < extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start

                                        frag_num += 1

                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= extend_base_threshold:
                                        break
                            elif longest_subject_start > longest_subject_end and ext_frag[2] > ext_frag[3]:
                                # reverse
                                if ext_frag[3] < longest_subject_end:
                                    # reverse extend
                                    if ext_frag[0] - longest_query_end < extend_base_threshold and longest_subject_end - \
                                            ext_frag[2] < extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start > \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start

                                        frag_num += 1

                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= extend_base_threshold:
                                        break
                    if cur_longest_query_len > 0:
                        longest_query_item = (
                            longest_query_start, longest_query_end, cur_longest_query_len, longest_subject_start,
                            longest_subject_end, frag_num)
                        if longest_query_item[2] <= 20000 and abs(longest_subject_end-longest_subject_start) <= 30000 and float(longest_query_item[2]) / query_len >= 0.9:
                            longest_queries.append(longest_query_item)
                            # is_replace = False
                            # for item in longest_queries.keys():
                            #     # subject direction should be same
                            #     direct = ''
                            #     if longest_query_item[3] < longest_query_item[4] and item[3] < item[4]:
                            #         direct = '+'
                            #     elif (longest_query_item[3] > longest_query_item[4] and item[3] > item[4]):
                            #         direct = '-'
                            #     if direct == '':
                            #         continue
                            #
                            #     # longer sequence will replace the shorter one
                            #     if abs(longest_query_item[4] - longest_query_item[3]) > abs(item[4] - item[3]):
                            #         if direct == '+':
                            #             # if overlap
                            #             if (longest_query_item[4] > item[3] and longest_query_item[4] < item[4]) \
                            #                     or (longest_query_item[3] > item[3] and longest_query_item[3] < item[4]):
                            #                 is_replace = True
                            #                 longest_queries[item] = -1
                            #                 # del longest_queries[item]
                            #                 # longest_queries[longest_query_item] = 1
                            #         elif direct == '-':
                            #             # if overlap
                            #             if (longest_query_item[3] > item[4] and longest_query_item[3] < item[3]) \
                            #                     or (longest_query_item[4] > item[4] and longest_query_item[4] < item[3]):
                            #                 is_replace = True
                            #                 longest_queries[item] = -1
                            #                 # del longest_queries[item]
                            #                 # longest_queries[longest_query_item] = 1
                            # new_longest_queries = {}
                            # for item in longest_queries.keys():
                            #     if longest_queries[item] == -1:
                            #         continue
                            #     new_longest_queries[item] = 1
                            # new_longest_queries[longest_query_item] = 1
                            # longest_queries = new_longest_queries
                #     if cur_longest_query_len > cluster_longest_query_len:
                #         cluster_longest_query_start = longest_query_start
                #         cluster_longest_query_end = longest_query_end
                #         cluster_longest_query_len = cur_longest_query_len
                # # keep this longest query
                # if cluster_longest_query_len != -1:
                #     longest_queries.append((cluster_longest_query_start, cluster_longest_query_end, cluster_longest_query_len))

        # # remove overlap copies
        # longest_queries_forward = []
        # longest_queries_reverse = []
        # for longest_query in longest_queries:
        #     if longest_query[3] < longest_query[4]:
        #         longest_queries_forward.append(longest_query)
        #     else:
        #         longest_queries_reverse.append(longest_query)
        #
        # longest_queries = set()
        # is_forward_overlap = [0 for i in range(len(longest_queries_forward))]
        # for i in range(len(longest_queries_forward)):
        #     cur_longest_query = longest_queries_forward[i]
        #     for j in range(len(longest_queries_forward)):
        #         if i == j or is_forward_overlap[i] == 1 or is_forward_overlap[j] == 1:
        #             continue
        #         next_longest_query = longest_queries_forward[j]
        #         # if overlap
        #         if (next_longest_query[4] > cur_longest_query[3] and next_longest_query[4] < cur_longest_query[4]) \
        #                 or (next_longest_query[3] > cur_longest_query[3] and next_longest_query[3] < cur_longest_query[4]):
        #             if cur_longest_query[4] - cur_longest_query[3] > next_longest_query[4] - next_longest_query[3]:
        #                 is_forward_overlap[j] = 1
        #             else:
        #                 is_forward_overlap[i] = 1
        #
        # is_reverse_overlap = [0 for i in range(len(longest_queries_reverse))]
        # for i in range(len(longest_queries_reverse)):
        #     cur_longest_query = longest_queries_reverse[i]
        #     for j in range(len(longest_queries_reverse)):
        #         if i == j or is_reverse_overlap[i] == 1 or is_reverse_overlap[j] == 1:
        #             continue
        #         next_longest_query = longest_queries_reverse[j]
        #         # if overlap
        #         if (next_longest_query[3] > cur_longest_query[4] and next_longest_query[3] < cur_longest_query[3]) \
        #                 or (next_longest_query[4] > cur_longest_query[4] and next_longest_query[4] < cur_longest_query[3]):
        #             if cur_longest_query[3] - cur_longest_query[4] > next_longest_query[3] - next_longest_query[4]:
        #                 is_reverse_overlap[j] = 1
        #             else:
        #                 is_reverse_overlap[i] = 1
        #
        # for k in range(len(longest_queries_forward)):
        #     if is_forward_overlap[k] == 1:
        #         continue
        #     longest_queries.add(longest_queries_forward[k])
        # for k in range(len(longest_queries_reverse)):
        #     if is_reverse_overlap[k] == 1:
        #         continue
        #     longest_queries.add(longest_queries_reverse[k])

        full_len_query_num = len(longest_queries)
        # print(full_len_query_num)
        # for longest_query in longest_queries.keys():
        #     if float(longest_query[2])/query_len >= 0.9:
        #         full_len_query_num += 1
        # if full_len_query_num < 10:
        #     false_positive_queries[query_name] = 1

        # full_len_records = full_len_records_dict[query_name]
        # full_connect_num = full_len_query_num - full_len_records
        #
        # total_records_num = total_records_num_dict[query_name]
        # true_connect_query_num = 0
        # for longest_query in longest_queries:
        #     if longest_query[5] >= 3:
        #         true_connect_query_num += 1

        # test_query_names = ['N_12004-len_871', 'N_14936-len_900', 'N_5106-len_742', 'N_34021-len_957', 'N_39546-len_572', 'N_35242-len_781']
        # test_query_name = 'N_12004-len_871'
        # if query_name in test_query_names:
        #     print('total_records_num: %d, full_len_query_num: %d, full_len_records: %d, full_connect_num: %d, true_connect_query_num: %d'
        #     %(total_records_num, full_len_query_num, full_len_records, full_connect_num, true_connect_query_num))
        #     #print(longest_queries)


        # if full_len_records < 4:
        #     false_positive_queries[query_name] = 1

        # if full_len_query_num < 2 or (query_len < 1000 and true_connect_query_num < 2) or full_len_records > 2*full_connect_num:
        if full_len_query_num < 3:
        #if true_connect_query_num < 2:
            false_positive_queries[query_name] = 1


        # if true_connect_query_num > 10 and false_positive_queries.__contains__(query_name):
        #     del false_positive_queries[query_name]
        # if full_len_query_num < 2 or (query_len < 1000 and full_len_records < 5):
        #     false_positive_queries[query_name] = 1

        # if false_positive_queries.__contains__(query_name) and full_len_query_num >= 20:
        #     del false_positive_queries[query_name]
        # if false_positive_queries.__contains__(query_name) and (full_len_query_num >= 20 or (full_len_query_num-full_len_records) >= full_len_records):
        #     del false_positive_queries[query_name]
        # if false_positive_queries.__contains__(query_name) and full_len_query_num >= 20:
        #     del false_positive_queries[query_name]

    #     if query_name in test_array:
    #         if not test_dict.__contains__(query_name):
    #             test_dict[query_name] = []
    #         array = test_dict[query_name]
    #         array.append(full_len_query_num)
    #         #print('query_name: %s, full_len_query_num: %d' %(query_name, full_len_query_num))
    # print(test_dict)
    # count = 0
    # true_positive_path = tmp_output_dir + '/true_positive.fa'
    # with open(true_positive_path, 'r') as f_r:
    #     for line in f_r:
    #         query_name = line.split('\t')[0]
    #         if false_positive_queries.__contains__(query_name):
    #             print(query_name)
    #             count += 1
    # print(count)

    seq_count = 0
    longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_duplication.fa'
    with open(longest_repeats_path, 'w') as f_save:
        for name in query_records.keys():
            if not false_positive_queries.__contains__(name):
                f_save.write('>' + name + '\n' + orig_contigs[name] + '\n')
                seq_count += 1
    print(seq_count)
    return longest_repeats_path


def get_REXdb_aligned_seq(blast_program_dir, protein_db_path, longest_repeats_path):
    (protein_dir, protein_filename) = os.path.split(protein_db_path)
    # parallel alignment
    tmp_blast_dir = tmp_output_dir + '/tmp_blast'
    os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)
    longestRepeatNames, longestRepeatContigs = read_fasta(longest_repeats_path)
    longest_repeats_files = []
    segments_cluster = divided_array(list(longestRepeatContigs.items()), threads)
    for partition_index, cur_segments in enumerate(segments_cluster):
        single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
        if not os.path.exists(single_tmp_dir):
            os.makedirs(single_tmp_dir)
        split_longest_repeat_file = single_tmp_dir + '/longest_repeats_split.fa'
        cur_contigs = {}
        for item in cur_segments:
            cur_contigs[item[0]] = item[1]
        store_fasta(cur_contigs, split_longest_repeat_file)
        os.system('cp ' + protein_db_path + ' ' + single_tmp_dir)
        longest_repeats_files.append((split_longest_repeat_file, single_tmp_dir + '/' + protein_filename,
                                      single_tmp_dir + '/longest_repeat.REXdb.out'))

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeats_files:
        job = ex.submit(multiple_alignment_blastx, file, blast_program_dir, tools_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    blastx2Results_path = longest_repeats_path + '.REXdb.tmpBlastXResults.out'
    if os.path.exists(blastx2Results_path):
        os.remove(blastx2Results_path)

    for job in as_completed(jobs):
        cur_blastx2Results_path = job.result()
        os.system('cat ' + cur_blastx2Results_path + ' >> ' + blastx2Results_path)
    return blastx2Results_path



def filter_derived_seq_v1(longest_repeats_path, blast_program_dir):
    tmp_blast_dir = tmp_output_dir + '/tmp_blast'
    os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    (repeat_dir, repeat_filename) = os.path.split(longest_repeats_path)

    orig_names, orig_contigs = read_fasta(longest_repeats_path)

    longest_repeat_files = []
    segments_cluster = divided_array(list(orig_contigs.items()), threads)
    for partition_index, cur_segments in enumerate(segments_cluster):
        single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
        if not os.path.exists(single_tmp_dir):
            os.makedirs(single_tmp_dir)
        split_repeat_file = single_tmp_dir + '/repeats_split.fa'
        cur_contigs = {}
        for item in cur_segments:
            cur_contigs[item[0]] = item[1]
        store_fasta(cur_contigs, split_repeat_file)
        os.system('cp ' + longest_repeats_path + ' ' + single_tmp_dir)
        longest_repeat_files.append((split_repeat_file, single_tmp_dir + '/' + repeat_filename,
                                     single_tmp_dir + '/longest_repeat.pairwise.out'))

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeat_files:
        job = ex.submit(multiple_alignment, file, blast_program_dir, tools_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    blastn2Results_path = tmp_output_dir + '/longest_repeat.pairwise.out'
    if os.path.exists(blastn2Results_path):
        os.remove(blastn2Results_path)

    for job in as_completed(jobs):
        cur_blastn2Results_path = job.result()
        os.system('cat ' + cur_blastn2Results_path + ' >> ' + blastn2Results_path)

    # filter derived sequences, which is caused by deletion of intact TE
    # if a sequence is contained in another sequence continiously or partly, it is a derived sequence
    # longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_duplication.fa'

    query_records = {}
    blastn2Results_path = tmp_output_dir + '/longest_repeat.pairwise.out'
    derived_queries = {}
    with open(blastn2Results_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            parts = line.split('\t')
            query_name = parts[0]
            query_len = int(query_name.split('-len_')[1])
            subject_name = parts[1]
            subject_len = int(subject_name.split('-len_')[1])
            identity = float(parts[2])
            alignment_len = int(parts[3])
            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            if query_name == subject_name or query_len > subject_len:
                continue

            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end, identity))

    for query_name in query_records.keys():
        query_len = int(query_name.split('-len_')[1])
        subject_dict = query_records[query_name]
        for subject_name in subject_dict.keys():
            query_seq = orig_contigs[query_name]
            query_seq_array = list(query_seq)

            subject_pos = subject_dict[subject_name]
            subject_pos.sort(key=lambda x: (x[0], x[1]))

            # As long as a sequence is included in another sequence, whether it is continuous or discrete in another sequence, then this sequence is redundant
            for i in range(len(subject_pos)):
                # keep a longest query start from each fragment
                origin_frag = subject_pos[i]
                longest_query_start = origin_frag[0]
                longest_query_end = origin_frag[1]
                for j in range(longest_query_start-1, longest_query_end):
                    query_seq_array[j] = 'X'

            masked_count = 0
            for k in range(len(query_seq_array)):
                if query_seq_array[k] == 'X':
                    masked_count += 1

            if float(masked_count) / query_len >= 0.99:
                derived_queries[query_name] = 1
                break

    print('derived seq size: %d' % (len(derived_queries)))

    longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_derived.fa'
    with open(longest_repeats_path, 'w') as f_save:
        for name in orig_names:
            if not derived_queries.__contains__(name):
                f_save.write('>' + name + '\n' + orig_contigs[name] + '\n')
    return longest_repeats_path


def filter_derived_seq(longest_repeats_path, blast_program_dir):
    tmp_blast_dir = tmp_output_dir + '/tmp_blast'
    os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    (repeat_dir, repeat_filename) = os.path.split(longest_repeats_path)

    longest_repeats_Names, longest_repeats_Contigs = read_fasta(longest_repeats_path)

    longest_repeat_files = []
    segments_cluster = divided_array(list(longest_repeats_Contigs.items()), threads)
    for partition_index, cur_segments in enumerate(segments_cluster):
        single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
        if not os.path.exists(single_tmp_dir):
            os.makedirs(single_tmp_dir)
        split_repeat_file = single_tmp_dir + '/repeats_split.fa'
        cur_contigs = {}
        for item in cur_segments:
            cur_contigs[item[0]] = item[1]
        store_fasta(cur_contigs, split_repeat_file)
        os.system('cp ' + longest_repeats_path + ' ' + single_tmp_dir)
        longest_repeat_files.append((split_repeat_file, single_tmp_dir + '/' + repeat_filename,
                                     single_tmp_dir + '/longest_repeat.pairwise.out'))

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeat_files:
        job = ex.submit(multiple_alignment, file, blast_program_dir, tools_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    blastn2Results_path = tmp_output_dir + '/longest_repeat.pairwise.out'
    if os.path.exists(blastn2Results_path):
        os.remove(blastn2Results_path)

    for job in as_completed(jobs):
        cur_blastn2Results_path = job.result()
        os.system('cat ' + cur_blastn2Results_path + ' >> ' + blastn2Results_path)

    # test_array = ['N_4243-len_7737', 'N_21637-len_426', 'N_21427-len_4703', 'N_3254-len_673', 'N_11132-len_392', 'N_38300-len_8735', 'N_12512-len_1050', 'N_9567-len_372', 'N_12025-len_216', 'N_33549-len_1106', 'N_37080-len_7122', 'N_26011-len_336', 'N_28977-len_608', 'N_10459-len_425', 'N_3483-len_4406', 'N_20566-len_282', 'N_21596-len_1053', 'N_36277-len_1459', 'N_21575-len_484', 'N_35361-len_7067', 'N_15501-len_524', 'N_11590-len_6097', 'N_12660-len_397', 'N_35581-len_8954', 'N_11473-len_4820', 'N_251-len_505', 'N_27450-len_419', 'N_4659-len_6533', 'N_16694-len_314', 'N_40226-len_395', 'N_18695-len_400', 'N_8745-len_324', 'N_16349-len_5098', 'N_33268-len_6941', 'N_827-len_202', 'N_29752-len_1448', 'N_16621-len_258', 'N_39743-len_1437', 'N_34549-len_622', 'N_29616-len_6281', 'N_12424-len_1417', 'N_11077-len_2066', 'N_20539-len_4537', 'N_30511-len_6541', 'N_17572-len_516', 'N_14999-len_1106', 'N_6991-len_4330', 'N_20492-len_4164', 'N_17718-len_1336', 'N_21010-len_427', 'N_7787-len_506', 'N_16585-len_250', 'N_32703-len_4516', 'N_24433-len_266', 'N_34623-len_1177', 'N_6905-len_378', 'N_32348-len_4225', 'N_26144-len_1989', 'N_22881-len_1138', 'N_6467-len_7421', 'N_11303-len_308', 'N_4416-len_410', 'N_17380-len_378', 'N_1371-len_7342', 'N_28590-len_5509', 'N_26728-len_9205', 'N_18057-len_478', 'N_22745-len_483', 'N_28719-len_5046', 'N_30556-len_4758', 'N_11507-len_454', 'N_7575-len_614', 'N_6510-len_223', 'N_3040-len_541', 'N_15470-len_626', 'N_6225-len_372', 'N_25293-len_336', 'N_39778-len_203', 'N_29727-len_2522', 'N_15671-len_5561', 'N_23910-len_677', 'N_30171-len_8342', 'N_21490-len_1591', 'N_41186-len_546', 'N_30214-len_8147', 'N_17719-len_609', 'N_40707-len_4340', 'N_5293-len_296', 'N_2835-len_7378', 'N_34932-len_472', 'N_31533-len_7701', 'N_11991-len_6803', 'N_8389-len_508', 'N_19235-len_1458', 'N_39879-len_2874', 'N_38749-len_425', 'N_15550-len_1462', 'N_13277-len_432', 'N_37880-len_7446', 'N_33736-len_269', 'N_20700-len_5123', 'N_9419-len_5034', 'N_27550-len_5989', 'N_37399-len_1936', 'N_24706-len_306', 'N_40969-len_7672', 'N_1811-len_983', 'N_31166-len_852', 'N_4778-len_409', 'N_30614-len_6003', 'N_13128-len_1389', 'N_428-len_7244', 'N_25001-len_337', 'N_15946-len_7409', 'N_39197-len_726', 'N_40860-len_4424', 'N_34000-len_7039', 'N_17615-len_8038', 'N_34386-len_4993', 'N_37447-len_4452', 'N_16678-len_2406', 'N_2109-len_1495', 'N_1358-len_414', 'N_27831-len_8215', 'N_15941-len_665', 'N_14688-len_1142', 'N_21849-len_1059', 'N_37913-len_1078', 'N_27779-len_470', 'N_23462-len_573', 'N_3132-len_7353', 'N_28491-len_1045', 'N_12920-len_1325', 'N_1845-len_415', 'N_13852-len_6730', 'N_32849-len_390', 'N_23621-len_6282', 'N_11081-len_425', 'N_39650-len_614', 'N_13729-len_5390', 'N_30572-len_7220', 'N_20491-len_4646', 'N_22459-len_1200', 'N_7191-len_384', 'N_17783-len_1028', 'N_6830-len_5896', 'N_21551-len_395', 'N_1026-len_5293', 'N_34380-len_1349', 'N_8521-len_5719', 'N_10773-len_8623', 'N_2833-len_8476', 'N_28398-len_352', 'N_24254-len_7308', 'N_26769-len_296', 'N_29339-len_357', 'N_36502-len_4732', 'N_30615-len_5565', 'N_9675-len_802', 'N_40794-len_4147', 'N_3481-len_1221', 'N_3532-len_6469', 'N_33829-len_492', 'N_6923-len_4674', 'N_34836-len_6769', 'N_22536-len_5831', 'N_1620-len_7196', 'N_21797-len_1603', 'N_3009-len_4500', 'N_8544-len_602', 'N_2122-len_1030', 'N_12426-len_983', 'N_13051-len_1148', 'N_39205-len_414', 'N_3673-len_344', 'N_11633-len_4120', 'N_24450-len_7452', 'N_22473-len_565', 'N_5518-len_4510', 'N_22456-len_6486', 'N_25263-len_361', 'N_7985-len_783', 'N_5914-len_4493', 'N_25491-len_4396', 'N_35208-len_7128', 'N_29377-len_539', 'N_25523-len_4440', 'N_3141-len_1077', 'N_12454-len_6681', 'N_19407-len_1495', 'N_8659-len_409', 'N_35360-len_7806', 'N_6364-len_5917', 'N_13113-len_324', 'N_25963-len_353', 'N_7891-len_1010', 'N_38053-len_315', 'N_34703-len_5112', 'N_32203-len_14076', 'N_8536-len_305', 'N_2518-len_5129', 'N_18977-len_507', 'N_13734-len_6167', 'N_19007-len_4651', 'N_40749-len_309', 'N_27567-len_1063', 'N_33751-len_458', 'N_6485-len_582', 'N_16635-len_429', 'N_13707-len_554', 'N_12541-len_1388', 'N_15504-len_406', 'N_19063-len_397', 'N_22152-len_601', 'N_7886-len_544', 'N_1272-len_4073', 'N_21035-len_5302', 'N_19999-len_334', 'N_812-len_424', 'N_37741-len_623', 'N_1076-len_525', 'N_26031-len_369', 'N_21018-len_5679', 'N_7247-len_3888', 'N_19608-len_613', 'N_30331-len_3719', 'N_6128-len_7764', 'N_11367-len_1637', 'N_7957-len_1715', 'N_22628-len_8326', 'N_16709-len_407', 'N_12165-len_7054', 'N_11494-len_6342', 'N_7633-len_7418', 'N_33585-len_1538', 'N_14063-len_7581', 'N_7729-len_545', 'N_3510-len_4146', 'N_22326-len_8605', 'N_781-len_5814', 'N_21561-len_1026', 'N_29943-len_5544', 'N_13827-len_381', 'N_20766-len_471', 'N_13106-len_217', 'N_16816-len_1351', 'N_23524-len_4111', 'N_28509-len_987', 'N_12534-len_614', 'N_38528-len_327', 'N_30616-len_5104', 'N_153-len_506', 'N_17638-len_6149', 'N_32691-len_308', 'N_22378-len_561', 'N_18235-len_288', 'N_8373-len_470', 'N_3059-len_4615', 'N_29550-len_621', 'N_32820-len_6045', 'N_8563-len_7546', 'N_40845-len_417', 'N_33975-len_520', 'N_37403-len_1029', 'N_10465-len_432', 'N_25950-len_4607', 'N_8303-len_378', 'N_39633-len_399', 'N_26058-len_469', 'N_17486-len_431', 'N_10467-len_4295', 'N_968-len_6356', 'N_32672-len_5191', 'N_34878-len_7372', 'N_13082-len_4670', 'N_30418-len_468', 'N_13023-len_245', 'N_38406-len_5059', 'N_32326-len_1372', 'N_30602-len_5691', 'N_24193-len_474', 'N_831-len_1381', 'N_34927-len_289', 'N_2541-len_7128']

    # filter derived sequences, which is caused by deletion of intact TE
    # if a sequence is contained in another sequence continiously or partly, it is a derived sequence
    longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_duplication.fa'
    orig_names, orig_contigs = read_fasta(longest_repeats_path)
    query_records = {}
    blastn2Results_path = tmp_output_dir + '/longest_repeat.pairwise.out'
    derived_queries = {}
    with open(blastn2Results_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            parts = line.split('\t')
            query_name = parts[0]
            query_len = int(query_name.split('-len_')[1])
            subject_name = parts[1]
            subject_len = int(subject_name.split('-len_')[1])
            identity = float(parts[2])
            alignment_len = int(parts[3])
            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            if query_name == subject_name or query_len > subject_len:
                continue

            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end, identity))

    for query_name in query_records.keys():
        # if query_name == 'N_27223-len_1686':
        #     print('here')
        query_len = int(query_name.split('-len_')[1])
        extend_base_threshold = int(query_len * extend_base_ratio)
        subject_dict = query_records[query_name]
        for subject_name in subject_dict.keys():
            longest_queries = []
            subject_pos = subject_dict[subject_name]
            subject_pos.sort(key=lambda x: (x[0], x[1]))

            # cluster_longest_query_start = -1
            # cluster_longest_query_end = -1
            # cluster_longest_query_len = -1

            # record visited fragments
            visited_frag = {}
            for i in range(len(subject_pos)):
                # keep a longest query start from each fragment
                origin_frag = subject_pos[i]
                if visited_frag.__contains__(origin_frag):
                    continue
                cur_frag_len = origin_frag[1] - origin_frag[0]
                cur_longest_query_len = cur_frag_len
                longest_query_start = origin_frag[0]
                longest_query_end = origin_frag[1]
                longest_subject_start = origin_frag[2]
                longest_subject_end = origin_frag[3]

                total_identity = origin_frag[4]
                identity_count = 1

                visited_frag[origin_frag] = 1
                # try to extend query
                for j in range(i + 1, len(subject_pos)):
                    ext_frag = subject_pos[j]
                    if visited_frag.__contains__(ext_frag):
                        continue
                    # could extend
                    # extend right
                    if ext_frag[1] > longest_query_end:
                        # judge subject direction
                        if longest_subject_start < longest_subject_end and ext_frag[2] < ext_frag[3]:
                            # +
                            if ext_frag[3] > longest_subject_end:
                                # forward extend
                                if ext_frag[0] - longest_query_end < extend_base_threshold and ext_frag[
                                    2] > longest_subject_end:
                                    # update the longest path
                                    longest_query_start = longest_query_start
                                    longest_query_end = ext_frag[1]
                                    longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                     ext_frag[2] else ext_frag[2]
                                    longest_subject_end = ext_frag[3]
                                    cur_longest_query_len = longest_query_end - longest_query_start

                                    total_identity += ext_frag[4]
                                    identity_count += 1

                                    visited_frag[ext_frag] = 1
                                elif ext_frag[0] - longest_query_end >= extend_base_threshold:
                                    break
                        elif longest_subject_start > longest_subject_end and ext_frag[2] > ext_frag[3]:
                            # reverse
                            if ext_frag[3] < longest_subject_end:
                                # reverse extend
                                if ext_frag[0] - longest_query_end < extend_base_threshold and longest_subject_end > \
                                        ext_frag[2]:
                                    # update the longest path
                                    longest_query_start = longest_query_start
                                    longest_query_end = ext_frag[1]
                                    longest_subject_start = longest_subject_start if longest_subject_start > \
                                                                                     ext_frag[2] else ext_frag[2]
                                    longest_subject_end = ext_frag[3]
                                    cur_longest_query_len = longest_query_end - longest_query_start

                                    total_identity += ext_frag[4]
                                    identity_count += 1

                                    visited_frag[ext_frag] = 1
                                elif ext_frag[0] - longest_query_end >= extend_base_threshold:
                                    break
                if cur_longest_query_len > 0:
                    longest_query_item = (longest_query_start, longest_query_end, cur_longest_query_len,
                                          float(total_identity) / identity_count)
                    # if query_name == 'N_31163-len_1026':
                    #     print(query_name, subject_name, longest_query_item)
                    if float(longest_query_item[2]) / query_len >= 0.99 and longest_query_item[3] >= 95:
                        longest_queries.append(longest_query_item)
            #     if cur_longest_query_len > cluster_longest_query_len:
            #         cluster_longest_query_start = longest_query_start
            #         cluster_longest_query_end = longest_query_end
            #         cluster_longest_query_len = cur_longest_query_len
            # # keep this longest query
            # if cluster_longest_query_len != -1:
            #     longest_queries.append((cluster_longest_query_start, cluster_longest_query_end, cluster_longest_query_len))

            # if query_name == 'N_5658-len_890':
            #     print(query_name, subject_name, longest_queries)

            if len(longest_queries) > 0:
                derived_queries[query_name] = 1
                break
    print('derived seq size: %d' % (len(derived_queries)))

    longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_derived.fa'
    with open(longest_repeats_path, 'w') as f_save:
        for name in orig_names:
            if not derived_queries.__contains__(name):
                f_save.write('>' + name + '\n' + orig_contigs[name] + '\n')
    return longest_repeats_path


def refine_LTR_sequences():
    longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_derived.fa'
    longest_repeats_multi_line_path = tmp_output_dir + '/longest_repeats.filter_derived.ml.fa'
    contigNames, contigs = read_fasta(longest_repeats_path)
    outfile = open(longest_repeats_multi_line_path, 'w')  # open outfile for writing
    for name in contigNames:
        print_seqs(name, contigs[name], 50, outfile)

    TRsearch_dir = '/home/hukang/repeat_detect_tools/REPET_linux-x64-3.0/bin'
    TRsearch_command = TRsearch_dir + '/TRsearch ' + longest_repeats_multi_line_path
    #os.system(TRsearch_command)

    TR_out = longest_repeats_multi_line_path + '.TR.set'
    # -----------------------------find sequences with obvious LTR structure--------------------------------------
    TE_merge_contigNames, TE_merge_contigs = read_fasta(longest_repeats_multi_line_path)
    type_sets = set()
    # query_structs = {query_name: {repeat_id: [r1, r2]} }
    LTR_query_structs = {}
    with open(TR_out, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            repeat_id = parts[0]
            query_name = parts[2]
            query_len = int(query_name.split('-len_')[1])
            start_pos = int(parts[3])
            end_pos = int(parts[4])
            struct_info = parts[1]
            struct_parts = struct_info.split('|')
            s_type = struct_parts[0].replace('Rep' + str(repeat_id), '')
            s_len = int(struct_parts[1].replace('len=', ''))
            s_identity = float(struct_parts[2].replace('id=', ''))
            TE_len = int(struct_parts[3].replace('lenTE=', ''))
            struct_tuple = (s_type, s_len, s_identity, TE_len, start_pos, end_pos)
            type_sets.add(s_type)

            if float(struct_tuple[3]) / query_len >= 0.95:
                if (struct_tuple[0] == 'LTR' or struct_tuple[0] == 'termLTR') \
                        and (struct_tuple[1] >= 100 and struct_tuple[1] <= 5000):
                    if not LTR_query_structs.__contains__(query_name):
                        LTR_query_structs[query_name] = {}
                    struct_records = LTR_query_structs[query_name]
                    if not struct_records.__contains__(repeat_id):
                        struct_records[repeat_id] = []
                    struct_tuples = struct_records[repeat_id]
                    struct_tuples.append(struct_tuple)

    # print(type_sets)
    # choose the most probably LTR (longer in LTR len)
    best_LTR_query_structs = {}
    for query_name in LTR_query_structs.keys():
        struct_records = LTR_query_structs[query_name]
        query_len = int(query_name.split('-len_')[1])
        if len(struct_records) <= 1:
            best_LTR_query_structs[query_name] = list(struct_records.values())[0]
        else:
            best_record = None
            for repeat_id in struct_records.keys():
                struct_tuples = struct_records[repeat_id]
                first_tuple = struct_tuples[0]
                second_tuple = struct_tuples[1]
                if best_record is None:
                    best_record = struct_tuples
                elif first_tuple[1] > best_record[0][1]:
                    best_record = struct_tuples
            best_LTR_query_structs[query_name] = best_record

    LTR_no_struct = {}
    clear_LTR_seqs = {}
    for query_name in best_LTR_query_structs:
        best_record = best_LTR_query_structs[query_name]
        first_tuple = best_record[0]
        second_tuple = best_record[1]

        leftLTR_start = first_tuple[4]
        leftLTR_end = first_tuple[5]

        ltr_len = int(first_tuple[1])
        ltr_identity = float(first_tuple[2])
        TE_len = int(first_tuple[3])

        rightLTR_start = second_tuple[4]
        rightLTR_end = second_tuple[5]
        query_seq = TE_merge_contigs[query_name]
        clear_LTR_seq = query_seq[leftLTR_start - 1: rightLTR_end]
        # print(, TE_len)
        # LTR_info = (left_ltr_start, left_ltr_end, right_ltr_start, right_ltr_end, complete_LTR-RT_len, dna_sequence)
        LTR_info = (1, ltr_len, len(clear_LTR_seq) - ltr_len + 1, len(clear_LTR_seq), len(clear_LTR_seq), clear_LTR_seq)
        # clear_LTR_seqs['N_'+str(node_index)+'-len_'+str(len(clear_LTR_seq))] = clear_LTR_seq
        clear_LTR_seqs[query_name] = clear_LTR_seq
        LTR_no_struct[query_name] = query_seq[leftLTR_end: rightLTR_start - 1]

    longest_repeats_path = tmp_output_dir + '/longest_repeats.refine.fa'
    refine_contigs = {}
    for query_name in TE_merge_contigNames:
        if clear_LTR_seqs.__contains__(query_name):
            refine_contigs[query_name] = clear_LTR_seqs[query_name]
        else:
            refine_contigs[query_name] = TE_merge_contigs[query_name]
    store_fasta(refine_contigs, longest_repeats_path)
    return longest_repeats_path, clear_LTR_seqs


def get_domain_sequences(longest_repeats_path, clear_LTR_seqs):
    TE_merge_classified_path = longest_repeats_path + '.final.classified'
    sample_name = 'dmel'
    TEClass_home = os.getcwd() + '/classification'
    TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
                      + ' --consensus ' + longest_repeats_path + ' --genome ' + reference \
                      + ' --thread_num ' + str(threads) + ' -o ' + tmp_output_dir
    #os.system(TEClass_command)

    query_sets = set()
    final_tmpBlastX_path = longest_repeats_path + '.tmpBlastXResults.out.bxsummary'
    with open(final_tmpBlastX_path, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            query_name = parts[0]
            # identity_parts = parts[6].split('/')
            # identity = float(identity_parts[0]) / int(identity_parts[1])
            # if identity >= 0.8:
            #     query_sets.add(query_name)
            query_sets.add(query_name)
    for query_name in clear_LTR_seqs.keys():
        query_sets.add(query_name)

    longest_repeats_contigNames, longest_repeats_contigs = read_fasta(longest_repeats_path)
    TE_domain_path = tmp_output_dir + '/TE.domain.fa'
    TE_domain_contigs = {}
    for query_name in query_sets:
        TE_domain_contigs[query_name] = longest_repeats_contigs[query_name]
    store_fasta(TE_domain_contigs, TE_domain_path)
    return TE_domain_path


def determine_struct_coding_region(longest_repeats_path, classified_path, final_tmpBlastX_path, final_tmpBlastn_path, RepLib_path, TR_out):
    # 1. get all LTR class sequences
    classified_contigNames, classified_contigs = read_fasta(classified_path)
    LTR_sequences = {}
    for name in classified_contigNames:
        class_name = name.split('#')[1]
        if class_name.__contains__('LTR'):
            LTR_sequences[name] = classified_contigs[name]

    query_classname = {}
    for name in classified_contigNames:
        parts = name.split('#')
        query_name = parts[0]
        classname = parts[1]
        query_classname[query_name] = classname

    # 2.同时考虑coding区域与结构区域。如果LTR序列具有LTR终端结构，则将终端区域和内部区域切分出来；如果LTR序列比对蛋白质有超过80%的identity，则将比对区域也取出来。
    # -----------------------------find sequences with obvious LTR structure--------------------------------------
    longest_repeats_contigNames, longest_repeats_contigs = read_fasta(longest_repeats_path)

    # LTR_query_structs = {query_name: {repeat_id: [r1, r2]} }
    LTR_query_structs = {}
    TIR_query_structs = {}
    with open(TR_out, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            repeat_id = parts[0]
            query_name = parts[2]
            query_len = int(query_name.split('-len_')[1])
            start_pos = int(parts[3])
            end_pos = int(parts[4])
            struct_info = parts[1]
            struct_parts = struct_info.split('|')
            s_type = struct_parts[0].replace('Rep' + str(repeat_id), '')
            s_len = int(struct_parts[1].replace('len=', ''))
            s_identity = float(struct_parts[2].replace('id=', ''))
            TE_len = int(struct_parts[3].replace('lenTE=', ''))
            struct_tuple = (s_type, s_len, s_identity, TE_len, start_pos, end_pos)

            if float(struct_tuple[3]) / query_len >= 0.85:
                if (struct_tuple[0] == 'LTR' or struct_tuple[0] == 'termLTR') \
                        and (struct_tuple[1] >= 100 and struct_tuple[1] <= 5000):
                    if not LTR_query_structs.__contains__(query_name):
                        LTR_query_structs[query_name] = {}
                    struct_records = LTR_query_structs[query_name]
                    if not struct_records.__contains__(repeat_id):
                        struct_records[repeat_id] = []
                    struct_tuples = struct_records[repeat_id]
                    struct_tuples.append(struct_tuple)
                elif (struct_tuple[0] == 'termTIR' or struct_tuple[0] == 'non-termTIR'):
                    if not TIR_query_structs.__contains__(query_name):
                        TIR_query_structs[query_name] = {}
                    struct_records = TIR_query_structs[query_name]
                    if not struct_records.__contains__(repeat_id):
                        struct_records[repeat_id] = []
                    struct_tuples = struct_records[repeat_id]
                    struct_tuples.append(struct_tuple)

    # choose the most probably LTR (longer in LTR len)
    best_LTR_query_structs = {}
    for query_name in LTR_query_structs.keys():
        struct_records = LTR_query_structs[query_name]
        if len(struct_records) <= 1:
            best_LTR_query_structs[query_name] = list(struct_records.values())[0]
        else:
            best_record = None
            for repeat_id in struct_records.keys():
                struct_tuples = struct_records[repeat_id]
                first_tuple = struct_tuples[0]
                second_tuple = struct_tuples[1]
                if best_record is None:
                    best_record = struct_tuples
                elif first_tuple[1] > best_record[0][1]:
                    best_record = struct_tuples
            best_LTR_query_structs[query_name] = best_record

    # choose the most probably TIR (longer in TE len)
    best_TIR_query_structs = {}
    for query_name in TIR_query_structs.keys():
        if best_LTR_query_structs.__contains__(query_name):
            continue
        struct_records = TIR_query_structs[query_name]
        query_len = int(query_name.split('-len_')[1])
        if len(struct_records) <= 1:
            best_TIR_query_structs[query_name] = list(struct_records.values())[0]
        else:
            best_record = None
            for repeat_id in struct_records.keys():
                struct_tuples = struct_records[repeat_id]
                first_tuple = struct_tuples[0]
                second_tuple = struct_tuples[1]
                if best_record is None:
                    best_record = struct_tuples
                elif first_tuple[3] > best_record[0][3]:
                    best_record = struct_tuples
            best_TIR_query_structs[query_name] = best_record

    # cut_LTR_sequences include LTR and LTR internal sequences
    keeped_sequences = {}
    query_class = {}
    # 1.align to known TEs
    rep_contigNames, rep_contigs = read_fasta(RepLib_path)
    queries = {}
    with open(final_tmpBlastn_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split(' ')
            #print(parts)
            score = int(parts[0])
            query_name = parts[4]
            q_start = int(parts[5])
            q_end = int(parts[6])
            if parts[8] == 'C':
                target_name = parts[9]
                target_start = int(parts[12])
                target_end = int(parts[11])
            else:
                target_name = parts[8]
                target_start = int(parts[9])
                target_end = int(parts[10])
            target_len = len(rep_contigs[target_name])
            if not queries.__contains__(query_name):
                queries[query_name] = []
            records = queries[query_name]
            records.append((score, q_start, q_end, target_name, target_start, target_end, target_len))

    for query_name in queries.keys():
        records = queries[query_name]
        records.sort(key=lambda x: -x[0])
        for record in records:
            if float(record[5] - record[4]) / record[6] >= 0.9:
                target_name = record[3]
                class_name = target_name.split('#')[1]
                seq = longest_repeats_contigs[query_name][record[1] - 1: record[2]]
                if keeped_sequences.__contains__(query_name):
                    query_name += '-random_' + str(random.randint(0, 10000))
                new_query_name = query_name+'-knownTEs'+'#'+class_name
                keeped_sequences[new_query_name] = seq

                # if not query_class.__contains__(new_query_name):
                #     query_class[new_query_name] = class_name
                # else:
                #     orig_class_name = query_class[new_query_name]
                #     if orig_class_name == 'confused':
                #         continue
                #     elif orig_class_name != class_name:
                #         class_name = 'confused'
                #         query_class[new_query_name] = class_name


    # 2.find LTRs
    for query_name in best_LTR_query_structs:
        if not query_classname[query_name].__contains__('LTR'):
            continue
        best_record = best_LTR_query_structs[query_name]
        first_tuple = best_record[0]
        second_tuple = best_record[1]

        leftLTR_start = first_tuple[4]
        leftLTR_end = first_tuple[5]

        ltr_len = int(first_tuple[1])
        ltr_identity = float(first_tuple[2])
        TE_len = int(first_tuple[3])

        rightLTR_start = second_tuple[4]
        rightLTR_end = second_tuple[5]
        query_seq = longest_repeats_contigs[query_name]
        left_ltr = query_seq[leftLTR_start - 1: leftLTR_end]
        right_ltr = query_seq[rightLTR_start - 1: rightLTR_end]
        internal_ltr = query_seq[leftLTR_end: rightLTR_start - 1]
        keeped_sequences[query_name+'-lLTR'+'#'+query_classname[query_name]] = left_ltr
        keeped_sequences[query_name+'-rLTR'+'#'+query_classname[query_name]] = right_ltr
        keeped_sequences[query_name+'-ILTR'+'#'+query_classname[query_name]] = internal_ltr

    for query_name in best_TIR_query_structs:
        if not query_classname[query_name].__contains__('DNA'):
            continue
        best_record = best_TIR_query_structs[query_name]
        first_tuple = best_record[0]
        second_tuple = best_record[1]

        leftTIR_start = first_tuple[4]
        leftTIR_end = first_tuple[5]

        tir_len = int(first_tuple[1])
        tir_identity = float(first_tuple[2])
        TE_len = int(first_tuple[3])

        rightTIR_start = second_tuple[4]
        rightTIR_end = second_tuple[5]
        query_seq = longest_repeats_contigs[query_name]
        left_tir = query_seq[leftTIR_start - 1: leftTIR_end]
        right_tir = query_seq[rightTIR_start - 1: rightTIR_end]
        internal_tir = query_seq[leftTIR_end: rightTIR_start - 1]
        keeped_sequences[query_name+'-TIR'+'#'+query_classname[query_name]] = left_tir+internal_tir+right_tir
        #keeped_sequences[query_name + '-TIR'] = query_seq

    # 3.如果序列比对到蛋白质有超过80%的identity，将序列取出；
    with open(final_tmpBlastX_path, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            query_name = parts[0]
            identity_parts = parts[6].split('/')
            identity = float(identity_parts[0]) / int(identity_parts[1])
            #if identity >= 0.8 and query_name not in keeped_query_names:
            if identity >= 0.8:
                seq = longest_repeats_contigs[query_name]
                keeped_sequences[query_name+'-knownProtein'+'#'+query_classname[query_name]] = seq


    # # 4.REXdb中超过80% identity, 80% cover的比对
    # protein_contigNames, protein_contigs = read_fasta(protein_db_path)
    # with open(blastx2Results_path, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.split('\t')
    #         query_name = parts[0]
    #         protein_name = parts[1]
    #         protein_len = len(protein_contigs[protein_name])
    #         protein_start = int(parts[8])
    #         protein_end = int(parts[9])
    #         identity = float(parts[2])
    #         if identity >= 0.8 and float(abs(protein_end-protein_start))/protein_len >= 0.8:
    #             query_sets.add(query_name)


    longest_repeats_path = tmp_output_dir + '/longest_repeats.confident.fa'
    store_fasta(keeped_sequences, longest_repeats_path)

    return longest_repeats_path

def run_TRsearch(TRsearch_dir, longest_repeats_multi_line_path, log):
    TRsearch_command = TRsearch_dir + '/TRsearch ' + longest_repeats_multi_line_path
    log.logger.debug(TRsearch_command)
    os.system(TRsearch_command)
    TR_out = longest_repeats_multi_line_path + '.TR.set'
    return TR_out


if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run kmerRepFinder...')
    parser.add_argument('G', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('-k', metavar='kmer size',
                        help='input kmer size, default = [ 31 ]')
    parser.add_argument('-t', metavar='thread num',
                        help='input thread num')
    parser.add_argument('a', metavar='alias name',
                        help='input alias name')
    parser.add_argument('-s', metavar='sensitive mode',
                        help='sensitive mode, default = [ 0 ]')
    parser.add_argument('--fault_tolerant_bases', metavar='fault_tolerant_bases',
                        help='the base number of fault tolerant in repeated kmers masking, default = [ 50 ]')
    parser.add_argument('-o', metavar='output dir',
                        help='output dir')
    parser.add_argument('--tandem_region_cutoff', metavar='tandem_region_cutoff',
                        help='Cutoff of the raw masked repeat regarded as tandem region, default = [ 0.8 ]')

    args = parser.parse_args()

    reference = args.G
    k_num = args.k
    threads = args.t
    alias = args.a
    output_dir = args.o
    sensitive_mode = args.s
    fault_tolerant_bases = args.fault_tolerant_bases
    tandem_region_cutoff = args.tandem_region_cutoff

    log = Logger('kmerRepFinder.log', level='debug')

    if reference is None:
        log.logger.error('\nreference path can not be empty')
        exit(-1)
    if output_dir is None:
        output_dir = os.getcwd() + '/output'
        log.logger.warning('\noutput directory path is empty, set to: ' + str(output_dir))

    if not os.path.isabs(reference):
        reference = os.path.abspath(reference)
    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(output_dir)

    if fault_tolerant_bases is not None:
        fault_tolerant_bases = int(fault_tolerant_bases)

    default_k_num = 31
    if k_num is None:
        k_num = int(default_k_num)
    else:
        k_num = int(k_num)

    default_threads = int(cpu_count())
    if threads is None:
        threads = int(default_threads)
    else:
        threads = int(threads)

    default_sensitive_mode = '0'
    if sensitive_mode is None:
        sensitive_mode = default_sensitive_mode
    is_sensitive = False
    if sensitive_mode == '1':
        is_sensitive = True

    default_fault_tolerant_bases = 100
    if fault_tolerant_bases is None:
        fault_tolerant_bases = default_fault_tolerant_bases


    default_tandem_region_cutoff = 0.5
    if tandem_region_cutoff is None:
        tandem_region_cutoff = default_tandem_region_cutoff

    flanking_len = 400
    max_single_repeat_len = 20000 # 20K bp
    extend_base_threshold = 200
    extend_base_ratio = 0.1
    skip_threshold = 500
    identity_threshold = 0.95
    length_similarity_cutoff = 0.95

    partitions_num = int(threads)

    if is_sensitive:
        alignment_tool = 'Blast'
        alignment_param = ''
        detect_mode = 'Sensitive'
    else:
        alignment_tool = 'Minimap2 + BWA'
        detect_mode = 'Normal'

    log.logger.info('\n-------------------------------------------------------------------------------------------\n'
                    'Copyright (C) 2022 Kang Hu ( kanghu@csu.edu.cn )\n'
                    'Hunan Provincial Key Lab on Bioinformatics, School of Computer Science and \n'
                    'Engineering, Central South University, Changsha 410083, P.R. China.\n'
                    '-------------------------------------------------------------------------------------------')

    log.logger.info('\nParameters configuration\n'
                    '====================================System settings========================================\n'
                    '  [Setting] The K-mer Size = [ ' + str(k_num) + 'bp]  Default( ' + str(default_k_num) + ' )\n'
                    '  [Setting] Threads = [ ' + str(threads) + ' ]  Default( ' + str(default_threads) + ' )\n'
                    '  [Setting] Detect Mode = [ ' + str(detect_mode) + ' ]  ( -s : 1 -> Sensitive, other -> Normal )\n'
                    '  [Setting] Reference sequences / assemblies path = [ ' + str(reference) + ' ]\n'
                    '  [Setting] Alias = [ ' + str(alias) + ' ]\n'
                    '  [Setting] Maximum base number of variants between repeats = [ ' + str(skip_threshold) + ' ]\n'
                    '  [Setting] Cutoff of the raw masked repeat regarded as tandem region = [ ' + str(tandem_region_cutoff) + ' ]\n'
                    '  [Setting] Output Directory = [' + str(output_dir) + ']'
                    )

    # preset steps, no needs to execute time-consuming job again
    # when the job is retried due to abnormal termination.


    # Step1. read configuration
    param_config_path = os.getcwd() + "/ParamConfig.json"
    # read param config
    with open(param_config_path, 'r') as load_f:
        param = json.load(load_f)
    load_f.close()

    chrom_seg_length = int(param['chrom_seg_length'])

    i = datetime.datetime.now()
    #tmp_output_dir = output_dir + '/CRD.' + str(i.date()) + '.' + str(i.hour) + '-' + str(i.minute) + '-' + str(i.second)
    tmp_output_dir = output_dir + '/CRD.2022-07-21.8-54-7'
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)

    total_starttime = time.time()
    tools_dir = os.getcwd() + '/tools'


    (ref_dir, ref_filename) = os.path.split(reference)
    (ref_name, ref_extension) = os.path.splitext(ref_filename)

    #os.system('cp ' + reference + ' ' + tmp_output_dir)
    reference = tmp_output_dir + '/' + ref_filename


    # --------------------------------------------------------------------------------------
    # background job: get LTR sequences from LTR_retriever for supplementary
    Genome_Tools_Home = param['Genome_Tools_Home']
    LTR_retriever_Home = param['LTR_retriever_Home']
    # run LTR_retriver background job for LTR supplementary
    # output of LTR_retriever
    # backjob = multiprocessing.Process(target=run_LTR_retriever_v1, args=(Genome_Tools_Home, LTR_retriever_Home, reference, tmp_output_dir, threads,))
    # backjob.start()


    # # # -------------------------------Stage01: this stage is used to generate kmer coverage repeats-------------------------------
    # # pipeline_starttime = time.time()
    # #
    # # starttime = time.time()
    # # # --------------------------------------------------------------------------------------
    # # # Step1. dsk get unique kmers, whose frequency >= 2
    # # freq_threshold = 2
    # # # ref_size = os.path.getsize(reference)
    # # # ref_size = ref_size / float(1024 * 1024)
    # # # if ref_size > 1024:
    # # #     #freq_threshold = 5
    # # #     log.logger.debug('warning: reference is larger than 1G, 切割reference为500M一个样本')
    # #
    # # # using multiple threads to gain speed
    # # reference_pre = convertToUpperCase_v1(reference)
    # # reference_tmp = multi_line(reference_pre, chrom_seg_length, k_num)
    # #
    # # cut_references = []
    # # cur_ref_contigs = {}
    # # cur_base_num = 0
    # # ref_index = 0
    # # with open(reference_tmp, 'r') as f_r:
    # #     for line in f_r:
    # #         line = line.replace('\n', '')
    # #         parts = line.split('\t')
    # #         ref_name = parts[0].replace('>', '')
    # #         start = parts[1]
    # #         seq = parts[2]
    # #         new_ref_name = ref_name + '$' + start
    # #         cur_ref_contigs[new_ref_name] = seq
    # #         cur_base_num += len(line)
    # #         if cur_base_num >= 500 * 1024 * 1024:
    # #             # store references
    # #             cur_ref_path = reference + '.cut'+str(ref_index)+'.fa'
    # #             store_fasta(cur_ref_contigs, cur_ref_path)
    # #             cut_references.append(cur_ref_path)
    # #             cur_ref_contigs = {}
    # #             cur_base_num = 0
    # #             ref_index += 1
    # #     if len(cur_ref_contigs) > 0:
    # #         cur_ref_path = reference + '.cut' + str(ref_index) + '.fa'
    # #         store_fasta(cur_ref_contigs, cur_ref_path)
    # #         cut_references.append(cur_ref_path)
    # #
    # # connected_repeats = {}
    # # for ref_index, cut_reference in enumerate(cut_references):
    # #     (cut_ref_dir, cut_ref_filename) = os.path.split(cut_reference)
    # #     (cut_ref_name, cut_ref_extension) = os.path.splitext(cut_ref_filename)
    # #
    # #     log.logger.debug('Start step1: get unique kmers')
    # #     dsk_h5_path = cut_ref_name + '.h5'
    # #     unique_kmer_path = tmp_output_dir + '/kmer.txt'
    # #     dsk_cmd1 = 'cd ' + ref_dir + ' && ' + tools_dir + '/dsk -file ' + cut_reference + ' -kmer-size ' + str(k_num) + ' -abundance-min ' + str(freq_threshold)
    # #     dsk_cmd2 = 'cd ' + ref_dir + ' && ' + tools_dir + '/dsk2ascii -file ' + dsk_h5_path + ' -out ' + unique_kmer_path
    # #     log.logger.debug(dsk_cmd1)
    # #     os.system(dsk_cmd1)
    # #     log.logger.debug(dsk_cmd2)
    # #     os.system(dsk_cmd2)
    # #
    # #     #tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-28.16-22-0'
    # #     # unique_kmer_path = tmp_output_dir + '/kmer.txt'
    # #     # --------------------------------------------------------------------------------------
    # #     # unique_kmer_map = {}
    # #     # with open(unique_kmer_path, 'r') as f_r:
    # #     #     for line in f_r:
    # #     #         line = line.replace('\n', '')
    # #     #         kmer = line.split(' ')[0]
    # #     #         r_kmer = getReverseSequence(kmer)
    # #     #         unique_key = kmer if kmer < r_kmer else r_kmer
    # #     #         if unique_key.__contains__('N'):
    # #     #             continue
    # #     #         unique_kmer_map[unique_key] = 1
    # #
    # #     # --------------------------------------------------------------------------------------
    # #     # Step2. each thread process a batch of kmers
    # #     log.logger.debug('Start step2: each thread process a batch of kmers')
    # #
    # #     # create thread pool, use multiple processes to execute
    # #     kmer_segments = []
    # #     with open(unique_kmer_path, 'r') as f_r:
    # #         for line in f_r:
    # #             line = line.replace('\n', '')
    # #             kmer_segments.append(line)
    # #     kmer_segments_cluster = split2cluster_normal(kmer_segments, partitions_num)
    # #
    # #     ex = ProcessPoolExecutor(partitions_num)
    # #     objs = []
    # #     for partiton_index in kmer_segments_cluster.keys():
    # #         cur_kmer_segments = kmer_segments_cluster[partiton_index]
    # #         obj = ex.submit(getUniqueKmer_v1, cur_kmer_segments, partiton_index)
    # #         objs.append(obj)
    # #     ex.shutdown(wait=True)
    # #     unique_kmer_map = {}
    # #     for obj in as_completed(objs):
    # #         for kmer in obj.result():
    # #             unique_kmer_map[kmer] = 1
    # #     log.logger.debug('generate unique_kmer_map finished')
    # #
    # #     # # store longest_repeats for testing
    # #     # unique_kmer_map_file = tmp_output_dir + '/unique_kmer_map.csv'
    # #     # with codecs.open(unique_kmer_map_file, 'w', encoding='utf-8') as f:
    # #     #     json.dump(unique_kmer_map, f)
    # #
    # #     # shared_unique_kmer_map = multiprocessing.Manager().dict()
    # #     # shared_unique_kmer_map.update(unique_kmer_map)
    # #
    # #     #reduce_partitions_num = judgeReduceThreads(unique_kmer_path, partitions_num, log)
    # #
    # #     # --------------------------------------------------------------------------------------
    # #     # Step3. get candidate repeats
    # #     log.logger.debug('Start step3: get candidate repeats from kmer coverage')
    # #     cur_connected_repeats = get_candidate_repeats(cut_reference, k_num, partitions_num, unique_kmer_map, fault_tolerant_bases, tmp_output_dir, log)
    # #
    # #     for name in cur_connected_repeats.keys():
    # #         repeat_list = cur_connected_repeats[name]
    # #         connected_repeats[name+'-'+str(ref_index)] = repeat_list
    # #
    # #
    # #     # # # single threads will ensure the accuracy
    # #     # # contigs = convertToUpperCase(reference)
    # #     # # repeat_dict, masked_ref = generate_candidate_repeats_v1(contigs, k_num, unique_kmer_map, fault_tolerant_bases)
    # #
    # # # generate repeats.fa and connected_regions
    # # repeats_path = tmp_output_dir + '/repeats.fa'
    # # node_index = 0
    # # with open(repeats_path, 'w') as f_save:
    # #     for ref_name in connected_repeats.keys():
    # #         repeat_list = connected_repeats[ref_name]
    # #         for repeat_item in repeat_list:
    # #             start_pos = repeat_item[0]
    # #             end_pos = repeat_item[1]
    # #             query_name = 'N' + str(node_index) + '-s_' + str(ref_name) + '-' + str(start_pos) + '-' + str(end_pos)
    # #             repeat = repeat_item[2]
    # #             if len(repeat) >= 80:
    # #                 f_save.write('>' + query_name + '\n' + repeat + '\n')
    # #                 node_index += 1
    # #
    # #
    # # # -------------------------------Stage02: this stage is used to do pairwise comparision, determine the repeat boundary-------------------------------
    # # blast_program_dir = param['RMBlast_Home']
    # # longest_repeats_path = determine_repeat_boundary(repeats_path, blast_program_dir)
    # #
    # #
    # # # --------------------------------------------consensus and filter tandem--------------------------------------------
    # # longest_repeats_path = tmp_output_dir + '/longest_repeats.fa'
    # # longest_repeats_consensus = tmp_output_dir + '/longest_repeats.cons.fa'
    # # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -aL ' + str(0.8) + ' -c ' + str(0.8) + ' -G 0 -g 1 -A 80 -i ' + longest_repeats_path + ' -o ' + longest_repeats_consensus + ' -T 0 -M 0'
    # # #cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -c ' + str(0.8) + ' -G 0 -g 1 -A 80 -i ' + longest_repeats_path + ' -o ' + longest_repeats_consensus + ' -T 0 -M 0'
    # # log.logger.debug(cd_hit_command)
    # # os.system(cd_hit_command)
    # #
    # # # --------------------------------------------------------------------------------------
    # # # Step6. filter low_complexity and tandem
    # # # >= tandem_region_cutoff region of the whole repeat region, then it should be filtered, since maybe false positive
    # # # keep sequences >= 80bp
    # # TRF_Path = param['TRF_Path']
    # # trf_dir = tmp_output_dir + '/trf_temp'
    # # os.system('rm -rf ' + trf_dir)
    # # if not os.path.exists(trf_dir):
    # #     os.makedirs(trf_dir)
    # # (repeat_dir, repeat_filename) = os.path.split(longest_repeats_consensus)
    # # (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)
    # # trf_command = 'cd ' + trf_dir + ' && ' + TRF_Path + ' ' + longest_repeats_consensus + ' 2 7 7 80 10 50 500 -f -d -m'
    # # log.logger.debug(trf_command)
    # # os.system(trf_command)
    # # trf_masked_repeats = trf_dir + '/' + repeat_filename + '.2.7.7.80.10.50.500.mask'
    # #
    # # trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    # # repeats_contigNames, repeats_contigs = read_fasta(longest_repeats_consensus)
    # # longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_tandem.fa'
    # # with open(longest_repeats_path, 'w') as f_save:
    # #     for name in trf_contigNames:
    # #         seq = trf_contigs[name]
    # #         if float(seq.count('N')) / len(seq) < tandem_region_cutoff and len(seq) >= 80:
    # #             f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')
    # #
    # # blast_program_dir = param['RMBlast_Home']
    # # # ---------------------------------filter segmental duplication---------------------------------------------
    # # #longest_repeats_path = tmp_output_dir + '/longest_repeats.fa'
    # # longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_tandem.fa'
    # # #longest_repeats_path = filter_duplication(longest_repeats_path, reference, blast_program_dir)
    # # #longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_duplication.fa'
    # #
    # # # ------------------------------------filter derived sequences--------------------------------------------------------
    # # # parallel
    # # #longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_duplication.fa'
    # # #longest_repeats_path = tmp_output_dir + '/TE.merge.fa'
    # # longest_repeats_path = filter_derived_seq_v1(longest_repeats_path, blast_program_dir)
    # longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_derived.fa'
    #
    # # #---------------------------------------HelitronScanner search Helitron structure--------------------------------------------------
    # # HelitronScanner_path = '/home/hukang/repeat_detect_tools/HelitronScanner/HelitronScanner.jar'
    # # TrainingSet_dir = '/home/hukang/repeat_detect_tools/TrainingSet'
    # # helitron_dir = tmp_output_dir + '/helitron'
    # # if not os.path.exists(helitron_dir):
    # #     os.makedirs(helitron_dir)
    # # # forward
    # # scanTail_command = 'java -jar ' + HelitronScanner_path + ' scanTail -lf ' + TrainingSet_dir + '/tail.lcvs -bs 0 -g ' + longest_repeats_path + ' -o ' + helitron_dir + '/tail.helitronscanner.out'
    # # scanHead_command = 'java -jar ' + HelitronScanner_path + ' scanHead -lf ' + TrainingSet_dir + '/head.lcvs -bs 0 -g ' + longest_repeats_path + ' -o ' + helitron_dir + '/head.helitronscanner.out'
    # # pairends_command = 'java -jar ' + HelitronScanner_path + ' pairends -head_score ' + helitron_dir + '/head.helitronscanner.out -tail_score ' + helitron_dir + '/tail.helitronscanner.out -output ' + helitron_dir + '/paired.helitrons'
    # # draw_command = 'java -jar ' + HelitronScanner_path + ' draw -p '+ helitron_dir + '/paired.helitrons -g ' + longest_repeats_path + ' -o ' + helitron_dir + '/draw_helitrons_pure --pure'
    # # os.system(scanTail_command)
    # # os.system(scanHead_command)
    # # os.system(pairends_command)
    # # os.system(draw_command)
    # #
    # # # reverse
    # # scanTail_command = 'java -jar ' + HelitronScanner_path + ' scanTail -lf ' + TrainingSet_dir + '/tail.lcvs -bs 0 -g ' + longest_repeats_path + ' --rc -o ' + helitron_dir + '/tail.helitronscanner.rc.out'
    # # scanHead_command = 'java -jar ' + HelitronScanner_path + ' scanHead -lf ' + TrainingSet_dir + '/head.lcvs -bs 0 -g ' + longest_repeats_path + ' --rc -o ' + helitron_dir + '/head.helitronscanner.rc.out'
    # # pairends_command = 'java -jar ' + HelitronScanner_path + ' pairends -head_score ' + helitron_dir + '/head.helitronscanner.rc.out -tail_score ' + helitron_dir + '/tail.helitronscanner.rc.out -output ' + helitron_dir + '/paired.rc.helitrons'
    # # draw_command = 'java -jar ' + HelitronScanner_path + ' draw -p ' + helitron_dir + '/paired.rc.helitrons -g ' + longest_repeats_path + ' -o ' + helitron_dir + '/draw_helitrons_pure.rc --pure'
    # # os.system(scanTail_command)
    # # os.system(scanHead_command)
    # # os.system(pairends_command)
    # # os.system(draw_command)
    # #
    # # merged_helitron_path = tmp_output_dir + '/helitron.fa'
    # # os.system('cat ' + helitron_dir + '/draw_helitrons_pure.hel.fa > ' + merged_helitron_path)
    # # os.system('cat ' + helitron_dir + '/draw_helitrons_pure.rc.hel.fa >> ' + merged_helitron_path)
    # #
    # # helitron_names, helitron_contigs = read_fasta(merged_helitron_path)
    # # # filter helitron names
    # # helitron_names_set = set()
    # # for name in helitron_names:
    # #     helitron_names_set.add(name.split('_#')[0])
    # #
    # # contigNames, contigs = read_fasta(longest_repeats_path)
    # # longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_helitron.fa'
    # # with open(longest_repeats_path, 'w') as f_save:
    # #     for name in contigNames:
    # #         if name not in helitron_names_set:
    # #             f_save.write('>'+name+'\n'+contigs[name]+'\n')
    #
    # #---------------------------------------classification and get--------------------------------------------------
    # ref_size = os.path.getsize(longest_repeats_path)
    # split_num = int(ref_size / float(0.5 * 1024 * 1024))
    # if split_num == 0:
    #     split_num = 1
    # classified_path = longest_repeats_path + '.final.classified'
    # sample_name = alias
    # TEClass_home = os.getcwd() + '/classification'
    # TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                   + ' --consensus ' + longest_repeats_path + ' --genome ' + reference \
    #                   + ' --thread_num ' + str(threads) + ' --split_num ' + str(split_num) + ' -o ' + tmp_output_dir
    # log.logger.debug(TEClass_command)
    # os.system(TEClass_command)
    # final_tmpBlastX_path = longest_repeats_path + '.tmpBlastXResults.out.bxsummary'
    # final_tmpBlastn_path = longest_repeats_path + '.tmpBlastnResults.out'
    #
    # # ---------------------------------------Terminal Repeat Search--------------------------------------------------
    # TRsearch_dir = '/home/hukang/repeat_detect_tools/REPET_linux-x64-3.0/bin'
    # contigNames, contigs = read_fasta(longest_repeats_path)
    # longest_repeats_cluster = split2cluster_normal(list(contigs.items()), partitions_num)
    #
    # longest_repeats_multi_line_dir = tmp_output_dir + '/multi_line'
    # if not os.path.exists(longest_repeats_multi_line_dir):
    #     os.makedirs(longest_repeats_multi_line_dir)
    #
    # pool = multiprocessing.Pool(processes=partitions_num)
    # for partition_index in longest_repeats_cluster.keys():
    #     cur_contigs = longest_repeats_cluster[partition_index]
    #     longest_repeats_multi_line_path = longest_repeats_multi_line_dir + '/'+str(partition_index)+'.fa'
    #
    #     outfile = open(longest_repeats_multi_line_path, 'w')  # open outfile for writing
    #     for item in cur_contigs:
    #         print_seqs(item[0], item[1], 50, outfile)
    #
    #     pool.apply_async(run_TRsearch, (TRsearch_dir, longest_repeats_multi_line_path, log,))
    # pool.close()
    # pool.join()
    #
    # # Step 3: merge final classified of each thread
    # longest_repeats_multi_line_path = tmp_output_dir + '/longest_repeats.multi_line.fa'
    # final_TR_out = longest_repeats_multi_line_path + '.TR.set'
    # if os.path.exists(final_TR_out):
    #     os.system('rm -f ' + final_TR_out)
    # for partition_index in range(partitions_num):
    #     cur_TR_out = longest_repeats_multi_line_dir + '/'+str(partition_index) + '.fa' + '.TR.set'
    #     merge_command = 'cat ' + cur_TR_out + ' >> ' + final_TR_out
    #     print(merge_command)
    #     os.system(merge_command)
    #
    #
    # # # 4.使用REXdb的蛋白质库比对比对，查找能够比对上的序列。
    # # blast_program_dir = param['RMBlast_Home']
    # # protein_db_path = os.getcwd() + '/REXdb/Viridiplantae_v3.0_ALL_protein-domains.fasta'
    # # #blastx2Results_path = get_REXdb_aligned_seq(blast_program_dir, protein_db_path, longest_repeats_path)
    # # blastx2Results_path = longest_repeats_path + '.REXdb.tmpBlastXResults.out'
    #
    # RepeatMasker_Home = '/home/hukang/repeat_detect_tools/RepeatMasker'
    # RepLib_path = RepeatMasker_Home + '/Libraries/RepeatMasker.lib'
    # # ---------------------------------------determine Structure/Coding region--------------------------------------------------
    # # longest_repeats_consensus = determine_struct_coding_region(longest_repeats_path, classified_path, final_tmpBlastX_path, TR_out, protein_db_path, blastx2Results_path)
    # longest_repeats_path = determine_struct_coding_region(longest_repeats_path, classified_path, final_tmpBlastX_path, final_tmpBlastn_path, RepLib_path, final_TR_out)
    #
    # # # add Helitron into longest_repeats_path
    # # node_index = 0
    # # with open(longest_repeats_path, 'a') as f_save:
    # #     for name in helitron_names:
    # #         seq = helitron_contigs[name]
    # #         f_save.write('>helitron_'+str(node_index)+'-len_'+str(len(seq))+'\n'+seq+'\n')
    # #         node_index += 1

    # 4. cd-hit-est 去冗余序列
    longest_repeats_path = tmp_output_dir + '/longest_repeats.confident.fa'
    longest_repeats_consensus = tmp_output_dir + '/longest_repeats.confident.cons.fa'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.95) + ' -G 0 -g 1 -A 80 -i ' + longest_repeats_path + ' -o ' + longest_repeats_consensus + ' -T 0 -M 0'
    #cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -aL ' + str(0.8) + ' -c ' + str(0.8) + ' -G 0 -g 1 -A 80 -i ' + longest_repeats_path + ' -o ' + longest_repeats_consensus + ' -T 0 -M 0'
    log.logger.debug(cd_hit_command)
    os.system(cd_hit_command)



    # # ------------------------------------use TRsearch to modify precise boundary of LTR--------------------------------------------------------
    # 1.将所有分类为LTR的序列拿出来
    # 2.根据TRsearch取出具有LTR结构的序列，并区分出LTR部分与Internal部分
    # 3.如果分类为DNA，根据TRsearch将TIR/Coding区域确定下来
    # 4.cd-hit-est 0.95生成非冗余序列



    # longest_repeats_path, clear_LTR_seqs = refine_LTR_sequences()
    # longest_repeats_path = tmp_output_dir + '/longest_repeats.refine.fa'
    #
    # # -----------------------------get TE sequences which have domain--------------------------------------
    # TE_domain_path = get_domain_sequences(longest_repeats_path, clear_LTR_seqs)
    # TE_domain_path = tmp_output_dir + '/TE.domain.fa'
    #
    # #backjob.join()
    # ltr_retriever_seq = tmp_output_dir + '/' + ref_filename + '.mod.LTRlib.fa'
    # ltr_contigNames, ltr_contigs = read_fasta(ltr_retriever_seq)
    #
    # merge_TE_path = tmp_output_dir + '/TE.merge.fa'
    # os.system('cat ' + TE_domain_path + ' > ' + merge_TE_path)
    #
    # node_index = 0
    # with open(merge_TE_path, 'a') as f_save:
    #     for name in ltr_contigNames:
    #         ltr_seq = ltr_contigs[name]
    #         f_save.write('>ltr_' + str(node_index) + '-len_' + str(len(ltr_seq)) + '\n' + ltr_seq + '\n')
    #         node_index += 1



















    # merge_TE_consensus = tmp_output_dir + '/TE.merge.cons.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.99) + ' -c ' + str(0.8) + ' -G 0 -g 1 -A 80 -i ' + merge_TE_path + ' -o ' + merge_TE_consensus + ' -T 0 -M 0'
    # log.logger.debug(cd_hit_command)
    # os.system(cd_hit_command)


    # # ------------------------------------filter false positive--------------------------------------------------------
    # longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_derived.fa'
    # # parallel alignment
    # tmp_blast_dir = tmp_output_dir + '/tmp_blast'
    # os.system('rm -rf ' + tmp_blast_dir)
    # if not os.path.exists(tmp_blast_dir):
    #     os.makedirs(tmp_blast_dir)
    # longestRepeatNames, longestRepeatContigs = read_fasta(longest_repeats_path)
    # longest_repeats_files = []
    # segments_cluster = divided_array(list(longestRepeatContigs.items()), threads)
    # for partition_index, cur_segments in enumerate(segments_cluster):
    #     single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
    #     if not os.path.exists(single_tmp_dir):
    #         os.makedirs(single_tmp_dir)
    #     split_longest_repeat_file = single_tmp_dir + '/longest_repeats_split.fa'
    #     cur_contigs = {}
    #     for item in cur_segments:
    #         cur_contigs[item[0]] = item[1]
    #     store_fasta(cur_contigs, split_longest_repeat_file)
    #     os.system('cp ' + reference + ' ' + single_tmp_dir)
    #     longest_repeats_files.append((split_longest_repeat_file, single_tmp_dir + '/' + ref_filename,
    #                                   single_tmp_dir + '/longest_repeat.ref.out'))
    #
    # ex = ProcessPoolExecutor(threads)
    # jobs = []
    # for file in longest_repeats_files:
    #     job = ex.submit(multiple_alignment, file, blast_program_dir, tools_dir)
    #     jobs.append(job)
    # ex.shutdown(wait=True)
    #
    # blastn2Results_path = tmp_output_dir + '/longest_repeat.ref.out'
    # if os.path.exists(blastn2Results_path):
    #     os.remove(blastn2Results_path)
    #
    # for job in as_completed(jobs):
    #     cur_blastn2Results_path = job.result()
    #     os.system('cat ' + cur_blastn2Results_path + ' >> ' + blastn2Results_path)



    # blastn2Results_path = tmp_output_dir + '/longest_repeat.ref.out'
    # filter_TE_path = tmp_output_dir + '/longest_repeats.cons.filter_tandem.fa'
    # repeats_contigNames, repeats_contigs = read_fasta(filter_TE_path)
    # query_records = {}
    # with open(blastn2Results_path, 'r') as f_r:
    #     for idx, line in enumerate(f_r):
    #         # print('current line idx: %d' % (idx))
    #         parts = line.split('\t')
    #         query_name = parts[0]
    #         subject_name = parts[1]
    #         identity = float(parts[2])
    #         alignment_len = int(parts[3])
    #         q_start = int(parts[6])
    #         q_end = int(parts[7])
    #         s_start = int(parts[8])
    #         s_end = int(parts[9])
    #         if identity < 80:
    #             continue
    #         if not query_records.__contains__(query_name):
    #             query_records[query_name] = []
    #         records = query_records[query_name]
    #         records.append((q_start, q_end, s_start, s_end))
    #
    # filter_fp_path = tmp_output_dir + '/repeats.filter_fp.fa'
    # with open(filter_fp_path, 'w') as f_save:
    #     for query_name in query_records.keys():
    #         records = query_records[query_name]
    #         if len(records) >= 80:
    #             f_save.write('>'+query_name+'\n'+repeats_contigs[query_name]+'\n')







    # #remove Chimerism
    # tmp_blast_dir = tmp_output_dir + '/tmp_blast'
    # os.system('rm -rf ' + tmp_blast_dir)
    # if not os.path.exists(tmp_blast_dir):
    #     os.makedirs(tmp_blast_dir)
    # blast_program_dir = param['RMBlast_Home']
    #
    # longest_repeats_path = tmp_output_dir + '/longest_repeats.fa'
    # longest_repeats_Names, longest_repeats_contigs = read_fasta(longest_repeats_path)
    # repeat_pool = longest_repeats_contigs
    #
    # max_iterate_num = 100
    # iterate_step = 200
    # cur_iterate_num = 1
    #
    # while cur_iterate_num <= max_iterate_num:
    #     pure_repeat_len = cur_iterate_num * iterate_step
    #     pure_repeat, repeat_pool = get_pure_repeat_from_pool(repeat_pool, pure_repeat_len)
    #
    #     log.logger.debug('iterate num: %d, repeat_pool size: %d' %(cur_iterate_num, len(repeat_pool)))
    #     if len(repeat_pool) <= 0:
    #         break
    #
    #     repeat_pool_files = []
    #     segments_cluster = divided_array(list(repeat_pool.items()), threads)
    #     for partition_index, cur_segments in enumerate(segments_cluster):
    #         single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
    #         if not os.path.exists(single_tmp_dir):
    #             os.makedirs(single_tmp_dir)
    #         split_repeat_pool_file = single_tmp_dir + '/repeat_pool_split.fa'
    #
    #         cur_contigs = {}
    #         for item in cur_segments:
    #             cur_contigs[item[0]] = item[1]
    #         if len(cur_contigs) <= 0:
    #             continue
    #         store_fasta(cur_contigs, split_repeat_pool_file)
    #
    #         pure_repeat_file = single_tmp_dir + '/pure_repeat.fa'
    #         store_fasta(pure_repeat, pure_repeat_file)
    #
    #         repeat_pool_files.append((pure_repeat_file, split_repeat_pool_file,
    #                                       single_tmp_dir + '/chimerism.out'))
    #
    #     ex = ProcessPoolExecutor(threads)
    #     jobs = []
    #     for partition_index, file in enumerate(repeat_pool_files):
    #         job = ex.submit(remove_chimerism, file, blast_program_dir, cur_iterate_num, partition_index)
    #         jobs.append(job)
    #     ex.shutdown(wait=True)
    #
    #     repeat_pool = {}
    #     node_index = 0
    #     for job in as_completed(jobs):
    #         cur_repeat_seqs = job.result()
    #         for seq in cur_repeat_seqs:
    #             query_name = 'N_' + str(node_index) + '-len_' + str(len(seq))
    #             repeat_pool[query_name] = seq
    #             node_index += 1
    #
    #     for name in pure_repeat.keys():
    #         seq = pure_repeat[name]
    #         query_name = 'N_' + str(node_index) + '-len_' + str(len(seq))
    #         repeat_pool[query_name] = seq
    #         node_index += 1
    #
    #     cur_iterate_num += 1
    #
    # pure_repeat_path = tmp_output_dir + '/pure_repeats.fa'
    # store_fasta(pure_repeat, pure_repeat_path)
    #
    # pure_repeats_consensus = tmp_output_dir + '/pure_repeats.cons.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -aL ' + str(0.8) + ' -c ' + str(0.8) + ' -G 0 -g 1 -A 80 -i ' + pure_repeat_path + ' -o ' + pure_repeats_consensus + ' -T 0 -M 0'
    # log.logger.debug(cd_hit_command)
    # os.system(cd_hit_command)



    # filter_TE_path = tmp_output_dir + '/longest_repeats.cons.filter_tandem.fa'
    # # --------------------------------------------------------------------------------------
    # # Step7. run TE classification to classify TE family
    # # problem may occur: RepeatClassifier: /usr/bin/perl^M: bad interpreter: No such file or directory
    # # solve: perl -i -pe 'y|\r||d' RepeatClassifier_path
    # starttime = time.time()
    # log.logger.debug('Start step8: get classified consensus sequence')
    # sample_name = alias
    # TEClass_home = os.getcwd() + '/classification'
    # TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                   + ' --consensus ' + filter_TE_path + ' --genome ' + reference \
    #                   + ' --thread_num ' + str(threads) + ' -o ' + tmp_output_dir
    # log.logger.debug(TEClass_command)
    # #os.system(TEClass_command)
    #
    # # --------------------------------------------------------------------------------------
    # # Step11. assign a family name for each classified TE consensus
    # classified_consensus_path = filter_TE_path + '.final.classified'
    # classified_contigNames, classified_contigs = read_fasta(classified_consensus_path)
    # sorted_classified_contigs = {k: v for k, v in sorted(classified_contigs.items(), key=lambda item: -len(item[1]))}
    # store_fasta(sorted_classified_contigs, classified_consensus_path)
    # with open(classified_consensus_path, 'w') as f_save:
    #     for f_id, name in enumerate(sorted_classified_contigs.keys()):
    #         sequence = sorted_classified_contigs[name]
    #         class_name = name.split('#')[1]
    #         if class_name == 'Unknown':
    #             continue
    #         f_save.write('>' + name + '\n' + sequence + '\n')


    # longest_repeats_path = tmp_output_dir + '/longest_repeats.fa'
    # blast_program_dir = param['RMBlast_Home']
    # repeats_path = tmp_output_dir + '/repeats.fa'
    #
    # # parallel
    # tmp_blast_dir = tmp_output_dir + '/tmp_blast_filter'
    # os.system('rm -rf ' + tmp_blast_dir)
    # if not os.path.exists(tmp_blast_dir):
    #     os.makedirs(tmp_blast_dir)
    #
    # (repeat_dir, repeat_filename) = os.path.split(repeats_path)
    # (repeat_name, repeat_extension) = os.path.splitext(repeats_path)
    #
    # longest_repeats_Names, longest_repeats_Contigs = read_fasta(longest_repeats_path)
    #
    # longest_repeat_files = []
    # segments_cluster = divided_array(list(longest_repeats_Contigs.items()), threads)
    # for partition_index, cur_segments in enumerate(segments_cluster):
    #     single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
    #     if not os.path.exists(single_tmp_dir):
    #         os.makedirs(single_tmp_dir)
    #     split_repeat_file = single_tmp_dir+'/repeats_split.fa'
    #     cur_contigs = {}
    #     for item in cur_segments:
    #         cur_contigs[item[0]] = item[1]
    #     store_fasta(cur_contigs, split_repeat_file)
    #     os.system('cp ' + repeats_path + ' ' + single_tmp_dir)
    #     longest_repeat_files.append((split_repeat_file, single_tmp_dir + '/' + repeat_filename,
    #                          single_tmp_dir + '/longest_repeat.pairwise.out'))
    #
    #
    # ex = ProcessPoolExecutor(threads)
    # longest_repeats_filter = {}
    # jobs = []
    # for file in longest_repeat_files:
    #     job = ex.submit(get_longest_repeats, file, blast_program_dir, extend_base_threshold)
    #     jobs.append(job)
    # ex.shutdown(wait=True)
    #
    # for job in as_completed(jobs):
    #     cur_longest_repeats = job.result()
    #     for query_name in cur_longest_repeats.keys():
    #         longest_repeats_filter[query_name] = cur_longest_repeats[query_name]
    #
    # longest_repeats_filter_path = tmp_output_dir + '/longest_repeats.filter.fa'
    # # with open(longest_repeats_filter_path, 'w') as f_save:
    # #     for query_name in longest_repeats_filter.keys():
    # #         item = longest_repeats_filter[query_name]
    # #         complete_query_num = item[5]
    # #         if complete_query_num >= 2:
    # #             longest_seq = longest_repeats_Contigs[query_name]
    # #             if len(longest_seq) >= 80:
    # #                 f_save.write('>'+query_name+'\n'+longest_seq+'\n')
    #
    # node_index = 0
    # with open(longest_repeats_filter_path, 'w') as f_save:
    #     for query_name in longest_repeats_filter.keys():
    #         item = longest_repeats_filter[query_name]
    #         start = item[0]
    #         end = item[1]
    #         confidence = item[5]
    #         longest_seq = longest_repeats_Contigs[query_name][start-1: end]
    #         if len(longest_seq) >= 80:
    #             f_save.write('>N_'+str(node_index)+'\tconfidence:'+str(confidence)+'\n'+longest_seq+'\n')
    #             node_index += 1
    #
    # longest_repeats_filter_consensus = tmp_output_dir + '/longest_repeats.filter.cons.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.95) + ' -G 0 -g 1 -A 80 -i ' + longest_repeats_filter_path + ' -o ' + longest_repeats_filter_consensus + ' -T 0 -M 0'
    # # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -c ' + str(0.95) + ' -G 0 -g 1 -A 80 -i ' + longest_repeats_path + ' -o ' + longest_repeats_consensus + ' -T 0 -M 0'
    # log.logger.debug(cd_hit_command)
    # os.system(cd_hit_command)







#     filter_TE_path = tmp_output_dir + '/repeats_connected.cons.filter_tandem.fa'
#     filter_TE_contigNames, filter_TE_contigs = read_fasta(filter_TE_path)
#     # --------------------------------------------------------------------------------------
#     # Step5. merge with LTR_retriever
#     merge_TE = tmp_output_dir + '/TE.merge.fa'
#     node_index = 0
#     with open(merge_TE, 'w') as f_save:
#         for name in filter_TE_contigNames:
#             f_save.write('>N_'+str(node_index)+'\n'+filter_TE_contigs[name]+'\n')
#             node_index += 1
#
# #backjob.join()
#     ltr_retriever_seq = tmp_output_dir + '/' + ref_filename + '.mod.LTRlib.fa'
#     ltr_contigNames, ltr_contigs = read_fasta(ltr_retriever_seq)
#     node_index = 0
#     with open(merge_TE, 'a') as f_save:
#         for name in ltr_contigNames:
#             f_save.write('>ltr_' + str(node_index) + '\n' + ltr_contigs[name] + '\n')
#             node_index += 1
#
#     merge_TE_consensus = tmp_output_dir + '/TE.merge.cons.fa'
#     cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -c ' + str(0.95) + ' -G 0 -g 1 -A 50 -i ' + merge_TE + ' -o ' + merge_TE_consensus + ' -T 0 -M 0'
#     log.logger.debug(cd_hit_command)
#     os.system(cd_hit_command)



    # # ------------------------------------ Remove nested TE insertion--------------------------------------------------
    # RepeatMasker_Home = param['RepeatMasker_Home']
    # rm_command = RepeatMasker_Home + '/RepeatMasker -pa ' + str(threads) + ' -q -no_is -norna -nolow -lib ' + merge_TE + ' ' + merge_TE
    # #os.system(rm_command)
    # rm_out = tmp_output_dir + '/nested.out'
    # os.system('cp ' + merge_TE + '.out ' + rm_out)
    # rm_out_tab = tmp_output_dir + '/nested.out.tab'
    # convert_rm2tab = 'cat ' + rm_out + ' | tr -s \' \' | sed \'s/^ *//g\' | tr \' \' \'\t\' > ' + rm_out_tab
    # os.system(convert_rm2tab)
    #
    # nested_TE_pool = {}
    # # parse out file
    # with open(rm_out_tab, 'r') as f_r:
    #     for i, line in enumerate(f_r):
    #         if i <= 2:
    #             continue
    #         parts = line.split('\t')
    #         query_name = parts[4]
    #         subject_name = parts[9]
    #         if query_name == subject_name:
    #             continue
    #         query_start = int(parts[5])
    #         query_end = int(parts[6])
    #         query_length = query_end + int(str(parts[7]).replace('(', '').replace(')', ''))
    #         direct = parts[8]
    #         if direct == '+':
    #             subject_start = int(parts[11])
    #             subject_end = int(parts[12])
    #             subject_length = subject_end + int(parts[13].replace('(', '').replace(')', ''))
    #         else:
    #             subject_start = int(parts[13])
    #             subject_end = int(parts[12])
    #             subject_length = subject_end + int(parts[11].replace('(', '').replace(')', ''))
    #         alignment_len = query_end - query_start
    #         # if query is a full query, which is included in subject, the subject is nested TE
    #         if float(alignment_len)/query_length >= 0.95 and float(alignment_len)/subject_length < 0.95:
    #             if not nested_TE_pool.__contains__(subject_name):
    #                 nested_TE_pool[subject_name] = []
    #             queries = nested_TE_pool[subject_name]
    #             queries.append((query_name, subject_start, subject_end))
    # print(nested_TE_pool)
    #
    #
    # TE_Names, TE_Contigs = read_fasta(merge_TE)
    # pure_TE_contigs = {}
    # nested_TE_contigs = {}
    # pure_TE = tmp_output_dir + '/pure_TE.fa'
    # nested_TE = tmp_output_dir + '/nested_TE.fa'
    # for name in TE_Names:
    #     nested_TE_names = nested_TE_pool.keys()
    #     if name in nested_TE_names:
    #         nested_TE_contigs[name] = TE_Contigs[name]
    #     else:
    #         pure_TE_contigs[name] = TE_Contigs[name]
    # store_fasta(pure_TE_contigs, pure_TE)
    # store_fasta(nested_TE_contigs, nested_TE)






















    # candidate_TE_fragments = tmp_output_dir + '/candidate_TE_fragments.fa'
    # node_index = 0
    # with open(candidate_TE_fragments, 'w') as f_save:
    #     for repeat_name in TE_frags.keys():
    #         frags = TE_frags[repeat_name]
    #         for frag in frags:
    #             f_save.write('>N_'+str(node_index)+'\n'+frag+'\n')
    #             node_index += 1
    #
    # candidate_TE_fragments_consensus = tmp_output_dir + '/candidate_TE_fragments.cons.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -c ' + str(0.8) + ' -A 80 -i ' + candidate_TE_fragments + ' -o ' + candidate_TE_fragments_consensus + ' -T 0 -M 0'
    # # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -c ' + str(0.8) + ' -G 0 -g 1 -A 80 -i ' + candidate_TE_fragments + ' -o ' + candidate_TE_fragments_consensus + ' -T 0 -M 0'
    # log.logger.debug(cd_hit_command)
    # #os.system(cd_hit_command)
    #

    #
    #
    # # remove complete TE freq <= 3
    # rm_output1 = tmp_output_dir + '/candidate_TE_fragments.cons.fa.out'
    # rm_output1_tab = tmp_output_dir + '/candidate_TE_fragments.cons.fa.out.tab'
    # convert_rm2tab = 'cat ' + rm_output1 + ' | tr -s \' \' | sed \'s/^ *//g\' | tr \' \' \'\t\' > ' + rm_output1_tab
    # os.system(convert_rm2tab)
    #
    # repeats_contigNames, repeats_contigs = read_fasta(candidate_TE_fragments_consensus)
    #
    # repeat_ref_pos = {}
    # with open(rm_output1_tab, 'r') as f_r:
    #     for i, line in enumerate(f_r):
    #         if i <= 2:
    #             continue
    #         parts = line.split('\t')
    #         # print(parts)
    #         chr_name = parts[4]
    #         chr_start = int(parts[5])
    #         chr_end = int(parts[6])
    #         direct = parts[8]
    #         repeat_name = parts[9]
    #         repeat_seq = repeats_contigs[repeat_name]
    #         if direct == '+':
    #             repeat_start = int(parts[11])
    #             repeat_end = int(parts[12])
    #         else:
    #             repeat_start = int(parts[13])
    #             repeat_end = int(parts[12])
    #
    #         if not repeat_ref_pos.__contains__(repeat_name):
    #             repeat_ref_pos[repeat_name] = []
    #         ref_pos = repeat_ref_pos[repeat_name]
    #         ref_pos.append((chr_name, chr_start, chr_end, repeat_start, repeat_end, len(repeat_seq)))
    #         frag_num += 1
    #
    # TE_frags = {}
    # for repeat_name in repeat_ref_pos.keys():
    #     ref_pos = repeat_ref_pos[repeat_name]
    #     complete_TE_num = 0
    #     for pos_item in ref_pos:
    #         repeat_start = pos_item[3]
    #         repeat_end = pos_item[4]
    #         repeat_len = pos_item[5]
    #
    #         chr_name = pos_item[0]
    #         chr_start = pos_item[1]
    #         chr_end = pos_item[2]
    #
    #         if float(repeat_end - repeat_start) / repeat_len >= 0.8:
    #             complete_TE_num += 1
    #
    #     if complete_TE_num > 10:
    #         # true TE fragment
    #         repeat_seq = repeats_contigs[repeat_name]
    #         TE_frags[repeat_name] = repeat_seq
    #
    # TE_frag_path = tmp_output_dir + '/TE_fragments.fa'
    # store_fasta(TE_frags, TE_frag_path)
    #
    # # --------------------------------------------------------------------------------------
    # # Step6. filter low_complexity and tandem
    # # >= tandem_region_cutoff region of the whole repeat region, then it should be filtered, since maybe false positive
    # # keep sequences >= 50bp
    # TRF_Path = param['TRF_Path']
    # trf_dir = tmp_output_dir + '/trf_temp'
    # os.system('rm -rf ' + trf_dir)
    # if not os.path.exists(trf_dir):
    #     os.makedirs(trf_dir)
    # (repeat_dir, repeat_filename) = os.path.split(TE_frag_path)
    # (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)
    # trf_command = 'cd ' + trf_dir + ' && ' + TRF_Path + ' ' + TE_frag_path + ' 2 7 7 80 10 50 500 -f -d -m'
    # log.logger.debug(trf_command)
    # os.system(trf_command)
    # trf_masked_repeats = trf_dir + '/' + repeat_filename + '.2.7.7.80.10.50.500.mask'
    #
    # trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    # repeats_contigNames, repeats_contigs = read_fasta(TE_frag_path)
    # filter_TE_path = tmp_output_dir + '/TE-filtered.fa'
    # with open(filter_TE_path, 'w') as f_save:
    #     for name in trf_contigNames:
    #         seq = trf_contigs[name]
    #         if float(seq.count('N')) / len(seq) < tandem_region_cutoff and len(seq) >= 50:
    #             f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')




    # # --------------------------------------------------------------------------------------
    # # Step4. determine repeats boundary
    # use_align_tools = 'bwa'
    # sam_path_bwa = run_alignment(repeats_path, reference, use_align_tools, threads, tools_dir)
    # cut_repeats_path = tmp_output_dir + '/repeats.cut.fa'
    # cut_repeat_v2(sam_path_bwa, repeats_path, cut_repeats_path)
    #
    # # merge redundant sequences
    # cut_repeats_consensus = tmp_output_dir + '/repeats.cut.cons.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(length_similarity_cutoff) + ' -c ' + str(identity_threshold) + ' -i ' + cut_repeats_path + ' -o ' + cut_repeats_consensus + ' -T 0 -M 0'
    # log.logger.debug(cd_hit_command)
    # os.system(cd_hit_command)
    #
    #
    # # --------------------------------------------------------------------------------------
    # # Step5. merge with LTR_retriever
    # merge_pure = tmp_output_dir + '/repeats.merge.fa'
    # merge_pure_consensus = tmp_output_dir + '/repeats.merge.consensus.fa'
    # os.system('cat ' + cut_repeats_consensus + ' > ' + merge_pure)
    # ltr_retriever_seq = tmp_output_dir + '/' + ref_filename + '.mod.LTRlib.fa'
    # backjob.join()
    # os.system('cat ' + ltr_retriever_seq + ' >> ' + merge_pure)
    # #cd_hit_command = tools_dir + '/cd-hit-est -s ' + str(length_similarity_cutoff) + ' -c ' + str(identity_threshold) + ' -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(length_similarity_cutoff) + ' -c ' + str(identity_threshold) + ' -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
    # log.logger.debug(cd_hit_command)
    # os.system(cd_hit_command)
    #
    # # --------------------------------------------------------------------------------------
    # # Step6. filter low_complexity and tandem
    # # script_path = tools_dir + '/filter-stage-1.prl'
    # # filter_repeats_path = tmp_output_dir + '/repeats-filtered.fa'
    # # filter_command = 'cat ' + merge_pure_consensus + ' | ' + script_path + ' > ' + filter_repeats_path
    # # os.system(filter_command)
    #
    # # >= tandem_region_cutoff region of the whole repeat region, then it should be filtered, since maybe false positive
    # # keep sequences >= 50bp
    # TRF_Path = param['TRF_Path']
    #
    # trf_dir = tmp_output_dir + '/trf_temp'
    # if not os.path.exists(trf_dir):
    #     os.makedirs(trf_dir)
    # (repeat_dir, repeat_filename) = os.path.split(merge_pure_consensus)
    # (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)
    # trf_command = 'cd ' + trf_dir + ' && ' + TRF_Path + ' ' + merge_pure_consensus + ' 2 7 7 80 10 50 500 -f -d -m'
    # log.logger.debug(trf_command)
    # os.system(trf_command)
    # trf_masked_repeats = trf_dir + '/' + repeat_filename + '.2.7.7.80.10.50.500.mask'
    #
    # trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    # repeats_contigNames, repeats_contigs = read_fasta(merge_pure_consensus)
    # filter_repeats_path = tmp_output_dir + '/repeats-filtered.fa'
    # with open(filter_repeats_path, 'w') as f_save:
    #     for name in trf_contigNames:
    #         seq = trf_contigs[name]
    #         if float(seq.count('N')) / len(seq) < tandem_region_cutoff and len(seq) >= 50:
    #             f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')
    #
    # # --------------------------------------------------------------------------------------
    # # Step7. run TE classification to classify TE family
    # starttime = time.time()
    # log.logger.debug('Start step8: get classified consensus sequence')
    # sample_name = alias
    # TEClass_home = os.getcwd() + '/classification'
    # TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                   + ' --consensus ' + filter_repeats_path + ' --genome ' + reference \
    #                   + ' --thread_num ' + str(threads) + ' -o ' + tmp_output_dir
    # log.logger.debug(TEClass_command)
    # os.system(TEClass_command)
    #
    # # --------------------------------------------------------------------------------------
    # # Step11. assign a family name for each classified TE consensus
    # classified_consensus_path = filter_repeats_path + '.final.classified'
    # classified_contigNames, classified_contigs = read_fasta(classified_consensus_path)
    # sorted_classified_contigs = {k: v for k, v in sorted(classified_contigs.items(), key=lambda item: -len(item[1]))}
    # family_path = tmp_output_dir + '/family_' + sample_name + '.fasta'
    # with open(family_path, 'w') as f_save:
    #     for f_id, name in enumerate(sorted_classified_contigs.keys()):
    #         sequence = sorted_classified_contigs[name]
    #         class_name = name.split('#')[1]
    #         # if len(sequence) < 80 and class_name == 'Unknown':
    #         #     continue
    #         f_save.write('>family-' + str(f_id) + '#' + class_name + '\n' + sequence + '\n')
    # endtime = time.time()
    # dtime = endtime - starttime
    # log.logger.debug("module8: get classified consensus sequence running time: %.8s s" % (dtime))
    #
    # pipeline_endtime = time.time()
    # pipeline_dtime = pipeline_endtime - pipeline_starttime
    # log.logger.debug("Total pipeline running time (no including RepeatMasker): %.8s s" % (pipeline_dtime))
    #
    # # --------------------------------------------------------------------------------------
    # # Step12. invoke RepeatMasker to align TE family to genome
    # starttime = time.time()
    # RepeatMasker_Home = param['RepeatMasker_Home']
    # RepeatMasker_output_dir = tmp_output_dir + '/' + sample_name
    # RepeatMasker_command = 'cd ' + tmp_output_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(threads) \
    #                        + ' -lib ' + family_path + ' -nolow -x -html -gff -dir ' + RepeatMasker_output_dir + ' ' + reference
    # os.system('rm -rf ' + RepeatMasker_output_dir)
    # log.logger.debug(RepeatMasker_command)
    # os.system(RepeatMasker_command)
    # endtime = time.time()
    # dtime = endtime - starttime
    # log.logger.debug("module9: invoke RepeatMasker to annotate genome running time: %.8s s" % (dtime))




