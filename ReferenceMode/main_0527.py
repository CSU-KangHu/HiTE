import argparse
import codecs
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
import pysam

#from command import run_bwa
from Util import convertToUpperCase, read_fasta, getReverseSequence, \
    Logger, split_repeats, compute_identity, run_alignment, multi_line, generate_blastlike_output, \
    get_multiple_alignment_repeat, split2cluster, cut_repeat_v1, judgeReduceThreads, get_ltr_suppl_from_ltrfinder, \
    store_fasta, printClass, parse_ref_blast_output, filter_LTR_high_similarity, get_alignment_info_v3, compare_seq, \
    getRegionCombination, getCombineFragments, convertToUpperCase_v1, generate_candidate_repeats_v2

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
        unique_key = kmer if kmer < r_kmer else r_kmer
        unique_kmers.append(unique_key)
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


def getLongestPath(grid, row, visited, skip_threshold, region_path):
    #longest_path = {}
    for i in range(row):
        for j in range(len(grid[i])):
            # start from each point in Matrix, try dfs to find a valid path
            path = {}
            dfs(grid, i, j, row, path, visited, skip_threshold)
            #print(path)
            region_path.append(path)
    #         if len(path) > len(longest_path):
    #             longest_path = path
    # print('longest path = ' + str(longest_path))


def dfs(grid, x, y, row, path, visited, skip_threshold):
    # if current node visited, return
    if visited[x][y]:
        return
    cur_node = grid[x][y]
    # add current node to path
    # if not path.__contains__(x):
    #     path[x] = []
    # p = path[x]
    # p.append(cur_node)
    path[x] = cur_node
    visited[x][y] = True
    # current node reach to the last row
    if x >= row-1:
        return
    # all cell in next row
    for j in range(len(grid[x+1])):
        # Pruning, nodes that have been visited do not need to be accessed repeatedly
        if visited[x+1][j]:
            continue
        child_node = grid[x+1][j]
        # child node is closed to current node, then there is a path between current node and child node
        if (child_node[0] - cur_node[1]) >= 0 and (child_node[0] - cur_node[1]) < skip_threshold:
            dfs(grid, x+1, j, row, path, visited, skip_threshold)


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

    default_fault_tolerant_bases = 0

    if fault_tolerant_bases is None:
        fault_tolerant_bases = default_fault_tolerant_bases


    default_tandem_region_cutoff = 0.8
    if tandem_region_cutoff is None:
        tandem_region_cutoff = default_tandem_region_cutoff


    skip_threshold = 50
    identity_threshold = 0.9
    length_similarity_cutoff = 0.9

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
    tmp_output_dir = output_dir + '/CRD.' + str(i.date()) + '.' + str(i.hour) + '-' + str(i.minute) + '-' + str(i.second)
#tmp_output_dir = output_dir + '/CRD.2022-05-12.23-7-32'
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)

    total_starttime = time.time()
    tools_dir = os.getcwd() + '/tools'

    (ref_dir, ref_filename) = os.path.split(reference)
    (ref_name, ref_extension) = os.path.splitext(ref_filename)

    starttime = time.time()
    # Step0. use RepeatMasker/trf to mask all low complexity/tandem repeats in raw repeat region
    # >= tandem_region_cutoff region of the whole repeat region, then it should be filtered, since maybe false positive
    TRF_Path = param['TRF_Path']

    trf_dir = tmp_output_dir + '/trf_temp'
    if not os.path.exists(trf_dir):
        os.makedirs(trf_dir)

    trf_command = 'cd ' + trf_dir + ' && ' + TRF_Path + ' ' + reference + ' 2 7 7 80 10 50 500 -f -d -m'
    log.logger.debug(trf_command)
    os.system(trf_command)
    trf_masked_repeats = trf_dir + '/' + ref_filename + '.2.7.7.80.10.50.500.mask'

    # trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    # repeats_contigNames, repeats_contigs = read_fasta(merge_pure_consensus)
    # repeats_path = tmp_output_dir + '/repeats.filter_tandem.fa'
    # with open(repeats_path, 'w') as f_save:
    #     for name in trf_contigNames:
    #         seq = trf_contigs[name]
    #         if float(seq.count('N')) / len(seq) < tandem_region_cutoff:
    #             f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')

    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("Step0: use trf to mask genome: %.8s s" % (dtime))

    # --------------------------------------------------------------------------------------
    # background job: get LTR sequences from LTR_retriever for supplementary
    Genome_Tools_Home = param['Genome_Tools_Home']
    LTR_retriever_Home = param['LTR_retriever_Home']
    # run LTR_retriver background job for LTR supplementary
    # output of LTR_retriever
    backjob = multiprocessing.Process(target=run_LTR_retriever_v1, args=(Genome_Tools_Home, LTR_retriever_Home, reference, tmp_output_dir, threads,))
    backjob.start()

    pipeline_starttime = time.time()

    starttime = time.time()
    # --------------------------------------------------------------------------------------
    # Step1. dsk get unique kmers, whose frequency >= 2
    freq_threshold = 2
    ref_size = os.path.getsize(reference)
    ref_size = ref_size / float(1024 * 1024)
    if ref_size > 1024:
        #freq_threshold = 5
        log.logger.debug('warning: reference is larger than 1G, increase kmer size to explicit the repeat boundary')
        k_num = 39
    log.logger.debug('Start step1: get unique kmers')
    dsk_h5_path = ref_name + '.h5'
    unique_kmer_path = tmp_output_dir + '/kmer.txt'
    dsk_cmd1 = 'cd ' + ref_dir + ' && ' + tools_dir + '/dsk -file ' + trf_masked_repeats + ' -kmer-size ' + str(k_num) + ' -abundance-min ' + str(freq_threshold)
    dsk_cmd2 = 'cd ' + ref_dir + ' && ' + tools_dir + '/dsk2ascii -file ' + dsk_h5_path + ' -out ' + unique_kmer_path
    log.logger.debug(dsk_cmd1)
    os.system(dsk_cmd1)
    log.logger.debug(dsk_cmd2)
    os.system(dsk_cmd2)

    # --------------------------------------------------------------------------------------
    # Step2. each thread process a batch of kmers
    log.logger.debug('Start step2: each thread process a batch of kmers')
    unique_kmer_map = {}
    with open(unique_kmer_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            kmer = line.split(' ')[0]
            r_kmer = getReverseSequence(kmer)
            unique_key = kmer if kmer < r_kmer else r_kmer
            if unique_key.__contains__('N'):
                continue
            unique_kmer_map[unique_key] = 1

    # --------------------------------------------------------------------------------------
    # Step3. split reference into segments
    log.logger.debug('Start step3: split reference into segments')
    reduce_partitions_num = judgeReduceThreads(unique_kmer_path, partitions_num, log)

    # using multiple threads to gain speed
    reference_pre = convertToUpperCase_v1(trf_masked_repeats)
    reference_tmp = multi_line(reference_pre, chrom_seg_length, k_num)

    segments = []
    with open(reference_tmp, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            segments.append(line)
    segments_cluster = split2cluster(segments, reduce_partitions_num)

    ex = ProcessPoolExecutor(reduce_partitions_num)
    repeat_dict = {}
    jobs = []
    for partiton_index in segments_cluster.keys():
        cur_segments = segments_cluster[partiton_index]
        job = ex.submit(generate_candidate_repeats_v2, cur_segments, k_num, unique_kmer_map, partiton_index, fault_tolerant_bases)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        cur_repeat_dict = job.result()
        for ref_name in cur_repeat_dict.keys():
            parts = ref_name.split('$')
            true_ref_name = parts[0]
            start_pos = int(parts[1])
            if not repeat_dict.__contains__(true_ref_name):
                repeat_dict[true_ref_name] = []
            new_repeat_list = repeat_dict[true_ref_name]
            cur_repeat_list = cur_repeat_dict[ref_name]
            for repeat_item in cur_repeat_list:
                new_repeat_item = (start_pos+repeat_item[0], start_pos+repeat_item[1], repeat_item[2])
                new_repeat_list.append(new_repeat_item)
    for ref_name in repeat_dict.keys():
        repeat_list = repeat_dict[ref_name]
        repeat_list.sort(key=lambda x: (x[1], x[2]))

    # # single threads will ensure the accuracy
    # contigs = convertToUpperCase(reference)
    # repeat_dict, masked_ref = generate_candidate_repeats_v1(contigs, k_num, unique_kmer_map, fault_tolerant_bases)

    # store repeat_dict for testing
    repeat_dict_file = tmp_output_dir + '/repeat_dict.csv'
    with codecs.open(repeat_dict_file, 'w', encoding='utf-8') as f:
        json.dump(repeat_dict, f)

    # new strategy by Kang Hu 2022/05/24
    # Step1: generate repeats.fa and connected_regions

    # connected_regions = {ref_name: {region_id: [(f1, start1, end1), (f2, start2, end2), (f3, start3, end3)], [(f4, start4, end4), (f5, start5, end5), (f6, start6, end6)]}}
    connected_regions = {}
    repeats_path = tmp_output_dir + '/repeats.fa'
    node_index = 0
    region_index = 0
    with open(repeats_path, 'w') as f_save:
        for ref_name in repeat_dict.keys():
            repeat_list = repeat_dict[ref_name]
            if not connected_regions.__contains__(ref_name):
                connected_regions[ref_name] = {}
            regions = connected_regions[ref_name]
            last_start_pos = -1
            last_end_pos = -1
            for repeat_item in repeat_list:
                start_pos = repeat_item[0]
                end_pos = repeat_item[1]
                query_name = 'N' + str(node_index) + '-s_' + str(ref_name) + '-' + str(start_pos) + '-' + str(end_pos)
                repeat = repeat_item[2]
                f_save.write('>' + query_name + '\n' + repeat + '\n')
                node_index += 1
                # generate connected_regions
                if last_start_pos == -1:
                    regions[region_index] = [(query_name, start_pos, end_pos)]
                else:
                    if (start_pos - last_end_pos) < skip_threshold:
                        # close to current region
                        cur_region = regions[region_index]
                        cur_region.append((query_name, start_pos, end_pos))
                        regions[region_index] = cur_region
                    else:
                        # far from current region, start a new region
                        region_index += 1
                        cur_region = []
                        cur_region.append((query_name, start_pos, end_pos))
                        regions[region_index] = cur_region
                last_start_pos = start_pos
                last_end_pos = end_pos
            connected_regions[ref_name] = regions

    # store connected_regions for testing
    connected_regions_file = tmp_output_dir + '/connected_regions.csv'
    with codecs.open(connected_regions_file, 'w', encoding='utf-8') as f:
        json.dump(connected_regions, f)

    # Step2: align repeat back to reference, generate frag_pos_dict
    repeat_contignames, repeat_contigs = read_fasta(repeats_path)
    use_align_tools = 'bwa'
    sam_path_bwa = run_alignment(repeats_path, trf_masked_repeats, use_align_tools, threads, tools_dir)
# sam_path_bwa = tmp_output_dir + '/repeats.sam'
    query_records = {}
    samfile = pysam.AlignmentFile(sam_path_bwa, "rb")
    for read in samfile.fetch():
        if read.is_unmapped:
            continue
        query_name = read.query_name
        reference_name = read.reference_name
        cigar = read.cigartuples
        cigarstr = read.cigarstring
        NM_tag = 0
        try:
            NM_tag = read.get_tag('NM')
        except KeyError:
            NM_tag = -1
        identity = compute_identity(cigarstr, NM_tag, 'BLAST')
        identity = float(identity) * 100
        is_reverse = read.is_reverse
        alignment_len = read.query_alignment_length
        # pos start from 1, change to 0
        t_start = int(read.reference_start)-1
        t_end = int(read.reference_end)-1
        if t_start > t_end:
            tmp = t_start
            t_start = t_end
            t_end = tmp

        if not query_records.__contains__(query_name):
            query_records[query_name] = []
        records = query_records[query_name]
        records.append((reference_name, alignment_len, identity, t_start, t_end))
        query_records[query_name] = records

    # frag_pos_dict = {f1: {ref1: [(start1, end1), (start2, end2), (start3, end3)]}}
    frag_pos_dict = {}
    for query_name in query_records.keys():
        complete_alignment_num = 0
        query_seq = repeat_contigs[query_name]
        query_len = len(query_seq)
        if not frag_pos_dict.__contains__(query_name):
            frag_pos_dict[query_name] = {}
        ref_pos = frag_pos_dict[query_name]
        # get fragments original position
        pos_parts = query_name.split('-s_')[1].split('-')
        original_ref = pos_parts[0]
        original_start = int(pos_parts[1])
        original_end = int(pos_parts[2])
        for i, record in enumerate(query_records[query_name]):
            reference_name = record[0]
            alignment_len = record[1]
            identity = record[2]
            t_start = record[3]
            t_end = record[4]
            if float(alignment_len) / query_len >= 0.95 and identity >= 95:
                complete_alignment_num += 1
            if not ref_pos.__contains__(reference_name):
                ref_pos[reference_name] = []
            same_chr_pos = ref_pos[reference_name]
            # alignment from original parts
            if original_ref == reference_name and abs(t_start - original_start) < 10 and abs(
                    t_end - original_end) < 10:
                continue
            same_chr_pos.append((t_start, t_end))
            ref_pos[reference_name] = same_chr_pos
        # sort pos in each ref_name
        for reference_name in ref_pos.keys():
            same_chr_pos = ref_pos[reference_name]
            same_chr_pos.sort(key=lambda x: (x[0], x[1]))
            ref_pos[reference_name] = same_chr_pos
        frag_pos_dict[query_name] = ref_pos

    # store frag_pos_dict for testing
    frag_pos_dict_file = tmp_output_dir + '/frag_pos_dict.csv'
    with codecs.open(frag_pos_dict_file, 'w', encoding='utf-8') as f:
        json.dump(frag_pos_dict, f)

    # Step3: generate pathMatrix
    # pathMatrix = {region_id: {ref_name: [[(start1, end1), (start2, end2), (start3, end3)], [(start1, end1), (start2, end2), (start3, end3)]]}}
    pathMatrix = {}
    for ref_name in connected_regions.keys():
        regions = connected_regions[ref_name]
        for region_index in regions.keys():
            if not pathMatrix.__contains__(region_index):
                pathMatrix[region_index] = {}
            cur_region_matrixs = pathMatrix[region_index]

            # get one region
            cur_region = regions[region_index]
            # get all ref_names for one region
            cur_ref_names_union = set()
            for i in range(len(cur_region)):
                frag_name = cur_region[i][0]
                if not frag_pos_dict.__contains__(frag_name):
                    frag_pos_dict[frag_name] = {}
                ref_pos = frag_pos_dict[frag_name]
                for ref in ref_pos.keys():
                    cur_ref_names_union.add(ref)

            for ref in cur_ref_names_union:
                if not cur_region_matrixs.__contains__(ref):
                    cur_region_matrixs[ref] = []
                cur_ref_matrix = cur_region_matrixs[ref]
                for i in range(len(cur_region)):
                    frag_name = cur_region[i][0]
                    if not frag_pos_dict.__contains__(frag_name):
                        frag_pos_dict[frag_name] = {}
                    ref_pos = frag_pos_dict[frag_name]
                    if not ref_pos.__contains__(ref):
                        ref_pos[ref] = []
                    same_chr_pos = ref_pos[ref]
                    cur_ref_matrix.append(same_chr_pos)
                cur_region_matrixs[ref] = cur_ref_matrix
            pathMatrix[region_index] = cur_region_matrixs

    # store pathMatrix for testing
    pathMatrix_file = tmp_output_dir + '/pathMatrix.csv'
    with codecs.open(pathMatrix_file, 'w', encoding='utf-8') as f:
        json.dump(pathMatrix, f)

    # Step1: According to fragment position matrix, get all valid paths
    region_paths = {}
    # go through each Matrix, compute the valid path
    for region_index in pathMatrix:
        cur_region_matrixs = pathMatrix[region_index]
        if not region_paths.__contains__(region_index):
            region_paths[region_index] = []
        region_path = region_paths[region_index]
        for ref in cur_region_matrixs.keys():
            cur_ref_matrix = cur_region_matrixs[ref]
            row = len(cur_ref_matrix)
            # define visited matrix, avoid duplicate visits
            visited = [[False for j in range(len(cur_ref_matrix[i]))] for i in range(row)]
            # print('region id=%s, ref=%s' %(region_index, ref))
            # core
            getLongestPath(cur_ref_matrix, row, visited, skip_threshold, region_path)
        region_path.sort(key=lambda x: (len(x)))
        region_paths[region_index] = region_path

    # store region_paths for testing
    region_paths_file = tmp_output_dir + '/region_paths.csv'
    with codecs.open(region_paths_file, 'w', encoding='utf-8') as f:
        json.dump(region_paths, f)

    # Step2: record all fragments in each region, including start and end position
    # generate region fragments
    # region_fragments = {region_id: {0: (f1,s1,e1), 1: (f2,s2,e2)}}
    region_fragments = {}
    # region_fragments1 = {region_id: {f1: (s1,e1), f2: (s2,e2)}}
    region_fragments1 = {}
    for ref in connected_regions.keys():
        region_dict = connected_regions[ref]
        for region_index in region_dict.keys():
            if not region_fragments.__contains__(region_index):
                region_fragments[region_index] = {}
            region_fragment = region_fragments[region_index]

            if not region_fragments1.__contains__(region_index):
                region_fragments1[region_index] = {}
            region_fragment1 = region_fragments1[region_index]

            region_frags = region_dict[region_index]
            for index, frag_item in enumerate(region_frags):
                region_fragment[index] = frag_item

            for index, frag_item in enumerate(region_frags):
                frag_name = frag_item[0]
                region_fragment1[frag_name] = (frag_item[1], frag_item[2])

    # store region_fragments for testing
    region_fragments_file = tmp_output_dir + '/region_fragments.csv'
    with codecs.open(region_fragments_file, 'w', encoding='utf-8') as f:
        json.dump(region_fragments, f)

    # store region_fragments1 for testing
    region_fragments1_file = tmp_output_dir + '/region_fragments1.csv'
    with codecs.open(region_fragments1_file, 'w', encoding='utf-8') as f:
        json.dump(region_fragments1, f)

    # Step3: According to valid paths, the fragments in the region are connected across gaps.
    # Fragments in region are connected according to the longer path in valid paths.

    # keeped_path = {region_id: [f1f2f3, f4, f5f6]}
    keeped_paths = {}
    for region_index in region_paths.keys():
        region_fragment = region_fragments[region_index]
        for path in region_paths[region_index]:
            isPathExist = True
            pathName = ''
            for i, frag_index in enumerate(path.keys()):
                if region_fragment.__contains__(frag_index):
                    frag_item = region_fragment[frag_index]
                    frag_name = frag_item[0]
                    if i == 0:
                        pathName += frag_name
                    else:
                        pathName += ',' + frag_name
                else:
                    isPathExist = False
                    break

            if isPathExist:
                if not keeped_paths.__contains__(region_index):
                    keeped_paths[region_index] = []
                paths = keeped_paths[region_index]
                paths.append(pathName)
                for frag_index in path.keys():
                    del region_fragment[frag_index]

    # store keeped_paths for testing
    keeped_paths_file = tmp_output_dir + '/keeped_paths.csv'
    with codecs.open(keeped_paths_file, 'w', encoding='utf-8') as f:
        json.dump(keeped_paths, f)

    # connected_frags = {region_id: [(f1,s1,e1), (f2f3f4,s2,e2), (f5,s3,e3)]}
    connected_frags = {}
    for region_index in keeped_paths.keys():
        keeped_frag_name = []
        paths = keeped_paths[region_index]
        region_fragment1 = region_fragments1[region_index]
        # connect path
        for path in paths:
            if len(path) <= 0:
                continue
            last_connected_frag_start = -1
            last_connected_frag_end = -1
            for i, frag_name in enumerate(path.split(',')):
                keeped_frag_name.append(frag_name)
                frag_item = region_fragment1[frag_name]
                if last_connected_frag_start == -1:
                    last_connected_frag_start = frag_item[0]
                last_connected_frag_end = frag_item[1]

            if not connected_frags.__contains__(region_index):
                connected_frags[region_index] = []
            connected_frag = connected_frags[region_index]
            connected_frag.append((path, last_connected_frag_start, last_connected_frag_end))

        for frag_name in region_fragments1[region_index].keys():
            frag_item = region_fragments1[region_index][frag_name]
            if frag_name not in keeped_frag_name:
                if not connected_frags.__contains__(region_index):
                    connected_frags[region_index] = []
                connected_frag = connected_frags[region_index]
                connected_frag.append((frag_name, frag_item[0], frag_item[1]))

        connected_frag = connected_frags[region_index]
        connected_frag.sort(key=lambda x: (x[1], x[2]))

    # store connected_frags for testing
    connected_frags_file = tmp_output_dir + '/connected_frags.csv'
    with codecs.open(connected_frags_file, 'w', encoding='utf-8') as f:
        json.dump(connected_frags, f)

    refNames, refContigs = read_fasta(trf_masked_repeats)
    repeats_connected_file = tmp_output_dir + '/repeats_connected.fa'
    repeats_connected = {}
    index = 0
    for region_index in connected_frags.keys():
        for connected_frag in connected_frags[region_index]:
            query_name = 'R' + str(index)
            frag_name = connected_frag[0].split(',')[0]
            ref_name = frag_name.split('-s_')[1].split('-')[0]
            seq = refContigs[ref_name][connected_frag[1]: connected_frag[2] + 1]
            index += 1
            repeats_connected[query_name] = seq
    sorted_repeats_connected = {k: v for k, v in sorted(repeats_connected.items(), key=lambda item: -len(item[1]))}
    store_fasta(sorted_repeats_connected, repeats_connected_file)


    # masked_path = tmp_output_dir + '/ref.masked.fa'
    # ref_index = 0
    # with open(masked_path, 'w') as f_save:
    #     for masked_sequence in masked_ref:
    #         f_save.write('>ref' + str(ref_index) + '\n' + "".join(masked_sequence) + '\n')
    #         ref_index += 1

    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("module1: Generate repeat file running time: %.8s s" % (dtime))


    # param_config_path = os.getcwd() + "/ParamConfig.json"
    # # read param config
    # with open(param_config_path, 'r') as load_f:
    #     param = json.load(load_f)
    # load_f.close()
    #
    # log = Logger('kmerRepFinder.log', level='debug')
    # tools_dir = os.getcwd() + '/tools'
    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-24.9-36-59'
    # output_dir = tmp_output_dir
    # reference = '/public/home/hpc194701009/Ref/dmel-all-chromosome-r5.43.fasta'
    # (ref_dir, ref_filename) = os.path.split(reference)
    # (ref_name, ref_extension) = os.path.splitext(ref_filename)
    # threads = 48
    # tandem_region_cutoff = 0.8
    # alias = 'dmel'
    # repeats_consensus = tmp_output_dir + '/repeats.consensus.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -s ' + str(length_similarity_cutoff) + ' -c ' + str(identity_threshold) + ' -i ' + repeats_path + ' -o ' + repeats_consensus + ' -T 0 -M 0'
    # log.logger.debug(cd_hit_command)
    # os.system(cd_hit_command)

    # newly strategy: 2022-04-29 by Kang Hu
    # 01: use bwa to get single mapped sequence
    candidate_repeats_path = repeats_connected_file
    blast_program_dir = param['RMBlast_Home']
    use_align_tools = 'bwa'
    sam_path_bwa = run_alignment(candidate_repeats_path, trf_masked_repeats, use_align_tools, threads, tools_dir)
    sam_paths = []
    sam_paths.append(sam_path_bwa)
    # unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds = get_alignment_info(sam_paths)
    new_mapping_repeatIds, query_position = get_alignment_info_v3(sam_paths, candidate_repeats_path)
    repeat_freq_path = tmp_output_dir + '/repeats.freq.fa'
    repeat_multiple_path = tmp_output_dir + '/repeats.multiple.fa'
    merge_repeat_contigNames, merge_repeat_contigs = read_fasta(candidate_repeats_path)
    # single repeat probably be Chimeric repeat
    with open(repeat_freq_path, 'w') as f_save:
        for repeat_id in new_mapping_repeatIds.keys():
            freq = new_mapping_repeatIds[repeat_id][0]
            seq = merge_repeat_contigs[repeat_id]
            if freq <= 1 or (freq < 5 and len(seq) < 80):
                continue
            f_save.write('>' + repeat_id + '\tcopies=' + str(freq) + '\n' + seq + '\n')

    with open(repeat_multiple_path, 'w') as f_save:
        for repeat_id in new_mapping_repeatIds.keys():
            freq = new_mapping_repeatIds[repeat_id][0]
            seq = merge_repeat_contigs[repeat_id]
            if freq <= 1:
                continue
            f_save.write('>' + repeat_id + '\n' + seq + '\n')

    # 06: merge
    merge_pure = tmp_output_dir + '/repeats.merge.pure.fa'
    merge_pure_consensus = tmp_output_dir + '/repeats.merge.pure.consensus.fa'
    os.system('cat ' + repeat_multiple_path + ' > ' + merge_pure)
    ltr_retriever_seq = tmp_output_dir + '/' + ref_filename + '.mod.LTRlib.fa'
    backjob.join()
    os.system('cat ' + ltr_retriever_seq + ' >> ' + merge_pure)
    #cd_hit_command = tools_dir + '/cd-hit-est -s ' + str(length_similarity_cutoff) + ' -c ' + str(identity_threshold) + ' -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -c ' + str(0.95) + ' -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
    log.logger.debug(cd_hit_command)
    os.system(cd_hit_command)

    # --------------------------------------------------------------------------------------
    # Step10. run TE classification to classify TE family
    starttime = time.time()
    log.logger.debug('Start step8: get classified consensus sequence')
    sample_name = alias
    TEClass_home = os.getcwd() + '/classification'
    TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
                      + ' --consensus ' + merge_pure_consensus + ' --genome ' + trf_masked_repeats \
                      + ' --thread_num ' + str(threads) + ' -o ' + tmp_output_dir
    log.logger.debug(TEClass_command)
    os.system(TEClass_command)

    # --------------------------------------------------------------------------------------
    # Step11. assign a family name for each classified TE consensus
    classified_consensus_path = merge_pure_consensus + '.final.classified'
    classified_contigNames, classified_contigs = read_fasta(classified_consensus_path)
    sorted_classified_contigs = {k: v for k, v in sorted(classified_contigs.items(), key=lambda item: -len(item[1]))}
    family_path = tmp_output_dir + '/family_' + sample_name + '.fasta'
    with open(family_path, 'w') as f_save:
        for f_id, name in enumerate(sorted_classified_contigs.keys()):
            sequence = sorted_classified_contigs[name]
            class_name = name.split('#')[1]
            if len(sequence) < 80 and class_name == 'Unknown':
                continue
            f_save.write('>family-' + str(f_id) + '#' + class_name + '\n' + sequence + '\n')
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("module8: get classified consensus sequence running time: %.8s s" % (dtime))

    pipeline_endtime = time.time()
    pipeline_dtime = pipeline_endtime - pipeline_starttime
    log.logger.debug("Total pipeline running time (no including RepeatMasker): %.8s s" % (pipeline_dtime))

    # --------------------------------------------------------------------------------------
    # Step12. invoke RepeatMasker to align TE family to genome
    starttime = time.time()
    RepeatMasker_Home = param['RepeatMasker_Home']
    RepeatMasker_output_dir = tmp_output_dir + '/' + sample_name
    RepeatMasker_command = 'cd ' + tmp_output_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(threads) \
                           + ' -lib ' + family_path + ' -nolow -x -html -gff -dir ' + RepeatMasker_output_dir + ' ' + reference
    os.system('rm -rf ' + RepeatMasker_output_dir)
    log.logger.debug(RepeatMasker_command)
    os.system(RepeatMasker_command)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("module9: invoke RepeatMasker to annotate genome running time: %.8s s" % (dtime))




