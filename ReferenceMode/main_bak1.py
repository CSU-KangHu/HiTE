import argparse
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
    store_fasta, printClass, parse_ref_blast_output, filter_LTR_high_similarity, get_alignment_info_v1, generate_candidate_repeats_v2, convertToUpperCase_v1

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
    #log.logger.debug('partition %d process: %d segments' %(partiton_index, len(cur_segments)))
    cur_lines = []
    cur_masked_segments = []
    for line in cur_segments:
        # parts = line.split("\t")
        # contigName = parts[0].split(" ")[0]
        # # remove ">"
        # start = int(parts[1])
        # line = parts[2]
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
            # kmer_num = three_bit_encode(kmer)
            # r_kmer_num = three_bit_encode(r_kmer)
            # unique_key = kmer_num if kmer_num < r_kmer_num else r_kmer_num
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
    #log.logger.debug('partition %d process finished' % (partiton_index))
    return repeat_list

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
    parser.add_argument('--min_ltr_complete_len', metavar='min_ltr_complete_len',
                        help='Minimum complete LTR length, default = [ 400 ]')
    parser.add_argument('--max_ltr_complete_len', metavar='max_ltr_complete_len',
                        help='Maximum complete LTR length, default = [ 22000 ]')
    parser.add_argument('--min_ltr_direct_repeat_len', metavar='min_ltr_direct_repeat_len',
                        help='Minimum LTR direct repeat length, default = [ 100 ]')
    parser.add_argument('--max_ltr_direct_repeat_len', metavar='max_ltr_direct_repeat_len',
                        help='Maximum LTR direct repeat length, default = [ 6000 ]')
    parser.add_argument('--min_tir_complete_len', metavar='min_tir_complete_len',
                        help='Minimum complete TIR length, default = [ 1000 ]')
    parser.add_argument('--max_tir_complete_len', metavar='max_tir_complete_len',
                        help='Maximum complete TIR length, default = [ 40000 ]')
    parser.add_argument('--min_tir_direct_repeat_len', metavar='min_tir_direct_repeat_len',
                        help='Minimum TIR direct repeat length, default = [ 0 ]')
    parser.add_argument('--max_tir_direct_repeat_len', metavar='max_tir_direct_repeat_len',
                        help='Maximum TIR direct repeat length, default = [ 1000 ]')
    parser.add_argument('--long_repeat_threshold', metavar='long_repeat_threshold',
                        help='Threshold of long repeat, default = [ 2000 ]')
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
    min_ltr_complete_len = args.min_ltr_complete_len
    max_ltr_complete_len = args.max_ltr_complete_len
    min_ltr_direct_repeat_len = args.min_ltr_direct_repeat_len
    max_ltr_direct_repeat_len = args.max_ltr_direct_repeat_len
    min_tir_complete_len = args.min_tir_complete_len
    max_tir_complete_len = args.max_tir_complete_len
    min_tir_direct_repeat_len = args.min_tir_direct_repeat_len
    max_tir_direct_repeat_len = args.max_tir_direct_repeat_len
    long_repeat_threshold = args.long_repeat_threshold
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

    default_fault_tolerant_bases = 10

    if fault_tolerant_bases is None:
        fault_tolerant_bases = default_fault_tolerant_bases

    # params get from LTRDetector
    default_min_ltr_complete_len = 400
    default_max_ltr_complete_len = 22000
    default_min_ltr_direct_repeat_len = 100
    default_max_ltr_direct_repeat_len = 6000

    # Class II transposons range in length from 1,000 to as many as 40,000 base pairs.
    # https://www.britannica.com/science/transposon
    default_min_tir_complete_len = 1000
    default_max_tir_complete_len = 40000
    default_min_tir_direct_repeat_len = 0
    default_max_tir_direct_repeat_len = 1000

    default_long_repeat_threshold = 2000

    if min_ltr_complete_len is None:
        min_ltr_complete_len = default_min_ltr_complete_len
    if max_ltr_complete_len is None:
        max_ltr_complete_len = default_max_ltr_complete_len
    if min_ltr_direct_repeat_len is None:
        min_ltr_direct_repeat_len = default_min_ltr_direct_repeat_len
    if max_ltr_direct_repeat_len is None:
        max_ltr_direct_repeat_len = default_max_ltr_direct_repeat_len
    if min_tir_complete_len is None:
        min_tir_complete_len = default_min_tir_complete_len
    if max_tir_complete_len is None:
        max_tir_complete_len = default_max_tir_complete_len
    if min_tir_direct_repeat_len is None:
        min_tir_direct_repeat_len = default_min_tir_direct_repeat_len
    if max_tir_direct_repeat_len is None:
        max_tir_direct_repeat_len = default_max_tir_direct_repeat_len
    if long_repeat_threshold is None:
        long_repeat_threshold = default_long_repeat_threshold
    if min_ltr_complete_len is None:
        min_ltr_complete_len = default_min_ltr_complete_len

    default_tandem_region_cutoff = 0.5
    if tandem_region_cutoff is None:
        tandem_region_cutoff = default_tandem_region_cutoff


    min_TE_len = 80

    domain_min_identity = 80
    domain_match_ratio = 0.8


    partitions_num = int(threads)

    if is_sensitive:
        alignment_tool = 'Blast'
        alignment_param = ''
        detect_mode = 'Sensitive'
    else:
        alignment_tool = 'Minimap2 + BWA'
        alignment_param = '\n  [Setting] Long Reads Threshold = [ ' + str(long_repeat_threshold) + ' bp]  Default( ' + str(default_long_repeat_threshold) + ' )'
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
                    '  [Setting] Base number of fault tolerant in repeated kmers masking = [ ' + str(fault_tolerant_bases) + ' ]\n'
                    '  [Setting] Cutoff of the raw masked repeat regarded as tandem region = [ ' + str(tandem_region_cutoff) + ' ]\n'
                    '  [Setting] Output Directory = [' + str(output_dir) + ']'
                    )

    log.logger.info('\n====================================Alignment settings========================================\n'
                    '  [Setting] Alignment Tool = [ ' + str(alignment_tool) + ' ]' + alignment_param)

    log.logger.info('\n====================================LTR settings========================================\n'
                    '  [Setting] Minimum complete LTR length  = [ ' + str(min_ltr_complete_len) + 'bp]  Default( ' + str(default_min_ltr_complete_len) + ' )\n'
                    '  [Setting] Maximum complete LTR length = [ ' + str(max_ltr_complete_len) + ' bp]  Default( ' + str(default_max_ltr_complete_len) + ' )\n'
                    '  [Setting] Minimum LTR direct repeat length = [ ' + str(min_ltr_direct_repeat_len) + ' bp]  Default( ' + str(default_min_ltr_direct_repeat_len) + ' )\n'
                    '  [Setting] Maximum LTR direct repeat length = [ ' + str(max_ltr_direct_repeat_len) + ' bp]  Default( ' + str(default_max_ltr_direct_repeat_len) + ' )'
                    )

    log.logger.info('\n====================================TIR settings========================================\n'
                    '  [Setting] Minimum complete TIR length  = [ ' + str(min_tir_complete_len) + 'bp]  Default( ' + str(default_min_tir_complete_len) + ' )\n'
                    '  [Setting] Maximum complete TIR length = [ ' + str(max_tir_complete_len) + ' bp]  Default( ' + str(default_max_tir_complete_len) + ' )\n'
                    '  [Setting] Minimum TIR direct repeat length = [ ' + str(min_tir_direct_repeat_len) + ' bp]  Default( ' + str(default_min_tir_direct_repeat_len) + ' )\n'
                    '  [Setting] Maximum TIR direct repeat length = [ ' + str(max_tir_direct_repeat_len) + ' bp]  Default( ' + str(default_max_tir_direct_repeat_len) + ' )'
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
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)

    total_starttime = time.time()
    tools_dir = os.getcwd() + '/tools'

    (ref_dir, ref_filename) = os.path.split(reference)
    (ref_name, ref_extension) = os.path.splitext(ref_filename)

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
        freq_threshold = 5
    log.logger.debug('Start step1: get unique kmers')
    log.logger.debug('warning: reference is larger than 1G, increase the frequency of unique kmers to gain speed')
    dsk_h5_path = ref_name + '.h5'
    unique_kmer_path = tmp_output_dir + '/kmer.txt'
    dsk_cmd1 = 'cd ' + ref_dir + ' && ' + tools_dir + '/dsk -file ' + reference + ' -kmer-size ' + str(k_num) + ' -abundance-min ' + str(freq_threshold)
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
    reference_pre = convertToUpperCase_v1(reference)
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
        job = ex.submit(generate_candidate_repeats_v2, cur_segments, k_num, unique_kmer_map, partiton_index,
                        fault_tolerant_bases)
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
                new_repeat_item = (start_pos + repeat_item[0], start_pos + repeat_item[1], repeat_item[2])
                new_repeat_list.append(new_repeat_item)
    for ref_name in repeat_dict.keys():
        repeat_list = repeat_dict[ref_name]
        repeat_list.sort(key=lambda x: (x[1], x[2]))

    repeats_path = tmp_output_dir + '/repeats.fa'
    node_index = 0
    with open(repeats_path, 'w') as f_save:
        for ref_name in repeat_dict.keys():
            repeat_list = repeat_dict[ref_name]
            for repeat_item in repeat_list:
                start_pos = repeat_item[0]
                end_pos = repeat_item[1]
                query_name = 'N' + str(node_index) + '-s_' + str(ref_name) + '-' + str(start_pos) + '-' + str(end_pos)
                repeat = repeat_item[2]
                f_save.write('>' + query_name + '\n' + repeat + '\n')
                node_index += 1

    # cur_segments = convertToUpperCase(reference)
    # repeat_segments = generate_candidate_repeats(cur_segments, k_num, unique_kmer_map, fault_tolerant_bases)
    #
    # repeats_path = tmp_output_dir + '/repeats.fa'
    # node_index = 0
    # with open(repeats_path, 'w') as f_save:
    #     for repeat in repeat_segments:
    #         f_save.write('>Node_'+str(node_index)+'-len_'+str(len(repeat))+'\n'+repeat+'\n')
    #         node_index += 1

    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("module1: Generate repeat file running time: %.8s s" % (dtime))

    starttime = time.time()
    # Step0. use RepeatMasker/trf to mask all low complexity/tandem repeats in raw repeat region
    # >= tandem_region_cutoff region of the whole repeat region, then it should be filtered, since maybe false positive
    TRF_Path = param['TRF_Path']
    RepeatMasker_Home = param['RepeatMasker_Home']
    trf_dir = tmp_output_dir + '/trf_temp'
    if not os.path.exists(trf_dir):
        os.makedirs(trf_dir)
    (repeat_dir, repeat_filename) = os.path.split(repeats_path)
    (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)
    trf_command = 'cd ' + trf_dir + ' && ' + TRF_Path + ' ' + repeats_path + ' 2 7 7 80 10 50 500 -f -d -m'
    log.logger.debug(trf_command)
    os.system(trf_command)
    trf_masked_repeats = trf_dir + '/' + repeat_filename + '.2.7.7.80.10.50.500.mask'

    # RepeatMasker_command = 'cd ' + trf_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(threads) + ' -noint -x -html -gff -dir ' + trf_dir + ' ' + repeats_path
    # log.logger.debug(RepeatMasker_command)
    # os.system(RepeatMasker_command)
    # trf_masked_repeats = trf_dir + '/' + repeat_filename + '.masked'

    trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    repeats_contigNames, repeats_contigs = read_fasta(repeats_path)
    repeats_path = tmp_output_dir + '/repeats.filter_tandem.fa'
    with open(repeats_path, 'w') as f_save:
        for name in trf_contigNames:
            seq = trf_contigs[name]
            if float(seq.count('N'))/len(seq) < tandem_region_cutoff and len(seq) >= 80:
                f_save.write('>'+name+'\n'+repeats_contigs[name]+'\n')

    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("Step0: use trf to mask genome: %.8s s" % (dtime))
    # mv_result_command = 'mv ' + trf_dir + '/' + ref_filename + '.2.7.7.80.10.50.500.mask ' + trf_dir + '/' + ref_filename + '.2.7.7.80.10.50.500.dat ' + tmp_output_dir
    # log.logger.debug(mv_result_command)
    # os.system(mv_result_command)
    # reference = tmp_output_dir + '/' + ref_filename + '.2.7.7.80.10.50.500.mask'

    # --------------------------------------------------------------------------------------
    # Step5. use blast to execute repeats self_alignment
    # starttime = time.time()
    # blastnResults_path = tmp_output_dir + '/tmpBlastResults.out'
    # candidate_repeats_path = tmp_output_dir + '/candidate_repeats.fa'
    # blast_program_dir = param['RMBlast_Home']
    # if is_sensitive:
    #     log.logger.debug('Start step5: use blast to execute repeats self_alignment')
    #     makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + repeats_path
    #     align_command = blast_program_dir + '/bin/blastn -db ' + repeats_path + ' -num_threads ' + str(threads) + ' -query ' + repeats_path + ' -outfmt 6 > ' + blastnResults_path
    #     log.logger.debug(makedb_command)
    #     os.system(makedb_command)
    #     log.logger.debug(align_command)
    #     os.system(align_command)
    #     # --------------------------------------------------------------------------------------
    #     # Step6. parse blast output, get candidate repeat sequence.
    #     # if a sequence align to multiple sequence, the match
    #     # segments with identity>=80% and length>80bp should be spliced.
    #     starttime1 = time.time()
    #     log.logger.debug('Start step6: splice candidate repeat sequence')
    #     parse_self_blast_output(blastnResults_path, repeats_path, candidate_repeats_path)
    #     endtime1 = time.time()
    #     dtime1 = endtime1 - starttime1
    #     log.logger.debug("module2: get candidate repeat sequence running time: %.8s s" % (dtime1))
    # else:
    #     log.logger.debug('Start step5: use minimap2+bwa to execute repeats alignment')
    #     # # cut raw repeat sequences
    #     # repeats_minimap2 = tmp_output_dir + '/repeats.stage1.minimap2.fa'
    #     # repeats_bwa = tmp_output_dir + '/repeats.stage1.bwa.fa'
    #     # split_repeats(repeats_path, long_repeat_threshold, repeats_minimap2, repeats_bwa)
    #     # use_align_tools = 'minimap2'
    #     # sam_path_minimap2 = run_alignment(repeats_minimap2, reference, use_align_tools, threads, tools_dir)
    #     # use_align_tools = 'bwa'
    #     # sam_path_bwa = run_alignment(repeats_bwa, reference, use_align_tools, threads, tools_dir)
    #     # sam_paths = []
    #     # sam_paths.append(sam_path_minimap2)
    #     # sam_paths.append(sam_path_bwa)
    #     # HS_gap = 0.10
    #     # ID_gap = 0.10
    #     # raw_cut_file = tmp_output_dir + '/repeats.cut.fasta'
    #     # cut_repeat_v1(sam_paths, HS_gap, ID_gap, repeats_path, raw_cut_file)
    #     raw_cut_file = repeats_path
    #
    #     # get multiple alignment
    #     repeats_minimap2 = tmp_output_dir + '/repeats.stage2.minimap2.fa'
    #     repeats_bwa = tmp_output_dir + '/repeats.stage2.bwa.fa'
    #     split_repeats(raw_cut_file, long_repeat_threshold, repeats_minimap2, repeats_bwa)
    #     use_align_tools = 'minimap2'
    #     sam_path_minimap2 = run_alignment(repeats_minimap2, reference, use_align_tools, threads, tools_dir)
    #     use_align_tools = 'bwa'
    #     sam_path_bwa = run_alignment(repeats_bwa, reference, use_align_tools, threads, tools_dir)
    #     sam_paths = []
    #     sam_paths.append((sam_path_minimap2, repeats_minimap2))
    #     sam_paths.append((sam_path_bwa, repeats_bwa))
    #     multi_mapping_repeatIds, not_multi_mapping_repeatIds = get_multiple_alignment_repeat(sam_paths)
    #     # filter not multiple alignment and generate blast like output
    #     generate_blastlike_output(sam_paths, blastnResults_path, not_multi_mapping_repeatIds)
    #     # --------------------------------------------------------------------------------------
    #     # Step6. parse blast output, get candidate repeat sequence.
    #     # if a sequence align to multiple sequence, the match
    #     # segments with identity>=80% and length>80bp should be spliced.
    #     starttime1 = time.time()
    #     log.logger.debug('Start step6: splice candidate repeat sequence')
    #     parse_ref_blast_output(blastnResults_path, reference, candidate_repeats_path)
    #     endtime1 = time.time()
    #     dtime1 = endtime1 - starttime1
    #     log.logger.debug("module2: get candidate repeat sequence running time: %.8s s" % (dtime1))
    # endtime = time.time()
    # dtime = endtime - starttime
    # log.logger.debug("module3: Blastn/(minimap2+bwa) get repeats alignment running time: %.8s s" % (dtime))

    # --------------------------------------------------------------------------------------
    # newly strategy: 2022-04-29 by Kang Hu
    # 01: use bwa to get single mapped sequence
    candidate_repeats_path = repeats_path
    blast_program_dir = param['RMBlast_Home']
    use_align_tools = 'bwa'
    sam_path_bwa = run_alignment(candidate_repeats_path, reference, use_align_tools, threads, tools_dir)
    sam_paths = []
    sam_paths.append(sam_path_bwa)
    HS_gap = 0.10
    ID_gap = 0.10
    # unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds = get_alignment_info(sam_paths)
    unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds, segmental_duplication_repeatIds = get_alignment_info_v1(sam_paths, candidate_repeats_path)
    # blastnResults_path = tmp_output_dir + '/tmpBlastResults.out'
    # makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + reference
    # align_command = blast_program_dir + '/bin/blastn -db ' + reference + ' -num_threads ' + str(threads) + ' -query ' + candidate_repeats_path + ' -outfmt 6 > ' + blastnResults_path
    # log.logger.debug(makedb_command)
    # os.system(makedb_command)
    # log.logger.debug(align_command)
    # os.system(align_command)
    # unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds = get_alignment_info_v2(blastnResults_path)
    single_mapped_path = tmp_output_dir + '/repeats.merge.consensus.single.fa'
    multiple_mapped_path = tmp_output_dir + '/repeats.merge.consensus.multiple.fa'
    segmental_duplication_path = tmp_output_dir + '/segmental_duplication.fa'
    merge_repeat_contigNames, merge_repeat_contigs = read_fasta(candidate_repeats_path)
    # single repeat probably be Chimeric repeat
    with open(single_mapped_path, 'w') as f_save:
        for repeat_id in single_mapped_repeatIds:
            f_save.write('>' + repeat_id + '\n' + merge_repeat_contigs[repeat_id] + '\n')

    with open(multiple_mapped_path, 'w') as f_save:
        for repeat_id in multi_mapping_repeatIds:
            seq = merge_repeat_contigs[repeat_id]
            #seq = seq.replace('N', '')
            f_save.write('>' + repeat_id + '\n' + seq + '\n')

    with open(segmental_duplication_path, 'w') as f_save:
        for repeat_id in segmental_duplication_repeatIds:
            seq = merge_repeat_contigs[repeat_id]
            #seq = seq.replace('N', '')
            f_save.write('>' + repeat_id + '\n' + seq + '\n')

    # 02: align single mapped sequence to reference
    # blastnResults_path = tmp_output_dir + '/repeats.merge.consensus.single.out'
    # makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + reference
    # align_command = blast_program_dir + '/bin/blastn -db ' + reference + ' -num_threads ' + str(
    #     threads) + ' -query ' + single_mapped_path + ' -outfmt 6 > ' + blastnResults_path
    # log.logger.debug(makedb_command)
    # os.system(makedb_command)
    # log.logger.debug(align_command)
    # os.system(align_command)

    # 03: cut segments from blast output
    # this method can be used to Chimeric removal

    # segments_path = tmp_output_dir + '/repeats.merge.consensus.single.segs.fa'
    # seg_threshold = 100
    # query_records = {}
    # single_mapped_segments = {}
    # with open(blastnResults_path, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.split('\t')
    #         query_name = parts[0]
    #         identity = float(parts[2])
    #         q_start = int(parts[6])
    #         q_end = int(parts[7])
    #         if identity < 80:
    #             continue
    #         if not query_records.__contains__(query_name):
    #             query_records[query_name] = []
    #         records = query_records[query_name]
    #         records.append((q_start, q_end))
    #         query_records[query_name] = records
    #
    # for query_name in query_records.keys():
    #     records = query_records[query_name]
    #     records.sort(key=lambda x: (x[0], x[1]))
    #     hash_records = {}
    #     for pos in records:
    #         if not hash_records.__contains__(pos[0]):
    #             hash_records[pos[0]] = []
    #         same_start_records = hash_records[pos[0]]
    #         same_start_records.append(pos)
    #
    #     if not single_mapped_segments.__contains__(query_name):
    #         single_mapped_segments[query_name] = set()
    #     segments_set = single_mapped_segments[query_name]
    #
    #     for i in range(len(records) - 1):
    #         for start_pos in hash_records.keys():
    #             if abs(records[i][0] - start_pos) >= seg_threshold:
    #                 break
    #             same_start_records = hash_records[start_pos]
    #             if same_start_records.__contains__(records[i]):
    #                 same_start_records.remove(records[i])
    #             for pos in same_start_records:
    #                 if abs(records[i][1] - pos[1]) < seg_threshold:
    #                     segments_set.add(records[i])
    #                     segments_set.add(pos)
    #                 else:
    #                     break
    #     single_mapped_segments[query_name] = segments_set
    #
    # with open(segments_path, 'w') as f_save:
    #     for query_name in single_mapped_segments.keys():
    #         for i, segment in enumerate(single_mapped_segments[query_name]):
    #             seg_seq = merge_repeat_contigs[query_name][segment[0] - 1: segment[1]]
    #             f_save.write('>' + query_name + '-seg_' + str(i) + '\n' + seg_seq + '\n')

    # # 04: find all multiple sequences
    # use_align_tools = 'bwa'
    # sam_path_bwa = run_alignment(segments_path, reference, use_align_tools, threads, tools_dir)
    # sam_paths = []
    # sam_paths.append(sam_path_bwa)
    # #unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds = get_alignment_info(sam_paths)
    # unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds = get_alignment_info_v1(sam_paths, segments_path, HS_gap, ID_gap)
    # segs_multiple_mapped_path = tmp_output_dir + '/repeats.merge.consensus.single.segs.consensus.multiple.fa'
    # segments_contigNames, segments_contigs = read_fasta(segments_path)
    # with open(segs_multiple_mapped_path, 'w') as f_save:
    #     for repeat_id in multi_mapping_repeatIds:
    #         seq = segments_contigs[repeat_id]
    #         #seq = seq.replace('N', '')
    #         f_save.write('>' + repeat_id + '\n' + seq + '\n')

    # 06: merge
    merge_pure = tmp_output_dir + '/repeats.merge.pure.fa'
    merge_pure_consensus = tmp_output_dir + '/repeats.merge.pure.consensus.fa'
    os.system('cat ' + multiple_mapped_path + ' >> ' + merge_pure)
    ltr_retriever_seq = tmp_output_dir + '/' + ref_filename + '.mod.LTRlib.fa'
    backjob.join()
    os.system('cat ' + ltr_retriever_seq + ' >> ' + merge_pure)
    #os.system('cat ' + segs_multiple_mapped_path + ' >> ' + merge_pure)
    cd_hit_command = tools_dir + '/cd-hit-est -s 0.8 -c 0.8 -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
    log.logger.debug(cd_hit_command)
    os.system(cd_hit_command)

    # --------------------------------------------------------------------------------------
    # Step7: get complete LTR and TIR sequences
    log.logger.debug('Start step7: get complete LTR and TIR sequences')
    starttime = time.time()
    sensitive_mode = 0
    if is_sensitive:
        sensitive_mode = 1
    TEFinder_cmd = 'python ' + os.getcwd() + '/TEFinder.v1.py --min_TE_len ' + str(min_TE_len) + \
                   ' --min_ltr_complete_len ' + str(min_ltr_complete_len) + ' --max_ltr_complete_len ' + \
                   str(max_ltr_complete_len) + ' --min_ltr_direct_repeat_len ' + str(min_ltr_direct_repeat_len) + \
                   ' --max_ltr_direct_repeat_len ' + str(max_ltr_direct_repeat_len) + ' --min_tir_complete_len ' + \
                   str(min_tir_complete_len) + ' --max_tir_complete_len ' + str(max_tir_complete_len) + \
                   ' --min_tir_direct_repeat_len ' + str(min_tir_direct_repeat_len) + ' --max_tir_direct_repeat_len ' + \
                   str(max_tir_direct_repeat_len) + ' --tmp_output_dir ' + tmp_output_dir + ' -R ' + reference + ' -t ' + \
                   str(threads) + ' --domain_min_identity ' + str(domain_min_identity) + ' --domain_match_ratio ' + str(domain_match_ratio) +\
                   ' -s ' + str(sensitive_mode) + ' --long_repeat_threshold ' + str(long_repeat_threshold) + ' --tools_dir ' + str(tools_dir) + ' --blast_program_dir ' + str(blast_program_dir)
    log.logger.debug(TEFinder_cmd)
    #os.system(TEFinder_cmd)
    # #output of TEFinder
    # ltr_repeats_path = tmp_output_dir + '/ltr_repeats.fa'
    # tir_repeats_path = tmp_output_dir + '/tir_repeats.fa'
    # merge_ltr_tir_path = tmp_output_dir + '/ltr_tir.merge.fa'
    # os.system('cat ' + ltr_repeats_path + ' >> ' + merge_ltr_tir_path)
    # os.system('cat ' + tir_repeats_path + ' >> ' + merge_ltr_tir_path)
    # merge_ltr_tir_consensus = tmp_output_dir + '/ltr_tir.merge.consensus.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS 0.8 -c 0.8 -g 1 -s 0.8 -G 1 -i ' + merge_ltr_tir_path + ' -o ' + merge_ltr_tir_consensus + ' -T 0 -M 0'
    # log.logger.debug(cd_hit_command)
    # os.system(cd_hit_command)
    # # classify
    # sample_name = alias
    # TEClass_home = os.getcwd() + '/classification'
    # TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                   + ' --consensus ' + merge_ltr_tir_consensus + ' --genome ' + reference \
    #                   + ' --thread_num ' + str(threads) + ' -o ' + tmp_output_dir
    # log.logger.debug(TEClass_command)
    # os.system(TEClass_command)
    # classified_merge_ltr_tir_consensus = merge_ltr_tir_consensus + '.final.classified'
    # classified_merge_contignames, classified_merge_contigs = read_fasta(classified_merge_ltr_tir_consensus)
    # # filter LTR not meet length requirement
    # LTR_set = {}
    # TIR_set = {}
    # for name in classified_merge_contignames:
    #     class_name = name.split('#')[1]
    #     sequence = classified_merge_contigs[name]
    #     seq_len = len(sequence)
    #     if class_name.__contains__('LTR'):
    #         if seq_len < min_ltr_complete_len or seq_len > max_ltr_complete_len:
    #             continue
    #         else:
    #             LTR_set[name] = sequence
    #     else:
    #         TIR_set[name] = sequence
    # store_fasta(LTR_set, ltr_repeats_path)
    # store_fasta(TIR_set, tir_repeats_path)
    #
    # endtime = time.time()
    # dtime = endtime - starttime
    # log.logger.debug("module5: get complete LTR and TIR sequences running time: %.8s s" % (dtime))


    merge_repeat_sequences = tmp_output_dir + '/repeats.merge.fa'
    merge_repeat_consensus = tmp_output_dir + '/repeats.merge.consensus.fa'
    # if os.path.exists(ltr_retriever_seq):
    #     # filter overlap sequence with LTR_retriever
    #     filter_ltr_repeats_path = tmp_output_dir + '/ltr_repeats.filter.fa'
    #     blastnResults_path = tmp_output_dir + '/tmp_ltr_blastn.out'
    #     makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + ltr_retriever_seq
    #     align_command = blast_program_dir + '/bin/blastn -db ' + ltr_retriever_seq + ' -num_threads ' + str(
    #         threads) + ' -query ' + ltr_repeats_path + ' -outfmt 6 > ' + blastnResults_path
    #     print(makedb_command)
    #     os.system(makedb_command)
    #     print(align_command)
    #     os.system(align_command)
    #     filter_LTR_high_similarity(blastnResults_path, ltr_retriever_seq, ltr_repeats_path, filter_ltr_repeats_path)
    #     os.system('cat ' + filter_ltr_repeats_path + ' >> ' + merge_repeat_sequences)
    # else:
    #     os.system('cat ' + ltr_repeats_path + ' >> ' + merge_repeat_sequences)
    # os.system('cat ' + tir_repeats_path + ' >> ' + merge_repeat_sequences)
    # os.system('cat ' + merge_pure_consensus + ' >> ' + merge_repeat_sequences)


    # starttime = time.time()
    # tools_dir = os.getcwd() + '/tools'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS 0.8 -c 0.8 -g 1 -s 0.8 -G 1 -i ' + merge_repeat_sequences + ' -o ' + merge_repeat_consensus + ' -T 0 -M 0'
    # log.logger.debug(cd_hit_command)
    # os.system(cd_hit_command)
    # endtime = time.time()
    # dtime = endtime - starttime
    # log.logger.debug("module7: get merge repeat consensus sequence running time: %.8s s" % (dtime))

    os.system('cat ' + merge_pure_consensus + ' >> ' + merge_repeat_consensus)
    # --------------------------------------------------------------------------------------
    # Step10. run TE classification to classify TE family
    starttime = time.time()
    log.logger.debug('Start step8: get classified consensus sequence')
    sample_name = alias
    TEClass_home = os.getcwd() + '/classification'
    TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
                      + ' --consensus ' + merge_repeat_consensus + ' --genome ' + reference \
                      + ' --thread_num ' + str(threads) + ' -o ' + tmp_output_dir
    log.logger.debug(TEClass_command)
    os.system(TEClass_command)

    # --------------------------------------------------------------------------------------
    # Step11. assign a family name for each classified TE consensus
    classified_consensus_path = merge_repeat_consensus + '.final.classified'
    classified_contignames, classified_contigs = read_fasta(classified_consensus_path)
    # filter LTR not meet length requirement
    repeat_set = {}
    for name in classified_contignames:
        class_name = name.split('#')[1]
        sequence = classified_contigs[name]
        seq_len = len(sequence)
        repeat_set[name] = sequence
        # if class_name.__contains__('LTR'):
        #     if seq_len < min_ltr_complete_len or seq_len > max_ltr_complete_len:
        #         continue
        #     else:
        #         repeat_set[name] = sequence
        # else:
        #     repeat_set[name] = sequence
    store_fasta(repeat_set, classified_consensus_path)

    # # filter unknown and short(<80bp)  TE
    # filter_classified_consensus_path = merge_repeat_consensus + '.final.filter_unknown.classified'
    # classified_contignames, classified_contigs = read_fasta(classified_consensus_path)
    # with open(filter_classified_consensus_path, 'w') as f_save:
    #     for name in classified_contignames:
    #         class_name = name.split('#')[1]
    #         if class_name == 'Unknown' or len(classified_contigs[name]) < 80:
    #             continue
    #         f_save.write('>' + name + '\n' + classified_contigs[name] + '\n')

    classified_contigNames, classified_contigs = read_fasta(classified_consensus_path)
    family_path = output_dir + '/family_' + sample_name + '.fasta'
    with open(family_path, 'w') as f_save:
        for f_id, name in enumerate(classified_contigNames):
            sequence = classified_contigs[name]
            if len(sequence) < 80:
                continue
            class_name = name.split('#')[1]
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
    RepeatMasker_output_dir = output_dir + '/' + sample_name
    RepeatMasker_command = 'cd ' + tmp_output_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(threads) \
                           + ' -lib ' + family_path + ' -nolow -x -html -gff -dir ' + RepeatMasker_output_dir + ' ' + reference
    os.system('rm -rf ' + RepeatMasker_output_dir)
    log.logger.debug(RepeatMasker_command)
    os.system(RepeatMasker_command)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("module9: invoke RepeatMasker to annotate genome running time: %.8s s" % (dtime))




