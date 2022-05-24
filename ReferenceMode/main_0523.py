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
    store_fasta, printClass, parse_ref_blast_output, filter_LTR_high_similarity, get_alignment_info_v3, compare_seq, getRegionCombination, getCombineFragments

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
    log.logger.debug('partition %d process finished' % (partiton_index))
    return repeat_list

def generate_candidate_repeats_v1(cur_segments, k_num, unique_kmer_map, fault_tolerant_bases):
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
    last_start_pos = -1
    last_end_pos = -1
    for seq_index, cur_masked_segment in enumerate(cur_masked_segments):
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
                    repeat_list.append((last_start_pos, cur_repeat_str))
                    cur_repeat_str = ''
                    try_connect_str = ''
                    last_start_pos = -1
        # keep last masked sequence
        if cur_repeat_str != '':
            repeat_list.append((last_start_pos, cur_repeat_str))
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


    skip_threshold = 200
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

    # create thread pool, use multiple processes to execute
    kmer_segments = []
    with open(unique_kmer_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            kmer_segments.append(line)
    kmer_segments_cluster = split2cluster(kmer_segments, partitions_num)

    ex = ProcessPoolExecutor(partitions_num)
    objs = []
    for partiton_index in kmer_segments_cluster.keys():
        cur_kmer_segments = kmer_segments_cluster[partiton_index]
        obj = ex.submit(getUniqueKmer_v1, cur_kmer_segments, partiton_index)
        objs.append(obj)
    ex.shutdown(wait=True)
    unique_kmer_map = {}
    for obj in as_completed(objs):
        for kmer in obj.result():
            unique_kmer_map[kmer] = 1
    log.logger.debug('generate unique_kmer_map finished')

    # --------------------------------------------------------------------------------------
    # Step3. split reference into segments
    log.logger.debug('Start step3: split reference into segments')
    reduce_partitions_num = judgeReduceThreads(unique_kmer_path, partitions_num, log)

    cur_segments = convertToUpperCase(reference)
    repeat_segments = generate_candidate_repeats(cur_segments, k_num, unique_kmer_map, fault_tolerant_bases)

    repeats_path = tmp_output_dir + '/repeats.fa'
    node_index = 0
    with open(repeats_path, 'w') as f_save:
        for repeat in repeat_segments:
            f_save.write('>N'+str(node_index)+'\n'+repeat+'\n')
            #f_save.write('>Node_' + str(node_index) + '-len_' + str(len(repeat)) + '\n' + repeat + '\n')
            node_index += 1

    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("module1: Generate repeat file running time: %.8s s" % (dtime))

    repeats_consensus = tmp_output_dir + '/repeats.consensus.fa'
    cd_hit_command = tools_dir + '/cd-hit-est -s ' + str(length_similarity_cutoff) + ' -c ' + str(identity_threshold) + ' -i ' + repeats_path + ' -o ' + repeats_consensus + ' -T 0 -M 0'
    log.logger.debug(cd_hit_command)
    os.system(cd_hit_command)

    # try to chain all fragments
    # --------------------------------------------------------------------------------------
    # newly strategy: 2022-04-29 by Kang Hu
    # 01: use bwa to get single mapped sequence
    candidate_repeats_path = repeats_consensus
    blast_program_dir = param['RMBlast_Home']
    use_align_tools = 'bwa'
    sam_path_bwa = run_alignment(candidate_repeats_path, reference, use_align_tools, threads, tools_dir)
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

    sam_path_bwa = run_alignment(repeat_multiple_path, reference, use_align_tools, threads, tools_dir)
    sam_paths = []
    sam_paths.append(sam_path_bwa)
    # unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds = get_alignment_info(sam_paths)
    new_mapping_repeatIds, query_position = get_alignment_info_v3(sam_paths, repeat_multiple_path)

    # --------------------------------------------------------------------------------------
    # skip variation between fragments: 2022-05-21 by Kang Hu
    # sort by start and end position
    sorted_fragments = {}
    for ref_name in query_position.keys():
        position_array = query_position[ref_name]
        position_array.sort(key=lambda x: (x[1], x[2]))
        sorted_fragments[ref_name] = position_array

    # find all regions could be connected, threshold = 200bp
    # region_list keeps all regions, which include all fragments can be connected
    # region_list = {
    #   R1: {
    #       ref_name: (F1, start, end)
    #   }
    # }
    region_list = {}
    region_index = 0
    for ref_name in sorted_fragments.keys():
        cur_ref_fragments = sorted_fragments[ref_name]
        last_end_pos = -1
        for item in cur_ref_fragments:
            query_name = item[0]
            start = item[1]
            end = item[2]
            region_id = 'R' + str(region_index)
            if not region_list.__contains__(region_id):
                region_list[region_id] = {}
            cur_region_dict = region_list[region_id]
            if last_end_pos == -1:
                # first fragment add into cur_region_list directly
                cur_region_dict[ref_name] = (query_name, start, end)
            else:
                # cur fragment close to last fragment
                if start - last_end_pos < skip_threshold:
                    last_frag = cur_region_dict[ref_name]
                    new_query_name = last_frag[0]+','+query_name
                    if end > last_frag[2]:
                        cur_region_dict[ref_name] = (new_query_name, last_frag[1], end)
                    else:
                        cur_region_dict[ref_name] = (new_query_name, last_frag[1], last_frag[2])
                else:
                    # cur fragment far from last fragment, start a new region
                    region_index += 1
                    region_id = 'R' + str(region_index)
                    if not region_list.__contains__(region_id):
                        region_list[region_id] = {}
                    cur_region_dict = region_list[region_id]
                    cur_region_dict[ref_name] = (query_name, start, end)
                    region_list[region_id] = cur_region_dict
            region_list[region_id] = cur_region_dict
            last_end_pos = end

    refNames, refContigs = read_fasta(reference)
    region_repeats = tmp_output_dir + '/repeats.region.fa'
    with open(region_repeats, 'w') as f_save:
        node_index = 0
        for region_id in region_list.keys():
            region_item = region_list[region_id]
            for ref_name in region_item.keys():
                frag_item = region_item[ref_name]
                seq = refContigs[ref_name][frag_item[1]: frag_item[2]]
                f_save.write('>R_' + str(node_index) + '\n' + seq + '\n')
                node_index += 1

    regionNames, regionContigs = read_fasta(region_repeats)


    query_records = {}

    # use_align_tools = 'bwa'
    # sam_path_bwa = run_alignment(region_repeats, reference, use_align_tools, threads, tools_dir)
    #
    # samfile = pysam.AlignmentFile(sam_path_bwa, "rb")
    # for read in samfile.fetch():
    #     if read.is_unmapped:
    #         continue
    #     query_name = read.query_name
    #     reference_name = read.reference_name
    #     cigar = read.cigartuples
    #     cigarstr = read.cigarstring
    #     NM_tag = 0
    #     try:
    #         NM_tag = read.get_tag('NM')
    #     except KeyError:
    #         NM_tag = -1
    #     identity = compute_identity(cigarstr, NM_tag, 'BLAST')
    #     identity = float(identity) * 100
    #     is_reverse = read.is_reverse
    #     alignment_len = read.query_alignment_length
    #     q_start = int(read.query_alignment_start)
    #     q_end = int(read.query_alignment_end)
    #     t_start = int(read.reference_start)
    #     t_end = int(read.reference_end)
    #     query_len = len(regionContigs[query_name])
    #
    #     if not query_records.__contains__(query_name):
    #         query_records[query_name] = []
    #     records = query_records[query_name]
    #     records.append((reference_name, identity, alignment_len, query_len, q_start, q_end, t_start, t_end))
    #     query_records[query_name] = records

    # blastn produce too many alignments
    blastnResults_path = tmp_output_dir + '/region.out'
    makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + reference
    align_command = blast_program_dir + '/bin/blastn -db ' + reference + ' -num_threads ' + str(threads) + ' -query ' + region_repeats + ' -word_size' + str(k_num) + ' -outfmt 6 > ' + blastnResults_path
    log.logger.debug(makedb_command)
    os.system(makedb_command)
    log.logger.debug(align_command)
    os.system(align_command)

    with open(blastnResults_path, 'r') as f_r:
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

            query_len = len(regionContigs[query_name])

            if not query_records.__contains__(query_name):
                query_records[query_name] = []
            records = query_records[query_name]
            records.append((target_name, identity, match_base, query_len, q_start, q_end, t_start, t_end))
            query_records[query_name] = records

    keep_repeats = {}
    for query_name in query_records.keys():
        complete_alignment_num = 0
        to_splice_seq = set()
        for record in query_records[query_name]:
            identity = record[1]
            match_base = record[2]
            query_len = record[3]
            q_start = record[4]
            q_end = record[5]
            if float(match_base) / query_len >= 0.95 and identity >= 95:
                complete_alignment_num += 1
            elif identity >= 95:
                to_splice_seq.add((q_start, q_end))
        if complete_alignment_num > 1:
            # keep as true repeat
            keep_repeats[query_name] = regionContigs[query_name]
        else:
            # keep cut sequence
            for i, seq in enumerate(to_splice_seq):
                new_query_name = query_name + '-s_' + str(i)
                cut_seq = regionContigs[query_name][seq[0]: seq[1]]
                keep_repeats[new_query_name] = seq

    region_cut_repeats = tmp_output_dir + '/repeats.region.cut.fa'
    region_cut_consensus = tmp_output_dir + '/repeats.region.cut.consensus.fa'
    store_fasta(keep_repeats, region_cut_repeats)
    # cd-hit remove redudant sequences
    cd_hit_command = tools_dir + '/cd-hit-est -s ' + str(length_similarity_cutoff) + ' -c ' + str(identity_threshold) + ' -i ' + region_cut_repeats + ' -o ' + region_cut_consensus + ' -T 0 -M 0'
    log.logger.debug(cd_hit_command)
    os.system(cd_hit_command)

    # bwa remove not multiple alignment repeats
    use_align_tools = 'bwa'
    sam_path_bwa = run_alignment(region_cut_consensus, reference, use_align_tools, threads, tools_dir)
    sam_paths = []
    sam_paths.append(sam_path_bwa)
    new_mapping_repeatIds, query_position = get_alignment_info_v3(sam_paths, region_cut_consensus)
    repeat_multiple_path = tmp_output_dir + '/repeats.region.cut.multiple.fa'
    region_cut_contigNames, region_cut_contigs = read_fasta(region_cut_consensus)
    with open(repeat_multiple_path, 'w') as f_save:
        for repeat_id in new_mapping_repeatIds.keys():
            freq = new_mapping_repeatIds[repeat_id][0]
            seq = region_cut_contigs[repeat_id]
            if freq <= 1:
                continue
            f_save.write('>' + repeat_id + '\n' + seq + '\n')

    # 06: merge
    merge_pure = tmp_output_dir + '/repeats.merge.pure.fa'
    merge_pure_consensus = tmp_output_dir + '/repeats.merge.pure.consensus.fa'
    os.system('cat ' + repeat_multiple_path + ' >> ' + merge_pure)
    ltr_retriever_seq = tmp_output_dir + '/' + ref_filename + '.mod.LTRlib.fa'
    backjob.join()
    os.system('cat ' + ltr_retriever_seq + ' >> ' + merge_pure)
    cd_hit_command = tools_dir + '/cd-hit-est -s ' + str(length_similarity_cutoff) + ' -c ' + str(identity_threshold) + ' -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
    log.logger.debug(cd_hit_command)
    os.system(cd_hit_command)

    starttime = time.time()
    # Step0. use RepeatMasker/trf to mask all low complexity/tandem repeats in raw repeat region
    # >= tandem_region_cutoff region of the whole repeat region, then it should be filtered, since maybe false positive
    TRF_Path = param['TRF_Path']

    trf_dir = tmp_output_dir + '/trf_temp'
    if not os.path.exists(trf_dir):
        os.makedirs(trf_dir)
    (repeat_dir, repeat_filename) = os.path.split(merge_pure_consensus)
    (repeat_name, repeat_extension) = os.path.splitext(repeat_filename)
    trf_command = 'cd ' + trf_dir + ' && ' + TRF_Path + ' ' + merge_pure_consensus + ' 2 7 7 80 10 50 500 -f -d -m'
    log.logger.debug(trf_command)
    os.system(trf_command)
    trf_masked_repeats = trf_dir + '/' + repeat_filename + '.2.7.7.80.10.50.500.mask'

    trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    repeats_contigNames, repeats_contigs = read_fasta(merge_pure_consensus)
    repeats_path = tmp_output_dir + '/repeats.filter_tandem.fa'
    with open(repeats_path, 'w') as f_save:
        for name in trf_contigNames:
            seq = trf_contigs[name]
            if float(seq.count('N')) / len(seq) < tandem_region_cutoff:
                f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')

    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("Step0: use trf to mask genome: %.8s s" % (dtime))

    # --------------------------------------------------------------------------------------
    # Step10. run TE classification to classify TE family
    starttime = time.time()
    log.logger.debug('Start step8: get classified consensus sequence')
    sample_name = alias
    TEClass_home = os.getcwd() + '/classification'
    TEClass_command = 'cd ' + TEClass_home + ' && python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
                      + ' --consensus ' + repeats_path + ' --genome ' + reference \
                      + ' --thread_num ' + str(threads) + ' -o ' + tmp_output_dir
    log.logger.debug(TEClass_command)
    os.system(TEClass_command)

    # --------------------------------------------------------------------------------------
    # Step11. assign a family name for each classified TE consensus
    classified_consensus_path = repeats_path + '.final.classified'
    classified_contigNames, classified_contigs = read_fasta(classified_consensus_path)
    family_path = tmp_output_dir + '/family_' + sample_name + '.fasta'
    with open(family_path, 'w') as f_save:
        for f_id, name in enumerate(classified_contigNames):
            sequence = classified_contigs[name]
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
    RepeatMasker_output_dir = output_dir + '/' + sample_name
    RepeatMasker_command = 'cd ' + tmp_output_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(threads) \
                           + ' -lib ' + family_path + ' -nolow -x -html -gff -dir ' + RepeatMasker_output_dir + ' ' + reference
    os.system('rm -rf ' + RepeatMasker_output_dir)
    log.logger.debug(RepeatMasker_command)
    os.system(RepeatMasker_command)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("module9: invoke RepeatMasker to annotate genome running time: %.8s s" % (dtime))




