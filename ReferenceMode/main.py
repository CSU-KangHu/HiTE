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
    store_fasta, printClass, parse_ref_blast_output, filter_LTR_high_similarity, get_alignment_info_v3

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


def compare_seq(self_info, other_info, identity_cutoff, length_similarity_cutoff,
                refContigs, output_dir, seq_idx, blast_program_dir):
    # self_info = (c, ref_name, combine_frag_start, combine_frag_end)
    ref_name = self_info[1]
    ref_seq = refContigs[ref_name]

    self_combine_frag_start = self_info[2]
    self_combine_frag_end = self_info[3]

    other_combine_frag_start = other_info[2]
    other_combine_frag_end = other_info[3]

    self_seq = ref_seq[self_combine_frag_start: self_combine_frag_end]
    self_contigs = {}
    self_contigs['self'] = self_seq

    other_seq = ref_seq[other_combine_frag_start: other_combine_frag_end]
    other_contigs = {}
    other_contigs['other'] = other_seq

    self_seq_path = output_dir + '/self_' + str(seq_idx) + '.fa'
    other_seq_path = output_dir + '/other_' + str(seq_idx) + '.fa'
    blastnResults_path = output_dir + '/blast_' + str(seq_idx) + '.out'
    store_fasta(self_contigs, self_seq_path)
    store_fasta(other_contigs, other_seq_path)

    makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + other_contigs
    align_command = blast_program_dir + '/bin/blastn -db ' + other_contigs + ' -query ' + self_contigs + ' -outfmt 6 > ' + blastnResults_path
    log.logger.debug(makedb_command)
    os.system(makedb_command)
    log.logger.debug(align_command)
    os.system(align_command)

    query_name_set = set()
    target_name_set = set()
    query_cluster = {}
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

            query_len = len(self_seq)
            target_len = len(other_seq)
            key = query_name + '$' +target_name
            if not query_cluster.__contains__(key):
                query_cluster[key] = ([], -1, -1)
            tuple = query_cluster[key]
            cluster = tuple[0]
            if identity >= identity_cutoff:
                cluster.append((q_start, q_end, t_start, t_end, identity))
            query_cluster[key] = (cluster, query_len, target_len)

    total_identity = 0
    avg_identity = 0
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

            identity = record[4]
            total_identity += identity
        avg_identity = float(total_identity)/len(tuple[0])
        for j in range(len(query_array)):
            if query_array[j] == 'X':
                query_masked_len += 1
        for j in range(len(target_array)):
            if target_array[j] == 'X':
                target_masked_len += 1
        if float(query_masked_len)/query_len >= length_similarity_cutoff and float(target_masked_len)/target_len >= length_similarity_cutoff:
            return avg_identity




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
    repeat_segments = generate_candidate_repeats_v1(cur_segments, k_num, unique_kmer_map, fault_tolerant_bases)

    repeats_path = tmp_output_dir + '/repeats.fa'
    node_index = 0
    with open(repeats_path, 'w') as f_save:
        for repeat in repeat_segments:
            f_save.write('>Node_'+str(node_index)+'-len_'+str(len(repeat))+'\n'+repeat+'\n')
            node_index += 1

    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("module1: Generate repeat file running time: %.8s s" % (dtime))

    repeats_consensus = tmp_output_dir + '/repeats.consensus.fa'
    cd_hit_command = tools_dir + '/cd-hit-est -s 0.95 -c 0.95 -i ' + repeats_path + ' -o ' + repeats_consensus + ' -T 0 -M 0'
    log.logger.debug(cd_hit_command)
    os.system(cd_hit_command)


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
    merge_repeat_contigNames, merge_repeat_contigs = read_fasta(candidate_repeats_path)
    # single repeat probably be Chimeric repeat
    with open(repeat_freq_path, 'w') as f_save:
        for repeat_id in new_mapping_repeatIds.keys():
            freq = new_mapping_repeatIds[repeat_id][0]
            seq = merge_repeat_contigs[repeat_id]
            if freq <= 1 or (freq < 5 and len(seq) < 80):
                continue
            f_save.write('>' + repeat_id + '\tcopies=' + str(freq+1) + '\n' + seq + '\n')

    # --------------------------------------------------------------------------------------
    # skip variation between fragments: 2022-05-19 by Kang Hu
    #sort by start and end position
    sorted_fragments = {k: v for k, v in sorted(query_position.items(), key=lambda item: (item[1][1], item[1][2]))}

    # find all regions could be connected, threshold = 200bp
    # region_list keeps all regions, which include all fragments can be connected
    # region_list = {
        # R1: {
            # ref_name: [(F1, start, end),(F2, start, end),(F3, start, end),(F4, start, end)]
        # }
    # }
    skip_threshold = 200
    region_list = {}
    last_end_pos = -1
    region_index = 0
    for ref_name in sorted_fragments.keys():
        cur_ref_fragments = sorted_fragments[ref_name]
        for item in cur_ref_fragments:
            query_name = item[0]
            start = item[1]
            end = item[2]
            region_id = 'R'+str(region_index)
            if not region_list.__contains__(region_id):
                region_list[region_id] = {}
            cur_region_dict = region_list[region_id]
            if not cur_region_dict.__contains__(ref_name):
                cur_region_dict[ref_name] = []
            cur_region_list = cur_region_dict[ref_name]
            if last_end_pos == -1:
                # first fragment add into cur_region_list directly
                cur_region_list.append((query_name, start, end))
            else:
                # cur fragment close to last fragment
                if start - last_end_pos < skip_threshold:
                    cur_region_list.append((query_name, start, end))
                else:
                    # cur fragment far from last fragment, start a new region
                    region_list[region_id] = cur_region_list
                    region_index += 1
                    region_id = 'R' + str(region_index)
                    if not region_list.__contains__(region_id):
                        region_list[region_id] = {}
                    cur_region_dict = region_list[region_id]
                    if not cur_region_dict.__contains__(ref_name):
                        cur_region_dict[ref_name] = []
                    cur_region_list = cur_region_dict[ref_name]
                    cur_region_list.append((query_name, start, end))
            last_end_pos = end

    # region_combination keeps all combination of fragment in one region
    # e.g., region_combination = {
        # R1: {
                # max_combination_len : 3,
                # combinations: {
                    # c=3: [F1F2F3],
                    # c=2: [F1F2, F2F3],
                    # c=1: [F1,F2,F3]
                # }
            # }
        # }

    # frag_hash keeps all fragment combination information: combination_len, reference, start, end
    # e.g., frag_hash = {
        # F1F2: {
            # R1: (c=2, ref_name, start, end)
            # R2: (c=2, ref_name, start, end)
        # }
    # }
    region_combination = {}
    frag_hash = {}
    for region_id in region_list.keys():
        cur_region_dict = region_list[region_id]
        for ref_name in cur_region_dict.keys():
            cur_region_list = cur_region_dict[ref_name]
            max_combination_len = len(cur_region_list)
            if not region_combination.__contains__(region_id):
                region_combination[region_id] = {}
            cur_region_combination = region_combination[region_id]
            cur_region_combination['max_combination_len'] = max_combination_len
            combinations = {}
            for c in range(1, max_combination_len+1):
                if not combinations.__contains__(c):
                    combinations[c] = []
                cur_combinations = combinations[c]
                for left in range(len(cur_region_list)-c):
                    # connect fragments with len=c
                    combine_name = ''
                    combine_frag_start = -1
                    combine_frag_end = -1
                    for l in range(c):
                        cur_frag = cur_region_list[left + l]
                        cur_frag_start = cur_frag[1]
                        cur_frag_end = cur_frag[2]
                        if combine_frag_start == -1:
                            combine_frag_start = cur_frag_start
                        combine_frag_end = cur_frag_end
                        if combine_name != '':
                            combine_name += ','
                        combine_name += cur_frag[0]
                    cur_combinations.append(combine_name)
                    if not frag_hash.__contains__(combine_name):
                        frag_hash[combine_name] = {}
                    cur_frag_dict = frag_hash[combine_name]
                    cur_frag_dict[region_id] = (c, ref_name, combine_frag_start, combine_frag_end)
                    frag_hash[combine_name] = cur_frag_dict

                    combine_name = ''
                    combine_frag_start = -1
                    combine_frag_end = -1
                combinations[c] = cur_combinations
            cur_region_combination['combinations'] = combinations

    refNames, refContigs = read_fasta(reference)
    # go through each region, find candidate combine fragments
    identity_threshold = 0.95
    length_similarity_cutoff = 0.95
    # regionContigs keeps all fragments in each region
    regionContigs = {}
    seq_idx = 0
    # go through each region
    for region_id in region_combination.keys():
        cur_region_combination = region_combination[region_id]
        combinations = cur_region_combination['combinations']
        max_combination_len = cur_region_combination['max_combination_len']
        find_best_combine = False
        if find_best_combine:
            continue
        # start from max length combination
        for c in range(max_combination_len, 0, -1):
            cur_combinations = combinations[c]
            max_identity = 0
            best_combine_name = None
            for combine_name in cur_combinations:
                # find in frag_hash
                frag_region_dict = frag_hash[combine_name]
                # self_info = (c, ref_name, combine_frag_start, combine_frag_end)
                self_info = frag_region_dict[region_id]
                for other_region_id in frag_region_dict.keys():
                    if region_id != other_region_id:
                        other_info = frag_region_dict[other_region_id]
                        # self info similar to other_info
                        identity = compare_seq(self_info, other_info, identity_threshold, length_similarity_cutoff, refContigs, output_dir, seq_idx, blast_program_dir)
                        seq_idx += 1
                        if identity > max_identity:
                            max_identity = identity
                            best_combine_name = combine_name
            # if current combination reach score threshold, then it can be used to replace the whole region
            if max_identity >= identity_threshold:
                if not regionContigs.__contains__(region_id):
                    regionContigs[region_id] = []
                final_frags = regionContigs[region_id]
                ref_name = self_info[1]
                # replace the whole region with best combine
                final_frags.append(best_combine_name, ref_name, self_info[2], self_info[3])
                # all_frags = [(F1, start, end),(F2, start, end),(F3, start, end),(F4, start, end)]
                all_frags = region_list[region_id][ref_name]
                best_frags = best_combine_name.split(',')
                # keep other fragments
                for frag in all_frags:
                    if frag[0] not in best_frags:
                        final_frags.append(frag[0], ref_name, frag[1], frag[2])
                regionContigs[region_id] = final_frags
                find_best_combine = True
                break

    connected_repeats = tmp_output_dir + '/repeats.connect.fa'
    with open(connected_repeats, 'w') as f_save:
        node_index = 0
        for region_id in regionContigs.keys():
            for frag in regionContigs[region_id]:
                seq = refContigs[frag[1]][frag[2]: frag[3]]
                f_save.write('>Node_'+node_index+'\n'+seq+'\n')


    # 06: merge
    merge_pure = tmp_output_dir + '/repeats.merge.pure.fa'
    merge_pure_consensus = tmp_output_dir + '/repeats.merge.pure.consensus.fa'
    os.system('cat ' + repeat_freq_path + ' >> ' + merge_pure)
    ltr_retriever_seq = tmp_output_dir + '/' + ref_filename + '.mod.LTRlib.fa'
    backjob.join()
    os.system('cat ' + ltr_retriever_seq + ' >> ' + merge_pure)
    cd_hit_command = tools_dir + '/cd-hit-est -s 0.95 -c 0.95 -i ' + merge_pure + ' -o ' + merge_pure_consensus + ' -T 0 -M 0'
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




