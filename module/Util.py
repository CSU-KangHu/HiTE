#!/usr/bin/env python
# The Util.py file contains many useful functions during development.
# Some of them are no longer in use, but they are kept for future reference and convenience.
import codecs
import json
import multiprocessing
import os
import sys
from datetime import datetime
import logging
import re
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Manager
from logging import handlers

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)

import numpy as np
import pandas as pd
from fuzzysearch import find_near_matches
from collections import Counter, defaultdict
import Levenshtein
import seaborn as sns
import subprocess
from PyPDF2 import PdfMerger
from matplotlib import pyplot as plt


class Logger(object):
    level_relations = {
        'debug':logging.DEBUG,
        'info':logging.INFO,
        'warning':logging.WARNING,
        'error':logging.ERROR,
        'crit':logging.CRITICAL
    }

    def __init__(self,filename,level='info',when='D',backCount=3,fmt='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'):
        self.logger = logging.getLogger(filename)
        format_str = logging.Formatter(fmt)
        self.logger.setLevel(self.level_relations.get(level))
        sh = logging.StreamHandler()
        sh.setFormatter(format_str)
        th = handlers.TimedRotatingFileHandler(filename=filename,when=when,backupCount=backCount,encoding='utf-8')
        th.setFormatter(format_str)
        self.logger.addHandler(sh)
        self.logger.addHandler(th)

if __name__ == '__main__':
    log = Logger('all.log',level='debug')
    log.logger.debug('debug')
    log.logger.info('info')
    Logger('error.log', level='error').logger.error('error')

def print_seqs(header, sequence, length, outfile):
    print('>' + header, file=outfile)
    while len(sequence) > 0:
        print(sequence[:length], file=outfile)
        sequence = sequence[length:]

def run_TRsearch(TRsearch_dir, longest_repeats_multi_line_path, longest_repeats_multi_line_dir):
    TRsearch_command = TRsearch_dir + '/TRsearch -i 0.75 ' + longest_repeats_multi_line_path
    os.system('cd ' + longest_repeats_multi_line_dir + ' && ' + TRsearch_command + '> /dev/null 2>&1')
    TR_out = longest_repeats_multi_line_path + '.TR.set'
    return TR_out

def run_command_with_timeout(command, timeout):
    process = subprocess.Popen(command, shell=True)
    start_time = time.time()
    while True:
        if process.poll() is not None:
            break
        if time.time() - start_time > timeout:
            process.terminate()
            process.wait()
            raise TimeoutError(f"Command '{command}' timed out after {timeout} seconds")
    return process.returncode

def run_HelitronScanner(sh_dir, temp_dir, cur_candidate_Helitrons_path, HSDIR, HSJAR, partition_index, debug):
    # cur_candidate_Helitrons_path = temp_dir + '/' + str(partition_index) + '.fa'
    # cur_candidate_Helitrons = {}
    # for item in cur_segments:
    #     query_name = item[0]
    #     orig_seq = item[1]
    #     cur_candidate_Helitrons[query_name] = orig_seq
    # store_fasta(cur_candidate_Helitrons, cur_candidate_Helitrons_path)

    HelitronScanner_command = 'cd ' + temp_dir + ' && ' + 'sh ' + sh_dir + '/run_helitron_scanner.sh ' \
                              + str(partition_index) + ' ' + cur_candidate_Helitrons_path + ' ' + HSDIR + ' ' + HSJAR + '> /dev/null 2>&1'
    # os.system(HelitronScanner_command + '> /dev/null 2>&1')
    # 在某些情况下，未知原因会导致HelitronScanner执行卡死，我们给每个进程限制最大的运行时间 5 min，如果还不结束就直接kill掉
    timeout = 300  # 5min
    try:
        return_code = run_command_with_timeout(HelitronScanner_command, timeout)
    except TimeoutError as e:
        print(e)
    if debug:
        print(HelitronScanner_command)

    cur_helitron_out = temp_dir + '/' + str(partition_index) + '.HelitronScanner.draw.hel.fa'
    cur_rc_helitron_out = temp_dir + '/' + str(partition_index) + '.HelitronScanner.draw.rc.hel.fa'
    cur_names, cur_contigs = read_fasta(cur_helitron_out)
    cur_rc_names, cur_rc_contigs = read_fasta(cur_rc_helitron_out)
    candidate_Helitrons = {}
    candidate_Helitrons.update(cur_contigs)
    candidate_Helitrons.update(cur_rc_contigs)
    return candidate_Helitrons

def run_HelitronScanner_v1(sh_dir, temp_dir, candidate_file, HSDIR, HSJAR, prefix):
    HelitronScanner_command = 'cd ' + temp_dir + ' && ' + 'sh ' + sh_dir + '/run_helitron_scanner.sh ' \
                              + str(prefix) + ' ' + candidate_file + ' ' + HSDIR + ' ' + HSJAR
    print(HelitronScanner_command)
    os.system(HelitronScanner_command + '> /dev/null 2>&1')

    cur_helitron_out = temp_dir + '/' + str(prefix) + '.HelitronScanner.draw.hel.fa'
    cur_rc_helitron_out = temp_dir + '/' + str(prefix) + '.HelitronScanner.draw.rc.hel.fa'
    cur_names, cur_contigs = read_fasta(cur_helitron_out)
    cur_rc_names, cur_rc_contigs = read_fasta(cur_rc_helitron_out)
    candidate_Helitrons = {}
    candidate_Helitrons.update(cur_contigs)
    candidate_Helitrons.update(cur_rc_contigs)
    return candidate_Helitrons, prefix

def run_EAHelitron(flanking_len, temp_dir, all_candidate_helitron_path, EAHelitron, partition_index):
    all_candidate_helitron_contigs = {}
    contigNames, contigs = read_fasta(all_candidate_helitron_path)
    for query_name in contigNames:
        seq = contigs[query_name]

        raw_start = flanking_len + 1
        raw_end = len(seq) - flanking_len

        if seq.__contains__('NNNNNNNNNN'):
            continue
        new_query_name = query_name + '-rawstart_' + str(raw_start) + '-rawend_' + str(raw_end)
        all_candidate_helitron_contigs[new_query_name] = seq
    store_fasta(all_candidate_helitron_contigs, all_candidate_helitron_path)
    EAHelitron_command = 'cd ' + temp_dir + ' && ' + 'perl ' + EAHelitron + '/EAHelitron -o ' + str(partition_index) + ' -u 20000 -T "ATC" -r 3 ' + all_candidate_helitron_path
    os.system(EAHelitron_command + '> /dev/null 2>&1')

    all_EAHelitron_res = temp_dir + '/' + str(partition_index) + '.5.fa'
    all_copies_out_names, all_copies_out_contigs = read_fasta_v1(all_EAHelitron_res)
    # group
    # group_copies_contigs -> {query_name: {name: seq}}
    group_copies_contigs = {}
    for cur_name in all_copies_out_contigs.keys():
        raw_name = cur_name.split(' ')[1]
        parts = raw_name.split(':')
        query_name = ':'.join(parts[:-1])
        if not group_copies_contigs.__contains__(query_name):
            group_copies_contigs[query_name] = {}
        cur_copies_out_contigs = group_copies_contigs[query_name]
        cur_copies_out_contigs[cur_name] = all_copies_out_contigs[cur_name]

    candidate_Helitrons = {}
    for query_name in group_copies_contigs.keys():
        cur_copies_out_contigs = group_copies_contigs[query_name]
        copies_candidate = {}
        # 2. Merge, select the sequence with the smallest distance (if same, get the longest)
        # as the representative sequence of this copy.
        # copies_candidate -> {copy_index: (min_distance_seq_name, min_distance, seq_len, first_6bp)}
        for name in cur_copies_out_contigs.keys():
            raw_name = name.split(' ')[1]
            parts = raw_name.split(':')
            raw_name = ':'.join(parts[:-1])
            pos_parts = parts[-1].split('..')
            cur_start = int(pos_parts[0])
            cur_end = int(pos_parts[1])
            if cur_start > cur_end:
                tmp = cur_start
                cur_start = cur_end
                cur_end = tmp
            cur_seq = cur_copies_out_contigs[name]
            # get the original boundary
            raw_start = int(raw_name.split('-rawstart_')[1].split('-')[0]) + 1
            raw_end = int(raw_name.split('-rawend_')[1])
            cur_distance = abs(cur_start - raw_start) + abs(cur_end - raw_end)
            seq_len = len(cur_seq)
            if not copies_candidate.__contains__(query_name):
                copies_candidate[query_name] = (name, cur_distance, seq_len)
            else:
                last_min_item = copies_candidate[query_name]
                if (cur_distance == last_min_item[1] and seq_len > last_min_item[2]) or cur_distance < last_min_item[1]:
                    copies_candidate[query_name] = (name, cur_distance, seq_len)
        for query_name in copies_candidate.keys():
            item = copies_candidate[query_name]
            candidate_Helitrons[query_name] = cur_copies_out_contigs[item[0]]

    return candidate_Helitrons


def run_polyATail(TRsearch_dir, longest_repeats_multi_line_path, longest_repeats_multi_line_dir):
    polyATail_command = TRsearch_dir + '/polyAtail -l 5 ' + longest_repeats_multi_line_path
    os.system('cd ' + longest_repeats_multi_line_dir + ' && ' + polyATail_command + '> /dev/null 2>&1')
    TR_out = longest_repeats_multi_line_path + '.polyA.set'
    return TR_out

def run_sinefinder(TRsearch_dir, input, input_dir):
    TRsearch_command = TRsearch_dir + '/itrsearch_bak -i 0.75 ' + input
    os.system('cd ' + input_dir + ' && ' +TRsearch_command + '> /dev/null 2>&1')
    TR_out = input + '.itr'
    return TR_out

def run_itrsearch_v1(TRsearch_dir, tir_seq):
    TRsearch_command = 'echo ' + tir_seq + ' | ' + TRsearch_dir + '/itrsearch_bak -i 0.8 -l 5 '
    #output = subprocess.check_output(["echo", tir_seq, "|", TRsearch_dir+"/itrsearch_bak", "-i", "0.7", "-l", "5"], shell=False)
    output = subprocess.check_output(TRsearch_command, shell=True)
    return output.decode('utf-8')

def run_itrsearch(TRsearch_dir, input, input_dir):
    TRsearch_command = TRsearch_dir + '/itrsearch -i 0.7 -l 7 ' + input
    #print(TRsearch_command + "> /dev/null 2>&1")
    if not os.path.exists(input_dir):
        os.makedirs(input_dir)
    TR_log = input + '.log'
    os.system('cd ' + input_dir + ' && ' +TRsearch_command + ' > ' + TR_log)
    TR_out = input + '.itr'
    return TR_out, TR_log

def run_ltrsearch(TRsearch_dir, input, input_dir):
    TRsearch_command = TRsearch_dir + '/ltrsearch -i 0.85 ' + input
    #print(TRsearch_command + "> /dev/null 2>&1")
    os.system('cd ' + input_dir + ' && ' +TRsearch_command + '> /dev/null 2>&1')
    TR_out = input + '.ltr'
    return TR_out

def multi_process_EAHelitron(longest_repeats_flanked_path, flanking_len, output, temp_dir, EAHelitron, threads):
    os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # After partitioning the files, perform parallel computations using multiple processes.
    fasta_file = longest_repeats_flanked_path
    subfile_size = 50000  # 50K
    output_dir = temp_dir
    split_fasta(fasta_file, subfile_size, output_dir)
    split_files = []
    for split_file_name in os.listdir(output_dir):
        split_file = output_dir + '/' + split_file_name
        split_files.append(split_file)

    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    for partition_index, split_file in enumerate(split_files):
        obj = ex.submit(run_EAHelitron, flanking_len, temp_dir, split_file, EAHelitron, partition_index)
        objs.append(obj)
        job_id += 1
    ex.shutdown(wait=True)
    candidate_Helitrons = {}
    for obj in as_completed(objs):
        cur_candidate_Helitrons = obj.result()
        candidate_Helitrons.update(cur_candidate_Helitrons)
    store_fasta(candidate_Helitrons, output)


def multi_process_helitronscanner(candidate_file, output, sh_dir, temp_dir, HSDIR, HSJAR, threads, debug):
    os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # After partitioning the files, perform parallel computations using multiple processes.
    fasta_file = candidate_file
    subfile_size = 50000  # 50K
    output_dir = temp_dir
    split_fasta(fasta_file, subfile_size, output_dir)
    split_files = []
    for split_file_name in os.listdir(output_dir):
        split_file = output_dir + '/' + split_file_name
        split_files.append(split_file)

    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    for partition_index, split_file in enumerate(split_files):
        obj = ex.submit(run_HelitronScanner, sh_dir, temp_dir, split_file, HSDIR, HSJAR, partition_index, debug)
        objs.append(obj)
        job_id += 1
    ex.shutdown(wait=True)
    candidate_Helitrons = {}
    for obj in as_completed(objs):
        cur_candidate_Helitrons = obj.result()
        candidate_Helitrons.update(cur_candidate_Helitrons)
    store_fasta(candidate_Helitrons, output)

def multi_process_TR(input, output, tmp_output_dir, TRsearch_dir):
    threads = 48
    # ---------------------------------------Terminal Repeat Search--------------------------------------------------
    contigNames, contigs = read_fasta(input)
    longest_repeats_cluster = split2cluster_normal(list(contigs.items()), threads)

    longest_repeats_multi_line_dir = tmp_output_dir + '/multi_line_test'
    os.system('rm -rf '+longest_repeats_multi_line_dir)
    if not os.path.exists(longest_repeats_multi_line_dir):
        os.makedirs(longest_repeats_multi_line_dir)

    pool = multiprocessing.Pool(processes=threads)
    for partition_index in longest_repeats_cluster.keys():
        cur_contigs = longest_repeats_cluster[partition_index]
        longest_repeats_multi_line_path = longest_repeats_multi_line_dir + '/'+str(partition_index)+'.fa'

        outfile = open(longest_repeats_multi_line_path, 'w')  # open outfile for writing
        for item in cur_contigs:
            print_seqs(item[0], item[1], 70, outfile)

        pool.apply_async(run_TRsearch, (TRsearch_dir, longest_repeats_multi_line_path, longest_repeats_multi_line_dir,))
    pool.close()
    pool.join()

    final_TR_out = output
    if os.path.exists(final_TR_out):
        os.system('rm -f ' + final_TR_out)
    for partition_index in range(threads):
        cur_TR_out = longest_repeats_multi_line_dir + '/'+str(partition_index) + '.fa' + '.TR.set'
        merge_command = 'cat ' + cur_TR_out + ' >> ' + final_TR_out
        os.system(merge_command)

def multi_process_polyATail(input, output, polyA_temp_dir, TRsearch_dir):
    threads = 48
    # ---------------------------------------Terminal Repeat Search--------------------------------------------------
    contigNames, contigs = read_fasta(input)
    longest_repeats_cluster = split2cluster_normal(list(contigs.items()), threads)

    os.system('rm -rf '+polyA_temp_dir)
    if not os.path.exists(polyA_temp_dir):
        os.makedirs(polyA_temp_dir)

    pool = multiprocessing.Pool(processes=threads)
    for partition_index in longest_repeats_cluster.keys():
        cur_contigs = longest_repeats_cluster[partition_index]
        longest_repeats_multi_line_path = polyA_temp_dir + '/'+str(partition_index)+'.fa'
        new_cur_contigs = {}
        for item in cur_contigs:
            new_cur_contigs[item[0]] = item[1]
        store_fasta(new_cur_contigs, longest_repeats_multi_line_path)
        pool.apply_async(run_polyATail, (TRsearch_dir, longest_repeats_multi_line_path, polyA_temp_dir,))
    pool.close()
    pool.join()

    final_TR_out = output
    if os.path.exists(final_TR_out):
        os.system('rm -f ' + final_TR_out)
    for partition_index in range(threads):
        cur_TR_out = polyA_temp_dir + '/'+str(partition_index) + '.fa' + '.polyA.set'
        merge_command = 'cat ' + cur_TR_out + ' >> ' + final_TR_out
        os.system(merge_command)


def multi_process_ltr(input, output, temp_dir, TRsearch_dir, threads=48):

    # ---------------------------------------Terminal Repeat Search--------------------------------------------------
    contigNames, contigs = read_fasta(input)
    longest_repeats_cluster = split2cluster_normal(list(contigs.items()), threads)

    os.system('rm -rf '+temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    pool = multiprocessing.Pool(processes=threads)
    for partition_index in longest_repeats_cluster.keys():
        cur_contigs = longest_repeats_cluster[partition_index]
        longest_repeats_multi_line_path = temp_dir + '/'+str(partition_index)+'.fa'

        outfile = open(longest_repeats_multi_line_path, 'w')  # open outfile for writing
        for item in cur_contigs:
            print_seqs(item[0], item[1], 70, outfile)

        pool.apply_async(run_ltrsearch, (TRsearch_dir, longest_repeats_multi_line_path, temp_dir,))
    pool.close()
    pool.join()

    final_ltr_out = output
    if os.path.exists(final_ltr_out):
        os.system('rm -f ' + final_ltr_out)
    for partition_index in range(threads):
        cur_TR_out = temp_dir + '/'+str(partition_index) + '.fa' + '.ltr'
        merge_command = 'cat ' + cur_TR_out + ' >> ' + final_ltr_out
        os.system(merge_command)

def multi_process_itr(input, output, longest_repeats_multi_line_dir, TRsearch_dir, threads = 48):
    # ---------------------------------------Terminal Repeat Search--------------------------------------------------
    contigNames, contigs = read_fasta(input)
    longest_repeats_cluster = split2cluster_normal(list(contigs.items()), threads)

    os.system('rm -rf '+longest_repeats_multi_line_dir)
    if not os.path.exists(longest_repeats_multi_line_dir):
        os.makedirs(longest_repeats_multi_line_dir)

    pool = multiprocessing.Pool(processes=threads)
    for partition_index in longest_repeats_cluster.keys():
        cur_contigs = longest_repeats_cluster[partition_index]
        longest_repeats_multi_line_path = longest_repeats_multi_line_dir + '/'+str(partition_index)+'.fa'

        outfile = open(longest_repeats_multi_line_path, 'w')  # open outfile for writing
        for item in cur_contigs:
            print_seqs(item[0], item[1], 70, outfile)

        pool.apply_async(run_itrsearch, (TRsearch_dir, longest_repeats_multi_line_path, longest_repeats_multi_line_dir,))
    pool.close()
    pool.join()

    final_itr_out = output
    if os.path.exists(final_itr_out):
        os.system('rm -f ' + final_itr_out)
    for partition_index in range(threads):
        cur_TR_out = longest_repeats_multi_line_dir + '/'+str(partition_index) + '.fa' + '.itr'
        merge_command = 'cat ' + cur_TR_out + ' >> ' + final_itr_out
        os.system(merge_command)

def multi_process_sinefinder(input, output, temp_dir, TRsearch_dir, threads = 48):
    # ---------------------------------------Terminal Repeat Search--------------------------------------------------
    contigNames, contigs = read_fasta(input)
    longest_repeats_cluster = split2cluster_normal(list(contigs.items()), threads)

    os.system('rm -rf '+temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    pool = multiprocessing.Pool(processes=threads)
    for partition_index in longest_repeats_cluster.keys():
        cur_contigs = longest_repeats_cluster[partition_index]
        longest_repeats_multi_line_path = temp_dir + '/'+str(partition_index)+'.fa'

        outfile = open(longest_repeats_multi_line_path, 'w')  # open outfile for writing
        for item in cur_contigs:
            print_seqs(item[0], item[1], 70, outfile)

        pool.apply_async(run_itrsearch, (TRsearch_dir, longest_repeats_multi_line_path, temp_dir,))
    pool.close()
    pool.join()

    final_itr_out = output
    if os.path.exists(final_itr_out):
        os.system('rm -f ' + final_itr_out)
    for partition_index in range(threads):
        cur_TR_out = temp_dir + '/'+str(partition_index) + '.fa' + '.itr'
        merge_command = 'cat ' + cur_TR_out + ' >> ' + final_itr_out
        os.system(merge_command)

def store_LTR_seq(ltrharvest_output, longest_repeats_path, confident_ltr_path, confident_ltr_cut_path):
    longest_repeats_contigNames, longest_repeats_contigs = read_fasta(longest_repeats_path)
    LTR_seqs = {}
    LTR_intact_seqs = {}
    node_index = 0
    with open(ltrharvest_output, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            if line.startswith('#'):
                continue
            parts = line.split('  ')
            lLTR_start = int(parts[3])
            lLTR_end = int(parts[4])
            lLTR_len = int(parts[5])
            left_tsd = parts[6]
            left_tsd_len = int(parts[7])
            rLTR_start = int(parts[9])
            rLTR_end = int(parts[10])
            rLTR_len = int(parts[11])
            right_tsd = parts[12]
            right_tsd_len = int(parts[13])
            LTR_similarity = float(parts[15])
            seq_id = int(parts[16])
            query_name = longest_repeats_contigNames[seq_id]
            query_seq = longest_repeats_contigs[query_name]
            if lLTR_len >= rLTR_len:
                LTR = query_seq[lLTR_start-1: lLTR_end]
            else:
                LTR = query_seq[rLTR_start - 1: rLTR_end]
            LTR_internal = query_seq[lLTR_end: rLTR_start - 1]

            LTR_query_name = 'N_' + str(node_index) + '-LTR' + \
                             '-lLTRStart_' + str(1) + '-lLTREnd_' + str(lLTR_end-lLTR_start+1) +\
                             '-rLTRStart_'+str(rLTR_start-lLTR_start+1)+'-rLTREnd_'+str(rLTR_end-lLTR_start+1) + '-tsd_' + left_tsd
            internal_query_name = 'N_' + str(node_index) + '-ILTR'
            LTR_seqs[LTR_query_name] = LTR
            LTR_seqs[internal_query_name] = LTR_internal
            LTR_intact_seqs[LTR_query_name] = query_seq[lLTR_start-1: rLTR_end]
            node_index += 1
    f_r.close()

    store_fasta(LTR_seqs, confident_ltr_cut_path)
    store_fasta(LTR_intact_seqs, confident_ltr_path)

def store_LTR_seq_v2(ltrharvest_output, longest_repeats_path, confident_ltr_path, confident_ltr_cut_path):
    longest_repeats_contigNames, longest_repeats_contigs = read_fasta(longest_repeats_path)
    LTR_seqs = {}
    LTR_intact_seqs = {}
    node_index = 0
    with open(ltrharvest_output, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            if line.startswith('#'):
                continue
            parts = line.split('  ')
            lLTR_start = int(parts[3])
            lLTR_end = int(parts[4])
            lLTR_len = int(parts[5])
            rLTR_start = int(parts[6])
            rLTR_end = int(parts[7])
            rLTR_len = int(parts[8])
            LTR_similarity = float(parts[9])
            seq_id = int(parts[10])
            query_name = longest_repeats_contigNames[seq_id]
            query_seq = longest_repeats_contigs[query_name]
            if lLTR_len >= rLTR_len:
                LTR = query_seq[lLTR_start-1: lLTR_end]
            else:
                LTR = query_seq[rLTR_start - 1: rLTR_end]
            LTR_internal = query_seq[lLTR_end: rLTR_start - 1]

            LTR_query_name = 'N_' + str(node_index) + '-LTR' + \
                             '-lLTRStart_' + str(1) + '-lLTREnd_' + str(lLTR_end-lLTR_start+1) +\
                             '-rLTRStart_'+str(rLTR_start-lLTR_start+1)+'-rLTREnd_'+str(rLTR_end-lLTR_start+1)
            internal_query_name = 'N_' + str(node_index) + '-ILTR'
            LTR_seqs[LTR_query_name] = LTR
            LTR_seqs[internal_query_name] = LTR_internal
            LTR_intact_seqs[LTR_query_name] = query_seq[lLTR_start-1: rLTR_end]
            node_index += 1
    f_r.close()

    store_fasta(LTR_seqs, confident_ltr_cut_path)
    store_fasta(LTR_intact_seqs, confident_ltr_path)

def store_LTR_seq_v1(ltrharvest_output, longest_repeats_path, confident_ltr_path, confident_ltr_cut_path):
    longest_repeats_contigNames, longest_repeats_contigs = read_fasta(longest_repeats_path)
    LTR_seqs = {}
    LTR_intact_seqs = {}
    node_index = 0
    with open(ltrharvest_output, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            if line.startswith('#'):
                continue
            parts = line.split(' ')
            lLTR_start = int(parts[3])
            lLTR_end = int(parts[4])
            lLTR_len = int(parts[5])
            rLTR_start = int(parts[6])
            rLTR_end = int(parts[7])
            rLTR_len = int(parts[8])
            seq_id = int(parts[10])
            query_name = longest_repeats_contigNames[seq_id]
            query_seq = longest_repeats_contigs[query_name]
            if lLTR_len >= rLTR_len:
                LTR = query_seq[lLTR_start-1: lLTR_end]
            else:
                LTR = query_seq[rLTR_start - 1: rLTR_end]
            LTR_internal = query_seq[lLTR_end: rLTR_start - 1]

            LTR_query_name = 'N_' + str(node_index) + '-LTR' + \
                             '-lLTRStart_' + str(1) + '-lLTREnd_' + str(lLTR_end-lLTR_start+1) +\
                             '-rLTRStart_'+str(rLTR_start-lLTR_start+1)+'-rLTREnd_'+str(rLTR_end-lLTR_start+1)
            internal_query_name = 'N_' + str(node_index) + '-ILTR'
            LTR_seqs[LTR_query_name] = LTR
            LTR_seqs[internal_query_name] = LTR_internal
            LTR_intact_seqs[LTR_query_name] = query_seq[lLTR_start-1: rLTR_end]
            node_index += 1
    f_r.close()

    store_fasta(LTR_seqs, confident_ltr_cut_path)
    store_fasta(LTR_intact_seqs, confident_ltr_path)

def run_LTR_detection(reference, tmp_output_dir, threads, LTR_harvest_parallel_Home, LTR_finder_parallel_Home, log):
    # 1.run LTR_harvest_parallel
    ltrharvest_output = reference + '.harvest.combine.scn'
    if os.path.isfile(ltrharvest_output):
        os.remove(ltrharvest_output)

    LTR_harvest_parallel_command = 'perl ' + LTR_harvest_parallel_Home + '/LTR_HARVEST_parallel -seq ' + reference + ' -threads ' + str(threads)
    log.logger.debug('cd ' + tmp_output_dir + ' && ' + LTR_harvest_parallel_command + ' > /dev/null 2>&1')
    os.system('cd ' + tmp_output_dir + ' && ' + LTR_harvest_parallel_command + ' > /dev/null 2>&1')

    # 1.run LTR_finder_parallel
    ltrfinder_output = reference + '.finder.combine.scn'
    if os.path.isfile(ltrfinder_output):
        os.remove(ltrfinder_output)
    LTR_finder_parallel_command = 'perl ' + LTR_finder_parallel_Home + '/LTR_FINDER_parallel -harvest_out -seq ' + reference + ' -threads ' + str(threads)
    log.logger.debug('cd ' + tmp_output_dir + ' && ' + LTR_finder_parallel_command + ' > /dev/null 2>&1')
    os.system('cd ' + tmp_output_dir + ' && ' + LTR_finder_parallel_command + ' > /dev/null 2>&1')

def run_LTR_harvest(reference, tmp_output_dir, threads, LTR_finder_parallel_Home, log):
    # starttime = time.time()
    # log.logger.debug('start LTR_harvest detection...')
    cut_references = []

    # If the genome exceeds 4G, cut the genome into 1G
    # to obtain the file size (in bytes).
    file_size = os.path.getsize(reference)
    one_g = 1024 ** 3
    file_size_gb = file_size / one_g
    if file_size_gb > 4:
        cur_ref_contigs = {}
        cur_base_num = 0
        ref_index = 0
        (ref_dir, ref_filename) = os.path.split(reference)
        (ref_name, ref_extension) = os.path.splitext(ref_filename)

        ref_names, ref_contigs = read_fasta(reference)
        for name in ref_names:
            seq = ref_contigs[name]
            cur_ref_contigs[name] = seq
            cur_base_num += len(seq)
            if cur_base_num >= one_g:
                # store references
                    cur_ref_path = tmp_output_dir + '/' +  ref_filename + '.ltr_cut' + str(ref_index) + '.fa'
                    store_fasta(cur_ref_contigs, cur_ref_path)
                    cut_references.append(cur_ref_path)
                    cur_ref_contigs = {}
                    cur_base_num = 0
                    ref_index += 1
        if len(cur_ref_contigs) > 0:
                cur_ref_path = cur_ref_path = tmp_output_dir + '/' +  ref_filename + '.ltr_cut' + str(ref_index) + '.fa'
                store_fasta(cur_ref_contigs, cur_ref_path)
                cut_references.append(cur_ref_path)
    else:
        cut_references.append(reference)

    if len(cut_references) < threads:
        thread_size = len(cut_references)
    else:
        thread_size = threads
    ex = ProcessPoolExecutor(thread_size)
    jobs = []
    for ref_index, cut_reference in enumerate(cut_references):
        job = ex.submit(run_LTR_harvest_single, cut_reference, tmp_output_dir, ref_index)
        jobs.append(job)
    ex.shutdown(wait=True)

    output = tmp_output_dir + '/genome_all.fa.harvest.scn'
    if os.path.isfile(output):
        os.remove(output)
    for job in as_completed(jobs):
        cur_output = job.result()
        os.system('cat ' + cur_output + ' >> ' + output)

    # run LTR_finder_parallel
    ltrfinder_output = reference + '.finder.combine.scn'
    if os.path.isfile(ltrfinder_output):
        os.remove(ltrfinder_output)
    for cut_reference in cut_references:
        LTR_finder_parallel_command = 'perl ' + LTR_finder_parallel_Home + '/LTR_FINDER_parallel -harvest_out -seq ' + cut_reference + ' -threads ' + str(threads)
        log.logger.debug('cd ' + tmp_output_dir + ' && ' + LTR_finder_parallel_command + ' > /dev/null 2>&1')
        os.system('cd ' + tmp_output_dir + ' && ' + LTR_finder_parallel_command + ' > /dev/null 2>&1')
        if len(cut_references) > 1:
            cur_ltrfinder_output = cut_reference + '.finder.combine.scn'
            os.system('cat ' + cur_ltrfinder_output + ' | grep -Ev \'^$|#\' >> ' + ltrfinder_output)

    # endtime = time.time()
    # dtime = endtime - starttime
    # log.logger.debug("LTR_harvest running time: %.8s s" % (dtime))

def run_LTR_harvest_single(reference, tmp_output_dir, ref_index):
    output = tmp_output_dir + '/genome_'+str(ref_index)+'.fa.harvest.scn'
    ltrharvest_command1 = 'gt suffixerator -db ' + reference + ' -indexname ' \
                          + reference + ' -tis -suf -lcp -des -ssp -sds -dna'
    ltrharvest_command2 = 'gt ltrharvest -index ' + reference \
                          + ' -seed 20 -minlenltr 100 -maxlenltr 7000 -similar 85 -motif TGCA -motifmis 1 -mintsd 4 -maxtsd 6 ' \
                            '-vic 10  > ' + output

    os.system(ltrharvest_command1)
    os.system(ltrharvest_command2)

    # Change the last column to chromosome name
    ref_names, ref_contigs = read_fasta(reference)
    new_lines = []
    with open(output, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            else:
                line = line.replace('\n', '')
                parts = line.split('  ')
                new_line = ''
                for i, p in enumerate(parts):
                    if i == len(parts)-1:
                        ref_index = int(p)
                        p = ref_names[ref_index]
                    else:
                        p += '  '
                    new_line += p
                new_lines.append(new_line)
    f_r.close()
    with open(output, 'w') as f_save:
        for line in new_lines:
            f_save.write(line+'\n')
    f_save.close
    return output

def run_LTR_retriever(reference, tmp_output_dir, threads, miu, log):
    starttime = time.time()
    log.logger.debug('start LTR_retriever detection...')
    LTR_retriever_command = 'cd ' + tmp_output_dir + ' && LTR_retriever -genome ' + reference \
                            + ' -inharvest ' + tmp_output_dir + '/genome_all.fa.rawLTR.scn -noanno -threads ' + str(threads) + ' -u ' + str(miu)
    log.logger.debug(LTR_retriever_command)
    os.system(LTR_retriever_command)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("LTR_retriever running time: %.8s s" % (dtime))

def run_GRF(GRF_Home, reference, tmp_output_dir, threads):
    grf_tir_command = GRF_Home + '/bin/grf-main -i ' + reference + ' -o ' + tmp_output_dir + ' -c 0 --min_tr 10 -t ' + str(threads)
    os.system(grf_tir_command)

    grf_mite_command = 'sh ' + GRF_Home + '/script/run_mite_detection.sh ' + reference + ' ' + tmp_output_dir + ' ' + str(threads)
    os.system(grf_mite_command)


def endWithPolyA(query_seq, query_start, query_end):
    prev_base = ''

    cur_polyA_len = 0
    cur_max_polyA_start = -1
    cur_max_polyA_end = -1

    max_polyA_len = 0
    max_polyA_start = -1
    max_polyA_end = -1
    for i in range(query_end, len(query_seq)):
        cur_base = query_seq[i]
        if cur_base == 'A':
            if prev_base != 'A':
                cur_max_polyA_start = i
            cur_max_polyA_end = i

            cur_polyA_len += 1
            if cur_polyA_len > max_polyA_len:
                max_polyA_len = cur_polyA_len
                max_polyA_start = cur_max_polyA_start
                max_polyA_end = cur_max_polyA_end
        else:
            cur_polyA_len = 0
            cur_max_polyA_start = -1
            cur_max_polyA_end = -1
        prev_base = cur_base

    if max_polyA_len >= 4:
        query_seq = query_seq[:max_polyA_end + 1]
    return query_seq

def get_longest_protein(LTR_out):
    query_records = {}
    with open(LTR_out, 'r') as f_r:
        for idx, line in enumerate(f_r):
            line = str(line).replace('\n', '')
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            identity = float(parts[2])
            alignment_len = int(parts[3])
            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            if query_name == subject_name or identity < 60:
                continue
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))
    f_r.close()

    fixed_extend_base_threshold = 200
    keep_longest_query = {}
    for idx, query_name in enumerate(query_records.keys()):
        subject_dict = query_records[query_name]

        longest_queries = []
        for subject_name in subject_dict.keys():
            subject_pos = subject_dict[subject_name]

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
                        if (frag[2] - exist_frag[3] < fixed_extend_base_threshold):
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
                        if (exist_frag[3] - frag[2] < fixed_extend_base_threshold):
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
                cur_cluster.sort(key=lambda x: (x[2], x[3]))

                cluster_longest_query_start = -1
                cluster_longest_query_end = -1
                cluster_longest_query_len = -1

                cluster_longest_subject_start = -1
                cluster_longest_subject_end = -1
                cluster_longest_subject_len = -1

                # record visited fragments
                visited_frag = {}
                for i in range(len(cur_cluster)):
                    # keep a longest query start from each fragment
                    origin_frag = cur_cluster[i]
                    if visited_frag.__contains__(origin_frag):
                        continue
                    cur_frag_len = abs(origin_frag[1] - origin_frag[0])
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
                        if ext_frag[3] > longest_subject_end:
                            # judge subject direction
                            if longest_query_start < longest_query_end and ext_frag[0] < ext_frag[1]:
                                # +
                                if ext_frag[1] > longest_query_end:
                                    # forward extend
                                    if ext_frag[0] - longest_query_end < fixed_extend_base_threshold and ext_frag[
                                        2] - longest_subject_end < fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start

                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= fixed_extend_base_threshold:
                                        break
                            elif longest_query_start > longest_query_end and ext_frag[0] > ext_frag[1]:
                                # reverse
                                if ext_frag[1] < longest_query_end:
                                    # reverse extend
                                    if longest_query_end - ext_frag[0] < fixed_extend_base_threshold and ext_frag[2] - longest_subject_end < fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_start - longest_query_end

                                        visited_frag[ext_frag] = 1
                                    elif longest_query_end - ext_frag[0] >= fixed_extend_base_threshold:
                                        break
                    if cur_longest_query_len > cluster_longest_query_len:
                        cluster_longest_query_start = longest_query_start
                        cluster_longest_query_end = longest_query_end
                        cluster_longest_query_len = cur_longest_query_len

                        cluster_longest_subject_start = longest_subject_start
                        cluster_longest_subject_end = longest_subject_end
                        cluster_longest_subject_len = longest_subject_end - longest_subject_start

                # keep this longest query
                if cluster_longest_query_len != -1:
                    longest_queries.append((cluster_longest_query_start, cluster_longest_query_end,
                                            cluster_longest_query_len, cluster_longest_subject_start,
                                            cluster_longest_subject_end, cluster_longest_subject_len, subject_name))

        longest_queries.sort(key=lambda x: -x[5])
        keep_longest_query[query_name] = longest_queries
    return keep_longest_query

def get_LINE_candidate(LINE_out, protein_path, orig_contigs):
    candidate_LINE = {}
    protein_names, protein_contigs = read_fasta(protein_path)
    keep_longest_query = get_longest_protein(LINE_out)
    node_index = 0
    for query_name in keep_longest_query.keys():
        protein_items = keep_longest_query[query_name]

        new_query_name = 'N_' + str(node_index)
        query_seq = orig_contigs[query_name]

        candidate_protein_items = {}
        last_query_pos = -1
        longest_pol_item = None
        for protein_item in protein_items:
            protein_len = protein_item[5]
            protein_name = protein_item[6]
            if str(protein_name).__contains__('gagpol'):
                type = 'gagpol'
            elif str(protein_name).__contains__('pol'):
                type = 'pol'
            elif str(protein_name).__contains__('gag'):
                type = 'gag'
            else:
                type = 'other'
            intact_protein_len = len(protein_contigs[protein_name])
            ratio_protein = float(protein_len) / intact_protein_len

            if ratio_protein >= 0.8:
                if not candidate_protein_items.__contains__(type):
                    candidate_protein_items[type] = (protein_len, protein_item)
                    if type == 'gagpol' or type == 'pol':
                        longest_pol_item = protein_item
                else:
                    orig_protein_len = candidate_protein_items[type][0]
                    if protein_len > orig_protein_len:
                        candidate_protein_items[type] = (protein_len, protein_item)
                        if type == 'gagpol' or type == 'pol':
                            longest_pol_item = protein_item

        if len(candidate_protein_items) > 0:
            # Judge the direction of the last protein
            if longest_pol_item is None:
                continue
            last_query_start = longest_pol_item[0]
            last_query_end = longest_pol_item[1]
            last_direct = '+'
            if last_query_start > last_query_end:
                last_direct = '-'
            if last_direct == '-':
                query_seq = getReverseSequence(query_seq)
                last_query_start = len(query_seq) - last_query_start + 1
                last_query_end = len(query_seq) - last_query_end + 1

            for type in candidate_protein_items.keys():
                protein_item = candidate_protein_items[type][1]
                protein_name = protein_item[6]

                query_start = protein_item[0]
                query_end = protein_item[1]

                direct = '+'
                if query_start > query_end:
                    direct = '-'

                # If the sequences have been reversely complementary, the original positions should also be interchanged.
                if last_direct == '-':
                    query_start = len(query_seq) - query_start + 1
                    query_end = len(query_seq) - query_end + 1

                new_query_name += ':' + protein_name + '_' + str(query_start) + '_' + str(query_end)

            query_seq = endWithPolyA(query_seq, last_query_start, last_query_end)
            new_query_name += '-len_'+str(len(query_seq))
            candidate_LINE[new_query_name] = query_seq
            node_index += 1
    return candidate_LINE

def multiple_alignment_blastx_v1(repeats_path, merge_distance):
    split_repeats_path = repeats_path[0]
    protein_db_path = repeats_path[1]
    blastx2Results_path = repeats_path[2]
    cur_table = repeats_path[3]
    align_command = 'blastx -db ' + protein_db_path + ' -num_threads ' \
                    + str(1) + ' -evalue 1e-20 -query ' + split_repeats_path + ' -outfmt 6 > ' + blastx2Results_path
    os.system(align_command)
    
    fixed_extend_base_threshold = merge_distance
    # Combine the segmented blastx alignments.
    query_names, query_contigs = read_fasta(split_repeats_path)

    # parse blastn output, determine the repeat boundary
    # query_records = {query_name: {subject_name: [(q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end)] }}
    query_records = {}
    with open(blastx2Results_path, 'r') as f_r:
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
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))
    f_r.close()

    keep_longest_query = {}
    longest_repeats = {}
    for idx, query_name in enumerate(query_records.keys()):
        query_len = len(query_contigs[query_name])
        #print('total query size: %d, current query name: %s, idx: %d' % (len(query_records), query_name, idx))

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
                if pos_item[0] > pos_item[1]:
                    reverse_pos.append(pos_item)
                else:
                    forward_pos.append(pos_item)
            forward_pos.sort(key=lambda x: (x[2], x[3]))
            reverse_pos.sort(key=lambda x: (-x[0], -x[1]))

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
                        if (frag[0] - exist_frag[1] < fixed_extend_base_threshold):
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
                        if (exist_frag[1] - frag[0] < fixed_extend_base_threshold):
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
                cur_cluster.sort(key=lambda x: (x[2], x[3]))

                cluster_longest_query_start = -1
                cluster_longest_query_end = -1
                cluster_longest_query_len = -1

                cluster_longest_subject_start = -1
                cluster_longest_subject_end = -1
                cluster_longest_subject_len = -1

                cluster_extend_num = 0

                # print('subject pos size: %d' %(len(cur_cluster)))
                # record visited fragments
                visited_frag = {}
                for i in range(len(cur_cluster)):
                    # keep a longest query start from each fragment
                    origin_frag = cur_cluster[i]
                    if visited_frag.__contains__(origin_frag):
                        continue
                    cur_frag_len = abs(origin_frag[1] - origin_frag[0])
                    cur_longest_query_len = cur_frag_len
                    longest_query_start = origin_frag[0]
                    longest_query_end = origin_frag[1]
                    longest_subject_start = origin_frag[2]
                    longest_subject_end = origin_frag[3]

                    cur_extend_num = 0

                    visited_frag[origin_frag] = 1
                    # try to extend query
                    for j in range(i + 1, len(cur_cluster)):
                        ext_frag = cur_cluster[j]
                        if visited_frag.__contains__(ext_frag):
                            continue

                        # could extend
                        # extend right
                        if ext_frag[3] > longest_subject_end:
                            # judge query direction
                            if longest_query_start < longest_query_end and ext_frag[0] < ext_frag[1]:
                                # +
                                if ext_frag[1] > longest_query_end:
                                    # forward extend
                                    if ext_frag[0] - longest_query_end < fixed_extend_base_threshold and ext_frag[
                                        2] - longest_subject_end < fixed_extend_base_threshold/3:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= fixed_extend_base_threshold:
                                        break
                            elif longest_query_start > longest_query_end and ext_frag[0] > ext_frag[1]:
                                # reverse
                                if ext_frag[1] < longest_query_end:
                                    # reverse extend
                                    if  longest_query_end - ext_frag[0] < fixed_extend_base_threshold and ext_frag[
                                        2] - longest_subject_end < fixed_extend_base_threshold/3:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_start - longest_query_end
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif longest_query_end -  ext_frag[0] >= fixed_extend_base_threshold:
                                        break
                    if cur_longest_query_len > cluster_longest_query_len:
                        cluster_longest_query_start = longest_query_start
                        cluster_longest_query_end = longest_query_end
                        cluster_longest_query_len = cur_longest_query_len

                        cluster_longest_subject_start = longest_subject_start
                        cluster_longest_subject_end = longest_subject_end
                        cluster_longest_subject_len = longest_subject_end - longest_subject_start

                        cluster_extend_num = cur_extend_num
                # keep this longest query
                if cluster_longest_query_len != -1:
                    longest_queries.append((cluster_longest_query_start, cluster_longest_query_end,
                                            cluster_longest_query_len, cluster_longest_subject_start,
                                            cluster_longest_subject_end, cluster_longest_subject_len, subject_name, cluster_extend_num))


        # we now consider, we should take some sequences from longest_queries to represent this query sequence.
        # we take the longest sequence by length, if the latter sequence overlap with the former sequence largely (50%),
        # continue find next sequence until the ratio of query sequence over 90% or no more sequences.
        longest_queries.sort(key=lambda x: -x[2])
        keep_longest_query[query_name] = longest_queries
    #print(keep_longest_query)
    with open(cur_table, 'w') as f_save:
        for query_name in keep_longest_query.keys():
            domain_array = keep_longest_query[query_name]
            # for domain_info in domain_array:
            #     f_save.write(query_name+'\t'+str(domain_info[6])+'\t'+str(domain_info[0])+'\t'+str(domain_info[1])+'\t'+str(domain_info[3])+'\t'+str(domain_info[4])+'\n')
            merge_domains = []
            domain_array.sort(key=lambda x: -x[2])
            for domain_info in domain_array:
                if len(merge_domains) == 0:
                    merge_domains.append(domain_info)
                else:
                    is_new_domain = True
                    for pre_domain in merge_domains:
                        pre_start = pre_domain[0]
                        pre_end = pre_domain[1]
                        #计算overlap
                        if pre_start > pre_end:
                            tmp = pre_start
                            pre_start = pre_end
                            pre_end = tmp
                        cur_start = domain_info[0]
                        cur_end = domain_info[1]
                        if cur_start > cur_end:
                            tmp = cur_start
                            cur_start = cur_end
                            cur_end = tmp
                        if cur_end >= pre_start and cur_end <= pre_end:
                            if cur_start <= pre_start:
                                overlap = cur_end - pre_start
                            else:
                                overlap = cur_end - cur_start
                        elif cur_end > pre_end:
                            if cur_start >= pre_start and cur_start <= pre_end:
                                overlap = pre_end - cur_start
                            else:
                                overlap = 0
                        else:
                            overlap = 0
                        
                        if float(overlap / domain_info[2]) > 0.5:
                            is_new_domain = False
                    if  is_new_domain:
                        merge_domains.append(domain_info)

            for domain_info in merge_domains:
                f_save.write(query_name+'\t'+str(domain_info[6])+'\t'+str(domain_info[0])+'\t'+str(domain_info[1])+'\t'+str(domain_info[3])+'\t'+str(domain_info[4])+'\n')

    f_save.close()
    return cur_table

def multiple_alignment_blastx(repeats_path, tools_dir):
    split_repeats_path = repeats_path[0]
    protein_db_path = repeats_path[1]
    temp_dir = repeats_path[2]
    blastx2Results_path = temp_dir + '/temp.out'
    align_command = 'blastx -db ' + protein_db_path + ' -num_threads ' \
                    + str(1) + ' -query ' + split_repeats_path + ' -word_size 7 -outfmt 6 > ' + blastx2Results_path
    os.system(align_command)

    # 1. First, align longest_repeats.fa to LINERep.fa to identify the complete pol.
    # 2. Determine the direction, if it's forward, search for polyA at the tail; if it's reverse, search for polyT at the head.
    # 3. After determining the tail, take a total of 28 kmers from 30-3bp after the tail. Then, extend the head (referring to the initial domain, if it's only pol, it's the pol position, if it's gag, it's the gag position) by 500bp.
    # 4. Then, align these 28 kmers to the extended sequence. From longest to shortest, take the best-matched position as the starting position of the TSD at the head.
    candidate_LINE = get_LINE_candidate(blastx2Results_path, protein_db_path, split_repeats_path)
    cur_line_path = temp_dir + '/LINE.fa'
    store_fasta(candidate_LINE, cur_line_path)

    return cur_line_path

def getUniqueKmer_v1(cur_segments, partiton_index):
    unique_kmers = []
    for line in cur_segments:
        kmer = line.split(' ')[0]
        r_kmer = getReverseSequence(kmer)
        #unique_key = kmer if kmer < r_kmer else r_kmer
        unique_kmers.append(kmer)
        unique_kmers.append(r_kmer)
    return unique_kmers

def getCombineFragments(region_combination_item, frag_hash, identity_threshold, length_similarity_cutoff, refContigs, region_dict, output_dir, blast_program_dir, partiton_index):
    # go through each region, find candidate combine fragments
    # regionContigs keeps all fragments in each region
    regionContigs = {}
    # go through each region
    region_id = region_combination_item[0]
    cur_region_combination = region_combination_item[1]
    combinations = cur_region_combination['combinations']
    max_combination_len = cur_region_combination['max_combination_len']
    # start from max length combination
    for c in range(max_combination_len, 0, -1):
        cur_combinations = combinations[str(c)]
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
                    identity = compare_seq(self_info, other_info, identity_threshold, length_similarity_cutoff,
                                           refContigs, output_dir, blast_program_dir, partiton_index)
                    if identity is not None and identity > max_identity:
                        max_identity = identity
                        best_combine_name = combine_name
        # if current combination reach score threshold, then it can be used to replace the whole region
        if max_identity >= identity_threshold:
            if not regionContigs.__contains__(region_id):
                regionContigs[region_id] = []
            final_frags = regionContigs[region_id]
            ref_name = self_info[1]
            # replace the whole region with best combine
            final_frags.append((best_combine_name, ref_name, self_info[2], self_info[3]))
            # all_frags = [(F1, start, end),(F2, start, end),(F3, start, end),(F4, start, end)]
            all_frags = region_dict[ref_name]
            best_frags = best_combine_name.split(',')
            # keep other fragments
            for frag in all_frags:
                if frag[0] not in best_frags:
                    final_frags.append((frag[0], ref_name, frag[1], frag[2]))
            regionContigs[region_id] = final_frags
            break
    return regionContigs

def getRegionCombination(region_item):
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

    region_id = region_item[0]
    cur_region_dict = region_item[1]
    for ref_name in cur_region_dict.keys():
        cur_region_list = cur_region_dict[ref_name]
        max_combination_len = len(cur_region_list)
        print('current region id: ' + str(region_id) + ', size: ' + str(max_combination_len))
        if not region_combination.__contains__(region_id):
            region_combination[region_id] = {}
        cur_region_combination = region_combination[region_id]
        cur_region_combination['max_combination_len'] = max_combination_len
        combinations = {}
        for c in range(1, max_combination_len + 1):
            if not combinations.__contains__(c):
                combinations[str(c)] = []
            cur_combinations = combinations[str(c)]
            for left in range(len(cur_region_list) - c + 1):
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
                    if cur_frag_end > combine_frag_end:
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
            combinations[str(c)] = cur_combinations
        cur_region_combination['combinations'] = combinations
    return region_combination, frag_hash


def compare_seq(self_info, other_info, identity_cutoff, length_similarity_cutoff,
                refContigs, output_dir, blast_program_dir, partiton_index):
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
    output_dir += '/blastn_tmp'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    self_seq_path = output_dir + '/self_' + str(partiton_index) + '.fa'
    other_seq_path = output_dir + '/other_' + str(partiton_index) + '.fa'
    blastnResults_path = output_dir + '/blast_' + str(partiton_index) + '.out'
    store_fasta(self_contigs, self_seq_path)
    store_fasta(other_contigs, other_seq_path)

    makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + other_seq_path
    align_command = blast_program_dir + '/bin/blastn -db ' + other_seq_path + ' -query ' + self_seq_path + ' -outfmt 6 > ' + blastnResults_path
    print(makedb_command)
    os.system(makedb_command)
    print(align_command)
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
            if identity >= identity_cutoff and \
                    float(match_base) / query_len >= length_similarity_cutoff and\
                    float(match_base) / target_len >= length_similarity_cutoff:
                return float(identity) / 100
    f_r.close()


def multi_line(fasta_path, line_len, k_num):
    k_num = int(k_num)
    tmp_fasta_path = fasta_path + ".tmp"
    contigNames, contigs = read_fasta(fasta_path)
    with open(tmp_fasta_path, 'w') as f_w:
        for contigName in contigNames:
            contig = contigs[contigName]
            # line = '>' + contigName + '\t' + contig + '\n'
            # f_w.write(line)
            start = 0
            end = len(contig) - k_num + 1
            while start < end:
                # add extra kmer length
                cur_end = start+line_len
                cur_end = cur_end if cur_end <= len(contig) else len(contig)
                seg = contig[start:cur_end]
                line = '>' + contigName + '\t' + str(start) + '\t' + seg + '\n'
                f_w.write(line)
                start += line_len
    f_w.close()
    return tmp_fasta_path

def convertToUpperCase(reference):
    cur_segments = []
    contigNames = []
    contigs = {}
    with open(reference, "r") as f_r:
        contigName = ''
        contigseq = ''
        for line in f_r:
            if line.startswith('>'):
                if contigName != '' and contigseq != '':
                    contigs[contigName] = contigseq
                    contigNames.append(contigName)
                    cur_segments.append(contigseq)
                contigName = line.strip()[1:].split(' ')[0]
                contigseq = ''
            else:
                contigseq += line.strip().upper()
        contigs[contigName] = contigseq
        contigNames.append(contigName)
        cur_segments.append(contigseq)
    f_r.close()
    return contigs

def convertToUpperCase_v1(reference):
    contigNames = []
    contigs = {}
    with open(reference, "r") as f_r:
        contigName = ''
        contigseq = ''
        for line in f_r:
            if line.startswith('>'):
                if contigName != '' and contigseq != '':
                    contigs[contigName] = contigseq
                    contigNames.append(contigName)
                contigName = line.strip()[1:].split(' ')[0]
                contigseq = ''
            else:
                contigseq += line.strip().upper()
        contigs[contigName] = contigseq
        contigNames.append(contigName)
    f_r.close()

    # (dir, filename) = os.path.split(reference)
    # (name, extension) = os.path.splitext(filename)
    # reference_pre = dir + '/' + name + '_preprocess' + extension
    with open(reference, "w") as f_save:
        for contigName in contigNames:
            contigseq = contigs[contigName]
            f_save.write(">" + contigName + '\n' + contigseq + '\n')
    f_save.close()
    return reference

def generate_candidate_repeats_v2(contigs, k_num, unique_kmer_map, partiton_index, fault_tolerant_bases):

    #print('partition_index: %d, total contigs: %d' %(partiton_index, len(contigs)))

    # file = open(unique_kmer_map_file, 'r')
    # js = file.read()
    # unique_kmer_map = json.loads(js)

    cur_masked_segments = {}
    for ref_name in contigs.keys():
        line = contigs[ref_name]
        masked_line = list(line)
        last_masked_pos = -1
        for i in range(len(line)-k_num+1):
            kmer = line[i: i+k_num]
            # get reverse complement kmer
            #r_kmer = getReverseSequence(kmer)
            # filter invalid kmer, contains 'N'
            if "N" in kmer:
                continue
            #unique_key = kmer if kmer < r_kmer else r_kmer

            if unique_kmer_map.__contains__(kmer):
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

    #print('partition_index: %d finish masking stage' % (partiton_index))

    repeat_dict = {}
    cur_repeat_str = ''
    try_connect_str = ''
    last_start_pos = -1
    last_end_pos = -1
    for seq_index, cur_masked_item in enumerate(cur_masked_segments.items()):
        ref_name = cur_masked_item[0]
        ref_seq = contigs[ref_name]
        #print('ref seq length: %d' %len(ref_seq))
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
                cur_repeat_str = cur_repeat_str + ref_seq[i]
                last_end_pos = i
            elif cur_repeat_str != '':
                # meet unmasked base
                if (i - last_end_pos) <= fault_tolerant_bases:
                    # skip gap
                    if try_connect_str == '':
                        try_connect_str = cur_repeat_str
                    try_connect_str = try_connect_str + ref_seq[i]
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
    return repeat_dict

def getReverseSequence(sequence):
    base_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    res = ''
    length = len(sequence)
    i = length - 1
    while i >= 0:
        base = sequence[i]
        if base not in base_map.keys():
            base = 'N'
        else:
            base = base_map[base]
        res += base
        i -= 1
    return res

def read_fasta(fasta_path):
    contignames = []
    contigs = {}
    if os.path.exists(fasta_path):
        with open(fasta_path, 'r') as rf:
            contigname = ''
            contigseq = ''
            for line in rf:
                if line.startswith('>'):
                    if contigname != '' and contigseq != '':
                        contigs[contigname] = contigseq
                        contignames.append(contigname)
                    contigname = line.strip()[1:].split(" ")[0].split('\t')[0]
                    contigseq = ''
                else:
                    contigseq += line.strip().upper()
            if contigname != '' and contigseq != '':
                contigs[contigname] = contigseq
                contignames.append(contigname)
        rf.close()
    return contignames, contigs

def read_fasta_v2(fasta_path):
    contignames = []
    contigs = {}
    if os.path.exists(fasta_path):
        with open(fasta_path, 'r') as rf:
            contigname = ''
            contigseq = ''
            for line in rf:
                if line.startswith('>'):
                    if contigname != '' and contigseq != '':
                        contigs[contigname] = contigseq
                        contignames.append(contigname)
                    contigname = line.strip()[1:].split(" ")[0].split('\t')[0]
                    contigseq = ''
                else:
                    contigseq += line.strip()
            if contigname != '' and contigseq != '':
                contigs[contigname] = contigseq
                contignames.append(contigname)
        rf.close()
    return contignames, contigs

def read_fasta_v1(fasta_path):
    contignames = []
    contigs = {}
    if os.path.exists(fasta_path):
        with open(fasta_path, 'r') as rf:
            contigname = ''
            contigseq = ''
            for line in rf:
                if line.startswith('>'):
                    if contigname != '' and contigseq != '':
                        contigs[contigname] = contigseq
                        contignames.append(contigname)
                    contigname = line.strip()[1:]
                    contigseq = ''
                else:
                    contigseq += line.strip().upper()
            if contigname != '' and contigseq != '':
                contigs[contigname] = contigseq
                contignames.append(contigname)
        rf.close()
    return contignames, contigs

def split_repeats(repeats_path, long_repeat_threshold, repeats_minimap2, repeats_bwa):
    r_contignames, r_contigs = read_fasta(repeats_path)
    long_seqs = {}
    short_seqs = {}
    for r_name in r_contignames:
        if len(r_contigs[r_name]) >= long_repeat_threshold:
            long_seqs[r_name] = r_contigs[r_name]
        else:
            short_seqs[r_name] = r_contigs[r_name]

    with open(repeats_minimap2, 'w') as f_save:
        for r_name in long_seqs.keys():
            f_save.write('>' + r_name + '\n' + long_seqs[r_name] + '\n')
    f_save.close()

    with open(repeats_bwa, 'w') as f_save:
        for r_name in short_seqs.keys():
            f_save.write('>' + r_name + '\n' + short_seqs[r_name] + '\n')
    f_save.close()



def compute_identity(cigar, NM_tag, method):
    n = int(NM_tag)
    if n == -1:
        return -1
    cigar = str(cigar)
    if method == "BLAST":
        l = 0
        it = re.finditer("(\d+)[MID]", cigar)
        for match in it:
            l += int(match.groups()[0])
        identity = float(l-n)/l
    elif method == "Gap-compressed":
        m = 0
        g = 0
        o = 0
        it = re.finditer("(\d+)M", cigar)
        for match in it:
            m += int(match.groups()[0])
        it = re.finditer("(\d+)[ID]", cigar)
        for match in it:
            g += int(match.groups()[0])
            o += 1
        identity = 1 - float(n-g+o) / (m+o)
    identity = format(identity, '.5f')
    return identity

def store2file(data_partition, cur_consensus_path):
    if len(data_partition) > 0:
        with open(cur_consensus_path, 'w') as f_save:
            for item in data_partition:
                f_save.write('>'+item[0]+'\n'+item[1]+'\n')
        f_save.close()

def PET(seq_item, partitions):
    # sort contigs by length
    original = seq_item
    original = sorted(original, key=lambda x: len(x[1]), reverse=True)
    return divided_array(original, partitions)

def divided_array(original_array, partitions):
    final_partitions = [[] for _ in range(partitions)]
    node_index = 0

    read_from_start = True
    read_from_end = False
    i = 0
    j = len(original_array) - 1
    while i <= j:
        # read from file start
        if read_from_start:
            final_partitions[node_index % partitions].append(original_array[i])
            i += 1
        if read_from_end:
            final_partitions[node_index % partitions].append(original_array[j])
            j -= 1
        node_index += 1
        if node_index % partitions == 0:
            # reverse
            read_from_end = bool(1 - read_from_end)
            read_from_start = bool(1 - read_from_start)
    return final_partitions


def multi_line(fasta_path, line_len):
    tmp_fasta_path = fasta_path + ".tmp"
    contigNames, contigs = read_fasta(fasta_path)
    with open(tmp_fasta_path, 'w') as f_w:
        for contigName in contigNames:
            contig = contigs[contigName]
            # line = '>' + contigName + '\t' + contig + '\n'
            # f_w.write(line)
            start = 0
            end = len(contig)
            while start < end:
                # add extra kmer length
                seg = contig[start:start+line_len]
                line = '>' + contigName + '\t' + str(start) + '\t' + seg + '\n'
                f_w.write(line)
                start += line_len
    f_w.close()
    return tmp_fasta_path

def split2cluster_normal(segments, partitions_num):
    avg_num = int(len(segments)/partitions_num)
    avg_num = len(segments) if avg_num == 0 else avg_num

    segments_cluster = {}
    cur_segment = []
    partition_index = 0
    last_index = -1
    for i in range(len(segments)):
        if i != 0 and i % avg_num == 0:
            segments_cluster[partition_index] = cur_segment
            cur_segment = []
            partition_index = partition_index + 1
            # last partition
            if partition_index == partitions_num-1:
                last_index = i
                break
        cur_segment.append(segments[i])
    # only one partition
    if len(cur_segment) > 0:
        segments_cluster[partition_index] = cur_segment
    else:
        if last_index != -1:
            for j in range(last_index, len(segments)):
                cur_segment.append(segments[j])
            segments_cluster[partition_index] = cur_segment
    return segments_cluster

def split2cluster(segments, partitions_num):
    avg_num = int(len(segments)/partitions_num)
    avg_num = len(segments) if avg_num == 0 else avg_num

    segments_cluster = {}
    cur_segment = {}
    partition_index = 0
    last_index = -1
    for i in range(len(segments)):
        if i != 0 and i % avg_num == 0:
            segments_cluster[partition_index] = cur_segment
            cur_segment = {}
            partition_index = partition_index + 1
            # last partition
            if partition_index == partitions_num-1:
                last_index = i
                break
        parts = segments[i].split('\t')
        ref_name = parts[0].replace('>', '')
        start = parts[1]
        seq = parts[2]
        new_ref_name = ref_name + '$' + start
        # seq = segments[i]
        # new_ref_name = 'ref$'+str(i)
        cur_segment[new_ref_name] = seq
    # only one partition
    if len(cur_segment) > 0:
        segments_cluster[partition_index] = cur_segment
    else:
        if last_index != -1:
            for j in range(last_index, len(segments)):
                parts = segments[j].split('\t')
                ref_name = parts[0].replace('>', '')
                start = parts[1]
                seq = parts[2]
                new_ref_name = ref_name + '$' + start
                cur_segment[new_ref_name] = seq
                # seq = segments[j]
                # new_ref_name = 'ref$' + str(j)
                # cur_segment[new_ref_name] = seq
            segments_cluster[partition_index] = cur_segment
    return segments_cluster

def filter_not_multi_mapping(cur_records, not_multi_mapping_repeatIds_dict, partiton_index, blast_records):
    log.logger.debug('partition %d process: %d records' % (partiton_index, len(cur_records)))
    multi_mapping_records = []
    for record in cur_records:
        query_name = record[0]
        # filter not multiple mapping repeat
        if not_multi_mapping_repeatIds_dict.__contains__(query_name):
            continue
        multi_mapping_records.append(record)
    blast_records[partiton_index] = multi_mapping_records

def get_alignment_info_v2(blastn_output):
    unmapped_repeatIds = []
    single_mapped_repeatIds = []
    multi_mapping_repeatIds = []

    query_records = {}
    with open(blastn_output, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            query_name = parts[0]
            target_name = parts[1]
            identity = float(parts[2])
            match_base = int(parts[3])
            query_length = int(query_name.split('-')[1].split('_')[1])

            if not query_records.__contains__(query_name):
                query_records[query_name] = []
            records = query_records[query_name]
            records.append((query_name, target_name, identity, match_base, query_length))
            query_records[query_name] = records
    f_r.close()

    for query_name in query_records.keys():
        complete_alignment_num = 0
        for record in query_records[query_name]:
            identity = record[2]
            match_base = record[3]
            query_len = record[4]
            # complete Match in cigar
            if float(match_base)/query_len >= 0.8 and identity >= 80:
                complete_alignment_num += 1

        if complete_alignment_num == 1:
            single_mapped_repeatIds.append(query_name)
        elif complete_alignment_num > 1:
            multi_mapping_repeatIds.append(query_name)
        else:
            unmapped_repeatIds.append(query_name)

    return unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds

def get_ltr_suppl_from_ltrfinder(merged_ltr, cluster_file, suppl_ltr_file):
    cluster_info = {}
    cluster_id = ''
    with open(cluster_file, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            if line.startswith('>Cluster'):
                cluster_id = line
                continue
            if cluster_id != '':
                if not cluster_info.__contains__(cluster_id):
                    cluster_info[cluster_id] = []
                cluster_records = cluster_info[cluster_id]
                cluster_records.append(line)
                cluster_info[cluster_id] = cluster_records
    f_r.close()

    keep_contigname = []
    for cluster_id in cluster_info.keys():
        ltr_retriever_count = 0
        contigname = ''
        for index, record in enumerate(cluster_info[cluster_id]):
            # representative record
            record = str(record)
            if record.endswith('... *') and record.__contains__('>Node_'):
                contigname = record.split('>')[1].replace('... *', '')
            if not record.__contains__('>Node_'):
                ltr_retriever_count += 1
        if ltr_retriever_count < 2 and contigname != '':
            keep_contigname.append(contigname)

    with open(suppl_ltr_file, 'w') as f_save:
        contignames, contigs = read_fasta(merged_ltr)
        for name in keep_contigname:
            for contigname in contignames:
                if contigname.__contains__(name):
                    f_save.write('>'+contigname+'\n'+contigs[contigname]+'\n')
    f_save.close()


def store_fasta(contigs, file_path):
    with open(file_path, 'w') as f_save:
        for name in contigs.keys():
            seq = contigs[name]
            f_save.write('>'+name+'\n'+seq+'\n')
    f_save.close()

def printClass(filepath, log):
    contignames, contigs = read_fasta(filepath)
    class_names = {}
    ltr_set = {}
    for name in contignames:
        class_name = name.split('#')[1]
        if class_name.__contains__('LTR'):
            ltr_set[name] = contigs[name]
        if not class_names.__contains__(class_name):
            class_names[class_name] = 0
        num = class_names[class_name]
        class_names[class_name] = num + 1
    log.logger.debug(class_names)
    return ltr_set

def parse_ref_blast_output(blastnResults_path, target_path, candidate_repeats_path):
    targetContigNames, targetContigs = read_fasta(target_path)
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
    f_r.close()
    #print(query_records)

    # strcuture: {'Node_1': {'Node_1': [(),(),], 'Node_2':[(),(),], ...}}
    # Node_0-len_5109 Node_0-len_5109 100.000 4651    0       0       459     5109    1       4651    0.0     8589
    # Node_0-len_5109 Node_30444-len_20481    100.000 217     0       0       1       217     20265   20481   1.37e-110       401

    # step2. splice sequence
    # a map is used to avoid adding redudant sequence
    candidate_family_repeat = []
    perfect_query = {}
    for query_name in query_records.keys():
        records = query_records[query_name]
        for target_name in records.keys():
            for record in records[target_name]:
                # identity < 80% should be neglected
                if record[0] < 80:
                    continue
                if record[0] >= 95:
                    if perfect_query.__contains__(query_name):
                        continue
                    else:
                        perfect_query[query_name] = 1
                t_start = record[5]
                t_end = record[6]
                if t_start > t_end:
                    t_tmp = t_start
                    t_start = t_end
                    t_end = t_tmp
                seg_seq = targetContigs[target_name][t_start: t_end]
                candidate_family_repeat.append(seg_seq)
                # if seg_seq.__contains__('N'):
                #     print((query_name, target_name, record))

    # step3. generate candidate repeats
    node_index = 0
    with open(candidate_repeats_path, 'w') as f_save:
        for sequence in candidate_family_repeat:
            f_save.write('>Node_'+str(node_index)+'-len_'+str(len(sequence))+'\n'+sequence+'\n')
            node_index += 1
    f_save.close()


def filter_LTR_high_similarity(blastnResults_path, target_path, query_path, filter_ltr_repeats_path):
    targetContigNames, targetContigs = read_fasta(target_path)
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
            same_target_records.append((identity, match_base, q_start, q_end, t_start, t_end))
            records[target_name] = same_target_records
            query_records[query_name] = records
    f_r.close()

    removed_names = set()
    for query_name in query_records.keys():
        records = query_records[query_name]
        for target_name in records.keys():
            target_seq = targetContigs[target_name]
            for record in records[target_name]:
                if record[0] >= 80 and float(record[1])/len(target_seq) >= 0.8:
                    removed_names.add(query_name)

    contignames, contigs = read_fasta(query_path)
    with open(filter_ltr_repeats_path, 'w') as f_save:
        for name in contignames:
            if name not in removed_names:
                f_save.write('>'+name+'\n'+contigs[name]+'\n')
    f_save.close()


def extract_tandem_from_trf(trf_data_path):
    tandem_elements = []
    with open(trf_data_path, 'r') as f_r:
        for line in f_r:
            parts = line.split(' ')
            if len(parts) == 15:
                tandem_elements.append(parts[13])
    f_r.close()
    return tandem_elements


def get_candidate_repeats(reference, k_num, reduce_partitions_num, unique_kmer_map, fault_tolerant_bases, tmp_output_dir, log):
    ref_name, ref_contigs = read_fasta(reference)
    segments = []
    for name in ref_name:
        parts = name.split('$')
        line = parts[0] + '\t' + parts[1] + '\t' + ref_contigs[name]
        segments.append(line)
    # with open(reference, 'r') as f_r:
    #     for line in f_r:
    #         if line.startswith('>'):
    #             continue
    #         line = line.replace('\n', '')
    #         segments.append(line)
    segments_cluster = split2cluster(segments, reduce_partitions_num)

    # partiton_index = 0
    # cur_segments = segments_cluster[partiton_index]
    # repeat_dict = generate_candidate_repeats_v2(cur_segments, k_num, unique_kmer_map, partiton_index, fault_tolerant_bases)

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

    # store connected_repeats for testing
    repeat_dict_file = tmp_output_dir + '/repeat_dict.csv'
    with codecs.open(repeat_dict_file, 'w', encoding='utf-8') as f:
        json.dump(repeat_dict, f)

    connected_repeats = {}
    for ref_name in repeat_dict.keys():
        repeat_list = repeat_dict[ref_name]
        repeat_list.sort(key=lambda x: (x[0], x[1]))

        if not connected_repeats.__contains__(ref_name):
            connected_repeats[ref_name] = []
        connected_repeat_list = connected_repeats[ref_name]

        # connect repeats
        last_start_pos = -1
        last_end_pos = -1
        last_repeat_str = ''
        for repeat_item in repeat_list:
            start_pos = repeat_item[0]
            end_pos = repeat_item[1]
            repeat_str = repeat_item[2]
            if last_start_pos != -1:
                if (start_pos - last_end_pos) == 1:
                    # connected repeat
                    last_end_pos = end_pos
                    last_repeat_str += repeat_str
                else:
                    # not connected repeat
                    # keep last connected repeat
                    connected_repeat_list.append((last_start_pos, last_end_pos, last_repeat_str))
                    last_start_pos = -1
                    last_end_pos = -1
                    last_repeat_str = ''
            if last_start_pos == -1:
                # start a new connected repeat
                last_start_pos = start_pos
                last_end_pos = end_pos
                last_repeat_str = repeat_str
        if last_start_pos != -1:
            connected_repeat_list.append((last_start_pos, last_end_pos, last_repeat_str))

    # store connected_repeats for testing
    connected_repeats_file = tmp_output_dir + '/connected_repeats.csv'
    with codecs.open(connected_repeats_file, 'w', encoding='utf-8') as f:
        json.dump(connected_repeats, f)

    return connected_repeats

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

def TSDsearch_v1(orig_seq, tir_start, tir_end):
    TIR_TSDs = [11, 10, 9, 8, 6, 5, 4, 3, 2]  # 4-> TTAA
    #TIR_seq = orig_seq[tir_start-1: tir_end]
    tsd_seq = ''
    for tsd_len in TIR_TSDs:
        if tir_start-1-tsd_len >= 0 and tir_end+tsd_len <= len(orig_seq):
            left_tsd = orig_seq[tir_start-1-tsd_len: tir_start-1]
            right_tsd = orig_seq[tir_end: tir_end+tsd_len]
            if left_tsd == right_tsd and tsd_len != 4:
                tsd_seq = left_tsd
                break
            elif tsd_len == 4 and left_tsd == 'TTAA':
                tsd_seq = left_tsd
                break
    return tsd_seq


def allow_mismatch(left_tsd, right_tsd, allow_mismatch_num):
    mismatch_num = 0
    for i in range(len(left_tsd)):
        if left_tsd[i] == right_tsd[i]:
            continue
        else:
            mismatch_num += 1
    if mismatch_num <= allow_mismatch_num:
        return True
    else:
        return False

def TSDsearch_v4(orig_seq, tir_start, tir_end):
    TIR_TSDs = [11, 10, 9, 8, 6, 5, 4, 3, 2]  # 4-> TTAA
    # TIR_seq = orig_seq[tir_start-1: tir_end]
    for tsd_len in TIR_TSDs:
        left_tsd_seq = ''
        right_tsd_seq = ''
        allow_mismatch_num = 1
        first_5bp = orig_seq[tir_start - 1: tir_start + 4]
        last_5bp = orig_seq[tir_end - 5: tir_end]
        if tir_start-1-tsd_len >= 0 and tir_end+tsd_len <= len(orig_seq):
            left_tsd = orig_seq[tir_start-1-tsd_len: tir_start-1]
            right_tsd = orig_seq[tir_end: tir_end+tsd_len]
            if left_tsd == right_tsd:
                if (tsd_len != 2 and tsd_len != 3 and tsd_len != 4) or (tsd_len == 4 and left_tsd == 'TTAA') \
                        or (tsd_len == 2 and (left_tsd == 'TA' or (first_5bp == 'CACTA' and last_5bp == 'TAGTG'))) \
                        or (tsd_len == 3 and (left_tsd == 'TAA' or left_tsd == 'TTA' or (first_5bp == 'CACTA' and last_5bp == 'TAGTG'))):
                    left_tsd_seq = left_tsd
                    right_tsd_seq = right_tsd
                    break
            elif tsd_len >= 8 and allow_mismatch(left_tsd, right_tsd, allow_mismatch_num):
                left_tsd_seq = left_tsd
                right_tsd_seq = right_tsd
                break

    return left_tsd_seq, right_tsd_seq

def TSDsearch_v3(orig_seq, tir_start, tir_end, tsd, plant):
    tsd_len = len(tsd)
    left_tsd_seq = ''
    right_tsd_seq = ''
    allow_mismatch_num = 1
    first_5bp = orig_seq[tir_start - 1: tir_start + 4]
    last_5bp = orig_seq[tir_end - 5: tir_end]
    first_3bp = orig_seq[tir_start - 1: tir_start + 2]
    last_3bp = orig_seq[tir_end - 3: tir_end]
    if tir_start-1-tsd_len >= 0 and tir_end+tsd_len <= len(orig_seq):
        left_tsd = orig_seq[tir_start-1-tsd_len: tir_start-1]
        right_tsd = orig_seq[tir_end: tir_end+tsd_len]
        if left_tsd == right_tsd:
            if (tsd_len != 2 and tsd_len != 3 and tsd_len != 4) or (tsd_len == 4 and left_tsd == 'TTAA') \
                    or (tsd_len == 2 and (left_tsd == 'TA' or (plant == 0 and first_3bp == 'CCC' and last_3bp == 'GGG'))) \
                    or (tsd_len == 3 and (left_tsd == 'TAA' or left_tsd == 'TTA'
                                          or (plant == 1 and ((first_5bp == 'CACTA' and last_5bp == 'TAGTG')
                                                              or (first_5bp == 'CACTG' and last_5bp == 'CAGTG'))))):
                left_tsd_seq = left_tsd
                right_tsd_seq = right_tsd

        elif tsd_len >= 8 and allow_mismatch(left_tsd, right_tsd, allow_mismatch_num):
            left_tsd_seq = left_tsd
            right_tsd_seq = right_tsd

    return left_tsd_seq, right_tsd_seq

def TSDsearch_v2_bak(orig_seq, tir_start, tir_end, TSD_set, plant):
    plant = int(plant)
    # 2->TA or 5'-CCC...GGG-3' in animal/fungi, 3-> 5'-CACT(A/G)...(C/T)AGTG-3' in plants or (TAA or TTA), 4-> TTAA,
    # 2->TA or 5'-CCC...GGG-3' in animal/fungi, 3-> 5'-CACT(A/G)...(C/T)AGTG-3' in plants or any 3bp (PIF-harbinger), 4-> TTAA,

    TIR_TSDs = [11, 10, 9, 8, 6, 5, 4, 3, 2]
    #TIR_seq = orig_seq[tir_start-1: tir_end]
    first_5bp = orig_seq[tir_start-1: tir_start+4]
    last_5bp = orig_seq[tir_end-5: tir_end]
    first_3bp = orig_seq[tir_start - 1: tir_start + 2]
    last_3bp = orig_seq[tir_end - 3: tir_end]
    tsd_seq = ''
    allow_mismatch_num = 1
    for tsd_len in TIR_TSDs:
        if tir_start - 1 - tsd_len >= 0 and tir_end + tsd_len <= len(orig_seq):
            left_tsd = orig_seq[tir_start - 1 - tsd_len: tir_start - 1]
            right_tsd = orig_seq[tir_end: tir_end + tsd_len]
            if left_tsd == right_tsd:
                if (tsd_len != 2 and tsd_len != 3 and tsd_len != 4) or (tsd_len == 4 and left_tsd == 'TTAA') \
                        or (tsd_len == 2 and (left_tsd == 'TA' or (plant == 0 and first_3bp == 'CCC' and last_3bp == 'GGG'))) \
                        or (tsd_len == 3 and (left_tsd == 'TAA' or left_tsd == 'TTA'
                                              or (plant == 1 and ((first_5bp == 'CACTA' and last_5bp == 'TAGTG')
                                                                  or (first_5bp == 'CACTG' and last_5bp == 'CAGTG'))))):
                    left_tsd_seq = left_tsd
                    right_tsd_seq = right_tsd
                    TSD_set.add((tir_start - tsd_len, tir_start - 1, left_tsd_seq, tir_end + 1, tir_end + tsd_len, right_tsd_seq,
                                tir_start, tir_end, tir_end - tir_start + 1))
            elif tsd_len >= 8 and allow_mismatch(left_tsd, right_tsd, allow_mismatch_num):
                left_tsd_seq = left_tsd
                right_tsd_seq = right_tsd
                TSD_set.add((tir_start - tsd_len, tir_start - 1, left_tsd_seq, tir_end + 1, tir_end + tsd_len, right_tsd_seq,
                     tir_start, tir_end, tir_end - tir_start + 1))

def TSDsearch_v2(orig_seq, tir_start, tir_end, TSD_set, plant):
    plant = int(plant)
    # 2->TA or in animal/fungi 5'-CCC...GGG-3', 3-> in plants 5'-CACT(A/G)...(C/T)AGTG-3' or any 3bp (PIF-harbinger), 4-> TTAA,
    TIR_TSDs = [11, 10, 9, 8, 6, 5, 4, 3, 2]
    #TIR_seq = orig_seq[tir_start-1: tir_end]
    first_5bp = orig_seq[tir_start-1: tir_start+4]
    last_5bp = orig_seq[tir_end-5: tir_end]
    first_3bp = orig_seq[tir_start - 1: tir_start + 2]
    last_3bp = orig_seq[tir_end - 3: tir_end]
    tsd_seq = ''
    allow_mismatch_num = 1
    for tsd_len in TIR_TSDs:
        if tir_start - 1 - tsd_len >= 0 and tir_end + tsd_len <= len(orig_seq):
            left_tsd = orig_seq[tir_start - 1 - tsd_len: tir_start - 1]
            right_tsd = orig_seq[tir_end: tir_end + tsd_len]
            if left_tsd == right_tsd:
                if (tsd_len != 2 and tsd_len != 4) or (tsd_len == 4 and left_tsd == 'TTAA') \
                        or (tsd_len == 2 and (left_tsd == 'TA' or (plant == 0 and first_3bp == 'CCC' and last_3bp == 'GGG'))):
                    left_tsd_seq = left_tsd
                    right_tsd_seq = right_tsd
                    TSD_set.add((tir_start - tsd_len, tir_start - 1, left_tsd_seq, tir_end + 1, tir_end + tsd_len, right_tsd_seq,
                                tir_start, tir_end, tir_end - tir_start + 1))
            elif tsd_len >= 8 and allow_mismatch(left_tsd, right_tsd, allow_mismatch_num):
                left_tsd_seq = left_tsd
                right_tsd_seq = right_tsd
                TSD_set.add((tir_start - tsd_len, tir_start - 1, left_tsd_seq, tir_end + 1, tir_end + tsd_len, right_tsd_seq,
                     tir_start, tir_end, tir_end - tir_start + 1))

def get_boundary_ungap_str(raw_align_seq, boundary_pos, search_len, direct):
    valid_count = 0
    col_index = boundary_pos
    ungap_str = ''
    if direct == 'right':
        while valid_count < search_len and col_index < len(raw_align_seq):
            cur_base = raw_align_seq[col_index]
            if cur_base == '-':
                col_index += 1
                continue
            ungap_str += cur_base
            col_index += 1
            valid_count += 1
    else:
        while valid_count < search_len and col_index >= 0:
            cur_base = raw_align_seq[col_index]
            if cur_base == '-':
                col_index -= 1
                continue
            ungap_str = cur_base + ungap_str
            col_index -= 1
            valid_count += 1
    return ungap_str

def TSDsearch_v5(raw_align_seq, cur_boundary_start, cur_boundary_end, plant):
    plant = int(plant)
    # 2->TA or animal/fungi 5'-CCC...GGG-3', 3-> plant 5'-CACT(A/G)...(C/T)AGTG-3' or （TAA或TTA）, 4-> TTAA,
    TIR_TSDs = [11, 10, 9, 8, 6, 5, 4, 3, 2]

    direct = 'right'
    first_5bp = get_boundary_ungap_str(raw_align_seq, cur_boundary_start, 5, direct)
    first_3bp = get_boundary_ungap_str(raw_align_seq, cur_boundary_start, 3, direct)
    direct = 'left'
    last_5bp = get_boundary_ungap_str(raw_align_seq, cur_boundary_end, 5, direct)
    last_3bp = get_boundary_ungap_str(raw_align_seq, cur_boundary_end, 3, direct)


    left_tsd_seq = ''
    right_tsd_seq = ''
    allow_mismatch_num = 1
    for tsd_len in TIR_TSDs:
        left_tsd = get_boundary_ungap_str(raw_align_seq, cur_boundary_start-1, tsd_len, 'left')
        right_tsd = get_boundary_ungap_str(raw_align_seq, cur_boundary_end+1, tsd_len, 'right')
        if len(left_tsd) != len(right_tsd) or len(left_tsd) != tsd_len:
            continue
        if left_tsd == right_tsd:
            if (tsd_len != 2 and tsd_len != 3 and tsd_len != 4) or (tsd_len == 4 and left_tsd == 'TTAA') \
                    or (tsd_len == 2 and (left_tsd == 'TA' or (plant == 0 and first_3bp == 'CCC' and last_3bp == 'GGG'))) \
                    or (tsd_len == 3 and (left_tsd == 'TAA' or left_tsd == 'TTA'
                                          or (plant == 1 and ((first_5bp == 'CACTA' and last_5bp == 'TAGTG')
                                                              or (first_5bp == 'CACTG' and last_5bp == 'CAGTG'))))):
                left_tsd_seq = left_tsd
                right_tsd_seq = right_tsd
        elif tsd_len >= 8 and allow_mismatch(left_tsd, right_tsd, allow_mismatch_num):
            left_tsd_seq = left_tsd
            right_tsd_seq = right_tsd
    return left_tsd_seq, right_tsd_seq



def TSDsearch_ltr(orig_seq, ltr_start, ltr_end, TSD_set):
    LTR_TSDs = [6, 5, 4]  # LTR:most of them is 5-'TG...CA-3'
    #LTR_seq = orig_seq[ltr_start-1: ltr_end]
    first_2bp = orig_seq[ltr_start-1: ltr_start+1]
    last_2bp = orig_seq[ltr_end-2: ltr_end]
    allow_mismatch_num = 1
    for tsd_len in LTR_TSDs:
        if ltr_start - 1 - tsd_len >= 0 and ltr_end + tsd_len <= len(orig_seq):
            left_tsd = orig_seq[ltr_start - 1 - tsd_len: ltr_start - 1]
            right_tsd = orig_seq[ltr_end: ltr_end + tsd_len]
            if left_tsd == right_tsd:
                left_tsd_seq = left_tsd
                right_tsd_seq = right_tsd
                TSD_set.add((ltr_start - tsd_len, ltr_start - 1, left_tsd_seq, ltr_end + 1, ltr_end + tsd_len, right_tsd_seq,
                            ltr_start, ltr_end, ltr_end - ltr_start + 1))
            elif first_2bp == 'TG' and last_2bp == 'CA' and allow_mismatch(left_tsd, right_tsd, allow_mismatch_num):
                left_tsd_seq = left_tsd
                right_tsd_seq = right_tsd
                TSD_set.add((ltr_start - tsd_len, ltr_start - 1, left_tsd_seq, ltr_end + 1, ltr_end + tsd_len, right_tsd_seq,
                     ltr_start, ltr_end, ltr_end - ltr_start + 1))

def TSDsearch_ltr_v1(orig_seq, ltr_start, ltr_end, tsd_len):
    left_tsd_seq = ''
    right_tsd_seq = ''
    allow_mismatch_num = 1
    first_2bp = orig_seq[ltr_start - 1: ltr_start + 1]
    last_2bp = orig_seq[ltr_end - 2: ltr_end]
    if ltr_start-1-tsd_len >= 0 and ltr_end+tsd_len <= len(orig_seq):
        left_tsd = orig_seq[ltr_start-1-tsd_len: ltr_start-1]
        right_tsd = orig_seq[ltr_end: ltr_end+tsd_len]
        if left_tsd == right_tsd:
            left_tsd_seq = left_tsd
            right_tsd_seq = right_tsd
        elif first_2bp == 'TG' and last_2bp == 'CA' and allow_mismatch(left_tsd, right_tsd, allow_mismatch_num):
            left_tsd_seq = left_tsd
            right_tsd_seq = right_tsd

    return left_tsd_seq, right_tsd_seq

def calculate_max_min(x, max, min):
    distance = max-min
    return round((x-min)/distance, 4)


def get_score(confident_TIR):
    # (tir_len, tsd, te_len, contigs[name])
    TE_len_list = []
    tir_len_list = []
    tsd_len_list = []
    for info in confident_TIR:
        TE_len_list.append(info[2])
        tir_len_list.append(info[0])
        tsd_len_list.append(len(info[1]))

    max_TE_len = max(TE_len_list)
    min_TE_len = min(TE_len_list)

    max_tir_len = max(tir_len_list)
    min_tir_len = min(tir_len_list)

    max_tsd_len = max(tsd_len_list)
    min_tsd_len = min(tsd_len_list)

    max_info = None
    max_score = -1
    for info in confident_TIR:
        if max_TE_len == min_TE_len:
            TE_len_normal = 0
        else:
            TE_len_normal = calculate_max_min(info[2], max_TE_len, min_TE_len)
        if max_tir_len == min_tir_len:
            tir_len_normal = 0
        else:
            tir_len_normal = calculate_max_min(info[0], max_tir_len, min_tir_len)
        if max_tsd_len == min_tsd_len:
            tsd_len_normal = 0
        else:
            tsd_len_normal = calculate_max_min(len(info[1]), max_tsd_len, min_tsd_len)
        score = 0.3*TE_len_normal+0.3*tir_len_normal+0.4*tsd_len_normal
        #score = 0.4 * TE_len_normal + 0.6 * tsd_len_normal
        if score > max_score:
            max_score = score
            max_info = info
    return max_info

def get_score_v1(confident_TIR):
    # (copy_num, tsd_len)
    copy_num_list = []
    tsd_len_list = []
    for info in confident_TIR:
        copy_num_list.append(info[0])
        tsd_len_list.append(info[1])

    max_copy_num = max(copy_num_list)
    min_copy_num = min(copy_num_list)

    max_tsd_len = max(tsd_len_list)
    min_tsd_len = min(tsd_len_list)

    max_info = None
    max_score = -1
    for info in confident_TIR:
        if max_copy_num == min_copy_num:
            TE_copy_num = 0
        else:
            TE_copy_num = calculate_max_min(info[0], max_copy_num, min_copy_num)
        if max_tsd_len == min_tsd_len:
            tsd_len_normal = 0
        else:
            tsd_len_normal = calculate_max_min(info[1], max_tsd_len, min_tsd_len)
        score = 0.6*TE_copy_num+0.4*tsd_len_normal
        if score > max_score:
            max_score = score
            max_info = info
    return max_info

def get_score_v2(confident_TIR):
    # (copy_num, tsd_len)
    tir_len_list = []
    tsd_len_list = []
    for info in confident_TIR:
        tir_len_list.append(info[0])
        tsd_len_list.append(info[1])

    max_tir_len = max(tir_len_list)
    min_tir_len = min(tir_len_list)

    max_tsd_len = max(tsd_len_list)
    min_tsd_len = min(tsd_len_list)

    max_info = None
    max_score = -1
    for info in confident_TIR:
        if max_tir_len == min_tir_len:
            TE_tir_len = 0
        else:
            TE_tir_len = calculate_max_min(info[0], max_tir_len, min_tir_len)
        if max_tsd_len == min_tsd_len:
            tsd_len_normal = 0
        else:
            tsd_len_normal = calculate_max_min(info[1], max_tsd_len, min_tsd_len)
        score = 0.3*TE_tir_len+0.7*tsd_len_normal
        if score > max_score:
            max_score = score
            max_info = info
    return max_info

def get_score_ltr(confident_LTR):
    # (tir_len, tsd, te_len, contigs[name])
    TE_len_list = []
    ltr_len_list = []
    tsd_len_list = []
    for info in confident_LTR:
        TE_len_list.append(info[2])
        ltr_len_list.append(info[0])
        tsd_len_list.append(len(info[1]))

    max_TE_len = max(TE_len_list)
    min_TE_len = min(TE_len_list)

    max_ltr_len = max(ltr_len_list)
    min_ltr_len = min(ltr_len_list)

    max_tsd_len = max(tsd_len_list)
    min_tsd_len = min(tsd_len_list)

    max_info = None
    max_score = -1
    for info in confident_LTR:
        if max_TE_len == min_TE_len:
            TE_len_normal = 0
        else:
            TE_len_normal = calculate_max_min(info[2], max_TE_len, min_TE_len)
        if max_ltr_len == min_ltr_len:
            ltr_len_normal = 0
        else:
            ltr_len_normal = calculate_max_min(info[0], max_ltr_len, min_ltr_len)
        if max_tsd_len == min_tsd_len:
            tsd_len_normal = 0
        else:
            tsd_len_normal = calculate_max_min(len(info[1]), max_tsd_len, min_tsd_len)
        score = 0.4*TE_len_normal+0.4*ltr_len_normal+0.2*tsd_len_normal
        if score > max_score:
            max_score = score
            max_info = info
    return max_info

def filter_dup_ltr(ltr_out, filter_dup_path):
    # Consideration is given to multiple factors such as TSD length, TIR length, and TE length.
    contignames, contigs = read_fasta_v1(ltr_out)
    filtered_contigs = {}
    for name in contignames:
        ltr_pos_infos = name.split(' ')[1].replace('LTR', '').replace('(', '').replace(')', '').split('..')
        left_ltr_pos = ltr_pos_infos[0].split(',')
        right_ltr_pos = ltr_pos_infos[1].split(',')

        left_ltr_start = left_ltr_pos[0]
        left_ltr_end = left_ltr_pos[1]
        right_ltr_start = right_ltr_pos[0]
        right_ltr_end = right_ltr_pos[1]

        ltr_len = int(name.split('Length ltr=')[1])
        parts = name.split('-C_')
        orig_query_name = parts[0]
        tsd = parts[1].split(' ')[0].split('-tsd_')[1]
        te_len = len(contigs[name])
        # te_len = int(parts[1].split('-')[1].split('_')[1])
        if not filtered_contigs.__contains__(orig_query_name):
            filtered_contigs[orig_query_name] = set()
        confident_TIR = filtered_contigs[orig_query_name]
        confident_TIR.add((ltr_len, tsd, te_len, contigs[name]))

    node_index = 0
    with open(filter_dup_path, 'w') as f_save:
        for name in filtered_contigs.keys():
            confident_LTR = filtered_contigs[name]
            highest_confident_LTR = get_score_ltr(confident_LTR)
            query_name = 'N_' + str(node_index) + '-ltrlen_' + str(highest_confident_LTR[0]) \
                         + '-TElen_' + str(highest_confident_LTR[2]) + '-lltr_' + str(left_ltr_start) + '..' + str(left_ltr_end)\
                         + '-rltr_' + str(right_ltr_start) + '..' + str(right_ltr_end) + '-tsd_' + str(highest_confident_LTR[1])
            f_save.write('>' + query_name + '\n' + highest_confident_LTR[3] + '\n')
            node_index += 1
    f_save.close()

def filter_dup_itr(tir_out, filter_dup_path):
    # Consideration is given to multiple factors such as TSD length, TIR length, and TE length.
    contignames, contigs = read_fasta_v1(tir_out)
    filtered_contigs = {}
    for name in contignames:
        tir_len = int(name.split('Length itr=')[1])
        parts = name.split('-C_')
        orig_query_name = parts[0]
        tsd = parts[1].split(' ')[0].split('-tsd_')[1]
        te_len = len(contigs[name])
        # te_len = int(parts[1].split('-')[1].split('_')[1])
        if not filtered_contigs.__contains__(orig_query_name):
            filtered_contigs[orig_query_name] = set()
        confident_TIR = filtered_contigs[orig_query_name]
        confident_TIR.add((tir_len, tsd, te_len, contigs[name]))

    node_index = 0
    with open(filter_dup_path, 'w') as f_save:
        for name in filtered_contigs.keys():
            confident_TIR = filtered_contigs[name]
            highest_confident_TIR = get_score(confident_TIR)
            query_name = 'N_' + str(node_index) + '-mintirlen_' + str(highest_confident_TIR[0]) + '-TElen_' + str(
                highest_confident_TIR[2]) + '-tsd_' + str(highest_confident_TIR[1])
            f_save.write('>' + query_name + '\n' + highest_confident_TIR[3] + '\n')
            node_index += 1
    f_save.close()

def filter_dup_itr_v1(cur_copies_out_contigs, seq_copynum):
    filtered_contigs = {}
    for name in cur_copies_out_contigs.keys():
        if seq_copynum.__contains__(name):
            copy_num = seq_copynum[name]
        else:
            copy_num = 0
        parts = name.split('-C_')
        orig_query_name = parts[0]
        tsd = parts[1].split('-tsd_')[1].split('-')[0]
        if not filtered_contigs.__contains__(orig_query_name):
            filtered_contigs[orig_query_name] = set()
        confident_TIR = filtered_contigs[orig_query_name]
        confident_TIR.add((copy_num, len(tsd), cur_copies_out_contigs[name]))

    res_contigs = {}
    for name in filtered_contigs.keys():
        confident_TIR = filtered_contigs[name]
        highest_confident_TIR = get_score_v1(confident_TIR)
        res_contigs[name] = highest_confident_TIR[2]
    return res_contigs

def filter_dup_itr_v2(cur_copies_out_contigs, TIR_len_dict):
    res_contigs = {}
    filtered_contigs = {}
    for name in cur_copies_out_contigs.keys():
        if not TIR_len_dict.__contains__(name):
            continue
        parts = name.split('-C_')
        orig_query_name = parts[0]
        tsd = parts[1].split('-tsd_')[1].split('-')[0]

        tir_len = TIR_len_dict[name]
        if not filtered_contigs.__contains__(orig_query_name):
            filtered_contigs[orig_query_name] = set()
        confident_TIR = filtered_contigs[orig_query_name]
        confident_TIR.add((tir_len, len(tsd), cur_copies_out_contigs[name]))
        
    for name in filtered_contigs.keys():
        confident_TIR = filtered_contigs[name]
        highest_confident_TIR = get_score_v2(confident_TIR)
        new_name = name + '-tir_' + str(highest_confident_TIR[0]) + '-tsd_' + str(highest_confident_TIR[1])
        res_contigs[new_name] = highest_confident_TIR[2]
    return res_contigs


def filter_dup_itr_v3(cur_copies_out_contigs, TIR_len_dict):
    res_contigs = {}
    min_distance = 100000
    min_distance_name = ''
    for name in cur_copies_out_contigs.keys():
        distance = int(name.split('-distance_')[1])
        if distance < min_distance:
            min_distance = distance
            min_distance_name = name

    parts = min_distance_name.split('-C_')
    orig_query_name = parts[0]
    tsd = parts[1].split('-tsd_')[1].split('-')[0]
    if not TIR_len_dict.__contains__(min_distance_name):
        tir_len = 0
    else:
        tir_len = TIR_len_dict[min_distance_name]
    new_name = orig_query_name + '-tir_' + str(tir_len) + '-tsd_' + str(tsd)
    # full-length TIR should shorter than 10K bp
    if len(cur_copies_out_contigs[min_distance_name]) < 10000:
        res_contigs[new_name] = cur_copies_out_contigs[min_distance_name]
    return res_contigs

def filter_dup_itr_v4(cur_copies_out_contigs):
    res_contigs = {}
    min_distance = 100000
    min_distance_name = ''
    for name in cur_copies_out_contigs.keys():
        distance = int(name.split('-distance_')[1])
        if distance < min_distance:
            min_distance = distance
            min_distance_name = name

    parts = min_distance_name.split('-C_')
    orig_query_name = parts[0]
    tsd = parts[1].split('-tsd_')[1].split('-')[0]
    new_name = orig_query_name + '-tsd_' + str(tsd)
    res_contigs[new_name] = cur_copies_out_contigs[min_distance_name]
    return res_contigs

def file_exist(resut_file):
    if os.path.isfile(resut_file):  # 输入是文件
        if os.path.getsize(resut_file) > 0:
            # 如果是FASTA文件类型
            if resut_file.endswith('.fa') or resut_file.endswith('.fasta'):
                names, contigs = read_fasta(resut_file)
                return len(contigs) > 0
            else:
                # 对非FASTA文件，检查是否包含非注释的有效内容
                with open(resut_file, 'r') as f_r:
                    for line in f_r:
                        if not line.startswith('#') and line.strip():
                            return True
                return False
        else:
            return False
    elif os.path.isdir(resut_file):  # 输入是目录
        # 检查目录是否非空
        return len(os.listdir(resut_file)) > 0
    else:
        # 输入既不是文件也不是目录
        return False


def run_remove_TR(target_file, trf_dir):
    if not os.path.exists(trf_dir):
        os.makedirs(trf_dir)
    # 1. Use TRF to search for tandem repeats in the target file.
    trf_command = 'cd ' + trf_dir + ' && trf ' + target_file + ' 2 7 7 80 10 50 500 -f -d -m'
    os.system(trf_command + ' > /dev/null 2>&1')

    (repeat_dir, repeat_filename) = os.path.split(target_file)
    trf_masked_repeats = trf_dir + '/' + repeat_filename + '.2.7.7.80.10.50.500.mask'
    # 2. Remove sequences that consist entirely of 'N's.
    trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    filter_contigs = {}
    for name in trf_contigNames:
        seq = trf_contigs[name]
        if all(c == "N" for c in seq):
            continue
        filter_contigs[name] = seq
    store_fasta(filter_contigs, target_file)

    return target_file

def run_TRF(input, input_dir, tandem_region_cutoff, TE_type):
    trf_dir = input_dir + '/trf'
    if not os.path.exists(trf_dir):
        os.makedirs(trf_dir)

    trf_command = 'cd ' + trf_dir + ' && trf ' + input + ' 2 7 7 80 10 50 500 -f -d -m'
    os.system(trf_command + ' > /dev/null 2>&1')

    (repeat_dir, repeat_filename) = os.path.split(input)
    trf_masked_repeats = trf_dir + '/' + repeat_filename + '.2.7.7.80.10.50.500.mask'

    trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    repeats_contigNames, repeats_contigs = read_fasta(input)
    repeats_path = input_dir + '/filter_tandem.fa'
    with open(repeats_path, 'w') as f_save:
        for name in trf_contigNames:
            seq = trf_contigs[name]
            if TE_type == 'ltr':
                ltr_type = name.split('-')[1]
                if ltr_type == 'LTR':
                    # Extract the 5' end of the sequence for 100bp, then determine if it is a tandem repeat.
                    start_seq = seq[0:100]
                    if float(start_seq.count('N')) / len(start_seq) >= tandem_region_cutoff:
                        continue
            elif TE_type == 'tir':
                # Extract 20bp from both the 5' and 3' ends of the sequence, then determine if it is a tandem repeat.
                start_seq = seq[0:20]
                end_seq = seq[-20:]
                if float(start_seq.count('N')) / len(start_seq) >= tandem_region_cutoff \
                        or float(end_seq.count('N')) / len(end_seq) >= tandem_region_cutoff:
                    continue
            if float(seq.count('N')) / len(seq) < tandem_region_cutoff and len(seq) >= 100:
                f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')
    f_save.close()

    return repeats_path

def multi_process_TRF(input, output, temp_dir, tandem_region_cutoff, threads = 48, TE_type = ''):
    contigNames, contigs = read_fasta(input)

    os.system('rm -rf '+temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    longest_repeat_files = []
    segments_cluster = divided_array(list(contigs.items()), threads)
    for partition_index, cur_segments in enumerate(segments_cluster):
        single_tmp_dir = temp_dir + '/' + str(partition_index)
        if not os.path.exists(single_tmp_dir):
            os.makedirs(single_tmp_dir)
        split_repeat_file = single_tmp_dir + '/repeats_split.fa'
        cur_contigs = {}
        for item in cur_segments:
            cur_contigs[item[0]] = item[1]
        if len(cur_contigs) > 0:
            store_fasta(cur_contigs, split_repeat_file)
            longest_repeat_files.append((split_repeat_file, single_tmp_dir))

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeat_files:
        input = file[0]
        input_dir = file[1]
        job = ex.submit(run_TRF, input, input_dir, tandem_region_cutoff, TE_type)
        jobs.append(job)
    ex.shutdown(wait=True)

    if os.path.exists(output):
        os.remove(output)
    for job in as_completed(jobs):
        cur_repeats_path = job.result()
        os.system('cat ' + cur_repeats_path + ' >> ' + output)

def multiple_alignment(repeats_path, blast_program_dir, tools_dir):
    split_repeats_path = repeats_path[0]
    original_repeats_path = repeats_path[1]
    blastn2Results_path = repeats_path[2]

    makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + original_repeats_path + ' > /dev/null 2>&1'
    align_command = blast_program_dir + '/bin/blastn -db ' + original_repeats_path + ' -num_threads ' \
                    + str(1) + ' -query ' + split_repeats_path + ' -outfmt 6 > ' + blastn2Results_path
    os.system(makedb_command)
    os.system(align_command)

    return blastn2Results_path

def get_longest_repeats_v1(repeats_path, blast_program_dir, fixed_extend_base_threshold, max_single_repeat_len, threads):
    split_repeats_path = repeats_path[0]
    original_repeats_path = repeats_path[1]
    blastn2Results_path = repeats_path[2]
    tmp_blast_dir = repeats_path[3]

    subject_tmp_dir = tmp_blast_dir + '/subject'
    for partition_index in range(threads):
        split_subject_file = subject_tmp_dir + '/' + str(partition_index) + '.fa'
        if not os.path.exists(split_subject_file):
            continue
        align_command = blast_program_dir + '/bin/blastn -db ' + split_subject_file + ' -num_threads ' \
                        + str(1) + ' -query ' + split_repeats_path
        if partition_index == 0:
            align_command1 = align_command + ' -outfmt 6 > ' + blastn2Results_path
            os.system(align_command1)
        else:
            align_command2 = align_command + ' -outfmt 6 >> ' + blastn2Results_path
            os.system(align_command2)

    query_names, query_contigs = read_fasta(split_repeats_path)

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
    f_r.close()

    keep_longest_query = {}
    longest_repeats = {}
    for idx, query_name in enumerate(query_records.keys()):
        query_len = len(query_contigs[query_name])
        #print('total query size: %d, current query name: %s, idx: %d' % (len(query_records), query_name, idx))

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
                        if (frag[2] - exist_frag[3] < fixed_extend_base_threshold):
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
                        if (exist_frag[3] - frag[2] < fixed_extend_base_threshold):
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

                cluster_longest_subject_start = -1
                cluster_longest_subject_end = -1
                cluster_longest_subject_len = -1

                cluster_extend_num = 0

                # print('subject pos size: %d' %(len(cur_cluster)))
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

                    cur_extend_num = 0

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
                                    if ext_frag[0] - longest_query_end < fixed_extend_base_threshold and ext_frag[
                                        2] - longest_subject_end < fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= fixed_extend_base_threshold:
                                        break
                            elif longest_subject_start > longest_subject_end and ext_frag[2] > ext_frag[3]:
                                # reverse
                                if ext_frag[3] < longest_subject_end:
                                    # reverse extend
                                    if ext_frag[
                                        0] - longest_query_end < fixed_extend_base_threshold and longest_subject_end - \
                                            ext_frag[2] < fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start > \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= fixed_extend_base_threshold:
                                        break
                    if cur_longest_query_len > cluster_longest_query_len:
                        cluster_longest_query_start = longest_query_start
                        cluster_longest_query_end = longest_query_end
                        cluster_longest_query_len = cur_longest_query_len

                        cluster_longest_subject_start = longest_subject_start
                        cluster_longest_subject_end = longest_subject_end
                        cluster_longest_subject_len = longest_subject_end - longest_subject_start

                        cluster_extend_num = cur_extend_num
                # keep this longest query
                if cluster_longest_query_len != -1:
                    longest_queries.append((cluster_longest_query_start, cluster_longest_query_end,
                                            cluster_longest_query_len, cluster_longest_subject_start,
                                            cluster_longest_subject_end, cluster_longest_subject_len, subject_name, cluster_extend_num))


        # we now consider, we should take some sequences from longest_queries to represent this query sequence.
        # we take the longest sequence by length, if the latter sequence overlap with the former sequence largely (50%),
        # continue find next sequence until the ratio of query sequence over 90% or no more sequences.
        longest_queries.sort(key=lambda x: -x[2])

        keep_longest_query[query_name] = longest_queries
        #parts = query_name.split('-s_')[1].split('-')
        #chr_name = parts[0]
        #chr_start = int(parts[2])
        #chr_end = int(parts[3])
        chr_name = ''
        chr_start = 0
        chr_end = 0

        is_flank = False
        query_seq = query_contigs[query_name]

        # Calculate the number of supported copies of all fragments, and keep only the longest boundary.
        copies = []
        for query in longest_queries:
            is_copy = False
            for i in range(len(copies)-1, -1, -1):
                copy = copies[i]
                if copy[0] <= query[1] and copy[0] >= query[0]:
                    overlap = query[1] - copy[0]
                elif query[0] >= copy[0] and query[1] <= copy[1]:
                    overlap = query[1] - query[0]
                elif copy[1] <= query[1] and copy[1] >= query[0]:
                    overlap = copy[1] - query[0]
                else:
                    overlap = 0
                if overlap > 0:
                    if float(overlap) / (query[1]-query[0]) >= 0.95:
                        if float(overlap) / (copy[1]-copy[0]) >= 0.95:
                            copies[i] = (min(copy[0], query[0]), max(copy[1], query[1]), copy[2] + 1)
                            is_copy = True
                            break
                        else:
                            is_copy = False
                            break
                    else:
                        is_copy = False
                        break
            if not is_copy and abs(query[1]-query[0]) <= max_single_repeat_len:
                copies.append((query[0], query[1], 1))

        copies.sort(key=lambda x: -(x[1]-x[0]))

        # Merge sequence fragments of TE.
        pure_copies = []
        for cur_copy in copies:
            # Remove copies where copy[2] < 2 and it is not a full-length duplicate, without other copies supporting it as a repeat.
            if cur_copy[2] < 2 and float(cur_copy[1] - cur_copy[0]) / query_len < 0.95:
                continue
            is_seg = False
            for i, pure_copy in enumerate(pure_copies):
                # Since the true TE ends are not repeats, try to take the longest boundary.
                if cur_copy[0] >= pure_copy[0] and abs(cur_copy[1]-pure_copy[1]) < 50:
                    pure_copies[i] = (pure_copy[0], max(cur_copy[1], pure_copy[1]), pure_copy[2])
                    is_seg = True
                    break
                elif cur_copy[1] <= pure_copy[1] and abs(cur_copy[0]-pure_copy[0]) < 50:
                    pure_copies[i] = (min(cur_copy[0], pure_copy[0]), pure_copy[1], pure_copy[2])
                    is_seg = True
                    break
            if not is_seg:
                pure_copies.append(cur_copy)

        intact_copies = []
        for copy in pure_copies:
            start_pos = copy[0] - 1
            end_pos = copy[1]
            ori_start_pos = start_pos
            ori_end_pos = end_pos

            copy_seq = query_seq[start_pos: end_pos]

            seq_ref_start = chr_start + ori_start_pos
            seq_ref_end = chr_start + ori_end_pos

            if len(copy_seq) >= 80:
                intact_copies.append((ori_start_pos, ori_end_pos, chr_name, seq_ref_start, seq_ref_end, copy_seq))
        longest_repeats[query_name] = intact_copies
    return longest_repeats, keep_longest_query

def pairwise_alignment(repeats_path, blast_program_dir):
    split_repeats_path = repeats_path[0]
    subject_path = repeats_path[1]
    blastn2Results_path = repeats_path[2]

    if not os.path.exists(subject_path):
        return
    align_command = blast_program_dir + '/bin/blastn -db ' + subject_path + ' -num_threads ' \
                    + str(1) + ' -query ' + split_repeats_path
    align_command2 = align_command + ' -outfmt 6 > ' + blastn2Results_path
    os.system(align_command2)
    return blastn2Results_path

def get_longest_repeats_v2(repeat_file, merged_output, fixed_extend_base_threshold, max_single_repeat_len):

    query_names, query_contigs = read_fasta(repeat_file)

    # parse blastn output, determine the repeat boundary
    # query_records = {query_name: {subject_name: [(q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end)] }}
    query_records = {}
    with open(merged_output, 'r') as f_r:
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
    f_r.close()

    keep_longest_query = {}
    longest_repeats = {}
    for idx, query_name in enumerate(query_records.keys()):
        query_len = len(query_contigs[query_name])
        #print('total query size: %d, current query name: %s, idx: %d' % (len(query_records), query_name, idx))

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
                        if (frag[2] - exist_frag[3] < fixed_extend_base_threshold):
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
                        if (exist_frag[3] - frag[2] < fixed_extend_base_threshold):
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

                cluster_longest_subject_start = -1
                cluster_longest_subject_end = -1
                cluster_longest_subject_len = -1

                cluster_extend_num = 0

                # print('subject pos size: %d' %(len(cur_cluster)))
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

                    cur_extend_num = 0

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
                                    if ext_frag[0] - longest_query_end < fixed_extend_base_threshold and ext_frag[
                                        2] - longest_subject_end < fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= fixed_extend_base_threshold:
                                        break
                            elif longest_subject_start > longest_subject_end and ext_frag[2] > ext_frag[3]:
                                # reverse
                                if ext_frag[3] < longest_subject_end:
                                    # reverse extend
                                    if ext_frag[
                                        0] - longest_query_end < fixed_extend_base_threshold and longest_subject_end - \
                                            ext_frag[2] < fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start > \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= fixed_extend_base_threshold:
                                        break
                    if cur_longest_query_len > cluster_longest_query_len:
                        cluster_longest_query_start = longest_query_start
                        cluster_longest_query_end = longest_query_end
                        cluster_longest_query_len = cur_longest_query_len

                        cluster_longest_subject_start = longest_subject_start
                        cluster_longest_subject_end = longest_subject_end
                        cluster_longest_subject_len = longest_subject_end - longest_subject_start

                        cluster_extend_num = cur_extend_num
                # keep this longest query
                if cluster_longest_query_len != -1:
                    longest_queries.append((cluster_longest_query_start, cluster_longest_query_end,
                                            cluster_longest_query_len, cluster_longest_subject_start,
                                            cluster_longest_subject_end, cluster_longest_subject_len, subject_name, cluster_extend_num))


        # we now consider, we should take some sequences from longest_queries to represent this query sequence.
        # we take the longest sequence by length, if the latter sequence overlap with the former sequence largely (50%),
        # continue find next sequence until the ratio of query sequence over 90% or no more sequences.
        longest_queries.sort(key=lambda x: -x[2])

        keep_longest_query[query_name] = longest_queries
        parts = query_name.split('$')
        chr_name = parts[0]
        chr_start = int(parts[1])
        #chr_end = int(parts[3])
        #ref_seq = ref_contigs[chr_name]

        is_flank = False
        flanking_len = 0

        query_seq = query_contigs[query_name]

        copies = []
        for query in longest_queries:
            is_copy = False
            for i in range(len(copies)-1, -1, -1):
                copy = copies[i]
                if copy[0] <= query[1] and copy[0] >= query[0]:
                    overlap = query[1] - copy[0]
                elif query[0] >= copy[0] and query[1] <= copy[1]:
                    overlap = query[1] - query[0]
                elif copy[1] <= query[1] and copy[1] >= query[0]:
                    overlap = copy[1] - query[0]
                else:
                    overlap = 0
                if overlap > 0:
                    if float(overlap) / (query[1]-query[0]) >= 0.95:
                        if float(overlap) / (copy[1]-copy[0]) >= 0.95:
                            copies[i] = (min(copy[0], query[0]), max(copy[1], query[1]), copy[2] + 1)
                            is_copy = True
                            break
                        else:
                            is_copy = False
                            break
                    else:
                        is_copy = False
                        break
            if not is_copy and abs(query[1]-query[0]) <= max_single_repeat_len:
                copies.append((query[0], query[1], 1))

        copies.sort(key=lambda x: -(x[1]-x[0]))

        pure_copies = []
        for cur_copy in copies:
            if cur_copy[2] < 2 and float(cur_copy[1] - cur_copy[0]) / query_len < 0.95:
                continue
            is_seg = False
            for i, pure_copy in enumerate(pure_copies):
                if cur_copy[0] >= pure_copy[0] and abs(cur_copy[1]-pure_copy[1]) < 50:
                    pure_copies[i] = (pure_copy[0], max(cur_copy[1], pure_copy[1]), pure_copy[2])
                    is_seg = True
                    break
                elif cur_copy[1] <= pure_copy[1] and abs(cur_copy[0]-pure_copy[0]) < 50:
                    pure_copies[i] = (min(cur_copy[0], pure_copy[0]), pure_copy[1], pure_copy[2])
                    is_seg = True
                    break
            if not is_seg:
                pure_copies.append(cur_copy)

        intact_copies = []
        for copy in pure_copies:
            start_pos = copy[0] - 1
            end_pos = copy[1]
            ori_start_pos = start_pos
            ori_end_pos = end_pos
            if is_flank:
                end_pos += 2*flanking_len
            copy_seq = query_seq[start_pos: end_pos]

            seq_ref_start = chr_start + ori_start_pos
            if is_flank:
                seq_ref_end = chr_start + ori_end_pos + 2*flanking_len
            else:
                seq_ref_end = chr_start + ori_end_pos

            if len(copy_seq) >= 80:
                intact_copies.append((ori_start_pos, ori_end_pos, chr_name, seq_ref_start, seq_ref_end, copy_seq))
        longest_repeats[query_name] = intact_copies

    # print(longest_repeats)
    return longest_repeats, keep_longest_query

def get_longest_repeats_v3(repeats_path, fixed_extend_base_threshold, max_single_repeat_len):
    split_repeats_path = repeats_path[0]
    target_files = repeats_path[1]
    blastn2Results_path = repeats_path[2]
    tmp_blast_dir = repeats_path[3]

    for i, cur_target in enumerate(target_files):
        align_command = 'blastn -db ' + cur_target + ' -num_threads ' \
                        + str(1) + ' -query ' + split_repeats_path + ' -outfmt 6'
        if i == 0:
            align_command += ' > ' + blastn2Results_path
        else:
            align_command += ' >> ' + blastn2Results_path
        #os.system(align_command)

    blastn2Results_path = '/home/hukang/HiTE/demo/test/total_blast_0.out'
    #print('coarse alignment -- alignment finished:' + str(split_repeats_path))

    query_names, query_contigs = read_fasta(split_repeats_path)

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
            if identity < 80 or (query_name==subject_name and q_start==s_start and q_end==s_end):
                continue
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))
    f_r.close()
    #print('coarse alignment -- file open finished:' + str(split_repeats_path))

    keep_longest_query = {}
    longest_repeats = {}
    for idx, query_name in enumerate(query_records.keys()):
        #query_len = len(query_contigs[query_name])
        #print('total query size: %d, current query name: %s, idx: %d' % (len(query_records), query_name, idx))

        subject_dict = query_records[query_name]

        # if there are more than one longest query overlap with the final longest query over 90%,
        # then it probably the true TE
        longest_queries = []
        for subject_name in subject_dict.keys():
            if query_name != 'chr_0$15000000' or subject_name != 'chr_0$15000000':
                continue
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
                        if (frag[2] - exist_frag[3] < fixed_extend_base_threshold):
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
                        if (exist_frag[3] - frag[2] < fixed_extend_base_threshold):
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

                cluster_longest_subject_start = -1
                cluster_longest_subject_end = -1
                cluster_longest_subject_len = -1

                cluster_extend_num = 0

                # print('subject pos size: %d' %(len(cur_cluster)))
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

                    # if longest_query_start == 456006 and longest_query_end == 457161:
                    #     print('here')

                    cur_extend_num = 0

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
                                    if ext_frag[0] - longest_query_end < fixed_extend_base_threshold and ext_frag[
                                        2] - longest_subject_end < fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= fixed_extend_base_threshold:
                                        break
                            elif longest_subject_start > longest_subject_end and ext_frag[2] > ext_frag[3]:
                                # reverse
                                if ext_frag[3] < longest_subject_end:
                                    # reverse extend
                                    if ext_frag[
                                        0] - longest_query_end < fixed_extend_base_threshold and longest_subject_end - \
                                            ext_frag[2] < fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start > \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= fixed_extend_base_threshold:
                                        break
                    if cur_longest_query_len > cluster_longest_query_len:
                        cluster_longest_query_start = longest_query_start
                        cluster_longest_query_end = longest_query_end
                        cluster_longest_query_len = cur_longest_query_len

                        cluster_longest_subject_start = longest_subject_start
                        cluster_longest_subject_end = longest_subject_end
                        cluster_longest_subject_len = longest_subject_end - longest_subject_start

                        cluster_extend_num = cur_extend_num
                # keep this longest query
                if cluster_longest_query_len != -1:
                    longest_queries.append((cluster_longest_query_start, cluster_longest_query_end,
                                            cluster_longest_query_len, cluster_longest_subject_start,
                                            cluster_longest_subject_end, cluster_longest_subject_len, subject_name, cluster_extend_num))


        # we now consider, we should take some sequences from longest_queries to represent this query sequence.
        # we take the longest sequence by length, if the latter sequence overlap with the former sequence largely (50%),
        # continue find next sequence until the ratio of query sequence over 90% or no more sequences.
        longest_queries.sort(key=lambda x: -x[2])
        keep_longest_query[query_name] = longest_queries
        parts = query_name.split('$')
        chr_name = parts[0]
        chr_start = int(parts[1])

        is_flank = False
        flanking_len = 0

        query_seq = query_contigs[query_name]

        # 计算所有片段的支持拷贝数，边界只保留最长边界
        copies = []
        for query in longest_queries:
            is_copy = False
            for i in range(len(copies)-1, -1, -1):
                copy = copies[i]
                #计算overlap
                if copy[0] <= query[1] and copy[0] >= query[0]:
                    overlap = query[1] - copy[0]
                elif query[0] >= copy[0] and query[1] <= copy[1]:
                    overlap = query[1] - query[0]
                elif copy[1] <= query[1] and copy[1] >= query[0]:
                    overlap = copy[1] - query[0]
                else:
                    overlap = 0
                if overlap > 0:
                    if float(overlap) / (query[1]-query[0]) >= 0.95:
                        if float(overlap) / (copy[1]-copy[0]) >= 0.95:
                            copies[i] = (min(copy[0], query[0]), max(copy[1], query[1]), copy[2] + 1)
                            is_copy = True
                            break
                        else:
                            is_copy = False
                            break
                    else:
                        is_copy = False
                        break
            if not is_copy and abs(query[1]-query[0]) <= max_single_repeat_len:
                copies.append((query[0], query[1], 1))

        copies.sort(key=lambda x: -(x[1]-x[0]))

        pure_copies = []
        for cur_copy in copies:
            # if cur_copy[2] < 2 and float(cur_copy[1] - cur_copy[0]) / query_len < 0.95:
            #     continue
            is_seg = False
            for i, pure_copy in enumerate(pure_copies):
                if cur_copy[0] >= pure_copy[0] and abs(cur_copy[1]-pure_copy[1]) < 50:
                    pure_copies[i] = (pure_copy[0], max(cur_copy[1], pure_copy[1]), pure_copy[2])
                    is_seg = True
                    break
                elif cur_copy[1] <= pure_copy[1] and abs(cur_copy[0]-pure_copy[0]) < 50:
                    pure_copies[i] = (min(cur_copy[0], pure_copy[0]), pure_copy[1], pure_copy[2])
                    is_seg = True
                    break
            if not is_seg:
                pure_copies.append(cur_copy)

        intact_copies = []
        for copy in pure_copies:
            start_pos = copy[0] - 1
            end_pos = copy[1]
            ori_start_pos = start_pos
            ori_end_pos = end_pos
            if is_flank:
                end_pos += 2*flanking_len
            copy_seq = query_seq[start_pos: end_pos]

            seq_ref_start = chr_start + ori_start_pos
            if is_flank:
                seq_ref_end = chr_start + ori_end_pos + 2*flanking_len
            else:
                seq_ref_end = chr_start + ori_end_pos

            if len(copy_seq) >= 80:
                intact_copies.append((ori_start_pos, ori_end_pos, chr_name, seq_ref_start, seq_ref_end, copy_seq))
        # if query_name == 'chr_0$15000000' or subject_name == 'chr_0$15000000':
        #     print(intact_copies[0])
        longest_repeats[query_name] = intact_copies

    #print('coarse alignment -- analyze finished:' + str(split_repeats_path))
    # print(longest_repeats)
    return longest_repeats, keep_longest_query, blastn2Results_path

def get_alignment_records(repeats_path):
    split_repeats_path = repeats_path[0]
    target_files = repeats_path[1]
    blastn2Results_path = repeats_path[2]

    for i, cur_target in enumerate(target_files):
        align_command = 'blastn -db ' + cur_target + ' -num_threads ' \
                        + str(1) + ' -evalue 1e-20 -query ' + split_repeats_path + ' -outfmt 6'
        if i == 0:
            align_command += ' > ' + blastn2Results_path
        else:
            align_command += ' >> ' + blastn2Results_path
        os.system(align_command)
    return blastn2Results_path

def get_longest_repeats_v5(query_records, fixed_extend_base_threshold, max_single_repeat_len, chr_pos_candidates, debug):
    # The purpose of this function is actually to find possible
    # candidate repeat sequences according to the comparison information.
    # longest_repeats -> {'chr1:100-1000': seq, }
    longest_repeats = {}
    for idx, query_name in enumerate(query_records.keys()):
        subject_dict = query_records[query_name]

        longest_queries = []
        for subject_name in subject_dict.keys():
            #judge whether the alignments on each subject can be merged.
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
                        if (frag[2] - exist_frag[3] < fixed_extend_base_threshold):
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
                        if (exist_frag[3] - frag[2] < fixed_extend_base_threshold):
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

                # print('subject pos size: %d' %(len(cur_cluster)))
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

                    cur_extend_num = 0

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
                                    if ext_frag[0] - longest_query_end < fixed_extend_base_threshold and ext_frag[
                                        2] - longest_subject_end < fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= fixed_extend_base_threshold:
                                        break
                            elif longest_subject_start > longest_subject_end and ext_frag[2] > ext_frag[3]:
                                # reverse
                                if ext_frag[3] < longest_subject_end:
                                    # reverse extend
                                    if ext_frag[
                                        0] - longest_query_end < fixed_extend_base_threshold and longest_subject_end - \
                                            ext_frag[2] < fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start > \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= fixed_extend_base_threshold:
                                        break
                    # keep this longest query
                    if cur_longest_query_len != -1:
                        longest_queries.append((longest_query_start, longest_query_end,
                                                cur_longest_query_len, longest_subject_start,
                                                longest_subject_end, longest_subject_end - longest_subject_start, subject_name, cur_extend_num))

        parts = query_name.split('$')
        query_chr_name = parts[0]
        query_chr_start = int(parts[1])

        query_repeatNames = []
        query_repeats = {}
        for repeat in longest_queries:
            subject_name = repeat[6]
            parts = subject_name.split('$')
            subject_chr_name = parts[0]
            subject_chr_start = int(parts[1])

            old_subject_start_pos = repeat[3] - 1
            old_subject_end_pos = repeat[4]
            subject_start_pos = subject_chr_start + old_subject_start_pos
            subject_end_pos = subject_chr_start + old_subject_end_pos
            subject_pos = subject_chr_name + ':' + str(subject_start_pos) + '-' + str(subject_end_pos)

            old_query_start_pos = repeat[0] - 1
            old_query_end_pos = repeat[1]

            cur_seq_len = old_query_end_pos - old_query_start_pos
            query_start_pos = query_chr_start + old_query_start_pos
            query_end_pos = query_chr_start + old_query_end_pos
            query_pos = query_chr_name + ':' + str(query_start_pos) + '-' + str(query_end_pos)

            if cur_seq_len >= 80 and cur_seq_len < max_single_repeat_len:
                if not chr_pos_candidates.__contains__(query_pos) and not chr_pos_candidates.__contains__(subject_pos):
                    query_repeatNames.append(query_pos)
                    query_repeats[query_pos] = (query_name, old_query_start_pos, old_query_end_pos)
                chr_pos_candidates[query_pos] = 1
                chr_pos_candidates[subject_pos] = 1

        merge_contigNames = process_all_seqs(query_repeatNames)
        for query_pos in merge_contigNames:
            longest_repeats[query_pos] = query_repeats[query_pos]
    return longest_repeats

def get_longest_repeats_v4(repeats_path, fixed_extend_base_threshold, max_single_repeat_len, debug):
    split_repeats_path = repeats_path[0]
    target_files = repeats_path[1]
    blastn2Results_path = repeats_path[2]

    subject_contigs = {}
    for i, cur_target in enumerate(target_files):
        align_command = 'blastn -db ' + cur_target + ' -num_threads ' \
                        + str(1) + ' -evalue 1e-20 -query ' + split_repeats_path + ' -outfmt 6'
        if i == 0:
            align_command += ' > ' + blastn2Results_path
        else:
            align_command += ' >> ' + blastn2Results_path
        os.system(align_command)
        cur_subject_names, cur_subject_contigs = read_fasta(cur_target)
        subject_contigs.update(cur_subject_contigs)

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
            if (query_name==subject_name and q_start==s_start and q_end==s_end):
                continue
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))
    f_r.close()

    skip_gap = fixed_extend_base_threshold
    keep_longest_query = {}
    # longest_repeats -> {'chr1:100-1000': seq, }
    longest_repeats = {}
    # Recording all possible integer positions, for example, chr1, start: 98, end: 995, will yield the following four coordinate sets:
    # chr_pos_candidates -> {'chr1:90-990': 1, 'chr1:90-1000': 1, 'chr1:100-990': 1, 'chr1:100-1000': 1}
    chr_pos_candidates = {}
    for idx, query_name in enumerate(query_records.keys()):
        subject_dict = query_records[query_name]

        # if there are more than one longest query overlap with the final longest query over 90%,
        # then it probably the true TE
        longest_queries = []
        for subject_name in subject_dict.keys():
            subject_pos = subject_dict[subject_name]

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
                        cur_subject_start = frag[2]
                        cur_query_end = frag[1]
                        prev_subject_end = exist_frag[3]
                        prev_query_end = exist_frag[1]
                        if (cur_subject_start - prev_subject_end < skip_gap and cur_query_end > prev_query_end):
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
                        cur_subject_start = frag[2]
                        cur_query_end = frag[1]
                        prev_subject_end = exist_frag[3]
                        prev_query_end = exist_frag[1]
                        if (prev_subject_end - cur_subject_start < skip_gap and cur_query_end > prev_query_end):
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

                # print('subject pos size: %d' %(len(cur_cluster)))
                # record visited fragments
                visited_frag = {}
                for i in range(len(cur_cluster)):
                    # keep a longest query start from each fragment
                    prev_frag = cur_cluster[i]
                    if visited_frag.__contains__(prev_frag):
                        continue
                    prev_query_start = prev_frag[0]
                    prev_query_end = prev_frag[1]
                    prev_subject_start = prev_frag[2]
                    prev_subject_end = prev_frag[3]
                    prev_query_seq = (min(prev_query_start, prev_query_end), max(prev_query_start, prev_query_end))
                    prev_subject_seq = (
                    min(prev_subject_start, prev_subject_end), max(prev_subject_start, prev_subject_end))
                    prev_query_len = abs(prev_query_end - prev_query_start)
                    prev_subject_len = abs(prev_subject_end - prev_subject_start)
                    cur_longest_query_len = prev_query_len

                    cur_extend_num = 0
                    visited_frag[prev_frag] = 1
                    # try to extend query
                    for j in range(i + 1, len(cur_cluster)):
                        cur_frag = cur_cluster[j]
                        if visited_frag.__contains__(cur_frag):
                            continue
                        cur_query_start = cur_frag[0]
                        cur_query_end = cur_frag[1]
                        cur_subject_start = cur_frag[2]
                        cur_subject_end = cur_frag[3]
                        cur_query_seq = (min(cur_query_start, cur_query_end), max(cur_query_start, cur_query_end))
                        cur_subject_seq = (
                        min(cur_subject_start, cur_subject_end), max(cur_subject_start, cur_subject_end))
                        cur_query_len = abs(cur_query_end - cur_query_start)
                        cur_subject_len = abs(cur_subject_end - cur_subject_start)

                        query_overlap_len = get_overlap_len(cur_query_seq, prev_query_seq)
                        is_same_query = float(query_overlap_len) / cur_query_len >= 0.5 or float(
                            query_overlap_len) / prev_query_len >= 0.5
                        subject_overlap_len = get_overlap_len(prev_subject_seq, cur_subject_seq)
                        is_same_subject = float(subject_overlap_len) / cur_subject_len >= 0.5 or float(
                            subject_overlap_len) / prev_subject_len >= 0.5

                        # could extend
                        # extend right
                        if cur_query_end > prev_query_end:
                            # judge subject direction
                            if prev_subject_start < prev_subject_end and cur_subject_start < cur_subject_end:
                                # +
                                if cur_subject_end > prev_subject_end:
                                    # forward extend
                                    if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                            and cur_subject_start - prev_subject_end < skip_gap:  # \
                                        # and not is_same_query and not is_same_subject:
                                        # update the longest path
                                        prev_query_start = prev_query_start
                                        prev_query_end = cur_query_end
                                        prev_subject_start = prev_subject_start if prev_subject_start < cur_subject_start else cur_subject_start
                                        prev_subject_end = cur_subject_end
                                        cur_longest_query_len = prev_query_end - prev_query_start
                                        cur_extend_num += 1
                                        visited_frag[cur_frag] = 1
                                    elif cur_query_start - prev_query_end >= skip_gap:
                                        break
                            elif prev_subject_start > prev_subject_end and cur_subject_start > cur_subject_end:
                                # reverse
                                if cur_subject_end < prev_subject_end:
                                    # reverse extend
                                    if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                            and prev_subject_end - cur_subject_start < skip_gap:  # \
                                        # and not is_same_query and not is_same_subject:
                                        # update the longest path
                                        prev_query_start = prev_query_start
                                        prev_query_end = cur_query_end
                                        prev_subject_start = prev_subject_start if prev_subject_start > cur_subject_start else cur_subject_start
                                        prev_subject_end = cur_subject_end
                                        cur_longest_query_len = prev_query_end - prev_query_start
                                        cur_extend_num += 1
                                        visited_frag[cur_frag] = 1
                                    elif cur_query_start - prev_query_end >= skip_gap:
                                        break
                    # keep this longest query
                    if cur_longest_query_len != -1:
                        longest_queries.append(
                            (prev_query_start, prev_query_end, cur_longest_query_len, prev_subject_start,
                             prev_subject_end, abs(prev_subject_end - prev_subject_start), subject_name,
                             cur_extend_num))

        # Retrieve all possible repetitive regions, label them by (chromosome:start-end), and record them.
        # For sequences that are repeated again, refrain from saving them again to reduce subsequent computational load.
        # In other words, each repetitive segment is kept only once.
        parts = query_name.split('$')
        query_chr_name = parts[0]
        query_chr_start = int(parts[1])

        query_repeatNames = []
        query_repeats = {}
        for repeat in longest_queries:
            # Subject
            subject_name = repeat[6]
            parts = subject_name.split('$')
            subject_chr_name = parts[0]
            subject_chr_start = int(parts[1])

            old_subject_start_pos = repeat[3] - 1
            old_subject_end_pos = repeat[4]
            subject_start_pos = subject_chr_start + old_subject_start_pos
            subject_end_pos = subject_chr_start + old_subject_end_pos

            # Rounding the positional coordinates to the nearest multiple of 10,
            # for example, 567, we obtain two numbers: 560 and 570.
            subject_start_pos1, subject_start_pos2 = get_integer_pos(subject_start_pos)
            subject_end_pos1, subject_end_pos2 = get_integer_pos(subject_end_pos)
            # Record the chromosomal coordinates corresponding to this segment.
            subject_pos = subject_chr_name + ':' + str(subject_start_pos) + '-' + str(subject_end_pos)
            subject_pos1 = subject_chr_name + ':' + str(subject_start_pos1) + '-' + str(subject_end_pos1)
            subject_pos2 = subject_chr_name + ':' + str(subject_start_pos1) + '-' + str(subject_end_pos2)
            subject_pos3 = subject_chr_name + ':' + str(subject_start_pos2) + '-' + str(subject_end_pos1)
            subject_pos4 = subject_chr_name + ':' + str(subject_start_pos2) + '-' + str(subject_end_pos2)

            # Query
            old_query_start_pos = repeat[0] - 1
            old_query_end_pos = repeat[1]
            query_start_pos = query_chr_start + old_query_start_pos
            query_end_pos = query_chr_start + old_query_end_pos
            cur_query_seq_len = abs(old_query_end_pos - old_query_start_pos)

            query_start_pos1, query_start_pos2 = get_integer_pos(query_start_pos)
            query_end_pos1, query_end_pos2 = get_integer_pos(query_end_pos)

            query_pos = query_chr_name + ':' + str(query_start_pos) + '-' + str(query_end_pos)
            query_pos1 = query_chr_name + ':' + str(query_start_pos1) + '-' + str(query_end_pos1)
            query_pos2 = query_chr_name + ':' + str(query_start_pos1) + '-' + str(query_end_pos2)
            query_pos3 = query_chr_name + ':' + str(query_start_pos2) + '-' + str(query_end_pos1)
            query_pos4 = query_chr_name + ':' + str(query_start_pos2) + '-' + str(query_end_pos2)

            # Determine if this sequence has already been saved. If not, record the sequence of this subject segment.
            if not chr_pos_candidates.__contains__(subject_pos1) \
                    and not chr_pos_candidates.__contains__(subject_pos2) \
                    and not chr_pos_candidates.__contains__(subject_pos3) \
                    and not chr_pos_candidates.__contains__(subject_pos4) \
                    and not chr_pos_candidates.__contains__(query_pos1) \
                    and not chr_pos_candidates.__contains__(query_pos2) \
                    and not chr_pos_candidates.__contains__(query_pos3) \
                    and not chr_pos_candidates.__contains__(query_pos4):

                if cur_query_seq_len >= 80 and cur_query_seq_len < max_single_repeat_len:
                    query_repeatNames.append(query_pos)
                    query_repeats[query_pos] = 1

            chr_pos_candidates[subject_pos1] = 1
            chr_pos_candidates[subject_pos2] = 1
            chr_pos_candidates[subject_pos3] = 1
            chr_pos_candidates[subject_pos4] = 1
            chr_pos_candidates[query_pos1] = 1
            chr_pos_candidates[query_pos2] = 1
            chr_pos_candidates[query_pos3] = 1
            chr_pos_candidates[query_pos4] = 1

        # Merging query repeats, if the overlap exceeds 95%, then discard.
        merge_contigNames = process_all_seqs(query_repeatNames)
        for name in merge_contigNames:
            longest_repeats[name] = query_repeats[name]
    return longest_repeats, keep_longest_query, blastn2Results_path

def get_overlap_len(seq1, seq2):
    """Calculate the overlap length between two sequences."""
    overlap_len = min(seq1[1], seq2[1]) - max(seq1[0], seq2[0])
    return overlap_len if overlap_len > 0 else 0

def get_gap_len(seq1, seq2):
    """Calculate the gap length between two sequences."""
    gap_len = max(seq1[0], seq2[0]) - min(seq1[1], seq2[1])
    return gap_len if gap_len > 0 else 0

def get_chr_pos(query_name, subject_name, query_frag, subject_frag):
    parts = query_name.split('$')
    query_chr_name = parts[0]
    query_chr_start = int(parts[1])
    old_query_start_pos = query_frag[0] - 1
    old_query_end_pos = query_frag[1]
    query_start_pos = query_chr_start + old_query_start_pos
    query_end_pos = query_chr_start + old_query_end_pos
    query_pos = query_chr_name + ':' + str(query_start_pos) + '-' + str(query_end_pos)
    # Record the chromosomal coordinates corresponding to the subject segment.
    parts = subject_name.split('$')
    subject_chr_name = parts[0]
    subject_chr_start = int(parts[1])
    old_subject_start_pos = subject_frag[0] - 1
    old_subject_end_pos = subject_frag[1]
    subject_start_pos = subject_chr_start + old_subject_start_pos
    subject_end_pos = subject_chr_start + old_subject_end_pos
    subject_pos = subject_chr_name + ':' + str(subject_start_pos) + '-' + str(subject_end_pos)
    return query_pos, subject_pos

def merge_same_fragments(all_fragments):
    combine_repeats = {}
    longest_repeats = {}
    merged_fragments = {}
    for query_name in all_fragments.keys():
        if not combine_repeats.__contains__(query_name):
            combine_repeats[query_name] = []
        cur_combine_repeats = combine_repeats[query_name]

        if not longest_repeats.__contains__(query_name):
            longest_repeats[query_name] = []
        cur_longest_repeats = longest_repeats[query_name]

        if not merged_fragments.__contains__(query_name):
            merged_fragments[query_name] = []
        cur_merge_fragments = merged_fragments[query_name]

        cur_query_fragments = list(all_fragments[query_name])

        cur_query_fragments.sort(key=lambda x: (x[0], x[1]))

        keep_fragments_clusters = []
        # Maintain the maximum end of a keep fragments.
        # If the start of the current fragment is greater than the current maximum end,
        # it means that the current fragment cannot overlap with the reserved fragment, so skip.
        keep_frag_max_end = 0
        for i, cur_frag in enumerate(cur_query_fragments):
            is_merge = False
            for j, cluster in reversed(list(enumerate(keep_fragments_clusters))):
                for prev_frag in cluster:
                    # judge the max end of the current fragment and the saved fragment,
                    # and when the start of the current fragment > max end, there is no overlap.
                    if cur_frag[0] > keep_frag_max_end:
                        break
                    overlap_len = get_overlap_len(prev_frag, cur_frag)
                    if overlap_len / abs(prev_frag[1] - prev_frag[0]) >= 0.95 or overlap_len / abs(cur_frag[1] - cur_frag[0]) >= 0.95:
                        is_merge = True
                        # Update the max end corresponding to the currently saved fragment.
                        if cur_frag[1] > keep_frag_max_end:
                            keep_frag_max_end = cur_frag[1]
                        break
                if is_merge:
                    # update the overlap_sum of all fragments in current cluster
                    total_consensus_score = 0
                    total_count = 0
                    for prev_frag in cluster:
                        overlap_len = get_overlap_len(prev_frag, cur_frag)
                        coverage_prev_frag = float(overlap_len) / abs(prev_frag[1] - prev_frag[0])
                        coverage_cur_frag = float(overlap_len) / abs(cur_frag[1] - cur_frag[0])
                        if coverage_prev_frag >= 0.95 and coverage_cur_frag >= 0.95:
                            consensus_score = 1 - abs(coverage_cur_frag - coverage_prev_frag)
                            prev_frag[2] += consensus_score
                            prev_frag[3] += 1
                            total_count += 1
                        else:
                            consensus_score = 0
                        total_consensus_score += consensus_score
                    # add into current cluster
                    cluster.append([cur_frag[0], cur_frag[1], total_consensus_score, total_count])
                    break
            if not is_merge:
                # Add the current fragment to the new cluster and set overlap_sum to 0.
                keep_fragments_clusters.append([[cur_frag[0], cur_frag[1], 0, 0]])
                # Save the max end corresponding to the current saved fragment.
                if cur_frag[1] > keep_frag_max_end:
                    keep_frag_max_end = cur_frag[1]
        # Since we cannot determine if the highest scoring sequence represents the true boundary or
        # if the longest sequence represents the true boundary, we retain both.
        for cluster in keep_fragments_clusters:
            #print(cluster)
            max_longest_length = -1
            longest_repeat_frag = None
            max_consensus_score = -1
            max_consensus_frag = None
            for frag in cluster:
                if frag[2] > max_consensus_score:
                    max_consensus_frag = frag
                    max_consensus_score = frag[2]
                frag_len = abs(frag[1]-frag[0])
                if frag_len > max_longest_length:
                    longest_repeat_frag = frag
                    max_longest_length = frag_len
            if max_consensus_frag is not None:
                cur_merge_fragments.append(max_consensus_frag)
            if longest_repeat_frag is not None:
                cur_longest_repeats.append(longest_repeat_frag)
            # If the coverage between max_consensus_frag and longest_repeat_frag does not exceed 95%, both are kept.
            overlap_len = get_overlap_len(max_consensus_frag, longest_repeat_frag)
            max_consensus_frag_coverage = float(overlap_len) / abs(max_consensus_frag[1] - max_consensus_frag[0])
            longest_repeat_frag_coverage = float(overlap_len) / abs(longest_repeat_frag[1] - longest_repeat_frag[0])
            if max_consensus_frag_coverage >= 0.95 and longest_repeat_frag_coverage >= 0.95:
                cur_combine_repeats.append(max_consensus_frag)
            else:
                cur_combine_repeats.append(max_consensus_frag)
                cur_combine_repeats.append(longest_repeat_frag)
    return merged_fragments, longest_repeats, combine_repeats

def process_seq_group(seq_group):
    seq_group.sort(key=lambda x: x[1] - x[0], reverse=True) # Sort by length.
    keep_seq = [True] * len(seq_group)  # Store a boolean list indicating whether each sequence should be retained.
    for i in range(len(seq_group)):
        if not keep_seq[i]:  # If the current sequence has been marked for deletion, skip it.
            continue
        seq1 = seq_group[i]
        for j in range(i + 1, len(seq_group)):
            if not keep_seq[j]:  # If the next sequence has been marked for deletion, skip it.
                continue
            seq2 = seq_group[j]
            # Check if the current sequence has over 95% overlap with the next sequence.
            overlap_len = get_overlap_len(seq2, seq1)
            if overlap_len / (seq2[1] - seq2[0]) >= 0.95:
                keep_seq[j] = False  # Mark the next sequence for deletion.
    # Return the retained sequences.
    # filter_res = [seq for i, seq in enumerate(seq_group) if keep_seq[i]]
    # filter_res.sort(key=lambda x: (x[0], x[1]))
    # res = [seq[2] for i, seq in enumerate(filter_res)]
    res = [seq[2] for i, seq in enumerate(seq_group) if keep_seq[i]]
    return res

def process_all_seqs(seq_list):
    seq_dict = {}
    for seq in seq_list:
        chrom, pos = seq.split(":")
        start, end = pos.split("-")
        start, end = int(start), int(end)
        if chrom not in seq_dict:
            seq_dict[chrom] = []
        seq_dict[chrom].append((start, end, seq))
    result = []
    for chrom, seq_group in seq_dict.items():
        result.extend(process_seq_group(seq_group))
    return result

# Take the position coordinate as an integer of 10, such as 567, and we get two numbers 560 and 570.
def get_integer_pos(x):
    y = (x // 10) * 10
    z = y + 10
    return y, z

def get_domain_info(cons, lib, output_table, threads, temp_dir):
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    blast_db_command = 'makeblastdb -dbtype prot -in ' + lib + ' > /dev/null 2>&1'
    os.system(blast_db_command)
    # 1. blastx -num_threads 1 -evalue 1e-20
    partitions_num = int(threads)
    consensus_contignames, consensus_contigs = read_fasta(cons)
    data_partitions = PET(consensus_contigs.items(), partitions_num)
    merge_distance = 100
    file_list = []
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for partition_index, data_partition in enumerate(data_partitions):
        if len(data_partition) <= 0:
            continue
        cur_consensus_path = temp_dir + '/'+str(partition_index)+'.fa'
        store2file(data_partition, cur_consensus_path)
        cur_output = temp_dir + '/'+str(partition_index)+'.out'
        cur_table = temp_dir + '/' + str(partition_index) + '.tbl'
        cur_file = (cur_consensus_path, lib, cur_output, cur_table)
        job = ex.submit(multiple_alignment_blastx_v1, cur_file, merge_distance)
        jobs.append(job)
    ex.shutdown(wait=True)

    # 2. generate table of query and domain
    os.system("echo 'TE_name\tdomain_name\tTE_start\tTE_end\tdomain_start\tdomain_end\n' > " + output_table)
    for job in as_completed(jobs):
        cur_table = job.result()
        os.system('cat ' + cur_table + ' >> ' + output_table)

def flanking_seq(longest_repeats_path, longest_repeats_flanked_path, reference, flanking_len):
    seq_names, seq_contigs = read_fasta(longest_repeats_path)
    ref_names, ref_contigs = read_fasta(reference)
    flanked_contigs = {}
    node_index = 0
    for name in seq_names:
        # chr_0:10023581-10023756
        ref_info = name.split(':')
        ref_name = ref_info[0]
        pos_info = ref_info[1].split('-')
        ref_start = int(pos_info[0]) + 1
        ref_end = int(pos_info[1])
        ref_seq = ref_contigs[ref_name]
        if ref_start - 1 - flanking_len < 0:
            ref_start = flanking_len + 1
        if ref_end + flanking_len > len(ref_seq):
            ref_end = len(ref_seq) - flanking_len
        flanked_seq = ref_seq[ref_start - 1 - flanking_len: ref_end + flanking_len]
        new_name = ref_name + ':' + str(ref_start-flanking_len) + '-' + str(ref_end + flanking_len)
        flanked_contigs[new_name] = flanked_seq
    store_fasta(flanked_contigs, longest_repeats_flanked_path)

def determine_repeat_boundary_v5(repeats_path, longest_repeats_path, prev_TE, fixed_extend_base_threshold,
                                 max_single_repeat_len, tmp_output_dir, threads, ref_index, reference, debug):
    repeatNames, repeatContigs = read_fasta(repeats_path)
    tmp_blast_dir = tmp_output_dir + '/trf_filter_'+str(ref_index)
    os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    # We found that the main reason for the memory spike in blastn alignment is the presence of numerous short tandem repeats, leading to complex alignments.
    # Since STRs are not our target, we should preprocess and filter out this part of the data.
    # Use TRF to remove tandem repeats from target files.
    # Sequence alignment consumes a significant amount of memory and disk space. Therefore, we also split the target sequences into individual sequences to reduce the memory required for each alignment, avoiding out of memory errors.
    # It is important to calculate the total number of bases in the sequences, and it must meet a sufficient threshold to increase CPU utilization.
    # 1. Split the sequences into 1Mb files to reduce TRF execution time.
    base_threshold = 1000000 #1Mb
    all_files = []
    file_index = 0
    base_count = 0
    cur_contigs = {}
    for name in repeatNames:
        cur_seq = repeatContigs[name]
        cur_contigs[name] = cur_seq
        base_count += len(cur_seq)
        if base_count >= base_threshold:
            cur_target = tmp_blast_dir + '/' + str(file_index) + '_target.fa'
            cur_trf = tmp_blast_dir + '/' + str(file_index) + '_trf'
            store_fasta(cur_contigs, cur_target)
            all_files.append((cur_target, cur_trf))
            cur_contigs = {}
            file_index += 1
            base_count = 0
    if len(cur_contigs) > 0:
        cur_target = tmp_blast_dir + '/' + str(file_index) + '_target.fa'
        cur_trf = tmp_blast_dir + '/' + str(file_index) + '_trf'
        store_fasta(cur_contigs, cur_target)
        all_files.append((cur_target, cur_trf))

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in all_files:
        cur_target = file[0]
        cur_trf = file[1]
        job = ex.submit(run_remove_TR, cur_target, cur_trf)
        jobs.append(job)
    ex.shutdown(wait=True)

    filter_tandem_file = tmp_output_dir + '/filter_tandem_' + str(ref_index) + '.fa'
    if os.path.exists(filter_tandem_file):
        os.system('rm -f ' + filter_tandem_file)
    for job in as_completed(jobs):
        cur_target_file = job.result()
        os.system('cat ' + cur_target_file + ' >> ' + filter_tandem_file)

    if debug != 1:
        os.system('rm -rf ' + tmp_blast_dir)

    # # 2. Mask the genome with the TEs identified in the previous rounds of HiTE.
    # RepeatMasker_command = 'cd ' + tmp_output_dir + ' && RepeatMasker -e ncbi -pa ' + str(threads) \
    #                        + ' -q -no_is -norna -nolow -div 20 -gff -lib ' + prev_TE + ' -cutoff 225 ' \
    #                        + filter_tandem_file
    # os.system(RepeatMasker_command)
    # masked_file_path = filter_tandem_file + '.masked'

    # 2. Mask the genome with the TEs identified in the previous rounds of HiTE.
    mask_genome_intactTE(prev_TE, filter_tandem_file, tmp_output_dir, threads, ref_index)
    masked_file_path = filter_tandem_file + '.masked'

    # 2.1. Remove sequences consisting entirely of 'N'.
    filter_tandem_file = tmp_output_dir + '/filter_TE_' + str(ref_index) + '.fa'
    masked_contigNames, masked_contigs = read_fasta(masked_file_path)
    filter_contigs = {}
    for name in masked_contigNames:
        seq = masked_contigs[name]
        if all(c == "N" for c in seq):
            continue
        filter_contigs[name] = seq
    store_fasta(filter_contigs, filter_tandem_file)

    # 3.Collect the filtered concatenated sequences and divide them into 10Mb files to enhance blastn CPU utilization.
    repeatNames, repeatContigs = read_fasta(filter_tandem_file)
    tmp_blast_dir = tmp_output_dir + '/longest_repeats_blast_' + str(ref_index)
    os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    # Sequence alignment consumes a significant amount of memory and disk space. Therefore, we also split the target sequences into individual sequences to reduce the memory required for each alignment, avoiding out of memory errors.
    # It is important to calculate the total number of bases in the sequences, and it must meet a sufficient threshold to increase CPU utilization.
    base_threshold = 1000000  # 1Mb
    target_files = []
    file_index = 0
    base_count = 0
    cur_contigs = {}
    for name in repeatNames:
        cur_seq = repeatContigs[name]
        cur_contigs[name] = cur_seq
        base_count += len(cur_seq)
        if base_count >= base_threshold:
            cur_target = tmp_blast_dir + '/' + str(file_index) + '_target.fa'
            store_fasta(cur_contigs, cur_target)
            target_files.append(cur_target)
            makedb_command = 'makeblastdb -dbtype nucl -in ' + cur_target + ' > /dev/null 2>&1'
            os.system(makedb_command)
            cur_contigs = {}
            file_index += 1
            base_count = 0
    if len(cur_contigs) > 0:
        cur_target = tmp_blast_dir + '/' + str(file_index) + '_target.fa'
        store_fasta(cur_contigs, cur_target)
        target_files.append(cur_target)
        makedb_command = 'makeblastdb -dbtype nucl -in ' + cur_target + ' > /dev/null 2>&1'
        os.system(makedb_command)

    repeat_files = []
    file_index = 0
    base_count = 0
    cur_contigs = {}
    for name in repeatNames:
        cur_seq = repeatContigs[name]
        cur_contigs[name] = cur_seq
        base_count += len(cur_seq)
        if base_count >= base_threshold:
            split_repeat_file = tmp_blast_dir + '/' + str(file_index) + '.fa'
            store_fasta(cur_contigs, split_repeat_file)
            output_file = tmp_blast_dir + '/' + str(file_index) + '.out'
            repeat_files.append((split_repeat_file, target_files, output_file))
            cur_contigs = {}
            file_index += 1
            base_count = 0
    if len(cur_contigs) > 0:
        split_repeat_file = tmp_blast_dir + '/' + str(file_index) + '.fa'
        store_fasta(cur_contigs, split_repeat_file)
        output_file = tmp_blast_dir + '/' + str(file_index) + '.out'
        repeat_files.append((split_repeat_file, target_files, output_file))

    ex = ProcessPoolExecutor(threads)

    jobs = []
    for file in repeat_files:
        job = ex.submit(get_longest_repeats_v4, file, fixed_extend_base_threshold,
                        max_single_repeat_len, debug)
        jobs.append(job)
    ex.shutdown(wait=True)

    total_out = tmp_output_dir + '/total_blast_' + str(ref_index) + '.out'
    os.system('rm -f ' + total_out)
    longest_repeats = {}
    for job in as_completed(jobs):
        cur_longest_repeats, cur_keep_longest_query, cur_blastn2Results_path = job.result()
        longest_repeats.update(cur_longest_repeats)
        if debug == 1:
            os.system('cat ' + cur_blastn2Results_path + ' >> ' + total_out)
        else:
            os.system('rm -f ' + cur_blastn2Results_path)

    # Currently, it is still in coordinate form. Convert it into a nucleotide sequence and store it.
    ref_contigNames, ref_contigs = read_fasta(reference)
    for frag_name in longest_repeats.keys():
        parts = frag_name.split(':')
        chr_name = parts[0]
        pos_info = parts[1].split('-')
        chr_start = int(pos_info[0])
        chr_end = int(pos_info[1])
        frag_seq = ref_contigs[chr_name][chr_start: chr_end]
        longest_repeats[frag_name] = frag_seq
    store_fasta(longest_repeats, longest_repeats_path)
    if debug != 1:
        os.system('rm -rf ' + tmp_blast_dir)
    return longest_repeats_path

def filter_out_by_category(TE_out, tmp_output_dir, category):
    tmp_out= tmp_output_dir + '/tmp.out'
    os.system('cp ' + TE_out + ' ' + tmp_out)
    if category == 'Total':
        return tmp_out
    else:
        print('Warning: you are using the parameter "--cat ' + str(category) + '". The program will only calculate transposons that include the label "' + str(category) + '". Please make sure to check the labels in "--standard_lib" and "--test_lib" to avoid miscalculations. For example, if you want to assess the performance of TIR transposons using "--cat DNA", please ensure that all types of TIR transposon labels follow the pattern "xxx#DNA/xxx", and that no non-TIR transposon labels contain "DNA".')
        lines = []
        with open(tmp_out, 'r') as f_r:
            for line in f_r:
                query_name = line.split('\t')[0]
                parts = query_name.split('#')
                type = parts[1]
                if category in type:
                    lines.append(line)
        filter_tmp_out = tmp_output_dir + '/tmp.filter.out'
        with open(filter_tmp_out, 'w') as f_save:
            for line in lines:
                f_save.write(line)
        return filter_tmp_out

def run_command(command):
    subprocess.run(command, check=True, shell=True)

def identify_terminals(split_file, output_dir, tool_dir):
    #ltr_log = split_file + '.ltr.log'
    tir_log = split_file + '.itr.log'
    #ltrsearch_command = 'cd ' + output_dir + ' && ' + tool_dir + '/ltrsearch -l 50 ' + split_file + ' > ' + ltr_log
    itrsearch_command = 'cd ' + output_dir + ' && ' + tool_dir + '/itrsearch -i 0.7 -l 7 ' + split_file+ ' > ' + tir_log
    #run_command(ltrsearch_command)
    run_command(itrsearch_command)
    ltr_file = split_file + '.ltr'
    tir_file = split_file + '.itr'

    # ltr_identity_dict = {}
    # sequence_id = None
    # with open(ltr_log, 'r') as f_r:
    #     for line in f_r:
    #         if line.startswith('load sequence'):
    #             sequence_id = line.split('\t')[0].split(' ')[3]
    #         elif line.__contains__('Match Percentage') and sequence_id is not None:
    #             identity = float(line.split(':')[1].strip().strip('%')) / 100
    #             ltr_identity_dict[sequence_id] = identity
    #             sequence_id = None

    tir_identity_dict = {}
    sequence_id = None
    with open(tir_log, 'r') as f_r:
        for line in f_r:
            if line.startswith('load sequence'):
                sequence_id = line.split('\t')[0].split(' ')[3]
            elif line.__contains__('Identity percentage') and sequence_id is not None:
                identity = float(line.split(':')[1].strip())
                tir_identity_dict[sequence_id] = identity
                sequence_id = None

    # Read ltr and itr files to get the start and end positions of ltr and itr.
    #ltr_names, ltr_contigs = read_fasta_v1(ltr_file)
    # LTR_info = {}
    # for i, ltr_name in enumerate(ltr_names):
    #     parts = ltr_name.split('\t')
    #     orig_name = parts[0].split(' ')[0]
    #     terminal_info = parts[-1]
    #     LTR_info_parts = terminal_info.split('LTR')[1].split(' ')[0].replace('(', '').replace(')', '').split('..')
    #     LTR_left_pos_parts = LTR_info_parts[0].split(',')
    #     LTR_right_pos_parts = LTR_info_parts[1].split(',')
    #     lLTR_start = int(LTR_left_pos_parts[0])
    #     lLTR_end = int(LTR_left_pos_parts[1])
    #     rLTR_start = int(LTR_right_pos_parts[1])
    #     rLTR_end = int(LTR_right_pos_parts[0])
    #     LTR_info[orig_name] = (lLTR_start, lLTR_end, rLTR_start, rLTR_end, ltr_identity_dict[orig_name])

    tir_names, tir_contigs = read_fasta_v1(tir_file)
    TIR_info = {}
    for i, tir_name in enumerate(tir_names):
        parts = tir_name.split('\t')
        orig_name = parts[0].split(' ')[0]
        terminal_info = parts[-1]
        TIR_info_parts = terminal_info.split('ITR')[1].split(' ')[0].replace('(', '').replace(')', '').split('..')
        TIR_left_pos_parts = TIR_info_parts[0].split(',')
        TIR_right_pos_parts = TIR_info_parts[1].split(',')
        lTIR_start = int(TIR_left_pos_parts[0])
        lTIR_end = int(TIR_left_pos_parts[1])
        rTIR_start = int(TIR_right_pos_parts[1])
        rTIR_end = int(TIR_right_pos_parts[0])
        TIR_info[orig_name] = (lTIR_start, lTIR_end, rTIR_start, rTIR_end, tir_identity_dict[orig_name])
    return TIR_info

def search_TSD_regular(motif, sequence):
    motif_length = len(motif)
    pattern = ''

    # Build a regular expression pattern based on motif length.
    if motif_length >= 8:
        for i in range(motif_length):
            pattern += f"{motif[:i]}[ACGT]{motif[i + 1:]}" if i < motif_length - 1 else motif[:i] + "[ACGT]"
            if i < motif_length - 1:
                pattern += "|"
    else:
        pattern = motif

    matches = re.finditer(pattern, sequence)

    found = False
    pos = None
    for match in matches:
        #print(f"Found motif at position {match.start()}: {match.group()}")
        found = True
        pos = match.start()
        break
    return found, pos

def search_confident_tsd(orig_seq, raw_tir_start, raw_tir_end, tsd_search_distance):
    # Change all coordinates to start from 0.
    raw_tir_start -= 1
    raw_tir_end -= 1

    orig_seq_len = len(orig_seq)
    # 1. First, take 2 * tsd_search_distance sequences near the start and end positions
    left_start = raw_tir_start - tsd_search_distance
    if left_start < 0:
        left_start = 0
    # We don’t search inwards here because we consider Repbase boundaries to be correct.
    # If we consider the boundaries to be incorrect, many abnormal TSDs may meet the requirements.
    # For simplicity, we assume that Repbase boundaries are correct.
    left_end = raw_tir_start
    left_round_seq = orig_seq[left_start: left_end]
    # Obtain the position offset of left_round_seq relative to the entire sequence to correct the subsequent TSD boundary positions.
    left_offset = left_start
    right_start = raw_tir_end + 1
    if right_start < 0:
        right_start = 0
    right_end = raw_tir_end + tsd_search_distance + 1
    right_round_seq = orig_seq[right_start: right_end]
    # Obtain the position offset of right_round_seq relative to the entire sequence to correct the subsequent TSD boundary positions.
    right_offset = right_start

    # 2. Split the left sequence into k-mers from large to small, then search for the right sequence with k-mers.
    # If found, record as a candidate TSD, and finally select the closest one to the original boundary as the TSD.
    TIR_TSDs = [15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2]
    # Record the position nearest to the original boundary.
    is_found = False
    tsd_set = []
    for k_num in TIR_TSDs:
        for i in range(len(left_round_seq) - k_num, -1, -1):
            left_kmer = left_round_seq[i: i + k_num]
            left_pos = left_offset + i + k_num
            if left_pos < 0 or left_pos > orig_seq_len-1:
                continue
            found_tsd, right_pos = search_TSD_regular(left_kmer, right_round_seq)
            if found_tsd and not left_kmer.__contains__('N'):
                right_pos = right_offset + right_pos - 1
                is_found = True
                # Calculate the distance from the original boundary.
                left_distance = abs(left_pos - raw_tir_start)
                right_distance = abs(right_pos - raw_tir_end)
                distance = left_distance + right_distance
                TSD_seq = left_kmer
                TSD_len = len(TSD_seq)
                tsd_set.append((distance, TSD_len, TSD_seq))
    tsd_set = sorted(tsd_set, key=lambda x: (x[0], -x[1]))

    if not is_found:
        TSD_seq = 'NA'
        TSD_len = 'NA'
        min_distance = -1
    else:
        TSD_seq = tsd_set[0][2]
        TSD_len = tsd_set[0][1]
        min_distance = tsd_set[0][0]
    return TSD_seq, TSD_len, min_distance

def find_nearest_polyA_v1(sequence, search_range=30, min_length=6):
    max_length = 0
    max_start = -1
    max_end = -1

    # 在序列开头处查找多聚A结构
    current_length = 0
    start = 0
    for i, base in enumerate(sequence):
        if i >= search_range:
            break
        if base == 'A':
            current_length += 1
            if current_length == 1:
                start = i
        else:
            if current_length >= min_length and current_length > max_length:
                max_length = current_length
                max_start = start
                max_end = i
            current_length = 0

    # 更新最长多聚A结构的起始和结束位置
    if current_length >= min_length and current_length > max_length:
        max_start = start
        max_end = len(sequence)
    seq1 = sequence[max_start:max_end]

    # 在序列结尾处查找多聚A结构
    current_length = 0
    start = 0
    for i in range(len(sequence) - 1, -1, -1):
        if len(sequence) - i >= search_range:
            break
        if sequence[i] == 'A':
            current_length += 1
            if current_length == 1:
                start = i
        else:
            if current_length >= min_length and current_length > max_length:
                max_length = current_length
                max_start = start
                max_end = i + 1
            current_length = 0
    seq2 = sequence[max_end: max_start+1]

    seq = seq1 if len(seq1) > len(seq2) else seq2
    return max_start, max_end, seq

def find_nearest_polyT_v1(sequence, search_range=30, min_length=6):
    max_length = 0
    max_start = -1
    max_end = -1

    # 在序列开头处查找多聚T结构
    current_length = 0
    start = 0
    for i, base in enumerate(sequence):
        if i >= search_range:
            break
        if base == 'T':
            current_length += 1
            if current_length == 1:
                start = i
        else:
            if current_length >= min_length and current_length > max_length:
                max_length = current_length
                max_start = start
                max_end = i
            current_length = 0

    # 更新最长多聚A结构的起始和结束位置
    if current_length >= min_length and current_length > max_length:
        max_start = start
        max_end = len(sequence)
    seq1 = sequence[max_start:max_end]

    # 在序列结尾处查找多聚A结构
    current_length = 0
    start = 0
    for i in range(len(sequence) - 1, -1, -1):
        if len(sequence) - i >= search_range:
            break
        if sequence[i] == 'T':
            current_length += 1
            if current_length == 1:
                start = i
        else:
            if current_length >= min_length and current_length > max_length:
                max_length = current_length
                max_start = start
                max_end = i + 1
            current_length = 0
    seq2 = sequence[max_end: max_start+1]

    seq = seq1 if len(seq1) > len(seq2) else seq2
    return max_start, max_end, seq

def run_EAHelitron_v1(temp_dir, all_candidate_helitron_path, EAHelitron, partition_index):
    # 输入是Helitron序列，输出是hairpin loop序列
    all_candidate_helitron_contigs = {}
    contigNames, contigs = read_fasta(all_candidate_helitron_path)
    for query_name in contigNames:
        seq = contigs[query_name]
        all_candidate_helitron_contigs[query_name] = seq
    store_fasta(all_candidate_helitron_contigs, all_candidate_helitron_path)
    EAHelitron_command = 'cd ' + temp_dir + ' && ' + 'perl ' + EAHelitron + '/EAHelitron -o ' + str(partition_index) + ' -u 20000 -T "TC" -r 3 ' + all_candidate_helitron_path
    os.system(EAHelitron_command + '> /dev/null 2>&1')

    all_EAHelitron_res = temp_dir + '/' + str(partition_index) + '.3.txt'
    all_copies_out_names, all_copies_out_contigs = read_fasta_v1(all_EAHelitron_res)
    # search for hairpin loop sequence
    copies_hairpin_loops = {}
    for cur_name in all_copies_out_contigs.keys():
        name_parts = cur_name.split(' ')
        raw_name = name_parts[1]
        parts = raw_name.split(':')
        query_name = ':'.join(parts[:-1])
        forward_loop = name_parts[3]
        mid_loop = name_parts[4]
        reverse_loop = getReverseSequence(forward_loop)
        hairpin_loop_seq = forward_loop + mid_loop + reverse_loop
        r_hairpin_loop_seq = getReverseSequence(hairpin_loop_seq)
        cur_tail_seq = all_copies_out_contigs[cur_name]
        if cur_tail_seq.__contains__(hairpin_loop_seq):
            final_hairpin_loop_seq = hairpin_loop_seq
        elif cur_tail_seq.__contains__(r_hairpin_loop_seq):
            final_hairpin_loop_seq = r_hairpin_loop_seq
        else:
            final_hairpin_loop_seq = 'None'
        copies_hairpin_loops[query_name] = final_hairpin_loop_seq
    return copies_hairpin_loops

def get_structure_info(input_file, query_name, query_copies, flank_query_copies, cluster_dir, search_struct, tools_dir):
    if str(query_name).__contains__('Helitron'):
        flanking_len = 5
    else:
        flanking_len = 50

    annotations = {}
    if search_struct:
        (file_dir, filename) = os.path.split(input_file)
        full_length_copies_file = input_file + '.copies.fa'
        store_fasta(query_copies, full_length_copies_file)
        flank_full_length_copies_file = input_file + '.flank.copies.fa'
        store_fasta(flank_query_copies, flank_full_length_copies_file)
        if not str(filename).__contains__('Helitron'):
            if str(filename).__contains__('TIR'):
                # get LTR/TIR length and identity for TIR transposons
                TIR_info = identify_terminals(full_length_copies_file, cluster_dir, tools_dir)
                for copy_name in query_copies:
                    TIR_str = 'tir='
                    if TIR_info.__contains__(copy_name):
                        lTIR_start, lTIR_end, rTIR_start, rTIR_end, identity = TIR_info[copy_name]
                        TIR_str += str(lTIR_start) + '-' + str(lTIR_end) + ',' + str(rTIR_start) + '-' + str(
                            rTIR_end) + ';tir_identity=' + str(identity)
                    else:
                        TIR_str += 'NA'
                    update_name = TIR_str

                    flank_seq = flank_query_copies[copy_name]
                    tir_start = flanking_len + 1
                    tir_end = len(flank_seq) - flanking_len
                    tsd_search_distance = flanking_len
                    cur_tsd, cur_tsd_len, min_distance = search_confident_tsd(flank_seq, tir_start, tir_end,
                                                                              tsd_search_distance)
                    update_name += ';tsd=' + cur_tsd + ';tsd_len=' + str(cur_tsd_len)
                    if not annotations.__contains__(query_name):
                        annotations[query_name] = []
                    annotation_list = annotations[query_name]
                    annotation_list.append((copy_name, update_name))
            elif str(filename).__contains__('Non_LTR'):
                # get TSD and polyA/T head or tail for non-ltr transposons
                for copy_name in query_copies:
                    sequence = query_copies[copy_name]
                    max_start, max_end, polyA = find_nearest_polyA_v1(sequence, min_length=6)
                    max_start, max_end, polyT = find_nearest_polyT_v1(sequence, min_length=6)
                    polyA_T = polyA if len(polyA) > len(polyT) else polyT
                    update_name = 'polya_t=' + polyA_T

                    flank_seq = flank_query_copies[copy_name]
                    tir_start = flanking_len + 1
                    tir_end = len(flank_seq) - flanking_len
                    tsd_search_distance = flanking_len
                    cur_tsd, cur_tsd_len, min_distance = search_confident_tsd(flank_seq, tir_start, tir_end,
                                                                              tsd_search_distance)
                    update_name += ';tsd=' + cur_tsd + ';tsd_len=' + str(cur_tsd_len)
                    if not annotations.__contains__(query_name):
                        annotations[query_name] = []
                    annotation_list = annotations[query_name]
                    annotation_list.append((copy_name, update_name))
        else:
            # search for hairpin loop
            EAHelitron = tools_dir + '/../bin/EAHelitron-master'
            copies_hairpin_loops = run_EAHelitron_v1(cluster_dir, flank_full_length_copies_file, EAHelitron, query_name)
            for copy_name in query_copies:
                if copies_hairpin_loops.__contains__(copy_name):
                    hairpin_loop = copies_hairpin_loops[copy_name]
                else:
                    hairpin_loop = 'NA'
                update_name = 'hairpin_loop=' + hairpin_loop
                if not annotations.__contains__(query_name):
                    annotations[query_name] = []
                annotation_list = annotations[query_name]
                annotation_list.append((copy_name, update_name))
    else:
        if not annotations.__contains__(query_name):
            annotations[query_name] = []
        annotation_list = annotations[query_name]
        for copy_name in query_copies.keys():
            annotation_list.append((copy_name, ''))
    return annotations

def get_full_length_copies_from_blastn(TE_lib, reference, blastn_out, tmp_output_dir, threads, divergence_threshold,
                                    full_length_threshold, search_struct, tools_dir):
    ref_names, ref_contigs = read_fasta(reference)

    query_names, query_contigs = read_fasta(TE_lib)
    new_query_contigs = {}
    for name in query_names:
        new_query_contigs[name.split('#')[0]] = query_contigs[name]
    query_contigs = new_query_contigs

    query_records = {}
    with open(blastn_out, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            info_parts = line.split('\t')
            query_name = info_parts[0].split('#')[0]
            subject_name = info_parts[1]
            q_start = int(info_parts[6])
            q_end = int(info_parts[7])
            s_start = int(info_parts[8])
            s_end = int(info_parts[9])
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))

    full_length_copies = {}
    flank_full_length_copies = {}
    copies_direct = {}
    for idx, query_name in enumerate(query_records.keys()):
        subject_dict = query_records[query_name]
        if query_name not in query_contigs:
            continue
        query_len = len(query_contigs[query_name])
        skip_gap = query_len * full_length_threshold

        if str(query_name).__contains__('Helitron'):
            flanking_len = 5
        else:
            flanking_len = 50

        # if there are more than one longest query overlap with the final longest query over 90%,
        # then it probably the true TE
        longest_queries = []
        for subject_name in subject_dict.keys():
            subject_pos = subject_dict[subject_name]

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
                        cur_subject_start = frag[2]
                        cur_query_end = frag[1]
                        prev_subject_end = exist_frag[3]
                        prev_query_end = exist_frag[1]
                        if (cur_subject_start - prev_subject_end < skip_gap and cur_query_end > prev_query_end):
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
                        cur_subject_start = frag[2]
                        cur_query_end = frag[1]
                        prev_subject_end = exist_frag[3]
                        prev_query_end = exist_frag[1]
                        if (prev_subject_end - cur_subject_start < skip_gap and cur_query_end > prev_query_end):
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

                # print('subject pos size: %d' %(len(cur_cluster)))
                # record visited fragments
                visited_frag = {}
                for i in range(len(cur_cluster)):
                    # keep a longest query start from each fragment
                    prev_frag = cur_cluster[i]
                    if visited_frag.__contains__(prev_frag):
                        continue
                    prev_query_start = prev_frag[0]
                    prev_query_end = prev_frag[1]
                    prev_subject_start = prev_frag[2]
                    prev_subject_end = prev_frag[3]
                    prev_query_seq = (min(prev_query_start, prev_query_end), max(prev_query_start, prev_query_end))
                    prev_subject_seq = (
                        min(prev_subject_start, prev_subject_end), max(prev_subject_start, prev_subject_end))
                    prev_query_len = abs(prev_query_end - prev_query_start)
                    prev_subject_len = abs(prev_subject_end - prev_subject_start)
                    cur_longest_query_len = prev_query_len

                    cur_extend_num = 0
                    visited_frag[prev_frag] = 1
                    # try to extend query
                    for j in range(i + 1, len(cur_cluster)):
                        cur_frag = cur_cluster[j]
                        if visited_frag.__contains__(cur_frag):
                            continue
                        cur_query_start = cur_frag[0]
                        cur_query_end = cur_frag[1]
                        cur_subject_start = cur_frag[2]
                        cur_subject_end = cur_frag[3]
                        cur_query_seq = (min(cur_query_start, cur_query_end), max(cur_query_start, cur_query_end))
                        cur_subject_seq = (min(cur_subject_start, cur_subject_end), max(cur_subject_start, cur_subject_end))

                        # could extend
                        # extend right
                        if cur_query_end > prev_query_end:
                            # judge subject direction
                            if prev_subject_start < prev_subject_end and cur_subject_start < cur_subject_end:
                                # +
                                if cur_subject_end > prev_subject_end:
                                    # forward extend
                                    if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                            and cur_subject_start - prev_subject_end < skip_gap:  # \
                                        # and not is_same_query and not is_same_subject:
                                        # update the longest path
                                        prev_query_start = prev_query_start
                                        prev_query_end = cur_query_end
                                        prev_subject_start = prev_subject_start if prev_subject_start < cur_subject_start else cur_subject_start
                                        prev_subject_end = cur_subject_end
                                        cur_longest_query_len = prev_query_end - prev_query_start
                                        cur_extend_num += 1
                                        visited_frag[cur_frag] = 1
                                    elif cur_query_start - prev_query_end >= skip_gap:
                                        break
                            elif prev_subject_start > prev_subject_end and cur_subject_start > cur_subject_end:
                                # reverse
                                if cur_subject_end < prev_subject_end:
                                    # reverse extend
                                    if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                            and prev_subject_end - cur_subject_start < skip_gap:  # \
                                        # and not is_same_query and not is_same_subject:
                                        # update the longest path
                                        prev_query_start = prev_query_start
                                        prev_query_end = cur_query_end
                                        prev_subject_start = prev_subject_start if prev_subject_start > cur_subject_start else cur_subject_start
                                        prev_subject_end = cur_subject_end
                                        cur_longest_query_len = prev_query_end - prev_query_start
                                        cur_extend_num += 1
                                        visited_frag[cur_frag] = 1
                                    elif cur_query_start - prev_query_end >= skip_gap:
                                        break
                    # keep this longest query
                    if cur_longest_query_len != -1:
                        longest_queries.append(
                            (prev_query_start, prev_query_end, cur_longest_query_len, prev_subject_start,
                             prev_subject_end, abs(prev_subject_end - prev_subject_start), subject_name,
                             cur_extend_num))

        # To determine whether each copy has a coverage exceeding the full_length_threshold with respect
        # to the consensus sequence, retaining full-length copies.
        query_copies = {}
        flank_query_copies = {}
        orig_query_len = len(query_contigs[query_name])
        # query_copies[query_name] = query_contigs[query_name]
        for repeat in longest_queries:
            if repeat[2] < full_length_threshold * query_len:
                continue
            # Subject
            subject_name = repeat[6]
            subject_chr_start = 0

            if repeat[3] > repeat[4]:
                direct = '-'
                old_subject_start_pos = repeat[4] - 1
                old_subject_end_pos = repeat[3]
            else:
                direct = '+'
                old_subject_start_pos = repeat[3] - 1
                old_subject_end_pos = repeat[4]
            subject_start_pos = subject_chr_start + old_subject_start_pos
            subject_end_pos = subject_chr_start + old_subject_end_pos

            subject_pos = subject_name + ':' + str(subject_start_pos) + '-' + str(subject_end_pos)
            subject_seq = ref_contigs[subject_name][subject_start_pos: subject_end_pos]

            flank_subject_seq = ref_contigs[subject_name][
                                subject_start_pos - flanking_len: subject_end_pos + flanking_len]
            copies_direct[subject_pos] = direct
            cur_query_len = repeat[2]
            coverage = float(cur_query_len) / orig_query_len
            if coverage >= full_length_threshold:
                query_copies[subject_pos] = subject_seq
                flank_query_copies[subject_pos] = flank_subject_seq
        full_length_copies[query_name] = query_copies
        flank_full_length_copies[query_name] = flank_query_copies

    # The candidate full-length copies and the consensus are then clustered using cd-hit-est,
    # retaining copies that belong to the same cluster as the consensus.
    split_files = []
    cluster_dir = tmp_output_dir + '/cluster'
    os.system('rm -rf ' + cluster_dir)
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)

    all_query_copies = {}
    for query_name in full_length_copies.keys():
        query_copies = full_length_copies[query_name]
        flank_query_copies = flank_full_length_copies[query_name]
        all_query_copies.update(query_copies)
        fc_path = cluster_dir + '/' + query_name + '.fa'
        store_fasta(query_copies, fc_path)
        split_files.append((fc_path, query_name, query_copies, flank_query_copies))

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for ref_index, cur_file in enumerate(split_files):
        input_file = cur_file[0]
        query_name = cur_file[1]
        query_copies = cur_file[2]
        flank_query_copies = cur_file[3]
        job = ex.submit(get_structure_info, input_file, query_name, query_copies,
                        flank_query_copies, cluster_dir, search_struct, tools_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    full_length_annotations = {}
    for job in as_completed(jobs):
        annotations = job.result()
        full_length_annotations.update(annotations)
    return full_length_annotations, copies_direct

def generate_full_length_out(BlastnOut, full_length_out, TE_lib, reference, tmp_output_dir, tools_dir, full_length_threshold, category):
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)
    filter_tmp_out = filter_out_by_category(BlastnOut, tmp_output_dir, category)

    threads = 1
    divergence_threshold = 20
    search_struct = False
    full_length_annotations, copies_direct = get_full_length_copies_from_blastn(TE_lib, reference, filter_tmp_out,
                                                                             tmp_output_dir, threads,
                                                                             divergence_threshold,
                                                                             full_length_threshold,
                                                                             search_struct, tools_dir)

    lines = set()
    for query_name in full_length_annotations.keys():
        query_name = str(query_name)
        for copy_annotation in full_length_annotations[query_name]:
            chr_pos = copy_annotation[0]
            annotation = copy_annotation[1]
            parts = chr_pos.split(':')
            chr_name = parts[0]
            chr_pos_parts = parts[1].split('-')
            chr_start = int(chr_pos_parts[0]) + 1
            chr_end = int(chr_pos_parts[1])
            new_line = (query_name, chr_name, chr_start, chr_end)
            lines.add(new_line)

    return lines

def mask_genome_intactTE(TE_lib, genome_path, work_dir, thread, ref_index):
    masked_genome_path = genome_path + '.masked'
    if file_exist(TE_lib):
        tmp_blast_dir = work_dir + '/mask_tmp_' + str(ref_index)
        lib_out = work_dir + '/prev_TE_'+str(ref_index)+'.out'
        multi_process_align(TE_lib, genome_path, lib_out, tmp_blast_dir, thread, is_removed_dir=True)
        full_length_out = work_dir + '/mask_full_length'+str(ref_index)+'.out'
        coverage_threshold = 0.8
        category = 'Total'
        sorted_lines = generate_full_length_out(lib_out, full_length_out, TE_lib, genome_path, tmp_blast_dir, '',
                                 coverage_threshold, category)

        ref_names, ref_contigs = read_fasta(genome_path)
        # mask genome
        for query_name, chr_name, chr_start, chr_end in sorted_lines:
            start = chr_start - 1
            end = chr_end
            ref_seq = ref_contigs[chr_name]
            mask_seq = ref_seq[:start] + 'N' * (end - start) + ref_seq[end:]
            ref_contigs[chr_name] = mask_seq

        store_fasta(ref_contigs, masked_genome_path)
    else:
        os.system('cp ' + genome_path + ' ' + masked_genome_path)


def get_TSD(all_copies, flanking_len):
    # tsd_info = {query_name: {copy1: tsd+','+seq}, {copy2: tsd+','+seq}, {total_copy_num:}, {tsd_copy_num:}}
    # copy: (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq)
    # (subject_name, subject_start, subject_end, query_len, direct)
    tsd_info = {}
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        tsd_copy_num = 0
        total_copy_len = 0
        for copy in copies:
            orig_seq = copy[4]
            total_copy_len += len(orig_seq)

            tir_start = flanking_len + 1
            tir_end = len(orig_seq)-flanking_len
            left_tsd_seq, right_tsd_seq = TSDsearch_v4(orig_seq, tir_start, tir_end)
            if left_tsd_seq != '':
                tsd_copy_num += 1
            copy_name = str(copy[0])+'-'+str(copy[1])+'-'+str(copy[2])+'-'+str(copy[3])
            if not tsd_info.__contains__(query_name):
                tsd_info[query_name] = {}
            info = tsd_info[query_name]
            info[copy_name] = left_tsd_seq+','+right_tsd_seq+','+orig_seq
        if not tsd_info.__contains__(query_name):
            tsd_info[query_name] = {}
        info = tsd_info[query_name]
        info['total_copy_num'] = len(copies)
        info['tsd_copy_num'] = tsd_copy_num
        info['total_copy_len'] = total_copy_len
    return tsd_info


def store_copies(tsd_info, copy_info_path):
    # tsd_info = {query_name: {copy1: tsd+','+seq}, {copy2: tsd+','+seq}, {total_copy_num:}, {tsd_copy_num:}}
    with open(copy_info_path, 'w') as f_save:
        for query_name in tsd_info.keys():
            f_save.write(query_name + '\n')
            info = tsd_info[query_name]
            for copy_name in info.keys():
                if copy_name != 'total_copy_num' and copy_name != 'tsd_copy_num' and copy_name != 'total_copy_len':
                    info_parts = info[copy_name].split(',')
                    left_tsd_seq = info_parts[0]
                    right_tsd_seq = info_parts[1]
                    copy_seq = info_parts[2]
                    f_save.write('\t' + str(copy_name) + '\tleft_tsd_seq: ' + str(left_tsd_seq) + '\tright_tsd_seq: ' + str(right_tsd_seq) + '\n')
                    f_save.write(copy_seq + '\n')
            total_copy_num = info['total_copy_num']
            tsd_copy_num = info['tsd_copy_num']
            f_save.write('\ttotal_copy_num: ' + str(total_copy_num) + ', tsd_copy_num: ' + str(tsd_copy_num) + '\n')
    f_save.close()

def store_copies_v1(copies, copy_info_path):
    # new_copies.append((ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq))
    with open(copy_info_path, 'w') as f_save:
        for query_name in copies.keys():
            f_save.write(query_name + '\n')
            copy_list = copies[query_name]
            for copy in copy_list:
                f_save.write('\t'+str(copy[0])+':'+str(copy[1])+'-'+str(copy[2])+'-'+str(copy[2]-copy[1]+1)+'\n')
                f_save.write(copy[4] + '\n')
    f_save.close()

def store_copies_seq(copies, copy_info_path):
    # new_copies.append((ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq))
    with open(copy_info_path, 'w') as f_save:
        for query_name in copies.keys():
            copy_list = copies[query_name]
            for i, copy in enumerate(copy_list):
                new_query_name = query_name + '-C_' + str(i)
                f_save.write('>'+new_query_name + '\n' + copy[4] + '\n')
    f_save.close()

def split_fasta(fasta_file, subfile_size, output_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Read the input fasta file and iterate over the records
    with open(fasta_file) as handle:
        current_subfile_size = 0
        current_subfile_index = 0
        current_subfile_name = os.path.join(output_dir, "subfile_"+str(current_subfile_index)+".fasta")
        current_subfile = open(current_subfile_name, "w")

        for line in handle:
            if line.startswith(">") and current_subfile_size >= subfile_size:
                # Close the current subfile and open a new one
                current_subfile.close()
                current_subfile_index += 1
                current_subfile_name = os.path.join(output_dir, "subfile_"+str(current_subfile_index)+".fasta")
                current_subfile = open(current_subfile_name, "w")
                current_subfile_size = 0

            # Write the current line to the current subfile
            current_subfile.write(line)
            current_subfile_size += len(line)

        # Close the last subfile
        current_subfile.close()

def search_confident_tir_batch_v1(split_file, flanking_len, tir_tsd_dir, TRsearch_dir, partition_index, plant):
    cur_contignames, cur_contigs = read_fasta(split_file)

    all_copies_itr_contigs = {}
    all_candidate_TIRs_path = tir_tsd_dir + '/' + str(partition_index) + '.fa'
    short_candidate_TIRs_path = tir_tsd_dir + '/' + str(partition_index) + '_s.fa'
    for query_name in cur_contignames:
        seq = cur_contigs[query_name]

        if seq.__contains__('NNNNNNNNNN'):
            continue

        tir_start = flanking_len + 1
        tir_end = len(seq) - flanking_len
        # Find all possible TSD sequences, calculate the distance between each sequence boundary
        # and the original boundary, and store it in the header
        tsd_search_distance = flanking_len
        cur_itr_contigs = search_confident_tir_v4(seq, tir_start, tir_end, tsd_search_distance, query_name, plant)
        all_copies_itr_contigs.update(cur_itr_contigs)

    # Save the sequence of short TIR and submit it to itrsearch to determine the length of TIR
    short_itr_contigs = get_short_tir_contigs(all_copies_itr_contigs, plant)
    # Since we're only interested in determining the presence of TIR structure, we don't need to consider the entire sequence.
    # Instead, we can simply select 40 base pairs from each end for analysis. This approach significantly speeds up the process.
    new_short_itr_contigs = {}
    for name in short_itr_contigs.keys():
        seq = short_itr_contigs[name]
        new_seq = seq[:40] + seq[-40:]
        new_short_itr_contigs[name] = new_seq
    store_fasta(new_short_itr_contigs, short_candidate_TIRs_path)
    short_copies_out, short_copies_log = run_itrsearch(TRsearch_dir, short_candidate_TIRs_path, tir_tsd_dir)
    raw_short_copies_out_name, raw_short_copies_out_contigs = read_fasta_v1(short_copies_out)

    # The remaining sequences are handed over to itrsearch to search for TIR structures
    for name in short_itr_contigs.keys():
        del all_copies_itr_contigs[name]
    new_all_copies_itr_contigs = {}
    for name in all_copies_itr_contigs.keys():
        seq = all_copies_itr_contigs[name]
        new_seq = seq[:40] + seq[-40:]
        new_all_copies_itr_contigs[name] = new_seq
    store_fasta(new_all_copies_itr_contigs, all_candidate_TIRs_path)
    all_copies_out, all_copies_log = run_itrsearch(TRsearch_dir, all_candidate_TIRs_path, tir_tsd_dir)

    # Filter out sequences with excessive differences in terminal TIR length
    # filter_large_gap_tirs(all_copies_out, all_copies_out)

    # Record the TIR length corresponding to each query
    TIR_len_dict = {}
    raw_all_copies_out_name, raw_all_copies_out_contigs = read_fasta_v1(all_copies_out)
    all_copies_out_name, all_copies_out_contigs = read_fasta(all_copies_out)
    for name in raw_all_copies_out_name:
        query_name = name.split(' ')[0]
        tir_len = int(name.split('Length itr=')[1])
        TIR_len_dict[query_name] = tir_len
    for name in raw_short_copies_out_name:
        query_name = name.split(' ')[0]
        tir_len = int(name.split('Length itr=')[1])
        TIR_len_dict[query_name] = tir_len

    # # Analyze the itrsearch log file and extract the sequence names of alignment offsets
    # fake_tirs = get_fake_tirs(all_copies_log)
    # for name in all_copies_out_name:
    #     seq = all_copies_itr_contigs[name]
    #     if name in fake_tirs:
    #         del all_copies_out_contigs[name]
    #     else:
    #         all_copies_out_contigs[name] = seq
    # all_copies_out_contigs.update(short_itr_contigs)

    for name in all_copies_out_name:
        seq = all_copies_itr_contigs[name]
        all_copies_out_contigs[name] = seq
    all_copies_out_contigs.update(short_itr_contigs)

    # To reduce the amount of subsequent computation, we will perform grouping here,
    # selecting the sequences with the longest TIR and TSD for each group.
    # Group by query_name, and only select one sequence for each group,
    # which is the one with the best combined copy number and TSD.
    # Group all_copies_out_contigs by query_name.
    # group_copies_contigs -> {query_name: {name: seq}}
    group_copies_contigs = {}
    for cur_name in all_copies_out_contigs.keys():
        query_name = cur_name.split('-C_')[0]
        if not group_copies_contigs.__contains__(query_name):
            group_copies_contigs[query_name] = {}
        cur_copies_out_contigs = group_copies_contigs[query_name]
        cur_copies_out_contigs[cur_name] = all_copies_out_contigs[cur_name]

    filter_dup_itr_contigs = {}
    for query_name in group_copies_contigs.keys():
        cur_copies_out_contigs = group_copies_contigs[query_name]
        # Select the one with the smallest distance
        cur_contigs = filter_dup_itr_v3(cur_copies_out_contigs, TIR_len_dict)
        filter_dup_itr_contigs.update(cur_contigs)
    return filter_dup_itr_contigs

def multi_process_tsd_v1(longest_repeats_flanked_path, tir_tsd_path, tir_tsd_dir, flanking_len, threads, TRsearch_dir, plant):
    os.system('rm -rf '+tir_tsd_dir)
    if not os.path.exists(tir_tsd_dir):
        os.makedirs(tir_tsd_dir)

    # After partitioning the files, perform parallel computations using multiple processes.
    fasta_file = longest_repeats_flanked_path
    subfile_size = 10000 #10K
    output_dir = tir_tsd_dir
    split_fasta(fasta_file, subfile_size, output_dir)
    split_files = []
    for split_file_name in os.listdir(output_dir):
        split_file = output_dir + '/' + split_file_name
        split_files.append(split_file)

    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    for partition_index, split_file in enumerate(split_files):
        obj = ex.submit(search_confident_tir_batch_v1, split_file, flanking_len, tir_tsd_dir, TRsearch_dir, partition_index, plant)
        objs.append(obj)
        job_id += 1
    ex.shutdown(wait=True)
    candidate_TIRs = {}
    for obj in as_completed(objs):
        cur_candidate_TIRs = obj.result()
        candidate_TIRs.update(cur_candidate_TIRs)
    store_fasta(candidate_TIRs, tir_tsd_path)
    os.system('rm -rf ' + tir_tsd_dir)

def multi_process_tsd(longest_repeats_flanked_path, tir_tsd_path, tir_tsd_dir, flanking_len, threads, TRsearch_dir, plant):
    os.system('rm -rf '+tir_tsd_dir)
    if not os.path.exists(tir_tsd_dir):
        os.makedirs(tir_tsd_dir)

    # After partitioning the files, perform parallel computations using multiple processes.
    fasta_file = longest_repeats_flanked_path
    subfile_size = 10000 #10K
    output_dir = tir_tsd_dir
    split_fasta(fasta_file, subfile_size, output_dir)
    split_files = []
    for split_file_name in os.listdir(output_dir):
        split_file = output_dir + '/' + split_file_name
        split_files.append(split_file)

    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    for partition_index, split_file in enumerate(split_files):
        obj = ex.submit(search_confident_tir_batch, split_file, flanking_len, tir_tsd_dir, TRsearch_dir, partition_index, plant)
        objs.append(obj)
        job_id += 1
    ex.shutdown(wait=True)
    candidate_TIRs = {}
    for obj in as_completed(objs):
        cur_candidate_TIRs = obj.result()
        candidate_TIRs.update(cur_candidate_TIRs)
    store_fasta(candidate_TIRs, tir_tsd_path)
    os.system('rm -rf ' + tir_tsd_dir)


def multi_process_ltr_tsd(raw_candidate_ltrs, ltr_tsd_path, cut_ltr_tsd_path, ltr_tsd_dir, flanking_len, threads, TRsearch_dir, plant):
    os.system('rm -rf '+ltr_tsd_dir)
    if not os.path.exists(ltr_tsd_dir):
        os.makedirs(ltr_tsd_dir)

    # (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq)
    segments_cluster = divided_array(list(raw_candidate_ltrs.items()), threads)

    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    for partition_index, cur_segments in enumerate(segments_cluster):
        obj = ex.submit(search_tsd_batch, cur_segments, flanking_len)
        objs.append(obj)
        job_id += 1
    ex.shutdown(wait=True)
    all_copies_ltr_contigs = {}
    for obj in as_completed(objs):
        cur_copies_ltr_contigs = obj.result()
        all_copies_ltr_contigs.update(cur_copies_ltr_contigs)

    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    cur_ltr_contigs = {}
    partition_index = 0
    for index, query_name in enumerate(all_copies_ltr_contigs.keys()):
        if index % 50 == 0:
            obj = ex.submit(search_tsd_ltr_batch, cur_ltr_contigs, ltr_tsd_dir, partition_index, TRsearch_dir)
            objs.append(obj)
            job_id += 1

            partition_index += 1
            cur_ltr_contigs = {}
        else:
            cur_ltr_contigs[query_name] = all_copies_ltr_contigs[query_name]

    if len(cur_ltr_contigs) > 0:
        obj = ex.submit(search_tsd_ltr_batch, cur_ltr_contigs, ltr_tsd_dir, partition_index, TRsearch_dir)
        objs.append(obj)
        job_id += 1
        partition_index += 1
    ex.shutdown(wait=True)
    all_copies_out_contigs = {}
    for obj in as_completed(objs):
        cur_copies_out_contigs = obj.result()
        all_copies_out_contigs.update(cur_copies_out_contigs)

    candidate_LTRs, candidate_cut_LTRs = search_candidate_ltr(all_copies_out_contigs)
    store_fasta(candidate_LTRs, ltr_tsd_path)
    store_fasta(candidate_cut_LTRs, cut_ltr_tsd_path)


def flanking_copies(all_copies, query_path, reference, flanking_len, copy_num=10, query_coverage=0.99):
    new_all_copies = {}
    query_names, query_contigs = read_fasta(query_path)
    ref_names, ref_contigs = read_fasta(reference)
    for query_name in all_copies.keys():
        query_seq = query_contigs[query_name]
        copies = all_copies[query_name]
        new_copies = []
        # 取最多copy_num条
        for i, copy in enumerate(copies):
            if copy_num != -1 and i >= copy_num:
                break
            ref_name = copy[0]
            copy_ref_start = int(copy[1])
            copy_ref_end = int(copy[2])
            if copy_ref_start > copy_ref_end:
                tmp = copy_ref_start
                copy_ref_start = copy_ref_end
                copy_ref_end = tmp
            copy_len = copy_ref_end - copy_ref_start + 1
            if copy_ref_start - 1 - flanking_len < 0 or copy_ref_end + flanking_len > len(ref_contigs[ref_name]):
                continue
            ref_seq = ref_contigs[ref_name]
            orig_copy_seq = ref_seq[copy_ref_start-1: copy_ref_end]
            copy_seq = ref_seq[copy_ref_start-1-flanking_len: copy_ref_end+flanking_len]
            # If full-length copy is taken, it is necessary to judge whether the first 5bp and the last 5bp
            # of the copy are highly similar to the original sequence, and only 1-bp mismatch is allowed.
            if query_coverage != 0.99 or \
                    (allow_mismatch(query_seq[0:5], orig_copy_seq[0:5], 1)
                     and allow_mismatch(query_seq[-5:], orig_copy_seq[-5:], 1)):
                new_copies.append((ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq))
        if len(new_copies) > 0:
            new_all_copies[query_name] = new_copies
    return new_all_copies

def flanking_copies_v2(all_copies, query_path, reference, flanking_len, copy_num=10, query_coverage=0.99):
    new_all_copies = {}
    query_names, query_contigs = read_fasta(query_path)
    ref_names, ref_contigs = read_fasta(reference)
    for query_name in all_copies.keys():
        query_seq = query_contigs[query_name]
        copies = all_copies[query_name]
        new_copies = []
        # get at most copy_num copies
        for i, copy in enumerate(copies):
            if copy_num != -1 and i >= copy_num:
                break
            ref_name = copy[0]
            copy_ref_start = int(copy[1])
            copy_ref_end = int(copy[2])
            if copy_ref_start > copy_ref_end:
                tmp = copy_ref_start
                copy_ref_start = copy_ref_end
                copy_ref_end = tmp
            copy_len = copy_ref_end - copy_ref_start + 1
            if copy_ref_start - 1 - flanking_len < 0 or copy_ref_end + flanking_len > len(ref_contigs[ref_name]):
                continue
            ref_seq = ref_contigs[ref_name]
            orig_copy_seq = ref_seq[copy_ref_start-1: copy_ref_end]
            copy_seq = ref_seq[copy_ref_start-1-flanking_len: copy_ref_end+flanking_len]
            new_copies.append((ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq))
        if len(new_copies) > 0:
            new_all_copies[query_name] = new_copies
    return new_all_copies

def flanking_copies_v1(all_copies, reference, flanking_len):
    new_all_copies = {}
    ref_names, ref_contigs = read_fasta(reference)
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        if len(copies) >= 1:
            copy = copies[0]
            ref_name = copy[0]
            #avg_identity = copy[3]
            copy_ref_start = int(copy[1])
            copy_ref_end = int(copy[2])
            direct = copy[4]
            copy_len = copy_ref_end - copy_ref_start + 1
            copy_seq = ref_contigs[ref_name][copy_ref_start-1-flanking_len: copy_ref_end+flanking_len]
            if direct == '-':
                copy_seq = getReverseSequence(copy_seq)
            new_all_copies[query_name] = (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq)
    return new_all_copies

def get_query_copies(cur_segments, query_contigs, subject_path, query_coverage, subject_coverage, query_fixed_extend_base_threshold=200, subject_fixed_extend_base_threshold=200, max_copy_num=100):
    all_copies = {}

    if subject_coverage > 0:
        subject_names, subject_contigs = read_fasta(subject_path)

    for item in cur_segments:
        query_name = item[0]
        subject_dict = item[1]

        longest_queries = []
        for subject_name in subject_dict.keys():
            subject_pos = subject_dict[subject_name]

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
                        if (frag[2] - exist_frag[3] < subject_fixed_extend_base_threshold and frag[1] > exist_frag[1]):
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
                        if (exist_frag[3] - frag[2] < subject_fixed_extend_base_threshold and frag[1] > exist_frag[1]):
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

                cluster_longest_subject_start = -1
                cluster_longest_subject_end = -1
                cluster_longest_subject_len = -1

                cluster_identity = 0
                cluster_extend_num = 0

                # record visited fragments
                visited_frag = {}
                for i in range(len(cur_cluster)):
                    # keep a longest query start from each fragment
                    origin_frag = cur_cluster[i]
                    if visited_frag.__contains__(origin_frag):
                        continue
                    cur_frag_len = origin_frag[1] - origin_frag[0] + 1
                    cur_longest_query_len = cur_frag_len
                    longest_query_start = origin_frag[0]
                    longest_query_end = origin_frag[1]
                    longest_subject_start = origin_frag[2]
                    longest_subject_end = origin_frag[3]

                    cur_identity = origin_frag[4]
                    cur_extend_num = 0

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
                                    if ext_frag[0] - longest_query_end < query_fixed_extend_base_threshold and ext_frag[
                                        2] - longest_subject_end < subject_fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start

                                        cur_identity += ext_frag[4]
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= query_fixed_extend_base_threshold:
                                        break
                            elif longest_subject_start > longest_subject_end and ext_frag[2] > ext_frag[3]:
                                # reverse
                                if ext_frag[3] < longest_subject_end:
                                    # reverse extend
                                    if ext_frag[
                                        0] - longest_query_end < query_fixed_extend_base_threshold and longest_subject_end - \
                                            ext_frag[2] < subject_fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start > \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start

                                        cur_identity += ext_frag[4]
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= query_fixed_extend_base_threshold:
                                        break
                    if cur_longest_query_len > cluster_longest_query_len:
                        cluster_longest_query_start = longest_query_start
                        cluster_longest_query_end = longest_query_end
                        cluster_longest_query_len = cur_longest_query_len

                        cluster_longest_subject_start = longest_subject_start
                        cluster_longest_subject_end = longest_subject_end
                        cluster_longest_subject_len = abs(longest_subject_end - longest_subject_start) + 1

                        cluster_identity = cur_identity
                        cluster_extend_num = cur_extend_num
                # keep this longest query
                if cluster_longest_query_len != -1:
                    longest_queries.append((cluster_longest_query_start, cluster_longest_query_end,
                                            cluster_longest_query_len, cluster_longest_subject_start,
                                            cluster_longest_subject_end, cluster_longest_subject_len, subject_name,
                                            cluster_extend_num, cluster_identity))

        longest_queries.sort(key=lambda x: -x[2])
        query_len = len(query_contigs[query_name])
        # query_len = int(query_name.split('-')[1].split('_')[1])
        copies = []
        keeped_copies = set()
        for query in longest_queries:
            if len(copies) > max_copy_num:
                break
            subject_name = query[6]
            subject_start = query[3]
            subject_end = query[4]
            direct = '+'
            if subject_start > subject_end:
                tmp = subject_start
                subject_start = subject_end
                subject_end = tmp
                direct = '-'
            item = (subject_name, subject_start, subject_end)
            if subject_coverage > 0:
                subject_len = len(subject_contigs[subject_name])
                cur_subject_coverage = float(query[5])/subject_len
                if float(query[2])/query_len >= query_coverage and cur_subject_coverage >= subject_coverage and item not in keeped_copies:
                    copies.append((subject_name, subject_start, subject_end, query[2], direct))
                    keeped_copies.add(item)
            else:
                if float(query[2]) / query_len >= query_coverage and item not in keeped_copies:
                    copies.append((subject_name, subject_start, subject_end, query[2], direct))
                    keeped_copies.add(item)
        #copies.sort(key=lambda x: abs(x[3]-(x[2]-x[1]+1)))
        all_copies[query_name] = copies
    return all_copies

def get_copies_v1(blastnResults_path, query_path, subject_path, query_coverage=0.95, subject_coverage=0):
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
            if query_name == subject_name:
                continue
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end, identity))
    f_r.close()

    query_names, query_contigs = read_fasta(query_path)
    cur_segments = list(query_records.items())
    all_copies = get_query_copies(cur_segments, query_contigs, subject_path, query_coverage, subject_coverage)

    return all_copies


def get_copies(blastnResults_path, query_path, subject_path, query_coverage=0.99, subject_coverage=0, threads=48):
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
            if query_name == subject_name:
                continue
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end, identity))
    f_r.close()

    query_names, query_contigs = read_fasta(query_path)
    segments_cluster = divided_array(list(query_records.items()), threads)

    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    for partition_index, cur_segments in enumerate(segments_cluster):
        obj = ex.submit(get_query_copies, cur_segments, query_contigs, subject_path, query_coverage, subject_coverage)
        objs.append(obj)
        job_id += 1
    ex.shutdown(wait=True)

    all_copies = {}
    for obj in as_completed(objs):
        cur_copies = obj.result()
        all_copies.update(cur_copies)
    return all_copies

def generate_candidate_ltrs(all_copies, reference, flanking_len):
    # Gather adjacent copies corresponding to the same query_name together, and then generate candidate LTR sequences.
    # copy_cluster -> {ref_name: [(left_ltr_ref_start, left_ltr_ref_end, right_ltr_ref_start, right_ltr_ref_end), ..., ]}
    copy_cluster = {}
    ref_names, ref_contigs = read_fasta(reference)
    for query_name in all_copies.keys():
        copies = list(all_copies[query_name])
        # Sort the replicas by ref_name, ref_start and ref_end.
        # copy -> (subject_name, subject_start, subject_end, query[2], direct)
        copies.sort(key=lambda x: (x[0], x[1], x[2]))

        last_copy = None
        for i, copy in enumerate(copies):
            ref_name = copy[0]
            copy_ref_start = int(copy[1])
            copy_ref_end = int(copy[2])
            if copy_ref_start > copy_ref_end:
                tmp = copy_ref_start
                copy_ref_start = copy_ref_end
                copy_ref_end = tmp
            copy_len = copy[3]
            direct = copy[4]
            copy = (ref_name, copy_ref_start, copy_ref_end, copy_len, direct)

            if last_copy is not None and ref_name == last_copy[0] and direct == last_copy[4]:
                # The length range of the sequence (100-7000) defines the starting position distance (1000,15000) between two copies of the same sequence.
                if (last_copy[3] >= 100 and last_copy[3] <= 7000) and (copy_len >= 100 and copy_len <= 7000):
                    # There is no Overlap between the current copy and the previous copy, and the distance between them satisfies 1000<= x <= 15000, so it may be an LTR element.
                    if copy_ref_start > last_copy[2] and (copy_ref_start - last_copy[1] >= 1000) and (
                            copy_ref_start - last_copy[1] <= 15000):
                        if not copy_cluster.__contains__(ref_name):
                            copy_cluster[ref_name] = []
                        candidate_ltrs = copy_cluster[ref_name]
                        candidate_ltrs.append((last_copy[1], last_copy[2], copy_ref_start, copy_ref_end))
            last_copy = copy

    flanked_ltr_candidates = {}
    # flanked_ltr_candidates -> {ref_name: [(left_ltr_ref_start, left_ltr_ref_end, right_ltr_ref_start, right_ltr_ref_end, flanked_seq), ..., ]}
    # Expand all candidate_ltr by 100bp to ensure that TSD can be searched.
    for ref_name in copy_cluster.keys():
        candidate_ltrs = copy_cluster[ref_name]
        # Sort by LTR_start, ltr_end, and filter out the repeated ltrs (with the starting position of the last ltr less than 10).
        candidate_ltrs.sort(key=lambda x: (x[0], x[3]))
        last_candidate_ltr = None
        for candidate_ltr in candidate_ltrs:
            ltr_ref_start = candidate_ltr[0]
            ltr_ref_end = candidate_ltr[3]
            if last_candidate_ltr is not None:
                if abs(ltr_ref_start - last_candidate_ltr[0]) < 10 and abs(ltr_ref_end - last_candidate_ltr[3]) < 10:
                    continue
            if ltr_ref_start - 1 - flanking_len < 0 or ltr_ref_end + flanking_len > len(ref_contigs[ref_name]):
                continue
            flanked_ltr_seq = ref_contigs[ref_name][ltr_ref_start - 1 - flanking_len: ltr_ref_end + flanking_len]

            if not flanked_ltr_candidates.__contains__(ref_name):
                flanked_ltr_candidates[ref_name] = []
            candidate_ltrs = flanked_ltr_candidates[ref_name]
            candidate_ltrs.append((candidate_ltr[0], candidate_ltr[1], candidate_ltr[2], candidate_ltr[3], flanked_ltr_seq))
            last_candidate_ltr = candidate_ltr
    return flanked_ltr_candidates

def multiple_alignment_blast(repeats_path, tools_dir):
    split_repeats_path = repeats_path[0]
    ref_db_path = repeats_path[1]
    blastn2Results_path = repeats_path[2]

    align_command = 'blastn -db ' + ref_db_path + ' -num_threads ' \
                    + str(1) + ' -query ' + split_repeats_path + ' -evalue 1e-20 -outfmt 6 > ' + blastn2Results_path
    os.system(align_command)

    return blastn2Results_path

def multiple_alignment_blast_and_get_copies_v1(repeats_path):
    split_repeats_path = repeats_path[0]
    split_ref_dir = repeats_path[1]
    blastn2Results_path = repeats_path[2]
    os.system('rm -f ' + blastn2Results_path)
    all_copies = {}
    repeat_names, repeat_contigs = read_fasta(split_repeats_path)
    remain_contigs = repeat_contigs
    for chr_name in os.listdir(split_ref_dir):
        if len(remain_contigs) > 0:
            if not str(chr_name).endswith('.fa'):
                continue
            chr_path = split_ref_dir + '/' + chr_name
            align_command = 'blastn -db ' + chr_path + ' -num_threads ' \
                            + str(1) + ' -query ' + split_repeats_path + ' -evalue 1e-20 -outfmt 6 > ' + blastn2Results_path
            os.system(align_command)
            # 由于我们只需要100个拷贝，因此如果有序列已经满足了，就不需要进行后续的比对了，这样在mouse这样的高拷贝大基因组上减少运行时间
            cur_all_copies = get_copies_v1(blastn2Results_path, split_repeats_path, '')
            for query_name in cur_all_copies.keys():
                copy_list = cur_all_copies[query_name]
                if query_name in all_copies:
                    prev_copy_list = all_copies[query_name]
                else:
                    prev_copy_list = []
                update_copy_list = prev_copy_list + copy_list
                all_copies[query_name] = update_copy_list
                if len(update_copy_list) >= 100:
                    del repeat_contigs[query_name]
            remain_contigs = repeat_contigs
            store_fasta(remain_contigs, split_repeats_path)

    # all_copies = get_copies_v1(blastn2Results_path, split_repeats_path, '')
    return all_copies

def multiple_alignment_blast_and_get_copies(repeats_path):
    split_repeats_path = repeats_path[0]
    ref_db = repeats_path[1]
    blastn2Results_path = repeats_path[2]
    os.system('rm -f ' + blastn2Results_path)
    all_copies = None
    repeat_names, repeat_contigs = read_fasta(split_repeats_path)
    if len(repeat_contigs) > 0:
        align_command = 'blastn -db ' + ref_db + ' -num_threads ' \
                        + str(1) + ' -query ' + split_repeats_path + ' -evalue 1e-20 -outfmt 6 >> ' + blastn2Results_path
        os.system(align_command)
        all_copies = get_copies_v1(blastn2Results_path, split_repeats_path, '')
    return all_copies


def judge_itr_structure(TSD_set, orig_seq, name, raw_tir_start, raw_tir_end, cur_candidate_TIRs_path, TRsearch_dir, tir_tsd_dir, plant):
    # TIR with short tir：(hAT, 5-27bp tir, 8bp tsd, len<4kb), (Mutator, long/short tir, 9-11bp tsd, ), (CACTA, 5bp tir, 2-3bp tsd, )

    # (left_tsd_start, left_tsd_end, left_tsd_seq, right_tsd_start, right_tsd_end, right_tsd_seq, tir_start, tir_end, tir_len)
    # Sort according to the distance between TIR_start and TIR_end and the original boundary, and the nearest one comes first.
    TSD_set = sorted(TSD_set, key=lambda x: abs(x[6] - raw_tir_start) + abs(x[7] - raw_tir_end))

    itr_contigs = {}
    cur_tsd_seq = ''
    for i, tsd_info in enumerate(TSD_set):
        left_tsd_start = tsd_info[0]
        left_tsd_end = tsd_info[1]
        left_tsd = tsd_info[2]
        cur_tsd_seq = left_tsd
        if left_tsd.__contains__('NN'):
            continue

        right_tsd_start = tsd_info[3]
        right_tsd_end = tsd_info[4]
        right_tsd = tsd_info[5]

        tir_start = tsd_info[6]
        tir_end = tsd_info[7]
        # tir_contig = {}
        tir_seq = orig_seq[tir_start - 1: tir_end]

        if len(tir_seq) < 100:
            continue

        tir_start = 1
        tir_end = len(tir_seq)
        first_5bp = tir_seq[tir_start - 1: tir_start + 4]
        last_5bp = getReverseSequence(tir_seq[tir_end - 5: tir_end])
        first_3bp = tir_seq[tir_start - 1: tir_start + 2]
        last_3bp = getReverseSequence(tir_seq[tir_end - 3: tir_end])
        tsd_len = len(left_tsd)

        if first_5bp == last_5bp:
            # hAT
            if tsd_len == 8 and len(tir_seq) < 4000:
                # name += ' Length itr=5'
                itr_contigs[name] = tir_seq
                break
            # Mutator
            elif tsd_len >= 9 and tsd_len <= 11:
                # name += ' Length itr=5'
                itr_contigs[name] = tir_seq
                break
            # CACTA (plant -> CACT[A/G], animal/fungi -> CCC)
            elif plant == 1 and tsd_len == 3 and (first_5bp == 'CACTA' or first_5bp == 'CACTG'):
                # name += ' Length itr=5'
                itr_contigs[name] = tir_seq
                break
        elif first_3bp == last_3bp and plant == 0 and tsd_len == 2 and (first_3bp == 'CCC'):
            # name += ' Length itr=3'
            itr_contigs[name] = tir_seq
            break
        else:
            # Determine if there is a TIR structure using itrsearch.
            # Search if there is a TIR structure in cur_candidate_TIRs. If it exists, we can skip evaluating other copies.
            output = run_itrsearch(TRsearch_dir, tir_seq)
            output = str(output)
            if output != '':
                seq = output.split('\n')[1]
                itr_contigs[name] = seq
                break
    return itr_contigs, cur_tsd_seq


def get_short_tir_contigs(cur_itr_contigs, plant):
    # Save sequences with short TIR structures.
    short_itr_contigs = {}
    for name in cur_itr_contigs.keys():
        tir_seq = cur_itr_contigs[name]
        # TIRs with short TIRs include: (hAT, 5-27bp TIR, 8bp TSD, length < 4kb), (Mutator, long/short TIR, 9-11bp TSD), (CACTA, 5bp TIR, 2-3bp TSD).
        # Determine if this TIR sequence has a TIR terminal structure.
        tir_start = 1
        tir_end = len(tir_seq)
        first_5bp = tir_seq[tir_start - 1: tir_start + 4]
        last_5bp = getReverseSequence(tir_seq[tir_end - 5: tir_end])
        first_3bp = tir_seq[tir_start - 1: tir_start + 2]
        last_3bp = getReverseSequence(tir_seq[tir_end - 3: tir_end])
        tsd_len = 0
        parts = name.split('-tsd_')
        if len(parts) > 1:
            sub_parts = parts[1].split('-')
            tsd_seq = sub_parts[0]
            tsd_len = len(tsd_seq)

        if first_5bp == last_5bp:
            # hAT
            if tsd_len == 8 and len(tir_seq) < 4000:
                # name += '-tsd_' + str(tsd_len)
                short_itr_contigs[name] = tir_seq
            # Mutator
            elif tsd_len >= 9 and tsd_len <= 11:
                # name += '-tsd_' + str(tsd_len)
                short_itr_contigs[name] = tir_seq
            # CACTA (plant -> CACT[A/G], animal/fungi -> CCC)
            elif plant == 1 and tsd_len == 3 and (first_5bp == 'CACTA' or first_5bp == 'CACTG'):
                # name += '-tsd_' + str(tsd_len)
                short_itr_contigs[name] = tir_seq
        elif first_3bp == last_3bp and plant == 0 and tsd_len == 2 and (first_3bp == 'CCC'):
            # name += '-tsd_' + str(tsd_len)
            short_itr_contigs[name] = tir_seq
    return short_itr_contigs


def filter_large_gap_tirs(input, output):
    contignames, contigs = read_fasta_v1(input)
    for name in contignames:
        parts = name.split(' ')
        # ITR(1,61)..(166,113)
        ITR_info = parts[1].replace('ITR', '').replace('(', '').replace(')', '')
        ITR_info_parts = ITR_info.split('..')
        ITR_left_pos_parts = ITR_info_parts[0].split(',')
        ITR_right_pos_parts = ITR_info_parts[1].split(',')
        lITR_start = int(ITR_left_pos_parts[0])
        lITR_end = int(ITR_left_pos_parts[1])
        lITR_len = lITR_end - lITR_start + 1
        rITR_start = int(ITR_right_pos_parts[0])
        rITR_end = int(ITR_right_pos_parts[1])
        rITR_len = rITR_start - rITR_end + 1
        if abs(rITR_len-lITR_len) > 2:
            del contigs[name]
    store_fasta(contigs, output)


def search_confident_tir_batch(split_file, flanking_len, tir_tsd_dir, TRsearch_dir, partition_index, plant):
    cur_contignames, cur_contigs = read_fasta(split_file)

    all_copies_itr_contigs = {}
    all_candidate_TIRs_path = tir_tsd_dir + '/' + str(partition_index) + '.fa'
    short_candidate_TIRs_path = tir_tsd_dir + '/' + str(partition_index) + '_s.fa'
    for query_name in cur_contignames:
        seq = cur_contigs[query_name]

        if seq.__contains__('NNNNNNNNNN'):
            continue

        tir_start = flanking_len + 1
        tir_end = len(seq) - flanking_len
        # Find all possible TSD sequences, calculate the distance between each sequence boundary
        # and the original boundary, and store it in the header
        tsd_search_distance = flanking_len
        cur_itr_contigs = search_confident_tir_v4(seq, tir_start, tir_end, tsd_search_distance, query_name, plant)
        all_copies_itr_contigs.update(cur_itr_contigs)

    # Save the sequence of short TIR and submit it toitrsearch to determine the length of TIR
    short_itr_contigs = get_short_tir_contigs(all_copies_itr_contigs, plant)
    store_fasta(short_itr_contigs, short_candidate_TIRs_path)
    short_copies_out, short_copies_log = run_itrsearch(TRsearch_dir, short_candidate_TIRs_path, tir_tsd_dir)
    raw_short_copies_out_name, raw_short_copies_out_contigs = read_fasta_v1(short_copies_out)

    # The remaining sequences are handed over to itrsearch to search for TIR structures
    for name in short_itr_contigs.keys():
        del all_copies_itr_contigs[name]
    store_fasta(all_copies_itr_contigs, all_candidate_TIRs_path)
    all_copies_out, all_copies_log = run_itrsearch(TRsearch_dir, all_candidate_TIRs_path, tir_tsd_dir)

    # Filter out sequences with excessive differences in terminal TIR length
    # filter_large_gap_tirs(all_copies_out, all_copies_out)

    # Record the TIR length corresponding to each query
    TIR_len_dict = {}
    raw_all_copies_out_name, raw_all_copies_out_contigs = read_fasta_v1(all_copies_out)
    all_copies_out_name, all_copies_out_contigs = read_fasta(all_copies_out)
    for name in raw_all_copies_out_name:
        query_name = name.split(' ')[0]
        tir_len = int(name.split('Length itr=')[1])
        TIR_len_dict[query_name] = tir_len
    for name in raw_short_copies_out_name:
        query_name = name.split(' ')[0]
        tir_len = int(name.split('Length itr=')[1])
        TIR_len_dict[query_name] = tir_len

    # Analyze theitrsearch log file and extract the sequence names of alignment offsets
    fake_tirs = get_fake_tirs(all_copies_log)
    for name in all_copies_out_name:
        if name in fake_tirs:
            del all_copies_out_contigs[name]
    all_copies_out_contigs.update(short_itr_contigs)

    # To reduce the amount of subsequent computation, we will perform grouping here,
    # selecting the sequences with the longest TIR and TSD for each group.
    # Group by query_name, and only select one sequence for each group,
    # which is the one with the best combined copy number and TSD.
    # Group all_copies_out_contigs by query_name.
    # group_copies_contigs -> {query_name: {name: seq}}
    group_copies_contigs = {}
    for cur_name in all_copies_out_contigs.keys():
        query_name = cur_name.split('-C_')[0]
        if not group_copies_contigs.__contains__(query_name):
            group_copies_contigs[query_name] = {}
        cur_copies_out_contigs = group_copies_contigs[query_name]
        cur_copies_out_contigs[cur_name] = all_copies_out_contigs[cur_name]

    filter_dup_itr_contigs = {}
    for query_name in group_copies_contigs.keys():
        cur_copies_out_contigs = group_copies_contigs[query_name]
        # Select the one with the smallest distance
        cur_contigs = filter_dup_itr_v3(cur_copies_out_contigs, TIR_len_dict)
        filter_dup_itr_contigs.update(cur_contigs)
    return filter_dup_itr_contigs

def get_fake_tirs(itrsearch_log):
    fake_tirs = set()
    alignments = {}
    line_count = 0
    query_name = ''
    with open(itrsearch_log, 'r') as f_r:
        for line in f_r:
            line_count += 1
            if line.startswith('load sequence'):
                parts = line.split('\t')
                query_name = parts[0].split(' ')[3]
                line_count = 0
            if line_count == 3 or line_count == 4 or line_count == 5:
                if line.strip() == '':
                    continue
                if query_name != '':
                    if not alignments.__contains__(query_name):
                        alignments[query_name] = []
                    details = alignments[query_name]
                    details.append(line)
    f_r.close()

    for query_name in alignments.keys():
        details = alignments[query_name]
        if len(details) != 3:
            continue
        query_seq = details[0]
        align_seq = details[1]
        target_seq = details[2]
        query_parts = query_seq.split(' ')
        target_parts = target_seq.split(' ')
        if len(query_parts) > 7 and len(target_parts) > 7:
            if query_seq[8] == '-' or target_seq[8] == '-' or (align_seq[8] != '|' and align_seq[9] != '|'):
                fake_tirs.add(query_name)
    return fake_tirs

def search_tsd_batch(cur_segments, flanking_len):
    all_copies_ltr_contigs = {}
    # {ref_name: [(left_ltr_start, left_ltr_end, right_ltr_start, right_ltr_end, flanked_ltr_seq)]}
    for item in cur_segments:
        ref_name = item[0]
        cur_candidate_LTRs = item[1]

        for copy_index, copy in enumerate(cur_candidate_LTRs):
            orig_seq = str(copy[4])

            if orig_seq.__contains__('NNNNNNNNNN'):
                continue

            ltr_len = copy[3]-copy[0]+1
            ltr_start = flanking_len + 1
            ltr_end = flanking_len + ltr_len

            tsd_search_distance = 50
            cur_ltr_contigs = search_confident_ltr(orig_seq, ltr_start, ltr_end, tsd_search_distance, ref_name, copy_index)
            all_copies_ltr_contigs.update(cur_ltr_contigs)
    return all_copies_ltr_contigs

def search_tsd_ltr_batch(all_copies_ltr_contigs, ltr_tsd_dir, partition_index, TRsearch_dir):
    all_candidate_LTRs_path = ltr_tsd_dir + '/' + str(partition_index) + '.fa'
    store_fasta(all_copies_ltr_contigs, all_candidate_LTRs_path)
    run_ltrsearch(TRsearch_dir, all_candidate_LTRs_path, ltr_tsd_dir)
    all_copies_out = all_candidate_LTRs_path + '.ltr'
    all_copies_out_name, all_copies_out_contigs = read_fasta_v1(all_copies_out)
    return all_copies_out_contigs

def rename_fasta(input, output, header='N'):
    names, contigs = read_fasta(input)
    node_index = 0
    with open(output, 'w') as f_save:
        for name in names:
            seq = contigs[name]
            f_save.write('>'+header+'_'+str(node_index)+'\n'+seq+'\n')
            node_index += 1
    f_save.close()

def rename_reference(input, output, chr_name_map):
    names, contigs = read_fasta(input)
    chr_name_dict = {}
    ref_index = 0
    with open(output, 'w') as f_save:
        for name in names:
            seq = contigs[name]
            new_name = 'chr_'+str(ref_index)
            f_save.write('>'+new_name+'\n'+seq+'\n')
            ref_index += 1
            chr_name_dict[new_name] = name
    f_save.close()
    with open(chr_name_map, 'w') as f_save:
        for new_name in chr_name_dict.keys():
            f_save.write(new_name+'\t'+chr_name_dict[new_name]+'\n')
    f_save.close()

def search_candidate_ltr(copies_out_contigs):
    all_copies_out_contigs = {}
    for name in copies_out_contigs.keys():
        seq = copies_out_contigs[name]

        parts = name.split(' ')
        orig_name = parts[0]
        # LTR(1,594)..(2154,2747)
        LTR_info = parts[1].replace('LTR', '').replace('(', '').replace(')', '')
        LTR_info_parts = LTR_info.split('..')
        LTR_left_pos_parts = LTR_info_parts[0].split(',')
        LTR_right_pos_parts = LTR_info_parts[1].split(',')
        lLTR_start = int(LTR_left_pos_parts[0])
        lLTR_end = int(LTR_left_pos_parts[1])
        rLTR_start = int(LTR_right_pos_parts[0])
        rLTR_end = int(LTR_right_pos_parts[1])
        new_query_name = orig_name + '-lLTRstart_' + str(lLTR_start) + '-lLTRend_' + str(lLTR_end) \
                         + '-rLTRstart_' + str(rLTR_start) + '-rLTRend_' + str(rLTR_end)
        all_copies_out_contigs[new_query_name] = seq


    candidate_LTRs = {}
    candidate_cut_LTRs = {}
    # group
    # group_copies_contigs -> {query_name: {name: seq}}
    group_copies_contigs = {}
    for cur_name in all_copies_out_contigs.keys():
        query_name = cur_name.split('-C_')[0]
        if not group_copies_contigs.__contains__(query_name):
            group_copies_contigs[query_name] = {}
        cur_copies_out_contigs = group_copies_contigs[query_name]
        cur_copies_out_contigs[cur_name] = all_copies_out_contigs[cur_name]

    for query_name in group_copies_contigs.keys():
        cur_copies_out_contigs = group_copies_contigs[query_name]
        # 1. Merge and select the sequence with the smallest distance (choose the longest in case of tie)
        # as the representative sequence for this copy.
        # copies_candidate -> {copy_index: (min_distance_seq_name, min_distance, tsd)}
        min_distance = 10000000
        min_distance_seq_len = 0
        min_distance_name = ''
        for name in cur_copies_out_contigs.keys():
            parts = name.split('-lLTRstart_')
            orig_name = parts[0]
            cur_distance = int(orig_name.split('-distance_')[1].split('-')[0])
            tsd_seq = orig_name.split('-tsd_')[1].split('-')[0]
            seq_len = len(cur_copies_out_contigs[name])
            if (cur_distance == min_distance and seq_len > min_distance_seq_len) or cur_distance < min_distance:
                min_distance = cur_distance
                min_distance_seq_len = seq_len
                min_distance_name = name
        if min_distance_name != '':
            seq = cur_copies_out_contigs[min_distance_name]

            parts = min_distance_name.split('-lLTRstart_')
            orig_name = parts[0]
            candidate_LTRs[orig_name] = seq

            # LTR(1,594)..(2154,2747)
            lLTR_start = int(min_distance_name.split('-lLTRstart_')[1].split('-')[0])
            lLTR_end = int(min_distance_name.split('-lLTRend_')[1].split('-')[0])
            rLTR_start = int(min_distance_name.split('-rLTRstart_')[1].split('-')[0])
            rLTR_end = int(min_distance_name.split('-rLTRend_')[1].split('-')[0])

            left_LTR = seq[lLTR_start - 1: lLTR_end]
            right_LTR = seq[rLTR_start - 1: rLTR_end]

            LTR_internal = seq[lLTR_end: rLTR_start - 1]

            left_LTR_name = orig_name + '-lLTR'
            right_LTR_name = orig_name + '-rLTR'
            internal_query_name = orig_name + '-ILTR'

            candidate_cut_LTRs[left_LTR_name] = left_LTR
            candidate_cut_LTRs[right_LTR_name] = right_LTR
            candidate_cut_LTRs[internal_query_name] = LTR_internal
    return candidate_LTRs, candidate_cut_LTRs

def search_confident_ltr(orig_seq, raw_ltr_start, raw_ltr_end, tsd_search_distance, ref_name, copy_index):
    ltr_contigs = {}
    orig_seq_len = len(orig_seq)
    ltr_starts = []
    ltr_ends = []

    for i in range(raw_ltr_start - tsd_search_distance, raw_ltr_start + tsd_search_distance + 1):
        if i >= 1 and i <= orig_seq_len:
            ltr_starts.append(i)

    for i in range(raw_ltr_end - tsd_search_distance, raw_ltr_end + tsd_search_distance + 1):
        if i >= 1 and i <= orig_seq_len:
            ltr_ends.append(i)

    TSD_set = set()
    for ltr_start in ltr_starts:
        for ltr_end in ltr_ends:
            TSDsearch_ltr(orig_seq, ltr_start, ltr_end, TSD_set)

    # (left_tsd_start, left_tsd_end, left_tsd_seq, right_tsd_start, right_tsd_end, right_tsd_seq, ltr_start, ltr_end, ltr_len)
    # ltr_start, ltr_end与原始边界的距离进行排序，越近的排在前面
    TSD_set = sorted(TSD_set, key=lambda x: abs(x[6] - raw_ltr_start) + abs(x[7] - raw_ltr_end))

    for i, tsd_info in enumerate(TSD_set):
        if i >= 20:
            break
        left_tsd_start = tsd_info[0]
        left_tsd_end = tsd_info[1]
        left_tsd = tsd_info[2]

        if left_tsd.__contains__('NN'):
            continue

        right_tsd_start = tsd_info[3]
        right_tsd_end = tsd_info[4]
        right_tsd = tsd_info[5]

        ltr_start = tsd_info[6]
        ltr_end = tsd_info[7]
        # tir_contig = {}
        ltr_seq = orig_seq[ltr_start - 1: ltr_end]

        distance = abs(ltr_start - raw_ltr_start) + abs(ltr_end - raw_ltr_end)

        if len(ltr_seq) < 100:
            continue

        new_query_name = ref_name + '-ltr_' + str(copy_index) + '-C_' + str(i) + '-tsd_' + left_tsd + '-distance_' + str(distance)
        ltr_contigs[new_query_name] = ltr_seq
    return ltr_contigs

def search_confident_tir(orig_seq, raw_tir_start, raw_tir_end, tsd_search_distance, query_name, plant):
    itr_contigs = {}
    orig_seq_len = len(orig_seq)
    tir_starts = []
    tir_ends = []

    for i in range(raw_tir_start - tsd_search_distance, raw_tir_start + tsd_search_distance + 1):
        if i >= 1 and i <= orig_seq_len:
            tir_starts.append(i)

    for i in range(raw_tir_end - tsd_search_distance, raw_tir_end + tsd_search_distance + 1):
        if i >= 1 and i <= orig_seq_len:
            tir_ends.append(i)

    TSD_set = set()
    for tir_start in tir_starts:
        for tir_end in tir_ends:
            TSDsearch_v2(orig_seq, tir_start, tir_end, TSD_set, plant)

    # (left_tsd_start, left_tsd_end, left_tsd_seq, right_tsd_start, right_tsd_end, right_tsd_seq, tir_start, tir_end, tir_len)
    TSD_set = sorted(TSD_set, key=lambda x: abs(x[6] - raw_tir_start) + abs(x[7] - raw_tir_end))

    query_dist = []
    for i, tsd_info in enumerate(TSD_set):
        left_tsd_start = tsd_info[0]
        left_tsd_end = tsd_info[1]
        left_tsd = tsd_info[2]

        if left_tsd.__contains__('NN'):
            continue

        right_tsd_start = tsd_info[3]
        right_tsd_end = tsd_info[4]
        right_tsd = tsd_info[5]

        tir_start = tsd_info[6]
        tir_end = tsd_info[7]

        # tir_contig = {}
        tir_seq = orig_seq[tir_start - 1: tir_end]
        distance = abs(tir_start - raw_tir_start) + abs(tir_end - raw_tir_end)

        if len(tir_seq) < 100:
            continue

        # new_query_name = query_name + '-C_' + str(copy_index) + '_' + str(i) + '-tsd_' + left_tsd + '-distance_' + str(distance)
        # itr_contigs[new_query_name] = tir_seq

        if tir_seq[0:2] == 'TG' and tir_seq[-2:] == 'CA':
            continue

        if str(tir_seq).startswith('TATATATA') or str(tir_seq).startswith('ATATATAT'):
            continue

        tir_start_5base = orig_seq[tir_start - 1: tir_start + 4]
        tir_end_5base = orig_seq[tir_end - 5: tir_end]
        if allow_mismatch(getReverseSequence(tir_start_5base), tir_end_5base, 1):
            new_query_name = query_name + '-C_' + str(i) + '-tsd_' + left_tsd + '-distance_'+ str(distance)
            itr_contigs[new_query_name] = tir_seq
            query_dist.append((new_query_name, distance))

    max_top_num = 10
    query_dist.sort(key=lambda x: x[1])
    top_itr_contigs = {}
    for i, item in enumerate(query_dist):
        if i >= max_top_num:
            break
        query_name = item[0]
        top_itr_contigs[query_name] = itr_contigs[query_name]
    return top_itr_contigs

def search_confident_tir_v4(orig_seq, raw_tir_start, raw_tir_end, tsd_search_distance, query_name, plant):
    # Change the coordinates to start with 0
    raw_tir_start -= 1
    raw_tir_end -= 1

    itr_contigs = {}
    orig_seq_len = len(orig_seq)
    TSD_set = set()
    # 1. take the sequence near the starting and ending positions
    left_start = raw_tir_start - tsd_search_distance
    if left_start < 0:
        left_start = 0
    left_end = raw_tir_start + tsd_search_distance + 1
    left_round_seq = orig_seq[left_start: left_end]
    # Get the position offset relative to the entire sequence, which is used to correct the position of the TSD boundary later
    left_offset = left_start
    right_start = raw_tir_end - tsd_search_distance
    if right_start < 0:
        right_start = 0
    right_end = raw_tir_end + tsd_search_distance + 1
    right_round_seq = orig_seq[right_start: right_end]
    right_offset = right_start

    # 2.Cut the sequences on the left and right into k-mers, store the k-mers on the left in a dict,
    # and then iterate through the k-mers on the right, determining whether they match the k-mers on the left.
    # If they match, they are a candidate TSD, and the position information is recorded.
    TIR_TSDs = [2, 3, 4, 5, 6, 8, 9, 10, 11]
    # exist_tsd -> {'TAA': {'left_pos': 100, 'right_pos': 200}}
    exist_tsd = {}
    for k_num in TIR_TSDs:
        for i in range(len(left_round_seq) - k_num + 1):
            left_kmer = left_round_seq[i: i + k_num]
            cur_pos = left_offset + i + k_num
            if cur_pos < 0 or cur_pos > orig_seq_len-1:
                continue
            if not exist_tsd.__contains__(left_kmer):
                exist_tsd[left_kmer] = {}
            pos_dict = exist_tsd[left_kmer]
            if not pos_dict.__contains__('left_pos'):
                pos_dict['left_pos'] = cur_pos
            else:
                prev_pos = pos_dict['left_pos']
                if abs(cur_pos-raw_tir_start) < abs(prev_pos-raw_tir_start):
                    pos_dict['left_pos'] = cur_pos
            exist_tsd[left_kmer] = pos_dict
    for k_num in TIR_TSDs:
        for i in range(len(right_round_seq) - k_num + 1):
            right_kmer = right_round_seq[i: i + k_num]
            cur_pos = right_offset + i - 1
            if cur_pos < 0 or cur_pos > orig_seq_len - 1:
                continue
            if exist_tsd.__contains__(right_kmer):
                # This is TSD
                pos_dict = exist_tsd[right_kmer]
                if not pos_dict.__contains__('right_pos'):
                    pos_dict['right_pos'] = cur_pos
                else:
                    prev_pos = pos_dict['right_pos']
                    # Determine who is closer to the original boundary
                    if abs(cur_pos - raw_tir_end) < abs(prev_pos - raw_tir_end):
                        pos_dict['right_pos'] = cur_pos
                exist_tsd[right_kmer] = pos_dict
                # Determine whether this TSD meets some basic requirements.
                tir_start = pos_dict['left_pos']
                tir_end = pos_dict['right_pos']
                first_3bp = orig_seq[tir_start: tir_start + 3]
                last_3bp = orig_seq[tir_end - 2: tir_end+1]
                if (k_num != 2 and k_num != 4) \
                        or (k_num == 4 and right_kmer == 'TTAA') \
                        or (k_num == 2 and (right_kmer == 'TA' or (plant == 0 and first_3bp == 'CCC' and last_3bp == 'GGG'))):
                    TSD_set.add((right_kmer, tir_start, tir_end))

    # sort by distance
    TSD_set = sorted(TSD_set, key=lambda x: abs(x[1] - raw_tir_start) + abs(x[2] - raw_tir_end))

    query_dist = []
    for i, tsd_info in enumerate(TSD_set):
        left_tsd = tsd_info[0]

        if left_tsd.__contains__('NN'):
            continue

        tir_start = tsd_info[1]
        tir_end = tsd_info[2]

        tir_seq = orig_seq[tir_start: tir_end+1]
        distance = abs(tir_start - raw_tir_start) + abs(tir_end - raw_tir_end)

        if len(tir_seq) < 100:
            continue

        # Filter out TIR sequences with the "TG..CA" motif, most of which should be false positives
        if tir_seq[0:2] == 'TG' and tir_seq[-2:] == 'CA':
            continue

        # Filter out TIR that starts and ends with "TATATATA"
        if str(tir_seq).startswith('TATATATA') or str(tir_seq).startswith('ATATATAT'):
            continue

        new_query_name = query_name + '-C_' + str(i) + '-tsd_' + left_tsd + '-distance_' + str(distance)
        itr_contigs[new_query_name] = tir_seq
        query_dist.append((new_query_name, distance))

    max_top_num = 100
    query_dist.sort(key=lambda x: x[1])
    top_itr_contigs = {}
    for i, item in enumerate(query_dist):
        if i >= max_top_num:
            break
        query_name = item[0]
        top_itr_contigs[query_name] = itr_contigs[query_name]
    return top_itr_contigs


def is_overlapped(s1, e1, s2, e2):
    if (s1 <= s2 and e1 >= s2) or (s1 >= s2 and e1 <= e2) or (s1 <= s2 and e1 >= e2) or (s1 <= e2 and e2 <= e1):
        return True
    else:
        return False


def overlap_with_boundary(q_start, q_end, s_start, s_end, flanking_len, flanking_region_distance, orig_query_len,
                          orig_subject_len):
    query_start_covered = False
    query_end_covered = False
    subject_start_covered = False
    subject_end_covered = False

    query_start_left = flanking_len + 1 - flanking_region_distance
    query_start_right = flanking_len + 1 + flanking_region_distance
    query_end_left = orig_query_len + flanking_len - flanking_region_distance
    query_end_right = orig_query_len + flanking_len + flanking_region_distance

    subject_start_left = flanking_len + 1 - flanking_region_distance
    subject_start_right = flanking_len + 1 + flanking_region_distance
    subject_end_left = orig_subject_len + flanking_len - flanking_region_distance
    subject_end_right = orig_subject_len + flanking_len + flanking_region_distance

    if is_overlapped(query_start_left, query_start_right, q_start, q_end):
        query_start_covered = True
    if is_overlapped(query_end_left, query_end_right, q_start, q_end):
        query_end_covered = True
    if is_overlapped(subject_start_left, subject_start_right, s_start, s_end):
        subject_start_covered = True
    if is_overlapped(subject_end_left, subject_end_right, s_start, s_end):
        subject_end_covered = True

    return query_start_covered, query_end_covered, subject_start_covered, subject_end_covered


def store_flank_align_groups(query_groups, flank_align_dir):
    for query_name in query_groups.keys():
        tmp_out = flank_align_dir + '/' + query_name + '.out'
        subject_groups = query_groups[query_name]
        with open(tmp_out, 'w') as f_save:
            for subject_name in subject_groups.keys():
                for item in subject_groups[subject_name]:
                    q_start = item[0]
                    q_end = item[1]
                    orig_query_len = item[2]
                    s_start = item[3]
                    s_end = item[4]
                    orig_subject_len = item[5]
                    flanking_len = item[6]
                    flanking_region_distance = item[7]
                    direct = item[8]
                    query_name = item[9]
                    subject_name = item[10]
                    f_save.write(query_name+'\t'+subject_name+'\t'+str(q_start)+'\t'+str(q_end)+'\t'+str(s_start)+'\t'+str(s_end)+'\t'+str(direct)+'\n')
        f_save.close()

def store_flank_align_groups_v1(groups, flank_align_dir):
    # groups -> {
    # orig_query_name: {
    #   query_name : {subject_name: []}
    #   }
    # }
    for orig_query_name in groups.keys():
        query_groups = groups[orig_query_name]
        for query_name in query_groups.keys():
            tmp_out = flank_align_dir + '/' + query_name + '.out'
            subject_groups = query_groups[query_name]
            with open(tmp_out, 'w') as f_save:
                for subject_name in subject_groups.keys():
                    for item in subject_groups[subject_name]:
                        q_start = item[0]
                        q_end = item[1]
                        orig_query_len = item[2]
                        s_start = item[3]
                        s_end = item[4]
                        orig_subject_len = item[5]
                        flanking_len = item[6]
                        flanking_region_distance = item[7]
                        direct = item[8]
                        query_name = item[9]
                        subject_name = item[10]
                        f_save.write(query_name+'\t'+subject_name+'\t'+str(q_start)+'\t'+str(q_end)+'\t'+str(s_start)+'\t'+str(s_end)+'\t'+str(direct)+'\n')
            f_save.close()

def flank_region_align_v5(candidate_sequence_path, real_TEs, flanking_len, reference, split_ref_dir,
                          TE_type, tmp_output_dir, threads, ref_index, log, subset_script_path, plant, debug,
                          iter_num, result_type='cons'):
    log.logger.info('------Determination of homology in regions outside the boundaries of ' + TE_type + ' copies')
    starttime = time.time()
    temp_dir = tmp_output_dir + '/' + TE_type + '_copies_' + str(ref_index) + '_' + str(iter_num)
    os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # We are considering that the current running time is too long, maybe it is related to submitting one sequence for Blastn alignment at a time.
    # We will try to combine 10 sequences together and run Blastn once.
    # To increase CPU utilization, we will submit one thread to process 10 sequences.
    batch_size = 10
    batch_id = 0
    names, contigs = read_fasta(candidate_sequence_path)
    total_names = set(names)
    split_files = []
    cur_contigs = {}
    for i, name in enumerate(names):
        cur_file = temp_dir + '/' + str(batch_id) + '.fa'
        cur_contigs[name] = contigs[name]
        if len(cur_contigs) == batch_size:
            store_fasta(cur_contigs, cur_file)
            split_files.append(cur_file)
            cur_contigs = {}
            batch_id += 1
    if len(cur_contigs) > 0:
        cur_file = temp_dir + '/' + str(batch_id) + '.fa'
        store_fasta(cur_contigs, cur_file)
        split_files.append(cur_file)
        batch_id += 1

    ref_contigs = {}
    # 遍历目录下的所有文件
    for filename in os.listdir(split_ref_dir):
        if filename.endswith('.fa'):
            file_path = os.path.join(split_ref_dir, filename)
            cur_names, cur_contigs = read_fasta(file_path)
            ref_contigs.update(cur_contigs)

    # ref_names, ref_contigs = read_fasta(reference)
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cur_split_files in split_files:
        job = ex.submit(get_full_length_copies, cur_split_files, split_ref_dir, debug)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_copies = {}
    for job in as_completed(jobs):
        cur_all_copies = job.result()
        all_copies.update(cur_all_copies)
    # extend copies
    batch_member_files = []
    new_all_copies = {}
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        for copy in copies:
            ref_name = copy[0]
            copy_ref_start = int(copy[1])
            copy_ref_end = int(copy[2])
            direct = copy[4]
            copy_len = copy_ref_end - copy_ref_start + 1
            if copy_ref_start - 1 - flanking_len < 0 or copy_ref_end + flanking_len > len(ref_contigs[ref_name]):
                continue
            copy_seq = ref_contigs[ref_name][copy_ref_start - 1 - flanking_len: copy_ref_end + flanking_len]
            if direct == '-':
                copy_seq = getReverseSequence(copy_seq)
            if len(copy_seq) < 100:
                continue
            new_name = ref_name + ':' + str(copy_ref_start) + '-' + str(copy_ref_end) + '(' + direct + ')'
            if not new_all_copies.__contains__(query_name):
                new_all_copies[query_name] = {}
            copy_contigs = new_all_copies[query_name]
            copy_contigs[new_name] = copy_seq
            new_all_copies[query_name] = copy_contigs
    for query_name in new_all_copies.keys():
        copy_contigs = new_all_copies[query_name]
        valid_query_filename = re.sub(r'[<>:"/\\|?*]', '-', query_name)
        cur_member_file = temp_dir + '/' + valid_query_filename + '.blast.bed.fa'
        store_fasta(copy_contigs, cur_member_file)
        query_seq = contigs[query_name]
        batch_member_files.append((query_name, query_seq, cur_member_file))

    # Determine whether the multiple sequence alignment of each copied file satisfies the homology rule
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for batch_member_file in batch_member_files:
        job = ex.submit(run_find_members_v8, batch_member_file, temp_dir, subset_script_path,
                        plant, TE_type, debug, result_type)
        jobs.append(job)
    ex.shutdown(wait=True)

    not_found_boundary = 0
    full_length1 = 0
    copy_nums = {}
    true_te_names = set()
    true_tes = {}
    low_copy_contigs = {}
    low_copy_path = temp_dir + '/low_copy_elements.fa'
    for job in as_completed(jobs):
        result_info = job.result()
        cur_name, cur_seq, info, copy_count = result_info
        if info == 'nb':
            not_found_boundary += 1
        elif info == 'fl1':
            full_length1 += 1
        elif info.startswith('copy_num:'):
            copy_num = int(info.split('copy_num:')[1])
            if not copy_nums.__contains__(copy_num):
                copy_nums[copy_num] = 0
            cur_copy_num = copy_nums[copy_num]
            cur_copy_num += 1
            copy_nums[copy_num] = cur_copy_num
        if cur_name is not None:
            if TE_type == 'tir':
                if cur_seq.startswith('TG') and cur_seq.endswith('CA'):
                    continue
                ############################
                # 对于拷贝数<=2的序列，大概率可能是假阳性，那我们判断它是否具有末端反向重复来过滤
                if copy_count <= 2:
                    low_copy_contigs[cur_name] = cur_seq
                else:
                    true_tes[cur_name] = cur_seq
                    true_te_names.add(cur_name)
                ############################
            else:
                true_tes[cur_name] = cur_seq
                true_te_names.add(cur_name)
    store_fasta(low_copy_contigs, low_copy_path)

    if TE_type == 'tir':
        # 找回具有TIR结构的低拷贝TIR
        cur_temp_dir = temp_dir + '/low_copy_itr'
        TRsearch_dir = cur_dir + '/tools'
        with_tir_path, no_tir_path = remove_no_tirs(low_copy_path, plant, TRsearch_dir, cur_temp_dir)
        with_tir_names, with_tir_contigs = read_fasta(with_tir_path)
        true_tes.update(with_tir_contigs)
        # print('recall by TIR structure: ' + str(len(with_tir_contigs)))
        # print(with_tir_names)

        # 找回具有完整domain的低拷贝TIR
        temp_dir = temp_dir + '/tir_domain'
        output_table = no_tir_path + '.tir_domain'
        tir_protein_db = cur_dir + '/library/TIRPeps.lib'
        get_domain_info(no_tir_path, tir_protein_db, output_table, threads, temp_dir)
        has_intact_protein_contigs = {}
        protein_names, protein_contigs = read_fasta(tir_protein_db)
        with open(output_table, 'r') as f_r:
            for i, line in enumerate(f_r):
                if i < 2:
                    continue
                parts = line.split('\t')
                te_name = parts[0]
                protein_name = parts[1]
                protein_start = int(parts[4])
                protein_end = int(parts[5])
                intact_protein_len = len(protein_contigs[protein_name])
                if float(abs(protein_end - protein_start)) / intact_protein_len >= 0.95:
                    has_intact_protein_contigs[te_name] = low_copy_contigs[te_name]
        true_tes.update(has_intact_protein_contigs)
        # print('recall by intact domain structure: ' + str(len(has_intact_protein_contigs)))
        # print(has_intact_protein_contigs.keys())
    store_fasta(true_tes, real_TEs)

    deleted_names = total_names.difference(true_te_names)
    if debug:
        print(deleted_names)
        print('deleted_names len: ' + str(len(deleted_names)))
        print('not found boundary num: ' + str(not_found_boundary) + ', full length 1: ' + str(full_length1))
        print(copy_nums)
    if debug != 1:
        os.system('rm -rf ' + temp_dir)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of determination of homology in regions outside the boundaries of  " + TE_type + " copies: %.8s s" % (dtime))

def multi_process_alignx(query_path, subject_path, blastnResults_path, blast_program_dir, tmp_output_dir, threads):
    tools_dir = ''

    tmp_blast_dir = tmp_output_dir + '/tmp_blast_test'
    os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    orig_names, orig_contigs = read_fasta(query_path)

    blast_db_command = blast_program_dir + '/bin/makeblastdb -dbtype prot -in ' + subject_path + ' > /dev/null 2>&1'
    os.system(blast_db_command)

    longest_repeat_files = []
    segments_cluster = divided_array(list(orig_contigs.items()), threads)
    for partition_index, cur_segments in enumerate(segments_cluster):
        single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
        print('current partition_index: ' + str(partition_index))
        if not os.path.exists(single_tmp_dir):
            os.makedirs(single_tmp_dir)
        split_repeat_file = single_tmp_dir + '/repeats_split.fa'
        cur_contigs = {}
        for item in cur_segments:
            cur_contigs[item[0]] = item[1]
        store_fasta(cur_contigs, split_repeat_file)
        repeats_path = (split_repeat_file, subject_path, single_tmp_dir + '/temp.out')
        longest_repeat_files.append(repeats_path)

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeat_files:
        job = ex.submit(multiple_alignment_blastx, file, blast_program_dir, tools_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    if os.path.exists(blastnResults_path):
        os.remove(blastnResults_path)

    for job in as_completed(jobs):
        cur_blastn2Results_path = job.result()
        os.system('cat ' + cur_blastn2Results_path + ' >> ' + blastnResults_path)

def multi_process_LINE(query_path, subject_path, candidate_LINE_path, blast_program_dir, tmp_blast_dir, threads, is_removed_dir=True):
    tools_dir = ''
    if is_removed_dir:
        os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    orig_names, orig_contigs = read_fasta(query_path)

    blast_db_command = blast_program_dir + '/bin/makeblastdb -dbtype prot -in ' + subject_path + ' > /dev/null 2>&1'
    os.system(blast_db_command)

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
        repeats_path = (split_repeat_file, subject_path, single_tmp_dir)
        longest_repeat_files.append(repeats_path)

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeat_files:
        job = ex.submit(multiple_alignment_blastx, file, blast_program_dir, tools_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    if os.path.exists(candidate_LINE_path):
        os.remove(candidate_LINE_path)
    for job in as_completed(jobs):
        cur_candidate_LINE_path = job.result()
        os.system('cat ' + cur_candidate_LINE_path + ' >> ' + candidate_LINE_path)

def multi_process_align_and_get_copies(query_path, subject_path, tmp_blast_dir, TE_type, threads, is_removed_dir=True, query_coverage=0.99, subject_coverage=0):
    tools_dir = ''
    if is_removed_dir:
        os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    orig_names, orig_contigs = read_fasta(query_path)

    blast_db_command = 'makeblastdb -dbtype nucl -in ' + subject_path + ' > /dev/null 2>&1'
    os.system(blast_db_command)

    longest_repeat_files = []
    file_index = 0
    cur_seq_index = 0
    cur_contigs = {}
    for name in orig_names:
        cur_contigs[name] = orig_contigs[name]
        cur_seq_index += 1
        if cur_seq_index >= 50:
            split_repeat_file = tmp_blast_dir + '/' + str(file_index) + '.fa'
            store_fasta(cur_contigs, split_repeat_file)
            output_file = tmp_blast_dir + '/' + str(file_index) + '.out'
            longest_repeat_files.append((split_repeat_file, subject_path, output_file))
            cur_contigs = {}
            file_index += 1
            cur_seq_index = 0
    if len(cur_contigs) > 0:
        split_repeat_file = tmp_blast_dir + '/' + str(file_index) + '.fa'
        store_fasta(cur_contigs, split_repeat_file)
        output_file = tmp_blast_dir + '/' + str(file_index) + '.out'
        longest_repeat_files.append((split_repeat_file, subject_path, output_file))

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeat_files:
        job = ex.submit(multiple_alignment_blast_and_get_copies, file)
        jobs.append(job)
    ex.shutdown(wait=True)

    all_copies = {}
    for job in as_completed(jobs):
        cur_copies = job.result()
        all_copies.update(cur_copies)

    return all_copies


def multi_process_align(query_path, subject_path, blastnResults_path, tmp_blast_dir, threads, is_removed_dir=True, is_remove_index=False):
    tools_dir = ''
    if is_removed_dir:
        os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    if os.path.exists(blastnResults_path):
        os.remove(blastnResults_path)

    orig_names, orig_contigs = read_fasta(query_path)

    blast_db_command = 'makeblastdb -dbtype nucl -in ' + subject_path + ' > /dev/null 2>&1'
    os.system(blast_db_command)

    longest_repeat_files = []
    segments_cluster = divided_array(list(orig_contigs.items()), threads)
    for partition_index, cur_segments in enumerate(segments_cluster):
        if len(cur_segments) <= 0:
            continue
        single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
        #print('current partition_index: ' + str(partition_index))
        if not os.path.exists(single_tmp_dir):
            os.makedirs(single_tmp_dir)
        split_repeat_file = single_tmp_dir + '/repeats_split.fa'
        cur_contigs = {}
        for item in cur_segments:
            cur_contigs[item[0]] = item[1]
        store_fasta(cur_contigs, split_repeat_file)
        repeats_path = (split_repeat_file, subject_path, single_tmp_dir + '/temp.out')
        longest_repeat_files.append(repeats_path)

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeat_files:
        job = ex.submit(multiple_alignment_blast, file, tools_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        cur_blastn2Results_path = job.result()
        os.system('cat ' + cur_blastn2Results_path + ' >> ' + blastnResults_path)

    if is_remove_index:
        os.system('rm -f ' + subject_path + '.*')

def remove_ltr_from_tir(confident_ltr_cut_path, confident_tir_path, threads, tmp_output_dir):
    subject_path = confident_ltr_cut_path
    query_path = confident_tir_path
    out_path = query_path + '.out'
    align_command = 'cd ' + tmp_output_dir + ' && RepeatMasker -lib ' + subject_path + ' -nolow -pa ' + str(threads) + ' ' + query_path + ' > /dev/null 2>&1'
    os.system(align_command)
    delete_tir_names = set()
    query_names, query_contigs = read_fasta(query_path)
    subject_names, subject_contigs = read_fasta(subject_path)
    if not os.path.exists(out_path) or os.path.getsize(out_path) <= 0:
        return
    with open(out_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = []
            for p in line.split(' '):
                if p.strip() != '':
                    parts.append(p)
            if len(parts) >= 15:
                direct = parts[8]
                query_name = parts[4]
                subject_name = parts[9]+'#'+parts[10]
                if direct == '+':
                    subject_start = int(parts[11])
                    subject_end = int(parts[12])
                    subject_len = subject_end-subject_start+1
                else:
                    subject_start = int(parts[13])
                    subject_end = int(parts[12])
                    subject_len = subject_end - subject_start + 1
                total_subject_len = len(subject_contigs[subject_name])
                if float(subject_len)/total_subject_len >= 0.95:
                    delete_tir_names.add(query_name)
    f_r.close()
    # print(delete_tir_names)
    # print(len(delete_tir_names))
    # remain_names = set(query_names).difference(delete_tir_names)
    # print(remain_names)
    # print(len(remain_names))
    for tir_name in delete_tir_names:
        del query_contigs[tir_name]
    store_fasta(query_contigs, query_path)


def search_boundary_homo_v4(valid_col_threshold, pos, matrix, row_num, col_num,
                            type, homo_threshold, int_homo_threshold, out_homo_threshold, debug, sliding_window_size):
    # We need a program that takes an alignment file 'align_file' and boundary positions 'start_pos' and 'end_pos' as inputs, and extracts effective 20 columns around the boundaries. It also checks if these 20 columns exhibit homology.
    # Key Definitions:
    # ① What is an effective column? A column that has at least half of the total copy count, i.e., at least total/2 non-empty bases.
    # ② How is homology calculated? If consistent bases exceed 80% of the total sequence count, the column is considered homologous; otherwise, it is not.
    # If there is homology in 10 out of 15bp outside the boundary, it is likely to be a false positive.

    # Functionality:
    # Given an alignment matrix and a starting column, search for effective columns, homologous columns (homologous columns are always effective columns) towards both ends, and count the number of homologous columns, continuous homologous columns, and continuous non-homologous columns.
    # If there are consecutive non-homologous columns within the boundary or consecutive homologous columns outside the boundary beyond the threshold, it is considered a false positive.
    # Record the base composition of each column.
    col_base_map = {}
    for col_index in range(col_num):
        if not col_base_map.__contains__(col_index):
            col_base_map[col_index] = {}
        base_map = col_base_map[col_index]
        # Calculate the base composition ratio in the current column.
        if len(base_map) == 0:
            for row in range(row_num):
                cur_base = matrix[row][col_index]
                if not base_map.__contains__(cur_base):
                    base_map[cur_base] = 0
                cur_count = base_map[cur_base]
                cur_count += 1
                base_map[cur_base] = cur_count
        if not base_map.__contains__('-'):
            base_map['-'] = 0

    search_len = 100
    if type == 'start':
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        col_index = pos
        homo_cols = []
        while valid_col_count < search_len and col_index < col_num / 2:
            # Starting from position 'pos', search for 15 effective columns to the right.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            no_gap_num = row_num - base_map['-']
            max_homo_ratio = 0
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                no_gap_num = row_num - base_map['-']
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / row_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo,
                     max_homo_ratio))
            col_index += 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align start right: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the left. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        cur_boundary = pos
        new_boundary_start = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            for item in window:
                cur_homo_ratio = item[5]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio) / sliding_window_size
            if avg_homo_ratio >= homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_start = window[0][0]
                break
        if new_boundary_start != pos and new_boundary_start != -1:
            if debug:
                print('align start right non-homology, new boundary: ' + str(new_boundary_start))
            cur_boundary = new_boundary_start

        col_index = cur_boundary - 1
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        homo_cols = []
        while valid_col_count < search_len and col_index >= 0:
            # Starting from position 'pos', search for 15 effective columns to the left.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            max_homo_base = None
            max_homo_ratio = 0
            no_gap_num = row_num - base_map['-']
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / row_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                        max_homo_base = base
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo, max_homo_ratio, max_homo_base))
            col_index -= 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align start left: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the left. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        homo_cols.reverse()
        new_boundary_start = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            for item in window:
                cur_homo_ratio = item[5]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio)/sliding_window_size

            if avg_homo_ratio >= homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_start = window[0][0]
                break
        if new_boundary_start != pos and new_boundary_start != -1:
            if debug:
                print('align start left homology, new boundary: ' + str(new_boundary_start))
            cur_boundary = new_boundary_start

        return True, cur_boundary
    else:
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        col_index = pos + 1
        homo_cols = []
        while valid_col_count < search_len and col_index < col_num:
            # Starting from position 'pos', search for 15 effective columns to the right.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            # If the number of non-empty rows exceeds the threshold, then it is an effective row.
            no_gap_num = row_num - base_map['-']
            max_homo_base = None
            max_homo_ratio = 0
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / row_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                        max_homo_base = base
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo, max_homo_ratio, max_homo_base))
            col_index += 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align end right: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the right. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        new_boundary_end = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            if i != 0:
                break
            window = homo_cols[i:i + sliding_window_size]

            avg_homo_ratio = 0
            for item in window:
                cur_homo_ratio = item[5]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio) / sliding_window_size

            if avg_homo_ratio >= out_homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_end = window[len(window)-1][0]
                break
        if new_boundary_end != pos and new_boundary_end != -1:
            if debug:
                print('align end right homology, new boundary: ' + str(new_boundary_end))
            return False, -1

        col_index = pos
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        homo_cols = []
        while valid_col_count < search_len and col_index >= col_num / 2:
            # Starting from position 'pos', search for 20 effective columns to the left.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            # If the number of non-empty rows exceeds the threshold, then it is an effective row.
            no_gap_num = row_num - base_map['-']
            max_homo_ratio = 0
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / row_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo, max_homo_ratio))
            col_index -= 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align end left: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the right. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        new_boundary_end = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            if i != 0:
                break
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            for item in window:
                cur_homo_ratio = item[5]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio) / sliding_window_size
            if avg_homo_ratio < int_homo_threshold:
                new_boundary_end = window[len(window)-1][0]
                break
        if new_boundary_end != pos and new_boundary_end != -1:
            if debug:
                print('align end left non-homology, new boundary: ' + str(new_boundary_end))
            return False, -1
        return True, pos

def search_boundary_homo_v5(valid_col_threshold, pos, matrix, row_num, col_num,
                            type, homo_threshold, int_homo_threshold, out_homo_threshold, debug, sliding_window_size):
    # We need a program that takes an alignment file 'align_file' and boundary positions 'start_pos' and 'end_pos' as inputs, and extracts effective 20 columns around the boundaries. It also checks if these 20 columns exhibit homology.
    # Key Definitions:
    # ① What is an effective column? A column that has at least half of the total copy count, i.e., at least total/2 non-empty bases.
    # ② How is homology calculated? If consistent bases exceed 80% of the total sequence count, the column is considered homologous; otherwise, it is not.
    # If there is homology in 10 out of 15bp outside the boundary, it is likely to be a false positive.

    # Functionality:
    # Given an alignment matrix and a starting column, search for effective columns, homologous columns (homologous columns are always effective columns) towards both ends, and count the number of homologous columns, continuous homologous columns, and continuous non-homologous columns.
    # If there are consecutive non-homologous columns within the boundary or consecutive homologous columns outside the boundary beyond the threshold, it is considered a false positive.
    # Record the base composition of each column.

    col_base_map = {}
    for col_index in range(col_num):
        if not col_base_map.__contains__(col_index):
            col_base_map[col_index] = {}
        base_map = col_base_map[col_index]
        # Calculate the base composition ratio in the current column.
        if len(base_map) == 0:
            for row in range(row_num):
                cur_base = matrix[row][col_index]
                if not base_map.__contains__(cur_base):
                    base_map[cur_base] = 0
                cur_count = base_map[cur_base]
                cur_count += 1
                base_map[cur_base] = cur_count
        if not base_map.__contains__('-'):
            base_map['-'] = 0

    search_len = 100
    if type == 'start':
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        col_index = pos
        homo_cols = []
        while valid_col_count < search_len and col_index < col_num / 2:
            # Starting from position 'pos', search for 15 effective columns to the right.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            no_gap_num = row_num - base_map['-']
            max_homo_ratio = 0
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                no_gap_num = row_num - base_map['-']
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / row_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo,
                     max_homo_ratio))
            col_index += 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align start right: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the left. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        cur_boundary = pos
        new_boundary_start = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            for item in window:
                cur_homo_ratio = item[5]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio) / sliding_window_size
            if avg_homo_ratio >= homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_start = window[0][0]
                break
        if new_boundary_start != pos and new_boundary_start != -1:
            if debug:
                print('align start right non-homology, new boundary: ' + str(new_boundary_start))
            cur_boundary = new_boundary_start

        col_index = cur_boundary - 1
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        homo_cols = []
        while valid_col_count < search_len and col_index >= 0:
            # Starting from position 'pos', search for 15 effective columns to the left.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            max_homo_base = None
            max_homo_ratio = 0
            no_gap_num = row_num - base_map['-']
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / row_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                        max_homo_base = base
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo, max_homo_ratio, max_homo_base))
            col_index -= 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align start left: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the left. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        homo_cols.reverse()
        new_boundary_start = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            for item in window:
                cur_homo_ratio = item[5]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio)/sliding_window_size

            if avg_homo_ratio >= homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_start = window[0][0]
                break
        if new_boundary_start != pos and new_boundary_start != -1:
            if debug:
                print('align start left homology, new boundary: ' + str(new_boundary_start))
            cur_boundary = new_boundary_start
        return True, cur_boundary
    else:
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        col_index = pos + 1
        homo_cols = []
        while valid_col_count < search_len and col_index < col_num:
            # Starting from position 'pos', search for 15 effective columns to the right.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            # If the number of non-empty rows exceeds the threshold, then it is an effective row.
            no_gap_num = row_num - base_map['-']
            max_homo_base = None
            max_homo_ratio = 0
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / row_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                        max_homo_base = base
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo, max_homo_ratio, max_homo_base))
            col_index += 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align end right: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the right. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        cur_boundary = pos
        homo_cols.reverse()
        new_boundary_end = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            if i != 0:
                break
            window = homo_cols[i:i + sliding_window_size]

            avg_homo_ratio = 0
            for item in window:
                cur_homo_ratio = item[5]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio) / sliding_window_size

            if avg_homo_ratio >= out_homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_end = window[0][0]
                break
        if new_boundary_end != pos and new_boundary_end != -1:
            if debug:
                print('align end right homology, new boundary: ' + str(new_boundary_end))
            cur_boundary = new_boundary_end
        return True, cur_boundary

        col_index = pos
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        homo_cols = []
        while valid_col_count < search_len and col_index >= col_num / 2:
            # Starting from position 'pos', search for 20 effective columns to the left.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            # If the number of non-empty rows exceeds the threshold, then it is an effective row.
            no_gap_num = row_num - base_map['-']
            max_homo_ratio = 0
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / row_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo, max_homo_ratio))
            col_index -= 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align end left: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the right. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        new_boundary_end = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            if i != 0:
                break
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            for item in window:
                cur_homo_ratio = item[5]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio) / sliding_window_size
            if avg_homo_ratio < int_homo_threshold:
                new_boundary_end = window[len(window)-1][0]
                break
        if new_boundary_end != pos and new_boundary_end != -1:
            if debug:
                print('align end left non-homology, new boundary: ' + str(new_boundary_end))
            cur_boundary = new_boundary_end
        return True, cur_boundary

def search_boundary_homo_v3(valid_col_threshold, pos, matrix, row_num, col_num,
                            type, homo_threshold, debug, sliding_window_size):
    # We need a program that takes an alignment file 'align_file' and boundary positions 'start_pos' and 'end_pos' as inputs, and extracts effective 20 columns around the boundaries. It also checks if these 20 columns exhibit homology.
    # Key Definitions:
    # ① What is an effective column? A column that has at least half of the total copy count, i.e., at least total/2 non-empty bases.
    # ② How is homology calculated? If consistent bases exceed 80% of the total sequence count, the column is considered homologous; otherwise, it is not.
    # If there is homology in 10 out of 15bp outside the boundary, it is likely to be a false positive.

    # Functionality:
    # Given an alignment matrix and a starting column, search for effective columns, homologous columns (homologous columns are always effective columns) towards both ends, and count the number of homologous columns, continuous homologous columns, and continuous non-homologous columns.
    # If there are consecutive non-homologous columns within the boundary or consecutive homologous columns outside the boundary beyond the threshold, it is considered a false positive.
    # Record the base composition of each column.
    col_base_map = {}
    for col_index in range(col_num):
        if not col_base_map.__contains__(col_index):
            col_base_map[col_index] = {}
        base_map = col_base_map[col_index]
        # Calculate the base composition ratio in the current column.
        if len(base_map) == 0:
            for row in range(row_num):
                cur_base = matrix[row][col_index]
                if not base_map.__contains__(cur_base):
                    base_map[cur_base] = 0
                cur_count = base_map[cur_base]
                cur_count += 1
                base_map[cur_base] = cur_count
        if not base_map.__contains__('-'):
            base_map['-'] = 0

    search_len = 100
    if type == 'start':
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        col_index = pos
        homo_cols = []
        while valid_col_count < search_len and col_index < col_num / 2:
            # Starting from position 'pos', search for 15 effective columns to the right.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            no_gap_num = row_num - base_map['-']
            max_homo_ratio = 0
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                no_gap_num = row_num - base_map['-']
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / row_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo,
                     max_homo_ratio))
            col_index += 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align start right: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the left. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        cur_boundary = pos
        new_boundary_start = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            first_candidate_boundary = -1
            for item in window:
                cur_homo_ratio = item[5]
                if cur_homo_ratio >= homo_threshold-0.1 and first_candidate_boundary == -1:
                    first_candidate_boundary = item[0]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio) / sliding_window_size
            if avg_homo_ratio >= homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_start = first_candidate_boundary
                break
        if new_boundary_start != cur_boundary and new_boundary_start != -1:
            if debug:
                print('align start right non-homology, new boundary: ' + str(new_boundary_start))
        cur_boundary = new_boundary_start

        col_index = cur_boundary
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        homo_cols = []
        while valid_col_count < search_len and col_index >= 0:
            # Starting from position 'pos', search for 15 effective columns to the left.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            max_homo_ratio = 0
            no_gap_num = row_num - base_map['-']
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / row_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo, max_homo_ratio))
            col_index -= 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align start left: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the left. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        homo_cols.reverse()
        new_boundary_start = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            first_candidate_boundary = -1
            for item in window:
                cur_homo_ratio = item[5]
                if cur_homo_ratio >= homo_threshold-0.1 and first_candidate_boundary == -1:
                    first_candidate_boundary = item[0]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio)/sliding_window_size
            if avg_homo_ratio >= homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_start = first_candidate_boundary
                break
        if new_boundary_start != cur_boundary and new_boundary_start != -1:
            if debug:
                print('align start left homology, new boundary: ' + str(new_boundary_start))
            cur_boundary = new_boundary_start

        return cur_boundary
    else:
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        col_index = pos
        homo_cols = []
        while valid_col_count < search_len and col_index < col_num:
            # Starting from position 'pos', search for 15 effective columns to the right.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            # If the number of non-empty rows exceeds the threshold, then it is an effective row.
            no_gap_num = row_num - base_map['-']
            max_homo_ratio = 0
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / row_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo, max_homo_ratio))
            col_index += 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align end right: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the right. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        cur_boundary = pos
        homo_cols.reverse()
        new_boundary_end = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            first_candidate_boundary = -1
            for item in window:
                cur_homo_ratio = item[5]
                if cur_homo_ratio >= homo_threshold-0.1 and first_candidate_boundary == -1:
                    first_candidate_boundary = item[0]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio) / sliding_window_size
            if avg_homo_ratio >= homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_end = first_candidate_boundary
                break
        if new_boundary_end != cur_boundary and new_boundary_end != -1:
            if debug:
                print('align end right homology, new boundary: ' + str(new_boundary_end))
            cur_boundary = new_boundary_end

        col_index = cur_boundary
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        homo_cols = []
        while valid_col_count < search_len and col_index >= col_num / 2:
            # Starting from position 'pos', search for 20 effective columns to the left.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            # If the number of non-empty rows exceeds the threshold, then it is an effective row.
            no_gap_num = row_num - base_map['-']
            max_homo_ratio = 0
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / row_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo, max_homo_ratio))
            col_index -= 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align end left: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the right. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        new_boundary_end = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            first_candidate_boundary = -1
            for item in window:
                cur_homo_ratio = item[5]
                if cur_homo_ratio >= homo_threshold-0.1 and first_candidate_boundary == -1:
                    first_candidate_boundary = item[0]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio) / sliding_window_size
            if avg_homo_ratio >= homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_end = first_candidate_boundary
                break
        if new_boundary_end != cur_boundary and new_boundary_end != -1:
            if debug:
                print('align end left non-homology, new boundary: ' + str(new_boundary_end))
        cur_boundary = new_boundary_end

        return cur_boundary


def search_boundary_homo_v6(valid_col_threshold, pos, matrix, row_num, col_num,
                            type, homo_threshold, debug, sliding_window_size):
    # We need a program that takes an alignment file 'align_file' and boundary positions 'start_pos' and 'end_pos' as inputs, and extracts effective 20 columns around the boundaries. It also checks if these 20 columns exhibit homology.
    # Key Definitions:
    # ① What is an effective column? A column that has at least half of the total copy count, i.e., at least total/2 non-empty bases.
    # ② How is homology calculated? If consistent bases exceed 80% of the total sequence count, the column is considered homologous; otherwise, it is not.
    # If there is homology in 10 out of 15bp outside the boundary, it is likely to be a false positive.

    # Functionality:
    # Given an alignment matrix and a starting column, search for effective columns, homologous columns (homologous columns are always effective columns) towards both ends, and count the number of homologous columns, continuous homologous columns, and continuous non-homologous columns.
    # If there are consecutive non-homologous columns within the boundary or consecutive homologous columns outside the boundary beyond the threshold, it is considered a false positive.
    # Record the base composition of each column.
    col_base_map = {}
    for col_index in range(col_num):
        if not col_base_map.__contains__(col_index):
            col_base_map[col_index] = {}
        base_map = col_base_map[col_index]
        # Calculate the base composition ratio in the current column.
        if len(base_map) == 0:
            for row in range(row_num):
                cur_base = matrix[row][col_index]
                if not base_map.__contains__(cur_base):
                    base_map[cur_base] = 0
                cur_count = base_map[cur_base]
                cur_count += 1
                base_map[cur_base] = cur_count
        if not base_map.__contains__('-'):
            base_map['-'] = 0

    search_len = 100
    if type == 'start':
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        col_index = pos
        homo_cols = []
        while valid_col_count < search_len and col_index < col_num / 2:
            # Starting from position 'pos', search for 15 effective columns to the right.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            no_gap_num = row_num - base_map['-']
            max_homo_ratio = 0
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                no_gap_num = row_num - base_map['-']
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / no_gap_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo,
                     max_homo_ratio))
            col_index += 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align start right: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the left. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        cur_boundary = pos
        new_boundary_start = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            first_candidate_boundary = -1
            for item in window:
                cur_homo_ratio = item[5]
                if cur_homo_ratio >= homo_threshold-0.1 and first_candidate_boundary == -1:
                    first_candidate_boundary = item[0]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio) / sliding_window_size
            if avg_homo_ratio >= homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_start = first_candidate_boundary
                break
        if new_boundary_start != cur_boundary and new_boundary_start != -1:
            if debug:
                print('align start right non-homology, new boundary: ' + str(new_boundary_start))
        cur_boundary = new_boundary_start

        col_index = cur_boundary
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        homo_cols = []
        while valid_col_count < search_len and col_index >= 0:
            # Starting from position 'pos', search for 15 effective columns to the left.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            max_homo_ratio = 0
            no_gap_num = row_num - base_map['-']
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / no_gap_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo, max_homo_ratio))
            col_index -= 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align start left: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the left. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        homo_cols.reverse()
        new_boundary_start = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            first_candidate_boundary = -1
            for item in window:
                cur_homo_ratio = item[5]
                if cur_homo_ratio >= homo_threshold-0.1 and first_candidate_boundary == -1:
                    first_candidate_boundary = item[0]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio)/sliding_window_size
            if avg_homo_ratio >= homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_start = first_candidate_boundary
                break
        if new_boundary_start != cur_boundary and new_boundary_start != -1:
            if debug:
                print('align start left homology, new boundary: ' + str(new_boundary_start))
            cur_boundary = new_boundary_start

        return cur_boundary
    else:
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        col_index = pos
        homo_cols = []
        while valid_col_count < search_len and col_index < col_num:
            # Starting from position 'pos', search for 15 effective columns to the right.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            # If the number of non-empty rows exceeds the threshold, then it is an effective row.
            no_gap_num = row_num - base_map['-']
            max_homo_ratio = 0
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / no_gap_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo, max_homo_ratio))
            col_index += 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align end right: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the right. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        cur_boundary = pos
        homo_cols.reverse()
        new_boundary_end = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            first_candidate_boundary = -1
            for item in window:
                cur_homo_ratio = item[5]
                if cur_homo_ratio >= homo_threshold-0.1 and first_candidate_boundary == -1:
                    first_candidate_boundary = item[0]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio) / sliding_window_size
            if avg_homo_ratio >= homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_end = first_candidate_boundary
                break
        if new_boundary_end != cur_boundary and new_boundary_end != -1:
            if debug:
                print('align end right homology, new boundary: ' + str(new_boundary_end))
            cur_boundary = new_boundary_end

        col_index = cur_boundary
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        homo_cols = []
        while valid_col_count < search_len and col_index >= col_num / 2:
            # Starting from position 'pos', search for 20 effective columns to the left.
            # Determine if the current column is effective.
            is_homo_col = False
            base_map = col_base_map[col_index]
            # If the number of non-empty rows exceeds the threshold, then it is an effective row.
            no_gap_num = row_num - base_map['-']
            max_homo_ratio = 0
            gap_num = base_map['-']
            # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
            if gap_num <= valid_col_threshold:
                valid_col_count += 1
                # Determine if the effective column is homologous.
                for base in base_map.keys():
                    if base == '-':
                        continue
                    cur_homo_ratio = float(base_map[base]) / no_gap_num
                    if cur_homo_ratio > max_homo_ratio:
                        max_homo_ratio = cur_homo_ratio
                    if cur_homo_ratio >= homo_threshold:
                        homo_col_count += 1
                        # Check for consecutive homologous columns.
                        if prev_homo:
                            con_homo += 1
                        is_homo_col = True
                        break
                if not is_homo_col:
                    max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
                    con_homo = 0

                    if prev_non_homo:
                        con_no_homo += 1
                    else:
                        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                        con_no_homo = 0
                    is_no_homo_col = True
                    prev_non_homo = True
                    prev_homo = False
                else:
                    max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
                    prev_homo = True
                    prev_non_homo = False
                    con_no_homo = 0
                    is_no_homo_col = False
                homo_cols.append(
                    (col_index, is_homo_col, con_homo, is_no_homo_col, con_no_homo, max_homo_ratio))
            col_index -= 1
        max_con_homo = con_homo if con_homo > max_con_homo else max_con_homo
        max_con_no_homo = con_no_homo if con_no_homo > max_con_no_homo else max_con_no_homo
        # if debug:
        #     print('align end left: ' + str(homo_col_count) + ', max continous homology bases: ' + str(max_con_homo)
        #           + ', max continous no-homology bases: ' + str(max_con_no_homo))
        #     print(homo_cols)

        # Use a sliding window to calculate the average homology of 10 consecutive bases starting from the right. Determine if it exceeds the threshold.
        # If it exceeds the threshold, obtain the first column with homology above the threshold within the 10bp, and consider it as the homologous boundary.
        new_boundary_end = -1
        for i in range(len(homo_cols) - sliding_window_size + 1):
            window = homo_cols[i:i + sliding_window_size]
            avg_homo_ratio = 0
            first_candidate_boundary = -1
            for item in window:
                cur_homo_ratio = item[5]
                if cur_homo_ratio >= homo_threshold-0.1 and first_candidate_boundary == -1:
                    first_candidate_boundary = item[0]
                avg_homo_ratio += cur_homo_ratio
            avg_homo_ratio = float(avg_homo_ratio) / sliding_window_size
            if avg_homo_ratio >= homo_threshold:
                # If homology in the sliding window exceeds the threshold, find the boundary.
                new_boundary_end = first_candidate_boundary
                break
        if new_boundary_end != cur_boundary and new_boundary_end != -1:
            if debug:
                print('align end left non-homology, new boundary: ' + str(new_boundary_end))
        cur_boundary = new_boundary_end

        return cur_boundary

def judge_boundary_v5(cur_seq, align_file, debug, TE_type, plant, result_type):
    # 1. Based on the 'remove gap' multi-alignment file, locate the position of the original sequence (anchor point).
    #     # Extend 20bp on both sides from the anchor point, extract the effective columns, and determine their homology.
    #     If it contradicts our rule, it is a false positive sequence.
    #     # --First, locate the TIR boundary position of the first sequence in the alignment file as the anchor point.
    #     # Take the first and last 20bp of the original sequence, and search on the aligned sequence without gaps.
    anchor_len = 20
    first_10bp = cur_seq[0:anchor_len]
    last_10bp = cur_seq[-anchor_len:]
    align_names, align_contigs = read_fasta(align_file)

    align_start = -1
    align_end = -1
    for name in align_names:
        raw_align_seq = align_contigs[name]
        align_seq = ''
        position_reflex = {}
        cur_align_index = 0
        for i, base in enumerate(raw_align_seq):
            if base == '-':
                continue
            else:
                align_seq += base
                position_reflex[cur_align_index] = i
                cur_align_index += 1

        start_dist = 2
        last_dist = 2
        first_matches = find_near_matches(first_10bp, align_seq, max_l_dist=start_dist)
        last_matches = find_near_matches(last_10bp, align_seq, max_l_dist=last_dist)
        last_matches = last_matches[::-1]
        if len(first_matches) > 0 and len(last_matches) > 0:
            align_no_gap_start = first_matches[0].start
            align_no_gap_end = last_matches[0].end - 1
            align_start = position_reflex[align_no_gap_start]
            align_end = position_reflex[align_no_gap_end]
            break
    if debug:
        print(align_file, align_start, align_end)
    if align_start == -1 or align_end == -1:
        if debug:
            print('not found boundary:' + align_file)
        return False, 'nb', '', 0
    align_names, align_contigs = read_fasta(align_file)
    if len(align_names) <= 0:
        if debug:
            print('align file size = 0, ' + align_file)
        return False, '', '', 0

    # 3. Take the full-length sequence to generate a consensus sequence.
    # There should be bases both up and down by 10bp at the anchor point.
    full_length_member_names = []
    full_length_member_contigs = {}
    anchor_len = 10
    for name in align_names:
        # 为了减少计算量，只取100条全长拷贝
        if len(full_length_member_names) > 100:
            break
        align_seq = align_contigs[name]
        if align_start - anchor_len >= 0:
            anchor_start = align_start - anchor_len
        else:
            anchor_start = 0
        anchor_start_seq = align_seq[anchor_start: align_start + anchor_len]
        if align_end + anchor_len < len(align_seq):
            anchor_end = align_end + anchor_len
        else:
            anchor_end = len(align_seq)
        anchor_end_seq = align_seq[align_end - anchor_len: anchor_end]

        if not all(c == '-' for c in list(anchor_start_seq)) and not all(c == '-' for c in list(anchor_end_seq)):
            full_length_member_names.append(name)
            full_length_member_contigs[name] = align_seq

    first_seq = full_length_member_contigs[full_length_member_names[0]]
    col_num = len(first_seq)
    row_num = len(full_length_member_names)
    if row_num <= 1:
        if debug:
            print('full length number = 1, ' + align_file)
        return False, 'fl1', '', 1
    matrix = [[''] * col_num for i in range(row_num)]
    for row, name in enumerate(full_length_member_names):
        seq = full_length_member_contigs[name]
        for col in range(len(seq)):
            matrix[row][col] = seq[col]

    # Starting from column 'align_start', search for 15 effective columns to the left.
    # Count the base composition of each column, in the format of {40: {A: 10, T: 5, C: 7, G: 9, '-': 20}},
    # which indicates the number of different bases in the current column.
    # Based on this, it is easy to determine whether the current column is effective and whether it is a homologous column.
    sliding_window_size = 10
    valid_col_threshold = int(row_num/2)

    if row_num <= 2:
        homo_threshold = 0.95
    elif row_num <= 5:
        homo_threshold = 0.9
    else:
        homo_threshold = 0.7

    ############################################
    # 我现在想先遍历一遍整个矩阵，找到合法边界的起始和终止位置，如果最终识别到的边界和合法边界有任何一边重叠了，这代表着无法找到同源和非同源边界，说明这是一个假阳性。

    # Record the base composition of each column.
    col_base_map = {}
    for col_index in range(col_num):
        if not col_base_map.__contains__(col_index):
            col_base_map[col_index] = {}
        base_map = col_base_map[col_index]
        # Calculate the base composition ratio in the current column.
        if len(base_map) == 0:
            for row in range(row_num):
                cur_base = matrix[row][col_index]
                if not base_map.__contains__(cur_base):
                    base_map[cur_base] = 0
                cur_count = base_map[cur_base]
                cur_count += 1
                base_map[cur_base] = cur_count
        if not base_map.__contains__('-'):
            base_map['-'] = 0

    # Initialize variables
    valid_left_boundary = -1
    valid_right_boundary = -1

    # First, search for the left boundary (valid_left_boundary)
    for left in range(col_num):
        base_map_left = col_base_map[left]
        gap_num_left = base_map_left.get('-', 0)
        copy_count_left = row_num  # Total number of rows is the total copy count
        if gap_num_left <= copy_count_left / 2:
            valid_left_boundary = left
            break  # Exit once the left boundary is found

    # If left boundary is found, search for the right boundary (valid_right_boundary)
    if valid_left_boundary != -1:
        for right in range(col_num - 1, -1, -1):
            base_map_right = col_base_map[right]
            gap_num_right = base_map_right.get('-', 0)
            copy_count_right = row_num  # Total number of rows is the total copy count
            if gap_num_right <= copy_count_right / 2:
                valid_right_boundary = right
                break  # Exit once the right boundary is found

    # Check if the left boundary is less than the right boundary
    if not (valid_left_boundary != -1 and valid_right_boundary != -1 and valid_left_boundary < valid_right_boundary):
        valid_left_boundary = -1
        valid_right_boundary = -1
    ##################################################

    homo_boundary_start = search_boundary_homo_v3(valid_col_threshold, align_start, matrix, row_num,
                                             col_num, 'start', homo_threshold, debug, sliding_window_size)
    if homo_boundary_start == -1:
        return False, '', '', row_num

    homo_boundary_end = search_boundary_homo_v3(valid_col_threshold, align_end, matrix, row_num,
                                           col_num, 'end', homo_threshold, debug, sliding_window_size)

    if homo_boundary_end == -1:
        return False, '', '', row_num


    # Generate a consensus sequence.
    model_seq = ''
    if result_type == 'cons':
        for col_index in range(homo_boundary_start, homo_boundary_end+1):
            base_map = col_base_map[col_index]
            # Identify the most frequent base that exceeds the threshold 'valid_col_threshold'.
            max_base_count = 0
            max_base = ''
            for cur_base in base_map.keys():
                cur_count = base_map[cur_base]
                if cur_count > max_base_count:
                    max_base_count = cur_count
                    max_base = cur_base
            if max_base_count >= int(row_num/2):
                if max_base != '-':
                    model_seq += max_base
                else:
                    continue
            else:
                max_base_count = 0
                max_base = ''
                for cur_base in base_map.keys():
                    if cur_base == '-':
                        continue
                    cur_count = base_map[cur_base]
                    if cur_count > max_base_count:
                        max_base_count = cur_count
                        max_base = cur_base
                model_seq += max_base
    else:
        for col_index in range(homo_boundary_start, homo_boundary_end + 1):
            cur_base = first_seq[col_index]
            if cur_base != '-':
                model_seq += cur_base
            else:
                continue

    # Determine the starting base of the homologous boundary.
    final_boundary_start = -1
    final_boundary_end = -1
    final_cons_seq = ''

    if homo_boundary_start <= valid_left_boundary or homo_boundary_end >= valid_right_boundary:
        final_cons_seq = ''
    else:
        if TE_type == 'tir':
            # (TA, TTA, TAA, TTAA)xxxxx...xxxxx(TA, TTA, TAA, TTAA) may lead to incorrect homologous boundaries.
            first_5bps = []
            end_5bps = []
            first_5bps.append((model_seq[0:5], 0))
            end_5bps.append((model_seq[-5:], 0))

            if model_seq.startswith('A'):
                first_5bps.append((model_seq[1:6], 1))
            if model_seq.startswith('AA') or model_seq.startswith('TA'):
                first_5bps.append((model_seq[2:7], 2))
            if model_seq.startswith('TAA') or model_seq.startswith('TTA'):
                first_5bps.append((model_seq[3:8], 3))
            if model_seq.startswith('TTAA'):
                first_5bps.append((model_seq[4:9], 4))

            if model_seq.endswith('T'):
                end_5bps.append((model_seq[len(model_seq)-6:len(model_seq)-1], 1))
            if model_seq.endswith('TT') or model_seq.endswith('TA'):
                end_5bps.append((model_seq[len(model_seq)-7:len(model_seq)-2], 2))
            if model_seq.endswith('TAA') or model_seq.endswith('TTA'):
                end_5bps.append((model_seq[len(model_seq)-8:len(model_seq)-3], 3))
            if model_seq.endswith('TTAA'):
                end_5bps.append((model_seq[len(model_seq)-9:len(model_seq)-4], 4))

            # Take all valid boundaries, then sort them based on the distance between 'first_5bp' and 'end_5bp' + the number of TSDs.
            # Smaller distances and more TSDs make it more likely to be a true boundary.
            all_boundaries = []
            for first_5bp in first_5bps:
                for end_5bp in end_5bps:
                    # Determine if there are two or more TSDs in the align file based on the boundary.
                    # If yes, it is a genuine boundary.
                    cur_boundary_start = homo_boundary_start + first_5bp[1]
                    cur_boundary_end = homo_boundary_end - end_5bp[1]
                    tsd_count = 0
                    for name in align_names:
                        # 1. Verify if there are bases at the boundary.
                        raw_align_seq = align_contigs[name]
                        boundary_start_base = raw_align_seq[cur_boundary_start]
                        boundary_end_base = raw_align_seq[cur_boundary_end]
                        if boundary_start_base == '-' or boundary_end_base == '-':
                            continue
                        # 2. Can TSDs be found at the boundary?
                        left_tsd_seq, right_tsd_seq = TSDsearch_v5(raw_align_seq, cur_boundary_start, cur_boundary_end, plant)
                        if left_tsd_seq != '':
                            tsd_count += 1
                    if tsd_count > 0:
                        edit_distance = Levenshtein.distance(getReverseSequence(first_5bp[0]), end_5bp[0])
                        all_boundaries.append((edit_distance, tsd_count, cur_boundary_start, cur_boundary_end, first_5bp[1], end_5bp[1]))
            all_boundaries.sort(key=lambda x: (x[0], -x[1]))
            if len(all_boundaries) > 0:
                boundary = all_boundaries[0]
                final_boundary_start = boundary[2]
                final_boundary_end = boundary[3]
                if boundary[5] != 0:
                    final_cons_seq = model_seq[boundary[4]: -boundary[5]]
                else:
                    final_cons_seq = model_seq[boundary[4]:]

    if final_cons_seq == '':
        is_TE = False
    else:
        is_TE = True

    if debug:
        print(align_file, is_TE, final_boundary_start, final_boundary_end)
    return is_TE, '', final_cons_seq, row_num


def judge_boundary_v9(cur_seq, align_file, debug, TE_type, plant, result_type):
    # 1. Based on the 'remove gap' multi-alignment file, locate the position of the original sequence (anchor point).
    #     # Extend 20bp on both sides from the anchor point, extract the effective columns, and determine their homology.
    #     If it contradicts our rule, it is a false positive sequence.
    #     # --First, locate the TIR boundary position of the first sequence in the alignment file as the anchor point.
    #     # Take the first and last 20bp of the original sequence, and search on the aligned sequence without gaps.

    anchor_len = 20
    first_10bp = cur_seq[0:anchor_len]
    last_10bp = cur_seq[-anchor_len:]
    align_names, align_contigs = read_fasta(align_file)
    align_start = -1
    align_end = -1
    for name in align_names:
        raw_align_seq = align_contigs[name]
        align_seq = ''
        position_reflex = {}
        cur_align_index = 0
        for i, base in enumerate(raw_align_seq):
            if base == '-':
                continue
            else:
                align_seq += base
                position_reflex[cur_align_index] = i
                cur_align_index += 1

        start_dist = 2
        last_dist = 2
        first_matches = find_near_matches(first_10bp, align_seq, max_l_dist=start_dist)
        last_matches = find_near_matches(last_10bp, align_seq, max_l_dist=last_dist)
        last_matches = last_matches[::-1]
        if len(first_matches) > 0 and len(last_matches) > 0:
            align_no_gap_start = first_matches[0].start
            align_no_gap_end = last_matches[0].end - 1
            align_start = position_reflex[align_no_gap_start]
            align_end = position_reflex[align_no_gap_end]
            break
    if debug:
        print(align_file, align_start, align_end)
    if align_start == -1 or align_end == -1:
        if debug:
            print('not found boundary:' + align_file)
        return False, 'nb', '', 0

    align_names, align_contigs = read_fasta(align_file)
    if len(align_names) <= 0:
        if debug:
            print('align file size = 0, ' + align_file)
        return False, '', '', 0

    # 3. Take the full-length sequence to generate a consensus sequence.
    # There should be bases both up and down by 10bp at the anchor point.
    full_length_member_names = []
    full_length_member_contigs = {}
    anchor_len = 10
    for name in align_names:
        # 为了减少计算量，只取100条全长拷贝
        if len(full_length_member_names) > 100:
            break
        align_seq = align_contigs[name]
        if align_start - anchor_len >= 0:
            anchor_start = align_start - anchor_len
        else:
            anchor_start = 0
        anchor_start_seq = align_seq[anchor_start: align_start + anchor_len]
        if align_end + anchor_len < len(align_seq):
            anchor_end = align_end + anchor_len
        else:
            anchor_end = len(align_seq)
        anchor_end_seq = align_seq[align_end - anchor_len: anchor_end]

        if not all(c == '-' for c in list(anchor_start_seq)) and not all(c == '-' for c in list(anchor_end_seq)):
            full_length_member_names.append(name)
            full_length_member_contigs[name] = align_seq

    first_seq = full_length_member_contigs[full_length_member_names[0]]
    col_num = len(first_seq)
    row_num = len(full_length_member_names)
    if row_num <= 1:
        if debug:
            print('full length number = 1, ' + align_file)
        return False, 'fl1', '', row_num
    matrix = [[''] * col_num for i in range(row_num)]
    for row, name in enumerate(full_length_member_names):
        seq = full_length_member_contigs[name]
        for col in range(len(seq)):
            matrix[row][col] = seq[col]

    # Starting from column 'align_start', search for 15 effective columns to the left.
    # Count the base composition of each column, in the format of {40: {A: 10, T: 5, C: 7, G: 9, '-': 20}},
    # which indicates the number of different bases in the current column.
    # Based on this, it is easy to determine whether the current column is effective and whether it is a homologous column.
    sliding_window_size = 10
    valid_col_threshold = int(row_num/2)

    if row_num <= 2:
        homo_threshold = 0.95
    elif row_num <= 5:
        homo_threshold = 0.9
    else:
        homo_threshold = 0.8

    homo_boundary_start = search_boundary_homo_v3(valid_col_threshold, align_start, matrix, row_num,
                                             col_num, 'start', homo_threshold, debug, sliding_window_size)
    if homo_boundary_start == -1:
        return False, '', '', row_num

    homo_boundary_end = search_boundary_homo_v3(valid_col_threshold, align_end, matrix, row_num,
                                                col_num, 'end', homo_threshold, debug, sliding_window_size)

    if homo_boundary_end == -1:
        return False, '', '', row_num

    # Iterate through each copy in the multiple sequence alignment, take 15-bp of bases with no gaps above
    # and below the homologous boundary, search for polyA/T within the window, locate the position of polyA/T,
    # and then search for TSDs of 8-bp or more in the upstream 30-bp of the other homologous boundary.

    # 迭代每个拷贝，从尾部搜索第一个polyA, 提取right TSD
    TSD_sizes = list(range(8, 21))
    end_5_window_size = 50
    cur_boundary_start = homo_boundary_start
    tsd_count = 0
    first_non_ltr_seq = ''
    for name in align_names:
        raw_align_seq = align_contigs[name]
        align_seq = ''
        gap_to_nogap = {}
        nogap_to_gap = {}
        cur_align_index = 0
        for i, base in enumerate(raw_align_seq):
            if base == '-':
                gap_to_nogap[i] = cur_align_index
                continue
            else:
                align_seq += base
                nogap_to_gap[cur_align_index] = i
                gap_to_nogap[i] = cur_align_index
                cur_align_index += 1

        end_5 = gap_to_nogap[cur_boundary_start]
        end_3, polyA_seq = find_tail_polyA(align_seq)
        # polyA的边界不能与同源边界相差太多
        homology_end_3 = gap_to_nogap[homo_boundary_end]
        if abs(end_3 - homology_end_3) > 10:
            continue

        found_TSD = False
        TSD_seq = ''
        # After locating the 3' end, attempt to search for a set of TSDs in the side wing (8-20) lengths.
        # Search for the corresponding length of TSD near the 5' end (30 bp), and once found, confirm the final 5' end.
        if end_3 != -1 and end_5 != -1:
            # Obtain all possible TSDs on the side wing of the 3' end.
            TSD_list = [(k, align_seq[end_3:end_3 + k]) for k in TSD_sizes]
            # Search for TSDs of various lengths near the 5' end (30 bp) (when TSD len >=8, allow 1bp mismatch).
            subsequence = align_seq[end_5 - end_5_window_size: end_5]
            for k, TSD in reversed(TSD_list):
                for i in range(0, len(subsequence) - k + 1):
                    kmer = subsequence[i:i + k]
                    dist = 1
                    if k == len(TSD) and k == len(kmer):
                        first_matches = find_near_matches(TSD, kmer, max_l_dist=dist)
                        if len(first_matches) > 0:
                            end_5 = end_5 - end_5_window_size + i + k
                            found_TSD = True
                            TSD_seq = TSD
                            break
                if found_TSD:
                    break
        if found_TSD:
            tsd_count += 1
            if first_non_ltr_seq == '':
                final_boundary_start = min(end_5, end_3)
                final_boundary_end = max(end_5, end_3)
                first_non_ltr_seq = align_seq[final_boundary_start: final_boundary_end]
                homo_boundary_start = nogap_to_gap[final_boundary_start]

    # Generate a consensus sequence.
    model_seq = ''
    if tsd_count >= 5 or tsd_count > row_num / 2:
        # Record the base composition of each column.
        col_base_map = {}
        for col_index in range(col_num):
            if not col_base_map.__contains__(col_index):
                col_base_map[col_index] = {}
            base_map = col_base_map[col_index]
            # Calculate the base composition ratio in the current column.
            if len(base_map) == 0:
                for row in range(row_num):
                    cur_base = matrix[row][col_index]
                    if not base_map.__contains__(cur_base):
                        base_map[cur_base] = 0
                    cur_count = base_map[cur_base]
                    cur_count += 1
                    base_map[cur_base] = cur_count
            if not base_map.__contains__('-'):
                base_map['-'] = 0
        for col_index in range(homo_boundary_start, homo_boundary_end+1):
            base_map = col_base_map[col_index]
            # Identify the most frequent base that exceeds the threshold 'valid_col_threshold'.
            max_base_count = 0
            max_base = ''
            for cur_base in base_map.keys():
                cur_count = base_map[cur_base]
                if cur_count > max_base_count:
                    max_base_count = cur_count
                    max_base = cur_base
            if max_base_count >= int(row_num/2):
                if max_base != '-':
                    model_seq += max_base
                else:
                    continue
            else:
                max_base_count = 0
                max_base = ''
                for cur_base in base_map.keys():
                    if cur_base == '-':
                        continue
                    cur_count = base_map[cur_base]
                    if cur_count > max_base_count:
                        max_base_count = cur_count
                        max_base = cur_base
                model_seq += max_base

    if model_seq == '' or len(model_seq) < 80:
        is_TE = False
    else:
        is_TE = True

    if debug:
        print(align_file, is_TE, homo_boundary_start, homo_boundary_end)
    return is_TE, '', model_seq, row_num

def most_common_element(arr):
    if len(arr) > 0:
        counter = Counter(arr)
        most_common = counter.most_common(1)
        if most_common[0][1] == 1:
            return arr[0]
        return most_common[0][0]
    else:
        return -1

def judge_boundary_v6(cur_seq, align_file, debug, TE_type, plant, result_type):
    # 1. Based on the 'remove gap' multi-alignment file, locate the position of the original sequence (anchor point).
    #     # Extend 20bp on both sides from the anchor point, extract the effective columns, and determine their homology.
    #     If it contradicts our rule, it is a false positive sequence.
    #     # --First, locate the TIR boundary position of the first sequence in the alignment file as the anchor point.
    #     # Take the first and last 20bp of the original sequence, and search on the aligned sequence without gaps.
    anchor_len = 20
    first_10bp = cur_seq[0:anchor_len]
    last_10bp = cur_seq[-anchor_len:]
    align_names, align_contigs = read_fasta(align_file)
    # Process the aligned sequence into a no-gap form: align_seq
    # position_reflex stores the mapping relationship between the processed and original positions.
    align_starts = []
    align_ends = []
    align_start = -1
    align_end = -1
    for name in align_names:
        raw_align_seq = align_contigs[name]
        align_seq = ''
        position_reflex = {}
        cur_align_index = 0
        for i, base in enumerate(raw_align_seq):
            if base == '-':
                continue
            else:
                align_seq += base
                position_reflex[cur_align_index] = i
                cur_align_index += 1
        # Locate the anchor point in the aligned sequence based on the original sequence.
        start_dist = 2
        last_dist = 2
        first_matches = find_near_matches(first_10bp, align_seq, max_l_dist=start_dist)
        last_matches = find_near_matches(last_10bp, align_seq, max_l_dist=last_dist)
        last_matches = last_matches[::-1]
        if len(first_matches) > 0 and len(last_matches) > 0:
            align_no_gap_start = first_matches[0].start
            align_no_gap_end = last_matches[0].end - 1
            # The position needs to be mapped back to the original aligned sequence.
            align_start = position_reflex[align_no_gap_start]
            align_end = position_reflex[align_no_gap_end]
            align_starts.append(align_start)
            align_ends.append(align_end)
            #break
    # Select the most frequently occurring position as the anchor point.
    align_start = most_common_element(align_starts)
    align_end = most_common_element(align_ends)

    if debug:
        print(align_file, align_start, align_end)
    if align_start == -1 or align_end == -1:
        if debug:
            print('not found boundary:' + align_file)
        return False, 'nb', '', 0
    align_names, align_contigs = read_fasta(align_file)
    if len(align_names) <= 0:
        if debug:
            print('align file size = 0, ' + align_file)
        return False, '', '', 0

    tail_motifs = ['CTAGT', 'CTAAT', 'CTGGT', 'CTGAT']
    # 2. For align_start and align_end, take all copies that are not empty at the boundaries.
    # Determine if the copy has high homology within the boundary and non-homology outside the boundary.
    # There should be bases both up and down by 2bp at the anchor point.
    start_member_names = []
    start_member_contigs = {}
    end_member_names = []
    end_member_contigs = {}
    search_len = 1
    for name in align_names:
        if len(start_member_contigs) > 100:
            break
        if len(end_member_contigs) > 100:
            break
        align_seq = align_contigs[name]
        if align_start - search_len >= 0:
            anchor_start = align_start - search_len
        else:
            anchor_start = 0
        anchor_start_seq = align_seq[anchor_start: align_start + search_len]
        if not all(c == '-' for c in list(anchor_start_seq)):
            start_member_names.append(name)
            start_member_contigs[name] = align_seq

        # Check if align_seq ends with CTRRT, and filter out those that do not.
        # Take the downstream 1bp and upstream 3bp of align_end, a total of 5bp.
        # Attempt to extend to the left by 3bp.
        tail_seq = align_seq[align_end]
        ext_len = 3
        col_index = align_end - 1
        ext_count = 0
        while ext_count < ext_len and col_index >= 0:
            cur_base = align_seq[col_index]
            if cur_base != '-':
                tail_seq = cur_base + tail_seq
                ext_count += 1
            col_index -= 1

        # Attempt to extend to the right by 1bp.
        ext_len = 1
        col_index = align_end + 1
        ext_count = 0
        while ext_count < ext_len and col_index < len(align_seq):
            cur_base = align_seq[col_index]
            if cur_base != '-':
                tail_seq += cur_base
                ext_count += 1
            col_index += 1

        if align_end + search_len < len(align_seq):
            anchor_end = align_end + search_len
        else:
            anchor_end = len(align_seq)
        anchor_end_seq = align_seq[align_end - search_len: anchor_end]
        if not all(c == '-' for c in list(anchor_end_seq)):
            end_member_names.append(name)
            end_member_contigs[name] = align_seq

    # Starting from the anchor point, extend 20bp on both sides and extract the effective columns. Determine their homology. If it contradicts our rule, it is a false positive sequence.
    # Store the sequence as a matrix for easy traversal and indexing.
    if len(end_member_names) <= 0:
        return False, '', '', 0
    first_seq = end_member_contigs[end_member_names[0]]
    col_num = len(first_seq)
    end_row_num = len(end_member_names)
    matrix = [[''] * col_num for i in range(end_row_num)]
    for row, name in enumerate(end_member_names):
        seq = end_member_contigs[name]
        for col in range(len(seq)):
            matrix[row][col] = seq[col]

    # Since Helitron copies are relatively rare, it is difficult to identify homologous boundaries using the same method as TIR.
    # Therefore, we use a method based on continuous homology to search.
    # Since the end of Helitron is fixed, we start searching from the end.
    # If the conditions are met and there is no change, it means the end has been correctly identified.
    sliding_window_size = 10
    valid_col_threshold = int(end_row_num/2)

    if end_row_num <= 2:
        homo_threshold = 0.95
        int_homo_threshold = 0.9
        out_homo_threshold = 0.95
    elif end_row_num <= 5:
        homo_threshold = 0.9
        int_homo_threshold = 0.85
        out_homo_threshold = 0.9
    else:
        homo_threshold = 0.7
        int_homo_threshold = 0.65
        out_homo_threshold = 0.7

    # First, check if the end is a genuine Helitron support.
    end_align_valid, homo_boundary_end = search_boundary_homo_v4(valid_col_threshold, align_end, matrix, end_row_num,
                                                 col_num, 'end', homo_threshold, int_homo_threshold, out_homo_threshold, debug, sliding_window_size)
    if not end_align_valid:
        if debug:
            print(align_file, end_align_valid)
        return end_align_valid, '', '', 0

    # Store the sequence as a matrix for easy traversal and indexing.
    first_seq = start_member_contigs[start_member_names[0]]
    col_num = len(first_seq)
    start_row_num = len(start_member_names)
    matrix = [[''] * col_num for i in range(start_row_num)]
    for row, name in enumerate(start_member_names):
        seq = start_member_contigs[name]
        for col in range(len(seq)):
            matrix[row][col] = seq[col]

    valid_col_threshold = int(start_row_num / 2)
    if start_row_num <= 2:
        homo_threshold = 0.95
        int_homo_threshold = 0.9
        out_homo_threshold = 0.95
    elif start_row_num <= 5:
        homo_threshold = 0.9
        int_homo_threshold = 0.85
        out_homo_threshold = 0.9
    else:
        homo_threshold = 0.7
        int_homo_threshold = 0.65
        out_homo_threshold = 0.7
    # Start searching for homologous column boundaries from the beginning.
    start_align_valid, homo_boundary_start = search_boundary_homo_v4(valid_col_threshold, align_start, matrix, start_row_num,
                                                                 col_num, 'start', homo_threshold, int_homo_threshold,
                                                                 out_homo_threshold, debug, sliding_window_size)

    # 3. Take the full-length sequence to generate a consensus sequence.
    # There should be bases both up and down by 1bp at the anchor point.
    full_length_member_names = []
    full_length_member_contigs = {}
    anchor_len = 1
    for name in align_names:
        # To reduce computational load, only take 100 full-length copies.
        if len(full_length_member_names) > 100:
            break
        align_seq = align_contigs[name]
        if homo_boundary_start - anchor_len >= 0:
            anchor_start = homo_boundary_start - anchor_len
        else:
            anchor_start = 0
        anchor_start_seq = align_seq[anchor_start: homo_boundary_start + anchor_len]
        if homo_boundary_end + anchor_len < len(align_seq):
            anchor_end = homo_boundary_end + anchor_len
        else:
            anchor_end = len(align_seq)
        anchor_end_seq = align_seq[homo_boundary_end - anchor_len: anchor_end]

        if not all(c == '-' for c in list(anchor_start_seq)) and not all(c == '-' for c in list(anchor_end_seq)):
            full_length_member_names.append(name)
            full_length_member_contigs[name] = align_seq

    if len(full_length_member_names) <= 0:
        if debug:
            print('full length align file size = 0, ' + align_file)
        return False, '', '', 0

    # Generate a consensus sequence using full-length copies.
    first_seq = full_length_member_contigs[full_length_member_names[0]]
    col_num = len(first_seq)
    row_num = len(full_length_member_names)
    matrix = [[''] * col_num for i in range(row_num)]
    for row, name in enumerate(full_length_member_names):
        seq = full_length_member_contigs[name]
        for col in range(len(seq)):
            matrix[row][col] = seq[col]

    # Record the base composition of each column.
    col_base_map = {}
    for col_index in range(col_num):
        if not col_base_map.__contains__(col_index):
            col_base_map[col_index] = {}
        base_map = col_base_map[col_index]
        # Calculate the base composition ratio in the current column.
        if len(base_map) == 0:
            for row in range(row_num):
                cur_base = matrix[row][col_index]
                if not base_map.__contains__(cur_base):
                    base_map[cur_base] = 0
                cur_count = base_map[cur_base]
                cur_count += 1
                base_map[cur_base] = cur_count
        if not base_map.__contains__('-'):
            base_map['-'] = 0

    # Generate a consensus sequence.
    model_seq = ''
    if result_type == 'cons':
        for col_index in range(homo_boundary_start, homo_boundary_end+1):
            base_map = col_base_map[col_index]
            # Identify the most frequent base that exceeds the threshold 'valid_col_threshold'.
            max_base_count = 0
            max_base = ''
            for cur_base in base_map.keys():
                cur_count = base_map[cur_base]
                if cur_count > max_base_count:
                    max_base_count = cur_count
                    max_base = cur_base
            if max_base_count >= int(row_num/2):
                if max_base != '-':
                    model_seq += max_base
                else:
                    continue
            else:
                model_seq += 'N'
    else:
        for col_index in range(homo_boundary_start, homo_boundary_end + 1):
            cur_base = first_seq[col_index]
            if cur_base != '-':
                model_seq += cur_base
            else:
                continue

    # Determine the starting base of the homologous boundary.
    # After determining homo_boundary_start at the starting end,
    # extend the model seq by 1bp and search for the position of 'ATC' in the model seq to determine the boundary.
    final_boundary_start = -1
    final_boundary_end = homo_boundary_end
    final_cons_seq = ''
    if TE_type == 'helitron':
        # For helitron transposons, you need to look for ATC...CTRR (R=A/G) T near the boundary.
        # Strategy 1: Extend 1bp outside the boundary.
        # If ATC or CTRRT can be found within 10bp of the boundary, it is considered a true Helitron.
        model_seq = ''
        for col_index in range(homo_boundary_start, homo_boundary_end + 1):
            base_map = col_base_map[col_index]
            # Identify the most frequent base that exceeds the threshold 'valid_col_threshold'.
            max_base_count = 0
            max_base = ''
            for cur_base in base_map.keys():
                cur_count = base_map[cur_base]
                if cur_count > max_base_count:
                    max_base_count = cur_count
                    max_base = cur_base
            if max_base_count >= int(row_num / 2):
                if max_base != '-':
                    model_seq += max_base
                else:
                    continue
            else:
                model_seq += 'N'
        # Attempt to extend by 1bp on both ends.
        ext_len = 1
        col_index = homo_boundary_start - 1
        ext_count = 0
        while ext_count < ext_len and col_index >= 0:
            base_map = col_base_map[col_index]
            # Identify the most frequent base that exceeds the threshold 'valid_col_threshold'.
            max_base_count = 0
            max_base = ''
            for cur_base in base_map.keys():
                cur_count = base_map[cur_base]
                if cur_count > max_base_count:
                    max_base_count = cur_count
                    max_base = cur_base
            if max_base_count >= int(row_num / 2) and max_base != '-':
                model_seq = max_base + model_seq
                ext_count += 1
            col_index -= 1

        col_index = homo_boundary_end + 1
        ext_count = 0
        while ext_count < ext_len and col_index < col_num:
            base_map = col_base_map[col_index]
            # Identify the most frequent base that exceeds the threshold 'valid_col_threshold'.
            max_base_count = 0
            max_base = ''
            for cur_base in base_map.keys():
                cur_count = base_map[cur_base]
                if cur_count > max_base_count:
                    max_base_count = cur_count
                    max_base = cur_base
            if max_base_count >= int(row_num / 2) and max_base != '-':
                model_seq += max_base
                ext_count += 1
            col_index += 1

        # Search for the position of 'ATC' in the model seq to determine the starting end boundary.
        # Find the last occurrence of a substring.
        # In the first 10bp and last 10bp, look for the first occurrence of aTC, CTRRt.
        cur_search_len = 10
        first_10bp = model_seq[0: cur_search_len]
        last_10bp = model_seq[-cur_search_len:]
        for tail_motif in tail_motifs:
            # Find the last occurrence of a substring.
            end_index = last_10bp.rfind(tail_motif)
            if end_index != -1:
                start_index = first_10bp.find('ATC')
                if start_index != -1:
                    final_boundary_start = homo_boundary_start - 1 + start_index + 1
                    final_boundary_end = homo_boundary_end + 1 - (cur_search_len - (end_index + 3) - 1)
                    if (cur_search_len - (end_index + 3) - 1) == 0:
                        final_cons_seq = model_seq[start_index + 1:]
                    else:
                        final_cons_seq = model_seq[start_index + 1:-(cur_search_len - (end_index + 3) - 1)]
                    break

    if final_cons_seq == '':
        is_TE = False
    else:
        is_TE = True

    if debug:
        print(align_file, is_TE, final_boundary_start, final_boundary_end)
    return is_TE, '', final_cons_seq, row_num

def judge_boundary_v7(cur_seq, align_file, debug, TE_type, plant, result_type):
    # 1. Based on the 'remove gap' multi-alignment file, locate the position of the original sequence (anchor point).
    #     # Extend 20bp on both sides from the anchor point, extract the effective columns, and determine their homology.
    #     If it contradicts our rule, it is a false positive sequence.
    #     # --First, locate the TIR boundary position of the first sequence in the alignment file as the anchor point.
    #     # Take the first and last 20bp of the original sequence, and search on the aligned sequence without gaps.
    anchor_len = 20
    first_10bp = cur_seq[0:anchor_len]
    last_10bp = cur_seq[-anchor_len:]
    align_names, align_contigs = read_fasta(align_file)
    # Process the aligned sequence into a no-gap form: align_seq
    # position_reflex stores the mapping relationship between the processed and original positions.
    align_starts = []
    align_ends = []
    align_start = -1
    align_end = -1
    for name in align_names:
        raw_align_seq = align_contigs[name]
        align_seq = ''
        position_reflex = {}
        cur_align_index = 0
        for i, base in enumerate(raw_align_seq):
            if base == '-':
                continue
            else:
                align_seq += base
                position_reflex[cur_align_index] = i
                cur_align_index += 1
        # Locate the anchor point in the aligned sequence based on the original sequence.
        start_dist = 2
        last_dist = 2
        first_matches = find_near_matches(first_10bp, align_seq, max_l_dist=start_dist)
        last_matches = find_near_matches(last_10bp, align_seq, max_l_dist=last_dist)
        last_matches = last_matches[::-1]
        if len(first_matches) > 0 and len(last_matches) > 0:
            align_no_gap_start = first_matches[0].start
            align_no_gap_end = last_matches[0].end - 1
            # The position needs to be mapped back to the original aligned sequence.
            align_start = position_reflex[align_no_gap_start]
            align_end = position_reflex[align_no_gap_end]
            align_starts.append(align_start)
            align_ends.append(align_end)
            #break
    # Select the most frequently occurring position as the anchor point.
    align_start = most_common_element(align_starts)
    align_end = most_common_element(align_ends)

    if debug:
        print(align_file, align_start, align_end)
    if align_start == -1 or align_end == -1:
        if debug:
            print('not found boundary:' + align_file)
        return False, 'nb', ''
    align_names, align_contigs = read_fasta(align_file)
    if len(align_names) <= 1:
        if debug:
            print('align file size <= 1, ' + align_file)
        return False, '', ''

    # 2. For align_start and align_end, take all copies that are not empty at the boundaries. Determine if the copy has high homology within the boundary and non-homology outside the boundary.
    # There should be bases both up and down by 2bp at the anchor point.
    start_member_names = []
    start_member_contigs = {}
    end_member_names = []
    end_member_contigs = {}
    search_len = 1
    for name in align_names:
        if len(start_member_contigs) > 100:
            break
        if len(end_member_contigs) > 100:
            break
        align_seq = align_contigs[name]
        if align_start - search_len >= 0:
            anchor_start = align_start - search_len
        else:
            anchor_start = 0
        anchor_start_seq = align_seq[anchor_start: align_start + search_len]
        if not all(c == '-' for c in list(anchor_start_seq)):
            start_member_names.append(name)
            start_member_contigs[name] = align_seq

        # Check if align_seq ends with CTRRT, and filter out those that do not.
        # Take the downstream 1bp and upstream 3bp of align_end, a total of 5bp.
        # Attempt to extend to the left by 3bp.
        tail_seq = align_seq[align_end]
        ext_len = 3
        col_index = align_end - 1
        ext_count = 0
        while ext_count < ext_len and col_index >= 0:
            cur_base = align_seq[col_index]
            if cur_base != '-':
                tail_seq = cur_base + tail_seq
                ext_count += 1
            col_index -= 1

        # Attempt to extend to the right by 1bp.
        ext_len = 1
        col_index = align_end + 1
        ext_count = 0
        while ext_count < ext_len and col_index < len(align_seq):
            cur_base = align_seq[col_index]
            if cur_base != '-':
                tail_seq += cur_base
                ext_count += 1
            col_index += 1

        if align_end + search_len < len(align_seq):
            anchor_end = align_end + search_len
        else:
            anchor_end = len(align_seq)
        anchor_end_seq = align_seq[align_end - search_len: anchor_end]
        if not all(c == '-' for c in list(anchor_end_seq)):
            end_member_names.append(name)
            end_member_contigs[name] = align_seq

    # Starting from the anchor point, extend 20bp on both sides and extract the effective columns. Determine their homology. If it contradicts our rule, it is a false positive sequence.
    # Store the sequence as a matrix for easy traversal and indexing.
    if len(end_member_names) <= 0:
        return False, '', ''
    first_seq = end_member_contigs[end_member_names[0]]
    col_num = len(first_seq)
    end_row_num = len(end_member_names)
    matrix = [[''] * col_num for i in range(end_row_num)]
    for row, name in enumerate(end_member_names):
        seq = end_member_contigs[name]
        for col in range(len(seq)):
            matrix[row][col] = seq[col]

    sliding_window_size = 10
    valid_col_threshold = int(end_row_num/2)

    if end_row_num <= 2:
        homo_threshold = 0.95
        int_homo_threshold = 0.9
        out_homo_threshold = 0.95
    elif end_row_num <= 5:
        homo_threshold = 0.9
        int_homo_threshold = 0.85
        out_homo_threshold = 0.9
    else:
        homo_threshold = 0.7
        int_homo_threshold = 0.65
        out_homo_threshold = 0.7

    end_align_valid, homo_boundary_end = search_boundary_homo_v5(valid_col_threshold, align_end, matrix, end_row_num,
                                                 col_num, 'end', homo_threshold, int_homo_threshold, out_homo_threshold, debug, sliding_window_size)
    if not end_align_valid:
        if debug:
            print(align_file, end_align_valid)
        return end_align_valid, '', ''

    # Store the sequence as a matrix for easy traversal and indexing.
    first_seq = start_member_contigs[start_member_names[0]]
    col_num = len(first_seq)
    start_row_num = len(start_member_names)
    matrix = [[''] * col_num for i in range(start_row_num)]
    for row, name in enumerate(start_member_names):
        seq = start_member_contigs[name]
        for col in range(len(seq)):
            matrix[row][col] = seq[col]

    valid_col_threshold = int(start_row_num / 2)
    if start_row_num <= 2:
        homo_threshold = 0.95
        int_homo_threshold = 0.9
        out_homo_threshold = 0.95
    elif start_row_num <= 5:
        homo_threshold = 0.9
        int_homo_threshold = 0.85
        out_homo_threshold = 0.9
    else:
        homo_threshold = 0.7
        int_homo_threshold = 0.65
        out_homo_threshold = 0.7
    # Start searching for homologous column boundaries from the beginning.
    start_align_valid, homo_boundary_start = search_boundary_homo_v5(valid_col_threshold, align_start, matrix, start_row_num,
                                                                 col_num, 'start', homo_threshold, int_homo_threshold,
                                                                 out_homo_threshold, debug, sliding_window_size)


    # Iterate through each copy in the multiple sequence alignment, take 15-bp of bases with no gaps above
    # and below the homologous boundary, search for polyA/T within the window, locate the position of polyA/T,
    # and then search for TSDs of 8-bp or more in the upstream 30-bp of the other homologous boundary.
    final_boundary_start = -1
    final_boundary_end = -1
    TSD_sizes = list(range(8, 31))
    end_5_window_size = 30
    cur_boundary_start = homo_boundary_start
    cur_boundary_end = homo_boundary_end
    tsd_count = 0
    first_non_ltr_seq = ''
    for name in align_names:
        raw_align_seq = align_contigs[name]
        align_seq = ''
        gap_to_nogap = {}
        nogap_to_gap = {}
        cur_align_index = 0
        for i, base in enumerate(raw_align_seq):
            if base == '-':
                gap_to_nogap[i] = cur_align_index
                continue
            else:
                align_seq += base
                nogap_to_gap[cur_align_index] = i
                gap_to_nogap[i] = cur_align_index
                cur_align_index += 1

        end_3 = -1
        end_5 = -1
        direct = ''
        prev_tail_len = 0
        # Obtain the polyA sequence near cur_boundary_start.
        start_pos, end_pos, polyA_seq = find_nearest_polyA(align_seq, gap_to_nogap[cur_boundary_start], window_size=15, min_length=5)
        if len(polyA_seq) > 0 and len(polyA_seq) > prev_tail_len:
            end_3 = start_pos
            end_5 = gap_to_nogap[cur_boundary_end]
            direct = '-'
            prev_tail_len = len(polyA_seq)
        # Obtain the polyT sequence near cur_boundary_start.
        start_pos, end_pos, polyT_seq = find_nearest_polyT(align_seq, gap_to_nogap[cur_boundary_start], window_size=15, min_length=5)
        if len(polyT_seq) > 0 and len(polyT_seq) > prev_tail_len:
            end_3 = start_pos
            end_5 = gap_to_nogap[cur_boundary_end]
            direct = '-'
            prev_tail_len = len(polyT_seq)
        # Obtain the polyA sequence near cur_boundary_end.
        start_pos, end_pos, polyA_seq = find_nearest_polyA(align_seq, gap_to_nogap[cur_boundary_end], window_size=15, min_length=5)
        if len(polyA_seq) > 0 and len(polyA_seq) > prev_tail_len:
            end_3 = end_pos
            end_5 = gap_to_nogap[cur_boundary_start]
            direct = '+'
            prev_tail_len = len(polyA_seq)
        # Obtain the polyT sequence near cur_boundary_end.
        start_pos, end_pos, polyT_seq = find_nearest_polyT(align_seq, gap_to_nogap[cur_boundary_end], window_size=15, min_length=5)
        if len(polyT_seq) > 0 and len(polyT_seq) > prev_tail_len:
            end_3 = end_pos
            end_5 = gap_to_nogap[cur_boundary_start]
            direct = '+'
            prev_tail_len = len(polyT_seq)

        found_TSD = False
        TSD_seq = ''
        # After locating the 3' end, attempt to search for a set of TSDs in the side wing (8-20) lengths.
        # Search for the corresponding length of TSD near the 5' end (30 bp), and once found, confirm the final 5' end.
        if end_3 != -1 and end_5 != -1 and direct is not None:
            if direct == '-':
                # Obtain all possible TSDs on the side wing of the 3' end.
                TSD_list = [(k, align_seq[end_3 - k:end_3]) for k in TSD_sizes]
                # Search for TSDs of various lengths near the 5' end (30 bp) (when TSD len >=8, allow 1bp mismatch).
                subsequence = align_seq[end_5:end_5 + end_5_window_size]
                for k, TSD in reversed(TSD_list):
                    for i in range(0, len(subsequence) - k + 1):
                        kmer = subsequence[i:i + k]
                        dist = 1
                        if k == len(TSD) and k == len(kmer):
                            first_matches = find_near_matches(TSD, kmer, max_l_dist=dist)
                            if len(first_matches) > 0:
                                end_5 = end_5 + i
                                # end_5 cannot be too far from the homologous boundary, < 5 bp.
                                if abs(end_5 - gap_to_nogap[cur_boundary_end]) <= 5:
                                    found_TSD = True
                                    TSD_seq = TSD
                                    break
                    if found_TSD:
                        break
            if direct == '+':
                # Obtain all possible TSDs on the side wing of the 3' end.
                TSD_list = [(k, align_seq[end_3:end_3 + k]) for k in TSD_sizes]
                # Search for TSDs of various lengths near the 5' end (30 bp) (when TSD len >=8, allow 1bp mismatch).
                subsequence = align_seq[end_5 - end_5_window_size:end_5]
                for k, TSD in reversed(TSD_list):
                    for i in range(0, len(subsequence) - k + 1):
                        kmer = subsequence[i:i + k]
                        dist = 1
                        if k == len(TSD) and k == len(kmer):
                            first_matches = find_near_matches(TSD, kmer, max_l_dist=dist)
                            if len(first_matches) > 0:
                                end_5 = end_5 - end_5_window_size + i + k
                                # end_5 cannot be too far from the homologous boundary, < 5 bp.
                                if abs(end_5 - gap_to_nogap[cur_boundary_start]) <= 5:
                                    found_TSD = True
                                    TSD_seq = TSD
                                    break
                    if found_TSD:
                        break
        if found_TSD:
            tsd_count += 1
            if first_non_ltr_seq == '':
                final_boundary_start = min(end_5, end_3)
                final_boundary_end = max(end_5, end_3)
                first_non_ltr_seq = align_seq[final_boundary_start: final_boundary_end]
                final_boundary_start = nogap_to_gap[final_boundary_start]
                final_boundary_end = nogap_to_gap[final_boundary_end]

    final_non_ltr_seq = ''
    if tsd_count > 1:
        final_non_ltr_seq = first_non_ltr_seq

    if final_non_ltr_seq == '':
        is_TE = False
    else:
        is_TE = True

    if debug:
        print(align_file, is_TE, final_boundary_start, final_boundary_end, len(align_names))
    return is_TE, '', final_non_ltr_seq

def judge_boundary_v8(align_file, debug):
    # 由于non-LTR没有清晰的结构特征，我们只能根据序列的重复特性进行识别。
    # 因此我们根据序列的多拷贝，生成一致性序列，如果一致性序列具有polyA tail，就认为是non-ltr。要求：生成一致性的要求严格一点，如某列必须 90% 相同的碱基。
    cons_threshold = 0.9
    align_names, align_contigs = read_fasta(align_file)
    first_seq = align_contigs[align_names[0]]
    col_num = len(first_seq)
    row_num = len(align_contigs)
    if row_num <= 1:
        if debug:
            print('full length number = 1, ' + align_file)
        return False, 'fl1', ''
    matrix = [[''] * col_num for i in range(row_num)]
    for row, name in enumerate(align_names):
        seq = align_contigs[name]
        for col in range(len(seq)):
            matrix[row][col] = seq[col]

    # Record the base composition of each column.
    col_base_map = {}
    for col_index in range(col_num):
        if not col_base_map.__contains__(col_index):
            col_base_map[col_index] = {}
        base_map = col_base_map[col_index]
        # Calculate the base composition ratio in the current column.
        if len(base_map) == 0:
            for row in range(row_num):
                cur_base = matrix[row][col_index]
                if not base_map.__contains__(cur_base):
                    base_map[cur_base] = 0
                cur_count = base_map[cur_base]
                cur_count += 1
                base_map[cur_base] = cur_count
        if not base_map.__contains__('-'):
            base_map['-'] = 0
    # Generate a consensus sequence.
    model_seq = ''
    for col_index in range(0, col_num):
        base_map = col_base_map[col_index]
        # Identify the most frequent base that exceeds the threshold 'valid_col_threshold'.
        max_base_count = 0
        max_base = ''
        non_empty_row_num = 0
        for cur_base in base_map.keys():
            if cur_base == '-':
                continue
            cur_count = base_map[cur_base]
            if cur_count > max_base_count:
                max_base_count = cur_count
                max_base = cur_base
            non_empty_row_num += cur_count
        if max_base_count >= 10 and max_base_count >= int(non_empty_row_num * cons_threshold):
            model_seq += max_base

    # 如果一致性序列以polyA结尾就认为是non-ltr
    if model_seq.endswith('AAAAAA'):
        return True, '', model_seq
    else:
        return False, '', model_seq


def get_full_length_member(query_path, reference, flanking_len):
    member_file = query_path + '.blast.bed.fa'
    blastn2Results_path = query_path + '.blast.out'
    align_command = 'blastn -db ' + reference + ' -query ' + query_path + ' -evalue 1e-20 -outfmt 6 > ' + blastn2Results_path
    os.system(align_command)
    all_copies = get_copies_v1(blastn2Results_path, query_path, reference, query_coverage=0.95)

    ref_names, ref_contigs = read_fasta(reference)
    new_all_copies = {}
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        for copy in copies:
            ref_name = copy[0]
            copy_ref_start = int(copy[1])
            copy_ref_end = int(copy[2])
            direct = copy[4]
            copy_len = copy_ref_end - copy_ref_start + 1
            if copy_ref_start - 1 - flanking_len < 0 or copy_ref_end + flanking_len > len(ref_contigs[ref_name]):
                continue
            copy_seq = ref_contigs[ref_name][copy_ref_start - 1 - flanking_len: copy_ref_end + flanking_len]
            if direct == '-':
                copy_seq = getReverseSequence(copy_seq)
            if len(copy_seq) < 100:
                continue
            new_name = ref_name + ':' + str(copy_ref_start) + '-' + str(copy_ref_end) + '(' + direct + ')'
            new_all_copies[new_name] = copy_seq
    store_fasta(new_all_copies, member_file)
    return member_file


def split_chromosomes(chromosomes_dict, max_length=200_000_000):
    """
    分割染色体序列，如果序列长度超过max_length，则将其分割成多个部分。

    参数:
    chromosomes_dict (dict): 一个字典，键为染色体名称，值为对应的DNA序列。
    max_length (int): 最大序列长度，超过此长度的序列将被分割。默认值为200 MB。

    返回:
    dict: 一个新的字典，包含分割后的染色体序列。
    """
    new_chromosomes_dict = {}
    for chrom, sequence in chromosomes_dict.items():
        if len(sequence) > max_length:
            num_parts = (len(sequence) + max_length - 1) // max_length  # 计算需要分割的部分数
            for i in range(num_parts):
                part_name = f"{chrom}_part{i + 1}"
                start = i * max_length
                end = min((i + 1) * max_length, len(sequence))
                new_chromosomes_dict[part_name] = sequence[start:end]
        else:
            new_chromosomes_dict[chrom] = sequence
    return new_chromosomes_dict

def split_dict_into_blocks(chromosomes_dict, threads, chunk_size):
    chromosomes_dict = split_chromosomes(chromosomes_dict, max_length=chunk_size)
    total_length = sum(len(seq) for seq in chromosomes_dict.values())
    target_length = total_length // threads

    blocks = []
    current_block = {}
    current_length = 0

    for chrom, seq in chromosomes_dict.items():
        current_block[chrom] = seq
        current_length += len(seq)

        if current_length >= target_length:
            blocks.append(current_block)
            current_block = {}
            current_length = 0

    if current_block:
        blocks.append(current_block)

    return blocks

def get_full_length_copies(query_path, split_ref_dir, debug):
    blastn2Results_path = query_path + '.blast.out'
    repeats_path = (query_path, split_ref_dir, blastn2Results_path)
    all_copies = multiple_alignment_blast_and_get_copies_v1(repeats_path)
    if debug != 1:
        os.remove(blastn2Results_path)
    return all_copies

def get_full_length_member_batch(query_path, reference, ref_contigs, temp_dir, flanking_len):
    batch_member_files = []
    blastn2Results_path = query_path + '.blast.out'
    repeats_path = (query_path, reference, blastn2Results_path)
    all_copies = multiple_alignment_blast_and_get_copies(repeats_path)
    print(all_copies)

    new_all_copies = {}
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        for copy in copies:
            ref_name = copy[0]
            copy_ref_start = int(copy[1])
            copy_ref_end = int(copy[2])
            direct = copy[4]
            copy_len = copy_ref_end - copy_ref_start + 1
            if copy_ref_start - 1 - flanking_len < 0 or copy_ref_end + flanking_len > len(ref_contigs[ref_name]):
                continue
            copy_seq = ref_contigs[ref_name][copy_ref_start - 1 - flanking_len: copy_ref_end + flanking_len]
            if direct == '-':
                copy_seq = getReverseSequence(copy_seq)
            if len(copy_seq) < 100:
                continue
            new_name = ref_name + ':' + str(copy_ref_start) + '-' + str(copy_ref_end) + '(' + direct + ')'
            if not new_all_copies.__contains__(query_name):
                new_all_copies[query_name] = {}
            copy_contigs = new_all_copies[query_name]
            copy_contigs[new_name] = copy_seq
            new_all_copies[query_name] = copy_contigs
    for query_name in new_all_copies.keys():
        copy_contigs = new_all_copies[query_name]
        cur_member_file = temp_dir + '/' + query_name + '.blast.bed.fa'
        store_fasta(copy_contigs, cur_member_file)
        batch_member_files.append((query_name, cur_member_file))
    return batch_member_files

def run_find_members_v8(batch_member_file, temp_dir, subset_script_path, plant, TE_type, debug, result_type):
    (query_name, cur_seq, member_file) = batch_member_file
    member_names, member_contigs = read_fasta(member_file)
    if len(member_names) > 100:
        sub_command = 'cd ' + temp_dir + ' && sh ' + subset_script_path + ' ' + member_file + ' 100 100 ' + ' > /dev/null 2>&1'
        os.system(sub_command)
        member_file += '.rdmSubset.fa'
    if not os.path.exists(member_file):
        return (None, None, '', 0)
    align_file = member_file + '.maf.fa'
    align_command = 'cd ' + temp_dir + ' && mafft --preservecase --quiet --thread 1 ' + member_file + ' > ' + align_file
    os.system(align_command)

    # Due to the characteristics of TIR, Helitron, and non-LTR elements,
    # their homology filtering methods have slight differences.
    is_TE = False
    cons_seq = ''
    info = None
    copy_num = 0
    if TE_type == 'tir':
        is_TE, info, cons_seq, copy_num = judge_boundary_v5(cur_seq, align_file, debug, TE_type, plant, result_type)
    elif TE_type == 'helitron':
        is_TE, info, cons_seq, copy_num = judge_boundary_v6(cur_seq, align_file, debug, TE_type, plant, result_type)
    elif TE_type == 'non_ltr':
        is_TE, info, cons_seq, copy_num = judge_boundary_v9(cur_seq, align_file, debug, TE_type, plant, result_type)

    if is_TE:
        return (query_name, cons_seq, info, copy_num)
    else:
        return (None, None, info, 0)


def FMEA(blastn2Results_path, fixed_extend_base_threshold):
    # parse blastn output, determine the repeat boundary
    # query_records = {query_name: {subject_name: [(q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end)] }}
    query_records = {}
    with open(blastn2Results_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            # print('current line idx: %d' % (idx))
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            identity = float(parts[2])
            alignment_len = int(parts[3])
            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            if query_name == subject_name and q_start == s_start and q_end == s_end:
                continue
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))
    f_r.close()

    skip_gap = fixed_extend_base_threshold
    longest_repeats = {}
    for idx, query_name in enumerate(query_records.keys()):
        subject_dict = query_records[query_name]

        # if there are more than one longest query overlap with the final longest query over 90%,
        # then it probably the true TE
        longest_queries = []
        for subject_name in subject_dict.keys():
            subject_pos = subject_dict[subject_name]

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
                        cur_subject_start = frag[2]
                        cur_query_end = frag[1]
                        prev_subject_end = exist_frag[3]
                        prev_query_end = exist_frag[1]
                        if (cur_subject_start - prev_subject_end < skip_gap and cur_query_end > prev_query_end):
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
                        cur_subject_start = frag[2]
                        cur_query_end = frag[1]
                        prev_subject_end = exist_frag[3]
                        prev_query_end = exist_frag[1]
                        if (prev_subject_end - cur_subject_start < skip_gap and cur_query_end > prev_query_end):
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

                # print('subject pos size: %d' %(len(cur_cluster)))
                # record visited fragments
                visited_frag = {}
                for i in range(len(cur_cluster)):
                    # keep a longest query start from each fragment
                    prev_frag = cur_cluster[i]
                    if visited_frag.__contains__(prev_frag):
                        continue
                    prev_query_start = prev_frag[0]
                    prev_query_end = prev_frag[1]
                    prev_subject_start = prev_frag[2]
                    prev_subject_end = prev_frag[3]
                    prev_query_seq = (min(prev_query_start, prev_query_end), max(prev_query_start, prev_query_end))
                    prev_subject_seq = (
                        min(prev_subject_start, prev_subject_end), max(prev_subject_start, prev_subject_end))
                    prev_query_len = abs(prev_query_end - prev_query_start)
                    prev_subject_len = abs(prev_subject_end - prev_subject_start)
                    cur_longest_query_len = prev_query_len

                    cur_extend_num = 0
                    visited_frag[prev_frag] = 1
                    # try to extend query
                    for j in range(i + 1, len(cur_cluster)):
                        cur_frag = cur_cluster[j]
                        if visited_frag.__contains__(cur_frag):
                            continue
                        cur_query_start = cur_frag[0]
                        cur_query_end = cur_frag[1]
                        cur_subject_start = cur_frag[2]
                        cur_subject_end = cur_frag[3]
                        cur_query_seq = (min(cur_query_start, cur_query_end), max(cur_query_start, cur_query_end))
                        cur_subject_seq = (
                            min(cur_subject_start, cur_subject_end), max(cur_subject_start, cur_subject_end))
                        cur_query_len = abs(cur_query_end - cur_query_start)
                        cur_subject_len = abs(cur_subject_end - cur_subject_start)

                        query_overlap_len = get_overlap_len(cur_query_seq, prev_query_seq)
                        is_same_query = float(query_overlap_len) / cur_query_len >= 0.5 or float(
                            query_overlap_len) / prev_query_len >= 0.5
                        subject_overlap_len = get_overlap_len(prev_subject_seq, cur_subject_seq)
                        is_same_subject = float(subject_overlap_len) / cur_subject_len >= 0.5 or float(
                            subject_overlap_len) / prev_subject_len >= 0.5

                        # could extend
                        # extend right
                        if cur_query_end > prev_query_end:
                            # judge subject direction
                            if prev_subject_start < prev_subject_end and cur_subject_start < cur_subject_end:
                                # +
                                if cur_subject_end > prev_subject_end:
                                    # forward extend
                                    if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                            and cur_subject_start - prev_subject_end < skip_gap:  # \
                                        # and not is_same_query and not is_same_subject:
                                        # update the longest path
                                        prev_query_start = prev_query_start
                                        prev_query_end = cur_query_end
                                        prev_subject_start = prev_subject_start if prev_subject_start < cur_subject_start else cur_subject_start
                                        prev_subject_end = cur_subject_end
                                        cur_longest_query_len = prev_query_end - prev_query_start
                                        cur_extend_num += 1
                                        visited_frag[cur_frag] = 1
                                    elif cur_query_start - prev_query_end >= skip_gap:
                                        break
                            elif prev_subject_start > prev_subject_end and cur_subject_start > cur_subject_end:
                                # reverse
                                if cur_subject_end < prev_subject_end:
                                    # reverse extend
                                    if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                            and prev_subject_end - cur_subject_start < skip_gap:  # \
                                        # and not is_same_query and not is_same_subject:
                                        # update the longest path
                                        prev_query_start = prev_query_start
                                        prev_query_end = cur_query_end
                                        prev_subject_start = prev_subject_start if prev_subject_start > cur_subject_start else cur_subject_start
                                        prev_subject_end = cur_subject_end
                                        cur_longest_query_len = prev_query_end - prev_query_start
                                        cur_extend_num += 1
                                        visited_frag[cur_frag] = 1
                                    elif cur_query_start - prev_query_end >= skip_gap:
                                        break
                    # keep this longest query
                    if cur_longest_query_len != -1:
                        longest_queries.append(
                            (prev_query_start, prev_query_end, cur_longest_query_len, prev_subject_start,
                             prev_subject_end, abs(prev_subject_end - prev_subject_start), subject_name,
                             cur_extend_num))
        if not longest_repeats.__contains__(query_name):
            longest_repeats[query_name] = []
        cur_longest_repeats = longest_repeats[query_name]
        for repeat in longest_queries:
            # Subject序列处理流程
            subject_name = repeat[6]
            old_subject_start_pos = repeat[3] - 1
            old_subject_end_pos = repeat[4]
            # Query序列处理流程
            old_query_start_pos = repeat[0] - 1
            old_query_end_pos = repeat[1]
            cur_query_seq_len = abs(old_query_end_pos - old_query_start_pos)
            cur_longest_repeats.append((query_name, old_query_start_pos, old_query_end_pos, subject_name, old_subject_start_pos, old_subject_end_pos))

    return longest_repeats

def cons_from_mafft(align_file):
    align_names, align_contigs = read_fasta(align_file)
    if len(align_names) <= 0:
        return None

    # Generate a consensus sequence using full-length copies.
    first_seq = align_contigs[align_names[0]]
    col_num = len(first_seq)
    row_num = len(align_names)
    matrix = [[''] * col_num for i in range(row_num)]
    for row, name in enumerate(align_names):
        seq = align_contigs[name]
        for col in range(len(seq)):
            matrix[row][col] = seq[col]
    # Record the base composition of each column.
    col_base_map = {}
    for col_index in range(col_num):
        if not col_base_map.__contains__(col_index):
            col_base_map[col_index] = {}
        base_map = col_base_map[col_index]
        # Calculate the percentage of each base in the current column.
        if len(base_map) == 0:
            for row in range(row_num):
                cur_base = matrix[row][col_index]
                if not base_map.__contains__(cur_base):
                    base_map[cur_base] = 0
                cur_count = base_map[cur_base]
                cur_count += 1
                base_map[cur_base] = cur_count
        if not base_map.__contains__('-'):
            base_map['-'] = 0

    ## Generate a consensus sequence.
    model_seq = ''
    for col_index in range(col_num):
        base_map = col_base_map[col_index]
        # Identify the most frequently occurring base if it exceeds the threshold valid_col_threshold.
        max_base_count = 0
        max_base = ''
        for cur_base in base_map.keys():
            # if cur_base == '-':
            #     continue
            cur_count = base_map[cur_base]
            if cur_count > max_base_count:
                max_base_count = cur_count
                max_base = cur_base
        if max_base_count >= int(row_num / 2):
            if max_base != '-':
                model_seq += max_base
            else:
                continue
        else:
            # Here, we do not use 'N' because it can make it difficult to find the boundary. Therefore, we take the base with the highest non-empty count.
            max_base_count = 0
            max_base = ''
            for cur_base in base_map.keys():
                if cur_base == '-':
                    continue
                cur_count = base_map[cur_base]
                if cur_count > max_base_count:
                    max_base_count = cur_count
                    max_base = cur_base
            model_seq += max_base
    return model_seq

def generate_cons(cluster_id, cur_cluster_path, cluster_dir):
    cur_cons = {}
    cur_contignames, cur_contigs = read_fasta(cur_cluster_path)
    align_file = cur_cluster_path + '.maf.fa'
    if len(cur_contignames) >= 1:
        if not file_exist(align_file):
            align_command = 'cd ' + cluster_dir + ' && mafft --preservecase --quiet --thread 1 ' + cur_cluster_path + ' > ' + align_file
            os.system(align_command)
        cons_seq = cons_from_mafft(align_file)
        cur_cons[cur_contignames[0]] = cons_seq
    return cur_cons

def deredundant_for_LTR(redundant_ltr, work_dir, threads):
    # We found that performing a direct mafft alignment on the redundant LTR library was too slow.
    # Therefore, we first need to use Blastn for alignment clustering, and then proceed with mafft processing.
    tmp_blast_dir = work_dir + '/LTR_blastn'
    blastnResults_path = work_dir + '/LTR_blastn.out'
    # 1. Start by performing an all-vs-all comparison using blastn.
    multi_process_align(redundant_ltr, redundant_ltr, blastnResults_path, tmp_blast_dir, threads, is_removed_dir=True)
    if not os.path.exists(blastnResults_path):
        return redundant_ltr
    # 2. Next, using the FMEA algorithm, bridge across the gaps and link together sequences that can be connected.
    fixed_extend_base_threshold = 1000
    longest_repeats = FMEA(blastnResults_path, fixed_extend_base_threshold)
    # 3. If the combined sequence length constitutes 95% or more of the original individual sequence lengths, we place these two sequences into a cluster.
    contigNames, contigs = read_fasta(redundant_ltr)
    keep_clusters = []
    relations = set()
    for query_name in longest_repeats.keys():
        longest_repeats_list = longest_repeats[query_name]
        for cur_longest_repeat in longest_repeats_list:
            query_name = cur_longest_repeat[0]
            query_len = len(contigs[query_name])
            q_len = abs(cur_longest_repeat[2] - cur_longest_repeat[1])
            subject_name = cur_longest_repeat[3]
            subject_len = len(contigs[subject_name])
            s_len = abs(cur_longest_repeat[4] - cur_longest_repeat[5])
            if float(q_len) / query_len >= 0.95 and float(s_len) / subject_len >= 0.95:
                # we consider the query and subject to be from the same family.
                if (query_name, subject_name) in relations:
                    continue
                relations.add((query_name, subject_name))
                relations.add((subject_name, query_name))
                is_new_cluster = True
                for cluster in keep_clusters:
                    if query_name in cluster or subject_name in cluster:
                        is_new_cluster = False
                        cluster.add(query_name)
                        cluster.add(subject_name)
                        break
                if is_new_cluster:
                    new_cluster = set()
                    new_cluster.add(query_name)
                    new_cluster.add(subject_name)
                    keep_clusters.append(new_cluster)
                    # print(keep_clusters)
    # Iterate through each cluster, if any element in the cluster overlaps with elements in other clusters, merge the clusters.
    merged_clusters = []
    while keep_clusters:
        current_cluster = keep_clusters.pop(0)
        for other_cluster in keep_clusters[:]:
            if current_cluster.intersection(other_cluster):
                current_cluster.update(other_cluster)
                keep_clusters.remove(other_cluster)
        merged_clusters.append(current_cluster)
    keep_clusters = merged_clusters
    # store cluster
    all_unique_name = set()
    raw_cluster_files = []
    cluster_dir = work_dir + '/raw_ltr_cluster'
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)
    for cluster_id, cur_cluster in enumerate(keep_clusters):
        cur_cluster_path = cluster_dir + '/' + str(cluster_id) + '.fa'
        cur_cluster_contigs = {}
        for ltr_name in cur_cluster:
            cur_cluster_contigs[ltr_name] = contigs[ltr_name]
            all_unique_name.add(ltr_name)
        store_fasta(cur_cluster_contigs, cur_cluster_path)
        raw_cluster_files.append((cluster_id, cur_cluster_path))
    # We save the sequences that did not appear in any clusters separately. These sequences do not require clustering.
    uncluster_path = work_dir + '/uncluster_ltr.fa'
    uncluster_contigs = {}
    for name in contigNames:
        if name not in all_unique_name:
            uncluster_contigs[name] = contigs[name]
    store_fasta(uncluster_contigs, uncluster_path)
    # 4. The final cluster should encompass all instances from the same family.
    # We then use the mafft+majority principle to generate a consensus sequence for each cluster.
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cluster_id, cur_cluster_path in raw_cluster_files:
        job = ex.submit(generate_cons, cluster_id, cur_cluster_path, cluster_dir)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_cons = {}
    for job in as_completed(jobs):
        cur_cons = job.result()
        all_cons.update(cur_cons)
    all_cons.update(uncluster_contigs)
    ltr_cons_path = redundant_ltr + '.cons'
    store_fasta(all_cons, ltr_cons_path)
    #rename_fasta(ltr_cons_path, ltr_cons_path, 'LTR')
    return ltr_cons_path

def find_tail_polyA(sequence, min_length=6):
    for i in range(len(sequence) - (min_length-1), -1, -1):
        six_mer = sequence[i:i + min_length]
        if six_mer == 'AAAAAA':
            return i + min_length, six_mer
    return -1, None

def find_nearest_polyA(sequence, position, window_size=25, min_length=6):
    subsequence = sequence[max(0, position - window_size):position + window_size]
    max_length = 0
    current_length = 0
    max_start = 0
    max_end = 0

    for i, base in enumerate(subsequence):
        if base == 'A':
            current_length += 1
            if current_length == 1:
                start = i
        else:
            if current_length >= min_length and current_length > max_length:
                max_length = current_length
                max_start = start
                max_end = i
            current_length = 0

    if current_length >= min_length and current_length > max_length:
        max_start = start
        max_end = len(subsequence)


    max_start = max(0, position - window_size + max_start)
    max_end = max(0, position - window_size + max_end)

    return max_start, max_end, sequence[max_start:max_end]

def find_nearest_polyT(sequence, position, window_size=25, min_length=6):
    subsequence = sequence[max(0, position - window_size):position + window_size]
    max_length = 0
    current_length = 0
    max_start = 0
    max_end = 0

    for i, base in enumerate(subsequence):
        if base == 'T':
            current_length += 1
            if current_length == 1:
                start = i
        else:
            if current_length >= min_length and current_length > max_length:
                max_length = current_length
                max_start = start
                max_end = i
            current_length = 0

    if current_length >= min_length and current_length > max_length:
        max_start = start
        max_end = len(subsequence)

    if max_start != 0:
        max_start = max(0, position - window_size + max_start)
        max_end = max(0, position - window_size + max_end)

    return max_start, max_end, sequence[max_start:max_end]

def search_polyA_TSD(seq, flanking_len, end_5_window_size, TSD_sizes):
    # 2. Identify sequences in longest_repeats.flanked.fa with polyA/T tail or polyA/T head.
    raw_start = flanking_len + 1
    raw_end = len(seq) - flanking_len
    end_3 = -1
    end_5 = -1
    direct = None
    prev_tail_len = 0
    # Obtain the polyA sequence near raw_end.
    start_pos, end_pos, polyA_seq = find_nearest_polyA(seq, raw_end)
    if len(polyA_seq) > 0 and len(polyA_seq) > prev_tail_len:
        end_3 = end_pos
        end_5 = raw_start
        direct = '+'
        prev_tail_len = len(polyA_seq)

    # # Obtain the polyA sequence near raw_start.
    # start_pos, end_pos, polyA_seq = find_nearest_polyA(seq, raw_start)
    # if len(polyA_seq) > 0 and len(polyA_seq) > prev_tail_len:
    #     end_3 = start_pos
    #     end_5 = raw_end
    #     direct = '-'
    #     prev_tail_len = len(polyA_seq)

    # Obtain the polyT sequence near raw_start.
    start_pos, end_pos, polyT_seq = find_nearest_polyT(seq, raw_start)
    if len(polyT_seq) > 0 and len(polyT_seq) > prev_tail_len:
        end_3 = start_pos
        end_5 = raw_end
        direct = '-'
        prev_tail_len = len(polyT_seq)

    # # Obtain the polyT sequence near raw_end.
    # start_pos, end_pos, polyT_seq = find_nearest_polyT(seq, raw_end)
    # if len(polyT_seq) > 0 and len(polyT_seq) > prev_tail_len:
    #     end_3 = end_pos
    #     end_5 = raw_start
    #     direct = '+'
    #     prev_tail_len = len(polyT_seq)

    found_TSD = False
    TSD_seq = ''
    # 3. After locating the 3' end, attempt to search for a set of TSDs with various lengths (up to 20bp)
    # on both sides within 50bp of the 5' end. Allow for 1bp mismatch when TSD length >=8.
    if end_3 != -1 and end_5 != -1 and direct is not None:
        if direct == '-':
            # Retrieve all possible TSDs on the 3' end side.
            TSD_list = [(k, seq[end_3 - k:end_3]) for k in TSD_sizes]
            # Search for TSDs of various lengths (up to 20bp) near the 5' end within 50bp. Allow for 1bp mismatch when TSD length >=8.
            subsequence = seq[max(0, end_5 - end_5_window_size):end_5 + end_5_window_size]
            for k, TSD in reversed(TSD_list):
                for i in range(0, len(subsequence) - k + 1):
                    kmer = subsequence[i:i + k]
                    if len(kmer) >= 8:
                        dist = 1
                    else:
                        dist = 0
                    if k == len(TSD) and k == len(kmer):
                        first_matches = find_near_matches(TSD, kmer, max_l_dist=dist)
                        if len(first_matches) > 0 and not TSD.__contains__('N'):
                            end_5 = max(0, end_5 - end_5_window_size) + i
                            found_TSD = True
                            TSD_seq = TSD
                            break
                if found_TSD:
                    break
        if direct == '+':
            # Retrieve all possible TSDs on the 5' end side.
            TSD_list = [(k, seq[end_3:end_3 + k]) for k in TSD_sizes]
            # Search for TSDs of various lengths (up to 20bp) near the 5' end within 50bp. Allow for 1bp mismatch when TSD length >=8.
            subsequence = seq[max(0, end_5 - end_5_window_size):end_5 + end_5_window_size]
            for k, TSD in reversed(TSD_list):
                for i in range(0, len(subsequence) - k + 1):
                    kmer = subsequence[i:i + k]
                    if len(kmer) >= 8:
                        dist = 1
                    else:
                        dist = 0
                    if k == len(TSD) and k == len(kmer):
                        first_matches = find_near_matches(TSD, kmer, max_l_dist=dist)
                        if len(first_matches) > 0 and not TSD.__contains__('N'):
                            end_5 = max(0, end_5 - end_5_window_size) + i + k
                            found_TSD = True
                            TSD_seq = TSD
                            break
                if found_TSD:
                    break
    if direct is None:
        non_ltr_seq = ''
    else:
        non_ltr_seq = seq[min(end_5, end_3): max(end_5, end_3)]
        if direct == '-':
            non_ltr_seq = getReverseSequence(non_ltr_seq)
    # print(end_5, end_3, found_TSD, TSD_seq)
    # print(non_ltr_seq)
    return found_TSD, TSD_seq, non_ltr_seq

def get_candidate_non_LTR(longest_repeats_flanked_path, flanking_len=50):
    TSD_sizes = list(range(8, 21))
    end_5_window_size = 25
    # 1. Determine the length of the sequence. If the sequence length is between 100-700, it is identified as a SINE.
    # If the sequence length is between 4-8 kb, it is identified as a LINE.
    contig_names, contigs = read_fasta(longest_repeats_flanked_path)
    raw_SINE_names = []
    raw_SINE_contigs = {}
    raw_LINE_names = []
    raw_LINE_contigs = {}
    for name in contig_names:
        seq_len = len(contigs[name]) - 2 * flanking_len
        if seq_len >= 100 and seq_len <= 700:
            raw_SINE_contigs[name] = contigs[name]
            raw_SINE_names.append(name)
        elif seq_len >= 4000 and seq_len <= 8000:
            raw_LINE_contigs[name] = contigs[name]
            raw_LINE_names.append(name)

    candidate_SINE_contigs = {}
    candidate_LINE_contigs = {}
    for name in raw_SINE_names:
        seq = raw_SINE_contigs[name]
        found_TSD, TSD_seq, non_ltr_seq = search_polyA_TSD(seq, flanking_len, end_5_window_size, TSD_sizes)
        seq_len = len(non_ltr_seq)
        # 4. If a TSD is found, it is retained. Otherwise, if it is a SINE, it is filtered out. If it is a LINE,
        # it is aligned with the LINE domain. If a protein sequence is present, it is retained. Otherwise, it is filtered out.
        if found_TSD and seq_len >= 100 and seq_len <= 700:
            new_name = name + '\tTSD:' + TSD_seq
            candidate_SINE_contigs[new_name] = non_ltr_seq
    for name in raw_LINE_names:
        seq = raw_LINE_contigs[name]
        found_TSD, TSD_seq, non_ltr_seq = search_polyA_TSD(seq, flanking_len, end_5_window_size, TSD_sizes)
        seq_len = len(non_ltr_seq)
        if found_TSD and seq_len >= 4000 and seq_len <= 8000:
            new_name = name + '\tTSD:' + TSD_seq
            candidate_LINE_contigs[new_name] = non_ltr_seq
    return candidate_SINE_contigs, candidate_LINE_contigs

def get_candidate_non_ltr_parallel(longest_repeats_flanked_path, work_dir, threads):
    tmp_blast_dir = work_dir + '/non_ltr_split'
    orig_names, orig_contigs = read_fasta(longest_repeats_flanked_path)
    longest_repeat_files = []
    segments_cluster = divided_array(list(orig_contigs.items()), threads)
    for partition_index, cur_segments in enumerate(segments_cluster):
        if len(cur_segments) <= 0:
            continue
        single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
        if not os.path.exists(single_tmp_dir):
            os.makedirs(single_tmp_dir)
        split_repeat_file = single_tmp_dir + '/repeats_split.fa'
        cur_contigs = {}
        for item in cur_segments:
            cur_contigs[item[0]] = item[1]
        store_fasta(cur_contigs, split_repeat_file)
        longest_repeat_files.append(split_repeat_file)

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeat_files:
        job = ex.submit(get_candidate_non_LTR, file)
        jobs.append(job)
    ex.shutdown(wait=True)
    candidate_SINE_contigs = {}
    candidate_LINE_contigs = {}
    for job in as_completed(jobs):
        cur_candidate_SINE_contigs, cur_candidate_LINE_contigs = job.result()
        candidate_SINE_contigs.update(cur_candidate_SINE_contigs)
        candidate_LINE_contigs.update(cur_candidate_LINE_contigs)
    candidate_SINE_path = work_dir + '/candidate_SINE.fa'
    candidate_LINE_path = work_dir + '/candidate_LINE.fa'
    store_fasta(candidate_SINE_contigs, candidate_SINE_path)
    store_fasta(candidate_LINE_contigs, candidate_LINE_path)
    return candidate_SINE_path, candidate_LINE_path


def get_full_length_copies_from_gff(TE_lib, reference, gff_path, tmp_output_dir, threads, divergence_threshold,
                                    full_length_threshold, search_struct, tools_dir):
    ref_names, ref_contigs = read_fasta(reference)

    query_names, query_contigs = read_fasta(TE_lib)
    new_query_contigs = {}
    for name in query_names:
        new_query_contigs[name.split('#')[0]] = query_contigs[name]
    query_contigs = new_query_contigs

    query_records = {}
    with open(gff_path, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            query_name = parts[8].split(' ')[1].replace('"', '').split(':')[1]
            subject_name = parts[0]
            info_parts = parts[8].split(' ')
            q_start = int(info_parts[2])
            q_end = int(info_parts[3])
            if parts[6] == '-':
                s_start = int(parts[4])
                s_end = int(parts[3])
            else:
                s_start = int(parts[3])
                s_end = int(parts[4])
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))

    full_length_copies = {}
    flank_full_length_copies = {}
    copies_direct = {}

    for idx, query_name in enumerate(query_records.keys()):
        subject_dict = query_records[query_name]
        query_len = len(query_contigs[query_name])
        skip_gap = query_len * full_length_threshold
        if str(query_name).__contains__('Helitron'):
            flanking_len = 5
        else:
            flanking_len = 50

        # if there are more than one longest query overlap with the final longest query over 90%,
        # then it probably the true TE
        longest_queries = []
        for subject_name in subject_dict.keys():
            subject_pos = subject_dict[subject_name]

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
                        cur_subject_start = frag[2]
                        cur_query_end = frag[1]
                        prev_subject_end = exist_frag[3]
                        prev_query_end = exist_frag[1]
                        if (cur_subject_start - prev_subject_end < skip_gap and cur_query_end > prev_query_end):
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
                        cur_subject_start = frag[2]
                        cur_query_end = frag[1]
                        prev_subject_end = exist_frag[3]
                        prev_query_end = exist_frag[1]
                        if (prev_subject_end - cur_subject_start < skip_gap and cur_query_end > prev_query_end):
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

                # print('subject pos size: %d' %(len(cur_cluster)))
                # record visited fragments
                visited_frag = {}
                for i in range(len(cur_cluster)):
                    # keep a longest query start from each fragment
                    prev_frag = cur_cluster[i]
                    if visited_frag.__contains__(prev_frag):
                        continue
                    prev_query_start = prev_frag[0]
                    prev_query_end = prev_frag[1]
                    prev_subject_start = prev_frag[2]
                    prev_subject_end = prev_frag[3]
                    prev_query_seq = (min(prev_query_start, prev_query_end), max(prev_query_start, prev_query_end))
                    prev_subject_seq = (
                        min(prev_subject_start, prev_subject_end), max(prev_subject_start, prev_subject_end))
                    prev_query_len = abs(prev_query_end - prev_query_start)
                    prev_subject_len = abs(prev_subject_end - prev_subject_start)
                    cur_longest_query_len = prev_query_len

                    cur_extend_num = 0
                    visited_frag[prev_frag] = 1
                    # try to extend query
                    for j in range(i + 1, len(cur_cluster)):
                        cur_frag = cur_cluster[j]
                        if visited_frag.__contains__(cur_frag):
                            continue
                        cur_query_start = cur_frag[0]
                        cur_query_end = cur_frag[1]
                        cur_subject_start = cur_frag[2]
                        cur_subject_end = cur_frag[3]
                        cur_query_seq = (min(cur_query_start, cur_query_end), max(cur_query_start, cur_query_end))
                        cur_subject_seq = (
                            min(cur_subject_start, cur_subject_end), max(cur_subject_start, cur_subject_end))
                        cur_query_len = abs(cur_query_end - cur_query_start)
                        cur_subject_len = abs(cur_subject_end - cur_subject_start)

                        query_overlap_len = get_overlap_len(cur_query_seq, prev_query_seq)
                        is_same_query = float(query_overlap_len) / cur_query_len >= 0.5 or float(
                            query_overlap_len) / prev_query_len >= 0.5
                        subject_overlap_len = get_overlap_len(prev_subject_seq, cur_subject_seq)
                        is_same_subject = float(subject_overlap_len) / cur_subject_len >= 0.5 or float(
                            subject_overlap_len) / prev_subject_len >= 0.5

                        # could extend
                        # extend right
                        if cur_query_end > prev_query_end:
                            # judge subject direction
                            if prev_subject_start < prev_subject_end and cur_subject_start < cur_subject_end:
                                # +
                                if cur_subject_end > prev_subject_end:
                                    # forward extend
                                    if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                            and cur_subject_start - prev_subject_end < skip_gap:  # \
                                        # and not is_same_query and not is_same_subject:
                                        # update the longest path
                                        prev_query_start = prev_query_start
                                        prev_query_end = cur_query_end
                                        prev_subject_start = prev_subject_start if prev_subject_start < cur_subject_start else cur_subject_start
                                        prev_subject_end = cur_subject_end
                                        cur_longest_query_len = prev_query_end - prev_query_start
                                        cur_extend_num += 1
                                        visited_frag[cur_frag] = 1
                                    elif cur_query_start - prev_query_end >= skip_gap:
                                        break
                            elif prev_subject_start > prev_subject_end and cur_subject_start > cur_subject_end:
                                # reverse
                                if cur_subject_end < prev_subject_end:
                                    # reverse extend
                                    if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                            and prev_subject_end - cur_subject_start < skip_gap:  # \
                                        # and not is_same_query and not is_same_subject:
                                        # update the longest path
                                        prev_query_start = prev_query_start
                                        prev_query_end = cur_query_end
                                        prev_subject_start = prev_subject_start if prev_subject_start > cur_subject_start else cur_subject_start
                                        prev_subject_end = cur_subject_end
                                        cur_longest_query_len = prev_query_end - prev_query_start
                                        cur_extend_num += 1
                                        visited_frag[cur_frag] = 1
                                    elif cur_query_start - prev_query_end >= skip_gap:
                                        break
                    # keep this longest query
                    if cur_longest_query_len != -1:
                        longest_queries.append(
                            (prev_query_start, prev_query_end, cur_longest_query_len, prev_subject_start,
                             prev_subject_end, abs(prev_subject_end - prev_subject_start), subject_name,
                             cur_extend_num))

        # To determine whether each copy has a coverage exceeding the full_length_threshold with respect
        # to the consensus sequence, retaining full-length copies.
        query_copies = {}
        flank_query_copies = {}
        # query_copies[query_name] = query_contigs[query_name]
        for repeat in longest_queries:
            if repeat[2] < full_length_threshold * query_len:
                continue
            # Subject
            subject_name = repeat[6]
            subject_chr_start = 0

            if repeat[3] > repeat[4]:
                direct = '-'
                old_subject_start_pos = repeat[4] - 1
                old_subject_end_pos = repeat[3]
            else:
                direct = '+'
                old_subject_start_pos = repeat[3] - 1
                old_subject_end_pos = repeat[4]
            subject_start_pos = subject_chr_start + old_subject_start_pos
            subject_end_pos = subject_chr_start + old_subject_end_pos

            subject_pos = subject_name + ':' + str(subject_start_pos) + '-' + str(subject_end_pos)
            subject_seq = ref_contigs[subject_name][subject_start_pos: subject_end_pos]

            flank_subject_seq = ref_contigs[subject_name][
                                subject_start_pos - flanking_len: subject_end_pos + flanking_len]
            copies_direct[subject_pos] = direct
            cur_query_len = repeat[2]
            cur_subject_len = repeat[5]
            min_cur_len = min(cur_query_len, cur_subject_len)
            max_cur_len = max(cur_query_len, cur_subject_len)
            coverage = float(min_cur_len) / max_cur_len
            if coverage >= full_length_threshold:
                query_copies[subject_pos] = subject_seq
                flank_query_copies[subject_pos] = flank_subject_seq
        full_length_copies[query_name] = query_copies
        flank_full_length_copies[query_name] = flank_query_copies

    # The candidate full-length copies and the consensus are then clustered using cd-hit-est,
    # retaining copies that belong to the same cluster as the consensus.
    split_files = []
    cluster_dir = tmp_output_dir + '/cluster'
    os.system('rm -rf ' + cluster_dir)
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)

    all_query_copies = {}
    for query_name in full_length_copies.keys():
        query_copies = full_length_copies[query_name]
        flank_query_copies = flank_full_length_copies[query_name]
        all_query_copies.update(query_copies)
        fc_path = cluster_dir + '/' + query_name + '.fa'
        fc_cons = cluster_dir + '/' + query_name + '.cons.fa'
        store_fasta(query_copies, fc_path)
        split_files.append((fc_path, query_name, query_copies, flank_query_copies))

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for ref_index, cur_file in enumerate(split_files):
        input_file = cur_file[0]
        query_name = cur_file[1]
        query_copies = cur_file[2]
        flank_query_copies = cur_file[3]
        job = ex.submit(get_structure_info, input_file, query_name, query_copies,
                        flank_query_copies, cluster_dir, search_struct, tools_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    full_length_annotations = {}
    for job in as_completed(jobs):
        annotations = job.result()
        full_length_annotations.update(annotations)
    return full_length_annotations, copies_direct

def get_full_length_copies_RM(TE_lib, reference, tmp_output_dir, threads, divergence_threshold, full_length_threshold,
                              search_struct, tools_dir):
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)
    tmp_TE_out = tmp_output_dir + '/TE_tmp.out'
    tmp_TE_gff = tmp_output_dir + '/TE_tmp.gff'

    RepeatMasker_command = 'cd ' + tmp_output_dir + ' && RepeatMasker -e ncbi -pa ' + str(threads) \
                           + ' -s -no_is -norna -nolow -div ' + str(divergence_threshold) \
                           + ' -gff -lib ' + TE_lib + ' -cutoff 225 ' + reference
    os.system(RepeatMasker_command + '> /dev/null 2>&1')

    mv_file_command = 'mv ' + reference + '.out ' + tmp_TE_out + ' && mv ' + reference + '.out.gff ' + tmp_TE_gff
    os.system(mv_file_command)

    full_length_annotations, copies_direct = get_full_length_copies_from_gff(TE_lib, reference, tmp_TE_gff,
                                                    tmp_output_dir, threads, divergence_threshold,
                                                    full_length_threshold, search_struct, tools_dir)
    return full_length_annotations, copies_direct

def multiple_alignment_blast_v1(repeats_path, tools_dir, coverage_threshold, category):
    split_repeats_path = repeats_path[0]
    target_files = repeats_path[1]
    blastn2Results_path = repeats_path[2]
    full_length_out = repeats_path[3]
    tmp_dir = repeats_path[4]
    genome_path = repeats_path[5]
    os.system('rm -f ' + blastn2Results_path)
    for target_file in target_files:
        align_command = 'blastn -db ' + target_file + ' -num_threads ' \
                        + str(1) + ' -query ' + split_repeats_path + ' -evalue 1e-20 -outfmt 6 >> ' + blastn2Results_path
        os.system(align_command)

    # invoke the function to retrieve the full-length copies.
    lines = generate_full_length_out(blastn2Results_path, full_length_out, split_repeats_path, genome_path, tmp_dir, tools_dir,
                             coverage_threshold, category)

    return lines

def multi_process_align_v1(query_path, subject_path, blastnResults_path,
                           tmp_blast_dir, threads, coverage_threshold,
                           category, is_removed_dir=True):
    tools_dir = ''
    if is_removed_dir:
        os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    if os.path.exists(blastnResults_path):
        os.remove(blastnResults_path)

    orig_names, orig_contigs = read_fasta(query_path)

    # blast_db_command = 'makeblastdb -dbtype nucl -in ' + subject_path + ' > /dev/null 2>&1'
    # os.system(blast_db_command)

    ref_names, ref_contigs = read_fasta(subject_path)
    # Sequence alignment consumes a significant amount of memory and disk space. Therefore, we also split the target sequences into individual sequences to reduce the memory required for each alignment, avoiding out of memory errors.
    # It is important to calculate the total number of bases in the sequences, and it must meet a sufficient threshold to increase CPU utilization.
    base_threshold = 10000000  # 10Mb
    target_files = []
    file_index = 0
    base_count = 0
    cur_contigs = {}
    for name in ref_names:
        cur_seq = ref_contigs[name]
        cur_contigs[name] = cur_seq
        base_count += len(cur_seq)
        if base_count >= base_threshold:
            cur_target = tmp_blast_dir + '/' + str(file_index) + '_target.fa'
            store_fasta(cur_contigs, cur_target)
            target_files.append(cur_target)
            makedb_command = 'makeblastdb -dbtype nucl -in ' + cur_target + ' > /dev/null 2>&1'
            os.system(makedb_command)
            cur_contigs = {}
            file_index += 1
            base_count = 0
    if len(cur_contigs) > 0:
        cur_target = tmp_blast_dir + '/' + str(file_index) + '_target.fa'
        store_fasta(cur_contigs, cur_target)
        target_files.append(cur_target)
        makedb_command = 'makeblastdb -dbtype nucl -in ' + cur_target + ' > /dev/null 2>&1'
        os.system(makedb_command)


    longest_repeat_files = []
    segments_cluster = divided_array(list(orig_contigs.items()), 2 * threads)
    for partition_index, cur_segments in enumerate(segments_cluster):
        if len(cur_segments) <= 0:
            continue
        single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
        #print('current partition_index: ' + str(partition_index))
        if not os.path.exists(single_tmp_dir):
            os.makedirs(single_tmp_dir)
        split_repeat_file = single_tmp_dir + '/repeats_split.fa'
        cur_contigs = {}
        for item in cur_segments:
            cur_contigs[item[0]] = item[1]
        store_fasta(cur_contigs, split_repeat_file)
        repeats_path = (split_repeat_file, target_files, single_tmp_dir + '/temp.out',
                        single_tmp_dir + '/full_length.out', single_tmp_dir + '/tmp',
                        subject_path)
        longest_repeat_files.append(repeats_path)

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeat_files:
        job = ex.submit(multiple_alignment_blast_v1, file, tools_dir, coverage_threshold, category)
        jobs.append(job)
    ex.shutdown(wait=True)

    lines = set()
    for job in as_completed(jobs):
        cur_lines = job.result()
        lines.update(cur_lines)
    lines = list(lines)
    sorted_lines = sorted(lines, key=lambda x: (x[1], x[2], x[3]))

    with open(blastnResults_path, 'w') as f_save:
        for line in sorted_lines:
            new_line = line[0] + '\t' + line[1] + '\t' + '-1' + '\t' + '-1' + '\t' + '-1' + '\t' + '-1' + '\t' + '-1' + '\t' + '-1' + '\t' + str(line[2]) + '\t' + str(line[3]) + '\t' + '-1' + '\t' + '-1' + '\n'
            f_save.write(new_line)

    if is_removed_dir:
        os.system('rm -rf ' + tmp_blast_dir)

def ReassignInconsistentLabels(TE_lib):
    # Set conflicting labels to "Unknown".
    TE_names, TE_contigs = read_fasta(TE_lib)
    new_TE_contigs = {}
    for name in TE_names:
        parts = name.split('#')
        raw_name = parts[0]
        label = parts[1]
        if ('LTR' in raw_name and 'LTR' not in label) \
                or ('TIR' in raw_name and 'DNA' not in label) \
                or ('Helitron' in raw_name and 'Helitron' not in label) \
                or ('Non_LTR' in raw_name and ('LINE' not in label or 'SINE' not in label)):
            new_label = 'Unknown'
            new_name = raw_name + '#' + new_label
            new_TE_contigs[new_name] = TE_contigs[name]
        else:
            new_TE_contigs[name] = TE_contigs[name]
    store_fasta(new_TE_contigs, TE_lib)

def lib_add_prefix(HiTE_lib, prefix):
    lib_names, lib_contigs = read_fasta(HiTE_lib)
    new_lib_contigs = {}
    for name in lib_names:
        new_name = prefix + '-' + name
        new_lib_contigs[new_name] = lib_contigs[name]
    store_fasta(new_lib_contigs, HiTE_lib)
    return HiTE_lib

def find_gene_relation_tes(genome_info_list, output_dir, recover, log):
    genome_num = len(genome_info_list)
    genome_num_threshold = int(0.8 * genome_num)

    te_gene_insertions_list = []
    for genome_info in genome_info_list:
        genome_name = genome_info["genome_name"]
        gene_gtf = genome_info["gene_gtf"]
        full_length_TE_gff = genome_info["full_length_TE_gff"]
        if full_length_TE_gff is not None \
                and os.path.exists(full_length_TE_gff) \
                and gene_gtf is not None \
                and os.path.exists(gene_gtf):
            output_file = output_dir + '/' + genome_name + '.te_gene_insertions.tsv'
            resut_file = output_file
            if not recover or not file_exist(resut_file):
                results = analyze_te_insertions(gene_gtf, full_length_TE_gff)
                save_results(results, resut_file)
            else:
                log.logger.info(resut_file + ' exists, skip...')
            te_gene_insertions_list.append((genome_name, resut_file))

    # 初始化一个字典，gene name为键，每个键的值是一个列表
    gene_te_associations = defaultdict(list)
    # 遍历te_gene_insertions_list
    for genome_name, resut_file in te_gene_insertions_list:
        results = load_results(resut_file)
        for _, row in results.iterrows():
            gene_name = row['Gene_name']
            gene_name_parts = str(gene_name).split('_')
            gene_name = gene_name_parts[-1]
            gene_te_associations[gene_name].append({
                'Genome_name': genome_name,
                'TE_name': row['TE_name'],
                'Chromosome': row['Chromosome'],
                'TE_start': row['TE_start'],
                'TE_end': row['TE_end'],
                'Gene_start': row['Gene_start'],
                'Gene_end': row['Gene_end'],
                'Position': row['Position']
            })
    output_file = output_dir + '/gene_te_associations.tsv'
    save_gene_te_associations_to_file(gene_te_associations, output_file)

    # 统计同一个 gene 下，相同的位置，相同的TE 出现了多少次
    gene_te_stats = {}
    for gene_name, te_list in gene_te_associations.items():
        # 初始化gene_te_stats的key值
        if gene_name not in gene_te_stats:
            gene_te_stats[gene_name] = defaultdict(lambda: {'count': 0, 'genomes': set()})
        for te_info in te_list:
            # 创建唯一标识符：TE的名称 + TE的位置 + TE的染色体
            te_key = (te_info['Position'], te_info['TE_name'])

            # 更新统计信息
            gene_te_stats[gene_name][te_key]['count'] += 1
            gene_te_stats[gene_name][te_key]['genomes'].add(te_info['Genome_name'])

    # 查看结果
    core_gene_te_list = []
    softcore_gene_te_list = []
    dispensable_gene_te_list = []
    private_gene_te_list = []
    unknown_gene_te_list = []
    for gene_name, stats  in gene_te_stats.items():
        for te_key, info in stats.items():
            position, te_name = te_key
            count = info['count']
            genome_set = info['genomes']
            if count == genome_num:
                core_gene_te_list.append((gene_name, te_name, position, count, genome_set))
            elif genome_num_threshold <= count < genome_num:
                softcore_gene_te_list.append((gene_name, te_name, position, count, genome_set))
            elif 2 <= count < genome_num_threshold:
                dispensable_gene_te_list.append((gene_name, te_name, position, count, genome_set))
            elif count == 1:
                private_gene_te_list.append((gene_name, te_name, position, count, genome_set))
            else:
                unknown_gene_te_list.append((gene_name, te_name, position, count, genome_set))
    output_file = output_dir + '/core_gene_te_relations.tsv' # 相同 gene 与 TE 的关系在 所有 基因组上出现
    save_gene_te_list(core_gene_te_list, output_file)
    output_file = output_dir + '/softcore_gene_te_relations.tsv'  # 相同 gene 与 TE 的关系在 绝大部分 基因组上出现
    save_gene_te_list(softcore_gene_te_list, output_file)
    output_file = output_dir + '/dispensable_gene_te_relations.tsv'  # 相同 gene 与 TE 的关系在 少部分 基因组上出现
    save_gene_te_list(dispensable_gene_te_list, output_file)
    output_file = output_dir + '/private_gene_te_relations.tsv'  # 相同 gene 与 TE 的关系在 单个 基因组上出现
    save_gene_te_list(private_gene_te_list, output_file)
    output_file = output_dir + '/unknown_gene_te_relations.tsv'  # 相同 gene 与 TE 的关系 没有在任一个  基因组上出现
    save_gene_te_list(unknown_gene_te_list, output_file)

def save_gene_te_list(core_gene_te_list, output_file):
    with open(output_file, 'w') as f_save:
        for item in core_gene_te_list:
            f_save.write(item[0] + '\t' + item[1] + '\t' + item[2] + '\t' + str(item[3]) + '\t' + str(item[4]) + '\n')


def analyze_te_insertions(gene_file, te_file, upstream_downstream=10_000):
    genes_by_chrom = read_gff(gene_file, feature_type='gene')
    tes_by_chrom = read_gff(te_file)

    results = []

    for chrom, genes in genes_by_chrom.items():
        if chrom in tes_by_chrom:
            tes = tes_by_chrom[chrom]
            te_index = 0

            for idx, gene in genes.iterrows():
                gene_name = extract_name(gene['attribute'], file_type='gtf')
                next_gene = genes.iloc[idx + 1] if idx + 1 < len(genes) else None

                # 剪枝：跳过TE的起始位置小于当前基因的上游范围的TE
                while te_index < len(tes) and tes.iloc[te_index]['end'] < gene['start'] - upstream_downstream:
                    te_index += 1

                # 从当前索引开始检查TE
                current_te_index = te_index

                while current_te_index < len(tes):
                    te = tes.iloc[current_te_index]
                    te_name = extract_name(te['attribute'], file_type='gff')

                    # 如果TE的起始位置超过了基因的下游范围，停止检查
                    if te['start'] > gene['end'] + upstream_downstream:
                        break

                    # 如果TE和下一个基因有交集了或超过了下一个基因，说明TE和当前基因没有关系了
                    if (next_gene is not None) and (te['end'] > next_gene['start']):
                        break

                    # 检查TE与基因的关系
                    position = check_te_in_gene(te, gene, upstream_downstream)
                    if position != 'None':
                        results.append({
                            'Gene_name': gene_name,
                            'TE_name': te_name,
                            'Chromosome': te['chromosome'],
                            'TE_start': int(te['start']),
                            'TE_end': int(te['end']),
                            'Gene_start': int(gene['start']),
                            'Gene_end': int(gene['end']),
                            'Position': position
                        })

                    current_te_index += 1

    return pd.DataFrame(results)


def save_results(df, output_file):
    df.to_csv(output_file, sep='\t', index=False)


def load_results(input_file):
    df = pd.read_csv(input_file, sep='\t')
    return df

# 读取gff文件，并只保留feature为gene的记录
def read_gff(file, feature_type=None):
    df = pd.read_csv(file, sep='\t', header=None,
                     names=['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    if feature_type:
        df = df[df['feature'] == feature_type]

    # 按染色体分组
    grouped = df.groupby('chromosome')
    return {chrom: grouped.get_group(chrom) for chrom in grouped.groups}

def extract_name(attribute, file_type='gff'):
    if file_type == 'gtf':
        match = re.search(r'gene_id\s+"([^";]+)"', attribute)
        if match:
            return match.group(1)
        else:
            return None
    else:
        match = re.search(r'(?:Name|name)=([^;]+)', attribute)
        if match:
            return match.group(1)
        else:
            return None

# 判断TE是否插入到gene的上下游或内部
def check_te_in_gene(te, gene, upstream_downstream=10_000):
    gene_start, gene_end, gene_strand = gene['start'], gene['end'], gene['strand']
    te_start, te_end = te['start'], te['end']

    if te_end < gene_start - upstream_downstream:
        return 'None'
    elif te_start >= gene_start and te_end <= gene_end:
        return 'Inside'
    elif te_start < gene_start and te_end >= gene_start - upstream_downstream:
        if gene_strand == '+':
            return 'Upstream'
        else:
            return 'Downstream'
    elif te_start <= gene_end + upstream_downstream and te_end > gene_end:
        if gene_strand == '+':
            return 'Downstream'
        else:
            return 'Upstream'
    else:
        return 'None'

# 将gene_te_associations转换为DataFrame
def save_gene_te_associations_to_file(gene_te_associations, output_file):
    rows = []
    for gene_name, associations in gene_te_associations.items():
        for association in associations:
            row = {
                'Gene_name': gene_name,
                'Genome_name': association['Genome_name'],
                'TE_name': association['TE_name'],
                'Chromosome': association['Chromosome'],
                'TE_start': association['TE_start'],
                'TE_end': association['TE_end'],
                'Gene_start': association['Gene_start'],
                'Gene_end': association['Gene_end'],
                'Position': association['Position']
            }
            rows.append(row)

    # 创建DataFrame
    df = pd.DataFrame(rows)

    # 保存为文件
    df.to_csv(output_file, sep='\t', index=False)

def deredundant_for_LTR_v5(redundant_ltr, work_dir, threads, type, coverage_threshold, debug):
    starttime = time.time()
    # We found that performing a direct mafft alignment on the redundant LTR library was too slow.
    # Therefore, we first need to use Blastn for alignment clustering, and then proceed with mafft processing.
    tmp_blast_dir = work_dir + '/LTR_blastn_' + str(type)
    blastnResults_path = work_dir + '/LTR_blastn_' + str(type) + '.out'
    # 1. Start by performing an all-vs-all comparison using blastn.
    multi_process_align(redundant_ltr, redundant_ltr, blastnResults_path, tmp_blast_dir, threads, is_removed_dir=True, is_remove_index=True)
    if not os.path.exists(blastnResults_path):
        return redundant_ltr
    # 2. Next, using the FMEA algorithm, bridge across the gaps and link together sequences that can be connected.
    query_records = {}
    with open(blastnResults_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            # print('current line idx: %d' % (idx))
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            if query_name == subject_name and q_start == s_start and q_end == s_end:
                continue
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))
    f_r.close()
    query_names, query_contigs = read_fasta(redundant_ltr)
    longest_repeats = FMEA_new(query_contigs, query_records, coverage_threshold)

    # 3. If the combined sequence length constitutes 95% or more of the original individual sequence lengths, we place these two sequences into a cluster.
    contigNames, contigs = read_fasta(redundant_ltr)
    keep_clusters = []
    redundant_ltr_names = set()
    for query_name in longest_repeats.keys():
        if query_name in redundant_ltr_names:
            continue
        longest_repeats_list = longest_repeats[query_name]
        cur_cluster = set()
        cur_cluster.add(query_name)
        for cur_longest_repeat in longest_repeats_list:
            query_name = cur_longest_repeat[0]
            query_len = len(contigs[query_name])
            q_len = abs(cur_longest_repeat[2] - cur_longest_repeat[1])
            subject_name = cur_longest_repeat[3]
            subject_len = len(contigs[subject_name])
            s_len = abs(cur_longest_repeat[4] - cur_longest_repeat[5])
            # 我们这里先将跨过 gap 之后的全长拷贝先聚类在一起，后续再使用 cd-hit 将碎片化合并到全长拷贝中
            if float(q_len) / query_len >= coverage_threshold or float(s_len) / subject_len >= coverage_threshold:
                # we consider the query and subject to be from the same family.
                cur_cluster.add(subject_name)
                redundant_ltr_names.add(subject_name)
        keep_clusters.append(cur_cluster)

    endtime = time.time()
    dtime = endtime - starttime
    # print("Running time of FMEA clustering: %.8s s" % (dtime))
    if not debug:
        os.system('rm -f ' + blastnResults_path)
        os.system('rm -rf ' + tmp_blast_dir)


    starttime = time.time()
    # store cluster
    all_unique_name = set()
    raw_cluster_files = []
    cluster_dir = work_dir + '/raw_ltr_cluster_' + str(type)
    os.system('rm -rf ' + cluster_dir)
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)
    for cluster_id, cur_cluster in enumerate(keep_clusters):
        cur_cluster_path = cluster_dir + '/' + str(cluster_id) + '.fa'
        cur_cluster_contigs = {}
        for ltr_name in cur_cluster:
            cur_cluster_contigs[ltr_name] = contigs[ltr_name]
            all_unique_name.add(ltr_name)
        store_fasta(cur_cluster_contigs, cur_cluster_path)
        raw_cluster_files.append((cluster_id, cur_cluster_path))
    # We save the sequences that did not appear in any clusters separately. These sequences do not require clustering.
    uncluster_path = work_dir + '/uncluster_ltr_' + str(type) + '.fa'
    uncluster_contigs = {}
    for name in contigNames:
        if name not in all_unique_name:
            uncluster_contigs[name] = contigs[name]
    store_fasta(uncluster_contigs, uncluster_path)

    # 4. The final cluster should encompass all instances from the same family.
    # We use Ninja to cluster families precisely, and
    # We then use the mafft+majority principle to generate a consensus sequence for each cluster.
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cluster_id, cur_cluster_path in raw_cluster_files:
        job = ex.submit(generate_cons_v1, cluster_id, cur_cluster_path, cluster_dir, 1)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_cons = {}
    for job in as_completed(jobs):
        cur_cons_contigs = job.result()
        all_cons.update(cur_cons_contigs)

    all_cons.update(uncluster_contigs)

    ltr_cons_path = redundant_ltr + '.tmp.cons'
    store_fasta(all_cons, ltr_cons_path)

    endtime = time.time()
    dtime = endtime - starttime
    # print("Running time of MSA cons: %.8s s" % (dtime))
    if not debug:
        os.system('rm -f ' + uncluster_path)
        os.system('rm -rf ' + cluster_dir)

    ltr_cons_cons = redundant_ltr + '.cons'
    # 调用 cd-hit-est 合并碎片化序列
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(coverage_threshold) \
                     + ' -G 0 -g 1 -A 80 -i ' + ltr_cons_path + ' -o ' + ltr_cons_cons + ' -T 0 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')

    #rename_fasta(ltr_cons_path, ltr_cons_path, 'LTR')
    return ltr_cons_path

def FMEA_new(query_contigs, query_records, full_length_threshold):
    # 我们现在尝试新的策略，直接在生成簇的时候进行扩展，同时新的比对片段和所有的扩展片段进行比较，判断是否可以扩展
    longest_repeats = {}
    for idx, query_name in enumerate(query_records.keys()):
        subject_dict = query_records[query_name]

        if query_name not in query_contigs:
            continue
        query_len = len(query_contigs[query_name])
        skip_gap = query_len * (1 - full_length_threshold)

        longest_queries = []
        for subject_name in subject_dict.keys():
            subject_pos = subject_dict[subject_name]

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

            forward_long_frags = {}
            frag_index_array = []
            frag_index = 0
            for k, frag in enumerate(forward_pos):
                is_update = False
                cur_subject_start = frag[2]
                cur_subject_end = frag[3]
                cur_query_start = frag[0]
                cur_query_end = frag[1]
                for cur_frag_index in reversed(frag_index_array):
                    cur_frag = forward_long_frags[cur_frag_index]
                    prev_subject_start = cur_frag[2]
                    prev_subject_end = cur_frag[3]
                    prev_query_start = cur_frag[0]
                    prev_query_end = cur_frag[1]

                    if cur_subject_start - prev_subject_end >= skip_gap:
                        break

                    if cur_subject_end > prev_subject_end:
                        # forward extend
                        if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                and cur_subject_start - prev_subject_end < skip_gap:  # \
                            # extend frag
                            prev_query_start = prev_query_start if prev_query_start < cur_query_start else cur_query_start
                            prev_query_end = cur_query_end
                            prev_subject_start = prev_subject_start if prev_subject_start < cur_subject_start else cur_subject_start
                            prev_subject_end = cur_subject_end
                            extend_frag = (prev_query_start, prev_query_end, prev_subject_start, prev_subject_end, subject_name)
                            forward_long_frags[cur_frag_index] = extend_frag
                            is_update = True
                if not is_update:
                    frag_index_array.append(frag_index)
                    forward_long_frags[frag_index] = (cur_query_start, cur_query_end, cur_subject_start, cur_subject_end, subject_name)
                    frag_index += 1
            longest_queries += list(forward_long_frags.values())

            reverse_long_frags = {}
            frag_index_array = []
            frag_index = 0
            for k, frag in enumerate(reverse_pos):
                is_update = False
                cur_subject_start = frag[2]
                cur_subject_end = frag[3]
                cur_query_start = frag[0]
                cur_query_end = frag[1]
                for cur_frag_index in reversed(frag_index_array):
                    cur_frag = reverse_long_frags[cur_frag_index]
                    prev_subject_start = cur_frag[2]
                    prev_subject_end = cur_frag[3]
                    prev_query_start = cur_frag[0]
                    prev_query_end = cur_frag[1]

                    if prev_subject_end - cur_subject_start >= skip_gap:
                        break

                    # reverse
                    if cur_subject_end < prev_subject_end:
                        # reverse extend
                        if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                and prev_subject_end - cur_subject_start < skip_gap:  # \
                            # extend frag
                            prev_query_start = prev_query_start
                            prev_query_end = cur_query_end
                            prev_subject_start = prev_subject_start if prev_subject_start > cur_subject_start else cur_subject_start
                            prev_subject_end = cur_subject_end
                            extend_frag = (prev_query_start, prev_query_end, prev_subject_start, prev_subject_end, subject_name)
                            reverse_long_frags[cur_frag_index] = extend_frag
                            is_update = True
                if not is_update:
                    frag_index_array.append(frag_index)
                    reverse_long_frags[frag_index] = (cur_query_start, cur_query_end, cur_subject_start, cur_subject_end, subject_name)
                    frag_index += 1
            longest_queries += list(reverse_long_frags.values())

        if not longest_repeats.__contains__(query_name):
            longest_repeats[query_name] = []
        cur_longest_repeats = longest_repeats[query_name]
        for repeat in longest_queries:
            # Subject序列处理流程
            subject_name = repeat[4]
            old_subject_start_pos = repeat[2] - 1
            old_subject_end_pos = repeat[3]
            # Query序列处理流程
            old_query_start_pos = repeat[0] - 1
            old_query_end_pos = repeat[1]
            cur_query_seq_len = abs(old_query_end_pos - old_query_start_pos)
            cur_longest_repeats.append((query_name, old_query_start_pos, old_query_end_pos, subject_name, old_subject_start_pos, old_subject_end_pos))

    return longest_repeats

def generate_cons_v1(cluster_id, cur_cluster_path, cluster_dir, threads):
    ltr_terminal_names, ltr_terminal_contigs = read_fasta(cur_cluster_path)
    temp_cluster_dir = cluster_dir
    cons_contigs = {}
    if len(ltr_terminal_contigs) >= 1:
        align_file = cur_cluster_path + '.maf.fa'
        align_command = 'cd ' + cluster_dir + ' && mafft --preservecase --quiet --thread ' + str(threads) + ' ' + cur_cluster_path + ' > ' + align_file
        # align_command = 'cd ' + cluster_dir + ' && famsa -t ' + str(threads) + ' -medoidtree ' + cur_cluster_path + ' ' + align_file + ' > /dev/null 2>&1'
        os.system(align_command)

        # 调用 Ninja 对多序列比对再次聚类
        cluster_file = align_file + '.dat'
        Ninja_command = 'Ninja --in ' + align_file + ' --out ' + cluster_file + ' --out_type c --corr_type m --cluster_cutoff 0.2 --threads ' + str(threads)
        os.system(Ninja_command + ' > /dev/null 2>&1')

        # 解析聚类文件，生成不同簇
        Ninja_cluster_dir = temp_cluster_dir + '/Ninja_' + str(cluster_id)
        if not os.path.exists(Ninja_cluster_dir):
            os.makedirs(Ninja_cluster_dir)
        clusters = read_Ninja_clusters(cluster_file)
        for cur_cluster_id in clusters.keys():
            cur_cluster_file = Ninja_cluster_dir + '/' + str(cur_cluster_id) + '.fa'
            cur_cluster_contigs = {}
            cur_ltr_name = ''
            for name in clusters[cur_cluster_id]:
                seq = ltr_terminal_contigs[name]
                cur_cluster_contigs[name] = seq
                cur_ltr_name = name
            store_fasta(cur_cluster_contigs, cur_cluster_file)

            cur_align_file = cur_cluster_file + '.maf.fa'
            if len(cur_cluster_contigs) >= 1:
                align_command = 'cd ' + Ninja_cluster_dir + ' && mafft --preservecase --quiet --thread ' + str(threads) + ' ' + cur_cluster_file + ' > ' + cur_align_file
                # align_command = 'cd ' + Ninja_cluster_dir + ' && famsa -t ' + str(threads) + ' -medoidtree ' + cur_cluster_file + ' ' + cur_align_file + ' > /dev/null 2>&1'
                os.system(align_command)
                cons_seq = cons_from_mafft_v1(cur_align_file)
                cons_contigs[cur_ltr_name] = cons_seq
    # 如果未能识别到可靠的一致性序列，则使用原始序列代替
    if len(cons_contigs) > 0:
        return cons_contigs
    else:
        return ltr_terminal_contigs

def read_Ninja_clusters(cluster_file):
    clusters = {}
    with open(cluster_file, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            cluster_id = int(parts[0])
            seq_name = parts[1]
            if cluster_id not in clusters:
                clusters[cluster_id] = []
            cur_cluster = clusters[cluster_id]
            cur_cluster.append(seq_name)

    return clusters

def cons_from_mafft_v1(align_file):
    align_names, align_contigs = read_fasta(align_file)
    if len(align_names) <= 0:
        return None

    # Generate a consensus sequence using full-length copies.
    first_seq = align_contigs[align_names[0]]
    col_num = len(first_seq)
    row_num = len(align_names)
    matrix = [[''] * col_num for i in range(row_num)]
    for row, name in enumerate(align_names):
        seq = align_contigs[name]
        for col in range(len(seq)):
            matrix[row][col] = seq[col]
    # Record the base composition of each column.
    col_base_map = {}
    for col_index in range(col_num):
        if not col_base_map.__contains__(col_index):
            col_base_map[col_index] = {}
        base_map = col_base_map[col_index]
        # Calculate the percentage of each base in the current column.
        if len(base_map) == 0:
            for row in range(row_num):
                cur_base = matrix[row][col_index]
                if not base_map.__contains__(cur_base):
                    base_map[cur_base] = 0
                cur_count = base_map[cur_base]
                cur_count += 1
                base_map[cur_base] = cur_count
        if not base_map.__contains__('-'):
            base_map['-'] = 0

    ## Generate a consensus sequence.
    model_seq = ''
    for col_index in range(col_num):
        base_map = col_base_map[col_index]
        # Identify the most frequently occurring base if it exceeds the threshold valid_col_threshold.
        max_base_count = 0
        max_base = ''
        for cur_base in base_map.keys():
            if cur_base == '-':
                continue
            cur_count = base_map[cur_base]
            if cur_count > max_base_count:
                max_base_count = cur_count
                max_base = cur_base
        if max_base_count > int(row_num / 2):
            if max_base != '-':
                model_seq += max_base
            else:
                continue
    return model_seq

def generate_bam_for_RNA_seq(genome_info_list, threads, recover, RNA_dir, log):
    new_batch_files = []
    for genome_info in genome_info_list:
        genome_name = genome_info["genome_name"]
        reference = genome_info["reference"]
        TE_gff = genome_info["TE_gff"]
        gene_gtf = genome_info["gene_gtf"]
        RNA_seq_dict = genome_info["RNA_seq"]
        full_length_TE_gff = genome_info["full_length_TE_gff"]

        if len(RNA_seq_dict) > 0:
            is_PE = RNA_seq_dict['is_PE']

            if is_PE:
                raw_RNA1 = os.path.join(RNA_dir, RNA_seq_dict['raw_RNA1'])
                raw_RNA2 = os.path.join(RNA_dir, RNA_seq_dict['raw_RNA2'])
                sorted_bam = RNA_dir + '/' + genome_name + '.output.sorted.bam'
                resut_file = sorted_bam
                if not recover or not os.path.exists(resut_file):
                    generate_bam(genome_path=reference, genome_annotation_file=gene_gtf,
                                             output_dir=RNA_dir, threads=threads, is_PE=True, raw_RNA1=raw_RNA1,
                                             raw_RNA2=raw_RNA2)
                else:
                    log.logger.info(resut_file + ' exists, skip...')

            else:
                raw_RNA = os.path.join(RNA_dir, RNA_seq_dict['raw_RNA'])
                sorted_bam = RNA_dir + '/' + genome_name + '.output.sorted.bam'
                resut_file = sorted_bam
                if not recover or not os.path.exists(resut_file):
                    generate_bam(genome_path=reference, output_dir=RNA_dir, threads=threads, is_PE=False, raw_RNA=raw_RNA)
                else:
                    log.logger.info(resut_file + ' exists, skip...')

            new_batch_files.append((genome_name, reference, TE_gff, full_length_TE_gff, gene_gtf, sorted_bam, is_PE))
    return new_batch_files

def generate_bam(genome_path, output_dir, threads, is_PE=True, **kwargs):
    # 2. 调用 hisat2 将RNA-seq比对到基因组上
    genome_dir = os.path.dirname(genome_path)
    genome_name = os.path.splitext(os.path.basename(genome_path))[0]
    output_sam = output_dir + '/' + genome_name + '.output.sam'
    sorted_bam = output_dir + '/' + genome_name + '.output.sorted.bam'
    hisat2_build = 'cd ' + genome_dir + ' && hisat2-build ' + genome_path + ' ' + genome_name

    conda_prefix = os.getenv("CONDA_PREFIX")  # 获取当前激活的Conda环境根目录
    if not conda_prefix:
        raise RuntimeError("No active Conda environment found.")

    # 1. 调用 trimmomatic 去掉低质量和测序adapter
    if is_PE:
        raw_RNA1 = kwargs['raw_RNA1']
        raw_RNA2 = kwargs['raw_RNA2']

        ILLUMINACLIP_path = os.path.join(conda_prefix, 'share/trimmomatic/adapters/TruSeq3-PE.fa')
        paired_trim_RNA1, paired_trim_RNA2, temp_files = PE_RNA_trim(raw_RNA1, raw_RNA2, ILLUMINACLIP_path, threads)

        os.system(hisat2_build)
        hisat2_align = 'cd ' + genome_dir + ' && hisat2 -x ' + genome_name + ' -1 ' + paired_trim_RNA1 + ' -2 ' + paired_trim_RNA2 + ' -S ' + output_sam + ' -p ' + str(threads)
        os.system(hisat2_align)
    else:
        raw_RNA = kwargs['raw_RNA']
        ILLUMINACLIP_path = os.path.join(conda_prefix, 'share/trimmomatic/adapters/TruSeq3-SE.fa')
        trim_RNA, temp_files = SE_RNA_trim(raw_RNA, ILLUMINACLIP_path, threads)

        os.system(hisat2_build)
        hisat2_align = 'cd ' + genome_dir + ' && hisat2 -x ' + genome_name + ' -U ' + trim_RNA + ' -S ' + output_sam + ' -p ' + str(threads)
        os.system(hisat2_align)

    for tmp_file in temp_files:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)

    # 3. 排序 sam 文件
    sorted_command = 'samtools sort -o ' + sorted_bam + ' ' + output_sam + ' && rm -f ' + output_sam
    os.system(sorted_command)
    return sorted_bam

def PE_RNA_trim(raw_RNA1, raw_RNA2, ILLUMINACLIP_path, threads):
    temp_files = []
    # 1. 先获取文件的名称，以方便后续构建临时文件
    raw_RNA1_dir = os.path.dirname(raw_RNA1)
    raw_RNA1_filename = os.path.basename(raw_RNA1)
    raw_RNA1_name = raw_RNA1_filename.split('.')[0]
    paired_trim_RNA1 = raw_RNA1_dir + '/' + raw_RNA1_name + '.paired_trim.fq.gz'
    unpaired_trim_RNA1 = raw_RNA1_dir + '/' + raw_RNA1_name + '.unpaired_trim.fq.gz'
    temp_files.append(paired_trim_RNA1)
    temp_files.append(unpaired_trim_RNA1)

    raw_RNA2_dir = os.path.dirname(raw_RNA2)
    raw_RNA2_filename = os.path.basename(raw_RNA2)
    raw_RNA2_name = raw_RNA2_filename.split('.')[0]
    paired_trim_RNA2 = raw_RNA2_dir + '/' + raw_RNA2_name + '.paired_trim.fq.gz'
    unpaired_trim_RNA2 = raw_RNA2_dir + '/' + raw_RNA2_name + '.unpaired_trim.fq.gz'
    temp_files.append(paired_trim_RNA2)
    temp_files.append(unpaired_trim_RNA2)

    # 2.调用 trimmomatic 去掉低质量和测序adapter
    trimmomatic_command = 'trimmomatic PE ' + raw_RNA1 + ' ' + raw_RNA2 + ' ' + paired_trim_RNA1 + ' ' + \
                          unpaired_trim_RNA1 + ' ' + paired_trim_RNA2 + ' ' + unpaired_trim_RNA2 + ' ' + \
                          'ILLUMINACLIP:' + ILLUMINACLIP_path + ':2:30:10' + \
                          ' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 ' + '-threads ' + str(threads)
    os.system(trimmomatic_command)
    return paired_trim_RNA1, paired_trim_RNA2, temp_files


def SE_RNA_trim(raw_RNA, ILLUMINACLIP_path, threads):
    temp_files = []
    # 1. 先获取文件的名称，以方便后续构建临时文件
    raw_RNA_dir = os.path.dirname(raw_RNA)
    raw_RNA_filename = os.path.basename(raw_RNA)
    raw_RNA_name = raw_RNA_filename.split('.')[0]
    trim_RNA = raw_RNA_dir + '/' + raw_RNA_name + '.trim.fq.gz'
    temp_files.append(trim_RNA)

    # 2.调用 trimmomatic 去掉低质量和测序adapter
    trimmomatic_command = 'trimmomatic SE ' + raw_RNA + ' ' + trim_RNA + ' ' + \
                          'ILLUMINACLIP:' + ILLUMINACLIP_path + ':2:30:10' + \
                          ' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 ' + '-threads ' + str(threads)
    os.system(trimmomatic_command)
    return trim_RNA, temp_files

def run_featurecounts(output_dir, RNA_tool_dir, sorted_bam, gene_gtf, genome_name, is_PE, recover, log):
    gene_express_count = output_dir + '/' + genome_name + '.count'
    resut_file = gene_express_count
    if not recover or not os.path.exists(resut_file):
        featurecounts_cmd = 'cd ' + output_dir + ' && Rscript ' + RNA_tool_dir + '/run-featurecounts.R' + ' -b ' + sorted_bam + ' -g ' + gene_gtf + ' -o ' + genome_name + \
                            ' --isPairedEnd ' + str(is_PE)
        log.logger.debug(featurecounts_cmd)
        os.system(featurecounts_cmd)
    else:
        log.logger.info(resut_file + ' exists, skip...')
    return gene_express_count

def quantitative_gene(new_batch_files, RNA_tool_dir, output_dir, threads, recover, log):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    for genome_name, reference, TE_gff, full_length_gff_path, gene_gtf, sorted_bam, is_PE in new_batch_files:
        obj = ex.submit(run_featurecounts, output_dir, RNA_tool_dir, sorted_bam, gene_gtf, genome_name, is_PE, recover, log)
        objs.append(obj)
        job_id += 1
    ex.shutdown(wait=True)
    gene_express_counts = []
    for obj in as_completed(objs):
        cur_gene_express_count = obj.result()
        gene_express_counts.append(cur_gene_express_count)

    output_table = output_dir + '/gene_express.table'
    merge_gene_express_table(gene_express_counts, output_table)
    return output_table

def merge_gene_express_table(gene_express_counts, output_table):
    TE_families = set()
    samples = set()
    TE_express = {}
    for cur_count_file in gene_express_counts:
        file_name = os.path.basename(cur_count_file)
        sample_name = file_name.replace('.count', '')
        samples.add(sample_name)
        with open(cur_count_file, 'r') as f_r:
            for i, line in enumerate(f_r):
                if i == 0:
                    continue
                parts = line.split('\t')
                te_name = parts[0]
                te_name = str(te_name).split('_')[-1]
                TE_families.add(te_name)
                counts = parts[1]
                fpkm = "{:.2f}".format(float(parts[2]))
                tpm = "{:.2f}".format(float(parts[3]))
                if te_name not in TE_express:
                    TE_express[te_name] = {}
                cur_TE_express = TE_express[te_name]
                cur_TE_express[sample_name] = counts + ',' + fpkm + ',' + tpm

    with open(output_table, 'w') as f_save:
        header = 'gene_id\t'
        for i, sample in enumerate(samples):
            if i != len(samples) - 1:
                header += sample + '\t'
            else:
                header += sample + '\n'
        f_save.write(header)

        for te_name in TE_express.keys():
            f_save.write(te_name + '\t')
            cur_TE_express = TE_express[te_name]
            for i, sample_name in enumerate(samples):
                if sample_name in cur_TE_express:
                    express_value = cur_TE_express[sample_name]
                else:
                    express_value = 'NA,NA,NA'
                if i != len(samples) - 1:
                    f_save.write(express_value + '\t')
                else:
                    f_save.write(express_value + '\n')

def summary_TEs(genome_info_list, panTE_lib, output_dir, recover, log):
    # 找到 core TEs (出现在100%的基因组)，softcore TEs (出现在80%以上基因组)，dispensable TEs (出现 2个-80%基因组)，private TEs (出现在1个基因组)
    genome_num = len(genome_info_list)
    genome_num_threshold = int(0.8 * genome_num)

    # 获取 TE_name 和 TE_class的对应关系
    te_classes, new_te_contigs = get_TE_class(panTE_lib)

    # 统计 TE 出现的拷贝次数 和 length coverage
    pan_te_fl_infos, pan_te_total_infos = get_panTE_info(genome_info_list, panTE_lib)

    # 生成一个 PAV.tsv 表格，行表示TE family，列表示基因组名称
    generate_panTE_PAV(new_te_contigs, pan_te_fl_infos, output_dir, log)

    # 获取 TE 出现在多少个不同的基因组
    te_fl_occur_genomes, te_occur_genomes = get_te_occur_genomes(new_te_contigs, pan_te_fl_infos, pan_te_total_infos)

    # 将 TE 分成 core, softcore, dispensable, private, unknown
    (core_fl_tes, softcore_fl_tes, dispensable_fl_tes, private_fl_tes, unknown_fl_tes,
     core_tes, softcore_tes, dispensable_tes, private_tes,
     unknown_tes) = get_core_softcore_dispensable_private_uknown_TEs(new_te_contigs, te_fl_occur_genomes,
                                                                     te_occur_genomes, genome_num, genome_num_threshold)

    genome_names = list(pan_te_total_infos.keys())
    all_class_names = list(set(te_classes.values()))

    pdf_files = []
    output_full_length_pdf = output_dir + '/TEs_ratio.full_length.pdf'
    output_pdf = output_dir + '/TEs_ratio.all.pdf'
    draw_four_types_TE_ratio(core_fl_tes, softcore_fl_tes, dispensable_fl_tes, private_fl_tes,
                             core_tes, softcore_tes, dispensable_tes, private_tes, output_full_length_pdf, output_pdf)
    pdf_files.append(output_full_length_pdf)
    pdf_files.append(output_pdf)


    output_full_length_pdf = output_dir + '/TEs_coverage.full_length.pdf'
    output_pdf = output_dir + '/TEs_coverage.all.pdf'
    draw_four_types_TE_coverage(genome_num, pan_te_fl_infos, core_fl_tes, softcore_fl_tes, dispensable_fl_tes,
                                private_fl_tes,
                                pan_te_total_infos, core_tes, softcore_tes, dispensable_tes, private_tes, output_full_length_pdf, output_pdf)
    pdf_files.append(output_full_length_pdf)
    pdf_files.append(output_pdf)


    output_full_length_pdf = output_dir + '/TEclasses_ratio.full_length.pdf'
    output_pdf = output_dir + '/TEclasses_ratio.all.pdf'
    draw_four_types_TE_class_ratio(te_classes, core_fl_tes, softcore_fl_tes, dispensable_fl_tes,
                                  private_fl_tes, unknown_fl_tes, core_tes, softcore_tes, dispensable_tes,
                                  private_tes, unknown_tes, output_full_length_pdf, output_pdf)
    pdf_files.append(output_full_length_pdf)
    pdf_files.append(output_pdf)


    output_full_length_pdf = output_dir + '/TEclasses_coverage.full_length.pdf'
    output_pdf = output_dir + '/TEclasses_coverage.all.pdf'
    draw_four_types_TE_class_coverage(genome_names, all_class_names, te_classes, pan_te_fl_infos, pan_te_total_infos, output_full_length_pdf, output_pdf)
    pdf_files.append(output_full_length_pdf)
    pdf_files.append(output_pdf)

    # 计算每个 genome 上的intact LTR插入时间
    output_pdf = output_dir + '/intact_LTR_insert_time.pdf'
    draw_intact_LTR_insert_time(genome_info_list, output_pdf)
    pdf_files.append(output_pdf)

    # 合并成一个Pdf
    TE_summary_pdf = output_dir + '/TE_summary.pdf'
    merger = PdfMerger()
    for pdf in pdf_files:
        merger.append(pdf)
    merger.write(TE_summary_pdf)
    merger.close()

    if os.path.exists(TE_summary_pdf):
        for pdf in pdf_files:
            os.remove(pdf)

def draw_four_types_TE_coverage(genome_num, pan_te_fl_infos, core_fl_tes, softcore_fl_tes, dispensable_fl_tes, private_fl_tes,
                                pan_te_total_infos, core_tes, softcore_tes, dispensable_tes, private_tes, output_full_length_pdf, output_pdf):
    # 计算四种类型的全长 TE 分别覆盖每个基因组的比例，并画出所有基因组上的饼图比例
    fl_genome_coverage_core = {}
    fl_genome_coverage_softcore = {}
    fl_genome_coverage_dispensable = {}
    fl_genome_coverage_private = {}
    for genome_name in pan_te_fl_infos.keys():
        if genome_name not in fl_genome_coverage_core:
            fl_genome_coverage_core[genome_name] = 0
        cur_fl_coverage_core = fl_genome_coverage_core[genome_name]
        if genome_name not in fl_genome_coverage_softcore:
            fl_genome_coverage_softcore[genome_name] = 0
        cur_fl_coverage_softcore = fl_genome_coverage_softcore[genome_name]
        if genome_name not in fl_genome_coverage_dispensable:
            fl_genome_coverage_dispensable[genome_name] = 0
        cur_fl_coverage_dispensable = fl_genome_coverage_dispensable[genome_name]
        if genome_name not in fl_genome_coverage_private:
            fl_genome_coverage_private[genome_name] = 0
        cur_fl_coverage_private = fl_genome_coverage_private[genome_name]

        for te_name in core_fl_tes.keys():
            te_fl_infos = pan_te_fl_infos[genome_name]
            if te_name in te_fl_infos:
                cur_copy_num, cur_fl_length = te_fl_infos[te_name]
                fl_genome_coverage_core[genome_name] = cur_fl_coverage_core + cur_fl_length

        for te_name in softcore_fl_tes.keys():
            te_fl_infos = pan_te_fl_infos[genome_name]
            if te_name in te_fl_infos:
                cur_copy_num, cur_fl_length = te_fl_infos[te_name]
                fl_genome_coverage_softcore[genome_name] = cur_fl_coverage_softcore + cur_fl_length

        for te_name in dispensable_fl_tes.keys():
            te_fl_infos = pan_te_fl_infos[genome_name]
            if te_name in te_fl_infos:
                cur_copy_num, cur_fl_length = te_fl_infos[te_name]
                fl_genome_coverage_dispensable[genome_name] = cur_fl_coverage_dispensable + cur_fl_length

        for te_name in private_fl_tes.keys():
            te_fl_infos = pan_te_fl_infos[genome_name]
            if te_name in te_fl_infos:
                cur_copy_num, cur_fl_length = te_fl_infos[te_name]
                fl_genome_coverage_private[genome_name] = cur_fl_coverage_private + cur_fl_length

    # print(fl_genome_coverage_core)
    # print(fl_genome_coverage_softcore)
    # print(fl_genome_coverage_dispensable)
    # print(fl_genome_coverage_private)

    sizes_list = []
    labels_list = []
    genome_names = []
    title = 'Full length TE Coverage'
    for genome_name in pan_te_fl_infos.keys():
        cur_sizes = []
        cur_labels = []
        coverage_core = fl_genome_coverage_core[genome_name]
        cur_sizes.append(coverage_core)
        cur_labels.append('Core TEs')
        coverage_softcore = fl_genome_coverage_softcore[genome_name]
        cur_sizes.append(coverage_softcore)
        cur_labels.append('Softcore TEs')
        coverage_dispensable = fl_genome_coverage_dispensable[genome_name]
        cur_sizes.append(coverage_dispensable)
        cur_labels.append('Dispensable TEs')
        coverage_private = fl_genome_coverage_private[genome_name]
        cur_sizes.append(coverage_private)
        cur_labels.append('Private TEs')
        sizes_list.append(cur_sizes)
        labels_list.append(cur_labels)
        genome_names.append(genome_name)
    draw_multiple_pie(genome_num, sizes_list, labels_list, genome_names, title, output_full_length_pdf)
    # =====================================================================
    # 计算四种类型的 TE 分别覆盖每个基因组的比例，并画出每个基因组的饼图
    genome_coverage_core = {}
    genome_coverage_softcore = {}
    genome_coverage_dispensable = {}
    genome_coverage_private = {}
    for genome_name in pan_te_total_infos.keys():
        if genome_name not in genome_coverage_core:
            genome_coverage_core[genome_name] = 0
        cur_coverage_core = genome_coverage_core[genome_name]
        if genome_name not in genome_coverage_softcore:
            genome_coverage_softcore[genome_name] = 0
        cur_coverage_softcore = genome_coverage_softcore[genome_name]
        if genome_name not in genome_coverage_dispensable:
            genome_coverage_dispensable[genome_name] = 0
        cur_coverage_dispensable = genome_coverage_dispensable[genome_name]
        if genome_name not in genome_coverage_private:
            genome_coverage_private[genome_name] = 0
        cur_coverage_private = genome_coverage_private[genome_name]

        for te_name in core_tes.keys():
            te_total_infos = pan_te_total_infos[genome_name]
            if te_name in te_total_infos:
                cur_copy_num, cur_length = te_total_infos[te_name]
                genome_coverage_core[genome_name] = cur_coverage_core + cur_length

        for te_name in softcore_tes.keys():
            te_total_infos = pan_te_total_infos[genome_name]
            if te_name in te_total_infos:
                cur_copy_num, cur_length = te_total_infos[te_name]
                genome_coverage_softcore[genome_name] = cur_coverage_softcore + cur_length

        for te_name in dispensable_tes.keys():
            te_total_infos = pan_te_total_infos[genome_name]
            if te_name in te_total_infos:
                cur_copy_num, cur_length = te_total_infos[te_name]
                genome_coverage_dispensable[genome_name] = cur_coverage_dispensable + cur_length

        for te_name in private_tes.keys():
            te_total_infos = pan_te_total_infos[genome_name]
            if te_name in te_total_infos:
                cur_copy_num, cur_length = te_total_infos[te_name]
                genome_coverage_private[genome_name] = cur_coverage_private + cur_length

    # print(genome_coverage_core)
    # print(genome_coverage_softcore)
    # print(genome_coverage_dispensable)
    # print(genome_coverage_private)

    sizes_list = []
    labels_list = []
    genome_names = []
    title = 'TE Coverage'
    for genome_name in pan_te_total_infos.keys():
        cur_sizes = []
        cur_labels = []
        coverage_core = genome_coverage_core[genome_name]
        cur_sizes.append(coverage_core)
        cur_labels.append('Core TEs')
        coverage_softcore = genome_coverage_softcore[genome_name]
        cur_sizes.append(coverage_softcore)
        cur_labels.append('Softcore TEs')
        coverage_dispensable = genome_coverage_dispensable[genome_name]
        cur_sizes.append(coverage_dispensable)
        cur_labels.append('Dispensable TEs')
        coverage_private = genome_coverage_private[genome_name]
        cur_sizes.append(coverage_private)
        cur_labels.append('Private TEs')
        sizes_list.append(cur_sizes)
        labels_list.append(cur_labels)
        genome_names.append(genome_name)
    draw_multiple_pie(genome_num, sizes_list, labels_list, genome_names, title, output_pdf)
    return genome_names

def draw_four_types_TE_class_ratio(te_classes, core_fl_tes, softcore_fl_tes, dispensable_fl_tes,
                                   private_fl_tes, unknown_fl_tes, core_tes, softcore_tes, dispensable_tes,
                                   private_tes, unknown_tes, output_full_length_pdf, output_pdf):
    # 统计全长 TEs 中，四种类型的TEs，每一种的 TE class 数量
    core_te_class_num = {}
    for te_name in core_fl_tes.keys():
        cur_class_name = te_classes[te_name]
        if cur_class_name not in core_te_class_num:
            core_te_class_num[cur_class_name] = 0
        cur_te_class_num = core_te_class_num[cur_class_name]
        core_te_class_num[cur_class_name] = cur_te_class_num + 1

    softcore_te_class_num = {}
    for te_name in softcore_fl_tes.keys():
        cur_class_name = te_classes[te_name]
        if cur_class_name not in softcore_te_class_num:
            softcore_te_class_num[cur_class_name] = 0
        cur_te_class_num = softcore_te_class_num[cur_class_name]
        softcore_te_class_num[cur_class_name] = cur_te_class_num + 1

    dispensable_te_class_num = {}
    for te_name in dispensable_fl_tes.keys():
        cur_class_name = te_classes[te_name]
        if cur_class_name not in dispensable_te_class_num:
            dispensable_te_class_num[cur_class_name] = 0
        cur_te_class_num = dispensable_te_class_num[cur_class_name]
        dispensable_te_class_num[cur_class_name] = cur_te_class_num + 1

    private_te_class_num = {}
    for te_name in private_fl_tes.keys():
        cur_class_name = te_classes[te_name]
        if cur_class_name not in private_te_class_num:
            private_te_class_num[cur_class_name] = 0
        cur_te_class_num = private_te_class_num[cur_class_name]
        private_te_class_num[cur_class_name] = cur_te_class_num + 1

    unknown_te_class_num = {}
    for te_name in unknown_fl_tes.keys():
        cur_class_name = te_classes[te_name]
        if cur_class_name not in unknown_te_class_num:
            unknown_te_class_num[cur_class_name] = 0
        cur_te_class_num = unknown_te_class_num[cur_class_name]
        unknown_te_class_num[cur_class_name] = cur_te_class_num + 1

    # print(core_te_class_num)
    # print(softcore_te_class_num)
    # print(dispensable_te_class_num)
    # print(private_te_class_num)
    # print(unknown_te_class_num)

    all_class_names = list(set(te_classes.values()))
    labels = ['Core TEs', 'Softcore TEs', 'Dispensable TEs', 'Private TEs']
    data_list = []
    for cur_class_name in all_class_names:
        cur_data_list = []
        if cur_class_name in core_te_class_num:
            core_num = core_te_class_num[cur_class_name]
        else:
            core_num = 0
        cur_data_list.append(core_num)
        if cur_class_name in softcore_te_class_num:
            softcore_num = softcore_te_class_num[cur_class_name]
        else:
            softcore_num = 0
        cur_data_list.append(softcore_num)
        if cur_class_name in dispensable_te_class_num:
            dispensable_num = dispensable_te_class_num[cur_class_name]
        else:
            dispensable_num = 0
        cur_data_list.append(dispensable_num)
        if cur_class_name in private_te_class_num:
            private_num = private_te_class_num[cur_class_name]
        else:
            private_num = 0
        cur_data_list.append(private_num)
        data_list.append(cur_data_list)
    draw_stacked_bar_chart(labels, data_list, all_class_names, 'Full length TE Classes Ratio', output_full_length_pdf)

    # 统计所有 TEs 中，四种类型的TEs，每一种的 TE class 数量
    core_te_class_num = {}
    for te_name in core_tes.keys():
        cur_class_name = te_classes[te_name]
        if cur_class_name not in core_te_class_num:
            core_te_class_num[cur_class_name] = 0
        cur_te_class_num = core_te_class_num[cur_class_name]
        core_te_class_num[cur_class_name] = cur_te_class_num + 1

    softcore_te_class_num = {}
    for te_name in softcore_tes.keys():
        cur_class_name = te_classes[te_name]
        if cur_class_name not in softcore_te_class_num:
            softcore_te_class_num[cur_class_name] = 0
        cur_te_class_num = softcore_te_class_num[cur_class_name]
        softcore_te_class_num[cur_class_name] = cur_te_class_num + 1

    dispensable_te_class_num = {}
    for te_name in dispensable_tes.keys():
        cur_class_name = te_classes[te_name]
        if cur_class_name not in dispensable_te_class_num:
            dispensable_te_class_num[cur_class_name] = 0
        cur_te_class_num = dispensable_te_class_num[cur_class_name]
        dispensable_te_class_num[cur_class_name] = cur_te_class_num + 1

    private_te_class_num = {}
    for te_name in private_tes.keys():
        cur_class_name = te_classes[te_name]
        if cur_class_name not in private_te_class_num:
            private_te_class_num[cur_class_name] = 0
        cur_te_class_num = private_te_class_num[cur_class_name]
        private_te_class_num[cur_class_name] = cur_te_class_num + 1

    unknown_te_class_num = {}
    for te_name in unknown_tes.keys():
        cur_class_name = te_classes[te_name]
        if cur_class_name not in unknown_te_class_num:
            unknown_te_class_num[cur_class_name] = 0
        cur_te_class_num = unknown_te_class_num[cur_class_name]
        unknown_te_class_num[cur_class_name] = cur_te_class_num + 1

    # print(core_te_class_num)
    # print(softcore_te_class_num)
    # print(dispensable_te_class_num)
    # print(private_te_class_num)
    # print(unknown_te_class_num)

    all_class_names = list(set(te_classes.values()))
    labels = ['Core TEs', 'Softcore TEs', 'Dispensable TEs', 'Private TEs']
    data_list = []
    for cur_class_name in all_class_names:
        cur_data_list = []
        if cur_class_name in core_te_class_num:
            core_num = core_te_class_num[cur_class_name]
        else:
            core_num = 0
        cur_data_list.append(core_num)
        if cur_class_name in softcore_te_class_num:
            softcore_num = softcore_te_class_num[cur_class_name]
        else:
            softcore_num = 0
        cur_data_list.append(softcore_num)
        if cur_class_name in dispensable_te_class_num:
            dispensable_num = dispensable_te_class_num[cur_class_name]
        else:
            dispensable_num = 0
        cur_data_list.append(dispensable_num)
        if cur_class_name in private_te_class_num:
            private_num = private_te_class_num[cur_class_name]
        else:
            private_num = 0
        cur_data_list.append(private_num)
        data_list.append(cur_data_list)
    draw_stacked_bar_chart(labels, data_list, all_class_names, 'TE Classes Ratio', output_pdf)
    return all_class_names

def draw_four_types_TE_class_coverage(genome_names, all_class_names, te_classes, pan_te_fl_infos, pan_te_total_infos, output_full_length_pdf, output_pdf):
    # 统计每个基因组的各种类型的全长TE的占比,绘制横向堆叠柱状图
    genome_fl_te_class_coverages = {}
    for genome_name in pan_te_fl_infos.keys():
        if genome_name not in genome_fl_te_class_coverages:
            genome_fl_te_class_coverages[genome_name] = {}
        cur_genome_fl_te_class_coverage = genome_fl_te_class_coverages[genome_name]
        te_fl_infos = pan_te_fl_infos[genome_name]
        for te_name in te_fl_infos.keys():
            cur_class_name = te_classes[te_name]
            cur_copy_num, cur_fl_length = te_fl_infos[te_name]

            if cur_class_name not in cur_genome_fl_te_class_coverage:
                cur_genome_fl_te_class_coverage[cur_class_name] = 0
            cur_class_coverage = cur_genome_fl_te_class_coverage[cur_class_name]
            cur_genome_fl_te_class_coverage[cur_class_name] = cur_class_coverage + cur_fl_length
    # print(genome_fl_te_class_coverages)
    labels = genome_names
    data_list = []
    for cur_class_name in all_class_names:
        cur_data_list = []
        for genome_name in genome_names:
            cur_genome_fl_te_class_coverage = genome_fl_te_class_coverages[genome_name]
            if cur_class_name in cur_genome_fl_te_class_coverage:
                cur_class_coverage = cur_genome_fl_te_class_coverage[cur_class_name]
            else:
                cur_class_coverage = 0
            cur_data_list.append(cur_class_coverage / 1000000)
        data_list.append(cur_data_list)
    draw_horizontal_stacked_bar_chart(labels, data_list, all_class_names, 'Full length TE Classes Coverage', output_full_length_pdf)

    # 统计每个基因组的各种类型的TE的占比,绘制横向堆叠柱状图
    genome_te_class_coverages = {}
    for genome_name in pan_te_total_infos.keys():
        if genome_name not in genome_te_class_coverages:
            genome_te_class_coverages[genome_name] = {}
        cur_genome_te_class_coverage = genome_te_class_coverages[genome_name]
        te_total_infos = pan_te_total_infos[genome_name]
        for te_name in te_total_infos.keys():
            cur_class_name = te_classes[te_name]
            cur_copy_num, cur_length = te_total_infos[te_name]

            if cur_class_name not in cur_genome_te_class_coverage:
                cur_genome_te_class_coverage[cur_class_name] = 0
            cur_class_coverage = cur_genome_te_class_coverage[cur_class_name]
            cur_genome_te_class_coverage[cur_class_name] = cur_class_coverage + cur_length
    # print(genome_te_class_coverages)
    labels = genome_names
    data_list = []
    for cur_class_name in all_class_names:
        cur_data_list = []
        for genome_name in genome_names:
            cur_genome_te_class_coverage = genome_te_class_coverages[genome_name]
            if cur_class_name in cur_genome_te_class_coverage:
                cur_class_coverage = cur_genome_te_class_coverage[cur_class_name]
            else:
                cur_class_coverage = 0
            cur_data_list.append(cur_class_coverage / 1000000)
        data_list.append(cur_data_list)
    draw_horizontal_stacked_bar_chart(labels, data_list, all_class_names, 'TE Classes Coverage', output_pdf)

def get_panTE_full_length_annotation(pan_te_full_length_annotations, core_fl_tes, softcore_fl_tes, dispensable_fl_tes, private_fl_tes, output_dir, recover, log):
    for genome_name in pan_te_full_length_annotations.keys():
        full_length_gff_path = output_dir + '/' + genome_name + '.full_length.gff'
        resut_file = full_length_gff_path
        if not recover or not os.path.exists(resut_file):
            full_length_annotations = pan_te_full_length_annotations[genome_name]
            TE_gff = output_dir + '/' + genome_name + '.temp.gff'
            intact_count = 0
            with open(TE_gff, 'w') as f_save:
                for query_name in full_length_annotations.keys():
                    if query_name in core_fl_tes:
                        type = 'Core_TEs'
                    elif query_name in softcore_fl_tes:
                        type = 'Softcore_TEs'
                    elif query_name in dispensable_fl_tes:
                        type = 'Dispensable_TEs'
                    elif query_name in private_fl_tes:
                        type = 'Private_TEs'
                    else:
                        type = 'Unknown_TEs'
                    for copy_annotation in full_length_annotations[query_name]:
                        classification = copy_annotation[0]
                        chr_name = copy_annotation[1]
                        chr_start = str(copy_annotation[2] + 1)
                        chr_end = str(copy_annotation[3])
                        intact_count += 1
                        update_annotation = 'id=te_intact_' + str(
                            intact_count) + ';name=' + query_name + ';classification=' + classification + ';type=' + type
                        f_save.write(
                            chr_name + '\t' + 'HiTE' + '\t' + classification + '\t' + chr_start + '\t' + chr_end + '\t' + '.\t' +
                            str(copy_annotation[4]) + '\t' + '.\t' + update_annotation + '\n')

            os.system('sort -k1,1 -k4n ' + TE_gff + ' > ' + full_length_gff_path)
            gff_lines = []
            with open(full_length_gff_path, 'r') as f_r:
                for line in f_r:
                    if line.startswith('#'):
                        continue
                    gff_lines.append(line)

            date = datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
            header = (
                "##gff-version 3\n"
                f"##date {date}\n"
            )
            with open(full_length_gff_path, "w") as gff_file:
                gff_file.write(header)
                for line in gff_lines:
                    gff_file.write(line)
            os.remove(TE_gff)
        else:
            if log is not None:
                log.logger.info(resut_file + ' exists, skip...')

def draw_intact_LTR_insert_time(genome_info_list, output_pdf):
    # 初始化列表，用于存储每个 genome_name 的插入时间数据
    all_data = []
    for genome_info in genome_info_list:
        genome_name = genome_info["genome_name"]
        intact_LTR_list = genome_info["intact_LTR_list"]

        # 读取 CSV 文件并跳过 # 注释行
        df = pd.read_csv(intact_LTR_list, sep='\t', comment='#',
                         names=['LTR_loc', 'Status', 'Motif', 'TSD', 'lTSD', 'rTSD', 'Internal', 'Identity', 'Strand', 'Undef1', 'Undef2', 'Classification', 'Insertion_Time'])

        df['Insertion_Time'] = pd.to_numeric(df['Insertion_Time'], errors='coerce')

        # 计算插入时间并转换为百万年，忽略 NaN
        df['Insertion_Time'] = df['Insertion_Time'] / 1_000_000

        # 为每个 genome_name 添加 Copia 和 Gypsy 的插入时间及分类信息
        for classification in ['LTR/Copia', 'LTR/Gypsy']:
            filtered_df = df[df['Classification'] == classification]
            for insertion_time in filtered_df['Insertion_Time']:
                all_data.append([genome_name, insertion_time, classification])

    # 将数据转换为 DataFrame
    all_data_df = pd.DataFrame(all_data, columns=['Genome', 'Insertion_Time', 'Classification'])

    # 绘制箱线图，使用 hue 区分 LTR/Copia 和 LTR/Gypsy
    plt.figure(figsize=(15, 15))
    medianprops = {'color': 'black', 'linewidth': 2.5}  # 中位数线样式
    ax = sns.boxplot(x='Genome', y='Insertion_Time', hue='Classification', data=all_data_df,
                     showfliers=False, medianprops=medianprops)

    # 图形调整
    plt.ylabel('Million years ago (Mya)')
    plt.title('Insertion Time for LTR/Copia and LTR/Gypsy across Genomes')
    plt.xticks(rotation=90)  # 如果 genome_name 太长，可以旋转 x 轴标签
    plt.tight_layout()

    # 保存图表
    plt.savefig(output_pdf)
    plt.close()

def get_TE_class(panTE_lib):
    te_classes = {}
    te_names, te_contigs = read_fasta(panTE_lib)
    new_te_contigs = {}
    for name in te_names:
        parts = name.split('#')
        raw_name = parts[0]
        class_name = parts[1]
        new_te_contigs[raw_name] = te_contigs[name]
        te_classes[raw_name] = class_name
    return te_classes, new_te_contigs

def get_panTE_info(genome_info_list, panTE_lib):
    pan_te_fl_infos = {}
    pan_te_total_infos = {}
    for genome_info in genome_info_list:
        genome_name = genome_info["genome_name"]
        TE_gff = genome_info["TE_gff"]
        full_length_copies_path = genome_info["full_length_copies"]
        te_fl_infos, te_total_infos = get_copy_and_length_from_gff(TE_gff, full_length_copies_path, panTE_lib)
        pan_te_fl_infos[genome_name] = te_fl_infos
        pan_te_total_infos[genome_name] = te_total_infos
    return pan_te_fl_infos, pan_te_total_infos

def get_te_occur_genomes(new_te_contigs, pan_te_fl_infos, pan_te_total_infos):
    te_fl_occur_genomes = {}
    te_occur_genomes = {}
    for te_name in new_te_contigs.keys():
        for genome_name in pan_te_total_infos.keys():
            te_fl_infos = pan_te_fl_infos[genome_name]
            te_total_infos = pan_te_total_infos[genome_name]
            if te_name in te_fl_infos:
                if te_name not in te_fl_occur_genomes:
                    te_fl_occur_genomes[te_name] = 0
                occur_genome_count = te_fl_occur_genomes[te_name]
                te_fl_occur_genomes[te_name] = occur_genome_count + 1

            if te_name in te_total_infos:
                if te_name not in te_occur_genomes:
                    te_occur_genomes[te_name] = 0
                occur_genome_count = te_occur_genomes[te_name]
                te_occur_genomes[te_name] = occur_genome_count + 1
    return te_fl_occur_genomes, te_occur_genomes

def get_core_softcore_dispensable_private_uknown_TEs(new_te_contigs, te_fl_occur_genomes, te_occur_genomes, genome_num, genome_num_threshold):
    core_fl_tes = {}
    softcore_fl_tes = {}
    dispensable_fl_tes = {}
    private_fl_tes = {}
    unknown_fl_tes = {}

    core_tes = {}
    softcore_tes = {}
    dispensable_tes = {}
    private_tes = {}
    unknown_tes = {}
    for te_name in new_te_contigs.keys():
        if te_name in te_fl_occur_genomes:
            cur_occur_num = te_fl_occur_genomes[te_name]
            if cur_occur_num == genome_num:
                core_fl_tes[te_name] = cur_occur_num
            elif genome_num_threshold <= cur_occur_num < genome_num:
                softcore_fl_tes[te_name] = cur_occur_num
            elif 2 <= cur_occur_num < genome_num_threshold:
                dispensable_fl_tes[te_name] = cur_occur_num
            elif cur_occur_num == 1:
                private_fl_tes[te_name] = cur_occur_num
            else:
                unknown_fl_tes[te_name] = cur_occur_num
        else:
            unknown_fl_tes[te_name] = 0

        if te_name in te_occur_genomes:
            cur_occur_num = te_occur_genomes[te_name]
            if cur_occur_num == genome_num:
                core_tes[te_name] = cur_occur_num
            elif genome_num_threshold <= cur_occur_num < genome_num:
                softcore_tes[te_name] = cur_occur_num
            elif 2 <= cur_occur_num < genome_num_threshold:
                dispensable_tes[te_name] = cur_occur_num
            elif cur_occur_num == 1:
                private_tes[te_name] = cur_occur_num
            else:
                unknown_tes[te_name] = cur_occur_num
        else:
            unknown_tes[te_name] = 0



    # print('core_fl_tes:' + str(core_fl_tes))
    # print('softcore_fl_tes:' + str(softcore_fl_tes))
    # print('dispensable_fl_tes:' + str(dispensable_fl_tes))
    # print('private_fl_tes:' + str(private_fl_tes))
    # print('unknown_fl_tes:' + str(unknown_fl_tes))
    # print('===========================================================================')
    # print('core_tes:' + str(core_tes))
    # print('softcore_tes:' + str(softcore_tes))
    # print('dispensable_tes:' + str(dispensable_tes))
    # print('private_tes:' + str(private_tes))
    # print('unknown_tes:' + str(unknown_tes))
    return (core_fl_tes, softcore_fl_tes, dispensable_fl_tes, private_fl_tes, unknown_fl_tes,
            core_tes, softcore_tes, dispensable_tes, private_tes, unknown_tes)

def draw_four_types_TE_ratio(core_fl_tes, softcore_fl_tes, dispensable_fl_tes, private_fl_tes,
                             core_tes, softcore_tes, dispensable_tes, private_tes, output_full_length_pdf, output_pdf):
    # 统计四种类型的全长TE数量 占 总体TE数量的比例，并画出饼图
    core_fl_tes_num = len(core_fl_tes)
    softcore_fl_tes_num = len(softcore_fl_tes)
    dispensable_fl_tes_num = len(dispensable_fl_tes)
    private_fl_tes_num = len(private_fl_tes)
    total_num = core_fl_tes_num + softcore_fl_tes_num + dispensable_fl_tes_num + private_fl_tes_num
    labels = ['Core TEs', 'Softcore TEs', 'Dispensable TEs', 'Private TEs']
    sizes = [float(core_fl_tes_num) * 100 / total_num,
             float(softcore_fl_tes_num) * 100 / total_num,
             float(dispensable_fl_tes_num) * 100 / total_num,
             float(private_fl_tes_num) * 100 / total_num]
    draw_pie(sizes, labels, output_full_length_pdf, title='Full length TEs Ratio', is_percent=True)

    core_tes_num = len(core_tes)
    softcore_tes_num = len(softcore_tes)
    dispensable_tes_num = len(dispensable_tes)
    private_tes_num = len(private_tes)
    total_num = core_tes_num + softcore_tes_num + dispensable_tes_num + private_tes_num
    labels = ['Core TEs', 'Softcore TEs', 'Dispensable TEs', 'Private TEs']
    sizes = [float(core_tes_num) * 100 / total_num,
             float(softcore_tes_num) * 100 / total_num,
             float(dispensable_tes_num) * 100 / total_num,
             float(private_tes_num) * 100 / total_num]
    draw_pie(sizes, labels, output_pdf, title='TEs Ratio', is_percent=True)

def draw_multiple_pie(sub_num, sizes_list, labels_list, genome_names, title, output_pdf):
    # colors = ['#D53E4F', '#FC8D59', '#3288BD', '#FEE08B', '#FFFFBF', '#E6F598', '#ABDDA4', '#66C2A5', '#3288BD',
    #           '#4A90E2', '#E94E77', '#F5A623', '#7ED321', '#BD10E0', '#003366', '#D0D0D0', '#00274D', '#F7E04C',
    #           '#6B8E23', '#0033A0', '#E50000', '#228B22', '#C2B280', '#5B3A29', '#1E90FF', '#FF4500']
    # colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue']

    colors = [
        '#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3',
        '#ff69b4', '#00fa9a', '#add8e6', '#ff6347', '#dcdcdc', '#ff4500', '#ffb6c1', '#87ceeb',
        '#b0e57c', '#ff8c00', '#f08080', '#ffb000', '#90ee90', '#20b2aa', '#b0c4de', '#d8bfd8',
        '#f5deb3', '#ff1493', '#40e0d0', '#ff6347', '#ffb6c1', '#ff4500'
    ]

    k = 3 # 一行最多显示3个子图
    num_rows = (sub_num + k - 1) // k  # 计算需要多少行

    # 创建图形和子图
    fig, axs = plt.subplots(num_rows, k, figsize=(5 * k, 5 * num_rows))
    axs = axs.flatten()  # 展平二维数组以方便访问

    for i in range(sub_num):
        wedges, texts, autotexts = axs[i].pie(sizes_list[i], colors=colors, autopct='%1.1f%%', startangle=140)
        axs[i].set_title(genome_names[i])
        axs[i].legend(wedges, labels_list[i], loc="best", fontsize='small')

    # 隐藏多余的子图
    for j in range(sub_num, len(axs)):
        axs[j].axis('off')

    fig.suptitle(title, y=1, fontweight='bold')
    # 调整布局
    plt.tight_layout(rect=[0, 0, 0.98, 0.95])
    plt.savefig(output_pdf)

    # 显示图形
    # plt.show()

def draw_stacked_bar_chart(labels, data_list, all_class_names, title, output_pdf):
    # colors = ['#D53E4F', '#FC8D59', '#3288BD', '#FEE08B', '#FFFFBF', '#E6F598', '#ABDDA4', '#66C2A5', '#3288BD',
    #           '#4A90E2', '#E94E77', '#F5A623', '#7ED321', '#BD10E0', '#003366', '#D0D0D0', '#00274D', '#F7E04C',
    #           '#6B8E23', '#0033A0', '#E50000', '#228B22', '#C2B280', '#5B3A29', '#1E90FF', '#FF4500']

    # colors = sns.color_palette("husl", 28)

    colors = [
        '#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3',
        '#ff69b4', '#00fa9a', '#add8e6', '#ff6347', '#dcdcdc', '#ff4500', '#ffb6c1', '#87ceeb',
        '#b0e57c', '#ff8c00', '#f08080', '#ffb000', '#90ee90', '#20b2aa', '#b0c4de', '#d8bfd8',
        '#f5deb3', '#ff1493', '#40e0d0', '#ff6347', '#ffb6c1', '#ff4500'
    ]

    # 设置柱子的宽度
    bar_width = 0.5

    # 设置x轴的位置
    x = np.arange(len(labels))

    # 绘制柱状图
    fig, ax = plt.subplots(figsize=(10, 7))

    bottom_data = np.zeros(len(x))
    for i in range(len(data_list)):
        # 绘制第一组数据
        bars1 = ax.bar(x, data_list[i], bar_width, bottom=bottom_data, label=all_class_names[i], color=colors[i])
        bottom_data += np.array(data_list[i])

    # 添加标签
    ax.set_ylabel('PanTE family number')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    plt.title(title, pad=20, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_pdf)
    # 显示图形
    # plt.show()

def draw_horizontal_stacked_bar_chart(labels, data_list, all_class_names, title, output_pdf):
    # colors = ['#D53E4F', '#FC8D59', '#3288BD', '#FEE08B', '#FFFFBF', '#E6F598', '#ABDDA4', '#66C2A5', '#3288BD',
    #           '#4A90E2', '#E94E77', '#F5A623', '#7ED321', '#BD10E0', '#003366', '#D0D0D0', '#00274D', '#F7E04C',
    #           '#6B8E23', '#0033A0', '#E50000', '#228B22', '#C2B280', '#5B3A29', '#1E90FF', '#FF4500']

    # colors = sns.color_palette("husl", 28)

    colors = [
        '#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3',
        '#ff69b4', '#00fa9a', '#add8e6', '#ff6347', '#dcdcdc', '#ff4500', '#ffb6c1', '#87ceeb',
        '#b0e57c', '#ff8c00', '#f08080', '#ffb000', '#90ee90', '#20b2aa', '#b0c4de', '#d8bfd8',
        '#f5deb3', '#ff1493', '#40e0d0', '#ff6347', '#ffb6c1', '#ff4500'
    ]

    # 示例数据
    y = np.arange(len(labels))  # y 轴位置（条形图的类别）

    fig, ax = plt.subplots(figsize=(15, 7))

    # 初始化 left_data 为与数据长度相同的零数组
    left_data = np.zeros(len(y))

    for i in range(len(data_list)):
        # 绘制每组数据的条形图
        bars = ax.barh(y, data_list[i], height=0.5, left=left_data, label=all_class_names[i], color=colors[i])
        # 更新 left_data，用于堆叠效果
        left_data += np.array(data_list[i])

    # 添加图例和标题
    ax.set_xlabel('Length (Mb)')
    ax.set_ylabel('Genomes')
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.legend(loc='best')
    plt.title(title, pad=20, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_pdf)
    # plt.show()


def get_full_length_copies_from_gff_v1(TE_lib, reference, gff_path, full_length_threshold, te_classes):
    ref_names, ref_contigs = read_fasta(reference)

    query_names, query_contigs = read_fasta(TE_lib)
    new_query_contigs = {}
    for name in query_names:
        new_query_contigs[name.split('#')[0]] = query_contigs[name]
    query_contigs = new_query_contigs

    query_records = {}
    with open(gff_path, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            query_name = parts[8].split(' ')[1].replace('"', '').split(':')[1]
            subject_name = parts[0]
            info_parts = parts[8].split(' ')
            q_start = int(info_parts[2])
            q_end = int(info_parts[3])
            s_start = int(parts[3])
            s_end = int(parts[4])
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))

    longest_repeats = FMEA_new(query_contigs, query_records, full_length_threshold)

    if str(query_name).__contains__('Helitron'):
        flanking_len = 5
    else:
        flanking_len = 50

    full_length_annotations = {}
    full_length_copies = {}
    flank_full_length_copies = {}
    for query_name in longest_repeats.keys():
        cur_longest_repeats = longest_repeats[query_name]
        query_len = len(query_contigs[query_name])
        annotations = []
        query_copies = {}
        flank_query_copies = {}
        for repeat in cur_longest_repeats:
            cur_query_len = abs(repeat[2] - repeat[1])
            if cur_query_len < full_length_threshold * query_len:
                continue
            # Subject
            subject_name = repeat[3]
            subject_chr_start = 0
            if repeat[4] > repeat[5]:
                direct = '-'
                old_subject_start_pos = repeat[5]
                old_subject_end_pos = repeat[4]
            else:
                direct = '+'
                old_subject_start_pos = repeat[4]
                old_subject_end_pos = repeat[5]
            subject_start_pos = subject_chr_start + old_subject_start_pos
            subject_end_pos = subject_chr_start + old_subject_end_pos
            subject_pos = subject_name + ':' + str(subject_start_pos) + '-' + str(subject_end_pos)
            subject_seq = ref_contigs[subject_name][subject_start_pos: subject_end_pos]
            flank_subject_seq = ref_contigs[subject_name][subject_start_pos - flanking_len: subject_end_pos + flanking_len]
            query_copies[subject_pos] = subject_seq
            flank_query_copies[subject_pos] = flank_subject_seq
            annotations.append((te_classes[query_name], subject_name, subject_start_pos, subject_end_pos, direct))
        full_length_copies[query_name] = query_copies
        flank_full_length_copies[query_name] = flank_query_copies
        full_length_annotations[query_name] = annotations
    return full_length_copies, flank_full_length_copies, full_length_annotations

def get_copy_and_length_from_gff(gff_path, full_length_copies_path, panTE_lib):
    with open(full_length_copies_path, 'r') as f:
        full_length_copies = json.load(f)

    te_names, te_contigs = read_fasta(panTE_lib)
    new_te_contigs = {}
    for name in te_names:
        raw_name = name.split('#')[0]
        new_te_contigs[raw_name] = te_contigs[name]

    te_fl_infos = {}
    for te_name in full_length_copies.keys():
        te_cons_len = len(new_te_contigs[te_name])
        query_copies = full_length_copies[te_name]
        for subject_pos in query_copies.keys():
            subject_seq = query_copies[subject_pos]
            if len(subject_seq) / te_cons_len > 0.95:
                if te_name not in te_fl_infos:
                    te_fl_infos[te_name] = (0, 0)
                cur_copy_num, cur_fl_length = te_fl_infos[te_name]
                cur_copy_num += 1
                cur_fl_length += len(subject_seq)
                te_fl_infos[te_name] = (cur_copy_num, cur_fl_length)

    te_total_infos = {}
    with open(gff_path, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            chr_start = int(parts[3])
            chr_end = int(parts[4])
            te_info = parts[8]
            te_parts = te_info.split(' ')
            te_name = te_parts[1].replace('"', '').split(':')[1]
            if te_name not in te_total_infos:
                te_total_infos[te_name] = (0, 0)
            cur_copy_num, cur_total_length = te_total_infos[te_name]
            cur_copy_num += 1
            cur_total_length += abs(chr_end - chr_start + 1)
            te_total_infos[te_name] = (cur_copy_num, cur_total_length)

    return te_fl_infos, te_total_infos

def draw_pie(sizes, labels, output_pdf, title, is_percent):
    fig = plt.figure(figsize=(10, 5))
    # 数据
    # colors = ['#D53E4F', '#FC8D59', '#3288BD', '#FEE08B', '#FFFFBF', '#E6F598', '#ABDDA4', '#66C2A5', '#3288BD',
    #           '#4A90E2', '#E94E77', '#F5A623', '#7ED321', '#BD10E0', '#003366', '#D0D0D0', '#00274D', '#F7E04C',
    #           '#6B8E23', '#0033A0', '#E50000', '#228B22', '#C2B280', '#5B3A29', '#1E90FF', '#FF4500']
    # colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue']
    colors = [
        '#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3',
        '#ff69b4', '#00fa9a', '#add8e6', '#ff6347', '#dcdcdc', '#ff4500', '#ffb6c1', '#87ceeb',
        '#b0e57c', '#ff8c00', '#f08080', '#ffb000', '#90ee90', '#20b2aa', '#b0c4de', '#d8bfd8',
        '#f5deb3', '#ff1493', '#40e0d0', '#ff6347', '#ffb6c1', '#ff4500'
    ]

    if is_percent:
        # 绘制饼图
        plt.pie(sizes, colors=colors, autopct='%1.1f%%', shadow=False, startangle=140)
    else:
        plt.pie(sizes, colors=colors, shadow=False, startangle=140)
    # 使饼图为圆形
    plt.axis('equal')
    plt.legend(labels, loc="best")
    plt.title(title, pad=20, fontweight='bold')
    plt.savefig(output_pdf)
    # 显示图表
    # plt.show()

def convertGeneAnnotation2GTF(genome_paths, script_dir, output_dir, log):
    log.logger.info("Start convert Gene annotation Files")
    gtf_output_dir = output_dir + '/gene_gtf_files'
    input_files = ''
    new_genome_paths = []
    for genome_name, reference, gene_gtf, RNA_seq_dict in genome_paths:
        if gene_gtf is not None:
            input_files += gene_gtf + ' '
            base_name = os.path.basename(gene_gtf)
            new_gene_gtf = os.path.join(gtf_output_dir, os.path.splitext(base_name)[0] + ".gtf")
            new_genome_paths.append((genome_name, reference, new_gene_gtf, RNA_seq_dict))
        else:
            new_genome_paths.append((genome_name, reference, gene_gtf, RNA_seq_dict))
    if input_files != '':
        command = 'python ' + script_dir + '/makeCleanGeneGTF.py --input_files ' + input_files + ' --output_dir ' + gtf_output_dir
        log.logger.debug(command)
        os.system(command)
    return new_genome_paths

def run_HybridLTR(reference, tmp_output_dir, HybridLTR_home, threads, miu, recover, debug, is_output_lib, log):
    starttime = time.time()
    log.logger.debug('start HybridLTR detection...')
    HybridLTR_command = 'cd ' + HybridLTR_home + ' && python main.py --genome ' + reference \
                            + ' --out_dir ' + tmp_output_dir + ' --thread ' + str(threads) + ' --miu ' + str(miu) \
                        + ' --recover ' + str(recover) + ' --debug ' + str(debug) + ' --is_output_lib ' + str(is_output_lib)
    log.logger.debug(HybridLTR_command)
    os.system(HybridLTR_command)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("HybridLTR running time: %.8s s" % (dtime))

def assign_label_to_lib(lib, intact_LTR_labels):
    ltr_names, ltr_contigs = read_fasta(lib)
    no_label_ltr_names = []
    no_label_ltr_contigs = {}
    for name in ltr_names:
        seq = ltr_contigs[name]
        name = name.split('#')[0].replace('-lLTR', '-LTR').replace('-rLTR', '-LTR')
        no_label_ltr_names.append(name)
        no_label_ltr_contigs[name] = seq

    # Assign identical IDs to the same LTR and INT.
    ltr_index = 0
    stored_names = set()
    confident_ltr_cut_contigs = {}
    for name in no_label_ltr_names:
        if name in stored_names:
            continue
        seq = no_label_ltr_contigs[name]
        intact_ltr_name = name[:-4]
        label = intact_LTR_labels[intact_ltr_name]
        ltr_type = name[-4:]
        new_name = 'LTR_' + str(ltr_index) + ltr_type + '#' + label
        confident_ltr_cut_contigs[new_name] = seq

        # find the other type
        if ltr_type == '_LTR':
            other_ltr_type = '_INT'
        else:
            other_ltr_type = '_LTR'
        other_name = intact_ltr_name + other_ltr_type
        if no_label_ltr_contigs.__contains__(other_name):
            new_name = 'LTR_' + str(ltr_index) + other_ltr_type + '#' + label
            confident_ltr_cut_contigs[new_name] = no_label_ltr_contigs[other_name]
            stored_names.add(other_name)

        ltr_index += 1
    store_fasta(confident_ltr_cut_contigs, lib)

def remove_no_tirs(confident_tir_path, plant, TRsearch_dir, temp_dir):
    tir_names, tir_contigs = read_fasta(confident_tir_path)
    # keep sequence with short tir
    short_itr_contigs = get_short_tir_contigs(tir_contigs, plant)

    # The remaining sequence is handed over to itrsearch for TIR structure searching.
    temp_tir_path = confident_tir_path + '.no_short_tir.fa'
    for name in short_itr_contigs.keys():
        del tir_contigs[name]

    store_fasta(tir_contigs, temp_tir_path)
    all_copies_out, all_copies_log = run_itrsearch(TRsearch_dir, temp_tir_path, temp_dir)
    all_copies_out_name, all_copies_out_contigs = read_fasta(all_copies_out)
    all_copies_out_contigs.update(short_itr_contigs)
    temp_tir_path = confident_tir_path + '.all_tir.fa'
    store_fasta(all_copies_out_contigs, temp_tir_path)

    no_tir_contigs = {}
    no_tir_path = confident_tir_path + '.no_tir.fa'
    for name in tir_names:
        if name not in all_copies_out_contigs:
            no_tir_contigs[name] = tir_contigs[name]
    store_fasta(no_tir_contigs, no_tir_path)
    return temp_tir_path, no_tir_path

def map_fragment(start, end, chunk_size):
    start_chunk = start // chunk_size
    end_chunk = end // chunk_size

    if start_chunk == end_chunk:
        return start_chunk
    elif abs(end_chunk * chunk_size - start) < abs(end - end_chunk * chunk_size):
        return end_chunk
    else:
        return start_chunk


def multi_process_align_v2(query_path, subject_path, blastnResults_path, tmp_blast_dir, threads, chrom_length, coverage_threshold, category, is_full_length, is_removed_dir=True):
    tools_dir = ''
    if is_removed_dir:
        os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    if os.path.exists(blastnResults_path):
        os.remove(blastnResults_path)

    orig_names, orig_contigs = read_fasta(query_path)

    # blast_db_command = 'makeblastdb -dbtype nucl -in ' + subject_path + ' > /dev/null 2>&1'
    # os.system(blast_db_command)

    ref_names, ref_contigs = read_fasta(subject_path)
    # Sequence alignment consumes a significant amount of memory and disk space. Therefore, we also split the target sequences into individual sequences to reduce the memory required for each alignment, avoiding out of memory errors.
    # It is important to calculate the total number of bases in the sequences, and it must meet a sufficient threshold to increase CPU utilization.
    base_threshold = 10000000  # 10Mb
    target_files = []
    file_index = 0
    base_count = 0
    cur_contigs = {}
    for name in ref_names:
        cur_seq = ref_contigs[name]
        cur_contigs[name] = cur_seq
        base_count += len(cur_seq)
        if base_count >= base_threshold:
            cur_target = tmp_blast_dir + '/' + str(file_index) + '_target.fa'
            store_fasta(cur_contigs, cur_target)
            target_files.append(cur_target)
            makedb_command = 'makeblastdb -dbtype nucl -in ' + cur_target + ' > /dev/null 2>&1'
            os.system(makedb_command)
            cur_contigs = {}
            file_index += 1
            base_count = 0
    if len(cur_contigs) > 0:
        cur_target = tmp_blast_dir + '/' + str(file_index) + '_target.fa'
        store_fasta(cur_contigs, cur_target)
        target_files.append(cur_target)
        makedb_command = 'makeblastdb -dbtype nucl -in ' + cur_target + ' > /dev/null 2>&1'
        os.system(makedb_command)


    longest_repeat_files = []
    # 为了保证处理大型library时，blastn比对结果不会过大，我们保证每个簇里的序列数量为固定值
    avg_cluster_size = 50
    cluster_num = int(len(orig_names) / avg_cluster_size) + 1
    segments_cluster = divided_array(list(orig_contigs.items()), cluster_num)
    for partition_index, cur_segments in enumerate(segments_cluster):
        if len(cur_segments) <= 0:
            continue
        single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
        #print('current partition_index: ' + str(partition_index))
        if not os.path.exists(single_tmp_dir):
            os.makedirs(single_tmp_dir)
        split_repeat_file = single_tmp_dir + '/repeats_split.fa'
        cur_contigs = {}
        for item in cur_segments:
            cur_contigs[item[0]] = item[1]
        store_fasta(cur_contigs, split_repeat_file)
        repeats_path = (split_repeat_file, target_files, single_tmp_dir + '/temp.out',
                        single_tmp_dir + '/full_length.out', single_tmp_dir + '/tmp',
                        subject_path)
        longest_repeat_files.append(repeats_path)

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeat_files:
        job = ex.submit(multiple_alignment_blast_v2, file, tools_dir, coverage_threshold, category, chrom_length, is_full_length)
        jobs.append(job)
    ex.shutdown(wait=True)

    # 合并所有进程的结果，总体去除冗余
    chr_segments_list = []
    for job in as_completed(jobs):
        cur_chr_segments = job.result()
        chr_segments_list.append(cur_chr_segments)

    # 由于可能会有多个序列比对到同一个位置，因此我们对于基因组上的某一个位置，我们只取一条比对
    segment_len = 100000  # 100K
    # chr_segments -> {chr1: {seg0: [(start, end, status)], seg1: []}}
    # Status: 0 indicates that the fragment is not marked as found, while 1 indicates that the fragment is marked as found.
    prev_chr_segments = {}
    total_chr_len = 0
    # Divide the chromosome evenly into N segments to store fragments in segments and reduce retrieval time.
    for chr_name in chrom_length.keys():
        chr_len = chrom_length[chr_name]
        total_chr_len += chr_len
        if not prev_chr_segments.__contains__(chr_name):
            prev_chr_segments[chr_name] = {}
        prev_chr_segment_list = prev_chr_segments[chr_name]
        num_segments = chr_len // segment_len
        if chr_len % segment_len != 0:
            num_segments += 1
        for i in range(num_segments):
            prev_chr_segment_list[i] = []

    for cur_chr_segments in chr_segments_list:
        # Map the fragments to the corresponding segment,
        # and check if there is an overlap of over 95% with the fragment in the segment.
        for chr_name in cur_chr_segments.keys():
            cur_chr_segment_dict = cur_chr_segments[chr_name]
            prev_chr_segment_list = prev_chr_segments[chr_name]
            for seg_index in cur_chr_segment_dict.keys():
                cur_segment_frags = cur_chr_segment_dict[seg_index]
                for cur_frag in cur_segment_frags:
                    start = cur_frag[0]
                    end = cur_frag[1]
                    seq_name = cur_frag[2]
                    coverage = cur_frag[3]
                    seg_index = map_fragment(start, end, segment_len)

                    prev_segment_frags = prev_chr_segment_list[seg_index]
                    # Check if there is an overlap of over 95% between the fragment in the segment and the test fragment.
                    is_found = False
                    for prev_frag in prev_segment_frags:
                        overlap_len = get_overlap_len(prev_frag, cur_frag)
                        if overlap_len / abs(prev_frag[1] - prev_frag[0]) >= coverage_threshold and overlap_len / abs(
                                end - start) >= coverage_threshold:
                            is_found = True
                            break
                    if not is_found:
                        prev_segment_frags.append([start, end, seq_name, coverage])

    with open(blastnResults_path, 'w') as f_save:
        for chr_name in prev_chr_segments.keys():
            cur_chr_segments = prev_chr_segments[chr_name]
            for seg_index in cur_chr_segments.keys():
                segment_frags = cur_chr_segments[seg_index]
                for frag in segment_frags:
                    new_line = frag[2] + '\t' + chr_name + '\t' + '-1' + '\t' + '-1' + '\t' + '-1' + '\t' + '-1' + '\t' + '-1' + '\t' + '-1' + '\t' + str(frag[0]) + '\t' + str(frag[1]) + '\t' + str(frag[3]) + '\t' + '-1' + '\n'
                    f_save.write(new_line)

    if is_removed_dir:
        os.system('rm -rf ' + tmp_blast_dir)


def multiple_alignment_blast_v2(repeats_path, tools_dir, coverage_threshold, category, chrom_length, is_full_length):
    split_repeats_path = repeats_path[0]
    target_files = repeats_path[1]
    blastn2Results_path = repeats_path[2]
    full_length_out = repeats_path[3]
    tmp_dir = repeats_path[4]
    genome_path = repeats_path[5]
    os.system('rm -f ' + blastn2Results_path)
    for target_file in target_files:
        align_command = 'blastn -db ' + target_file + ' -num_threads ' \
                        + str(1) + ' -query ' + split_repeats_path + ' -evalue 1e-20 -outfmt 6 >> ' + blastn2Results_path
        os.system(align_command)

    # invoke the function to retrieve the full-length copies.
    lines = generate_full_length_out_v2(blastn2Results_path, full_length_out, split_repeats_path, genome_path, tmp_dir, tools_dir,
                             coverage_threshold, category, is_full_length)

    # 去除冗余的影响
    lines = list(lines)
    sorted_lines = sorted(lines, key=lambda x: (x[1], x[2], x[3]))
    test_fragments = {}
    for line in sorted_lines:
        seq_name = line[0]
        chr_name = line[1]
        chr_start = line[2]
        chr_end = line[3]
        coverage = line[4]
        if chr_name not in test_fragments:
            test_fragments[chr_name] = []
        fragments = test_fragments[chr_name]
        fragments.append((chr_start, chr_end, seq_name, coverage))

    # 由于可能会有多个序列比对到同一个位置，因此我们对于基因组上的某一个位置，我们只取一条比对
    segment_len = 100000  # 100K
    # chr_segments -> {chr1: {seg0: [(start, end, status)], seg1: []}}
    # Status: 0 indicates that the fragment is not marked as found, while 1 indicates that the fragment is marked as found.
    chr_segments = {}
    total_chr_len = 0
    # Divide the chromosome evenly into N segments to store fragments in segments and reduce retrieval time.
    for chr_name in chrom_length.keys():
        chr_len = chrom_length[chr_name]
        total_chr_len += chr_len
        if not chr_segments.__contains__(chr_name):
            chr_segments[chr_name] = {}
        cur_chr_segments = chr_segments[chr_name]
        num_segments = chr_len // segment_len
        if chr_len % segment_len != 0:
            num_segments += 1
        for i in range(num_segments):
            cur_chr_segments[i] = []

    # Map the fragments to the corresponding segment,
    # and check if there is an overlap of over 95% with the fragment in the segment.
    for chr_name in test_fragments.keys():
        fragments = test_fragments[chr_name]
        cur_chr_segments = chr_segments[chr_name]
        for cur_frag in fragments:
            start = cur_frag[0]
            end = cur_frag[1]
            seq_name = cur_frag[2]

            # if seq_name == 'chr_11_15708136-15719185-lLTR' and chr_name == 'Chr14' and start == 21886838 and end == 21890006:
            #     print('h')

            coverage = cur_frag[3]
            seg_index = map_fragment(start, end, segment_len)
            segment_frags = cur_chr_segments[seg_index]
            # Check if there is an overlap of over 95% between the fragment in the segment and the test fragment.
            is_found = False
            for prev_frag in segment_frags:
                overlap_len = get_overlap_len(prev_frag, cur_frag)
                if overlap_len / abs(prev_frag[1] - prev_frag[0]) >= coverage_threshold and overlap_len / abs(
                        end - start) >= coverage_threshold:
                    is_found = True
                    break
            if not is_found:
                segment_frags.append([start, end, seq_name, coverage])

    return chr_segments


def generate_full_length_out_v2(BlastnOut, full_length_out, TE_lib, reference, tmp_output_dir, tools_dir, full_length_threshold, category, is_full_length):
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)
    filter_tmp_out = filter_out_by_category(BlastnOut, tmp_output_dir, category)

    threads = 1
    divergence_threshold = 20
    search_struct = False
    full_length_annotations, copies_direct, all_query_copies = get_full_length_copies_from_blastn_v2(TE_lib, reference, filter_tmp_out,
                                                                             tmp_output_dir, threads,
                                                                             divergence_threshold,
                                                                             full_length_threshold,
                                                                             search_struct, tools_dir)
    lines = set()
    if not is_full_length:
        for query_name in all_query_copies.keys():
            query_copies = all_query_copies[query_name]
            for subject_pos in query_copies.keys():
                chr_name, chr_start, chr_end, coverage = query_copies[subject_pos]
                new_line = (query_name, chr_name, chr_start, chr_end, coverage)
                lines.add(new_line)
    else:
        for query_name in full_length_annotations.keys():
            query_name = str(query_name)
            for copy_annotation in full_length_annotations[query_name]:
                chr_pos = copy_annotation[0]
                annotation = copy_annotation[1]
                parts = chr_pos.split(':')
                chr_name = parts[0]
                chr_pos_parts = parts[1].split('-')
                chr_start = int(chr_pos_parts[0]) + 1
                chr_end = int(chr_pos_parts[1])
                new_line = (query_name, chr_name, chr_start, chr_end, -1)
                lines.add(new_line)

    return lines


def get_full_length_copies_from_blastn_v2(TE_lib, reference, blastn_out, tmp_output_dir, threads, divergence_threshold,
                                    full_length_threshold, search_struct, tools_dir):
    ref_names, ref_contigs = read_fasta(reference)

    query_names, query_contigs = read_fasta(TE_lib)
    new_query_contigs = {}
    for name in query_names:
        new_query_contigs[name.split('#')[0]] = query_contigs[name]
    query_contigs = new_query_contigs

    query_records = {}
    with open(blastn_out, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            info_parts = line.split('\t')
            query_name = info_parts[0].split('#')[0]
            subject_name = info_parts[1]
            q_start = int(info_parts[6])
            q_end = int(info_parts[7])
            s_start = int(info_parts[8])
            s_end = int(info_parts[9])
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))

    all_query_copies = {}
    full_length_copies = {}
    flank_full_length_copies = {}
    copies_direct = {}
    for idx, query_name in enumerate(query_records.keys()):
        subject_dict = query_records[query_name]
        if query_name not in query_contigs:
            continue
        query_len = len(query_contigs[query_name])
        skip_gap = query_len * (1 - 0.95)

        if str(query_name).__contains__('Helitron'):
            flanking_len = 5
        else:
            flanking_len = 50

        # if there are more than one longest query overlap with the final longest query over 90%,
        # then it probably the true TE
        longest_queries = []
        for subject_name in subject_dict.keys():
            subject_pos = subject_dict[subject_name]

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
                        cur_subject_start = frag[2]
                        cur_query_end = frag[1]
                        prev_subject_end = exist_frag[3]
                        prev_query_end = exist_frag[1]
                        if (cur_subject_start - prev_subject_end < skip_gap and cur_query_end > prev_query_end):
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
                        cur_subject_start = frag[2]
                        cur_query_end = frag[1]
                        prev_subject_end = exist_frag[3]
                        prev_query_end = exist_frag[1]
                        if (prev_subject_end - cur_subject_start < skip_gap and cur_query_end > prev_query_end):
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

                # print('subject pos size: %d' %(len(cur_cluster)))
                # record visited fragments
                visited_frag = {}
                for i in range(len(cur_cluster)):
                    # keep a longest query start from each fragment
                    prev_frag = cur_cluster[i]
                    if visited_frag.__contains__(prev_frag):
                        continue
                    prev_query_start = prev_frag[0]
                    prev_query_end = prev_frag[1]
                    prev_subject_start = prev_frag[2]
                    prev_subject_end = prev_frag[3]
                    prev_query_seq = (min(prev_query_start, prev_query_end), max(prev_query_start, prev_query_end))
                    prev_subject_seq = (
                        min(prev_subject_start, prev_subject_end), max(prev_subject_start, prev_subject_end))
                    prev_query_len = abs(prev_query_end - prev_query_start)
                    prev_subject_len = abs(prev_subject_end - prev_subject_start)
                    cur_longest_query_len = prev_query_len

                    cur_extend_num = 0
                    visited_frag[prev_frag] = 1
                    # try to extend query
                    for j in range(i + 1, len(cur_cluster)):
                        cur_frag = cur_cluster[j]
                        if visited_frag.__contains__(cur_frag):
                            continue
                        cur_query_start = cur_frag[0]
                        cur_query_end = cur_frag[1]
                        cur_subject_start = cur_frag[2]
                        cur_subject_end = cur_frag[3]
                        cur_query_seq = (min(cur_query_start, cur_query_end), max(cur_query_start, cur_query_end))
                        cur_subject_seq = (min(cur_subject_start, cur_subject_end), max(cur_subject_start, cur_subject_end))

                        # could extend
                        # extend right
                        if cur_query_end > prev_query_end:
                            # judge subject direction
                            if prev_subject_start < prev_subject_end and cur_subject_start < cur_subject_end:
                                # +
                                if cur_subject_end > prev_subject_end:
                                    # forward extend
                                    if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                            and cur_subject_start - prev_subject_end < skip_gap:  # \
                                        # and not is_same_query and not is_same_subject:
                                        # update the longest path
                                        prev_query_start = prev_query_start
                                        prev_query_end = cur_query_end
                                        prev_subject_start = prev_subject_start if prev_subject_start < cur_subject_start else cur_subject_start
                                        prev_subject_end = cur_subject_end
                                        cur_longest_query_len = prev_query_end - prev_query_start
                                        cur_extend_num += 1
                                        visited_frag[cur_frag] = 1
                                    elif cur_query_start - prev_query_end >= skip_gap:
                                        break
                            elif prev_subject_start > prev_subject_end and cur_subject_start > cur_subject_end:
                                # reverse
                                if cur_subject_end < prev_subject_end:
                                    # reverse extend
                                    if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                            and prev_subject_end - cur_subject_start < skip_gap:  # \
                                        # and not is_same_query and not is_same_subject:
                                        # update the longest path
                                        prev_query_start = prev_query_start
                                        prev_query_end = cur_query_end
                                        prev_subject_start = prev_subject_start if prev_subject_start > cur_subject_start else cur_subject_start
                                        prev_subject_end = cur_subject_end
                                        cur_longest_query_len = prev_query_end - prev_query_start
                                        cur_extend_num += 1
                                        visited_frag[cur_frag] = 1
                                    elif cur_query_start - prev_query_end >= skip_gap:
                                        break
                    # keep this longest query
                    if cur_longest_query_len != -1:
                        longest_queries.append(
                            (prev_query_start, prev_query_end, cur_longest_query_len, prev_subject_start,
                             prev_subject_end, abs(prev_subject_end - prev_subject_start), subject_name,
                             cur_extend_num))

        # To determine whether each copy has a coverage exceeding the full_length_threshold with respect
        # to the consensus sequence, retaining full-length copies.
        full_length_query_copies = {}
        full_length_flank_query_copies = {}
        query_copies = {}
        orig_query_len = len(query_contigs[query_name])
        for repeat in longest_queries:
            # Subject
            subject_name = repeat[6]
            subject_chr_start = 0

            if repeat[3] > repeat[4]:
                direct = '-'
                old_subject_start_pos = repeat[4] - 1
                old_subject_end_pos = repeat[3]
            else:
                direct = '+'
                old_subject_start_pos = repeat[3] - 1
                old_subject_end_pos = repeat[4]
            subject_start_pos = subject_chr_start + old_subject_start_pos
            subject_end_pos = subject_chr_start + old_subject_end_pos

            subject_pos = subject_name + ':' + str(subject_start_pos) + '-' + str(subject_end_pos)
            subject_seq = ref_contigs[subject_name][subject_start_pos: subject_end_pos]

            flank_subject_seq = ref_contigs[subject_name][
                                subject_start_pos - flanking_len: subject_end_pos + flanking_len]
            copies_direct[subject_pos] = direct
            cur_query_len = repeat[2]
            coverage = float(cur_query_len) / orig_query_len
            if coverage >= full_length_threshold:
                full_length_query_copies[subject_pos] = subject_seq
                full_length_flank_query_copies[subject_pos] = flank_subject_seq
            query_copies[subject_pos] = (subject_name, subject_start_pos, subject_end_pos, coverage)
        full_length_copies[query_name] = full_length_query_copies
        flank_full_length_copies[query_name] = full_length_flank_query_copies
        all_query_copies[query_name] = query_copies

    # The candidate full-length copies and the consensus are then clustered using cd-hit-est,
    # retaining copies that belong to the same cluster as the consensus.
    split_files = []
    cluster_dir = tmp_output_dir + '/cluster'
    os.system('rm -rf ' + cluster_dir)
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)

    for query_name in full_length_copies.keys():
        query_copies = full_length_copies[query_name]
        flank_query_copies = flank_full_length_copies[query_name]
        fc_path = cluster_dir + '/' + query_name + '.fa'
        store_fasta(query_copies, fc_path)
        split_files.append((fc_path, query_name, query_copies, flank_query_copies))

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for ref_index, cur_file in enumerate(split_files):
        input_file = cur_file[0]
        query_name = cur_file[1]
        query_copies = cur_file[2]
        flank_query_copies = cur_file[3]
        job = ex.submit(get_structure_info, input_file, query_name, query_copies,
                        flank_query_copies, cluster_dir, search_struct, tools_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    full_length_annotations = {}
    for job in as_completed(jobs):
        annotations = job.result()
        full_length_annotations.update(annotations)
    return full_length_annotations, copies_direct, all_query_copies

def generate_panTE_PAV(new_te_contigs, pan_te_fl_infos, output_dir, log):
    pav_table = output_dir + '/panHiTE_PAV.tsv'
    lines = []
    first_line = 'TE_families\t'
    genome_names = pan_te_fl_infos.keys()
    for i, genome_name in enumerate(genome_names):
        first_line += genome_name
        if i != len(genome_names) - 1:
            first_line += '\t'
        else:
            first_line += '\n'
    lines.append(first_line)

    te_fl_occur_genomes = {}
    for te_name in new_te_contigs.keys():
        cur_line = te_name + '\t'
        for i, genome_name in enumerate(genome_names):
            te_fl_infos = pan_te_fl_infos[genome_name]
            if te_name in te_fl_infos:
                cur_copy_num, cur_fl_length = te_fl_infos[te_name]
            else:
                cur_copy_num = 0
            cur_line += str(cur_copy_num)
            if i != len(genome_names) - 1:
                cur_line += '\t'
            else:
                cur_line += '\n'
        lines.append(cur_line)

    with open(pav_table, 'w') as f_save:
        for line in lines:
            f_save.write(line)

    # 调用 drawCorePanPAV.R 生成TE饱和曲线图
    script_path = cur_dir + '/RNA_seq/drawCorePanPAV.R'
    cmd = 'cd ' + output_dir + ' && Rscript ' + script_path + ' ' + pav_table + ' 500 panHiTE'
    log.logger.debug(cmd)
    os.system(cmd)

    return te_fl_occur_genomes



