import codecs
import json
import multiprocessing
import os

import logging
import re
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from logging import handlers

import subprocess

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

def run_HelitronScanner(sh_dir, temp_dir, cur_segments, HSDIR, HSJAR, partition_index):
    candidate_Helitrons = {}
    cur_candidate_Helitrons_path = temp_dir + '/' + str(partition_index) + '.fa'
    cur_candidate_Helitrons = {}
    for item in cur_segments:
        query_name = item[0]
        copies = item[1]
        if len(copies) > 0:
            copy = copies[0]
            orig_seq = str(copy[4])
            # 如果有连续10个以上的N或者搜索的TSD有连续>=2个N，就过滤掉
            if orig_seq.__contains__('NNNNNNNNNN'):
                continue
            cur_candidate_Helitrons[query_name] = orig_seq
    store_fasta(cur_candidate_Helitrons, cur_candidate_Helitrons_path)
    HelitronScanner_command = 'cd ' + temp_dir + ' && ' + 'sh ' + sh_dir + '/run_helitron_scanner.sh ' \
                              + str(partition_index) + ' ' + cur_candidate_Helitrons_path + ' ' + HSDIR + ' ' + HSJAR
    os.system(HelitronScanner_command + '> /dev/null 2>&1')

    cur_helitron_out = temp_dir + '/' + str(partition_index) + '.HelitronScanner.draw.hel.fa'
    cur_rc_helitron_out = temp_dir + '/' + str(partition_index) + '.HelitronScanner.draw.rc.hel.fa'
    cur_names, cur_contigs = read_fasta(cur_helitron_out)
    cur_rc_names, cur_rc_contigs = read_fasta(cur_rc_helitron_out)
    candidate_Helitrons.update(cur_contigs)
    candidate_Helitrons.update(cur_rc_contigs)
    return candidate_Helitrons



def run_EAHelitron(flanking_len, temp_dir, cur_segments, EAHelitron, partition_index):
    all_candidate_helitron_contigs = {}
    all_candidate_helitron_path = temp_dir + '/' + str(partition_index) + '.fa'
    for item in cur_segments:
        query_name = item[0]
        seq = item[1]

        raw_start = flanking_len + 1
        raw_end = len(seq) - flanking_len

        # 如果有连续10个以上的N或者搜索的TSD有连续>=2个N，就过滤掉
        if seq.__contains__('NNNNNNNNNN'):
            continue
        new_query_name = query_name + '-rawstart_' + str(raw_start) + '-rawend_' + str(raw_end)
        all_candidate_helitron_contigs[new_query_name] = seq

    store_fasta(all_candidate_helitron_contigs, all_candidate_helitron_path)
    EAHelitron_command = 'cd ' + temp_dir + ' && ' + 'perl ' + EAHelitron + '/EAHelitron -o ' + str(partition_index) + ' -u 20000 -T "ATC" -r 3 ' + all_candidate_helitron_path
    os.system(EAHelitron_command + '> /dev/null 2>&1')

    all_EAHelitron_res = temp_dir + '/' + str(partition_index) + '.5.fa'
    all_copies_out_names, all_copies_out_contigs = read_fasta_v1(all_EAHelitron_res)
    # 对all_copies_out_contigs按照name进行分组
    # group_copies_contigs -> {query_name: {name: seq}}
    group_copies_contigs = {}
    for cur_name in all_copies_out_contigs.keys():
        raw_name = cur_name.split(' ')[1]
        parts = raw_name.split(':')
        query_name = parts[0]
        if not group_copies_contigs.__contains__(query_name):
            group_copies_contigs[query_name] = {}
        cur_copies_out_contigs = group_copies_contigs[query_name]
        cur_copies_out_contigs[cur_name] = all_copies_out_contigs[cur_name]

    candidate_Helitrons = {}
    for query_name in group_copies_contigs.keys():
        cur_copies_out_contigs = group_copies_contigs[query_name]
        #记录同一个query_name下不同的copy_index对应的值
        copies_candidate = {}
        # 2. 合并，选择distance最小（相同则取最长）的那条序列,当做这条拷贝的代表序列
        # 根据name中的copy_index进行分组，copy_index = name.split('-C_')[1].split('-')[0]。
        # copies_candidate -> {copy_index: (min_distance_seq_name, min_distance, seq_len, first_6bp)}
        for name in cur_copies_out_contigs.keys():
            raw_name = name.split(' ')[1]
            parts = raw_name.split(':')
            raw_name = parts[0]
            pos_parts = parts[1].split('..')
            cur_start = int(pos_parts[0])
            cur_end = int(pos_parts[1])
            if cur_start > cur_end:
                tmp = cur_start
                cur_start = cur_end
                cur_end = tmp
            cur_seq = cur_copies_out_contigs[name]
            # 取原始边界
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
    seq_names, seq_contigs = read_fasta(longest_repeats_flanked_path)
    # ---------------------------------------Terminal Repeat Search--------------------------------------------------
    segments_cluster = divided_array(list(seq_contigs.items()), threads)

    os.system('rm -rf '+temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    for partition_index, cur_segments in enumerate(segments_cluster):
        obj = ex.submit(run_EAHelitron, flanking_len, temp_dir, cur_segments, EAHelitron, partition_index)
        objs.append(obj)
        job_id += 1
    ex.shutdown(wait=True)
    candidate_Helitrons = {}
    for obj in as_completed(objs):
        cur_candidate_Helitrons = obj.result()
        candidate_Helitrons.update(cur_candidate_Helitrons)
    store_fasta(candidate_Helitrons, output)


def multi_process_helitronscanner(all_copies, output, sh_dir, temp_dir, HSDIR, HSJAR, threads):
    segments_cluster = divided_array(list(all_copies.items()), threads)

    os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    for partition_index, cur_segments in enumerate(segments_cluster):
        obj = ex.submit(run_HelitronScanner, sh_dir, temp_dir, cur_segments, HSDIR, HSJAR, partition_index)
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

def run_LTR_harvest(reference, tmp_output_dir, log):
    starttime = time.time()
    log.logger.debug('start LTR_harvest detection...')
    ltrharvest_command1 = 'gt suffixerator -db ' + reference + ' -indexname ' \
                          + reference + ' -tis -suf -lcp -des -ssp -sds -dna'
    ltrharvest_command2 = 'gt ltrharvest -index ' + reference \
                          + ' -seed 20 -minlenltr 100 -maxlenltr 7000 -similar 85 -motif TGCA -motifmis 1 -mintsd 4 -maxtsd 6 ' \
                            '-vic 10 -seqids yes > ' + tmp_output_dir + '/genome_all.fa.harvest.scn'

    os.system(ltrharvest_command1)
    os.system(ltrharvest_command2)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.debug("LTR_harvest running time: %.8s s" % (dtime))

def run_LTR_retriever(reference, tmp_output_dir, threads, log):
    starttime = time.time()
    log.logger.debug('start LTR_retriever detection...')
    LTR_retriever_command = 'cd ' + tmp_output_dir + ' && LTR_retriever -genome ' + reference \
                            + ' -inharvest ' + tmp_output_dir + '/genome_all.fa.rawLTR.scn -noanno -threads ' + str(threads)
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
            #判断最后一个protein的比对方向
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

                #如果序列已经反向互补了，那么原本的位置也应该互换
                if last_direct == '-':
                    query_start = len(query_seq) - query_start + 1
                    query_end = len(query_seq) - query_end + 1

                new_query_name += ':' + protein_name + '_' + str(query_start) + '_' + str(query_end)

            query_seq = endWithPolyA(query_seq, last_query_start, last_query_end)
            new_query_name += '-len_'+str(len(query_seq))
            candidate_LINE[new_query_name] = query_seq
            node_index += 1
    return candidate_LINE

def multiple_alignment_blastx(repeats_path, blast_program_dir, tools_dir):
    split_repeats_path = repeats_path[0]
    protein_db_path = repeats_path[1]
    temp_dir = repeats_path[2]
    blastx2Results_path = temp_dir + '/temp.out'
    align_command = blast_program_dir + '/bin/blastx -db ' + protein_db_path + ' -num_threads ' \
                    + str(1) + ' -query ' + split_repeats_path + ' -word_size 7 -outfmt 6 > ' + blastx2Results_path
    os.system(align_command)

    # 1.首先把longest_repeats.fa比对到LINERep.fa上识别完整的pol。
    # 2.判断方向，正向就去尾部寻找polyA，反向就头部寻找polyT。
    # 3.确定尾巴之后，取尾巴后的30-3bp一共28个kmer，然后再延伸头部（指的是最开始domain，如果只有pol就是pol位置，如果是gag就是gag位置）500bp，
    # 4.再将这28个kmer比对到延伸的序列，从长到短取最佳匹配的位置当做头部的TSD起始位置。
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

# def cut_repeat_v2(sam_path_bwa, repeats_path, cut_repeats_path):
#     query_records = {}
#     samfile = pysam.AlignmentFile(sam_path_bwa, "rb")
#     for read in samfile.fetch():
#         if read.is_unmapped:
#             continue
#         query_name = read.query_name
#         reference_name = read.reference_name
#         cigar = read.cigartuples
#         cigarstr = read.cigarstring
#         NM_tag = 0
#         try:
#             NM_tag = read.get_tag('NM')
#         except KeyError:
#             NM_tag = -1
#         identity = compute_identity(cigarstr, NM_tag, 'BLAST')
#         identity = float(identity) * 100
#         is_reverse = read.is_reverse
#         alignment_len = read.query_alignment_length
#         # pos start from 1, change to 0
#         q_start = int(read.query_alignment_start)  # [q_start, q_end)
#         q_end = int(read.query_alignment_end)
#         if q_start > q_end:
#             tmp = q_start
#             q_start = q_end
#             q_end = tmp
#         if not query_records.__contains__(query_name):
#             query_records[query_name] = []
#         records = query_records[query_name]
#         records.append((reference_name, alignment_len, identity, q_start, q_end))
#         query_records[query_name] = records
#
#     repeat_contignames, repeat_contigs = read_fasta(repeats_path)
#     cut_repeats = {}
#     for query_name in query_records.keys():
#         query_seq = repeat_contigs[query_name]
#         query_len = len(query_seq)
#         records = query_records[query_name]
#         for i, record in enumerate(records):
#             # filter first alignment
#             if i == 0:
#                 continue
#             identity = record[2]
#             q_start = record[3]
#             q_end = record[4]
#             if identity < 95:
#                 continue
#             # get repeats boundary by getting all alignment sequences
#             new_seq = query_seq[q_start: q_end]
#             new_query_name = query_name + '-p_' + str(i) + '-len_' + str(len(new_seq))
#             cut_repeats[new_query_name] = new_seq
#     store_fasta(cut_repeats, cut_repeats_path)

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

# remove one perfect match record, exclude a situation
# that sequence has only one perfect match position,
# and many clips position
# def generate_blastlike_output(sam_paths, blastnResults_path, not_multi_mapping_repeatIds):
#     not_multi_mapping_repeatIds_dict = {}
#     for repeat_id in not_multi_mapping_repeatIds:
#         not_multi_mapping_repeatIds_dict[repeat_id] = 1
#     query_names = {}
#     with open(blastnResults_path, 'w') as f_save:
#         for item in sam_paths:
#             sam_path = item[0]
#             repeats_path = item[1]
#             contignames, contigs = read_fasta(repeats_path)
#             samfile = pysam.AlignmentFile(sam_path, "rb")
#             for read in samfile.fetch():
#                 if read.is_unmapped:
#                     continue
#                 query_name = read.query_name
#                 # filter not multiple mapping repeat
#                 if not_multi_mapping_repeatIds_dict.__contains__(query_name):
#                     continue
#                 target_name = read.reference_name
#                 cigar = read.cigarstring
#                 NM_tag = 0
#                 try:
#                     NM_tag = read.get_tag('NM')
#                 except KeyError:
#                     NM_tag = -1
#                 identity = compute_identity(cigar, NM_tag, 'BLAST')
#                 identity = float(identity) * 100
#                 match_base = int(read.query_alignment_length)
#                 q_start = int(read.query_alignment_start)
#                 q_end = int(read.query_alignment_end)
#                 t_start = int(read.reference_start)
#                 t_end = int(read.reference_end)
#                 query_length = int(len(contigs[query_name]))
#
#                 if not query_names.__contains__(query_name) and (identity > 95 and float(match_base)/query_length > 0.95):
#                     query_names[query_name] = 1
#                     continue
#
#                 if read.is_reverse:
#                     temp = t_start
#                     t_start = t_end
#                     t_end = temp
#                     strand = '+' if t_end >= t_start else '-'
#
#                 f_save.write(str(query_name) + '\t' + str(target_name) + '\t' + str(identity) + '\t' + str(match_base)
#                              + '\t' + str('X') + '\t' + str('X') + '\t' + str(q_start) + '\t' + str(q_end)
#                              + '\t' + str(t_start) + '\t' + str(t_end) + '\n')
#
# def generate_blastlike_output_parallel(sam_paths, blastnResults_path, not_multi_mapping_repeatIds, partitions_num):
#     not_multi_mapping_repeatIds_dict = {}
#     for repeat_id in not_multi_mapping_repeatIds:
#         not_multi_mapping_repeatIds_dict[repeat_id] = 1
#
#     parse_records = []
#     for sam_path in sam_paths:
#         samfile = pysam.AlignmentFile(sam_path, "rb")
#         for read in samfile.fetch():
#             if read.is_unmapped:
#                 continue
#             query_name = read.query_name
#             target_name = read.reference_name
#             cigar = read.cigarstring
#             try:
#                 NM_tag = read.get_tag('NM')
#             except KeyError:
#                 NM_tag = -1
#             identity = compute_identity(cigar, NM_tag, 'BLAST')
#             identity = float(identity) * 100
#             match_base = int(read.query_alignment_length)
#             q_start = int(read.query_alignment_start)
#             q_end = int(read.query_alignment_end)
#             t_start = int(read.reference_start)
#             t_end = int(read.reference_end)
#
#             if read.is_reverse:
#                 temp = t_start
#                 t_start = t_end
#                 t_end = temp
#                 strand = '+' if t_end >= t_start else '-'
#             parse_records.append((query_name, target_name, identity, match_base, q_start, q_end, t_start, t_end))
#
#     records_cluster = split2cluster(parse_records, partitions_num)
#
#     blast_records = multiprocessing.Manager().dict()
#     pool = multiprocessing.Pool(processes=partitions_num)
#     for partiton_index in records_cluster.keys():
#         cur_records = records_cluster[partiton_index]
#         pool.apply_async(filter_not_multi_mapping, (cur_records, not_multi_mapping_repeatIds_dict, partiton_index, blast_records,))
#     pool.close()
#     pool.join()
#
#     with open(blastnResults_path, 'w') as f_save:
#         for partiton_index in blast_records.keys():
#             for record in blast_records[partiton_index]:
#                 f_save.write(str(record[0]) + '\t' + str(record[1]) + '\t' + str(record[2]) + '\t' + str(record[3])
#                              + '\t' + str('X') + '\t' + str('X') + '\t' + str(record[4]) + '\t' + str(record[5])
#                              + '\t' + str(record[6]) + '\t' + str(record[7]) + '\n')
#
# def cut_repeat_v1(sam_paths, HS_gap, ID_gap, repeats_file, raw_cut_file):
#     repeat_contignames, repeat_contigs = read_fasta(repeats_file)
#     all_fragments = {}
#
#     # if a repeat can only align to a position complete,
#     # and other records contain large H/S or I/D, it should be spliced
#     query_records = {}
#     for sam_path in sam_paths:
#         samfile = pysam.AlignmentFile(sam_path, "rb")
#         for read in samfile.fetch():
#             if read.is_unmapped:
#                 continue
#             query_name = read.query_name
#             reference_name = read.reference_name
#             cigar = read.cigartuples
#             cigarstr = read.cigarstring
#             is_reverse = read.is_reverse
#
#             if not query_records.__contains__(query_name):
#                 query_records[query_name] = []
#             records = query_records[query_name]
#             records.append((query_name, reference_name, cigar, cigarstr, is_reverse))
#             query_records[query_name] = records
#
#     # delete Chimerism to avoid fragments
#     # for query_name in query_records.keys():
#     #     is_chimerism = False
#     #     for record in query_records[query_name]:
#     #         query_seq = repeat_contigs[query_name]
#     #         query_len = len(query_seq)
#     #         float(record.query_alignment_length)/query_len > 90
#     #     if is_chimerism:
#     #         del query_records['query_name']
#
#     pattern = r'[0-9]\d+M'
#     repeats_tobe_spliced = {}
#     for query_name in query_records.keys():
#         complete_alignment_num = 0
#         for record in query_records[query_name]:
#             cigar = record[2]
#             cigarstr = str(record[3])
#             # complete Match in cigar
#             if re.match(pattern, cigarstr) is not None:
#                 complete_alignment_num += 1
#             else:
#                 # if cigar contains small I/D or small H/S, it can be seen as a complete record
#                 query_seq = repeat_contigs[query_name]
#                 query_len = len(query_seq)
#                 is_complete = True
#                 for c in cigar:
#                     if (c[0] == 4 and c[1] >= HS_gap * query_len) \
#                             or (c[0] == 5 and c[1] >= HS_gap * query_len) \
#                             or (c[0] == 1 and c[1] >= ID_gap * query_len) \
#                             or (c[0] == 2 and c[1] >= ID_gap * query_len):
#                         is_complete = False
#                         break
#                 if is_complete:
#                     complete_alignment_num += 1
#         if complete_alignment_num == 1:
#             repeats_tobe_spliced[query_name] = 1
#
#     #print(repeats_tobe_spliced)
#
#     for query_name in repeats_tobe_spliced.keys():
#         for record in query_records[query_name]:
#             reference_name = record[1]
#             cigar = record[2]
#             cigarstr = record[3]
#             is_reverse = record[4]
#             query_seq = repeat_contigs[query_name]
#             query_len = len(query_seq)
#             if is_reverse:
#                 query_seq = getReverseSequence(query_seq)
#             if not all_fragments.__contains__(query_name):
#                 all_fragments[query_name] = []
#             fragments = []
#
#             # parse cigar
#             repeat_index = 0
#             last_cigar = -1
#             frag = ''
#             frag_start_pos = 0
#             is_split = False
#             for c in cigar:
#                 if last_cigar != c[0]:
#                     last_cigar = c[0]
#                     # large gap, split repeat
#                     if ((c[0] == 4 and c[1] >= HS_gap * query_len)
#                                        or (c[0] == 5 and c[1] >= HS_gap * query_len)
#                                        or (c[0] == 1 and c[1] >= ID_gap * query_len)
#                                        or (c[0] == 2 and c[1] >= ID_gap * query_len)):
#                         is_split = True
#                         if frag != '':
#                             fragments.append((frag_start_pos, len(frag), frag))
#                             frag_start_pos += len(frag)
#                             if c[0] != 2:
#                                 frag_start_pos += c[1]
#                             frag = ''
#
#                 if (c[0] == 4 or c[0] == 5) and c[1] >= HS_gap * query_len:
#                     # if cigar is large H/S, repeat index increment
#                     repeat_index += c[1]
#                     continue
#                 elif c[0] == 2 or c[0] == 3:
#                     # if cigar is D/N, repeat index should stay
#                     continue
#                 elif c[0] == 1 and c[1] >= ID_gap * query_len:
#                     # if cigar is large I, repeat index increment
#                     repeat_index += c[1]
#                     continue
#                 else:
#                     # if cigar is M or small I/D or small H/S, store sequence and repeat index increment
#                     frag += query_seq[repeat_index:repeat_index+c[1]]
#                     repeat_index += c[1]
#             if frag != '':
#                 fragments.append((frag_start_pos, len(frag), frag))
#             old_fragments = all_fragments[query_name]
#             all_fragments[query_name] = old_fragments + fragments
#
#     # if keep original repeat
#     # de-duplicate reverse-complementarty repeat
#     all_unique_fragments = {}
#     for query_name in all_fragments.keys():
#         frag_seq_set = []
#         frag_set = []
#
#         for frag in all_fragments[query_name]:
#             if frag[2] not in frag_seq_set and getReverseSequence(frag[2]) not in frag_seq_set:
#                 frag_seq_set.append(frag[2])
#             frag_set.append(frag)
#         all_unique_fragments[query_name] = frag_set
#
#     node_index = 0
#     with open(raw_cut_file, 'w') as f_save:
#         for query_name in all_unique_fragments.keys():
#             for unique_frag in all_unique_fragments[query_name]:
#                 f_save.write('>Node_' + str(node_index) + '-len_' + str(len(unique_frag[2]))
#                              + '\n' + unique_frag[2] + '\n')
#                 node_index += 1
#
#         for query_name in query_records.keys():
#             if not repeats_tobe_spliced.__contains__(query_name):
#                 seq = repeat_contigs[query_name]
#                 f_save.write('>Node_' + str(node_index) + '-len_' + str(len(seq))
#                              + '\n' + seq + '\n')
#                 node_index += 1
#
#
# def cut_repeat(sam_paths, HS_gap, ID_gap, repeats_file):
#     repeat_contignames, repeat_contigs = read_fasta(repeats_file)
#     all_fragments = {}
#     #original_repeat_occurrences = {}
#     for sam_path in sam_paths:
#         samfile = pysam.AlignmentFile(sam_path, "rb")
#         for read in samfile.fetch():
#             if read.is_unmapped:
#                 continue
#             query_name = read.query_name
#             reference_name = read.reference_name
#             cigar = read.cigartuples
#             cigarstr = read.cigarstring
#             is_reverse = read.is_reverse
#             query_seq = repeat_contigs[query_name]
#             query_len = len(query_seq)
#             if is_reverse:
#                 query_seq = getReverseSequence(query_seq)
#             if not all_fragments.__contains__(query_name):
#                 all_fragments[query_name] = []
#             fragments = []
#
#             # parse cigar
#             repeat_index = 0
#             last_cigar = -1
#             frag = ''
#             frag_start_pos = 0
#             is_split = False
#             for c in cigar:
#                 if last_cigar != c[0]:
#                     last_cigar = c[0]
#                     # large gap, split repeat
#                     if ((c[0] == 4 and c[1] >= HS_gap)
#                                        or (c[0] == 5 and c[1] >= HS_gap)
#                                        or (c[0] == 1 and c[1] >= ID_gap * query_len)
#                                        or (c[0] == 2 and c[1] >= ID_gap * query_len)):
#                         is_split = True
#                         if frag != '':
#                             fragments.append((frag_start_pos, len(frag), frag))
#                             frag_start_pos += len(frag)
#                             if c[0] != 2:
#                                 frag_start_pos += c[1]
#                             frag = ''
#
#                 if (c[0] == 4 or c[0] == 5) and c[1] >= HS_gap:
#                     # if cigar is large H/S, repeat index increment
#                     repeat_index += c[1]
#                     continue
#                 elif c[0] == 2 or c[0] == 3:
#                     # if cigar is D/N, repeat index should stay
#                     continue
#                 elif c[0] == 1 and c[1] >= ID_gap * query_len:
#                     # if cigar is large I, repeat index increment
#                     repeat_index += c[1]
#                     continue
#                 else:
#                     # if cigar is M or small I/D or small H/S, store sequence and repeat index increment
#                     frag += query_seq[repeat_index:repeat_index+c[1]]
#                     repeat_index += c[1]
#             if frag != '':
#                 fragments.append((frag_start_pos, len(frag), frag))
#             old_fragments = all_fragments[query_name]
#             all_fragments[query_name] = old_fragments + fragments
#         samfile.close()
#
#     # if keep original repeat
#     # de-duplicate reverse-complementarty repeat
#     all_unique_fragments = {}
#     for query_name in all_fragments.keys():
#         frag_seq_set = []
#         frag_set = []
#
#         for frag in all_fragments[query_name]:
#             if frag[2] not in frag_seq_set and getReverseSequence(frag[2]) not in frag_seq_set:
#                 frag_seq_set.append(frag[2])
#             frag_set.append(frag)
#         all_unique_fragments[query_name] = frag_set
#
#     return all_unique_fragments
#
# def get_multiple_alignment_repeat(sam_paths):
#     not_multi_mapping_repeatIds = []
#     multi_mapping_repeatIds = []
#     mapping_repeatIds = {}
#     for item in sam_paths:
#         sam_path = item[0]
#         repeats_path = item[1]
#         samfile = pysam.AlignmentFile(sam_path, "rb")
#         for read in samfile.fetch():
#             query_name = read.query_name
#             if not mapping_repeatIds.__contains__(query_name):
#                 mapping_repeatIds[query_name] = 0
#
#             if read.is_unmapped:
#                 continue
#             else:
#                 count = mapping_repeatIds[query_name]
#                 count += 1
#                 mapping_repeatIds[query_name] = count
#
#     for query_name in mapping_repeatIds.keys():
#         count = mapping_repeatIds[query_name]
#         if count > 1:
#             multi_mapping_repeatIds.append(query_name)
#         else:
#             not_multi_mapping_repeatIds.append(query_name)
#
#     return multi_mapping_repeatIds, not_multi_mapping_repeatIds

# def judgeReduceThreads(unique_kmer_path, threads, log):
#     file_size = os.path.getsize(unique_kmer_path) / (1024 * 1024 * 1024)
#     mem = psutil.virtual_memory()
#     free_memory = float(mem.free) / (1024 * 1024 * 1024)
#
#     if free_memory / 4 < file_size * threads:
#         reduce_threads = int(free_memory / (4 * file_size))
#         log.logger.debug('----Warning:\nDetect the free memory of your machine is %f GB.\n'
#               'The kmer.txt file is %f GB, which will be used in each thread.\n'
#               'The number of thread you set is %d. According to our experience,\n'
#               'To avoid the risk of out of memory, free memory should be more than 4 times\n'
#               'higher than thread_num*kmer_size. Thus, reduce the number of thread to %d.\n' % (
#               free_memory, file_size, threads, reduce_threads))
#         return reduce_threads
#     else:
#         return threads

# def get_alignment_info(sam_paths):
#     unmapped_repeatIds = []
#     single_mapped_repeatIds = []
#     multi_mapping_repeatIds = []
#     mapping_repeatIds = {}
#
#     for item in sam_paths:
#         sam_path = item[0]
#         repeats_path = item[1]
#         samfile = pysam.AlignmentFile(sam_path, "rb")
#         for read in samfile.fetch():
#             query_name = read.query_name
#             if not mapping_repeatIds.__contains__(query_name):
#                 mapping_repeatIds[query_name] = 0
#
#             if read.is_unmapped:
#                 continue
#             else:
#                 count = mapping_repeatIds[query_name]
#                 count += 1
#                 mapping_repeatIds[query_name] = count
#
#     for query_name in mapping_repeatIds.keys():
#         count = mapping_repeatIds[query_name]
#         if count <= 0:
#             unmapped_repeatIds.append(query_name)
#         elif count == 1:
#             single_mapped_repeatIds.append(query_name)
#         else:
#             multi_mapping_repeatIds.append(query_name)
#
#     return unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds
#
#
# def get_alignment_info_v4(sam_paths, repeats_file):
#     repeat_contignames, repeat_contigs = read_fasta(repeats_file)
#
#     unmapped_repeatIds = []
#     single_mapped_repeatIds = []
#     multi_mapping_repeatIds = []
#     segmental_duplication_repeatIds = []
#
#     query_records = {}
#     for sam_path in sam_paths:
#         samfile = pysam.AlignmentFile(sam_path, "rb")
#         for read in samfile.fetch():
#             if read.is_unmapped:
#                 continue
#             query_name = read.query_name
#             reference_name = read.reference_name
#             cigar = read.cigartuples
#             cigarstr = read.cigarstring
#             NM_tag = 0
#             try:
#                 NM_tag = read.get_tag('NM')
#             except KeyError:
#                 NM_tag = -1
#             identity = compute_identity(cigarstr, NM_tag, 'BLAST')
#             identity = float(identity) * 100
#             is_reverse = read.is_reverse
#             alignment_len = read.query_alignment_length
#
#             if not query_records.__contains__(query_name):
#                 query_records[query_name] = []
#             records = query_records[query_name]
#             records.append((query_name, reference_name, cigar, cigarstr, is_reverse, alignment_len, identity))
#             query_records[query_name] = records
#
#     for query_name in query_records.keys():
#         records = query_records[query_name]
#         freq = len(records)
#         if freq == 1:
#             single_mapped_repeatIds.append(query_name)
#         elif freq <= 4:
#             segmental_duplication_repeatIds.append(query_name)
#         else:
#             multi_mapping_repeatIds.append((query_name, freq))
#     return single_mapped_repeatIds, multi_mapping_repeatIds, segmental_duplication_repeatIds
#
# def get_alignment_info_v1(sam_paths, repeats_file):
#     repeat_contignames, repeat_contigs = read_fasta(repeats_file)
#
#     unmapped_repeatIds = []
#     single_mapped_repeatIds = []
#     multi_mapping_repeatIds = []
#     segmental_duplication_repeatIds = []
#
#     query_records = {}
#     for sam_path in sam_paths:
#         samfile = pysam.AlignmentFile(sam_path, "rb")
#         for read in samfile.fetch():
#             if read.is_unmapped:
#                 continue
#             query_name = read.query_name
#             reference_name = read.reference_name
#             cigar = read.cigartuples
#             cigarstr = read.cigarstring
#             NM_tag = 0
#             try:
#                 NM_tag = read.get_tag('NM')
#             except KeyError:
#                 NM_tag = -1
#             identity = compute_identity(cigarstr, NM_tag, 'BLAST')
#             identity = float(identity) * 100
#             is_reverse = read.is_reverse
#             alignment_len = read.query_alignment_length
#
#             if not query_records.__contains__(query_name):
#                 query_records[query_name] = []
#             records = query_records[query_name]
#             records.append((query_name, reference_name, cigar, cigarstr, is_reverse, alignment_len, identity))
#             query_records[query_name] = records
#
#     for query_name in query_records.keys():
#         complete_alignment_num = 0
#         high_identity_num = 0
#         # other cigars all the same (except first one) are regarded as segmental duplication(LCR)
#         # is_special_lcr = True
#         # last_cigarstr = ''
#         # first_cigarstr = ''
#         for i, record in enumerate(query_records[query_name]):
#             if i == 0:
#                 continue
#             cigar = record[2]
#             cigarstr = str(record[3])
#             alignment_len = record[5]
#             identity = record[6]
#             query_seq = repeat_contigs[query_name]
#             query_len = len(query_seq)
#             if float(alignment_len) / query_len >= 0.9 and identity >= 80:
#                 complete_alignment_num += 1
#                 if identity >= 90:
#                     high_identity_num += 1
#
#         if complete_alignment_num == 0:
#             single_mapped_repeatIds.append(query_name)
#         elif complete_alignment_num > 0:
#             # low copy number and all of them are high identicial, they are LCR
#             # if complete_alignment_num < 4 and (high_identity_num >= len(query_records[query_name])-1):
#             if complete_alignment_num < 5 and high_identity_num >= 1:
#                 segmental_duplication_repeatIds.append(query_name)
#             else:
#                 multi_mapping_repeatIds.append(query_name)
#         else:
#             unmapped_repeatIds.append(query_name)
#
#     return unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds, segmental_duplication_repeatIds
#
# def get_alignment_info_v3(sam_paths, repeats_file):
#     repeat_contignames, repeat_contigs = read_fasta(repeats_file)
#     mapping_repeatIds = {}
#     query_records = {}
#     for sam_path in sam_paths:
#         samfile = pysam.AlignmentFile(sam_path, "rb")
#         for read in samfile.fetch():
#             if read.is_unmapped:
#                 continue
#             query_name = read.query_name
#             reference_name = read.reference_name
#             cigar = read.cigartuples
#             cigarstr = read.cigarstring
#             NM_tag = 0
#             try:
#                 NM_tag = read.get_tag('NM')
#             except KeyError:
#                 NM_tag = -1
#             identity = compute_identity(cigarstr, NM_tag, 'BLAST')
#             identity = float(identity) * 100
#             is_reverse = read.is_reverse
#             alignment_len = read.query_alignment_length
#             q_start = int(read.query_alignment_start)
#             q_end = int(read.query_alignment_end)
#             t_start = int(read.reference_start)
#             t_end = int(read.reference_end)
#
#             if not query_records.__contains__(query_name):
#                 query_records[query_name] = []
#             records = query_records[query_name]
#             records.append((query_name, reference_name, cigar, cigarstr, is_reverse, alignment_len, identity, t_start, t_end))
#             query_records[query_name] = records
#
#     query_position = {}
#     for query_name in query_records.keys():
#         complete_alignment_num = 0
#         high_identity_num = 0
#         query_seq = repeat_contigs[query_name]
#         query_len = len(query_seq)
#         for i, record in enumerate(query_records[query_name]):
#             reference_name = record[1]
#             cigar = record[2]
#             cigarstr = str(record[3])
#             alignment_len = record[5]
#             identity = record[6]
#             t_start = record[7]
#             t_end = record[8]
#             if float(alignment_len) / query_len >= 0.8 and identity >= 80:
#                 complete_alignment_num += 1
#                 if identity >= 90:
#                     high_identity_num += 1
#             if t_start > t_end:
#                 tmp = t_end
#                 t_end = t_start
#                 t_start = tmp
#             if not query_position.__contains__(reference_name):
#                 query_position[reference_name] = []
#             same_chr_seq = query_position[reference_name]
#             same_chr_seq.append((query_name, t_start, t_end))
#             query_position[reference_name] = same_chr_seq
#         mapping_repeatIds[query_name] = (complete_alignment_num, query_len)
#     new_mapping_repeatIds = {k: v for k, v in sorted(mapping_repeatIds.items(), key=lambda item: (-item[1][1], -item[1][0]))}
#
#     return new_mapping_repeatIds, query_position

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

    # jobs = []
    # pool = multiprocessing.Pool(processes=reduce_partitions_num)
    # for partiton_index in segments_cluster.keys():
    #     cur_segments = segments_cluster[partiton_index]
    #     #future = pool.map_async(generate_candidate_repeats_v2, iterable=args)
    #     res = pool.apply_async(generate_candidate_repeats_v2, (cur_segments, k_num, unique_kmer_map, partiton_index, fault_tolerant_bases,))
    #     jobs.append(res)
    # pool.close()
    # pool.join()
    #
    # repeat_dict = {}
    # for res in jobs:
    #     cur_repeat_dict = res.get()
    #     for ref_name in cur_repeat_dict.keys():
    #         parts = ref_name.split('$')
    #         true_ref_name = parts[0]
    #         start_pos = int(parts[1])
    #         if not repeat_dict.__contains__(true_ref_name):
    #             repeat_dict[true_ref_name] = []
    #         new_repeat_list = repeat_dict[true_ref_name]
    #         cur_repeat_list = cur_repeat_dict[ref_name]
    #         for repeat_item in cur_repeat_list:
    #             new_repeat_item = (start_pos + repeat_item[0], start_pos + repeat_item[1], repeat_item[2])
    #             new_repeat_list.append(new_repeat_item)


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

def TSDsearch_v2(orig_seq, tir_start, tir_end, TSD_set, plant):
    plant = int(plant)
    # 2->TA或者animal/fungi中的5'-CCC...GGG-3', 3-> plant中的5'-CACT(A/G)...(C/T)AGTG-3' 或者是 （TAA或TTA）, 4-> TTAA,
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


def TSDsearch_ltr(orig_seq, ltr_start, ltr_end, TSD_set):
    LTR_TSDs = [6, 5, 4]  # LTR:绝大部分是 5-'TG...CA-3'
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
            #如果是TG..CA，允许TSD有1bp误差
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
        # 如果是TG..CA，允许TSD有1bp误差
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
    # 综合考虑TSD len, tir len, TE len等多个因素
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
    # 综合考虑TSD len, tir len, TE len等多个因素
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
    # 综合考虑拷贝数, TSD len等多个因素,占比(6/4)
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

def file_exist(resut_file):
    if os.path.exists(resut_file) and os.path.getsize(resut_file) > 0:
        if resut_file.endswith('.fa') or resut_file.endswith('.fasta'):
            names, contigs = read_fasta(resut_file)
            if len(contigs) > 0:
                return True
            else:
                return False
        else:
            line_count = 0
            with open(resut_file, 'r') as f_r:
                for line in f_r:
                    if line.startswith('#'):
                        continue
                    line_count += 1
                    if line_count > 0:
                        return True
            f_r.close()
            return False
    else:
        return False


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
                    # 提取序列的5'的100bp，然后判断是否是串联重复
                    start_seq = seq[0:100]
                    if float(start_seq.count('N')) / len(start_seq) >= tandem_region_cutoff:
                        continue
            elif TE_type == 'tir':
                #提取序列的5'和3'的20bp，然后判断是否是串联重复
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

        #合并序列片段TE
        pure_copies = []
        for cur_copy in copies:
            # 去掉copy[2] < 2 且不是全长拷贝的copies，没有其余的拷贝支持它是一个重复。
            if cur_copy[2] < 2 and float(cur_copy[1] - cur_copy[0]) / query_len < 0.95:
                continue
            is_seg = False
            for i, pure_copy in enumerate(pure_copies):
                # 因为真实的TE两端不会是重复，因此尽量取最长的边界
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

        #合并序列片段TE
        pure_copies = []
        for cur_copy in copies:
            # 去掉copy[2] < 2 且不是全长拷贝的copies，没有其余的拷贝支持它是一个重复。
            if cur_copy[2] < 2 and float(cur_copy[1] - cur_copy[0]) / query_len < 0.95:
                continue
            is_seg = False
            for i, pure_copy in enumerate(pure_copies):
                # 因为真实的TE两端不会是重复，因此尽量取最长的边界
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

def get_longest_repeats_v3(repeats_path, fixed_extend_base_threshold, max_single_repeat_len, threads):
    split_repeats_path = repeats_path[0]
    original_repeats_path = repeats_path[1]
    blastn2Results_path = repeats_path[2]
    tmp_blast_dir = repeats_path[3]


    align_command = 'blastn -db ' + original_repeats_path + ' -num_threads ' \
                    + str(1) + ' -query ' + split_repeats_path + ' -outfmt 6 > ' + blastn2Results_path
    os.system(align_command)

    print('coarse alignment -- alignment finished:' + str(split_repeats_path))

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
            if identity < 80:
                continue
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end))
    f_r.close()
    print('coarse alignment -- file open finished:' + str(split_repeats_path))

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

        #合并序列片段TE
        pure_copies = []
        for cur_copy in copies:
            # 去掉copy[2] < 2 且不是全长拷贝的copies，没有其余的拷贝支持它是一个重复。
            if cur_copy[2] < 2 and float(cur_copy[1] - cur_copy[0]) / query_len < 0.95:
                continue
            is_seg = False
            for i, pure_copy in enumerate(pure_copies):
                # 因为真实的TE两端不会是重复，因此尽量取最长的边界
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

    print('coarse alignment -- analyze finished:' + str(split_repeats_path))
    # print(longest_repeats)
    return longest_repeats, keep_longest_query

def flanking_seq(longest_repeats_path, longest_repeats_flanked_path, reference, flanking_len):
    seq_names, seq_contigs = read_fasta(longest_repeats_path)
    ref_names, ref_contigs = read_fasta(reference)
    flanked_contigs = {}
    node_index = 0
    for name in seq_names:
        # N_126-len_955-ref_NC_029267.1_7568408_7569363
        ref_info = name.split('-ref_')[1].split('-')
        ref_name = ref_info[0]
        ref_start = int(ref_info[1]) + 1
        ref_end = int(ref_info[2])
        ref_seq = ref_contigs[ref_name]
        if ref_start - 1 - flanking_len < 0 or ref_end + flanking_len > len(ref_seq):
            continue
        flanked_seq = ref_seq[ref_start - 1 - flanking_len: ref_end + flanking_len]
        flanked_contigs[name] = flanked_seq
    store_fasta(flanked_contigs, longest_repeats_flanked_path)

def determine_repeat_boundary_v1(repeats_path, longest_repeats_path, blast_program_dir,
                                 fixed_extend_base_threshold, max_single_repeat_len, tmp_output_dir, threads):
    repeatNames, repeatContigs = read_fasta(repeats_path)
    # parallel
    tmp_blast_dir = tmp_output_dir + '/longest_repeats_blast'
    os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    #(repeat_dir, repeat_filename) = os.path.split(repeats_path)

    # makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + repeats_path + ' > /dev/null 2>&1'
    # os.system(makedb_command)

    ## before method
    # repeat_files = []
    # data_partitions = PET(list(repeatContigs.items()), threads)
    # for partition_index, data_partition in enumerate(data_partitions):
    #     single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
    #     if not os.path.exists(single_tmp_dir):
    #         os.makedirs(single_tmp_dir)
    #     split_repeat_file = single_tmp_dir + '/repeats_split.fa'
    #     store2file(data_partition, split_repeat_file)
    #     repeat_files.append((split_repeat_file, repeats_path,
    #                          single_tmp_dir + '/repeat.pairwise.out'))

    # 2022-12-22 method
    # repeats.fa通常会有几十上百M，只对query切割进行比对还是会有某个query的输出占据几个G的情况，有可能会导致内存溢出。
    # 因此我们对query和subject同时进行切分比对，完成后将query对应的output合并
    data_partitions = PET(list(repeatContigs.items()), threads)
    for partition_index, data_partition in enumerate(data_partitions):
        subject_tmp_dir = tmp_blast_dir + '/subject'
        if not os.path.exists(subject_tmp_dir):
            os.makedirs(subject_tmp_dir)
        split_subject_file = subject_tmp_dir + '/' + str(partition_index) + '.fa'
        store2file(data_partition, split_subject_file)

        makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + split_subject_file + ' > /dev/null 2>&1'
        os.system(makedb_command)


    repeat_files = []
    file_index = 0
    cur_seq_index = 0
    cur_contigs = {}
    for name in repeatNames:
        cur_contigs[name] = repeatContigs[name]
        cur_seq_index += 1
        # 获取query_name中包含的染色体名称
        # parts = name.split('-s_')[1].split('-')
        # chr_name = parts[0]

        if cur_seq_index >= 10:
            split_repeat_file = tmp_blast_dir + '/' + str(file_index) + '.fa'
            store_fasta(cur_contigs, split_repeat_file)
            output_file = tmp_blast_dir + '/' + str(file_index) + '.out'
            repeat_files.append((split_repeat_file, repeats_path, output_file, tmp_blast_dir))
            cur_contigs = {}
            file_index += 1
            cur_seq_index = 0
    if len(cur_contigs) > 0:
        split_repeat_file = tmp_blast_dir + '/' + str(file_index) + '.fa'
        store_fasta(cur_contigs, split_repeat_file)
        output_file = tmp_blast_dir + '/' + str(file_index) + '.out'
        repeat_files.append((split_repeat_file, repeats_path, output_file, tmp_blast_dir))

    ex = ProcessPoolExecutor(threads)
    longest_repeats = {}
    keep_longest_query = {}
    jobs = []
    for file in repeat_files:
        #为了减少内存，只传递需要的reference sequence
        job = ex.submit(get_longest_repeats_v1, file, blast_program_dir, fixed_extend_base_threshold,
                        max_single_repeat_len, threads)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        cur_longest_repeats, cur_keep_longest_query = job.result()
        for query_name in cur_longest_repeats.keys():
            longest_repeats[query_name] = cur_longest_repeats[query_name]
        for query_name in cur_keep_longest_query.keys():
            keep_longest_query[query_name] = cur_keep_longest_query[query_name]

    # longest_repeats_file = tmp_output_dir + '/longest_repeats.csv'
    # with open(longest_repeats_file, 'w') as f_save:
    #     for query_name in keep_longest_query.keys():
    #         longest_queries = keep_longest_query[query_name]
    #         for query in longest_queries:
    #             f_save.write(query_name+'\t'+str(query[0])+'\t'+str(query[1])+'\t'+str(query[2])+'\t'+str(query[3])+'\t'+str(query[4])+'\n')

    # # store longest_repeats for testing
    # longest_repeats_file = tmp_output_dir + '/longest_repeats.csv'
    # with codecs.open(longest_repeats_file, 'w', encoding='utf-8') as f:
    #     json.dump(longest_repeats, f)

    node_index = 0
    with open(longest_repeats_path, 'w') as f_save:
        for query_name in longest_repeats.keys():
            seqs_tuples = longest_repeats[query_name]
            for longest_seq in seqs_tuples:
                # orig_start_pos = longest_seq[0]
                # orig_end_pos = longest_seq[1]
                # chr_name = longest_seq[2]
                # seq_ref_start = longest_seq[3]
                # seq_ref_end = longest_seq[4]
                seq = longest_seq[5]
                # 如果有连续10个以上的N过滤掉
                if len(seq) >= 100 and not seq.__contains__('NNNNNNNNNN'):
                    f_save.write('>N_' + str(node_index) + '-len_' + str(len(seq)) + '\n' + seq + '\n')
                    # f_save.write('>N_' + str(node_index) + '-len_' + str(len(seq)) +
                    #              '-ref_' + chr_name + '-' + str(seq_ref_start) + '-' + str(seq_ref_end) + '\n' + seq + '\n')
                    node_index += 1
    f_save.close()
    return longest_repeats_path

def determine_repeat_boundary_v2(repeats_path, longest_repeats_path, blast_program_dir,
                                 fixed_extend_base_threshold, max_single_repeat_len, tmp_output_dir, threads):
    repeatNames, repeatContigs = read_fasta(repeats_path)
    # parallel
    tmp_blast_dir = tmp_output_dir + '/longest_repeats_blast'
    os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    # 2022-12-22 method
    # repeats.fa通常会有几十上百M，只对query切割进行比对还是会有某个query的输出占据几个G的情况，有可能会导致内存溢出。
    # 因此我们对query和subject同时进行切分比对，完成后将query对应的output合并
    subject_list = []
    subject_tmp_dir = tmp_blast_dir + '/subject'
    if not os.path.exists(subject_tmp_dir):
        os.makedirs(subject_tmp_dir)
    data_partitions = PET(list(repeatContigs.items()), threads)
    for partition_index, data_partition in enumerate(data_partitions):
        split_subject_file = subject_tmp_dir + '/' + str(partition_index) + '.fa'
        store2file(data_partition, split_subject_file)
        subject_list.append(split_subject_file)
        makedb_command = blast_program_dir + '/bin/makeblastdb -dbtype nucl -in ' + split_subject_file + ' > /dev/null 2>&1'
        os.system(makedb_command)


    repeat_files = []
    file_index = 0
    cur_seq_index = 0
    cur_contigs = {}
    for name in repeatNames:
        cur_contigs[name] = repeatContigs[name]
        cur_seq_index += 1

        if cur_seq_index >= 10:
            split_repeat_file = tmp_blast_dir + '/' + str(file_index) + '.fa'
            store_fasta(cur_contigs, split_repeat_file)
            repeat_files.append(split_repeat_file)
            #output_file = tmp_blast_dir + '/' + str(file_index) + '.out'
            #repeat_files.append((split_repeat_file, subject_list, output_file, tmp_blast_dir))
            cur_contigs = {}
            file_index += 1
            cur_seq_index = 0
    if len(cur_contigs) > 0:
        split_repeat_file = tmp_blast_dir + '/' + str(file_index) + '.fa'
        store_fasta(cur_contigs, split_repeat_file)
        repeat_files.append(split_repeat_file)
        #output_file = tmp_blast_dir + '/' + str(file_index) + '.out'
        #repeat_files.append((split_repeat_file, subject_list, output_file, tmp_blast_dir))

    #记录repeat_file对应的所有输出文件.out,比对完之后需要合并然后分析
    output_tmp_dir = tmp_blast_dir + '/out'
    if not os.path.exists(output_tmp_dir):
        os.makedirs(output_tmp_dir)
    out_file_index = 0
    repeat_outs = {}
    repeat_combination = []
    for repeat_file in repeat_files:
        if not repeat_outs.__contains__(repeat_file):
            repeat_outs[repeat_file] = []
        repeat_out_list = repeat_outs[repeat_file]
        for split_subject_file in subject_list:
            output_file = output_tmp_dir + '/' + str(out_file_index) + '.out'
            repeat_combination.append(repeat_file, split_subject_file, output_file)
            repeat_out_list.append(output_file)

    #进行pairwise比对
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in repeat_combination:
        # 为了减少内存，只传递需要的reference sequence
        job = ex.submit(pairwise_alignment, file, blast_program_dir)
        jobs.append(job)
    ex.shutdown(wait=True)
    finished_outs = []
    for job in as_completed(jobs):
        cur_blastn2Results_path = job.result()
        finished_outs.append(cur_blastn2Results_path)

    #合并repeat对应的out文件
    merged_outputs = {}
    for repeat_file in repeat_outs.keys():
        repeat_out_list = repeat_outs[repeat_file]
        merged_output = repeat_file + '.out'
        for repeat_out in repeat_out_list:
            if not os.path.exists(repeat_out):
                continue
            os.system('cat '+repeat_out+' >> '+merged_output)
        merged_outputs[repeat_file] = merged_output

    #将合并的output交给get_longest_repeats_v2函数处理
    ex = ProcessPoolExecutor(threads)
    longest_repeats = {}
    keep_longest_query = {}
    jobs = []
    for repeat_file in merged_outputs.keys():
        merged_output = merged_outputs[repeat_file]
        job = ex.submit(get_longest_repeats_v2, repeat_file, merged_output, fixed_extend_base_threshold,
                        max_single_repeat_len)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        cur_longest_repeats, cur_keep_longest_query = job.result()
        for query_name in cur_longest_repeats.keys():
            longest_repeats[query_name] = cur_longest_repeats[query_name]
        for query_name in cur_keep_longest_query.keys():
            keep_longest_query[query_name] = cur_keep_longest_query[query_name]

    # longest_repeats_file = tmp_output_dir + '/longest_repeats.csv'
    # with open(longest_repeats_file, 'w') as f_save:
    #     for query_name in keep_longest_query.keys():
    #         longest_queries = keep_longest_query[query_name]
    #         for query in longest_queries:
    #             f_save.write(query_name+'\t'+str(query[0])+'\t'+str(query[1])+'\t'+str(query[2])+'\t'+str(query[3])+'\t'+str(query[4])+'\n')

    # # store longest_repeats for testing
    # longest_repeats_file = tmp_output_dir + '/longest_repeats.csv'
    # with codecs.open(longest_repeats_file, 'w', encoding='utf-8') as f:
    #     json.dump(longest_repeats, f)

    node_index = 0
    with open(longest_repeats_path, 'w') as f_save:
        for query_name in longest_repeats.keys():
            seqs_tuples = longest_repeats[query_name]
            for longest_seq in seqs_tuples:
                orig_start_pos = longest_seq[0]
                orig_end_pos = longest_seq[1]
                chr_name = longest_seq[2]
                seq_ref_start = longest_seq[3]
                seq_ref_end = longest_seq[4]
                seq = longest_seq[5]
                # 如果有连续10个以上的N过滤掉
                if len(seq) >= 100 and not seq.__contains__('NNNNNNNNNN'):
                    f_save.write('>N_' + str(node_index) + '-len_' + str(len(seq)) +
                                 '-ref_' + chr_name + '-' + str(seq_ref_start) + '-' + str(seq_ref_end) + '\n' + seq + '\n')
                    node_index += 1
    f_save.close()
    return longest_repeats_path

def determine_repeat_boundary_v3(repeats_path, longest_repeats_path, fixed_extend_base_threshold, max_single_repeat_len, tmp_output_dir, threads, ref_index, log):
    repeatNames, repeatContigs = read_fasta(repeats_path)
    # parallel
    tmp_blast_dir = tmp_output_dir + '/longest_repeats_blast_'+str(ref_index)
    os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    makedb_command = 'makeblastdb -dbtype nucl -in ' + repeats_path + ' > /dev/null 2>&1'
    os.system(makedb_command)

    seq_num = 1
    repeat_files = []
    file_index = 0
    cur_seq_index = 0
    cur_contigs = {}
    for name in repeatNames:
        cur_contigs[name] = repeatContigs[name]
        cur_seq_index += 1
        if cur_seq_index >= seq_num:
            split_repeat_file = tmp_blast_dir + '/' + str(file_index) + '.fa'
            store_fasta(cur_contigs, split_repeat_file)
            output_file = tmp_blast_dir + '/' + str(file_index) + '.out'
            repeat_files.append((split_repeat_file, repeats_path, output_file, tmp_blast_dir))
            cur_contigs = {}
            file_index += 1
            cur_seq_index = 0
    if len(cur_contigs) > 0:
        split_repeat_file = tmp_blast_dir + '/' + str(file_index) + '.fa'
        store_fasta(cur_contigs, split_repeat_file)
        output_file = tmp_blast_dir + '/' + str(file_index) + '.out'
        repeat_files.append((split_repeat_file, repeats_path, output_file, tmp_blast_dir))

    log.logger.debug('coarse alignment -- total file size:'+str(len(repeat_files)))
    ex = ProcessPoolExecutor(threads)
    longest_repeats = {}
    keep_longest_query = {}
    jobs = []
    for file in repeat_files:
        job = ex.submit(get_longest_repeats_v3, file, fixed_extend_base_threshold,
                        max_single_repeat_len, threads)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        cur_longest_repeats, cur_keep_longest_query = job.result()
        for query_name in cur_longest_repeats.keys():
            longest_repeats[query_name] = cur_longest_repeats[query_name]
        for query_name in cur_keep_longest_query.keys():
            keep_longest_query[query_name] = cur_keep_longest_query[query_name]

    # longest_repeats_file = tmp_output_dir + '/longest_repeats.csv'
    # with open(longest_repeats_file, 'w') as f_save:
    #     for query_name in keep_longest_query.keys():
    #         longest_queries = keep_longest_query[query_name]
    #         for query in longest_queries:
    #             f_save.write(query_name+'\t'+str(query[0])+'\t'+str(query[1])+'\t'+str(query[2])+'\t'+str(query[3])+'\t'+str(query[4])+'\n')

    # # store longest_repeats for testing
    # longest_repeats_file = tmp_output_dir + '/longest_repeats.csv'
    # with codecs.open(longest_repeats_file, 'w', encoding='utf-8') as f:
    #     json.dump(longest_repeats, f)

    node_index = 0
    with open(longest_repeats_path, 'w') as f_save:
        for query_name in longest_repeats.keys():
            seqs_tuples = longest_repeats[query_name]
            for longest_seq in seqs_tuples:
                orig_start_pos = longest_seq[0]
                orig_end_pos = longest_seq[1]
                chr_name = longest_seq[2]
                seq_ref_start = longest_seq[3]
                seq_ref_end = longest_seq[4]
                seq = longest_seq[5]
                # 如果有连续10个以上的N过滤掉
                if len(seq) >= 100 and not seq.__contains__('NNNNNNNNNN'):
                    f_save.write('>N_' + str(node_index) + '-len_' + str(len(seq)) +
                                 '-ref_' + chr_name + '-' + str(seq_ref_start) + '-' + str(seq_ref_end) + '\n' + seq + '\n')
                    node_index += 1
    f_save.close()
    return longest_repeats_path


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

            tir_start = flanking_len + 1  # (坐标都是以1开始的，start和end都是取到的)
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

def multi_process_tsd(longest_repeats_flanked_path, tir_tsd_path, tir_tsd_filter_dup_path, tir_tsd_dir, flanking_len, threads, TRsearch_dir, plant, reference):
    os.system('rm -rf '+tir_tsd_dir)
    if not os.path.exists(tir_tsd_dir):
        os.makedirs(tir_tsd_dir)

    seq_names, seq_contigs = read_fasta(longest_repeats_flanked_path)

    # (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq)
    segments_cluster = divided_array(list(seq_contigs.items()), threads)

    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    for partition_index, cur_segments in enumerate(segments_cluster):
        obj = ex.submit(search_confident_tir_batch, cur_segments, flanking_len, tir_tsd_dir, TRsearch_dir, partition_index, plant)
        objs.append(obj)
        job_id += 1
    ex.shutdown(wait=True)
    candidate_TIRs = {}
    for obj in as_completed(objs):
        cur_candidate_TIRs = obj.result()
        candidate_TIRs.update(cur_candidate_TIRs)
    store_fasta(candidate_TIRs, tir_tsd_path)

    # 重新比对到基因组，获取候选TIR拷贝数
    seq_copynum = {}
    temp_dir = tir_tsd_dir + '/tir_temp'
    all_copies = multi_process_align_and_get_copies(tir_tsd_path, reference,
                                                    temp_dir, 'tir', threads)
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        seq_copynum[query_name] = len(copies)

    # 按照query_name进行分组，同一组里只取一条序列，即拷贝数和TSD综合最优的那一条
    # 对all_copies_out_contigs按照query_name进行分组
    # group_copies_contigs -> {query_name: {name: seq}}
    group_copies_contigs = {}
    for cur_name in candidate_TIRs.keys():
        query_name = cur_name.split('-C_')[0]
        if not group_copies_contigs.__contains__(query_name):
            group_copies_contigs[query_name] = {}
        cur_copies_out_contigs = group_copies_contigs[query_name]
        cur_copies_out_contigs[cur_name] = candidate_TIRs[cur_name]

    filter_dup_itr_contigs = {}
    for query_name in group_copies_contigs.keys():
        cur_copies_out_contigs = group_copies_contigs[query_name]
        # 选择拷贝数和TSD综合最优的那一条
        cur_contigs = filter_dup_itr_v1(cur_copies_out_contigs, seq_copynum)
        filter_dup_itr_contigs.update(cur_contigs)

    store_fasta(filter_dup_itr_contigs, tir_tsd_filter_dup_path)

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

def multi_process_tsd_v1(all_copies, flanking_len, threads, plant):
    # (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq)
    segments_cluster = divided_array(list(all_copies.items()), threads)

    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    for partition_index, cur_segments in enumerate(segments_cluster):
        obj = ex.submit(search_confident_tir_batch_v1, cur_segments, flanking_len, plant)
        objs.append(obj)
        job_id += 1
    ex.shutdown(wait=True)
    confident_copies = {}
    for obj in as_completed(objs):
        cur_confident_copies = obj.result()
        confident_copies.update(cur_confident_copies)
    return confident_copies

def search_confident_tir_batch_v1(cur_segments, flanking_len, plant):
    tsd_search_distance = flanking_len
    confident_copies = {}
    for item in cur_segments:
        query_name = item[0]
        copies = item[1]
        for i, copy in enumerate(copies):
            orig_seq = copy[4]
            tsd = copy[5]
            tir_start = flanking_len + 1
            tir_end = flanking_len + copy[3]
            copy_ref_name = copy[0]
            copy_ref_start = copy[1]
            copy_ref_end = copy[2]
            copy_ref_info = str(copy_ref_name) + ':' + str(copy_ref_start) + '_' + str(copy_ref_end)
            new_copy = search_confident_tir_v1(orig_seq, tir_start, tir_end, copy_ref_name, copy_ref_start, copy_ref_end, tsd_search_distance, tsd, plant)
            # copy -> (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq, tsd)
            if new_copy is not None:
                if not confident_copies.__contains__(query_name):
                    confident_copies[query_name] = []
                copy_list = confident_copies[query_name]
                copy_list.append(new_copy)
    return confident_copies

def multi_process_tsd_v2(all_copies, flanking_len, threads, plant):
    # (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq)
    segments_cluster = divided_array(list(all_copies.items()), threads)

    job_id = 0
    ex = ProcessPoolExecutor(threads)
    objs = []
    for partition_index, cur_segments in enumerate(segments_cluster):
        obj = ex.submit(search_confident_tir_batch_v2, cur_segments, flanking_len, plant)
        objs.append(obj)
        job_id += 1
    ex.shutdown(wait=True)
    confident_copies = {}
    for obj in as_completed(objs):
        cur_confident_copies = obj.result()
        confident_copies.update(cur_confident_copies)
    return confident_copies

def search_confident_tir_batch_v2(cur_segments, flanking_len, plant):
    tsd_search_distance = flanking_len
    confident_copies = {}
    new_confident_copies = {}
    for item in cur_segments:
        query_name = item[0]
        copies = item[1]
        #记录TSD长度出现的次数
        TSD_len_count = {}
        for i, copy in enumerate(copies):
            orig_seq = copy[4]
            tsd = copy[5]
            tir_start = flanking_len + 1
            tir_end = flanking_len + copy[3]
            copy_ref_name = copy[0]
            copy_ref_start = copy[1]
            copy_ref_end = copy[2]
            copy_ref_info = str(copy_ref_name) + ':' + str(copy_ref_start) + '_' + str(copy_ref_end)

            new_copy = search_confident_tir_v2(orig_seq, tir_start, tir_end, copy_ref_name, copy_ref_start, copy_ref_end, tsd_search_distance, tsd, plant)
            # copy -> (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq, tsd)
            if new_copy is not None:
                if not confident_copies.__contains__(query_name):
                    confident_copies[query_name] = []
                copy_list = confident_copies[query_name]
                copy_list.append(new_copy)
                new_copy_tsd_len = len(new_copy[5])
                if not TSD_len_count.__contains__(new_copy_tsd_len):
                    TSD_len_count[new_copy_tsd_len] = 0
                count = TSD_len_count[new_copy_tsd_len]
                count += 1
                TSD_len_count[new_copy_tsd_len] = count

        # 判断当前拷贝中哪一个长度的TSD出现次数最多，那么该TSD对应的拷贝，最有可能是真实的TIR
        max_count = 0
        max_count_tsd_len = -1
        for tsd_len in TSD_len_count.keys():
            count = TSD_len_count[tsd_len]
            if count > max_count:
                max_count = count
                max_count_tsd_len = tsd_len
        # 只取具有最多次数TSD长度的拷贝
        if max_count_tsd_len != -1:
            copy_list = confident_copies[query_name]
            for copy in copy_list:
                cur_tsd_len = len(copy[5])
                #如果最高次数TSD长度是2，则要求拷贝中的TSD

                if cur_tsd_len == max_count_tsd_len:
                    if not new_confident_copies.__contains__(query_name):
                        new_confident_copies[query_name] = []
                    copy_list = new_confident_copies[query_name]
                    copy_list.append(copy)
    return new_confident_copies

def search_confident_tir_v1(orig_seq, raw_tir_start, raw_tir_end, copy_ref_name, copy_ref_start, copy_ref_end, tsd_search_distance, tsd, plant):
    orig_seq_len = len(orig_seq)
    tir_starts = []
    tir_ends = []

    for i in range(raw_tir_start - tsd_search_distance, raw_tir_start + tsd_search_distance + 1):
        if i >= 1 and i <= orig_seq_len:
            tir_starts.append(i)

    for i in range(raw_tir_end - tsd_search_distance, raw_tir_end + tsd_search_distance + 1):
        if i >= 1 and i <= orig_seq_len:
            tir_ends.append(i)

    boundary_set = set()
    for tir_start in tir_starts:
        for tir_end in tir_ends:
            left_tsd_seq, right_tsd_seq = TSDsearch_v3(orig_seq, tir_start, tir_end, tsd, plant)
            if left_tsd_seq != '':
                boundary_set.add((tir_start, tir_end, left_tsd_seq))

    #取最靠近原始边界，且具有相同长度TSD的拷贝当做可靠拷贝
    boundary_set = sorted(boundary_set, key=lambda x: abs(x[0] - raw_tir_start) + abs(x[1] - raw_tir_end))
    new_copy = None
    if len(boundary_set) > 0:
        new_boundary = boundary_set[0]
        tsd = new_boundary[2]
        tir_start = new_boundary[0]
        tir_end = new_boundary[1]
        copy_len = tir_end - tir_start + 1
        copy_seq = orig_seq[tir_start-1: tir_end]
        copy_ref_start += tir_start - raw_tir_start
        copy_ref_end += tir_end - raw_tir_end
        new_copy = (copy_ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq, tsd)

    return new_copy

def search_confident_tir_v2(orig_seq, raw_tir_start, raw_tir_end, copy_ref_name, copy_ref_start, copy_ref_end, tsd_search_distance, tsd, plant):
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
    #boundary_set = set()
    for tir_start in tir_starts:
        for tir_end in tir_ends:
            TSDsearch_v2(orig_seq, tir_start, tir_end, TSD_set, plant)
            # left_tsd_seq, right_tsd_seq = TSDsearch_v3(orig_seq, tir_start, tir_end, tsd, plant)
            # if left_tsd_seq != '':
            #     boundary_set.add((tir_start, tir_end, left_tsd_seq))

    #取最靠近原始边界，且具有TSD的拷贝
    # (left_tsd_start, left_tsd_end, left_tsd_seq, right_tsd_start, right_tsd_end, right_tsd_seq, tir_start, tir_end, tir_len)
    # 按照tir_start, tir_end与原始边界的距离进行排序，越近的排在前面
    TSD_set = sorted(TSD_set, key=lambda x: abs(x[6] - raw_tir_start) + abs(x[7] - raw_tir_end))
    # boundary_set = sorted(boundary_set, key=lambda x: abs(x[0] - raw_tir_start) + abs(x[1] - raw_tir_end))
    new_copy = None
    if len(TSD_set) > 0:
        new_boundary = TSD_set[0]
        tsd = new_boundary[2]
        tir_start = new_boundary[6]
        tir_end = new_boundary[7]
        copy_len = tir_end - tir_start + 1
        copy_seq = orig_seq[tir_start-1: tir_end]
        copy_ref_start += tir_start - raw_tir_start
        copy_ref_end += tir_end - raw_tir_end
        new_copy = (copy_ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq, tsd)

    return new_copy

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
            #如果是取全长拷贝，则需要判断拷贝的前5bp和后5bp是否与原始序列高度相似，只允许1 bp mismatch
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

def get_query_copies(cur_segments, query_contigs, subject_path, query_coverage, subject_coverage, query_fixed_extend_base_threshold, subject_fixed_extend_base_threshold):
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
                        if (frag[2] - exist_frag[3] < subject_fixed_extend_base_threshold):
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
                        if (exist_frag[3] - frag[2] < subject_fixed_extend_base_threshold):
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

        # if query_name.__contains__('N_21202-len_1496-ref_NC_029260.1-22739190-22740686-C_29-tsd_TTTTTTCAT-distance_18'):
        #     print('here')


        # we now consider, we should take some sequences from longest_queries to represent this query sequence.
        # we take the longest sequence by length, if the latter sequence overlap with the former sequence largely (50%),
        # continue find next sequence until the ratio of query sequence over 90% or no more sequences.
        longest_queries.sort(key=lambda x: -x[2])
        query_len = len(query_contigs[query_name])
        # query_len = int(query_name.split('-')[1].split('_')[1])
        copies = []
        keeped_copies = set()
        for query in longest_queries:
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
        copies.sort(key=lambda x: abs(x[3]-(x[2]-x[1]+1)))
        all_copies[query_name] = copies
    return all_copies

def get_copies_v1(blastnResults_path, query_path, subject_path, query_coverage=0.99, subject_coverage=0,
               query_fixed_extend_base_threshold=1000, subject_fixed_extend_base_threshold=1000):
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
    all_copies = get_query_copies(cur_segments, query_contigs, subject_path, query_coverage, subject_coverage, query_fixed_extend_base_threshold, subject_fixed_extend_base_threshold)

    return all_copies


def get_copies(blastnResults_path, query_path, subject_path, query_coverage=0.99, subject_coverage=0,
               query_fixed_extend_base_threshold=1000, subject_fixed_extend_base_threshold=1000, threads=48):
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
        obj = ex.submit(get_query_copies, cur_segments, query_contigs, subject_path, query_coverage, subject_coverage, query_fixed_extend_base_threshold, subject_fixed_extend_base_threshold)
        objs.append(obj)
        job_id += 1
    ex.shutdown(wait=True)

    all_copies = {}
    for obj in as_completed(objs):
        cur_copies = obj.result()
        all_copies.update(cur_copies)
    return all_copies

def generate_candidate_ltrs(all_copies, reference, flanking_len):
    # 将同一个query_name对应的相邻拷贝聚在一起，然后生成候选LTR序列
    # copy_cluster -> {ref_name: [(left_ltr_ref_start, left_ltr_ref_end, right_ltr_ref_start, right_ltr_ref_end), ..., ]}
    copy_cluster = {}
    ref_names, ref_contigs = read_fasta(reference)
    for query_name in all_copies.keys():
        copies = list(all_copies[query_name])
        # 将copies按照ref_name, ref_start, ref_end排序
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
                # 序列的长度范围（100-7000），定义同一个序列的两个拷贝之间的起始位置距离（ 1000，15000）
                if (last_copy[3] >= 100 and last_copy[3] <= 7000) and (copy_len >= 100 and copy_len <= 7000):
                    # 当前拷贝与前一个拷贝是没有Overlap的，且之间的距离满足1000<= x <= 15000，则可能是一个LTR element
                    if copy_ref_start > last_copy[2] and (copy_ref_start - last_copy[1] >= 1000) and (
                            copy_ref_start - last_copy[1] <= 15000):
                        if not copy_cluster.__contains__(ref_name):
                            copy_cluster[ref_name] = []
                        candidate_ltrs = copy_cluster[ref_name]
                        candidate_ltrs.append((last_copy[1], last_copy[2], copy_ref_start, copy_ref_end))
            last_copy = copy

    flanked_ltr_candidates = {}
    # flanked_ltr_candidates -> {ref_name: [(left_ltr_ref_start, left_ltr_ref_end, right_ltr_ref_start, right_ltr_ref_end, flanked_seq), ..., ]}
    # 将所有的candidate_ltr扩展100bp，以确保能搜索到TSD
    for ref_name in copy_cluster.keys():
        candidate_ltrs = copy_cluster[ref_name]
        # 按照ltr_start, ltr_end进行排序，过滤掉重复（与上一个ltr的起始位置均小于10）的LTR。
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
                    + str(1) + ' -query ' + split_repeats_path + ' -outfmt 6 > ' + blastn2Results_path
    os.system(align_command)

    return blastn2Results_path

def multiple_alignment_blast_and_get_copies(repeats_path, tools_dir, query_coverage, TE_type, subject_coverage=0):
    split_repeats_path = repeats_path[0]
    ref_db_path = repeats_path[1]
    blastn2Results_path = repeats_path[2]
    all_copies = None
    repeat_names, repeat_contigs = read_fasta(split_repeats_path)
    if len(repeat_contigs) > 0:

        align_command = 'blastn -db ' + ref_db_path + ' -num_threads ' \
                        + str(1) + ' -query ' + split_repeats_path + ' -outfmt 6 > ' + blastn2Results_path
        os.system(align_command)

        # if TE_type == 'ltr':
        #     all_copies = get_copies_v1(blastn2Results_path, split_repeats_path, ref_db_path,
        #                                query_coverage=query_coverage, query_fixed_extend_base_threshold=10000,
        #                                subject_fixed_extend_base_threshold=10000)
        # else:

        all_copies = get_copies_v1(blastn2Results_path, split_repeats_path, ref_db_path, query_coverage=query_coverage, subject_coverage=subject_coverage)

    return all_copies

def run_blast_align(query_path, subject_path, output, flanking_len, flanking_region_distance, flank_align_dir):
    if file_exist(query_path) and file_exist(subject_path):
        blast_db_command = 'makeblastdb -dbtype nucl -in ' + subject_path + ' > /dev/null 2>&1'
        align_command = 'blastn -db ' + subject_path + ' -num_threads ' \
                        + str(1) + ' -query ' + query_path + ' -outfmt 6 > ' + output
        os.system(blast_db_command)
        os.system(align_command)

    delete_names, appeared_names = judge_flank_align(flanking_region_distance, output, flanking_len, flank_align_dir)
    return delete_names, appeared_names


def judge_itr_structure(TSD_set, orig_seq, name, raw_tir_start, raw_tir_end, cur_candidate_TIRs_path, TRsearch_dir, tir_tsd_dir, plant):
    # 具有短tir的TIR有：(hAT, 5-27bp tir, 8bp tsd, len<4kb), (Mutator, long/short tir, 9-11bp tsd, ), (CACTA, 5bp tir, 2-3bp tsd, )

    # (left_tsd_start, left_tsd_end, left_tsd_seq, right_tsd_start, right_tsd_end, right_tsd_seq, tir_start, tir_end, tir_len)
    # 按照tir_start, tir_end与原始边界的距离进行排序，越近的排在前面
    TSD_set = sorted(TSD_set, key=lambda x: abs(x[6] - raw_tir_start) + abs(x[7] - raw_tir_end))

    itr_contigs = {}
    cur_tsd_seq = ''
    # 遍历所有的候选TSD，只取一条
    for i, tsd_info in enumerate(TSD_set):
        left_tsd_start = tsd_info[0]
        left_tsd_end = tsd_info[1]
        left_tsd = tsd_info[2]
        cur_tsd_seq = left_tsd
        # 如果有连续10个以上的N或者搜索的TSD有连续>=2个N，就过滤掉
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

        # 判断这条tir_seq是否具有TIR终端结构
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
            # 通过itrsearch判断是否有TIR结构
            # 搜索 cur_candidate_TIRs 中是否存在tir结构，如果存在就可以跳过判断其他拷贝了
            output = run_itrsearch(TRsearch_dir, tir_seq)
            output = str(output)
            if output != '':
                seq = output.split('\n')[1]
                itr_contigs[name] = seq
                break
    return itr_contigs, cur_tsd_seq


def get_short_tir_contigs(cur_itr_contigs, plant):
    # 保存那些具有短tir结构的序列
    short_itr_contigs = {}
    for name in cur_itr_contigs.keys():
        tir_seq = cur_itr_contigs[name]
        # 具有短tir的TIR有：(hAT, 5-27bp tir, 8bp tsd, len<4kb), (Mutator, long/short tir, 9-11bp tsd, ), (CACTA, 5bp tir, 2-3bp tsd, )
        # 判断这条tir_seq是否具有TIR终端结构
        tir_start = 1
        tir_end = len(tir_seq)
        first_5bp = tir_seq[tir_start - 1: tir_start + 4]
        last_5bp = getReverseSequence(tir_seq[tir_end - 5: tir_end])
        first_3bp = tir_seq[tir_start - 1: tir_start + 2]
        last_3bp = getReverseSequence(tir_seq[tir_end - 3: tir_end])
        tsd_seq = name.split('-tsd_')[1].split('-')[0]
        tsd_len = len(tsd_seq)

        if first_5bp == last_5bp:
            # hAT
            if tsd_len == 8 and len(tir_seq) < 4000:
                # name += ' Length itr=5'
                short_itr_contigs[name] = tir_seq
            # Mutator
            elif tsd_len >= 9 and tsd_len <= 11:
                # name += ' Length itr=5'
                short_itr_contigs[name] = tir_seq
            # CACTA (plant -> CACT[A/G], animal/fungi -> CCC)
            elif plant == 1 and tsd_len == 3 and (first_5bp == 'CACTA' or first_5bp == 'CACTG'):
                # name += ' Length itr=5'
                short_itr_contigs[name] = tir_seq
        elif first_3bp == last_3bp and plant == 0 and tsd_len == 2 and (first_3bp == 'CCC'):
            # name += ' Length itr=3'
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


def search_confident_tir_batch(cur_segments, flanking_len, tir_tsd_dir, TRsearch_dir, partition_index, plant):
    all_copies_itr_contigs = {}
    all_candidate_TIRs_path = tir_tsd_dir + '/' + str(partition_index) + '.fa'
    for item in cur_segments:
        query_name = item[0]
        seq = item[1]

        #如果有连续10个以上的N或者搜索的TSD有连续>=2个N，就过滤掉
        if seq.__contains__('NNNNNNNNNN'):
            continue

        tir_start = flanking_len + 1
        tir_end = len(seq) - flanking_len
        # 寻找所有可能的TSD序列，计算每条序列的边界与原始边界的距离，并存到header里
        tsd_search_distance = flanking_len
        cur_itr_contigs = search_confident_tir_v3(seq, tir_start, tir_end, tsd_search_distance, query_name, plant)
        all_copies_itr_contigs.update(cur_itr_contigs)

    #保存短tir的序列
    short_itr_contigs = get_short_tir_contigs(all_copies_itr_contigs, plant)

    #剩下的序列交由itrsearch去搜索TIR结构
    for name in short_itr_contigs.keys():
        del all_copies_itr_contigs[name]
    store_fasta(all_copies_itr_contigs, all_candidate_TIRs_path)
    all_copies_out, all_copies_log = run_itrsearch(TRsearch_dir, all_candidate_TIRs_path, tir_tsd_dir)
    # 过滤掉终端TIR长度差距过大的序列
    # filter_large_gap_tirs(all_copies_out, all_copies_out)
    all_copies_out_name, all_copies_out_contigs = read_fasta(all_copies_out)
    # 解析itrsearch log文件，提取比对偏移的序列名称
    fake_tirs = get_fake_tirs(all_copies_log)
    #过滤掉可能是fake tir的序列
    for name in all_copies_out_name:
        if name in fake_tirs:
            del all_copies_out_contigs[name]

    all_copies_out_contigs.update(short_itr_contigs)
    return all_copies_out_contigs

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
        #print(query_name)
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

            #如果有连续10个以上的N或者搜索的TSD有连续>=2个N，就过滤掉
            if orig_seq.__contains__('NNNNNNNNNN'):
                continue

            ltr_len = copy[3]-copy[0]+1
            ltr_start = flanking_len + 1
            ltr_end = flanking_len + ltr_len
            # 寻找所有可能的TSD序列，计算每条序列的边界与原始边界的距离，并存到header里
            tsd_search_distance = 50
            cur_ltr_contigs = search_confident_ltr(orig_seq, ltr_start, ltr_end, tsd_search_distance, ref_name, copy_index)
            all_copies_ltr_contigs.update(cur_ltr_contigs)
    return all_copies_ltr_contigs

def search_tsd_ltr_batch(all_copies_ltr_contigs, ltr_tsd_dir, partition_index, TRsearch_dir):
    all_candidate_LTRs_path = ltr_tsd_dir + '/' + str(partition_index) + '.fa'
    #序列交由ltrsearch去搜索TIR结构
    store_fasta(all_copies_ltr_contigs, all_candidate_LTRs_path)
    run_ltrsearch(TRsearch_dir, all_candidate_LTRs_path, ltr_tsd_dir)
    all_copies_out = all_candidate_LTRs_path + '.ltr'
    all_copies_out_name, all_copies_out_contigs = read_fasta_v1(all_copies_out)
    return all_copies_out_contigs

def rename_fasta(input, output):
    names, contigs = read_fasta(input)
    node_index = 0
    with open(output, 'w') as f_save:
        for name in names:
            seq = contigs[name]
            f_save.write('>N_'+str(node_index)+'\n'+seq+'\n')
            node_index += 1
    f_save.close()

def rename_reference(input, output):
    names, contigs = read_fasta(input)
    ref_index = 0
    with open(output, 'w') as f_save:
        for name in names:
            seq = contigs[name]
            new_name = 'chr_'+str(ref_index)
            print(new_name)
            f_save.write('>'+new_name+'\n'+seq+'\n')
            ref_index += 1
    f_save.close()

def search_candidate_ltr(copies_out_contigs):
    #对all_copies_out_contigs进行rename
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
    #对all_copies_out_contigs按照query_name进行分组
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
        # 1. 合并，选择distance最小（相同则取最长）的那条序列,当做这条拷贝的代表序列
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

    # 遍历所有的候选TSD，控制只选择前20条序列
    for i, tsd_info in enumerate(TSD_set):
        if i >= 20:
            break
        left_tsd_start = tsd_info[0]
        left_tsd_end = tsd_info[1]
        left_tsd = tsd_info[2]

        # 如果有连续10个以上的N或者搜索的TSD有连续>=2个N，就过滤掉
        if left_tsd.__contains__('NN'):
            continue

        right_tsd_start = tsd_info[3]
        right_tsd_end = tsd_info[4]
        right_tsd = tsd_info[5]

        ltr_start = tsd_info[6]
        ltr_end = tsd_info[7]
        # tir_contig = {}
        ltr_seq = orig_seq[ltr_start - 1: ltr_end]
        #计算与原始边界的距离
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
    # 按照tir_start, tir_end与原始边界的距离进行排序，越近的排在前面
    TSD_set = sorted(TSD_set, key=lambda x: abs(x[6] - raw_tir_start) + abs(x[7] - raw_tir_end))

    query_dist = []
    # 遍历所有的候选TSD，控制只选择前20条序列
    for i, tsd_info in enumerate(TSD_set):
        left_tsd_start = tsd_info[0]
        left_tsd_end = tsd_info[1]
        left_tsd = tsd_info[2]

        # 如果有连续10个以上的N或者搜索的TSD有连续>=2个N，就过滤掉
        if left_tsd.__contains__('NN'):
            continue

        right_tsd_start = tsd_info[3]
        right_tsd_end = tsd_info[4]
        right_tsd = tsd_info[5]

        tir_start = tsd_info[6]
        tir_end = tsd_info[7]

        # tir_contig = {}
        tir_seq = orig_seq[tir_start - 1: tir_end]
        #计算与原始边界的距离
        distance = abs(tir_start - raw_tir_start) + abs(tir_end - raw_tir_end)

        if len(tir_seq) < 100:
            continue

        # new_query_name = query_name + '-C_' + str(copy_index) + '_' + str(i) + '-tsd_' + left_tsd + '-distance_' + str(distance)
        # itr_contigs[new_query_name] = tir_seq

        #过滤掉具有TG..CA motif的TIR序列，绝大多数应该是假阳性
        if tir_seq[0:2] == 'TG' and tir_seq[-2:] == 'CA':
            continue

        # 过滤掉以TATATATA开头和结束的TIR
        if str(tir_seq).startswith('TATATATA') or str(tir_seq).startswith('ATATATAT'):
            continue

        # 如果以候选TSD定位边界，且tir的起始和结束5bp满足高相似性，则大概率这是一条真实的具有TSD+TIR结构的序列
        tir_start_5base = orig_seq[tir_start - 1: tir_start + 4]
        tir_end_5base = orig_seq[tir_end - 5: tir_end]
        if allow_mismatch(getReverseSequence(tir_start_5base), tir_end_5base, 1):
            new_query_name = query_name + '-C_' + str(i) + '-tsd_' + left_tsd + '-distance_'+ str(distance)
            itr_contigs[new_query_name] = tir_seq
            query_dist.append((new_query_name, distance))

        # # 我们使用itrsearch的比对信息排除了比对偏移的情况，所以不需要5bp限制
        # new_query_name = query_name + '-C_' + str(i) + '-tsd_' + left_tsd + '-distance_' + str(distance)
        # itr_contigs[new_query_name] = tir_seq
        # query_dist.append((new_query_name, distance))

    #取distance最小的top 10
    max_top_num = 10
    query_dist.sort(key=lambda x: x[1])
    top_itr_contigs = {}
    for i, item in enumerate(query_dist):
        if i >= max_top_num:
            break
        query_name = item[0]
        top_itr_contigs[query_name] = itr_contigs[query_name]
    return top_itr_contigs

def search_confident_tir_v3(orig_seq, raw_tir_start, raw_tir_end, tsd_search_distance, query_name, plant):
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
    # 按照tir_start, tir_end与原始边界的距离进行排序，越近的排在前面
    TSD_set = sorted(TSD_set, key=lambda x: abs(x[6] - raw_tir_start) + abs(x[7] - raw_tir_end))

    query_dist = []
    # 遍历所有的候选TSD，控制只选择前20条序列
    for i, tsd_info in enumerate(TSD_set):
        left_tsd_start = tsd_info[0]
        left_tsd_end = tsd_info[1]
        left_tsd = tsd_info[2]

        # 如果有连续10个以上的N或者搜索的TSD有连续>=2个N，就过滤掉
        if left_tsd.__contains__('NN'):
            continue

        right_tsd_start = tsd_info[3]
        right_tsd_end = tsd_info[4]
        right_tsd = tsd_info[5]

        tir_start = tsd_info[6]
        tir_end = tsd_info[7]

        # tir_contig = {}
        tir_seq = orig_seq[tir_start - 1: tir_end]
        #计算与原始边界的距离
        distance = abs(tir_start - raw_tir_start) + abs(tir_end - raw_tir_end)

        if len(tir_seq) < 100:
            continue

        # new_query_name = query_name + '-C_' + str(copy_index) + '_' + str(i) + '-tsd_' + left_tsd + '-distance_' + str(distance)
        # itr_contigs[new_query_name] = tir_seq

        #过滤掉具有TG..CA motif的TIR序列，绝大多数应该是假阳性
        if tir_seq[0:2] == 'TG' and tir_seq[-2:] == 'CA':
            continue

        # 过滤掉以TATATATA开头和结束的TIR
        if str(tir_seq).startswith('TATATATA') or str(tir_seq).startswith('ATATATAT'):
            continue

        # 如果以候选TSD定位边界，且tir的起始和结束5bp满足高相似性，则大概率这是一条真实的具有TSD+TIR结构的序列
        tir_start_5base = orig_seq[tir_start - 1: tir_start + 4]
        tir_end_5base = orig_seq[tir_end - 5: tir_end]
        if allow_mismatch(getReverseSequence(tir_start_5base), tir_end_5base, 1):
            new_query_name = query_name + '-C_' + str(i) + '-tsd_' + left_tsd + '-distance_'+ str(distance)
            itr_contigs[new_query_name] = tir_seq
            query_dist.append((new_query_name, distance))

    top_itr_contigs = {}
    for i, item in enumerate(query_dist):
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

def judge_flank_align(flanking_region_distance, output, flanking_len, flank_align_dir):
    # 按照 orig_query_name 进行分组，每个 orig_query_name 为一个单独的单元
    # query_groups -> {query_name: {subject_name: []}}
    query_groups = {}
    with open(output, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]

            orig_query_name = query_name.split('-c_')[0]
            orig_subject_name = subject_name.split('-c_')[0]

            if orig_query_name != orig_subject_name:
                continue

            # 按照subject_name进行分组；按照分组去确定该subject_name是否该过滤。
            if not query_groups.__contains__(orig_query_name):
                query_groups[orig_query_name] = {}
            subject_groups = query_groups[orig_query_name]

            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            direct = '+'
            if s_start > s_end:
                direct = '-'
            if query_name == subject_name or not query_name.__contains__('-c_0'):
                continue

            query_parts = query_name.split('-c_')
            orig_query_len = int(query_parts[1].split('-')[1])

            subject_parts = subject_name.split('-c_')
            orig_subject_len = int(subject_parts[1].split('-')[1])

            if not subject_groups.__contains__(subject_name):
                subject_groups[subject_name] = []
            group = subject_groups[subject_name]
            group.append((q_start, q_end, orig_query_len, s_start, s_end, orig_subject_len, flanking_len, flanking_region_distance, direct, query_name, subject_name))
    f_r.close()

            # # 判断C0与C1,C2...等拷贝的比对情况，如果有flanking区域包含在比对区域内，那么这条拷贝应该被抛弃，如果所有拷贝被抛弃，则该条序列应该是假阳性。
            # if direct == '+' and ((q_start < (flanking_len + 1 - flanking_region_distance) and s_start < (flanking_len + 1 - flanking_region_distance))
            #                       or (q_end > (orig_query_len + flanking_len + flanking_region_distance) and s_end > (orig_subject_len + flanking_len + flanking_region_distance))):
            #     deleted_subject_names.add(subject_name)
            # elif direct == '-' and ((q_start < (flanking_len + 1 - flanking_region_distance) and s_start > (orig_subject_len + flanking_len + flanking_region_distance))
            #                       or (q_end > (orig_query_len + flanking_len + flanking_region_distance) and s_end < (flanking_len + 1 - flanking_region_distance))):
            #     deleted_subject_names.add(subject_name)

    #存储query_groups中的比对情况，方便我们后期调试
    store_flank_align_groups(query_groups, flank_align_dir)


    delete_names = set()
    appeared_names = set()
    for query_name in query_groups.keys():
        # if query_name.__contains__('N_54503-len_5175-ref_NC_029264.1-11302200-11307375-C_27-tsd_TAA-distance_9'):
        #     print('here')

        # 统计整个flanking区域都能比对的次数，超过1次就过滤掉
        complete_false_num = 0
        # 统计明显超过边界的次数，超过5次就过滤掉
        obvious_false_num = 0
        appeared_names.add(query_name)
        subject_groups = query_groups[query_name]
        if len(subject_groups) == 0:
            delete_names.add(query_name)
        total_subject_names = set()
        deleted_subject_names = set()
        for subject_name in subject_groups.keys():
            total_subject_names.add(subject_name)
            group = subject_groups[subject_name]
            is_query_start_covered = False
            is_query_end_covered = False
            is_subject_start_covered = False
            is_subject_end_covered = False
            for item in group:
                q_start = item[0]
                q_end = item[1]
                orig_query_len = item[2]
                s_start = item[3]
                s_end = item[4]
                orig_subject_len = item[5]
                flanking_len = item[6]
                flanking_region_distance = item[7]
                direct = item[8]


                query_start_covered, query_end_covered, \
                subject_start_covered, subject_end_covered = overlap_with_boundary(q_start, q_end, s_start, s_end,
                                                                                         flanking_len,
                                                                                         flanking_region_distance,
                                                                                         orig_query_len, orig_subject_len)
                is_query_start_covered = is_query_start_covered or query_start_covered
                is_query_end_covered = is_query_end_covered or query_end_covered
                is_subject_start_covered = is_subject_start_covered or subject_start_covered
                is_subject_end_covered = is_subject_end_covered or subject_end_covered

                #如果边界位置向里向外各(20/40)bp的序列具有同源性，说明边界周围的序列具有同源性
                if direct == '+' and ((q_start < (flanking_len+1-20) and q_end > (flanking_len+1+20)) or
                                      (q_start < (orig_query_len+flanking_len-20) and q_end > (orig_query_len+flanking_len+20)) or
                                      (s_start < (flanking_len+1-20) and s_end > (flanking_len+1+20)) or
                                      (s_start < (orig_subject_len+flanking_len-20) and s_end > (orig_subject_len+flanking_len+20))):
                    obvious_false_num += 1
                elif direct == '-' and ((q_start < (flanking_len+1-20) and q_end > (flanking_len+1+20)) or
                                        (q_start < (orig_query_len+flanking_len-20) and q_end > (orig_query_len+flanking_len+20)) or
                                        (s_start > (flanking_len+1+20) and s_end < (flanking_len+1-20)) or
                                        (s_start > (orig_subject_len+flanking_len+20) and s_end < (orig_subject_len+flanking_len-20))):
                    obvious_false_num += 1

                if direct == '+' and ((q_start < (flanking_len+1-40) and q_end > (flanking_len+1+40)) or
                                      (q_start < (orig_query_len+flanking_len-40) and q_end > (orig_query_len+flanking_len+40)) or
                                      (s_start < (flanking_len+1-40) and s_end > (flanking_len+1+40)) or
                                      (s_start < (orig_subject_len+flanking_len-40) and s_end > (orig_subject_len+flanking_len+40))):
                    complete_false_num += 1
                elif direct == '-' and ((q_start < (flanking_len+1-40) and q_end > (flanking_len+1+40)) or
                                        (q_start < (orig_query_len+flanking_len-40) and q_end > (orig_query_len+flanking_len+40)) or
                                        (s_start > (flanking_len+1+40) and s_end < (flanking_len+1-40)) or
                                        (s_start > (orig_subject_len+flanking_len+40) and s_end < (orig_subject_len+flanking_len-40))):
                    complete_false_num += 1



                # 判断C0与C1,C2...等拷贝的比对情况，如果有边界外的区域能和拷贝高同源性，那么这条拷贝应该被抛弃，如果所有拷贝被抛弃，则该条序列应该是假阳性。
                if direct == '+' and \
                        (q_start < (flanking_len + 1 - flanking_region_distance) and q_end > (flanking_len + 1)) or \
                        (q_start < (orig_query_len + flanking_len) and q_end > (orig_query_len + flanking_len + flanking_region_distance)) or \
                        (s_start < (flanking_len + 1 - flanking_region_distance) and s_end > (flanking_len + 1)) or \
                        (s_start < (orig_subject_len + flanking_len) and s_end > (orig_subject_len + flanking_len + flanking_region_distance)):
                    deleted_subject_names.add(subject_name)
                    break
                elif direct == '-' and \
                        (q_start < (flanking_len + 1 - flanking_region_distance) and q_end > (flanking_len + 1)) or \
                        (q_start < (orig_query_len + flanking_len) and q_end > (orig_query_len + flanking_len + flanking_region_distance)) or \
                        (s_end < (flanking_len + 1 - flanking_region_distance) and s_start > (flanking_len + 1)) or \
                        (s_end < (orig_subject_len + flanking_len) and s_start > (orig_subject_len + flanking_len + flanking_region_distance)):
                    deleted_subject_names.add(subject_name)
                    break

            if not is_query_start_covered or not is_query_end_covered or\
                    not is_subject_start_covered or not is_subject_end_covered:
                deleted_subject_names.add(subject_name)

        # 如果50%拷贝的subject都有flanking区域能够比对，则这条序列应该被丢弃
        if complete_false_num >= 1 or obvious_false_num >= 5 or (len(total_subject_names) > 0 and float(len(deleted_subject_names))/len(total_subject_names) >= 0.5):
            delete_names.add(query_name)
    return delete_names, appeared_names

    # 如果有3条以上拷贝的subject都有flanking区域能够比对，则这条序列应该被丢弃
    # if len(total_subject_names) > 0 and len(deleted_subject_names) >= 3:
    #     return cur_orig_query_name
    # else:
    #     return ''

def flank_region_align_v1(candidate_sequence_path, flanking_len, similar_ratio, reference, TE_type, tmp_output_dir, threads, ref_index, log):
    log.logger.info('------generating candidate ' + TE_type + ' copies')
    starttime = time.time()
    tir_tsd_temp_dir = tmp_output_dir + '/' + TE_type + '_blast'
    all_copies = multi_process_align_and_get_copies(candidate_sequence_path, reference, tir_tsd_temp_dir, TE_type, threads)

    # multi_process_align(candidate_sequence_path, reference, blastnResults_path, blast_program_dir, tir_tsd_temp_dir, threads)
    # if not os.path.exists(blastnResults_path):
    #     return
    # all_copies = get_copies(blastnResults_path, candidate_sequence_path, reference, threads=threads)

    # 过滤掉拷贝数小于2, flanking copies
    ref_names, ref_contigs = read_fasta(reference)
    new_all_copies = {}
    for query_name in all_copies.keys():
        if TE_type == 'tir':
            # tsd = query_name.split('-tsd_')[1]
            tsd = ''
        else:
            tsd = ''
        copies = all_copies[query_name]
        if TE_type != 'ltr' and len(copies) < 2:
            continue
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
            if not new_all_copies.__contains__(query_name):
                new_all_copies[query_name] = []
            copy_list = new_all_copies[query_name]
            copy_list.append((ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq, tsd))

    # store confident_copies for testing
    confident_copies_file = tmp_output_dir + '/'+TE_type+'_copies.info'
    with open(confident_copies_file, 'w') as f_save:
        f_save.write('# all copies have been flanked ' + str(flanking_len) +' bp at both ends\n')
        for orig_query_name in new_all_copies.keys():
            f_save.write(orig_query_name + '\n')
            for copy in new_all_copies[orig_query_name]:
                f_save.write('\tfrom:' + str(copy[0]) + '_' + str(copy[1]) + '_' + str(copy[2]) + '_' + str(
                    copy[2] - copy[1] + 1) + '\n')
                f_save.write(copy[4] + '\n')
    f_save.close()

    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of generating candidate " + TE_type + " copies: %.8s s" % (dtime))

    log.logger.info('------flanking ' + TE_type + ' alignment')
    starttime = time.time()
    if ref_index == -1:
        flank_align_dir = tmp_output_dir + '/flank_' + TE_type + '_align'
    else:
        flank_align_dir = tmp_output_dir + '/flank_'+TE_type+'_align_' + str(ref_index)
    os.system('rm -rf '+flank_align_dir)
    if not os.path.exists(flank_align_dir):
        os.makedirs(flank_align_dir)

    # 保留单拷贝的LTR转座子
    single_copy_LTR = set()

    flanking_region_distance = int(flanking_len * similar_ratio)
    #flanking_region_distance = 10
    ex = ProcessPoolExecutor(threads)
    jobs = []
    #一条orig_query_name调用一次比对，CPU利用率太低，我们尝试50条序列调用一次
    batch_size = 50
    cur_query_copies = {}
    cur_subject_copies = {}
    batch_num = 0
    for index, orig_query_name in enumerate(new_all_copies.keys()):
        if index % batch_size == 0 and len(cur_query_copies) > 0:
            cur_query_copy_path = flank_align_dir + '/' + str(batch_num) + '_query.fa'
            store_fasta(cur_query_copies, cur_query_copy_path)
            cur_subject_copy_path = flank_align_dir + '/' + str(batch_num) + '_subject.fa'
            store_fasta(cur_subject_copies, cur_subject_copy_path)
            output_path = flank_align_dir + '/' + str(batch_num) + '.out'
            job = ex.submit(run_blast_align, cur_query_copy_path, cur_subject_copy_path, output_path,
                            flanking_len, flanking_region_distance, flank_align_dir)
            jobs.append(job)
            batch_num += 1
            cur_query_copies = {}
            cur_subject_copies = {}

        copy_list = new_all_copies[orig_query_name]

        if TE_type == 'ltr' and len(copy_list) <= 1:
            single_copy_LTR.add(orig_query_name)
            continue

        for i, copy in enumerate(copy_list):
            # copyt -> (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq, tsd)
            ref_name = copy[0]
            ref_start = int(copy[1])
            ref_end = int(copy[2])
            ref_flank_seq = copy[4]
            if len(ref_flank_seq) < 100:
                continue
            if i == 0:
                cur_query_copies[orig_query_name + '-c_' + str(i) + '-' + str(ref_end - ref_start + 1)] = ref_flank_seq
            else:
                cur_subject_copies[orig_query_name + '-c_' + str(i) + '-' + str(ref_end - ref_start + 1)] = ref_flank_seq

    if len(cur_query_copies) > 0:
        cur_query_copy_path = flank_align_dir + '/' + str(batch_num) + '_query.fa'
        store_fasta(cur_query_copies, cur_query_copy_path)
        cur_subject_copy_path = flank_align_dir + '/' + str(batch_num) + '_subject.fa'
        store_fasta(cur_subject_copies, cur_subject_copy_path)
        output_path = flank_align_dir + '/' + str(batch_num) + '.out'
        job = ex.submit(run_blast_align, cur_query_copy_path, cur_subject_copy_path, output_path,
                        flanking_len, flanking_region_distance, flank_align_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    deleted_names = set()
    appeared_names = set()
    for job in as_completed(jobs):
        cur_delete_names, cur_appeared_names = job.result()
        deleted_names.update(cur_delete_names)
        appeared_names.update(cur_appeared_names)

    for name in new_all_copies.keys():
        if name not in appeared_names:
            deleted_names.add(name)

    for cur_delete_name in deleted_names:
        if cur_delete_name != '' and cur_delete_name not in single_copy_LTR:
            del new_all_copies[cur_delete_name]

    print('deleted_names len: ' + str(len(deleted_names)))
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of flanking " + TE_type + " alignment: %.8s s" % (dtime))
    return new_all_copies

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

    #这里是不是可以考虑把序列划分的更细，以此来减少任务的不均衡
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

    # longest_repeat_files = []
    # segments_cluster = divided_array(list(orig_contigs.items()), threads)
    # for partition_index, cur_segments in enumerate(segments_cluster):
    #     if len(cur_segments) <= 0:
    #         continue
    #     single_tmp_dir = tmp_blast_dir + '/' + str(partition_index)
    #     if not os.path.exists(single_tmp_dir):
    #         os.makedirs(single_tmp_dir)
    #     split_repeat_file = single_tmp_dir + '/repeats_split.fa'
    #     cur_contigs = {}
    #     for item in cur_segments:
    #         cur_contigs[item[0]] = item[1]
    #     store_fasta(cur_contigs, split_repeat_file)
    #     repeats_path = (split_repeat_file, subject_path, single_tmp_dir + '/temp.out')
    #     longest_repeat_files.append(repeats_path)

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for file in longest_repeat_files:
        job = ex.submit(multiple_alignment_blast_and_get_copies, file, tools_dir, query_coverage, TE_type, subject_coverage=subject_coverage)
        jobs.append(job)
    ex.shutdown(wait=True)

    all_copies = {}
    for job in as_completed(jobs):
        cur_copies = job.result()
        all_copies.update(cur_copies)

    return all_copies


def multi_process_align(query_path, subject_path, blastnResults_path, tmp_blast_dir, threads, is_removed_dir=True):
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

def remove_ltr_from_tir(confident_ltr_cut_path, confident_tir_path, all_copies):
    # copy -> (subject_name, subject_start, subject_end, query[2], direct)
    query_names, query_contigs = read_fasta(confident_ltr_cut_path)
    subject_names, subject_contigs = read_fasta(confident_tir_path)
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        for copy in copies:
            subject_name = copy[0]
            if subject_contigs.__contains__(subject_name):
                del subject_contigs[subject_name]
            # subject_len = len(subject_contigs[subject_name])
            # alignment_len = copy[2] - copy[1] + 1
            # # 这是一条LTR序列，应该被过滤掉
            # if float(alignment_len) / subject_len >= 0.95:
            #     del query_contigs[query_name]
            #     break
    store_fasta(subject_contigs, confident_tir_path)


if __name__ == '__main__':
    trf_data_path = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/ReferenceMode/output/CRD.2022-05-04.9-30-26/trf/dmel-all-chromosome-r5.43.fasta.2.7.7.80.10.50.500.dat'
    tandem_elements = extract_tandem_from_trf(trf_data_path)
    tandem_path = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/ReferenceMode/output/CRD.2022-05-04.9-30-26/trf/tandem.fa'

    with open(tandem_path, 'w') as f_save:
        for index, elem in enumerate(tandem_elements):
            f_save.write('>Node_'+str(index)+'\n'+elem+'\n')
    f_save.close()