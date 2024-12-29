import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
import logging
from logging import handlers
from fuzzysearch import find_near_matches
import Levenshtein
from concurrent.futures import ProcessPoolExecutor, as_completed

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


def save_dict_to_fasta(data_dict, fasta_file):
    with open(fasta_file, 'w') as file:
        for identifier, sequence in data_dict.items():
            # 写入FASTA头部
            file.write(f">{identifier}\n")

            # 分割序列以确保每行不超过70个字符
            seq_lines = [sequence[i:i + 70] for i in range(0, len(sequence), 70)]

            # 写入序列数据
            for line in seq_lines:
                file.write(line + "\n")

def map_chr_position(LTR_detector_scn_file, scn_file, position_map):
    ltr_candidates, ltr_lines = read_scn(LTR_detector_scn_file, log=None)
    with open(scn_file, 'w') as f_save:
        f_save.write('# LtrDetector\n')
        for candidate_index in ltr_candidates.keys():
            line = ltr_lines[candidate_index]
            parts = line.split(' ')
            ltr_start = int(parts[0])
            ltr_end = int(parts[1])
            chr_name = parts[11]
            left_ltr_start = int(parts[3])
            left_ltr_end = int(parts[4])
            right_ltr_start = int(parts[6])
            right_ltr_end = int(parts[7])
            raw_chr_name, start, end = position_map[chr_name]
            parts[0] = str(start + ltr_start)
            parts[1] = str(start + ltr_end)
            parts[11] = raw_chr_name
            parts[3] = str(start + left_ltr_start)
            parts[4] = str(start + left_ltr_end)
            parts[6] = str(start + right_ltr_start)
            parts[7] = str(start + right_ltr_end)
            new_line = ' '.join(parts)
            f_save.write(new_line + '\n')



def split_chromosomes(chromosomes_dict, max_length=400_000_000):
    """
    分割染色体序列，如果序列长度超过max_length，则将其分割成多个部分。

    参数:
    chromosomes_dict (dict): 一个字典，键为染色体名称，值为对应的DNA序列。
    max_length (int): 最大序列长度，超过此长度的序列将被分割。默认值为200 MB。

    返回:
    dict: 一个新的字典，包含分割后的染色体序列。
    """
    new_chromosomes_dict = {}
    position_map = {}
    for chrom, sequence in chromosomes_dict.items():
        if len(sequence) > max_length:
            num_parts = (len(sequence) + max_length - 1) // max_length  # 计算需要分割的部分数
            for i in range(num_parts):
                part_name = f"{chrom}_part{i + 1}"
                start = i * max_length
                end = min((i + 1) * max_length, len(sequence))
                new_chromosomes_dict[part_name] = sequence[start:end]
                position_map[part_name] = (chrom, start, end)
        else:
            new_chromosomes_dict[chrom] = sequence
            position_map[chrom] = (chrom, 0, len(sequence))
    return new_chromosomes_dict, position_map

def split_dict_into_blocks(chromosomes_dict, threads):
    # chromosomes_dict = split_chromosomes(chromosomes_dict, max_length=chunk_size)
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

def read_scn(scn_file, log, remove_dup=False):
    ltr_candidates = {}
    ltr_lines = {}
    candidate_index = 0
    existing_records = set()
    remove_count = 0
    total_lines = 0
    with open(scn_file, 'r') as f_r:
        for line in f_r:
            if line.startswith('#') or line.strip() == '':
                continue
            line = line.replace('\n', '')
            parts = line.split(' ')
            ltr_start = int(parts[0])
            ltr_end = int(parts[1])
            chr_name = parts[11]
            total_lines += 1
            if remove_dup:
                cur_record = (ltr_start, ltr_end, chr_name)
                # 过滤掉冗余的记录
                if cur_record in existing_records:
                    remove_count += 1
                    continue
                existing_records.add(cur_record)
            left_ltr_start = int(parts[3])
            left_ltr_end = int(parts[4])
            right_ltr_start = int(parts[6])
            right_ltr_end = int(parts[7])

            ltr_candidates[candidate_index] = (chr_name, left_ltr_start, left_ltr_end, right_ltr_start, right_ltr_end)
            ltr_lines[candidate_index] = line
            candidate_index += 1
    if remove_dup and log is not None:
        log.logger.debug('Total LTR num: ' + str(total_lines))
        log.logger.debug('Remove ' + str(remove_count) + ' replicate LTR in scn file, remaining LTR num: ' + str(len(ltr_candidates)))
    return ltr_candidates, ltr_lines

def store_scn(confident_lines, confident_scn):
    with open(confident_scn, 'w') as f_save:
        for line in confident_lines:
            f_save.write(line + '\n')

def rename_fasta(input, output, header='N'):
    names, contigs = read_fasta(input)
    node_index = 0
    with open(output, 'w') as f_save:
        for name in names:
            seq = contigs[name]
            f_save.write('>'+header+'_'+str(node_index)+'\n'+seq+'\n')
            node_index += 1
    f_save.close()

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

def store_fasta(contigs, file_path):
    with open(file_path, 'w') as f_save:
        for name in contigs.keys():
            seq = contigs[name]
            f_save.write('>'+name+'\n'+seq+'\n')
    f_save.close()

def get_full_length_copies_batch(confident_ltr_internal, split_ref_dir, threads, temp_dir):
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    batch_size = 10
    batch_id = 0
    names, contigs = read_fasta(confident_ltr_internal)
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

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cur_split_files in split_files:
        job = ex.submit(get_full_length_copies, cur_split_files, split_ref_dir, debug=0)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_copies = {}
    for job in as_completed(jobs):
        cur_all_copies = job.result()
        all_copies.update(cur_all_copies)
    return all_copies

def get_full_length_copies_batch_v1(confident_ltr_internal, split_ref_dir, threads, temp_dir, max_copy_num, full_length_threshold):
    os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    batch_size = 10
    batch_id = 0
    names, contigs = read_fasta(confident_ltr_internal)
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

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cur_split_files in split_files:
        job = ex.submit(get_full_length_copies_v1, cur_split_files, split_ref_dir, max_copy_num, full_length_threshold, debug=1)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_copies = {}
    for job in as_completed(jobs):
        cur_all_copies = job.result()
        all_copies.update(cur_all_copies)
    return all_copies

def get_full_length_copies_v1(query_path, split_ref_dir, max_copy_num, full_length_threshold, debug):
    blastn2Results_path = query_path + '.blast.out'
    repeats_path = (query_path, split_ref_dir, blastn2Results_path)
    all_copies = multiple_alignment_blast_and_get_copies_v2(repeats_path, max_copy_num, full_length_threshold)
    if debug != 1:
        os.remove(blastn2Results_path)
    return all_copies

def multiple_alignment_blast_and_get_copies_v2(repeats_path, max_copy_num, full_length_threshold):
    split_repeats_path = repeats_path[0]
    split_ref_dir = repeats_path[1]
    # raw_blastn2Results_path = repeats_path[2]
    # os.system('rm -f ' + raw_blastn2Results_path)
    blastn2Results_path = repeats_path[2]
    os.system('rm -f ' + blastn2Results_path)
    all_copies = {}
    repeat_names, repeat_contigs = read_fasta(split_repeats_path)
    remain_contigs = repeat_contigs
    for chr_name in os.listdir(split_ref_dir):
        # blastn2Results_path = raw_blastn2Results_path + '_' + str(chr_name)
        if len(remain_contigs) > 0:
            if not str(chr_name).endswith('.fa'):
                continue
            chr_path = split_ref_dir + '/' + chr_name
            align_command = 'blastn -db ' + chr_path + ' -num_threads ' \
                            + str(1) + ' -query ' + split_repeats_path + ' -evalue 1e-20 -outfmt 6 > ' + blastn2Results_path
            os.system(align_command)
            # 由于我们只需要100个拷贝，因此如果有序列已经满足了，就不需要进行后续的比对了，这样在mouse这样的高拷贝大基因组上减少运行时间
            cur_all_copies = get_copies_v2(blastn2Results_path, split_repeats_path, max_copy_num, full_length_threshold)
            for query_name in cur_all_copies.keys():
                copy_list = cur_all_copies[query_name]
                if query_name in all_copies:
                    prev_copy_list = all_copies[query_name]
                else:
                    prev_copy_list = []
                update_copy_list = prev_copy_list + copy_list
                all_copies[query_name] = update_copy_list
                if len(update_copy_list) >= max_copy_num:
                    del repeat_contigs[query_name]
            remain_contigs = repeat_contigs
            store_fasta(remain_contigs, split_repeats_path)

    # all_copies = get_copies_v1(blastn2Results_path, split_repeats_path, '')
    return all_copies

def get_copies_v2(blastnResults_path, query_path, max_copy_num, full_length_threshold=0.95):
    query_names, query_contigs = read_fasta(query_path)
    longest_repeats = FMEA_new(query_path, blastnResults_path, full_length_threshold)
    # (query_name, old_query_start_pos, old_query_end_pos, subject_name, old_subject_start_pos, old_subject_end_pos)
    all_copies = {}
    for query_name in longest_repeats.keys():
        longest_queries = longest_repeats[query_name]
        longest_queries.sort(key=lambda x: -x[2])
        query_len = len(query_contigs[query_name])
        copies = []
        keeped_copies = set()
        for query in longest_queries:
            if len(copies) > max_copy_num:
                break
            subject_name = query[3]
            subject_start = query[4]
            subject_end = query[5]
            cur_query_len = abs(query[2] - query[1])
            direct = '+'
            if subject_start > subject_end:
                tmp = subject_start
                subject_start = subject_end
                subject_end = tmp
                direct = '-'
            item = (subject_name, subject_start, subject_end)

            if float(cur_query_len) / query_len >= full_length_threshold and item not in keeped_copies:
                copies.append((subject_name, subject_start, subject_end, cur_query_len, direct))
                keeped_copies.add(item)
        # copies.sort(key=lambda x: abs(x[3]-(x[2]-x[1]+1)))
        all_copies[query_name] = copies

    return all_copies

def flank_region_align_v5(candidate_sequence_path, flanking_len, reference, split_ref_dir, TE_type, tmp_output_dir, threads, ref_index, log, subset_script_path, plant, debug, iter_num, result_type='cons'):
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

    ref_names, ref_contigs = read_fasta(reference)
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
        cur_member_file = temp_dir + '/' + query_name + '.blast.bed.fa'
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

    true_tes = {}
    for job in as_completed(jobs):
        cur_name, is_TE = job.result()
        true_tes[cur_name] = is_TE

    if debug != 1:
        os.system('rm -rf ' + temp_dir)
    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of determination of homology in regions outside the boundaries of  " + TE_type + " copies: %.8s s" % (dtime))
    return true_tes

def get_full_length_copies(query_path, split_ref_dir, debug):
    blastn2Results_path = query_path + '.blast.out'
    repeats_path = (query_path, split_ref_dir, blastn2Results_path)
    all_copies = multiple_alignment_blast_and_get_copies_v1(repeats_path)
    if debug != 1:
        os.remove(blastn2Results_path)
    return all_copies

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

def get_query_copies(cur_segments, query_contigs, subject_path, query_coverage, subject_coverage, query_fixed_extend_base_threshold=1000, subject_fixed_extend_base_threshold=1000, max_copy_num=100):
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


def generate_left_frame_from_seq(candidate_sequence_path, reference, threads, temp_dir, output_dir, split_ref_dir):
    debug = 0
    flanking_len = 50
    starttime = time.time()
    if os.path.exists(temp_dir):
        os.system('rm -rf ' + temp_dir)
    os.makedirs(temp_dir)
    if os.path.exists(output_dir):
        os.system('rm -rf ' + output_dir)
    os.makedirs(output_dir)

    # We are considering that the current running time is too long, maybe it is related to submitting one sequence for Blastn alignment at a time.
    # We will try to combine 10 sequences together and run Blastn once.
    # To increase CPU utilization, we will submit one thread to process 10 sequences.
    batch_size = 10
    batch_id = 0
    names, contigs = read_fasta(candidate_sequence_path)
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

    ref_names, ref_contigs = read_fasta(reference)
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
        cur_member_file = temp_dir + '/' + query_name + '.blast.bed.fa'
        store_fasta(copy_contigs, cur_member_file)
        query_seq = contigs[query_name]
        batch_member_files.append((query_name, query_seq, cur_member_file))


    subset_script_path = os.getcwd() + '/tools/ready_for_MSA.sh'
    # Determine whether the multiple sequence alignment of each copied file satisfies the homology rule
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for batch_member_file in batch_member_files:
        job = ex.submit(generate_msa, batch_member_file, temp_dir, output_dir, subset_script_path, debug)
        jobs.append(job)
    ex.shutdown(wait=True)

    all_left_frames = []
    for job in as_completed(jobs):
        left_frame_path = job.result()
        all_left_frames.append(left_frame_path)

    endtime = time.time()
    dtime = endtime - starttime
    return all_left_frames

def generate_both_ends_frame_from_seq(candidate_sequence_path, reference, flanking_len,
                                      threads, temp_dir, output_dir, full_length_output_dir, split_ref_dir,
                                      max_copy_num, full_length_threshold):
    debug = 0
    starttime = time.time()
    if os.path.exists(temp_dir):
        os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    if os.path.exists(output_dir):
        os.system('rm -rf ' + output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if os.path.exists(full_length_output_dir):
        os.system('rm -rf ' + full_length_output_dir)
    if not os.path.exists(full_length_output_dir):
        os.makedirs(full_length_output_dir)

    # We are considering that the current running time is too long, maybe it is related to submitting one sequence for Blastn alignment at a time.
    # We will try to combine 10 sequences together and run Blastn once.
    # To increase CPU utilization, we will submit one thread to process 10 sequences.
    batch_size = 10
    batch_id = 0
    names, contigs = read_fasta(candidate_sequence_path)
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

    ref_names, ref_contigs = read_fasta(reference)
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cur_split_files in split_files:
        job = ex.submit(get_full_length_copies_v1, cur_split_files, split_ref_dir, max_copy_num, full_length_threshold, debug)
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
        cur_member_file = temp_dir + '/' + query_name + '.blast.bed.fa'
        store_fasta(copy_contigs, cur_member_file)
        batch_member_files.append((query_name, cur_member_file))

    endtime = time.time()
    dtime = endtime - starttime
    print("Running time of get copies: %.8s s" % (dtime))

    starttime = time.time()
    # subset_script_path = config.project_dir + '/tools/ready_for_MSA.sh'
    # Determine whether the multiple sequence alignment of each copied file satisfies the homology rule
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for batch_member_file in batch_member_files:
        job = ex.submit(generate_msa, batch_member_file, temp_dir, output_dir, full_length_output_dir, flanking_len, debug)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        both_end_frame_paths = job.result()

    if not debug:
        os.system('rm -rf ' + temp_dir)

    endtime = time.time()
    dtime = endtime - starttime
    print("Running time of MSA: %.8s s" % (dtime))


def generate_both_ends_frame_for_intactLTR(candidate_sequence_path, reference, flanking_len, threads, temp_dir,
                                           output_dir, full_length_output_dir, split_ref_dir, log):
    debug = 0
    # flanking_len = 100
    starttime = time.time()
    if os.path.exists(temp_dir):
        os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    if os.path.exists(output_dir):
        os.system('rm -rf ' + output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if os.path.exists(full_length_output_dir):
        os.system('rm -rf ' + full_length_output_dir)
    if not os.path.exists(full_length_output_dir):
        os.makedirs(full_length_output_dir)

    # We are considering that the current running time is too long, maybe it is related to submitting one sequence for Blastn alignment at a time.
    # We will try to combine 10 sequences together and run Blastn once.
    # To increase CPU utilization, we will submit one thread to process 10 sequences.
    batch_size = 10
    batch_id = 0
    names, contigs = read_fasta(candidate_sequence_path)
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

    ref_names, ref_contigs = read_fasta(reference)
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

    # 判断LTR序列是否在整个基因组上存在至少一个以上的全长拷贝
    single_copy_names = []
    filtered_intact_count = 0
    for name in all_copies.keys():
        if len(all_copies[name]) < 2:
            single_copy_names.append(name)
            filtered_intact_count += 1
    if log is not None:
        log.logger.info('Filter the number of intact LTR <= 1 full-length copy: ' + str(filtered_intact_count))

    # extend copies
    batch_member_files = []
    new_all_copies = {}
    for query_name in all_copies.keys():
        if query_name in single_copy_names:
            continue
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
            # 由于我们只是判断它的侧翼区域是否同源，因此我们不需要太长的序列
            left_end = 2 * flanking_len if 2 * flanking_len < len(copy_seq) else len(copy_seq)
            copy_seq = copy_seq[0: left_end] + copy_seq[-left_end: ]

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

    # subset_script_path = config.project_dir + '/tools/ready_for_MSA.sh'
    # Determine whether the multiple sequence alignment of each copied file satisfies the homology rule
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for batch_member_file in batch_member_files:
        job = ex.submit(generate_msa, batch_member_file, temp_dir, output_dir, full_length_output_dir, flanking_len, debug)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        left_frame_path, full_length_align_file = job.result()

    endtime = time.time()
    dtime = endtime - starttime

def merge_terminals(confident_ltr_terminal_cons, threads):
    blastn_command = 'blastn -query ' + confident_ltr_terminal_cons + ' -subject ' + confident_ltr_terminal_cons + ' -num_threads ' + str(threads) + ' -outfmt 6 '
    contignames, contigs = read_fasta(confident_ltr_terminal_cons)

    result = subprocess.run(blastn_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            executable='/bin/bash')

    remove_frag_ltr = set()
    duplicate_records = set()
    if result.returncode == 0:
        lines = result.stdout.split('\n')
        for line in lines:
            parts = line.split('\t')
            if len(parts) != 12:
                continue
            query_name = parts[0]
            subject_name = parts[1]
            if query_name == subject_name:
                continue
            query_len = len(contigs[query_name])
            subject_len = len(contigs[subject_name])

            query_start = int(parts[6])
            query_end = int(parts[7])
            if abs(query_end - query_start) / query_len >= 0.95:
                if (query_name, subject_name) not in duplicate_records and (subject_name, query_name) not in duplicate_records:
                    if query_len < subject_len:
                        remove_frag_ltr.add(query_name)
                    else:
                        remove_frag_ltr.add(subject_name)
                    duplicate_records.add((query_name, subject_name))
                    duplicate_records.add((subject_name, query_name))

    print('remove fragmented LTR num: ' + str(len(remove_frag_ltr)))

    for ltr_name in remove_frag_ltr:
        del contigs[ltr_name]
    store_fasta(contigs, confident_ltr_terminal_cons)



def extract_copies(member_file, max_num):
    contignames, contigs = read_fasta(member_file)
    new_contigs = {}
    for name in contignames[:max_num]:
        new_contigs[name] = contigs[name]
    # sorted_contigs = dict(sorted(contigs.items(), key=lambda item: len(item[1]))[:max_num])
    member_file += '.rdmSubset.fa'
    store_fasta(new_contigs, member_file)
    return member_file

def generate_msa(batch_member_file, temp_dir, output_dir, full_length_output_dir, flanking_len, debug):
    (query_name, member_file) = batch_member_file

    member_names, member_contigs = read_fasta(member_file)
    if len(member_names) > 100:
        # 抽取100条最长的 拷贝
        max_num = 100
        member_file = extract_copies(member_file, max_num)
    if not os.path.exists(member_file):
        return None, None
    align_file = member_file + '.maf.fa'
    align_command = 'cd ' + temp_dir + ' && mafft --preservecase --quiet --thread 1 ' + member_file + ' > ' + align_file
    os.system(align_command)
    # left_frame_path = get_left_frame(query_name, cur_seq, align_file, output_dir, debug)
    if len(member_names) >= 1:
        cur_seq = member_contigs[member_names[0]][flanking_len:-flanking_len]
        both_end_frame_path, full_length_align_file  = get_both_ends_frame(query_name, cur_seq, align_file, output_dir, full_length_output_dir, flanking_len, debug)
    else:
        both_end_frame_path = ''
    return both_end_frame_path


def remove_sparse_col_in_align_file(align_file, align_start, align_end):
    align_names, align_contigs = read_fasta(align_file)
    first_seq = align_contigs[align_names[0]]
    col_num = len(first_seq)
    row_num = len(align_names)

    matrix = [[''] * col_num for i in range(row_num)]
    for row, name in enumerate(align_names):
        seq = align_contigs[name]
        for col in range(len(seq)):
            matrix[row][col] = seq[col]

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

    new_align_start = -1
    new_align_end = -1
    valid_col = 0
    valid_col_threshold = row_num / 2
    sparse_cols = []
    for col_index in range(col_num):
        base_map = col_base_map[col_index]
        gap_num = base_map['-']
        if col_index == align_start:
            new_align_start = valid_col
        elif col_index == align_end:
            new_align_end = valid_col
        # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
        elif gap_num > valid_col_threshold:
            sparse_cols.append(col_index)
            continue
        valid_col += 1

    clean_align_file = align_file + '.clean.fa'
    with open(clean_align_file, 'w') as f_save:
        for name in align_names:
            seq = align_contigs[name]
            new_seq = ''
            for i in range(len(seq)):
                if i in sparse_cols:
                    continue
                else:
                    new_seq += seq[i]
            f_save.write('>' + name + '\n' + new_seq + '\n')
    return clean_align_file, new_align_start, new_align_end

def get_both_ends_frame(query_name, cur_seq, align_file, output_dir, full_length_output_dir, flanking_len, debug):
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
        return None, None
    # print(align_start, align_end)
    # 清理 align file 中的稀疏列，并生成新的 align_start, align_end 位置
    align_file, align_start, align_end = remove_sparse_col_in_align_file(align_file, align_start, align_end)
    # print(align_start, align_end)
    align_names, align_contigs = read_fasta(align_file)

    # 取 align_start 的外侧 100 bp
    out_threshold = flanking_len
    start_align_file = output_dir + '/' + query_name + '.matrix'
    with open(start_align_file, 'w') as f_save:
        for name in align_names:
            raw_align_seq = align_contigs[name]
            # 取左侧 100 bp frame
            if align_start - out_threshold < 0:
                start_pos = 0
            else:
                start_pos = align_start - out_threshold
            start_seq = raw_align_seq[start_pos: align_start]
            # 补全不足 100 bp的部分
            seq1 = '-' * (out_threshold - len(start_seq)) + start_seq

            # 取右侧 100 bp frame
            if align_end + out_threshold > len(raw_align_seq):
                end_pos = len(raw_align_seq)
            else:
                end_pos = align_end + out_threshold
            end_seq = raw_align_seq[align_end: end_pos]
            # 补全不足 100 bp的部分
            seq2 = end_seq + '-' * (out_threshold - len(end_seq))

            f_save.write(seq1+'\t'+seq2+'\n')

    full_length_align_file = full_length_output_dir + '/' + query_name + '.matrix'
    with open(full_length_align_file, 'w') as f_save:
        for name in align_names:
            raw_align_seq = align_contigs[name]
            # 取左侧 100 bp frame
            if align_start - out_threshold < 0:
                start_pos = 0
            else:
                start_pos = align_start - out_threshold
            start_seq = raw_align_seq[start_pos: align_start]
            # 补全不足 100 bp的部分
            seq1 = '-' * (out_threshold - len(start_seq)) + start_seq

            # 取中间的序列
            middle_seq = raw_align_seq[align_start: align_end]

            # 取右侧侧 100 bp frame
            if align_end + out_threshold > len(raw_align_seq):
                end_pos = len(raw_align_seq)
            else:
                end_pos = align_end + out_threshold
            end_seq = raw_align_seq[align_end: end_pos]
            # 补全不足 100 bp的部分
            seq2 = end_seq + '-' * (out_threshold - len(end_seq))

            f_save.write(seq1 + middle_seq + seq2 + '\n')

    return start_align_file, full_length_align_file


def get_both_ends_frame_v2(query_name, cur_seq, align_file, full_length_output_dir, flanking_len, debug):
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
        return None, None
    # print(align_start, align_end)
    # 清理 align file 中的稀疏列，并生成新的 align_start, align_end 位置
    align_file, align_start, align_end = remove_sparse_col_in_align_file(align_file, align_start, align_end)
    # print(align_start, align_end)
    align_names, align_contigs = read_fasta(align_file)

    # 取 align_start 的外侧 100 bp
    out_threshold = flanking_len
    full_length_align_file = full_length_output_dir + '/' + query_name + '.matrix'
    with open(full_length_align_file, 'w') as f_save:
        for name in align_names:
            raw_align_seq = align_contigs[name]
            # 取左侧 100 bp frame
            if align_start - out_threshold < 0:
                start_pos = 0
            else:
                start_pos = align_start - out_threshold
            start_seq = raw_align_seq[start_pos: align_start]
            # 补全不足 100 bp的部分
            seq1 = '-' * (out_threshold - len(start_seq)) + start_seq

            # 取中间的序列
            middle_seq = raw_align_seq[align_start: align_end]

            # 取右侧侧 100 bp frame
            if align_end + out_threshold > len(raw_align_seq):
                end_pos = len(raw_align_seq)
            else:
                end_pos = align_end + out_threshold
            end_seq = raw_align_seq[align_end: end_pos]
            # 补全不足 100 bp的部分
            seq2 = end_seq + '-' * (out_threshold - len(end_seq))

            f_save.write(seq1 + middle_seq + seq2 + '\n')

    return full_length_align_file

def get_both_ends_frame_v1(query_name, cur_seq, align_file, output_dir, flanking_len, debug):
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
        return None, None
    # print(align_start, align_end)
    # 清理 align file 中的稀疏列，并生成新的 align_start, align_end 位置
    align_file, align_start, align_end = remove_sparse_col_in_align_file(align_file, align_start, align_end)
    # print(align_start, align_end)
    align_names, align_contigs = read_fasta(align_file)

    # 取 align_start 的外侧 100 bp
    out_threshold = flanking_len
    start_align_file = output_dir + '/' + query_name + '.matrix'
    with open(start_align_file, 'w') as f_save:
        for name in align_names:
            raw_align_seq = align_contigs[name]
            # 取左侧 100 bp frame
            if align_start - out_threshold < 0:
                start_pos = 0
            else:
                start_pos = align_start - out_threshold
            start_seq = raw_align_seq[start_pos: align_start]
            # 补全不足 100 bp的部分
            seq1 = '-' * (out_threshold - len(start_seq)) + start_seq

            # 取右侧 100 bp frame
            if align_end + out_threshold > len(raw_align_seq):
                end_pos = len(raw_align_seq)
            else:
                end_pos = align_end + out_threshold
            end_seq = raw_align_seq[align_end: end_pos]
            # 补全不足 100 bp的部分
            seq2 = end_seq + '-' * (out_threshold - len(end_seq))

            f_save.write(seq1+'\t'+seq2+'\n')

    return start_align_file

def get_left_frame(query_name, cur_seq, align_file, output_dir, debug):
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
        return False

    align_names, align_contigs = read_fasta(align_file)
    if len(align_names) <= 0:
        if debug:
            print('align file size = 0, ' + align_file)
        return False

    # 取 align_start 的外侧 100 bp
    out_threshold = 100
    start_align_file = output_dir + '/' + query_name + '.matrix'
    with open(start_align_file, 'w') as f_save:
        for name in align_names:
            raw_align_seq = align_contigs[name]
            if align_start - out_threshold < 0:
                start_pos = 0
            else:
                start_pos = align_start - out_threshold
            start_seq = raw_align_seq[start_pos: align_start]
            # 补全不足 100 bp的部分
            seq = '-' * (out_threshold - len(start_seq)) + start_seq
            f_save.write(seq+'\n')
    return start_align_file





def run_find_members_v8(batch_member_file, temp_dir, subset_script_path, plant, TE_type, debug, result_type):
    (query_name, cur_seq, member_file) = batch_member_file

    member_names, member_contigs = read_fasta(member_file)
    if len(member_names) > 100:
        sub_command = 'cd ' + temp_dir + ' && sh ' + subset_script_path + ' ' + member_file + ' 100 100 ' + ' > /dev/null 2>&1'
        os.system(sub_command)
        member_file += '.rdmSubset.fa'
    if not os.path.exists(member_file):
        return (query_name, False)
    align_file = member_file + '.maf.fa'
    align_command = 'cd ' + temp_dir + ' && mafft --preservecase --quiet --thread 1 ' + member_file + ' > ' + align_file
    os.system(align_command)

    is_TE = judge_boundary_v9(cur_seq, align_file, debug, TE_type, plant, result_type)
    return query_name, is_TE


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

def TSDsearch_v5(raw_align_seq, cur_boundary_start, cur_boundary_end):
    # 2->TA or animal/fungi 5'-CCC...GGG-3', 3-> plant 5'-CACT(A/G)...(C/T)AGTG-3' or （TAA或TTA）, 4-> TTAA,
    TIR_TSDs = [11, 10, 9, 8, 7, 6, 5, 4, 3, 2]

    left_tsd_seq = ''
    right_tsd_seq = ''
    allow_mismatch_num = 1
    for tsd_len in TIR_TSDs:
        left_tsd = get_boundary_ungap_str(raw_align_seq, cur_boundary_start-1, tsd_len, 'left')
        right_tsd = get_boundary_ungap_str(raw_align_seq, cur_boundary_end+1, tsd_len, 'right')
        if len(left_tsd) != len(right_tsd) or len(left_tsd) != tsd_len:
            continue
        if left_tsd == right_tsd:
            left_tsd_seq = left_tsd
            right_tsd_seq = right_tsd
        elif tsd_len >= 8 and allow_mismatch(left_tsd, right_tsd, allow_mismatch_num):
            left_tsd_seq = left_tsd
            right_tsd_seq = right_tsd
    return left_tsd_seq, right_tsd_seq


def read_matrix_file(matrix_file):
    align_names = []
    align_contigs = {}
    index_num = 0
    with open(matrix_file, 'r') as f_r:
        for line in f_r:
            name = 'seq_' + str(index_num)
            align_names.append(name)
            seq = line.replace('\n', '')
            align_contigs[name] = seq
            index_num += 1
    return align_names, align_contigs


def judge_boundary_v5_batch(job_list, full_length_output_dir, TE_type, flanking_len):
    results = []
    for ltr_name in job_list:
        cur_matrix_file = full_length_output_dir + '/' + ltr_name + '.matrix'
        if not os.path.exists(cur_matrix_file):
            continue
        ltr_name, is_TE, info, final_cons_seq = judge_boundary_v5(cur_matrix_file, ltr_name, TE_type, flanking_len)
        results.append((ltr_name, is_TE, info, final_cons_seq))
    return results

def judge_boundary_v5(matrix_file, ltr_name, TE_type, flanking_len):
    align_names, align_contigs = read_matrix_file(matrix_file)

    debug = 0
    col_num = -1
    row_num = 0
    lines = []
    no_empty_row = 0
    with open(matrix_file, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '').split('\t')[0]
            col_num = len(line)
            row_num += 1
            lines.append(line)
            if line != '-' * col_num:
                no_empty_row += 1

    # 过滤掉单拷贝的LTR，因为我们没办法判断它是否是LTR
    if row_num <= 1:
        return ltr_name, False, '', ''

    matrix = [[''] * col_num for i in range(row_num)]
    for row, line in enumerate(lines):
        for col in range(len(line)):
            matrix[row][col] = line[col]

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

    # Starting from column 'align_start', search for 15 effective columns to the left.
    # Count the base composition of each column, in the format of {40: {A: 10, T: 5, C: 7, G: 9, '-': 20}},
    # which indicates the number of different bases in the current column.
    # Based on this, it is easy to determine whether the current column is effective and whether it is a homologous column.
    sliding_window_size = 10
    valid_col_threshold = int(row_num/2)

    if row_num <= 5:
        homo_threshold = 0.95
    elif row_num <= 10:
        homo_threshold = 0.9
    elif row_num <= 50:
        homo_threshold = 0.8
    else:
        homo_threshold = 0.75

    align_start = flanking_len
    align_end = col_num - flanking_len

    homo_boundary_start = search_boundary_homo_v3(valid_col_threshold, align_start, matrix, row_num,
                                             col_num, 'start', homo_threshold, debug, sliding_window_size, flanking_len)
    if homo_boundary_start == -1:
        return ltr_name, False, '', ''

    homo_boundary_end = search_boundary_homo_v3(valid_col_threshold, align_end, matrix, row_num,
                                           col_num, 'end', homo_threshold, debug, sliding_window_size, flanking_len)

    if homo_boundary_end == -1:
        return ltr_name, False, '', ''

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

    # Determine the starting base of the homologous boundary.
    final_boundary_start = -1
    final_boundary_end = -1
    final_cons_seq = ''
    if TE_type == 'tir':
        # Generate a consensus sequence.
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
                    left_tsd_seq, right_tsd_seq = TSDsearch_v5(raw_align_seq, cur_boundary_start, cur_boundary_end)
                    if left_tsd_seq != '':
                        tsd_count += 1
                # TSD的数量超过拷贝数的 1/4，认为是TIR转座子
                # if tsd_count > 10:
                if tsd_count > 10 or float(tsd_count) / row_num >= 0.25:
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
    elif TE_type == 'helitron':
        # 我们对识别到的同源边界两侧扩展 10 bp，然后生成的一致性序列 model_seq
        # Generate a consensus sequence.
        extend_len = 10
        model_seq = ''
        if homo_boundary_start - extend_len >= 0 and homo_boundary_end + 1 + extend_len <= col_num:
            for col_index in range(homo_boundary_start - extend_len, homo_boundary_end + 1 + extend_len):
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
        final_cons_seq = model_seq

    if final_cons_seq == '':
        is_TE = False
    else:
        is_TE = True

    if debug:
        print(matrix_file, is_TE, final_boundary_start, final_boundary_end)
    return ltr_name, is_TE, '', final_cons_seq



def judge_boundary_v6_batch(job_list, full_length_output_dir, flanking_len):
    results = []
    for ltr_name in job_list:
        cur_matrix_file = full_length_output_dir + '/' + ltr_name + '.matrix'
        if not os.path.exists(cur_matrix_file):
            continue
        ltr_name, is_TE, info, final_cons_seq = judge_boundary_v6(cur_matrix_file, ltr_name, flanking_len)
        results.append((ltr_name, is_TE, info, final_cons_seq))
    return results

def judge_boundary_v6(matrix_file, ltr_name, flanking_len):
    align_names, align_contigs = read_matrix_file(matrix_file)

    debug = 0
    col_num = -1
    row_num = 0
    lines = []
    no_empty_row = 0
    with open(matrix_file, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '').split('\t')[0]
            col_num = len(line)
            row_num += 1
            lines.append(line)
            if line != '-' * col_num:
                no_empty_row += 1

    # 过滤掉单拷贝的LTR，因为我们没办法判断它是否是LTR
    if row_num <= 1:
        return ltr_name, False, '', ''

    matrix = [[''] * col_num for i in range(row_num)]
    for row, line in enumerate(lines):
        for col in range(len(line)):
            matrix[row][col] = line[col]

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

    # Starting from column 'align_start', search for 15 effective columns to the left.
    # Count the base composition of each column, in the format of {40: {A: 10, T: 5, C: 7, G: 9, '-': 20}},
    # which indicates the number of different bases in the current column.
    # Based on this, it is easy to determine whether the current column is effective and whether it is a homologous column.
    sliding_window_size = 10
    valid_col_threshold = int(row_num/2)

    if row_num <= 5:
        homo_threshold = 0.95
    elif row_num <= 10:
        homo_threshold = 0.9
    elif row_num <= 50:
        homo_threshold = 0.8
    else:
        homo_threshold = 0.75

    align_start = flanking_len
    align_end = col_num - flanking_len

    homo_boundary_start = search_boundary_homo_v3(valid_col_threshold, align_start, matrix, row_num,
                                             col_num, 'start', homo_threshold, debug, sliding_window_size, flanking_len)
    if homo_boundary_start == -1:
        return ltr_name, False, '', ''

    homo_boundary_end = search_boundary_homo_v3(valid_col_threshold, align_end, matrix, row_num,
                                           col_num, 'end', homo_threshold, debug, sliding_window_size, flanking_len)
    if homo_boundary_end == -1:
        return ltr_name, False, '', ''


    # Iterate through each copy in the multiple sequence alignment, take 15-bp of bases with no gaps above
    # and below the homologous boundary, search for polyA/T within the window, locate the position of polyA/T,
    # and then search for TSDs of 8-bp or more in the upstream 30-bp of the other homologous boundary.

    # 迭代每个拷贝，从尾部搜索第一个polyA/polyT, 提取right TSD
    TSD_sizes = list(range(8, 21))
    end_5_window_size = 20
    end_3_window_size = 20
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

        homology_end_3 = gap_to_nogap[homo_boundary_end]
        end_5 = gap_to_nogap[cur_boundary_start]
        end_3, polyA_seq = find_tail_polyA(align_seq)
        # polyA的边界不能与同源边界相差太多
        if end_3 == -1 or abs(end_3 - homology_end_3) > 10:
            end_3, polyT_seq = find_tail_polyT(align_seq)
        if end_3 == -1 or abs(end_3 - homology_end_3) > 10:
            end_3 = homology_end_3

        found_TSD = False
        TSD_seq = ''
        # After locating the 3' end, attempt to search for a set of TSDs in the side wing (8-20) lengths.
        # Search for the corresponding length of TSD near the 5' end (30 bp), and once found, confirm the final 5' end.
        if end_3 != -1 and end_5 != -1:
            # Obtain all possible TSDs on the side wing of the 3' end.
            # 允许 max_offset = 5
            TSD_list = [(k, align_seq[end_3: end_3 + k]) for k in reversed(TSD_sizes)]
            # max_offset = 1
            # TSD_list = []
            # for i in range(max_offset):
            #     end_3 += i
            #     for k in reversed(TSD_sizes):
            #         TSD_list.append((k, align_seq[end_3: end_3 + k]))

            # Search for TSDs of various lengths near the 5' end (30 bp) (when TSD len >=8, allow 1bp mismatch).
            subsequence = align_seq[end_5 - end_5_window_size: end_5]
            for k, TSD in TSD_list:
                for i in range(0, len(subsequence) - k + 1):
                    kmer = subsequence[i:i + k]
                    dist = 1
                    if k == len(TSD) and k == len(kmer):
                        first_matches = find_near_matches(TSD, kmer, max_l_dist=dist)
                        if len(first_matches) > 0:
                            end_5 = end_5 - end_5_window_size + i + k
                            found_TSD = True
                            TSD_seq = TSD
                            # print(name, TSD_seq)
                            break
                if found_TSD:
                    break
        if found_TSD:
            tsd_count += 1
    # print(tsd_count)
    # Generate a consensus sequence.
    model_seq = ''
    # if tsd_count > 10:
    if tsd_count > 10 or float(tsd_count) / row_num >= 0.25:
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
        model_seq = ''
    else:
        is_TE = True

    if debug:
        print(matrix_file, is_TE, homo_boundary_start, homo_boundary_end)
    return ltr_name, is_TE, '', model_seq


def judge_boundary_v7(matrix_file, ltr_name, TE_type, flanking_len):
    align_names, align_contigs = read_matrix_file(matrix_file)

    debug = 0
    col_num = -1
    row_num = 0
    lines = []
    no_empty_row = 0
    with open(matrix_file, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '').split('\t')[0]
            col_num = len(line)
            row_num += 1
            lines.append(line)
            if line != '-' * col_num:
                no_empty_row += 1

    # 过滤掉单拷贝的LTR，因为我们没办法判断它是否是LTR
    if row_num <= 1:
        return ltr_name, False, '', ''

    matrix = [[''] * col_num for i in range(row_num)]
    for row, line in enumerate(lines):
        for col in range(len(line)):
            matrix[row][col] = line[col]

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

    # Starting from column 'align_start', search for 15 effective columns to the left.
    # Count the base composition of each column, in the format of {40: {A: 10, T: 5, C: 7, G: 9, '-': 20}},
    # which indicates the number of different bases in the current column.
    # Based on this, it is easy to determine whether the current column is effective and whether it is a homologous column.
    sliding_window_size = 10
    valid_col_threshold = int(row_num/2)

    if row_num <= 5:
        homo_threshold = 0.95
    elif row_num <= 10:
        homo_threshold = 0.9
    elif row_num <= 50:
        homo_threshold = 0.8
    else:
        homo_threshold = 0.75

    align_start = flanking_len
    align_end = col_num - flanking_len

    homo_boundary_start = search_boundary_homo_v3(valid_col_threshold, align_start, matrix, row_num,
                                             col_num, 'start', homo_threshold, debug, sliding_window_size)
    if homo_boundary_start == -1:
        return ltr_name, False, '', ''

    homo_boundary_end = search_boundary_homo_v3(valid_col_threshold, align_end, matrix, row_num,
                                           col_num, 'end', homo_threshold, debug, sliding_window_size)

    if homo_boundary_end == -1:
        return ltr_name, False, '', ''

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

    # Determine the starting base of the homologous boundary.
    final_boundary_start = -1
    final_boundary_end = -1
    final_cons_seq = ''
    if TE_type == 'tir':
        TSD_sizes = list(range(2, 12))
        end_5_window_size = 20
        end_5_positions = {}
        end_3_positions = {}
        # Determine if there are two or more TSDs in the align file based on the boundary.
        # If yes, it is a genuine boundary.
        tsd_count = 0
        for name in align_names:
            # 1. Verify if there are bases at the boundary.
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

            end_5 = gap_to_nogap[homo_boundary_start]
            end_3 = gap_to_nogap[homo_boundary_end]

            found_TSD = False
            # After locating the 3' end, attempt to search for a set of TSDs in the side wing (8-20) lengths.
            # Search for the corresponding length of TSD near the 5' end (30 bp), and once found, confirm the final 5' end.
            if end_3 != -1 and end_5 != -1:
                # Obtain all possible TSDs on the side wing of the 3' end.
                # 允许 max_offset = 5
                max_offset = 5
                TSD_list = []
                subsequence = align_seq[end_3 - max_offset: end_3 + end_5_window_size]
                for k in reversed(TSD_sizes):
                    for i in range(0, len(subsequence) - k + 1):
                        kmer = subsequence[i:i + k]
                        cur_end_3 = end_3 - max_offset + i
                        TSD_list.append((k, kmer, cur_end_3))

                all_valid_tsds = []
                # Search for TSDs of various lengths near the 5' end (30 bp) (when TSD len >=8, allow 1bp mismatch).
                subsequence = align_seq[end_5 - end_5_window_size: end_5 + max_offset]
                for k, TSD, cur_end_3 in TSD_list:
                    for i in range(0, len(subsequence) - k + 1):
                        kmer = subsequence[i:i + k]
                        if k >= 8:
                            dist = 1
                        else:
                            dist = 0
                        if k == len(TSD) and k == len(kmer):
                            # first_matches = find_near_matches(TSD, kmer, max_l_dist=dist)
                            # if len(first_matches) > 0:
                            if allow_mismatch(TSD, kmer, allow_mismatch_num=dist):
                                cur_end_5 = end_5 - end_5_window_size + i + k
                                distance = abs(cur_end_5-end_5) + abs(cur_end_3-end_3)
                                all_valid_tsds.append((cur_end_5, cur_end_3, distance, k, kmer, TSD))
                if len(all_valid_tsds) > 0:
                    cur_end_5, cur_end_3, distance, k, left_tsd, right_tsd = sorted(all_valid_tsds, key=lambda x: x[2])[0]

                    gap_end_5 = nogap_to_gap[cur_end_5]
                    if gap_end_5 not in end_5_positions:
                        end_5_positions[gap_end_5] = 0
                    cur_count = end_5_positions[gap_end_5]
                    end_5_positions[gap_end_5] = cur_count + 1

                    gap_end_3 = nogap_to_gap[cur_end_3]
                    if gap_end_3 not in end_3_positions:
                        end_3_positions[gap_end_3] = 0
                    cur_count = end_3_positions[gap_end_3]
                    end_3_positions[gap_end_3] = cur_count + 1

                    tsd_count += 1
        # 取出现次数最多的end_5当作最后的end_5
        homo_boundary_start = sorted(end_5_positions.items(), key=lambda item: item[1], reverse=True)[0][0]
        homo_boundary_end = sorted(end_3_positions.items(), key=lambda item: item[1], reverse=True)[0][0]

        # Generate a consensus sequence.
        model_seq = ''
        # TSD的数量超过拷贝数的 1/4，认为是TIR转座子
        # if tsd_count > 10:
        if tsd_count > 10 or float(tsd_count) / row_num >= 0.25:
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
        final_cons_seq = model_seq
    elif TE_type == 'helitron':
        # 我们对识别到的同源边界两侧扩展 10 bp，然后生成的一致性序列 model_seq
        # Generate a consensus sequence.
        extend_len = 10
        model_seq = ''
        if homo_boundary_start - extend_len >= 0 and homo_boundary_end + 1 + extend_len <= col_num:
            for col_index in range(homo_boundary_start - extend_len, homo_boundary_end + 1 + extend_len):
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
        final_cons_seq = model_seq

    if final_cons_seq == '':
        is_TE = False
    else:
        is_TE = True

    if debug:
        print(matrix_file, is_TE, final_boundary_start, final_boundary_end)
    return ltr_name, is_TE, '', final_cons_seq

def judge_boundary_v8(matrix_file, ltr_name, flanking_len):
    align_names, align_contigs = read_matrix_file(matrix_file)

    debug = 0
    col_num = -1
    row_num = 0
    lines = []
    no_empty_row = 0
    with open(matrix_file, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '').split('\t')[0]
            col_num = len(line)
            row_num += 1
            lines.append(line)
            if line != '-' * col_num:
                no_empty_row += 1

    # 过滤掉单拷贝的LTR，因为我们没办法判断它是否是LTR
    if row_num <= 1:
        return ltr_name, False, '', ''

    matrix = [[''] * col_num for i in range(row_num)]
    for row, line in enumerate(lines):
        for col in range(len(line)):
            matrix[row][col] = line[col]

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

    # Starting from column 'align_start', search for 15 effective columns to the left.
    # Count the base composition of each column, in the format of {40: {A: 10, T: 5, C: 7, G: 9, '-': 20}},
    # which indicates the number of different bases in the current column.
    # Based on this, it is easy to determine whether the current column is effective and whether it is a homologous column.
    sliding_window_size = 10
    valid_col_threshold = int(row_num/2)

    if row_num <= 5:
        homo_threshold = 0.95
    elif row_num <= 10:
        homo_threshold = 0.9
    elif row_num <= 50:
        homo_threshold = 0.8
    else:
        homo_threshold = 0.75

    align_start = flanking_len
    align_end = col_num - flanking_len

    homo_boundary_start = search_boundary_homo_v3(valid_col_threshold, align_start, matrix, row_num,
                                             col_num, 'start', homo_threshold, debug, sliding_window_size)
    if homo_boundary_start == -1:
        return ltr_name, False, '', ''

    homo_boundary_end = search_boundary_homo_v3(valid_col_threshold, align_end, matrix, row_num,
                                           col_num, 'end', homo_threshold, debug, sliding_window_size)
    if homo_boundary_end == -1:
        return ltr_name, False, '', ''


    # Iterate through each copy in the multiple sequence alignment, take 15-bp of bases with no gaps above
    # and below the homologous boundary, search for polyA/T within the window, locate the position of polyA/T,
    # and then search for TSDs of 8-bp or more in the upstream 30-bp of the other homologous boundary.

    # 迭代每个拷贝，从尾部搜索第一个polyA/polyT, 提取right TSD
    TSD_sizes = list(range(8, 21))
    end_5_window_size = 20
    cur_boundary_start = homo_boundary_start
    tsd_count = 0
    end_5_positions = {}
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

        homology_end_3 = gap_to_nogap[homo_boundary_end]
        end_5 = gap_to_nogap[cur_boundary_start]
        end_3, polyA_seq = find_tail_polyA(align_seq)
        # polyA的边界不能与同源边界相差太多
        if end_3 == -1 or abs(end_3 - homology_end_3) > 10:
            end_3, polyT_seq = find_tail_polyT(align_seq)
        if end_3 == -1 or abs(end_3 - homology_end_3) > 10:
            end_3 = homology_end_3

        found_TSD = False
        # After locating the 3' end, attempt to search for a set of TSDs in the side wing (8-20) lengths.
        # Search for the corresponding length of TSD near the 5' end (30 bp), and once found, confirm the final 5' end.
        if end_3 != -1 and end_5 != -1:
            # Obtain all possible TSDs on the side wing of the 3' end.
            # 允许 max_offset = 5
            TSD_list = [(k, align_seq[end_3: end_3 + k]) for k in reversed(TSD_sizes)]

            # Search for TSDs of various lengths near the 5' end (30 bp) (when TSD len >=8, allow 1bp mismatch).
            max_offset = 5
            subsequence = align_seq[end_5 - end_5_window_size: end_5 + max_offset]
            for k, TSD in TSD_list:
                for i in range(0, len(subsequence) - k + 1):
                    kmer = subsequence[i:i + k]
                    dist = 1
                    if k == len(TSD) and k == len(kmer):
                        first_matches = find_near_matches(TSD, kmer, max_l_dist=dist)
                        if len(first_matches) > 0:
                            end_5 = end_5 - end_5_window_size + i + k
                            gap_end_5 = nogap_to_gap[end_5]
                            if gap_end_5 not in end_5_positions:
                                end_5_positions[gap_end_5] = 0
                            cur_count = end_5_positions[gap_end_5]
                            end_5_positions[gap_end_5] = cur_count + 1
                            found_TSD = True
                            break
                if found_TSD:
                    break
        if found_TSD:
            tsd_count += 1
    # 取出现次数最多的end_5当作最后的end_5
    homo_boundary_start = sorted(end_5_positions.items(), key=lambda item: item[1], reverse=True)[0][0]

    # Generate a consensus sequence.
    model_seq = ''
    # if tsd_count > 10:
    if tsd_count > 10 or float(tsd_count) / row_num >= 0.25:
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
        model_seq = ''
    else:
        is_TE = True

    if debug:
        print(matrix_file, is_TE, homo_boundary_start, homo_boundary_end)
    return ltr_name, is_TE, '', model_seq

def find_tail_polyA(sequence, min_length=6):
    for i in range(len(sequence) - (min_length-1), -1, -1):
        six_mer = sequence[i:i + min_length]
        if six_mer == 'AAAAAA':
            return i + min_length, six_mer
    return -1, None

def find_tail_polyT(sequence, min_length=6):
    for i in range(len(sequence) - (min_length-1), -1, -1):
        six_mer = sequence[i:i + min_length]
        if six_mer == 'TTTTTT':
            return i + min_length, six_mer
    return -1, None

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
        return False, 'nb', ''

    align_names, align_contigs = read_fasta(align_file)
    if len(align_names) <= 0:
        if debug:
            print('align file size = 0, ' + align_file)
        return False, '', ''

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
        return False, 'fl1', ''
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
        return False, '', ''

    homo_boundary_end = search_boundary_homo_v3(valid_col_threshold, align_end, matrix, row_num,
                                                col_num, 'end', homo_threshold, debug, sliding_window_size)

    if homo_boundary_end == -1:
        return False, '', ''

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
    return is_TE, '', model_seq

def judge_boundary_v10(cur_seq, align_file, debug, TE_type, plant, result_type):
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
        return False

    align_names, align_contigs = read_fasta(align_file)
    if len(align_names) <= 0:
        if debug:
            print('align file size = 0, ' + align_file)
        return False

    # 3. Take the full-length sequence to generate a consensus sequence.
    # There should be bases both up and down by 10bp at the anchor point.
    full_length_member_names = []
    full_length_member_contigs = {}
    anchor_len = 10
    for name in align_names:
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
        return False
    matrix = [[''] * col_num for i in range(row_num)]
    for row, name in enumerate(full_length_member_names):
        seq = full_length_member_contigs[name]
        for col in range(len(seq)):
            matrix[row][col] = seq[col]

    # Starting from column 'align_start', search for 15 effective columns to the left.
    # Count the base composition of each column, in the format of {40: {A: 10, T: 5, C: 7, G: 9, '-': 20}},
    # which indicates the number of different bases in the current column.
    # Based on this, it is easy to determine whether the current column is effective and whether it is a homologous column.
    sliding_window_size = 100
    valid_col_threshold = int(row_num/2)

    if row_num <= 2:
        homo_threshold = 0.95
    elif row_num <= 5:
        homo_threshold = 0.9
    else:
        homo_threshold = 0.8

    homology_boundary_shift_threshold = 10

    is_TE = (search_boundary_homo_v4(valid_col_threshold, align_start, matrix, row_num, col_num, 'start',
                                    homo_threshold, debug, sliding_window_size)
             and search_boundary_homo_v4(valid_col_threshold, align_end, matrix, row_num, col_num, 'end',
                                         homo_threshold, debug, sliding_window_size))



    # Generate a consensus sequence.
    model_seq = ''
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
    for col_index in range(align_start, align_end+1):
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


    if debug:
        print(align_file, is_TE, align_start, align_end)
    return is_TE


def judge_boundary_v11(cur_seq, align_file, debug, TE_type, plant, result_type):
    # 我们取 align_start 和 align_end 的外侧 100 bp，然后使用 Ninja 进行聚类
    # 假阳性的聚类数量应该比较少，而真实LTR 聚类数量会较多

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
        return False

    align_names, align_contigs = read_fasta(align_file)
    if len(align_names) <= 0:
        if debug:
            print('align file size = 0, ' + align_file)
        return False

    is_TE = True
    # 取 align_start 和 align_end 的外侧 100 bp
    out_threshold = 100

    start_align_file = align_file + '.start'
    start_contigs = {}
    for name in align_names:
        raw_align_seq = align_contigs[name]
        if align_start - out_threshold < 0:
            start_pos = 0
        else:
            start_pos = align_start - out_threshold
        start_seq = raw_align_seq[start_pos: align_start]
        start_contigs[name] = start_seq
    save_dict_to_fasta(start_contigs, start_align_file)

    # 调用 Ninja 对多序列比对再次聚类
    # align_start_out = '/home/hukang/LTR_Benchmarking/LTR_libraries/LtrHomo/dmel_ltr_detector/tmp_dir/LTR_copies_0_0/Chr1:63812-65704.start.maf.fa'
    cluster_file = align_file + '.dat'
    Ninja_command = 'Ninja --in ' + align_file + ' --out ' + cluster_file + ' --out_type c --corr_type m --cluster_cutoff 0.2 --threads 1'
    #os.system(Ninja_command + ' > /dev/null 2>&1')
    # os.system(Ninja_command)
    print(Ninja_command)
    run_command(Ninja_command)

    # 解析聚类文件，生成不同簇
    clusters = read_Ninja_clusters(cluster_file)
    print(len(clusters))

    if debug:
        print(align_file, is_TE, align_start, align_end)
    return is_TE

def search_boundary_homo_v3(valid_col_threshold, pos, matrix, row_num, col_num,
                            type, homo_threshold, debug, sliding_window_size, flanking_len):
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

    search_len = flanking_len
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
                for base in base_map.keys():
                    if base == '-':
                        continue
                    # 修正bug，row_num 替换成 no_gap_num
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
                    # 修正bug，row_num 替换成 no_gap_num
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
                    # 修正bug，row_num 替换成 no_gap_num
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
                    # 修正bug，row_num 替换成 no_gap_num
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

def search_boundary_homo_v4(valid_col_threshold, pos, matrix, row_num, col_num,
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

    search_len = 120
    if type == 'start':
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        is_TE = True

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
                for base in base_map.keys():
                    if base == '-':
                        continue
                    # 修正bug，row_num 替换成 no_gap_num
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
        window = homo_cols[0:0 + sliding_window_size]
        if len(window) > 0:
            # 计算窗口中的同源列的比例
            homo_col_num = 0
            for item in window:
                if item[1]:
                    homo_col_num += 1
            if float(homo_col_num)/len(window) >= homo_threshold:
                is_TE &= True
            else:
                is_TE &= False
        else:
            is_TE &= False

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
                    # 修正bug，row_num 替换成 no_gap_num
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
        window = homo_cols[0:0 + sliding_window_size]
        if len(window) > 0:
            # 计算窗口中的同源列的比例
            homo_col_num = 0
            for item in window:
                if item[1]:
                    homo_col_num += 1
            if float(homo_col_num) / len(window) >= homo_threshold:
                is_TE &= False
            else:
                is_TE &= True
        else:
            is_TE &= False

        return is_TE
    else:
        valid_col_count = 0
        homo_col_count = 0

        max_con_homo = 0
        con_homo = 0
        prev_homo = False

        max_con_no_homo = 0
        con_no_homo = 0
        prev_non_homo = False

        is_TE = True

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
                    # 修正bug，row_num 替换成 no_gap_num
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
        window = homo_cols[0:0 + sliding_window_size]
        if len(window) > 0:
            # 计算窗口中的同源列的比例
            homo_col_num = 0
            for item in window:
                if item[1]:
                    homo_col_num += 1
            if float(homo_col_num) / len(window) >= homo_threshold:
                is_TE &= False
            else:
                is_TE &= True
        else:
            is_TE &= False

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
                    # 修正bug，row_num 替换成 no_gap_num
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
        window = homo_cols[0:0 + sliding_window_size]
        if len(window) > 0:
            # 计算窗口中的同源列的比例
            homo_col_num = 0
            for item in window:
                if item[1]:
                    homo_col_num += 1
            if float(homo_col_num) / len(window) >= homo_threshold:
                is_TE &= True
            else:
                is_TE &= False
        else:
            is_TE &= False

        return is_TE

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

def rename_reference(input, output, chr_name_map):
    names, contigs = read_fasta(input)
    chr_name_dict = {}
    ref_index = 0
    with open(output, 'w') as f_save:
        for name in names:
            seq = contigs[name]
            new_name = 'Chr'+str(ref_index)
            f_save.write('>'+new_name+'\n'+seq+'\n')
            ref_index += 1
            chr_name_dict[new_name] = name
    f_save.close()
    with open(chr_name_map, 'w') as f_save:
        for new_name in chr_name_dict.keys():
            f_save.write(new_name+'\t'+chr_name_dict[new_name]+'\n')
    f_save.close()

def parse_ltr_log(ltr_log):
    ltr_pos = {}
    cur_seq_name = None
    with open(ltr_log, 'r') as f_r:
        for line in f_r:
            if line.startswith('load sequence'):
                seq_name = line.split(' ')[-1].split('\t')[0]
            elif 'Length ltr=' in line:
                parts = line.strip().split(' ')[0].split('..')
                lLTR_info = parts[0].replace('(', '').replace(')', '').split(',')
                rLTR_info = parts[1].replace('(', '').replace(')', '').split(',')
                lLTR_start = lLTR_info[0]
                lLTR_end = lLTR_info[1]
                rLTR_start = rLTR_info[0]
                rLTR_end = rLTR_info[1]
                pos = (lLTR_start, lLTR_end, rLTR_start, rLTR_end)
                if seq_name is not None:
                    ltr_pos[seq_name] = pos
                    seq_name = None
    return ltr_pos

def identify_LTR_terminals(split_file, output_dir, tool_dir):
    ltr_file = split_file + '.ltr'
    ltr_log = ltr_file + '.log'
    tir_file = split_file + '.itr'
    tir_log = tir_file + '.log'

    ltrsearch_command = 'cd ' + output_dir + ' && ' + tool_dir + '/ltrsearch -l 100 ' + split_file + ' > ' + ltr_log
    run_command(ltrsearch_command)

    # 解析.log 文件，读取 'load sequence' 和 'Length ltr=' 标志
    LTR_info = parse_ltr_log(ltr_log)

    # Update the header of the split_file, adding two columns LTR:1-206,4552-4757 TIR:1-33,3869-3836.
    update_split_file = split_file + '.updated'
    update_contigs = {}
    names, contigs = read_fasta_v1(split_file)
    for name in names:
        orig_name = name.split('\t')[0]
        LTR_str = 'LTR:'
        lLTR_start, lLTR_end, rLTR_start, rLTR_end = LTR_info[orig_name]
        LTR_str += str(lLTR_start) + '-' + str(lLTR_end) + ',' + str(rLTR_start) + '-' + str(rLTR_end)
        update_name = name + '\t' + LTR_str
        update_contigs[update_name] = contigs[name]
    store_fasta(update_contigs, update_split_file)
    return update_split_file

def split_fasta(cur_path, output_dir, num_chunks):
    split_files = []

    if os.path.exists(output_dir):
        os.system('rm -rf ' + output_dir)
    os.makedirs(output_dir)

    names, contigs = read_fasta_v1(cur_path)
    num_names = len(names)
    chunk_size = num_names // num_chunks

    for i in range(num_chunks):
        chunk_start = i * chunk_size
        chunk_end = chunk_start + chunk_size if i < num_chunks - 1 else num_names
        chunk = names[chunk_start:chunk_end]
        output_path = output_dir + '/out_' + str(i) + '.fa'
        with open(output_path, 'w') as out_file:
            for name in chunk:
                seq = contigs[name]
                out_file.write('>'+name+'\n'+seq+'\n')
        split_files.append(output_path)
    return split_files

def generate_LTR_terminal_info(data_path, work_dir, tool_dir, threads):
    output_dir = work_dir + '/temp'
    # Split the file into threads blocks.
    split_files = split_fasta(data_path, output_dir, threads)

    # Parallelize the identification of LTR and TIR.
    cur_update_path = data_path + '.update'
    os.system('rm -f ' + cur_update_path)
    with ProcessPoolExecutor(threads) as executor:
        futures = []
        for split_file in split_files:
            future = executor.submit(identify_LTR_terminals, split_file, output_dir, tool_dir)
            futures.append(future)
        executor.shutdown(wait=True)

        is_exit = False
        for future in as_completed(futures):
            update_split_file = future.result()
            os.system('cat ' + update_split_file + ' >> ' + cur_update_path)

    return cur_update_path

# Parse output of LTR_harvest
def get_LTR_seq_from_scn(genome, scn_path, ltr_terminal, ltr_internal, ltr_intact, dirty_dicts, ltr_intact_list, miu, coverage_threshold=0.95):
    ref_names, ref_contigs = read_fasta(genome)
    ltr_lines = []
    with open(scn_path, 'r') as f_r:
        for i, line in enumerate(f_r):
            if line.startswith('#'):
                continue
            else:
                line = line.replace('\n', '')
                ltr_lines.append(line)
    f_r.close()

    deredundant_lines = []
    redundant_index = set()
    # 去除坐标重复的LTR
    for i in range(len(ltr_lines) - 1):
        if i in redundant_index:
            continue
        cur_line = ltr_lines[i]
        parts = cur_line.split(' ')
        chr_name = parts[11]
        chr_start = int(parts[3]) - 1
        chr_end = int(parts[7])
        cur_frag = (chr_start, chr_end)
        for j in range(i+1, len(ltr_lines)):
            next_line = ltr_lines[j]
            parts = next_line.split(' ')
            next_chr_name = parts[11]
            next_chr_start = int(parts[3]) - 1
            next_chr_end = int(parts[7])
            next_frag = (next_chr_start, next_chr_end)
            overlap_len = get_overlap_len(next_frag, cur_frag)
            if (next_chr_name == chr_name
                    and overlap_len / abs(next_frag[1] - next_frag[0]) >= coverage_threshold
                    and overlap_len / abs(cur_frag[1] - cur_frag[0]) >= coverage_threshold):
                redundant_index.add(j)
            else:
                break
        deredundant_lines.append(cur_line)
    # print(len(ltr_lines), len(deredundant_lines))

    LTR_terminals = {}
    LTR_ints = {}
    LTR_intacts = {}
    for line in deredundant_lines:
        parts = line.split(' ')
        LTR_start = int(parts[0])
        LTR_end = int(parts[1])
        chr_name = parts[11]
        lLTR_start = int(parts[3])
        lLTR_end = int(parts[4])
        rLTR_start = int(parts[6])
        rLTR_end = int(parts[7])
        lLTR_seq = ref_contigs[chr_name][lLTR_start - 1: lLTR_end]
        rLTR_seq = ref_contigs[chr_name][rLTR_start - 1: rLTR_end]
        LTR_int_seq = ref_contigs[chr_name][lLTR_end: rLTR_start - 1]
        intact_seq = ref_contigs[chr_name][lLTR_start - 1: rLTR_end]

        cur_name = chr_name + '-' + str(lLTR_start) + '-' + str(lLTR_end)
        if (len(lLTR_seq) > 0 and len(rLTR_seq) > 0 and len(LTR_int_seq) > 0
                and 'NNNNNNNNNN' not in lLTR_seq and 'NNNNNNNNNN' not in rLTR_seq and 'NNNNNNNNNN' not in LTR_int_seq):
            lLTR_name = chr_name + ':' + str(LTR_start) + '..' + str(LTR_end) + '-lLTR' + '#LTR'
            LTR_terminals[lLTR_name] = lLTR_seq
            rLTR_name = chr_name + ':' + str(LTR_start) + '..' + str(LTR_end) + '-rLTR' + '#LTR'
            LTR_terminals[rLTR_name] = rLTR_seq
            intact_name = chr_name + ':' + str(LTR_start) + '..' + str(LTR_end)
            LTR_intacts[intact_name] = intact_seq
            # 如果当前是dirty (内部序列包含其他full-length LTR)，就过滤该内部序列
            if cur_name not in dirty_dicts:
                LTR_int_name = chr_name + ':' + str(LTR_start) + '..' + str(LTR_end) + '-int' + '#LTR'
                LTR_ints[LTR_int_name] = LTR_int_seq
    store_fasta(LTR_ints, ltr_internal)
    store_fasta(LTR_terminals, ltr_terminal)
    store_fasta(LTR_intacts, ltr_intact)

    # 生成LTR_retriever .pass.list格式
    # #LTR_loc        Category        Motif   TSD     5_TSD 3_TSD       Internal        Identity      Strand  SuperFamily  TE_type     Insertion_Time
    # chr_1:146286..139837    pass    motif:TGCA      TSD:GTATA       139832..139836  146287..146291  IN:140816..145316       1.0000  -       Copia   LTR     0
    with open(ltr_intact_list, 'w') as f_save:
        # header = '#LTR_loc\tMotif\tTSD\tInternal\tIdentity\tInsertion_Time'
        # f_save.write(header + '\n')
        for line in deredundant_lines:
            parts = line.split(' ')
            LTR_start = int(parts[0])
            LTR_end = int(parts[1])
            chr_name = parts[11]
            ref_seq = ref_contigs[chr_name]
            lLTR_start = int(parts[3])
            lLTR_end = int(parts[4])
            rLTR_start = int(parts[6])
            rLTR_end = int(parts[7])
            identity = float(parts[9]) / 100
            insertion_time = int(estimate_insert_time(identity, float(miu)))
            ltr_name = chr_name + ':' + str(LTR_start) + '..' + str(LTR_end)
            ltr_name, has_structure, tsd_seq = is_ltr_has_structure(ltr_name, line, ref_seq)
            motif = ref_seq[lLTR_start - 1: lLTR_start - 1 + 2] + ref_seq[rLTR_end - 2: rLTR_end]
            if tsd_seq == '':
                tsd_seq = 'NA'
            cur_record = ltr_name + '\t' + 'pass' + '\t' + 'motif:' + motif + '\t' + 'TSD:' + tsd_seq + '\t' \
                         + 'NA..NA\tNA..NA\t' + 'IN:' + str(lLTR_end) + '..' + str(rLTR_start - 1) + '\t' \
                         + str(identity) + '\t.\tNA\tNA' + '\t' + str(insertion_time)
            f_save.write(cur_record + '\n')

def estimate_insert_time(identity, miu):
    # 估计序列的插入时间
    d = 1 - identity
    K=  -3/4*math.log(1-d*4/3)
    T = K/(2*miu)
    return T

def merge_intervals(intervals):
    # 如果列表为空，直接返回空列表
    if not intervals:
        return []
    # 按照 start 值对 intervals 进行排序
    intervals.sort(key=lambda x: x[0])
    # 初始化合并后的区间列表
    merged = [intervals[0]]
    for current in intervals:
        # 获取上一个合并后的区间
        last = merged[-1]
        # 如果当前区间与上一个区间重叠或相邻，合并它们
        if current[0] <= last[1]:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            # 否则，添加当前区间
            merged.append(current)
    return merged

def get_overlap_len(seq1, seq2):
    """Calculate the overlap length between two sequences."""
    overlap_len = min(seq1[1], seq2[1]) - max(seq1[0], seq2[0])
    return overlap_len if overlap_len > 0 else 0

def merge_overlap_seq(seq1, seq2):
    return (min(seq1[0], seq2[0]), max(seq1[1], seq2[1]))

def multi_process_align_v1(query_path, subject_path, blastnResults_path, tmp_blast_dir, threads, chrom_length, coverage_threshold, category, is_full_length, is_removed_dir=True):
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
        job = ex.submit(multiple_alignment_blast_v1, file, tools_dir, coverage_threshold, category, chrom_length, is_full_length)
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

                    # if seq_name == 'chr_11_15708136-15719185-lLTR' and chr_name == 'Chr14' and start == 21886838 and end == 21890006:
                    #     print('i')

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


def multi_process_align_v3(query_path, subject_path, blastnResults_path, tmp_blast_dir, threads, coverage_threshold, is_removed_dir=True):
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
        job = ex.submit(multiple_alignment_blast_v3, file, tools_dir, coverage_threshold)
        jobs.append(job)
    ex.shutdown(wait=True)

    # 合并所有进程的结果，总体去除冗余
    full_length_annotation = {}
    for job in as_completed(jobs):
        cur_full_length_annotation = job.result()
        full_length_annotation.update(cur_full_length_annotation)

    if is_removed_dir:
        os.system('rm -rf ' + tmp_blast_dir)
    return full_length_annotation


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
        skip_gap = query_len * (1 - full_length_threshold)
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

def multi_process_align_v2(query_path, subject_path, blastnResults_path, tmp_blast_dir, threads, chrom_length, coverage_threshold, category, TRsearch_dir, is_removed_dir=True):
    if is_removed_dir:
        os.system('rm -rf ' + tmp_blast_dir)

    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    if os.path.exists(blastnResults_path):
        os.remove(blastnResults_path)


    # 由于 blastn 未能将一些差异性的一致性序列比对到应有的位置，因此我们调用 RepeatMasker 来进行比对
    intact_dir = tmp_blast_dir + '/intact_tmp'
    divergence_threshold = 20
    full_length_threshold = 0.8
    search_struct = False
    full_length_annotations, copies_direct = get_full_length_copies_RM(query_path, subject_path, intact_dir, threads,
                                                                       divergence_threshold,
                                                                       full_length_threshold, search_struct,
                                                                       TRsearch_dir)
    lines = []
    for seq_name in full_length_annotations.keys():
        for copy in full_length_annotations[seq_name]:
            parts = copy[0].split(':')
            chr_name = parts[0]
            pos_parts = parts[1].split('-')
            chr_start = int(pos_parts[0]) + 1
            chr_end = int(pos_parts[1])
            lines.append((seq_name, chr_name, chr_start, chr_end))


    lines = list(lines)
    sorted_lines = sorted(lines, key=lambda x: (x[1], x[2], x[3]))
    test_fragments = {}
    for line in sorted_lines:
        seq_name = line[0]
        chr_name = line[1]
        chr_start = line[2]
        chr_end = line[3]
        if chr_name not in test_fragments:
            test_fragments[chr_name] = []
        fragments = test_fragments[chr_name]
        fragments.append((chr_start, chr_end, seq_name))

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
                segment_frags.append([start, end, seq_name])

    with open(blastnResults_path, 'w') as f_save:
        for chr_name in chr_segments.keys():
            cur_chr_segments = chr_segments[chr_name]
            for seg_index in cur_chr_segments.keys():
                segment_frags = cur_chr_segments[seg_index]
                for frag in segment_frags:
                    new_line = frag[2] + '\t' + chr_name + '\t' + '-1' + '\t' + '-1' + '\t' + '-1' + '\t' + '-1' + '\t' + '-1' + '\t' + '-1' + '\t' + str(frag[0]) + '\t' + str(frag[1]) + '\t' + '-1' + '\t' + '-1' + '\n'
                    f_save.write(new_line)

    if is_removed_dir:
        os.system('rm -rf ' + tmp_blast_dir)

def map_fragment(start, end, chunk_size):
    start_chunk = start // chunk_size
    end_chunk = end // chunk_size

    if start_chunk == end_chunk:
        return start_chunk
    elif abs(end_chunk * chunk_size - start) < abs(end - end_chunk * chunk_size):
        return end_chunk
    else:
        return start_chunk

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

def multiple_alignment_blast_v1(repeats_path, tools_dir, coverage_threshold, category, chrom_length, is_full_length):
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


def multiple_alignment_blast_v3(repeats_path, tools_dir, coverage_threshold):
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
    full_length_annotations = generate_full_length_out_v3(blastn2Results_path, full_length_out, split_repeats_path, genome_path, tmp_dir, tools_dir,
                             coverage_threshold)

    return full_length_annotations

def generate_full_length_out(BlastnOut, full_length_out, TE_lib, reference, tmp_output_dir, tools_dir, full_length_threshold, category, is_full_length):
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)
    filter_tmp_out = filter_out_by_category(BlastnOut, tmp_output_dir, category)

    threads = 1
    divergence_threshold = 20
    search_struct = False
    full_length_annotations, copies_direct, all_query_copies = get_full_length_copies_from_blastn(TE_lib, reference, filter_tmp_out,
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

def generate_full_length_out_v3(BlastnOut, full_length_out, TE_lib, reference, tmp_output_dir, tools_dir, full_length_threshold):
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)

    threads = 1
    divergence_threshold = 20
    search_struct = False
    full_length_annotations, copies_direct = get_full_length_copies_from_blastn_v3(TE_lib, reference, BlastnOut,
                                                                             tmp_output_dir, threads,
                                                                             divergence_threshold,
                                                                             full_length_threshold,
                                                                             search_struct, tools_dir)

    return full_length_annotations

def filter_out_by_category(TE_out, tmp_output_dir, category):
    tmp_out= tmp_output_dir + '/tmp.out'
    os.system('cp ' + TE_out + ' ' + tmp_out)
    if category == 'Total':
        return tmp_out
    else:
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

def get_full_length_copies_from_blastn_v3(TE_lib, reference, blastn_out, tmp_output_dir, threads, divergence_threshold,
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
        skip_gap = query_len * (1 - full_length_threshold)

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
            cur_subject_len = repeat[5]
            query_coverage = float(cur_query_len) / orig_query_len
            query_subject_ratio = float(min(cur_query_len, cur_subject_len))/ max(cur_query_len, cur_subject_len)
            if query_coverage >= full_length_threshold and query_subject_ratio >= full_length_threshold:
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
            EAHelitron = os.getcwd() + '/../bin/EAHelitron-master'
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

def identify_terminals(split_file, output_dir, tool_dir):
    #ltr_log = split_file + '.ltr.log'
    tir_log = split_file + '.itr.log'
    #ltrsearch_command = 'cd ' + output_dir + ' && ' + tool_dir + '/ltrsearch -l 50 ' + split_file + ' > ' + ltr_log
    itrsearch_command = 'cd ' + output_dir + ' && ' + tool_dir + '/itrsearch -i 0.7 -l 7 ' + split_file+ ' > ' + tir_log
    #run_command(ltrsearch_command)
    run_command(itrsearch_command)
    ltr_file = split_file + '.ltr'
    tir_file = split_file + '.itr'

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

def run_command(command):
    subprocess.run(command, check=True, shell=True)

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

def is_recombination(query_seq, subject_seq, candidate_index):
    query_len = len(query_seq)

    # 创建临时文件
    query_file = tempfile.mkstemp()[1]
    subject_file = tempfile.mkstemp()[1]
    output_file = tempfile.mkstemp()[1]

    # 将查询序列和参考序列写入临时文件
    with open(query_file, 'w') as f:
        f.write(query_seq)
    with open(subject_file, 'w') as f:
        f.write(subject_seq)

    # 运行 BLAST，将输出写入文件
    exec_command = [
        "blastn",
        "-subject", subject_file,
        "-query", query_file,
        "-outfmt", "6",
        "-out", output_file
    ]
    process = subprocess.run(exec_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # 检查 BLAST 是否成功运行
    if process.returncode != 0:
        print("BLAST error:", process.stderr)
        # 清理临时文件
        os.remove(query_file)
        os.remove(subject_file)
        os.remove(output_file)
        return False, candidate_index

    # 从输出文件中读取结果
    with open(output_file, 'r') as f:
        lines = f.readlines()

    # 清理临时文件
    os.remove(query_file)
    os.remove(subject_file)
    os.remove(output_file)

    # 处理结果
    for line in lines:
        parts = line.strip().split('\t')
        if len(parts) != 12:
            continue
        query_start = int(parts[6])
        query_end = int(parts[7])
        alignment_length = abs(query_end - query_start) + 1
        if alignment_length / query_len >= 0.95:
            return True, candidate_index
    return False, candidate_index

def filter_single_ltr(output_path, intact_output_path, leftLtr2Candidates, ltr_lines, reference, flanking_len, tmp_output_dir, split_ref_dir, threads, log):
    true_ltrs = {}
    with open(output_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ltr_name = parts[0]
            is_ltr = int(parts[1])
            if is_ltr:
                true_ltrs[ltr_name] = is_ltr

    confident_lines = []
    for name in leftLtr2Candidates.keys():
        # if name not in FP_ltrs:
        if name in true_ltrs:
            candidate_index = leftLtr2Candidates[name]
            line = ltr_lines[candidate_index]
            confident_lines.append((name, line))

    intact_ltrs = {}
    intact_ltr_path = tmp_output_dir + '/intact_ltr.fa'
    ref_names, ref_contigs = read_fasta(reference)
    for name, line in confident_lines:
        parts = line.split(' ')
        ltr_start = int(parts[0])
        ltr_end = int(parts[1])
        chr_name = parts[11]
        left_ltr_start = int(parts[3])
        left_ltr_end = int(parts[4])
        right_ltr_start = int(parts[6])
        right_ltr_end = int(parts[7])

        intact_ltr_seq =  ref_contigs[chr_name][left_ltr_start-1: right_ltr_end]
        if len(intact_ltr_seq) > 0:
            # intactLTR_name = chr_name + '_' + str(left_ltr_start) + '-' + str(right_ltr_end) + '-intactLTR' + '#LTR'
            intact_ltrs[name] = intact_ltr_seq
    store_fasta(intact_ltrs, intact_ltr_path)

    # 对于全长拷贝数 >=2 的LTR序列，检测是否由于输入基因组包含多个冗余 contigs 造成
    # 即，对全长LTR生成窗口，然后调用规则方法判断窗口两侧是否存在同源性
    temp_dir = tmp_output_dir + '/intact_ltr_filter'
    output_dir = tmp_output_dir + '/intact_ltr_both_frames'
    full_length_output_dir = tmp_output_dir + '/intact_ltr_full_length_frames'
    generate_both_ends_frame_for_intactLTR(intact_ltr_path, reference, flanking_len, threads, temp_dir, output_dir, full_length_output_dir,
                                      split_ref_dir, log)
    sliding_window_size = 20
    type = 'intact copy'
    judge_ltr_from_both_ends_frame_for_intactLTR(output_dir, intact_output_path, threads, type, log, sliding_window_size=sliding_window_size)


def filter_single_copy_ltr(output_path, single_copy_internals_file, ltr_copies, internal_ltrs,
                           ltr_protein_db, other_protein_db, tmp_output_dir, threads,
                           reference, leftLtr2Candidates, ltr_lines, log, debug):
    # # 我们只保留具有 TSD + TG...CA 结构的单拷贝 或者 内部具有完整protein的。
    # 我们只保留内部具有完整或 >=100 bp protein

    # 判断单拷贝的内部序列是否有完整的蛋白质
    single_copy_internals = {}
    for name in ltr_copies.keys():
        if len(ltr_copies[name]) <= 1:
            single_copy_internals[name] = internal_ltrs[name]
    store_fasta(single_copy_internals, single_copy_internals_file)

    # temp_dir = tmp_output_dir + '/ltr_domain'
    # query_protein_types = get_domain_info_v2(single_copy_internals_file, protein_db, threads, temp_dir, tool_dir)

    temp_dir = tmp_output_dir + '/ltr_domain'
    output_table = single_copy_internals_file + '.ltr_domain'
    get_domain_info(single_copy_internals_file, ltr_protein_db, output_table, threads, temp_dir)
    is_single_ltr_has_intact_protein = {}
    protein_names, protein_contigs = read_fasta(ltr_protein_db)
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
                is_single_ltr_has_intact_protein[te_name] = True

    if not debug:
        os.system('rm -rf ' + temp_dir)
        os.system('rm -f ' + output_table)

    temp_dir = tmp_output_dir + '/other_domain'
    output_table = single_copy_internals_file + '.other_domain'
    get_domain_info(single_copy_internals_file, other_protein_db, output_table, threads, temp_dir)
    is_single_ltr_has_intact_other_protein = {}
    protein_names, protein_contigs = read_fasta(other_protein_db)
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
                is_single_ltr_has_intact_other_protein[te_name] = True

    if not debug:
        os.system('rm -rf ' + temp_dir)
        os.system('rm -f ' + output_table)

    # 判断单拷贝LTR是否有TSD结构
    is_single_ltr_has_structure = {}
    no_structure_single_count = 0
    ref_names, ref_contigs = read_fasta(reference)
    for ltr_name in single_copy_internals.keys():
        candidate_index = leftLtr2Candidates[ltr_name]
        line = ltr_lines[candidate_index]
        parts = line.split(' ')
        chr_name = parts[11]
        ref_seq = ref_contigs[chr_name]
        ltr_name, has_structure, tsd_seq = is_ltr_has_structure(ltr_name, line, ref_seq)
        is_single_ltr_has_structure[ltr_name] = has_structure
        if not has_structure:
            no_structure_single_count += 1
    if log is not None:
        log.logger.info('Filter the number of no structure single copy LTR: ' + str(no_structure_single_count) + ', remaining single copy LTR num: ' + str(len(single_copy_internals)-no_structure_single_count))


    remain_intact_count = 0
    filtered_intact_count = 0
    with open(output_path, 'w') as f_save:
        for name in ltr_copies.keys():
            if len(ltr_copies[name]) >= 2:
                cur_is_ltr = 1
                remain_intact_count += 1
            else:
                if name in is_single_ltr_has_intact_other_protein and is_single_ltr_has_intact_other_protein[name]:
                    cur_is_ltr = 0
                    filtered_intact_count += 1
                else:
                    if (name in is_single_ltr_has_intact_protein and is_single_ltr_has_intact_protein[name]) \
                            and (name in is_single_ltr_has_structure and is_single_ltr_has_structure[name]):
                        cur_is_ltr = 1
                        remain_intact_count += 1
                    else:
                        cur_is_ltr = 0
                        filtered_intact_count += 1
            f_save.write(name + '\t' + str(cur_is_ltr) + '\n')
    if log is not None:
        log.logger.info('Filter the number of non-intact LTR: ' + str(filtered_intact_count) + ', remaining LTR num: ' + str(remain_intact_count))

    if not debug:
        os.system('rm -f ' + single_copy_internals_file)

def filter_ltr_by_copy_num(output_path, leftLtr2Candidates, ltr_lines, reference, tmp_output_dir, split_ref_dir, threads, full_length_threshold, debug):
    true_ltrs = {}
    with open(output_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ltr_name = parts[0]
            is_ltr = int(parts[1])
            if is_ltr:
                true_ltrs[ltr_name] = is_ltr

    confident_lines = []
    for name in leftLtr2Candidates.keys():
        # if name not in FP_ltrs:
        if name in true_ltrs:
            candidate_index = leftLtr2Candidates[name]
            line = ltr_lines[candidate_index]
            confident_lines.append((name, line))

    internal_ltrs = {}
    intact_ltrs = {}
    intact_ltr_path = tmp_output_dir + '/intact_ltr.fa'
    ref_names, ref_contigs = read_fasta(reference)
    for name, line in confident_lines:
        parts = line.split(' ')
        chr_name = parts[11]
        left_ltr_start = int(parts[3])
        left_ltr_end = int(parts[4])
        right_ltr_start = int(parts[6])
        right_ltr_end = int(parts[7])

        ref_seq = ref_contigs[chr_name]

        intact_ltr_seq = ref_seq[left_ltr_start-1: right_ltr_end]
        internal_ltr_seq = ref_seq[left_ltr_end: right_ltr_start - 1]
        internal_ltrs[name] = internal_ltr_seq
        if len(intact_ltr_seq) > 0:
            intact_ltrs[name] = intact_ltr_seq
    store_fasta(intact_ltrs, intact_ltr_path)

    temp_dir = tmp_output_dir + '/intact_ltr_filter'
    ltr_copies = filter_ltr_by_copy_num_sub(intact_ltr_path, threads, temp_dir, split_ref_dir, full_length_threshold, max_copy_num=10)

    if not debug:
        os.system('rm -rf ' + temp_dir)
        os.system('rm -f ' + intact_ltr_path)

    # 我们获取拷贝之后再和原始结果计算overlap，如果overlap超过 95% 才算一个真的全长拷贝，否则不算。
    # 经常会出现获取了两个拷贝，但是实际上都是同一个拷贝(坐标overlap或者来自冗余contig)，因此我们要对拷贝去冗余
    intact_ltr_copies = get_intact_ltr_copies(ltr_copies, ltr_lines, full_length_threshold, reference)

    # 过滤来自冗余contig的拷贝，即取拷贝的左侧100bp+右侧100bp的序列，任意一侧能够很好比对就说明这两个全长拷贝是来自冗余contig造成的
    temp_dir = tmp_output_dir + '/intact_ltr_deredundant'
    # intact_ltr_copies = remove_copies_from_redundant_contig(intact_ltr_copies, reference, temp_dir, threads)
    intact_ltr_copies = remove_copies_from_redundant_contig_v1(intact_ltr_copies, reference, temp_dir, threads)

    if not debug:
        os.system('rm -rf ' + temp_dir)

    return intact_ltr_copies, internal_ltrs

def remove_copies_from_redundant_contig_v1(intact_ltr_copies, reference, temp_dir, threads, flanking_len=100):
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    ref_names, ref_contigs = read_fasta(reference)

    all_list = []
    for ltr_name in intact_ltr_copies.keys():
        cur_left_copy_path = temp_dir + '/' + ltr_name + '.left_copies'
        cur_right_copy_path = temp_dir + '/' + ltr_name + '.right_copies'
        copies = intact_ltr_copies[ltr_name]
        cur_left_copy_contigs = {}
        cur_right_copy_contigs = {}
        copy_name_2_copy_dict = {}
        for cur_copy in copies:
            # 将所有拷贝都存成文件，然后
            cur_copy_name = str(cur_copy[0]) + '-' + str(cur_copy[1]) + '-' + str(cur_copy[2])
            ref_seq = ref_contigs[cur_copy[0]]
            cur_copy_left_flank_seq = ref_seq[cur_copy[1] - flanking_len: cur_copy[1]]
            cur_copy_right_flank_seq = ref_seq[cur_copy[2]: cur_copy[2] + flanking_len]
            cur_left_copy_contigs[cur_copy_name] = cur_copy_left_flank_seq
            cur_right_copy_contigs[cur_copy_name] = cur_copy_right_flank_seq
            copy_name_2_copy_dict[cur_copy_name] = cur_copy
        store_fasta(cur_left_copy_contigs, cur_left_copy_path)
        store_fasta(cur_right_copy_contigs, cur_right_copy_path)
        all_list.append((ltr_name, cur_left_copy_contigs, cur_left_copy_path, cur_right_copy_contigs, cur_right_copy_path, copy_name_2_copy_dict))

    part_size = len(all_list) // threads
    divided_job_list = []
    # 划分前 n-1 部分
    for i in range(threads - 1):
        divided_job_list.append(all_list[i * part_size: (i + 1) * part_size])
    # 最后一部分包含剩余的所有元素
    divided_job_list.append(all_list[(threads - 1) * part_size:])

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for job_list in divided_job_list:
        job = ex.submit(get_non_redundant_copies_v1_batch, job_list)
        jobs.append(job)
    ex.shutdown(wait=True)
    intact_ltr_copies = {}
    for job in as_completed(jobs):
        results = job.result()
        for ltr_name, cur_intact_ltr_copies in results:
            intact_ltr_copies[ltr_name] = cur_intact_ltr_copies
    return intact_ltr_copies


def get_non_redundant_copies_v1_batch(job_list):
    results = []
    for ltr_name, cur_left_copy_contigs, cur_left_copy_path, cur_right_copy_contigs, cur_right_copy_path, copy_name_2_copy_dict in job_list:
        ltr_name, intact_copies = get_non_redundant_copies_v1(ltr_name, cur_left_copy_contigs, cur_left_copy_path, cur_right_copy_contigs,
                                    cur_right_copy_path, copy_name_2_copy_dict)
        results.append((ltr_name, intact_copies))
    return results

def get_non_redundant_copies_v1(ltr_name, cur_left_copy_contigs, cur_left_copy_path, cur_right_copy_contigs, cur_right_copy_path, copy_name_2_copy_dict):
    blastn_command = 'blastn -query ' + cur_left_copy_path + ' -subject ' + cur_left_copy_path + ' -num_threads ' + str(1) + ' -outfmt 6 '
    result = subprocess.run(blastn_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            executable='/bin/bash')
    redundant_copies = set()
    if result.returncode == 0:
        lines = result.stdout.split('\n')
        for line in lines:
            parts = line.split('\t')
            if len(parts) != 12:
                continue
            query_name = parts[0]
            subject_name = parts[1]
            if query_name == subject_name:
                continue
            query_len = len(cur_left_copy_contigs[query_name])
            subject_len = len(cur_left_copy_contigs[subject_name])

            query_start = int(parts[6])
            query_end = int(parts[7])
            subject_start = int(parts[8])
            subject_end = int(parts[9])
            if abs(query_end - query_start) / query_len >= 0.95 and abs(
                    subject_end - subject_start) / subject_len >= 0.95:
                redundant_copies.add(query_name)
                redundant_copies.add(subject_name)
    blastn_command = 'blastn -query ' + cur_right_copy_path + ' -subject ' + cur_right_copy_path + ' -num_threads ' + str(1) + ' -outfmt 6 '
    result = subprocess.run(blastn_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            executable='/bin/bash')
    if result.returncode == 0:
        lines = result.stdout.split('\n')
        for line in lines:
            parts = line.split('\t')
            if len(parts) != 12:
                continue
            query_name = parts[0]
            subject_name = parts[1]
            if query_name == subject_name:
                continue
            query_len = len(cur_right_copy_contigs[query_name])
            subject_len = len(cur_right_copy_contigs[subject_name])

            query_start = int(parts[6])
            query_end = int(parts[7])
            subject_start = int(parts[8])
            subject_end = int(parts[9])
            if abs(query_end - query_start) / query_len >= 0.95 and abs(
                    subject_end - subject_start) / subject_len >= 0.95:
                redundant_copies.add(query_name)
                redundant_copies.add(subject_name)

    # 遍历 cur_copy_contigs，取所有非冗余的拷贝
    intact_copies = []
    for cur_copy_name in cur_left_copy_contigs.keys():
        if cur_copy_name not in redundant_copies:
            intact_copies.append(copy_name_2_copy_dict[cur_copy_name])
    # 在冗余拷贝中任取一个加入到拷贝中
    if len(redundant_copies) > 0:
        intact_copies.append(copy_name_2_copy_dict[list(redundant_copies)[0]])
    return ltr_name, intact_copies

def remove_copies_from_redundant_contig(intact_ltr_copies, reference, temp_dir, threads, flanking_len=100):
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    ref_names, ref_contigs = read_fasta(reference)

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for ltr_name in intact_ltr_copies.keys():
        cur_copy_path = temp_dir + '/' + ltr_name + '.copies'
        copies = intact_ltr_copies[ltr_name]
        cur_copy_contigs = {}
        copy_name_2_copy_dict = {}
        for cur_copy in copies:
            # 将所有拷贝都存成文件，然后
            cur_copy_name = str(cur_copy[0]) + '-' + str(cur_copy[1]) + '-' + str(cur_copy[2])
            ref_seq = ref_contigs[cur_copy[0]]
            cur_copy_flank_seq = ref_seq[cur_copy[1] - flanking_len: cur_copy[1]] + ref_seq[cur_copy[2]: cur_copy[2] + flanking_len]
            cur_copy_contigs[cur_copy_name] = cur_copy_flank_seq
            copy_name_2_copy_dict[cur_copy_name] = cur_copy
        store_fasta(cur_copy_contigs, cur_copy_path)

        job = ex.submit(get_non_redundant_copies, ltr_name, cur_copy_contigs, cur_copy_path, copy_name_2_copy_dict)
        jobs.append(job)
    ex.shutdown(wait=True)
    intact_ltr_copies = {}
    for job in as_completed(jobs):
        ltr_name, cur_intact_ltr_copies = job.result()
        intact_ltr_copies[ltr_name] = cur_intact_ltr_copies
    return intact_ltr_copies

def get_non_redundant_copies(ltr_name, cur_copy_contigs, cur_copy_path, copy_name_2_copy_dict):
    blastn_command = 'blastn -query ' + cur_copy_path + ' -subject ' + cur_copy_path + ' -num_threads ' + str(1) + ' -outfmt 6 '
    result = subprocess.run(blastn_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            executable='/bin/bash')
    redundant_copies = set()
    if result.returncode == 0:
        lines = result.stdout.split('\n')
        for line in lines:
            parts = line.split('\t')
            if len(parts) != 12:
                continue
            query_name = parts[0]
            subject_name = parts[1]
            if query_name == subject_name:
                continue
            query_len = len(cur_copy_contigs[query_name])
            subject_len = len(cur_copy_contigs[subject_name])

            query_start = int(parts[6])
            query_end = int(parts[7])
            subject_start = int(parts[8])
            subject_end = int(parts[9])
            if abs(query_end - query_start) / query_len >= 0.95 and abs(
                    subject_end - subject_start) / subject_len >= 0.95:
                redundant_copies.add(query_name)
                redundant_copies.add(subject_name)
    # 遍历 cur_copy_contigs，取所有非冗余的拷贝
    intact_copies = []
    for cur_copy_name in cur_copy_contigs.keys():
        if cur_copy_name not in redundant_copies:
            intact_copies.append(copy_name_2_copy_dict[cur_copy_name])
    # 在冗余拷贝中任取一个加入到拷贝中
    if len(redundant_copies) > 0:
        intact_copies.append(copy_name_2_copy_dict[list(redundant_copies)[0]])
    return ltr_name, intact_copies

def filter_ltr_by_copy_num_sub(candidate_sequence_path, threads, temp_dir, split_ref_dir, full_length_threshold, max_copy_num=10):
    debug = 0
    # flanking_len = 100
    if os.path.exists(temp_dir):
        os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # We are considering that the current running time is too long, maybe it is related to submitting one sequence for Blastn alignment at a time.
    # We will try to combine 10 sequences together and run Blastn once.
    # To increase CPU utilization, we will submit one thread to process 10 sequences.
    batch_size = 1
    batch_id = 0
    names, contigs = read_fasta(candidate_sequence_path)
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


    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cur_split_files in split_files:
        job = ex.submit(get_full_length_copies_v1, cur_split_files, split_ref_dir, max_copy_num, full_length_threshold,
                        debug)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_copies = {}
    for job in as_completed(jobs):
        cur_all_copies = job.result()
        all_copies.update(cur_all_copies)

    return all_copies

def get_domain_info(cons, lib, output_table, threads, temp_dir):
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    db_prefix = os.path.basename(lib)
    db_files = [os.path.join(os.path.dirname(lib), f"{db_prefix}.phr"),
                os.path.join(os.path.dirname(lib), f"{db_prefix}.pin"),
                os.path.join(os.path.dirname(lib), f"{db_prefix}.psq")]

    if all(os.path.exists(f) for f in db_files):
        print(f"BLAST database exist, skip creating：{db_prefix}")
    else:
        blast_db_command = f"makeblastdb -dbtype prot -in {lib} > /dev/null 2>&1"
        print(f"Creating BLAST database：{db_prefix}")
        os.system(blast_db_command)

    # 1. blastx -num_threads 1 -evalue 1e-20
    partitions_num = int(threads)
    split_files = split_fasta(cons, temp_dir, partitions_num)

    merge_distance = 300
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for partition_index, cur_consensus_path in enumerate(split_files):
        if not file_exist(cur_consensus_path):
            continue
        cur_output = temp_dir + '/' + str(partition_index) + '.out'
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
        # print('total query size: %d, current query name: %s, idx: %d' % (len(query_records), query_name, idx))

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
                                        2] - longest_subject_end < fixed_extend_base_threshold / 3:
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
                                    if longest_query_end - ext_frag[0] < fixed_extend_base_threshold and ext_frag[
                                        2] - longest_subject_end < fixed_extend_base_threshold / 3:
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
                                    elif longest_query_end - ext_frag[0] >= fixed_extend_base_threshold:
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
                                            cluster_longest_subject_end, cluster_longest_subject_len, subject_name,
                                            cluster_extend_num))

        # we now consider, we should take some sequences from longest_queries to represent this query sequence.
        # we take the longest sequence by length, if the latter sequence overlap with the former sequence largely (50%),
        # continue find next sequence until the ratio of query sequence over 90% or no more sequences.
        longest_queries.sort(key=lambda x: -x[2])
        keep_longest_query[query_name] = longest_queries
    # print(keep_longest_query)
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
                        # 计算overlap
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
                    if is_new_domain:
                        merge_domains.append(domain_info)

            for domain_info in merge_domains:
                f_save.write(query_name + '\t' + str(domain_info[6]) + '\t' + str(domain_info[0]) + '\t' + str(
                    domain_info[1]) + '\t' + str(domain_info[3]) + '\t' + str(domain_info[4]) + '\n')

    f_save.close()
    return cur_table

def get_domain_info_v1(cons, lib, threads, temp_dir):
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    db_prefix = os.path.basename(lib)
    db_files = [os.path.join(os.path.dirname(lib), f"{db_prefix}.phr"),
                os.path.join(os.path.dirname(lib), f"{db_prefix}.pin"),
                os.path.join(os.path.dirname(lib), f"{db_prefix}.psq")]

    if all(os.path.exists(f) for f in db_files):
        print(f"BLAST database exist, skip creating：{db_prefix}")
    else:
        blast_db_command = f"makeblastdb -dbtype prot -in {lib} > /dev/null 2>&1"
        print(f"Creating BLAST database：{db_prefix}")
        os.system(blast_db_command)

    partitions_num = int(threads)
    split_files = split_fasta(cons, temp_dir, partitions_num)
    merge_distance = 100
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for partition_index, cur_consensus_path in enumerate(split_files):
        cur_output = temp_dir + '/'+str(partition_index)+'.out'
        cur_table = temp_dir + '/' + str(partition_index) + '.tbl'
        cur_file = (cur_consensus_path, lib, cur_output, cur_table)
        job = ex.submit(multiple_alignment_blastx_v2, cur_file)
        jobs.append(job)
    ex.shutdown(wait=True)

    query_protein_types = {}
    for job in as_completed(jobs):
        cur_query_protein_types = job.result()
        query_protein_types.update(cur_query_protein_types)
    return query_protein_types

def multiple_alignment_blastx_v2(repeats_path):
    split_repeats_path = repeats_path[0]
    protein_db_path = repeats_path[1]
    blastx2Results_path = repeats_path[2]
    cur_table = repeats_path[3]
    orf_path = split_repeats_path + '.TE.orfs'
    getorf_command = 'getorf -sequence ' + split_repeats_path + ' --outseq ' + orf_path + ' -minsize 400 -reverse ' + ' > /dev/null 2>&1'
    os.system(getorf_command)
    align_command = 'blastp -db ' + protein_db_path + ' -num_threads ' \
                    + str(1) + ' -query ' + orf_path + ' -outfmt 6 | sort -k1,1 -k12,12nr | sort -u -k1,1 > ' + blastx2Results_path
    os.system(align_command)

    query_protein_types = {}
    with open(blastx2Results_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            parts = line.split('\t')
            query_name = parts[0].rpartition('_')[0]
            subject_name = parts[1]
            protein_type = subject_name.split('#')[1]
            if not query_protein_types.__contains__(query_name):
                query_protein_types[query_name] = set()
            protein_types = query_protein_types[query_name]
            protein_types.add(protein_type)
    f_r.close()
    return query_protein_types

def get_domain_info_v2(cons, lib, threads, temp_dir, tool_dir):
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    db_prefix = os.path.basename(lib)
    db_files = [os.path.join(os.path.dirname(lib), f"{db_prefix}.phr"),
                os.path.join(os.path.dirname(lib), f"{db_prefix}.pin"),
                os.path.join(os.path.dirname(lib), f"{db_prefix}.psq")]

    if all(os.path.exists(f) for f in db_files):
        print(f"BLAST database exist, skip creating：{db_prefix}")
    else:
        blast_db_command = f"makeblastdb -dbtype prot -in {lib} > /dev/null 2>&1"
        print(f"Creating BLAST database：{db_prefix}")
        os.system(blast_db_command)

    cons_names, cons_contigs = read_fasta(cons)
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for name in cons_names:
        cur_contigs = {}
        cur_contigs[name] = cons_contigs[name]
        cur_path = temp_dir + '/'+str(name)+'.fa'
        store_fasta(cur_contigs, cur_path)
        cur_file = (name, cur_path, lib)
        job = ex.submit(multiple_alignment_blastx_v3, cur_file, tool_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    query_protein_types = {}
    for job in as_completed(jobs):
        cur_query_protein_types = job.result()
        query_protein_types.update(cur_query_protein_types)
    return query_protein_types

def multiple_alignment_blastx_v3(repeats_path, tool_dir):
    query_name = repeats_path[0]
    split_repeats_path = repeats_path[1]
    protein_db_path = repeats_path[2]

    get_domain_command = 'sh ' + tool_dir + '/get_domain_table.sh ' + split_repeats_path + ' ' + protein_db_path + ' > /dev/null 2>&1'
    os.system(get_domain_command)

    orf_table = split_repeats_path + '.orftetable.clean'
    query_protein_types = {}
    if os.path.exists(orf_table):
        with open(orf_table, 'r') as f_r:
            for idx, line in enumerate(f_r):
                parts = line.split('\t')
                if len(parts) != 9:
                    continue
                subject_name = parts[2]
                protein_type = parts[3]
                orf_start = int(parts[0])
                orf_end = int(parts[1])
                orf_len = abs(orf_end-orf_start)
                protein_len = int(parts[8])
                if float(protein_len)/orf_len < 0.5:
                    continue
                if query_name not in query_protein_types:
                    query_protein_types[query_name] = set()
                protein_types = query_protein_types[query_name]
                protein_types.add(protein_type)
        f_r.close()
    return query_protein_types

def get_recombination_ltr(ltr_candidates, ref_contigs, threads, log):
    log.logger.info('Start get recombination ltr')
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for candidate_index in ltr_candidates.keys():
        (chr_name, left_ltr_start, left_ltr_end, right_ltr_start, right_ltr_end) = ltr_candidates[candidate_index]
        if chr_name not in ref_contigs:
            log.logger.error(
                'Error: Chromosome names in the SCN file do not match the input genome names. Please correct this and rerun.')
            sys.exit(-1)
        ref_seq = ref_contigs[chr_name]
        left_ltr_name = chr_name + ':' + str(left_ltr_start) + '-' + str(left_ltr_end)
        left_ltr_seq = ref_seq[left_ltr_start - 1: left_ltr_end]

        int_ltr_name = chr_name + ':' + str(left_ltr_end) + '-' + str(right_ltr_start)
        int_ltr_seq = ref_seq[left_ltr_end: right_ltr_start - 1]
        job = ex.submit(is_recombination, left_ltr_seq, int_ltr_seq, candidate_index)
        jobs.append(job)
    ex.shutdown(wait=True)

    recombination_candidates = []
    for job in as_completed(jobs):
        cur_is_recombination, cur_candidate_index = job.result()
        if cur_is_recombination:
            recombination_candidates.append(cur_candidate_index)
    return recombination_candidates

def remove_dirty_LTR(confident_lines, log):
    log.logger.info('Start remove dirty LTR')
    new_confident_lines = []
    dirty_lines = []
    dirty_dicts = {}
    for i, cur_line in enumerate(confident_lines):
        parts = cur_line.split(' ')
        cur_chr_name = parts[11]
        cur_left_ltr_start = int(parts[3])
        cur_left_ltr_end = int(parts[4])
        cur_right_ltr_start = int(parts[6])
        cur_right_ltr_end = int(parts[7])
        is_dirty = False
        for j in range(i+1, len(confident_lines)):
            next_line = confident_lines[j]
            parts = next_line.split(' ')
            next_chr_name = parts[11]
            next_left_ltr_start = int(parts[3])
            next_left_ltr_end = int(parts[4])
            next_right_ltr_start = int(parts[6])
            next_right_ltr_end = int(parts[7])

            if cur_chr_name != next_chr_name:
                break

            # 如果 next_left_ltr_start 和 next_right_ltr_end 在 当前LTR序列的内部，则认为当前LTR不干净
            if cur_left_ltr_end < next_left_ltr_start < cur_right_ltr_start and cur_left_ltr_end < next_right_ltr_end < cur_right_ltr_start:
                is_dirty = True
                break
            if next_left_ltr_start > cur_right_ltr_end:
                break

        if not is_dirty:
            new_confident_lines.append(cur_line)
        else:
            dirty_lines.append(cur_line)
            cur_name = cur_chr_name + '-' + str(cur_left_ltr_start) + '-' + str(cur_left_ltr_end)
            dirty_dicts[cur_name] = 1
    # log.logger.debug('Remove dirty LTR: ' + str(len(dirty_lines)) + ', remaining LTR num: ' + str(len(new_confident_lines)))
    # print(dirty_lines)
    return dirty_dicts

def deredundant_for_LTR(redundant_ltr, work_dir, threads, type, coverage_threshold):
    starttime = time.time()
    # We found that performing a direct mafft alignment on the redundant LTR library was too slow.
    # Therefore, we first need to use Blastn for alignment clustering, and then proceed with mafft processing.
    tmp_blast_dir = work_dir + '/LTR_blastn_' + str(type)
    blastnResults_path = work_dir + '/LTR_blastn_' + str(type) + '.out'
    # 1. Start by performing an all-vs-all comparison using blastn.
    multi_process_align(redundant_ltr, redundant_ltr, blastnResults_path, tmp_blast_dir, threads, is_removed_dir=True)
    if not os.path.exists(blastnResults_path):
        return redundant_ltr
    # 2. Next, using the FMEA algorithm, bridge across the gaps and link together sequences that can be connected.
    full_length_threshold = 0.8
    longest_repeats = FMEA_new(redundant_ltr, blastnResults_path, full_length_threshold)

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
            # 我们这里先将跨过 gap 之后的全长拷贝先聚类在一起，后续再使用 cd-hit 将碎片化合并到全长拷贝中
            if float(q_len) / query_len >= coverage_threshold or float(s_len) / subject_len >= coverage_threshold:
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

    # # Iterate through each cluster, if any element in the cluster overlaps with elements in other clusters, merge the clusters.
    # merged_clusters = []
    # while keep_clusters:
    #     current_cluster = keep_clusters.pop(0)
    #     for other_cluster in keep_clusters[:]:
    #         if current_cluster.intersection(other_cluster):
    #             current_cluster.update(other_cluster)
    #             keep_clusters.remove(other_cluster)
    #     merged_clusters.append(current_cluster)
    # keep_clusters = merged_clusters

    endtime = time.time()
    dtime = endtime - starttime
    print("Running time of FMEA clustering: %.8s s" % (dtime))

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
        job = ex.submit(generate_cons, cluster_id, cur_cluster_path, cluster_dir, 1)
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
    print("Running time of MSA cons: %.8s s" % (dtime))

    ltr_cons_cons = redundant_ltr + '.cons'
    # 调用 cd-hit-est 合并碎片化序列
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + ltr_cons_path + ' -o ' + ltr_cons_cons + ' -T 0 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')

    #rename_fasta(ltr_cons_path, ltr_cons_path, 'LTR')
    return ltr_cons_path

def deredundant_for_LTR_v5(redundant_ltr, work_dir, threads, type, coverage_threshold, debug):
    starttime = time.time()
    # We found that performing a direct mafft alignment on the redundant LTR library was too slow.
    # Therefore, we first need to use Blastn for alignment clustering, and then proceed with mafft processing.
    tmp_blast_dir = work_dir + '/LTR_blastn_' + str(type)
    blastnResults_path = work_dir + '/LTR_blastn_' + str(type) + '.out'
    # 1. Start by performing an all-vs-all comparison using blastn.
    multi_process_align(redundant_ltr, redundant_ltr, blastnResults_path, tmp_blast_dir, threads, is_removed_dir=True)
    if not os.path.exists(blastnResults_path):
        return redundant_ltr
    # 2. Next, using the FMEA algorithm, bridge across the gaps and link together sequences that can be connected.
    longest_repeats = FMEA_new(redundant_ltr, blastnResults_path, coverage_threshold)

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
        job = ex.submit(generate_cons, cluster_id, cur_cluster_path, cluster_dir, 1)
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

def deredundant_for_LTR_v2(redundant_ltr, work_dir, threads, type, coverage_threshold):
    # We found that performing a direct mafft alignment on the redundant LTR library was too slow.
    # Therefore, we first need to use Blastn for alignment clustering, and then proceed with mafft processing.
    tmp_blast_dir = work_dir + '/LTR_blastn_' + str(type)
    blastnResults_path = work_dir + '/LTR_blastn_' + str(type) + '.out'
    # 1. Start by performing an all-vs-all comparison using blastn.
    multi_process_align(redundant_ltr, redundant_ltr, blastnResults_path, tmp_blast_dir, threads, is_removed_dir=True)
    if not os.path.exists(blastnResults_path):
        return redundant_ltr
    # 2. Next, using the FMEA algorithm, bridge across the gaps and link together sequences that can be connected.
    longest_repeats = FMEA_new(redundant_ltr, blastnResults_path, coverage_threshold)
    # 3. If the combined sequence length constitutes 95% or more of the original individual sequence lengths, we place these two sequences into a cluster.
    contigNames, contigs = read_fasta(redundant_ltr)
    keep_clusters = []
    for query_name in longest_repeats.keys():
        longest_repeats_list = longest_repeats[query_name]
        cur_cluster = set()
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
                cur_cluster.add(query_name)
                cur_cluster.add(subject_name)
        keep_clusters.append(cur_cluster)
    # store cluster
    all_unique_name = set()
    raw_cluster_files = []
    cluster_dir = work_dir + '/raw_ltr_cluster_' + str(type)
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)
    for cluster_id, cur_cluster in enumerate(keep_clusters):
        if len(cur_cluster) > 1:
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
        job = ex.submit(generate_cons, cluster_id, cur_cluster_path, cluster_dir)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_cons = {}
    for job in as_completed(jobs):
        cur_cons_contigs = job.result()
        all_cons.update(cur_cons_contigs)
    all_cons.update(uncluster_contigs)
    ltr_cons_path = redundant_ltr + '.tmp.cons'
    store_fasta(all_cons, ltr_cons_path)

    ltr_cons_cons = redundant_ltr + '.cons'
    # 调用 cd-hit-est 合并碎片化序列
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + ltr_cons_path + ' -o ' + ltr_cons_cons + ' -T 0 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')

    #rename_fasta(ltr_cons_path, ltr_cons_path, 'LTR')
    return ltr_cons_path

def deredundant_for_LTR_v3(redundant_ltr, work_dir, threads, type, coverage_threshold, split_ref_dir, reference):
    # We found that performing a direct mafft alignment on the redundant LTR library was too slow.
    # Therefore, we first need to use Blastn for alignment clustering, and then proceed with mafft processing.
    tmp_blast_dir = work_dir + '/LTR_blastn_' + str(type)
    blastnResults_path = work_dir + '/LTR_blastn_' + str(type) + '.out'
    # 1. Start by performing an all-vs-all comparison using blastn.
    multi_process_align(redundant_ltr, redundant_ltr, blastnResults_path, tmp_blast_dir, threads, is_removed_dir=True)
    if not os.path.exists(blastnResults_path):
        return redundant_ltr
    # 2. Next, using the FMEA algorithm, bridge across the gaps and link together sequences that can be connected.
    longest_repeats = FMEA_new(redundant_ltr, blastnResults_path, coverage_threshold)
    # 3. If the combined sequence length constitutes 95% or more of the original individual sequence lengths, we place these two sequences into a cluster.
    contigNames, contigs = read_fasta(redundant_ltr)
    represent_ltr_names = set()
    redundant_ltr_names = set()
    for query_name in longest_repeats.keys():
        longest_repeats_list = longest_repeats[query_name]
        for cur_longest_repeat in longest_repeats_list:
            query_name = cur_longest_repeat[0]
            query_len = len(contigs[query_name])
            q_len = abs(cur_longest_repeat[2] - cur_longest_repeat[1])
            subject_name = cur_longest_repeat[3]
            subject_len = len(contigs[subject_name])
            s_len = abs(cur_longest_repeat[4] - cur_longest_repeat[5])
            # 我们这里先将跨过 gap 之后的全长拷贝先聚类在一起，后续再使用 cd-hit 将碎片化合并到全长拷贝中
            if float(q_len) / query_len >= coverage_threshold and float(s_len) / subject_len >= coverage_threshold:
                # we consider the query and subject to be from the same family.
                if query_name not in redundant_ltr_names:
                    represent_ltr_names.add(query_name)
                redundant_ltr_names.add(subject_name)

    # 将在 represent_ltr_names 和 (不在 represent_ltr_names 和 redundant_ltr_names)的ltr存储成文件，然后获取它们的拷贝，进行MSA
    represent_ltr1 = redundant_ltr + '.rep1'
    represent_ltr_contigs = {}
    for ltr_name in contigNames:
        if ltr_name in represent_ltr_names:
            represent_ltr_contigs[ltr_name] = contigs[ltr_name]
        elif ltr_name not in redundant_ltr_names:
            represent_ltr_contigs[ltr_name] = contigs[ltr_name]
    store_fasta(represent_ltr_contigs, represent_ltr1)
    # 调用cd-hit-est 再次压缩，获取非冗余
    represent_ltr2 = represent_ltr1 + '.rep2'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + represent_ltr1 + ' -o ' + represent_ltr2 + ' -T 0 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')
    represent_ltr_names, represent_ltr_contigs = read_fasta(represent_ltr2)

    temp_dir = work_dir + '/ltr_copies_' + str(type)
    max_copy_num = 10
    all_copies = get_full_length_copies_batch_v1(represent_ltr2, split_ref_dir, threads, temp_dir, max_copy_num, coverage_threshold)

    raw_copy_cluster_files = []
    ref_names, ref_contigs = read_fasta(reference)
    # 获取每个序列对应的全长拷贝，放在一个文件中
    copy_cluster_dir = work_dir + '/raw_ltr_copies_cluster_' + str(type)
    no_copy_path = work_dir + '/no_copy_ltr_' + str(type) + '.fa'
    no_copy_contigs = {}
    if not os.path.exists(copy_cluster_dir):
        os.makedirs(copy_cluster_dir)
    for ltr_name in represent_ltr_names:
        cur_copy_cluster_path = copy_cluster_dir + '/' + str(ltr_name) + '.fa'
        cur_copy_cluster_contigs = {}
        if ltr_name in all_copies:
            for copy in all_copies[ltr_name]:
                chr_name = copy[0]
                chr_start = int(copy[1])
                chr_end = int(copy[2])
                copy_name = chr_name + '-' + str(chr_start) + '-' + str(chr_end)
                te_seq = ref_contigs[chr_name][chr_start: chr_end]
                cur_copy_cluster_contigs[copy_name] = te_seq
        else:
            no_copy_contigs[ltr_name] = represent_ltr_contigs[ltr_name]
        store_fasta(cur_copy_cluster_contigs, cur_copy_cluster_path)
        raw_copy_cluster_files.append(cur_copy_cluster_path)
    store_fasta(no_copy_contigs, no_copy_path)

    # 4. The final cluster should encompass all instances from the same family.
    # We use Ninja to cluster families precisely, and
    # We then use the mafft+majority principle to generate a consensus sequence for each cluster.
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cluster_id, cur_copy_cluster_path in enumerate(raw_copy_cluster_files):
        job = ex.submit(generate_cons, cluster_id, cur_copy_cluster_path, copy_cluster_dir)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_cons = {}
    for job in as_completed(jobs):
        cur_cons_contigs = job.result()
        all_cons.update(cur_cons_contigs)
    all_cons.update(no_copy_contigs)
    ltr_cons_path = redundant_ltr + '.tmp.cons'
    store_fasta(all_cons, ltr_cons_path)

    ltr_cons_cons = redundant_ltr + '.cons'
    # 调用 cd-hit-est 合并碎片化序列
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + ltr_cons_path + ' -o ' + ltr_cons_cons + ' -T 0 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')
    return ltr_cons_path

def deredundant_for_LTR_v4(redundant_ltr, work_dir, threads, type, coverage_threshold, split_ref_dir, reference):
    # We found that performing a direct mafft alignment on the redundant LTR library was too slow.
    # Therefore, we first need to use Blastn for alignment clustering, and then proceed with mafft processing.
    tmp_blast_dir = work_dir + '/LTR_blastn_' + str(type)
    blastnResults_path = work_dir + '/LTR_blastn_' + str(type) + '.out'
    # 1. Start by performing an all-vs-all comparison using blastn.
    multi_process_align(redundant_ltr, redundant_ltr, blastnResults_path, tmp_blast_dir, threads, is_removed_dir=True)
    if not os.path.exists(blastnResults_path):
        return redundant_ltr
    # 2. Next, using the FMEA algorithm, bridge across the gaps and link together sequences that can be connected.
    longest_repeats = FMEA_new(redundant_ltr, blastnResults_path, coverage_threshold)
    # 3. If the combined sequence length constitutes 95% or more of the original individual sequence lengths, we place these two sequences into a cluster.
    contigNames, contigs = read_fasta(redundant_ltr)
    represent_ltr_names = {}
    redundant_ltr_names = set()
    for query_name in longest_repeats.keys():
        longest_repeats_list = longest_repeats[query_name]
        for cur_longest_repeat in longest_repeats_list:
            query_name = cur_longest_repeat[0]
            query_len = len(contigs[query_name])
            q_len = abs(cur_longest_repeat[2] - cur_longest_repeat[1])
            subject_name = cur_longest_repeat[3]
            subject_len = len(contigs[subject_name])
            s_len = abs(cur_longest_repeat[4] - cur_longest_repeat[5])
            # 我们这里先将跨过 gap 之后的全长拷贝先聚类在一起，后续再使用 cd-hit 将碎片化合并到全长拷贝中
            if float(q_len) / query_len >= coverage_threshold or float(s_len) / subject_len >= coverage_threshold:
                # 只保留长的那一条当作代表
                if q_len > s_len:
                    represent_ltr_name = query_name
                    redundant_ltr_name = subject_name
                else:
                    represent_ltr_name = subject_name
                    redundant_ltr_name = query_name

                redundant_ltr_names.add(redundant_ltr_name)
                if redundant_ltr_name in represent_ltr_names:
                    del represent_ltr_names[redundant_ltr_name]
                if represent_ltr_name not in redundant_ltr_names:
                    represent_ltr_names[represent_ltr_name] = 1

    # 将在 represent_ltr_names 和 (不在 represent_ltr_names 和 redundant_ltr_names)的ltr存储成文件，然后获取它们的拷贝，进行MSA
    represent_ltr1 = redundant_ltr + '.rep1'
    represent_ltr_contigs = {}
    for ltr_name in contigNames:
        if ltr_name in represent_ltr_names:
            represent_ltr_contigs[ltr_name] = contigs[ltr_name]
        elif ltr_name not in redundant_ltr_names:
            represent_ltr_contigs[ltr_name] = contigs[ltr_name]
    store_fasta(represent_ltr_contigs, represent_ltr1)
    # 调用cd-hit-est 再次压缩，获取非冗余
    represent_ltr2 = represent_ltr1 + '.rep2'
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + represent_ltr1 + ' -o ' + represent_ltr2 + ' -T 0 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')
    represent_ltr_names, represent_ltr_contigs = read_fasta(represent_ltr2)

    temp_dir = work_dir + '/ltr_copies_' + str(type)
    max_copy_num = 10
    all_copies = get_full_length_copies_batch_v1(represent_ltr2, split_ref_dir, threads, temp_dir, max_copy_num, coverage_threshold)

    raw_copy_cluster_files = []
    ref_names, ref_contigs = read_fasta(reference)
    # 获取每个序列对应的全长拷贝，放在一个文件中
    copy_cluster_dir = work_dir + '/raw_ltr_copies_cluster_' + str(type)
    no_copy_path = work_dir + '/no_copy_ltr_' + str(type) + '.fa'
    no_copy_contigs = {}
    if not os.path.exists(copy_cluster_dir):
        os.makedirs(copy_cluster_dir)
    for ltr_name in represent_ltr_names:
        cur_copy_cluster_path = copy_cluster_dir + '/' + str(ltr_name) + '.fa'
        cur_copy_cluster_contigs = {}
        if ltr_name in all_copies:
            for copy in all_copies[ltr_name]:
                chr_name = copy[0]
                chr_start = int(copy[1])
                chr_end = int(copy[2])
                copy_name = chr_name + '-' + str(chr_start) + '-' + str(chr_end)
                te_seq = ref_contigs[chr_name][chr_start: chr_end]
                cur_copy_cluster_contigs[copy_name] = te_seq
        else:
            no_copy_contigs[ltr_name] = represent_ltr_contigs[ltr_name]
        store_fasta(cur_copy_cluster_contigs, cur_copy_cluster_path)
        raw_copy_cluster_files.append(cur_copy_cluster_path)
    store_fasta(no_copy_contigs, no_copy_path)

    # 4. The final cluster should encompass all instances from the same family.
    # We use Ninja to cluster families precisely, and
    # We then use the mafft+majority principle to generate a consensus sequence for each cluster.
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cluster_id, cur_copy_cluster_path in enumerate(raw_copy_cluster_files):
        job = ex.submit(generate_cons, cluster_id, cur_copy_cluster_path, copy_cluster_dir)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_cons = {}
    for job in as_completed(jobs):
        cur_cons_contigs = job.result()
        all_cons.update(cur_cons_contigs)
    all_cons.update(no_copy_contigs)
    ltr_cons_path = redundant_ltr + '.tmp.cons'
    store_fasta(all_cons, ltr_cons_path)

    ltr_cons_cons = redundant_ltr + '.cons'
    # 调用 cd-hit-est 合并碎片化序列
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + ltr_cons_path + ' -o ' + ltr_cons_cons + ' -T 0 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')
    return ltr_cons_path

def deredundant_for_LTR_v1(redundant_ltr, work_dir, reference, split_ref_dir, threads, coverage_threshold, type):
    # We found that performing a direct mafft alignment on the redundant LTR library was too slow.
    # Therefore, we first need to use Blastn for alignment clustering, and then proceed with mafft processing.
    tmp_blast_dir = work_dir + '/LTR_blastn_' + str(type)
    blastnResults_path = work_dir + '/LTR_blastn_' + str(type) + '.out'
    # 1. Start by performing an all-vs-all comparison using blastn.
    multi_process_align(redundant_ltr, redundant_ltr, blastnResults_path, tmp_blast_dir, threads, is_removed_dir=True)
    if not os.path.exists(blastnResults_path):
        return redundant_ltr
    # 2. Next, using the FMEA algorithm, bridge across the gaps and link together sequences that can be connected.
    full_length_threshold = coverage_threshold
    longest_repeats = FMEA_new(redundant_ltr, blastnResults_path, full_length_threshold)
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
            # 我们这里先将跨过 gap 之后的全长拷贝先聚类在一起
            if float(q_len) / query_len >= coverage_threshold or float(s_len) / subject_len >= coverage_threshold:
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
    cluster_dir = work_dir + '/raw_ltr_cluster_' + str(type)
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

    # 获取冗余LTR的所有全长拷贝
    temp_dir = work_dir + '/redundant_ltr_copy_' + str(type)
    max_copy_num = 10
    all_copies = get_full_length_copies_batch_v1(redundant_ltr, split_ref_dir, threads, temp_dir, max_copy_num, coverage_threshold)

    raw_copy_cluster_files = []
    ref_names, ref_contigs = read_fasta(reference)
    # 获取每个簇序列对应的全长拷贝，放在一个文件中
    copy_cluster_dir = work_dir + '/raw_ltr_copies_cluster_' + str(type)
    if not os.path.exists(copy_cluster_dir):
        os.makedirs(copy_cluster_dir)
    for cluster_id, cur_cluster_path in raw_cluster_files:
        cluster_names, cluster_contigs = read_fasta(cur_cluster_path)
        cur_copy_cluster_path = copy_cluster_dir + '/' + str(cluster_id) + '.fa'
        cur_copy_cluster_contigs = {}
        deredundant_copy_set = {}
        for query_name in cluster_names:
            if query_name in all_copies:
                for copy in all_copies[query_name]:
                    chr_name = copy[0]
                    chr_start = int(copy[1])
                    chr_end = int(copy[2])
                    is_redundant = False
                    # 判断当前拷贝是否已存储
                    if chr_name not in deredundant_copy_set:
                        deredundant_copy_set[chr_name] = []
                    copy_list = deredundant_copy_set[chr_name]
                    cur_copy = (chr_start, chr_end)
                    for exist_copy in copy_list:
                        overlap_len = get_overlap_len(exist_copy, cur_copy)
                        if overlap_len / abs(exist_copy[1] - exist_copy[0]) >= coverage_threshold and overlap_len / abs(cur_copy[1] - cur_copy[0]) >= coverage_threshold:
                            is_redundant = True
                            break
                    if not is_redundant:
                        copy_list.append(cur_copy)
                        deredundant_copy_set[chr_name] = copy_list

                        copy_name = chr_name + '-' + str(chr_start) + '-' + str(chr_end)
                        te_seq = ref_contigs[chr_name][chr_start: chr_end]
                        cur_copy_cluster_contigs[copy_name] = te_seq
            else:
                # 对于未能找到全长拷贝的序列，我们使用header中所代表的序列表示
                parts = query_name.split('#')[0].split('-')
                chr_name = parts[0]
                chr_start = int(parts[1])
                chr_end = int(parts[2])
                is_redundant = False
                # 判断当前拷贝是否已存储
                if chr_name not in deredundant_copy_set:
                    deredundant_copy_set[chr_name] = []
                copy_list = deredundant_copy_set[chr_name]
                cur_copy = (chr_start, chr_end)
                for exist_copy in copy_list:
                    overlap_len = get_overlap_len(exist_copy, cur_copy)
                    if overlap_len / abs(exist_copy[1] - exist_copy[0]) >= coverage_threshold and overlap_len / abs(
                            cur_copy[1] - cur_copy[0]) >= coverage_threshold:
                        is_redundant = True
                        break
                if not is_redundant:
                    copy_list.append(cur_copy)
                    deredundant_copy_set[chr_name] = copy_list

                    copy_name = chr_name + '-' + str(chr_start) + '-' + str(chr_end)
                    te_seq = ref_contigs[chr_name][chr_start: chr_end]
                    cur_copy_cluster_contigs[copy_name] = te_seq
        store_fasta(cur_copy_cluster_contigs, cur_copy_cluster_path)
        raw_copy_cluster_files.append((cluster_id, cur_copy_cluster_path, cluster_contigs))


    # 4. The final cluster should encompass all instances from the same family.
    # We use Ninja to cluster families precisely, and
    # We then use the mafft+majority principle to generate a consensus sequence for each cluster.
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cluster_id, cur_copy_cluster_path, cluster_contigs in raw_copy_cluster_files:
        job = ex.submit(generate_cons_v2, cluster_id, cur_copy_cluster_path, copy_cluster_dir, cluster_contigs)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_cons = {}
    for job in as_completed(jobs):
        cur_cons_contigs = job.result()
        all_cons.update(cur_cons_contigs)
    all_cons.update(uncluster_contigs)
    ltr_cons_path = redundant_ltr + '.tmp.cons'
    store_fasta(all_cons, ltr_cons_path)

    ltr_cons_cons = redundant_ltr + '.cons'
    # 调用 cd-hit-est 合并碎片化序列
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + ltr_cons_path + ' -o ' + ltr_cons_cons + ' -T 0 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')

    #rename_fasta(ltr_cons_path, ltr_cons_path, 'LTR')
    return ltr_cons_path

def generate_cons_v2(cluster_id, cur_cluster_path, cluster_dir, cluster_contigs):
    ltr_internal_names, ltr_internal_contigs = read_fasta(cur_cluster_path)
    temp_cluster_dir = cluster_dir
    cons_contigs = {}
    # 使用大多数规则，我们认为只有超过5个拷贝支持的簇，才能恢复碎片LTR内部序列。如果没有找到，则使用原始序列
    is_found_cons = False
    if len(ltr_internal_contigs) >= 1:
        align_file = cur_cluster_path + '.maf.fa'
        align_command = 'cd ' + cluster_dir + ' && mafft --preservecase --quiet --thread -1 ' + cur_cluster_path + ' > ' + align_file
        os.system(align_command)

        # 调用 Ninja 对多序列比对再次聚类
        cluster_file = align_file + '.dat'
        Ninja_command = 'Ninja --in ' + align_file + ' --out ' + cluster_file + ' --out_type c --corr_type m --cluster_cutoff 0.2 --threads 1'
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
                seq = ltr_internal_contigs[name]
                cur_cluster_contigs[name] = seq
                cur_ltr_name = name
            store_fasta(cur_cluster_contigs, cur_cluster_file)

            cur_align_file = cur_cluster_file + '.maf.fa'
            if len(cur_cluster_contigs) >= 5:
                align_command = 'cd ' + Ninja_cluster_dir + ' && mafft --preservecase --quiet --thread -1 ' + cur_cluster_file + ' > ' + cur_align_file
                os.system(align_command)
                cons_seq = cons_from_mafft(cur_align_file)
                cons_contigs[cur_ltr_name] = cons_seq
                is_found_cons = True
    if is_found_cons:
        return cons_contigs
    else:
        return cluster_contigs

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
        # else:
        #     # Here, we do not use 'N' because it can make it difficult to find the boundary. Therefore, we take the base with the highest non-empty count.
        #     max_base_count = 0
        #     max_base = ''
        #     for cur_base in base_map.keys():
        #         if cur_base == '-':
        #             continue
        #         cur_count = base_map[cur_base]
        #         if cur_count > max_base_count:
        #             max_base_count = cur_count
        #             max_base = cur_base
        #     model_seq += max_base
    return model_seq


def generate_cons(cluster_id, cur_cluster_path, cluster_dir, threads):
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
                cons_seq = cons_from_mafft(cur_align_file)
                cons_contigs[cur_ltr_name] = cons_seq
    # 如果未能识别到可靠的一致性序列，则使用原始序列代替
    if len(cons_contigs) > 0:
        return cons_contigs
    else:
        return ltr_terminal_contigs

def generate_cons_v1(cluster_id, cur_cluster_path, cluster_dir, threads):
    ltr_terminal_names, ltr_terminal_contigs = read_fasta(cur_cluster_path)
    temp_cluster_dir = cluster_dir
    cons_contigs = {}
    if len(ltr_terminal_contigs) >= 1:
        align_file = cur_cluster_path + '.maf.fa'
        align_command = 'cd ' + cluster_dir + ' && famsa -t ' + str(threads) + ' ' + cur_cluster_path + ' ' + align_file + ' > /dev/null 2>&1'
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
                cons_seq = cons_from_mafft(cur_align_file)
                cons_contigs[cur_ltr_name] = cons_seq
    # 如果未能识别到可靠的一致性序列，则使用原始序列代替
    if len(cons_contigs) > 0:
        return cons_contigs
    else:
        return ltr_terminal_contigs

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

def FMEA_LTR(blastn2Results_path, fixed_extend_base_threshold):
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
            if query_name == subject_name:
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

def FMEA_new(query_path, blastn2Results_path, full_length_threshold):
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

    query_names, query_contigs = read_fasta(query_path)

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

def multi_process_align_blastx(query_path, subject_path, blastxResults_path, tmp_blast_dir, threads, is_removed_dir=True):
    tools_dir = ''
    if is_removed_dir:
        os.system('rm -rf ' + tmp_blast_dir)
    if not os.path.exists(tmp_blast_dir):
        os.makedirs(tmp_blast_dir)

    if os.path.exists(blastxResults_path):
        os.remove(blastxResults_path)

    orig_names, orig_contigs = read_fasta(query_path)

    blast_db_command = 'makeblastdb -dbtype prot -in ' + subject_path + ' > /dev/null 2>&1'
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
        job = ex.submit(multiple_alignment_blastx, file, tools_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        cur_blastn2Results_path = job.result()
        os.system('cat ' + cur_blastn2Results_path + ' >> ' + blastxResults_path)

def multiple_alignment_blast(repeats_path, tools_dir):
    split_repeats_path = repeats_path[0]
    ref_db_path = repeats_path[1]
    blastn2Results_path = repeats_path[2]

    align_command = 'blastn -db ' + ref_db_path + ' -num_threads ' \
                    + str(1) + ' -query ' + split_repeats_path + ' -evalue 1e-20 -outfmt 6 > ' + blastn2Results_path
    os.system(align_command)

    return blastn2Results_path

def multiple_alignment_blastx(repeats_path, tools_dir):
    split_repeats_path = repeats_path[0]
    ref_db_path = repeats_path[1]
    blastn2Results_path = repeats_path[2]

    align_command = 'blastx -word_size 3 -max_target_seqs 10 -db ' + ref_db_path + ' -num_threads ' \
                    + str(1) + ' -query ' + split_repeats_path + ' -evalue 1e-3 -outfmt 6 > ' + blastn2Results_path
    os.system(align_command)

    return blastn2Results_path

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

def convert_LtrDetector_scn(LtrDetector_output, scn_file):
    with open(scn_file, 'w') as f_save:
        f_save.write('# LtrDetector \n')
        with open(LtrDetector_output, 'r') as f_r:
            for i, line in enumerate(f_r):
                line = str(line).replace('\n', '')
                parts = line.split('\t')
                # 由于多线程运行的LtrDetector结果会出现多个header，因此我们排除掉这些无用的信息
                if len(parts) != 18 or parts[0] == '':
                    continue
                ltr_start = int(parts[1]) + 1
                ltr_end = int(parts[2])
                ltr_len = ltr_end - ltr_start + 1
                left_ltr_start = int(parts[3]) + 1
                left_ltr_end = int(parts[4])
                left_ltr_len = left_ltr_end - left_ltr_start + 1
                right_ltr_start = int(parts[5]) + 1
                right_ltr_end = int(parts[6])
                right_ltr_len = right_ltr_end - right_ltr_start + 1
                ltr_identity = float(parts[7])
                seq_id = 'NA'
                chr_name = parts[0]
                new_line = str(ltr_start) + ' ' + str(ltr_end) + ' ' + str(ltr_len) + ' ' + \
                           str(left_ltr_start) + ' ' + str(left_ltr_end) + ' ' + str(left_ltr_len) + ' ' + \
                           str(right_ltr_start) + ' ' + str(right_ltr_end) + ' ' + str(right_ltr_len) + ' ' + \
                           str(ltr_identity) + ' ' + str(seq_id) + ' ' + chr_name + '\n'
                f_save.write(new_line)

def judge_both_ends_frame_v1(maxtrix_file, debug=1):
    # 我现在想的假阳性过滤方法：
    # 1. 对matrix file 搜索同源边界，如果不存在，则说明是真实LTR，否则为假阳性

    is_ltr = True
    seq_name = os.path.basename(maxtrix_file).split('.')[0]
    # Step3. 对候选随机序列 搜索同源边界。我们将窗口设置为40.
    is_left_ltr, new_boundary_start = judge_left_frame_LTR(maxtrix_file)
    if debug:
        print(maxtrix_file, is_left_ltr, new_boundary_start)
    is_right_ltr, new_boundary_end = judge_right_frame_LTR(maxtrix_file)
    if debug:
        print(maxtrix_file, is_right_ltr, new_boundary_end)
    is_ltr &= is_left_ltr and is_right_ltr
    return seq_name, is_ltr

def judge_right_frame_LTR(matrix_file, flanking_len, sliding_window_size=20):
    pos = 0
    debug = 0
    col_num = -1
    row_num = 0
    lines = []
    no_empty_row = 0
    with open(matrix_file, 'r') as f_r:
        for line in f_r:
            parts = line.replace('\n', '').split('\t')
            line = parts[1]
            col_num = len(line)
            row_num += 1
            lines.append(line)
            if line != '-' * col_num:
                no_empty_row += 1

    # 过滤掉单拷贝的LTR，因为我们没办法判断它是否是LTR
    if row_num <= 1:
        return True, -1

    matrix = [[''] * col_num for i in range(row_num)]
    for row, line in enumerate(lines):
        for col in range(len(line)):
            matrix[row][col] = line[col]

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

    search_len = flanking_len
    valid_col_threshold = int(row_num / 2)

    if row_num <= 5:
        homo_threshold = 0.95
    elif row_num <= 10:
        homo_threshold = 0.9
    elif row_num <= 50:
        homo_threshold = 0.85
    else:
        homo_threshold = 0.85
    # homo_threshold = 0.9

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
        if no_gap_num <= 1:
            col_index += 1
            continue
        max_homo_ratio = 0
        gap_num = base_map['-']
        # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
        if gap_num <= valid_col_threshold:
            valid_col_count += 1
            # Determine if the effective column is homologous.
            for base in base_map.keys():
                if base == '-':
                    continue
                # 修正bug，row_num 替换成 no_gap_num
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
            if cur_homo_ratio >= homo_threshold - 0.1 and first_candidate_boundary == -1:
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

    if new_boundary_end != -1 and abs(new_boundary_end - pos) > 20:
        return False, new_boundary_end
    else:
        return True, new_boundary_end

def judge_left_frame_LTR(matrix_file, flanking_len, sliding_window_size=20):
    pos = flanking_len - 1
    debug = 0
    col_num = -1
    row_num = 0
    lines = []
    no_empty_row = 0
    with open(matrix_file, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '').split('\t')[0]
            col_num = len(line)
            row_num += 1
            lines.append(line)
            if line != '-' * col_num:
                no_empty_row += 1

    # 过滤掉单拷贝的LTR，因为我们没办法判断它是否是LTR
    if row_num <= 1:
        return True, -1

    matrix = [[''] * col_num for i in range(row_num)]
    for row, line in enumerate(lines):
        for col in range(len(line)):
            matrix[row][col] = line[col]

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

    search_len = flanking_len
    valid_col_threshold = int(row_num / 2)

    if row_num <= 5:
        homo_threshold = 0.95
    elif row_num <= 10:
        homo_threshold = 0.9
    elif row_num <= 50:
        homo_threshold = 0.85
    else:
        homo_threshold = 0.85

    # homo_threshold = 0.85

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
    while valid_col_count < search_len and col_index >= 0:
        # Starting from position 'pos', search for 15 effective columns to the left.
        # Determine if the current column is effective.
        is_homo_col = False
        base_map = col_base_map[col_index]
        max_homo_ratio = 0
        no_gap_num = row_num - base_map['-']
        if no_gap_num <= 1:
            col_index -= 1
            continue
        gap_num = base_map['-']
        # If the number of gaps in the current column is <= half of the copy count, then it is an effective column.
        if gap_num <= valid_col_threshold:
            valid_col_count += 1
            # Determine if the effective column is homologous.
            for base in base_map.keys():
                if base == '-':
                    continue
                # 修正bug，row_num 替换成 no_gap_num
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
    if new_boundary_start != pos and new_boundary_start != -1:
        if debug:
            print('align start left homology, new boundary: ' + str(new_boundary_start))
        cur_boundary = new_boundary_start

    if new_boundary_start != -1 and abs(new_boundary_start - pos) > 5:
        return False, new_boundary_start
    else:
        return True, new_boundary_start

def filter_ltr_by_homo(dl_output_path, homo_output_path, matrix_dir, threads):
    true_ltr_names = []
    ltr_dict = {}
    with open(dl_output_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ltr_name = parts[0]
            is_ltr = int(parts[1])
            ltr_dict[ltr_name] = is_ltr
            if is_ltr:
                true_ltr_names.append(ltr_name)

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for ltr_name in true_ltr_names:
        cur_matrix_file = matrix_dir + '/' + ltr_name + '.matrix'
        job = ex.submit(judge_both_ends_frame_v1, cur_matrix_file, debug=0)
        jobs.append(job)
    ex.shutdown(wait=True)
    filter_ltr_names = []
    for job in as_completed(jobs):
        cur_seq_name, cur_is_ltr = job.result()
        if cur_is_ltr:
            cur_is_ltr = 1
        else:
            cur_is_ltr = 0
            filter_ltr_names.append(cur_seq_name)
        ltr_dict[cur_seq_name] = cur_is_ltr
    print('Deep Learning LTR num: ' + str(len(true_ltr_names)) + ', Homology filter LTR num: ' + str(len(filter_ltr_names)))
    print(filter_ltr_names)
    with open(homo_output_path, 'w') as f_save:
        for ltr_name in ltr_dict.keys():
            f_save.write(ltr_name+'\t'+str(ltr_dict[ltr_name])+'\n')


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

def get_non_empty_seq(raw_align_seq):
    align_seq = ''
    for i, base in enumerate(raw_align_seq):
        if base == '-':
            continue
        else:
            align_seq += base
    return align_seq


def is_rich_in_ta(sequence, threshold=0.8):
    """
    判断给定的DNA序列是否富含TA
    :param sequence: DNA序列，只包含ATCG
    :param threshold: 富集阈值，默认为0.8（即TA的比例超过80%）
    :return: 如果序列富含TA，返回True；否则返回False
    """
    # 统计TA的数量
    ta_count = sequence.count('T') + sequence.count('A')

    # 计算TA的比例
    ta_ratio = ta_count / len(sequence)

    # 判断是否富含TA
    return ta_ratio > threshold

# 判断窗口的拷贝是否具有TSD特征
def is_TIR_frame(matrix_file, ltr_name, debug=1):
    # 我们只过滤 4-6bp 以外的tsd
    TIR_TSDs = [15, 14, 13, 12, 11, 10, 9, 8, 7, 3, 2]
    has_tsd_copy_count = 0
    copy_count = 0
    with open(matrix_file, 'r') as f_r:
        for i, line in enumerate(f_r):
            parts = line.split('\t')
            left_line = parts[0]
            right_line = parts[1]
            copy_count += 1
            left_non_empty_seq = get_non_empty_seq(left_line)
            right_non_empty_seq = get_non_empty_seq(right_line)
            for tsd_len in TIR_TSDs:
                left_tsd = left_non_empty_seq[-tsd_len:]
                right_tsd = right_non_empty_seq[:tsd_len]
                if len(left_tsd) != tsd_len or len(right_tsd) != tsd_len:
                    continue
                allow_mismatch_num = 0
                if len(left_tsd) >= 8:
                    allow_mismatch_num = 1
                if allow_mismatch(left_tsd, right_tsd, allow_mismatch_num):
                    # 规定了几种特定的短TSD
                    if (tsd_len == 2 and left_tsd != 'TA') \
                            or (tsd_len == 3 and (left_tsd != 'TTA' or left_tsd != 'TAA')) \
                            or (tsd_len == 4 and left_tsd != 'TTAA'):
                        continue
                    # if is_rich_in_ta(left_tsd):
                    #     continue
                    has_tsd_copy_count += 1
                    # print(i, left_tsd, right_tsd)
                    break
    # print(has_tsd_copy_count)
    is_TIR = False
    if has_tsd_copy_count >= 10:
        is_TIR = True
    return ltr_name, is_TIR


def get_confident_TIR(candidate_tir_path, tool_dir):
    output_dir = os.path.dirname(candidate_tir_path)
    itrsearch_command = 'cd ' + output_dir + ' && ' + tool_dir + '/itrsearch -i 0.7 -l 7 ' + candidate_tir_path + ' > /dev/null 2>&1'
    run_command(itrsearch_command)
    tir_file = candidate_tir_path + '.itr'
    tir_names, tir_contigs = read_fasta_v1(tir_file)
    TIR_info = {}
    for tir_name in tir_names:
        parts = tir_name.split('\t')
        orig_name = parts[0]
        terminal_info = parts[-1]
        TIR_info_parts = terminal_info.split('ITR')[1].split(' ')[0].replace('(', '').replace(')', '').split('..')
        TIR_left_pos_parts = TIR_info_parts[0].split(',')
        TIR_right_pos_parts = TIR_info_parts[1].split(',')
        lTIR_start = int(TIR_left_pos_parts[0])
        lTIR_end = int(TIR_left_pos_parts[1])
        rTIR_start = int(TIR_right_pos_parts[1])
        rTIR_end = int(TIR_right_pos_parts[0])
        TIR_info[orig_name.split(' ')[0]] = (lTIR_start, lTIR_end, rTIR_start, rTIR_end)
    return TIR_info

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
    HelitronScanner_command = 'cd ' + temp_dir + ' && ' + 'sh ' + sh_dir + '/run_helitron_scanner.sh ' \
                              + str(partition_index) + ' ' + cur_candidate_Helitrons_path + ' ' + HSDIR + ' ' + HSJAR + ' > /dev/null 2>&1'
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
    EAHelitron_command = 'cd ' + temp_dir + ' && ' + 'perl ' + EAHelitron + '/EAHelitron -o ' + str(partition_index) + ' -u 20000 -T "ATC" -r 3 ' + all_candidate_helitron_path + '> /dev/null 2>&1'
    os.system(EAHelitron_command)

    all_EAHelitron_res = temp_dir + '/' + str(partition_index) + '.5.fa'
    all_copies_out_names, all_copies_out_contigs = read_fasta_v1(all_EAHelitron_res)

    candidate_Helitrons = {}
    # 如果识别到的边界离原始边界过远，则丢弃
    for cur_name in all_copies_out_names:
        raw_name = cur_name.split(' ')[1]
        parts = raw_name.split(':')
        raw_name = ':'.join(parts[:-1])
        pos_parts = parts[-1].split('..')
        cur_start = int(pos_parts[0])
        cur_end = int(pos_parts[1])
        if cur_start > cur_end:
            tmp = cur_start
            cur_start = cur_end
            cur_end = tmp
        cur_seq = all_candidate_helitron_contigs[raw_name]

        out_threshold = 15
        if cur_start > out_threshold or len(cur_seq) - cur_end > out_threshold:
            continue
        candidate_Helitrons[raw_name] = all_copies_out_contigs[cur_name]

    return candidate_Helitrons

def multi_process_helitronscanner(candidate_file, output, sh_dir, temp_dir, HSDIR, HSJAR, threads, debug):
    os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # After partitioning the files, perform parallel computations using multiple processes.
    fasta_file = candidate_file
    subfile_size = 50000  # 50K
    output_dir = temp_dir
    split_fasta_v1(fasta_file, subfile_size, output_dir)
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

def split_fasta_v1(fasta_file, subfile_size, output_dir):
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

def multi_process_EAHelitron(longest_repeats_flanked_path, flanking_len, output, temp_dir, EAHelitron, threads):
    os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # After partitioning the files, perform parallel computations using multiple processes.
    fasta_file = longest_repeats_flanked_path
    subfile_size = 50000  # 50K
    output_dir = temp_dir
    split_fasta_v1(fasta_file, subfile_size, output_dir)
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

def get_confident_Helitron(candidate_helitron_path, confident_helitron_path, project_dir, HS_temp_dir, EA_temp_dir, flanking_len, threads):
    sh_dir = project_dir + '/bin'
    HSDIR = project_dir + '/bin/HelitronScanner/TrainingSet'
    HSJAR = project_dir + '/bin/HelitronScanner/HelitronScanner.jar'
    debug = 0
    # run HelitronScanner
    multi_process_helitronscanner(candidate_helitron_path, confident_helitron_path, sh_dir, HS_temp_dir,
                                  HSDIR, HSJAR, threads, debug)

    helitron_names, helitron_contigs = read_fasta(confident_helitron_path)
    Helitron_info = {}
    for helitron_name in helitron_names:
        parts = helitron_name.split('_#')
        orig_name = parts[0]
        Helitron_info[orig_name] = 1
    # print('HelitronScanner:' + str(Helitron_info)+ ', ' + str(len(Helitron_info)))

    # run EAHelitron
    EAHelitron = project_dir + '/bin/EAHelitron-master'
    multi_process_EAHelitron(candidate_helitron_path, flanking_len, confident_helitron_path, EA_temp_dir,
                             EAHelitron, threads)
    helitron_names, helitron_contigs = read_fasta(confident_helitron_path)
    for helitron_name in helitron_names:
        parts = helitron_name.split('-rawstart_')
        orig_name = parts[0]
        Helitron_info[orig_name] = 1
    # print('EAHelitron:' + str(Helitron_info) + ', ' + str(len(Helitron_info)))
    return Helitron_info

def search_ltr_structure(ltr_name, left_seq, right_seq):
    # # 搜索左右两侧是否存在TG...CA 或 TSD
    # has_tgca = False
    # if 'TG' in left_seq and 'CA' in right_seq:
    #     has_tgca = True

    has_tsd = False
    tsd_lens = [6, 5, 4]
    tsd_seq = ''
    exist_tsd = set()
    for k_num in tsd_lens:
        for i in range(len(left_seq) - k_num + 1):
            left_kmer = left_seq[i: i + k_num]
            if 'N' not in left_kmer:
                exist_tsd.add(left_kmer)
    for k_num in tsd_lens:
        if has_tsd:
            break
        for i in range(len(right_seq) - k_num + 1):
            right_kmer = right_seq[i: i + k_num]
            if 'N' not in right_kmer and right_kmer in exist_tsd:
                has_tsd = True
                tsd_seq = right_kmer
                break
    # print(ltr_name, left_seq, right_seq, has_tsd, tsd_seq)

    return ltr_name, has_tsd, tsd_seq

def has_consecutive_repeats_near_tail(sequence, min_repeat_length=2, min_repeats=3, max_distance_from_tail=5):
    """
    检查序列尾部是否存在至少 min_repeats 次串联重复子序列，
    并且这些重复子序列离序列尾部的距离不超过 max_distance_from_tail 个碱基对。

    :param sequence: 要检查的序列（字符串或列表）
    :param min_repeat_length: 最小重复单元长度，默认为2
    :param min_repeats: 至少重复次数，默认为3
    :param max_distance_from_tail: 串联重复离尾部的最大距离，默认为5
    :return: 如果尾部存在至少 min_repeats 次串联重复且距离不超过 max_distance_from_tail 返回True，否则返回False
    """
    sequence_length = len(sequence)

    # 检查尾部区域
    for repeat_length in range(min_repeat_length, sequence_length // min_repeats + 1):
        for i in range(max(0, sequence_length - repeat_length * min_repeats - max_distance_from_tail),
                       sequence_length - repeat_length * min_repeats + 1):
            repeat = sequence[i:i + repeat_length]
            match = True
            for j in range(1, min_repeats):
                if sequence[i + j * repeat_length:i + (j + 1) * repeat_length] != repeat:
                    match = False
                    break
            if match:
                return True
    return False


def has_polyA_or_polyT_near_tail(sequence, poly_length=6, max_mismatches=1, max_distance_from_tail=5):
    """
    检查序列尾部是否存在 polyA 或 polyT 序列，并且这些序列离序列尾部的距离不超过 max_distance_from_tail，
    容许最多 max_mismatches 个错配。

    :param sequence: 要检查的序列（字符串）
    :param poly_length: polyA 或 polyT 序列的长度，默认为5
    :param max_mismatches: 允许的最大错配数，默认为1
    :param max_distance_from_tail: polyA 或 polyT 离尾部的最大距离，默认为5
    :return: 如果尾部存在符合条件的 polyA 或 polyT 序列返回True，否则返回False
    """
    sequence_length = len(sequence)
    tail_start = max(0, sequence_length - poly_length - max_distance_from_tail)
    polyA = 'A' * poly_length
    polyT = 'T' * poly_length

    for i in range(tail_start, sequence_length - poly_length + 1):
        segment = sequence[i:i + poly_length]
        mismatches_A = sum(1 for a, b in zip(segment, polyA) if a != b)
        mismatches_T = sum(1 for a, b in zip(segment, polyT) if a != b)
        if mismatches_A <= max_mismatches or mismatches_T <= max_mismatches:
            return True
    return False

# 识别 polyA 和 simple repeat
def identify_SINE_tail(sequence, tail_length=20):
    tail_sequence = sequence[-tail_length:]
    if has_polyA_or_polyT_near_tail(tail_sequence) or has_consecutive_repeats_near_tail(tail_sequence):
        return True
    else:
        return False

def filter_ltr_by_structure(output_path, structure_output_path, leftLtr2Candidates, ltr_lines, reference, threads, log):
    ref_names, ref_contigs = read_fasta(reference)

    true_ltrs = {}
    ltr_dict = {}
    with open(output_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ltr_name = parts[0]
            is_ltr = int(parts[1])
            if is_ltr:
                true_ltrs[ltr_name] = is_ltr
            ltr_dict[ltr_name] = is_ltr

    confident_lines = []
    for name in leftLtr2Candidates.keys():
        # if name not in FP_ltrs:
        if name in true_ltrs:
            candidate_index = leftLtr2Candidates[name]
            line = ltr_lines[candidate_index]
            confident_lines.append((name, line))

    left_size = 8
    internal_size = 3
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for ltr_name, line in confident_lines:
        parts = line.split(' ')
        chr_name = parts[11]
        ref_seq = ref_contigs[chr_name]

        lLTR_start = int(parts[3])
        lLTR_end = int(parts[4])
        rLTR_start = int(parts[6])
        rLTR_end = int(parts[7])

        # 取左/右侧 8bp + 3bp
        # 计算左LTR的切片索引，并确保它们在范围内
        left_start = max(lLTR_start - 1 - left_size, 0)
        left_end = min(lLTR_start + internal_size - 1, len(ref_seq))
        left_seq = ref_seq[left_start: left_end]

        # 计算右LTR的切片索引，并确保它们在范围内
        right_start = max(rLTR_end - internal_size, 0)
        right_end = min(rLTR_end + left_size, len(ref_seq))
        right_seq = ref_seq[right_start: right_end]

        job = ex.submit(search_ltr_structure, ltr_name, left_seq, right_seq)
        jobs.append(job)
    ex.shutdown(wait=True)

    FP_ltrs = {}
    for job in as_completed(jobs):
        cur_seq_name, cur_is_tp = job.result()
        if not cur_is_tp:
            FP_ltrs[cur_seq_name] = 1

    log.logger.debug('LTR num: ' + str(len(true_ltrs)) + ', LTR structure filter num: ' + str(len(FP_ltrs)) + ', remaining LTR num: ' + str(len(true_ltrs) - len(FP_ltrs)))
    with open(structure_output_path, 'w') as f_save:
        for ltr_name in ltr_dict.keys():
            if ltr_name in FP_ltrs:
                ltr_dict[ltr_name] = 0
            f_save.write(ltr_name + '\t' + str(ltr_dict[ltr_name]) + '\n')

def is_ltr_has_structure(ltr_name, line, ref_seq):
    left_size = 8
    internal_size = 3
    parts = line.split(' ')
    chr_name = parts[11]
    lLTR_start = int(parts[3])
    lLTR_end = int(parts[4])
    rLTR_start = int(parts[6])
    rLTR_end = int(parts[7])

    # 取左/右侧 8bp + 3bp
    # 计算左LTR的切片索引，并确保它们在范围内
    left_start = max(lLTR_start - 1 - left_size, 0)
    left_end = min(lLTR_start + internal_size - 1, len(ref_seq))
    left_seq = ref_seq[left_start: left_end]

    # 计算右LTR的切片索引，并确保它们在范围内
    right_start = max(rLTR_end - internal_size, 0)
    right_end = min(rLTR_end + left_size, len(ref_seq))
    right_seq = ref_seq[right_start: right_end]
    cur_seq_name, cur_is_tp, tsd_seq = search_ltr_structure(ltr_name, left_seq, right_seq)
    # cur_seq_name, cur_is_tp = search_ltr_structure_low_copy(ltr_name, left_seq, right_seq)
    return cur_seq_name, cur_is_tp, tsd_seq

def filter_tir(output_path, tir_output_path, full_length_output_dir, threads, left_LTR_contigs, tmp_output_dir, tool_dir, flanking_len, log, debug):
    true_ltr_names = []
    ltr_dict = {}
    with open(output_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ltr_name = parts[0]
            is_ltr = int(parts[1])
            ltr_dict[ltr_name] = is_ltr
            if is_ltr:
                true_ltr_names.append(ltr_name)

    part_size = len(true_ltr_names) // threads
    divided_job_list = []
    # 划分前 n-1 部分
    for i in range(threads - 1):
        divided_job_list.append(true_ltr_names[i * part_size: (i + 1) * part_size])
    # 最后一部分包含剩余的所有元素
    divided_job_list.append(true_ltr_names[(threads - 1) * part_size:])

    # 1. 先在左右窗口中识别同源边界，并搜索 TSD 结构
    # 2. 如果存在，则生成一致性序列，并调用itrsearch判断是否具有终端反向重复结构
    # 检查窗口，是否具有可靠TSD的TIR
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for job_list in divided_job_list:
        TE_type = 'tir'
        job = ex.submit(judge_boundary_v5_batch, job_list, full_length_output_dir, TE_type, flanking_len)
        jobs.append(job)
    ex.shutdown(wait=True)
    # 收集候选的 TIR 序列
    candidate_tirs = {}
    for job in as_completed(jobs):
        results = job.result()
        for cur_seq_name, is_tir, info, cons_seq in results:
            if len(cons_seq) > 0:
                candidate_tirs[cur_seq_name] = cons_seq

    candidate_tir_path = tmp_output_dir + '/candidate_tir.fa'
    store_fasta(candidate_tirs, candidate_tir_path)

    # 检查 候选TIR是否具有 terminal inverted repeats 结构。至此，我们获得了有TSD+TIR结构支持的TIR转座子
    all_confident_tirs = get_confident_TIR(candidate_tir_path, tool_dir)

    # 剩余的序列都当作候选的LTR转座子
    candidate_ltrs = {}
    confident_tirs = {}
    for ltr_name in true_ltr_names:
        if ltr_name not in all_confident_tirs:
            if ltr_name in left_LTR_contigs:
                candidate_ltrs[ltr_name] = left_LTR_contigs[ltr_name]
        else:
            if ltr_name in candidate_tirs:
                confident_tirs[ltr_name] = candidate_tirs[ltr_name]
    candidate_ltr_path = tmp_output_dir + '/candidate_ltr.fa'
    store_fasta(candidate_ltrs, candidate_ltr_path)
    confident_tir_path = tmp_output_dir + '/confident_tir.fa'
    store_fasta(confident_tirs, confident_tir_path)

    # 我们再基于同源搜索比对的方法，过滤掉所有与TIR转座子高度同源的候选LTR序列
    if file_exist(candidate_ltr_path) and file_exist(confident_tir_path):
        # 如果候选LTR中存在TIR，说明是假阳性，应过滤
        blastnResults_path = tmp_output_dir + '/rm_tir.out'
        temp_blast_dir = tmp_output_dir + '/rm_tir_blast'
        multi_process_align(candidate_ltr_path, confident_tir_path, blastnResults_path, temp_blast_dir, threads)
        remain_candidate_tirs = find_tir_in_ltr(blastnResults_path, candidate_ltr_path, confident_tir_path)
        confident_tirs.update(remain_candidate_tirs)

        if not debug:
            os.system('rm -f ' + blastnResults_path)
            os.system('rm -rf ' + temp_blast_dir)

    filter_ltr_names = []
    for ltr_name in true_ltr_names:
        if ltr_name in confident_tirs:
            cur_is_ltr = 0
            filter_ltr_names.append(ltr_name)
        else:
            cur_is_ltr = 1
        ltr_dict[ltr_name] = cur_is_ltr

    log.logger.debug('LTR num: ' + str(len(true_ltr_names)) + ', TIR filter LTR num: ' + str(len(filter_ltr_names)) + ', remaining LTR num: ' + str(len(true_ltr_names)-len(filter_ltr_names)))
    # log.logger.debug(filter_ltr_names)
    with open(tir_output_path, 'w') as f_save:
        for ltr_name in ltr_dict.keys():
            f_save.write(ltr_name+'\t'+str(ltr_dict[ltr_name])+'\n')

    if not debug:
        os.system('rm -f ' + candidate_tir_path)
        os.system('rm -f ' + candidate_ltr_path)


def filter_helitron(output_path, helitron_output_path, full_length_output_dir, threads, left_LTR_contigs, tmp_output_dir, project_dir, flanking_len, log, debug):
    true_ltr_names = []
    ltr_dict = {}
    with open(output_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ltr_name = parts[0]
            is_ltr = int(parts[1])
            ltr_dict[ltr_name] = is_ltr
            if is_ltr:
                true_ltr_names.append(ltr_name)

    part_size = len(true_ltr_names) // threads
    divided_job_list = []
    # 划分前 n-1 部分
    for i in range(threads - 1):
        divided_job_list.append(true_ltr_names[i * part_size: (i + 1) * part_size])
    # 最后一部分包含剩余的所有元素
    divided_job_list.append(true_ltr_names[(threads - 1) * part_size:])

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for job_list in divided_job_list:
        TE_type = 'helitron'
        job = ex.submit(judge_boundary_v5_batch, job_list, full_length_output_dir, TE_type, flanking_len)
        jobs.append(job)
    ex.shutdown(wait=True)
    # 收集候选的 TIR 序列
    candidate_helitrons = {}
    for job in as_completed(jobs):
        results = job.result()
        for cur_seq_name, is_tir, info, cons_seq in results:
            if len(cons_seq) > 0:
                candidate_helitrons[cur_seq_name] = cons_seq

    candidate_helitron_path = tmp_output_dir + '/candidate_helitron.fa'
    store_fasta(candidate_helitrons, candidate_helitron_path)

    tmp_helitron_path = tmp_output_dir + '/tmp_helitron.fa'
    HS_temp_dir = tmp_output_dir + '/HS_temp'
    EA_temp_dir = tmp_output_dir + '/EA_temp'
    # 检查 候选 Helitron 是否具有 Helitron发卡环等 结构。至此，我们获得了有结构支持的Helitron转座子
    all_confident_helitrons = get_confident_Helitron(candidate_helitron_path, tmp_helitron_path, project_dir, HS_temp_dir, EA_temp_dir, flanking_len, threads)

    if not debug:
        os.system('rm -f ' + tmp_helitron_path)
        os.system('rm -rf ' + HS_temp_dir)
        os.system('rm -rf ' + EA_temp_dir)

    # 剩余的序列都当作候选的LTR转座子
    candidate_ltrs = {}
    confident_helitrons = {}
    for ltr_name in true_ltr_names:
        if ltr_name not in all_confident_helitrons:
            candidate_ltrs[ltr_name] = left_LTR_contigs[ltr_name]
        else:
            confident_helitrons[ltr_name] = candidate_helitrons[ltr_name]
    candidate_ltr_path = tmp_output_dir + '/candidate_ltr.fa'
    store_fasta(candidate_ltrs, candidate_ltr_path)
    confident_helitron_path = tmp_output_dir + '/confident_helitron.fa'
    store_fasta(confident_helitrons, confident_helitron_path)

    # 我们再基于同源搜索比对的方法，过滤掉所有与TIR转座子高度同源的候选LTR序列
    if file_exist(candidate_ltr_path) and file_exist(confident_helitron_path):
        # 如果候选LTR中存在Helitron，说明是假阳性，应过滤
        blastnResults_path = tmp_output_dir + '/rm_helitron.out'
        temp_blast_dir = tmp_output_dir + '/rm_helitron_blast'
        multi_process_align(candidate_ltr_path, confident_helitron_path, blastnResults_path, temp_blast_dir, threads)
        remain_candidate_helitrons = find_tir_in_ltr(blastnResults_path, candidate_ltr_path, confident_helitron_path)
        confident_helitrons.update(remain_candidate_helitrons)

        if not debug:
            os.system('rm -f ' + blastnResults_path)
            os.system('rm -rf ' + temp_blast_dir)

    filter_ltr_names = []
    for ltr_name in true_ltr_names:
        if ltr_name in confident_helitrons:
            cur_is_ltr = 0
            filter_ltr_names.append(ltr_name)
        else:
            cur_is_ltr = 1
        ltr_dict[ltr_name] = cur_is_ltr

    if log is not None:
        log.logger.debug('LTR num: ' + str(len(true_ltr_names)) + ', Helitron filter LTR num: ' + str(len(filter_ltr_names)) + ', remaining LTR num: ' + str(len(true_ltr_names)-len(filter_ltr_names)))
        # log.logger.debug(filter_ltr_names)

    with open(helitron_output_path, 'w') as f_save:
        for ltr_name in ltr_dict.keys():
            f_save.write(ltr_name+'\t'+str(ltr_dict[ltr_name])+'\n')

    if not debug:
        os.system('rm -f ' + candidate_helitron_path)
        os.system('rm -f ' + candidate_ltr_path)

def filter_sine(output_path, sine_output_path, full_length_output_dir, threads, left_LTR_contigs, tmp_output_dir, flanking_len, log, debug):
    true_ltr_names = []
    ltr_dict = {}
    with open(output_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ltr_name = parts[0]
            is_ltr = int(parts[1])
            ltr_dict[ltr_name] = is_ltr
            if is_ltr:
                true_ltr_names.append(ltr_name)

    part_size = len(true_ltr_names) // threads
    divided_job_list = []
    # 划分前 n-1 部分
    for i in range(threads - 1):
        divided_job_list.append(true_ltr_names[i * part_size: (i + 1) * part_size])
    # 最后一部分包含剩余的所有元素
    divided_job_list.append(true_ltr_names[(threads - 1) * part_size:])

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for job_list in divided_job_list:
        job = ex.submit(judge_boundary_v6_batch, job_list, full_length_output_dir, flanking_len)
        jobs.append(job)
    ex.shutdown(wait=True)
    # 收集候选的 TIR 序列
    confident_sines = {}
    for job in as_completed(jobs):
        results = job.result()
        for cur_seq_name, is_sine, info, cons_seq in results:
            if len(cons_seq) > 0:
                confident_sines[cur_seq_name] = cons_seq
    confident_sine_path = tmp_output_dir + '/confident_sine.fa'
    store_fasta(confident_sines, confident_sine_path)

    # 剩余的序列都当作候选的LTR转座子
    candidate_ltrs = {}
    for ltr_name in true_ltr_names:
        if ltr_name not in confident_sines:
            candidate_ltrs[ltr_name] = left_LTR_contigs[ltr_name]
    candidate_ltr_path = tmp_output_dir + '/candidate_ltr.fa'
    store_fasta(candidate_ltrs, candidate_ltr_path)

    # 我们再基于同源搜索比对的方法，过滤掉所有与TIR转座子高度同源的候选LTR序列
    if file_exist(candidate_ltr_path) and file_exist(confident_sine_path):
        # 如果候选LTR中存在TIR，说明是假阳性，应过滤
        blastnResults_path = tmp_output_dir + '/rm_sine.out'
        temp_blast_dir = tmp_output_dir + '/rm_sine_blast'
        multi_process_align(candidate_ltr_path, confident_sine_path, blastnResults_path, temp_blast_dir, threads)
        remain_candidate_tirs = find_tir_in_ltr(blastnResults_path, candidate_ltr_path, confident_sine_path)
        confident_sines.update(remain_candidate_tirs)

        if not debug:
            os.system('rm -f ' + blastnResults_path)
            os.system('rm -rf ' + temp_blast_dir)

    filter_ltr_names = []
    for ltr_name in true_ltr_names:
        if ltr_name in confident_sines:
            cur_is_ltr = 0
            filter_ltr_names.append(ltr_name)
        else:
            cur_is_ltr = 1
        ltr_dict[ltr_name] = cur_is_ltr

    log.logger.debug('LTR num: ' + str(len(true_ltr_names)) + ', SINE filter LTR num: ' + str(len(filter_ltr_names)) + ', remaining LTR num: ' + str(len(true_ltr_names)-len(filter_ltr_names)))
    # log.logger.debug(filter_ltr_names)
    with open(sine_output_path, 'w') as f_save:
        for ltr_name in ltr_dict.keys():
            f_save.write(ltr_name+'\t'+str(ltr_dict[ltr_name])+'\n')

    if not debug:
        os.system('rm -f ' + candidate_ltr_path)

def filter_tir_by_tsd(dl_output_path, tsd_output_path, matrix_dir, threads, left_LTR_contigs, tmp_output_dir, tool_dir,
                      log):
    true_ltr_names = []
    ltr_dict = {}
    with open(dl_output_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ltr_name = parts[0]
            is_ltr = int(parts[1])
            ltr_dict[ltr_name] = is_ltr
            if is_ltr:
                true_ltr_names.append(ltr_name)

    # 检查窗口，是否具有可靠TSD的TIR
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for ltr_name in true_ltr_names:
        cur_matrix_file = matrix_dir + '/' + ltr_name + '.matrix'
        job = ex.submit(is_TIR_frame, cur_matrix_file, ltr_name)
        jobs.append(job)
    ex.shutdown(wait=True)
    # 收集候选的 TIR 序列
    candidate_tirs = {}
    candidate_ltrs = {}
    for job in as_completed(jobs):
        cur_seq_name, cur_is_tir = job.result()
        if cur_is_tir:
            candidate_tirs[cur_seq_name] = left_LTR_contigs[cur_seq_name]
        else:
            candidate_ltrs[cur_seq_name] = left_LTR_contigs[cur_seq_name]

    candidate_ltr_path = tmp_output_dir + '/candidate_ltr.fa'
    store_fasta(candidate_ltrs, candidate_ltr_path)

    # 判断序列中是否有 terminal inverted repeats 来判断是否有 TIR
    candidate_tir_path = tmp_output_dir + '/candidate_tir.fa'
    store_fasta(candidate_tirs, candidate_tir_path)
    confident_tir_path = tmp_output_dir + '/confident_tir.fa'
    confident_tirs = get_confident_TIR(candidate_tir_path, tool_dir)
    confident_tir_contigs = {}
    for name in candidate_tirs.keys():
        if name in confident_tirs:
            confident_tir_contigs[name] = candidate_tirs[name]
    store_fasta(confident_tir_contigs, confident_tir_path)

    # 如果候选LTR中存在TIR，说明是假阳性，应过滤
    blastnResults_path = tmp_output_dir + '/rm_tir.out'
    temp_blast_dir = tmp_output_dir + '/rm_tir_blast'
    multi_process_align(candidate_ltr_path, confident_tir_path, blastnResults_path, temp_blast_dir, threads)
    candidate_tirs = find_tir_in_ltr(blastnResults_path, candidate_ltr_path, confident_tir_path)
    confident_tirs.update(candidate_tirs)

    filter_ltr_names = []
    for ltr_name in true_ltr_names:
        if ltr_name in confident_tirs:
            cur_is_ltr = 0
            filter_ltr_names.append(ltr_name)
        else:
            cur_is_ltr = 1
        ltr_dict[ltr_name] = cur_is_ltr

    log.logger.debug('LTR num: ' + str(len(true_ltr_names)) + ', TIR TSD filter LTR num: ' + str(
        len(filter_ltr_names)) + ', remaining LTR num: ' + str(len(true_ltr_names) - len(filter_ltr_names)))
    log.logger.debug(filter_ltr_names)
    with open(tsd_output_path, 'w') as f_save:
        for ltr_name in ltr_dict.keys():
            f_save.write(ltr_name + '\t' + str(ltr_dict[ltr_name]) + '\n')

def find_tir_in_ltr(blastnResults_path, query_path, subject_path, coverage = 0.95):
    full_length_threshold = coverage
    longest_repeats = FMEA_new(query_path, blastnResults_path, full_length_threshold)

    query_names, query_contigs = read_fasta(query_path)
    subject_names, subject_contigs = read_fasta(subject_path)

    candidate_tirs = {}
    for query_name in longest_repeats.keys():
        for item in longest_repeats[query_name]:
            query_name = item[0]
            subject_name = item[3]

            query_len = len(query_contigs[query_name])
            subject_len = len(subject_contigs[subject_name])
            q_start = int(item[1])
            q_end = int(item[2])
            q_len = abs(q_end - q_start)
            s_start = int(item[4])
            s_end = int(item[5])
            s_len = abs(s_end - s_start)
            # 我们认为 比对上的 query 和 subject 部分应该占它自身的95%以上
            if float(q_len) / query_len >= coverage and float(s_len) / subject_len >= coverage:
                candidate_tirs[query_name] = 1
                break
    return candidate_tirs

def get_low_copy_LTR_sub(job_list, low_copy_output_dir, copy_num_threshold):
    finished_list = []
    for matrix_file in job_list:
        cur_copy_num = 0
        with open(matrix_file, 'r') as f_r:
            for line in f_r:
                cur_copy_num += 1
        # 我们增加一个限制，LTR终端的拷贝应该要超过3个
        if cur_copy_num <= copy_num_threshold:
            os.system('cp ' + matrix_file + ' ' + low_copy_output_dir)
        finished_list.append(matrix_file)


def create_or_clear_directory(directory_path):
    # 如果目录已存在，则删除该目录及其中的所有内容
    if os.path.exists(directory_path):
        shutil.rmtree(directory_path)
    # 创建新的目录
    os.makedirs(directory_path)

def get_low_copy_LTR(output_dir, low_copy_output_dir, threads, copy_num_threshold=3):
    create_or_clear_directory(low_copy_output_dir)

    all_matrix_files = []
    for name in os.listdir(output_dir):
        matrix_file = output_dir + '/' + name
        all_matrix_files.append(matrix_file)

    part_size = len(all_matrix_files) // threads
    divided_job_list = []
    # 划分前 n-1 部分
    for i in range(threads - 1):
        divided_job_list.append(all_matrix_files[i * part_size: (i + 1) * part_size])
    # 最后一部分包含剩余的所有元素
    divided_job_list.append(all_matrix_files[(threads - 1) * part_size:])

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for job_list in divided_job_list:
        job = ex.submit(get_low_copy_LTR_sub, job_list, low_copy_output_dir, copy_num_threshold)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        finished_list = job.result()

def get_high_copy_LTR_sub(job_list, high_copy_output_dir, copy_num_threshold):
    finished_list = []
    for matrix_file in job_list:
        cur_copy_num = 0
        with open(matrix_file, 'r') as f_r:
            for line in f_r:
                cur_copy_num += 1
        # 我们增加一个限制，LTR终端的拷贝应该要超过3个
        if cur_copy_num > copy_num_threshold:
            os.system('cp ' + matrix_file + ' ' + high_copy_output_dir)
        finished_list.append(matrix_file)

def get_high_copy_LTR(output_dir, high_copy_output_dir, threads, copy_num_threshold=3):
    create_or_clear_directory(high_copy_output_dir)

    all_matrix_files = []
    for name in os.listdir(output_dir):
        matrix_file = output_dir + '/' + name
        all_matrix_files.append(matrix_file)

    part_size = len(all_matrix_files) // threads
    divided_job_list = []
    # 划分前 n-1 部分
    for i in range(threads - 1):
        divided_job_list.append(all_matrix_files[i * part_size: (i + 1) * part_size])
    # 最后一部分包含剩余的所有元素
    divided_job_list.append(all_matrix_files[(threads - 1) * part_size:])

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for job_list in divided_job_list:
        job = ex.submit(get_high_copy_LTR_sub, job_list, high_copy_output_dir, copy_num_threshold)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        finished_list = job.result()


def random_downsample(data_list, num_samples):
    """
    随机下采样列表，以减少元素数量到指定的num_samples。

    :param data_list: 原始列表。
    :param num_samples: 需要下采样到的元素数量。
    :return: 下采样后的列表。
    """
    # 确保请求的样本数不大于原始列表的元素数量
    if num_samples > len(data_list):
        raise ValueError("样本数不能大于原始列表的元素数量")

    # 使用random.sample进行下采样
    sampled_list = random.sample(data_list, num_samples)

    return sampled_list

def find_files_recursively(root_dir, file_extension=''):
    """
    递归搜索指定目录及其子目录中的文件，并返回文件路径列表。

    参数:
    root_dir (str): 根目录路径。
    file_extension (str): 可选的文件扩展名，例如 '.txt'。如果不指定，则搜索所有文件。

    返回:
    files (list): 匹配的文件路径列表。
    """
    files = []

    # 遍历目录中的每一个条目
    for root, dirs, file_names in os.walk(root_dir):
        # 过滤出具有指定扩展名的文件
        if file_extension:
            filtered_file_names = [f for f in file_names if f.endswith(file_extension)]
        else:
            filtered_file_names = file_names

            # 为每个匹配的文件构建完整的文件路径，并添加到列表中
        for file_name in filtered_file_names:
            files.append(os.path.join(root, file_name))

    return files

def judge_ltr_from_both_ends_frame(output_dir, output_path, threads, type, flanking_len, log, sliding_window_size=20):
    file_extension = '.matrix'
    all_matrix_files = find_files_recursively(output_dir, file_extension)

    part_size = len(all_matrix_files) // threads
    divided_job_list = []
    # 划分前 n-1 部分
    for i in range(threads - 1):
        divided_job_list.append(all_matrix_files[i * part_size: (i + 1) * part_size])
    # 最后一部分包含剩余的所有元素
    divided_job_list.append(all_matrix_files[(threads - 1) * part_size:])

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for job_list in divided_job_list:
        job = ex.submit(judge_both_ends_frame_batch, job_list, sliding_window_size, flanking_len, debug=1)
        jobs.append(job)
    ex.shutdown(wait=True)

    FP_num = 0
    true_ltrs = {}
    for job in as_completed(jobs):
        results = job.result()
        for cur_name, cur_is_ltr in results:
            if cur_is_ltr:
                cur_is_ltr = 1
            else:
                cur_is_ltr = 0
                FP_num += 1
            true_ltrs[cur_name] = cur_is_ltr

    if log is not None:
        log.logger.debug(type + ' LTR num: ' + str(len(all_matrix_files)) + ', LTR Homo filter LTR num: ' + str(FP_num) + ', remain LTR num: ' + str(len(all_matrix_files) - FP_num))

    with open(output_path, 'w') as f_save:
        for cur_name in true_ltrs.keys():
            f_save.write(cur_name + '\t' + str(true_ltrs[cur_name]) + '\n')

def judge_ltr_from_both_ends_frame_for_intactLTR(output_dir, output_path, threads, type, log, sliding_window_size=10):
    file_extension = '.matrix'
    all_matrix_files = find_files_recursively(output_dir, file_extension)

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for matrix_file in all_matrix_files:
        job = ex.submit(judge_both_ends_frame, matrix_file, sliding_window_size, debug=1)
        jobs.append(job)
    ex.shutdown(wait=True)

    FP_num = 0
    true_ltrs = {}
    for job in as_completed(jobs):
        cur_name, cur_is_ltr = job.result()
        # print(cur_name, cur_is_ltr)
        cur_is_ltr = True
        if cur_is_ltr:
            cur_is_ltr = 1
        else:
            cur_is_ltr = 0
            FP_num += 1
        true_ltrs[cur_name] = cur_is_ltr

    if log is not None:
        log.logger.debug(type + ' LTR num: ' + str(len(all_matrix_files)) + ', LTR Homo filter LTR num: ' + str(FP_num) + ', remain LTR num: ' + str(len(all_matrix_files) - FP_num))

    with open(output_path, 'w') as f_save:
        for cur_name in true_ltrs.keys():
            f_save.write(cur_name + '\t' + str(true_ltrs[cur_name]) + '\n')

def search_ltr_structure_low_copy(ltr_name, left_seq, right_seq):
    has_structure = False
    tsd_lens = [6, 5, 4]
    tsd_seq = ''
    exist_tsd = set()
    tsd_motif_distance = 3
    for k_num in tsd_lens:
        for i in range(len(left_seq) - k_num + 1):
            left_kmer = left_seq[i: i + k_num]
            tsd_right_seq = left_seq[i + k_num: i + k_num + tsd_motif_distance]
            if 'N' not in left_kmer and len(set(left_kmer)) != 1 and 'TG' in tsd_right_seq:
                exist_tsd.add(left_kmer)
    for k_num in tsd_lens:
        if has_structure:
            break
        for i in range(len(right_seq) - k_num + 1):
            right_kmer = right_seq[i: i + k_num]
            start_pos = max(i - tsd_motif_distance, 0)
            tsd_left_seq = right_seq[start_pos: i]
            if 'N' not in right_kmer and right_kmer in exist_tsd and 'CA' in tsd_left_seq:
                has_structure = True
                tsd_seq = right_kmer
                break
    # print(ltr_name, left_seq, right_seq, has_tsd, tsd_seq)
    return ltr_name, has_structure

def filter_ltr_by_structure_low_copy(output_path, structure_output_path, leftLtr2Candidates, ltr_lines, reference, threads, log):
    ref_names, ref_contigs = read_fasta(reference)

    true_ltrs = {}
    ltr_dict = {}
    with open(output_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ltr_name = parts[0]
            is_ltr = int(parts[1])
            if is_ltr:
                true_ltrs[ltr_name] = is_ltr
            ltr_dict[ltr_name] = is_ltr

    confident_lines = []
    for name in leftLtr2Candidates.keys():
        # if name not in FP_ltrs:
        if name in true_ltrs:
            candidate_index = leftLtr2Candidates[name]
            line = ltr_lines[candidate_index]
            confident_lines.append((name, line))

    left_size = 6
    internal_size = 0
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for ltr_name, line in confident_lines:
        parts = line.split(' ')
        chr_name = parts[11]
        ref_seq = ref_contigs[chr_name]

        lLTR_start = int(parts[3])
        lLTR_end = int(parts[4])
        rLTR_start = int(parts[6])
        rLTR_end = int(parts[7])

        # 取左/右侧 8bp + 3bp
        # 计算左LTR的切片索引，并确保它们在范围内
        left_start = max(lLTR_start - 1 - left_size, 0)
        left_end = min(lLTR_start + internal_size - 1, len(ref_seq))
        left_seq = ref_seq[left_start: left_end]

        # 计算右LTR的切片索引，并确保它们在范围内
        right_start = max(rLTR_end - internal_size, 0)
        right_end = min(rLTR_end + left_size, len(ref_seq))
        right_seq = ref_seq[right_start: right_end]

        job = ex.submit(search_ltr_structure_low_copy, ltr_name, left_seq, right_seq)
        jobs.append(job)
    ex.shutdown(wait=True)

    FP_ltrs = {}
    for job in as_completed(jobs):
        cur_seq_name, cur_is_tp = job.result()
        if not cur_is_tp:
            FP_ltrs[cur_seq_name] = 1

    log.logger.debug('LTR num: ' + str(len(true_ltrs)) + ', LTR low copy structure filter num: ' + str(len(FP_ltrs)) + ', remaining LTR num: ' + str(len(true_ltrs) - len(FP_ltrs)))
    with open(structure_output_path, 'w') as f_save:
        for ltr_name in ltr_dict.keys():
            if ltr_name in FP_ltrs:
                ltr_dict[ltr_name] = 0
            f_save.write(ltr_name + '\t' + str(ltr_dict[ltr_name]) + '\n')

def judge_ltr_has_structure(lc_output_path, structure_output_path, leftLtr2Candidates, ltr_lines, reference, log):
    ref_names, ref_contigs = read_fasta(reference)

    true_ltr_names = []
    ltr_dict = {}
    with open(lc_output_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            ltr_name = parts[0]
            is_ltr = int(parts[1])
            ltr_dict[ltr_name] = is_ltr
            if is_ltr:
                true_ltr_names.append(ltr_name)

    filter_ltr_names = []
    for ltr_name in true_ltr_names:
        candidate_index = leftLtr2Candidates[ltr_name]
        ltr_line = ltr_lines[candidate_index]
        parts = ltr_line.split(' ')
        chr_name = parts[11]
        ref_seq = ref_contigs[chr_name]
        left_ltr_start = int(parts[3])
        left_ltr_end = int(parts[4])
        right_ltr_start = int(parts[6])
        right_ltr_end = int(parts[7])

        left_ltr_seq = ref_seq[left_ltr_start - 1: left_ltr_end]
        right_ltr_seq = ref_seq[right_ltr_start - 1: right_ltr_end]

        tsd_lens = [6, 5, 4]
        allow_mismatch_num = 0
        has_structure = False
        for tsd_len in tsd_lens:
            left_tsd = ref_seq[left_ltr_start - 1 - tsd_len: left_ltr_start - 1]
            right_tsd = ref_seq[right_ltr_end: right_ltr_end + tsd_len]
            if allow_mismatch(left_tsd, right_tsd, allow_mismatch_num):
                has_structure = True
                break
        if has_structure:
            ltr_dict[ltr_name] = 1
        else:
            ltr_dict[ltr_name] = 0
            filter_ltr_names.append(ltr_name)

    log.logger.debug('Low copy LTR after HomoFilter num: ' + str(len(true_ltr_names)) + ', LTR structure filter LTR num: ' + str(len(filter_ltr_names)) + ', remain LTR num: ' + str(len(true_ltr_names)-len(filter_ltr_names)))

    with open(structure_output_path, 'w') as f_save:
        for ltr_name in ltr_dict.keys():
            f_save.write(ltr_name + '\t' + str(ltr_dict[ltr_name]) + '\n')

def judge_both_ends_frame_batch(job_list, sliding_window_size, flanking_len, debug=1):
    results = []
    for maxtrix_file in job_list:
        seq_name, is_ltr = judge_both_ends_frame(maxtrix_file, sliding_window_size, flanking_len, debug=1)
        results.append((seq_name, is_ltr))
    return results


def judge_both_ends_frame(maxtrix_file, sliding_window_size, flanking_len, debug=1):
    # 我现在想的假阳性过滤方法：
    # 1. 对matrix file 搜索同源边界，如果不存在，则说明是真实LTR，否则为假阳性
    is_ltr = True
    seq_name = os.path.basename(maxtrix_file).split('.')[0]
    # Step3. 对候选随机序列 搜索同源边界。我们将窗口设置为20.
    is_left_ltr, new_boundary_start = judge_left_frame_LTR(maxtrix_file, flanking_len, sliding_window_size=sliding_window_size)
    # if debug:
    #     print('left', maxtrix_file, is_left_ltr, new_boundary_start)
    is_right_ltr, new_boundary_end = judge_right_frame_LTR(maxtrix_file, flanking_len, sliding_window_size=sliding_window_size)
    # if debug:
    #     print('right', maxtrix_file, is_right_ltr, new_boundary_end)
    is_ltr &= is_left_ltr and is_right_ltr
    return seq_name, is_ltr

def alter_deep_learning_results(dl_output_path, hc_output_path, alter_dl_output_path, high_copy_output_dir, log):
    ltr_dict1 = {}
    if file_exist(dl_output_path):
        with open(dl_output_path, 'r') as f_r:
            for line in f_r:
                line = line.replace('\n', '')
                parts = line.split('\t')
                ltr_name = parts[0]
                is_ltr = int(parts[1])
                ltr_dict1[ltr_name] = is_ltr

    ltr_dict2 = {}
    if file_exist(hc_output_path):
        with open(hc_output_path, 'r') as f_r:
            for line in f_r:
                line = line.replace('\n', '')
                parts = line.split('\t')
                ltr_name = parts[0]
                is_ltr = int(parts[1])
                ltr_dict2[ltr_name] = is_ltr

    if len(ltr_dict1) == 0 and len(ltr_dict2) == 0:
        # 如果同源搜索模块和深度学习模块同时不执行
        file_extension = '.matrix'
        all_matrix_files = find_files_recursively(high_copy_output_dir, file_extension)
        ltr_dict = {}
        for matrix_file in all_matrix_files:
            seq_name = os.path.basename(matrix_file).split('.')[0]
            ltr_dict[seq_name] = 1
    elif len(ltr_dict1) == 0:
        ltr_dict = ltr_dict2
        log.logger.debug('No deep learning prediction is found, use the homology rule prediction result.')
    elif len(ltr_dict2) == 0:
        ltr_dict = ltr_dict1
        log.logger.debug('No homology rule prediction is found, use the deep learning prediction result.')
    else:
        for ltr_name in ltr_dict2.keys():
            is_ltr = ltr_dict2[ltr_name]
            # 当同源方法预测为 0，即非LTR时，以其预测结果为准
            if not is_ltr:
                ltr_dict1[ltr_name] = is_ltr
        ltr_dict = ltr_dict1
        log.logger.debug('Adjust the deep learning prediction result using the homology rule prediction result.')

    with open(alter_dl_output_path, 'w') as f_save:
        for ltr_name in ltr_dict.keys():
            f_save.write(ltr_name + '\t' + str(ltr_dict[ltr_name]) + '\n')

def judge_scn_line_by_flank_seq_v2(job_list):
    results = {}
    for cur_job in job_list:
        candidate_index, left_ltr, right_ltr, lLTR_len, rLTR_len, lLTR_start, lLTR_end, rLTR_start, rLTR_end = cur_job
        extend_len = 50
        error_region_len = 10

        # 创建临时文件
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file, \
             tempfile.NamedTemporaryFile(mode='w', delete=False) as subject_file, \
             tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:

            # 将序列写入临时文件
            query_file.write(left_ltr)
            subject_file.write(right_ltr)

            # 运行 BLAST，将输出写入文件
            blastn_command = [
                "blastn",
                "-subject", subject_file.name,
                "-query", query_file.name,
                "-outfmt", "6",
                "-num_threads", "1",
                "-out", output_file.name
            ]
            result = subprocess.run(blastn_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            # 读取 BLAST 结果
            adjust_boundary = None
            q_start_offset = 0
            q_end_offset = 0
            s_start_offset = 0
            s_end_offset = 0
            precise_boundary_alignments = []
            if result.returncode == 0:
                with open(output_file.name, 'r') as f:
                    blastn_lines = f.readlines()
                for blastn_line in blastn_lines:
                    if blastn_line.strip():
                        parts = blastn_line.strip().split('\t')
                        identity = float(parts[2])
                        q_start = int(parts[6])
                        q_end = int(parts[7])
                        s_start = int(parts[8])
                        s_end = int(parts[9])
                        # 比对落在了左侧边界上
                        if (((extend_len - error_region_len) <= q_start <= (extend_len + error_region_len)
                                and q_end <= (extend_len + lLTR_len + error_region_len)
                                and (extend_len - error_region_len) <= s_start <= (extend_len + error_region_len)
                                and s_end <= (extend_len + rLTR_len + error_region_len))):
                            q_start_offset = q_start - extend_len
                            s_start_offset = s_start - extend_len
                            precise_boundary_alignments.append((q_start, q_end))
                        # 比对落在了右侧边界上
                        if (((extend_len - error_region_len) <= q_start
                                and (extend_len + lLTR_len - error_region_len) <= q_end <= (extend_len + lLTR_len + error_region_len)
                                and (extend_len - error_region_len) <= s_start
                                and (extend_len + rLTR_len - error_region_len) <= s_end <= (extend_len + rLTR_len + error_region_len))):
                            q_end_offset = q_end + 1 - (lLTR_len + extend_len)
                            s_end_offset = s_end + 1 - (rLTR_len + extend_len)
                            precise_boundary_alignments.append((q_start, q_end))

            # 清理临时文件
            os.remove(query_file.name)
            os.remove(subject_file.name)
            os.remove(output_file.name)

            # 合并 precise_boundary_alignments 坐标后，判断是否覆盖 lLTR_len 的 95%
            merged_alignments = merge_intervals(precise_boundary_alignments)
            cover_len = sum(end - start for start, end in merged_alignments)
            if float(cover_len) / lLTR_len >= 0.9:
                is_FP = False
                lLTR_start += q_start_offset
                lLTR_end += q_end_offset
                rLTR_start += s_start_offset
                rLTR_end += s_end_offset
                if q_start_offset != 0 or q_end_offset != 0 or s_start_offset != 0 or s_end_offset != 0:
                    adjust_boundary = (lLTR_start, lLTR_end, rLTR_start, rLTR_end)
            else:
                is_FP = True

            if is_FP:
                is_ltr = 0
            else:
                is_ltr = 1

            results[candidate_index] = (is_ltr, adjust_boundary)
    return results

def filter_ltr_by_flank_seq_v2(scn_file, filter_scn, reference, threads, log):
    ref_names, ref_contigs = read_fasta(reference)
    ltr_candidates, ltr_lines = read_scn(scn_file, log)

    job_list = []
    for candidate_index in ltr_lines.keys():
        line = ltr_lines[candidate_index]
        parts = line.split(' ')
        chr_name = parts[11]
        ref_seq = ref_contigs[chr_name]
        extend_len = 50
        lLTR_start = int(parts[3])
        lLTR_end = int(parts[4])
        lLTR_len = int(parts[5])
        rLTR_start = int(parts[6])
        rLTR_end = int(parts[7])
        rLTR_len = int(parts[8])
        left_ltr = ref_seq[lLTR_start - extend_len: lLTR_end + extend_len]
        right_ltr = ref_seq[rLTR_start - extend_len: rLTR_end + extend_len]
        job_list.append((candidate_index, left_ltr, right_ltr, lLTR_len, rLTR_len, lLTR_start, lLTR_end, rLTR_start, rLTR_end))

    part_size = len(job_list) // threads
    divided_job_list = []
    # 划分前 n-1 部分
    for i in range(threads - 1):
        divided_job_list.append(job_list[i * part_size: (i + 1) * part_size])
    # 最后一部分包含剩余的所有元素
    divided_job_list.append(job_list[(threads - 1) * part_size:])

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cur_job_list in divided_job_list:
        job = ex.submit(judge_scn_line_by_flank_seq_v2, cur_job_list)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_results = {}
    for job in as_completed(jobs):
        cur_results = job.result()
        all_results.update(cur_results)

    adjust_boundaries = {}
    is_ltrs = {}
    fp_count = 0
    for candidate_index in all_results.keys():
        is_ltr, adjust_boundary = all_results[candidate_index]
        if is_ltr == 0:
            fp_count += 1
        is_ltrs[candidate_index] = is_ltr
        adjust_boundaries[candidate_index] = adjust_boundary

    confident_lines = []
    for candidate_index in ltr_lines.keys():
        is_ltr = is_ltrs[candidate_index]
        if is_ltr == 1:
            line = ltr_lines[candidate_index]
            confident_lines.append(line)
            # 生成新的调整边界后的记录
            adjust_boundary = adjust_boundaries[candidate_index]
            if adjust_boundary is not None:
                parts = line.split(' ')
                parts[3] = str(adjust_boundary[0])
                parts[4] = str(adjust_boundary[1])
                parts[5] = str(adjust_boundary[1] - adjust_boundary[0] + 1)
                parts[6] = str(adjust_boundary[2])
                parts[7] = str(adjust_boundary[3])
                parts[8] = str(adjust_boundary[3] - adjust_boundary[2] + 1)
                new_line = ' '.join(parts)
                confident_lines.append(new_line)
    store_scn(confident_lines, filter_scn)
    if log is not None:
        log.logger.debug('Remove False Positive LTR terminal: ' + str(fp_count) + ', remaining LTR num: ' + str(len(confident_lines)))


def filter_ltr_by_flanking_cluster_sub_batch(job_list, target_dir, temp_dir):
    results = []
    for matrix_file in job_list:
        matrix_file_name = os.path.basename(matrix_file)
        cur_seq_file = temp_dir + '/' + matrix_file_name + '.fa'
        cur_seq_cons = cur_seq_file + '.cons'
        cur_seq_name, cur_is_ltr = filter_ltr_by_flanking_cluster_sub(matrix_file, target_dir, cur_seq_file, cur_seq_cons)
        results.append((cur_seq_name, cur_is_ltr))
    return results

def filter_ltr_by_flanking_cluster_sub(matrix_file, target_dir, cur_seq_file, cur_seq_cons):
    cur_seq_name = os.path.basename(matrix_file).split('.')[0]
    contigs = {}
    with open(matrix_file, 'r') as f_r:
        for i, line in enumerate(f_r):
            seq = line.replace('-', '').replace('\t', '').replace('\n', '')
            seq_name = 'seq_' + str(i)
            contigs[seq_name] = seq
    store_fasta(contigs, cur_seq_file)
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.95) \
                     + ' -G 0 -g 1 -A 50 -i ' + cur_seq_file + ' -o ' + cur_seq_cons + ' -T 1 -M 0'
    os.system(cd_hit_command + ' > /dev/null 2>&1')
    # 解析聚类文件中的单拷贝
    cluster_file = cur_seq_cons + '.clstr'
    cluster_idx = -1
    clusters = {}
    with open(cluster_file, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            if line.startswith('>'):
                cluster_idx = line.split(' ')[1]
            else:
                if not clusters.__contains__(cluster_idx):
                    clusters[cluster_idx] = []
                cur_cluster = clusters[cluster_idx]
                name = line.split(',')[1].split(' ')[1].strip()[1:]
                name = name[0: len(name) - 3]
                cur_cluster.append(name)
    # 计算随机簇序列数量
    homo_seq_num = 0
    random_seq_num = 0
    for cluster_idx in clusters.keys():
        cur_cluster = clusters[cluster_idx]
        if len(cur_cluster) > 1:
            homo_seq_num += len(cur_cluster)
        else:
            random_seq_num += len(cur_cluster)

    cur_is_ltr = True
    if homo_seq_num >= random_seq_num:
        cur_is_ltr = False
    if cur_is_ltr:
        cp_command = 'cp -f ' + matrix_file + ' ' + target_dir
        # print(cp_command)
        os.system(cp_command)
    return cur_seq_name, cur_is_ltr

def filter_ltr_by_flanking_cluster(output_dir, target_dir, temp_dir, threads, log):
    # 根据窗口两侧 50 bp序列连在一起，使用cd-hit-est聚类，真正 LTR 的绝大多数拷贝应该不能够聚类在一起（考虑到LTR插入到其他TE中，而整个 TE 转座子几次，导致会出现少数一致的侧翼区域）
    file_extension = '.matrix'
    all_matrix_files = find_files_recursively(output_dir, file_extension)
    os.system('rm -rf ' + target_dir)
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    part_size = len(all_matrix_files) // threads
    divided_job_list = []
    # 划分前 n-1 部分
    for i in range(threads - 1):
        divided_job_list.append(all_matrix_files[i * part_size: (i + 1) * part_size])
    # 最后一部分包含剩余的所有元素
    divided_job_list.append(all_matrix_files[(threads - 1) * part_size:])

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for job_list in divided_job_list:
        job = ex.submit(filter_ltr_by_flanking_cluster_sub_batch, job_list, target_dir, temp_dir)
        jobs.append(job)
    ex.shutdown(wait=True)

    input_num = 0
    keep_num = 0
    for job in as_completed(jobs):
        cur_results = job.result()
        for cur_name, cur_is_ltr in cur_results:
            if cur_is_ltr:
                keep_num += 1
            input_num += 1
    if log is not None:
        log.logger.debug('Input LTR num: ' + str(input_num) + ', keep LTR num: ' + str(keep_num))

def get_intact_ltr_copies(ltr_copies, ltr_lines, full_length_threshold, reference):
    std_fragments = {}
    for candidate_index in ltr_lines.keys():
        line = ltr_lines[candidate_index]
        parts = line.split(' ')
        chr_name = parts[11]
        chr_start = int(parts[3]) - 1
        chr_end = int(parts[7])
        seq_name = chr_name + '_' + str(chr_start) + '-' + str(chr_end)
        if chr_name not in std_fragments:
            std_fragments[chr_name] = []
        fragments = std_fragments[chr_name]
        fragments.append((chr_start, chr_end, seq_name))

    coverage_threshold = full_length_threshold
    # 由于比对的问题，获取的拷贝并不总是带终端的，因此实际上不算全长拷贝。考虑到intact LTR 总是带终端的，它的拷贝应该也能被原始的程序识别到。
    # 因此我们获取拷贝之后再和原始结果计算overlap，如果overlap超过 95% 才算一个真的全长拷贝，否则不算。
    test_fragments = {}
    for ltr_name in ltr_copies:
        for copy in ltr_copies[ltr_name]:
            chr_name = copy[0]
            chr_start = copy[1]
            chr_end = copy[2]
            seq_name = chr_name + '_' + str(chr_start) + '-' + str(chr_end)
            if chr_name not in test_fragments:
                test_fragments[chr_name] = []
            fragments = test_fragments[chr_name]
            fragments.append((chr_start, chr_end, seq_name, ltr_name))

    names, contigs = read_fasta(reference)
    chrom_length = {}
    for i, name in enumerate(names):
        chr_len = len(contigs[name])
        chrom_length[name] = chr_len

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
    for chr_name in std_fragments.keys():
        fragments = std_fragments[chr_name]
        cur_chr_segments = chr_segments[chr_name]
        for cur_frag in fragments:
            start = cur_frag[0]
            end = cur_frag[1]
            seq_name = cur_frag[2]
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
                segment_frags.append([start, end, seq_name])

    intact_ltr_copies = {}
    # Map the fragments to the corresponding segment,
    # and check if there is an overlap of over 95% with the fragment in the segment.
    for chr_name in test_fragments.keys():
        fragments = test_fragments[chr_name]
        cur_chr_segments = chr_segments[chr_name]
        for cur_frag in fragments:
            start = cur_frag[0]
            end = cur_frag[1]
            seq_name = cur_frag[2]
            ltr_name = cur_frag[3]
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
            if is_found:
                if ltr_name not in intact_ltr_copies:
                    intact_ltr_copies[ltr_name] = []
                copies = intact_ltr_copies[ltr_name]
                copies.append((chr_name, start, end))

    # 去除全长拷贝中的冗余拷贝，即短的那一条包含在长的里面
    for ltr_name in intact_ltr_copies.keys():
        copies = intact_ltr_copies[ltr_name]
        copies.sort(key=lambda x: (x[2]-x[1]))
        filtered_copies = []
        for i in range(len(copies)):
            chr_name_i, start_i, end_i = copies[i]
            length_i = end_i - start_i + 1
            is_redundant = False
            for j in range(i + 1, len(copies)):
                chr_name_j, start_j, end_j = copies[j]
                length_j = end_j - start_j + 1
                if chr_name_i == chr_name_j:
                    overlap_start = max(start_i, start_j)
                    overlap_end = min(end_i, end_j)
                    overlap_length = overlap_end - overlap_start + 1
                    # 判断是否冗余（即 95% 的拷贝 i 被包含在拷贝 j 中）
                    if overlap_length >= 0.95 * length_i:
                        is_redundant = True
                        break
            if not is_redundant:
                filtered_copies.append((chr_name_i, start_i, end_i))
        intact_ltr_copies[ltr_name] = filtered_copies
    return intact_ltr_copies

def get_ltr_from_line(cur_internal_seqs):
    all_lines = {}

    for candidate_index, lLTR_end, int_seq, chr_name, seq_id, line in cur_internal_seqs:
        # 创建临时文件
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file, \
             tempfile.NamedTemporaryFile(mode='w', delete=False) as subject_file, \
             tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:

            # 将序列写入临时文件
            query_file.write(int_seq)
            subject_file.write(int_seq)

            # 运行 BLAST，将输出写入文件
            blastn_command = [
                "blastn",
                "-subject", subject_file.name,
                "-query", query_file.name,
                "-outfmt", "6",
                "-evalue", "1e-20",
                "-num_threads", "1",
                "-out", output_file.name
            ]
            result = subprocess.run(blastn_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            # 读取 BLAST 结果
            new_lines = []
            new_lines.append(line.split(' '))
            if result.returncode == 0:
                with open(output_file.name, 'r') as f:
                    blastn_lines = f.readlines()
                for blastn_line in blastn_lines:
                    if blastn_line.strip():
                        parts = blastn_line.strip().split('\t')
                        q_start = int(parts[6])
                        q_end = int(parts[7])
                        q_len = abs(q_end - q_start)
                        s_start = int(parts[8])
                        s_end = int(parts[9])
                        s_len = abs(s_end - s_start)
                        identity = float(parts[2])
                        if identity > 95 and q_len > 300 and s_len > 300 and q_end > q_start and s_end > s_start and s_start - q_end > 500:
                            new_lLTR_start = lLTR_end + q_start
                            new_lLTR_end = lLTR_end + q_end
                            new_rLTR_start = lLTR_end + s_start
                            new_rLTR_end = lLTR_end + s_end
                            new_line = [
                                new_lLTR_start, new_rLTR_end, new_rLTR_end - new_lLTR_start + 1,
                                new_lLTR_start, new_lLTR_end, new_lLTR_end - new_lLTR_start + 1,
                                new_rLTR_start, new_rLTR_end, new_rLTR_end - new_rLTR_start + 1,
                                identity, seq_id, chr_name
                            ]
                            new_lines.append(new_line)

            # 清理临时文件
            os.remove(query_file.name)
            os.remove(subject_file.name)
            os.remove(output_file.name)

            # 过滤冗余行
            new_lines.sort(key=lambda x: (int(x[0]), -int(x[1])))
            filtered_new_lines = []
            for i in range(len(new_lines) - 1, -1, -1):
                cur_item = new_lines[i]
                cur_lLTR_start = int(cur_item[0])
                cur_rLTR_end = int(cur_item[1])
                cur_LTR_len = int(cur_item[2])
                is_redundant = False
                for j in range(i - 1, -1, -1):
                    next_item = new_lines[j]
                    next_lLTR_start = int(next_item[0])
                    next_rLTR_end = int(next_item[1])
                    next_LTR_len = int(next_item[2])
                    overlap_start = max(cur_lLTR_start, next_lLTR_start)
                    overlap_end = min(cur_rLTR_end, next_rLTR_end)
                    overlap_length = overlap_end - overlap_start + 1
                    if overlap_length >= 0.95 * cur_LTR_len and overlap_length >= 0.95 * next_LTR_len:
                        is_redundant = True
                        break
                if not is_redundant:
                    filtered_new_lines.append(cur_item)

            # 将结果转换为字符串
            new_line_strs = [' '.join(map(str, line)) for line in reversed(filtered_new_lines)]
            all_lines[candidate_index] = new_line_strs

    return all_lines

def get_all_potential_ltr_lines(confident_lines, reference, threads, temp_path, log):
    log.logger.info('Start get all potential ltr lines')
    ref_names, ref_contigs = read_fasta(reference)

    part_size = len(confident_lines) // threads
    result = []
    # 划分前 n-1 部分
    for i in range(threads - 1):
        result.append(confident_lines[i * part_size: (i + 1) * part_size])
    # 最后一部分包含剩余的所有元素
    result.append(confident_lines[(threads - 1) * part_size:])

    ex = ProcessPoolExecutor(threads)
    jobs = []
    candidate_index = 0
    for cur_lines in result:
        cur_internal_seqs = []
        for cur_line in cur_lines:
            parts = cur_line.split(' ')
            chr_name = parts[11]
            ref_seq = ref_contigs[chr_name]
            seq_id = parts[10]
            lLTR_start = int(parts[3])
            lLTR_end = int(parts[4])
            rLTR_start = int(parts[6])
            rLTR_end = int(parts[7])
            int_seq = ref_seq[lLTR_end: rLTR_start]
            cur_internal_seqs.append((candidate_index, lLTR_end, int_seq, chr_name, seq_id, cur_line))
            candidate_index += 1

        job = ex.submit(get_ltr_from_line, cur_internal_seqs)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_lines = {}
    temp_lines = {}
    for job in as_completed(jobs):
        cur_all_lines = job.result()
        all_lines.update(cur_all_lines)

    for candidate_index in all_lines.keys():
        new_lines = all_lines[candidate_index]
        if len(new_lines) > 1:
            temp_lines[candidate_index] = new_lines

    # 把ltr_copies存成文件
    with open(temp_path, 'w', encoding='utf-8') as f:
        json.dump(temp_lines, f, ensure_ascii=False, indent=4)

    new_confident_lines = []
    for i in range(len(confident_lines)):
        new_lines = all_lines[i]
        for cur_line in new_lines:
            new_confident_lines.append(cur_line)

    if log is not None:
        log.logger.info('Find all potential LTR, raw ltr num:' + str(len(confident_lines)) + ', current ltr num:' + str(len(new_confident_lines)))
    return new_confident_lines

