import codecs
import json
import multiprocessing
import os

import logging
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from logging import handlers

import psutil
import pysam

from command import run_bwa, run_minimap2, run_bowtie2

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

    (dir, filename) = os.path.split(reference)
    (name, extension) = os.path.splitext(filename)
    reference_pre = dir + '/' + name + '_preprocess' + extension
    with open(reference_pre, "w") as f_save:
        for contigName in contigNames:
            contigseq = contigs[contigName]
            f_save.write(">" + contigName + '\n' + contigseq + '\n')
    f_save.close()
    return reference_pre

def generate_candidate_repeats_v2(contigs, k_num, unique_kmer_map, partiton_index, fault_tolerant_bases):
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
            elif cur_repeat_str != '':
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
    return repeat_dict

def cut_repeat_v2(sam_path_bwa, repeats_path, cut_repeats_path):
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
        q_start = int(read.query_alignment_start)  # [q_start, q_end)
        q_end = int(read.query_alignment_end)
        if q_start > q_end:
            tmp = q_start
            q_start = q_end
            q_end = tmp
        if not query_records.__contains__(query_name):
            query_records[query_name] = []
        records = query_records[query_name]
        records.append((reference_name, alignment_len, identity, q_start, q_end))
        query_records[query_name] = records

    repeat_contignames, repeat_contigs = read_fasta(repeats_path)
    cut_repeats = {}
    for query_name in query_records.keys():
        query_seq = repeat_contigs[query_name]
        query_len = len(query_seq)
        records = query_records[query_name]
        for i, record in enumerate(records):
            # filter first alignment
            if i == 0:
                continue
            identity = record[2]
            q_start = record[3]
            q_end = record[4]
            if identity < 95:
                continue
            # get repeats boundary by getting all alignment sequences
            new_seq = query_seq[q_start: q_end]
            new_query_name = query_name + '-p_' + str(i) + '-len_' + str(len(new_seq))
            cut_repeats[new_query_name] = new_seq
    store_fasta(cut_repeats, cut_repeats_path)

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
    with open(fasta_path, 'r') as rf:
        contigname = ''
        contigseq = ''
        for line in rf:
            if line.startswith('>'):
                if contigname != '' and contigseq != '':
                    contigs[contigname] = contigseq
                    contignames.append(contigname)
                contigname = line.strip()[1:].split(" ")[0]
                contigseq = ''
            else:
                contigseq += line.strip().upper()
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

    with open(repeats_bwa, 'w') as f_save:
        for r_name in short_seqs.keys():
            f_save.write('>' + r_name + '\n' + short_seqs[r_name] + '\n')

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

def run_alignment(repeat_contig_path, reference_path, use_align_tools, thread_num, tools_dir):
    sam_path = ''
    if use_align_tools == 'bwa':
        sam_path = run_bwa(repeat_contig_path, reference_path, thread_num, tools_dir)
    elif use_align_tools == 'minimap2':
        sam_path = run_minimap2(repeat_contig_path, reference_path, thread_num, tools_dir)
    elif use_align_tools == 'bowtie2':
        sam_path = run_bowtie2(repeat_contig_path, reference_path, thread_num, tools_dir)
    return sam_path

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
            end = len(contig)
            while start < end:
                # add extra kmer length
                seg = contig[start:start+line_len]
                line = '>' + contigName + '\t' + str(start) + '\t' + seg + '\n'
                f_w.write(line)
                start += line_len
    f_w.close()
    return tmp_fasta_path

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
        cur_segment[new_ref_name] = seq
    # only one partition
    if len(cur_segment) > 0:
        segments_cluster[partition_index] = cur_segment
    else:
        if last_index != -1:
            for j in range(last_index, len(segments)):
                parts = segments[i].split('\t')
                ref_name = parts[0].replace('>', '')
                start = parts[1]
                seq = parts[2]
                new_ref_name = ref_name + '$' + start
                cur_segment[new_ref_name] = seq
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
def generate_blastlike_output(sam_paths, blastnResults_path, not_multi_mapping_repeatIds):
    not_multi_mapping_repeatIds_dict = {}
    for repeat_id in not_multi_mapping_repeatIds:
        not_multi_mapping_repeatIds_dict[repeat_id] = 1
    query_names = {}
    with open(blastnResults_path, 'w') as f_save:
        for item in sam_paths:
            sam_path = item[0]
            repeats_path = item[1]
            contignames, contigs = read_fasta(repeats_path)
            samfile = pysam.AlignmentFile(sam_path, "rb")
            for read in samfile.fetch():
                if read.is_unmapped:
                    continue
                query_name = read.query_name
                # filter not multiple mapping repeat
                if not_multi_mapping_repeatIds_dict.__contains__(query_name):
                    continue
                target_name = read.reference_name
                cigar = read.cigarstring
                NM_tag = 0
                try:
                    NM_tag = read.get_tag('NM')
                except KeyError:
                    NM_tag = -1
                identity = compute_identity(cigar, NM_tag, 'BLAST')
                identity = float(identity) * 100
                match_base = int(read.query_alignment_length)
                q_start = int(read.query_alignment_start)
                q_end = int(read.query_alignment_end)
                t_start = int(read.reference_start)
                t_end = int(read.reference_end)
                query_length = int(len(contigs[query_name]))

                if not query_names.__contains__(query_name) and (identity > 95 and float(match_base)/query_length > 0.95):
                    query_names[query_name] = 1
                    continue

                if read.is_reverse:
                    temp = t_start
                    t_start = t_end
                    t_end = temp
                    strand = '+' if t_end >= t_start else '-'

                f_save.write(str(query_name) + '\t' + str(target_name) + '\t' + str(identity) + '\t' + str(match_base)
                             + '\t' + str('X') + '\t' + str('X') + '\t' + str(q_start) + '\t' + str(q_end)
                             + '\t' + str(t_start) + '\t' + str(t_end) + '\n')

def generate_blastlike_output_parallel(sam_paths, blastnResults_path, not_multi_mapping_repeatIds, partitions_num):
    not_multi_mapping_repeatIds_dict = {}
    for repeat_id in not_multi_mapping_repeatIds:
        not_multi_mapping_repeatIds_dict[repeat_id] = 1

    parse_records = []
    for sam_path in sam_paths:
        samfile = pysam.AlignmentFile(sam_path, "rb")
        for read in samfile.fetch():
            if read.is_unmapped:
                continue
            query_name = read.query_name
            target_name = read.reference_name
            cigar = read.cigarstring
            try:
                NM_tag = read.get_tag('NM')
            except KeyError:
                NM_tag = -1
            identity = compute_identity(cigar, NM_tag, 'BLAST')
            identity = float(identity) * 100
            match_base = int(read.query_alignment_length)
            q_start = int(read.query_alignment_start)
            q_end = int(read.query_alignment_end)
            t_start = int(read.reference_start)
            t_end = int(read.reference_end)

            if read.is_reverse:
                temp = t_start
                t_start = t_end
                t_end = temp
                strand = '+' if t_end >= t_start else '-'
            parse_records.append((query_name, target_name, identity, match_base, q_start, q_end, t_start, t_end))

    records_cluster = split2cluster(parse_records, partitions_num)

    blast_records = multiprocessing.Manager().dict()
    pool = multiprocessing.Pool(processes=partitions_num)
    for partiton_index in records_cluster.keys():
        cur_records = records_cluster[partiton_index]
        pool.apply_async(filter_not_multi_mapping, (cur_records, not_multi_mapping_repeatIds_dict, partiton_index, blast_records,))
    pool.close()
    pool.join()

    with open(blastnResults_path, 'w') as f_save:
        for partiton_index in blast_records.keys():
            for record in blast_records[partiton_index]:
                f_save.write(str(record[0]) + '\t' + str(record[1]) + '\t' + str(record[2]) + '\t' + str(record[3])
                             + '\t' + str('X') + '\t' + str('X') + '\t' + str(record[4]) + '\t' + str(record[5])
                             + '\t' + str(record[6]) + '\t' + str(record[7]) + '\n')

def cut_repeat_v1(sam_paths, HS_gap, ID_gap, repeats_file, raw_cut_file):
    repeat_contignames, repeat_contigs = read_fasta(repeats_file)
    all_fragments = {}

    # if a repeat can only align to a position complete,
    # and other records contain large H/S or I/D, it should be spliced
    query_records = {}
    for sam_path in sam_paths:
        samfile = pysam.AlignmentFile(sam_path, "rb")
        for read in samfile.fetch():
            if read.is_unmapped:
                continue
            query_name = read.query_name
            reference_name = read.reference_name
            cigar = read.cigartuples
            cigarstr = read.cigarstring
            is_reverse = read.is_reverse

            if not query_records.__contains__(query_name):
                query_records[query_name] = []
            records = query_records[query_name]
            records.append((query_name, reference_name, cigar, cigarstr, is_reverse))
            query_records[query_name] = records

    # delete Chimerism to avoid fragments
    # for query_name in query_records.keys():
    #     is_chimerism = False
    #     for record in query_records[query_name]:
    #         query_seq = repeat_contigs[query_name]
    #         query_len = len(query_seq)
    #         float(record.query_alignment_length)/query_len > 90
    #     if is_chimerism:
    #         del query_records['query_name']

    pattern = r'[0-9]\d+M'
    repeats_tobe_spliced = {}
    for query_name in query_records.keys():
        complete_alignment_num = 0
        for record in query_records[query_name]:
            cigar = record[2]
            cigarstr = str(record[3])
            # complete Match in cigar
            if re.match(pattern, cigarstr) is not None:
                complete_alignment_num += 1
            else:
                # if cigar contains small I/D or small H/S, it can be seen as a complete record
                query_seq = repeat_contigs[query_name]
                query_len = len(query_seq)
                is_complete = True
                for c in cigar:
                    if (c[0] == 4 and c[1] >= HS_gap * query_len) \
                            or (c[0] == 5 and c[1] >= HS_gap * query_len) \
                            or (c[0] == 1 and c[1] >= ID_gap * query_len) \
                            or (c[0] == 2 and c[1] >= ID_gap * query_len):
                        is_complete = False
                        break
                if is_complete:
                    complete_alignment_num += 1
        if complete_alignment_num == 1:
            repeats_tobe_spliced[query_name] = 1

    #print(repeats_tobe_spliced)

    for query_name in repeats_tobe_spliced.keys():
        for record in query_records[query_name]:
            reference_name = record[1]
            cigar = record[2]
            cigarstr = record[3]
            is_reverse = record[4]
            query_seq = repeat_contigs[query_name]
            query_len = len(query_seq)
            if is_reverse:
                query_seq = getReverseSequence(query_seq)
            if not all_fragments.__contains__(query_name):
                all_fragments[query_name] = []
            fragments = []

            # parse cigar
            repeat_index = 0
            last_cigar = -1
            frag = ''
            frag_start_pos = 0
            is_split = False
            for c in cigar:
                if last_cigar != c[0]:
                    last_cigar = c[0]
                    # large gap, split repeat
                    if ((c[0] == 4 and c[1] >= HS_gap * query_len)
                                       or (c[0] == 5 and c[1] >= HS_gap * query_len)
                                       or (c[0] == 1 and c[1] >= ID_gap * query_len)
                                       or (c[0] == 2 and c[1] >= ID_gap * query_len)):
                        is_split = True
                        if frag != '':
                            fragments.append((frag_start_pos, len(frag), frag))
                            frag_start_pos += len(frag)
                            if c[0] != 2:
                                frag_start_pos += c[1]
                            frag = ''

                if (c[0] == 4 or c[0] == 5) and c[1] >= HS_gap * query_len:
                    # if cigar is large H/S, repeat index increment
                    repeat_index += c[1]
                    continue
                elif c[0] == 2 or c[0] == 3:
                    # if cigar is D/N, repeat index should stay
                    continue
                elif c[0] == 1 and c[1] >= ID_gap * query_len:
                    # if cigar is large I, repeat index increment
                    repeat_index += c[1]
                    continue
                else:
                    # if cigar is M or small I/D or small H/S, store sequence and repeat index increment
                    frag += query_seq[repeat_index:repeat_index+c[1]]
                    repeat_index += c[1]
            if frag != '':
                fragments.append((frag_start_pos, len(frag), frag))
            old_fragments = all_fragments[query_name]
            all_fragments[query_name] = old_fragments + fragments

    # if keep original repeat
    # de-duplicate reverse-complementarty repeat
    all_unique_fragments = {}
    for query_name in all_fragments.keys():
        frag_seq_set = []
        frag_set = []

        for frag in all_fragments[query_name]:
            if frag[2] not in frag_seq_set and getReverseSequence(frag[2]) not in frag_seq_set:
                frag_seq_set.append(frag[2])
            frag_set.append(frag)
        all_unique_fragments[query_name] = frag_set

    node_index = 0
    with open(raw_cut_file, 'w') as f_save:
        for query_name in all_unique_fragments.keys():
            for unique_frag in all_unique_fragments[query_name]:
                f_save.write('>Node_' + str(node_index) + '-len_' + str(len(unique_frag[2]))
                             + '\n' + unique_frag[2] + '\n')
                node_index += 1

        for query_name in query_records.keys():
            if not repeats_tobe_spliced.__contains__(query_name):
                seq = repeat_contigs[query_name]
                f_save.write('>Node_' + str(node_index) + '-len_' + str(len(seq))
                             + '\n' + seq + '\n')
                node_index += 1


def cut_repeat(sam_paths, HS_gap, ID_gap, repeats_file):
    repeat_contignames, repeat_contigs = read_fasta(repeats_file)
    all_fragments = {}
    #original_repeat_occurrences = {}
    for sam_path in sam_paths:
        samfile = pysam.AlignmentFile(sam_path, "rb")
        for read in samfile.fetch():
            if read.is_unmapped:
                continue
            query_name = read.query_name
            reference_name = read.reference_name
            cigar = read.cigartuples
            cigarstr = read.cigarstring
            is_reverse = read.is_reverse
            query_seq = repeat_contigs[query_name]
            query_len = len(query_seq)
            if is_reverse:
                query_seq = getReverseSequence(query_seq)
            if not all_fragments.__contains__(query_name):
                all_fragments[query_name] = []
            fragments = []

            # parse cigar
            repeat_index = 0
            last_cigar = -1
            frag = ''
            frag_start_pos = 0
            is_split = False
            for c in cigar:
                if last_cigar != c[0]:
                    last_cigar = c[0]
                    # large gap, split repeat
                    if ((c[0] == 4 and c[1] >= HS_gap)
                                       or (c[0] == 5 and c[1] >= HS_gap)
                                       or (c[0] == 1 and c[1] >= ID_gap * query_len)
                                       or (c[0] == 2 and c[1] >= ID_gap * query_len)):
                        is_split = True
                        if frag != '':
                            fragments.append((frag_start_pos, len(frag), frag))
                            frag_start_pos += len(frag)
                            if c[0] != 2:
                                frag_start_pos += c[1]
                            frag = ''

                if (c[0] == 4 or c[0] == 5) and c[1] >= HS_gap:
                    # if cigar is large H/S, repeat index increment
                    repeat_index += c[1]
                    continue
                elif c[0] == 2 or c[0] == 3:
                    # if cigar is D/N, repeat index should stay
                    continue
                elif c[0] == 1 and c[1] >= ID_gap * query_len:
                    # if cigar is large I, repeat index increment
                    repeat_index += c[1]
                    continue
                else:
                    # if cigar is M or small I/D or small H/S, store sequence and repeat index increment
                    frag += query_seq[repeat_index:repeat_index+c[1]]
                    repeat_index += c[1]
            if frag != '':
                fragments.append((frag_start_pos, len(frag), frag))
            old_fragments = all_fragments[query_name]
            all_fragments[query_name] = old_fragments + fragments
        samfile.close()

    # if keep original repeat
    # de-duplicate reverse-complementarty repeat
    all_unique_fragments = {}
    for query_name in all_fragments.keys():
        frag_seq_set = []
        frag_set = []

        for frag in all_fragments[query_name]:
            if frag[2] not in frag_seq_set and getReverseSequence(frag[2]) not in frag_seq_set:
                frag_seq_set.append(frag[2])
            frag_set.append(frag)
        all_unique_fragments[query_name] = frag_set

    return all_unique_fragments

def get_multiple_alignment_repeat(sam_paths):
    not_multi_mapping_repeatIds = []
    multi_mapping_repeatIds = []
    mapping_repeatIds = {}
    for item in sam_paths:
        sam_path = item[0]
        repeats_path = item[1]
        samfile = pysam.AlignmentFile(sam_path, "rb")
        for read in samfile.fetch():
            query_name = read.query_name
            if not mapping_repeatIds.__contains__(query_name):
                mapping_repeatIds[query_name] = 0

            if read.is_unmapped:
                continue
            else:
                count = mapping_repeatIds[query_name]
                count += 1
                mapping_repeatIds[query_name] = count

    for query_name in mapping_repeatIds.keys():
        count = mapping_repeatIds[query_name]
        if count > 1:
            multi_mapping_repeatIds.append(query_name)
        else:
            not_multi_mapping_repeatIds.append(query_name)

    return multi_mapping_repeatIds, not_multi_mapping_repeatIds

def judgeReduceThreads(unique_kmer_path, threads, log):
    file_size = os.path.getsize(unique_kmer_path) / (1024 * 1024 * 1024)
    mem = psutil.virtual_memory()
    free_memory = float(mem.free) / (1024 * 1024 * 1024)

    if free_memory / 4 < file_size * threads:
        reduce_threads = int(free_memory / (4 * file_size))
        log.logger.debug('----Warning:\nDetect the free memory of your machine is %f GB.\n'
              'The kmer.txt file is %f GB, which will be used in each thread.\n'
              'The number of thread you set is %d. According to our experience,\n'
              'To avoid the risk of out of memory, free memory should be more than 4 times\n'
              'higher than thread_num*kmer_size. Thus, reduce the number of thread to %d.\n' % (
              free_memory, file_size, threads, reduce_threads))
        return reduce_threads
    else:
        return threads

def get_alignment_info(sam_paths):
    unmapped_repeatIds = []
    single_mapped_repeatIds = []
    multi_mapping_repeatIds = []
    mapping_repeatIds = {}

    for item in sam_paths:
        sam_path = item[0]
        repeats_path = item[1]
        samfile = pysam.AlignmentFile(sam_path, "rb")
        for read in samfile.fetch():
            query_name = read.query_name
            if not mapping_repeatIds.__contains__(query_name):
                mapping_repeatIds[query_name] = 0

            if read.is_unmapped:
                continue
            else:
                count = mapping_repeatIds[query_name]
                count += 1
                mapping_repeatIds[query_name] = count

    for query_name in mapping_repeatIds.keys():
        count = mapping_repeatIds[query_name]
        if count <= 0:
            unmapped_repeatIds.append(query_name)
        elif count == 1:
            single_mapped_repeatIds.append(query_name)
        else:
            multi_mapping_repeatIds.append(query_name)

    return unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds


def get_alignment_info_v4(sam_paths, repeats_file):
    repeat_contignames, repeat_contigs = read_fasta(repeats_file)

    unmapped_repeatIds = []
    single_mapped_repeatIds = []
    multi_mapping_repeatIds = []
    segmental_duplication_repeatIds = []

    query_records = {}
    for sam_path in sam_paths:
        samfile = pysam.AlignmentFile(sam_path, "rb")
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

            if not query_records.__contains__(query_name):
                query_records[query_name] = []
            records = query_records[query_name]
            records.append((query_name, reference_name, cigar, cigarstr, is_reverse, alignment_len, identity))
            query_records[query_name] = records

    for query_name in query_records.keys():
        records = query_records[query_name]
        freq = len(records)
        if freq == 1:
            single_mapped_repeatIds.append(query_name)
        elif freq <= 4:
            segmental_duplication_repeatIds.append(query_name)
        else:
            multi_mapping_repeatIds.append((query_name, freq))
    return single_mapped_repeatIds, multi_mapping_repeatIds, segmental_duplication_repeatIds

def get_alignment_info_v1(sam_paths, repeats_file):
    repeat_contignames, repeat_contigs = read_fasta(repeats_file)

    unmapped_repeatIds = []
    single_mapped_repeatIds = []
    multi_mapping_repeatIds = []
    segmental_duplication_repeatIds = []

    query_records = {}
    for sam_path in sam_paths:
        samfile = pysam.AlignmentFile(sam_path, "rb")
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

            if not query_records.__contains__(query_name):
                query_records[query_name] = []
            records = query_records[query_name]
            records.append((query_name, reference_name, cigar, cigarstr, is_reverse, alignment_len, identity))
            query_records[query_name] = records

    for query_name in query_records.keys():
        complete_alignment_num = 0
        high_identity_num = 0
        # other cigars all the same (except first one) are regarded as segmental duplication(LCR)
        # is_special_lcr = True
        # last_cigarstr = ''
        # first_cigarstr = ''
        for i, record in enumerate(query_records[query_name]):
            if i == 0:
                continue
            cigar = record[2]
            cigarstr = str(record[3])
            alignment_len = record[5]
            identity = record[6]
            query_seq = repeat_contigs[query_name]
            query_len = len(query_seq)
            if float(alignment_len) / query_len >= 0.9 and identity >= 80:
                complete_alignment_num += 1
                if identity >= 90:
                    high_identity_num += 1

        if complete_alignment_num == 0:
            single_mapped_repeatIds.append(query_name)
        elif complete_alignment_num > 0:
            # low copy number and all of them are high identicial, they are LCR
            # if complete_alignment_num < 4 and (high_identity_num >= len(query_records[query_name])-1):
            if complete_alignment_num < 5 and high_identity_num >= 1:
                segmental_duplication_repeatIds.append(query_name)
            else:
                multi_mapping_repeatIds.append(query_name)
        else:
            unmapped_repeatIds.append(query_name)

    return unmapped_repeatIds, single_mapped_repeatIds, multi_mapping_repeatIds, segmental_duplication_repeatIds

def get_alignment_info_v3(sam_paths, repeats_file):
    repeat_contignames, repeat_contigs = read_fasta(repeats_file)
    mapping_repeatIds = {}
    query_records = {}
    for sam_path in sam_paths:
        samfile = pysam.AlignmentFile(sam_path, "rb")
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
            q_start = int(read.query_alignment_start)
            q_end = int(read.query_alignment_end)
            t_start = int(read.reference_start)
            t_end = int(read.reference_end)

            if not query_records.__contains__(query_name):
                query_records[query_name] = []
            records = query_records[query_name]
            records.append((query_name, reference_name, cigar, cigarstr, is_reverse, alignment_len, identity, t_start, t_end))
            query_records[query_name] = records

    query_position = {}
    for query_name in query_records.keys():
        complete_alignment_num = 0
        high_identity_num = 0
        query_seq = repeat_contigs[query_name]
        query_len = len(query_seq)
        for i, record in enumerate(query_records[query_name]):
            reference_name = record[1]
            cigar = record[2]
            cigarstr = str(record[3])
            alignment_len = record[5]
            identity = record[6]
            t_start = record[7]
            t_end = record[8]
            if float(alignment_len) / query_len >= 0.8 and identity >= 80:
                complete_alignment_num += 1
                if identity >= 90:
                    high_identity_num += 1
            if t_start > t_end:
                tmp = t_end
                t_end = t_start
                t_start = tmp
            if not query_position.__contains__(reference_name):
                query_position[reference_name] = []
            same_chr_seq = query_position[reference_name]
            same_chr_seq.append((query_name, t_start, t_end))
            query_position[reference_name] = same_chr_seq
        mapping_repeatIds[query_name] = (complete_alignment_num, query_len)
    new_mapping_repeatIds = {k: v for k, v in sorted(mapping_repeatIds.items(), key=lambda item: (-item[1][1], -item[1][0]))}

    return new_mapping_repeatIds, query_position

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


def store_fasta(contigs, file_path):
    with open(file_path, 'w') as f_save:
        for name in contigs.keys():
            seq = contigs[name]
            f_save.write('>'+name+'\n'+seq+'\n')


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
    #print(query_records)

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


def extract_tandem_from_trf(trf_data_path):
    tandem_elements = []
    with open(trf_data_path, 'r') as f_r:
        for line in f_r:
            parts = line.split(' ')
            if len(parts) == 15:
                tandem_elements.append(parts[13])
    return tandem_elements


def get_candidate_repeats(reference, chrom_seg_length, k_num, reduce_partitions_num, unique_kmer_map, fault_tolerant_bases, tmp_output_dir):
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

def connect_frags(connected_regions, repeats_path, reference, threads, tools_dir, tmp_output_dir, skip_threshold):
    # Step2: align repeat back to reference, generate frag_pos_dict
    repeat_contignames, repeat_contigs = read_fasta(repeats_path)
    use_align_tools = 'bwa'
    sam_path_bwa = run_alignment(repeats_path, reference, use_align_tools, threads, tools_dir)
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
        t_start = int(read.reference_start) - 1
        t_end = int(read.reference_end) - 1
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

    return connected_frags





if __name__ == '__main__':
    trf_data_path = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/ReferenceMode/output/CRD.2022-05-04.9-30-26/trf/dmel-all-chromosome-r5.43.fasta.2.7.7.80.10.50.500.dat'
    tandem_elements = extract_tandem_from_trf(trf_data_path)
    tandem_path = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/ReferenceMode/output/CRD.2022-05-04.9-30-26/trf/tandem.fa'

    with open(tandem_path, 'w') as f_save:
        for index, elem in enumerate(tandem_elements):
            f_save.write('>Node_'+str(index)+'\n'+elem+'\n')