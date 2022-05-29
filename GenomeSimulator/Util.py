import multiprocessing
import os

import logging
import re
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
                seg = contig[start:start+line_len+(k_num-1)]
                line = '>' + contigName + '\t' + str(start) + '\t' + seg + '\n'
                f_w.write(line)
                start += line_len
    f_w.close()
    return tmp_fasta_path

def convertToUpperCase(reference):
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
                contigname = contigname.split("\t")[0]
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
