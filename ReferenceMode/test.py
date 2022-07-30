import codecs
import json
import os
import sys
import time

import pysam
from concurrent.futures import ProcessPoolExecutor, as_completed

from Util import convertToUpperCase, read_fasta, getReverseSequence, \
    Logger, split_repeats, compute_identity, run_alignment, multi_line, generate_blastlike_output, \
    get_multiple_alignment_repeat, split2cluster, cut_repeat_v1, judgeReduceThreads, get_ltr_suppl_from_ltrfinder, \
    store_fasta, printClass, parse_ref_blast_output, filter_LTR_high_similarity, get_alignment_info_v3, compare_seq, \
    getRegionCombination, getCombineFragments, convertToUpperCase_v1, generate_candidate_repeats_v2

def test(num):
    for i in range(10000):
        num += num
    return num

# def cut_repeat(sam_path_bwa, repeats_path, cut_repeats_path):
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

def get_longest_repeats(repeats_path, blast_program_dir, tools_dir, extend_base_threshold):
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
            subject_pos.sort(key=lambda x: (x[2], x[3]))

            #cluster all closed fragments
            clusters = {}
            cluster_index = 0
            last_frag = None
            last_direct = '+'
            for k, frag in enumerate(subject_pos):
                if k != 0:
                    if frag[2] > frag[3]:
                        direct = '-'
                    else:
                        direct = '+'
                    if direct != last_direct \
                        or (direct == '+' and frag[2]-last_frag[3] >= extend_base_threshold) \
                        or (direct == '-' and frag[3]-last_frag[2] >= extend_base_threshold):
                        cluster_index += 1
                if not clusters.__contains__(cluster_index):
                    clusters[cluster_index] = []
                cur_cluster = clusters[cluster_index]
                cur_cluster.append(frag)
                last_frag = frag
                # judge subject direction
                if last_frag[2] > last_frag[3]:
                    last_direct = '-'
                else:
                    last_direct = '+'
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
                            if longest_subject_start < longest_subject_end:
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
                            else:
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
                                        cur_longest_query_len = longest_query_start - longest_query_end

                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= extend_base_threshold:
                                        break
                    if cur_longest_query_len > cluster_longest_query_len:
                        cluster_longest_query_start = longest_query_start
                        cluster_longest_query_end = longest_query_end
                        cluster_longest_query_len = cur_longest_query_len
                # keep this longest query
                longest_queries.append((cluster_longest_query_start, cluster_longest_query_end, cluster_longest_query_len))
        # generate fasta, use cd-hit-est to cluster sequences
        local_longest_query_file = single_tmp_dir + '/local_longest_query_' +str(idx) + '.fa'
        l_idx = 0
        with open(local_longest_query_file, 'w') as f_save:
            for item in longest_queries:
                local_longest_seq = split_repeats_contigs[query_name][item[0]-1: item[1]]
                f_save.write('>L_'+str(l_idx)+'\n'+local_longest_seq+'\n')
                l_idx += 1
        local_longest_query_consensus = single_tmp_dir + '/local_longest_query_' +str(idx) + '.cons.fa'
        cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -c ' + str(0.95) + ' -G 0 -g 1 -d 0 -A 80 -i ' + local_longest_query_file + ' -o ' + local_longest_query_consensus
        os.system(cd_hit_command)

        cluster_file = local_longest_query_consensus + '.clstr'
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
                    name = name[0: len(name)-3]
                    cur_cluster.append(name)
                    if line.endswith('*'):
                        clusters['rep_' + str(cluster_idx)] = name

        local_longest_names, local_longest_contigs = read_fasta(local_longest_query_file)
        new_longest_queries = []
        for cluster_idx in clusters.keys():
            if cluster_idx.isdecimal():
                cluster_size = len(clusters[cluster_idx])
                if cluster_size <= 1:
                    continue
                name = clusters['rep_' + str(cluster_idx)]
                new_longest_queries.append((name, cluster_size))
        new_longest_queries.sort(key=lambda x: -x[1])

        sorted_longest_queries = []
        for item in new_longest_queries:
            name = item[0]
            seq = local_longest_contigs[name]
            sorted_longest_queries.append((item[1], seq))

        longest_repeats[query_name] = sorted_longest_queries
    # print(longest_repeats)
    return longest_repeats


def filter_derived_seq(longest_repeats_path, blast_program_dir):
    extend_base_ratio = 0.1

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
                    if query_name == 'N_12064-len_756':
                        print(query_name, subject_name, longest_query_item)
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

def filter_derived_seq_v1(longest_repeats_path, blast_program_dir):

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


if __name__ == '__main__':
    rm2_dir = '/home/hukang/rm2_test/rice/test2'
    perfect = rm2_dir + '/perfect.families'
    perfect_names = set()
    with open(perfect, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            perfect_names.add(line)

    repbase_node_names = {}
    summary_dir = rm2_dir + '/summary_files'
    for name in perfect_names:
        file = summary_dir + '/' + name + '.summary'
        query_records = {}
        with open(file, 'r') as f_r:
            for line in f_r:
                parts = line.split('\t')
                repbase_len = int(parts[1])
                repbase_start = int(parts[2])
                repbase_end = int(parts[3])
                node_name = parts[4]
                if not query_records.__contains__(node_name):
                    query_records[node_name] = []
                records = query_records[node_name]
                records.append((name, repbase_len, repbase_start, repbase_end, node_name))

        perfect_node_found = False
        for node_name in query_records.keys():
            records = query_records[node_name]
            repbase_len = records[0][1]
            repbase_array = ['' for i in range(repbase_len)]
            for record in records:
                repbase_start = record[2]
                repbase_end = record[3]
                #print(query_records[node_name])
                for k in range(repbase_start-1, repbase_end):
                    repbase_array[k] = 'X'

            count = 0
            for i in range(repbase_len):
                if repbase_array[i] == 'X':
                    count += 1
            if float(count)/repbase_len >= 0.95:
                if not repbase_node_names.__contains__(name):
                    repbase_node_names[name] = set()
                node_names = repbase_node_names[name]
                node_names.add(node_name)

                perfect_node_found = True
        if not perfect_node_found:
            print('no node cover %s' %(name))
    #print(repbase_node_names)
    print(len(repbase_node_names))

    #-----------------------------------------------------------

    rm2_dir = '/home/hukang/rm2_test/rice/test'
    perfect = rm2_dir + '/perfect.families'
    perfect_names = set()
    with open(perfect, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            perfect_names.add(line)

    repbase_node_names1 = {}
    summary_dir = rm2_dir + '/summary_files'
    for name in perfect_names:
        file = summary_dir + '/' + name + '.summary'
        query_records = {}
        with open(file, 'r') as f_r:
            for line in f_r:
                parts = line.split('\t')
                repbase_len = int(parts[1])
                repbase_start = int(parts[2])
                repbase_end = int(parts[3])
                node_name = parts[4]
                if not query_records.__contains__(node_name):
                    query_records[node_name] = []
                records = query_records[node_name]
                records.append((name, repbase_len, repbase_start, repbase_end, node_name))

        perfect_node_found = False
        for node_name in query_records.keys():
            records = query_records[node_name]
            repbase_len = records[0][1]
            repbase_array = ['' for i in range(repbase_len)]
            for record in records:
                repbase_start = record[2]
                repbase_end = record[3]
                # print(query_records[node_name])
                for k in range(repbase_start - 1, repbase_end):
                    repbase_array[k] = 'X'

            count = 0
            for i in range(repbase_len):
                if repbase_array[i] == 'X':
                    count += 1
            if float(count) / repbase_len >= 0.95:
                if not repbase_node_names1.__contains__(name):
                    repbase_node_names1[name] = set()
                node_names = repbase_node_names1[name]
                node_names.add(node_name)

                perfect_node_found = True
        if not perfect_node_found:
            print('no node cover %s' % (name))
    #print(repbase_node_names1)
    print(len(repbase_node_names1))

    #--------------------find difference-------------------------
    diff_set = set(repbase_node_names.keys()) - set(repbase_node_names1.keys())
    #print(diff_set)
    diff_names = {}
    for repbase_name in diff_set:
        if not diff_names.__contains__(repbase_name):
            diff_names[repbase_name] = set()
        node_names = diff_names[repbase_name]
        for node_name in repbase_node_names[repbase_name]:
            node_names.add(node_name)
    print(diff_names)

    super_count = {}
    for repbase_name in diff_names.keys():
        super_class = repbase_name.split('-')[0]
        if not super_count.__contains__(super_class):
            super_count[super_class] = 0
        count = super_count[super_class]
        count += 1
        super_count[super_class] = count
    print(super_count)



    # tmp_output_dir = '/home/hukang/KRF_output/dmel/CRD.2022-07-09.16-8-30'
    # blast_program_dir = '/public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2'
    # # ------------------------------------filter derived sequences--------------------------------------------------------
    # # parallel
    # longest_repeats_path = tmp_output_dir + '/longest_repeats.filter_duplication.fa'
    # # longest_repeats_path = tmp_output_dir + '/TE.merge.fa'
    # longest_repeats_path = filter_derived_seq_v1(longest_repeats_path, blast_program_dir)


    # tmp_output_dir = '/home/hukang/rm2_test/test'
    #
    # unique_name = set()
    #
    # perfect_file = tmp_output_dir + '/perfect.families'
    # perfect_set = set()
    # line_count = 0
    # with open(perfect_file, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         perfect_set.add(line)
    #         unique_name.add(line)
    #         line_count += 1
    # print('line count: %d, perfect size: %d' %(line_count, len(perfect_set)))
    #
    # good_file = tmp_output_dir + '/good.families'
    # good_set = set()
    # line_count = 0
    # with open(good_file, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         good_set.add(line)
    #         unique_name.add(line)
    #         line_count += 1
    # print('line count: %d, good size: %d' % (line_count, len(good_set)))
    #
    # present_file = tmp_output_dir + '/present.all.families'
    # present_set = set()
    # line_count = 0
    # with open(present_file, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         present_set.add(line)
    #         unique_name.add(line)
    #         line_count += 1
    # print('line count: %d, present size: %d' % (line_count, len(present_set)))
    #
    # print('unique size: %d' %len(unique_name))


    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-07-01.11-20-22'
    # target_query_name = 'N_5431-len_237'
    # blastn2Results_path = tmp_output_dir + '/longest_repeat.pairwise.out'
    # with open(blastn2Results_path, 'r') as f_r:
    #     for idx, line in enumerate(f_r):
    #         line = line.replace('\n', '')
    #         parts = line.split('\t')
    #         query_name = parts[0]
    #         query_len = int(query_name.split('-len_')[1])
    #         subject_name = parts[1]
    #         subject_len = int(subject_name.split('-len_')[1])
    #         identity = float(parts[2])
    #         alignment_len = int(parts[3])
    #         q_start = int(parts[6])
    #         q_end = int(parts[7])
    #         s_start = int(parts[8])
    #         s_end = int(parts[9])
    #         if query_name == target_query_name:
    #             print(line)



    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-06-30.11-42-3'
    # step1_true_positive = tmp_output_dir + '/true_positive.fa.step1'
    # step2_true_positive = tmp_output_dir + '/true_positive.fa.step2'
    # tp1 = []
    # with open(step1_true_positive, 'r') as f_r:
    #     for line in f_r:
    #         query_name = line.split('\t')[0]
    #         tp1.append(query_name)
    # tp2 = []
    # with open(step2_true_positive, 'r') as f_r:
    #     for line in f_r:
    #         query_name = line.split('\t')[0]
    #         tp2.append(query_name)
    # print(len(tp1), len(tp2), len(tp1)-len(tp2), len(set(tp1) - set(tp2)))
    # print(set(tp1) - set(tp2))

    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-06-09.9-10-58'
    # repeats_path = tmp_output_dir + '/longest_repeats.cons.fa'
    # repeatNames, repeatContigs = read_fasta(repeats_path)
    # sorted_repeatContigs = {k: v for k, v in sorted(repeatContigs.items(), key=lambda item: -len(item[1]))}
    #
    # print(len(list(sorted_repeatContigs.items())[0][1]))

    # repeat_dict_file = tmp_output_dir + '/repeat_dict.csv'
    # file = open(repeat_dict_file, 'r')
    # js = file.read()
    # repeat_dict = json.loads(js)
    #
    # max_len = 0
    # for ref_name in repeat_dict.keys():
    #     repeat_list = repeat_dict[ref_name]
    #     for repeat_item in repeat_list:
    #         start_pos = repeat_item[0]
    #         end_pos = repeat_item[1]
    #         repeat_str = repeat_item[2]
    #         cur_len = len(repeat_str)
    #         if cur_len > max_len:
    #             max_len = cur_len
    # print(max_len)


    # identity_threshold = 0.95
    # length_threshold = 0.95
    # out_file = tmp_dir + '/longest_repeats.out'
    # query_records = {}
    # with open(out_file, 'r') as f_r:
    #     for idx, line in enumerate(f_r):
    #         parts = line.split('\t')
    #         query_name = parts[0]
    #         query_len = int(query_name.split('-')[1][4:])
    #         subject_name = parts[1]
    #         identity = float(parts[2])
    #         alignment_len = int(parts[3])
    #         q_start = int(parts[6])
    #         q_end = int(parts[7])
    #         s_start = int(parts[8])
    #         s_end = int(parts[9])
    #         if not query_records.__contains__(query_name):
    #             query_records[query_name] = []
    #         records = query_records[query_name]
    #         records.append((query_len, alignment_len, identity))
    #
    # false_positive_names = []
    # false_positive_file = tmp_dir + '/false_positive.fa'
    # with open(false_positive_file, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         false_positive_names.append(line)
    #
    # SD_names = []
    # for query_name in query_records.keys():
    #     is_SD = False
    #     for record in records:
    #         if float(record[1])/record[0] >= length_threshold and record[2] >= identity_threshold:
    #             is_SD = False
    #     if is_SD:
    #         SD_names.append(query_name)
    #
    # false_positive_filter_SD = set(false_positive_names)-set(SD_names)
    # print('original false_positive size: %d, segmental duplication size: %d,filter SD false_positive size: %d'
    #       %(len(false_positive_names), len(SD_names), len(false_positive_filter_SD)))



    # single_tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-06-09.9-10-58/tmp_blast/0'
    # repeats_path = (single_tmp_dir+'/repeats_split.fa', single_tmp_dir + '/repeats.fa', single_tmp_dir + '/repeat.pairwise.out')
    # blast_program_dir = '/public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2'
    # tools_dir = os.getcwd() + '/tools'
    # extend_base_threshold = 100
    # longest_repeats = get_longest_repeats(repeats_path, blast_program_dir, tools_dir, extend_base_threshold)
    # print(longest_repeats)


    # merged_out = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-06-09.9-10-58/repeat.pairwise.out'
    # test_out = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-06-09.9-10-58/test.out'
    # os.system('rm -f ' + merged_out)
    # tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-06-09.9-10-58/tmp_blast'
    # for i in range(48):
    #     single_dir = tmp_dir + '/' + str(i)
    #     os.system('cat ' + single_dir + '/repeat.pairwise.out >> ' +merged_out)

    # query_records = []
    # with open(merged_out, 'r') as f_r:
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
    #         if query_name == 'N16925-s_2R-2179831-2199855':
    #             query_records.append(line)
    # with open(test_out, 'w') as f_save:
    #     for line in query_records:
    #         f_save.write(line + '\n')

    # A_set = {}
    # item = (10, 20, 10, 5)
    # item1 = (10, 20, 10, 5)
    # A_set[item] = 1
    # print(A_set.__contains__(item1))


    # reference = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/repbase/zebrep.ref'
    # refNames, refContigs = read_fasta(reference)
    # refContigs = {k: v for k, v in
    #                          sorted(refContigs.items(), key=lambda item: -len(item[1]))}
    #
    # sorted_reference = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/repbase/zebrep.sorted.ref'
    # with open(sorted_reference, 'w') as f_save:
    #     for name in refContigs.keys():
    #         ref_seq = refContigs[name]
    #         f_save.write('>'+name+'\n'+ref_seq+'\n')


    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-21.0-9-18/region_combination_tmp'
    # file = tmp_output_dir + '/rc_0.csv'
    # region_combination_size = float(sys.getsizeof(region_combination))

    # ex = ProcessPoolExecutor(48)
    # bigint = 1024*1024*1024
    # big_list = list(x for x in range(bigint))
    #
    # MAX_JOBS_IN_QUEUE = 500
    # jobs_left = len(big_list)
    # jobs_iter = iter(big_list)
    # jobs = {}
    # while jobs_left:
    #     for num in jobs_iter:
    #         job = ex.submit(test, num)
    #         jobs[job] = 1
    #         if len(jobs) > MAX_JOBS_IN_QUEUE:
    #             break  # limit the job submission for now job
    #
    #     for job in as_completed(jobs):
    #         jobs_left -= 1
    #         num = job.result()
    #         del jobs[job]
    #         break
    # ex.shutdown(wait=True)
    # param_config_path = os.getcwd() + "/ParamConfig.json"
    # # read param config
    # with open(param_config_path, 'r') as load_f:
    #     param = json.load(load_f)
    # load_f.close()

    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/dmel/CRD.2022-05-26.10-28-8'
    # reference = '/public/home/hpc194701009/Ref/dmel-all-chromosome-r5.43.fasta'
    #
    # (ref_dir, ref_filename) = os.path.split(reference)
    # (ref_name, ref_extension) = os.path.splitext(ref_filename)
    # output_dir = tmp_output_dir
    # repeats_consensus = tmp_output_dir + '/repeats.fa'
    # unique_kmer_path = tmp_output_dir + '/kmer.txt'
    # skip_threshold = 200
    # identity_threshold = 0.90
    # length_similarity_cutoff = 0.90
    # tandem_region_cutoff = 0.5
    # k_num = 31
    # threads = 48
    # partitions_num = threads
    # tools_dir = os.getcwd() + '/tools'
    # alias = 'dmel'
    # # chrom_seg_length = int(param['chrom_seg_length'])
    # # fault_tolerant_bases = 100
    #
    # reference = '/public/home/hpc194701009/KmerRepFinder_test/genome/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna'
    # repbase = '/public/home/hpc194701009/repeat_detect_tools/RepeatMasker-4.1.2/RepBase26.05.fasta/drorep.ref'
    # use_align_tools = 'bwa'
    # sam_path_bwa = run_alignment(repbase, reference, use_align_tools, threads, tools_dir)

#
#
#     # Step1: generate candidate repeat regions
#     unique_kmer_map = {}
#     with open(unique_kmer_path, 'r') as f_r:
#         for line in f_r:
#             line = line.replace('\n', '')
#             kmer = line.split(' ')[0]
#             r_kmer = getReverseSequence(kmer)
#             unique_key = kmer if kmer < r_kmer else r_kmer
#             if unique_key.__contains__('N'):
#                 continue
#             unique_kmer_map[unique_key] = 1
#
#     # using multiple threads to gain speed
#     reference_pre = convertToUpperCase_v1(reference)
#     reference_tmp = multi_line(reference_pre, chrom_seg_length, k_num)
#
#     segments = []
#     with open(reference_tmp, 'r') as f_r:
#         for line in f_r:
#             line = line.replace('\n', '')
#             segments.append(line)
#     segments_cluster = split2cluster(segments, partitions_num)
#
#     ex = ProcessPoolExecutor(partitions_num)
#     repeat_dict = {}
#     jobs = []
#     for partiton_index in segments_cluster.keys():
#         cur_segments = segments_cluster[partiton_index]
#         job = ex.submit(generate_candidate_repeats_v2, cur_segments, k_num, unique_kmer_map, partiton_index,
#                         fault_tolerant_bases)
#         jobs.append(job)
#     ex.shutdown(wait=True)
#
#     for job in as_completed(jobs):
#         cur_repeat_dict = job.result()
#         for ref_name in cur_repeat_dict.keys():
#             parts = ref_name.split('$')
#             true_ref_name = parts[0]
#             start_pos = int(parts[1])
#             if not repeat_dict.__contains__(true_ref_name):
#                 repeat_dict[true_ref_name] = []
#             new_repeat_list = repeat_dict[true_ref_name]
#             cur_repeat_list = cur_repeat_dict[ref_name]
#             for repeat_item in cur_repeat_list:
#                 new_repeat_item = (start_pos + repeat_item[0], start_pos + repeat_item[1], repeat_item[2])
#                 new_repeat_list.append(new_repeat_item)
#     for ref_name in repeat_dict.keys():
#         repeat_list = repeat_dict[ref_name]
#         repeat_list.sort(key=lambda x: (x[1], x[2]))
#
#     repeats_path = tmp_output_dir + '/repeats.fa'
#     node_index = 0
#     with open(repeats_path, 'w') as f_save:
#         for ref_name in repeat_dict.keys():
#             repeat_list = repeat_dict[ref_name]
#             for repeat_item in repeat_list:
#                 start_pos = repeat_item[0]
#                 end_pos = repeat_item[1]
#                 query_name = 'N' + str(node_index) + '-s_' + str(ref_name) + '-' + str(start_pos) + '-' + str(end_pos)
#                 repeat = repeat_item[2]
#                 f_save.write('>' + query_name + '\n' + repeat + '\n')
#                 node_index += 1
#
#      # store repeat_dict for testing
#     repeat_dict_file = tmp_output_dir + '/repeat_dict.csv'
#     with codecs.open(repeat_dict_file, 'w', encoding='utf-8') as f:
#         json.dump(repeat_dict, f)
#
#     # Step2: ensure repeats boundary
#     # round-1
#     repeats_path = tmp_output_dir + '/repeats.fa'
#     use_align_tools = 'bwa'
#     sam_path_bwa = run_alignment(repeats_path, reference, use_align_tools, threads, tools_dir)
# #sam_path_bwa = tmp_output_dir + '/repeats.sam'
#     cut_repeats_path = tmp_output_dir + '/repeats.cut.fa'
#     cut_repeat(sam_path_bwa, repeats_path, cut_repeats_path)
#
#     # Step3: merge redundant sequences
#     cut_repeats_consensus = tmp_output_dir + '/repeats.cut.cons.fa'
#     cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -c ' + str(0.95) + ' -i ' + cut_repeats_path + ' -o ' + cut_repeats_consensus + ' -T 0 -M 0'
# #log.logger.debug(cd_hit_command)
#     os.system(cd_hit_command)
#
#
#     # Step4: get multiple alignment sequences
#     cut_repeats_consensus = tmp_output_dir + '/repeats.cut.cons.fa'
#     use_align_tools = 'bwa'
#     sam_path_bwa = run_alignment(cut_repeats_consensus, reference, use_align_tools, threads, tools_dir)
# #sam_path_bwa = tmp_output_dir + '/repeats.cut.cons.sam'
#     sam_paths = []
#     sam_paths.append(sam_path_bwa)
#     new_mapping_repeatIds, query_position = get_alignment_info_v3(sam_paths, cut_repeats_consensus)
#     repeat_multiple_path = tmp_output_dir + '/repeats.cut.cons.multiple.fa'
#     cut_repeat_contigNames, cut_repeat_contigs = read_fasta(cut_repeats_consensus)
#     node_index = 0
#     with open(repeat_multiple_path, 'w') as f_save:
#         for repeat_id in new_mapping_repeatIds.keys():
#             freq = new_mapping_repeatIds[repeat_id][0]
#             seq = cut_repeat_contigs[repeat_id]
#             #f_save.write('>' + repeat_id + '\n' + seq + '\n')
#             #f_save.write('>' + repeat_id + '\tcopies=' + str(freq) + '\n' + seq + '\n')
#             f_save.write('>N' + str(node_index) + '\n' + seq + '\n')
#             node_index += 1

#     repeat_multiple_path = tmp_output_dir + '/repeats.cut.cons.multiple.fa'
#     starttime = time.time()
#     # Step0. use RepeatMasker/trf to mask all low complexity/tandem repeats in raw repeat region
#     # >= tandem_region_cutoff region of the whole repeat region, then it should be filtered, since maybe false positive
#     RepeatMasker_Home = param['RepeatMasker_Home']
#     RepeatMasker_output_dir = tmp_output_dir + '/noint'
#     RepeatMasker_command = 'cd ' + tmp_output_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(threads) \
#                            + ' -noint -x -dir ' + RepeatMasker_output_dir + ' ' + repeat_multiple_path
# #os.system('rm -rf ' + RepeatMasker_output_dir)
# #log.logger.debug(RepeatMasker_command)
#     #os.system(RepeatMasker_command)
#
#     (repeat_multiple_dir, repeat_multiple_filename) = os.path.split(repeat_multiple_path)
#     (repeat_multiple_name, repeat_multiple_extension) = os.path.splitext(repeat_multiple_filename)
#     trf_masked_repeats = RepeatMasker_output_dir + '/' + repeat_multiple_filename + '.masked'
#     trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
#     repeats_contigNames, repeats_contigs = read_fasta(repeat_multiple_path)
#     repeats_path = tmp_output_dir + '/repeats.filter_tandem.fa'
#     with open(repeats_path, 'w') as f_save:
#         for name in trf_contigNames:
#             seq = trf_contigs[name]
#             if float(seq.count('X')) / len(seq) < tandem_region_cutoff:
#                 f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')
#     RepeatMasker_Home = param['RepeatMasker_Home']
#     RepeatMasker_output_dir = tmp_output_dir + '/noint'
#     RepeatMasker_command = 'cd ' + tmp_output_dir + ' && ' + RepeatMasker_Home + '/RepeatMasker -parallel ' + str(
#         threads) + ' -lib ' + repeats_path + ' -dir ' + RepeatMasker_output_dir + ' ' + reference
#     os.system(RepeatMasker_command)

#     TRF_Path = param['TRF_Path']
#
#     trf_dir = tmp_output_dir + '/trf_temp'
#     if not os.path.exists(trf_dir):
#         os.makedirs(trf_dir)
#
#     trf_command = 'cd ' + trf_dir + ' && ' + TRF_Path + ' ' + reference + ' 2 7 7 80 10 50 500 -f -d -m'
# #log.logger.debug(trf_command)
#     os.system(trf_command)
#     trf_masked_repeats = trf_dir + '/' + ref_filename + '.2.7.7.80.10.50.500.mask'

    # trf_contigNames, trf_contigs = read_fasta(trf_masked_repeats)
    # repeats_contigNames, repeats_contigs = read_fasta(merge_pure_consensus)
    # repeats_path = tmp_output_dir + '/repeats.filter_tandem.fa'
    # with open(repeats_path, 'w') as f_save:
    #     for name in trf_contigNames:
    #         seq = trf_contigs[name]
    #         if float(seq.count('N')) / len(seq) < tandem_region_cutoff:
    #             f_save.write('>' + name + '\n' + repeats_contigs[name] + '\n')

    # endtime = time.time()
    # dtime = endtime - starttime
#log.logger.debug("Step0: use trf to mask genome: %.8s s" % (dtime))











    #tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_git/KmerRepFinder/GenomeSimulator/10M_low_freq_out/krf_output/CRD.2022-05-25.20-12-12'

    # Step3: According to valid paths, the fragments in the region are connected across gaps.
    # Fragments in region are connected according to the longer path in valid paths.


    # connected_frags_file = tmp_output_dir + '/connected_frags.csv'
    # file = open(connected_frags_file, 'r')
    # js = file.read()
    # connected_frags = json.loads(js)
    #
    # refNames, refContigs = read_fasta(reference)
    # repeats_connected_file = tmp_output_dir + '/repeats_connected.fa'
    # repeats_connected = {}
    # index = 0
    # for region_index in connected_frags.keys():
    #     for connected_frag in connected_frags[region_index]:
    #         frag_name = connected_frag[0].split(',')[0]
    #         ref_name = frag_name.split('-s_')[1].split('-')[0]
    #         query_name = 'R' + str(index) + '-' + frag_name
    #         seq = refContigs[ref_name][connected_frag[1]: connected_frag[2] + 1]
    #         index += 1
    #         repeats_connected[query_name] = seq
    # sorted_repeats_connected = {k: v for k, v in sorted(repeats_connected.items(), key=lambda item: -len(item[1]))}
    # store_fasta(sorted_repeats_connected, repeats_connected_file)