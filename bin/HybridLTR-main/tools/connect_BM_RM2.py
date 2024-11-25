import argparse
import os

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

def FMEA_new(query_path, file_work_with_file, full_length_threshold):
    # parse file_work_with_file output, determine the repeat boundary
    # query_records = {query_name: {subject_name: [(q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end)] }}
    query_records = {}
    with open(file_work_with_file, 'r') as f_r:
        for idx, line in enumerate(f_r):
            # print('current line idx: %d' % (idx))
            parts = line.split('\t')
            query_name = parts[4]
            subject_name = parts[9]
            divergence = float(parts[1])/100
            q_start = int(parts[5])
            q_end = int(parts[6])
            direct = parts[8]
            if direct == '+':
                s_start = int(parts[11])
                s_end = int(parts[12])
            else:
                s_start = int(parts[12])
                s_end = int(parts[13])
            if query_name == subject_name and q_start == s_start and q_end == s_end:
                continue
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end, divergence))
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
            frag_index = 0
            for k, frag in enumerate(forward_pos):
                is_update = False
                cur_subject_start = frag[2]
                cur_subject_end = frag[3]
                cur_query_start = frag[0]
                cur_query_end = frag[1]
                cur_divergence = frag[4]
                cur_divergence_bases = abs(cur_query_end - cur_query_start) * cur_divergence
                for cur_frag_index in forward_long_frags.keys():
                    cur_frag = forward_long_frags[cur_frag_index]
                    prev_subject_start = cur_frag[2]
                    prev_subject_end = cur_frag[3]
                    prev_query_start = cur_frag[0]
                    prev_query_end = cur_frag[1]
                    prev_divergence = cur_frag[5]
                    prev_divergence_bases = abs(prev_query_end - prev_query_start) * prev_divergence
                    if cur_subject_end > prev_subject_end:
                        # forward extend
                        if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                and cur_subject_start - prev_subject_end < skip_gap:  # \
                            # extend frag
                            prev_query_start = prev_query_start
                            prev_query_end = cur_query_end
                            prev_subject_start = prev_subject_start if prev_subject_start < cur_subject_start else cur_subject_start
                            prev_subject_end = cur_subject_end
                            # 更新divergence
                            update_divergence = float(cur_divergence_bases + prev_divergence_bases)/abs(prev_query_end - prev_query_start)
                            extend_frag = (prev_query_start, prev_query_end, prev_subject_start, prev_subject_end, subject_name, update_divergence)
                            forward_long_frags[cur_frag_index] = extend_frag
                            is_update = True
                if not is_update:
                    forward_long_frags[frag_index] = (cur_query_start, cur_query_end, cur_subject_start, cur_subject_end, subject_name, cur_divergence)
                    frag_index += 1
            longest_queries += list(forward_long_frags.values())

            reverse_long_frags = {}
            frag_index = 0
            for k, frag in enumerate(reverse_pos):
                is_update = False
                cur_subject_start = frag[2]
                cur_subject_end = frag[3]
                cur_query_start = frag[0]
                cur_query_end = frag[1]
                cur_divergence = frag[4]
                cur_divergence_bases = abs(cur_query_end - cur_query_start) * cur_divergence
                for cur_frag_index in reverse_long_frags.keys():
                    cur_frag = reverse_long_frags[cur_frag_index]
                    prev_subject_start = cur_frag[2]
                    prev_subject_end = cur_frag[3]
                    prev_query_start = cur_frag[0]
                    prev_query_end = cur_frag[1]
                    prev_divergence = cur_frag[5]
                    prev_divergence_bases = abs(prev_query_end - prev_query_start) * prev_divergence
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
                            # 更新divergence
                            update_divergence = float(cur_divergence_bases + prev_divergence_bases) / abs(prev_query_end - prev_query_start)
                            extend_frag = (prev_query_start, prev_query_end, prev_subject_start, prev_subject_end, subject_name, update_divergence)
                            reverse_long_frags[cur_frag_index] = extend_frag
                            is_update = True
                if not is_update:
                    reverse_long_frags[frag_index] = (cur_query_start, cur_query_end, cur_subject_start, cur_subject_end, subject_name, cur_divergence)
                    frag_index += 1
            longest_queries += list(reverse_long_frags.values())

        if not longest_repeats.__contains__(query_name):
            longest_repeats[query_name] = []
        cur_longest_repeats = longest_repeats[query_name]
        for repeat in longest_queries:
            # Subject序列处理流程
            divergence = repeat[5] * 100
            subject_name = repeat[4]
            old_subject_start_pos = repeat[2]
            old_subject_end_pos = repeat[3]
            # Query序列处理流程
            old_query_start_pos = repeat[0]
            old_query_end_pos = repeat[1]
            cur_query_seq_len = abs(old_query_end_pos - old_query_start_pos)
            cur_longest_repeats.append((query_name, old_query_start_pos, old_query_end_pos, subject_name, old_subject_start_pos, old_subject_end_pos, divergence))

    return longest_repeats


if __name__ == '__main__':
    describe_info = ''
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--in_file', required=True, metavar='in_file',
                        help='Input file of BM_RM2, name: file_work_with.txt')
    parser.add_argument('--std_lib', required=True, metavar='std_lib',
                        help='standard library')
    parser.add_argument('--test_lib', required=True, metavar='test_lib',
                        help='test library')

    args = parser.parse_args()

    in_file = args.in_file
    std_lib = args.std_lib
    test_lib = args.test_lib

    full_length_threshold = 0.95

    longest_repeats = FMEA_new(test_lib, in_file, full_length_threshold)

    test_names, test_contigs = read_fasta(test_lib)
    std_names, std_contigs = read_fasta(std_lib)
    std_classifications = {}
    for name in std_names:
        parts = name.split('#')
        std_classifications[parts[0]] = parts[1]


    perfect_families = set()
    for query_name in longest_repeats.keys():
        for item in longest_repeats[query_name]:
            query_name = item[0]
            subject_name = item[3]
            label = std_classifications[subject_name]
            query_len = len(test_contigs[query_name])
            subject_len = len(std_contigs[subject_name + '#' + label])
            q_start = item[1]
            q_end = item[2]
            q_len = abs(q_end - q_start)
            query_remains = query_len - q_end

            subject_start = item[4]
            subject_end = item[5]
            divergence = round(item[6], 1)
            s_len = abs(subject_end - subject_start)

            direct = 'C' if subject_start > subject_end else '+'

            if float(q_len) / query_len >= 0.95 and float(s_len) / subject_len >= 0.95:
                if float(min(q_len, s_len)) / max(q_len, s_len) >= 0.95:
                    # print('full_length: ' + query_name)
                    perfect_families.add(subject_name+'#'+label)
                # else:
                #     print('Insertions or Deletions in ' + query_name)
    print('Perfect families: ' + str(len(perfect_families)))
    # print(perfect_families)


            # # 我们连接完之后，不能继续交给 BM_RM2 去获得 perfect families，因为它无法处理嵌合的情况
            # new_line = 'NA' + '\t' + str(divergence) + '\t' + 'NA' + '\t' + 'NA' + '\t' + query_name + \
            #            '\t' + str(q_start) + '\t' + str(q_end) + '\t' + str(query_remains) + \
            #            '\t' + direct + '\t' + subject_name + '\t' + label
            # if direct == '+':
            #     subject_remains = subject_len - subject_end
            #     new_line += '\t' + str(subject_start) + '\t' + str(subject_end) + '\t' + str(subject_remains)
            # else:
            #     subject_remains = subject_len - subject_start
            #     new_line += '\t' + str(subject_remains) + '\t' + str(subject_start) + '\t' + str(subject_end)
            # new_line += '\tNA\tNA'
            # print(new_line)



