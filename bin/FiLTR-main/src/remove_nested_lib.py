import argparse
import os
import sys

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)

from Util import read_fasta, store_fasta, multi_process_align, get_overlap_len, merge_overlap_seq


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
                for cur_frag_index in forward_long_frags.keys():
                    cur_frag = forward_long_frags[cur_frag_index]
                    prev_subject_start = cur_frag[2]
                    prev_subject_end = cur_frag[3]
                    prev_query_start = cur_frag[0]
                    prev_query_end = cur_frag[1]

                    if cur_subject_end > prev_subject_end:
                        # forward extend
                        if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
                                and cur_subject_start - prev_subject_end < skip_gap:  # \
                            # extend frag
                            prev_query_start = prev_query_start
                            prev_query_end = cur_query_end
                            prev_subject_start = prev_subject_start if prev_subject_start < cur_subject_start else cur_subject_start
                            prev_subject_end = cur_subject_end
                            extend_frag = (prev_query_start, prev_query_end, prev_subject_start, prev_subject_end, subject_name)
                            forward_long_frags[cur_frag_index] = extend_frag
                            is_update = True
                if not is_update:
                    forward_long_frags[frag_index] = (cur_query_start, cur_query_end, cur_subject_start, cur_subject_end, subject_name)
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
                for cur_frag_index in reverse_long_frags.keys():
                    cur_frag = reverse_long_frags[cur_frag_index]
                    prev_subject_start = cur_frag[2]
                    prev_subject_end = cur_frag[3]
                    prev_query_start = cur_frag[0]
                    prev_query_end = cur_frag[1]

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

def remove_self_alignment(blastnResults_path, new_blastnResults_path):
    records = []
    with open(blastnResults_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            if query_name == subject_name:
                continue
            else:
                records.append(line)
    f_r.close()

    with open(new_blastnResults_path, 'w') as f_save:
        for line in records:
            f_save.write(line)
    f_save.close()

def remove_nest(blastnResults_path, query_path, subject_path, output, iter_num, coverage = 0.95):
    full_length_threshold = coverage
    longest_repeats = FMEA_new(query_path, blastnResults_path, full_length_threshold)

    minlen = 500
    new_subject_contigs = {}
    query_names, query_contigs = read_fasta(query_path)
    subject_names, subject_contigs = read_fasta(subject_path)

    nested_alignments = {}

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
            if s_start > s_end:
                tmp = s_start
                s_start = s_end
                s_end = tmp

            # 我们认为 比对上的 query 部分应该占它自身的95%以上，同时 比对上的 subject 部分应该不超过自身的 80%，此时我们认为query是subject的一个嵌合元素
            if float(q_len) / query_len >= coverage and float(s_len) / subject_len < 0.8:
                # 保留subject上所有的比对位置，如果有重叠就合并重叠
                if subject_name not in nested_alignments:
                    nested_alignments[subject_name] = []
                keep_alignments = nested_alignments[subject_name]
                is_overlap = False
                cur_frag = (s_start, s_end, [query_name])
                for i, align in enumerate(keep_alignments):
                    prev_frag = align
                    overlap_len = get_overlap_len(prev_frag, cur_frag)
                    if overlap_len > 0:
                        is_overlap = True
                        merge_frag = merge_overlap_seq(prev_frag, cur_frag)
                        keep_alignments[i] = (merge_frag[0], merge_frag[1], prev_frag[2] + cur_frag[2])
                        break
                if not is_overlap:
                    # if subject_name == 'Chr451_8898275-8906382-int#LTR':
                    #     print('h')
                    keep_alignments.append(cur_frag)

    # 现在我们开始去除嵌合元素
    for subject_name in nested_alignments.keys():
        subject_seq = subject_contigs[subject_name]
        align_pos = nested_alignments[subject_name]
        sorted_align_pos = sorted(align_pos, key=lambda x: (x[0], x[1]))
        # 合并有交集的位置
        merge_sorted_align_pos = []
        merge_frag_index = set()
        for i in range(len(sorted_align_pos)):
            cur_frag = sorted_align_pos[i]
            if i in merge_frag_index:
                continue
            for j in range(i+1, len(sorted_align_pos)):
                prev_frag = sorted_align_pos[j]
                overlap_len = get_overlap_len(prev_frag, cur_frag)
                if overlap_len > 0:
                    merge_frag = merge_overlap_seq(prev_frag, cur_frag)
                    cur_frag = (merge_frag[0], merge_frag[1], prev_frag[2] + cur_frag[2])
                    merge_frag_index.add(j)
                else:
                    break
            merge_sorted_align_pos.append(cur_frag)

        # print(subject_name, merge_sorted_align_pos)
        subject_new_seq = ''
        start_pos = 0
        for i, cur_frag in enumerate(merge_sorted_align_pos):
            s_start = cur_frag[0]
            s_end = cur_frag[1]
            subject_new_seq += subject_seq[start_pos: s_start]
            start_pos = s_end
            if i == len(merge_sorted_align_pos) - 1:
                subject_new_seq += subject_seq[s_end:]
        subject_new_len = len(subject_new_seq)
        if subject_new_len >= minlen:
            new_subject_contigs[subject_name] = subject_new_seq

    for subject_name in subject_contigs.keys():
        if subject_name in new_subject_contigs:
            subject_contigs[subject_name] = new_subject_contigs[subject_name]
    store_fasta(subject_contigs, output)


if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run remove nested TEs')
    parser.add_argument('-t', metavar='threads number',
                        help='Input threads number.')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='Please enter the directory for output. Use an absolute path.')
    parser.add_argument('--max_iter_num', metavar='max_iter_num',
                        help='Maximum number of iterations to resolve nested insertions.')
    parser.add_argument('--input1', metavar='input1',
                        help='The path of input1')
    parser.add_argument('--input2', metavar='input2',
                        help='The path of input2')
    parser.add_argument('--output', metavar='output',
                        help='The path of output')

    args = parser.parse_args()

    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    max_iter_num = args.max_iter_num
    input1 = args.input1
    input2 = args.input2
    output = args.output

    default_max_iter_num = 3
    if max_iter_num is None:
        max_iter_num = default_max_iter_num
    else:
        max_iter_num = int(max_iter_num)

    iter_num = 0
    while iter_num < max_iter_num:
        # remove nested TE
        blastnResults_path = tmp_output_dir + '/rm_nested.self.out'
        confident_TE_blast_dir = tmp_output_dir + '/rm_nested_blast'
        multi_process_align(input1, input2, blastnResults_path, confident_TE_blast_dir, threads)
        clean_output = output
        remove_nest(blastnResults_path, input1, input2, clean_output, iter_num, coverage=0.95)
        input1 = clean_output
        input2 = clean_output
        iter_num += 1



