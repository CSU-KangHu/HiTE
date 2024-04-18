#-- coding: UTF-8 --
# This script is used to assess the performance of different TE libraries,
# incorporating both the full-length characteristics from RepeatModeler2 and specific metrics from EDTA such as TP, FP, and FN.
import argparse
import os

from Util import read_fasta, get_overlap_len, get_gap_len, merge_same_fragments, get_chr_pos, store_fasta, get_full_length_copies_from_blastn, multi_process_align, generate_full_length_out, \
    multi_process_align_v1


def get_chr_fragments(BlastnOut):
    chr_fragments = {}
    with open(BlastnOut, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            chr_name = parts[1]
            chr_start = int(parts[8])
            chr_end = int(parts[9])
            if not chr_fragments.__contains__(chr_name):
                chr_fragments[chr_name] = []
            fragments = chr_fragments[chr_name]
            fragments.append((chr_start, chr_end))
    return chr_fragments

def get_FN_evaluation(repbase_BlastnOut, test_BlastnOut, FN_BlastnOut, FP_BlastnOut, chrom_length, coverage_threshold):
    repbase_fragments = get_chr_fragments(repbase_BlastnOut)
    test_fragments = get_chr_fragments(test_BlastnOut)
    # Divide the chromosome into segments of 10k each. Map the fragments to a specific segment based on their start and end positions.
    # For example, if a fragment has start and end positions of 9549 and 9617 respectively,
    # and both values modulo 10k result in 0, then it is mapped to the 0th segment.
    # If a fragment falls between the start and end of two consecutive segments,
    # for instance, with start and end positions of 9645 and 15966 respectively,
    # and their modulo 10k values are 0 and 1, then it could be mapped to either the 0th or 1st segment.
    # In this case, we calculate 10k-9645=355 and 15966-10k=5966.
    # Since 5966 is greater than 355, this sequence should be mapped to the 1st segment.
    segment_len = 100000 # 100K
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
    # Map the fragments from Repbase to the corresponding segment for storage.
    for chr_name in repbase_fragments.keys():
        cur_chr_segments = chr_segments[chr_name]
        for frag in repbase_fragments[chr_name]:
            start = frag[0]
            end = frag[1]
            seg_index = map_fragment(start, end, segment_len)
            segment_frags = cur_chr_segments[seg_index]
            segment_frags.append([frag[0], frag[1], 0])
    # Map the fragments from the test set to the corresponding segment,
    # and check if there is an overlap of over 95% with the fragment in the segment.
    TP = 0
    FP = 0
    FN = 0
    target_len = 0
    FP_set = set()
    for chr_name in test_fragments.keys():
        fragments = test_fragments[chr_name]
        # chr_name = genome_name_dict[chr_name]
        cur_chr_segments = chr_segments[chr_name]
        for cur_frag in fragments:
            start = cur_frag[0]
            end = cur_frag[1]
            seg_index = map_fragment(start, end, segment_len)
            segment_frags = cur_chr_segments[seg_index]
            # Check if there is an overlap of over 95% between the fragment in the segment and the test fragment.
            is_found = False
            for prev_frag in segment_frags:
                overlap_len = get_overlap_len(prev_frag, cur_frag)
                if overlap_len / abs(prev_frag[1] - prev_frag[0]) >= coverage_threshold and overlap_len / abs(cur_frag[1] - cur_frag[0]) >= coverage_threshold:
                    # Change the status of prev_frag to 1.
                    is_found = True
                    prev_frag[2] = 1
                    TP += abs(cur_frag[1] - cur_frag[0])
            if not is_found:
                FP += abs(cur_frag[1] - cur_frag[0])
                FP_set.add((chr_name, start, end))
    FN_set = set()
    for chr_name in chr_segments.keys():
        cur_chr_segments = chr_segments[chr_name]
        for seg_index in cur_chr_segments.keys():
            segment_frags = cur_chr_segments[seg_index]
            for frag in segment_frags:
                if frag[2] == 0:
                    FN += abs(frag[1] - frag[0])
                    FN_set.add((chr_name, frag[0], frag[1]))
                target_len += abs(frag[1] - frag[0])
    TN = total_chr_len - target_len

    # Iterate through repbase_BlastnOut, locate and save the FN column.
    FN_lines = []
    with open(repbase_BlastnOut, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            query_name = parts[0]
            chr_name = parts[1]
            chr_start = int(parts[8])
            chr_end = int(parts[9])
            item = (chr_name, chr_start, chr_end)
            if item in FN_set:
                new_line = query_name + '\t' + chr_name + '\t' + str(chr_start) + '\t' +  str(chr_end) + '\n'
                FN_lines.append(new_line)

    with open(FN_BlastnOut, 'w') as f_save:
        for line in FN_lines:
            f_save.write(line)

    # # Calculate the top 5 TEs with the highest FN counts.
    # top_FN = {}
    # with open(FN_BlastnOut, 'w') as f_save:
    #     for line in FN_lines:
    #         parts = line.split('\t')
    #         TE_name = parts[1]
    #         TE_start = int(parts[7])
    #         TE_end = int(parts[8])
    #         FN_base = abs(TE_end-TE_start)
    #         if not top_FN.__contains__(TE_name):
    #             top_FN[TE_name] = 0
    #         cur_FN = top_FN[TE_name]
    #         cur_FN += FN_base
    #         top_FN[TE_name] = cur_FN
    #         f_save.write(line)
    # sorted_items = sorted(top_FN.items(), key=lambda x: x[1], reverse=True)[:20]
    # # print(sorted_items)

    # Iterate through test_BlastnOut, locate and save the FP column.
    FP_lines = []
    with open(test_BlastnOut, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            query_name = parts[0]
            chr_name = parts[1]
            chr_start = int(parts[8])
            chr_end = int(parts[9])
            item = (chr_name, chr_start, chr_end)
            if item in FP_set:
                new_line = query_name + '\t' + chr_name + '\t' + str(chr_start) + '\t' + str(chr_end) + '\n'
                FP_lines.append(new_line)

    with open(FP_BlastnOut, 'w') as f_save:
        for line in FP_lines:
            f_save.write(line)

    sensitivity = round(float(TP) / (TP + FN), 4)
    # specificity = round(float(TN) / (TN + FP), 4)
    # accuracy = round(float(TP + TN) / (TP + TN + FP + FN), 4)
    precision = round(float(TP) / (TP + FP), 4)
    F1 = round(float(2 * TP) / (2 * TP + FP + FN), 4)
    FDR = round(float(FP) / (TP + FP), 4)
    print('TP:', TP)
    print('FP:', FP)
    # print('TN:', TN)
    print('FN:', FN)
    print('sensitivity:', sensitivity)
    # print('specificity:', specificity)
    # print('accuracy:', accuracy)
    print('precision:', precision)
    print('FDR:', FDR)
    print('F1:', F1)

def get_evaluation_sample(genome_path, standard_lib_out, test_lib_out, work_dir, coverage_threshold):
    # Step 0. 获取基因组的长度
    names, contigs = read_fasta(genome_path)
    chrom_length = {}
    for i, name in enumerate(names):
        chr_len = len(contigs[name])
        chrom_length[name] = chr_len

    print('test libray:')
    FN_BlastnOut = work_dir + '/FN.blastn.out'
    FP_BlastnOut = work_dir + '/FP.blastn.out'
    get_FN_evaluation(standard_lib_out, test_lib_out, FN_BlastnOut, FP_BlastnOut, chrom_length, coverage_threshold)

def map_fragment(start, end, chunk_size):
    start_chunk = start // chunk_size
    end_chunk = end // chunk_size

    if start_chunk == end_chunk:
        return start_chunk
    elif abs(end_chunk * chunk_size - start) < abs(end - end_chunk * chunk_size):
        return end_chunk
    else:
        return start_chunk

def sort_blastn_out(lib_out):
    chr_list = {}
    with open(lib_out, 'r') as f_r:
        for line in f_r:
            cur_tuple = tuple(line.split('\t'))
            chr_name = cur_tuple[1]
            if chr_name not in chr_list:
                chr_list[chr_name] = []
            pos_list = chr_list[chr_name]
            pos_list.append(cur_tuple)
    sorted_lib_out = lib_out + '.sorted'
    with open(sorted_lib_out, 'w') as f_save:
        for chr_name in chr_list.keys():
            pos_list = chr_list[chr_name]
            sorted_pos_list = sorted(pos_list, key=lambda x: (int(x[8]), int(x[9])))
            for item in sorted_pos_list:
                f_save.write('\t'.join(item))
    return sorted_lib_out

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE library evaluation...')
    parser.add_argument('-g', metavar='genome_path',
                        help='e.g., ')
    parser.add_argument('--standard_lib', metavar='standard_lib',
                        help='Path of standard library')
    parser.add_argument('--test_lib', metavar='test_lib',
                        help='Path of test library')
    parser.add_argument('--work_dir', metavar='work_dir',
                        help='work directory')
    parser.add_argument('--thread', metavar='thread_num',
                        help='Input thread num.')
    parser.add_argument('--skip_blastn', metavar='skip_blastn',
                        help='Should the blastn alignment process be skipped? If yes, please ensure that the blastn output file already exists.')
    parser.add_argument('--coverage_threshold', metavar='coverage_threshold',
                        help='Finding thresholds for full-length copies. Three common thresholds: 0.80, 0.95, and 0.99.')
    parser.add_argument('--cat', metavar='TE category',
                        help='TE category: LTR|nonLTR|LINE|SINE|TIR|Helitron|Total')

    args = parser.parse_args()
    genome_path = args.g
    standard_lib = args.standard_lib
    work_dir = args.work_dir
    thread = int(args.thread)
    test_lib = args.test_lib
    coverage_threshold = args.coverage_threshold
    category = args.cat
    skip_blastn = args.skip_blastn

    default_skip_blastn = 0
    default_coverage_threshold = 0.95

    if skip_blastn is not None:
        skip_blastn = int(skip_blastn)
    else:
        skip_blastn = default_skip_blastn

    if coverage_threshold is not None:
        coverage_threshold = float(coverage_threshold)
    else:
        coverage_threshold = default_coverage_threshold

    std_tmp_blast_dir = work_dir + '/std_blastn'
    standard_lib_out = work_dir + '/std_full_length.out'

    test_tmp_blast_dir = work_dir + '/test_blastn'
    test_lib_out = work_dir + '/test_full_length.out'


    if skip_blastn != 1:
        # 由于RepeatMasker会跨过一些大gap，导致比对不准确，影响我们的评价结果，因此我们使用blastn进行比对
        # 此外，我们发现会出现 test library 过大，导致 blastn 比对 out 过大的问题，我们先在单个线程中调用获取全长拷贝的函数，然后再合并所有的全长拷贝
        multi_process_align_v1(standard_lib, genome_path, standard_lib_out, std_tmp_blast_dir, thread, coverage_threshold, category, is_removed_dir=True)
        multi_process_align_v1(test_lib, genome_path, test_lib_out, test_tmp_blast_dir, thread, coverage_threshold, category, is_removed_dir=True)

    tools_dir = os.getcwd() + '/tools'
    get_evaluation_sample(genome_path, standard_lib_out, test_lib_out, work_dir, coverage_threshold)
