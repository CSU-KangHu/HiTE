#-- coding: UTF-8 --
# This script is used to assess the performance of different TE libraries,
# incorporating both the full-length characteristics from RepeatModeler2 and specific metrics from EDTA such as TP, FP, and FN.
import argparse
import os

from Util import read_fasta, get_overlap_len, get_gap_len, merge_same_fragments, get_chr_pos, store_fasta

def transfer_RMOut2BlastnOut(RMOut, BlastnOut, consensus_path, tools_dir, coverage_threshold):
    cons_names, cons_contigs = read_fasta(consensus_path)
    cons_len = {}
    for name in cons_names:
        new_name = name.split('#')[0]
        cons_len[new_name] = len(cons_contigs[name])
    # 1. Convert the .out file to a .bed file.
    convert2bed_command = 'perl ' + tools_dir + '/RMout_to_bed.pl ' + RMOut + ' base1'
    print(convert2bed_command)
    os.system(convert2bed_command)
    bed_file = RMOut + '.bed'
    lines = []
    with open(bed_file, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            query_name = parts[0]
            q_start = parts[1]
            q_end = parts[2]
            subject_info = parts[3]
            subject_pats = subject_info.split(';')
            direction = subject_pats[8]
            subject_name = subject_pats[9]
            if direction == '+':
                s_start = subject_pats[11]
                s_end = subject_pats[12]
            else:
                s_start = subject_pats[12]
                s_end = subject_pats[13]
            # Retrieve the full-length copy.
            if float(abs(int(s_end)-int(s_start)))/cons_len[subject_name] >= coverage_threshold:
                new_line = query_name+'\t'+subject_name+'\t'+'-1'+'\t'+'-1'+'\t'+'-1'+'\t'+'-1'+'\t'+q_start+'\t'+q_end+'\t'+s_start+'\t'+s_end+'\t'+'-1'+'\t'+'-1'+'\n'
                lines.append(new_line)
    with open(BlastnOut, 'w') as f_save:
        for line in lines:
            f_save.write(line)


def get_chr_fragments(BlastnOut):
    chr_fragments = {}
    with open(BlastnOut, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            chr_name = parts[0]
            chr_start = int(parts[6])
            chr_end = int(parts[7])
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
        cur_chr_segment = chr_segments[chr_name]
        num_segments = chr_len // segment_len
        if chr_len % segment_len != 0:
            num_segments += 1
        for i in range(num_segments):
            cur_chr_segment[i] = []
    # Map the fragments from Repbase to the corresponding segment for storage.
    for chr_name in repbase_fragments.keys():
        cur_chr_segment = chr_segments[chr_name]
        for frag in repbase_fragments[chr_name]:
            start = frag[0]
            end = frag[1]
            seg_index = map_fragment(start, end, segment_len)
            segment_frags = cur_chr_segment[seg_index]
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
        cur_chr_segment = chr_segments[chr_name]
        for cur_frag in fragments:
            start = cur_frag[0]
            end = cur_frag[1]
            seg_index = map_fragment(start, end, segment_len)
            segment_frags = cur_chr_segment[seg_index]
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
        cur_chr_segment = chr_segments[chr_name]
        for seg_index in cur_chr_segment.keys():
            segment_frags = cur_chr_segment[seg_index]
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
            chr_name = parts[0]
            chr_start = int(parts[6])
            chr_end = int(parts[7])
            TE_name = parts[1]
            TE_start = int(parts[8])
            TE_end = int(parts[9])
            item = (chr_name, chr_start, chr_end)
            if item in FN_set:
                FN_lines.append(line)

    # Calculate the top 5 TEs with the highest FN counts.
    top_FN = {}
    with open(FN_BlastnOut, 'w') as f_save:
        for line in FN_lines:
            parts = line.split('\t')
            TE_name = parts[1]
            TE_start = int(parts[7])
            TE_end = int(parts[8])
            FN_base = abs(TE_end-TE_start)
            if not top_FN.__contains__(TE_name):
                top_FN[TE_name] = 0
            cur_FN = top_FN[TE_name]
            cur_FN += FN_base
            top_FN[TE_name] = cur_FN
            f_save.write(line)
    sorted_items = sorted(top_FN.items(), key=lambda x: x[1], reverse=True)[:5]
    #print(sorted_items)

    # Iterate through test_BlastnOut, locate and save the FP column.
    FP_lines = []
    with open(test_BlastnOut, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            chr_name = parts[0]
            chr_start = int(parts[6])
            chr_end = int(parts[7])
            TE_name = parts[1]
            TE_start = int(parts[8])
            TE_end = int(parts[9])
            item = (chr_name, chr_start, chr_end)
            if item in FP_set:
                FP_lines.append(line)

    with open(FP_BlastnOut, 'w') as f_save:
        for line in FP_lines:
            f_save.write(line)

    sensitivity = round(float(TP) / (TP + FN), 4)
    specificity = round(float(TN) / (TN + FP), 4)
    accuracy = round(float(TP + TN) / (TP + TN + FP + FN), 4)
    precision = round(float(TP) / (TP + FP), 4)
    F1 = round(float(2 * TP) / (2 * TP + FP + FN), 4)
    FDR = round(float(FP) / (TP + FP), 4)
    print('TP:', TP)
    print('FP:', FP)
    print('TN:', TN)
    print('FN:', FN)
    print('sensitivity:', sensitivity)
    print('specificity:', specificity)
    print('accuracy:', accuracy)
    print('precision:', precision)
    print('FDR:', FDR)
    print('F1:', F1)

def get_evaluation_sample(genome_path, repbase_path, repbase_RMOut, test_path, test_RMOut, work_dir, tools_dir, coverage_threshold):
    # Step 0. 获取基因组的长度
    names, contigs = read_fasta(genome_path)
    chrom_length = {}
    for i, name in enumerate(names):
        chr_len = len(contigs[name])
        chrom_length[name] = chr_len

    # Step 1. transform RepeatMasker out format to blastn format
    repbase_BlastnOut = work_dir + '/repbase.blastn.out'
    transfer_RMOut2BlastnOut(repbase_RMOut, repbase_BlastnOut, repbase_path, tools_dir, coverage_threshold)

    test_BlastnOut = work_dir + '/test.blastn.out'
    transfer_RMOut2BlastnOut(test_RMOut, test_BlastnOut, test_path, tools_dir, coverage_threshold)
    print('test libray:')

    FN_BlastnOut = work_dir + '/FN.blastn.out'
    FP_BlastnOut = work_dir + '/FP.blastn.out'
    get_FN_evaluation(repbase_BlastnOut, test_BlastnOut, FN_BlastnOut, FP_BlastnOut, chrom_length, coverage_threshold)

def map_fragment(start, end, chunk_size):
    start_chunk = start // chunk_size
    end_chunk = end // chunk_size

    if start_chunk == end_chunk:
        return start_chunk
    elif abs(end_chunk * chunk_size - start) < abs(end - end_chunk * chunk_size):
        return end_chunk
    else:
        return start_chunk

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE library evaluation...')
    parser.add_argument('-g', metavar='genome_path',
                        help='e.g., ')
    parser.add_argument('--standard_lib', metavar='standard_lib',
                        help='Path of standard library')
    parser.add_argument('--standard_lib_out', metavar='standard_lib_out',
                        help='e.g., Path of standard library .out file')
    parser.add_argument('--test_lib', metavar='test_lib',
                        help='Path of test library')
    parser.add_argument('--test_lib_out', metavar='test_lib_out',
                        help='e.g., Path of test library .out file')
    parser.add_argument('--work_dir', metavar='work_dir',
                        help='work directory')
    parser.add_argument('--coverage_threshold', metavar='coverage_threshold',
                        help='coverage threshold')

    args = parser.parse_args()
    genome_path = args.g
    standard_lib = args.standard_lib
    standard_lib_out = args.standard_lib_out
    work_dir = args.work_dir
    test_lib = args.test_lib
    test_lib_out = args.test_lib_out
    coverage_threshold = args.coverage_threshold

    default_coverage_threshold = 0.95

    if coverage_threshold is not None:
        coverage_threshold = float(coverage_threshold)
    else:
        coverage_threshold = default_coverage_threshold

    tools_dir = os.getcwd() + '/tools'
    get_evaluation_sample(genome_path, standard_lib, standard_lib_out, test_lib, test_lib_out, work_dir, tools_dir, coverage_threshold)
