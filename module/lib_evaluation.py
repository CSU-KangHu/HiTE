#-- coding: UTF-8 --
# 这个脚本是用来评估不同的 TE library 的性能，它结合了RepeatModeler2的序列全长特征以及EDTA的具体细节指标例如TP,FP,FN。
# 与 lib_evaluation_bak.py不同的是，我们这里对.out文件进行过滤，只保留完整全长的拷贝
import argparse
import os

from Util import read_fasta, get_overlap_len, get_gap_len, merge_same_fragments, get_chr_pos, store_fasta

def transfer_RMOut2BlastnOut(RMOut, BlastnOut, consensus_path, tools_dir, coverage_threshold):
    cons_names, cons_contigs = read_fasta(consensus_path)
    cons_len = {}
    for name in cons_names:
        new_name = name.split('#')[0]
        cons_len[new_name] = len(cons_contigs[name])

    # 1. 将.out文件转.bed文件
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
            # 取全长拷贝
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
    # 将染色体划分成 10k 一段，根据碎片的起始和结束位置将碎片映射到某一段中。例如碎片的起始和结束位置分别是9549和9617，它们对10k取模都是0，因此映射到第0段。
    # 如果碎片映射到了前一段和后一段中间，例如某个碎片的起始和结束位置分别为9645和15966，因为它们对10k取模分别是0和1，因此它们要么映射到第0段，要么第1段。
    # 此时我们计算10k-9645=355,15966-10k=5966,5966>355，因此应该把这条序列映射到第1段。
    segment_len = 100000 # 100K
    # chr_segments -> {chr1: {seg0: [(start, end, status)], seg1: []}} status:0, 1 0代表片段未被标记为found，1表示片段被标记为found
    chr_segments = {}
    total_chr_len = 0
    # 根据染色体的长度将其均分成N段，这样做是将fragment分段存储，减少检索时间
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
    # 将repbase的fragment映射到对应的segment中存储
    for chr_name in repbase_fragments.keys():
        cur_chr_segment = chr_segments[chr_name]
        for frag in repbase_fragments[chr_name]:
            start = frag[0]
            end = frag[1]
            seg_index = map_fragment(start, end, segment_len)
            segment_frags = cur_chr_segment[seg_index]
            segment_frags.append([frag[0], frag[1], 0])
    # 将test的fragment映射到对应的segment中，判断segment中是否有fragment和test fragment有超过95%的overlap
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
            # 判断segment中是否有fragment和test fragment有超过95%的overlap
            is_found = False
            for prev_frag in segment_frags:
                overlap_len = get_overlap_len(prev_frag, cur_frag)
                if overlap_len / abs(prev_frag[1] - prev_frag[0]) >= coverage_threshold and overlap_len / abs(cur_frag[1] - cur_frag[0]) >= coverage_threshold:
                    # 将prev_frag的status改为1
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

    # 遍历 repbase_BlastnOut ,寻找并保存 FN 栏目
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

    # 统计一下当前FN最多的top 5 的TE
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

    # 遍历 test_BlastnOut ,寻找并保存 FP 栏目
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

def evaluation(repbase_BlastnOut, test_BlastnOut, chrom_length):
    repbase_fragments = get_chr_fragments(repbase_BlastnOut)
    test_fragments = get_chr_fragments(test_BlastnOut)
    # 将染色体划分成 10k 一段，根据碎片的起始和结束位置将碎片映射到某一段中。例如碎片的起始和结束位置分别是9549和9617，它们对10k取模都是0，因此映射到第0段。
    # 如果碎片映射到了前一段和后一段中间，例如某个碎片的起始和结束位置分别为9645和15966，因为它们对10k取模分别是0和1，因此它们要么映射到第0段，要么第1段。
    # 此时我们计算10k-9645=355,15966-10k=5966,5966>355，因此应该把这条序列映射到第1段。
    segment_len = 100000 # 100K
    # chr_segments -> {chr1: {seg0: [(start, end, status)], seg1: []}} status:0, 1 0代表片段未被标记为found，1表示片段被标记为found
    chr_segments = {}
    total_chr_len = 0
    # 根据染色体的长度将其均分成N段，这样做是将fragment分段存储，减少检索时间
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
    # 将repbase的fragment映射到对应的segment中存储
    for chr_name in repbase_fragments.keys():
        cur_chr_segment = chr_segments[chr_name]
        for frag in repbase_fragments[chr_name]:
            start = frag[0]
            end = frag[1]
            seg_index = map_fragment(start, end, segment_len)
            segment_frags = cur_chr_segment[seg_index]
            segment_frags.append([frag[0], frag[1], 0])
    # 将test的fragment映射到对应的segment中，判断segment中是否有fragment和test fragment有超过95%的overlap
    TP = 0
    FP = 0
    FN = 0
    target_len = 0
    for chr_name in test_fragments.keys():
        cur_chr_segment = chr_segments[chr_name]
        for cur_frag in test_fragments[chr_name]:
            start = cur_frag[0]
            end = cur_frag[1]
            seg_index = map_fragment(start, end, segment_len)
            segment_frags = cur_chr_segment[seg_index]
            # 判断segment中是否有fragment和test fragment有超过95%的overlap
            is_found = False
            for prev_frag in segment_frags:
                overlap_len = get_overlap_len(prev_frag, cur_frag)
                if overlap_len / abs(prev_frag[1] - prev_frag[0]) >= 0.95 and overlap_len / abs(cur_frag[1] - cur_frag[0]) >= 0.95:
                    # 将prev_frag的status改为1
                    is_found = True
                    prev_frag[2] = 1
                    TP += abs(cur_frag[1] - cur_frag[0])
            if not is_found:
                FP += abs(cur_frag[1] - cur_frag[0])
    for chr_name in chr_segments.keys():
        cur_chr_segment = chr_segments[chr_name]
        for seg_index in cur_chr_segment.keys():
            segment_frags = cur_chr_segment[seg_index]
            for frag in segment_frags:
                if frag[2] == 0:
                    FN += abs(frag[1] - frag[0])
                target_len += abs(frag[1] - frag[0])
    TN = total_chr_len - target_len

    sensitivity = round(float(TP) / (TP + FN), 4)
    specificity = round(float(TN) / (TN + FP), 4)
    accuracy = round(float(TP + TN) / (TP + TN + FP + FN), 4)
    precision = round(float(TP) / (TP + FP), 4)
    F1 = round(float(2 * TP) / (2 * TP + FP + FN), 4)
    FDR = round(float(FP) / (TP + FP), 4)
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

def evaluation_sample(genome_path, repbase_path, repbase_RMOut, test_path, test_RMOut, work_dir, tools_dir):
    # Step 0. 获取基因组的长度
    names, contigs = read_fasta(genome_path)
    chrom_length = {}
    for name in names:
        chr_len = len(contigs[name])
        chrom_length[name] = chr_len

    # Step 1. transform RepeatMasker out format to blastn format
    repbase_BlastnOut = work_dir + '/repbase.blastn.out'
    transfer_RMOut2BlastnOut(repbase_RMOut, repbase_BlastnOut, repbase_path, tools_dir)

    test_BlastnOut = work_dir + '/test.blastn.out'
    transfer_RMOut2BlastnOut(test_RMOut, test_BlastnOut, test_path, tools_dir)
    print('test libray:')
    evaluation(repbase_BlastnOut, test_BlastnOut, chrom_length)

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
