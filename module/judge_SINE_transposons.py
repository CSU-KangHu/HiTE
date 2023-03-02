import argparse
import codecs
import os
import sys

import json

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta
from main import get_copies, multi_process_align, flank_region_align, flank_region_align_v1, multi_process_LINE

def cut_reference(fasta_path, line_len):
    tmp_fasta_path = fasta_path + ".helitron.tmp"
    contigNames, contigs = read_fasta(fasta_path)
    with open(tmp_fasta_path, 'w') as f_w:
        for contigName in contigNames:
            contig = contigs[contigName]
            start = 0
            end = len(contig)
            while start < end:
                seg = contig[start:start+line_len]
                line = '>' + contigName + '-' + str(start) + '\n' + seg + '\n'
                f_w.write(line)
                start += line_len
    f_w.close()
    return tmp_fasta_path

if __name__ == '__main__':
    # 1.parse args
    parser = argparse.ArgumentParser(description='run HiTE...')
    parser.add_argument('-i', metavar='input',
                        help='input longest_repeats.fa')
    parser.add_argument('-g', metavar='Genome assembly',
                        help='input genome assembly path')
    parser.add_argument('-t', metavar='threads number',
                        help='input threads number')
    parser.add_argument('--tmp_output_dir', metavar='tmp_output_dir',
                        help='e.g., /public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/dmel')
    parser.add_argument('--blast_program_dir', metavar='blast_program_dir',
                        help='e.g., /public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2')
    parser.add_argument('--ref_index', metavar='ref_index',
                        help='e.g., 0')

    args = parser.parse_args()

    longest_repeats_path = args.i
    reference = args.g
    threads = int(args.t)
    tmp_output_dir = args.tmp_output_dir
    blast_program_dir = args.blast_program_dir
    ref_index = args.ref_index

    candidate_SINE_path = tmp_output_dir + '/candidate_sine_'+str(ref_index)+'.fa'
    blastnResults_path = tmp_output_dir + '/candidate_sine_blast.out'
    temp_dir = tmp_output_dir + '/sine_tmp'
    multi_process_LINE(longest_repeats_path, library_path, candidate_LINE_path, blastnResults_path, blast_program_dir, temp_dir, threads)

    # candidate_LINE_path = tmp_output_dir + '/candidate_line_'+str(ref_index)+'.fa'
    # contignames, contigs = read_fasta(candidate_LINE_path)
    #
    # confident_copies = {}
    # if len(contigs) > 0:
    #     candidate_line_rename = tmp_output_dir + '/candidate_line_'+str(ref_index)+'.rename.fa'
    #     node_index = 0
    #     with open(candidate_line_rename, 'w') as f_save:
    #         for name in contignames:
    #             seq = contigs[name]
    #             new_name = 'LINE_'+str(node_index)
    #             f_save.write('>'+new_name+'\n'+seq+'\n')
    #             node_index += 1
    #
    #     candidate_line_names, candidate_line_contigs = read_fasta(candidate_line_rename)
    #     blastnResults_path = tmp_output_dir + '/line.ref.out'
    #     line_temp_dir = tmp_output_dir + '/line_blast'
    #     multi_process_align(candidate_line_rename, reference, blastnResults_path, blast_program_dir, line_temp_dir, threads)
    #     all_copies = get_copies(blastnResults_path, candidate_line_rename)
    #
    #     flanking_len = 50
    #     # query_name: [(subject_name, subject_start, subject_end, query[2], direct)]
    #     #copy -> (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq, tsd)
    #     ref_names, ref_contigs = read_fasta(reference)
    #     for name in all_copies.keys():
    #         copies = all_copies[name]
    #         if not confident_copies.__contains__(name):
    #             confident_copies[name] = []
    #         copy_list = confident_copies[name]
    #         for copy in copies:
    #             ref_name = copy[0]
    #             copy_ref_start = int(copy[1])
    #             copy_ref_end = int(copy[2])
    #             copy_len = copy_ref_end - copy_ref_start + 1
    #             copy_seq = ref_contigs[ref_name][copy_ref_start - 1 - flanking_len: copy_ref_end + flanking_len]
    #             copy_list.append((copy[0], copy[1], copy[2], copy_len, copy_seq, ''))
    #
    #     # 我们将copies扩展50bp，一个orig_query_name对应一个文件，然后做自比对。
    #     # 解析每个自比对文件，判断C0与C1,C2...等拷贝的比对情况，如果有flanking区域包含在比对区域内，那么这条拷贝应该被抛弃，如果所有拷贝被抛弃，则该条序列应该是假阳性。
    #     flank_align_dir = tmp_output_dir + '/line_flank_align'
    #     confident_copies = flank_region_align_v1(confident_copies, flanking_len, reference, flank_align_dir, blast_program_dir, threads)
    #
    # confident_line_path = tmp_output_dir + '/confident_line_'+str(ref_index)+'.fa'
    # confident_line = {}
    # for name in confident_copies.keys():
    #     copy_list = confident_copies[name]
    #     if len(copy_list) >= 2:
    #         confident_line[name] = candidate_line_contigs[name]
    # store_fasta(confident_line, confident_line_path)




