import argparse
import codecs
import json
import os
import sys


cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, read_fasta_v1, multi_process_align_and_get_copies, flanking_copies


def flanking_and_getcopies(input):
    threads = 48
    reference = tmp_output_dir + '/GCF_001433935.1_IRGSP-1.0_genomic.fna'
    temp_dir = tmp_output_dir + '/temp'
    blast_program_dir = '/public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2'
    all_copies = multi_process_align_and_get_copies(input, reference, blast_program_dir,
                                                    temp_dir, 'tir', threads)

    # 在copies的两端 flanking 20bp的序列
    flanking_len = 20
    all_copies = flanking_copies(all_copies, input, reference, flanking_len, copy_num=1)
    ltr_copies = tmp_output_dir + '/ltr_copies.csv'
    # 存储all copies
    with codecs.open(ltr_copies, 'w', encoding='utf-8') as f:
        json.dump(all_copies, f)

    return all_copies


def get_logo_seq(ltr_copies):
    # (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq)
    start_logos = {}
    tail_logos = {}
    for name in ltr_copies.keys():
        copies = ltr_copies[name]
        copy = copies[0]
        copy_seq = copy[4]
        head_seq = copy_seq[0:50]
        tail_seq = copy_seq[-60:-10]
        start_logos[name] = head_seq
        tail_logos[name] = tail_seq
    return start_logos, tail_logos




if __name__ == '__main__':
    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa'
    # confident_ltr_terminal_path = tmp_output_dir + '/confident_ltr_cut.terminal.fa'
    # ltr_copies = flanking_and_getcopies(confident_ltr_terminal_path)

    confident_helitron_path = tmp_output_dir + '/confident_helitron_0.rename.cons.fa'
    helitron_copies = flanking_and_getcopies(confident_helitron_path)

    # enspm_path = tmp_output_dir + '/enspm.fa'
    # confident_TE_path = tmp_output_dir + '/confident_TE.cons.fa.final.classified'
    # contignames, contigs = read_fasta(confident_TE_path)
    # node_index = 0
    # Enspm_contigs = {}
    # for name in contignames:
    #     if name.__contains__('#DNA/CMC-EnSpm'):
    #         seq = contigs[name]
    #         Enspm_contigs[name] = seq
    # store_fasta(Enspm_contigs, enspm_path)
    #
    # Enspm_copies = flanking_and_getcopies(enspm_path)

    #取序列头部前20bp至头部后20bp
    #ltr_start_logos, ltr_tail_logos = get_logo_seq(ltr_copies)
    helitron_start_logos, helitron_tail_logos = get_logo_seq(helitron_copies)
    #Enspm_start_logos, Enspm_tail_logos = get_logo_seq(Enspm_copies)

    # ltr_terminal_start_path = tmp_output_dir + '/ltr_terminal_start.fa'
    # ltr_terminal_end_path = tmp_output_dir + '/ltr_terminal_end.fa'
    # store_fasta(ltr_start_logos, ltr_terminal_start_path)
    # store_fasta(ltr_tail_logos, ltr_terminal_end_path)

    helitron_terminal_start_path = tmp_output_dir + '/helitron_terminal_start.fa'
    helitron_terminal_end_path = tmp_output_dir + '/helitron_terminal_end.fa'
    store_fasta(helitron_start_logos, helitron_terminal_start_path)
    store_fasta(helitron_tail_logos, helitron_terminal_end_path)

    # EnSpm_terminal_start_path = tmp_output_dir + '/EnSpm_terminal_start.fa'
    # EnSpm_terminal_end_path = tmp_output_dir + '/EnSpm_terminal_end.fa'
    # store_fasta(Enspm_start_logos, EnSpm_terminal_start_path)
    # store_fasta(Enspm_tail_logos, EnSpm_terminal_end_path)
