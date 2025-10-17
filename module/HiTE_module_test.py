#!/usr/bin/env python
import argparse
import os
import re
import shutil
import sys

import codecs

import json
import tempfile
import time
import math
import uuid

#import regex
from fuzzysearch import find_near_matches
import Levenshtein
from matplotlib_venn import venn3
from numpy import median
import numpy as np
import matplotlib

from judge_TIR_transposons import is_transposons
from module.benchmarking import run_benchmarking

matplotlib.use('pdf')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from openpyxl.utils import get_column_letter
from pandas import ExcelWriter
import matplotlib.ticker as ticker
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from collections import Counter, defaultdict

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, read_fasta_v1, rename_fasta, getReverseSequence, allow_mismatch, \
    run_itrsearch, multi_process_itr, filter_large_gap_tirs, multi_process_align_and_get_copies, \
    store_copies_v1, get_TSD, store_copies, store_LTR_seq_v1, store_LTR_seq, store_LTR_seq_v2, rename_reference, \
    run_LTR_harvest, run_LTR_retriever, multi_process_align, \
    get_copies, TSDsearch_v4, overlap_with_boundary, get_copies_v1, convertToUpperCase_v1, \
    search_confident_tir, store_copies_seq, PET, multiple_alignment_blastx_v1, store2file, \
    TSDsearch_v2, remove_ltr_from_tir, \
    multi_process_tsd, get_domain_info, run_HelitronScanner, run_HelitronScanner_v1, get_longest_repeats_v3, \
    flanking_seq, multi_process_helitronscanner, split_fasta, get_longest_repeats_v4, \
    process_all_seqs, get_short_tir_contigs, multi_process_EAHelitron, search_confident_tir_v4, \
    multiple_alignment_blast_and_get_copies, split_dict_into_blocks, file_exist, flank_region_align_v5, \
    multi_process_align_v1, FMEA, judge_boundary_v8, judge_boundary_v9, run_find_members_v8, deredundant_for_LTR_v5, \
    get_candidate_non_LTR, get_candidate_non_ltr_parallel, find_nearest_polyA, find_nearest_tandem, search_polyA_TSD, \
    split_internal_out, create_or_clear_directory, copy_files, filter_short_contigs_in_genome, parse_clstr_file, \
    save_to_file, te_genome_align_minimap2, convert_paf_to_blast6_format, generate_full_length_out_v1, load_from_file, \
    multi_process_tsd_v1, judge_boundary_v5, get_genome_size, get_fixed_extend_base_threshold, \
    remove_sparse_col_in_align_file, split_and_store_sequences, get_full_length_copies_minimap2, \
    sequence2sequenceBlastn, is_TE_from_align_file
from judge_TIR_transposons import run_TIR_detection
from pan_recover_low_copy_TEs import get_pan_genome_copies, filter_true_TEs


def filter_repbase_nonTE():
    # 水稻
    repbase_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/repbase'
    TE_repbase_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/only_TE/repbase'
    if not os.path.exists(TE_repbase_dir):
        os.makedirs(TE_repbase_dir)
    repbase_names = ['athrep.ref', 'cbrrep.ref', 'drorep.ref', 'maize.ref', 'oryrep.ref', 'zebrep.ref', 'edcotrep.ref']
    nonTE_tags = ['Satellite', 'SAT', 'MSAT', 'Multicopy gene', 'rRNA', 'tRNA', 'snRNA', 'Integrated Virus', 'DNA Virus', 'Caulimoviridae']
    for repbase_name in repbase_names:
        repbase_path = repbase_dir+'/'+repbase_name
        repbase_names, repbase_contigs = read_fasta_v1(repbase_path)

        TE_contigs = {}
        for name in repbase_names:
            tag = name.split('\t')[1]
            if nonTE_tags.__contains__(tag):
                continue
            if (repbase_name == 'edcotrep.ref' and not name.__contains__('Solanum tuberosum')):
                continue
            TE_contigs[name] = repbase_contigs[name]
        store_fasta(TE_contigs, TE_repbase_dir+'/'+repbase_name)

def generate_repbases():
    # 水稻
    repbase_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/Repbase_28.06/dmel'
    repbase_path = repbase_dir + '/drorep.ref'
    repbase_names, repbase_contigs = read_fasta_v1(repbase_path)
    tmp_out_dir = repbase_dir + '/dmel'
    # tags = set()
    # for name in repbase_names:
    #     if not name.__contains__('Solanum tuberosum'):
    #         continue
    #     tag = name.split('\t')[1]
    #     tags.add(tag)
    # print(tags)
    # print(len(tags))

    ltr_tags = ['Gypsy', 'Copia', 'LTR Retrotransposon', 'BEL', 'LTR', 'Endogenous Retrovirus', 'Caulimoviridae']
    tir_tags = ['Mariner/Tc1', 'DNA transposon', 'EnSpm/CACTA', 'MuDR', 'hAT', 'Harbinger', 'Transib', 'piggyBac', 'P', 'DNA', 'Sola2', 'Kolobok', ]
    helitron_tags = ['Helitron', 'MINIME_DN']
    non_ltr_tags = ['L1', 'SINE2/tRNA', 'Non-LTR Retrotransposon', 'SINE', 'R1', 'Jockey', 'CR1', 'R2', 'RTEX', 'Hero', 'RTE', 'tRNA']

    if not os.path.exists(tmp_out_dir):
        os.makedirs(tmp_out_dir)
    ltr_repbase_path = tmp_out_dir + '/ltr.repbase.ref'
    tir_repbase_path = tmp_out_dir + '/tir.repbase.ref'
    helitron_repbase_path = tmp_out_dir + '/helitron.repbase.ref'
    non_ltr_repbase_path = tmp_out_dir + '/non_ltr.repbase.ref'
    all_repbase_path = tmp_out_dir + '/all.repbase.ref'

    all_contigs = {}
    ltr_contigs = {}
    tir_contigs = {}
    helitron_contigs = {}
    non_ltr_contigs = {}
    for name in repbase_names:
        # if not name.__contains__('Solanum tuberosum'):
        #     continue
        tag = name.split('\t')[1]
        if tag in ltr_tags:
            ltr_contigs[name] = repbase_contigs[name]
        elif tag in tir_tags:
            tir_contigs[name] = repbase_contigs[name]
        elif tag in helitron_tags:
            helitron_contigs[name] = repbase_contigs[name]
        elif tag in non_ltr_tags:
            non_ltr_contigs[name] = repbase_contigs[name]
        else:
            print(tag)
        all_contigs[name] = repbase_contigs[name]
    store_fasta(ltr_contigs, ltr_repbase_path)
    store_fasta(tir_contigs, tir_repbase_path)
    store_fasta(helitron_contigs, helitron_repbase_path)
    store_fasta(non_ltr_contigs, non_ltr_repbase_path)
    store_fasta(all_contigs, all_repbase_path)

def generate_rm2():
    # 水稻
    repbase_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/Repbase_28.06/zebrafish'
    repbase_path = repbase_dir + '/zebrafish.lib'
    tmp_out_dir = repbase_dir + '/zebrafish'
    repbase_names, repbase_contigs = read_fasta(repbase_path)
    # repbase中的文件可能有重复，我们根据序列，去掉完全一致的重复
    seq_set = set()
    for name in repbase_names:
        seq = repbase_contigs[name]
        if seq in seq_set:
            del repbase_contigs[name]
        else:
            seq_set.add(seq)
    store_fasta(repbase_contigs, repbase_path)
    repbase_names, repbase_contigs = read_fasta(repbase_path)

    tags = set()
    for name in repbase_names:
        tag = name.split('#')[1]
        tags.add(tag)
    print(tags)
    print(len(tags))

    # ltr_tags = ['LTR/Gypsy', 'LTR/Copia', 'LTR/Pao', 'LTR/Cassandra', 'LTR', 'LTR/ERVK', 'LTR/ERV1', 'LTR/Unknown', 'LTR/Caulimovirus']
    # tir_tags = ['Mariner/Tc1', 'DNA transposon', 'DNA/TcMar-Stowaway', 'DNA/TcMar', 'DNA/Maverick', 'DNA/TcMar-Pogo', 'DNA/hAT-Charlie', 'DNA/CMC-EnSpm', 'DNA/CMC', 'DNA/MULE-MuDR', 'DNA/hAT-Tag1', 'DNA/hAT-Ac', 'DNA/hAT-Tip100', 'DNA/PIF-Harbinger', 'Transib', 'piggyBac', 'DNA/P', 'DNA', 'Sola2', 'Kolobok', ]
    # helitron_tags = ['RC/Helitron', 'MINIME_DN']
    # non_ltr_tags = ['LINE/L1', 'LINE/RTE-BovB', 'Retroposon', 'Retroposon/L1-derived', 'SINE/tRNA', 'SINE/tRNA-RTE', 'SINE/ID', 'LINE/Rex-Babar', 'SINE', 'R1', 'Jockey', 'CR1', 'R2', 'RTEX', 'Hero', 'RTE']
    # unknown_tags = ['Unknown']

    if not os.path.exists(tmp_out_dir):
        os.makedirs(tmp_out_dir)
    ltr_repbase_path = tmp_out_dir + '/ltr.ref'
    tir_repbase_path = tmp_out_dir + '/tir.ref'
    helitron_repbase_path = tmp_out_dir + '/helitron.ref'
    non_ltr_repbase_path = tmp_out_dir + '/non_ltr.ref'
    unknown_repbase_path = tmp_out_dir + '/unknown.ref'

    ltr_contigs = {}
    tir_contigs = {}
    helitron_contigs = {}
    non_ltr_contigs = {}
    unknown_contigs = {}
    for name in repbase_names:
        tag = name.split('#')[1]
        if 'LTR' in tag:
            ltr_contigs[name] = repbase_contigs[name]
        elif 'DNA' in tag:
            tir_contigs[name] = repbase_contigs[name]
        elif 'Helitron' in tag:
            helitron_contigs[name] = repbase_contigs[name]
        elif 'LINE' in tag or 'SINE' in tag:
            non_ltr_contigs[name] = repbase_contigs[name]
        elif 'Unknown' in tag:
            unknown_contigs[name] = repbase_contigs[name]
        else:
            print(tag)
    store_fasta(ltr_contigs, ltr_repbase_path)
    store_fasta(tir_contigs, tir_repbase_path)
    store_fasta(helitron_contigs, helitron_repbase_path)
    store_fasta(non_ltr_contigs, non_ltr_repbase_path)
    store_fasta(unknown_contigs, unknown_repbase_path)

def generate_HiTE():
    # 水稻
    repbase_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/HiTE_lib/potato'
    repbase_path = repbase_dir + '/confident_TE.cons.fa.classified'
    repbase_names, repbase_contigs = read_fasta(repbase_path)
    tags = set()
    for name in repbase_names:
        tag = name.split('#')[1]
        tags.add(tag)
    print(tags)
    print(len(tags))

    ltr_tags = ['LTR/Gypsy', 'LTR/Copia', 'LTR/Pao', 'LTR/Cassandra', 'LTR', 'LTR/ERVK', 'LTR/ERV1', 'LTR/Unknown', 'LTR/DIRS', 'LTR/Caulimovirus']
    tir_tags = ['Mariner/Tc1', 'DNA transposon', 'DNA/TcMar-Stowaway', 'DNA/hAT-Charlie', 'DNA/TcMar-Pogo', 'DNA/CMC-EnSpm', 'DNA/MULE-MuDR', 'DNA/Crypton-H', 'DNA/hAT-Tag1', 'DNA/hAT-Ac', 'DNA/hAT-Tip100', 'DNA/hAT', 'DNA/PIF-Harbinger', 'DNA/TcMar-Tigger', 'piggyBac', 'DNA/P', 'DNA', 'Sola2', 'Kolobok', ]
    helitron_tags = ['RC/Helitron', 'MINIME_DN']
    non_ltr_tags = ['LINE/L1', 'LINE/RTE-BovB', 'Retroposon', 'Retroposon/L1-derived', 'Retroposon/RTE-derived', 'SINE/tRNA', 'SINE/tRNA-RTE', 'LINE/Rex-Babar', 'SINE', 'R1', 'Jockey', 'CR1', 'R2', 'RTEX', 'Hero', 'RTE']
    unknown_tags = ['Unknown']

    tmp_out_dir = repbase_dir + '/potato'
    if not os.path.exists(tmp_out_dir):
        os.makedirs(tmp_out_dir)
    ltr_repbase_path = tmp_out_dir + '/ltr.HiTE.ref'
    tir_repbase_path = tmp_out_dir + '/tir.HiTE.ref'
    helitron_repbase_path = tmp_out_dir + '/helitron.HiTE.ref'
    non_ltr_repbase_path = tmp_out_dir + '/non_ltr.HiTE.ref'
    unknown_repbase_path = tmp_out_dir + '/unknown.HiTE.ref'

    ltr_contigs = {}
    tir_contigs = {}
    helitron_contigs = {}
    non_ltr_contigs = {}
    unknown_contigs = {}
    for name in repbase_names:
        tag = name.split('#')[1]
        if tag in ltr_tags:
            ltr_contigs[name] = repbase_contigs[name]
        elif tag in tir_tags:
            tir_contigs[name] = repbase_contigs[name]
        elif tag in helitron_tags:
            helitron_contigs[name] = repbase_contigs[name]
        elif tag in non_ltr_tags:
            non_ltr_contigs[name] = repbase_contigs[name]
        elif tag in unknown_tags:
            unknown_contigs[name] = repbase_contigs[name]
    store_fasta(ltr_contigs, ltr_repbase_path)
    store_fasta(tir_contigs, tir_repbase_path)
    store_fasta(helitron_contigs, helitron_repbase_path)
    store_fasta(non_ltr_contigs, non_ltr_repbase_path)
    store_fasta(unknown_contigs, unknown_repbase_path)

def generate_EDTA():
    # 水稻
    repbase_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/EDTA_lib/potato_reference'
    repbase_path = repbase_dir + '/GCF_000226075.1_SolTub_3.0_genomic.fna.mod.EDTA.TElib.fa'
    repbase_names, repbase_contigs = read_fasta(repbase_path)
    tags = set()
    for name in repbase_names:
        tag = name.split('#')[1]
        tags.add(tag)
    print(tags)
    print(len(tags))

    ltr_tags = ['LTR/Gypsy', 'LTR/Copia', 'LTR/Pao', 'LTR/Cassandra', 'LTR', 'LTR/ERVK', 'LTR/ERV1', 'LTR/Unknown', 'LTR/unknown', 'LTR/DIRS']
    tir_tags = ['MITE/DTA', 'MITE/DTC', 'MITE/DTM', 'MITE/DTH', 'MITE/DTT', 'TIR/PIF_Harbinger', 'DNA/DTH', 'DNA/DTC', 'DNA/DTM', 'DNA/DTT', 'DNA/DTA', 'TIR/MuDR_Mutator', 'Mariner/Tc1', 'DNA transposon', 'DNA/TcMar-Stowaway', 'DNA/hAT-Charlie', 'DNA/CMC-EnSpm', 'DNA/MULE-MuDR', 'DNA/hAT-Tag1', 'DNA/hAT-Ac', 'DNA/hAT-Tip100', 'DNA/hAT', 'DNA/PIF-Harbinger', 'DNA/TcMar-Tigger', 'piggyBac', 'DNA/P', 'DNA', 'Sola2', 'Kolobok', ]
    helitron_tags = ['RC/Helitron', 'MINIME_DN', 'DNA/Helitron']
    non_ltr_tags = ['LINE/L1', 'LINE/RTE-BovB', 'LINE/unknown', 'Retroposon', 'Retroposon/L1-derived', 'SINE/tRNA', 'LINE/Rex-Babar', 'SINE', 'R1', 'Jockey', 'CR1', 'R2', 'RTEX', 'Hero', 'RTE']
    unknown_tags = ['Unknown']
    
    tmp_out_dir = repbase_dir + '/potato'
    if not os.path.exists(tmp_out_dir):
        os.makedirs(tmp_out_dir)
    ltr_repbase_path = tmp_out_dir + '/ltr.EDTA.ref'
    tir_repbase_path = tmp_out_dir + '/tir.EDTA.ref'
    helitron_repbase_path = tmp_out_dir + '/helitron.EDTA.ref'
    non_ltr_repbase_path = tmp_out_dir + '/non_ltr.EDTA.ref'
    unknown_repbase_path = tmp_out_dir + '/unknown.EDTA.ref'

    ltr_contigs = {}
    tir_contigs = {}
    helitron_contigs = {}
    non_ltr_contigs = {}
    unknown_contigs = {}
    for name in repbase_names:
        raw_name = name.split('#')[0]
        tag = name.split('#')[1]
        if tag in ltr_tags:
            ltr_contigs[name] = repbase_contigs[name]
        elif tag in tir_tags:
            tir_contigs[raw_name] = repbase_contigs[name]
        elif tag in helitron_tags:
            helitron_contigs[raw_name] = repbase_contigs[name]
        elif tag in non_ltr_tags:
            non_ltr_contigs[name] = repbase_contigs[name]
        elif tag in unknown_tags:
            unknown_contigs[name] = repbase_contigs[name]
    store_fasta(ltr_contigs, ltr_repbase_path)
    store_fasta(tir_contigs, tir_repbase_path)
    store_fasta(helitron_contigs, helitron_repbase_path)
    store_fasta(non_ltr_contigs, non_ltr_repbase_path)
    store_fasta(unknown_contigs, unknown_repbase_path)

def generate_zebrafish_repbases():
    # 水稻
    repbase_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/repbase'
    repbase_path = repbase_dir + '/zebrep.ref.final.classified'
    repbase_names, repbase_contigs = read_fasta(repbase_path)

    tags = set()
    for name in repbase_names:
        tag = name.split('#')[1]
        tags.add(tag)
    print(tags)
    print(len(tags))

    tmp_out_dir = repbase_dir + '/drerio'
    ltr_repbase_path = tmp_out_dir + '/ltr.repbase.ref'
    tir_repbase_path = tmp_out_dir + '/tir.repbase.ref'
    helitron_repbase_path = tmp_out_dir + '/helitron.repbase.ref'
    non_ltr_repbase_path = tmp_out_dir + '/non_ltr.repbase.ref'

    ltr_contigs = {}
    tir_contigs = {}
    helitron_contigs = {}
    non_ltr_contigs = {}
    for name in repbase_names:
        tag = name.split('#')[1]
        if tag.__contains__('LTR'):
            ltr_contigs[name] = repbase_contigs[name]
        elif tag.__contains__('DNA') and not tag.__contains__('Crypton') and not tag.__contains__('Maverick'):
            tir_contigs[name] = repbase_contigs[name]
        elif tag.__contains__('Helitron'):
            helitron_contigs[name] = repbase_contigs[name]
        elif tag.__contains__('LINE') or tag.__contains__('SINE'):
            non_ltr_contigs[name] = repbase_contigs[name]
    store_fasta(ltr_contigs, ltr_repbase_path)
    store_fasta(tir_contigs, tir_repbase_path)
    store_fasta(helitron_contigs, helitron_repbase_path)
    store_fasta(non_ltr_contigs, non_ltr_repbase_path)

def rename_TE_file():
    tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/dmel'
    candidate_tir_path = tmp_dir + '/confident_tir.fa'
    candidate_tir_rename_path = tmp_dir + '/confident_tir.rename.fa'
    rename_fasta(candidate_tir_path, candidate_tir_rename_path)

    candidate_tir_rename_consensus = tmp_dir + '/confident_tir.rename.cons.fa'
    tools_dir = os.getcwd() + '/../tools'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + candidate_tir_rename_path + ' -o ' + candidate_tir_rename_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    candidate_helitron_path = tmp_dir + '/confident_helitron_0.fa'
    candidate_helitron_rename_path = tmp_dir + '/confident_helitron_0.rename.fa'
    rename_fasta(candidate_helitron_path, candidate_helitron_rename_path)

    candidate_helitron_rename_consensus = tmp_dir + '/confident_helitron_0.rename.cons.fa'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + candidate_helitron_rename_path + ' -o ' + candidate_helitron_rename_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    candidate_other_path = tmp_dir + '/confident_other0.fa'
    candidate_other_rename_path = tmp_dir + '/confident_other0.rename.fa'
    rename_fasta(candidate_other_path, candidate_other_rename_path)

    candidate_other_rename_consensus = tmp_dir + '/confident_other0.rename.cons.fa'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + candidate_other_rename_path + ' -o ' + candidate_other_rename_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

def rename_EDTA_file():
    tmp_dir = '/public/home/hpc194701009/repeat_detect_tools/EDTA-master/genome_test/oryza_sativa'
    raw_EDTA_tir = tmp_dir + '/EDTA_TIR/GCF_001433935.1_IRGSP-1.0_genomic.fna.mod.TIR.raw.fa'
    EDTA_tir = tmp_dir + '/EDTA_TIR/EDTA_tir.fa'
    rename_fasta(raw_EDTA_tir, EDTA_tir)

    raw_EDTA_helitron = tmp_dir + '/EDTA_Helitron/GCF_001433935.1_IRGSP-1.0_genomic.fna.mod.Helitron.raw.fa'
    EDTA_helitron = tmp_dir + '/EDTA_Helitron/EDTA_helitron.fa'
    rename_fasta(raw_EDTA_helitron, EDTA_helitron)

    raw_EDTA_ltr = tmp_dir + '/EDTA_LTR/GCF_001433935.1_IRGSP-1.0_genomic.fna.mod.LTR.raw.fa'
    EDTA_ltr = tmp_dir + '/EDTA_LTR/EDTA_ltr.fa'
    rename_fasta(raw_EDTA_ltr, EDTA_ltr)


def get_repbase_copies(tir_repbase_path, copy_info_path):
    threads = 48

    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/dmel'
    reference = tmp_output_dir + '/dmel-all-chromosome-r5.43.fasta'
    temp_dir = tmp_output_dir + '/temp'
    blast_program_dir = '/public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2'
    all_copies = multi_process_align_and_get_copies(tir_repbase_path, reference, blast_program_dir,
                                                    temp_dir, 'tir', threads, query_coverage=0.99)

    # 在copies的两端 flanking 20bp的序列
    flanking_len = 20
    all_copies = flanking_copies(all_copies, tir_repbase_path, reference, flanking_len, copy_num=10)

    # 判断copy的TSD信息
    # tsd_info = {query_name: {copy1: tsd+','+seq}, {copy2: tsd+','+seq}, {total_copy_num:}, {tsd_copy_num:}}
    tsd_info = get_TSD(all_copies, flanking_len)

    copies_list = []
    for query_name in tsd_info.keys():
        info = tsd_info[query_name]
        total_copy_num = info['total_copy_num']
        total_copy_len = info['total_copy_len']
        copies_list.append((query_name, total_copy_num, total_copy_len))
    copies_list.sort(key=lambda x: -x[2])

    store_copies(tsd_info, copy_info_path)

    return copies_list

def summary_TGCA_motif():
    tmp_out_dir = repbase_dir + '/rice'
    tir_repbase_path = tmp_out_dir + '/tir.repbase.ref'

    #统计tir repbase中有多少序列是有TGCA motif
    tir_names, tir_contigs = read_fasta(tir_repbase_path)
    count = 0
    for name in tir_names:
        seq = tir_contigs[name]
        if seq[0:2] == 'TG' and seq[-2:] == 'CA':
            count += 1
            print(name)
    print(count)

def summary_tir_5bp(tir_repbase_path):
    tir_names, tir_contigs = read_fasta(tir_repbase_path)
    diversity_set = set()
    allow_mismatch_num = 2
    count = 0
    for name in tir_names:
        seq = tir_contigs[name]
        first_5bp = seq[0:5]
        end_5bp = getReverseSequence(seq[-5:])
        if allow_mismatch(first_5bp, end_5bp, allow_mismatch_num):
            count += 1
        else:
            diversity_set.add(name)
    print(count)
    print(len(tir_names))
    print(diversity_set)

def summary_tir_len(tir_repbase_path):
    # 统计一下tir_repbase的tir有多长
    TRsearch_dir = '/public/home/hpc194701009/HiTE/ReferenceMode/tools/'
    run_itrsearch(TRsearch_dir, tir_repbase_path, tmp_out_dir)
    tir_repbase_path = tmp_out_dir + '/tir.repbase.ref.itr'
    tir_names, tir_contigs = read_fasta_v1(tir_repbase_path)
    count = 0
    for name in tir_names:
        parts = name.split('ITR')
        # ITR(1,61)..(166,113)
        ITR_info = parts[1].split(' ')[0].replace('(', '').replace(')', '')

        ITR_info_parts = ITR_info.split('..')

        ITR_left_pos_parts = ITR_info_parts[0].split(',')
        ITR_right_pos_parts = ITR_info_parts[1].split(',')
        lITR_start = int(ITR_left_pos_parts[0])
        lITR_end = int(ITR_left_pos_parts[1])
        lITR_len = lITR_end - lITR_start + 1
        rITR_start = int(ITR_right_pos_parts[0])
        rITR_end = int(ITR_right_pos_parts[1])
        rITR_len = rITR_start - rITR_end + 1
        print(lITR_len, rITR_len)
        if abs(rITR_len - lITR_len) > 2:
            count += 1
    print(count, len(tir_names))

    min_len = 10000000
    max_len = 0
    tir_len_list = []
    less_15bp = []
    for name in tir_names:
        tir_len = int(name.split('Length itr=')[1]) + 1
        tir_len_list.append(tir_len)
        if tir_len < min_len:
            min_len = tir_len
        if tir_len > max_len:
            max_len = tir_len
        if tir_len < 11:
            less_15bp.append(tir_len)
    avg_tir_len = sum(tir_len_list)/len(tir_len_list)
    print(avg_tir_len)
    print(min_len)
    print(max_len)
    print(len(less_15bp))
    print(len(tir_names))

def test_filter_large_gap_tirs():
    # 测试把confident_tir.fa的差距过大的terminal过滤掉是什么结果
    TRsearch_dir = '/public/home/hpc194701009/HiTE/ReferenceMode/tools/'
    tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa'
    confident_tir_path = tmp_dir + '/confident_tir.fa'
    confident_tir_out = tmp_dir + '/confident_tir.fa.itr'
    temp_dir = tmp_dir + '/tir_test'
    multi_process_itr(confident_tir_path, confident_tir_out, temp_dir, TRsearch_dir, threads=48)
    filter_large_gap_tirs(confident_tir_out, confident_tir_out)

def test_filter_diff_5bp_tirs():
    # 测试过滤掉confident_tir.fa中开始5bp，结束5bp不一致的序列
    TRsearch_dir = '/public/home/hpc194701009/HiTE/ReferenceMode/tools/'
    tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa'
    confident_tir_path = tmp_dir + '/confident_tir.fa'
    confident_tir_filter_path = tmp_dir + '/confident_tir.filter5bp.fa'
    tir_names, tir_contigs = read_fasta(confident_tir_path)
    for name in tir_names:
        seq = tir_contigs[name]
        first_5bp = seq[0:5]
        end_5bp = seq[-5:]
        if first_5bp != getReverseSequence(end_5bp):
            del tir_contigs[name]
    store_fasta(tir_contigs, confident_tir_filter_path)

def summary_not_perfect_repbase():
    #我们观察，有哪些序列是没有出现在比对中
    align_file = '/homeb/hukang/KmerRepFinder_test/library/get_family_summary_test/perfect.families'
    query_names = set()
    with open(align_file, 'r') as f_r:
        for line in f_r:
            query_name = line.replace('\n', '')
            query_names.add(query_name)

    tmp_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/repbase/ath'
    names, contigs = read_fasta(tmp_dir + '/tir.repbase.ref')
    names = set(names)
    diff_set = names.difference(query_names)
    print(diff_set)
    print('not appear size: ' + str(len(diff_set)))

    lost_tirs = '/homeb/hukang/KmerRepFinder_test/library/tir_test/lost_tir.fa'
    lost_contigs = {}
    for name in diff_set:
        lost_contigs[name] = contigs[name]
    store_fasta(lost_contigs, lost_tirs)

def summary_not_appear_repbase():
    #我们观察，有哪些序列是没有出现在比对中
    align_file = '/homeb/hukang/KmerRepFinder_test/library/get_family_summary_test/file_final.0.1.txt'
    query_names = set()
    with open(align_file, 'r') as f_r:
        for line in f_r:
            query_name = line.split('\t')[0]
            query_names.add(query_name)

    tmp_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/repbase/ath'
    names, contigs = read_fasta(tmp_dir + '/tir.repbase.ref')
    names = set(names)
    diff_set = names.difference(query_names)
    print(diff_set)
    print('not appear size: ' + str(len(diff_set)))

    lost_tirs = '/homeb/hukang/KmerRepFinder_test/library/tir_test/lost_tir.fa'
    lost_contigs = {}
    for name in diff_set:
        lost_contigs[name] = contigs[name]
    store_fasta(lost_contigs, lost_tirs)

    # diff_set1 = set()
    # #去除掉TGCA后还有多少条
    # for name in diff_set:
    #     seq = contigs[name]
    #     if seq[0:2] == 'TG' and seq[-2:] == 'CA':
    #         continue
    #     else:
    #         diff_set1.add(name)
    # print(diff_set1)
    # print('not appear size: ' + str(len(diff_set1)))

def not_found_repbase(tir_repbase_path, copies_list):
    all_names, contigs = read_fasta(tir_repbase_path)
    present_file = '/public/home/hpc194701009/KmerRepFinder_test/library/get_family_summary_test/present.all.families'
    present_repbase_names = set()
    with open(present_file, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            present_repbase_names.add(line)
    all_names = set(all_names)
    diff_names = all_names.difference(present_repbase_names)
    print(diff_names)
    print(len(diff_names))
    #输出未找到的repbase，按照总长度排序
    not_found_copies = []
    for copy in copies_list:
        name = copy[0]
        if name in diff_names:
            not_found_copies.append(copy)
    print(not_found_copies)


def test_LTR_finder():
    ref_index = 0
    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/dmel'
    tools_dir = os.getcwd() + '/../tools'
    candidate_ltr_finder_cut_path = tmp_output_dir + '/candidate_ltr_finder_cut_' + str(ref_index) + '.fa'
    candidate_ltr_finder_cut_rename_path = tmp_output_dir + '/candidate_ltr_finder_cut_' + str(ref_index) + '.rename.fa'
    rename_fasta(candidate_ltr_finder_cut_path, candidate_ltr_finder_cut_rename_path)

    candidate_ltr_finder_cut_consensus = tmp_output_dir + '/candidate_ltr_finder_cut_' + str(ref_index) + '.rename.cons.fa'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + candidate_ltr_finder_cut_rename_path + ' -o ' + candidate_ltr_finder_cut_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)


def test_LTR_harvest():
    ref_index = 0
    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/dmel'
    tools_dir = os.getcwd() + '/../tools'
    candidate_ltr_harvest_cut_path = tmp_output_dir + '/candidate_ltr_harvest_cut_'+str(ref_index)+'.fa'
    candidate_ltr_harvest_cut_rename = tmp_output_dir + '/candidate_ltr_harvest_cut_'+str(ref_index)+'.rename.fa'
    rename_fasta(candidate_ltr_harvest_cut_path, candidate_ltr_harvest_cut_rename)

    candidate_ltr_harvest_cut_consensus = tmp_output_dir + '/candidate_ltr_harvest_cut_'+str(ref_index)+'.rename.cons.fa'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + candidate_ltr_harvest_cut_rename + ' -o ' + candidate_ltr_harvest_cut_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)


def get_LTR_seqs():
    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/dmel'
    ltr_harvest_output = tmp_output_dir + '/genome_0.fa.harvest.scn'
    ltr_finder_output = tmp_output_dir + '/dmel-all-chromosome-r5.43.fasta.cut0.rename.fa.finder.combine.scn'
    ref_rename_path = tmp_output_dir + '/dmel-all-chromosome-r5.43.fasta.cut0.rename.fa'
    candidate_ltr_path = tmp_output_dir + '/candidate_ltr_finder_0.fa'
    candidate_ltr_cut_path = tmp_output_dir + '/candidate_ltr_finder_cut_0.fa'
    store_LTR_seq_v1(ltr_finder_output, ref_rename_path, candidate_ltr_path, candidate_ltr_cut_path)

    cut_reference = tmp_output_dir + '/dmel-all-chromosome-r5.43.fasta.cut0.fa'
    candidate_ltr_path = tmp_output_dir + '/candidate_ltr_harvest_0.fa'
    candidate_ltr_cut_path = tmp_output_dir + '/candidate_ltr_harvest_cut_0.fa'
    store_LTR_seq_v2(ltr_harvest_output, cut_reference, candidate_ltr_path, candidate_ltr_cut_path)

def identify_new_TIR(tmp_output_dir):
    # 用cd-hit-est，使用-aS 0.8 –aL 0.8 –c 0.8进行聚类，然后分析类中没有curated library出现的转座子为新的TIR。
    confident_tir_path = tmp_output_dir + '/confident_tir.fa'
    repbase_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/repbase/rice'
    tir_repbase_path = repbase_dir + '/tir.repbase.ref'

    total_tir_path = tmp_output_dir + '/total_tir.fa'
    total_tir_consensus = tmp_output_dir + '/total_tir.cons.fa'
    os.system('cat '+tir_repbase_path+' '+confident_tir_path+' > '+total_tir_path)

    confident_tir_names, confident_tir_contigs = read_fasta(confident_tir_path)

    tools_dir = os.getcwd() + '/../tools'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -aL ' + str(0.8) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + total_tir_path + ' -o ' + total_tir_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    cluster_file = total_tir_consensus + '.clstr'
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
                name = name[0: len(name) - 3]
                cur_cluster.append(name)
                if line.endswith('*'):
                    clusters['rep_' + str(cluster_idx)] = name
    #print(clusters)
    new_tir_contigs = {}
    exist_tir_names = set()
    for key in clusters.keys():
        if not key.startswith('rep_'):
            cluster = clusters[key]
            has_repbase = False
            for name in cluster:
                if not name.startswith('TIR_'):
                    has_repbase = True
                    break
            #没有repbase序列，代表这是我们新识别的TIR
            if not has_repbase:
                represent_name = clusters['rep_'+str(key)]
                new_tir_contigs[represent_name] = confident_tir_contigs[represent_name]
            else:
                exist_tir_names.add(clusters['rep_'+str(key)])
    print('novel TIR:')
    #print(new_tir_contigs.keys())
    print(len(new_tir_contigs))

    novel_tir_consensus = tmp_output_dir + '/novel_tir.fa'
    store_fasta(new_tir_contigs, novel_tir_consensus)

    # print(exist_tir_names)
    # print(len(exist_tir_names))

    # 1.用itrseach识别repbase_TIR以及新TIR中的TIRs
    TRsearch_dir = '/home/hukang/HiTE/tools'
    run_itrsearch(TRsearch_dir, tir_repbase_path, repbase_dir)
    tir_repbase_out = tir_repbase_path + '.itr'
    repbase_itr_names, repbase_itr_contigs = read_fasta_v1(tir_repbase_out)

    run_itrsearch(TRsearch_dir, novel_tir_consensus, tmp_output_dir)
    novel_tir_out = novel_tir_consensus + '.itr'
    novel_itr_names, novel_itr_contigs = read_fasta_v1(novel_tir_out)

    repbase_tir_contigs = {}
    for name in repbase_itr_names:
        orig_name = name.split('\t')[0]
        seq = repbase_itr_contigs[name]
        LTR_info_parts = name.split('ITR')[1].split(' ')[0].replace('(', '').replace(')', '').split('..')
        LTR_left_pos_parts = LTR_info_parts[0].split(',')
        LTR_right_pos_parts = LTR_info_parts[1].split(',')
        lLTR_start = int(LTR_left_pos_parts[0])
        lLTR_end = int(LTR_left_pos_parts[1])
        rLTR_start = int(LTR_right_pos_parts[1])
        rLTR_end = int(LTR_right_pos_parts[0])

        left_LTR = seq[lLTR_start - 1: lLTR_end]
        right_LTR = seq[rLTR_start - 1: rLTR_end]
        repbase_tir_contigs[orig_name+'-lTIR'] = left_LTR
        repbase_tir_contigs[orig_name + '-rTIR'] = right_LTR

    novel_tir_contigs = {}
    for name in novel_itr_names:
        orig_name = name.split(' ')[0]
        seq = novel_itr_contigs[name]
        LTR_info_parts = name.split('ITR')[1].split(' ')[0].replace('(', '').replace(')', '').split('..')
        LTR_left_pos_parts = LTR_info_parts[0].split(',')
        LTR_right_pos_parts = LTR_info_parts[1].split(',')
        lLTR_start = int(LTR_left_pos_parts[0])
        lLTR_end = int(LTR_left_pos_parts[1])
        rLTR_start = int(LTR_right_pos_parts[1])
        rLTR_end = int(LTR_right_pos_parts[0])


        left_LTR = seq[lLTR_start - 1: lLTR_end]
        right_LTR = seq[rLTR_start - 1: rLTR_end]
        novel_tir_contigs[orig_name+'-lTIR'] = left_LTR
        novel_tir_contigs[orig_name + '-rTIR'] = right_LTR
    #print(len(repbase_tir_contigs), len(novel_tir_contigs))

    repbase_tirs = tmp_output_dir + '/repbase.itr.fa'
    novel_tirs = tmp_output_dir + '/novel.itr.fa'
    store_fasta(repbase_tir_contigs, repbase_tirs)
    store_fasta(novel_tir_contigs, novel_tirs)
    total_tirs = tmp_output_dir + '/total.itr.fa'
    total_tirs_consensus = tmp_output_dir + '/total.itr.cons.fa'
    os.system('cat '+repbase_tirs+' '+novel_tirs+' > '+total_tirs)
    # 2.用 cd-hit-est -aS 0.8 -aL 0.8 -c 0.8 对新的TIRs和repbase_TIR的TIRs进行聚类。
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + total_tirs + ' -o ' + total_tirs_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    # 3.分析类中没有repbase_TIR的序列，为新的TIRs；有repbase_TIR的序列，为已知的TIRs。
    ##对于具有已知TIRs的，novel TIR，很可能是non-autonomous TIR。
    cluster_file = total_tirs_consensus + '.clstr'
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
                name = name[0: len(name) - 3]
                cur_cluster.append(name)
                if line.endswith('*'):
                    clusters['rep_' + str(cluster_idx)] = name
    #print(clusters)
    new_tir_names = set()

    for key in clusters.keys():
        if not key.startswith('rep_'):
            cluster = clusters[key]
            has_repbase = False
            for name in cluster:
                if not name.startswith('TIR_'):
                    has_repbase = True
                    break
            #没有repbase序列，代表这是我们新识别的TIR
            if not has_repbase:
                for name in cluster:
                    new_tir_names.add(name.split('-')[0])
    print('novel TIR with new TIR terminals:')
    #print(new_tir_names)
    print(len(new_tir_names))
    return new_tir_names

    # #统计novel TIR的拷贝数量
    # threads = 40
    # reference = tmp_output_dir + '/GCF_001433935.1_IRGSP-1.0_genomic.fna'
    # temp_dir = tmp_output_dir + '/temp'
    # blast_program_dir = '/home/hukang/repeat_detect_tools/rmblast-2.9.0-p2'
    # all_copies = multi_process_align_and_get_copies(novel_tir_consensus, reference, blast_program_dir,
    #                                                 temp_dir, 'tir', threads)
    #
    # # 在copies的两端 flanking 20bp的序列
    # flanking_len = 20
    # all_copies, tsd_info = flanking_copies(all_copies, novel_tir_consensus, reference, flanking_len, copy_num=-1)
    #
    # multicopy_novel_tirs_num = 0
    # # 统计每个query_name拷贝的数量
    # # query_name -> (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq)
    # query_copy_num = {}
    # for query_name in all_copies.keys():
    #     copies = all_copies[query_name]
    #     query_copy_num[query_name] = len(copies)
    #     if len(copies) >= 2:
    #         multicopy_novel_tirs_num += 1
    # #print(query_copy_num)
    # print('novel TIR with copy number >= 2:')
    # print(multicopy_novel_tirs_num)
    #
    # query_copy_num_path = tmp_output_dir + '/novel_tir_copies_num.csv'
    # # 存储query_copy_num
    # with codecs.open(query_copy_num_path, 'w', encoding='utf-8') as f:
    #     json.dump(query_copy_num, f)
    #
    # novel_tir_copies = tmp_output_dir + '/novel_tir_copies.csv'
    # # 存储all copies
    # with codecs.open(novel_tir_copies, 'w', encoding='utf-8') as f:
    #     json.dump(all_copies, f)
    #
    # # 判断copy的TSD信息
    # # tsd_info = {query_name: {copy1: tsd+','+seq}, {copy2: tsd+','+seq}, {total_copy_num:}, {tsd_copy_num:}}
    # tsd_info = get_TSD(all_copies, flanking_len)
    #
    # copies_list = []
    # for query_name in tsd_info.keys():
    #     info = tsd_info[query_name]
    #     total_copy_num = info['total_copy_num']
    #     total_copy_len = info['total_copy_len']
    #     copies_list.append((query_name, total_copy_num, total_copy_len))
    # copies_list.sort(key=lambda x: -x[2])
    #
    # copy_info_path = tmp_output_dir + '/novel_tir.copies.info'
    # store_copies(tsd_info, copy_info_path)
    #
    # #计算Venn图指标
    # #1.novel tir copy > 3 and novel tirs with new terminals
    # #2.novel tir copy > 3 and novel tirs with known terminals
    # #3.novel tir copy <=3 and novel tirs with new terminals
    # #4.novel tir copy <=3 and novel tirs with known terminals
    # query_copy_over_3 = set()
    # query_copy_less_3 = set()
    # for query_name in query_copy_num.keys():
    #     copy_num = query_copy_num[query_name]
    #     if copy_num >= 2:
    #         query_copy_over_3.add(query_name)
    #     else:
    #         query_copy_less_3.add(query_name)
    #
    # c1 = 0
    # c2 = 0
    # c3 = 0
    # c4 = 0
    # for query_name in new_tir_contigs.keys():
    #     if query_name in query_copy_over_3 and query_name in new_tir_names:
    #         c1 += 1
    #     elif query_name in query_copy_over_3 and query_name not in new_tir_names:
    #         c2 += 1
    #     elif query_name in query_copy_less_3 and query_name in new_tir_names:
    #         c3 += 1
    #     elif query_name in query_copy_less_3 and query_name not in new_tir_names:
    #         c4 += 1
    #
    # print(c1, c2, c3, c4)
    # print(new_tir_names)
    # print(new_tir_contigs.keys())
    # print(query_copy_over_3)



def draw_dist(input_file):
    #tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa'
    #tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library/RepeatMasking_test/rice_no_kmer'
    #query_copy_num_path = tmp_output_dir + '/novel_tir_copies_num.csv'
    # query_copy_num_path = input_file
    # file = open(query_copy_num_path, 'r')
    # js = file.read()
    # query_copy_num = json.loads(js)
    # print(query_copy_num)
    # count = 0
    # for name in query_copy_num.keys():
    #     num = query_copy_num[name]
    #     if num == 1:
    #         count += 1
    # print(count)

    #画小提琴图
    # output_path = tmp_output_dir + '/novel_tir_copy_dist.txt'
    # with open(output_path, 'w') as f_save:
    #     f_save.write('copy number\tmethods\n')
    #     label = 'novel TIR elements'
    #     for name in query_copy_num.keys():
    #         num = query_copy_num[name]
    #         f_save.write(str(num) + '\t' + str(label) + '\n')
    #
    # df = pd.read_csv(output_path, sep='\t', encoding='utf-8')
    # print(df)
    #sns.boxplot(x=df["methods"], y=df["copy number"])
    #sns.violinplot(x=df["methods"], y=df["copy number"])
    # my_pal = {"HiTE-Helitron": "#4497B1", "HiTE-Helitron-NoFiltering": "#F7B92E"}
    # sns.violinplot(x=df["methods"], y=df["length"], palette=my_pal)
   # plt.show()
    # # Calculate number of obs per group & median to position labels
    # medians = df.groupby(['methods'])['copy number'].median().values
    # print(medians)

    query_copy_num = []
    with open(input_file, 'r') as f_r:
        for i,line in enumerate(f_r):
            line = line.replace('\n', '')
            if i == 0 or line.strip()== '""':
                continue
            query_copy_num.append(float(line))

    y = list(query_copy_num)
    print(len(y))
    x = pd.Series(y, name="copy number")
    sns.set_theme(style="ticks", font='Times New Roman', font_scale=1.4)
    sns.set_context("paper")
    ax = sns.distplot(x, kde=True, color='green')
    #ax.set_xlim(5, )
    plt.show()
    # plt.savefig(tmp_output_dir + "/copy_num.eps", format='eps', dpi=1000)
    # plt.savefig(tmp_output_dir + "/copy_num.svg", format='svg', dpi=150, bbox_inches='tight')


def test_EAHelitron():
    # EAHelitron_command = 'cd ' + temp_dir + ' && ' + 'perl ' + EAHelitron + '/EAHelitron -o ' + str(
    #     partition_index) + ' -u 20000 -T "ATC" -r 3 ' + all_candidate_helitron_path
    # os.system(EAHelitron_command + '> /dev/null 2>&1')

    EAHelitron_output = '/home/hukang/EDTA/krf_test/rice/EAHelitron_output/EAHelitron.5.fa'
    EAHelitron_output_clean = '/home/hukang/EDTA/krf_test/rice/EAHelitron_output/EAHelitron.5.clean.fa'
    # names, contigs = read_fasta_v1(EAHelitron_output)
    # new_contigs = {}
    # for name in names:
    #     raw_name = name.split(' ')[1]
    #     parts = raw_name.split(':')
    #     raw_name = parts[0]
    #     seq = contigs[name]
    #     if not new_contigs.__contains__(raw_name):
    #         new_contigs[raw_name] = seq
    #     else:
    #         old_seq = new_contigs[raw_name]
    #         if len(seq) < len(old_seq):
    #             new_contigs[raw_name] = seq
    # store_fasta(new_contigs, EAHelitron_output_clean)

    #统计其平均长度
    EAHelitron_output_clean = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa/EAHelitron_output/EAHelitron.5.clean.fa'
    names, contigs = read_fasta_v1(EAHelitron_output_clean)
    print(len(contigs))
    total_len = 0
    for name in contigs.keys():
        seq = contigs[name]
        total_len += len(seq)
    avg_len = float(total_len) / len(contigs)
    print(avg_len)


def reduce_library_size():
    # 1.去掉nested TE, 让所有模型都保持 full-length
    threads = 48
    blast_program_dir = '/public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2'
    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/drerio'
    tir_path = tmp_output_dir + '/confident_tir.rename.cons.fa'
    clean_tir_path = tmp_output_dir + '/confident_tir.clean.fa'
    remove_nested_command = 'python remove_nested_lib.py ' \
                            + ' -t ' + str(threads) + ' --blast_program_dir ' + blast_program_dir \
                            + ' --tmp_output_dir ' + tmp_output_dir + ' --max_iter_num ' + str(3) \
                            + ' --input1 ' + tir_path \
                            + ' --input2 ' + tir_path \
                            + ' --output ' + clean_tir_path
    #os.system(remove_nested_command)

    # 2.cd-hit -aS 0.95 -c 0.8看能否
    tools_dir = os.getcwd() + '/../tools'
    tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/drerio'
    clean_tir_consensus = tmp_output_dir + '/confident_tir.clean.cons.fa'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + clean_tir_path + ' -o ' + clean_tir_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)


def run_LTR_test():
    log = Logger('HiTE.log', level='debug')
    param_config_path = os.getcwd() + "/../ParamConfig.json"
    # read param config
    with open(param_config_path, 'r') as load_f:
        param = json.load(load_f)
    load_f.close()
    Genome_Tools_Home = param['Genome_Tools_Home']
    LTR_finder_parallel_Home = param['LTR_finder_parallel_Home']
    LTR_retriever_Home = param['LTR_retriever_Home']

    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/maize/LTR_test'
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)
    ref_rename_path = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Zea_mays/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.rename.fna'
    ref_index = 'all'
    threads = 48

    # run_LTR_harvest(Genome_Tools_Home, ref_rename_path, tmp_output_dir, ref_index, log)
    # LTR_finder_parallel_command = 'perl ' + LTR_finder_parallel_Home + '/LTR_FINDER_parallel -harvest_out -seq ' + ref_rename_path + ' -threads ' + str(
    #     threads)
    # os.system('cd ' + tmp_output_dir + ' && ' + LTR_finder_parallel_command + ' > /dev/null 2>&1')
    ltrharvest_output = tmp_output_dir + '/genome_' + str(ref_index) + '.fa.harvest.scn'
    ltrfinder_output = tmp_output_dir + '/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.rename.fna.finder.combine.scn'
    ltr_output = tmp_output_dir + '/genome_' + str(ref_index) + '.fa.rawLTR.scn'
    os.system('cat ' + ltrharvest_output + ' ' + ltrfinder_output + ' > ' + ltr_output)
    run_LTR_retriever(LTR_retriever_Home, ref_rename_path, tmp_output_dir, ref_index, threads, log)


def test_no_RepeatMasking_time():
    log = Logger('HiTE.log', level='debug')

    tmp_dir = '/homeb/hukang/KmerRepFinder_test/library/all_tools_run_lib'
    species_name = ['cb', 'dmel', 'rice', 'drerio']
    ref_indexs_dict = {'cb': [0], 'dmel': [0], 'rice': [0], 'drerio': [0, 1, 2, 3, 4, 5, 6, 7]}

    references = ['GCF_000004555.2_CB4_genomic.rename.fna', 'dmel-all-chromosome-r5.43.rename.fasta',
                  'GCF_001433935.1_IRGSP-1.0_genomic.rename.fna', 'GCF_000002035.6_GRCz11_genomic.rename.fnaa']

    blast_program_dir = '/home/hukang/repeat_detect_tools/rmblast-2.9.0-p2'
    fixed_extend_base_threshold = 1000
    max_repeat_len = 30000
    threads = 40
    debug = 1
    for index, species in enumerate(species_name):
        tmp_output_dir = tmp_dir + '/' + species + '/HiTE'
        ref_indexs = ref_indexs_dict[species]
        ref_name = references[index]
        for ref_index in ref_indexs:
            starttime = time.time()
            log.logger.info('Species: ' + species)
            log.logger.info('Start 2.2: Coarse-grained boundary mapping')
            log.logger.info('------generate longest_repeats.fa')
            repeats_path = tmp_output_dir + '/' + ref_name + '.cut' + str(ref_index) + '.fa'
            longest_repeats_path = tmp_output_dir + '/genome_longest_repeats_' + str(ref_index) + '.fa'
            log.logger.debug(repeats_path)
            log.logger.debug(longest_repeats_path)
            # -------------------------------Stage02: this stage is used to do pairwise comparision, determine the repeat boundary-------------------------------
            determine_repeat_boundary_v2(repeats_path, longest_repeats_path, blast_program_dir,
                                         fixed_extend_base_threshold, max_repeat_len, tmp_output_dir, threads)

            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of generating longest_repeats.fa: %.8s s" % (dtime))
        break

def test_connect_RepeatMasking_results(tmp_output_dir):
    log = Logger('HiTE.log', level='debug')
    ref_index = 0

    longest_repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.fa'
    repeats_path = tmp_output_dir + '/repeats_' + str(ref_index) + '.fa'

    #连接repeats序列
    avg_repeat_len = 500000
    repeat_names, repeat_contigs = read_fasta(repeats_path)

    connected_repeats = {}
    node_index = 0
    cur_repeat_len = 0
    cur_repeat = ''
    for name in repeat_names:
        seq = repeat_contigs[name]
        cur_repeat_len += len(seq)
        cur_repeat += seq
        if cur_repeat_len >= avg_repeat_len:
            connected_repeats['N_'+str(node_index)] = cur_repeat
            node_index += 1
            cur_repeat_len = 0
            cur_repeat = ''
    if cur_repeat != '':
        connected_repeats['N_' + str(node_index)] = cur_repeat

    repeats_path = tmp_output_dir + '/connected_repeats_' + str(ref_index) + '.fa'
    store_fasta(connected_repeats, repeats_path)

    blast_program_dir = '/home/hukang/repeat_detect_tools/rmblast-2.9.0-p2'
    fixed_extend_base_threshold = 1000
    max_repeat_len = 30000
    threads = 40
    starttime = time.time()
    log.logger.info('Start 2.2: Coarse-grained boundary mapping')
    log.logger.info('------generate longest_repeats.fa')

    # -------------------------------Stage02: this stage is used to do pairwise comparision, determine the repeat boundary-------------------------------
    determine_repeat_boundary_v1(repeats_path, longest_repeats_path, blast_program_dir,
                                 fixed_extend_base_threshold, max_repeat_len, tmp_output_dir,
                                 threads)

    endtime = time.time()
    dtime = endtime - starttime
    log.logger.info("Running time of generating longest_repeats.fa: %.8s s" % (dtime))

def flanking_copies(all_copies, query_path, reference, flanking_len, copy_num=10, query_coverage=0.99):
    new_all_copies = {}
    tsd_info = {}
    query_names, query_contigs = read_fasta(query_path)
    ref_names, ref_contigs = read_fasta(reference)
    for query_name in all_copies.keys():
        tsd_copy_num = 0
        total_copy_len = 0
        query_seq = query_contigs[query_name]
        copies = all_copies[query_name]
        new_copies = []
        # 取最多copy_num条
        for i, copy in enumerate(copies):
            if copy_num != -1 and i >= copy_num:
                break
            ref_name = copy[0]
            copy_ref_start = int(copy[1])
            copy_ref_end = int(copy[2])
            if copy_ref_start > copy_ref_end:
                tmp = copy_ref_start
                copy_ref_start = copy_ref_end
                copy_ref_end = tmp
            copy_len = copy_ref_end - copy_ref_start + 1
            if copy_ref_start - 1 - flanking_len < 0 or copy_ref_end + flanking_len > len(ref_contigs[ref_name]):
                continue
            ref_seq = ref_contigs[ref_name]
            orig_copy_seq = ref_seq[copy_ref_start-1: copy_ref_end]
            copy_seq = ref_seq[copy_ref_start-1-flanking_len: copy_ref_end+flanking_len]
            #如果是取全长拷贝，则需要判断拷贝的前5bp和后5bp是否与原始序列高度相似，只允许1 bp mismatch
            if query_coverage != 0.99 or \
                    (allow_mismatch(query_seq[0:5], orig_copy_seq[0:5], 1)
                     and allow_mismatch(query_seq[-5:], orig_copy_seq[-5:], 1)):
                new_copies.append((ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq))

            total_copy_len += len(copy_seq)
            tir_start = flanking_len + 1  # (坐标都是以1开始的，start和end都是取到的)
            tir_end = len(copy_seq) - flanking_len
            left_tsd_seq, right_tsd_seq = TSDsearch_v4(copy_seq, tir_start, tir_end)
            if left_tsd_seq != '':
                tsd_copy_num += 1
            copy_name = str(copy[0]) + '-' + str(copy[1]) + '-' + str(copy[2]) + '-' + str(copy_ref_end-copy_ref_start+1)
            if not tsd_info.__contains__(query_name):
                tsd_info[query_name] = {}
            info = tsd_info[query_name]
            info[copy_name] = left_tsd_seq + ',' + right_tsd_seq + ',' + copy_seq
        if not tsd_info.__contains__(query_name):
            tsd_info[query_name] = {}
        info = tsd_info[query_name]
        info['total_copy_num'] = len(copies)
        info['tsd_copy_num'] = tsd_copy_num
        info['total_copy_len'] = total_copy_len

        if len(new_copies) > 0:
            new_all_copies[query_name] = new_copies
    return new_all_copies, tsd_info

def get_seq_copies(input, tmp_output_dir):
    threads = 40
    reference = tmp_output_dir + '/all.chrs.con'

    blast_program_dir = '/home/hukang/repeat_detect_tools/rmblast-2.9.0-p2'
    temp_dir = tmp_output_dir + '/repbase_blast'
    blastnResults_path = temp_dir + '/temp.out'
    multi_process_align(input, reference, blastnResults_path, blast_program_dir, temp_dir, threads)
    all_copies = get_copies(blastnResults_path, input, reference, threads=threads)

    # 在copies的两端 flanking 20bp的序列
    flanking_len = 0
    all_copies, tsd_info = flanking_copies(all_copies, input, reference, flanking_len, copy_num=-1)

    copy_info_path = tmp_output_dir + '/sMITE.copies.info'
    store_copies(tsd_info, copy_info_path)

    copy_info_path = tmp_output_dir + '/sMITE.copies.fa'
    store_copies_seq(all_copies, copy_info_path)


    # 输出无拷贝和单拷贝序列个数
    delete_repbase_names = set()
    no_copy_num = 0
    single_copy_num = 0
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        if len(copies) == 1:
            single_copy_num += 1
            delete_repbase_names.add(query_name)
        elif len(copies) == 0:
            no_copy_num += 1
            delete_repbase_names.add(query_name)
    print('no_copy_num: ' + str(no_copy_num) + ', single_copy_num: ' + str(single_copy_num))


def lost_TIRs(tmp_output_dir):
    tools_dir = os.getcwd() + '/../tools'
    raw_tirs = tmp_output_dir + '/tir_tsd_0.filter_tandem.fa'
    raw_tir_consensus = tmp_output_dir + '/tir_tsd_0.cons.fa'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + raw_tirs + ' -o ' + raw_tir_consensus + ' -T 0 -M 0'
    #os.system(cd_hit_command)

    rename_tir_consensus = tmp_output_dir + '/tir_tsd_0.rename.cons.fa'
    #rename_fasta(raw_tir_consensus, rename_tir_consensus)

    repbase_path = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/repbase/rice/tir.repbase.ref'
    repbase_names, repbase_contigs = read_fasta(repbase_path)

    #1.提取未过滤多出的100个perfect序列
    raw_perfect_path = '/homeb/hukang/KmerRepFinder_test/library/raw_TIRs/perfect.families'
    filter_perfect_path = '/homeb/hukang/KmerRepFinder_test/library/filtered_TIRs/perfect.families'

    raw_perfect = set()
    filter_perfect = set()
    with open(raw_perfect_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            raw_perfect.add(line)
    with open(filter_perfect_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            filter_perfect.add(line)
    lost_tirs = raw_perfect.difference(filter_perfect)
    print('lost tir size: ' + str(len(lost_tirs)))

    lost_tirs_path = tmp_output_dir + '/lost_tir.fa'
    with open(lost_tirs_path, 'w') as f_save:
        for name in lost_tirs:
            f_save.write('>'+name+'\n'+repbase_contigs[name]+'\n')

    #2.将序列比对到参考上，获取拷贝
    get_seq_copies(lost_tirs_path, tmp_output_dir)


def get_fake_tirs(itrsearch_log):
    fake_tirs = set()
    alignments = {}
    line_count = 0
    query_name = ''
    with open(itrsearch_log, 'r') as f_r:
        for line in f_r:
            line_count += 1
            if line.startswith('load sequence'):
                parts = line.split('\t')
                query_name = parts[0].split(' ')[3]
                line_count = 0
            if line_count == 3 or line_count == 4 or line_count == 5:
                if line.strip() == '':
                    continue
                if query_name != '':
                    if not alignments.__contains__(query_name):
                        alignments[query_name] = []
                    details = alignments[query_name]
                    details.append(line)

    for query_name in alignments.keys():
        details = alignments[query_name]
        if len(details) != 3:
            continue
        query_seq = details[0]
        align_seq = details[1]
        target_seq = details[2]
        print(query_name)
        query_parts = query_seq.split(' ')
        target_parts = target_seq.split(' ')
        if len(query_parts) > 7 and len(target_parts) > 7:
            if query_seq[8] == '-' or target_seq[8] == '-' or (align_seq[8] != '|' and align_seq[9] != '|'):
                fake_tirs.add(query_name)
                print('Fake')
            print(query_seq)
            print(align_seq)
            print(target_seq)
            print(query_seq[8])
            print(align_seq[8])
            print(target_seq[8])
    return fake_tirs


def parse_novel_tir_locations(tmp_output_dir):
    # 解析出novel tir对应的染色体坐标
    novel_tir_copies = tmp_output_dir + '/novel_tir_copies.csv'
    file = open(novel_tir_copies, 'r')
    js = file.read()
    all_copies = json.loads(js)

    novel_tir_locations = tmp_output_dir + '/novel_tir_locations.txt'
    chr_names = set()
    #{query_name: set()}
    matrix = {}
    query_names = all_copies.keys()
    with open(novel_tir_locations, 'w') as f_save:
        for query_name in query_names:
            # 取所有的chr name
            new_copies = all_copies[query_name]
            for item in new_copies:
                chr_names.add(item[0].strip())

        chr_names = list(chr_names)
        header = 'id\t'
        for chr_name in chr_names:
            header += chr_name+'\t'
        f_save.write(header + '\n')

        for query_name in query_names:
            # 取query_name <-> chr_name的矩阵
            new_copies = all_copies[query_name]
            if not matrix.__contains__(query_name):
                matrix[query_name] = set()
            chr_set = matrix[query_name]
            for item in new_copies:
                chr_set.add(item[0])

        for query_name in query_names:
            f_save.write(query_name+'\t')
            chr_set = matrix[query_name]
            for chr_name in chr_names:
                if chr_name in chr_set:
                    f_save.write(str(1) + '\t')
                else:
                    f_save.write(str(0) + '\t')
            f_save.write('\n')
    # new_copies.append((ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq))


def get_length_dist(paths, labels, my_pal, output_path):
    with open(output_path, 'w') as f_save:
        f_save.write('length\tmethods\n')
        for i, path in enumerate(paths):
            label = labels[i]
            names, contigs = read_fasta(path)
            for name in names:
                seq_len = len(contigs[name])
                # if seq_len >= 5000:
                #     continue
                f_save.write(str(seq_len)+'\t'+str(label)+'\n')
    df = pd.read_csv(output_path, sep='\t', encoding='utf-8')
    print(df)
    sns.set_context("paper", rc={"font.size":20,"axes.titlesize":8,"axes.labelsize":20})
    b = sns.violinplot(x=df["methods"], y=df["length"], palette=my_pal)
    b.set_xlabel("", fontsize=20)
    b.set_ylabel("Length", fontsize=20)
    plt.tick_params(axis='x', which='major', labelsize=15)
    #plt.show()


    # Calculate number of obs per group & median to position labels
    medians = df.groupby(['methods'])['length'].median().values
    print(medians)

def get_length_dist_boxplot(paths, labels, my_pal, output_path, output_fig):
    with open(output_path, 'w') as f_save:
        f_save.write('length\tmethods\n')
        for i, path in enumerate(paths):
            label = labels[i]
            names, contigs = read_fasta(path)
            for name in names:
                seq_len = len(contigs[name])
                f_save.write(str(seq_len)+'\t'+str(label)+'\n')
    df = pd.read_csv(output_path, sep='\t', encoding='utf-8')
    print(df)

    # 自定义颜色列表
    my_colors = ["#293991", "#9DCFEF"]
    # 设置自定义颜色列表
    sns.set_palette(my_colors)

    # 创建子图布局
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    # 绘制 HiTE-TIR 方法的长度分布
    sns.distplot(df[df['methods'] == labels[0]]['length'], kde=False, bins=5, ax=ax[0], rug=True, hist_kws={'edgecolor': 'black'})
    ax[0].set_title(labels[0]+' Length Distribution')
    ax[0].set_xlabel('Length')
    ax[0].set_ylabel('Count')
    # 绘制 Maize Library 方法的长度分布
    sns.distplot(df[df['methods'] == labels[1]]['length'], kde=False, bins=5, ax=ax[1], rug=True, hist_kws={'edgecolor': 'black'})
    ax[1].set_title(labels[1]+' Length Distribution')
    ax[1].set_xlabel('Length')
    ax[1].set_ylabel('Count')
    # 调整子图之间的间距
    plt.tight_layout()
    # 显示图形
    plt.savefig(output_fig, format='png')

    # y = list(query_copy_num)
    # x = pd.Series(y, name="copy number")
    # sns.set_theme(style="ticks", font='Times New Roman', font_scale=1.4)
    # sns.set_context("paper")
    # ax = sns.distplot(x, kde=True, color='green')
    # # ax.set_xlim(5, )
    # plt.show()

    # sns.set_context("paper", rc={"font.size":20,"axes.titlesize":8,"axes.labelsize":20})
    # #b = sns.violinplot(x=df["methods"], y=df["length"], palette=my_pal)
    # #b = sns.boxplot(x=df["methods"], y=df["length"], palette=my_pal)
    # b = sns.barplot(x=df["methods"], y=df["length"])
    # b.set_xlabel("", fontsize=20)
    # b.set_ylabel("Length", fontsize=20)
    # plt.tick_params(axis='x', which='major', labelsize=15)
    # #plt.show()
    # plt.tight_layout()
    # plt.savefig(output_fig, format='png')

    # Calculate number of obs per group & median to position labels
    medians = df.groupby(['methods'])['length'].median().values
    print(medians)

def get_copies_dist_boxplot(paths, copies_dirs, labels, output_path, output_fig):
    with open(output_path, 'w') as f_save:
        f_save.write('copies\tmethods\n')
        for i, path in enumerate(paths):
            contigNames, contigs = read_fasta(path)
            for file_name in contigNames:
                file_path = copies_dirs[i] + '/' + file_name + '.blast.bed.fa'
                names, contigs = read_fasta(file_path)
                f_save.write(str(len(names))+'\t'+str(labels[i])+'\n')
    df = pd.read_csv(output_path, sep='\t', encoding='utf-8')
    print(df)

    # 创建子图布局
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    # 绘制 HiTE-TIR 方法的长度分布
    sns.distplot(df[df['methods'] == labels[0]]['copies'], kde=False, bins=5, ax=ax[0], rug=True,
                 hist_kws={'edgecolor': 'black'})
    ax[0].set_title(labels[0] + ' Copies Distribution')
    ax[0].set_xlabel('Copies')
    ax[0].set_ylabel('Count')
    # 绘制 Maize Library 方法的长度分布
    sns.distplot(df[df['methods'] == labels[1]]['copies'], kde=False, bins=5, ax=ax[1], rug=True,
                 hist_kws={'edgecolor': 'black'})
    ax[1].set_title(labels[1] + ' Copies Distribution')
    ax[1].set_xlabel('Copies')
    ax[1].set_ylabel('Count')
    # 调整子图之间的间距
    plt.tight_layout()
    # 显示图形
    plt.savefig(output_fig, format='png')

    # Calculate number of obs per group & median to position labels
    medians = df.groupby(['methods'])['copies'].median().values
    print(medians)

def get_insert_time_dist_boxplot(output_path, output_fig):
    df = pd.read_csv(output_path, sep='\t', encoding='utf-8')
    print(df)

    # 自定义颜色列表
    my_colors = ["#293991", "#9DCFEF"]
    # 设置自定义颜色列表
    sns.set_palette(my_colors)

    # 创建子图布局
    #fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    # 绘制 HiTE-TIR 方法的长度分布
    ax = sns.distplot(df[df['methods'] == 'HiTE-TIR']['insert_time'], kde=False, bins=50, rug=False,
                 hist_kws={'edgecolor': 'black'})
    ax.set_title('HiTE-TIR' + ' Insertion Time Distribution')
    ax.set_xlabel('Mya (million years ago)')
    ax.set_ylabel('Count')
    # 调整子图之间的间距
    #plt.tight_layout()
    # 显示图形
    plt.savefig(output_fig, format='png')

def test_remove_copy1_tirs(tmp_output_dir):
    query_copy_num_path = tmp_output_dir + '/novel_tir_copies_num.csv'
    file = open(query_copy_num_path, 'r')
    js = file.read()
    query_copy_num = json.loads(js)
    count = 0
    for name in query_copy_num.keys():
        num = query_copy_num[name]
        if num == 1:
            count += 1
    print(count)

    confident_tir_path = tmp_output_dir + '/confident_tir.rename.cons.fa'
    contignames, contigs = read_fasta(confident_tir_path)
    for name in query_copy_num.keys():
        num = query_copy_num[name]
        if num == 1:
            del contigs[name]

    new_confident_tir_path = tmp_output_dir + '/confident_tir.copy2.fa'
    store_fasta(contigs, new_confident_tir_path)


def generate_seq_logos(tmp_output_dir):
    labels = ['LTR', 'Helitron', 'EnSpm']
    for label in labels:
        output_path = tmp_output_dir + '/logos_'+label+'.txt'
        logos = {}
        logo_len = 30
        if label == 'LTR':
            TE_path = tmp_output_dir + '/confident_ltr_cut.terminal.fa'
        elif label == 'Helitron':
            TE_path = tmp_output_dir + '/confident_helitron_0.rename.cons.fa'
            helitron_names, helitron_contigs = read_fasta(TE_path)
            new_helitron_contigs = {}
            for name in helitron_names:
                seq = helitron_contigs[name]
                new_helitron_contigs[name+'#Helitron'] = seq
            store_fasta(new_helitron_contigs, TE_path)
        else:
            TE_path = tmp_output_dir + '/confident_TE.cons.fa.final.classified'
        names, contigs = read_fasta(TE_path)
        for name in names:
            class_name = name.split('#')[1]
            if class_name.__contains__(label):
                #取序列的开始和结束序列
                seq = contigs[name]
                seq_start = seq[0:logo_len]
                seq_end = seq[-logo_len:]
                if not logos.__contains__('5\'-'+label):
                    logos['5\'-'+label] = []
                seqs = logos['5\'-'+label]
                seqs.append(seq_start)

                if not logos.__contains__('3\'-'+label):
                    logos['3\'-'+label] = []
                seqs = logos['3\'-'+label]
                seqs.append(seq_end)
        title_names = logos.keys()
        with open(output_path, 'w') as f_save:
            for logo in title_names:
                f_save.write(logo + '\t')
            f_save.write('\n')
            line_num = 0
            stop = False
            while(not stop):
                col_finish_num = 0
                for name in title_names:
                    seqs = logos[name]
                    if line_num >= len(seqs):
                        col_finish_num += 1
                        f_save.write(' \t')
                    else:
                        f_save.write(seqs[line_num]+'\t')
                if col_finish_num >= len(title_names):
                    stop = True
                f_save.write('\n')
                line_num += 1


def test_split_genome():
    tmp_output_dir = '/home/hukang/HiTE-2.0.1/demo'
    reference = tmp_output_dir + '/genome.fa'
    chunk_size = 4

    reference_pre = convertToUpperCase_v1(reference)
    ref_names, ref_contigs = read_fasta(reference_pre)

    cut_references = []
    cur_ref_contigs = {}
    ref_index = 0
    single_batch_size = chunk_size * 1024 * 1024
    start = 0
    cur_seq = ''
    for ref_name in ref_names:
        seq = ref_contigs[ref_name]
        cur_ref_name = ref_name + '$' + str(start)
        while len(seq) > 0:
            seq_len = len(seq)
            if len(cur_seq) + seq_len <= single_batch_size:
                cur_seq += seq
                cur_ref_contigs[cur_ref_name] = cur_seq
                cur_seq = ''
                seq = ''
                start = 0
            else:
                remain_size = single_batch_size - len(cur_seq)
                remain_seq = seq[0: remain_size]
                cur_ref_contigs[cur_ref_name] = remain_seq
                start += remain_size
                cur_seq = ''
                # store references
                cur_ref_path = reference + '.cut' + str(ref_index) + '.fa'
                store_fasta(cur_ref_contigs, cur_ref_path)
                cut_references.append(cur_ref_path)
                cur_ref_contigs = {}
                ref_index += 1

                cur_ref_name = ref_name + '$' + str(start)
                seq = seq[remain_size:]

    if len(cur_ref_contigs) > 0:
        cur_ref_path = reference + '.cut' + str(ref_index) + '.fa'
        store_fasta(cur_ref_contigs, cur_ref_path)
        cut_references.append(cur_ref_path)

    blast_program_dir = ''
    if blast_program_dir == '':
        (status, blast_program_path) = subprocess.getstatusoutput('which makeblastdb')
        blast_program_dir = os.path.dirname(os.path.dirname(blast_program_path))

    fixed_extend_base_threshold = 1000
    max_repeat_len = 30000
    threads = 40
    for ref_index, cut_reference in enumerate(cut_references):
        print(cut_reference)
        (cut_ref_dir, cut_ref_filename) = os.path.split(cut_reference)
        (cut_ref_name, cut_ref_extension) = os.path.splitext(cut_ref_filename)


        longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.flanked.fa'
        longest_repeats_path = tmp_output_dir + '/longest_repeats_' + str(ref_index) + '.fa'
        resut_file = longest_repeats_path

        repeats_path = cut_reference
        # -------------------------------Stage02: this stage is used to do pairwise comparision, determine the repeat boundary-------------------------------
        determine_repeat_boundary_v3(repeats_path, longest_repeats_path, blast_program_dir,
                                     fixed_extend_base_threshold, max_repeat_len, tmp_output_dir, threads)


def generate_MITE_identity_dist(sMITE_path, Hi_TIR_Ghd2_path, tmp_output_dir, dist_path):
    TRsearch_dir = '/home/hukang/HiTE-2.0.1/ReferenceMode/tools'
    sMITE_TR_out, sMITE_TR_log = run_itrsearch(TRsearch_dir, sMITE_path, tmp_output_dir)
    Hi_TIR_Ghd2_TR_out, Hi_TIR_Ghd2_TR_log = run_itrsearch(TRsearch_dir, Hi_TIR_Ghd2_path, tmp_output_dir)

    sMITE_identities = []
    Hi_TIR_Ghd2_identities = []

    with open(sMITE_TR_log, 'r') as f_r:
        for line in f_r:
            if line.__contains__('Identity percentage : '):
                identity = float(line.split('Identity percentage :')[1].strip())
                sMITE_identities.append(identity)

    with open(Hi_TIR_Ghd2_TR_log, 'r') as f_r:
        for line in f_r:
            if line.__contains__('Identity percentage : '):
                identity = float(line.split('Identity percentage :')[1].strip())
                Hi_TIR_Ghd2_identities.append(identity)

    print(len(sMITE_identities), len(Hi_TIR_Ghd2_identities))
    with open(dist_path, 'w') as f_save:
        f_save.write('identity' + '\t' + 'Type\n')
        for identity in sMITE_identities:
            f_save.write(str(identity)+'\t'+'sMITE\n')
        for identity in Hi_TIR_Ghd2_identities:
            f_save.write(str(identity)+'\t'+'Hi_TIR_Ghd2\n')


def draw_violin(dist_path, my_pal):
    df = pd.read_csv(dist_path, sep='\t', encoding='utf-8')
    print(df)
    sns.violinplot(x=df["Type"], y=df["identity"], palette=my_pal)
    #plt.show()
    plt.savefig('/home/hukang/Novel_TIR_Ghd2.png', format='png')

    # Calculate number of obs per group & median to position labels
    medians = df.groupby(['Type'])['identity'].median().values
    print(medians)


def draw_stripplot():
    fig, ax = plt.subplots()
    #设置风格
    #sns.set_style('whitegrid')
    sns.set(context="notebook", style='whitegrid', font_scale=1)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    # 构建数据
    tips=pd.read_csv('/home/hukang/nextflow_runtime.csv')
    print(tips)
    """
    案例2：
    根据x的类别进行分组统计
    """

    sns.stripplot(x="process",y="Execution time (minutes)", data=tips, ax=ax, s=15, linewidth=1, palette=sns.hls_palette(8))
    ax.grid(True)
    plt.xticks(rotation=60)
    plt.tight_layout()
    plt.show()
    #plt.savefig('/home/hukang/nextflow_runtime.png', format='png')

def darw_barplot(input):
    #sns.set(context="notebook", style='whitegrid', font_scale=1)
    sns.set_style('whitegrid')
    result=pd.read_csv(input)
    #print(result)
    bar_plot = sns.barplot(y = result['time'].unique(), x = -result[result['species'] == 'Rice']['number'], color = "DarkSalmon",
                       data = result, order = result['time'].unique()[::-1],)
    bar_plot = sns.barplot(y = result['time'].unique(), x = result[result['species'] == 'Maize']['number'], color = "CadetBlue",
                        data = result, order = result['time'].unique()[::-1],)
    #plt.xticks([-650,-400,-200,0,200,400,650],[650,200,100,0,100,200,250])
    # plt.rcParams['font.sans-serif'] = ['SimHei']
    #plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
    plt.rcParams['axes.unicode_minus'] = True
    bar_plot.set(xlabel="Number of LTR (log)", ylabel="Mya", title = "")
    #plt.show()
    plt.savefig('/home/hukang/LTR_insert_time.png', format='png')

def generate_insert_time(ltr_file, type):
    #将LTR_retriever的插入时间统计到文件
    times = []
    with open(ltr_file, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            if line.startswith('#'):
                continue
            elif parts[9] == type:
                time = int(parts[11])
                times.append(time)
    f_r.close()
    #print(max(times))

    #按照20W年分组
    u = 2
    bin = u*100000

    #1.先确定最大年龄
    max_time = int(max(times)/bin) + 1
    #print(max_time)
    #2.遍历分组
    time_group = {}
    for t in times:
        g_num = float(int(t/bin)*u/10)
        if not time_group.__contains__(g_num):
            time_group[g_num] = 0
        cur_num = time_group[g_num]
        time_group[g_num] = cur_num + 1
    return time_group

def generate_insertion_time(type):
    output = '/home/hukang/insert_time.csv'
    ltr_file = '/homeb/hukang/KmerRepFinder_test/library/all_tools_run_lib/rice_v7/HiTE/all.chrs.rename.fa.pass.list'
    speices1 = 'Rice'
    time_group1 = generate_insert_time(ltr_file, type)
    #print(time_group1)

    ltr_file = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test2/maize/genome.rename.fa.pass.list_3.3e-08'
    speices2 = 'Maize'
    time_group2 = generate_insert_time(ltr_file, type)
    #print(time_group2)

    max_insert_time = 3.6
    line1 = []
    with open(output, 'w') as f_save:
        f_save.write('time,species,number\n')
        #对time进行填充，补充那些为空的数据
        times = [round(num,1) for num in np.arange(0,max_insert_time+0.1, 0.2)]
        for t in times:
            if not time_group1.__contains__(t):
                time_group1[t] = 0
            if not time_group2.__contains__(t):
                time_group2[t] = 0

        keys = sorted(time_group1.keys(), reverse=False)
        out_max_insert_count1 = 0
        for i, g_num in enumerate(keys):
            if g_num >= max_insert_time:
                out_max_insert_count1 += time_group1[g_num]
            else:
                if i < len(keys)-1:
                    next_g_num = str(keys[i+1])
                else:
                    next_g_num = ''
                num = time_group1[g_num]
                # 转化为Log
                if num >= 1:
                    num = math.log(num)
                line = str(g_num)+'-'+str(next_g_num)+','+speices1+','+str(num)+'\n'
                line1.append(line)
                f_save.write(line)
        # 转化为Log
        if out_max_insert_count1 >= 1:
            out_max_insert_count1 = math.log(out_max_insert_count1)
        line = str(max_insert_time) + '-' + str('') + ',' + speices1 + ',' + str(out_max_insert_count1) + '\n'
        line1.append(line)
        f_save.write(line)

        line2 = []
        out_max_insert_count2 = 0
        keys = sorted(time_group2.keys(), reverse=False)
        for i, g_num in enumerate(keys):
            if g_num >= max_insert_time:
                out_max_insert_count2 += time_group2[g_num]
            else:
                if i < len(keys)-1 and g_num<max_insert_time:
                    next_g_num = str(keys[i+1])
                else:
                    next_g_num = ''
                num = time_group2[g_num]
                # 转化为Log
                if num >= 1:
                    num = math.log(num)
                line = str(g_num)+'-'+str(next_g_num)+','+speices2+','+str(num)+'\n'
                line2.append(line)
                f_save.write(line)
        if out_max_insert_count2 >= 1:
            out_max_insert_count2 = math.log(out_max_insert_count2)
        line = str(max_insert_time) + '-' + str('') + ',' + speices2 + ',' + str(out_max_insert_count2) + '\n'
        line2.append(line)
        f_save.write(line)
    print(len(line1), len(line2))
    return output

def get_cd_hit_cluster(cluster_file):
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
                name = name[0: len(name) - 3]
                seq_len = int(line.split(',')[0].split('\t')[1].replace('nt', ''))
                cur_cluster.append((seq_len, name))
                if line.endswith('*'):
                    clusters['rep_' + str(cluster_idx)] = name
    f_r.close()
    return clusters

def clean_nested_rep(raw_input, temp_dir, threads, test_home, TRsearch_dir, cur_round):
    #去除嵌合步骤：
        # 1. 找到簇中长度最短的序列x，遍历类簇中的每一条序列,如果序列长度超过2x, 则归为潜在嵌合TE，否则为普通TE
        # 2. 如果簇中存在嵌合TE，则将嵌合TE收集到nested数组中；如果簇中不存在嵌合，则保存Rep序列为代表性序列，其余的序列收集到delete数组中
        # 3. 对原始输入fasta文件去除nested+delete元素。
        # 4. 我们得到nested序列，Rep序列以及剩余序列remain.fa（指的是可能插入到其他TE中的TE）。将nested和Rep序列合并到final nested和final Rep中。
        # 5. 剩余序列做为输入，重新聚类-aS 0.8，进行步骤1-4,直到剩余文件remain.fa为空或者达到最大迭代次数5。
        # 6. 我们最后得到final nested和final Rep文件。将两个文件输入到remove nested中进行去除嵌合。
        # 7. 观察去除嵌合后的序列是否仍包含TIR结构，如果不包含，过滤掉去除嵌合序列；如果包含则加入到Rep中。
        # 8. 对新的Rep序列进行聚类, 迭代1-7
    if cur_round == 0:
        os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    short_coverage = 0.95
    similarity = 0.95
    cons_input = temp_dir + '/raw.cons.fa'
    cd_hit_command = 'cd-hit-est -aS ' + str(short_coverage) + ' -c ' + str(similarity) \
                     + ' -d 0 -G 0 -g 1 -A 80 -i ' + raw_input + ' -o ' + cons_input + ' -T 0 -M 0 ' + ' > /dev/null 2>&1'
    os.system(cd_hit_command)

    max_round = 5
    cur_round = 0
    names, contigs = read_fasta(raw_input)
    final_Rep_seqs = []
    final_nested_seqs = []
    remain_contigs = contigs
    while cur_round < max_round and len(remain_contigs) > 0:
        cluster_file = cons_input + '.clstr'
        clusters = get_cd_hit_cluster(cluster_file)

        Rep_seqs = []
        nested_seqs = []
        delete_seqs = []
        # step: 1-2
        for key in clusters.keys():
            if not key.startswith('rep_'):
                cluster = clusters[key]
                #取cluster中长度最短的序列
                cluster.sort(key=lambda x: x[0])
                min_len = cluster[0][0]
                c_putative_nested = []
                c_normal = []
                for seq_item in cluster:
                    if seq_item[0] >= 2*min_len:
                        c_putative_nested.append(seq_item[1])
                    else:
                        c_normal.append(seq_item[1])
                if len(c_putative_nested) == 0:
                    rep_seq = clusters['rep_'+key]
                    Rep_seqs.append(rep_seq)
                    delete_seqs.extend(c_normal)
                else:
                    nested_seqs.extend(c_putative_nested)
        # step: 3
        #print(len(contigs))
        remove_seqs = nested_seqs + delete_seqs
        cur_contigs = remain_contigs
        for seq_name in remove_seqs:
            del cur_contigs[seq_name]
        #print(Rep_seqs[0:1])
        #print(nested_seqs[0:1])
        #print(len(contigs))
        final_Rep_seqs.extend(Rep_seqs)
        final_nested_seqs.extend(nested_seqs)
        # step: 4-5
        remain_file = temp_dir + '/remain_' + str(cur_round) + '.fa'
        remain_cons = temp_dir + '/remain_' + str(cur_round) + '.cons.fa'
        store_fasta(cur_contigs, remain_file)
        remain_contigs = cur_contigs
        cd_hit_command = 'cd-hit-est -aS ' + str(short_coverage) + ' -c ' + str(similarity) \
                     + ' -d 0 -G 0 -g 1 -A 80 -i ' + remain_file + ' -o ' + remain_cons + ' -T 0 -M 0 ' + ' > /dev/null 2>&1'
        os.system(cd_hit_command)
        cons_input = remain_cons
        cur_round += 1
    # print(len(final_Rep_seqs))
    # print(len(final_nested_seqs))
    names, contigs = read_fasta(raw_input)
    pure_rep = temp_dir + '/pure_tir_rep.fa'
    with open(pure_rep, 'w') as f_save:
        for rep_name in final_Rep_seqs:
            rep_seq = contigs[rep_name]
            f_save.write('>'+rep_name+'\n'+rep_seq+'\n')
    f_save.close()
    nested_file = temp_dir + '/nested_tir.fa'
    with open(nested_file, 'w') as f_save:
        for nested_name in final_nested_seqs:
            nested_seq = contigs[nested_name]
            f_save.write('>'+nested_name+'\n'+nested_seq+'\n')
    f_save.close()

    # step: 6
    clean_nested_path = temp_dir + '/tir_nested.clean.fa'
    remove_nested_command = 'python3 ' + test_home + '/remove_nested_lib.py ' \
                            + ' -t ' + str(threads) \
                            + ' --tmp_output_dir ' + temp_dir + ' --max_iter_num ' + str(5) \
                            + ' --input1 ' + pure_rep \
                            + ' --input2 ' + nested_file \
                            + ' --output ' + clean_nested_path
    print(remove_nested_command)
    os.system(remove_nested_command)
    # names, contigs = read_fasta(clean_nested_path)
    # print(len(contigs))

    # step: 7
    clean_nested_out, clean_nested_log = run_itrsearch(TRsearch_dir, clean_nested_path, temp_dir)      
    # names, contigs = read_fasta(clean_nested_out)
    # print(len(contigs))

    pure_with_clean_file = temp_dir + '/pure_with_clean_tir.fa'
    os.system('cat ' + clean_nested_out + ' > ' + pure_with_clean_file)
    os.system('cat ' + pure_rep + ' >> ' + pure_with_clean_file)

    return pure_with_clean_file

def run_BM_RM2(TE_path, res_out, temp_dir, rm2_script, lib_path):
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    threads = 30
    rm2_command = 'RepeatMasker -lib '+lib_path+' -nolow -pa '+str(threads)+' '+TE_path
    rm2_res_command = 'cd '+temp_dir+' && rm -rf * && sh '+rm2_script+ ' ' + TE_path + '.out > ' + res_out
    os.system(rm2_command)
    print(rm2_res_command)
    os.system(rm2_res_command)

def run_EDTA(out_dir, TE_path, reference, EDTA_home, TE_out, std_TE_out):
    command1 = 'cd '+out_dir+' && RepeatMasker -e ncbi -pa 40 -q -no_is -norna -nolow -div 40 -lib '+TE_path+' -cutoff 225 '+reference+' && mv '+reference+'.out '+TE_out
    command2 = 'cd '+out_dir+' && perl '+EDTA_home+'/lib-test.pl -genome '+reference+' -std '+std_TE_out+' -tst '+TE_out+' -cat Total'
    os.system(command1)
    os.system(command2)

def run_find_members(cur_file, reference, temp_dir, script_path):
    copy_command = 'cd ' + temp_dir + ' && sh ' + script_path + ' ' + reference + ' ' + cur_file + ' 0 0 ' + ' > /dev/null 2>&1'
    os.system(copy_command)
    output_file = cur_file + '.blast.bed.fa'
    names, contigs = read_fasta(cur_file)
    names1, contigs1 = read_fasta(output_file)
    if len(names1) > 1:
        return (names[0], output_file)
    else:
        return (None, None)

def get_mask_regions(std_out):
    chr_masks = {}
    te_masks = {}
    with open(std_out, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            chr_name = parts[0]
            chr_start = int(parts[1])
            chr_end = int(parts[2])
            te_name = parts[3].split(';')[9]
            if chr_name not in chr_masks:
                chr_masks[chr_name] = []
            segs = chr_masks[chr_name]
            segs.append((chr_start, chr_end, te_name))

            if te_name not in te_masks:
                te_masks[te_name] = {}
            chr_dict = te_masks[te_name]
            if chr_name not in chr_dict:
                chr_dict[chr_name] = []
            chr_segs = chr_dict[chr_name]
            chr_segs.append((chr_start, chr_end))
    return chr_masks, te_masks

def BM_EDTA():
    # 自己编写程序，尝试获得和BM_EDTA一致的结果
    genome_dir = '/home/hukang/test/HiTE/demo/ath/BM_EDTA'
    work_dir = '/home/hukang/test/HiTE/demo/ath/BM_EDTA/01_out'
    genome = genome_dir + '/01.col.fa'
    std_out = work_dir + '/repbase.edta.out'
    test_out = work_dir + '/HiTE.edta.out'
    # Convert .out file to .bed file
    convert2bed_command = 'perl ' + cur_dir + '/tools/RMout_to_bed.pl ' + std_out + ' base1'
    os.system(convert2bed_command)
    convert2bed_command = 'perl ' + cur_dir + '/tools/RMout_to_bed.pl ' + test_out + ' base1'
    os.system(convert2bed_command)
    std_out += '.bed'
    test_out += '.bed'
    std_chr_masks, std_te_masks = get_mask_regions(std_out)
    test_chr_masks, test_te_masks = get_mask_regions(test_out)
    # print(std_chr_masks)
    # print(std_te_masks)
    # print(test_chr_masks)
    # print(test_te_masks)

    # 计算TP,FP,FN,TN等指标
    # 我们获取基因组的base，然后分别把比对mask上去，那我们就可以知道哪些碱基是TP和FP了。
    # 由于我们有碱基与 te_name 的对应关系，我们可以获得每个 te 对应的TP,FP数量
    ref_names, ref_contigs = read_fasta(genome)
    base_masked = {}
    for name in ref_names:
        ref_seq = ref_contigs[name]
        ref_arr = [0] * len(ref_seq)
        base_masked[name] = ref_arr
    # 将金标准覆盖碱基，将覆盖区域置为1
    # 将测试覆盖碱基，将有>=1的区域置为2，否则置为-1
    # 目前有 0， 1， 2，-1 四种组成，TP: 2; TN: 0; FP: -1; FN: 1
    for chr_name in std_chr_masks.keys():
        ref_arr = base_masked[chr_name]
        for seg in std_chr_masks[chr_name]:
            chr_start = seg[0]
            chr_end = seg[1]
            for i in range(chr_start - 1, chr_end):
                ref_arr[i] = 1
    for chr_name in test_chr_masks.keys():
        ref_arr = base_masked[chr_name]
        for seg in test_chr_masks[chr_name]:
            chr_start = seg[0]
            chr_end = seg[1]
            for i in range(chr_start - 1, chr_end):
                if ref_arr[i] >= 1:
                    ref_arr[i] = 2
                else:
                    ref_arr[i] = -1
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    for chr_name in base_masked.keys():
        ref_arr = base_masked[chr_name]
        for i in range(len(ref_arr)):
            if ref_arr[i] == 2:
                TP += 1
            elif ref_arr[i] == 0:
                TN += 1
            elif ref_arr[i] == -1:
                FP += 1
            elif ref_arr[i] == 1:
                FN += 1
    print('TP:' + str(TP))
    print('FN:' + str(FN))
    print('TN:' + str(TN))
    print('FP:' + str(FP))

    sens = float(TP) / (TP + FN)
    spec = float(TN) / (FP + TN)
    accu = float(TP + TN) / (TP + TN + FP + FN)
    prec = float(TP) / (TP + FP)
    FDR = float(FP) / (TP + FP)
    F1 = float(2 * TP) / (2 * TP + FP + FN)
    print('Sensitivity: ' + str(sens))
    print('Specificity: ' + str(spec))
    print('Accuracy: ' + str(accu))
    print('Precision: ' + str(prec))
    print('FDR: ' + str(FDR))
    print('F1 measure: ' + str(F1))

    # 计算每个TE 对应的 TP, FP, FN
    FP_ind = {}
    for te_name in test_te_masks.keys():
        TP = 0
        FP = 0
        FN = 0
        chr_dict = test_te_masks[te_name]
        for chr_name in chr_dict.keys():
            ref_arr = base_masked[chr_name]
            for seg in chr_dict[chr_name]:
                chr_start = seg[0]
                chr_end = seg[1]
                for i in range(chr_start - 1, chr_end):
                    if ref_arr[i] == 2:
                        TP += 1
                    elif ref_arr[i] == -1:
                        FP += 1
                    elif ref_arr[i] == 1:
                        FN += 1
        FP_ind[te_name] = (TP, FP, FN)
    sorted_FP_ind = sorted(FP_ind.items(), key=lambda item: item[1][1], reverse=True)

    # 序列化字典并保存到文件
    FP_ind_json = work_dir + '/FP_ind.json'
    with open(FP_ind_json, 'w', encoding='utf-8') as f:
        json.dump(sorted_FP_ind, f, ensure_ascii=False, indent=4)

    # 计算每个 TE 对应的 TP, FP, FN
    FN_ind = {}
    for te_name in std_te_masks.keys():
        TP = 0
        FP = 0
        FN = 0
        chr_dict = std_te_masks[te_name]
        for chr_name in chr_dict.keys():
            ref_arr = base_masked[chr_name]
            for seg in chr_dict[chr_name]:
                chr_start = seg[0]
                chr_end = seg[1]
                for i in range(chr_start - 1, chr_end):
                    if ref_arr[i] == 2:
                        TP += 1
                    elif ref_arr[i] == -1:
                        FP += 1
                    elif ref_arr[i] == 1:
                        FN += 1
        FN_ind[te_name] = (TP, FP, FN)
    sorted_FN_ind = sorted(FN_ind.items(), key=lambda item: item[1][2], reverse=True)

    # 序列化字典并保存到文件
    FN_ind_json = work_dir + '/FN_ind.json'
    with open(FN_ind_json, 'w', encoding='utf-8') as f:
        json.dump(sorted_FN_ind, f, ensure_ascii=False, indent=4)

    # with open(te_ind_json, 'r', encoding='utf-8') as f:
    #     sorted_te_ind = json.load(f)
    # print(sorted_te_ind)

def BM_HiTE():
    # 统计一下哪些序列贡献的FP最多
    work_dir = '/public/home/hpc194701009/ath_pan_genome/pan_genome/ath/panHiTE_serial_output'
    FP_file = work_dir + '/FP.blastn.out'
    FP_count = {}
    with open(FP_file, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            seq_name = parts[0]
            chr_start = int(parts[2])
            chr_end = int(parts[3])
            cover_len = abs(chr_end-chr_start)
            if seq_name not in FP_count:
                FP_count[seq_name] = 0
            cur_cover_len = FP_count[seq_name]
            cur_cover_len += cover_len
            FP_count[seq_name] = cur_cover_len
    sorted_FP_count = sorted(FP_count.items(), key=lambda item: -item[1])

    file_path = work_dir + '/sorted_FP_count.txt'
    # 将排序后的结果写入文件
    with open(file_path, 'w') as file:
        for key, value in sorted_FP_count:
            file.write(f"{key}: {value}\n")

    FN_file = work_dir + '/FN.blastn.out'
    FN_count = {}
    with open(FN_file, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            seq_name = parts[0]
            chr_start = int(parts[2])
            chr_end = int(parts[3])
            cover_len = abs(chr_end - chr_start)
            if seq_name not in FN_count:
                FN_count[seq_name] = 0
            cur_cover_len = FN_count[seq_name]
            cur_cover_len += cover_len
            FN_count[seq_name] = cur_cover_len
    sorted_FN_count = sorted(FN_count.items(), key=lambda item: -item[1])

    file_path = work_dir + '/sorted_FN_count.txt'
    # 将排序后的结果写入文件
    with open(file_path, 'w') as file:
        for key, value in sorted_FN_count:
            file.write(f"{key}: {value}\n")

def is_only_chars(s, target_char):
    return all(c == target_char for c in s)

def run_find_members_v1(cur_file, reference, temp_dir, member_script_path, subset_script_path, TRsearch_dir):
    # step1: 搜索是否具有清晰边界
    # 1.找members，扩展members 10bp,进行多序列比对
    copy_command = 'cd ' + temp_dir + ' && sh ' + member_script_path + ' ' + reference + ' ' + cur_file + ' 0 10 ' + ' > /dev/null 2>&1'
    os.system(copy_command)
    member_file = cur_file + '.blast.bed.fa'
    extend_member_file = cur_file + '.ext.blast.bed.fa'
    os.system('mv ' + member_file + ' ' + extend_member_file)
    member_names, member_contigs = read_fasta(extend_member_file)
    if len(member_names) > 100:
        sub_command = 'cd ' + temp_dir + ' && sh ' + subset_script_path + ' ' + extend_member_file + ' 100 50 ' + ' > /dev/null 2>&1'
        os.system(sub_command)
        extend_member_file +=  '.rdmSubset.fa'
    elif len(member_names) <= 1:
        return None, None
    align_file = cur_file + '.ext.maf.fa'
    align_command = 'cd ' + temp_dir + ' && mafft --quiet --thread 1 ' + extend_member_file + ' > ' + align_file
    os.system(align_command)
    cons_file = cur_file + '.ext.cons.fa'
    # 2.使用默认参数的一致性序列命令生成一致性序列
    cons_command = 'cons -sequence ' + align_file + ' -outseq ' + cons_file + ' > /dev/null 2>&1'
    os.system(cons_command)
    # 3.小写转大写
    cons_names, cons_contigs = read_fasta(cons_file)
    p_cons_contigs = {}
    for name in cons_names:
        seq = cons_contigs[name]
        seq = seq.upper()
        p_cons_contigs[name] = seq
    store_fasta(p_cons_contigs, cons_file)
    cons_seq = p_cons_contigs[cons_names[0]]
    # 4.取原始序列的首尾10bp，分别在一致性序列上搜索，获取在一致性序列上的边界。（迭代搜索）
    cur_names, cur_contigs = read_fasta(cur_file)
    cur_seq = cur_contigs[cur_names[0]]
    first_10bp = cur_seq[0:10]
    last_10bp = cur_seq[-10:]
    #对应序列从0开始的索引
    start_pos = -1
    end_pos = -1
    start_dist = 3
    last_dist = 3
    first_matches = find_near_matches(first_10bp, cons_seq, max_l_dist=start_dist)
    while len(first_matches) == 0:
        start_dist += 1
        first_matches = find_near_matches(first_10bp, cons_seq, max_l_dist=start_dist)
    for f_m in first_matches:
        if f_m.start >= 10:
            start_pos = f_m.start
            break
    last_matches = find_near_matches(last_10bp, cons_seq, max_l_dist=last_dist)
    while len(last_matches) == 0:
        last_dist += 1
        last_matches = find_near_matches(last_10bp, cons_seq, max_l_dist=last_dist)
    for l_m in reversed(last_matches):
        if l_m.end <= len(cons_seq)-10:
            end_pos = l_m.end-1
            break
    if start_pos == -1 or end_pos == -1:
        #没有找到边界，异常情况
        return None, None
    # 5.如果边界前的序列为N、2bp TA、4bp TTAA，则序列为真实TIR边界，再经由下面的流程判断未扩展的一致性序列是否具有TIR结构
    # 取一致性序列的前后2bp, 4bp
    # 5. 取边界前后的10bp，将10b划分成前6bp和后4bp；真实的TIR前6bp至多有1个非N，后4bp有几种情况：全N(允许1bp非N)，NNTA-TANN，NTAA-TAAN, NTTA-TTAN, NTNA-TNAN, TTAA
    cons_first_6bp = cons_seq[start_pos-10: start_pos-4] 
    cons_first_4bp = cons_seq[start_pos-4: start_pos]
    cons_last_6bp = cons_seq[end_pos+5: end_pos+11] 
    cons_last_4bp = cons_seq[end_pos+1: end_pos+5]
    if allow_mismatch(cons_first_6bp, 'NNNNNN', 1) and allow_mismatch(cons_last_6bp, 'NNNNNN', 1):
        if (allow_mismatch(cons_first_4bp, 'NNNN', 1) and allow_mismatch(cons_last_4bp, 'NNNN', 1)) or \
            (allow_mismatch(cons_first_4bp, 'NNTA', 1) and allow_mismatch(cons_last_4bp, 'TANN', 1)) or \
            (allow_mismatch(cons_first_4bp, 'NTAA', 1) and allow_mismatch(cons_last_4bp, 'TAAN', 1)) or \
            (allow_mismatch(cons_first_4bp, 'NTTA', 1) and allow_mismatch(cons_last_4bp, 'TTAN', 1)) or \
            (allow_mismatch(cons_first_4bp, 'NTNA', 1) and allow_mismatch(cons_last_4bp, 'TNAN', 1)) or \
            (allow_mismatch(cons_first_4bp, 'TTAA', 1) and allow_mismatch(cons_last_4bp, 'TTAA', 1)):

            #具有真实TIR边界
            # step2: 搜索一致性序列是否具有TIR结构
            # 1.找members，对members进行多比对
            copy_command = 'cd ' + temp_dir + ' && sh ' + member_script_path + ' ' + reference + ' ' + cur_file + ' 0 0 ' + ' > /dev/null 2>&1'
            os.system(copy_command)
            member_file = cur_file + '.blast.bed.fa'
            member_names, member_contigs = read_fasta(member_file)
            if len(member_names) > 100:
                sub_command = 'cd ' + temp_dir + ' && sh ' + subset_script_path + ' ' + member_file + ' 100 50 ' + ' > /dev/null 2>&1'
                os.system(sub_command)
                member_file +=  '.rdmSubset.fa'
            align_file = cur_file + '.maf.fa'
            align_command = 'cd ' + temp_dir + ' && mafft --quiet --thread 1 ' + member_file + ' > ' + align_file
            os.system(align_command)
            cons_file = cur_file + '.cons.fa'
            # 2.生成一致性序列，每个位点保留出现2次以上的碱基
            cons_command = 'cons -sequence ' + align_file + ' -outseq ' + cons_file + ' -plurality 2 ' + ' > /dev/null 2>&1'
            os.system(cons_command)
            # 3.小写转大写，并去掉一致性序列中的N
            final_cons_file = cur_file + '.final.cons.fa'
            cons_names, cons_contigs = read_fasta(cons_file)
            p_cons_contigs = {}
            for name in cons_names:
                seq = cons_contigs[name]
                seq = seq.upper().replace("N", "")
                p_cons_contigs[name] = seq
            store_fasta(p_cons_contigs, final_cons_file)
            # 4.判断生成的一致性序列是否有TIR结构
            names, contigs = read_fasta(cur_file)
            final_cons_out, final_cons_log = run_itrsearch(TRsearch_dir, final_cons_file, temp_dir)
            final_cons_names, final_cons_contigs = read_fasta(final_cons_out)
            # 返回序列名称和一致性序列
            if len(final_cons_names) >= 1:
                query_name = names[0]
                return query_name, final_cons_contigs[final_cons_names[0]]
            else:
                return None, None
        else:
            return None, None
    else:
        return None, None

def get_no_gap_seq(seq, pos, direct, sub_len):
    count = 0
    res = ''
    i = pos
    char = seq[i]
    while count < sub_len:
        if direct == 'right' and i < len(seq):
            if seq[i] != '-':
                res += seq[i]
                count += 1
            i += 1
        elif direct == 'left' and i >= 0:
            if seq[i] != '-':
                res = seq[i] + res
                count += 1
            i -= 1
    return res


def run_find_members_v2(cur_file, reference, temp_dir, member_script_path, plant, TE_type):
    # if cur_file.__contains__('N_89008'):
    #     print('here')

    # 1.找members，扩展members 20bp,进行多序列比对
    extend_len = 20
    copy_command = 'cd ' + temp_dir + ' && sh ' + member_script_path + ' ' + reference + ' ' + cur_file + ' 0 ' + str(extend_len) + ' > /dev/null 2>&1'
    #print(copy_command)
    os.system(copy_command)
    member_file = cur_file + '.blast.bed.fa'
    extend_member_file = cur_file + '.ext.blast.bed.fa'
    os.system('mv ' + member_file + ' ' + extend_member_file)

    align_file = cur_file + '.ext.maf.fa'
    align_command = 'cd ' + temp_dir + ' && mafft --quiet --thread 1 ' + extend_member_file + ' > ' + align_file
    os.system(align_command)

    # 2. 根据比对文件取全长的拷贝。
    # 2.1 先定位比对文件中第一条序列的TIR边界位置，作为锚点
    # 取原始序列的首尾20bp，分别在对齐序列上搜索
    # 因为比对文件存在gap，所以我们先在没有gap的序列上搜索位置，然后映射到有gap序列的位置
    cur_names, cur_contigs = read_fasta(cur_file)
    cur_seq = cur_contigs[cur_names[0]]
    first_10bp = cur_seq[0:20]
    last_10bp = cur_seq[-20:]
    align_names, align_contigs = read_fasta(align_file)
    align_start = -1
    align_end = -1
    for name in align_names:
        raw_align_seq = align_contigs[name]
        align_seq = ''
        position_reflex = {}
        cur_align_index = 0
        for i, base in enumerate(raw_align_seq):
            if base == '-':
                continue
            else:
                align_seq += base
                position_reflex[cur_align_index] = i
                cur_align_index += 1

        start_dist = 1
        last_dist = 1
        first_matches = find_near_matches(first_10bp, align_seq, max_l_dist=start_dist)
        last_matches = find_near_matches(last_10bp, align_seq, max_l_dist=last_dist)
        last_matches = last_matches[::-1]
        if len(first_matches) > 0 and len(last_matches) > 0:
            align_no_gap_start = first_matches[0].start
            align_no_gap_end = last_matches[0].end - 1
            #需要将位置映射到原始比对序列上
            align_start = position_reflex[align_no_gap_start]
            align_end = position_reflex[align_no_gap_end]
            break
    # 2.2 其余全长序列在锚点位置上下2bp应该有碱基
    #重新遍历，找到全长拷贝
    full_member_contigs = {}
    for name in align_names:
        # 为了减少计算量，只取100条全长拷贝
        if len(full_member_contigs) > 100:
            break
        align_seq = align_contigs[name]
        if align_start-2 >= 0 and align_end + 3 < len(align_seq):
            anchor_start_seq = align_seq[align_start-2: align_start+3]
            anchor_end_seq = align_seq[align_end - 2: align_end + 3]
            if anchor_start_seq != '-----' and anchor_end_seq != '-----':
                full_member_contigs[name] = align_seq
    full_align_file = cur_file + '.full.maf.fa'
    store_fasta(full_member_contigs, full_align_file)

    # 3. 根据全长拷贝来确定序列是否为真实TIR或Helitron
    # 全长拷贝为1的Helitron不敢认为是真实的Helitron，因为数量太多了
    if len(full_member_contigs) < 1:
        return None, None
    elif len(full_member_contigs) == 1:
        if TE_type == 'TIR':
            # 如果只有一条全长拷贝，但是满足一些固定特征的序列，我认为是真实TIR
            # 包括TSD>=6 或者
            # animal/fungi中的5'-CCC...GGG-3', 3-> plant中的5'-CACT(A/G)...(C/T)AGTG-3'
            query_name = cur_names[0]
            align_seq = full_member_contigs.get(next(iter(full_member_contigs)))
            tsd_len = int(query_name.split('-tsd_')[1].split('-')[0])
            left_tsd = get_no_gap_seq(align_seq, align_start - 1, 'left', tsd_len)
            first_3bp = get_no_gap_seq(align_seq, align_start, 'right', 3)
            last_3bp = get_no_gap_seq(align_seq, align_end, 'left', 3)
            first_5bp = get_no_gap_seq(align_seq, align_start, 'right', 5)
            last_5bp = get_no_gap_seq(align_seq, align_end, 'left', 5)
            if tsd_len >= 6 or (plant == 0 and first_3bp == 'CCC' and last_3bp == 'GGG') \
                    or (plant == 1 and ((first_5bp == 'CACTA' and last_5bp == 'TAGTG')
                                        or (first_5bp == 'CACTG' and last_5bp == 'CAGTG'))):
                return cur_names[0], cur_seq
            else:
                return None, None
        elif TE_type == 'Helitron':
            return None, None
    elif len(full_member_contigs) < 10:
        #全长拷贝较少的序列，挨个获取边界外10bp非空区域。循环遍历，任意两条序列，只要在边界外出现了高同源性，则认为是假阳性序列。
        all_first_10bp_seqs = []
        all_last_10bp_seqs = []
        for name in full_member_contigs.keys():
            seq = full_member_contigs[name]
            #往前遍历，取10个非空组成的序列
            count = 0
            first_10bp_no_gap = ""
            i = align_start-1
            while count < 10 and i >= 0:
                if seq[i] != '-':
                    first_10bp_no_gap = seq[i] + first_10bp_no_gap
                    count += 1
                i -= 1
            all_first_10bp_seqs.append(first_10bp_no_gap)
            #往后遍历，取10个非空组成的序列
            count = 0
            last_10bp_no_gap = ""
            i = align_end + 1
            while count < 10 and i < len(seq):
                if seq[i] != '-':
                    last_10bp_no_gap += seq[i]
                    count += 1
                i += 1
            all_last_10bp_seqs.append(last_10bp_no_gap)
        #循环遍历all_first_10bp_seqs和all_last_10bp_seqs
        for i in range(len(all_first_10bp_seqs)):
            cur_first = all_first_10bp_seqs[i]
            cur_last = all_last_10bp_seqs[i]
            for j in range(i + 1, len(all_first_10bp_seqs)):
                next_first = all_first_10bp_seqs[j]
                next_last = all_last_10bp_seqs[j]
                #是否具有高同源性
                start_dist = 1
                last_dist = 1
                if cur_first != '' and next_first != '':
                    first_matches = find_near_matches(cur_first, next_first, max_l_dist=start_dist)
                    if len(first_matches) > 0:
                        return None, None
                if cur_last != '' and next_last != '':
                    last_matches = find_near_matches(cur_last, next_last, max_l_dist=last_dist)
                    if len(last_matches) > 0:
                        return None, None
        #所有拷贝都没有同源性，说明TIR边界很清晰
        return cur_names[0], cur_seq
    else:
        # 1. 使用cialign去除比对中间的gaps
        rm_gap_align_file = full_align_file + '_cleaned.fasta'
        rm_gap_command = 'cd ' + temp_dir + ' && CIAlign --infile ' + full_align_file + ' --outfile_stem ' + full_align_file + ' --remove_insertions ' + ' > /dev/null 2>&1'
        os.system(rm_gap_command)

        # 2. 定位remove gap之后的边界位置
        align_names, align_contigs = read_fasta(rm_gap_align_file)
        align_start = -1
        align_end = -1
        for name in align_names:
            align_seq = align_contigs[name]
            start_dist = 1
            last_dist = 1
            first_matches = find_near_matches(first_10bp, align_seq, max_l_dist=start_dist)
            last_matches = find_near_matches(last_10bp, align_seq, max_l_dist=last_dist)
            last_matches = last_matches[::-1]
            if len(first_matches) > 0 and len(last_matches) > 0:
                align_start = first_matches[0].start
                align_end = last_matches[0].end - 1
                break
        if align_start == -1 or align_end == -1:
            # 没有找到边界，异常情况
            return None, None

        cons_file = cur_file + '.ext.cons.fa'
        # 3.使用默认参数的一致性序列命令生成一致性序列
        cons_command = 'cons -sequence ' + rm_gap_align_file + ' -outseq ' + cons_file + ' > /dev/null 2>&1'
        os.system(cons_command)
        # 4.小写转大写
        cons_names, cons_contigs = read_fasta(cons_file)
        p_cons_contigs = {}
        for name in cons_names:
            seq = cons_contigs[name]
            seq = seq.upper()
            p_cons_contigs[name] = seq
        store_fasta(p_cons_contigs, cons_file)
        if len(cons_names) <= 0:
            # 没有生成一致性序列
            return None, None
        cons_seq = p_cons_contigs[cons_names[0]]

        # 5.根据一致性序列过滤掉具有较高数量的全长拷贝的假阳性序列，边界外有高同源性的序列。
        # 如果边界前20bp内有5个以上的N序列，不存在连续5个非N，边界以内的序列全都是高一致性序列（不存在连续2个N）
        # 说明序列是真实的TE
        if align_start-extend_len < 0:
            cur_start_pos = 0
        else:
            cur_start_pos = align_start-extend_len
        if align_end+extend_len+1 > len(cons_seq):
            cur_end_pos = len(cons_seq)
        else:
            cur_end_pos = align_end+extend_len+1
        cons_first_20bp = cons_seq[cur_start_pos: align_start]
        cons_last_20bp = cons_seq[align_end + 1: cur_end_pos]
        #连续5个非N
        first_contig_5 = re.search("[^N]{5,}", cons_first_20bp)
        last_contig_5 = re.search("[^N]{5,}", cons_last_20bp)
        #边界向内取20bp序列，判断是否有连续2个N
        cons_first_internal_20bp = cons_seq[align_start: align_start+20]
        cons_last_internal_20bp = cons_seq[align_end-19: align_end+1]
        first_contig_2 = re.search("[N]{2,}", cons_first_internal_20bp)
        last_contig_2 = re.search("[N]{2,}", cons_last_internal_20bp)
        if cons_first_20bp.count("N") > 5 and cons_last_20bp.count("N") > 5 \
                and not first_contig_5 and not last_contig_5\
                and not first_contig_2 and not last_contig_2:
            #取一致性序列
            return cur_names[0], cons_seq[align_start: align_end+1]
        else:
            return None, None


def keep_multi_copy(raw_input, output, reference, temp_dir, script_path, threads):
    os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    names, contigs = read_fasta(raw_input)
    split_files = []
    for i, name in enumerate(names):
        cur_file = temp_dir + '/' + str(i) + '.fa'
        cur_contigs = {}
        cur_contigs[name] = contigs[name]
        store_fasta(cur_contigs, cur_file)
        split_files.append(cur_file)

    ex = ProcessPoolExecutor(threads)
    jobs = []
    for ref_index, cur_file in enumerate(split_files):
        job = ex.submit(run_find_members, cur_file, reference, temp_dir, script_path)
        jobs.append(job)
    ex.shutdown(wait=True)

    multi_contigs = {}
    for job in as_completed(jobs):
        cur_name, output_file = job.result()
        if cur_name is not None:
            multi_contigs[cur_name] = contigs[cur_name]
    store_fasta(multi_contigs, output)

def build_lib(raw_input, reference, threads, temp_copies_dir, temp_nested_dir, test_home, TRsearch_dir, script_path, log):
    # 1. 过滤掉单拷贝的TIR序列
    multi = temp_copies_dir + '/multi.fa'
    log.logger.debug('Start filtering single copy TE instances')
    keep_multi_copy(raw_input, multi, reference, temp_copies_dir, script_path, threads)
    log.logger.debug('Finish filtering single copy TE instances')

    log.logger.debug('Start unwrapping nested TE')
    # 2. 使用cd-hit按照80-80-80 rule 聚类，并迭代解开嵌合序列
    max_round = 5
    cur_round = 0
    stable = False
    cur_input = multi
    cur_output = None
    while cur_round < max_round and not stable:
        log.logger.debug('current round: ' + str(cur_round))
        cur_output = clean_nested_rep(cur_input, temp_nested_dir, threads, test_home, TRsearch_dir, cur_round)
        names1, contigs1 = read_fasta(cur_input)
        names2, contigs2 = read_fasta(cur_output)
        if len(names1) == len(names2):
            stable= True
        cur_input = cur_output
        cur_round += 1
    log.logger.debug('Finish unwrapping nested TE')
    return cur_output


def change_LTR_insertion_time(orig_miu, miu):
    ltr_file = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test2/maize/genome.rename.fa.pass.list'
    new_ltr_file = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test2/maize/genome.rename.fa.pass.list_'+str(miu)
    lines = []
    with open(ltr_file, 'r') as f_r:
        for i, line in enumerate(f_r):
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            # d = 1 - float(parts[7])
            # K=  -3/4*math.log(1-d*4/3)
            # T = K/(2*orig_miu)
            # print(line)
            # print(d, K, T)

            #根据插入时间反推一致性序列
            T = float(parts[11])
            K = T*2*orig_miu
            # (-4*K)/3 = ln(1-d*4/3)
            # math.exp((-4*K)/3) = 1-d*4/3
            identity = 1 - (1-math.exp((-4*K)/3))*3/4
            print(identity)
            # 重新计算时间
            d = 1 - identity
            K =  -3/4*math.log(1-d*4/3)
            T = int(round(K/(2*miu),1))
            new_line = ''
            for j, p in enumerate(parts):
                if j != len(parts)-1:
                    new_line += p + '\t'
                else:
                    new_line += str(T) + '\n'
            lines.append(new_line)
            # if i > 10:
            #     break
    print(lines)
    with open(new_ltr_file, 'w') as f_w:
        for line in lines:
            f_w.write(line)
    f_w.close()


def analyze_new_TIRs(tmp_output_dir):
    new_tir_names = identify_new_TIR(tmp_output_dir)
    novel_tir_path = tmp_output_dir + '/novel_tir.fa'
    novel_tir_names, novel_tir_contigs = read_fasta(novel_tir_path)
    novel_tir_names.sort(key=lambda x: int(x.split('_')[1]))

    #获取TIR长度
    novel_tir_len_path = tmp_output_dir + '/novel.itr.fa'
    tir_len_names, tir_len_contigs = read_fasta(novel_tir_len_path)

    #获取novel TIR的拷贝数和多序列比对文件
    plant = 1
    member_script_path = '/home/hukang/HiTE/tools/make_fasta_from_blast.sh'
    subset_script_path = '/home/hukang/HiTE/tools/ready_for_MSA.sh'
    reference = '/home/hukang/EDTA/krf_test/rice/GCF_001433935.1_IRGSP-1.0_genomic.rename.fna'
    threads = 40
    flanking_len = 50
    similar_ratio = 0.2
    TE_type = 'tir'
    ref_index = 0
    iter_num = 0
    log = Logger(tmp_output_dir + '/HiTE.log', level='debug')
    debug = 1
    output = tmp_output_dir + '/confident_tir_' + str(ref_index) + '.fa'
    flank_region_align_v4(novel_tir_path, output, flanking_len, similar_ratio, reference, TE_type, tmp_output_dir, threads,
                          ref_index, log, member_script_path, subset_script_path, plant, debug, iter_num, 'cons')
    temp_dir = tmp_output_dir + '/' + TE_type + '_copies_' + str(ref_index) + '_' + str(iter_num)

    msa_dir = tmp_output_dir + '/msa'
    if not os.path.exists(msa_dir):
        os.makedirs(msa_dir)

    # 获取novel TIR的蛋白质结构组成
    protein_path = '/home/hukang/HiTE/library/RepeatPeps.lib'
    output_table = novel_tir_path + '.domain'
    domain_temp_dir = tmp_output_dir + '/domain'
    get_domain_info(novel_tir_path, protein_path, output_table, threads, domain_temp_dir)
    domain_info = {}
    with open(output_table, 'r') as f_r:
        for i, line in enumerate(f_r):
            if i <= 1:
                continue
            parts = line.split('\t')
            tir_name = parts[0]
            domain_name = parts[1]
            TE_start = parts[2]
            TE_end = parts[3]
            domain_start = parts[4]
            domain_end = parts[5]
            if not domain_info.__contains__(tir_name):
                domain_info[tir_name] = []
            info_array = domain_info[tir_name]
            info_array.append((domain_name, TE_start, TE_end, domain_start, domain_end))


    #存储所有的novel tir
    #统计新的TIR和终端有多少个
    known_terminal_count = 0
    novel_terminal_count = 0
    data = {}
    data_names = []
    data_tir_lens = []
    data_tir_types = []
    date_copy_nums = []
    data_msa_files = []
    data_domain_names = []
    data_domain_TE_starts = []
    data_domain_TE_ends = []
    data_domain_starts = []
    data_domain_ends = []
    for name in novel_tir_names:
        #获取TIR长度
        if tir_len_contigs.__contains__(name + '-lTIR'):
            tir_len = len(tir_len_contigs[name + '-lTIR'])
        else:
            continue
        #获取TIR类型
        if name in new_tir_names:
            type = 'novel_terminal'
            novel_terminal_count += 1
        else:
            type = 'known_terminal'
            known_terminal_count += 1
        #获取TIR的拷贝数
        member_file = temp_dir + '/' + name + '.blast.bed.fa'
        member_names, member_contigs = read_fasta(member_file)
        #获取TIR的MSA文件
        file_name = name + '.blast.bed.fa.maf.fa'
        file_path = temp_dir + '/' + file_name
        if not os.path.exists(file_path):
            print('mas not exist: ' + file_path)
        os.system('cp ' + file_path + ' ' + msa_dir)

        #获取domain信息
        if not domain_info.__contains__(name):
            data_names.append(name)
            data_tir_lens.append(tir_len)
            data_tir_types.append(type)
            date_copy_nums.append(len(member_names))
            data_msa_files.append(file_name)
            data_domain_names.append('')
            data_domain_TE_starts.append('')
            data_domain_TE_ends.append('')
            data_domain_starts.append('')
            data_domain_ends.append('')
        else:
            info_array = domain_info[name]
            for j, info in enumerate(info_array):
                if j == 0:
                    data_names.append(name)
                    data_tir_lens.append(tir_len)
                    data_tir_types.append(type)
                    date_copy_nums.append(len(member_names))
                    data_msa_files.append(file_name)
                    data_domain_names.append(info[0])
                    data_domain_TE_starts.append(info[1])
                    data_domain_TE_ends.append(info[2])
                    data_domain_starts.append(info[3])
                    data_domain_ends.append(info[4])
                else:
                    data_names.append('')
                    data_tir_lens.append('')
                    data_tir_types.append('')
                    date_copy_nums.append('')
                    data_msa_files.append('')
                    data_domain_names.append(info[0])
                    data_domain_TE_starts.append(info[1])
                    data_domain_TE_ends.append(info[2])
                    data_domain_starts.append(info[3])
                    data_domain_ends.append(info[4])
    data['name'] = data_names
    data['terminal tir len'] = data_tir_lens
    data['terminal type'] = data_tir_types
    data['copy num'] = date_copy_nums
    data['msa file'] = data_msa_files
    data['domain name'] = data_domain_names
    data['TE start'] = data_domain_TE_starts
    data['TE end'] = data_domain_TE_ends
    data['domain start'] = data_domain_starts
    data['domain end'] = data_domain_ends
    #print(data)
    print(novel_terminal_count, known_terminal_count)

    df = pd.DataFrame(data)

    # 将 DataFrame 存储到 Excel 文件中
    with pd.ExcelWriter(tmp_output_dir + '/data.xlsx', engine="openpyxl") as writer:
        to_excel_auto_column_weight(df, writer, f'novel TIR information')


def to_excel_auto_column_weight(df: pd.DataFrame, writer: ExcelWriter, sheet_name="Shee1"):
    """DataFrame保存为excel并自动设置列宽"""
    # 数据 to 写入器，并指定sheet名称
    df.to_excel(writer, sheet_name=sheet_name, index=False)
    #  计算每列表头的字符宽度
    column_widths = (
        df.columns.to_series().apply(lambda x: len(str(x).encode('gbk'))).values
    )
    #  计算每列的最大字符宽度
    max_widths = (
        df.astype(str).applymap(lambda x: len(str(x).encode('gbk'))).agg(max).values
    )
    # 取前两者中每列的最大宽度
    widths = np.max([column_widths, max_widths], axis=0)
    # 指定sheet，设置该sheet的每列列宽
    worksheet = writer.sheets[sheet_name]
    for i, width in enumerate(widths, 1):
        # openpyxl引擎设置字符宽度时会缩水0.5左右个字符，所以干脆+2使左右都空出一个字宽。
        worksheet.column_dimensions[get_column_letter(i)].width = width + 2


def analyze_potato_libs():
    #这个实验分析流程：
    # 1.应该可以拿到5个不同的TE库，分别是来自于RepBase，RepeatMasker，RepeatModeler2，EDTA，以及HiTE。
    # 2.有如下的情况可以分析：a) HiTE独有，其它4个库没有的TE序列。 b) HiTE与其它任意一个库共有，但是剩余的库没有的TE序列。 c）HiTE没有，但是其它库有的TE序列。
    # 3.我们统计上述三种情况的TE，然后可以拿出一些具体的样例进行分析，目的需要证明两个点：a) HiTE找到的序列 (无论独有或者是共有)，是真实的TE序列。b）HiTE没有找到的序列，不是TE或者不是全长TE（缺乏TE结构）。
    #work_dir = '/homeb/hukang/KmerRepFinder_test/library/ath_test'
    work_dir = '/home/hukang/HiTE_test/ath_test'
    lib_paths = []
    HiTE_lib = work_dir + '/HiTE_lib.fa'
    repbase_lib = work_dir + '/repbase.ref'
    #rm_lib = work_dir + '/ath_repeatmasker.ref'
    rm2_lib = work_dir + '/rm2_lib.fa'
    edta_lib = work_dir + '/EDTA_lib.fa'
    lib_paths.append((repbase_lib, 'Repbase'))
    #lib_paths.append((rm_lib, 'RM'))
    lib_paths.append((HiTE_lib, 'HiTE'))
    lib_paths.append((rm2_lib, 'RM2'))
    lib_paths.append((edta_lib, 'EDTA'))
    analyze_lib_name = 'HiTE'

    # repbase1_lib = work_dir + '/rm2_lib1.fa'
    # repbase_names, repbase_contigs = read_fasta_v1(rm2_lib)
    # store_fasta(repbase_contigs, repbase1_lib)

    # 0.对library重命名，添加工具名称，方便后续分别；合并所有的library
    merge_lib = work_dir + '/merge.fa'
    merge_lib_contigs = {}
    for i, path_item in enumerate(lib_paths):
        path = path_item[0]
        tool_name = path_item[1]
        lib_names, lib_contigs = read_fasta(path)
        for name in lib_names:
            new_name = tool_name + '-' + name
            merge_lib_contigs[new_name] = lib_contigs[name]
    store_fasta(merge_lib_contigs, merge_lib)

    merge_lib_cons = work_dir + '/merge.cons.fa'
    # 1.使用cd-hit-est聚类
    tools_dir = os.getcwd() + '/../tools'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -aL ' + str(0.8) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + merge_lib + ' -o ' + merge_lib_cons + ' -T 0 -M 0'
    #os.system(cd_hit_command)

    # 2. 分析聚类文件
    cluster_file = merge_lib_cons + '.clstr'
    cluster_idx = -1
    clusters = {}
    cluster_rep = {}
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
                name = name[0: len(name) - 3]
                cur_cluster.append(name)
                if line.endswith('*'):
                    cluster_rep['rep_' + str(cluster_idx)] = name
    print(len(clusters))
    #使用一个dict记录每个簇，包含的库的类别，format: {Rep_seq: {HiTE: 1, Repbase: 1, RM: 1, RM2: 1, EDTA: 0}}
    represent_dict = {}
    for cluster_rep_name in cluster_rep.keys():
        rep_name = cluster_rep[cluster_rep_name]
        cluster_idx = cluster_rep_name.replace('rep_', '')
        cluster = clusters[cluster_idx]
        has_HiTE = 0
        has_Repbase = 0
        has_RM = 0
        has_RM2 = 0
        has_EDTA = 0
        for name in cluster:
            if name.startswith('HiTE-'):
                has_HiTE = 1
            elif name.startswith('Repbase-'):
                has_Repbase = 1
            elif name.startswith('RM-'):
                has_RM = 1
            elif name.startswith('RM2-'):
                has_RM2 = 1
            elif name.startswith('EDTA-'):
                has_EDTA = 1
        record = {}
        record['HiTE'] = has_HiTE
        record['Repbase'] = has_Repbase
        record['RM'] = has_RM
        record['RM2'] = has_RM2
        record['EDTA'] = has_EDTA
        represent_dict[rep_name] = record

    #print(represent_dict)
    #print(len(represent_dict))

    # 统计待分析Library锁独有的，共有的，缺失的序列

    #pd.set_option('display.max_rows', None)
    #转dataframe
    df = pd.DataFrame(represent_dict).T
    print('Total clusters:')
    print(df)

    # 使用条件筛选和逻辑运算符进行统计
    filtered_df = df[(df['HiTE'] == 1) & (df['Repbase'] == 0) & (df['RM'] == 0) & (df['RM2'] == 0) & (df['EDTA'] == 0)]
    print('HiTE unique clusters:')
    print(filtered_df)

    filtered_df = df[(df['HiTE'] == 1) & ((df['Repbase'] == 1) | (df['RM'] == 1) | (df['RM2'] == 1) | (df['EDTA'] == 1))]
    print('HiTE share clusters:')
    print(filtered_df)

    filtered_df = df[(df['HiTE'] == 0) & ((df['Repbase'] == 1) | (df['RM'] == 1)) | ((df['RM2'] == 1) | (df['EDTA'] == 1))]
    print('without HiTE clusters:')
    print(filtered_df)
    rows = filtered_df.shape[0]
    print(rows)

    # # 排除掉EDTA库
    # # 使用条件筛选和逻辑运算符进行统计
    # filtered_df = df[(df['HiTE'] == 1) & (df['Repbase'] == 0) & (df['RM'] == 0) & (df['RM2'] == 0)]
    # print('HiTE unique clusters:')
    # print(filtered_df)
    #
    # filtered_df = df[
    #     (df['HiTE'] == 1) & ((df['Repbase'] == 1) | (df['RM'] == 1) | (df['RM2'] == 1))]
    # print('HiTE share clusters:')
    # print(filtered_df)
    #
    # filtered_df = df[
    #     (df['HiTE'] == 0) & ((df['Repbase'] == 1) | (df['RM'] == 1) | (df['RM2'] == 1))]
    # print('without HiTE clusters:')
    # print(filtered_df)

    # # 排除掉RepeatMasker库
    # # 使用条件筛选和逻辑运算符进行统计
    # filtered_df = df[(df['HiTE'] == 1) & (df['Repbase'] == 0)]
    # print('HiTE unique clusters:')
    # print(filtered_df)
    #
    # filtered_df = df[
    #     (df['HiTE'] == 1) & ((df['Repbase'] == 1))]
    # print('HiTE share clusters:')
    # print(filtered_df)
    #
    # filtered_df = df[
    #     (df['HiTE'] == 0) & ((df['Repbase'] == 1))]
    # print('without HiTE clusters:')
    # print(filtered_df)


def filter_last_round_tir_with_itrsearch(tmp_output_dir, reference, ref_index, TRsearch_dir, member_script_path, subset_script_path):
    i = 2
    result_type = 'cons'
    input_file = tmp_output_dir + '/confident_tir_' + str(ref_index) + '.r' + str(i-1) + '.fa'
    output_file = tmp_output_dir + '/confident_tir_' + str(ref_index) + '.r' + str(i) + '.fa'
    flank_region_align_v3(input_file, output_file, flanking_len, similar_ratio, reference, TE_type, tmp_output_dir,
                          threads, ref_index, log, member_script_path, subset_script_path, plant, debug, i, result_type)

    confident_tir_path = tmp_output_dir + '/confident_tir_' + str(ref_index) + '.r' + str(i) + '.fa'
    tir_names, tir_contigs = read_fasta(confident_tir_path)
    # 保存短tir的序列
    short_itr_contigs = get_short_tir_contigs(tir_contigs, plant)

    # 剩下的序列交由itrsearch去搜索TIR结构
    confident_tir_path = tmp_output_dir + '/confident_tir_' + str(ref_index) + '.r' + str(i) + '.no_short_tir.fa'
    for name in short_itr_contigs.keys():
        del tir_contigs[name]
    store_fasta(tir_contigs, confident_tir_path)
    all_copies_out, all_copies_log = run_itrsearch(TRsearch_dir, confident_tir_path, tmp_output_dir)
    all_copies_out_name, all_copies_out_contigs = read_fasta(all_copies_out)
    all_copies_out_contigs.update(short_itr_contigs)
    confident_tir_path = tmp_output_dir + '/confident_tir_' + str(ref_index) + '.r' + str(i) + '.all_tir.fa'
    store_fasta(all_copies_out_contigs, confident_tir_path)
    rename_fasta(confident_tir_path, confident_tir_path, 'TIR')

    confident_tir_cons = confident_tir_path + '.cons.fa'
    # 生成一致性序列
    cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + confident_tir_path + ' -o ' + confident_tir_cons + ' -T 0 -M 0'
    os.system(cd_hit_command)

    confident_tir_path = tmp_output_dir + '/confident_tir_' + str(ref_index) + '.fa'
    rename_fasta(confident_tir_cons, confident_tir_path, 'TIR')

def split_module_test():
    EAHelitron = '/public/home/hpc194701009/HiTE/bin/EAHelitron-master'
    work_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/HiTE_lib/module_test'
    # 1.生成未过滤文件
    unfilter_files = ['rice/tir_tsd_0.fa']
    for unfilter_file in unfilter_files:
        unfilter_path = work_dir + '/' + unfilter_file
        unfilter_cons = unfilter_path + '.cons'
        # 生成一致性序列
        cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
                        + ' -G 0 -g 1 -A 80 -i ' + unfilter_path + ' -o ' + unfilter_cons + ' -T 0 -M 0'
        os.system(cd_hit_command)

        longest_repeats_flanked_path = work_dir + '/rice/longest_repeats_0.flanked.fa'
        flanking_len = 50
        candidate_helitron_path = work_dir + '/rice/EAHelitron.fa'
        temp_dir = work_dir + '/rice/helitron_tmp_'+str(ref_index)
        multi_process_EAHelitron(longest_repeats_flanked_path, flanking_len, candidate_helitron_path, temp_dir, EAHelitron, threads)


def maize_centromeres():
    # 1.将玉米中心的区域抽取出来，形成maize_centromeres.fna
    reference = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Zea_mays_Mol7/Zm-Mo17-REFERENCE-CAU-2.0.fa'
    refer_centromeres = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Zea_mays_Mol7/Zm-Mo17-centromeres.fa'
    # os.system('makeblastdb -in ' + refer_centromeres + ' -dbtype nucl')
    # centromeres_regions = [('chr1', 137230000, 138850000), ('chr2', 95730000, 98660000),
    #                        ('chr3', 87930000, 89870000), ('chr4', 110310000, 112800000),
    #                        ('chr5', 102910000, 105670000), ('chr6', 78590000, 80440000),
    #                        ('chr7', 54930000, 57020000), ('chr8', 49960000, 52215000),
    #                        ('chr9', 59700000, 61890000), ('chr10', 50280000, 52320000)]
    # ref_names, ref_contigs = read_fasta(reference)
    # centromeres_contigs = {}
    # for cr in centromeres_regions:
    #     centromeres_seq = ref_contigs[cr[0]][cr[1]-1: cr[2]]
    #     name = cr[0] + ':' + str(cr[1]) + '-' + str(cr[2])
    #     centromeres_contigs[name] = centromeres_seq
    # store_fasta(centromeres_contigs, refer_centromeres)

    # 2.获取HiTE TIR和Helitron在maize_centromeres.fna上的全长拷贝（通过我们的获取拷贝的方法）
    # 3.收集在maize_centromeres.fna上具有全长拷贝的TIR和Helitron
    HiTE_lib = '/public/home/hpc194701009/KmerRepFinder_test/library/HiTE_lib/maize_mol7'
    tir_path = HiTE_lib + '/confident_tir.fa'
    tir_names, tir_contigs = read_fasta(tir_path)
    blastn2Results_path = tir_path + '.blast.out'
    repeats_path = (tir_path, refer_centromeres, blastn2Results_path)
    all_copies = multiple_alignment_blast_and_get_copies(repeats_path)
    print(all_copies)
    keep_tir_names = []
    for tir_name in all_copies.keys():
        if len(all_copies[tir_name]) > 0:
            keep_tir_names.append(tir_name)
    keep_tir_path = HiTE_lib + '/centromeres_tir.fa'
    with open(keep_tir_path, 'w') as f_save:
        for tir_name in keep_tir_names:
            tir_seq = tir_contigs[tir_name]
            f_save.write('>'+tir_name+'\n'+tir_seq+'\n')

    # 4.用cd-hit-est判断有多少新的TIR和Helitron
    # 用cd-hit-est，使用-aS 0.8 –aL 0.8 –c 0.8进行聚类，然后分析类中没有curated library出现的转座子为新的TIR。
    repbase_path = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/only_TE/repbase/maize.ref'
    repeatmaseker_path = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/maize_repeatmasker.ref'
    rm2_path = '/public/home/hpc194701009/KmerRepFinder_test/library/rm2_run_lib/maize_mol7/maize_mol7-families.fa'
    maize_lib_path = '/public/home/hpc194701009/repeat_detect_tools/EDTA-master/genome_test/maize_mol7/maizeTE02052020.fa'
    total_tir_path = HiTE_lib + '/total_tir.fa'
    total_tir_consensus = HiTE_lib + '/total_tir.cons.fa'
    os.system('cat ' + repbase_path + ' ' + repeatmaseker_path + ' ' + rm2_path + ' ' + maize_lib_path + ' ' + keep_tir_path + ' > ' + total_tir_path)

    confident_tir_names, confident_tir_contigs = read_fasta(keep_tir_path)

    tools_dir = os.getcwd() + '/../tools'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -aL ' + str(0.8) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + total_tir_path + ' -o ' + total_tir_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    cluster_file = total_tir_consensus + '.clstr'
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
                name = name[0: len(name) - 3]
                cur_cluster.append(name)
                if line.endswith('*'):
                    clusters['rep_' + str(cluster_idx)] = name
    # print(clusters)
    new_tir_contigs = {}
    exist_tir_names = set()
    for key in clusters.keys():
        if not key.startswith('rep_'):
            cluster = clusters[key]
            has_repbase = False
            for name in cluster:
                if not name.startswith('TIR_') and not name.startswith('Helitron_'):
                    has_repbase = True
                    break
            # 没有repbase序列，代表这是我们新识别的TIR
            if not has_repbase:
                represent_name = clusters['rep_' + str(key)]
                new_tir_contigs[represent_name] = confident_tir_contigs[represent_name]
            else:
                exist_tir_names.add(clusters['rep_' + str(key)])
    print('novel TIR:')
    print(len(new_tir_contigs))
    novel_tir_consensus = HiTE_lib + '/novel_tir.fa'
    store_fasta(new_tir_contigs, novel_tir_consensus)

    #将novel_tir的拷贝存成excel文件
    data = {}
    data_tir_names = []
    data_chr_name = []
    data_chr_start = []
    date_chr_end = []
    data_length = []
    for name in new_tir_contigs.keys():
        copies = all_copies[name]
        for i,copy in enumerate(copies):
            chr_info = copy[0]
            chr_parts = chr_info.split(':')
            chr_name = chr_parts[0]
            chr_start_base = int(chr_parts[1].split('-')[0])
            chr_start = chr_start_base + int(copy[1])
            chr_end = chr_start_base + int(copy[2])
            length = abs(chr_end-chr_start)
            if i==0:
                data_tir_names.append(name)
                data_chr_name.append(chr_name)
                data_chr_start.append(chr_start)
                date_chr_end.append(chr_end)
                data_length.append(length)
            else:
                data_tir_names.append('')
                data_chr_name.append(chr_name)
                data_chr_start.append(chr_start)
                date_chr_end.append(chr_end)
                data_length.append(length)
    data['TIR Name'] = data_tir_names
    data['Chr'] = data_chr_name
    data['Start (bp)'] = data_chr_start
    data['End (bp)'] = date_chr_end
    data['Size (bp)'] = data_length
    df = pd.DataFrame(data)
    # 将 DataFrame 存储到 Excel 文件中
    with pd.ExcelWriter(HiTE_lib + '/centromeres_tir.xlsx', engine="openpyxl") as writer:
        to_excel_auto_column_weight(df, writer, f'Centromeres novel TIRs')


    helitron_path = HiTE_lib + '/confident_helitron.fa'
    helitron_names, helitron_contigs = read_fasta(helitron_path)
    blastn2Results_path = helitron_path + '.blast.out'
    repeats_path = (helitron_path, refer_centromeres, blastn2Results_path)
    all_copies = multiple_alignment_blast_and_get_copies(repeats_path)
    print(all_copies)
    keep_helitron_names = []
    for helitron_name in all_copies.keys():
        if len(all_copies[helitron_name]) > 0:
            keep_helitron_names.append(helitron_name)
    keep_helitron_path = HiTE_lib + '/centromeres_helitron.fa'
    with open(keep_helitron_path, 'w') as f_save:
        for helitron_name in keep_helitron_names:
            helitron_seq = helitron_contigs[helitron_name]
            f_save.write('>' + helitron_name + '\n' + helitron_seq + '\n')
    # 用cd-hit-est，使用-aS 0.8 –aL 0.8 –c 0.8进行聚类，然后分析类中没有curated library出现的转座子为新的TIR。
    total_helitron_path = HiTE_lib + '/total_helitron.fa'
    total_helitron_consensus = HiTE_lib + '/total_helitron.cons.fa'
    os.system('cat ' + repbase_path + ' ' + repeatmaseker_path + ' ' + rm2_path + ' ' + maize_lib_path + ' ' + keep_helitron_path + ' > ' + total_helitron_path)

    confident_helitron_names, confident_helitron_contigs = read_fasta(keep_helitron_path)

    tools_dir = os.getcwd() + '/../tools'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -aL ' + str(0.8) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + total_helitron_path + ' -o ' + total_helitron_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    cluster_file = total_helitron_consensus + '.clstr'
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
                name = name[0: len(name) - 3]
                cur_cluster.append(name)
                if line.endswith('*'):
                    clusters['rep_' + str(cluster_idx)] = name
    # print(clusters)
    new_helitron_contigs = {}
    exist_helitron_names = set()
    for key in clusters.keys():
        if not key.startswith('rep_'):
            cluster = clusters[key]
            has_repbase = False
            for name in cluster:
                if not name.startswith('TIR_') and not name.startswith('Helitron_'):
                    has_repbase = True
                    break
            # 没有repbase序列，代表这是我们新识别的TIR
            if not has_repbase:
                represent_name = clusters['rep_' + str(key)]
                new_helitron_contigs[represent_name] = confident_helitron_contigs[represent_name]
            else:
                exist_helitron_names.add(clusters['rep_' + str(key)])
    print('novel Helitron:')
    print(len(new_helitron_contigs))
    novel_helitron_consensus = HiTE_lib + '/novel_helitron.fa'
    store_fasta(new_helitron_contigs, novel_helitron_consensus)

    # 将novel_helitron的拷贝存成excel文件
    data = {}
    data_helitron_names = []
    data_chr_name = []
    data_chr_start = []
    date_chr_end = []
    data_length = []
    for name in new_helitron_contigs.keys():
        copies = all_copies[name]
        for i, copy in enumerate(copies):
            chr_info = copy[0]
            chr_parts = chr_info.split(':')
            chr_name = chr_parts[0]
            chr_start_base = int(chr_parts[1].split('-')[0])
            chr_start = chr_start_base + int(copy[1])
            chr_end = chr_start_base + int(copy[2])
            length = abs(chr_end - chr_start)
            if i == 0:
                data_helitron_names.append(name)
                data_chr_name.append(chr_name)
                data_chr_start.append(chr_start)
                date_chr_end.append(chr_end)
                data_length.append(length)
            else:
                data_helitron_names.append('')
                data_chr_name.append(chr_name)
                data_chr_start.append(chr_start)
                date_chr_end.append(chr_end)
                data_length.append(length)
    data['Helitron Name'] = data_helitron_names
    data['Chr'] = data_chr_name
    data['Start (bp)'] = data_chr_start
    data['End (bp)'] = date_chr_end
    data['Size (bp)'] = data_length
    df = pd.DataFrame(data)
    # 将 DataFrame 存储到 Excel 文件中
    with pd.ExcelWriter(HiTE_lib + '/centromeres_helitron.xlsx', engine="openpyxl") as writer:
        to_excel_auto_column_weight(df, writer, f'Centromeres novel Helitrons')


def get_full_length_TE(tir_path, refer_centromeres):
    tir_names, tir_contigs = read_fasta(tir_path)
    blastn2Results_path = tir_path + '.blast.out'
    repeats_path = (tir_path, refer_centromeres, blastn2Results_path)
    all_copies = multiple_alignment_blast_and_get_copies(repeats_path)
    print(all_copies)
    keep_tir_names = []
    for tir_name in all_copies.keys():
        if len(all_copies[tir_name]) > 0:
            keep_tir_names.append(tir_name)
    keep_tir_path = tir_path + '_centromeres_tir.fa'
    with open(keep_tir_path, 'w') as f_save:
        for tir_name in keep_tir_names:
            tir_seq = tir_contigs[tir_name]
            f_save.write('>' + tir_name + '\n' + tir_seq + '\n')
    return keep_tir_path, all_copies


def get_overlap_library(HiTE_keep_tir_path, maize_lib_keep_tir_path, total_tir_path):
    total_tir_consensus = total_tir_path + '.cons.fa'
    os.system('cat ' + HiTE_keep_tir_path + ' ' + maize_lib_keep_tir_path + ' > ' + total_tir_path)

    HiTE_tir_names, HiTE_tir_contigs = read_fasta(HiTE_keep_tir_path)
    maize_lib_tir_names, maize_lib_contigs = read_fasta(maize_lib_keep_tir_path)

    tools_dir = os.getcwd() + '/../tools'
    cd_hit_command = tools_dir + '/cd-hit-est -d 0 -aS ' + str(0.8) + ' -aL ' + str(0.8) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + total_tir_path + ' -o ' + total_tir_consensus + ' -T 0 -M 0'
    os.system(cd_hit_command)

    cluster_file = total_tir_consensus + '.clstr'
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
                name = name[0: len(name) - 3]
                cur_cluster.append(name)
                if line.endswith('*'):
                    clusters['rep_' + str(cluster_idx)] = name
    # print(clusters)
    new_tir_contigs = {}
    share_tir_contigs = {}
    repbase_tir_contigs = {}
    for key in clusters.keys():
        if not key.startswith('rep_'):
            cluster = clusters[key]
            has_repbase = False
            has_HiTE = False
            for name in cluster:
                if not name.startswith('TIR_') and not name.startswith('Helitron_'):
                    has_repbase = True
                else:
                    has_HiTE = True
            represent_name = clusters['rep_' + str(key)]
            # 没有repbase序列，代表这是我们新识别的TIR
            if not has_repbase:
                new_tir_contigs[represent_name] = HiTE_tir_contigs[represent_name]
            elif has_HiTE:
                if HiTE_tir_contigs.__contains__(represent_name):
                    seq = HiTE_tir_contigs[represent_name]
                elif maize_lib_contigs.__contains__(represent_name):
                    seq = maize_lib_contigs[represent_name]
                share_tir_contigs[represent_name] = seq
            else:
                repbase_tir_contigs[represent_name] = maize_lib_contigs[represent_name]
    print('HiTE unique:' + str(len(new_tir_contigs)) + ', Library unique:' + str(len(repbase_tir_contigs)) + ', share:' + str(len(share_tir_contigs)))
    return new_tir_contigs, repbase_tir_contigs, share_tir_contigs


def maize_centromeres_v1():
    #我们分析在中心粒区域，玉米注释库与HiTE TIR独有，共享的数量分别有多少

    # 1.将玉米中心的区域抽取出来，形成maize_centromeres.fna
    reference = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Zea_mays_Mol7/Zm-Mo17-REFERENCE-CAU-2.0.fa'
    refer_centromeres = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Zea_mays_Mol7/Zm-Mo17-centromeres.fa'
    # os.system('makeblastdb -in ' + refer_centromeres + ' -dbtype nucl')
    # centromeres_regions = [('chr1', 137230000, 138850000), ('chr2', 95730000, 98660000),
    #                        ('chr3', 87930000, 89870000), ('chr4', 110310000, 112800000),
    #                        ('chr5', 102910000, 105670000), ('chr6', 78590000, 80440000),
    #                        ('chr7', 54930000, 57020000), ('chr8', 49960000, 52215000),
    #                        ('chr9', 59700000, 61890000), ('chr10', 50280000, 52320000)]
    # ref_names, ref_contigs = read_fasta(reference)
    # centromeres_contigs = {}
    # for cr in centromeres_regions:
    #     centromeres_seq = ref_contigs[cr[0]][cr[1]-1: cr[2]]
    #     name = cr[0] + ':' + str(cr[1]) + '-' + str(cr[2])
    #     centromeres_contigs[name] = centromeres_seq
    # store_fasta(centromeres_contigs, refer_centromeres)


    #我们将maize_lib中的不同类型的转座子区分开
    maize_lib_path = '/public/home/hpc194701009/repeat_detect_tools/EDTA-master/genome_test/maize_mol7/maizeTE02052020.fa'
    maize_lib_dir = '/public/home/hpc194701009/repeat_detect_tools/EDTA-master/genome_test/maize_mol7/maize_lib'
    if not os.path.exists(maize_lib_dir):
        os.makedirs(maize_lib_dir)
    ltr_path = maize_lib_dir + '/ltr.fa'
    tir_path = maize_lib_dir + '/tir.fa'
    helitron_path = maize_lib_dir + '/helitron.fa'
    non_ltr_path = maize_lib_dir + '/non_ltr.fa'
    satellite_path = maize_lib_dir + '/satellite_path.fa'
    # contigNames, contigs = read_fasta(maize_lib_path)
    # class_set = set()
    # ltr_contigs = {}
    # tir_contigs = {}
    # helitron_contigs = {}
    # non_ltr_contigs = {}
    # satellite_contigs = {}
    # for name in contigNames:
    #     class_name = name.split('#')[1]
    #     if class_name.__contains__('LTR'):
    #         ltr_contigs[name] = contigs[name]
    #     elif class_name.__contains__('Helitron'):
    #         helitron_contigs[name] = contigs[name]
    #     elif class_name.__contains__('MITE') or class_name.__contains__('DNA'):
    #         tir_contigs[name] = contigs[name]
    #     elif class_name.__contains__('LINE'):
    #         non_ltr_contigs[name] = contigs[name]
    #     else:
    #         satellite_contigs[name] = contigs[name]
    #     class_set.add(class_name)
    # print(class_set)
    # #根据class_name将注释库进行分类
    # store_fasta(ltr_contigs, ltr_path)
    # store_fasta(tir_contigs, tir_path)
    # store_fasta(helitron_contigs, helitron_path)
    # store_fasta(non_ltr_contigs, non_ltr_path)
    # store_fasta(satellite_contigs, satellite_path)


    # 2.获取HiTE TIR和Helitron在maize_centromeres.fna上的全长拷贝（通过我们的获取拷贝的方法）
    # 3.收集在maize_centromeres.fna上具有全长拷贝的TIR和Helitron
    HiTE_lib = '/public/home/hpc194701009/KmerRepFinder_test/library/HiTE_lib/maize_mol7'
    HiTE_tir_path = HiTE_lib + '/confident_tir.fa'
    HiTE_keep_tir_path, HiTE_all_copies = get_full_length_TE(HiTE_tir_path, refer_centromeres)

    maize_lib_keep_tir_path, maize_lib_all_copies = get_full_length_TE(tir_path, refer_centromeres)

    # 4.用cd-hit-est判断HiTE, maize lib分别独有，共享的有多少
    # 用cd-hit-est，使用-aS 0.8 –aL 0.8 –c 0.8进行聚类，然后分析类中没有curated library出现的转座子为新的TIR。
    total_tir_path = HiTE_lib + '/total_tir.fa'
    new_tir_contigs, repbase_tir_contigs, share_tir_contigs = get_overlap_library(HiTE_keep_tir_path, maize_lib_keep_tir_path, total_tir_path)

    HiTE_novel_tir_consensus = HiTE_lib + '/HiTE_novel_tir.fa'
    store_fasta(new_tir_contigs, HiTE_novel_tir_consensus)
    Lib_novel_tir_consensus = HiTE_lib + '/library_novel_tir.fa'
    store_fasta(repbase_tir_contigs, Lib_novel_tir_consensus)
    share_tir_consensus = HiTE_lib + '/share_tir.fa'
    store_fasta(share_tir_contigs, share_tir_consensus)

    # 将novel_tir的拷贝存成excel文件
    data = {}
    data_tir_names = []
    data_chr_name = []
    data_chr_start = []
    date_chr_end = []
    data_length = []
    for name in new_tir_contigs.keys():
        copies = HiTE_all_copies[name]
        for i, copy in enumerate(copies):
            chr_info = copy[0]
            chr_parts = chr_info.split(':')
            chr_name = chr_parts[0]
            chr_start_base = int(chr_parts[1].split('-')[0])
            chr_start = chr_start_base + int(copy[1])
            chr_end = chr_start_base + int(copy[2])
            length = abs(chr_end - chr_start)
            if i == 0:
                data_tir_names.append(name)
                data_chr_name.append(chr_name)
                data_chr_start.append(chr_start)
                date_chr_end.append(chr_end)
                data_length.append(length)
            else:
                data_tir_names.append('')
                data_chr_name.append(chr_name)
                data_chr_start.append(chr_start)
                date_chr_end.append(chr_end)
                data_length.append(length)
    data['TIR Name'] = data_tir_names
    data['Chr'] = data_chr_name
    data['Start (bp)'] = data_chr_start
    data['End (bp)'] = date_chr_end
    data['Size (bp)'] = data_length
    df = pd.DataFrame(data)
    # 将 DataFrame 存储到 Excel 文件中
    with pd.ExcelWriter(HiTE_lib + '/centromeres_tir.xlsx', engine="openpyxl") as writer:
        to_excel_auto_column_weight(df, writer, f'Centromeres novel TIRs')


    HiTE_helitron_path = HiTE_lib + '/confident_helitron.fa'

    HiTE_keep_helitron_path, HiTE_all_copies = get_full_length_TE(HiTE_helitron_path, refer_centromeres)

    maize_lib_keep_helitron_path, maize_lib_all_copies = get_full_length_TE(helitron_path, refer_centromeres)

    total_helitron_path = HiTE_lib + '/total_helitron.fa'
    # 4.用cd-hit-est判断HiTE, maize lib分别独有，共享的有多少
    # 用cd-hit-est，使用-aS 0.8 –aL 0.8 –c 0.8进行聚类，然后分析类中没有curated library出现的转座子为新的TIR。
    new_helitron_contigs, repbase_helitron_contigs, share_helitron_contigs = get_overlap_library(HiTE_keep_helitron_path,
                                                                                                 maize_lib_keep_helitron_path,
                                                                                                 total_helitron_path)
    HiTE_novel_helitron_consensus = HiTE_lib + '/HiTE_novel_helitron.fa'
    store_fasta(new_helitron_contigs, HiTE_novel_helitron_consensus)
    Lib_novel_helitron_consensus = HiTE_lib + '/library_novel_helitron.fa'
    store_fasta(repbase_helitron_contigs, Lib_novel_helitron_consensus)
    share_helitron_consensus = HiTE_lib + '/share_helitron.fa'
    store_fasta(share_helitron_contigs, share_helitron_consensus)

    # 将novel_helitron的拷贝存成excel文件
    data = {}
    data_helitron_names = []
    data_chr_name = []
    data_chr_start = []
    date_chr_end = []
    data_length = []
    for name in new_helitron_contigs.keys():
        copies = HiTE_all_copies[name]
        for i, copy in enumerate(copies):
            chr_info = copy[0]
            chr_parts = chr_info.split(':')
            chr_name = chr_parts[0]
            chr_start_base = int(chr_parts[1].split('-')[0])
            chr_start = chr_start_base + int(copy[1])
            chr_end = chr_start_base + int(copy[2])
            length = abs(chr_end - chr_start)
            if i == 0:
                data_helitron_names.append(name)
                data_chr_name.append(chr_name)
                data_chr_start.append(chr_start)
                date_chr_end.append(chr_end)
                data_length.append(length)
            else:
                data_helitron_names.append('')
                data_chr_name.append(chr_name)
                data_chr_start.append(chr_start)
                date_chr_end.append(chr_end)
                data_length.append(length)
    data['Helitron Name'] = data_helitron_names
    data['Chr'] = data_chr_name
    data['Start (bp)'] = data_chr_start
    data['End (bp)'] = date_chr_end
    data['Size (bp)'] = data_length
    df = pd.DataFrame(data)
    # 将 DataFrame 存储到 Excel 文件中
    with pd.ExcelWriter(HiTE_lib + '/centromeres_helitron.xlsx', engine="openpyxl") as writer:
        to_excel_auto_column_weight(df, writer, f'Centromeres novel Helitrons')


def maize_Mol7():
    #我们分析整个基因组，玉米注释库与HiTE TIR独有，共享的数量分别有多少
    # 1.将玉米中心的区域抽取出来，形成maize_centromeres.fna
    reference = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Zea_mays_Mol7/Zm-Mo17-REFERENCE-CAU-2.0.fa'
    blast_db_command = 'makeblastdb -dbtype nucl -in ' + reference + ' > /dev/null 2>&1'
    os.system(blast_db_command)

    #我们将maize_lib中的不同类型的转座子区分开
    maize_lib_dir = '/public/home/hpc194701009/repeat_detect_tools/EDTA-master/genome_test/maize_mol7/maize_lib'
    tir_path = maize_lib_dir + '/tir.fa'
    helitron_path = maize_lib_dir + '/helitron.fa'

    # 2.获取HiTE TIR和Helitron在maize_centromeres.fna上的全长拷贝（通过我们的获取拷贝的方法）
    # 3.收集在maize_centromeres.fna上具有全长拷贝的TIR和Helitron
    HiTE_lib = '/public/home/hpc194701009/KmerRepFinder_test/library/HiTE_lib/maize_mol7'
    HiTE_tir_path = HiTE_lib + '/confident_tir.fa'
    HiTE_keep_tir_path, HiTE_all_copies = get_full_length_TE(HiTE_tir_path, reference)

    maize_lib_keep_tir_path, maize_lib_all_copies = get_full_length_TE(tir_path, reference)

    # 4.用cd-hit-est判断HiTE, maize lib分别独有，共享的有多少
    # 用cd-hit-est，使用-aS 0.8 –aL 0.8 –c 0.8进行聚类，然后分析类中没有curated library出现的转座子为新的TIR。
    total_tir_path = HiTE_lib + '/total_tir.fa'
    new_tir_contigs, repbase_tir_contigs, share_tir_contigs = get_overlap_library(HiTE_keep_tir_path, maize_lib_keep_tir_path, total_tir_path)

    HiTE_novel_tir_consensus = HiTE_lib + '/HiTE_novel_tir.fa'
    store_fasta(new_tir_contigs, HiTE_novel_tir_consensus)
    Lib_novel_tir_consensus = HiTE_lib + '/library_novel_tir.fa'
    store_fasta(repbase_tir_contigs, Lib_novel_tir_consensus)
    share_tir_consensus = HiTE_lib + '/share_tir.fa'
    store_fasta(share_tir_contigs, share_tir_consensus)

    # 将novel_tir的拷贝存成excel文件
    data = {}
    data_tir_names = []
    data_chr_name = []
    data_chr_start = []
    date_chr_end = []
    data_length = []
    for name in new_tir_contigs.keys():
        copies = HiTE_all_copies[name]
        for i, copy in enumerate(copies):
            chr_name = copy[0]
            chr_start = int(copy[1])
            chr_end = int(copy[2])
            length = abs(chr_end - chr_start)
            if i == 0:
                data_tir_names.append(name)
                data_chr_name.append(chr_name)
                data_chr_start.append(chr_start)
                date_chr_end.append(chr_end)
                data_length.append(length)
            else:
                data_tir_names.append('')
                data_chr_name.append(chr_name)
                data_chr_start.append(chr_start)
                date_chr_end.append(chr_end)
                data_length.append(length)
    data['TIR Name'] = data_tir_names
    data['Chr'] = data_chr_name
    data['Start (bp)'] = data_chr_start
    data['End (bp)'] = date_chr_end
    data['Size (bp)'] = data_length
    df = pd.DataFrame(data)
    # 将 DataFrame 存储到 Excel 文件中
    with pd.ExcelWriter(HiTE_lib + '/genome_tir.xlsx', engine="openpyxl") as writer:
        to_excel_auto_column_weight(df, writer, f'Genome novel TIRs')


    HiTE_helitron_path = HiTE_lib + '/confident_helitron.fa'

    HiTE_keep_helitron_path, HiTE_all_copies = get_full_length_TE(HiTE_helitron_path, reference)

    maize_lib_keep_helitron_path, maize_lib_all_copies = get_full_length_TE(helitron_path, reference)

    total_helitron_path = HiTE_lib + '/total_helitron.fa'
    # 4.用cd-hit-est判断HiTE, maize lib分别独有，共享的有多少
    # 用cd-hit-est，使用-aS 0.8 –aL 0.8 –c 0.8进行聚类，然后分析类中没有curated library出现的转座子为新的TIR。
    new_helitron_contigs, repbase_helitron_contigs, share_helitron_contigs = get_overlap_library(HiTE_keep_helitron_path,
                                                                                                 maize_lib_keep_helitron_path,
                                                                                                 total_helitron_path)
    HiTE_novel_helitron_consensus = HiTE_lib + '/HiTE_novel_helitron.fa'
    store_fasta(new_helitron_contigs, HiTE_novel_helitron_consensus)
    Lib_novel_helitron_consensus = HiTE_lib + '/library_novel_helitron.fa'
    store_fasta(repbase_helitron_contigs, Lib_novel_helitron_consensus)
    share_helitron_consensus = HiTE_lib + '/share_helitron.fa'
    store_fasta(share_helitron_contigs, share_helitron_consensus)

    # 将novel_helitron的拷贝存成excel文件
    data = {}
    data_helitron_names = []
    data_chr_name = []
    data_chr_start = []
    date_chr_end = []
    data_length = []
    for name in new_helitron_contigs.keys():
        copies = HiTE_all_copies[name]
        for i, copy in enumerate(copies):
            chr_name = copy[0]
            chr_start = int(copy[1])
            chr_end = int(copy[2])
            length = abs(chr_end - chr_start)
            if i == 0:
                data_helitron_names.append(name)
                data_chr_name.append(chr_name)
                data_chr_start.append(chr_start)
                date_chr_end.append(chr_end)
                data_length.append(length)
            else:
                data_helitron_names.append('')
                data_chr_name.append(chr_name)
                data_chr_start.append(chr_start)
                date_chr_end.append(chr_end)
                data_length.append(length)
    data['Helitron Name'] = data_helitron_names
    data['Chr'] = data_chr_name
    data['Start (bp)'] = data_chr_start
    data['End (bp)'] = date_chr_end
    data['Size (bp)'] = data_length
    df = pd.DataFrame(data)
    # 将 DataFrame 存储到 Excel 文件中
    with pd.ExcelWriter(HiTE_lib + '/genome_helitron.xlsx', engine="openpyxl") as writer:
        to_excel_auto_column_weight(df, writer, f'Genome novel Helitrons')


def maize_length_compare():
    HiTE_lib = '/public/home/hpc194701009/KmerRepFinder_test/library/HiTE_lib/maize_mol7/genome'
    #HiTE_lib = '/public/home/hpc194701009/KmerRepFinder_test/library/HiTE_lib/maize_mol7/centromeres'
    # paths = [HiTE_lib + '/HiTE_novel_tir.fa', HiTE_lib + '/library_novel_tir.fa']
    # labels = ['HiTE-TIR', 'Maize Library']

    # paths = [HiTE_lib + '/HiTE_novel_tir.internal.fa', HiTE_lib + '/library_novel_tir.internal.fa']
    # labels = ['HiTE-TIR-internal', 'Maize Library-internal']

    HiTE_lib = '/public/home/hpc194701009/KmerRepFinder_test/library/tir_test'
    paths = [HiTE_lib + '/real_tir.fa', HiTE_lib + '/real_helitron.fa']
    labels = ['Maize Library TIR', 'Maize Library Helitron']

    my_pal = {"Maize Library TIR": "#4497B1", "Maize Library Helitron": "#F7B92E"}
    output_path = HiTE_lib + '/tir_length_dist.txt'
    output_fig = HiTE_lib + '/tir_length_dist.png'

    # paths = [HiTE_lib + '/HiTE_novel_helitron.fa', HiTE_lib + '/library_novel_helitron.fa']
    # labels = ['HiTE-Helitron', 'Maize Library']
    #
    # my_pal = {"HiTE-TIR": "#4497B1", "Maize Library": "#F7B92E"}
    # output_path = HiTE_lib + '/helitron_length_dist.txt'
    # output_fig = HiTE_lib + '/helitron_length_dist.png'
    get_length_dist_boxplot(paths, labels, my_pal, output_path, output_fig)


def get_TIR_internal():
    # 1. 获取itr文件，解析itr长度，获取Internal序列，并存成文件
    HiTE_lib = '/public/home/hpc194701009/KmerRepFinder_test/library/HiTE_lib/maize_mol7/genome'
    paths = [HiTE_lib + '/HiTE_novel_tir.fa.itr', HiTE_lib + '/library_novel_tir.fa.itr']
    internal_paths = [HiTE_lib + '/HiTE_novel_tir.internal.fa', HiTE_lib + '/library_novel_tir.internal.fa']

    for i, path in enumerate(paths):
        internal_contigs = {}
        contigNames, contigs = read_fasta_v1(path)
        for name in contigNames:
            itr_len = int(name.split('Length itr=')[1])
            seq = contigs[name]
            internal_seq = seq[itr_len: len(seq)-itr_len]
            internal_contigs[name.split(' ')[0]] = internal_seq
        store_fasta(internal_contigs, internal_paths[i])

def merge_overlapping_intervals(intervals):
    # 按照起始位置对基因组位置数组进行排序
    sorted_intervals = sorted(intervals, key=lambda x: x[0])

    merged_intervals = []
    current_interval = sorted_intervals[0]

    # 遍历排序后的基因组位置数组
    for interval in sorted_intervals[1:]:
        # 如果当前基因组位置与下一个基因组位置有重叠
        if current_interval[1] >= interval[0]:
            # 更新当前基因组位置的结束位置
            current_interval = (current_interval[0], max(current_interval[1], interval[1]))
        else:
            # 如果没有重叠，则将当前基因组位置添加到结果数组中，并更新当前基因组位置为下一个基因组位置
            merged_intervals.append(current_interval)
            current_interval = interval

    # 将最后一个基因组位置添加到结果数组中
    merged_intervals.append(current_interval)

    return merged_intervals

def analyze_EDTA_false_positives():
    # 1.将EDTA.out转换成可分析格式，读取文件EDTA.gff
    # 2、获取每个element对应的染色体以及染色体上的起始、终止位置，位置可能会有overlap，合并有overlap的坐标，此时得到每个element对应的不重叠的位置
    # 3.对所有分类为DNA且不是DNA/Helitron的元素的坐标进行合并，此时得到所有DNA转座子的不重叠的位置
    work_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/EDTA_lib/potato_reference'
    edta_gff_path = work_dir + '/EDTA.gff'
    # te_chrs -> {te_name: {chr_name: []}}
    te_chrs = {}
    #class_chrs = -> {class_name: {chr_name: []}}
    class_chrs = {}
    # main_class_chrs = -> {class_name: {chr_name: []}}
    main_class_chrs = {}

    te_class_dict = {}

    #保存每个te对应的碱基数量
    te_bases = {}
    #保存每个 class 对应的碱基数量
    class_bases = {}
    # 保存每个 main_class 对应的碱基数量
    main_class_bases = {}
    with open(edta_gff_path, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            chr_name = parts[0]
            chr_start = int(parts[3])
            chr_end = int(parts[4])
            info = parts[8]
            info_parts = info.split(';')
            te_name = info_parts[1].split('=')[1]
            te_class = info_parts[2].split('=')[1]
            if not te_chrs.__contains__(te_name):
                te_chrs[te_name] = {}
            chr_dict = te_chrs[te_name]
            if not chr_dict.__contains__(chr_name):
                chr_dict[chr_name] = []
            chr_pos = chr_dict[chr_name]
            chr_pos.append((chr_start, chr_end))

            if not class_chrs.__contains__(te_class):
                class_chrs[te_class] = {}
            chr_dict = class_chrs[te_class]
            if not chr_dict.__contains__(chr_name):
                chr_dict[chr_name] = []
            chr_pos = chr_dict[chr_name]
            chr_pos.append((chr_start, chr_end))

            if te_class.__contains__('LTR'):
                main_class = 'LTR'
            elif te_class.__contains__('Helitron'):
                main_class = 'Helitron'
            elif te_class.__contains__('DNA'):
                main_class = 'DNA'
            elif te_class.__contains__('MITE'):
                main_class = 'MITE'

            if not main_class_chrs.__contains__(main_class):
                main_class_chrs[main_class] = {}
            chr_dict = main_class_chrs[main_class]
            if not chr_dict.__contains__(chr_name):
                chr_dict[chr_name] = []
            chr_pos = chr_dict[chr_name]
            chr_pos.append((chr_start, chr_end))

            te_class_dict[te_name] = te_class
            #print(chr_name, chr_start, chr_end, te_name, te_class)
            # if len(te_chrs) > 50:
            #     break
    for te_name in te_chrs.keys():
        chr_dict = te_chrs[te_name]
        total_bases = 0
        for chr_name in chr_dict.keys():
            chr_pos = chr_dict[chr_name]
            merged_intervals = merge_overlapping_intervals(chr_pos)
            for segment in merged_intervals:
                total_bases += abs(segment[1] - segment[0])
        te_bases[te_name] = total_bases
    #print(te_bases)

    #统计占比最高的前10个DNA转座子
    te_base_list = list(te_bases.items())
    te_base_list.sort(key=lambda x: -x[1])
    top_dna = []
    for item in te_base_list:
        te_name = item[0]
        class_name = te_class_dict[te_name]
        if not class_name.__contains__('Helitron') and (class_name.__contains__('DNA') or class_name.__contains__('MITE')):
            top_dna.append((te_name, class_name, item[1]))

    for class_name in class_chrs.keys():
        chr_dict = class_chrs[class_name]
        total_bases = 0
        for chr_name in chr_dict.keys():
            chr_pos = chr_dict[chr_name]
            merged_intervals = merge_overlapping_intervals(chr_pos)
            for segment in merged_intervals:
                total_bases += abs(segment[1] - segment[0])
        class_bases[class_name] = total_bases
    #print(class_bases)

    for class_name in main_class_chrs.keys():
        chr_dict = main_class_chrs[class_name]
        total_bases = 0
        for chr_name in chr_dict.keys():
            chr_pos = chr_dict[chr_name]
            merged_intervals = merge_overlapping_intervals(chr_pos)
            for segment in merged_intervals:
                total_bases += abs(segment[1] - segment[0])
        main_class_bases[class_name] = total_bases

    # store for testing
    te_base_file = work_dir + '/te_base.json'
    with codecs.open(te_base_file, 'w', encoding='utf-8') as f:
        json.dump(te_bases, f)
    class_base_file = work_dir + '/class_base.json'
    with codecs.open(class_base_file, 'w', encoding='utf-8') as f:
        json.dump(class_bases, f)
    main_class_base_file = work_dir + '/main_class_base.json'
    with codecs.open(main_class_base_file, 'w', encoding='utf-8') as f:
        json.dump(main_class_bases, f)
    top_dna_file = work_dir + '/top_dna.json'
    with codecs.open(top_dna_file, 'w', encoding='utf-8') as f:
        json.dump(top_dna, f)

def estimate_insert_time(identity, miu):
    # 估计序列的插入时间
    d = 1 - identity
    K=  -3/4*math.log(1-d*4/3)
    T = K/(2*miu)
    return T

def analyz_maize_TIR_insert_time():
    work_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/HiTE_lib/maize_mol7'
    tir_path = work_dir + '/confident_tir.fa'
    log_path = tir_path + '.log'
    identity_list = []
    with open(log_path, 'r') as f_r:
        for line in f_r:
            if line.__contains__('Identity percentage'):
                identity = float(line.split(':')[1].strip())
                identity_list.append(identity)

    zero_count = 0
    miu = 3.3e-8
    insert_time_list = []
    for identity in identity_list:
        T = estimate_insert_time(identity, miu)
        if T == 0:
            zero_count += 1
        insert_time_list.append(T/1000000)
    print(zero_count)
    output_path = work_dir + '/tir_insert_time.txt'
    output_fig = work_dir + '/tir_insert_time.png'
    with open(output_path, 'w') as f_save:
        f_save.write('insert_time\tmethods\n')
        for insert_time in insert_time_list:
            f_save.write(str(insert_time)+'\t'+'HiTE-TIR'+'\n')
    get_insert_time_dist_boxplot(output_path, output_fig)

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


def find_longest_tandem_repeat_tail(sequence, tail_length=20, min_repeats=2, min_unit_length=2, max_unit_length=6):
    """
    查找序列末尾最长的串联重复及其终止索引。
    :param sequence: 输入的序列。
    :param tail_length: 检测序列末尾的长度。
    :param min_repeats: 最小重复次数。
    :param min_unit_length: 最小重复单元长度。
    :param max_unit_length: 最大重复单元长度。
    :return: 如果找到串联重复，返回 (重复单元, 终止索引)；否则返回 (None, -1)。
    """
    # 截取序列末尾的指定长度
    tail_seq = sequence[-tail_length:] if len(sequence) >= tail_length else sequence

    longest_unit = None
    longest_end_index = -1
    max_repeat_length = 0  # 用于记录最大的 motif 长度 × 重复次数

    # 遍历可能的重复单元长度
    for unit_length in range(min_unit_length, max_unit_length + 1):
        # 从末尾开始，检查是否有足够的重复单元
        for start in range(len(tail_seq) - unit_length * min_repeats, -1, -1):
            # 提取候选重复单元
            unit = tail_seq[start:start + unit_length]
            # 检查是否满足最小重复次数
            repeat_count = 1
            for i in range(1, (len(tail_seq) - start) // unit_length):
                if tail_seq[start + i * unit_length:start + (i + 1) * unit_length] != unit:
                    break
                repeat_count += 1
            # 如果满足最小重复次数，计算总长度
            if repeat_count >= min_repeats:
                total_length = unit_length * repeat_count
                # 更新最长的串联重复
                if total_length > max_repeat_length:
                    max_repeat_length = total_length
                    longest_unit = unit
                    longest_end_index = len(sequence) - len(tail_seq) + start + total_length

    return longest_unit, longest_end_index

def draw_intact_LTR_insert_time(intact_ltr_paths, output_pdf):
    # 初始化列表，用于存储每个 genome_name 的插入时间数据
    all_data = []

    for genome_name, file in intact_ltr_paths:

        # 读取 CSV 文件并跳过 # 注释行
        df = pd.read_csv(file, sep='\t', comment='#',
                         names=['LTR_loc', 'Status', 'Motif', 'TSD', 'UN1', 'UN2', 'Internal', 'Identity', 'UN3', 'UN4', 'UN5', 'Classification', 'Insertion_Time'])

        df['Insertion_Time'] = pd.to_numeric(df['Insertion_Time'], errors='coerce')

        # 计算插入时间并转换为百万年，忽略 NaN
        df['Insertion_Time'] = df['Insertion_Time'] / 1_000_000

        # 为每个 genome_name 添加 Copia 和 Gypsy 的插入时间及分类信息
        for classification in ['LTR/Copia', 'LTR/Gypsy']:
            filtered_df = df[df['Classification'] == classification]
            for insertion_time in filtered_df['Insertion_Time']:
                all_data.append([genome_name, insertion_time, classification])

    # 将数据转换为 DataFrame
    all_data_df = pd.DataFrame(all_data, columns=['Genome', 'Insertion_Time', 'Classification'])

    # 定义时间区间（根据实际数据调整）
    bins = np.linspace(0, 10, 41)  # 将插入时间划分为10个区间
    bin_labels = [f'{bins[i]}-{bins[i + 1]}' for i in range(len(bins) - 1)]

    # 创建一个DataFrame来存储每个基因组在各时间区间内的LTR插入频率
    genome_ids = all_data_df['Genome'].unique()  # 获取基因组ID
    heatmap_df = pd.DataFrame(0, index=genome_ids, columns=bin_labels)  # 初始化为0

    # 按基因组统计各时间区间的插入次数
    for genome in genome_ids:
        genome_data = all_data_df[all_data_df['Genome'] == genome]['Insertion_Time']
        counts, _ = np.histogram(genome_data, bins=bins)  # 计算每个区间的插入次数
        heatmap_df.loc[genome] = counts  # 将结果保存到对应基因组的行

    # 绘制热图
    plt.figure(figsize=(12, 8))
    sns.heatmap(heatmap_df, cmap='coolwarm', annot=True, fmt="d", cbar_kws={'label': 'Frequency of Insertions'})

    # 添加标题和轴标签
    # plt.title(f'Density of LTR Insertion Times', fontsize=16)
    plt.xlabel('Insertion Time (Million years ago)', fontsize=12)
    # plt.ylabel('Genomes', fontsize=12)

    # 显示图形
    plt.tight_layout()
    # plt.show()
    plt.savefig(output_pdf, format='pdf', bbox_inches='tight')

# work_dir = '/home/hukang/test/HiTE/demo'
# log = Logger(work_dir + '/HiTE.log', level='debug')

def parse_coordinates(header):
    """解析header中的坐标信息"""
    parts = header.split('-')
    chr_name = parts[0]
    start = int(parts[1])
    end = int(parts[2])
    return (chr_name, start, end)

def calculate_overlap(coord1, coord2):
    """计算两条序列的重叠比例"""
    chr1, start1, end1 = coord1
    chr2, start2, end2 = coord2

    # 如果染色体不同，直接返回0
    if chr1 != chr2:
        return 0

    # 计算重叠区域
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap_length = max(0, overlap_end - overlap_start)

    # 计算每条序列的长度
    length1 = end1 - start1
    length2 = end2 - start2

    # 计算重叠比例
    overlap_ratio1 = overlap_length / length1
    overlap_ratio2 = overlap_length / length2

    # 返回较大的重叠比例
    return min(overlap_ratio1, overlap_ratio2)

def remove_duplicates(input_fasta, output_fasta, threshold=0.95):
    """去除重叠比例超过阈值的重复序列"""
    input_names, input_contigs = read_fasta(input_fasta)
    unique_records = []
    for record in input_names:
        coords = parse_coordinates(record)
        is_duplicate = False
        # 检查是否与已保留的序列重叠超过阈值
        for unique_record in unique_records:
            unique_coords = parse_coordinates(unique_record)
            overlap_ratio = calculate_overlap(coords, unique_coords)
            if overlap_ratio >= threshold:
                is_duplicate = True
                break
        # 如果不是重复序列，则保留
        if not is_duplicate:
            unique_records.append(record)
    unique_contigs = {}
    for name in unique_records:
        unique_contigs[name] = input_contigs[name]
    # 将去重后的序列写入新的FASTA文件
    store_fasta(unique_contigs, output_fasta)

# from matplotlib_venn import venn3
# from itertools import combinations

def run_command(species):
    cur_dir = os.path.join(work_dir, species)
    intact_LTR_path = os.path.join(cur_dir, 'intact_LTR.fa')
    intact_LTR_list = os.path.join(cur_dir, 'intact_LTR.list')
    cmd = f'python {tool_dir}/recompute_LTR_insertion_time.py --intact_LTR_path {intact_LTR_path} --intact_LTR_list {intact_LTR_list} --miu {miu}'
    print(f"Running command for {species}: {cmd}")
    os.system(cmd)


def parse_panHiTE_miss_bed(bed_file):
    """
    解析panHiTE_miss.out.bed文件，提取TE拷贝信息
    返回: TE_copies, TE_full_copies
    """
    TE_copies = defaultdict(list)
    TE_full_copies = defaultdict(list)

    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue

            # 解析TE信息
            info_parts = parts[3].split(';')
            if len(info_parts) < 14:
                continue

            TE_name = info_parts[9]
            direction = info_parts[8]
            div = float(info_parts[1])
            # 解析位置信息
            try:
                if direction == '+':
                    TE_start = int(info_parts[11])
                    TE_end = int(info_parts[12])
                    TE_remain = int(info_parts[13].replace('(', '').replace(')', ''))
                else:
                    TE_start = int(info_parts[13])
                    TE_end = int(info_parts[12])
                    TE_remain = int(info_parts[11].replace('(', '').replace(')', ''))

                # 计算TE长度和完整性
                total_TE_length = TE_end + TE_remain
                TE_length = TE_end - TE_start
                completeness_ratio = TE_length / total_TE_length if total_TE_length > 0 else 0

                # 存储拷贝信息
                copy_info = (TE_start, TE_end, TE_remain, div)
                TE_copies[TE_name].append(copy_info)

                # 判断是否为全长拷贝（完整性≥95%）
                if completeness_ratio >= 0.95:
                    TE_full_copies[TE_name].append(copy_info)

            except (ValueError, IndexError) as e:
                print(f"警告: 解析行时出错: {line.strip()}, 错误: {e}")
                continue

    return TE_copies, TE_full_copies


def load_te_metadata(repbase_path, panEDTA_unique_file, total_repbase_header):
    """
    加载TE元数据：类型和长度信息
    """
    # 从Repbase加载所有TE信息
    repbase_names, repbase_contigs = read_fasta_v1(repbase_path)

    TE_types = {}
    TE_lengths = {}

    class_map = {'LTR Retrotransposon': 'LTR', 'Non-LTR Retrotransposon': 'Non-LTR', 'DNA transposon': 'DNA', 'EnSpm/CACTA': 'CACTA'}

    for name in repbase_names:
        parts = name.split('\t')
        if len(parts) >= 2:
            TE_name = parts[0]
            classification = parts[1]
            if classification in class_map:
                classification = class_map[classification]
            TE_types[TE_name] = classification
            TE_lengths[TE_name] = len(repbase_contigs[name])

    # 从panEDTA文件加载特定的TE子集
    panEDTA_types = {}
    panEDTA_lengths = {}
    with open(panEDTA_unique_file, 'r') as f:
        for line in f:
            TE_name = line.strip()
            if TE_name in TE_types:
                panEDTA_types[TE_name] = TE_types[TE_name]
                panEDTA_lengths[TE_name] = TE_lengths[TE_name]
            else:
                print(f"警告: TE名称 {TE_name} 在Repbase中未找到")

    # 其余来自Repbase的TE
    other_types = {}
    other_lengths = {}
    with open(total_repbase_header, 'r') as f:
        for line in f:
            TE_name = line.strip()
            if TE_name not in panEDTA_types:
                other_types[TE_name] = TE_types[TE_name]
                other_lengths[TE_name] = TE_lengths[TE_name]

    return panEDTA_types, panEDTA_lengths, other_types, other_lengths, TE_types, TE_lengths


def plot_TE_distribution(TE_types, TE_lengths, TE_full_copies, output_dir):
    """
    绘制TE分布图：全长拷贝分布、类型数量、长度分布、平均identity分布
    """
    # 设置全局字体为DejaVu Sans
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans']

    # 统计每个类型的TE数量
    type_counter = Counter(TE_types.values())

    # 按类型分组数据
    type_lengths = defaultdict(list)
    type_copy_counts = defaultdict(list)  # 存储每个TE的拷贝数量
    type_avg_identities = defaultdict(list)  # 存储每个TE的平均identity

    for te_name, te_type in TE_types.items():
        # 长度数据
        if te_name in TE_lengths:
            type_lengths[te_type].append(TE_lengths[te_name])

        # 拷贝数量和identity数据
        if te_name in TE_full_copies:
            copies_list = TE_full_copies[te_name]
            copy_count = len(copies_list)
            type_copy_counts[te_type].append(copy_count)

            # 计算该TE所有拷贝的平均identity
            identities = []
            for copy_info in copies_list:
                if len(copy_info) >= 4:  # 确保有div值
                    div = copy_info[3]
                    identity = float(100 - div) /100
                    identities.append(identity)

            if identities:  # 如果有有效的identity数据
                avg_identity = np.mean(identities)
                type_avg_identities[te_type].append(avg_identity)


    # 计算统计信息
    type_avg_length = {te_type: np.mean(lengths) if lengths else 0
                       for te_type, lengths in type_lengths.items()}
    type_median_length = {te_type: np.median(lengths) if lengths else 0
                          for te_type, lengths in type_lengths.items()}
    type_avg_copies = {te_type: np.mean(copies) if copies else 0
                       for te_type, copies in type_copy_counts.items()}
    type_median_copies = {te_type: np.median(copies) if copies else 0
                          for te_type, copies in type_copy_counts.items()}
    type_avg_of_avg_identity = {te_type: np.mean(identities) if identities else 0
                                for te_type, identities in type_avg_identities.items()}
    type_median_of_avg_identity = {te_type: np.median(identities) if identities else 0
                                   for te_type, identities in type_avg_identities.items()}

    # 按数量排序类型
    types_sorted = [k for k, v in sorted(type_counter.items(),
                                         key=lambda x: x[1], reverse=True)]
    counts_sorted = [type_counter[t] for t in types_sorted]

    # 创建图形 - 2行2列布局
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))

    # 图1: TE全长拷贝数量分布（左上）
    copy_bins = [0, 1, 5, 10, 20, 50, 100, float('inf')]
    copy_labels = ['1', '2-5', '6-10', '11-20', '21-50', '51-100', '>100']

    copy_distribution = {}
    for te_type in types_sorted:
        copy_distribution[te_type] = [0] * len(copy_labels)
        for copy_count in type_copy_counts[te_type]:
            for i, (start, end) in enumerate(zip(copy_bins[:-1], copy_bins[1:])):
                if start <= copy_count < end:
                    copy_distribution[te_type][i] += 1
                    break

    bottom = np.zeros(len(copy_labels))
    colors_viridis = plt.cm.viridis(np.linspace(0, 1, len(types_sorted)))

    for i, te_type in enumerate(types_sorted):
        ax1.bar(copy_labels, copy_distribution[te_type], bottom=bottom,
                label=te_type, color=colors_viridis[i], alpha=0.8)
        bottom += copy_distribution[te_type]

    ax1.set_title('(a) TE Full-length Copy Distribution by Type',
                  fontsize=18, fontweight='bold', fontname='DejaVu Sans')
    ax1.set_xlabel('Copy Number Range', fontsize=18, fontname='DejaVu Sans')
    ax1.set_ylabel('Number of TEs', fontsize=18, fontname='DejaVu Sans')
    ax1.tick_params(axis='x', rotation=45, labelsize=15)
    ax1.tick_params(axis='y', labelsize=15)
    ax1.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # 图2: TE类型数量分布（右上）
    colors_set3 = plt.cm.Set3(np.linspace(0, 1, len(types_sorted)))
    bars2 = ax2.bar(types_sorted, counts_sorted, color=colors_set3)
    ax2.set_title('(b) TE Classification Distribution',
                  fontsize=18, fontweight='bold', fontname='DejaVu Sans')
    ax2.set_xlabel('TE Classification', fontsize=18, fontname='DejaVu Sans')
    ax2.set_ylabel('Number of TEs', fontsize=18, fontname='DejaVu Sans')
    ax2.tick_params(axis='x', rotation=45, labelsize=15)
    ax2.tick_params(axis='y', labelsize=15)
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # 添加数量标签
    for bar, count in zip(bars2, counts_sorted):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1,
                 f'{count}', ha='center', va='bottom', fontsize=9, fontname='DejaVu Sans')

    # 图3: TE长度分布（左下）
    length_bins = [0, 500, 1000, 2000, 5000, 10000, float('inf')]
    bin_labels = ['<500', '500-1k', '1k-2k', '2k-5k', '5k-10k', '>10k']

    length_distribution = {}
    for te_type in types_sorted:
        length_distribution[te_type] = [0] * len(bin_labels)
        for length in type_lengths[te_type]:
            for i, (start, end) in enumerate(zip(length_bins[:-1], length_bins[1:])):
                if start <= length < end:
                    length_distribution[te_type][i] += 1
                    break

    bottom = np.zeros(len(bin_labels))
    for i, te_type in enumerate(types_sorted):
        ax3.bar(bin_labels, length_distribution[te_type], bottom=bottom,
                label=te_type, color=colors_viridis[i], alpha=0.8)
        bottom += length_distribution[te_type]

    ax3.set_title('(c) TE Length Distribution by Type',
                  fontsize=18, fontweight='bold', fontname='DejaVu Sans')
    ax3.set_xlabel('Length Range (bp)', fontsize=18, fontname='DejaVu Sans')
    ax3.set_ylabel('Number of TEs', fontsize=18, fontname='DejaVu Sans')
    ax3.tick_params(axis='x', rotation=45, labelsize=15)
    ax3.tick_params(axis='y', labelsize=15)
    ax3.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # 图4: TE平均Identity分布（右下）- 修改为按TE平均值
    identity_bins = [0, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 1.0]
    identity_labels = ['<0.7', '0.7-0.8', '0.8-0.9', '0.9-0.95', '0.95-0.98', '0.98-0.99', '≥0.99']

    identity_distribution = {}
    for te_type in types_sorted:
        identity_distribution[te_type] = [0] * len(identity_labels)
        for avg_identity in type_avg_identities[te_type]:
            for i, (start, end) in enumerate(zip(identity_bins[:-1], identity_bins[1:])):
                if start <= avg_identity < end:
                    identity_distribution[te_type][i] += 1
                    break
            else:
                # avg_identity == 1.0 的情况
                identity_distribution[te_type][-1] += 1

    bottom = np.zeros(len(identity_labels))
    for i, te_type in enumerate(types_sorted):
        ax4.bar(identity_labels, identity_distribution[te_type], bottom=bottom,
                label=te_type, color=colors_viridis[i], alpha=0.8)
        bottom += identity_distribution[te_type]

    ax4.set_title('(d) TE Average Identity Distribution by Type',
                  fontsize=18, fontweight='bold', fontname='DejaVu Sans')
    ax4.set_xlabel('Average Identity Range', fontsize=18, fontname='DejaVu Sans')
    ax4.set_ylabel('Number of TEs', fontsize=18, fontname='DejaVu Sans')  # 注意：这里是TE数量
    ax4.tick_params(axis='x', rotation=45, labelsize=15)
    ax4.tick_params(axis='y', labelsize=15)
    ax4.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # 添加图例
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=18, prop={'family': 'DejaVu Sans'})
    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=18, prop={'family': 'DejaVu Sans'})
    ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=18, prop={'family': 'DejaVu Sans'})

    # 调整布局并保存
    plt.tight_layout()
    output_path = os.path.join(output_dir, 'TE_comprehensive_distribution.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()

    # 打印统计摘要
    print_summary(types_sorted, counts_sorted, type_avg_length,
                  type_median_length, type_avg_copies, type_median_copies,
                  type_avg_of_avg_identity, type_median_of_avg_identity)


def print_summary(types_sorted, counts_sorted, avg_lengths,
                  median_lengths, avg_copies, median_copies,
                  avg_identities, median_identities):
    """打印统计摘要"""
    print("\n" + "=" * 50)
    print("TE Distribution Summary")
    print("=" * 50)
    print(f"Total number of TEs: {sum(counts_sorted)}")
    print(f"Number of TE types: {len(types_sorted)}")

    print("\nTop 5 TE types by count:")
    for i in range(min(5, len(types_sorted))):
        te_type = types_sorted[i]
        print(f"  {te_type}: {counts_sorted[i]} TEs")

    print("\nLength statistics (top 5 types):")
    for te_type in types_sorted[:5]:
        if te_type in avg_lengths:
            print(f"  {te_type}: avg={avg_lengths[te_type]:.1f}bp, "
                  f"median={median_lengths[te_type]:.1f}bp")

    print("\nCopy number statistics (top 5 types):")
    for te_type in types_sorted[:5]:
        if te_type in avg_copies:
            print(f"  {te_type}: avg={avg_copies[te_type]:.1f}, "
                  f"median={median_copies[te_type]:.1f}")

    print("\nAverage identity statistics (top 5 types):")
    for te_type in types_sorted[:5]:
        if te_type in avg_identities and avg_identities[te_type] > 0:
            print(f"  {te_type}: avg={avg_identities[te_type]:.3f}, "
                  f"median={median_identities[te_type]:.3f}")

if __name__ == '__main__':
    # # 画插入时间图
    # data_dir = '/home/hukang/test/data'
    # output_pdf = '/home/hukang/test/intact_LTR_insert_time.pdf'
    # intact_ltr_paths = []
    # for genome_name in os.listdir(data_dir):
    #     input_dir = data_dir + '/' + genome_name
    #     cur_list = input_dir + '/intact_LTR.list'
    #     intact_ltr_paths.append((genome_name, cur_list))
    # draw_intact_LTR_insert_time(intact_ltr_paths, output_pdf)



    """主函数"""
    # 设置目录和文件路径
    cur_dir = '/home/hukang/panHiTE_reviewer_test/panHiTE_uncover_TEs_analysis'

    # 文件路径
    panHiTE_miss_bed = os.path.join(cur_dir, 'panHiTE_miss.out.bed')
    total_repbase_header = os.path.join(cur_dir, 'panHiTE_miss.txt')
    panEDTA_unique_file = os.path.join(cur_dir, 'panEDTA_vs_Repbase.header.txt')
    repbase_path = os.path.join(cur_dir, 'athrep.ref')
    copies_output_file = os.path.join(cur_dir, 'panHiTE_miss.full_length.txt')

    # 1. 加载TE元数据
    print("Loading TE metadata...")
    panEDTA_types, panEDTA_lengths, other_types, other_lengths, all_TE_types, all_TE_lengths = load_te_metadata(
        repbase_path, panEDTA_unique_file, total_repbase_header)

    # 2. 解析panHiTE结果
    print("Parsing panHiTE miss results...")
    TE_copies, TE_full_copies = parse_panHiTE_miss_bed(panHiTE_miss_bed)

    # 3. 保存全长拷贝信息
    print("Saving full-length copy information...")
    with open(copies_output_file, 'w') as f:
        for te_name, copies_list in TE_full_copies.items():
            f.write(f"{te_name}\t{len(copies_list)}\t{copies_list}\n")

    # 4. 筛选panEDTA相关的全长拷贝
    panEDTA_copies = {te_name: copies for te_name, copies in TE_copies.items()
                           if te_name in panEDTA_types}
    panEDTA_full_copies = {te_name: copies for te_name, copies in TE_full_copies.items()
                           if te_name in panEDTA_types}

    print(f"Total TEs in panEDTA: {len(panEDTA_types)}")
    print(f"Total TEs in Other Repbase: {len(other_types)}")
    print(f"TEs with full-length copies: {len(panEDTA_full_copies)}")

    # 5. 筛选其他Repbase相关的全长拷贝
    other_copies = {te_name: copies for te_name, copies in TE_copies.items()
                         if te_name not in panEDTA_types}
    other_full_copies = {te_name: copies for te_name, copies in TE_full_copies.items()
                         if te_name not in panEDTA_types}

    # 6. 绘制分布图
    print("Plotting TE distributions...")
    # plot_TE_distribution(panEDTA_types, panEDTA_lengths, panEDTA_full_copies, cur_dir)
    # print(len(panEDTA_copies), len(panEDTA_full_copies))
    plot_TE_distribution(other_types, other_lengths, other_full_copies, cur_dir)
    print(len(other_copies), len(other_full_copies))


    # cur_dir = '/public/home/hpc194701009/ath_pan_genome/pan_genome/ath/32_ath/panHiTE2_32_ath'
    # original_candidates = cur_dir + '/865.fasta'
    # keep_tir_low_copy = cur_dir + '/865_tir_low.fasta'
    # keep_helitron_low_copy = cur_dir + '/865_helitron_low.fasta'
    # keep_non_ltr_low_copy = cur_dir + '/865_non_ltr_low.fasta'
    # original_candidates = cur_dir + '/514.fasta'
    # keep_tir_low_copy = cur_dir + '/514_tir_low.fasta'
    # keep_helitron_low_copy = cur_dir + '/514_helitron_low.fasta'
    # keep_non_ltr_low_copy = cur_dir + '/514_non_ltr_low.fasta'
    # original_candidates = cur_dir + '/415.fasta'
    # keep_tir_low_copy = cur_dir + '/415_tir_low.fasta'
    # keep_helitron_low_copy = cur_dir + '/415_helitron_low.fasta'
    # keep_non_ltr_low_copy = cur_dir + '/415_non_ltr_low.fasta'
    # os.system('cp ' + original_candidates + ' ' + keep_tir_low_copy)
    # os.system('cp ' + original_candidates + ' ' + keep_helitron_low_copy)
    # os.system('cp ' + original_candidates + ' ' + keep_non_ltr_low_copy)
    # genome_dir = '/public/home/hpc194701009/ath_pan_genome/pan_genome/ath/32_ath/genomes/'
    # genome_path = '/public/home/hpc194701009/ath_pan_genome/pan_genome/ath/32_ath/genome_list_bak'
    # genome_paths = []
    # with open(genome_path, 'r') as f_r:
    #     for line in f_r:
    #         genome_name = line.replace('\n', '')
    #         genome_paths.append((genome_name, genome_dir + '/' + genome_name))
    # # temp_dir = cur_dir + '/514_copies'
    # # temp_dir = cur_dir + '/865_copies'
    # temp_dir = cur_dir + '/415_copies'
    # threads = 48
    # log = Logger(cur_dir + '/pan_recover_low_copy_TEs.log', level='debug')
    # get_copies_dir = os.path.join(temp_dir, 'get_copies')
    # pan_genome_copies_dict = get_pan_genome_copies(keep_tir_low_copy, keep_helitron_low_copy, keep_non_ltr_low_copy, genome_paths, threads, get_copies_dir, log)
    #
    # # store connected_repeats for testing
    # pan_genome_copies_dict_file = temp_dir + '/pan_genome_copies_dict.csv'
    # with codecs.open(pan_genome_copies_dict_file, 'w', encoding='utf-8') as f:
    #     json.dump(pan_genome_copies_dict, f)
    #
    # project_dir = '/public/home/hpc194701009/test/test/HiTE'
    # subset_script_path = project_dir + '/tools/ready_for_MSA.sh'
    # plant = 0
    # debug = 0
    # result_type = 'cons'
    # # 依次对tir, helitron, non_ltr 类型的TE调用动态边界调整算法，过滤出真实的TE
    # for TE_type in pan_genome_copies_dict.keys():
    #     cur_temp_dir, batch_member_files = pan_genome_copies_dict[TE_type]
    #     real_TEs = os.path.join(temp_dir, 'real_' + TE_type + '.fa')
    #     for batch_member_file in batch_member_files:
    #         (query_name, cur_seq, trunc_member_file, extend_member_file) = batch_member_file
    #         is_TE_from_align_file(extend_member_file, query_name, cur_seq, temp_dir, subset_script_path, plant, TE_type, debug, result_type)


    # te_names, te_contigs = read_fasta(original_candidates)
    # not_found_names = []
    # msa_dir = cur_dir + '/panHiTE_unique_msa_files'
    # for file_name in os.listdir(msa_dir):
    #     te_name = file_name.split('.')[0]
    #     del te_contigs[te_name]
    # print(te_contigs.keys())
    # print(len(te_contigs.keys()))
    #
    # remain_candidates = cur_dir + '/371.fasta'
    # store_fasta(te_contigs, remain_candidates)





    # max_copy_num = 100
    # threads = 40
    # temp_dir = cur_dir + '/get_copies'
    # for ref_name, reference in genomes:
    #     ref_names, ref_contigs = read_fasta(reference)
    #     all_copies = get_full_length_copies_minimap2(panHiTE_miss, reference, temp_dir, max_copy_num, threads)
    #     store_dir = cur_dir + '/' + ref_name
    #     if not os.path.exists(store_dir):
    #         os.makedirs(store_dir)
    #     copies_file = store_dir + '/TE_copy_numbers.txt'
    #     with open(copies_file, 'w') as f_save:
    #         for te_name, copies_list in all_copies.items():
    #             copy_count = len(copies_list)
    #             f_save.write(f"{te_name}\t{copy_count}\n")
    #     for query_name in all_copies:
    #         query_file = store_dir + '/' + query_name + '.fa'
    #         query_contigs = {}
    #         query_contigs[query_name] = panHiTE_miss_contigs[query_name]
    #         store_fasta(query_contigs, query_file)
    #
    #         query_copy_file = store_dir + '/' + query_name + '.copies.fa'
    #         copy_contigs = {}
    #         for copy in all_copies[query_name]:
    #             subject_name = copy[0]
    #             subject_start = copy[1]
    #             subject_end = copy[2]
    #             copy_seq = ref_contigs[subject_name][subject_start: subject_end]
    #             copy_name = subject_name + ':' + str(subject_start) + '-' + str(subject_end)
    #             copy_contigs[copy_name] = copy_seq
    #         store_fasta(copy_contigs, query_copy_file)
    #         blastn2Results_path = sequence2sequenceBlastn(query_file, query_copy_file)
    #
    #         avg_identity = 0
    #         line_num = 0
    #         get_longest_repeats_v4, blastn2Results_path, fixed_extend_base_threshold, max_single_repeat_len, debug




    # cur_dir = '/home/hukang/test/HiTE/demo/out29/NeuralTE_all'
    # classified_path = os.path.join(cur_dir, 'classified_TE.fa')
    # panTE_path = os.path.join(cur_dir, 'panTE.merge_recover.fa')
    # classified_names, classified_contigs = read_fasta(classified_path)
    # panTE_names, panTE_contigs = read_fasta(panTE_path)
    # classified_dict = {}
    # for name in classified_names:
    #     parts = name.split('#')
    #     raw_name = parts[0]
    #     class_name = parts[1]
    #     classified_dict[raw_name] = class_name
    #
    # new_panTE_contigs = {}
    # for name in panTE_names:
    #     parts = name.split('#')
    #     raw_name = parts[0]
    #     class_name = parts[1]
    #     if raw_name in classified_dict:
    #         cur_class_name = classified_dict[raw_name]
    #         new_name = raw_name+'#'+cur_class_name
    #     else:
    #         new_name = name
    #     new_panTE_contigs[new_name] = panTE_contigs[name]
    # store_fasta(new_panTE_contigs, panTE_path)

    # sf_dict = {'RLC': 'LTR/Copia', 'RLG': 'LTR/Gypsy', 'RLB': 'LTR/Pao',
    #            'RLR': 'LTR/ERV', 'RLE': 'LTR/ERV', 'RYD': 'LTR/DIRS',
    #            'RYN': 'LTR/Ngaro', 'RYV': 'LTR/Viper', 'RPP': 'LINE/Penelope',
    #            'RIR': 'LINE/R2', 'RIT': 'LINE/RTE-RTE', 'RIJ': 'LINE/I-Jockey',
    #            'RIL': 'LINE/L1', 'RII': 'LINE/I', 'RST': 'SINE/tRNA',
    #            'RSL': 'SINE/7SL', 'RSS': 'SINE/5S', 'DTT': 'DNA/TcMar',
    #            'DTA': 'DNA/hAT', 'DTM': 'DNA/MULE', 'DTE': 'DNA/Merlin',
    #            'DTR': 'DNA/CMC-Transib', 'DTP': 'DNA/P', 'DTB': 'DNA/PiggyBac',
    #            'DTH': 'DNA/PIF-Harbinger', 'DTC': 'DNA/CMC-EnSpm', 'DYC': 'DNA/Crypton',
    #            'DHH': 'RC/Helitron', 'DMM': 'DNA/Maverick', 'DTX': 'DNA',
    #            'XXX': 'Unknown', 'RIX': 'LINE', 'RLX': 'LTR', 'SIX': 'SINE', 'DXX': 'DNA'
    #            }
    # clariTeRep_path = os.path.join(cur_dir, 'clariTeRep.fna')
    # clariTeRep_names, clariTeRep_contigs = read_fasta_v1(clariTeRep_path)
    # new_clariTeRep_contigs = {}
    # for name in clariTeRep_names:
    #     parts = name.split(' ')
    #     raw_name = parts[0]
    #     super_family = parts[2].split('=')[1]
    #     class_name = sf_dict[super_family]
    #     new_name = raw_name + '#' + class_name
    #     new_clariTeRep_contigs[new_name] = clariTeRep_contigs[name]
    # new_clariTeRep_path = os.path.join(cur_dir, 'clariTeRep.rename.fna')
    # store_fasta(new_clariTeRep_contigs, new_clariTeRep_path)


    # cur_dir = '/home/hukang/test/HiTE/demo/out43'
    # tbl_file = cur_dir + '/wheat_A.fasta.sorted.gff.tbl'
    # rawname_value = {}
    # with open(tbl_file, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         elements = line.split()
    #         rawname_value[elements[0]] = elements[3]
    # lib_path = cur_dir + '/panTE.merge_recover.fa'
    # names, contigs = read_fasta(lib_path)
    #
    # target_rawname_value = {}
    # for name in names:
    #     parts = name.split('#')
    #     class_name = parts[1]
    #     raw_name = parts[0]
    #     # if 'MULE' in class_name or 'Harbinger' in class_name or 'TcMar' in class_name or 'Unknown' in class_name:
    #     if 'MULE' in class_name or 'Harbinger' in class_name or 'TcMar' in class_name:
    #         target_rawname_value[raw_name] = rawname_value[raw_name]
    #
    # # 排序（去掉百分号并转成 float）
    # sorted_data = dict(
    #     sorted(target_rawname_value.items(), key=lambda x: float(x[1].strip('%')), reverse=True)
    # )
    # print(sorted_data)


    # lib = cur_dir + '/panTE.merge_recover.fa'
    # org_lib = cur_dir + '/panTE.merge_recover.original.fa'
    # total_lib = cur_dir + '/total_TE.cons.fa'
    # names, contigs = read_fasta(lib)
    # total_names, total_contigs = read_fasta(total_lib)
    # total_dict = {}
    # for name in total_names:
    #     parts = name.split('#')
    #     raw_name = parts[0]
    #     class_name = parts[1]
    #     total_dict[raw_name] = class_name
    #
    # new_lib_contigs = {}
    # for name in names:
    #     seq = contigs[name]
    #     parts = name.split('#')
    #     raw_name = parts[0]
    #     class_name = parts[1]
    #     if class_name == 'Unknown':
    #         if raw_name in total_dict:
    #             new_class_name = total_dict[raw_name]
    #             new_name = raw_name + '#' + new_class_name
    #         else:
    #             new_name = name
    #     else:
    #         new_name = name
    #     new_lib_contigs[new_name] = seq
    # store_fasta(new_lib_contigs, org_lib)

    # tmp_output_dir = cur_dir
    # ref_index = 0
    # tmp_blast_dir = os.path.join(tmp_output_dir, f'longest_repeats_blast_{ref_index}')
    # if not os.path.exists(tmp_blast_dir):
    #     os.makedirs(tmp_blast_dir)

    # filter_tandem_file = tmp_output_dir + '/genome.cut0.fa'
    # # 读取 FASTA 文件
    # repeat_names, repeat_contigs = read_fasta(filter_tandem_file)
    #
    # # 分割序列并存储
    # target_tuples = split_and_store_sequences(repeat_names, repeat_contigs, tmp_blast_dir, base_threshold=1_000_000)





    # file_name = 'Chr1-722405358-722410410-tir_12-tsd_ATG.blast.bed.fa'
    # file_path = cur_dir + file_name
    # names, contigs = read_fasta(file_path)
    # rec_file_path = cur_dir + file_name + '.rec.fa'
    # rec_contigs = {}
    # for name in names:
    #     seq = contigs[name]
    #     if len(seq) > 1000:
    #         rec_seq = seq[0:500] + seq[-500:]
    #         rec_contigs[name] = rec_seq
    # store_fasta(rec_contigs, rec_file_path)


    # cmc_path = os.path.join(cur_dir, 'CMC-EnSpm.fa')
    # cmc_names, cmc_contigs = read_fasta(cmc_path)
    # len_dict = {}
    # for name in cmc_names:
    #     parts = name.split('#')
    #     raw_name = parts[0]
    #     class_name = parts[1]
    #     len_dict[raw_name] = len(cmc_contigs[name])
    # TE_annotation_path = os.path.join(cur_dir, 'TE_v1.1.gff3')
    # cmc_copies = {}
    # with open(TE_annotation_path, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.split('\t')
    #         raw_name = parts[8].split('Motif:')[1].split(';')[0]
    #         if raw_name in len_dict:
    #             pos = parts[8].split('Motif:')[1].split(';')[-1].strip()
    #             if pos == '':
    #                 start = int(parts[9])
    #                 end = int(parts[10])
    #             else:
    #                 pos_parts = pos.split(' ')
    #                 start = int(pos_parts[0])
    #                 end = int(pos_parts[1])
    #             cur_len = abs(end-start)
    #             if float(cur_len)/len_dict[raw_name] >= 0.99:
    #                 if raw_name not in cmc_copies:
    #                     cmc_copies[raw_name] = []
    #                 copies = cmc_copies[raw_name]
    #                 copies.append((line, cur_len))
    # cmc_copies_path = os.path.join(cur_dir, 'cmc.copies')
    # # save_to_file(cmc_copies, cmc_copies_path)
    #
    # cmc_copies = load_from_file(cmc_copies_path)
    # copy_lens = {}
    # for name in cmc_copies.keys():
    #     total_len = 0
    #     for line, cur_len in cmc_copies[name]:
    #         total_len += cur_len
    #     copy_lens[name] = total_len
    # sorted_copy_lens = sorted(copy_lens.items(), key=lambda x: -x[1])
    # print(sorted_copy_lens, len(sorted_copy_lens))
    # # for copy in cmc_copies['TREP1802']:
    # #     print(copy[0])

    # [('rph7_rep_0240', 7577758), ('rph7_rep_0208', 5236507), ('ctgF_rep_0341', 3825249), ('ctgG_rep_0154', 3431853), ('TREP1802', 3180454), ('TREP3415', 2791784), ('sr2_rep_0101', 2780524), ('ctgE_rep_0053', 2682375), ('ctgL_rep_0051', 2546533), ('ctgF_rep_0445', 2491249), ('ctgD_rep_0144', 2235961), ('ctgH_rep_0316', 2178090), ('ctgA_rep_0463', 2041075), ('rph7_rep_0583', 1921262), ('ctgF_rep_0343', 1800030), ('TREP1302', 1530764), ('ctgB_rep_0167', 1468240), ('TREP1803', 1409710), ('ctgC_rep_0034', 1195118), ('ctgH_rep_0033', 1158818), ('TREP764', 1126303), ('ctgH_rep_0299', 988485), ('ctgF_rep_0158', 878678), ('ctgF_rep_0024', 857972), ('TREP1304', 840452), ('ctgE_rep_0141', 799408), ('TREP2249', 716391), ('TREP3416', 697394), ('sr2_rep_0031', 678385), ('ctgF_rep_0212', 660717), ('TREP773', 598807), ('TREP1303', 596905), ('ctgF_rep_0275', 586565), ('sr2_rep_0097', 575566), ('ctgL_rep_0039', 545074), ('TREP780', 487336), ('TREP747', 487308), ('ctgF_rep_0435', 447786), ('ctgL_rep_0098', 427496), ('ctgL_rep_0056', 418757), ('ctgF_rep_0313', 417936), ('ctgA_rep_0070', 397124), ('ctgB_rep_0250', 392253), ('ctgF_rep_0301', 362369), ('ctgE_rep_0123', 346704), ('TREP2205', 346311), ('ctgH_rep_0260', 339801), ('rph7_rep_0509', 339420), ('ctgI_rep_0115', 335765), ('ctgH_rep_0281', 321755), ('rph7_rep_0420', 315874), ('TREP1318', 261999), ('ctgB_rep_0036', 251966), ('rph7_rep_0466', 247057), ('rph7_rep_0236', 234041), ('rph7_rep_0012', 220665), ('TREP783', 202520), ('ctgF_rep_0047', 193851), ('ctgA_rep_0202', 184812), ('TREP757', 177729), ('TREP3029', 175141), ('ctgL_rep_0035', 171747), ('ctgF_rep_0297', 164101), ('ctgE_rep_0185', 160637), ('ctgB_rep_0224', 125622), ('sr2_rep_0100', 125185), ('ctgA_rep_0292', 109636), ('ctgA_rep_0262', 105224), ('TREP1477', 103221), ('TREP3021', 102663), ('TREP775', 99009), ('ctgF_rep_0314', 86999), ('rph7_rep_0358', 76602), ('TREP3018', 71810), ('sr2_rep_0138', 71469), ('TREP746', 71277), ('ctgF_rep_0151', 70184), ('rph7_rep_01895', 68039), ('ctgH_rep_00165', 65016), ('TREP3005', 64112), ('rph7_rep_0588', 62950), ('ctgL_rep_0062', 60464), ('rph7_rep_0431', 58360), ('TREP2065', 56800), ('rph7_rep_0229', 55677), ('ctgD_rep_0052', 55536), ('TREP3026', 54436), ('rph7_rep_0475', 52451), ('ctgL_rep_0023', 49425), ('ctgD_rep_0065', 48659), ('ctgG_rep_0005', 43883), ('ctgE_rep_0083', 39423), ('rph7_rep_0306', 34228), ('ctgE_rep_0109', 34128), ('TREP2221', 33341), ('ctgA_rep_0232', 32852), ('ctgD_rep_0097', 32592), ('rph7_rep_0433', 30955), ('TREP3022', 30239), ('sr2_rep_0022', 30192), ('ctgF_rep_0130', 28969), ('rph7_rep_04301', 28685), ('sr2_rep_0038', 28091), ('TREP3020', 27398), ('TREP1241', 26998), ('rph7_rep_0009', 26311), ('TREP2064', 24868), ('TREP779', 24090), ('TREP1704', 24050), ('ctgF_rep_0370', 23368), ('TREP3009', 22682), ('ctgF_rep_0300', 21746), ('TREP3008', 21260), ('ctgE_rep_0134', 20972), ('TREP784', 20788), ('TREP2243', 19697), ('rph7_rep_0242', 17499), ('rph7_rep_0173', 16985), ('ctgB_rep_0105', 16254), ('ctgF_rep_0304', 15000), ('rph7_rep_0232', 14089), ('TREP3006', 14070), ('ctgF_rep_0291', 12296), ('TREP787', 9098), ('TREP3030', 6596), ('ctgF_rep_0265', 6468), ('TREP1478', 4009), ('TREP763', 3717), ('rph7_rep_0015', 3552), ('sr2_rep_0137', 2837), ('TREP1245', 2510), ('sr2_rep_0194', 2495), ('rph7_rep_0297', 2450), ('sr2_rep_0168', 2404), ('ctgD_rep_0253', 2193), ('ctgC_rep_0031', 1440), ('ctgH_rep_0090', 889), ('ctgC_rep_0088', 449), ('ctgF_rep_0398', 112)]
    # for name in cmc_copies.keys():
    #     print(name, len(cmc_copies[name]))

    # def is_overlapping(start1, end1, start2, end2):
    #     """
    #     判断两个基因组坐标是否重叠
    #     :param start1: 第一个区间的起始点
    #     :param end1: 第一个区间的结束点
    #     :param start2: 第二个区间的起始点
    #     :param end2: 第二个区间的结束点
    #     :return: 如果重叠返回 True，否则返回 False
    #     """
    #     return start1 <= end2 and end1 >= start2
    #
    # blastn_out = os.path.join(cur_dir, '46_target.fa_vs_24_target.fa.out')
    # q_loc = (854330, 867169)
    # with open(blastn_out, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.split('\t')
    #         q_start = int(parts[6])
    #         q_end = int(parts[7])
    #         s_start = int(parts[8])
    #         s_end = int(parts[9])
    #         if is_overlapping(q_start, q_end, q_loc[0], q_loc[1]):
    #             print(line)

    # # 取这两个拷贝出来
    # TE1 = ('Chr1A', 35568574, 35581458)
    # TE2 = ('Chr1A', 297140742, 297153200)
    # TE3 = ('Chr1A', 57854340, 57867216)
    # genome_path = cur_dir + '/Chr1A.fa'
    # ref_names, ref_contigs = read_fasta(genome_path)
    # ref_seq = ref_contigs['Chr1A']
    # TE3_seq = ref_seq[TE3[1]: TE3[2]]
    # print(TE3_seq)
    # # TE1_seq = ref_seq[TE1[1]: TE1[2]]
    # # TE2_seq = ref_seq[TE2[1]: TE2[2]]
    # # print(TE1_seq)
    # # print(TE2_seq)

    # blastn2Results_path = cur_dir + '/test.out'
    # fixed_extend_base_threshold = 4000
    # max_single_repeat_len = 30_000
    # debug = 0
    # a = get_longest_repeats_v4(blastn2Results_path, fixed_extend_base_threshold, max_single_repeat_len, debug)
    # l = load_from_file(a)
    # print(l)

    # tmp_output_dir = '/home/hukang/test/HiTE/demo/test'
    # longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_0.flanked.fa'
    # reference = '/home/hukang/test/HiTE/demo/GCF_000004555.2_CB4_genomic.rename.fna'
    # prev_TE = None
    # flanking_len = 50
    # threads = 1
    # debug = 0
    # split_ref_dir = None
    # all_low_copy_tir = tmp_output_dir + '/all_low_copy_tir.fa'
    # plant = 0
    # ref_index = 0
    # min_TE_len = 80
    # is_recover = 1
    # log = Logger(tmp_output_dir+'/HiTE_tir.log', level='debug')
    # run_TIR_detection(tmp_output_dir, longest_repeats_flanked_path, reference, prev_TE, flanking_len, threads,
    #                   debug, split_ref_dir, all_low_copy_tir, plant, ref_index, min_TE_len, is_recover, log)

    # project_dir = '/home/hukang/test/HiTE'
    # TRsearch_dir = project_dir + '/tools'
    # subset_script_path = project_dir + '/tools/ready_for_MSA.sh'
    # tmp_output_dir = '/home/hukang/test/HiTE/demo/out29/NeuralTE_all/tir_out'
    # ref_index = 0
    # longest_repeats_flanked_path = '/home/hukang/test/HiTE/demo/out29/NeuralTE_all/test.fa'
    # flanking_len = 50
    # threads = 40
    # plant = 1
    # #
    # # tir_tsd_path = tmp_output_dir + '/tir_tsd_' + str(ref_index) + '.fa'
    # # resut_file = tir_tsd_path
    # # tir_tsd_dir = tmp_output_dir + '/tir_tsd_temp_' + str(ref_index)
    # # multi_process_tsd_v1(longest_repeats_flanked_path, tir_tsd_path, tir_tsd_dir, flanking_len, threads,
    # #                      TRsearch_dir, plant)
    #
    #
    # tir_tsd_cons = '/home/hukang/test/HiTE/demo/out29/NeuralTE_all/tir_out/tir_tsd_0.cons.fa'
    # reference = cur_dir + '/Chr1A.fa'
    # debug = 1
    # split_ref_dir = '/home/hukang/test/HiTE/demo/out29/NeuralTE_all/ref_chr'
    # all_low_copy_tir = '/home/hukang/test/HiTE/demo/out29/NeuralTE_all/tir_low_copy.fa'
    # min_TE_len = 80
    # is_recover = 0
    # # is_transposons(tir_tsd_cons, reference, threads, tmp_output_dir, ref_index, log,
    # #                subset_script_path, plant, debug, TRsearch_dir, split_ref_dir, all_low_copy_tir, min_TE_len,
    # #                is_recover)
    #
    # # cur_seq = 'CACTACAAAAAAAAATACACTTCCGTGATGATACATGTTTGTCATAGTAGGTCGCTTTTTTTGTCATGCATGTACATCCATGACAAATTTATGACAGAATCAAGATAGTCATACATGTGCTGTCGTAGAAGTGTTCCATGACATTACCAAAATTATCATCATGGAAGTGTCTAGTTCCATGACGATAAATCGTGCGTCACAGAAGTGCTTTCATCAAGGGTGACCAACACGTGGCATCCACCGTAACGGAACGCCGTTAAGCTATCGGGTTGGGTTTTGGATCCGATAACCCGTTAACAGCCCCGACCAATGGGAAATTTCCACGTGTAAAATTCTTATTGGCCGGACAAAACACGTGTCAGCTCGTCAGTGGGTCAGATAGGCGCCTATGATATGTCGACACGTGGCACGGCCCAACAGTGGCCCATTTAGCTTACAAAGCCGGCCCGTTTAACTTGGTCAAAAGTTAACGGGCTGTCCCATGAAAAGCCTGTTAACGGTCTCTTTGCAAATAGCCCATTTTATGGCCCGTTAACTCACGGCCCATTAGGCCCTAAAGGAAATCGGCCCAACAACATCATGTGGGTTGTCGAATATAATACCAACCCATTTCACTTTTGGCCCATGTATGGCCCACGACGTCTTTCAGCCCATATGAGGCCTTCGTATCTTTCGGCCATTACAGGCCCACGGTTAGTCTAGCCCATAATGAACAGGAATTTTCTTTATACCCATTAACGGCCCGTGATTTACATGGGCCGTTTCCAGCCCGTGTTAGCTGTCGGCCTATTGACGACCCATACGTTCTTGGGCTCCTTTTAGGCCCTTGATTACTTCTGGCCCGTTACTGGCCTGTTCCACTAATGGGCCAAATTCGGCCCATGGCAAGAGTCGGCCCGTTACTGGCCTATTACCCTAATGGGCCAAAATCGTCCCATGGCAAGAGTCGGCCCGTTTCTGGCATGTTAACCTGTTGCCCCGTTTCTAGCCCATCCTATATTCTGGCCCATTAACGACCCATTATGCCTGTGATAGAATTAACCCTTTGTTGTCCTACGGCCTATTAACCTCTTGTTACCGTGCTAGGCCGATACCAATTACGCCCTATTATGGCCCATGTAGACACATTTATCCGACGACCCGAGGCCCACCGATTACAGGCCCATTTATCGACGGCCCGTAGGTGACCCATGGATCCTACGACCCGTATATGGCCCATGGTAGTTGCGGCCACTAGCAAACCAGGAAAAAAGAAGACTAGGAAATAAATAAGGCCGAAACTAACGCTAGGCTATGAAGGCGATTGCACAGCTTACATCCACACGGCATCAAAGATCGCCACCAGTGCAAATATAGGGAACACCCTACACTATAGAAAAATGGCTTGCTGTTTTCTTTAGGCGGTGGCTGCACGTTAAGAATACATTTTGTATCGCACCAAAACAAATTACAACTATATAACAAAAGGAAGGTTGGCATAAAACTTAACAGTCTAGCAGATAAATGCACCACCAGAAGTTCAGAAGCTCAGTTCATTTGACCTGACTGTTACGTTTTCATATTGTACTGTCCAATGGAGCACCATAACCCTCGCCTTCATTTCCCCTAGGCCCCTGAACAACATAGCAAATGTCTGATAGTCTGGTTTCAAAACCATGCATATCTTGACGACATTGAGCAATGTTTCTCCTTGTTTCACAAAGGGTCTCTGTTAGGGAATAGACTTGTGTCTGGAGCGCAGATACAACATTGCTTTGAGCTTACATATCTTGATCATGGGGCGGTTTGGACGATATTGCCTTAACAACCAACCTAGTATTACACAAGAACGTGCTTTTGGCACTGTTAGTGGACAGGTACTGACCCACTGCAGCAAGAGCTGACATTGCGGTTATATTTGCCTCGCCACCTTCAGAAGGTGGCGGTTCCACCATTTTTTCCATAGCTTGCTGAAAAGTTGAGGTTTTACATTAGTAGTACTAGAACAGAAGATATTAAGGAGGCAGCAAAACCAAACAAAATAGCTTTGGTCTATATAGAACTTATTCTGGTGAACCCATGTAAAATATTACTTGACTTTTATTGCAAAGTACATAATAAAACCAGGTTCCAAACATATGACCATGTACCATCGCTAATGTATTGTCTATTTTTTCACATTGTGTTGAGGGCTAAGGATGAAGAGACGAAGTTTATAATACGGTAAGCATAGTATGACGTCATTCATATCACAATACAACTGAGAACAGACAAACATGCAATTGTAATGTTGTGTGGGCAAGTCCAGCACGGAACTAGATTTCGGGTCAAGATGAAAGGAACATCACTCCTCTCTTTTAAGCCTTGTAAGGTGCAATGATACGCGAAGAAAATTGTGTATAATATTGAAAAGAGTAATAGACATTGGCTGCAAACCATGCAGTTCAAGGCACAAGTCAGAGTAACTAGTAGTTAAAAGGCATGAGGATTAAGGACTTACAACAGCGGCTTTGACTGGTGTAGTCATGCCCTTCTTCTTGCTGGTGTGACAGTCCTTGAAGATTTCAACAACATTTGGTTCAGGTACTTTTTGGTCCTTGCGGGCTTTCCTCTGCATTTGGAACAAGTCAATATAGATCTGATAATATGGTACAGCGAAAAGAGTGAAACAGCAGCTAGTTACAAGAGCCTCGCAGTGTGCAATATAGCTACGACATCCTGTCGTCTGTTGGAATTTCACTTTCGTACGGTTGGCCTTGTTCTTTGAACAGTTCACCTACAAGAAGATAGATTTATATGGCACATGCGTAAATAGGACTTCAACATGTGATATGATGACATATGTATAATGTACCTGATACTTCGGATCGGACCAGTGTTTAACGAGGCCTCTCCAGTCTTCATCAGATATATTTTCCACTGGAGATGTTTGGGCAAGTTCACTGTTAGCCTTGCCTTCAAAGTGAGTTTTCCTCAAGTTATACCGATACTGTCACAGAGCAGACATGAAAACATGGGTGCAAGCTTGTCTGGTTGCATCATCTTGGCTATCCAACTTGAACCTCATCTGTTGAAAATTATGGGAGTCTCATTATCAGCAAGCTAGAAAGTATTGGACAGAGCAAGACGATATTACTAGTTTCCATACATAACCATACTTACAGATAAATGGTCGAGGAATGTGTTGAACTGGGTTTCGTCTTTGTCATTCCTGTACTGAATCCATGTTGGGAGGATACATACATGACACGTAACAGCAACAGGTGCCTCTGATACTAACTTGGCTGACTCTGTAGCATCACGTGGCCTCTTTAAACCTGCCTCAAAATGGATCTCCATTCTTCCTCCTCTAGATTTAGTTAACCTATCGAGCATTATCCCTGATGTTTGTTTCCTCTTGCGCCTAGGTATATGAGGAGACTGTGAGAGTACATATAAAGAAGTGTTAGTGTACAGCAAACTATGAGACTACATATAAATAAGCATGCAACTGTAACTTGACTCCAGCAACTACCTTCTTGCTGTTGAAGCTCACAAGGTTCATCTGATATTGCCAACTCGTCTTGTGTCAAGAGTGCTTGGGAAGTAGTGTCAGGTGGCAAGGGAGCTTGGGAACGTGCTGCCAGCTCAGCTGACGCACGAGTACGTCGTGCAGCTGGTAGTACTGTGGTAGGTGTTGCAAGAACTAGGGTTGTCTCTTCTTGTTTGGTCGAAGCTATTTGTTTCTATTTGTCTTTTTTTGCAACTTGTGTAGCATGGGGTGACGTGTTTGTAGAATTTTCTGCTGGACTTTGACCTTCTATTGCACCATCAGTACCAGCAGAATCTGCAGGATTCATCCGCTTCTCATTCCCTGCCATTCTGTGTTTCTTCCTCTTTCTCTGTGCCATGTTAAGTCAAGCCAACAAGAAAAACGAATAACAGTTTAGTTCTTTAGAACAAATTGATGTGTGATGCAGACGAACTGAACTTTGATGTTAGACAACCACATGATAGGAGGATGTAATATGTATGAAAAACAGGACGGAAGAGATAACCAGAACTATGATGTGTAAACTAAGAATAAGACATTTCACCAAATTAGAAGCAATGCAATGGAAGGTAGACACCATATGAAAGATGCTACAAATATCGCTGTGCCAAGTTAACAGTGTCTCATCATAAATGGCATTTAAACAAATGTGAAGCAGTACAAGAAATAAATCATATGCAAGATGCAAACAACATCACACTACAATGTTAGAGGTCTCTCGTCATTAGTGCCATTGCATATTCACAACATTGCATATTCACAGCTTATGATGTGTAAACTAAGGAGAAGGCATTTCAACAAATTAGAAGCAATGAAAGCTAGACACCATATGAAAGATGCAACGAATATCACTCTACCAAGTTAAAACGGTGTCGTCATAAGTGCCATTTCAACAAATTAGTAGTACAAGCTAGAAAAGAATGTGAGATGTAACCTATATCACACTCCGAAGTCAAACGTGTCTTATGATAAGTGATTGTTGCAGGTTACACAACAGTAGTAATATGATGTAAGAACATGGTGTGATAGACAGATATTGTATGTTCTAGAAACAGGAACACGTGTCTTACTCAAGTATAATGTGTAAAGTAAGCAGATGGAATTCAACAAAATTAGAAGCAATACAAGCTAGACATCACACATAACATGCAACCACATATAACAGTGCCAAGTTAGAGGGAGCACCGGTGATGCTGCACTAGGGGAGACATGCTCATCATAAACACAGCTTTGTTGAGTGTCTATGCATGTGTTGTCTTGGAAATCATCTTGGGCGTGTGGAGGAGTAGAAGCTGTAGAGGCCCTCCGTTCGGAAGGAGCACCTTTATCCGCTGCCTTGCTAACTTCCTTCCGCTTTATTGATTTTGAATTGTCAAGTAAAACCGACAAGAAAAAGCTATAAAATTGTCTCTTCAGAACTCATTCATGCATGATACTGACAATCACTACTTTGATGTAAGGCAAACACATGATACCAGGATGGAATATGTACGAAAAACAGGACGACATGGCATGTTCACAACTGTTTGTGTAAACTAAGTAGAAACCATTCCAGCAAATAAGAAGCAATGCAAGCTAGACACCACATGAAAGATGCAACCAATATGACTCTGCCAAGTTAAAAGCCTCCCATCATAAATGGAATTTCAACAAATTAGAAGCAGCGCATGCTATAAATCATATAAAAGATGCAACTAATATCGCACTGCATATACTAAAGGTGTCTCTTCATGTTGATTGATGCGAGATACAAAAAAAACAATAGTTTTATGTAAGGACATGCTTAATACACAGATGTACTATGTACTAAATACAGGAAGGCATGTCACATTAACAGGTATAATGTATAAGCTAAGCAGACTAGCACATGCAGTTTAACCAGATTGGAAGAATTACAATCTAGACACCATATGCAACATGCAACCACATATCACGCTGCCAAGTTAGAGCAAGCACCTGCGATGCAAGTGTAGGGGAAATGCTGTAATCAACTACAGAGCTACGTTCACTATCTTCATGTAGCATTTCTGTTTCTTGAGAATCATGTCGGGAGGCTGACGCTGCATTCTTTAAATGAGGAGTAGTCCCCTTGCATGCCTGCCTCGTCATAAGTGATTGATGAGGGATACACAAGAACATTAGTTTGATGTAAGAACATGGTTAATACACAGATGTACTAGTATGTACCAAAAACAGAAAGGCATGTCATATTCAAAGGTATAATGTGTAAGCTAAGTAGATGCAATTTAACCAAATTGGAAGCAATACAAGCTAGACATTCGGCTGGAGTGGTGAGCATGATCATCTAGGACCCGAGAGCCTGGTGGGAGGCGGGCTGCTTCGCCACCATTATAGCTTTTTAGTTGGCTTCCAGGTGATGCCTATCTTTGCTGGGCTTGTCACTGACTCGGCCAGTGCCCTGACTAGCTATTTGTCTTCTAGTGCTCACGGGATGGTGGTTCATGTAAAAAAATATCTTTCCTCATCTTATCGACTTTGGCCTTTAGAATACAAGCTATACACCATATGGAACATGCAACCACATATGACACCGCGAAGTTAAAGCAAGCACCTGGGATGGTGGTTCATTGCAAGCGAGGTCATCAACTCTAGAGCTACGTTCCCCGTCTCCATCTGCTGAGACTGCACCCATTAAAGCATTGAAACGTTTCTTTCCTAGCCGCACAACAAACTGGAGCGACTTCTCTCTTTTCTTTTTGGCTTCTCTCATGGCAGCAGAGTCTGAAATACCCTTGCATAAAAACACAACAGTGAGAGTAAGGGTGAAGTGTGTGCAGAACTAGTGCACTGCAAAAGTAGCAGATAGCAGGTGGCTGAAAAATAAATGACAAGAGACATATGATAGCTCTACAATATTGCATGGCAAACAAAATTAGATTCTACTTCGTATGATTTTGTAATATATGAGCACAGGAAATGCAGGTACTCCTCCAGCACACCACCACATGTGGTCATATCTAAGATAACAGTCGACTGAAAATAAATAGGTCAGTCGATTTTTTTAATGAACCGTCCGATCTATAAGGAACGGGTAGATGGCCTTTGTCTACCTCCCGCCAGTCACTGTTGACGATCTTTCTCAAACCGCCAGCGACCACCGGCCCAGACCCCTGCCTCCGCCGCCTCGCCCTCCCATCGACCTCACACCACCGCCGTGCTTAGCAATGCCCCCAGGCCATCCCTTTAATCTCGACGGACCTCCAGATCGGGCGCCAGTAGCTCTAGGCTGCACACGCCCCTATGATAAACACCAAGGCGCTGCTTCCCCAGCGTAATAGGCGTACTAGCGCCTGTTACGCCGGGGAATGCTTTCGTTAATCGGGAATCAAGTGGCCAGAAATTGAAATTCAGGGGGGCGATGGACCTCTAGCTAAGGTGAAATAGTATATGAGGAGAAACCTTGTGCATCGTGAGTTACTAGGAACATCTAGTATCAAGAAGAAGCGGTCAACTAGTGTTTTCTGTGGGATGAATACCGTCCAATCAGAAATGTAAGATTTGCATGAATAAGGGAAGAGGAGTACCTTTTTCGGTTGCCGGTGTGGTCACCCCCACATGTCGCGCCGATCTTCTTGTTTTTACCGGGGAAACCGATGTTCTAGAGGGTGCATGGTTGGACGAACCCACGGCCAAATTGTTGACTGTGTTGCCGTGGCCGTATGTCGGCGGAGGGGAAGCAGATCCGAGGCGGCGAAGGACGGCATAGCAGGACAGCTCCGGTGTGGTGGAGGGTGGCGAAGTTGGAGCGGCAGCGGTGAGGCGTAGGACGGCGTGGTTGAATTGGCTCGGGTACGGACGGCTCGGCCTCGGCTTCAACTACTAGCCATGGAGGGCGTTGCGGTTGAGGTTGCGGTGGATGATGACATTGGGGTATCGGCTCCGGCGTGATGGAGGAGGGCGATGGTTGATGGGTGGGATGGGGTGGCGGACGGAGGACGCCGGGGTTTCACGGCTAGTGGAGAGGTAGTTTTGGCGCCCACAGAGTACGAATGGGGAATCGGGCGGGTGGGGGGAACCATGTTTCGGTCTGGGACGCGCTTGTTTTTGGATTTTTTGAACAGAAGATGGGACGCGCTTGTTTGAAATTTGGGGAAAGTACAAACTTTGCCTCCCTCTTAAATTTGGGACCTACTGCGGGTCGGGGGTAGGATGGTAATCCCAGGTGTCCCAAATAGTGGGCGGGAGCGATTTCGGCCGCGCGCATGTGTGCTTAGGCGCCATTGTGTATGAATGTTTTTCGGATGCAGTTGCGTGGCGTCTAGATGAAGCTCCAAACCTACCCCGGTTTAGACCTTTGGACGTCGGGCCAGTAGGTTTCACGGGTCGGAGTTATTAGAAAAGTAATTTCCTTGTACTACCAAACATTGGTCTCACGCGGTTTCGGGCGAGTAGAGGCTCTTTTGGCGCTTCGTTTGAAATTTTGAAACACAAGCATTCACGTTTAGAGTACTCTTGAACTATGAGGATTGACCTAAATATACAACGTTATATTTGACCTTGTATGTTTGAAAAACTCATTAAATTATCTGCTTCTTTTTGAAATTTTAACACTTGAAAATCCCATAACTACGATTATTTTGGAATGGAGGAGCTAAATTTGAATCTATTTTGTACACGACTACATATGCTATTATGTACAAAATTAGAGCAAAGATTGGTGCAAATTCAAATAATTGAATGGACTTCATGTTGAATTCAATTTTTTTAAACCACCTTAAGTTGAATTCAAATGTATGCTAGTTCATAGGTAGTAGCTATTTATAATGTCCATTGTGTAACTTTCGTATGTAGAATGTGCTTTACACACACAACATGTACTATACTCATGTTTATATTCGATTGCTAAAGCGATGCATTGTCTGGAGTTCAAAATACTAACTATGTGGATTAATTCTCGAATAATAACATAGACATGCTCAAATTCGATAGAATCTAGATGTCTGATGGATCAATGACCGCTCCTTACGCGAATTGCAAATATCATGATCCCATCCTAATATTTGGATACAATTTATATCTTTAATCAATGAGTAGAGCAGTCTTCTACATGTAGTTTCACTCTCCATCGGTGGTCTCTCACACATGCTCACACACTCTCTCCCGAGGTGTCGTTCTCGCACACATGCTCTAGATATAACCCTCCCTCTCTCTCGTAACCATTTACAACTGTTTTCACTAACCTTTATCACAAGCACTCTTCCTCTCGCAAACATACACACACACAATCCTCACCAACCCTTTCTATATACATGGACCCGCCTACGTGTCTCGTGTACGCACATTATCCTTATGCCTGTCTCTCATGCATTGTCTCACACCTCTCTTCCTAATCGACGAGTATGATTGTAGCAGAGCATTCTGTAGGATGCCCCTTAAAGTTAGGTTGGATGAATAGTAATATCAGAGATAAAGGGGTTACCTTAGTTTGAGAACCTAGGGTTACAAATCCTAGCTAGCTTTGTCAGTGGACACAGAGTCCTTCACAATTCTGCTATCTTTGTCTTGTACGACTAGTAGTCCATGGGTCTTCCTAAATGCGTCACGGGCCAACCAATGTGGCGCAGGACGGAACAGAATGAAGGGCCTCACCTCGACCCACCGGACCTCATGCGGTGACACACCCTCTGCCCCCACGTCTCATTATCTTCTCACACCCACAAATACCAACACAAACACCACACACATAGAAAATTTCTCCTGGTGCCTCTATTCCCTCCCCCCTTGCATATACACACTCAACAGAATGGAATCTCTTGTCGTGACATCATATTCGTCTCTCTCTCTCTCTAGACCACACTCATTTCTCGTGCTATCTCTCCCACGCCCCCCCCGTTGGCCCCTCTATCCATGCCTCTCTCGAATACGTAATACACACACACATACCTCTCTAGGCCCCGCCAAATGCATCAACTACACATTCACACATGTTACATCACGTACCCCATTGATCTAAAAAGCCCTCTCTTCATGTTTATTTCGCATACAACATTCCTTTTGAATAAAGTGTGGCCAAAAAAATATTTTAGAGTTTCTACATATATACAACATTACATAGTTATATGGAGGAGTACGAACGCTAATTTGCATTGCATGGTGCCCACACTGATATTTATTCCGCACTGTGAATTTGTATATTTTTATACTATATACAAAGTTAGTCATAATAAAGTGAAACGAAGAGCATCTCTCTCTCCATCGCATACCCCCCCTGTACGCCCTTTCATGTATGCATACAAAACTATCTCCAGGTCCTCTAAAAGTAAGCCCCCCAGAAAATTCAGGCATGTTTCATCTTTCTCTCGCGATATATGCAATGACGATCATTGTATTTGAAGCGAAAGCATGATACACACATTGTACTTCAAAATCCACTCATTACGGTCACCATTTATGTTTAGTGGTTATTATGCAGTCACGGGCTCACCGTACAACTCTATATTTGCTCATGACATGACTAGTTGACTAGTTAAATGCTCCATGTGGCACCTCTAGTGTTGCATCACATTATCTCACTACCAACATGCACACCTGTCGAGTTGTATCTACTCTACATCTACACTCTTATAAAATACATAAAATACACTCTTATAAAAACAGAGTTGGTGATGATGGTGTGCCTGCCATCCTGCAATATAGGCCGTCCGATTTATATCTGACGAATAGGAAAGAAACTATGGCAATTTTGCAAAAAGGTACCCACACACCTCTCCACATTTGCAAATAAGGCCTTCCCTCAGTCATCCTTTTCTCTCACAAGATAAACAAGTCATACAAATGCATCTTGATGTTACGTGCAACGCACGGGCATCTTGCTAGTAAGAATGAAAACAAACGAGAAAAAATATATTTGGACGGTGGTTCAAAGTCATGACCGCAGGCACAGCGGACAAGGTGGCCTTGTATTGGTATGCTGAGCCATCTGTTTTAACACACAGGAGCTGGTGTGGTGATTATACACACACACACGCGCATGCACACGAACTGCATCCACGCAACACTCACACAAACACATGCACCCACGCACGCACGCAACACTCACACGTACACACGCAGCCACGCACGCACGCAACACTCACACGCACACACGCACGCATGCATGCACACACCAGTACACCACACGACCAACACACACTCTCACGTTCAAGTGCTCAACCTACACGTGCATGCGCACACATCAGGCCCTGACAGACACATACACAACTAGGTGATTCCTGCGCACGGGTTGGAGCTTTGTCGCGCCTGCGTAAATTTGTGTTTATCTGACCTTTTGCCTATGCGAACTAACAGGTGGGACCTACAACGAACGTGGTCAATTTAACTAGTAAACATGGCGCCACACAGACTATTTGGCCGGTCAACAGAGTGCAATCCGATGCTGAGCTGTGAGCAAGTGACTGTATAATCACCGCAAAAGGTATATAGTGACCGTAATCAGCTTATTCTGAAGTTCGGTGACCGTATTACACTTTTCTCTCTGTAGGCCCCCCAAAACTACATCTCTACACATTCTCGCATGTTACATCTAGCAACTATGCAATACAGCACTACTACTACACTCTAACACAATAAATCATATGTGTTTTTAGTTTTACTAGATATCATATTGTTACTTATAATTATTCAGTTTGCATACATTAAAATTGCCTTTGAAATGTATGGCCAGTATCATCTACCGCACGGGAAAAATCCCGCCCTTTCCTGTTTAATCGAGTTTGTATCGATGGATCGTGCGAAACAACAAAACTATACTAGTACTATAAAGTACAAGGGACAGTTCACTGGGCCGTCGTGTCATCTACTCCTCCGTCGTAAGAGTGGAGCGCTCGCTTTCATAAATATTCTAGTCCCTTTCACCGACATGTGGGACATCAAGCTACCGGGTCCATGTGCCATACCAGAATACAGTGCAAAGCAGCTGGGGTGAAGCACCCCAGTCATGGCGCCGGCGTAGTAATGAGCAATGAGGCATCAAATCACAAAGCTGCACCTTTCCCGCAGTTAAGTTGGAGGGACGAAGCAACCATCATTTCACTCGTTGCTGAATCCGCTTCTCCCTCATCCTCTTCAACTTAACCACGTGTTGCTTCTGCAGTTGTTCTGCGTCATGTGGGACGCGGATCGGCTTGGGAAAGCACTGACACGTGGGACCGCATGTTAGCAAACTGCCAGGACAATGTCCTCCGAAAATAAGAAGCCTAATCCTAGACTTACAAAGGGCTACGTAGCTGCTACGTATACTTTCCATCTAATCTATATATGCCCCCCTTGATTTTCAGGTGGGGTGGGCCTAATCCTCTCCTTTAATCCAATCAAATAGTCCCACATCACATCACTTTGTAACTTTACGTAAACCATTACGTGGATGTAGCATTGCTAAAATAAGAAACCTCCCCCGCCATCCGCGCAGCAAAATCACCTCCTTTTTACTAATCTCGATTCTATAATTTTTTAGATCAATATGATTGAGATTTGTATAGATAGCACACAAGAGCAAAACCAAGTGGTACGGTTAAGGGCATCCACAATGTGTAGCCAAATGTGAAAATAGCTTTTGGCGTCGTGGGAAGCATTGTGGATGGTGATGCCAAAGCGAAAAAGATGGAACCGAAGCTCCATGATTGATTGATTGAAGGAATAGAAAAATGCAACGTCTCTCCCCTTTCTCCCATGCTTGGACGCATGTAGTACTAGAGTATAGAGCACAGAGTCTACTCCCCTGTCCCTTCTATCACAAATCACAAAGAAGTGAGGCGTCAAATCACAAAAGCATGGACGAAGCAACCATCATTTCACTGTTCGCTGACTCCGCTTCTCCATCGTCTTCTTCAAGAAGCCATCTCACCGCATTAGTTACTCCTAGTAGTACTAGTAAAGGACTTCCTGGATGCAAATGGTGTTCTCCGTCGAACGCACTTCGCCAAACCATCAAAAAAAATATCGTCGACGTCATCGCAATGGGGAAGGAAGTCAAATATAGTTACTACCTCCATTTTTTTGTACACAAGGCCAGTATATCAAAATTACAATTTGGAAAGGGCCACTAACACTAATCGAAGCAAAATCGATGGGGTTAGCGACCAAGGACATTAATATCCCCTGCATGCATGCGAAAGTGAGAAGGTTGGTTTACATGGCGCTGCATTAATCAAGCGATGAGGAGTGAGATGGTTTGGTCTTTGGAGAAGAAAATACGCATTAATTGACACGTGAAACATGGATTTGTATCGATTAAATCCCGCGAAGGAAAACCCAGCCTTTCTTTTTTAACACTAGTACTCCTAGTATGTACTAAAAAGAAACAGAGGGAGTACGTACTTATCCTGTTTGGGTCTAGGCCTTGCATTGCCAAACTTCGCCACACATTTTTGCCAAGCTTGCCTAAAGTTTTCAAAATGAGAAGGAGGTACACTTGCGCACGGCTATCGCGGTCGCACCGCACCCTCACACGGTCGCTCCTGCACGCCGCGATCGCACCGCACCCTCACACGGTCGCTCCTGCACGCCGCAAACCCAACCAAGGGAACGCACCCTCACTAACACGTGCTATCGCATGTGAACACACCAAACTCAGCCATCCTATCGCTGACACATAATTTTCGTTCATAGAGTATGTTTATCCAACCGTATGTTGGGGAGTATCTGTACGGCCCGTGCGGTGGTCTGTTGGGGAAGAAGAAGAGGTGAGGAAGAAGGAAAGAAGGGTTAGGAGACGTAGGATATATGTTTTTACACAGCAAAAGATCCATGAAAACAACAACATAAAGAAAAACTATCACTGCTAGCTGAAAGAGACGGCCTCGTCGTCTTTGGCTGCCGGTGCGAACCAGCCGTCGCTGCCGACGTGTGTCTCCGCACCGTGGTTGTCCTCGGCCTCCGGCTACTCGTTCAAGCGGAGCGTGCAGGCGCTCTGCACGTATAAGAGCACAGCCCGGTTCAAGTCGAGGATGTCATGTTCCGCCTGCAGGAGGGCCTCGATGGCAGCGGCATCCGCGCGCTCTGCGTCGAGTCAAGCCTCCGCCGCCTGGTTCAAGTCCAGGATGTCTGGTTCCGCCTGCAGGAGGGCGTCGATGAGAGCGGCATCCGCGCGCTCTGCGTTGAGTCGAGCCTCTGCCGCCCGGTTCAAGTCCAGGACGTCCGGTTCCGCCTGCAGGAGGGCCTTAATGGCAGCGCCATCCGCGCGCTCTGCCTCAAGTTGAGCCTCTGCCGCCCGGTTCACATCCACGACATCCGGTTCCGCCTGCAGGAGAGCCTCAATGGCATCGGCATGCGCGTGCTCTGCGTCGAGTTGAGCATCTGCCTCCCGGTCCTCCTTTAGTATCGCCATGTTGAACCGCTGGTCCACCTGTACCATGCGGGGCTCGGACGCCGGCACGAAGTCGACGTCCATTATCTCATAACCATCGGGGGCCTTCTGCGCCGCGACCTCTGCCATGAACATTGGATCGGCTTCCTCCTCCTCGTAGCTTGGATCCAGAGAACTCGGGGATTGGCCCGCCGCAAAGCGAGCTTGGGAACGGAGGTCTGTGGCGCCGACCGCCACCGCGATTTCTGCGTGGCGCTCCGGCGTCAGCATCTTGTAGTAAGTGATCTTGCAGCGCACGATGGCTTTCCGGCGAACTAGGGGACGCTTTATGCGGGCTAGGGGATAGAAAGGTGGGGGAATGGGCTGGAGATTTGGTGTGGTTTTAGGGAAGTGGCTGGCGTAGGCCGAGCAGCGAAGGTTCTGGTGTGAATAGTGGTTCCGGTGACGACCGGTCAATCGGTTTTGACTGTACATCGCTCTTGTCATGTTCGCGTTGCAGCACGGCACGTTGGACCCTCCACGTCGCCCACCGAATGGAACGCGCTTCGTCTTGCGTGTTCGCTTGCTTGTGGGACTGCAGTTGAGAGCAAACCGACGTCCCTCCCCCTCCCCCGCCCGCATCATTTATTATACCACCACTCCCTTCGTTTTATGCAGAGTAAATAAAACCAGAGGGAGTATTACTACGGATGAGGATCCCCTGACGTGGCCGAACACATTGAGTAACTGACATGCGGGGCCTAGCAAGCATGGGGTCCGCCAGTAAGTCACCGAAAGGGAGGGTAAGGCAGGGATACCTACGTCAGAGGATCCACTTCCACTATACTACTTTGACCAAATATTAGAGTAATAATATATGGCATGCAACTTACATAAATGTTAGAGTAATAATATATTACATGCAACTTACACAAAGCATACAGGGTCTAATGTGTTTTTCGAGGCTAACTTTGACCAAATGTTAGAGTAATAATATATGACATGCAACTTACACAAAGCATATGGTCAAATTTGTATGTGAAAGGAGTTTCCAATGATATAATTTTCACATTATACATCTCATGTACTATTAATCTTGTCAATAGTCAAATTGGAAACCATCTCGAGTACAAATCTAGGAGTACGTACGAATCTAACAATATTGATTGGGAGTTGTTGAATGTTGATACCTTTTTTTATCAAAGTAGATATACTTTGACTACAGATAAAACTAATATGTAAACTACTCCCTCCTTCCATCTGTATAAGGCCTAATGCGTTTTCGAGGCTAACTTTGACCAAAAGGTAGAGCAATAATATATGACATGCAACTTAAATAAAGCATACCATCAAATTCATATGTGAAAGGAGCTTCCAATGATATAATTTTTGTATTATACATGTCATGTACTATTAATCTTGTCAATAGTCAAAGGCGGTCTTGAAAAATGCATTAGGCCCAATATAGATGGAAGGAGGGAGTAGAAAGGATAGGAGGGAGTATATTGTACGTTCAAAAATGAAACGAACATTGAATACCCTTTATATATGATTTGCGGATGCTATTATCACCGGTTCAAAATTTTGAATCCGCCTTAGGTATCACGTGCTCCCACTCATTCCTTCTTGATTTCATCTCGGGGTCCAGAGATCTATTGAAAGAGGACTAGGGTGGGTACATGCGAATGATACTCGGCGGAATGTTCAAATGAGGCATGCGGGAGGATGGACTCGACTGGAGGGCGACATGATTGATGGAGAAGAGGGTTGGAGAGGAATGAGAAATAAGCAAATGCCACATGATTATACTGGAACGGCTCAAATAAACAATAAGTTGTCTCGCTTGGCCTACATAGTAAAAGGCCTAATGCGCAACGCGAAGCACTAATGCTCGAATATGGCCGGGAAAACAGGGAAACAGACATGTGGGGCACACCCGTCGGGACGACGTAGCACAGACTAGTCCAATGGAGGAGCTGGCCCACGACCTCCTCGCCACTGACAACTGTTCATGTGGGGCGTTGGGGCCAACAACGGGGCGCAGCTCCAGCACCGCCTCTAGACACGGCGACCGCGGTGAGCGACACGCTCATATCGTCGCTAGACACACTCGGTTTCAGTTTGCCAAGGGTTGCGTTAGTGCTGTGGCACGACCAACGGTATGTTTTGTTTGTGCCCAAGGTTGCCCCATCAAAGCATTGGGTAGCCAAAATTTTGGTTGAGGTATTGGTTGCCCATGATTTGGCCGACATTGGCAAGAAAAATGAACTTGAGTTGGCTAGAGTTCATTGGCATGCCAAAGAAATGGCAACCATCCAAACAAAGACCAATCTTTGGGTCATGACCAAATTTTTGGTTGAGGTATTGGTTGCCCATGATTTGGCCGACATTAGCAAGAAAAATGAACTAGAGTTGGCTAGAGTTCATTGGCATGCCAAAAAAATGGCAACCATCCAAACAAAGACCAATCTTTTGGTCATGACCAAAATTTTGGCAAGGTGCACTTTGGCCACAATCCAAACACACCCCAACACCCGACTCATCGAGGATGCACGGGGCACCGCGACATGTCCGCGGAGTGCTGTAGCGCGGCCAACTGCATCCTCTGGACCTGAGACCATGCCGGGTCATCTTGCTTTTTTAATGACATGCGGGGATCGCATGTGAGCAAAGGGCATCTCCAACGTTGCCACGACCGCAGCTGACTACCCTGGCACACCAAAACCCTAGTCCCTCGTGCCATCGCGCGGTGCGTTCCCAGCCAGCGCCAGCTGCCTGCATTCATGCCGTCCGTTCGTATGTCGCCATTAAGAGGCTCGGTGCCCGAGGAACAACTCCGGCGACTTAATGACACACCCGGACGCTTGGCCTCACCGGAATGCGCTACTTAAAGCGGCCACGCCAAGCTAACCTCCACACCATAGCACATCGTCATCCTACCTGGAACCACATCCGCCATGAACGCAAGCGTGAACGCTCTGCAGGAGAGCTTGTCTTCGGGGATGGAGCTCGAGGTCGCCGCCCTCGCTGAAGTCGCCGCACAAGCGGTAGCCACGCTGGCGGCCGATGACGGGAGCACCGTCAGCCACGCGTCCTACTTGCCGCCGTCCAACGTCCGGCACCGCCGAGGTCGGGGACTCCTCCAAGGATGAGTAGGGCATGGGAGGCGGCAGGGCATGGTGTGCCAACTGGCCGGTCCGTCTCCCGCATGTTCTACTTTTCCTCGCCGGAGACCGCACCTTCACTTCACGGGTCTTGAACCCCGCCCCCGTACATGGAGATCAAAGATGGAGGACGTCGGTCGCCGCCGTAGAATAGGTTTAGGGTCGAGTTTATTTTATTTTGTTGTCCTATAAAAGTTTGAAATGTAATGAAAATCTGCCGTGTTTGCGTTAATCTCGGCCAGTGTATATGAACTTTCACAATAATATACATATTGTAGGAGTACTAGCAAGATGCCCGTGCATTGCACGGAAGATCAAGATCTTGTGGGAGAAAAGGATGAACGAGGGAAAGCCTTATCTGCAAATGTGGAGATGAGTGCGGGTAAATTGTCATAGTTTCCTTCCTATCCGTTAGATATAGATCGGACGGCCTATATTGCAGGATCAGAAACATCTTTCTACATCTACGGGAACCTACGAAAGGGTTAACGTAAATTGCCTCTCTCTCTCCTCCCCTGATTTTCATTGTGATGGGCCCCTCCCTCCCCCCAATCTACAATCAACAGCTCACACATTTCACATTATGTTAATTACGTAATTGGGATTACGTAGGTGTAGCATTACTCGTCCTAAAAAACCCAACTAAGCCAACCTCCCGGCATATCACGCAGTTAAAGACCCGGCCGCGTCGTTTTAGTCGGGATTACCGGATTACTAGTAGTAGTATCGATGGATCGCGCGAAGAAACAGATATTTTTTTTATGGGGGCAAAGCACCTAGTTTAAAAGGAAAAAATGCGGTACACTCACAGCACTCCTGTCGCGGTCGCATTGCACCCCACACACGTTCGGTCCTGCACACGGCTCACGCAAACGGAACGAGCTTCGGCTTGCCTCTTCACTGACACGTGGGACTGCACTTGGAGAAACCGACCATGCACCAAACTCCCCCCGCCCACCGCGCGAAAATCCTCCCCCCCACACCCATAAATTCCCCCCAAATCCGATTCAACCGCCCTTCGGCCAACTAGCCAATAATATTCCCCAACCCCCCCCCCCCCCTGTCCCCCTCCACTCCATTCACTCTGCCAAAGCCGCCCTCCCCGCTGCGCCGATACGTCGGCACCATGGTTCGGCGATCCTCTCGCTCCCACAGTACGCTCTCGTAGACTGCCGTCCTGTAGCCGCGCCGCCACCTCTGGCTACGTTCGCGTGGAGACGCACGACAGACAGCTTCTCACACGGTACGCGCATTCTAGACTGAGGCCCCGTTCATGTCGGTGTAGCAATGAATGTTAGCTGATGGACGCTGTTAATCTCACTGTTTGCCAAAGTAGTTTACTGTCAATATGCTCTTCTTAAGAAGCTTCCTTGCCCTGTTCTGCACTCAATCTGCTTAGCTGAGATGGCATGCCAAGGCCGATTAGTTGCTAAAATATCCCGCCATATATTTATTGTGACGTAGATTGCCATGGTCTAGTAGCCAATCTGTTAGGTGAAATAGCTCTACACCGTTTTATACTCGATCTTTGTTGCTGCGATGGCCTTCGATAGTTCGATCGCTGTGTGCTCGGCCTCATTGATTGCTAAAACAGCCTGCCATAATTTGGCACTCAAATCTGTTTGCTGAGTTAGCTTTACATATATTTCGTTACTTAGATCTGCTCGTCCAAATATCTTGCCATAGCTTTAGACATCGACCTAGGTATTTTTTCTGGACTTCATAAAAAAGCATATATGTTTTTCCTATAATATCTAGTAGAAGAACTAATTTGTATGTCTAACAGATGGACAGGACTTGGATAACTTCTGCTCAAAGATTCTCCGCTGCATATGTCGAGGGGGTTGAAAACTTCATGAAATTTATCAGAGCTGAGTACGGTGGTCCGAAATCAGATGTGCTCTGCCCATGTAGCAGTTGTATGAATTCAGTTACAAGACCCCAGTCAACTGTGCAAAATCATCTACACTTGTATGGGATGTCGGTCACATATTCTAGGTGGGTTCATCATGGTGACGCTGTGAACGTCAATGTTATTGACTACGTGGAAGCAGCAGATCACCATCTTGGTCTGCGTGATGATCAGGTGGAAGAGGAGGAGGAGGTGGTGGTGGCGGAGCCAGTGAGTTTGACCAACATTGAAACAATGCTATGAAATGCTCGTGCATTCCGTGAACTTTCACCTGCAGAAGAAAAACGGTGGGCCCGCATGTTGGAACAATGCAACATTGATGTCACCCCAGGAAATAAGCTGTCAGTATTCTCAGCTATGGTCACCTTTCTTCAGGTGAAGACATCTGAGCGCATGACCAACAATCGATGCGATGTTGGCTGCTTTCCGCGAATCTTTCCCAGATGCGTCTGAGCTGCCACACACCTACAGTAAAATGAAGAATTTCCTTCGTGCATTTGGAATTGGATATGATATGATCCATGTTTGTAAGAATAATTGTGTTCTGTTCCGGAAGGATTATGCCAACTTAAGTGAATGCCCGAAATGCAAATCATCAAGATGGAAAGATGGCGATGCTGTGAAGAGGATTCCTCATAATGTTCTCAAACATTTTCCAATTACACCAAGATTGCAGAGGTTGTTTCATGATGCTGAAACAAGAGAGGATGTACTGTGGCATTCTAGGAACCAGGAGTATAGAGATCAAAATGTAATGAGCCATCCATCTCATGGTAGTGAGTGGAAAAGTATCAATGATAAACACAAAGAGTTTGCTGCTGACCCAAGAAACATTAGACTTGGCTTAGCTTCAGATGGATTTAACCCATTCGGCCACCAGAGCGCCACATATAGCATGTGGCTAGTGCTTGCTATCCCTTACAACATGCCTCCAAATGTCCGCACCAAAGAATCAAACTACATGACAGCCTTGCTCATCCCAAGTCCAAAAAGTCCTGGAAAGGATTTTGATTTGTTCATGGAGCCTCTTGTGGAGGAACTTTAACAGCTATGGAAGGGTGTTCTCACTCGAGACCTATATAGCAGCCCACCAGCTGATTTCTTTCTGCGTGCTGTTATAATTTGGTGCATCCATGATTATCCGGCTTTGGGCACTATGTCAGGGCGAACGACACATGGTTACAATCCATGTGTTCGCTGTGATAGGAATCCGCTGTCATATGCAATACTTAGCAAGATCTGTTACATTAGACACCGCTGTTTCCTTGCCAAGGACAAGCCGCTTCCTAGAAAATACCAAATACATGTGTTCAATGCAAAGCATGAAAACCGTGATGCGCCAAAGAGGCTCACCGCCGATGAGTTGCAAGTGGAATTAGAGAAGGTCAGGCATATTACACTAGGAAACCATCCTGATAATGGTAGCGGGAAAAGGAAGCGTGGCCGGGCAGAAGAGAGATTATTGTTTAACCGCAGGTCCACTTTGTGGGACTTGGAGTATTGGAAAGATTTGGATCTGCGGCATAATCTTGATGTGATGCACATCGAGAAAAATATATGTGATTGCATTATCGGCACACTTCTTAATATTGAAGGCAAGACGAAAGATACCTTAAAATCTAGGATTGATTTGACACACCTGGGTATCAGACAGGATTTGCTGGTGCAAGATGAAGGCAAACCACGGGATATGGCACCAGCTGTGTACGTCTTGGACAAGGTAAAAAGAAAAGAATTCTGTGAGGTCCTGTCACGTGTGAGATTCCCACATGGATTTGCTTCCAACCCTGAAAGGAGAGTCAGTGCAGATGGAAACAAGGTACAAGGGTTGAAAACTCATGACTGCCACGTCCTACTTGAAAGGGTTTTACCTGCTATCCTTAGAGGATTGGGCTGCCCTGAGTTATACAGAGCAGCTGCAGAGTTGGGACAATTCTTTAGGGAACTCTGCAGTAGAAATATCAGGATAGATGCTTTGGAGCGTGTTAGAGACAAGATACAACCTATCCTATGCGACCTTGAGAAGATATATCCTCCAACCTTTTTGATGTGATGGTGCATTTGGCTATTCATCTACCTGATGAGGCACTACCTAGAGGTCCATTACAGTATGGCTGGATGTACCCTATTGAAAGGCGGCTAGGCACTTTCAAGGGCTATGTTAGGAACAGAGCTAGACCCGAGGGTTCCATTGCAGAGGCCTACATTGCTACAGAAGCATTGATATTCATCTCAAAATACATTGAAACAGCTGATCAGGTTAGCAAAGAGGTGGGTGAAGACAATCCCGGGCTCAATGTTTTCGATTATTCTGTTCGAGTTACAGGGAAGAGTCAACAAGAGGACAAACCTAAAGATTTGGACAAAATGGGTTGGTATGTGTTGAATAACTGTCCTAAGATACTACCTTATATCAAGTAAGTGCAGTACTGCAGCTTATACATCGTAATCTATCAAACTTGCAGCACATTCTTATATATTTAGATGTTTGACGTGATCACTATGTTGTGCAGCATCTACAAAGAGGAGTTACTGCCGCAAAATCCAAGAAACATTGACAAACTGGTTATGGCAGGATTTGCGAAATGGTTCAAGAACCATGTAAGCTTTTGAAGTGCACTGTCAATTTTTTTAATATATTGAGTACATAAGATGATCAGCTTGTCTATTTTCTAACCATAGGTTAAGAAGATGCGGGAGGATGGGCAGGCAGTTGATGATGCCCTTTACTCACTAGCAATGGGTCCTGATACTCGGGTAAGACATTATGAATCTTGCGTTGTTGGAGATGTGCGCTACAACACCCTTGCACAAGACGAAGGCAGGAAGACACAAAATAGTGCCATCATGAGCACAGATACGTATGACAAAGAGACAACTGAAATGTATGCTAACATAACAGACATTGTTCAGTTGCAGTATATCTCCAGTTTCGAGGATCATCGGTGTGTGGTTCTGTTGTGTTGTCGTTGGTATAACCTGTTTTCCAGGATCGCAAAACCCAGAGCTGATGATTATTTCAAATCCATCAATGTCAAGGCGGCGTACCAGACAATGATCCTTTTATTTTGGCAAATCAAGCAACACAGATATTTTTCTTGGAAGACACATTTGCACGTAGCAATGACTGGAGAGTATTGCAAAGGTTTGAGCAGAGGAATTCGTTTAATGAAGTTGCACAACAAGATGATGCTTACACTTCTCCTGATGTACAAGATAACACAGATGTTCCTAATATCTTTGAGAACTATCATGTCAATGACGCCGGTGAAAACATTTCCTGTCGTGCTGTGGACATACAAGAGTTGATCAAGAAGAAGCCAACGTTCGAGGACGTTGAGGACGAAGAAGAAGATGACACCGTGGGGAATTACGATTCAGACTAATACACATGGCGAAGATGTTGACGCTGCTGCTGCGGATGATGATTACATTGCTTTTGTTGTATTTGCGTGTCAGAAGACGTTTTTTTATCTACTATGTGAAATTGTGTTTGCTGTGACAACTGAAACTTGATGAACATGTTGTTTAGTATTTTGCTGTGAGAACATCACTTGTTATGTATGCTGATTTTTGTTTGGTGATTGTATGACAACTAAAACATGATGACATTGTTGTTTAGTTGTAATCATATTTGGCTGTGAAACTTGAATTTGATGACATAGTTTGCTGTTGGACTGCAATTTGCTCGACGTCTTCGAATGAGAACAAAGCTGACCCGACAGGGCCTCTGTTATTTATTTATTTATATGAGCAAAAGGGGATACCCCCTGATTTCCGTTAATAGAAATCATTAGATGTTCACAACACTACGCCCAACCTGCTACACAGGTTTCATTCATTAGTAGTTTCAGCAAATTGCTCATGTCGTAGCCACTGAATACATCACGAGTGCTACATACACTACAAATAAAGAGACAATCATCCTACAGAGCACATCGATTTTGCCCATTATCTTCACAAGCTCTTTCATCTCCTCATTGCTGTCAAGGCCGTTCAGCAGCTGTTGCTTGGCCATGACCAGAGGAACTGCATCTTTAACCTTCTGCTCTGCATCTTCCTCTTCTGTCGTTTCTTTAGTTACTTGTCCAACACCCCAACCTGCTGGTATTGGGATTCCTCGGCTAACCAAATAATCAGCGTAGTTGCATTCCCCCACATCAGCAACTTCCCAGAAATACAAATCACACGAATCTTCCCTTTTCTGCACGAATCAAATTGGTAGATCATCAACCAAACTATTAAAAAAATCAGCAACTGAAAGCAGCAGAACAACACTTAAACATATTATTACCCCATGATTCGGGCATTTGTAGAAGACTCAGCCAGGGTTGCGGATTGTGGTTGAGACGCGACGGATGACTCTGCGTGATTTGCAGTAGTGGCACCACACCAACGGCAGCAGGCTGCGATGAGCAATCGGCTGCGGTGCGCGTGCCGGCGAGGGTGTGATGAGACCCATAGGTGATGGTACGGGTTGCGACTTGGCGCATATCTGGTTGTTGCCTGCTAAAGAGCTGGTGGACGAGTAGCCAGCGGACATGGTTGCTGCGACGAGAGCTTCTGATGCTAGGCAGAGTAAGTAGAGAAGAGGAGAAGCGAACGTTGGATATGTCAGTTTCATTACTTGCAGAGTAGCATAGCACCATATATAGTCATAGTCTAGGTTTGGATTCGAACCAGCTACAGGAGCACGTCTTTATTTTTTCTAAAACTTGGCAACGTCTCTTTAGTGGTACACTTGCTTGTTCTGTACGCCTTCCCATGTCTTTGATTTTCTTTCCCTTTTTCTTTGCGAGGACCTGATAGCTGGGACCCACCGGAAGGGCCTCTGTATTTCACGAAAAAACATTCCCCCCGCTGACATGTCAGACCCACCAGCTATATCTTCGCACGCAAGGAAGTGCCTCCTTATTACGCACAAAAAAATGAATACCCCCCTGCTAGCTGGGACCCACCATAGTGGGAGGATGACTTGTGGGCCAACAAAGTTGACGGGGACAGAGGGCTTTGTCAACTTAGTCAATATGAACGATTCTAGCTCCAGTGACCGTACGATGTCCATCCAACGGGCGTAGTGCTTCTTCAACCTCTGGTCTTCTTGCAACAGCCGCCCAAATCAGCGCCGGTCGTGCCGCATGCTCCTGCCTCCTATGGCCGGCTGTGCTGCTGTAGAGGCCTCACCGCCCCTACTATTCCCACCGCTGGCCAGGCCCTGCGGCGACGGCAGCCGCATCCGAACCAGTGAACCCTCGTACTCCTCTCCGCACGGGCTTCCACTGCCGCGTCTTCCCCGGCTCCGCGTCGTCCCCTTCCTAGGCCTCGCCGTCGTCCACCGCCGTGGTTCTCTCCGTGCGGCGTGGTGAACGTGGTCAAGGAACGACTTTCATCGGAAGAGTACTGTATGTGGAGAGGCTGACAGCTGGGTCCACGGCCATAGCCCAGTTTTTTGTGATTTGCCAAGTAAGTCGCTTTGTCAGGCCTATTGGGCTGCAAATCTTCCAAGACGAGGAGAGCTTTCATTCAGGTGGCCGAGAAAATGGCCCATCAGTAATGAGAAATGGGCTGTACATTTTTAAAACACATCAAACCGGCAATTAGTTTCAAATATCTTTTTTTTCATTTCAAGGTTTTAAATTACATTAATTTTTATGCGTGGAGAATTTGTTGGATTTTATATTGATATACATTTATTTTTAAAATCAGTTTGAATGTGAGTCGAAATTTCGGGATTAAAAACAGTTCGGACCGCACCGAAATATGCAAAATTTTGTATAACTTTTTAACCGTGGCCACAATTTGGGCTATAATGCTAACAAAAAGACTATGGGCTCCAAAAAAAATCTTAAGAATTAATAAATGGGCTGAAAATTATTAGAAATAATGGCAGATGGGCTACATGTTGTTTTCCACAGATTTGAGGCTTTCCTAAAAAAAGGTTGACACACAAGCAGTGACTGTTGGATGTCCATCCAACGGCCGTCATGCTTCTTCAATCTCTTCTCTTCCTGCTCCAGCCGCTCAAACAGGCACCGGCGGGACTGCCTGCTCCCTCCTCCCCGCGGCCGGCTCTGCTGCCGCGCGGGCCTCACCGCCCCACTCTACTCCCATCGCTGGCCTAGCCATCCCTCTACTCACCCACACCTGCTGTTATTCTCCAGCGACGGCAAACGAACTAGTAAACCCTCGTACAGTCGTACTCCCCTCCGCATGGGAAACAACTATCGAGTTTTCCCTGCCTTTGTGTCATTCCCTTCCTAGGCCTCGTCGTCGTCCACCGCTCTGGTGCTCTCGGCGCGGCCTGGTCAACGTGGTCAACGACCGACATGCATCTGAAGTGGACTGTACGAGGAGAGGCCGACAGCTGGGTCCACGGCCGCACGGAAGGAAATGCCTCCTTATTACGTGCAAAATAATGATTCCTCCACCTGACATTAGGGACCCACCGAAAGGGCCTCTGTATTTCGCGAAAAAAACGTTACCGCCGCTGACAGCTCGGACCCACCAGCTGTATCTTGGCACGCAAGGAAGTGCCTCCTTATTACGCACAAAAAAATGAATACTCCCCCTGTTAGCTGGGACCCAGTATACTGGCAGGCTAACTTGTGGGCCTACTAAGTTGACGGGGACGGAGGGCTTTGTCAACTTAGTCAATAATGCATGATTCTAGCTCTAGTGACCGTATGATGTCCATCCAACGGTCGTAATGCTTCTTCAACCTCTGGTCTTCTTGCTGTAGCCGCCCAAACCAGCGCCGGTCGTGCCTCGTGCTCCTACCTCCCGTGGCCGGCTGCGATGCGGCGGAGGCCTCACTGCCCCCTACTACTCCCACCGCTGGCCAGGCCATCCCTCTACTCACCCACACCCCCTGTTACTCTGCGGCGACGGCAGCCTCACACTGCAGCGAACCAGTGCCCTCGTACTCCTCTATTTGTGGGCATCCACTGCCGCGTCTTCCCCGGCTCTGCTTCGTCCCCTTCCCAGGCCTCGCCGTCGTCCACCACCCTGGTTCGCTCGGCGTGGCGTGGTCAACATGGTCAAGGAACGACTTCCATCGGACGTGGACTGTACGTGGAGAGGCTGACAGCTGGGTCCACGGCCGCAGCAAGGAAGTGCCTCCATATTACATGCAAAATAATTATTCCTCCACCTGACAGCAGGGACCCACCGGACGGGCCACCAGTATTTTGCAAAAAAACCGTTTCCCCCTGACTGCTGGGACCCACCAGACGGGCCACCGTATTTCATGAAAAACACGTTCCCCCCGCTGTCAGCTTGGACCCACCGGAAGTGCCTCCTTATTACGCAGAAAAAAATGAATACTCCTCCGGCTAGCTGGGATCCACCTTGGTGGGAGGCTGACTTGTGGGCCTACTAAGTTGACGGGGATGAAGGGCTTTGTCAACTTAGTCAATATGAACGATTCTAGCTCCATTGACCGTACGATGTCCATCCAACGGTCGTAGTGCTTCTTCAACCTCTGGTGTTCTTGCTCCAGCCGCCCAAAGCAGCGCCGGTCATGCCGCATGCTCCTGCCTCCCGTGGCCGGCCGTGCTGCCGTGGAGGCCTCACCGCCCCTACTATTCCCACCGCTGGCCAGGCCCTGCGGCATTGACAGCCTCACACCGCAGCCGAACCAGTGAATCCCTCGTACTCCTCTCCGCGCACGCGCAGGCTTCCACTGCCGCGTCTTCCCCGGCTCCCCGTCGTCCCCTTCCTAGGCCTCACCGTCGTCCACCGCCCTGGTGCTCTCGGTGCAGCGTTGTCAATGTTGTCAAGGAACGACTTCCATCGGACATGGACTATACGTGGAGAGGATGACAGCTGGGTCCATGGCGGCAGCAAGGAAGTGCCTCCTTATTATGCGGAAAATAATGATTCCTCCACCTGACAGCAGGGACCCACCGGACGGGCCATCGTATTTTGCGAAAAAACGTTTCCCCCTTGACTGCTGGGACCCACCAGCTACACCTTCACACGCAAGGAAGTGCCTGACAGTCGGGACCCACCTGGTCGAAGCGTATGTAGCATTGTCATTCTGGTCGTGAACATGAATGTACATATATACTGGTGGATGTAGAGGCGCGCACGTGTCGTAGTAGAGGCGGCCAGGGTGCAAGAAAGAAAATACGGCCACGTATGTGTACATACGGGCGGGGTCTCGAACGCCTACTCGCGCATACATACGGCGAGGGCTCATGTACATGGCTGGGTCGGAACGGAGAAACAACGTCGTCGTCGTGTTTATGGGTCGGAACGGAATGCGTGGTCGTGTTCATCGGGAGGGCTTGGACGGAAGAGGCGATGGAAACGAGGCCTGGCGTACCGCACAACGGAGGAAACGGACCTCCTACGGTCGAAACGGGGGTCCTGTTGATCGGGAGTGGTGTGGCGTACCGCAAAACGGAGGAAACGGACCTCCTACGGTCGAAACGGGGGTCCTGTTGATAGGGAGGGGTGTGGCGTACCGCAAAACGGAGGAAATGGACCTCCTACGGTCGAAATGGGGGTCCTGTTGATCGGGAGGGATGTGGCGTACCGCAAAACAGAGGAAACGGACTTGTGTTGGAGCGCTACGGTCGAAACAGGGGTCCTGTTCATCGGGAGGGGTGTGGCGTACCGCAAAACGGGACTCCACGGGATACTGTTCGTCTCCACCGTCAACCCCCTCCAGCCTCCACGAGCTACTGTTCATCCACCATCGACCTCCTCCAGCCTCCACGGGCTCCTATTCATCCAGCCTCCACTGGGCGCTACTCCACCGTCTACTGTTCAACCACCCCTCCACACCACGGGGTCCTGTTCAACCACCCCTCCACGGGCACCCCTCCATCGTCTACTGTTCATCCTGCCCTCCACACCACGCGGTCCTGTTCATCCACCCCTCCACGGGCACCCCTCCACCGTCTACTGTTCATCCAGCCCTCCACCACACGATGGGGTCCTGTTCATCCAAAGGCAACGCCACCGCTCACTGTTCATCCAACACCCCCTCCCCCTCCCCCCCGCGCAACACTCACTGTTCATCCCATCGATCGGCTTCAGTTAGCAGCAGTAGGCAAGGAATCGCTCGATTGGGTTCAGTTAACAGCCATCGATCAATCGCTCGGGTTCAGTAATGCGTAGCCTGCAGCGCAATCGCTCGGGTTTAGTTAGAGCCAACGCCTTGCACACATGCGCGTACGTGTACGAGAGAAACGCGCATCGCTCGGCCCCCGACCTCCCACCGTAACCGGGAACTCCCCGAAATTTTCCTCCCCCTCACTTCTACCACGGTTTTTTCCGTCATGGACAGCCCAAAGAATGTCATGCAGCTGCGTCTCCGTCTCGCCCAGGACGAAAAGCCCATTTTCTGTCATGATTTCTTTGTCATAGAAGTAGGAGCCCACCACATCTATGATGATACCATATTTTGTCACAATTATCATCATAGAAGTGTCATATGTATGACAGAAAAAAATTTCATTCGGCCCAAAATGTCACGGATGTGTCTTTTTTTGTAGTGTGTGCT'
    # # align_file = '/home/hukang/test/HiTE/demo/out29/NeuralTE_all/tir_work/judge_TIR_transposons_aed87012-16ca-40bc-a65a-090961d8724f/tir_copies_0_0/Chr1A-84162232-84190733-tir_20-tsd_GCT.blast.bed.fa.maf.fa.clean.fa'
    # # TE_type = 'tir'
    # # result_type = 'cons'
    # # is_TE, info, cons_seq, copy_num = judge_boundary_v5(cur_seq, align_file, debug, TE_type, plant, result_type)
    # # print(is_TE, info, cons_seq, copy_num)




    # output_dir = '/public/home/hpc194701009/ath_pan_genome/pan_genome/test'
    # merge_te_file = output_dir + '/pan_te.tmp.fa'
    #
    # unique_id = uuid.uuid4()
    # temp_dir = '/tmp/HiTE_module_test_' + str(unique_id)
    # create_or_clear_directory(temp_dir)
    # os.system('cp ' + merge_te_file + ' ' + temp_dir)
    # merge_te_file = temp_dir + '/pan_te.tmp.fa'
    # other_path, internal_path = split_internal_out(merge_te_file, temp_dir)
    #
    # terminal_coverage_threshold = 0.95
    # threads = 40
    # deredundant_for_LTR_v5(other_path, temp_dir, threads, 'terminal', terminal_coverage_threshold, debug=1)
    #
    # internal_coverage_threshold = 0.8
    # threads = 40
    # deredundant_for_LTR_v5(internal_path, temp_dir, threads, 'internal', internal_coverage_threshold, debug=1)
    #
    # # 计算完之后将结果拷贝回输出目录
    # copy_files(temp_dir, output_dir)
    #
    # # 删除临时目录
    # if os.path.exists(temp_dir):
    #     shutil.rmtree(temp_dir)


    # # panHiTE, panEDTA 和 Repbase 库三者的交集关系
    # # cd-hit-est不准确，很多序列其实和Repbase重叠，但是cd-hit-est无法聚类。比如：panHiTE_39-TIR_51#DNA/MULE能和MuDR-N3_AT比对上。
    # # 尝试使用BM_RM2的结果调整clusters
    # # cluster_file = '/home/hukang/test/demo/three_lib_compare/three_panTE.lib.cons.clstr'
    # cluster_file = '/home/hukang/test/demo/maize_lib_compare/maize_merge.fa.cons.clstr'
    # clusters = parse_clstr_file(cluster_file)
    # new_clusters = defaultdict()
    # # 遍历簇，记录sequence_name对应cluster_id，便于快速查询
    # seq_name2cluster_id = defaultdict()
    # for cluster, sequences in clusters.items():
    #     names = []
    #     for seq in sequences:
    #         seq = seq.split('#')[0]
    #         names.append(seq)
    #         seq_name2cluster_id[seq] = cluster
    #     new_clusters[cluster] = names
    # clusters = new_clusters
    #
    # # file_path = '/home/hukang/test/demo/three_lib_compare/repbase_panHiTE_rm2_test/file_final.0.1.txt'
    # # with open(file_path, 'r') as f:
    # #     for line in f:
    # #         parts = line.strip().split()
    # #         seq1 = parts[0].split(' ')[0]  # 第1列：序列1
    # #         seq2 = parts[4].split('#')[0]  # 第5列：序列2
    # #         cov1 = float(parts[3])  # 第4列：序列1对序列2的覆盖度
    # #         cov2 = float(parts[6])  # 第7列：序列2对序列1的覆盖度
    # #         similarity = (100 - float(parts[7])) / 100
    # #         if seq1 != seq2 and cov1 > 0.8 and cov2 > 0.8 and similarity > 0.8:
    # #             # 读取clusters，将panHiTE的cluster加入到Repbase cluster，并将原本的panHiTE的cluster删除掉
    # #             repbase_cluster_id = seq_name2cluster_id[seq1]
    # #             panHiTE_cluster_id = seq_name2cluster_id[seq2]
    # #             panHiTE_cluster = clusters.get(panHiTE_cluster_id, [])
    # #             repbase_cluster = clusters.get(repbase_cluster_id, [])
    # #             clusters[repbase_cluster_id] = repbase_cluster + panHiTE_cluster
    # #             if panHiTE_cluster_id in clusters:
    # #                 del clusters[panHiTE_cluster_id]
    # #
    # #
    # # file_path = '/home/hukang/test/demo/three_lib_compare/repbase_panEDTA_rm2_test/file_final.0.1.txt'
    # # with open(file_path, 'r') as f:
    # #     for line in f:
    # #         parts = line.strip().split()
    # #         seq1 = parts[0].split(' ')[0]  # 第1列：序列1
    # #         seq2 = parts[4].split('#')[0]  # 第5列：序列2
    # #         cov1 = float(parts[3])  # 第4列：序列1对序列2的覆盖度
    # #         cov2 = float(parts[6])  # 第7列：序列2对序列1的覆盖度
    # #         similarity = (100 - float(parts[7])) / 100
    # #         if seq1 != seq2 and cov1 > 0.8 and cov2 > 0.8 and similarity > 0.8:
    # #             # 读取clusters，将panHiTE的cluster加入到Repbase cluster，并将原本的panHiTE的cluster删除掉
    # #             repbase_cluster_id = seq_name2cluster_id[seq1]
    # #             panHiTE_cluster_id = seq_name2cluster_id[seq2]
    # #             panHiTE_cluster = clusters.get(panHiTE_cluster_id, [])
    # #             repbase_cluster = clusters.get(repbase_cluster_id, [])
    # #             clusters[repbase_cluster_id] = repbase_cluster + panHiTE_cluster
    # #             if panHiTE_cluster_id in clusters:
    # #                 del clusters[panHiTE_cluster_id]
    # #
    # # file_path = '/home/hukang/test/demo/three_lib_compare/panEDTA_panHiTE_rm2_test/file_final.0.1.txt'
    # # with open(file_path, 'r') as f:
    # #     for line in f:
    # #         parts = line.strip().split()
    # #         seq1 = parts[0].split(' ')[0]  # 第1列：序列1
    # #         seq2 = parts[4].split('#')[0]  # 第5列：序列2
    # #         cov1 = float(parts[3])  # 第4列：序列1对序列2的覆盖度
    # #         cov2 = float(parts[6])  # 第7列：序列2对序列1的覆盖度
    # #         similarity = (100 - float(parts[7])) / 100
    # #         if seq1 != seq2 and cov1 > 0.8 and cov2 > 0.8 and similarity > 0.8:
    # #             # 读取clusters，将panHiTE的cluster加入到Repbase cluster，并将原本的panHiTE的cluster删除掉
    # #             repbase_cluster_id = seq_name2cluster_id[seq1]
    # #             panHiTE_cluster_id = seq_name2cluster_id[seq2]
    # #             panHiTE_cluster = clusters.get(panHiTE_cluster_id, [])
    # #             repbase_cluster = clusters.get(repbase_cluster_id, [])
    # #             clusters[repbase_cluster_id] = repbase_cluster + panHiTE_cluster
    # #             if panHiTE_cluster_id in clusters:
    # #                 del clusters[panHiTE_cluster_id]
    #
    # total_panHiTE = []
    # total_panEDTA = []
    # total_Repbase = []
    #
    # only_panHiTE = []
    # only_panEDTA = []
    # only_Repbase = []
    #
    # panHiTE_vs_panEDTA = []
    # panHiTE_vs_Repbase = []
    # panEDTA_vs_Repbase = []
    # panHiTE_vs_panEDTA_vs_Repbase = []
    #
    # for cluster, sequences in clusters.items():
    #     total_sources = set()
    #     unique_sequences = []
    #     for seq in sequences:
    #         source = seq.split('_')[0]  # 获取序列的来源
    #         if source not in total_sources:
    #             unique_sequences.append(seq)
    #             if source == 'panHiTE':
    #                 total_panHiTE.append(seq)
    #             elif source == 'panEDTA':
    #                 total_panEDTA.append(seq)
    #             elif source == 'Repbase':
    #                 total_Repbase.append(seq)
    #         total_sources.add(source)
    #     if len(total_sources) == 1:
    #         seq = sequences[0]
    #         source = seq.split('_')[0]  # 获取序列的来源
    #         if source == 'panHiTE':
    #             only_panHiTE.append(seq)
    #         elif source == 'panEDTA':
    #             only_panEDTA.append(seq)
    #         elif source == 'Repbase':
    #             only_Repbase.append(seq)
    #     elif len(total_sources) == 2:
    #         if 'panHiTE' in total_sources and 'panEDTA' in total_sources:
    #             panHiTE_vs_panEDTA.append(unique_sequences)
    #         elif 'panHiTE' in total_sources and 'Repbase' in total_sources:
    #             panHiTE_vs_Repbase.append(unique_sequences)
    #         elif 'panEDTA' in total_sources and 'Repbase' in total_sources:
    #             panEDTA_vs_Repbase.append(unique_sequences)
    #     elif len(total_sources) == 3:
    #         panHiTE_vs_panEDTA_vs_Repbase.append(unique_sequences)
    #
    # # 绘制韦恩图
    # plt.figure(figsize=(6, 6))
    # venn = venn3(
    #     subsets=(
    #         len(only_panEDTA),
    #         len(only_panHiTE),
    #         len(panHiTE_vs_panEDTA),
    #         len(only_Repbase),
    #         len(panEDTA_vs_Repbase),
    #         len(panHiTE_vs_Repbase),
    #         len(panHiTE_vs_panEDTA_vs_Repbase)
    #     ),
    #     set_labels=('panEDTA', 'panHiTE', 'Repbase')
    # )
    #
    # # plt.title("Venn Diagram of Sequence Intersections")
    # # plt.savefig('/home/hukang/test/demo/three_lib_compare/output_plot.png')
    # plt.savefig('/home/hukang/test/demo/maize_lib_compare/output_plot.png')
    #
    # # 定义存储文件的路径
    # # output_dir = "/home/hukang/test/demo/three_lib_compare/"
    # output_dir = "/home/hukang/test/demo/maize_lib_compare"
    #
    # # 需要保存的集合及对应的文件名
    # output_files = {
    #     "total_panHiTE.txt": total_panHiTE,
    #     "total_panEDTA.txt": total_panEDTA,
    #     "total_Repbase.txt": total_Repbase,
    #     "only_panHiTE.txt": only_panHiTE,
    #     "only_panEDTA.txt": only_panEDTA,
    #     "only_Repbase.txt": only_Repbase,
    #     "panHiTE_vs_panEDTA.txt": panHiTE_vs_panEDTA,
    #     "panHiTE_vs_Repbase.txt": panHiTE_vs_Repbase,
    #     "panEDTA_vs_Repbase.txt": panEDTA_vs_Repbase,
    #     "panHiTE_vs_panEDTA_vs_Repbase.txt": panHiTE_vs_panEDTA_vs_Repbase
    # }
    #
    # # 遍历字典，写入文件
    # for filename, data_set in output_files.items():
    #     file_path = output_dir + filename
    #     with open(file_path, "w") as f:
    #         for item in data_set:
    #             f.write(str(item) + "\n")  # 每个序列单独占一行


    # # 我们想知道 panHiTE unique 里除去-LTR之外，有多少包含完整domain
    # tmp_output_dir = '/home/hukang/test/HiTE/demo'
    # threads = 48
    # protein_path = '/home/hukang/test/HiTE/library/RepeatPeps.lib'
    # confident_TE_consensus = '/home/hukang/test/HiTE/demo/865.fasta'
    # output_table = confident_TE_consensus + '.domain'
    # temp_dir = tmp_output_dir + '/domain'
    # get_domain_info(confident_TE_consensus, protein_path, output_table, threads, temp_dir)


    # # 重新对32个拟南芥的插入时间进行计算
    # tool_dir = '/public/home/hpc194701009/test/test/HiTE/tools'
    # # miu = str(7e-9)
    # # work_dir = '/public/home/hpc194701009/ath_pan_genome/pan_genome/rice/16_rice/panHiTE1/LTR_insertion_times'
    # # miu = str(1.3e-8)
    # work_dir = '/public/home/hpc194701009/ath_pan_genome/pan_genome/maize/panHiTE1/LTR_insertion_times'
    # miu = str(3e-8)
    #
    # # 获取所有物种目录
    # species_list = os.listdir(work_dir)
    #
    # # 设置并行进程池的大小（可以根据你的机器配置调整）
    # max_workers = 40  # 默认使用所有可用的 CPU 核心
    # print(f"Using {max_workers} processes in parallel.")
    #
    # # 使用 ProcessPoolExecutor 并行执行
    # with ProcessPoolExecutor(max_workers=max_workers) as executor:
    #     executor.map(run_command, species_list)

    # import pysam
    # def get_copies_minimap2(query_fasta_path, sam_path, full_length_coverage_threshold=0.95,
    #                         full_length_identity_threshold=0.8):
    #     """
    #     根据覆盖度和一致性阈值筛选比对结果。
    #
    #     参数:
    #         query_fasta_path (str): 查询序列的 FASTA 文件路径。
    #         sam_path (str): SAM 文件路径。
    #         full_length_coverage_threshold (float): 覆盖度阈值（默认 0.95）。
    #         full_length_identity_threshold (float): 一致性阈值（默认 0.8）。
    #
    #     返回:
    #         dict: 包含所有符合条件比对的字典，格式为 {query_name: [(subject_name, subject_start, subject_end, aligned_length, direct), ...]}。
    #     """
    #     # 读取 FASTA 文件
    #     query_names, query_contigs = read_fasta(query_fasta_path)
    #
    #     # 打开 SAM 文件
    #     samfile = pysam.AlignmentFile(sam_path, "r")
    #
    #     # 存储所有符合条件的比对
    #     all_copies = {}
    #
    #     # 遍历 SAM 文件中的每条比对记录
    #     for read in samfile:
    #         if not read.is_unmapped:  # 只处理比对上的 reads
    #             query_name = read.query_name
    #             query_length = len(query_contigs[query_name])  # query 的总长度
    #             aligned_length = read.query_alignment_length  # 比对上的 query 长度
    #             nm = read.get_tag("NM")  # 编辑距离
    #             matches = aligned_length - nm  # 匹配的碱基数
    #             identity = matches / aligned_length  # 一致性
    #             coverage = aligned_length / query_length  # 覆盖度
    #
    #             # 提取 subject 信息
    #             subject_name = read.reference_name
    #             subject_start = read.reference_start + 1  # SAM 坐标从 0 开始，转换为 1-based
    #             subject_end = read.reference_end  # pysam 自动计算结束位置
    #             direct = "-" if read.is_reverse else "+"  # 比对方向
    #
    #             # 如果覆盖度和一致性满足阈值，则保存比对信息
    #             if coverage >= full_length_coverage_threshold and identity >= full_length_identity_threshold:
    #                 if query_name not in all_copies:
    #                     all_copies[query_name] = []
    #                 all_copies[query_name].append((subject_name, subject_start, subject_end, aligned_length, direct))
    #
    #     # 关闭 SAM 文件
    #     samfile.close()
    #
    #     return all_copies
    #
    #
    # candidate_sequence_path = '/public/home/hpc194701009/ath_pan_genome/pan_genome/wheat/intact_ltr.fa'
    # sam_path = '/public/home/hpc194701009/ath_pan_genome/pan_genome/wheat/intact_ltr.asm20.sam'
    # all_copies_path = '/public/home/hpc194701009/ath_pan_genome/pan_genome/wheat/intact_ltr.asm20.copies'
    # all_copies = get_copies_minimap2(candidate_sequence_path, sam_path)
    # # 将字典存储为 JSON 文件
    # with open(all_copies_path, "w", encoding="utf-8") as file:
    #     json.dump(all_copies, file, indent=4, ensure_ascii=False)

    # import pysam
    #
    #
    # def get_copies_minimap2(query_fasta_path, sam_path, full_length_coverage_threshold=0.95):
    #     """
    #     根据覆盖度和一致性阈值筛选比对结果。
    #
    #     参数:
    #         query_fasta_path (str): 查询序列的 FASTA 文件路径。
    #         sam_path (str): SAM 文件路径。
    #         full_length_coverage_threshold (float): 覆盖度阈值（默认 0.95）。
    #
    #     返回:
    #         dict: 包含所有符合条件比对的字典，格式为 {query_name: [(subject_name, subject_start, subject_end, aligned_length, direct), ...]}。
    #     """
    #     # 读取 FASTA 文件
    #     query_names, query_contigs = read_fasta(query_fasta_path)
    #
    #     # 打开 SAM 文件
    #     samfile = pysam.AlignmentFile(sam_path, "r")
    #
    #     # 存储所有符合条件的比对
    #     all_copies = {}
    #
    #     # 遍历 SAM 文件中的每条比对记录
    #     for read in samfile:
    #         if not read.is_unmapped:  # 只处理比对上的 reads
    #             query_name = read.query_name
    #             query_length = len(query_contigs[query_name])  # query 的总长度
    #             aligned_length = read.query_alignment_length  # 比对上的 query 长度
    #             # nm = read.get_tag("NM")  # 编辑距离
    #             # matches = aligned_length - nm  # 匹配的碱基数
    #             query_coverage = aligned_length / query_length  # 覆盖度
    #             total_target_length = 0
    #             cigar = read.cigartuples
    #             for op, length in cigar:
    #                 if op == 0 or op == 2:  # M (match) 或 D (deletion)
    #                     total_target_length += length
    #             target_coverage = aligned_length / total_target_length  # 覆盖度
    #
    #             # 提取 subject 信息
    #             subject_name = read.reference_name
    #             subject_start = read.reference_start + 1  # SAM 坐标从 0 开始，转换为 1-based
    #             subject_end = read.reference_end  # pysam 自动计算结束位置
    #             direct = "-" if read.is_reverse else "+"  # 比对方向
    #
    #             # 如果覆盖度和一致性满足阈值，则保存比对信息
    #             if query_coverage >= full_length_coverage_threshold and target_coverage >= full_length_coverage_threshold:
    #                 if query_name not in all_copies:
    #                     all_copies[query_name] = []
    #                 all_copies[query_name].append((subject_name, subject_start, subject_end, aligned_length, direct))
    #
    #     # 关闭 SAM 文件
    #     samfile.close()
    #
    #     return all_copies
    #
    #
    # def filter_ltr_by_copy_num_minimap(candidate_sequence_path, threads, temp_dir, reference, multi_records=100,
    #                                    full_length_coverage=0.95):
    #     if os.path.exists(temp_dir):
    #         os.system('rm -rf ' + temp_dir)
    #     os.makedirs(temp_dir, exist_ok=True)
    #
    #     sam_path = os.path.join(temp_dir, 'intact_LTR.sam')
    #     # 1. 使用minimap2进行比对
    #     build_index_command = ['minimap2', '-I', '20G', '-d', reference + '.mmi', reference]
    #     # 执行命令并捕获输出
    #     try:
    #         result = subprocess.run(build_index_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    #                                 check=True)
    #     except subprocess.CalledProcessError as e:
    #         print(f"minimap2 build index failed with error code {e.returncode}.")
    #
    #     align_command = [
    #         'minimap2',
    #         '-ax', 'asm20',
    #         '-N', str(multi_records),
    #         '-p', '0.2',
    #         "-t", str(threads),
    #         reference + '.mmi',
    #         candidate_sequence_path
    #     ]
    #
    #     # 执行命令并重定向输出到文件
    #     try:
    #         with open(sam_path, "w") as sam_file:
    #             result = subprocess.run(align_command, stdout=sam_file, stderr=subprocess.DEVNULL, check=True)
    #     except subprocess.CalledProcessError as e:
    #         print(f"minimap2 alignment failed with error code {e.returncode}.")
    #     except Exception as e:
    #         print(f"An error occurred: {e}")
    #
    #     minimap_single_copy = os.path.join(temp_dir, 'minimap_single_copy.fa')
    #     all_ltr_names, all_ltr_contigs = read_fasta(candidate_sequence_path)
    #
    #     # 解析 sam 文件，获取拷贝
    #     all_copies = get_copies_minimap2(candidate_sequence_path, sam_path,
    #                                      full_length_coverage_threshold=full_length_coverage)
    #
    #     minimap_single_copy_contigs = {}
    #     for ltr_name in all_copies:
    #         if len(all_copies[ltr_name]) <= 1:
    #             minimap_single_copy_contigs[ltr_name] = all_ltr_contigs[ltr_name]
    #     store_fasta(minimap_single_copy_contigs, minimap_single_copy)
    #
    #     return all_copies, minimap_single_copy
    #
    # import pickle
    #
    # def filter_ltr_by_copy_num_blastn(candidate_sequence_path, threads, temp_dir, split_ref_dir, full_length_threshold,
    #                                max_copy_num=10):
    #     debug = 0
    #     if os.path.exists(temp_dir):
    #         os.system('rm -rf ' + temp_dir)
    #     os.makedirs(temp_dir, exist_ok=True)
    #
    #     batch_size = 1
    #     batch_id = 0
    #     names, contigs = read_fasta(candidate_sequence_path)
    #     split_files = []
    #     cur_contigs = {}
    #
    #     for i, name in enumerate(names):
    #         cur_file = os.path.join(temp_dir, f"{batch_id}.fa")
    #         cur_contigs[name] = contigs[name]
    #         if len(cur_contigs) == batch_size:
    #             store_fasta(cur_contigs, cur_file)
    #             split_files.append(cur_file)
    #             cur_contigs = {}
    #             batch_id += 1
    #
    #     if cur_contigs:
    #         cur_file = os.path.join(temp_dir, f"{batch_id}.fa")
    #         store_fasta(cur_contigs, cur_file)
    #         split_files.append(cur_file)
    #         batch_id += 1
    #
    #     ex = ProcessPoolExecutor(threads)
    #     jobs = []
    #     result_files = []
    #
    #     for idx, cur_split_file in enumerate(split_files):
    #         result_file = os.path.join(temp_dir, f"result_{idx}.pkl")
    #         result_files.append(result_file)
    #         job = ex.submit(run_and_save_result, cur_split_file, split_ref_dir, max_copy_num, full_length_threshold,
    #                         debug, result_file)
    #         jobs.append(job)
    #
    #     ex.shutdown(wait=True)
    #
    #     # 合并所有结果
    #     all_copies = {}
    #     for result_file in result_files:
    #         if os.path.exists(result_file):
    #             with open(result_file, "rb") as f:
    #                 cur_all_copies = pickle.load(f)
    #                 all_copies.update(cur_all_copies)
    #
    #     return all_copies
    #
    #
    # def run_and_save_result(split_file, split_ref_dir, max_copy_num, full_length_threshold, debug, result_file):
    #     """运行 get_full_length_copies_v1，并将结果存入文件"""
    #     cur_all_copies = get_full_length_copies_v1(split_file, split_ref_dir, max_copy_num, full_length_threshold,
    #                                                debug)
    #     with open(result_file, "wb") as f:
    #         pickle.dump(cur_all_copies, f)
    #
    #
    # def get_full_length_copies_v1(query_path, split_ref_dir, max_copy_num, full_length_threshold, debug):
    #     blastn2Results_path = query_path + '.blast.out'
    #     repeats_path = (query_path, split_ref_dir, blastn2Results_path)
    #     all_copies = multiple_alignment_blast_and_get_copies_v2(repeats_path, max_copy_num, full_length_threshold)
    #     if debug != 1:
    #         if os.path.exists(blastn2Results_path):
    #             os.remove(blastn2Results_path)
    #     return all_copies
    #
    #
    # def multiple_alignment_blast_and_get_copies_v2(repeats_path, max_copy_num, full_length_threshold):
    #     split_repeats_path = repeats_path[0]
    #     split_ref_dir = repeats_path[1]
    #     # raw_blastn2Results_path = repeats_path[2]
    #     # os.system('rm -f ' + raw_blastn2Results_path)
    #     blastn2Results_path = repeats_path[2]
    #     os.system('rm -f ' + blastn2Results_path)
    #     all_copies = {}
    #     repeat_names, repeat_contigs = read_fasta(split_repeats_path)
    #     remain_contigs = repeat_contigs
    #     for chr_name in os.listdir(split_ref_dir):
    #         # blastn2Results_path = raw_blastn2Results_path + '_' + str(chr_name)
    #         if len(remain_contigs) > 0:
    #             if not str(chr_name).endswith('.fa'):
    #                 continue
    #             chr_path = split_ref_dir + '/' + chr_name
    #             align_command = 'blastn -db ' + chr_path + ' -num_threads ' \
    #                             + str(
    #                 1) + ' -query ' + split_repeats_path + ' -evalue 1e-20 -outfmt 6 > ' + blastn2Results_path
    #             os.system(align_command)
    #             # 由于我们只需要100个拷贝，因此如果有序列已经满足了，就不需要进行后续的比对了，这样在mouse这样的高拷贝大基因组上减少运行时间
    #             cur_all_copies = get_copies_v2(blastn2Results_path, split_repeats_path, max_copy_num,
    #                                            full_length_threshold)
    #             for query_name in cur_all_copies.keys():
    #                 copy_list = cur_all_copies[query_name]
    #                 if query_name in all_copies:
    #                     prev_copy_list = all_copies[query_name]
    #                 else:
    #                     prev_copy_list = []
    #                 update_copy_list = prev_copy_list + copy_list
    #                 all_copies[query_name] = update_copy_list
    #                 if len(update_copy_list) >= max_copy_num:
    #                     del repeat_contigs[query_name]
    #             remain_contigs = repeat_contigs
    #             store_fasta(remain_contigs, split_repeats_path)
    #
    #     # all_copies = get_copies_v1(blastn2Results_path, split_repeats_path, '')
    #     return all_copies
    #
    #
    # def get_copies_v2(blastnResults_path, query_path, max_copy_num, full_length_threshold=0.95):
    #     query_names, query_contigs = read_fasta(query_path)
    #     longest_repeats = FMEA_new(query_path, blastnResults_path, full_length_threshold)
    #     # (query_name, old_query_start_pos, old_query_end_pos, subject_name, old_subject_start_pos, old_subject_end_pos)
    #     all_copies = {}
    #     for query_name in longest_repeats.keys():
    #         longest_queries = longest_repeats[query_name]
    #         longest_queries.sort(key=lambda x: -x[2])
    #         query_len = len(query_contigs[query_name])
    #         copies = []
    #         keeped_copies = set()
    #         for query in longest_queries:
    #             if len(copies) > max_copy_num:
    #                 break
    #             subject_name = query[3]
    #             subject_start = query[4]
    #             subject_end = query[5]
    #             cur_query_len = abs(query[2] - query[1])
    #             direct = '+'
    #             if subject_start > subject_end:
    #                 tmp = subject_start
    #                 subject_start = subject_end
    #                 subject_end = tmp
    #                 direct = '-'
    #             item = (subject_name, subject_start, subject_end)
    #
    #             if float(cur_query_len) / query_len >= full_length_threshold and item not in keeped_copies:
    #                 copies.append((subject_name, subject_start, subject_end, cur_query_len, direct))
    #                 keeped_copies.add(item)
    #         # copies.sort(key=lambda x: abs(x[3]-(x[2]-x[1]+1)))
    #         all_copies[query_name] = copies
    #
    #     return all_copies
    #
    #
    # def FMEA_new(query_path, blastn2Results_path, full_length_threshold):
    #     longest_repeats = {}
    #     if not os.path.exists(blastn2Results_path):
    #         return longest_repeats
    #     # parse blastn output, determine the repeat boundary
    #     # query_records = {query_name: {subject_name: [(q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end)] }}
    #     query_records = {}
    #     with open(blastn2Results_path, 'r') as f_r:
    #         for idx, line in enumerate(f_r):
    #             # print('current line idx: %d' % (idx))
    #             parts = line.split('\t')
    #             query_name = parts[0]
    #             subject_name = parts[1]
    #             identity = float(parts[2])
    #             alignment_len = int(parts[3])
    #             q_start = int(parts[6])
    #             q_end = int(parts[7])
    #             s_start = int(parts[8])
    #             s_end = int(parts[9])
    #             if query_name == subject_name and q_start == s_start and q_end == s_end:
    #                 continue
    #             if not query_records.__contains__(query_name):
    #                 query_records[query_name] = {}
    #             subject_dict = query_records[query_name]
    #
    #             if not subject_dict.__contains__(subject_name):
    #                 subject_dict[subject_name] = []
    #             subject_pos = subject_dict[subject_name]
    #             subject_pos.append((q_start, q_end, s_start, s_end))
    #     f_r.close()
    #
    #     query_names, query_contigs = read_fasta(query_path)
    #
    #     # 我们现在尝试新的策略，直接在生成簇的时候进行扩展，同时新的比对片段和所有的扩展片段进行比较，判断是否可以扩展
    #     for idx, query_name in enumerate(query_records.keys()):
    #         subject_dict = query_records[query_name]
    #
    #         if query_name not in query_contigs:
    #             continue
    #         query_len = len(query_contigs[query_name])
    #         skip_gap = query_len * (1 - full_length_threshold)
    #
    #         longest_queries = []
    #         for subject_name in subject_dict.keys():
    #             subject_pos = subject_dict[subject_name]
    #
    #             # cluster all closed fragments, split forward and reverse records
    #             forward_pos = []
    #             reverse_pos = []
    #             for pos_item in subject_pos:
    #                 if pos_item[2] > pos_item[3]:
    #                     reverse_pos.append(pos_item)
    #                 else:
    #                     forward_pos.append(pos_item)
    #             forward_pos.sort(key=lambda x: (x[2], x[3]))
    #             reverse_pos.sort(key=lambda x: (-x[2], -x[3]))
    #
    #             forward_long_frags = {}
    #             frag_index_array = []
    #             frag_index = 0
    #             for k, frag in enumerate(forward_pos):
    #                 is_update = False
    #                 cur_subject_start = frag[2]
    #                 cur_subject_end = frag[3]
    #                 cur_query_start = frag[0]
    #                 cur_query_end = frag[1]
    #                 for cur_frag_index in reversed(frag_index_array):
    #                     cur_frag = forward_long_frags[cur_frag_index]
    #                     prev_subject_start = cur_frag[2]
    #                     prev_subject_end = cur_frag[3]
    #                     prev_query_start = cur_frag[0]
    #                     prev_query_end = cur_frag[1]
    #
    #                     if cur_subject_start - prev_subject_end >= skip_gap:
    #                         break
    #
    #                     if cur_subject_end > prev_subject_end:
    #                         # forward extend
    #                         if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
    #                                 and cur_subject_start - prev_subject_end < skip_gap:  # \
    #                             # extend frag
    #                             prev_query_start = prev_query_start if prev_query_start < cur_query_start else cur_query_start
    #                             prev_query_end = cur_query_end
    #                             prev_subject_start = prev_subject_start if prev_subject_start < cur_subject_start else cur_subject_start
    #                             prev_subject_end = cur_subject_end
    #                             extend_frag = (
    #                             prev_query_start, prev_query_end, prev_subject_start, prev_subject_end, subject_name)
    #                             forward_long_frags[cur_frag_index] = extend_frag
    #                             is_update = True
    #                 if not is_update:
    #                     frag_index_array.append(frag_index)
    #                     forward_long_frags[frag_index] = (
    #                     cur_query_start, cur_query_end, cur_subject_start, cur_subject_end, subject_name)
    #                     frag_index += 1
    #             longest_queries += list(forward_long_frags.values())
    #
    #             reverse_long_frags = {}
    #             frag_index_array = []
    #             frag_index = 0
    #             for k, frag in enumerate(reverse_pos):
    #                 is_update = False
    #                 cur_subject_start = frag[2]
    #                 cur_subject_end = frag[3]
    #                 cur_query_start = frag[0]
    #                 cur_query_end = frag[1]
    #                 for cur_frag_index in reversed(frag_index_array):
    #                     cur_frag = reverse_long_frags[cur_frag_index]
    #                     prev_subject_start = cur_frag[2]
    #                     prev_subject_end = cur_frag[3]
    #                     prev_query_start = cur_frag[0]
    #                     prev_query_end = cur_frag[1]
    #
    #                     if prev_subject_end - cur_subject_start >= skip_gap:
    #                         break
    #
    #                     # reverse
    #                     if cur_subject_end < prev_subject_end:
    #                         # reverse extend
    #                         if cur_query_start - prev_query_end < skip_gap and cur_query_end > prev_query_end \
    #                                 and prev_subject_end - cur_subject_start < skip_gap:  # \
    #                             # extend frag
    #                             prev_query_start = prev_query_start
    #                             prev_query_end = cur_query_end
    #                             prev_subject_start = prev_subject_start if prev_subject_start > cur_subject_start else cur_subject_start
    #                             prev_subject_end = cur_subject_end
    #                             extend_frag = (
    #                             prev_query_start, prev_query_end, prev_subject_start, prev_subject_end, subject_name)
    #                             reverse_long_frags[cur_frag_index] = extend_frag
    #                             is_update = True
    #                 if not is_update:
    #                     frag_index_array.append(frag_index)
    #                     reverse_long_frags[frag_index] = (
    #                     cur_query_start, cur_query_end, cur_subject_start, cur_subject_end, subject_name)
    #                     frag_index += 1
    #             longest_queries += list(reverse_long_frags.values())
    #
    #         if not longest_repeats.__contains__(query_name):
    #             longest_repeats[query_name] = []
    #         cur_longest_repeats = longest_repeats[query_name]
    #         for repeat in longest_queries:
    #             # Subject序列处理流程
    #             subject_name = repeat[4]
    #             old_subject_start_pos = repeat[2] - 1
    #             old_subject_end_pos = repeat[3]
    #             # Query序列处理流程
    #             old_query_start_pos = repeat[0] - 1
    #             old_query_end_pos = repeat[1]
    #             cur_query_seq_len = abs(old_query_end_pos - old_query_start_pos)
    #             cur_longest_repeats.append((query_name, old_query_start_pos, old_query_end_pos, subject_name,
    #                                         old_subject_start_pos, old_subject_end_pos))
    #
    #     return longest_repeats


    # work_dir = '/public/home/hpc194701009/ath_pan_genome/pan_genome/ath/32_ath/panHiTE2_32_ath/pan_run_hite_single/01.col.fa'
    # candidate_sequence_path = work_dir + '/intact_LTR.fa'
    # reference = '/public/home/hpc194701009/ath_pan_genome/pan_genome/ath/32_ath/genomes/01.col.fa'
    # threads = 40
    # temp_dir = work_dir + '/copies_temp_bak'
    # full_length_threshold = 0.95
    #
    # # split_ref_dir = work_dir + '/ref_chr'
    # # test_home = '/public/home/hpc194701009/test/test/HiTE/bin/FiLTR-main/src'
    # # split_genome_command = [
    # #     "cd", test_home, "&&", "python3", f"{test_home}/split_genome_chunks.py",
    # #     "-g", reference, "--tmp_output_dir", work_dir
    # # ]
    # # split_genome_command = " ".join(split_genome_command)
    # #
    # # try:
    # #     result = subprocess.run(split_genome_command, shell=True, check=True)
    # # except subprocess.CalledProcessError as e:
    # #     print(f"BLAT alignment failed with error code {e.returncode}.")
    #
    # all_copies, _ = filter_ltr_by_copy_num_minimap(candidate_sequence_path, threads, temp_dir, reference, multi_records=20)
    # all_copies_file = temp_dir + '/all_copies.minimap2.asm20.json'
    # # all_copies = filter_ltr_by_copy_num_sub(candidate_sequence_path, threads, temp_dir, split_ref_dir, full_length_threshold, max_copy_num=20)
    # # all_copies_file = temp_dir + '/all_copies.blastn.json'
    # with open(all_copies_file, "w", encoding="utf-8") as file:
    #     json.dump(all_copies, file, indent=4, ensure_ascii=False)











    # def filter_intact_ltr(intact_ltr_path, tmp_output_dir, threads):
    #     """
    #     使用minimap2对LTR序列进行自比对，去除冗余序列，减少后续计算量。
    #
    #     参数:
    #         intact_ltr_path (str): 完整LTR序列的输入文件路径。
    #         tmp_output_dir (str): 临时输出目录路径。
    #         threads (int): 运行minimap2时使用的线程数。
    #         intact_ltrs (dict): 完整LTR序列的字典，键为LTR名称，值为序列。
    #
    #     返回:
    #         str: 去除冗余后的完整LTR序列文件路径。
    #     """
    #     # 创建minimap2比对的临时目录
    #     self_minimap_temp_dir = os.path.join(tmp_output_dir, 'intact_ltr_self_minimap_filter')
    #     os.makedirs(self_minimap_temp_dir, exist_ok=True)
    #
    #     # 使用minimap2比对去除冗余
    #     all_copies, _ = filter_ltr_by_copy_num_minimap(intact_ltr_path, threads, self_minimap_temp_dir, intact_ltr_path,
    #                                                    multi_records=100)
    #
    #     # 初始化冗余和非冗余LTR名称集合
    #     redundant_ltr_names = set()
    #     unique_ltr_names = set()
    #
    #     # 遍历所有比对结果，筛选出非冗余的LTR名称
    #     for ltr_name in all_copies:
    #         if ltr_name not in redundant_ltr_names:
    #             unique_ltr_names.add(ltr_name)
    #             for copy in all_copies[ltr_name]:
    #                 copy_name = copy[0]
    #                 redundant_ltr_names.add(copy_name)
    #
    #     intact_ltr_names, intact_ltrs = read_fasta(intact_ltr_path)
    #     # 构建去除冗余后的完整LTR序列字典
    #     clean_intact_ltr_path = intact_ltr_path + '.clean'
    #     intact_ltrs_clean = {ltr_name: intact_ltrs[ltr_name] for ltr_name in unique_ltr_names}
    #
    #     # 将去除冗余后的完整LTR序列存储为FASTA文件
    #     store_fasta(intact_ltrs_clean, clean_intact_ltr_path)
    #
    #     shutil.rmtree(self_minimap_temp_dir)
    #
    #     return clean_intact_ltr_path
    #
    #
    # def run_test(temp_dir, intact_ltr_path):
    #     for _ in range(3):
    #         intact_ltr_path = filter_intact_ltr(intact_ltr_path, temp_dir, threads)
    #
    #     # 首先用minimap2 asm5 快速过滤出单拷贝的LTR (可能会误将一些低拷贝的认为是单拷贝)；然后用 blastn 再次对单拷贝的LTR进行比对，准确地获取单拷贝
    #     minimap_temp_dir = os.path.join(temp_dir, 'intact_ltr_minimap_filter')
    #     ltr_copies, minimap_single_copy = filter_ltr_by_copy_num_minimap(intact_ltr_path, threads, minimap_temp_dir,
    #                                                                      reference, multi_records=50)
    #
    #     # split_ref_dir = temp_dir + '/ref_chr'
    #     # test_home = '/public/home/hpc194701009/test/test/HiTE/bin/FiLTR-main/src'
    #     # split_genome_command = [
    #     #     "cd", test_home, "&&", "python3", f"{test_home}/split_genome_chunks.py",
    #     #     "-g", reference, "--tmp_output_dir", temp_dir
    #     # ]
    #     # split_genome_command = " ".join(split_genome_command)
    #     # try:
    #     #     result = subprocess.run(split_genome_command, shell=True, check=True)
    #     # except subprocess.CalledProcessError as e:
    #     #     print(f"BLAT alignment failed with error code {e.returncode}.")
    #     #
    #     # blastn_temp_dir = os.path.join(temp_dir, 'intact_ltr_blastn_filter')
    #     # minimap_single_copy_copies = filter_ltr_by_copy_num_blastn(minimap_single_copy, threads, blastn_temp_dir,
    #     #                                                            split_ref_dir, full_length_threshold,
    #     #                                                            max_copy_num=10)
    #     # # 对 ltr_copies 进行更新
    #     # for ltr_name in minimap_single_copy_copies:
    #     #     ltr_copies[ltr_name] = minimap_single_copy_copies[ltr_name]
    #
    #     all_intact_ltr_copies = os.path.join(temp_dir, 'all_intact_ltr.copies')
    #     save_to_file(ltr_copies, all_intact_ltr_copies)
    #     print('all intact ltr num: ' + str(len(ltr_copies)))
    #
    #     # 去掉单拷贝
    #     multi_copies_path = os.path.join(temp_dir, 'multi_intact_ltr.copies')
    #     multi_copies = {}
    #     for ltr_name in ltr_copies:
    #         if len(ltr_copies[ltr_name]) > 1:
    #             multi_copies[ltr_name] = ltr_copies[ltr_name]
    #     save_to_file(multi_copies, multi_copies_path)
    #     print('multiple intact ltr num: ' + str(len(multi_copies)))
    #
    #
    # tmp_output_dir = '/public/home/hpc194701009/ath_pan_genome/pan_genome/wheat/test'
    # intact_ltr_path = tmp_output_dir + '/intact_ltr.fa'
    # reference = '/public/home/hpc194701009/ath_pan_genome/pan_genome/wheat/genomes/CM42.fa'
    # threads = 48
    # full_length_threshold = 0.95
    #
    # # 创建本地临时目录，存储计算结果
    # unique_id = uuid.uuid4()
    # temp_dir = '/tmp/test' + str(unique_id)
    # try:
    #     create_or_clear_directory(temp_dir)
    #     # 运行计算任务
    #     run_test(temp_dir, intact_ltr_path)
    #
    #     # 计算完之后将结果拷贝回输出目录
    #     copy_files(temp_dir, tmp_output_dir)
    #
    # except Exception as e:
    #     # 如果出现异常，打印错误信息并删除临时目录
    #     print(f"An error occurred: {e}")
    #     if os.path.exists(temp_dir):
    #         shutil.rmtree(temp_dir)
    #     raise  # 重新抛出异常，以便上层代码可以处理
    #
    # else:
    #     # 如果没有异常，删除临时目录
    #     if os.path.exists(temp_dir):
    #         shutil.rmtree(temp_dir)


    # def genome_genome_align_minimap2(query_fasta, target_fasta, output_paf, threads, multi_records=100):
    #     """
    #     使用 Minimap2 对两个 FASTA 文件进行比对，并将结果保存为 PAF 格式。
    #     :param query_fasta: 查询序列文件路径
    #     :param target_fasta: 目标序列文件路径
    #     :param output_paf: 输出的 PAF 文件路径
    #     """
    #     command = [
    #         "minimap2",
    #         "-x", "asm20",  # 使用适合组装的参数
    #         # "-N", str(multi_records),
    #         # "-p", "0.2",
    #         "-P",
    #         "-g 500",
    #         "--max-chain-skip 5",
    #         "-f 0",
    #         "--no-long-join",
    #         "-c",
    #         "-t", str(threads),
    #         target_fasta,
    #         query_fasta,
    #         "-o", output_paf
    #     ]
    #     # with open(output_sam, "w") as sam_file:
    #     #     result = subprocess.run(command, stdout=sam_file, stderr=subprocess.DEVNULL, check=True)
    #     subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    # # # 测试一下粗比对是否可以用minimap2替代
    # work_dir = '/home/hukang/test/HiTE/demo/out4'
    # query_path = os.path.join(work_dir, '9_target.fa')
    # target_path = os.path.join(work_dir, '0_target.fa')
    # thread = 1
    # lib_out = work_dir + '/test.out'
    # lib_paf = work_dir + '/test.paf'
    # if os.path.exists(lib_paf):
    #     os.remove(lib_paf)
    # # target_dir = os.path.join(work_dir, 'longest_repeats_blast_25')
    # # for target_name in os.listdir(target_dir):
    # #     if target_name.endswith('.fa'):
    # #         target_path = os.path.join(target_dir, target_name)
    # #         cur_lib_paf = target_name + '.paf'
    # #         te_genome_align_minimap2(query_path, target_path, cur_lib_paf, thread, multi_records=100)
    # #         os.system('cat ' + cur_lib_paf + ' >> ' + lib_paf)
    #
    # genome_genome_align_minimap2(query_path, target_path, lib_paf, thread, multi_records=100)
    # # 将 PAF 文件转换为类似 BLAST outfmt 6 的格式
    # convert_paf_to_blast6_format(lib_paf, lib_out)

    # lib_out = '/tmp/coarse_boundary_083091c4-a07e-4b36-8776-2298bab374f7/prev_TE_0.out'
    # TE_lib = '/home/hukang/test/HiTE/demo/out14/prev_TE.fa'
    # genome_path = '/tmp/coarse_boundary_083091c4-a07e-4b36-8776-2298bab374f7/filter_tandem_0.fa'
    # tmp_blast_dir = '/tmp/coarse_boundary_083091c4-a07e-4b36-8776-2298bab374f7/'
    # coverage_threshold = 0.8
    # category = 'Total'
    # debug = 1
    # output_files = generate_full_length_out_v1(lib_out, TE_lib, genome_path, tmp_blast_dir, '',
    #                                            coverage_threshold, category, debug=debug)
    #
    # masked_genome_path = genome_path + '.masked'
    # ref_names, ref_contigs = read_fasta(genome_path)
    # for output_file in output_files:
    #     sorted_lines = load_from_file(output_file)
    #     # mask genome
    #     for query_name, chr_name, chr_start, chr_end in sorted_lines:
    #         print(query_name, chr_name, chr_start, chr_end)
    #         start = chr_start - 1
    #         end = chr_end
    #         ref_seq = ref_contigs[chr_name]
    #         mask_seq = ref_seq[:start] + 'N' * (end - start) + ref_seq[end:]
    #         ref_contigs[chr_name] = mask_seq
    #     del sorted_lines
    # store_fasta(ref_contigs, masked_genome_path)
    # del ref_contigs

    # candidate_sequence_path = '/home/hukang/test/HiTE/demo/out4/test.fa'
    # full_length_coverage = 0.98
    # sam_path = '/home/hukang/test/HiTE/demo/out4/intact_LTR.sam'
    # all_copies = get_copies_minimap2(candidate_sequence_path, sam_path,
    #                                  full_length_coverage_threshold=full_length_coverage)
    # print(all_copies)


    # # 我们现在测试一下minimap2和blastn分别在判断一条假阳性的内部序列时，会得到几条拷贝
    # candidate_sequence_path = '/public/home/hpc194701009/LTR_libraries/Ours/zebrafish/HiTE_LTR_latest1/test.fa'
    # threads = 40
    # full_length_threshold = 0.95
    # temp_dir = '/public/home/hpc194701009/LTR_libraries/Ours/zebrafish/HiTE_LTR_latest1/temp_dir'
    # reference = '/public/home/hpc194701009/LTR_libraries/Ours/zebrafish/HiTE_LTR_latest1/genome.rename.fa'
    # all_copies, minimap_single_copy = filter_ltr_by_copy_num_minimap(candidate_sequence_path, threads, temp_dir, reference, multi_records=100,
    #                                full_length_coverage=0.95)
    # print('minimap2:')
    # print(all_copies)
    #
    # tmp_output_dir = '/public/home/hpc194701009/LTR_libraries/Ours/zebrafish/HiTE_LTR_latest1'
    # blastn_temp_dir = os.path.join(tmp_output_dir, 'intact_ltr_blastn_filter')
    # split_ref_dir = tmp_output_dir + '/ref_chr'
    # test_home = '/public/home/hpc194701009/test/test/HiTE/bin/FiLTR-main/src'
    # split_genome_command = [
    #     "cd", test_home, "&&", "python3", f"{test_home}/split_genome_chunks.py",
    #     "-g", reference, "--tmp_output_dir", tmp_output_dir
    # ]
    # split_genome_command = " ".join(split_genome_command)
    # try:
    #     result = subprocess.run(split_genome_command, shell=True, check=True)
    # except subprocess.CalledProcessError as e:
    #     print(f"BLAT alignment failed with error code {e.returncode}.")
    # ltr_copies = filter_ltr_by_copy_num_blastn(candidate_sequence_path, threads, blastn_temp_dir, split_ref_dir, full_length_threshold, max_copy_num=10)
    # print('blastn:')
    # print(ltr_copies)




    # temp_dir = '/public/home/hpc194701009/LTR_libraries/Ours/zebrafish/HiTE_LTR_latest3'
    # reference = temp_dir + '/genome.rename.fa'
    # curated_lib = '/public/home/hpc194701009/test/test/HiTE/library/zebrafish.ltr.ref'
    # TE_lib = temp_dir + '/confident_TE.cons.fa'
    # log = Logger(temp_dir + '/benchmarking.log', level='debug')
    # BM_RM2 = 1
    # BM_EDTA = 1
    # BM_HiTE = 1
    # threads = 48
    # coverage_threshold = 0.95
    # recover = 0
    # EDTA_home = '/public/home/hpc194701009/repeat_detect_tools/EDTA'
    # run_benchmarking(temp_dir, reference, curated_lib, TE_lib, BM_RM2, BM_EDTA, BM_HiTE, threads, coverage_threshold,
    #                  recover, EDTA_home, log)

    # # 测试test中是否包含蛋白质
    # tmp_output_dir = '/public/home/hpc194701009/LTR_libraries/Ours/zebrafish/HiTE_LTR_latest3'
    # single_copy_internals_file = tmp_output_dir + '/test.fa'
    # temp_dir = tmp_output_dir + '/ltr_domain'
    # output_table = single_copy_internals_file + '.ltr_domain'
    # ltr_protein_db = '/public/home/hpc194701009/test/test/HiTE/bin/FiLTR-main/databases/LTRPeps.lib'
    # threads = 48
    # get_domain_info(single_copy_internals_file, ltr_protein_db, output_table, threads, temp_dir)
    # is_single_ltr_has_intact_protein = {}
    # protein_names, protein_contigs = read_fasta(ltr_protein_db)
    # with open(output_table, 'r') as f_r:
    #     for i, line in enumerate(f_r):
    #         if i < 2:
    #             continue
    #         parts = line.split('\t')
    #         te_name = parts[0]
    #         protein_name = parts[1]
    #         protein_start = int(parts[4])
    #         protein_end = int(parts[5])
    #         intact_protein_len = len(protein_contigs[protein_name])
    #         if float(abs(protein_end - protein_start)) / intact_protein_len >= 0.95:
    #             is_single_ltr_has_intact_protein[te_name] = True
    # print(is_single_ltr_has_intact_protein)






    # BM_EDTA()
    # BM_HiTE()

    # # 对TIR protein 进行筛选，去掉非法的
    # tir_protein_db = cur_dir + '/library/TIRPeps.lib'
    # protein_names, protein_contigs = read_fasta(tir_protein_db)
    # store_fasta(protein_contigs, tir_protein_db)


    # # 将RepeatModeler2中的non-ltr抽取出来
    # work_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/Repbase_28.06/human'
    # te_path = work_dir + '/human.lib'
    # names, contigs = read_fasta(te_path)
    # label_set = set()
    # for name in names:
    #     label = name.split('#')[1]
    #     label_set.add(label)
    # print(label_set)
    # non_ltr_path = work_dir + '/non_ltr.fa'
    # non_ltr_contigs = {}
    # for name in names:
    #     if 'LINE' in name or 'SINE' in name:
    #         non_ltr_contigs[name] = contigs[name]
    # store_fasta(non_ltr_contigs, non_ltr_path)

    # align_file = '/home/hukang/test/HiTE/demo/test1/non_ltr_copies_0_0/chr_0:13880942-13881327.blast.bed.fa.maf.fa'
    # debug = 1
    # #is_TE, info, cons_seq = judge_boundary_v8(align_file, debug)
    # cur_seq = 'TCCGCATTTGTTTGTCCCAAATTTGATGGAACAGTGTTCTGAAGATTACTGTAGATTACTGTAGAAAATCTCACAGTTTGCAGTTTTCTTTTGTCAACGATTTTTATCCACATTTTGACGAAAAAACTACAAAACGAAAATGGCTTCTACAGTAATCTAAACACCATCCCATCAATTTTCAGACACAAAATTACCTGGAAAAATCGGAAAAGAATACTTCAAAAATAGCTGAAAACCCGAAAAATGTCTGAAAATTTTGCCGAAAAAATCGGATAAAAAA'
    # TE_type = 'non-ltr'
    # plant = 0
    # result_type = 'cons'
    # is_TE, info, cons_seq = judge_boundary_v9(cur_seq, align_file, debug, TE_type, plant, result_type)
    # print(is_TE, info, cons_seq)

    # # work_dir = '/tmp/run_hite_main'
    # # # candidate_tir_path = work_dir + '/more_tir.fa'
    # # # candidate_tir_path = work_dir + '/confident_helitron.fa'
    # # # candidate_tir_path = work_dir + '/confident_tir.fa'
    # # candidate_tir_path = work_dir + '/test.fa'
    # # confident_tir_path = work_dir + '/test.final.fa'
    # # flanking_len = 50
    # # # reference = '/homeb/hukang/KmerRepFinder_test/genome/ucsc_genome/rice.fa'
    # # # reference = '/home/hukang/test/test/HiTE/demo/32_ath/genomes/39.ah_7.fa'
    # # # reference = '/home/hukang/test/test/HiTE/demo/32_ath/genomes/03.yilong.fa'
    # # # reference = '/home/hukang/test/test/HiTE/demo/32_ath/genomes/01.col.fa'
    # # # reference = '/home/hukang/test/test/HiTE/demo/32_ath/genomes/23.sij_1.fa'
    # # # reference = '/home/hukang/test/test/HiTE/demo/32_ath/genomes/41.sorbo.fa'
    # # # reference = '/home/hukang/test/test/HiTE/demo/32_ath/genomes/04.bor_1.fa'
    # # # # reference = '/homeb/hukang/KmerRepFinder_test/genome/04.bor_1.fa'
    # # # # reference = '/homeb/hukang/KmerRepFinder_test/genome/02.tibet.fa'
    # reference = '/home/hukang/repeat_detect_tools/TEtrimmer-main/liuyuyang/athaliana.fna'
    # # test_home = '/home/hukang/test/HiTE/module'
    # tmp_output_dir = '/home/hukang/repeat_detect_tools/TEtrimmer-main/liuyuyang/tirs'
    # work_dir = tmp_output_dir
    # if not os.path.exists(tmp_output_dir):
    #     os.makedirs(tmp_output_dir)
    # chrom_seg_length = 1_000_000
    # chunk_size = 400
    # split_genome_command = 'split_genome_chunks.py -g ' \
    #                        + reference + ' --tmp_output_dir ' + tmp_output_dir \
    #                        + ' --chrom_seg_length ' + str(chrom_seg_length) + ' --chunk_size ' + str(chunk_size)
    # # os.system(split_genome_command)
    #

    # reference = '/home/hukang/test/test/HiTE/demo/32_ath/genomes/13.got_22.fa'
    # tmp_output_dir = '/home/hukang/test/HiTE/demo/ath/BM_EDTA'
    # work_dir = tmp_output_dir
    # TE_type = 'tir'
    # split_ref_dir = tmp_output_dir + '/ref_chr'
    # threads = 40
    # ref_index = 0
    # subset_script_path = '/home/hukang/test/HiTE/tools/ready_for_MSA.sh'
    # plant = 0
    # debug = 1
    # # candidate_tir_path = '/home/hukang/repeat_detect_tools/TEtrimmer-main/liuyuyang/result.fa'
    # candidate_tir_path = tmp_output_dir + '/test.fa'
    # contig_names, contigs = read_fasta(candidate_tir_path)
    # #candidate_tir_path_rename = '/home/hukang/repeat_detect_tools/TEtrimmer-main/liuyuyang/result.rename.fa'
    # candidate_tir_path_rename = tmp_output_dir + '/test.fa'
    # store_fasta(contigs, candidate_tir_path_rename)
    # confident_tir_path = tmp_output_dir + '/confident.fa'
    # all_low_copy = tmp_output_dir + '/all_low_copy.fa'
    # log = Logger(tmp_output_dir + '/test.log', level='debug')
    # flanking_len = 200
    # flank_region_align_v5(candidate_tir_path_rename, confident_tir_path, flanking_len, reference, split_ref_dir,
    #                       TE_type, work_dir, threads, ref_index, log, subset_script_path,
    #                       plant, debug, 1, all_low_copy, result_type='cons')











    # seq = 'TAAAACGACTTAAATTTGGCCCAAAATAGCCGAAAAGGACGCCGGTTAAAGTCGACGGTCCGTCTCTTTGCGGATGACACCATCATCTACCTGACTGTAACAAGTGAAGAAAGTTGCAAACTACTGCAAAAAGACCTTGATGTTCTTCAGGAATGGGGTAAAACCTGGAAAATGGAATTCAATACCACTAAATGCGAAGTGCTGCGAATTAAACCCAAACGGCAATTCCAACACACTTACTACCTTCATGGCGTTCCACTTAAGGAAGTCAAACAGAGCAAATACCTGGGCGTACAAATATCACACGACTTAAAATGGAATAAACATATAGACGATATAACAGCGAAATCAAGCAGAACTCTTGGCTTCCTTAAACGGAACCTGAAGATAAATTCCCCAACACTTAAAGCAAAAGCCTACCAAGGTCTGGTTCGGCCGAAATTAGAATACTGTGCTAGTGTGTGGGATCCCCGTTCACGTGTTGAAAATAATGGATCGTACAAGATAGAGATGATACAACGCCGTGCTGCGCGGTGGTGCCTACGGAGATACGGAAGAACAGATAGTGTTACTGCAATGCTTGATGCCCTAGGATGGCGGACTTTGGAGCAGCGAAGAGCAGATATCCGACTGACAACTCTGTTTAAAATCACCAGAGGACTACTCTCAGTCAATACTCAAGGTCTCCTACGTCCCGTAACCCGCAAAACCCGTCGTAGTCATCCAGAATCGTTTATCCCGCTAAGTACTTCTTCGGCATCTGAACACCTCTCTTTTTTCTCAAGAACAATAATCCAATGGAACAATCTGCCTGCCCTCTTGTTTAACGAACCGTGCTCTTTACCTATTTTTAAAGACAAAATAAGCCGTTTAACTCACCTTCCTGTAAATTAAGCTATAATCTGACCCCCCCAGCCAACAATCCTCAAATGTGAGGGGCAGGCTGTAATTGGTAAGTAAGTAAGTAAGAGGACGCCGGTTAAAAAACTCACCATCAACACCAAATTCGTCTCTGAGCAAATTGC'
    # # seq = 'TTAAAACGACTTAAATTTGGCCCAAAATAGCCGAAAAGGACGCCGGTTAAAGTCGACGGTCCGTCTCTTTGCGGATGACACCATCATCTACCTGACTGTAACAAGTGAAGAAAGTTGCAAACTACTGCAAAAAGACCTTGATGTTCTTCAGGAATGGGGTAAAACCTGGAAAATGGAATTCAATACCACTAAATGCGAAGTGCTGCGAATTAAACCCAAACGGCAATTCCAACACACTTACTACCTTCATGGCGTTCCACTTAAGGAAGTCAAACAGAGCAAATACCTGGGCGTACAAATATCACACGACTTAAAATGGAATAAACATATAGACGATATAACAGCGAAATCAAGCAGAACTCTTGGCTTCCTTAAACGGAACCTGAAGATAAATTCCCCAACACTTAAAGCAAAAGCCTACCAAGGTCTGGTTCGGCCGAAATTAGAATACTGTGCTAGTGTGTGGGATCCCCGTTCACGTGTTGAAAATAATGGATCGTACAAGATAGAGATGATACAACGCCGTGCTGCGCGGTGGTGCCTACGGAGATACGGAAGAACAGATAGTGTTACTGCAATGCTTGATGCCCTAGGATGGCGGACTTTGGAGCAGCGAAGAGCAGATATCCGACTGACAACTCTGTTTAAAATCACCAGAGGACTACTCTCAGTCAATACTCAAGGTCTCCTACGTCCCGTAACCCGCAAAACCCGTCGTAGTCATCCAGAATCGTTTATCCCGCTAAGTACTTCTTCGGCATCTGAACACCTCTCTTTTTTCTCAAGAACAATAATCCAATGGAACAATCTGCCTGCCCTCTTGTTTAACGAACCGTGCTCTTTACCTATTTTTAAAGACAAAATAAGCCGTTTAACTCACCTTCCTGTAAATTAAGCTATAATCTGACCCCCCCAGCCAACAATCCTCAAATGTGAGGGGCAGGCTGTAATTGGTAATGTGGACGCCGGTTAAAAAACTCACCATCAACACCAAATTCGTCTCTGAGCA'
    # flanking_len = 50
    # raw_end = len(seq) - flanking_len
    #
    # end_5_window_size = 25
    # TSD_sizes = list(range(8, 21))
    # found_TSD, TSD_seq, non_ltr_seq = search_polyA_TSD(seq, flanking_len, end_5_window_size, TSD_sizes)
    # print(found_TSD, TSD_seq, non_ltr_seq)




    # work_dir = '/tmp/run_hite_main_genome'
    # longest_repeats_flanked_path = '/tmp/run_hite_main_genome/test.fa'
    # threads = 10
    # get_candidate_non_ltr_parallel(longest_repeats_flanked_path, work_dir, threads)









    # # 将candidate_non_ltr.fa进行分类，抽取出LINE/SINE标签
    # work_dir = '/homeb/hukang/KmerRepFinder_test/library/HiTE_lib/mouse_bak'
    # ref_index = 0
    # use_NeuralTE = 0
    # is_wicker = 0
    # #candidate_non_ltr_path = work_dir + '/candidate_non_ltr.fa'
    # candidate_non_ltr_path = work_dir + '/confident_non_ltr.fa'
    # confident_non_ltr_path = candidate_non_ltr_path + '.filter'
    # tmp_output_dir = work_dir
    # threads = 40
    #
    # # candidate_non_ltr_cons = candidate_non_ltr_path + '.cons'
    # # cd_hit_command = 'cd-hit-est -aS ' + str(0.8) + ' -c ' + str(0.8) \
    # #                  + ' -G 0 -g 1 -A 80 -i ' + candidate_non_ltr_path + ' -o ' + candidate_non_ltr_cons + ' -T 0 -M 0'
    # # os.system(cd_hit_command + ' > /dev/null 2>&1')
    #
    # redundant_seq = candidate_non_ltr_path
    # # 使用 blastn 比对， 去掉冗余
    # tmp_blast_dir = work_dir + '/non_ltr_blastn'
    # blastnResults_path = work_dir + '/non_ltr_blastn.out'
    # #blastnResults_path = work_dir + '/test.out'
    # # 1. Start by performing an all-vs-all comparison using blastn.
    # #multi_process_align(redundant_seq, redundant_seq, blastnResults_path, tmp_blast_dir, threads, is_removed_dir=True)
    # # 2. Next, using the FMEA algorithm, bridge across the gaps and link together sequences that can be connected.
    # fixed_extend_base_threshold = 1000
    # longest_repeats = FMEA(blastnResults_path, fixed_extend_base_threshold)
    # # 3. If the combined sequence length constitutes 95% or more of the original individual sequence lengths, we place these two sequences into a cluster.
    # # 3. 如果两条序列是从某一端开始比对的，且overlap占短的序列的60%以上，我们认为短序列应该是长序列的一个truncated copy
    # contigNames, contigs = read_fasta(redundant_seq)
    # keep_clusters = []
    # relations = set()
    # close_end_threshold = 10
    # for query_name in longest_repeats.keys():
    #     longest_repeats_list = longest_repeats[query_name]
    #     for cur_longest_repeat in longest_repeats_list:
    #         close_q_start = False
    #         close_q_end = False
    #         close_s_start = False
    #         close_s_end = False
    #         query_name = cur_longest_repeat[0]
    #         query_len = len(contigs[query_name])
    #         q_len = abs(cur_longest_repeat[2] - cur_longest_repeat[1])
    #         q_start = min(cur_longest_repeat[1], cur_longest_repeat[2])
    #         q_end = max(cur_longest_repeat[1], cur_longest_repeat[2])
    #         if abs(q_start - 0) < close_end_threshold:
    #             close_q_start = True
    #         if abs(q_end - query_len) < close_end_threshold:
    #             close_q_end = True
    #         subject_name = cur_longest_repeat[3]
    #         subject_len = len(contigs[subject_name])
    #         s_len = abs(cur_longest_repeat[4] - cur_longest_repeat[5])
    #         s_start = min(cur_longest_repeat[4], cur_longest_repeat[5])
    #         s_end = max(cur_longest_repeat[4], cur_longest_repeat[5])
    #         if abs(s_start - 0) < close_end_threshold:
    #             close_s_start = True
    #         if abs(s_end - subject_len) < close_end_threshold:
    #             close_s_end = True
    #
    #         is_full_length_copy = False
    #         if float(q_len) / query_len >= 0.95 and float(s_len) / subject_len >= 0.95:
    #             is_full_length_copy = True
    #         # is_truncated_copy = False
    #         # if ((close_q_start or close_q_end) and (close_s_start or close_s_end)
    #         #         and float(q_len) / query_len >= 0.8 or float(s_len) / subject_len >= 0.8):
    #         #     is_truncated_copy = True
    #         if is_full_length_copy:# or is_truncated_copy:
    #             # we consider the query and subject to be from the same family.
    #             if (query_name, subject_name) in relations:
    #                 continue
    #             relations.add((query_name, subject_name))
    #             relations.add((subject_name, query_name))
    #             is_new_cluster = True
    #             for cluster in keep_clusters:
    #                 if query_name in cluster or subject_name in cluster:
    #                     is_new_cluster = False
    #                     cluster.add(query_name)
    #                     cluster.add(subject_name)
    #                     break
    #             if is_new_cluster:
    #                 new_cluster = set()
    #                 new_cluster.add(query_name)
    #                 new_cluster.add(subject_name)
    #                 keep_clusters.append(new_cluster)
    #                 # print(keep_clusters)
    # # Iterate through each cluster, if any element in the cluster overlaps with elements in other clusters, merge the clusters.
    # merged_clusters = []
    # while keep_clusters:
    #     current_cluster = keep_clusters.pop(0)
    #     for other_cluster in keep_clusters[:]:
    #         if current_cluster.intersection(other_cluster):
    #             current_cluster.update(other_cluster)
    #             keep_clusters.remove(other_cluster)
    #     merged_clusters.append(current_cluster)
    # keep_clusters = merged_clusters
    # # store cluster
    # # 挑选cluster中最长的序列当作代表性序列
    # representative_non_ltr = {}
    # representative_non_ltr_path = work_dir + '/representative_non_ltr.fa'
    # all_unique_name = set()
    # raw_cluster_files = []
    # cluster_dir = work_dir + '/raw_ltr_cluster'
    # if not os.path.exists(cluster_dir):
    #     os.makedirs(cluster_dir)
    # for cluster_id, cur_cluster in enumerate(keep_clusters):
    #     cur_cluster_path = cluster_dir + '/' + str(cluster_id) + '.fa'
    #     max_len = 80
    #     longest_non_ltr = None
    #     cur_cluster_contigs = {}
    #     for non_ltr_name in cur_cluster:
    #         seq = contigs[non_ltr_name]
    #         cur_cluster_contigs[non_ltr_name] = seq
    #         all_unique_name.add(non_ltr_name)
    #         if len(seq) > max_len:
    #             longest_non_ltr = non_ltr_name
    #     if longest_non_ltr is not None:
    #         representative_non_ltr[longest_non_ltr] = contigs[longest_non_ltr]
    #     store_fasta(cur_cluster_contigs, cur_cluster_path)
    #     raw_cluster_files.append((cluster_id, cur_cluster_path))
    # # We save the sequences that did not appear in any clusters separately. These sequences do not require clustering.
    # uncluster_path = work_dir + '/uncluster_ltr.fa'
    # uncluster_contigs = {}
    # for name in contigNames:
    #     if name not in all_unique_name:
    #         seq = contigs[name]
    #         uncluster_contigs[name] = seq
    #         if len(seq) > 80:
    #             representative_non_ltr[name] = seq
    # store_fasta(uncluster_contigs, uncluster_path)
    # store_fasta(representative_non_ltr, representative_non_ltr_path)


    # # classify intact-LTRs
    # if use_NeuralTE:
    #     # classify LTR using NeuralTE
    #     NeuralTE_home = '/home/hukang/test/HiTE/bin/NeuralTE'
    #     NeuralTE_output_dir = tmp_output_dir + '/NeuralTE'
    #     if not os.path.exists(NeuralTE_output_dir):
    #         os.makedirs(NeuralTE_output_dir)
    #     NeuralTE_command = 'python ' + NeuralTE_home + '/src/Classifier.py --data ' + candidate_non_ltr_path \
    #                        + ' --use_TSD 0 --model_path ' \
    #                        + NeuralTE_home + '/models/NeuralTE_model.h5 --out_dir ' \
    #                        + NeuralTE_output_dir + ' --thread ' + str(threads) + ' --is_wicker ' + str(is_wicker)
    #     os.system(NeuralTE_command + ' > /dev/null 2>&1')
    #     classified_TE_path = NeuralTE_output_dir + '/classified_TE.fa'
    # else:
    #     TEClass_home = '/home/hukang/test/HiTE/classification'
    #     # classify LTR using RepeatClassifier
    #     sample_name = 'test'
    #     TEClass_command = 'python ' + TEClass_home + '/TEClass_parallel.py --sample_name ' + sample_name \
    #                       + ' --consensus ' + candidate_non_ltr_path + ' --genome 1' \
    #                       + ' --thread_num ' + str(threads) + ' --split_num ' + str(threads) + ' -o ' + tmp_output_dir
    #     os.system(TEClass_command)
    #     classified_TE_path = candidate_non_ltr_path + '.classified'
    # # merge classified LTRs and TEs
    # os.system('cp ' + classified_TE_path + ' ' + confident_non_ltr_path)
    # contignames, contigs = read_fasta(confident_non_ltr_path)
    # non_ltr_contigs = {}
    # for name in contignames:
    #     label = name.split('#')[1]
    #     if 'LINE' in label or 'SINE' in label:
    #         non_ltr_contigs[name] = contigs[name]
    # store_fasta(non_ltr_contigs, confident_non_ltr_path)



    # # 获取zebrafinch lib
    # work_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/Repbase_28.06/zebrafinch'
    # lib_path = work_dir + '/zebrafinch.lib'
    # lib_contigs = {}
    # for filename in os.listdir(work_dir):
    #     file_path = work_dir + '/' + filename
    #     if os.path.isfile(file_path):
    #         contignames, contigs = read_fasta_v1(file_path)
    #         for name in contignames:
    #             if 'Taeniopygia' in name:
    #                 lib_contigs[name] = contigs[name]
    # store_fasta(lib_contigs, lib_path)

    # # 将zebrafinch .out文件中的 Unspecified 换成具体类别
    # labels = set()
    # names, contigs = read_fasta_v1(lib_path)
    # te_labels = {}
    # for name in names:
    #     parts = name.split('\t')
    #     raw_name = parts[0]
    #     label = parts[1]
    #     labels.add(label)
    #     if label.__contains__('ERV'):
    #         label = 'LTR/ERV'
    #     elif label.__contains__('R2'):
    #         label = 'LINE/R2'
    #     else:
    #         label = 'Unknown'
    #     te_labels[raw_name] = label
    # print(labels)
    # out_path = work_dir + '/zebrafinch.out'
    # new_out_path = work_dir + '/zebrafinch.new.out'
    # with open(new_out_path, 'w') as f_save:
    #     with open(out_path, 'r') as f_r:
    #         for line in f_r:
    #             if 'Unspecified' in line:
    #                 parts = line.split(' ')
    #                 te_name = None
    #                 for p in parts:
    #                     if p.strip() == '':
    #                         continue
    #                     if p == 'Unspecified':
    #                         break
    #                     te_name = p
    #                 print(te_name)
    #                 label = te_labels[te_name]
    #                 line = line.replace('Unspecified', label)
    #             f_save.write(line)







    # 获取maize T2T中全长的CRM注释
    # work_dir = '/home/hukang/test/HiTE/library/test'
    # lib_path = work_dir + '/CRM.ref'
    # gff_file = work_dir + '/Zm-Mo17-REFERENCE-CAU-2.0.fa.RMout.gff'
    # contigNames, contigs = read_fasta(lib_path)
    # with open(gff_file, 'r') as f_r:
    #     for line in f_r:
    #         if 'CRM' not in line:
    #             continue
    #         info_parts = line.split('\t')[8].split(' ')
    #         te_name = info_parts[1].replace('"', '').split(':')[1]
    #         te_start = int(info_parts[2])
    #         te_end = int(info_parts[3])
    #         cur_te_len = abs(te_end-te_start)
    #         new_te_name = te_name + '#LTR/CRM'
    #         te_len = len(contigs[new_te_name])
    #         if cur_te_len/te_len >= 0.95:
    #             print(line)

    # # 将chr_1转为chr1
    # file_path = work_dir + '/genome.rename.fa.pass.list.gff3'
    # new_file_path = work_dir + '/genome.rename.fa.pass.list.new.gff3'
    # with open(new_file_path, 'w') as f_save:
    #     with open(file_path, 'r') as f_r:
    #         for line in f_r:
    #             if not line.startswith('#'):
    #                 parts = line.split('\t')
    #                 raw_chr_name = parts[0]
    #                 chr_name = 'chr'+ str(int(parts[0].split('_')[1])+1)
    #                 line = line.replace(raw_chr_name, chr_name)
    #             f_save.write(line)


    # # 输出染色体长度
    # reference = '/home/hukang/EDTA/krf_test/maize/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.rename.fna'
    # contignames, contigs = read_fasta(reference)
    # for name in contignames:
    #     print(name+','+str(len(contigs[name])))

    # work_dir = '/home/hukang/test/HiTE/demo/test2'
    # ltr_list = work_dir + '/genome.rename.fa.pass.list'
    # genome_path = work_dir + '/genome.rename.fa'
    # intact_LTR_path = work_dir + '/intact_LTR.fa'
    #
    # segmentLTRs = {}
    # ltr_contigs = {}
    # ref_contigNames, ref_contigs = read_fasta(genome_path)
    # ltr_index = 1
    # with open(ltr_list, 'r') as f_r:
    #     for line in f_r:
    #         if line.startswith('#'):
    #             continue
    #         direct = line.split('\t')[8]
    #         ltr_name = line.split('\t')[0]
    #         internal_info = line.split('\t')[6]
    #         int_start = int(internal_info.split(':')[1].split('..')[0])
    #         int_end = int(internal_info.split(':')[1].split('..')[1])
    #         parts = ltr_name.split(':')
    #         chr_name = parts[0]
    #         chr_start = int(parts[1].split('..')[0])
    #         chr_end = int(parts[1].split('..')[1])
    #         if direct == '-':
    #             temp = chr_end
    #             chr_end = chr_start
    #             chr_start = temp
    #         chr_seq = ref_contigs[chr_name]
    #         ltr_seq = chr_seq[chr_start-1: chr_end]
    #
    #         terminal_seq1 = chr_seq[chr_start-1: int_start-1]
    #         internal_seq = chr_seq[int_start-1: int_end]
    #         terminal_seq2 = chr_seq[int_end: chr_end]
    #         current_name = 'ltr_'+str(ltr_index)
    #         segmentLTRs[current_name+'_LTR'] = terminal_seq1
    #         segmentLTRs[current_name + '_INT'] = internal_seq
    #         ltr_contigs[current_name] = ltr_seq
    #         ltr_index += 1
    #         # print(terminal_seq1)
    #         # print(internal_seq)
    #         # print(terminal_seq2)
    #         # print(ltr_seq)
    #         break
    # store_fasta(ltr_contigs, intact_LTR_path)




    # repbase_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/repbase'
    # tmp_out_dir = repbase_dir + '/rice'
    # ltr_repbase_path = tmp_out_dir + '/ltr.repbase.ref'
    # tir_repbase_path = tmp_out_dir + '/tir.repbase.ref'
    #log = Logger(tmp_output_dir+'/HiTE.log', level='debug')

    #将repbase中的TIR序列，用最新的过滤方法，看它认为哪些是假阳性
    # plant = 1
    #TE_type = 'TIR'
    # tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/tir_test'
    #raw_input = tmp_dir + '/tir.repbase.ref'
    #output = tmp_dir + '/real_lost_tirs.fa'
    # output = tmp_dir + '/real_tir.fa'
    #output = tmp_dir + '/real_helitron.fa'
    # member_script_path = '/home/hukang/TE_ManAnnot/bin/make_fasta_from_blast.sh'
    # subset_script_path = '/home/hukang/TE_ManAnnot/bin/ready_for_MSA.sh'
    # reference = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test2/rice/genome.rename.fa'
    # member_script_path = '/public/home/hpc194701009/HiTE/tools/make_fasta_from_blast.sh'
    # subset_script_path = '/public/home/hpc194701009/HiTE/tools/ready_for_MSA.sh'
    #reference = '/home/hukang/EDTA/krf_test/dmel/dmel-all-chromosome-r5.43.rename.fasta'
    #reference = '/home/hukang/EDTA/krf_test/drerio/GCF_000002035.6_GRCz11_genomic.rename.fna'
    #reference = '/home/hukang/EDTA/krf_test/maize/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.rename.fna'
    #reference = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test2/rice/GCF_001433935.1_IRGSP-1.0_genomic.fna'
    #reference = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test2/rice_v7/all.chrs.con'
    #reference = '/home/hukang/EDTA/krf_test/ath/GCF_000001735.4_TAIR10.1_genomic.rename.fna'
    #reference = '/home/hukang/HiTE/demo/genome.fa'
    #reference = '/home/hukang/EDTA/krf_test/Zea_mays_Mol7/Zm-Mo17-REFERENCE-CAU-2.0.fa'
    #reference = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Zea_mays_Mol7/Zm-Mo17-REFERENCE-CAU-2.0.fa'
    # reference = '/public/home/hpc194701009/KmerRepFinder_test/genome/GCF_000226075.1_SolTub_3.0_genomic.fna'
    # temp_dir = tmp_dir + '/copies'
    # threads = 40
    # filter_boundary_homo(raw_input, output, reference, member_script_path, subset_script_path, temp_dir, threads, plant,
    #                     TE_type)

    # 我想尝试一下把获得拷贝的方法换成member_script，过滤方法还是老的过滤方法
    #raw_input = tmp_dir + '/fake_helitron.fa'
    #raw_input = tmp_dir + '/lost_tir.fa'
    # raw_input = tmp_dir + '/test.fa'
    #raw_input = tmp_dir + '/test_helitron.fa'
    #raw_input = tmp_dir + '/helitron.repbase.ref'
    #raw_input = tmp_dir + '/tir.repbase.ref'
    #raw_input = tmp_dir + '/novel_tir_repbase.fa'
    #raw_input = tmp_dir + '/library_novel_tir.fa.itr'
    #raw_input = tmp_dir + '/library_novel_helitron.fa'
    # #去掉其中的分类
    # rename_contigs = {}
    # contigNames, contigs = read_fasta(raw_input)
    # for name in contigNames:
    #     new_name = name.split('#')[0]
    #     rename_contigs[new_name] = contigs[name]
    # store_fasta(rename_contigs, raw_input)

    # flanking_len = 50
    # similar_ratio = 0.2
    # TE_type = 'tir'
    #TE_type = 'helitron'
    # ref_index = 0
    #log = Logger(tmp_dir+'/HiTE.log', level='debug')
    #confident_copies = flank_region_align_v2(raw_input, flanking_len, similar_ratio, reference, TE_type, tmp_dir, threads, ref_index, log, member_script_path, subset_script_path, plant)
    # debug = 1
    #output = tmp_dir + '/real_tir_'+str(ref_index)+'.fa'
    #output = tmp_dir + '/confident_helitron_' + str(ref_index) + '.r1.fa'
    #output1 = tmp_dir + '/confident_helitron_' + str(ref_index) + '.fa'
    # result_type = 'cons'
    #flank_region_align_v4(raw_input, output, flanking_len, similar_ratio, reference, TE_type, tmp_dir, threads, ref_index, log, member_script_path, subset_script_path, plant, debug, 0, result_type)


    # confident_helitron_path = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test2/rice/candidate_helitron_0.cons.fa'
    # rename_fasta(confident_helitron_path, confident_helitron_path, 'Helitron')

    # #从present.all.families和perfect.families中，提取非perfect序列，然后看这些序列是否是真的未找到。
    # work_dir = '/homeb/hukang/KmerRepFinder_test/library/get_family_summary_test'
    # all_present = '/homeb/hukang/KmerRepFinder_test/library/tir_test/real_test.fa.itr'
    # perfect = work_dir + '/perfect.families'
    #
    # perfect_set = set()
    # contignames, contigs = read_fasta(all_present)
    # all_present_set = set(contignames)
    # with open(perfect, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         perfect_set.add(line)
    # lost_tirs = all_present_set.difference(perfect_set)
    # print(lost_tirs)
    # print(len(lost_tirs))


    # #读取novel TIR的拷贝数列，统计拷贝数分布
    # # 读取 Excel 表格
    # df = pd.read_excel('/homeb/hukang/KmerRepFinder_test/library/all_tools_run_lib2/rice/novel_tir/data.xlsx', engine='openpyxl')
    #
    # # 统计 "terminal type" 列中每个不同值的数量
    # counts = df['terminal type'].value_counts()
    #
    # # 输出 "novel_terminal" 和 "known_terminal" 的数量
    # print('novel_terminal:', counts['novel_terminal'])
    # print('known_terminal:', counts['known_terminal'])
    #
    # #使用列名访问数据
    # column_data = df['copy num']
    #
    # column_data.to_csv('/homeb/hukang/KmerRepFinder_test/library/nextflow_test2/rice/novel_tir/data.csv', index=False)
    #
    # draw_dist('/homeb/hukang/KmerRepFinder_test/library/nextflow_test2/rice/novel_tir/data.csv')

    # 获取新的TIR转座子，得到它们的多序列比对，蛋白质结构信息
    # tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test4/rice'
    #analyze_new_TIRs(tmp_output_dir)


    # #目前已经能够包含绝大多数的真实TIR，接下来我们需要过滤掉绝大多数的假阳性
    # #我们将confident_tir_0.fa评测RM2结果，大致知道哪些是假阳性的；然后将其用我们的过滤方法，看能否过滤掉这些假阳性
    # #1.收集假阳性序列
    # raw_tir_path = '/homeb/hukang/KmerRepFinder_test/library/helitron_test/confident_helitron_0.fa'
    # file_path = '/homeb/hukang/KmerRepFinder_test/library/get_family_summary_test/file_final.0.1.txt'
    # real_tir_names = set()
    #
    # #统计那些完全没有在金标准中出现的序列
    # with open(file_path, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.split('\t')
    #         query_name = parts[4]
    #         # c1 = float(parts[3])
    #         # c2 = float(parts[6])
    #         # if c1 >= 0.95 and c2 >= 0.95:
    #         #     real_tir_names.add(query_name)
    #         real_tir_names.add(query_name)
    #         #total_names.add(query_name)
    # raw_names, raw_contigs = read_fasta(raw_tir_path)
    # total_names = set(raw_names)
    # fake_tir_names = total_names.difference(real_tir_names)
    # print(fake_tir_names)
    # print(len(fake_tir_names))
    #
    # fake_tir_path = tmp_dir + '/fake_helitron.fa'
    # fake_tirs = {}
    # for name in fake_tir_names:
    #     seq = raw_contigs[name]
    #     fake_tirs[name] = seq
    # store_fasta(fake_tirs, fake_tir_path)






    # lib = '/home/hukang/anaconda2/envs/HiTE1/share/RepeatMasker/Libraries/RepeatPeps.lib'
    # output_table = tmp_dir + '/test.domain.info'
    # thread = 10
    # temp_dir = tmp_dir + '/domain'
    # get_domain_info(cons, lib, output_table, threads, temp_dir)



   
    #分组散点图
    #draw_stripplot()
    # type = 'Copia'
    # output_path = generate_insertion_time(type)
    #
    # #金字塔图
    # darw_barplot(output_path)

    # orig_miu = 1.3e-8
    # miu = 3.3e-8
    # change_LTR_insertion_time(orig_miu, miu)

    #tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library/all_tools_run_lib/rice_v7/HiTE'
    #generate_zebrafish_repbases()
    #generate_repbases()
    # generate_rm2()
    #generate_EDTA()
    # generate_HiTE()
    # input = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Solanum_tuberosum/GCF_000226075.1_SolTub_3.0_genomic.fna'
    # output = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Solanum_tuberosum/GCF_000226075.1_SolTub_3.0_genomic.rename.fna'
    # rename_reference(input, output)

    #测试LTR_finder结果
    # test_LTR_finder()
    # test_LTR_harvest()

    #get_LTR_seqs()

    #identify_new_TIR(tmp_output_dir)

    #parse_novel_tir_locations(tmp_output_dir)

    #draw_dist()

    #reduce_library_size()

    #generate_seq_logos(tmp_output_dir)

    # tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library/all_tools_run_lib2/rice'
    # tmp_output_dir1 = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test2/rice'
    # paths = [tmp_output_dir + '/confident_other.fa', '/home/hukang/HiTE/library/non_LTR.lib']
    # labels = ['HiTE-Non-LTR', 'Non-LTR']
    # my_pal = {"HiTE-Non-LTR": "#4497B1", "Non-LTR": "#F7B92E"}
    # output_path = tmp_output_dir + '/non_ltr_length_dist.txt'

    # paths = [tmp_output_dir+'/confident_tir.fa', tmp_output_dir1+'/tir_tsd_0.cons.fa']
    # labels = ['HiTE-TIR', 'HiTE-TIR-NoFiltering']
    # my_pal = {"HiTE-TIR": "#4497B1", "HiTE-TIR-NoFiltering": "#F7B92E"}
    # output_path = tmp_output_dir + '/tir_length_dist.txt'

    # paths = [tmp_output_dir + '/confident_TE.cons.fa.classified', tmp_output_dir1 + '/longest_repeats_0.cons.rename.fa']
    # labels = ['HiTE', 'HiTE-FMEA']
    # my_pal = {"HiTE": "#4497B1", "HiTE-FMEA": "#F7B92E"}
    # output_path = tmp_output_dir + '/TE_length_dist.txt'

    # paths = [tmp_output_dir + '/confident_helitron.fa', tmp_output_dir + '/candidate_helitron_0.fa']
    # labels = ['HiTE-Helitron', 'HiTE-Helitron-NoFiltering']
    # my_pal = {"HiTE-Helitron": "#4497B1", "HiTE-Helitron-NoFiltering": "#F7B92E"}
    # output_path = tmp_output_dir + '/helitron_length_dist.txt'
    # get_length_dist(paths, labels, my_pal, output_path)

    #run_LTR_test()

    # lost_tirs_path = tmp_output_dir + '/test.fa'
    # get_seq_copies(lost_tirs_path, tmp_output_dir)
    # tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library/all_tools_run_lib/rice_v7/HiTE'
    # sMITE_path = tmp_output_dir + '/sMITE.copies.fa'
    # Novel_TIR_Ghd2_path = tmp_output_dir + '/Novel_TIR_Ghd2.copies.fa'
    # dist_path = tmp_output_dir + '/MITE_dist.txt'
    # generate_MITE_identity_dist(sMITE_path, Novel_TIR_Ghd2_path, tmp_output_dir, dist_path)
    # dist_path = tmp_output_dir + '/MITE_dist.txt'
    # my_pal = {"sMITE": "#16499D", "Novel_TIR_Ghd2": "#E71F19"}
    # draw_violin(dist_path, my_pal)

    # fixed_extend_base_threshold = 1000
    # max_single_repeat_len = 3000000000
    # debug = 1
    # output_dir = '/home/hukang/HiTE/demo/test'
    # repeats_path = (output_dir+'/genome.cut0.fa', [output_dir+'/genome.cut0.fa'], output_dir+'/genome.cut0.fa', '')
    # get_longest_repeats_v4(repeats_path, fixed_extend_base_threshold, max_single_repeat_len, debug)

    # candidate_helitron_path = tmp_output_dir + '/candidate_helitron_' + str(ref_index) + '.fa'
    # resut_file = candidate_helitron_path
    # longest_repeats_flanked_path = tmp_output_dir + '/longest_repeats_0.flanked.fa'
    # sh_dir = '/home/hukang/HiTE/module'
    # HSDIR = '/home/hukang/repeat_detect_tools/TrainingSet'
    # HSJAR = '/home/hukang/repeat_detect_tools/HelitronScanner/HelitronScanner.jar'
    # # 运行helitronscanner
    # HS_temp_dir = tmp_output_dir + '/HS_temp'
    # if not os.path.exists(HS_temp_dir):
    #     os.makedirs(HS_temp_dir)
    # candidate_helitron_path = tmp_output_dir + '/candidate_helitron_' + str(ref_index) + '.fa'
    # multi_process_helitronscanner(longest_repeats_flanked_path, candidate_helitron_path, sh_dir, HS_temp_dir, HSDIR,
    #                               HSJAR, threads)
    # candidate_helitron_contignames, candidate_helitron_contigs = read_fasta(candidate_helitron_path)
    # os.system('rm -rf ' + HS_temp_dir)

    # # 看一下能否过滤或者合并longest_repeats
    # max_single_repeat_len = 30000
    # output_dir = '/homeb/hukang/KmerRepFinder_test/library/all_tools_run_lib1/rice'
    # longest_repeats_path = output_dir + '/longest_repeats_0.fa'
    # contigName, contigs = read_fasta(longest_repeats_path)
    # r1_longest_repeats_path = output_dir + '/longest_repeats_0.r1.fa'
    # with open(r1_longest_repeats_path, 'w') as f_save:
    #     for name in contigName:
    #         seq = contigs[name]
    #         if len(seq) < max_single_repeat_len:
    #             f_save.write('>'+name+'\n'+seq+'\n')

    #尝试合并
    # output_dir = '/home/hukang/HiTE/demo/test'
    # max_single_repeat_len = 30000
    # r1_longest_repeats_path = output_dir + '/longest_repeats_0.fa'
    # contigNames, contigs = read_fasta(r1_longest_repeats_path)
    # merge_contigNames = process_all_seqs(contigNames)
    # print(len(contigNames), len(merge_contigNames))
    # r2_longest_repeats_path = output_dir + '/longest_repeats_0.r2.fa'
    # with open(r2_longest_repeats_path, 'w') as f_save:
    #     for name in merge_contigNames:
    #         seq = contigs[name]
    #         if len(seq) < max_single_repeat_len:
    #             f_save.write('>'+name+'\n'+seq+'\n')

    #summary_not_perfect_repbase()

    #analyze_potato_libs()

    #filter_repbase_nonTE()

    #分模块测试LTR、TIR、Helitron、non-LTR
    #split_module_test()

    # #下面这段代码是把老版本的TIR结果用itrsearch+short_itr过滤一遍
    # tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test2/potato'
    # log = Logger(tmp_output_dir + '/HiTE.log', level='debug')
    # reference = tmp_output_dir + '/C514_hifiasm_ctg.rename.fasta'
    # ref_index = 1
    # TRsearch_dir = '/home/hukang/HiTE/tools'
    # member_script_path = '/home/hukang/HiTE/tools/make_fasta_from_blast.sh'
    # subset_script_path = '/home/hukang/HiTE/tools/ready_for_MSA.sh'
    # filter_last_round_tir_with_itrsearch(tmp_output_dir, reference, ref_index, TRsearch_dir, member_script_path, subset_script_path)

    # #测试get_tsd函数是否正确工作
    # seq = 'GAAGCAGCGCAAGACCAAGGCCAGACGACGAGAAGCTAAGCGAAATAATAACTTAGGGGGTGTTTGGTTGAGCTGTGGTTTTTAGAAAAGCTGCTGTGAGCTGTGGCCTGTGGAAAAGCAGCTGTTGAAAATCAGCTGTGAAAAAAGCAGAAGCCCGTTTGGTTAGAGCAGCTGTGAAACTGTAGGTTGTAGAAGAAATACCTGTAATACCCCTAATGGCATGTGGGACTGTTTATATGCAAATTACTCATACAAACCTAATATGTTCAATATTTTGTATTTGAATGCCCAATTTTTAAAGAAAAAATGTATAATTTTTTCATCAATCAAAGTTTTTCATCAACATATATTAGTTTTTTATCCATCAAAGTGATATGGTCCACATAAACAGAGCAAAAGATAAATGGTTCATACAATGATGCCATATCTAAAAATAACTATCCTCCTCTTGCACTAATCAAAGCATCAGCTATCCTAGTACGAATAGTATTCATGGTATCATCATTCTCTGCATCTTCACTCTCATCGTCTGACGTGTCACTTGTATCTTCTATTAAGTAGTCCTCATGAGCATCACATCTCTCAAACTCCTTATCACGTAAATTGCTATCACGGACGAAGTTATGCATAGCAAGACAAGCCATGATGATATGCTTCTGAGTAGTAGGTGAGAAACTTGGGATGCCCTTCAAAATACGCCACTTTTGCTTCAAGACTCCAAAAGACCTTTCAATAACATTTCGAAGGGAGGAATGCAAGAAATTGAACACCTCATACTTTCCTTGAGGAGGTCCCGAACGATGCCGAAATTCTGGTATATGGTATGTGGTTCCTTTGAAAGGTGCAGGATACCCTATTCGGTTTGGATAACCTGAGTCGACGAGATAATATTTTCCTATAAAACTAAGTTAGTAATGCAAATAATAATACTAGTGTACAAATATTAAAATTAATGAAGGAACCATTTATACCTGTAGGAGGCATGGGAAATGAAGGAAAATTAGCCAATGCATGACTGAGGATGCGCGAGTCATGTGCAGAACCTGGCCAACCAGCAACCGCAAAGATAAACCTCATGTCAAAGTCACAAATAGCTAGCACATTCTGAGATGTATATCCATGACGACATGTGTGGCTAATTGTTTCATCCAATGGCACTACCACTGAAATATGTGAACCATCTATTGCTCCAATCGCATCTTTAAAGTATGGCCAAAATCGATCCTCTCTCAACCTTCCATGCTCCGTAGAAAAAGTAGGATCTCTAGGAGTAATATTGTCTTTAGCTAACTTGCGCAAACACTTCAACACTTCATGAAACTTGGTGTGAACCATCCAAAGTGACCGTGTGAAACGACTTTCAACTTGTGAAAATGATTGTGGCCCTCCGACGATCCACAAAAACATCGTCAATGACTCAATGGATGAAACATTGTTGCTCGAGCTCAATCCGTAGGTAGAGACTAGCAAGTCATGAAGTGCATGAAAAATCTCAGGGCTCATCCGAAACATCTTATAAAAATACCTTGGACGATCCATCAACCTCATCACCCATTGAATTCCAGTTTCTGGCGTTTGTCTGTATGCACCTTTTTGCAAATACCGCTCAGTGTACAAATTGGCCATTTCTCCCTTAATAGAAGCATGTTGGGCAAGGTCAGCTAGCATCAATGATCCCTTTTGCCCTTCTTTGGCTACATCTACTAATTTTCCATCCATCTACAGAAGATATAGTAACAATTAACGGCAAACAAGTACAACAACAAGTACAGAACATAACATCACATAGGGAAGTAAATATTAGCTCACAACATGGAAATAACGGTTTTACATGTCCAAGCAAATATTGGTTTACAACATAGAAATAAATAACTCAAAGTGTTCTCACATAGGGAAGCAAATATCACATAAACAAGGCAACAATCGATAAAAATAAAATTTTAATGTAGGGAATCACTAACGCTTCATCATCTCCTTTTCATGCTCCCTCTCGATCAAATCGAGTCTTCCTTCATTTGTTTCTAGGTAGCTAAAGACCTCCCTAAATTCAGGCTTCATAATCAAAGTTGTGGCTGTGTGCATTAGAGCAGTTTTTTCTTCTACCCCACAATCTTTCATCATCTTCATAGCCTCTTTAATGGTAGGGACTACATTGCATGGAGAAGCTGACGATGCTTCAACACTTGTACAAATCTTTTCTGCCGCCGTTTCTATCTTCAAACAAGTGTTCTTGTACATACGAGCGAATGGGCTCTTCTCATCCTTATCTCGAACAATAGTGGACGACCCCTTGCGTTTATTTTTTTCGGGGTTCACTGTCTCTGCAGGTTTTGGTACCTCAGCCACATCATCATCACTTGCCTCATCGGAGGATACATTTCCAGGACATGATGCTGAAGCCCCACTAACATGTGACTTGCCAAAAATGATATGCAAATCATCTAGGAATTTTGGTCCTTGTTTGCGAAACTTCACATGATCGCATTTGATTCCTTTCTCACGATTGTTGCATCTCTGAAATTAAATAAATAAAGTTAGCAAGGCACACTACAACAAAGAATAGAAAAAAAATATCACAGGGTGAAAGAATAGTACTTACATCAAGGTGTTCTTCCCACCATTCCGGCGAAAAGACAATAGTTTGTTTTGCCTCATCCCATCCAAGACCAGTGGCACAATTCTTGAACTCCATAAACATTGTATATTCCTTTTTCAACACATCTAACTTGTTCTTCAATTGTTTTTTGGTTCGCTCGTCACCCGATTTTTCTGCAAACTTGGAAACAACATTTTTCCAACCCGTAGTGGTGAAAATTCCCATGGGTCTATTCCCAGCTTCTATCTCTTCTTTGCAAACATCAAACAAATGCCTTACAAAACCATCACCCCAATCTGCCTTGTCACCCATTCTTGAACCAAATAAATTTTTGTTCAATCAATTTGCACGCCAAAAATCAATTTGCTGTAACATTCTAGGGGCAACACATGCATCAACACAAGAGCAGGAGCATTCGGCCTCACCTTTCATGCCGGCCGGAGCAGTGGCGAAGGAGTCGTAGTAGCGTCGGGCGGCGTGGGGCGGCGTCGGGCAGGGGTGGCTTTGGGCAGCCTGGGGCAGCGTCGGCCGGCGTCGGGGATCATCGGCCGGCATGGGGCGGCGTCGGCCGGCGTCGGGCAGGCAACGAGGAGACGGCCGTCGGCCGTCGGGCGACCTGGGGCGACGTCGGGTGGCGTGGGGGCGGGCGACGGGGAGACGGGCGATGGAGCAAGCCGGCACCGGGAGATGCGGCGCGCGACGGCGGGACCACAGGCGGGAGCACCCGGCGGTGGGAGCACGGCGGCGCGCGCTGGCGGCGGCAAACCCTAGCGCCGCGTCGAACCCATGGCCAGCGCTTGGCCTCGTGAAAGGGGAACGGGAAAGAAAGAAAGAGCGAATCGTTTTCTTCCGTAATGGCATATGTGGGTAATTTTTTCCCACAAAAGCCAGTGAAAGCAGGGGTAGGGGCTGCTTTTGGATTTGTACTACACTAAAAGCAGCTTTGGAGCAAAAGCACTTGTGGCTTTTGGGCCCTTTGGTTAGCTTTTGGCTTTTGCAAAAGCAAAAGCAGGTTGAAAAGCCCAACCAAACACACCCTTAATAAGTGAGCTTAGCGAAGCGACATAGGACGAAAGTCGGACGACGTGAAGCTAGGTGAAACAACGGACTGAACTACTGAAGCACCAGGTGGCGAAAGCCCAGCAACGCAAGGCGAAGCAGGCTGCAGCCAGTGGCGAATCTAGAAAATCGA'
    # flanking_len = 150
    # tir_start = flanking_len + 1
    # tir_end = len(seq) - flanking_len
    # tsd_search_distance = flanking_len
    # query_name = 'Chr2:4549470-4553257'
    # plant = 1
    # cur_itr_contigs = search_confident_tir_v4(seq, tir_start, tir_end, tsd_search_distance, query_name, plant)
    # print(cur_itr_contigs)

    #Maize Mol7 centromeres experiment
    #maize_centromeres()
    #maize_centromeres_v1()

    #maize_Mol7()
    #get_TIR_internal()
    #maize_length_compare()

    # #获取拷贝目录下，每条序列的拷贝数量分布
    # work_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/tir_test'
    # paths = [work_dir + '/real_tir.fa', work_dir + '/real_helitron.fa']
    #
    # copies_dirs = [work_dir + '/tir_copies_0_0', work_dir+'/helitron_copies_0_0']
    # labels = ['Maize Library Confident TIR', 'Maize Library Confident Helitron']
    # output_path = '/public/home/hpc194701009/KmerRepFinder_test/library/tir_test/tir_copies.txt'
    # output_fig = '/public/home/hpc194701009/KmerRepFinder_test/library/tir_test/tir_copies.png'
    # get_copies_dist_boxplot(paths, copies_dirs, labels, output_path, output_fig)

    # # 估计序列的插入时间
    # tir_identity = 0.79
    # miu = 1.3e-8
    # d = 1 - tir_identity
    # K=  -3/4*math.log(1-d*4/3)
    # T = K/(2*miu)
    # print(d, K, T)


    # reference = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Zea_mays_Mol7/Zm-Mo17-REFERENCE-CAU-2.0.fa'
    # refer_centromeres = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Zea_mays_Mol7/Zm00014ba364380.fa'
    # centromeres_regions = [('chr8', 51478257, 51481204)]
    # ref_names, ref_contigs = read_fasta(reference)
    # centromeres_contigs = {}
    # for cr in centromeres_regions:
    #     centromeres_seq = ref_contigs[cr[0]][cr[1]-1: cr[2]]
    #     name = cr[0] + ':' + str(cr[1]) + '-' + str(cr[2])
    #     centromeres_contigs[name] = centromeres_seq
    # store_fasta(centromeres_contigs, refer_centromeres)


    # #分析EDTA中的TIR假阳性
    # #analyze_EDTA_false_positives()
    # #用我们的过滤方法，看去掉假阳性后，是否能获得一个接近HiTE和Repbase的比例
    # #tir_tsd_path = tmp_dir + '/tir.EDTA.ref'
    # candidate_helitron_path = tmp_dir + '/helitron.EDTA.ref'
    # TRsearch_dir = '/public/home/hpc194701009/HiTE/tools'
    # ref_names, ref_contigs = read_fasta(reference)
    # split_ref_dir = tmp_dir + '/ref_chr'
    # os.system('rm -rf ' + split_ref_dir)
    # if not os.path.exists(split_ref_dir):
    #     os.makedirs(split_ref_dir)
    # ref_blocks = split_dict_into_blocks(ref_contigs, 100)
    # for i, block in enumerate(ref_blocks):
    #     chr_path = split_ref_dir + '/ref_block_' + str(i) + '.fa'
    #     store_fasta(block, chr_path)
    #     os.system('makeblastdb -in ' + chr_path + ' -dbtype nucl')
    # is_recover = 1
    # flanking_len = 50
    # similar_ratio = 0.2
    # TE_type = 'tir'
    #
    # # 使用多比对序列+滑动窗口模式过滤掉假阳性序列
    # # 多轮迭代是为了找到更加准确的边界
    # iter_num = 3
    # input_file = tir_tsd_path
    # for i in range(iter_num):
    #     result_type = 'cons'
    #     output_file = tmp_dir + '/confident_tir_' + str(ref_index) + '.r' + str(i) + '.fa'
    #     resut_file = output_file
    #     if not is_recover or not file_exist(resut_file):
    #         flank_region_align_v5(input_file, output_file, flanking_len, reference, split_ref_dir,
    #                               TE_type, tmp_dir, threads,
    #                               ref_index, log, member_script_path, subset_script_path, plant, debug, i, result_type)
    #     input_file = output_file
    #
    # confident_tir_path = tmp_dir + '/confident_tir_' + str(ref_index) + '.r' + str(iter_num - 1) + '.fa'
    # tir_names, tir_contigs = read_fasta(confident_tir_path)
    # for name in tir_names:
    #     seq = tir_contigs[name]
    #     # 去掉序列首尾的TA或AT
    #     while seq.startswith("TATA") or seq.startswith("ATAT"):
    #         seq = seq[4:]
    #     while seq.endswith("TATA") or seq.endswith("ATAT"):
    #         seq = seq[:-4]
    #     tir_contigs[name] = seq
    #
    # store_fasta(tir_contigs, confident_tir_path)
    # all_copies_out, all_copies_log = run_itrsearch(TRsearch_dir, confident_tir_path, tmp_dir)
    # all_copies_out_name, all_copies_out_contigs = read_fasta(all_copies_out)
    # confident_tir_path = tmp_dir + '/confident_tir_' + str(ref_index) + '.r' + str(iter_num - 1) + '.all_tir.fa'
    # store_fasta(all_copies_out_contigs, confident_tir_path)
    # rename_fasta(confident_tir_path, confident_tir_path, 'TIR')
    #
    # confident_tir_cons = confident_tir_path + '.cons.fa'
    # # 生成一致性序列
    # cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + confident_tir_path + ' -o ' + confident_tir_cons + ' -T 0 -M 0'
    # os.system(cd_hit_command)
    #
    # confident_tir_path = tmp_dir + '/confident_tir_' + str(ref_index) + '.fa'
    # rename_fasta(confident_tir_cons, confident_tir_path, 'TIR')

    # # 对HelitronScanner识别的结果使用同源性过滤方法过滤
    # confident_helitron_path = tmp_dir + '/confident_helitron_' + str(ref_index) + '.fa'
    # resut_file = confident_helitron_path
    # if not is_recover or not file_exist(resut_file):
    #     flanking_len = 50
    #     similar_ratio = 0.2
    #     TE_type = 'helitron'
    #     # 多轮迭代是为了找到更加准确的边界
    #     iter_num = 3
    #     input_file = candidate_helitron_path
    #     for i in range(iter_num):
    #         result_type = 'cons'
    #         output_file = tmp_dir + '/confident_helitron_' + str(ref_index) + '.r' + str(i) + '.fa'
    #         resut_file = output_file
    #         if not is_recover or not file_exist(resut_file):
    #             flank_region_align_v5(input_file, output_file, flanking_len, reference, split_ref_dir,
    #                                   TE_type,
    #                                   tmp_dir, threads,
    #                                   ref_index, log, member_script_path, subset_script_path, 1, debug, i, result_type)
    #         input_file = output_file
    #     cur_confident_helitron_path = tmp_dir + '/confident_helitron_' + str(ref_index) + '.r' + str(
    #         iter_num - 1) + '.fa'
    #     cur_confident_helitron_cons = tmp_dir + '/confident_helitron_' + str(ref_index) + '.r' + str(
    #         iter_num - 1) + '.cons.fa'
    #     # 生成一致性序列
    #     cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                      + ' -G 0 -g 1 -A 80 -i ' + cur_confident_helitron_path + ' -o ' + cur_confident_helitron_cons + ' -T 0 -M 0'
    #     os.system(cd_hit_command)
    #     rename_fasta(cur_confident_helitron_cons, confident_helitron_path, 'Helitron')


    #统计玉米中所有TIR转座子的插入时间
    #analyz_maize_TIR_insert_time()

    # tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa'
    # # confident_ltr_terminal_path = tmp_output_dir + '/confident_ltr_cut.terminal.fa'
    # # ltr_copies = flanking_and_getcopies(confident_ltr_terminal_path)
    #
    # confident_helitron_path = tmp_output_dir + '/confident_helitron_0.rename.cons.fa'
    # helitron_copies = flanking_and_getcopies(confident_helitron_path)
    #
    # # enspm_path = tmp_output_dir + '/enspm.fa'
    # # confident_TE_path = tmp_output_dir + '/confident_TE.cons.fa.final.classified'
    # # contignames, contigs = read_fasta(confident_TE_path)
    # # node_index = 0
    # # Enspm_contigs = {}
    # # for name in contignames:
    # #     if name.__contains__('#DNA/CMC-EnSpm'):
    # #         seq = contigs[name]
    # #         Enspm_contigs[name] = seq
    # # store_fasta(Enspm_contigs, enspm_path)
    # #
    # # Enspm_copies = flanking_and_getcopies(enspm_path)
    #
    # # 取序列头部前20bp至头部后20bp
    # # ltr_start_logos, ltr_tail_logos = get_logo_seq(ltr_copies)
    # helitron_start_logos, helitron_tail_logos = get_logo_seq(helitron_copies)
    # # Enspm_start_logos, Enspm_tail_logos = get_logo_seq(Enspm_copies)
    #
    # # ltr_terminal_start_path = tmp_output_dir + '/ltr_terminal_start.fa'
    # # ltr_terminal_end_path = tmp_output_dir + '/ltr_terminal_end.fa'
    # # store_fasta(ltr_start_logos, ltr_terminal_start_path)
    # # store_fasta(ltr_tail_logos, ltr_terminal_end_path)
    #
    # helitron_terminal_start_path = tmp_output_dir + '/helitron_terminal_start.fa'
    # helitron_terminal_end_path = tmp_output_dir + '/helitron_terminal_end.fa'
    # store_fasta(helitron_start_logos, helitron_terminal_start_path)
    # store_fasta(helitron_tail_logos, helitron_terminal_end_path)
    #
    # # EnSpm_terminal_start_path = tmp_output_dir + '/EnSpm_terminal_start.fa'
    # # EnSpm_terminal_end_path = tmp_output_dir + '/EnSpm_terminal_end.fa'
    # # store_fasta(Enspm_start_logos, EnSpm_terminal_start_path)
    # # store_fasta(Enspm_tail_logos, EnSpm_terminal_end_path)

