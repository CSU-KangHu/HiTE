#-- coding: UTF-8 --
import argparse
import os
import re
import sys

import codecs

import json
import time

#import regex
#from fuzzysearch import find_near_matches

#import numpy as np
# from matplotlib import pyplot as plt
# import seaborn as sns
# import pandas as pd

import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, read_fasta_v1, rename_fasta, getReverseSequence, allow_mismatch, \
    run_itrsearch, multi_process_itr, filter_large_gap_tirs, multi_process_align_and_get_copies, \
    store_copies_v1, get_TSD, store_copies, store_LTR_seq_v1, store_LTR_seq, store_LTR_seq_v2, rename_reference, \
    run_LTR_harvest, run_LTR_retriever, determine_repeat_boundary_v2, determine_repeat_boundary_v1, multi_process_align, \
    get_copies, TSDsearch_v4, overlap_with_boundary, judge_flank_align, get_copies_v1, convertToUpperCase_v1, \
    determine_repeat_boundary_v3, search_confident_tir, store_copies_seq, PET, multiple_alignment_blastx_v1, store2file


def generate_repbases():
    # 水稻
    repbase_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/repbase'
    repbase_path = repbase_dir + '/zebrep.ref'
    repbase_names, repbase_contigs = read_fasta_v1(repbase_path)
    tags = set()
    for name in repbase_names:
        tag = name.split('\t')[1]
        tags.add(tag)
    print(tags)
    print(len(tags))

    ltr_tags = ['Gypsy', 'Copia', 'LTR Retrotransposon', 'BEL', 'LTR', 'Endogenous Retrovirus']
    tir_tags = ['Mariner/Tc1', 'DNA transposon', 'EnSpm/CACTA', 'MuDR', 'hAT', 'Harbinger', 'Transib', 'piggyBac', 'P', 'DNA', 'Sola2', 'Kolobok', ]
    helitron_tags = ['Helitron', 'MINIME_DN']
    non_ltr_tags = ['L1', 'SINE2/tRNA', 'Non-LTR Retrotransposon', 'SINE', 'R1', 'Jockey', 'CR1', 'R2', 'RTEX', 'Hero', 'RTE']
    tmp_out_dir = repbase_dir + '/dmel'
    ltr_repbase_path = tmp_out_dir + '/ltr.repbase.ref'
    tir_repbase_path = tmp_out_dir + '/tir.repbase.ref'
    helitron_repbase_path = tmp_out_dir + '/helitron.repbase.ref'
    non_ltr_repbase_path = tmp_out_dir + '/non_ltr.repbase.ref'

    ltr_contigs = {}
    tir_contigs = {}
    helitron_contigs = {}
    non_ltr_contigs = {}
    for name in repbase_names:
        tag = name.split('\t')[1]
        if tag in ltr_tags:
            ltr_contigs[name] = repbase_contigs[name]
        elif tag in tir_tags:
            tir_contigs[name] = repbase_contigs[name]
        elif tag in helitron_tags:
            helitron_contigs[name] = repbase_contigs[name]
        elif tag in non_ltr_tags:
            non_ltr_contigs[name] = repbase_contigs[name]
    store_fasta(ltr_contigs, ltr_repbase_path)
    store_fasta(tir_contigs, tir_repbase_path)
    store_fasta(helitron_contigs, helitron_repbase_path)
    store_fasta(non_ltr_contigs, non_ltr_repbase_path)

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

def summary_not_appear_repbase():
    #我们观察，有哪些序列是没有出现在比对中
    align_file = '/public/home/hpc194701009/KmerRepFinder_test/library/get_family_summary_test/file_final.0.1.txt'
    query_names = set()
    with open(align_file, 'r') as f_r:
        for line in f_r:
            query_name = line.split('\t')[4]
            query_names.add(query_name)

    tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/drerio'
    names, contigs = read_fasta(tmp_dir + '/confident_tir.rename.cons.fa')
    names = set(names)
    diff_set = names.difference(query_names)
    print(diff_set)
    print('not appear size: ' + str(len(diff_set)))

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
    confident_tir_path = tmp_output_dir + '/confident_tir.rename.cons.fa'
    repbase_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/repbase/rice'
    tir_repbase_path = repbase_dir + '/tir.repbase.ref'
    confident_tir_consensus = tmp_output_dir + '/confident_tir.copy2.fa'
    total_tir_path = tmp_output_dir + '/total_tir.fa'
    total_tir_consensus = tmp_output_dir + '/total_tir.cons.fa'
    os.system('cat '+tir_repbase_path+' '+confident_tir_consensus+' > '+total_tir_path)

    confident_tir_names, confident_tir_contigs = read_fasta(confident_tir_consensus)

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
                if not name.startswith('N_'):
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
    TRsearch_dir = '/home/hukang/HiTE-2.0.1/ReferenceMode/tools'
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
                if not name.startswith('N_'):
                    has_repbase = True
                    break
            #没有repbase序列，代表这是我们新识别的TIR
            if not has_repbase:
                for name in cluster:
                    new_tir_names.add(name.split('-')[0])
    print('novel TIR with new TIR terminals:')
    #print(new_tir_names)
    print(len(new_tir_names))

    #统计novel TIR的拷贝数量
    threads = 40
    reference = tmp_output_dir + '/GCF_001433935.1_IRGSP-1.0_genomic.fna'
    temp_dir = tmp_output_dir + '/temp'
    blast_program_dir = '/home/hukang/repeat_detect_tools/rmblast-2.9.0-p2'
    all_copies = multi_process_align_and_get_copies(novel_tir_consensus, reference, blast_program_dir,
                                                    temp_dir, 'tir', threads)

    # 在copies的两端 flanking 20bp的序列
    flanking_len = 20
    all_copies, tsd_info = flanking_copies(all_copies, novel_tir_consensus, reference, flanking_len, copy_num=-1)

    multicopy_novel_tirs_num = 0
    # 统计每个query_name拷贝的数量
    # query_name -> (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq)
    query_copy_num = {}
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        query_copy_num[query_name] = len(copies)
        if len(copies) >= 2:
            multicopy_novel_tirs_num += 1
    #print(query_copy_num)
    print('novel TIR with copy number >= 2:')
    print(multicopy_novel_tirs_num)

    query_copy_num_path = tmp_output_dir + '/novel_tir_copies_num.csv'
    # 存储query_copy_num
    with codecs.open(query_copy_num_path, 'w', encoding='utf-8') as f:
        json.dump(query_copy_num, f)

    novel_tir_copies = tmp_output_dir + '/novel_tir_copies.csv'
    # 存储all copies
    with codecs.open(novel_tir_copies, 'w', encoding='utf-8') as f:
        json.dump(all_copies, f)

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

    copy_info_path = tmp_output_dir + '/novel_tir.copies.info'
    store_copies(tsd_info, copy_info_path)

    #计算Venn图指标
    #1.novel tir copy > 3 and novel tirs with new terminals
    #2.novel tir copy > 3 and novel tirs with known terminals
    #3.novel tir copy <=3 and novel tirs with new terminals
    #4.novel tir copy <=3 and novel tirs with known terminals
    query_copy_over_3 = set()
    query_copy_less_3 = set()
    for query_name in query_copy_num.keys():
        copy_num = query_copy_num[query_name]
        if copy_num >= 2:
            query_copy_over_3.add(query_name)
        else:
            query_copy_less_3.add(query_name)

    c1 = 0
    c2 = 0
    c3 = 0
    c4 = 0
    for query_name in new_tir_contigs.keys():
        if query_name in query_copy_over_3 and query_name in new_tir_names:
            c1 += 1
        elif query_name in query_copy_over_3 and query_name not in new_tir_names:
            c2 += 1
        elif query_name in query_copy_less_3 and query_name in new_tir_names:
            c3 += 1
        elif query_name in query_copy_less_3 and query_name not in new_tir_names:
            c4 += 1

    print(c1, c2, c3, c4)
    print(new_tir_names)
    print(new_tir_contigs.keys())
    print(query_copy_over_3)



def draw_dist():
    #tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa'
    tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library/RepeatMasking_test/rice_no_kmer'
    query_copy_num_path = tmp_output_dir + '/novel_tir_copies_num.csv'
    file = open(query_copy_num_path, 'r')
    js = file.read()
    query_copy_num = json.loads(js)
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

    y = list(query_copy_num.values())
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
                if seq_len >= 5000:
                    continue
                f_save.write(str(seq_len)+'\t'+str(label)+'\n')
    df = pd.read_csv(output_path, sep='\t', encoding='utf-8')
    print(df)
    sns.set_context("paper", rc={"font.size":20,"axes.titlesize":8,"axes.labelsize":20})
    b = sns.violinplot(x=df["methods"], y=df["length"], palette=my_pal)
    b.set_xlabel("", fontsize=20)
    b.set_ylabel("Length", fontsize=20)
    plt.tick_params(axis='x', which='major', labelsize=15)
    plt.show()


    # Calculate number of obs per group & median to position labels
    medians = df.groupby(['methods'])['length'].median().values
    print(medians)


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
    plt.show()

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
    sns.stripplot(x="process",y="Execution time (minutes)", data=tips, ax=ax, s=10)
    ax.grid(True)
    plt.xticks(rotation=60)
    plt.tight_layout()
    plt.savefig('/home/hukang/nextflow_runtime.png', format='png')

def darw_barplot(input):
    #sns.set(context="notebook", style='whitegrid', font_scale=1)
    sns.set_style('whitegrid')
    result=pd.read_csv(input)
    #print(result)
    bar_plot = sns.barplot(y = result['time'].unique(), x = -result['number'], color = "DarkSalmon", 
                       data = result, order = result['time'].unique()[::-1],)
    bar_plot = sns.barplot(y = result['time'].unique(), x = result['number'], color = "CadetBlue",
                        data = result, order = result['time'].unique()[::-1],)
    #plt.xticks([-650,-400,-200,0,200,400,650],[650,200,100,0,100,200,250])
    # plt.rcParams['font.sans-serif'] = ['SimHei']
    #plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
    plt.rcParams['axes.unicode_minus'] = True
    bar_plot.set(xlabel="Number of LTR", ylabel="Mya", title = "")
    plt.savefig('/home/hukang/LTR_insert_time.png', format='png')

def generate_insert_time(ltr_file):
    #将LTR_retriever的插入时间统计到文件
    times = []
    with open(ltr_file, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            else:
                time = int(line.split('\t')[11])
                times.append(time)
    f_r.close()
    #print(max(times))

    #按照10W年分组
    u = 1
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

def generate_insertion_time():
    output = '/home/hukang/insert_time.csv'
    ltr_file = '/homeb/hukang/KmerRepFinder_test/library/all_tools_run_lib/rice_v7/HiTE/all.chrs.rename.fa.pass.list'
    speices1 = 'Rice'
    time_group1 = generate_insert_time(ltr_file)
    #print(time_group1)

    ltr_file = '/homeb/hukang/KmerRepFinder_test/library/all_tools_run_lib/rice_v7/HiTE/all.chrs.rename.fa.pass.list'
    speices2 = 'P.abies'
    time_group2 = generate_insert_time(ltr_file)
    #print(time_group2)

    with open(output, 'w') as f_save:
        f_save.write('time,species,number\n')
        keys = sorted(time_group1.keys(), reverse=False)
        for i, g_num in enumerate(keys):
            if i < len(keys)-1:
                next_g_num = str(keys[i+1])
            else:
                next_g_num = ''
            num = time_group1[g_num]
            f_save.write(str(g_num)+'-'+str(next_g_num)+','+speices1+','+str(num)+'\n')
        keys = sorted(time_group2.keys(), reverse=False)
        for i, g_num in enumerate(keys):
            if i < len(keys)-1:
                next_g_num = str(keys[i+1])
            else:
                next_g_num = ''
            num = time_group2[g_num]
            f_save.write(str(g_num)+'-'+str(next_g_num)+','+speices2+','+str(num)+'\n')

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

def run_BM_RM2(TE_path, res_out, tmp_dir):
    threads = 30
    rm2_script = '/homeb/hukang/KmerRepFinder_test/library/get_family_summary_paper.sh'
    lib_path = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/repbase/oryrep.ref'
    rm2_command = 'RepeatMasker -lib '+lib_path+' -nolow -pa '+str(threads)+' '+TE_path
    rm2_res_command = 'cd '+tmp_dir+'/rm2_test && rm -rf * && sh '+rm2_script+ ' ' + TE_path + '.out > ' + res_out
    os.system(rm2_command)
    print(rm2_res_command)
    os.system(rm2_res_command)


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
    # 5.如果边界前的序列为N、2bp TA、4bp TTAA，则序列为真实TIR边界，再经由下面的流程判断未扩展的一致性序列是否具有TIR结构
    # 取一致性序列的前后2bp, 4bp
    cons_first_2bp = cons_seq[start_pos-2: start_pos] 
    cons_last_2bp = cons_seq[end_pos+1: end_pos+3]
    cons_first_4bp = cons_seq[start_pos-4: start_pos] 
    cons_last_4bp = cons_seq[end_pos+1: end_pos+5]
    cons_before_boundary = cons_seq[0: start_pos]
    cons_after_boundary = cons_seq[end_pos+1: ]
    cons_before_2bp_boundary = cons_seq[0: start_pos-2]
    cons_after_2bp_boundary = cons_seq[end_pos+3: ]
    cons_before_4bp_boundary = cons_seq[0: start_pos-4]
    cons_after_4bp_boundary = cons_seq[end_pos+5: ]
    if (is_only_chars(cons_before_boundary, 'N') and is_only_chars(cons_after_boundary, 'N')) or \
        (cons_first_2bp == 'TA' and cons_last_2bp == 'TA' and is_only_chars(cons_before_2bp_boundary, 'N') and is_only_chars(cons_after_2bp_boundary, 'N')) or \
        (cons_first_4bp == 'TTAA' and cons_last_4bp == 'TTAA' and is_only_chars(cons_before_4bp_boundary, 'N') and is_only_chars(cons_after_4bp_boundary, 'N')):
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


def get_domain_info(cons, lib, output_table, threads, temp_dir):
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    # 1. 将cons进行划分，每一块利用blastx -num_threads 1 -evalue 1e-20 进行cons与domain进行比对
    partitions_num = int(threads)
    consensus_contignames, consensus_contigs = read_fasta(cons)
    data_partitions = PET(consensus_contigs.items(), partitions_num)
    merge_distance = 100
    file_list = []
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for partition_index, data_partition in enumerate(data_partitions):
        if len(data_partition) <= 0:
            continue
        cur_consensus_path = temp_dir + '/'+str(partition_index)+'.fa'
        store2file(data_partition, cur_consensus_path)
        cur_output = temp_dir + '/'+str(partition_index)+'.out'
        cur_file = (cur_consensus_path, lib, cur_output)
        job = ex.submit(multiple_alignment_blastx_v1, cur_file, merge_distance)
        jobs.append(job)
    ex.shutdown(wait=True)

    # 2. 生成一个query与domain的最佳比对表
    if os.path.exists(output_table):
        os.remove(output_table)
    for job in as_completed(jobs):
        cur_output = job.result()
        os.system('cat ' + cur_output + ' >> ' + output_table)

if __name__ == '__main__':
    repbase_dir = '/homeb/hukang/KmerRepFinder_test/library/curated_lib/repbase'
    tmp_out_dir = repbase_dir + '/rice'
    ltr_repbase_path = tmp_out_dir + '/ltr.repbase.ref'
    tir_repbase_path = tmp_out_dir + '/tir.repbase.ref'
    tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library/RepeatMasking_test/rice_no_kmer'
    log = Logger(tmp_output_dir+'/HiTE.log', level='debug')
    # #分析
    tmp_dir = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test1/rice'
    # names, contigs = read_fasta(tmp_dir+'/confident_tir_0.fa')
    # new_contigs = {}
    # node_index = 0
    # for name in names:
    #     parts = name.split('-tir_')[1].split('-tsd_')
    #     tir_len = parts[0]
    #     tsd_len = parts[1]
    #     new_name = 'N_'+str(node_index) + '-tir_'+tir_len+'-tsd_'+tsd_len
    #     node_index += 1
    #     new_contigs[new_name] = contigs[name]
    # store_fasta(new_contigs, tmp_dir+'/confident_tir_0.fa')

    # 分析cd-hit聚类，去除嵌合
    #clean_nested_rep()
    # TE_path = tmp_dir + '/confident_tir_0.multi.cons.fa'
    # res_out = tmp_dir + '/multi_res.log'
    # TE_path = tmp_dir + '/pure_tir_rep.fa'
    # res_out = tmp_dir + '/pure_res.log'
    # TE_path = tmp_dir + '/pure_with_clean_tir.fa'
    # res_out = tmp_dir + '/pure_with_clean_res.log'
    # run_BM_RM2(TE_path, res_out, tmp_dir)

    # 将我们的结果使用make_fasta_from_blast.sh，过滤掉单拷贝的序列
    #keep_multi_copy()

    tmp_dir = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test1/rice'
    raw_input = tmp_dir + '/real_tirs.fa'
    temp_copies_dir = tmp_dir + '/copies'
    member_script_path = '/home/hukang/TE_ManAnnot/bin/make_fasta_from_blast.sh'
    subset_script_path = '/home/hukang/TE_ManAnnot/bin/ready_for_MSA.sh'
    reference = tmp_dir + '/GCF_001433935.1_IRGSP-1.0_genomic.rename.fa'
    threads = 40
    temp_nested_dir = tmp_dir+'/nested'
    test_home = '/home/hukang/HiTE/module'
    TRsearch_dir = '/home/hukang/HiTE/tools'
    #cur_output = build_lib(raw_input, reference, threads, temp_copies_dir, temp_nested_dir, test_home, TRsearch_dir, script_path, log)
    #print(cur_output)

    TE_path = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test1/rice/nested/pure_with_clean_tir.fa'
    res_out = tmp_dir + '/pure_with_clean_res.log'
    #run_BM_RM2(TE_path, res_out, tmp_dir)

    # # 1.将具有TIR+TSD结构的序列比对到reference上，取序列成员
    # # 2.遍历成员序列，如果有超过3条成员序列具有TIR结构则认定为真实TIR
    # #raw_input = tmp_dir + '/test.fa'
    # raw_input = tmp_dir + '/tir_tsd_0.filter_tandem.fa'
    # real_tirs = tmp_dir + '/real_tirs.fa'
    # temp_dir = tmp_dir + '/copies'
    # os.system('rm -rf ' + temp_dir)
    # if not os.path.exists(temp_dir):
    #     os.makedirs(temp_dir)
    # names, contigs = read_fasta(raw_input)
    # split_files = []
    # for i, name in enumerate(names):
    #     cur_file = temp_dir + '/' + str(name) + '.fa'
    #     cur_contigs = {}
    #     cur_contigs[name] = contigs[name]
    #     store_fasta(cur_contigs, cur_file)
    #     split_files.append(cur_file)

    # ex = ProcessPoolExecutor(threads)
    # jobs = []
    # for ref_index, cur_file in enumerate(split_files):
    #     job = ex.submit(run_find_members_v1, cur_file, reference, temp_dir, member_script_path, subset_script_path, TRsearch_dir)
    #     jobs.append(job)
    # ex.shutdown(wait=True)

    # true_tirs = {}
    # for job in as_completed(jobs):
    #     cur_name, cur_seq = job.result()
    #     if cur_name is not None:
    #         true_tirs[cur_name] = cur_seq
    # store_fasta(true_tirs, real_tirs)

    # real_tirs_rename = tmp_dir + '/real_tirs.rename.fa'
    # rename_fasta(real_tirs, real_tirs_rename)

    # res_out = tmp_dir + '/real_tirs_res.log'
    # run_BM_RM2(real_tirs_rename, res_out, tmp_dir)

    cons = tmp_dir + '/test.fa'
    # names,contigs = read_fasta(cons)
    # seq = contigs[names[0]]
    # r_seq = getReverseSequence(seq)
    # contigs['N_2'] = r_seq
    # store_fasta(contigs, cons)

    lib = '/home/hukang/anaconda2/envs/HiTE1/share/RepeatMasker/Libraries/RepeatPeps.lib'
    output_table = tmp_dir + '/test.domain.info'
    thread = 10
    temp_dir = tmp_dir + '/domain'
    get_domain_info(cons, lib, output_table, threads, temp_dir)


    # 定义要搜索的字符串
    # string_to_search = ''
    # STR = "nnTCCACAATAGTCCCnGGATGCAATAAAAACCGGGGCnnnnnnnnnTAnnnnnAnnnnnAGATnnnnnnnnnnnnnnnnnGATCTTTnAGTCCCGGTTnCAAAAGGGTAAnCGGGCATATTTGATCnTTTAGnnnnnnnnTCCCGGTTTGTGTnnnnnnnnnnnnnnnnnnnnnnnnnnnnnTACCAACCGGGnnnnnnnnnnnACTAAAGnnnATCnnnnnnnnnnnnnnnnnnAnTnnCTTTAGTCnnnnnnnnnnCCGGTTCGAnATGCTGTCAnnnGGGCCTGTCAGGnnnnnnnCCCCnCnCCnnnnnnnnnGGnnnnnnnnnnnnnnnnnGATnnnnnnnnnnCTTnnnnnnnnTAGTCCCGGnTTGGnnnnnTAnnCACCnAACCGnGGACTAAAGATCnnnnnnnnnnnnnn"
    # # 定义要查找的模式
    # pattern = "CCAACCGGGAC"
    # a = find_near_matches(pattern, STR, max_l_dist=2)
    # for b in reversed(a):
    #     print(b)
    # 使用 fuzzysearch 函数进行模糊匹配
    # match = regex.fuzzysearch(pattern, string_to_search)
    
    # # 输出匹配结果
    # if match:
    #     print(f"找到匹配项：{match.group()}")
    # else:
    #     print("未找到匹配项")



    # #恢复real_tirs.fa
    # real_tirs = tmp_dir + '/real_tirs.fa'
    # cur_dir = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test1/rice/copies'
    # #regex = re.compile(r"-tsd_\d+\.fa$")
    # regex = re.compile(r".final.cons.fa$")
    # if os.path.exists(real_tirs):
    #     os.remove(real_tirs)
    # # 遍历目录中的所有文件名
    # new_contigs = {}
    # for filename in os.listdir(cur_dir):
    #     # 如果文件名符合正则表达式，则输出
    #     if regex.search(filename):
    #         names, contigs = read_fasta(cur_dir+'/'+filename)
    #         if len(names) == 0:
    #             continue
    #         new_contigs[filename] = contigs[names[0]]
    #         print(filename+'\n'+contigs[names[0]]+'\n')
    # store_fasta(new_contigs, real_tirs)
    #         #print(filename)




    # #分析我们丢失了哪些perfect，是不是sub_families
    # dir1 = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test1/rice/rm2_test'
    # dir2 = '/homeb/hukang/KmerRepFinder_test/library/nextflow_test1/rice/rm2_test1'
    # perfect1 = set()
    # perfect2 = set()
    # with open(dir1+'/perfect.families', 'r') as f_r:
    #     for line in f_r:
    #         perfect1.add(line.replace('\n', ''))
    
    # with open(dir2+'/perfect.families', 'r') as f_r:
    #     for line in f_r:
    #         perfect2.add(line.replace('\n', ''))
    
    # s = perfect1.difference(perfect2)
    # print(s)
    # print(len(s))
    # sub_count = 0
    # sub_names = []
    # for name in s:
    #     parts = name.split('-')
    #     if len(parts) > 1:
    #         sub_name = parts[1]
    #         if sub_name.startswith('N'):
    #             sub_count += 1
    #             sub_names.append(name)
    # print(sub_count)
    # print(sub_names)


   
    #分组散点图
    #draw_stripplot()
    
    #generate_insertion_time()
            

    #金字塔图
    #darw_barplot(input)

    #tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library/all_tools_run_lib/rice_v7/HiTE'
    #generate_zebrafish_repbases()

    # copy_info_path = tmp_out_dir + '/tir.repbase.copies.info'
    # copies_list = get_repbase_copies(tir_repbase_path, copy_info_path)
    # print(copies_list)

    #rename_TE_file()

    #rename_EDTA_file()
    #summary_not_appear_repbase()
    #not_found_repbase(tir_repbase_path, copies_list)

    #测试LTR_finder结果
    # test_LTR_finder()
    # test_LTR_harvest()

    #get_LTR_seqs()

    #identify_new_TIR(tmp_output_dir)

    #parse_novel_tir_locations(tmp_output_dir)

    #draw_dist()

    #test_remove_copy1_tirs(tmp_output_dir)

    #test_EAHelitron()

    #reduce_library_size()

    #generate_seq_logos(tmp_output_dir)

    #test_split_genome()
    #summary_tir_5bp(tir_repbase_path)


    # plant = 1
    # query_name = 'N_136771-len_520-ref_Chr12-22058590-22059110'
    # seq = 'TCAGGCAACGCGGTAGCTCCACGTTACTAATTAGCAAAACTAACACGACCAACAAACGAAAACTAAATAGGTACTTTCTCCATTTCACAATGTAAGTTATTCTAGCATTTTGCACATTTATAACAATGTTAATGAATCTAGATAAATATATATGTATTAATATCTAGATTCATTAACATCAATATAAATGTGGGAAATACTAGAATGACTTACATTGTGAAACGGAGGGAGTACATCATTGCTCAGATGACTTGCCTATTATATTTTCTCACTTCATATGGAAGTCCATATAAACTCTCTAAATTGTTACTCCCTCCTTCCCTAAATATTTGACACCGTTAACTTTTTAAATATGTTTGACCGTTCGTTTTATTCAAAAACTTTTGTGATATGTGTAAAACTATATGTATACATAAAAGTATATTTAACAATAAATCAAATGATAGAAAAAGAATTAACAATTACTTAAATTTTTTGAATAAGACGAACGGTCAAACATTTTTAAAAAAGTCAACGGCATCAAACATTTTGGGATGGAGGTAGTATGTACTTGTGTTTTTTCTCCAAAAACTTTACGTGCCAGCATATTGCTTAGATGACTTGTCTATTCTATTCTTTTA'
    # flanking_len = 50
    # tir_start = flanking_len + 1
    # tir_end = len(seq) - flanking_len
    # # 寻找所有可能的TSD序列，计算每条序列的边界与原始边界的距离，并存到header里
    # tsd_search_distance = 50
    # top_itr_contigs = search_confident_tir(seq, tir_start, tir_end, tsd_search_distance, query_name, plant)
    # print(top_itr_contigs)

    paths = [tmp_output_dir + '/confident_other_0.fa', '/home/hukang/HiTE-2.0.1/ReferenceMode/library/non_LTR.lib']
    labels = ['HiTE-Non-LTR', 'Non-LTR']
    my_pal = {"HiTE-Non-LTR": "#4497B1", "Non-LTR": "#F7B92E"}
    output_path = tmp_output_dir + '/non_ltr_length_dist.txt'

    # paths = [tmp_output_dir+'/confident_tir.rename.cons.fa', tmp_output_dir+'/tir_tsd_0.cons.rename.fa']
    # labels = ['HiTE-TIR', 'HiTE-TIR-NoFiltering']
    # my_pal = {"HiTE-TIR": "#4497B1", "HiTE-TIR-NoFiltering": "#F7B92E"}
    # output_path = tmp_output_dir + '/tir_length_dist.txt'

    # paths = [tmp_output_dir + '/confident_TE.cons.fa.final.classified', tmp_output_dir + '/longest_repeats_0.cons.rename.fa']
    # labels = ['HiTE', 'HiTE-FMEA']
    # my_pal = {"HiTE": "#4497B1", "HiTE-FMEA": "#F7B92E"}
    # output_path = tmp_output_dir + '/TE_length_dist.txt'

    # paths = [tmp_output_dir + '/confident_helitron_0.rename.cons.fa', tmp_output_dir + '/candidate_helitron_0.cons.rename.fa']
    # labels = ['HiTE-Helitron', 'HiTE-Helitron-NoFiltering']
    # my_pal = {"HiTE-Helitron": "#4497B1", "HiTE-Helitron-NoFiltering": "#F7B92E"}
    # output_path = tmp_output_dir + '/helitron_length_dist.txt'
    #get_length_dist(paths, labels, my_pal, output_path)

    # ref_dir = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Zea_mays'
    # reference = ref_dir + '/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna'
    # rename_ref = ref_dir + '/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.rename.fna'
    # rename_reference(reference, rename_ref)

    #run_LTR_test()

    tools_dir = os.getcwd() + '/../tools'
    # tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa'
    # repeats = tmp_output_dir + '/candidate_helitron_0.fa'
    # repeats_cons = tmp_output_dir + '/candidate_helitron_0.cons.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + repeats + ' -o ' + repeats_cons + ' -T 0 -M 0'
    # os.system(cd_hit_command)
    #
    # repeats_rename_cons = tmp_output_dir + '/candidate_helitron_0.cons.rename.fa'
    # rename_fasta(repeats_cons, repeats_rename_cons)

    # repeats = tmp_output_dir + '/tir_tsd_0.filter_tandem.fa'
    # repeats_cons = tmp_output_dir + '/tir_tsd_0.cons.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + repeats + ' -o ' + repeats_cons + ' -T 0 -M 0'
    # os.system(cd_hit_command)
    #
    # repeats_rename_cons = tmp_output_dir + '/tir_tsd_0.cons.rename.fa'
    # rename_fasta(repeats_cons, repeats_rename_cons)

    # repeats = tmp_output_dir + '/longest_repeats_0.fa'
    # repeats_cons = tmp_output_dir + '/longest_repeats_0.cons.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + repeats + ' -o ' + repeats_cons + ' -T 0 -M 0'
    # os.system(cd_hit_command)
    #
    # repeats_rename_cons = tmp_output_dir + '/longest_repeats_0.cons.rename.fa'
    # rename_fasta(repeats_cons, repeats_rename_cons)


    #
    # repeats_cons = tmp_dir + '/tir_tsd_0.cons.fa'
    # repeats_rename_cons = tmp_dir + '/tir_tsd_0.cons.rename.fa'
    # rename_fasta(repeats_cons, repeats_rename_cons)

    # test_no_RepeatMasking_time()

    #tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library/RepeatMasking_test/rice_kmer'
    #tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library/recover_test/cb_super'
    #test_connect_RepeatMasking_results(tmp_output_dir)


    # repeats = tmp_output_dir + '/confident_tir_0.fa'
    # repeats_rename = tmp_output_dir + '/confident_tir_0.rename.fa'
    # rename_fasta(repeats, repeats_rename)
    #
    # repeats_rename = tmp_output_dir + '/confident_tir_0.rename.fa'
    # repeats_rename_cons = tmp_output_dir + '/confident_tir_0.cons.rename.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + repeats_rename + ' -o ' + repeats_rename_cons + ' -T 0 -M 0'
    # os.system(cd_hit_command)


    #lost_TIRs(tmp_output_dir)

    # lost_tirs_path = tmp_output_dir + '/test.fa'
    # get_seq_copies(lost_tirs_path, tmp_output_dir)

    # sMITE_path = tmp_output_dir + '/sMITE.copies.fa'
    # Hi_TIR_Ghd2_path = tmp_output_dir + '/Hi_TIR_Ghd2.copies.fa'
    # dist_path = tmp_output_dir + '/MITE_dist.txt'
    # generate_MITE_identity_dist(sMITE_path, Hi_TIR_Ghd2_path, tmp_output_dir, dist_path)
    dist_path = tmp_output_dir + '/MITE_dist.txt'
    my_pal = {"sMITE": "#16499D", "Hi_TIR_Ghd2": "#E71F19"}
    #draw_violin(dist_path, my_pal)




    # repeats = tmp_output_dir + '/tir_tsd_0.filter_tandem.fa'
    #repeats_cons = tmp_output_dir + '/tir_tsd_0.cons.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + repeats + ' -o ' + repeats_cons + ' -T 0 -M 0'
    # os.system(cd_hit_command)

    # repeats_cons_rename = tmp_output_dir + '/tir_tsd_0.rename.cons.fa'
    # rename_fasta(repeats_cons, repeats_cons_rename)

    # flanking_region_distance = 10
    # flanking_len = 50
    # flank_align_dir = tmp_output_dir + '/flank_tir_align_0'
    # mark = 'N_21202-len_1496-ref_NC_029260.1-22739190-22740686-C_29-tsd_TTTTTTCAT-distance_18'
    # for name in os.listdir(flank_align_dir):
    #     if name.endswith('.out'):
    #         with open(flank_align_dir+'/'+name) as f_r:
    #             for line in f_r:
    #                 if line.__contains__(mark):
    #                     print(name)

    TE_type = 'tir'
    tir_tsd_temp_dir = tmp_output_dir + '/' + TE_type + '_blast'
    # mark = 'N_21202-len_1496-ref_NC_029260.1-22739190-22740686-C_29-tsd_TTTTTTCAT-distance_18'
    # for name in os.listdir(tir_tsd_temp_dir):
    #     if name.endswith('.out'):
    #         with open(tir_tsd_temp_dir+'/'+name) as f_r:
    #             for line in f_r:
    #                 if line.__contains__(mark):
    #                     print(name)

    # # output = flank_align_dir + '/N_54503-len_5175-ref_NC_029264.1-11302200-11307375-C_27-tsd_TAA-distance_9.out'
    # output = flank_align_dir + '/1952.out'
    # judge_flank_align(flanking_region_distance, output, flanking_len, flank_align_dir)


    # split_repeats_path = tir_tsd_temp_dir + '/790.fa'
    # ref_db_path = tir_tsd_temp_dir + '/790.fa'
    # blastnResults_path = tir_tsd_temp_dir + '/790.out'
    # ref_db_path = tmp_output_dir + '/GCF_001433935.1_IRGSP-1.0_genomic.fna.cut0.fa'
    # query_coverage = 0.99
    # subject_coverage = 0
    # all_copies = get_copies_v1(blastnResults_path, split_repeats_path, ref_db_path, query_coverage=query_coverage, subject_coverage=subject_coverage)

    #解析itrsearch log文件，提取比对偏移的序列名称
    tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library'
    # for name in os.listdir(tmp_output_dir):
    #     if name.endswith('.log'):
    #         file = tmp_output_dir+'/'+name
    #         print('file:' + file)
    #         get_fake_tirs(file)

    #统计下开头就是mismatch或indel的repbase序列有多少个
    # itrsearch_log = tmp_output_dir + '/curated_lib/repbase/rice/tir.repbase.ref.log'
    # get_fake_tirs(itrsearch_log)

    # tmp_output_dir = '/home/hukang/HiTE/demo/test'
    # ref_name = 'genome.fa'
    # keep_files = ['genome_all.fa.harvest.scn', ref_name + '.rename.fa' + '.finder.combine.scn',
    #               ref_name + '.rename.fa' + '.LTRlib.fa', 'confident_TE.cons.fa',
    #               'confident_TE.cons.fa.classified', 'longest_repeats_*.flanked.fa', 'longest_repeats_*.fa',
    #               'confident_tir_*.fa', 'confident_helitron_*.fa', 'confident_other_*.fa']
    #
    # all_files = os.listdir(tmp_output_dir)
    # for filename in all_files:
    #     for keep_file in keep_files:
    #         if re.search('r'+keep_file, filename) is None:
    #             #os.system('rm -rf ' + tmp_output_dir + '/' + filename)
    #             #print('rm -rf ' + filename)
    #             l = 1
    #         else:
    #             print(keep_file)
    #             print(filename)
    #             print('here')

    # max_ram_usage = 0
    # memlog = '/home/hukang/memlog.txt'
    # with open(memlog, 'r') as f_r:
    #     for line in f_r:
    #         ram_usage = float(line.split(':')[1].split('GB')[0].strip())
    #         #print(ram_usage)
    #         if ram_usage > max_ram_usage:
    #             max_ram_usage = ram_usage
    # print(max_ram_usage)
