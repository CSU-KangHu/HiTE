#-- coding: UTF-8 --
import argparse
import os
import sys

import codecs

import json
import time

# import numpy as np
# from matplotlib import pyplot as plt
# import seaborn as sns
# import pandas as pd
import subprocess

cur_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(cur_dir)
from Util import read_fasta, store_fasta, Logger, read_fasta_v1, rename_fasta, getReverseSequence, allow_mismatch, \
    run_itrsearch, multi_process_itr, filter_large_gap_tirs, multi_process_align_and_get_copies, flanking_copies, \
    store_copies_v1, get_TSD, store_copies, store_LTR_seq_v1, store_LTR_seq, store_LTR_seq_v2, rename_reference, \
    run_LTR_harvest, run_LTR_retriever, determine_repeat_boundary_v2


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
    allow_mismatch_num = 0
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

def identify_new_TIR():
    # 用cd-hit-est，使用-aS 0.8 –aL 0.8 –c 0.8进行聚类，然后分析类中没有curated library出现的转座子为新的TIR。
    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa'
    confident_tir_path = tmp_output_dir + '/confident_tir.rename.cons.fa'
    repbase_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/repbase/rice'
    tir_repbase_path = repbase_dir + '/tir.repbase.ref'
    confident_tir_consensus = tmp_output_dir + '/confident_tir.rename.cons.fa'
    total_tir_path = tmp_output_dir + '/total_tir.fa'
    total_tir_consensus = tmp_output_dir + '/total_tir.cons.fa'
    os.system('cat '+tir_repbase_path+' '+confident_tir_consensus+' > '+total_tir_path)

    confident_tir_names, confident_tir_contigs = read_fasta(confident_tir_consensus)

    tools_dir = os.getcwd() + '/../tools'
    cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.8) + ' -aL ' + str(0.8) + ' -c ' + str(0.8) \
                     + ' -G 0 -g 1 -A 80 -i ' + total_tir_path + ' -o ' + total_tir_consensus + ' -T 0 -M 0'
    #os.system(cd_hit_command)

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
    TRsearch_dir = '/public/home/hpc194701009/repeat_detect_tools/REPET_linux-x64-3.0/bin'
    #run_itrsearch(TRsearch_dir, tir_repbase_path, repbase_dir)
    tir_repbase_out = tir_repbase_path + '.itr'
    repbase_itr_names, repbase_itr_contigs = read_fasta_v1(tir_repbase_out)

    #run_itrsearch(TRsearch_dir, novel_tir_consensus, tmp_output_dir)
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
    #os.system(cd_hit_command)

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
    print('novel TIR with novel TIR terminals:')
    print(new_tir_names)
    print(len(new_tir_names))

    #统计novel TIR的拷贝数量
    threads = 48
    reference = tmp_output_dir + '/GCF_001433935.1_IRGSP-1.0_genomic.fna'
    temp_dir = tmp_output_dir + '/temp'
    blast_program_dir = '/public/home/hpc194701009/repeat_detect_tools/rmblast-2.9.0-p2'
    all_copies = multi_process_align_and_get_copies(novel_tir_consensus, reference, blast_program_dir,
                                                    temp_dir, 'tir', threads)

    # 在copies的两端 flanking 20bp的序列
    flanking_len = 20
    all_copies = flanking_copies(all_copies, novel_tir_consensus, reference, flanking_len, copy_num=-1)

    multicopy_novel_tirs_num = 0
    # 统计每个query_name拷贝的数量
    # query_name -> (ref_name, copy_ref_start, copy_ref_end, copy_len, copy_seq)
    query_copy_num = {}
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        query_copy_num[query_name] = len(copies)
        if len(copies) > 3:
            multicopy_novel_tirs_num += 1
    print(query_copy_num)
    print('novel TIR with copy number > 3:')
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


def draw_dist():
    tmp_output_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa'
    query_copy_num_path = tmp_output_dir + '/novel_tir_copies_num.csv'
    file = open(query_copy_num_path, 'r')
    js = file.read()
    query_copy_num = json.loads(js)
    print(query_copy_num)

    y = list(query_copy_num.values())
    x = pd.Series(y, name="copy number")
    sns.set_theme(style="whitegrid", font='Times New Roman', font_scale=1.4)
    ax = sns.distplot(x, kde=True)
    #plt.show()
    #plt.savefig(tmp_output_dir + "/copy_num.eps", format='eps', dpi=1000)

    plt.savefig(tmp_output_dir + "/copy_num.svg", format='svg', dpi=150, bbox_inches='tight')


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
                                         fixed_extend_base_threshold, max_repeat_len, tmp_output_dir, debug, threads)

            endtime = time.time()
            dtime = endtime - starttime
            log.logger.info("Running time of generating longest_repeats.fa: %.8s s" % (dtime))
        break


if __name__ == '__main__':
    repbase_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/repbase'
    tmp_out_dir = repbase_dir + '/dmel'
    ltr_repbase_path = tmp_out_dir + '/ltr.repbase.ref'
    tir_repbase_path = tmp_out_dir + '/tir.repbase.ref'

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

    #identify_new_TIR()

    #draw_dist()

    #test_EAHelitron()

    #reduce_library_size()

    # ref_dir = '/public/home/hpc194701009/WebTE_Lib/New_cash_crops/Zea_mays'
    # reference = ref_dir + '/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna'
    # rename_ref = ref_dir + '/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.rename.fna'
    # rename_reference(reference, rename_ref)

    #run_LTR_test()

    # tools_dir = os.getcwd() + '/../tools'
    tmp_dir = '/public/home/hpc194701009/KmerRepFinder_test/library/KmerRepFinder_lib/test_2022_0914/oryza_sativa'
    # repeats = tmp_dir + '/candidate_helitron_0.fa'
    # repeats_cons = tmp_dir + '/candidate_helitron_0.cons.fa'
    # cd_hit_command = tools_dir + '/cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + repeats + ' -o ' + repeats_cons + ' -T 0 -M 0'
    # os.system(cd_hit_command)


    # repeats_cons = tmp_dir + '/candidate_helitron_0.cons.fa'
    # repeats_rename_cons = tmp_dir + '/candidate_helitron_0.cons.rename.fa'
    # rename_fasta(repeats_cons, repeats_rename_cons)
    #
    # repeats_cons = tmp_dir + '/tir_tsd_0.cons.fa'
    # repeats_rename_cons = tmp_dir + '/tir_tsd_0.cons.rename.fa'
    # rename_fasta(repeats_cons, repeats_rename_cons)

    #test_no_RepeatMasking_time()
    tmp_output_dir = '/homeb/hukang/KmerRepFinder_test/library/recover_test/cb'
    ref_name = 'GCF_000004555.2_CB4_genomic.rename'
    keep_files_temp = ['longest_repeats_*.flanked.fa', 'longest_repeats_*.fa',
                       'confident_tir_*.fa', 'confident_helitron_*.fa', 'confident_other_*.fa']
    keep_files = ['genome_all.fa.harvest.scn', ref_name + '.rename.fa' + '.finder.combine.scn',
                  ref_name + '.rename.fa' + '.LTRlib.fa', 'confident_TE.cons.fa',
                  'confident_TE.cons.fa.final.classified']

    cut_references = ['0']
    for ref_index, cut_reference in enumerate(cut_references):
        for filename in keep_files_temp:
            keep_files.append(filename.replace('*', str(ref_index)))

    all_files = os.listdir(tmp_output_dir)
    for filename in all_files:
        if filename not in keep_files:
            os.system('rm -rf ' + tmp_output_dir + '/' + filename)




