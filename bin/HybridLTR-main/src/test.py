#-- coding: UTF-8 --
import base64
import heapq
import json
import os
import random
import re
import subprocess
import sys
from collections import defaultdict
import math
from itertools import product

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn import datasets
from sklearn.manifold import TSNE
from concurrent.futures import ProcessPoolExecutor, as_completed


current_folder = os.path.dirname(os.path.abspath(__file__))
# 添加 configs 文件夹的路径到 Python 路径
configs_folder = os.path.join(current_folder, "..")  # 需要根据实际目录结构调整
sys.path.append(configs_folder)

from configs import config
from utils.data_util import read_fasta, store_fasta, read_fasta_v1, replace_non_atcg, get_flanking_copies, \
    get_copies_TSD_info, search_TSD_regular, extract_non_autonomous, run_command, \
    transfer_RMOut2Bed, generate_random_sequences, generate_random_sequence, word_seq, generate_kmer_list, \
    cluster_sequences, generate_random_sequence_v1, generate_expand_matrix, find_files_recursively, \
    split_list_into_groups, expand_matrix_dir, get_matrix_feature_v3, generate_kmer_dic, expand_matrix_dir_v1, \
    generate_sort_matrix
from utils.evaluate_util import generate_TERL_dataset, generate_ClassifyTE_dataset, evaluate_RepeatClassifier, \
    evaluate_TERL, evaluate_DeepTE, transform_DeepTE_to_fasta, add_ClassifyTE_classification, evaluate_ClassifyTE, \
    evaluate_TEsorter, merge_overlapping_intervals, get_metrics_by_label, plot_3D_param, analyze_class_ratio_gff
from Util import judge_both_ends_frame_v1, filter_ltr_by_homo, multi_process_align_v2, get_full_length_copies_RM, \
    read_scn, get_recombination_ltr, merge_terminals, deredundant_for_LTR, cons_from_mafft, generate_full_length_out, \
    is_TIR_frame, get_both_ends_frame, generate_both_ends_frame_from_seq, rename_reference, get_overlap_len, \
    merge_overlap_seq, judge_left_frame_LTR, get_LTR_seq_from_scn, search_ltr_structure, get_low_copy_LTR, \
    get_high_copy_LTR, judge_right_frame_LTR, judge_boundary_v5, judge_boundary_v6, remove_sparse_col_in_align_file, \
    random_downsample, find_tir_in_ltr, get_full_length_copies_batch, multi_process_align_v3, \
    judge_ltr_from_both_ends_frame, generate_both_ends_frame_for_intactLTR, multi_process_align, judge_both_ends_frame, \
    get_full_length_copies, get_full_length_copies_v1, map_fragment, get_intact_ltr_copies, \
    get_copies_v2, get_domain_info_v1, get_domain_info_v2, remove_copies_from_redundant_contig, \
    remove_copies_from_redundant_contig_v1, is_ltr_has_structure, \
    deredundant_for_LTR_v5, get_ltr_from_line, get_all_potential_ltr_lines, \
    multi_process_align_v1, filter_ltr_by_copy_num_sub, Logger, alter_deep_learning_results
from clean_LTR_internal import purge_internal_seq, purge_internal_seq_by_table
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score, classification_report

def connect_LTR(repbase_path):
    # Preprocess. connect LTR and LTR_internal
    # considering reverse complementary sequence
    raw_names, raw_contigs = read_fasta(repbase_path)
    label_names, label_contigs = read_fasta_v1(repbase_path)
    # store repbase name and label
    repbase_labels = {}
    for name in label_names:
        parts = name.split('\t')
        repbase_name = parts[0]
        classification = parts[1]
        species_name = parts[2]
        repbase_labels[repbase_name] = (classification, species_name)

    # 获取所有LTR序列
    LTR_names = set()
    for name in raw_names:
        # find LTR internal
        parts = name.split('-LTR')
        if len(parts) > 1:
            suffix = parts[1]
            prefix = parts[0]
            # find the LTR seq
            ltr_name = prefix + '-LTR' + suffix
            internal_name1 = prefix + '-I' + suffix
            internal_name2 = prefix + '-INT' + suffix
            LTR_names.add(ltr_name)
            LTR_names.add(internal_name1)
            LTR_names.add(internal_name2)

    # 存储分段的LTR与完整LTR的对应关系
    SegLTR2intactLTR = {}
    new_names = []
    new_contigs = {}
    for name in raw_names:
        if name in LTR_names:
            parts = name.split('-LTR')
            if len(parts) > 1:
                # 为LTR终端序列
                suffix = parts[1]
                prefix = parts[0]
                # find the LTR seq
                ltr_name = prefix + '-LTR' + suffix
                internal_name1 = prefix + '-I' + suffix
                internal_name2 = prefix + '-INT' + suffix
                if raw_contigs.__contains__(ltr_name):
                    if raw_contigs.__contains__(internal_name1):
                        internal_name = internal_name1
                        internal_seq = raw_contigs[internal_name1]
                    elif raw_contigs.__contains__(internal_name2):
                        internal_name = internal_name2
                        internal_seq = raw_contigs[internal_name2]
                    else:
                        internal_name = None
                        internal_seq = None
                    if internal_seq is not None:
                        intact_ltr_name = prefix + '-intactLTR' + suffix
                        intact_ltr_seq = raw_contigs[ltr_name] + internal_seq + raw_contigs[ltr_name]
                        new_names.append(intact_ltr_name)
                        new_contigs[intact_ltr_name] = intact_ltr_seq
                        repbase_labels[intact_ltr_name] = repbase_labels[ltr_name]
                        SegLTR2intactLTR[ltr_name] = intact_ltr_name
                        SegLTR2intactLTR[internal_name] = intact_ltr_name
                    else:
                        new_names.append(name)
                        new_contigs[name] = raw_contigs[name]
        else:
            new_names.append(name)
            new_contigs[name] = raw_contigs[name]

    # Step4. store Repbase sequence with classification, species_name, and TSD sequence
    # 去掉processed_TE_path中存在重复的LTR，例如Copia-1_AA-intactLTR1和Copia-1_AA-intactLTR2，取其中具有合法TSD那个。两个都有，则随机去一个；两个都没有，优先取有TSD那个，否则随机取一个。
    # get all classification
    all_classification = set()
    final_repbase_contigs = {}
    duplicate_ltr = set()
    for query_name in new_names:
        label_item = repbase_labels[query_name]
        cur_prefix = query_name.split('-intactLTR')[0]
        if not duplicate_ltr.__contains__(cur_prefix):
            new_name = query_name + '\t' + label_item[0] + '\t' + label_item[1]
            duplicate_ltr.add(cur_prefix)
            final_repbase_contigs[new_name] = new_contigs[query_name]
            all_classification.add(label_item[0])
    store_fasta(final_repbase_contigs, repbase_path)

    # 存储分段的LTR与完整LTR的对应关系
    SegLTR2intactLTRMap = config.work_dir + '/segLTR2intactLTR.map'
    with open(SegLTR2intactLTRMap, 'a+') as f_save:
        for name in SegLTR2intactLTR.keys():
            intact_ltr_name = SegLTR2intactLTR[name]
            f_save.write(name + '\t' + intact_ltr_name + '\n')
    return repbase_path, repbase_labels

def get_all_files(directory):
    all_files = []
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            all_files.append(file_path)
    return all_files


def generate_predict_gff(rice_gff, rice_repbase_gff, label_dict, chromosomes):
    with open(rice_repbase_gff, 'w') as f_save:
        with open(rice_gff, 'r') as f_r:
            for line in f_r:
                if not line.startswith('#'):
                    # line = 'Chr' + line
                    parts = line.split('\t')
                    chr_name = parts[0]
                    if chr_name not in chromosomes:
                        continue
                    seq_name = parts[8].split(' ')[1].replace('\"', '').split(':')[1]
                    label = label_dict[seq_name]
                    new_line = ''
                    for i in range(len(parts)):
                        if i == 2:
                            new_line += label
                        else:
                            new_line += parts[i]
                        if i != len(parts) - 1:
                            new_line += '\t'
                else:
                    new_line = line
                f_save.write(new_line)

def generate_seq_logos(TE_path, tmp_output_dir):
    labels = config.all_wicker_class.keys()
    label_terminals = {}
    logo_len = 30
    label_max_num = 500
    names, contigs = read_fasta_v1(TE_path)
    for name in names:
        label_name = name.split('\t')[1]
        if labels.__contains__(label_name):
            if not label_terminals.__contains__(label_name):
                label_terminals[label_name] = {}
            label_logos = label_terminals[label_name]
            # if len(label_logos) > 2*label_max_num:
            #     continue
            #取序列的开始和结束序列
            seq = contigs[name]
            seq_start = seq[0:logo_len]
            seq_end = seq[-logo_len:]
            if not label_logos.__contains__('5\'-'+label_name):
                label_logos['5\'-'+label_name] = []
            seqs = label_logos['5\'-'+label_name]
            seqs.append(seq_start)

            if not label_logos.__contains__('3\'-'+label_name):
                label_logos['3\'-'+label_name] = []
            seqs = label_logos['3\'-'+label_name]
            seqs.append(seq_end)

    for label_name in label_terminals.keys():
        label_logos = label_terminals[label_name]
        output_path = tmp_output_dir + '/logos_' + label_name + '.txt'
        title_names = label_logos.keys()
        with open(output_path, 'w') as f_save:
            # for logo in title_names:
            #     f_save.write(logo + '\t')
            # f_save.write('\n')
            line_num = 0
            stop = False
            while(not stop):
                col_finish_num = 0
                for name in title_names:
                    seqs = label_logos[name]
                    if line_num >= len(seqs):
                        col_finish_num += 1
                        f_save.write(' \t')
                    else:
                        f_save.write(seqs[line_num]+'\t')
                if col_finish_num >= len(title_names):
                    stop = True
                f_save.write('\n')
                line_num += 1

def identify_new_TIR(confident_tir_path, tir_repbase_path, tmp_output_dir):
    # 用cd-hit-est，使用-aS 0.8 –aL 0.8 –c 0.8进行聚类，然后分析类中没有curated library出现的转座子为新的TIR。
    total_tir_path = tmp_output_dir + '/total_tir.fa'
    total_tir_consensus = tmp_output_dir + '/total_tir.cons.fa'
    os.system('cat '+tir_repbase_path+' '+confident_tir_path+' > '+total_tir_path)

    confident_tir_names, confident_tir_contigs = read_fasta(confident_tir_path)

    cd_hit_command = 'cd-hit-est -aS ' + str(0.8) + ' -aL ' + str(0.8) + ' -c ' + str(0.8) \
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

def analyz_TIR_insert_time(tir_path, work_dir, miu, type, color):
    tool_dir = config.project_dir + '/tools'
    log_path = tir_path + '.log'
    itrsearch_command = 'cd ' + work_dir + ' && ' + tool_dir + '/itrsearch -i 0.7 -l 7 ' + tir_path + ' > ' + log_path
    run_command(itrsearch_command)

    identity_list = []
    with open(log_path, 'r') as f_r:
        for line in f_r:
            if line.__contains__('Identity percentage'):
                identity = float(line.split(':')[1].strip())
                identity_list.append(identity)

    zero_count = 0
    insert_time_list = []
    for identity in identity_list:
        T = estimate_insert_time(identity, miu)
        if T == 0:
            zero_count += 1
        insert_time_list.append(T/1000000)
    print(zero_count)
    output_path = work_dir + '/'+type+'_insert_time.txt'
    output_fig = work_dir + '/'+type+'_insert_time.png'
    with open(output_path, 'w') as f_save:
        f_save.write('insert_time\tmethods\n')
        for insert_time in insert_time_list:
            f_save.write(str(insert_time)+'\t'+'HiTE-TIR'+'\n')
    get_insert_time_dist_boxplot(output_path, output_fig, type, color)

def estimate_insert_time(identity, miu):
    # 估计序列的插入时间
    d = 1 - identity
    K=  -3/4*math.log(1-d*4/3)
    T = K/(2*miu)
    return T

def get_insert_time_dist_boxplot(output_path, output_fig, type, color):
    df = pd.read_csv(output_path, sep='\t', encoding='utf-8')
    print(df)
    # sns.set(font_scale=1.4)
    # sns.set_style("whitegrid", {'axes.facecolor': 'none'})
    # 自定义颜色列表
    # my_colors = ["#293991", "#9DCFEF"]
    my_colors = []
    my_colors.append(color)
    # 设置自定义颜色列表
    sns.set_palette(my_colors)

    # 创建子图布局
    #fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    # 绘制 HiTE-TIR 方法的长度分布
    ax = sns.displot(df[df['methods'] == 'HiTE-TIR']['insert_time'], kde=False, bins=50, rug=False)
    #ax.set_title(type + ' Insertion Time Distribution')
    # ax.set_xlabel('Mya (million years ago)')
    # ax.set_ylabel('Count')
    ax.set(title=type + ' Insertion Time Distribution')
    # ax.set(xlabel='Mya (million years ago)')
    # ax.set(ylabel='Count')
    plt.xlabel('Mya (million years ago)', fontsize=14)
    plt.ylabel('Count', fontsize=14)
    plt.title(type + ' Insertion Time Distribution', fontsize=16)

    # 调整子图之间的间距
    plt.tight_layout()
    # 显示图形
    #plt.show()
    plt.savefig(output_fig, format='png')
    plt.clf()


def connect_segments(segs):
    connected_segs = []
    visited_index = set()
    for i in range(len(segs)):
        if i in visited_index:
            continue
        cur_seg = segs[i]
        visited_index.add(i)
        for j in range(i+1, len(segs)):
            if j in visited_index:
                continue
            next_seg = segs[j]
            overlap_len = get_overlap_len(cur_seg, next_seg)
            if overlap_len > 0:
                merge_seg = merge_overlap_seq(cur_seg, next_seg)
                cur_seg = merge_seg
                visited_index.add(j)
            else:
                break
        connected_segs.append(cur_seg)
    return connected_segs


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
    work_dir = '/home/hukang/HybridLTR/demo/zebrafish'
    genome = work_dir + '/genome.rename.fa'
    std_out = work_dir + '/repbase.edta.out'
    test_out = work_dir + '/HiTE.edta.out'
    # Convert .out file to .bed file
    convert2bed_command = 'perl ' + config.project_dir + '/tools/RMout_to_bed.pl ' + std_out + ' base1'
    os.system(convert2bed_command)
    convert2bed_command = 'perl ' + config.project_dir + '/tools/RMout_to_bed.pl ' + test_out + ' base1'
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
    work_dir = '/home/hukang/HybridLTR/demo/rice_high_FP'
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



# from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
# from cryptography.hazmat.primitives import padding
# from cryptography.hazmat.backends import default_backend
import os

# 固定密钥和初始化向量
KEY = b'This is a key123'  # 必须是16, 24, 或 32 字节长 (128, 192, or 256 bits)
IV = b'This is an IV456'  # 必须是16字节长 (128 bits)


# 编码（加密）
def encode(plaintext, key, iv):
    cipher = Cipher(algorithms.AES(key), modes.CBC(iv), backend=default_backend())
    encryptor = cipher.encryptor()

    padder = padding.PKCS7(algorithms.AES.block_size).padder()
    padded_data = padder.update(plaintext.encode()) + padder.finalize()

    ciphertext = encryptor.update(padded_data) + encryptor.finalize()
    return base64.b64encode(ciphertext).decode('utf-8')


# 解码（解密）
def decode(ciphertext, key, iv):
    cipher = Cipher(algorithms.AES(key), modes.CBC(iv), backend=default_backend())
    decryptor = cipher.decryptor()

    padded_plaintext = decryptor.update(base64.b64decode(ciphertext)) + decryptor.finalize()

    unpadder = padding.PKCS7(algorithms.AES.block_size).unpadder()
    plaintext = unpadder.update(padded_plaintext) + unpadder.finalize()
    return plaintext.decode()


def matrix2align(matrix_file, align_file):
    with open(align_file, 'w') as f_save:
        seq_id = 0
        with open(matrix_file, 'r') as f_r:
            for line in f_r:
                f_save.write('>seq_'+str(seq_id) + '\n' + line)
                seq_id += 1

work_dir = '/home/hukang/test/HiTE/demo'
log = Logger(work_dir + '/HiTE.log', level='debug')

if __name__ == '__main__':
    # # 测试一下 deep learning模块的性能
    # work_dir = '/home/hukang/left_LTR_real_dataset/five_species_high_copy_bak/Arabidopsis_thaliana'
    # dl_output_path = work_dir + '/out/is_LTR_deep.txt'
    # true_labels = []
    # predict_labels = []
    # with open(dl_output_path, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.replace('\n', '').split('\t')
    #         seq_name = parts[0]
    #         predict_label = int(parts[1])
    #         if seq_name.startswith('chr_') or seq_name.startswith('Chr'):
    #             true_label = 0
    #         else:
    #             true_label = 1
    #         true_labels.append(true_label)
    #         predict_labels.append(predict_label)
    #
    # # 计算 Precision
    # precision = precision_score(true_labels, predict_labels, average='macro')
    # print(f"Precision: {precision:.4f}")
    # # 计算 Recall
    # recall = recall_score(true_labels, predict_labels, average='macro')
    # print(f"Recall: {recall:.4f}")
    # # 计算 F1 Score
    # f1 = f1_score(true_labels, predict_labels, average='macro')
    # print(f"F1 Score: {f1:.4f}")
    # # 计算 Accuracy
    # accuracy = accuracy_score(true_labels, predict_labels)
    # print(f"Accuracy: {accuracy:.4f}")
    # # 生成分类报告
    # report = classification_report(true_labels, predict_labels)
    # print("Classification Report:")
    # print(report)
    #
    # # 测试一下同源算法的性能
    # threads = 40
    # flanking_len = 100
    # positive_hc_output_path = os.path.join(work_dir, 'is_LTR_homo.positive.txt')
    # type = 'High copy'
    # high_copy_output_dir = os.path.join(work_dir, 'positive')
    # judge_ltr_from_both_ends_frame(high_copy_output_dir, positive_hc_output_path, threads, type, flanking_len, log)
    #
    # negative_hc_output_path = os.path.join(work_dir, 'is_LTR_homo.negative.txt')
    # type = 'High copy'
    # high_copy_output_dir = os.path.join(work_dir, 'negative')
    # judge_ltr_from_both_ends_frame(high_copy_output_dir, negative_hc_output_path, threads, type, flanking_len, log)
    #
    # hc_output_path = os.path.join(work_dir, 'is_LTR_homo.txt')
    # os.system('cat ' + positive_hc_output_path + ' ' + negative_hc_output_path + ' > ' + hc_output_path)
    #
    # alter_dl_output_path = os.path.join(work_dir, 'is_LTR_deep.alter.txt')
    # alter_deep_learning_results(dl_output_path, hc_output_path, alter_dl_output_path, high_copy_output_dir, log)
    #
    # true_labels = []
    # predict_labels = []
    # with open(alter_dl_output_path, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.replace('\n', '').split('\t')
    #         seq_name = parts[0]
    #         predict_label = int(parts[1])
    #         if seq_name.startswith('chr_') or seq_name.startswith('Chr'):
    #             true_label = 0
    #         else:
    #             true_label = 1
    #         true_labels.append(true_label)
    #         predict_labels.append(predict_label)
    #
    # # 计算 Precision
    # precision = precision_score(true_labels, predict_labels, average='macro')
    # print(f"Precision: {precision:.4f}")
    # # 计算 Recall
    # recall = recall_score(true_labels, predict_labels, average='macro')
    # print(f"Recall: {recall:.4f}")
    # # 计算 F1 Score
    # f1 = f1_score(true_labels, predict_labels, average='macro')
    # print(f"F1 Score: {f1:.4f}")
    # # 计算 Accuracy
    # accuracy = accuracy_score(true_labels, predict_labels)
    # print(f"Accuracy: {accuracy:.4f}")
    # # 生成分类报告
    # report = classification_report(true_labels, predict_labels)
    # print("Classification Report:")
    # print(report)

    # 测试一下在大规模基因组上获取拷贝是否有问题
    threads = 48
    full_length_threshold = 0.95
    tmp_output_dir = '/tmp/judge_LTR_transposons_3d88e304-8a84-4c49-9cf8-bb47526822f2/HybridLTR_output'
    split_ref_dir = tmp_output_dir + '/ref_chr'
    intact_ltr_path = tmp_output_dir + '/intact_ltr.fa'
    temp_dir = tmp_output_dir + '/intact_ltr_filter'
    reference = '/tmp/judge_LTR_transposons_3d88e304-8a84-4c49-9cf8-bb47526822f2/genome.rename.fa'
    ltr_copies = filter_ltr_by_copy_num_sub(intact_ltr_path, threads, temp_dir, split_ref_dir, full_length_threshold,
                                            max_copy_num=10)



    # # 对真实数据集进行拷贝数扩展
    # threads = 40
    # source_dir = '/home/hukang/left_LTR_real_dataset/raw_data/both_ends_frames_remove_gap/positive'
    # target_dir = '/home/hukang/left_LTR_real_dataset/raw_data/both_ends_frames_remove_gap_sort/positive'
    # min_raw_copy_num = 0
    # expand_matrix_dir_v1(source_dir, target_dir, threads, min_raw_copy_num)
    #
    # source_dir = '/home/hukang/left_LTR_real_dataset/raw_data/both_ends_frames_remove_gap/negative'
    # target_dir = '/home/hukang/left_LTR_real_dataset/raw_data/both_ends_frames_remove_gap_sort/negative'
    # min_raw_copy_num = 0
    # expand_matrix_dir_v1(source_dir, target_dir, threads, min_raw_copy_num)

    # source_dir = '/home/hukang/left_LTR_real_dataset/raw_data/both_ends_frames_other_species/non_ltr_tir_negative'
    # target_dir = '/home/hukang/left_LTR_real_dataset/raw_data/both_ends_frames_other_species/non_ltr_tir_negative_sort'
    # min_raw_copy_num = 0
    # expand_matrix_dir_v1(source_dir, target_dir, threads, min_raw_copy_num)

    # source_dir = '/home/hukang/LTR_Benchmarking/LTR_libraries/NeuralLTR/rice_MSUv7/LTR_finder/ltr_left_frames'
    # target_dir = '/home/hukang/LTR_Benchmarking/LTR_libraries/NeuralLTR/rice_MSUv7/LTR_finder/ltr_left_frames_expand'
    # min_raw_copy_num = 0
    # expand_matrix_dir(source_dir, target_dir, threads, min_raw_copy_num)

    # # 过滤掉数据集中拷贝数 <= 5 的数据
    # threshold = 5
    # source_dir = '/home/hukang/left_LTR_real_dataset/raw_data/both_ends_frames_other_species'
    # target_dir = '/home/hukang/left_LTR_real_dataset/raw_data/both_ends_frames_other_species_bak'
    # file_extension = '.matrix'
    # all_positive_matrix_files = find_files_recursively(source_dir, file_extension)
    # for matrix_file in all_positive_matrix_files:
    #     row_num = 0
    #     with open(matrix_file, 'r') as f_r:
    #         for line in f_r:
    #             row_num += 1
    #     if row_num > threshold:
    #         new_matrix_file = matrix_file.replace('both_ends_frames_other_species', 'both_ends_frames_other_species_bak')
    #         target_dir = os.path.dirname(new_matrix_file)
    #         if not os.path.exists(target_dir):
    #             os.makedirs(target_dir)
    #         os.system('cp ' + matrix_file + ' ' + new_matrix_file)

    # tmp_output_dir = '/home/hukang/NeuralLTR/demo/test1'
    # internal_seq_filter = tmp_output_dir + '/confident_ltr.internal.fa.filter_tandem'
    # project_dir = config.project_dir
    # lib_dir = project_dir + '/databases'
    # line_db = lib_dir + '/Tpases020812LINE'
    # threads = 40
    #
    # clean_internal_seq = internal_seq_filter + '.clean'
    # output_table = internal_seq_filter + '.domain'
    # temp_dir = tmp_output_dir + '/domain'
    # # get_domain_info(internal_seq_filter, line_db, output_table, threads, temp_dir)
    # purge_internal_seq_by_table(internal_seq_filter, line_db, clean_internal_seq, output_table)

    # BM_EDTA()
    # BM_HiTE()

    # identity = 98.0 / 100
    # miu = float(str(1.3e-8))
    # time = int(estimate_insert_time(identity, miu))
    # print(time)


    # reference = '/home/hukang/HybridLTR/demo/GCF_001433935.1_IRGSP-1.0_genomic.rename.fna'
    # tmp_output_dir = '/home/hukang/HybridLTR/demo/rice_high_FP'
    # threads = 40
    # split_ref_dir = tmp_output_dir + '/ref_chr'
    # full_length_threshold = 0.95
    # confident_lines = [('chr_0:1524180-1531026', '1524180 1531026 6847 1524180 1524621 442 1530585 1531026 442 98.9 12 chr_0')]
    #
    #
    # internal_ltrs = {}
    # intact_ltrs = {}
    # intact_ltr_path = tmp_output_dir + '/intact_ltr.fa'
    # ref_names, ref_contigs = read_fasta(reference)
    # for name, line in confident_lines:
    #     parts = line.split(' ')
    #     chr_name = parts[11]
    #     left_ltr_start = int(parts[3])
    #     left_ltr_end = int(parts[4])
    #     right_ltr_start = int(parts[6])
    #     right_ltr_end = int(parts[7])
    #
    #     ref_seq = ref_contigs[chr_name]
    #
    #     intact_ltr_seq = ref_seq[left_ltr_start - 1: right_ltr_end]
    #     internal_ltr_seq = ref_seq[left_ltr_end: right_ltr_start - 1]
    #     internal_ltrs[name] = internal_ltr_seq
    #     if len(intact_ltr_seq) > 0:
    #         intact_ltrs[name] = intact_ltr_seq
    # store_fasta(intact_ltrs, intact_ltr_path)
    #
    # temp_dir = tmp_output_dir + '/intact_ltr_filter'
    # ltr_copies = filter_ltr_by_copy_num_sub(intact_ltr_path, threads, temp_dir, split_ref_dir, full_length_threshold,
    #                                         max_copy_num=10)
    # print(ltr_copies)


    # work_dir = '/home/hukang/HybridLTR/demo/rice_high_FP'
    # genome_path = work_dir + '/genome.rename.fa'
    # test_lib = work_dir + '/test.fa'
    # # test_lib = work_dir + '/confident_ltr.fa'
    # test_tmp_blast_dir = work_dir + '/test_blastn'
    # test_lib_out = work_dir + '/test_full_length.out'
    # thread = 40
    #
    # names, contigs = read_fasta(genome_path)
    # chrom_length = {}
    # for i, name in enumerate(names):
    #     chr_len = len(contigs[name])
    #     chrom_length[name] = chr_len
    #
    # coverage_threshold = 0.95
    # category = 'Total'
    # is_full_length = 1
    # multi_process_align_v1(test_lib, genome_path, test_lib_out, test_tmp_blast_dir, thread, chrom_length,
    #                        coverage_threshold, category, is_full_length, is_removed_dir=False)



    # # # 测试内部序列去冗余能否解决玉米的问题
    # tmp_output_dir = '/home/hukang/HybridLTR/demo/maize'
    # reference = tmp_output_dir + '/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.rename.fna'
    # scn_file = tmp_output_dir + '/confident_ltr.scn'
    # dirty_dicts = {}
    #
    # confident_ltr_terminal = tmp_output_dir + '/confident_ltr.terminal.fa'
    # confident_ltr_internal = tmp_output_dir + '/confident_ltr.internal.fa'
    # # get_LTR_seq_from_scn(reference, scn_file, confident_ltr_terminal, confident_ltr_internal, dirty_dicts)
    #
    # type = 'internal'
    # threads = 40
    # internal_coverage_threshold = 0.8
    # deredundant_for_LTR_v5(confident_ltr_internal, tmp_output_dir, threads, type, internal_coverage_threshold)

    # test1_out = tmp_output_dir + '/chr_7_83450278-83457709-int.out'
    # test2_out = tmp_output_dir + '/chr_9_55186310-55193461-int.out'
    # test1_set = set()
    # test2_set = set()
    # with open(test1_out, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.split('\t')
    #         query_name = parts[0]
    #         subject_name = parts[1]
    #         if query_name == subject_name:
    #             continue
    #         test1_set.add(query_name)
    #         test1_set.add(subject_name)
    #
    # with open(test2_out, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.split('\t')
    #         query_name = parts[0]
    #         subject_name = parts[1]
    #         if query_name == subject_name:
    #             continue
    #         test2_set.add(query_name)
    #         test2_set.add(subject_name)
    #
    # intersection = test1_set.intersection(test2_set)
    # print(intersection)

    # work_dir = '/home/hukang/HybridLTR/demo/rice_high_FP'
    # reference = '/home/hukang/HybridLTR/demo/GCF_001433935.1_IRGSP-1.0_genomic.rename.fna'
    # ref_names, ref_contigs = read_fasta(reference)

    # scn_file = work_dir + '/remove_recomb.scn'
    # # scn_file = work_dir + '/test.scn'
    # ltr_candidates, ltr_lines = read_scn(scn_file, log=None, remove_dup=True)
    #
    # confident_lines = []
    # for candidate_index in ltr_lines.keys():
    #     confident_lines.append(ltr_lines[candidate_index])
    #
    # threads = 40
    # temp_path = work_dir + '/all_potential_ltr.json'
    # new_confident_lines = get_all_potential_ltr_lines(confident_lines, reference, threads, temp_path, log=None)
    # # print(new_confident_lines)


    # candidate_index = 0
    # line = '1258001 1264025 6025 1258001 1258445 445 1263582 1264025 444 97.0 NA chr_0'
    # #line = '64181422 64192039 10618 64181422 64181595 174 64191866 64192039 174 98.0 NA Chr426 '
    # parts = line.split(' ')
    # chr_name = parts[11]
    # ref_seq = ref_contigs[chr_name]
    # extend_len = 50
    # left_ltr = ref_seq[1258001 - extend_len: 1258445 + extend_len]
    # right_ltr = ref_seq[1263582 - extend_len: 1264025 + extend_len]
    # job_list = [(candidate_index, left_ltr, right_ltr, 445, 444, 1258001, 1258445, 1263582, 1264025)]
    # print(results)


    # # 过滤来自冗余contig的拷贝，即取拷贝的左侧100bp+右侧100bp组成的序列仍然能够很好比对
    # work_dir = '/home/hukang/HybridLTR-main/demo/zebrafish'
    # reference = work_dir + '/genome.rename.fa'
    # flanking_len = 100
    # ref_names, ref_contigs = read_fasta(reference)
    #
    # intact_ltr_copies = {}
    # # copy1 = ('Chr307', 880067, 891620)
    # # copy2 = ('Chr361', 259235, 270788)
    # # copies = [copy1, copy2]
    # # intact_ltr_copies['chr_12-880068-880443'] = copies
    #
    # # copy1 = ('Chr349', 199084, 220908)
    # # copy2 = ('Chr431', 66061484, 66083654)
    # # copies = [copy1, copy2]
    # # intact_ltr_copies['chr_1246-199085-199697'] = copies
    #
    # copy1 = ('Chr358', 152302, 173068)
    # copy2 = ('Chr307', 8124643, 8145856)
    # copies = [copy1, copy2]
    # intact_ltr_copies['chr_1395-152303-152654'] = copies
    #
    # temp_dir = work_dir + '/intact_ltr_deredundant'
    # threads= 40
    # intact_ltr_copies = remove_copies_from_redundant_contig_v1(intact_ltr_copies, reference, temp_dir, threads)
    # print(intact_ltr_copies)


    # copy1_flank_seq = ref_contigs[copy1[0]][copy1[1]- flanking_len: copy1[2] + flanking_len]
    # copy2_flank_seq = ref_contigs[copy2[0]][copy2[1] - flanking_len: copy2[2] + flanking_len]
    # print(copy1_flank_seq)
    # print(copy2_flank_seq)


    # work_dir = '/home/hukang/HybridLTR-main/demo/test_ath30'
    # matrix_file = work_dir + '/ltr_both_frames/chr_1-4149396-4150946.matrix'
    # keep_matrix_file = work_dir + '/ltr_both_frames_sort/chr_1-4149396-4150946.matrix'
    # temp_seq = work_dir + '/chr_1-4149396-4150946.fa'
    # temp_cons = work_dir + '/chr_1-4149396-4150946.cons'
    # is_keep = generate_sort_matrix(matrix_file, keep_matrix_file, temp_seq, temp_cons)
    # print(is_keep)

    # ltr_copies_path = '/home/hukang/HybridLTR-main/demo/rice/ltr_copies.json'
    # with open(ltr_copies_path, 'r', encoding='utf-8') as f:
    #     intact_ltr_copies = json.load(f)
    #
    # for ltr_name in intact_ltr_copies.keys():
    #     copies = intact_ltr_copies[ltr_name]
    #     copies.sort(key=lambda x: (x[2]-x[1]))
    #     filtered_copies = []
    #     for i in range(len(copies)):
    #         chr_name_i, start_i, end_i = copies[i]
    #         length_i = end_i - start_i + 1
    #         is_redundant = False
    #         for j in range(i + 1, len(copies)):
    #             chr_name_j, start_j, end_j = copies[j]
    #             length_j = end_j - start_j + 1
    #             if chr_name_i == chr_name_j:
    #                 overlap_start = max(start_i, start_j)
    #                 overlap_end = min(end_i, end_j)
    #                 overlap_length = overlap_end - overlap_start + 1
    #                 # 判断是否冗余（即 95% 的拷贝 i 被包含在拷贝 j 中）
    #                 if overlap_length >= 0.95 * length_i:
    #                     is_redundant = True
    #                     break
    #         if not is_redundant:
    #             filtered_copies.append((chr_name_i, start_i, end_i))
    #     if len(filtered_copies) != len(copies):
    #         print(ltr_name)
    #     intact_ltr_copies[ltr_name] = filtered_copies

    # work_dir = '/home/hukang/HybridLTR/demo/test_ath3'
    # single_copy_internals_file = work_dir + '/chr_2-13561769-13562356.fa'
    # protein_db = '/home/hukang/HybridLTR/databases/RepeatPeps.lib'
    # threads = 40
    # temp_dir = work_dir + '/domain'
    # tool_dir = '/home/hukang/HybridLTR/tools'
    # query_protein_types = get_domain_info_v2(single_copy_internals_file, protein_db, threads, temp_dir, tool_dir)
    # print(query_protein_types)
    #
    # name = 'chr_2-13561769-13562356'
    # # 如果单拷贝LTR的内部序列有蛋白质、且蛋白质类型单一就认为是真实的LTR
    # # if name in is_single_ltr_has_intact_protein and is_single_ltr_has_intact_protein[name]:
    # if name in query_protein_types:
    #     protein_types = query_protein_types[name]
    #     first_protein = next(iter(protein_types))
    #     if len(protein_types) == 1 and 'LTR' in first_protein:
    #         cur_is_ltr = 1
    #     else:
    #         cur_is_ltr = 0
    # else:
    #     cur_is_ltr = 0
    # print(cur_is_ltr)

    # work_dir = '/home/hukang/HybridLTR-main/demo/test_ath27'
    # single_copy_internals_file = work_dir + '/single_copy_internal.fa'
    # protein_db = '/home/hukang/HybridLTR-main/databases/RepeatPeps.lib'
    # threads = 40
    # temp_dir = work_dir + '/domain'
    # query_protein_types = get_domain_info_v1(single_copy_internals_file, protein_db, threads, temp_dir)
    # print(query_protein_types)

    # # 将scn_line写成fasta格式文件
    # # 199083 220916 21834 199085 199697 613 220297 220908 612 89.0 NA chr_1246
    # # 199083 220916 21834 199085 199696 612 220297 220909 613 89.0 NA chr_1246
    # lines = ['199083 220916 21834 199085 199697 613 220297 220908 612 89.0 NA chr_1246', '199083 220916 21834 199085 199696 612 220297 220909 613 89.0 NA chr_1246']
    # genome = '/home/hukang/LTR_Benchmarking/LTR_libraries/LTR_detector/zebrafish/GCF_000002035.6_GRCz11_genomic.rename.fna'
    # work_dir = '/home/hukang/HybridLTR-main/demo/test_zebrafish'
    # test_path = work_dir + '/test.fa'
    # ref_names, ref_contigs = read_fasta(genome)
    # contigs = {}
    # for line in lines:
    #     parts = line.split(' ')
    #     chr_name = parts[11]
    #     ref_seq = ref_contigs[chr_name]
    #
    #     lLTR_start = int(parts[3])
    #     lLTR_end = int(parts[4])
    #     rLTR_start = int(parts[6])
    #     rLTR_end = int(parts[7])
    #     intact_LTR = ref_seq[lLTR_start - 1: rLTR_end]
    #     internal_LTR = ref_seq[lLTR_end: rLTR_start]
    #     seq_name = str(chr_name) + '-' + str(lLTR_start) + '-' + str(rLTR_end)
    #     contigs[seq_name] = intact_LTR
    # store_fasta(contigs, test_path)
    #
    #
    # work_dir = '/home/hukang/HybridLTR-main/demo/zebrafish'
    # test_path = work_dir + '/chr_1246_199083-220916-int.fa'
    # split_ref_dir = work_dir + '/ref_chr'
    # threads = 40
    # max_copy_num = 100
    # full_length_threshold = 0.95
    # ltr_copies = get_full_length_copies_v1(test_path, split_ref_dir, max_copy_num, full_length_threshold, debug=1)
    # print(ltr_copies)
    # print(len(ltr_copies))





    # # 将 RepeatPeps.lib 拆成 LTRPeps.lib 和 OtherPeps.lib
    # work_dir = '/home/hukang/HybridLTR-main/databases'
    # RepeatPeps_lib = work_dir + '/RepeatPeps.lib'
    # LTR_lib = work_dir + '/LTRPeps.lib'
    # Other_lib = work_dir + '/OtherPeps.lib'
    # names, contigs = read_fasta(RepeatPeps_lib)
    # LTR_contigs = {}
    # Other_contigs = {}
    # TE_types = set()
    # for name in names:
    #     parts = name.split('#')
    #     raw_name = parts[0]
    #     TE_type = parts[1]
    #     if 'LTR' in TE_type:
    #         LTR_contigs[name] = contigs[name]
    #     else:
    #         Other_contigs[name] = contigs[name]
    #     TE_types.add(TE_type)
    # # print(TE_types)
    # store_fasta(LTR_contigs, LTR_lib)
    # store_fasta(Other_contigs, Other_lib)





    # # # 根据LTR终端的比对信息确定是否是LTR
    # #genome = '/home/hukang/HybridLTR-main/demo/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.rename.fna'
    # genome = '/home/hukang/HybridLTR-main/demo/GCF_000001735.4_TAIR10.1_genomic.rename.fna'
    # ref_names, ref_contigs = read_fasta(genome)
    # scn_file = '/home/hukang/HybridLTR-main/demo/test_ath6/test.scn'
    # filter_scn = '/home/hukang/HybridLTR-main/demo/test_ath6/filter_test.scn'
    #
    # # scn_file = '/home/hukang/HybridLTR-main/demo/test_dmel6/remove_recomb.scn'
    # # filter_scn = '/home/hukang/HybridLTR-main/demo/test_dmel6/filter_terminal_align.scn'
    # threads = 40
    # filter_ltr_by_flank_seq_v1(scn_file, filter_scn, genome, threads, log=None)

    # work_dir = '/home/hukang/HybridLTR-main/demo/rice'
    # split_ref_dir = work_dir + '/ref_chr'
    # threads = 40
    # confident_ltr_internal = work_dir + '/test.fa'
    # temp_dir = work_dir + '/internal_blast_temp'
    # max_copy_num = 100
    # full_length_threshold = 0.95
    # ltr_copies = get_full_length_copies_v1(confident_ltr_internal, split_ref_dir, max_copy_num, full_length_threshold, debug=1)
    # print(ltr_copies)
    # print(len(ltr_copies))
    #
    # scn_file = work_dir + '/filter_terminal_align.scn'
    # ltr_candidates, ltr_lines = read_scn(scn_file, log=None)
    # reference = '/home/hukang/HybridLTR-main/demo/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.rename.fna'
    # intact_ltr_copies = get_intact_ltr_copies(ltr_copies, ltr_lines, full_length_threshold, reference)
    # print(intact_ltr_copies)









    # # 去掉斑马鱼金标准中的DIRS转座子
    # work_dir = '/home/hukang/HybridLTR/library'
    # repbase_path = work_dir + '/zebrafish.ltr.ref'
    # new_repbase_path = work_dir + '/zebrafish.ltr.ref.bak'
    # rep_names, rep_contigs = read_fasta(repbase_path)
    # new_rep_contigs = {}
    # for name in rep_names:
    #     if 'DIRS' in name:
    #         continue
    #     new_rep_contigs[name] = rep_contigs[name]
    # store_fasta(new_rep_contigs, new_repbase_path)


    # tmp_output_dir = '/home/hukang/LTR_Benchmarking/LTR_libraries/Ours/maize/test'
    # confident_ltr_tmp = tmp_output_dir + '/confident_ltr.tmp.fa'
    # input1 = confident_ltr_tmp
    # confident_ltr_internal = tmp_output_dir + '/test.fa'
    # input2 = confident_ltr_internal
    # threads = 40
    # iter_num = 1
    # blastnResults_path = tmp_output_dir + '/rm_nested.self.out'
    # confident_TE_blast_dir = tmp_output_dir + '/rm_nested_blast'
    # output = confident_ltr_internal + '.clean_nested'
    # # multi_process_align(input1, input2, blastnResults_path, confident_TE_blast_dir, threads)
    # clean_output = output
    # remove_nest(blastnResults_path, input1, input2, clean_output, iter_num, coverage=0.95)

    # work_dir = '/home/hukang/HybridLTR-main/demo/zebrafish'
    # matrix_file = work_dir + '/chr_13-9928724-9929264.matrix'
    # align_file = work_dir + '/chr_13-9928724-9929264.matrix.maf.fa'
    # matrix2align(matrix_file, align_file)

    # clean_align_file = work_dir + '/chr_4:115052633-115053548.blast.bed.fa.rdmSubset.fa.clean.maf.fa'
    # remove_sparse_col_in_align_file(align_file, clean_align_file)

    # ltr_name = 'test'
    # ltr_name, is_TE, info, model_seq = judge_boundary_v6(matrix_file, ltr_name)
    # print(ltr_name, is_TE, info, model_seq)

    # confident_ltr = work_dir + '/test.fa'
    # reference = '/home/hukang/LTR_Benchmarking/LTR_libraries/LTR_detector/zebrafish/GCF_000002035.6_GRCz11_genomic.rename.fna'
    # ref_names, ref_contigs = read_fasta(genome)
    # chr_name = 'chr_24'
    # start = 24664670
    # end = 24668730
    # flank_len = 100
    # seq = ref_contigs[chr_name][start-1-flank_len: end+flank_len]
    #print(seq)

    # work_dir = '/home/hukang/LTR_Benchmarking/LTR_libraries/Ours/maize/test'
    # reference = '/home/hukang/LTR_Benchmarking/LTR_libraries/LTR_detector/maize/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.rename.fna'
    # tmp_output_dir = work_dir
    # test_home = os.getcwd()
    # split_genome_command = 'cd ' + test_home + ' && python3 ' + test_home + '/split_genome_chunks.py -g ' \
    #                        + reference + ' --tmp_output_dir ' + tmp_output_dir
    # os.system(split_genome_command)
    #
    # split_ref_dir = work_dir + '/ref_chr'
    # threads = 40
    # confident_ltr_internal = work_dir + '/test.fa'
    # temp_dir = work_dir + '/internal_blast_temp'
    # all_copies = get_full_length_copies(confident_ltr_internal, split_ref_dir, debug=1)
    # print(all_copies)
    # print(len(all_copies))

    # # 对于全长拷贝数 >=2 的LTR序列，检测是否由于输入基因组包含多个冗余 contigs 造成
    # # 即，对全长LTR生成窗口，然后调用规则方法判断窗口两侧是否存在同源性
    # tmp_output_dir = work_dir
    # split_ref_dir = work_dir + '/../ref_chr'
    # threads = 40
    # flanking_len = 100
    # temp_dir = tmp_output_dir + '/intact_ltr_filter'
    # output_dir = tmp_output_dir + '/intact_ltr_both_frames'
    # full_length_output_dir = tmp_output_dir + '/intact_ltr_full_length_frames'
    # generate_both_ends_frame_for_intactLTR(confident_ltr, reference, flanking_len, threads, temp_dir, output_dir,
    #                                        full_length_output_dir,
    #                                        split_ref_dir, log=None)
    #
    # type = 'High copy'
    # lc_output_path = tmp_output_dir + '/intact_LTR_homo.txt'
    # judge_ltr_from_both_ends_frame(output_dir, lc_output_path, threads, type, log=None)

    # blastnResults_path = work_dir + '/test.out'
    # tmp_blast_dir = work_dir + '/internal_blast_temp'
    # coverage_threshold = 0.8
    # full_length_annotation = multi_process_align_v3(confident_ltr_internal, genome, blastnResults_path, tmp_blast_dir, threads,
    #                        coverage_threshold, is_removed_dir=True)
    # print(full_length_annotation)

    # blastnResults_path = work_dir + '/test.out'
    # candidate_ltr_path = work_dir + '/test.fa'
    # confident_sine_path = work_dir + '/confident_tir.fa'
    # remain_candidate_tirs = find_tir_in_ltr(blastnResults_path, candidate_ltr_path, confident_sine_path)
    # print(remain_candidate_tirs)

    # query_name = 'chr_4_1150'
    # cur_seq = 'CCGAGCCCTGGGGTCGGGCGAAGCGGAGTTTCGTCGTCTTCCGGGTGCTAGCCCGAGTCCGAGCCCTGGGGTCGGGCGGAGCGGAGTTCGCCGTCTTCCGGGTCTTAGCCCGAGTCCGAGCCCTGGGGTCGGGCGGAGCGGAGTTCGCCGTCTTCCGGGTCTTAGCCCGAGTCCGAGCCCTGGGGTCGGGCGGAGCGGAGTTCGCCGTCTTCCGGGTCTTAGCCCGAGTCCGAGCCCTGGGGTCGGGCGGAGCGGAGTTCGCCGTCTTCCGGGTCTTAGCCCGAGTCCGAGCCCTGGGGTCGGGCGGAGCGGAGTTCGCCGTCTTCCGGGTCTTAGCCCGAGTCCGAGCCCTGGGGTCGGGCGGAGCGGAGTTCGCCGTCTTCCGGGTCTTAGCCCGAGTCCGAGCCCTGGGGTCGGGCGGAGCGGAGTTCGCCGTGGCGCCTTTGGCAAGGCCTGACTGCCTGTCAGACTCACTCTGTCGAGCGGCACTGCAGTCGGAGTGGCGCAGGCGGCGCTGTCCTTCTGTCAGACTGGCCAGTGGAGCAGTGGAGTGACGGCGGTCACTTCGGCTCTGCCGGGGGCGCGTGTCAGGATAGAGGTGTCAGGCCATCTTTGCGTTAAATGCCCCTACAATTTGGTCAGTCGGTGCGGCGATTTAGTCAAGGTTGCTTCTGAGTGAAGCCAAGGCCTCGGGCAAGCCGGTGATGTGTCCGCCATAAAAAGGGGGCCTCGGGCGAGACGGAAGTCTCTCGAGGTCGGCTGCCTTTGGCCGAGGCTAGGCTCGGGTGAAGCGTGATCGAGTCACTCGTGTGGACCGATCCCTGACTTAATCGTACCCATCAGGCCTTTGCAGCTTTATGCTGATGGGGGTTACCAGCTGAGAATTAGGCGTCTTGAGGGTACCCCTAATTATGGTCCCCGAC'
    # output_dir = work_dir + '/output'
    # full_length_output_dir = work_dir + '/full_length'
    # debug = 1
    # both_end_frame_path, full_length_align_file = get_both_ends_frame(query_name, cur_seq, align_file, output_dir,
    #                                                                   full_length_output_dir, debug)


    # # 删除正样本中多余的文件
    # work_dir = '/home/hukang/left_LTR_real_dataset/raw_data/both_ends_frames_remove_gap'
    # positive_dir = work_dir + '/positive'
    # for name in os.listdir(positive_dir):
    #     # 删除full_length目录
    #     species_dir = positive_dir + '/' + name
    #     full_length_dir = species_dir + '/full_length'
    #     os.system('rm -rf ' + full_length_dir)
    #
    #     # 将positive目录下的文件都移动到 species目录下
    #     cur_positive_dir = species_dir + '/positive'
    #     os.system('mv ' + cur_positive_dir + '/* ' + species_dir)
    #     os.system('rm -rf ' + cur_positive_dir)


    # # 过滤正负样本中的拷贝数 > 5
    # work_dir = '/home/hukang/left_LTR_real_dataset/raw_data/both_ends_frames_remove_gap'
    # threshold = 5
    # file_extension = '.matrix'
    # all_positive_matrix_files = find_files_recursively(work_dir, file_extension)
    # for matrix_file in all_positive_matrix_files:
    #     row_num = 0
    #     with open(matrix_file, 'r') as f_r:
    #         for line in f_r:
    #             row_num += 1
    #     if row_num <= threshold:
    #         os.remove(matrix_file)

    # # 对负样本进行下采样
    # negative_dir = '/home/hukang/left_LTR_real_dataset/raw_data/both_ends_frames_remove_gap/negative'
    # file_extension = '.matrix'
    # all_matrix_files = find_files_recursively(negative_dir, file_extension)
    # keep_files = random_downsample(all_matrix_files, 9021)
    # # 遍历原始列表，删除未被选中的文件
    # for file_path in all_matrix_files:
    #     if file_path not in keep_files:
    #         if os.path.exists(file_path):
    #             os.remove(file_path)
    #             print(f"Deleted: {file_path}")


    # sequence = 'CAGTGGCACATTAGGTAGTGCTATCGCCTCACAGCAAAAAGGTCGCTGGTTAGAGCCTCGGCTGGATCAGTTGGCGTTTCTGTGTGAAGTTTGCATGTTCTTCCTGCGTTCGTGTGATTCTCCTCCGGGTGCTTCGGTTTCCCCCACAGTCCAAAGACATGCGGTACAGGTGAATTGAGTAGGCTAAATTGTCTGTAGTGTATGAGTGTATGTGTGAATGAGCTTGTGTGGATGTTTCCCAGAGATGGGTTGCGGCGGAAAGGGCATAAGCTGCATAAAAAAAAACGTGCTGGATAAGTTGGTGGTTTATTTCACTGTGGCGACACCGGATTAAAAAAGGGACTAAGCTGAAAATAAAATGAATGAATGAATGAA'
    # has_sine_tail = identify_SINE_tail(sequence)
    # print(has_sine_tail)

    # # 遍历所有的正样本LTR序列，看看有多少 LTR 满足sine tail
    # repbase_path = '/home/hukang/NeuralTE_dataset/raw_Repbase/all_repbase.ref'
    # repbase_names, repbase_contigs = read_fasta(repbase_path)
    # ltr_num = 0
    # sine_num = 0
    # for name in repbase_names:
    #     if '-LTR' in name or '_LTR' in name:
    #         ltr_num += 1
    #         seq = repbase_contigs[name]
    #         has_sine_tail = identify_SINE_tail(seq)
    #         if has_sine_tail and not seq.startswith('TG') and not seq.endswith('CA'):
    #             sine_num += 1
    #             print(name)
    # print(ltr_num, sine_num)

    # repbase_names, repbase_contigs = read_fasta_v1(repbase_path)
    # sine_num = 0
    # sine_tail_num = 0
    # for name in repbase_names:
    #     parts = name.split('\t')
    #     if len(parts) != 3:
    #         continue
    #     raw_name = parts[0]
    #     label = parts[1]
    #     if 'SINE' in label:
    #         sine_num += 1
    #         seq = repbase_contigs[name]
    #         has_sine_tail = identify_SINE_tail(seq)
    #         if has_sine_tail:
    #             sine_tail_num += 1
    #         else:
    #             print(name)
    # print(sine_num, sine_tail_num)

    # genome = '/home/hukang/LTR_Benchmarking/LTR_libraries/LTR_detector/maize/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.rename.fna'
    # work_dir = '/home/hukang/LTR_Benchmarking/LTR_libraries/Ours/maize/test'
    # scn_path = work_dir + '/test.scn'
    # ltr_terminal = work_dir + '/terminal.fa'
    # ltr_internal = work_dir + '/internal.fa'
    # tool_dir = '/home/hukang/HybridLTR/tools'
    # threads = 40
    # tmp_output_dir = work_dir
    # get_LTR_seq_from_scn(genome, scn_path, ltr_terminal, ltr_internal, tmp_output_dir, tool_dir, threads)

    # line = '6399305 6419978 20674 6399305 6400644 1340 6418639 6419978 1340 98.0 NA Chr1'
    # genome = '/home/hukang/HybridLTR-main/demo/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.rename.fna'
    # # line = '44156089 44175327 19239 44156089 44156442 354 44174966 44175327 362 85.0 NA chr_16'
    # # line = '12521271 12537542 16272 12521281 12521656 376 12537167 12537533 367 86.0 NA chr_8'
    # # genome = '/home/hukang/LTR_Benchmarking/LTR_libraries/LTR_detector/zebrafish/GCF_000002035.6_GRCz11_genomic.rename.fna'
    # ref_names, ref_contigs = read_fasta(genome)
    # left_size = 20
    # internal_size = 3
    #
    # parts = line.split(' ')
    # chr_name = parts[11]
    # ref_seq = ref_contigs[chr_name]
    #
    # LTR_start = int(parts[0])
    # LTR_end = int(parts[1])
    # lLTR_start = int(parts[3])
    # lLTR_end = int(parts[4])
    # rLTR_start = int(parts[6])
    # rLTR_end = int(parts[7])
    # intact_LTR = ref_seq[lLTR_start-1: rLTR_end]
    # print(intact_LTR)
    # print(len(intact_LTR))

    # left_LTR_seq = ref_seq[lLTR_start: lLTR_end]
    # print(left_LTR_seq)
    # print(len(left_LTR_seq))
    #
    # # 取左/右侧 8bp + 3bp
    # # 计算左LTR的切片索引，并确保它们在范围内
    # left_start = max(LTR_start - 1 - left_size, 0)
    # left_end = min(LTR_start + internal_size, len(ref_seq))
    # left_seq = ref_seq[left_start: left_end]
    #
    # # 计算右LTR的切片索引，并确保它们在范围内
    # right_start = max(LTR_end - 1 - internal_size, 0)
    # right_end = min(LTR_end + left_size, len(ref_seq))
    # right_seq = ref_seq[right_start: right_end]
    #
    # ltr_name = 'chr_2_13128730-13141096'
    # ltr_name, has_structure = search_ltr_structure(ltr_name, left_seq, right_seq, left_LTR_seq)
    # print(ltr_name, has_structure)


    # tmp_output_dir = '/home/hukang/NeuralLTR/demo/test4'
    # reference = '/home/hukang/NeuralLTR/demo/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.rename.fna'
    # confident_ltr_terminal = tmp_output_dir + '/confident_ltr.terminal.fa'
    # confident_ltr_internal = tmp_output_dir + '/confident_ltr.internal.fa'
    # scn_file = tmp_output_dir + '/confident_ltr.scn'
    # get_LTR_seq_from_scn(reference, scn_file, confident_ltr_terminal, confident_ltr_internal)

    # work_dir = '/home/hukang/LTR_Benchmarking/LTR_libraries/Ours/zebrafish'
    # matrix_file = work_dir + '/high_copy_frames/chr_22:40137057-40138267.matrix'
    # is_left_ltr, new_boundary_start = judge_left_frame_LTR(matrix_file)
    # debug = 1
    # if debug:
    #     print('left', matrix_file, is_left_ltr, new_boundary_start)
    # is_right_ltr, new_boundary_end = judge_right_frame_LTR(matrix_file)
    # if debug:
    #     print('right', matrix_file, is_right_ltr, new_boundary_end)

    # cur_matrix_file = '/home/hukang/HybridLTR-main/demo/test_zebrafish/chr_1395_152303-152654.matrix'
    # ltr_name = ''
    # TE_type = 'tir'
    # flanking_len = 500
    # ltr_name, is_TE, info, final_cons_seq = judge_boundary_v5(cur_matrix_file, ltr_name, TE_type, flanking_len)
    # print(ltr_name, is_TE, info, final_cons_seq)

    # candidate_tir_path = '/home/hukang/HybridLTR-main/demo/test_ath19/candidate_tir.fa'
    # threads = 40
    # print(all_confident_tirs)
    # print(len(all_confident_tirs))

    # output_dir = '/home/hukang/LTR_Benchmarking/LTR_libraries/Ours/zebrafish/test'
    # expand_output_dir = '/home/hukang/LTR_Benchmarking/LTR_libraries/Ours/zebrafish/test_expand'
    # threads = 40
    # min_raw_copy_num = 0
    # expand_matrix_dir_v1(output_dir, expand_output_dir, threads, min_raw_copy_num)

    # 提取Repbase中的SINE元素


    # # 提取Repbase中的Helitron元素，并生成加密library
    # work_dir = '/home/hukang/NeuralTE_dataset/raw_Repbase'
    # repbase_path = work_dir + '/all_repbase.ref'
    # helitron_lib = work_dir + '/all_helitron.ref'
    # repbase_names, repbase_contigs = read_fasta_v1(repbase_path)
    # helitron_contigs = {}
    # for name in repbase_names:
    #     parts = name.split('\t')
    #     if len(parts) != 3:
    #         continue
    #     label = parts[1]
    #     if label == 'Helitron':
    #         seq = repbase_contigs[name]
    #         # encoded_name = encode(name, KEY, IV)
    #         # encoded_text = encode(seq, KEY, IV)
    #         # helitron_contigs[encoded_name] = encoded_text
    #         helitron_contigs[name] = seq
    # store_fasta(helitron_contigs, helitron_lib)


    # original_text = "This is a secret message."
    #
    # # 编码
    # encoded_text = encode(original_text, KEY, IV)
    # print(f"Encoded text: {encoded_text}")
    #
    # # 解码
    # decoded_text = decode(encoded_text, KEY, IV)
    # print(f"Decoded text: {decoded_text}")


    # work_dir = '/home/hukang/NeuralLTR/demo/test3'
    # high_copy_output_dir = work_dir + '/high_copy_frames'
    # expand_high_copy_output_dir = high_copy_output_dir + '_expand'
    # threads = 4
    # min_raw_copy_num = 0
    # expand_matrix_dir_v1(high_copy_output_dir, expand_high_copy_output_dir, threads, min_raw_copy_num)


    # work_dir = '/home/hukang/NeuralLTR/demo/test1'
    # query_path = work_dir + '/confident_ltr.internal.fa'
    # blast_other = work_dir + '/test.out'
    # purge_internal_seq(query_path, blast_other)

    # genome = '/home/hukang/NeuralLTR/demo/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.rename.fna'
    # ref_names, ref_contigs = read_fasta(genome)
    # seq_name = 'Chr450_26134790-26154779-int'
    # seq = ref_contigs['Chr450'][26135247: 26154322]
    # print(seq)
    # # Chr450_26134790-26154779-int


    # # 将玉米的参考基因组变成Chr+number格式
    # tmp_output_dir = '/home/hukang/LTR_Benchmarking/LTR_libraries/LTR_finder/Zea_mays'
    # reference = tmp_output_dir + '/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.filtered'
    # # rename reference
    # ref_rename_path = tmp_output_dir + '/genome.rename.fa'
    # chr_name_map = tmp_output_dir + '/chr_name.map'
    # # rename_reference(reference, ref_rename_path, chr_name_map)
    # chr_names_dicts = {}
    # with open(chr_name_map, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         parts = line.split('\t')
    #         chr_names_dicts[parts[1]] = parts[0]
    # reference = ref_rename_path
    # scn_file = tmp_output_dir + '/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.filtered.finder.combine.scn'
    # rename_scn_file = tmp_output_dir + '/LTR_finder.rename.scn'
    # lines = []
    # with open(scn_file, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         if line.startswith('#'):
    #             lines.append(line)
    #         else:
    #             parts = line.split(' ')
    #             parts[-1] = chr_names_dicts[parts[-1]]
    #             newline = ' '.join(parts)
    #             lines.append(newline)
    #
    # with open(rename_scn_file, 'w') as f_save:
    #     for line in lines:
    #         f_save.write(line + '\n')

    # confident_ltr_terminal = '/home/hukang/LTR_Benchmarking/LTR_libraries/NeuralLTR/Zea_mays/LTR_finder/confident_ltr.terminal.fa'
    # confident_ltr_internal_cons = '/home/hukang/LTR_Benchmarking/LTR_libraries/NeuralLTR/Zea_mays/LTR_finder/confident_ltr.terminal.fa.cdhit.cons.onlyshort'
    # # tmp_output_dir = '/home/hukang/LTR_Benchmarking/LTR_libraries/NeuralLTR/Zea_mays/LTR_finder'
    # # threads = 40
    # # deredundant_for_LTR(confident_ltr_terminal, tmp_output_dir, threads)
    #
    # cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 80 -i ' + confident_ltr_terminal + ' -o ' + confident_ltr_internal_cons + ' -T 0 -M 0'
    # os.system(cd_hit_command + ' > /dev/null 2>&1')

    # cur_align_file = '/home/hukang/HybridLTR-main/demo/test_dmel10/raw_ltr_copies_cluster_internal/Ninja_13/0.fa.maf.fa'
    # cons_seq = cons_from_mafft(cur_align_file)
    # print(cons_seq)
    # print(len(cons_seq))

    # # 输入k-mer特征的分布特点
    # matrix_file = '/home/hukang/left_LTR_real_dataset/clean_data/both_ends_frames_clean/negative/Oryza_sativa/Chr8_seg_7_12763883-12775122-lLTR.matrix'
    # # matrix_file = '/home/hukang/left_LTR_real_dataset/clean_data/both_ends_frames_clean/positive/Oryza_sativa/Gypsy-197_OS-LTR.matrix'
    # # matrix_file = '/home/hukang/NeuralLTR/demo/test1/high_copy_frames/Chr449:1527981-1529068.matrix'
    # sample_list = [(matrix_file, 0)]
    #
    # kmer_sizes = config.kmer_sizes
    # for kmer_size in kmer_sizes:
    #     config.kmer_sizes = [kmer_size]
    #     (feature_num_list1, feature_num_list2, label, row_num, matrix_file) = get_matrix_feature_v3(sample_list)[0]
    #
    #     feature_num_list1 = feature_num_list1[:-1]
    #     feature_num_list2 = feature_num_list2[:-1]
    #
    #     # 提取k-mer和对应的频次
    #     k_mers1 = list(range(len(feature_num_list1)))
    #     # k_mers1 = list(generate_kmer_dic(kmer_size).keys())
    #     frequencies1 = list(feature_num_list1)
    #     combined_dict1 = dict(zip(k_mers1, frequencies1))
    #     combined_dict1 = {k: v for k, v in combined_dict1.items() if v != 0}
    #     print(combined_dict1)
    #     top_ten_items = heapq.nlargest(50, combined_dict1.items(), key=lambda item: item[1])
    #     print(top_ten_items)
    #
    #     # 提取k-mer和对应的频次
    #     k_mers2 = list(range(len(feature_num_list2)))
    #     # k_mers2 = list(generate_kmer_dic(kmer_size).keys())
    #     frequencies2 = list(feature_num_list2)
    #     combined_dict2 = dict(zip(k_mers2, frequencies2))
    #     combined_dict2 = {k: v for k, v in combined_dict2.items() if v != 0}
    #     print(combined_dict2)
    #     top_ten_items = heapq.nlargest(50, combined_dict2.items(), key=lambda item: item[1])
    #     print(top_ten_items)
    #
    #     fig, axs = plt.subplots(2, 1, figsize=(8, 10))  # 2行1列的子图布局
    #
    #     # 在第一个子图上绘制k_mer_freq1的散点图
    #     axs[0].scatter(k_mers1, frequencies1, alpha=0.7, color='blue')
    #     axs[0].set_title('k-mer Frequencies - Dataset 1')
    #     axs[0].set_xlabel('k-mers')
    #     axs[0].set_ylabel('Frequencies')
    #
    #     # 在第二个子图上绘制k_mer_freq2的散点图
    #     axs[1].scatter(k_mers2, frequencies2, alpha=0.7, color='green')
    #     axs[1].set_title('k-mer Frequencies - Dataset 2')
    #     axs[1].set_xlabel('k-mers')
    #     axs[1].set_ylabel('Frequencies')
    #
    #     # 调整子图间距
    #     plt.tight_layout()
    #
    #     # 显示图形
    #     plt.show()




    # tmp_output_dir = '/home/hukang/NeuralLTR/demo/test1'
    # # Unpack nested TEs within TEs
    # input1_path = tmp_output_dir + '/confident_ltr.tmp.fa'
    # input2_path = tmp_output_dir + '/confident_ltr.internal.fa'
    # clean_LTR_path = input2_path + '.clean_nested'
    # threads = 40
    # remove_nested_command = 'python3 /home/hukang/NeuralLTR/src/remove_nested_lib.py ' \
    #                         + ' -t ' + str(threads) \
    #                         + ' --tmp_output_dir ' + tmp_output_dir + ' --max_iter_num ' + str(1) \
    #                         + ' --input1 ' + input1_path \
    #                         + ' --input2 ' + input2_path \
    #                         + ' --output ' + clean_LTR_path
    # os.system(remove_nested_command)

    # work_dir = '/home/hukang/LTR_Benchmarking/LTR_libraries/NeuralLTR/rice_MSUv7/LTR_finder'
    # in_file = work_dir + '/rm2_test_redundant/file_work_with1.txt'
    # std_lib = '/home/hukang/NeuralLTR/library/rice.ltr.ref'
    # test_lib = work_dir + '/confident_ltr.redundant.fa'
    # command = 'python /home/hukang/NeuralLTR/tools/connect_BM_RM2.py --in_file '+in_file+' --std_lib '+std_lib+' --test_lib ' + test_lib
    # os.system(command)

    # set1 = {'Gypsy-8B_OS-I#LTR/Gypsy', 'Copia-14B_OS-I#LTR/Copia', 'Gypsy-231_OS-LTR#LTR/Gypsy', 'Gypsy-116_OS-LTR#LTR/Gypsy', 'SZ-10_LTR#LTR/Gypsy', 'CPSC4A_I#LTR/Copia', 'Gypsy-84_OS-I#LTR/Gypsy', 'Gypsy-249_OS-LTR#LTR/Gypsy', 'Gypsy-138_OS-LTR#LTR/Gypsy', 'Copia-108_OS-LTR#LTR/Copia', 'LTR-18F_OS-LTR#LTR/Gypsy', 'Gypsy-130C_OS-LTR#LTR/Gypsy', 'Copia-104_OS-I#LTR/Copia', 'Gypsy-115_OS-LTR#LTR/Gypsy', 'RETRO2_LTR#LTR/Gypsy', 'Copia-96_OS-I#LTR/Copia', 'RETRO2B-I#LTR/Gypsy', 'RETROSAT-3B_LTR#LTR/Gypsy', 'Copia-72_OS-I#LTR/Copia', 'Gypsy-101_OS-LTR#LTR/Gypsy', 'Copia-75_OS-I#LTR/Copia', 'Copia-136_OS-I#LTR/Copia', 'RETRO3D_LTR#LTR/Gypsy', 'RETROSAT4_LTR#LTR/Gypsy', 'COPIA2-I_OS#LTR/Copia', 'OSCOPIA2_I#LTR/Copia', 'Gypsy-46_OS-LTR#LTR/Gypsy', 'Copia-47_OS-I#LTR/Copia', 'COPIA3-LTR_OS#LTR/Copia', 'RETRO3_LTR#LTR/Gypsy', 'Copia-39_OS-LTR#LTR/Copia', 'Copia-118_OS-LTR#LTR/Copia', 'Gypsy-86_OS-LTR#LTR/Gypsy', 'Copia-120_OS-I#LTR/Copia', 'Gypsy-72_OS-LTR#LTR/Gypsy', 'Gypsy-106_OS-I#LTR/Gypsy', 'Gypsy-228_OS-I#LTR/Gypsy', 'SZLTR#LTR/Gypsy', 'RETROFIT4_I#LTR/Copia', 'Copia-26_OS-LTR#LTR/Copia', 'SZ-54C_LTR#LTR/Gypsy', 'Copia-108_OS-I#LTR/Copia', 'Gypsy-25C_OS-LTR#LTR/Gypsy', 'SZ-64_LTR#LTR/Gypsy', 'Copia-18_OS-I#LTR/Copia', 'Copia-110_OS-LTR#LTR/Copia', 'Gypsy-217_OS-I#LTR/Gypsy', 'Gypsy-160_OS-LTR#LTR/Gypsy', 'SZ-22_I#LTR/Gypsy', 'Gypsy-61_OS-I#LTR/Gypsy', 'SZ-8_I#LTR/Gypsy', 'SZ-4_I#LTR/Gypsy', 'SZ-26_LTR#LTR/Gypsy', 'SC-4_I#LTR/Copia', 'Gypsy-249_OS-I#LTR/Gypsy', 'Gypsy-76_OS-I#LTR/Gypsy', 'SZ-57_LTR#LTR', 'Gypsy-46_OS-I#LTR/Gypsy', 'OSCOPIA1_LTR#LTR/Copia', 'SZ-7B-I#LTR/Gypsy', 'BAJIEIN#LTR/Gypsy', 'Gypsy-72_OS-I#LTR/Gypsy', 'Copia-27_OS-I#LTR/Copia', 'Gypsy-248_OS-LTR#LTR/Gypsy', 'RIRE5-LTR_OS#LTR/Copia', 'Copia-76_OS-LTR#LTR/Copia', 'Gypsy-6B_OS-LTR#LTR/Gypsy', 'Gypsy-205_OS-I#LTR/Gypsy', 'CPSC4A_LTR#LTR/Copia', 'Gypsy-248_OS-I#LTR/Gypsy', 'RETROFIT2_I#LTR/Copia', 'SZ-54B_LTR#LTR/Gypsy', 'Copia-52_OS-I#LTR/Copia', 'OSR39_I#LTR/Gypsy', 'RETROFIT6_LTR#LTR/Copia', 'SZ-55_I#LTR/Copia', 'Gypsy-71_OS-LTR#LTR/Gypsy', 'SZ-22_LTR#LTR/Gypsy', 'LTR-11_OS-LTR#LTR', 'RETROSAT2LTRA#LTR/Gypsy', 'Gypsy-165_OS-LTR#LTR/Gypsy', 'RIRE5-I_OS#LTR/Copia', 'CPSC3_LTR#LTR/Copia', 'Copia-44_OS-LTR#LTR/Copia', 'Gypsy-80B_OS-LTR#LTR/Gypsy', 'SZ-25LTR#LTR/Copia', 'Gypsy-250_OS-LTR#LTR/Gypsy', 'Gypsy-86_OS-I#LTR/Gypsy', 'Gypsy-67_OS-I#LTR/Gypsy', 'Copia-114_OS-LTR#LTR/Copia', 'SZ-10_I#LTR/Gypsy', 'RETRO2B-LTR#LTR/Gypsy', 'GYPSO_LTR#LTR/Gypsy', 'Gypsy-115_OS-I#LTR/Gypsy', 'Copia-73_OS-I#LTR/Copia', 'Copia-19_OS-I#LTR/Copia', 'Gypsy-67_OS-LTR#LTR/Gypsy', 'SZ-56_LTR#LTR/Gypsy', 'Copia-103_OS-I#LTR/Copia', 'Gypsy-104_OS-I#LTR/Gypsy', 'GYPSY-B_LTR#LTR/Gypsy', 'Copia-135_OS-LTR#LTR/Copia', 'Gypsy-40_OS-LTR#LTR/Gypsy', 'RETROSOR2_I#LTR/Gypsy', 'Copia-43_OS-I#LTR/Copia', 'RIREXC_I#LTR/Gypsy', 'Copia-97_OS-I#LTR/Copia', 'SC-9_LTR#LTR/Copia', 'Copia-46_OS-LTR#LTR/Copia', 'Copia-49_OS-I#LTR/Copia', 'Copia-13B_OS-LTR#LTR/Copia', 'RETROSOR2_LTR#LTR/Gypsy', 'RETROSAT-3B_I#LTR/Gypsy', 'SZ-50_LTR#LTR/Gypsy', 'Copia-37_OS-LTR#LTR/Copia', 'Gypsy-25_OS-LTR#LTR/Gypsy', 'SZ-66B-LTR#LTR/Gypsy', 'LTR-18K_OS-LTR#LTR', 'Copia-5_OS-LTR#LTR/Copia', 'RETRO3_I#LTR/Gypsy', 'Copia-91_OS-I#LTR/Copia', 'Gypsy-25_OS-I#LTR/Gypsy', 'Gypsy-218_OS-LTR#LTR/Gypsy', 'Gypsy-92_OS-I#LTR/Gypsy', 'Copia-38_OS-LTR#LTR/Copia', 'Gypsy-56_OS-I#LTR/Gypsy', 'Gypsy-76_OS-LTR#LTR/Gypsy', 'Gypsy-91_OS-LTR#LTR/Gypsy', 'TRUNCATOR#LTR/Gypsy', 'SZ-7C-I#LTR/Gypsy', 'COPIA1-I_OS#LTR/Copia', 'Copia-67_OS-LTR#LTR/Copia', 'Copia-55_OS-LTR#LTR/Copia', 'Gypsy-75_OS-I#LTR/Gypsy', 'LTR-18E_OS-LTR#LTR', 'Gypsy-5_OS-I#LTR/Gypsy', 'SZ-59_I#LTR/Gypsy', 'Copia-122_OS-LTR#LTR/Copia', 'Gypsy-154_OS-I#LTR/Gypsy', 'COPIA3-I_OS#LTR/Copia', 'Gypsy-209_OS-LTR#LTR/Gypsy', 'CRM-LTR_OS#LTR', 'Gypsy-87_OS-LTR#LTR/Gypsy', 'Gypsy-135_OS-I#LTR/Gypsy', 'CRM-I_OS#LTR/Gypsy', 'Copia-14_OS-I#LTR/Copia', 'RETROFIT6_I#LTR/Copia', 'SC-6_I#LTR/Copia', 'Gypsy-154_OS-LTR#LTR/Gypsy', 'Gypsy-83_OS-I#LTR/Gypsy', 'RIREXH_I#LTR/Gypsy', 'Copia-87_OS-LTR#LTR/Copia', 'Gypsy-84_OS-LTR#LTR/Gypsy', 'SZ-54D_LTR#LTR/Gypsy', 'RETRO2_I#LTR/Gypsy', 'RIRE3_LTR#LTR/Gypsy', 'Gypsy-7_OS-LTR#LTR/Gypsy', 'Gypsy-212_OS-LTR#LTR/Gypsy', 'Copia-65_OS-I#LTR/Copia', 'Gypsy-96_OS-LTR#LTR/Gypsy', 'OSTONOR1_LTR#LTR/Copia', 'Gypsy-33C_OS-LTR#LTR/Gypsy', 'RIREXF_I#LTR/Gypsy', 'SZ-54A_LTR#LTR/Gypsy', 'Copia-134_OS-LTR#LTR/Copia', 'Gypsy-93_OS-LTR#LTR/Gypsy', 'RIREXG_I#LTR/Gypsy', 'Gypsy-3B_OS-I#LTR/Gypsy', 'Copia-46_OS-I#LTR/Copia', 'RIREX_LTR#LTR/Gypsy', 'Gypsy-225_OS-LTR#LTR/Gypsy', 'Gypsy-105_OS-LTR#LTR/Gypsy', 'Gypsy-182_OS-LTR#LTR/Gypsy', 'Gypsy-66_OS-I#LTR/Gypsy', 'Gypsy-81_OS-I#LTR/Gypsy', 'Copia-13_OS-LTR#LTR/Copia', 'Gypsy-224_OS-LTR#LTR/Gypsy', 'RIREXE_I#LTR/Gypsy', 'SZ-54_LTR#LTR/Gypsy', 'Gypsy-140_OS-I#LTR/Gypsy', 'Gypsy-206_OS-LTR#LTR/Gypsy', 'SZ-55B_LTR#LTR/Copia', 'Gypsy-145_OS-LTR#LTR/Gypsy', 'Copia-62_OS-I#LTR/Copia', 'Gypsy-25G_OS-LTR#LTR/Gypsy', 'Copia-47_OS-LTR#LTR/Copia', 'SC-10B_LTR#LTR/Copia', 'Gypsy-144_OS-I#LTR/Gypsy', 'Gypsy-35_OS-LTR#LTR/Gypsy', 'Gypsy-1_OSJ-I#LTR/Gypsy', 'LTR-33_OS-I#LTR', 'Copia-31_OS-LTR#LTR/Copia', 'SZ-54B_I#LTR/Gypsy', 'Gypsy-175_OS-I#LTR/Gypsy', 'SZ-37_LTR#LTR/Copia', 'Gypsy-225_OS-I#LTR/Gypsy', 'Gypsy-20C_OS-I#LTR/Gypsy', 'Gypsy-199_OS-I#LTR/Gypsy', 'Gypsy-171_OS-I#LTR/Gypsy', 'Copia-93_OS-LTR#LTR/Copia', 'RIRE8A_LTR#LTR/Gypsy', 'Gypsy-160_OS-I#LTR/Gypsy', 'Gypsy-247_OS-LTR#LTR/Gypsy', 'Gypsy-219_OS-I#LTR/Gypsy', 'ATLANTYS-I_OS#LTR/Gypsy', 'RIRE7_I#LTR/Gypsy', 'RETROFIT3_LTR#LTR/Copia', 'Gypsy-209C_OS-I#LTR/Gypsy', 'Copia-9_OS-I#LTR/Copia', 'SC-7_I#LTR/Copia', 'Copia-104_OS-LTR#LTR/Copia', 'Gypsy-6_OS-I#LTR/Gypsy', 'Gypsy-138_OS-I#LTR/Gypsy', 'RETROSAT6_LTR#LTR/Gypsy', 'SZ-33LTR#LTR/Gypsy', 'Copia-57_OS-LTR#LTR/Copia', 'SC-1_LTR#LTR/Copia', 'RETROSAT4_I#LTR/Gypsy', 'Gypsy-172B_OS-LTR#LTR/Gypsy', 'Gypsy-10B_OS-LTR#LTR/Gypsy', 'Gypsy-163_OS-LTR#LTR/Gypsy', 'Gypsy-8B_OS-LTR#LTR/Gypsy', 'Gypsy-118_OS-LTR#LTR/Gypsy', 'OSR39_LTR#LTR/Gypsy', 'GYPSI_I#LTR/Gypsy', 'RETROSAT3_I#LTR/Gypsy', 'Gypsy-217_OS-LTR#LTR/Gypsy', 'Copia-51_OS-LTR#LTR/Copia', 'Copia-16_OS-LTR#LTR/Copia', 'Copia-133_OS-I#LTR/Copia', 'SZ-56A_I#LTR/Gypsy', 'LTR-15_OS-I#LTR', 'GYPSY-A_I#LTR/Gypsy', 'Gypsy-74_OS-LTR#LTR/Gypsy', 'Copia-36_OS-LTR#LTR/Copia', 'Gypsy-113_OS-LTR#LTR/Gypsy', 'Copia-125_OS-LTR#LTR/Copia', 'Copia-76_OS-I#LTR/Copia', 'Gypsy-212_OS-I#LTR/Gypsy', 'Gypsy-5B_OS-I#LTR/Gypsy', 'Gypsy-44_OS-LTR#LTR/Gypsy', 'Copia-58_OS-LTR#LTR/Copia', 'LTR-15B_OS-I#LTR/Copia', 'Copia-51_OS-I#LTR/Copia', 'Copia-72_OS-LTR#LTR/Copia', 'RIRE3_I#LTR/Gypsy', 'Gypsy-127_OS-LTR#LTR/Gypsy', 'Gypsy-251_OS-LTR#LTR/Gypsy', 'SZ-54D_I#LTR/Gypsy', 'Copia-52_OS-LTR#LTR/Copia', 'Gypsy-101_OS-I#LTR/Gypsy', 'Copia-54_OS-I#LTR/Copia', 'OSR42_I#LTR/Gypsy', 'Gypsy-73_OS-I#LTR/Gypsy', 'RIREXB_I#LTR/Gypsy', 'COPI2_LTR#LTR/Copia', 'Gypsy-66_OS-LTR#LTR/Gypsy', 'Gypsy-176_OS-LTR#LTR/Gypsy', 'RIRE8B_I#LTR/Gypsy', 'Gypsy-208_OS-I#LTR/Gypsy', 'Copia-115_OS-I#LTR/Copia', 'SC9A_LTR#LTR/Copia', 'Gypsy-163_OS-I#LTR/Gypsy', 'LTR-18F_OS-I#LTR/Gypsy', 'Gypsy-87_OS-I#LTR/Gypsy', 'Gypsy-264_OS-I#LTR/Gypsy', 'LTR10_OS#LTR', 'OSR38_I#LTR/Gypsy', 'Copia-43_OS-LTR#LTR/Copia', 'Gypsy-146B_OS-LTR#LTR/Gypsy', 'Gypsy-20_OS-LTR#LTR/Gypsy', 'Copia-132_OS-I#LTR/Copia', 'RETROSAT2LTR#LTR/Gypsy', 'RIREXD_I#LTR/Gypsy', 'RIRE8D_I#LTR/Gypsy', 'RETROSAT6_I#LTR/Gypsy', 'Gypsy-81_OS-LTR#LTR/Gypsy', 'Gypsy-145_OS-I#LTR/Gypsy', 'Gypsy-284_OS-I#LTR/Gypsy', 'SC-9_I#LTR/Copia', 'TRUNCATOR2_LTR#LTR/Gypsy', 'Gypsy-80_OS-LTR#LTR/Gypsy', 'Gypsy-213_OS-LTR#LTR/Gypsy', 'Gypsy-276_OS-I#LTR/Gypsy', 'Gypsy-56_OS-LTR#LTR/Gypsy', 'Gypsy-61_OS-LTR#LTR/Gypsy', 'SZ-54C_I#LTR/Gypsy', 'RETROFIT7_I#LTR/Copia', 'Gypsy-188_OS-I#LTR/Gypsy', 'Gypsy-33_OS-LTR#LTR/Gypsy', 'Gypsy-169_OS-LTR#LTR/Gypsy', 'Gypsy-25D_OS-LTR#LTR/Gypsy', 'Gypsy-102_OS-LTR#LTR/Gypsy', 'SZ-7A_LTR#LTR/Gypsy', 'Copia-38_OS-I#LTR/Copia', 'SC-3_I#LTR/Copia', 'COPI2_I#LTR/Copia', 'COPIA2-LTR_OS#LTR/Copia', 'SZ-56A_LTR#LTR/Gypsy', 'Gypsy-3_OS-I#LTR/Gypsy', 'Gypsy-65_OS-LTR#LTR/Gypsy', 'Gypsy-168_OS-LTR#LTR/Gypsy', 'Copia-61_OS-LTR#LTR/Copia', 'Gypsy-125_OS-LTR#LTR/Gypsy', 'Copia-96_OS-LTR#LTR/Copia', 'Copia-90_OS-I#LTR/Copia', 'Copia-90_OS-LTR#LTR/Copia', 'RETROSAT3_LTR#LTR/Gypsy', 'Gypsy-126_OS-LTR#LTR/Gypsy', 'GYPSO_I#LTR/Gypsy', 'Gypsy-90_OS-LTR#LTR/Gypsy', 'Copia-101_OS-I#LTR/Copia', 'rn_179-105_IR#LTR/Copia', 'Gypsy-173_OS-I#LTR/Gypsy', 'RETROSAT2_I#LTR/Gypsy', 'Copia-118_OS-I#LTR/Copia', 'LTR-26_OS-I#LTR/Gypsy', 'OSTONOR-1B_LTR#LTR/Copia', 'RIRE8A_I#LTR/Gypsy', 'Gypsy-102_OS-I#LTR/Gypsy', 'SZ-55C_LTR#LTR/Copia', 'SZ-64_I#LTR/Gypsy', 'LTR-26_OS-LTR#LTR', 'Copia-57_OS-I#LTR/Copia', 'SZ-50_I#LTR/Gypsy', 'Gypsy-89_OS-LTR#LTR/Gypsy', 'Gypsy-179B_OS-LTR#LTR/Gypsy', 'Copia-101_OS-LTR#LTR/Copia', 'Copia-41_OS-LTR#LTR/Copia', 'Copia-19_OS-LTR#LTR/Copia', 'Gypsy-255_OS-LTR#LTR/Gypsy', 'COPI1_I#LTR/Copia', 'Copia-131_OS-LTR#LTR/Copia', 'Gypsy-7B_OS-I#LTR/Gypsy', 'SZ-7B-LTR#LTR/Gypsy', 'Gypsy-25F_OS-LTR#LTR/Gypsy', 'LTR-36_OS-LTR#LTR', 'Copia-65_OS-LTR#LTR/Copia', 'Gypsy-180_OS-LTR#LTR/Gypsy', 'Copia-7_OS-LTR#LTR/Copia', 'Gypsy-126_OS-I#LTR/Gypsy', 'Copia-55_OS-I#LTR/Copia', 'Gypsy-250_OS-I#LTR/Gypsy', 'RIREXJ_I#LTR/Gypsy', 'Gypsy-103_OS-LTR#LTR/Gypsy', 'SZ-63_I#LTR/Gypsy', 'Gypsy-74_OS-I#LTR/Gypsy', 'Gypsy-246_OS-I#LTR/Gypsy', 'Gypsy-118_OS-I#LTR/Gypsy', 'Copia-132_OS-LTR#LTR/Copia', 'Gypsy-254_OS-LTR#LTR/Gypsy', 'Gypsy-95_OS-I#LTR/Gypsy', 'Copia-107_OS-LTR#LTR/Copia', 'Gypsy-47_OS-LTR#LTR/Gypsy', 'Gypsy-111_OS-LTR#LTR/Gypsy', 'Gypsy-254_OS-I#LTR/Gypsy', 'Gypsy-119_OS-LTR#LTR/Gypsy', 'Gypsy-135_OS-LTR#LTR/Gypsy', 'Copia-31_OS-I#LTR/Copia', 'OSTONOR-1B_I#LTR/Copia', 'SZ-31B-LTR#LTR/Gypsy', 'COPIO_I#LTR/Copia', 'SZ-53_LTR#LTR/Gypsy', 'Copia-49_OS-LTR#LTR/Copia', 'Copia-135_OS-I#LTR/Copia', 'Gypsy-8_OS-LTR#LTR/Gypsy', 'Copia-107_OS-I#LTR/Copia', 'SZ-7_I#LTR/Gypsy', 'RETROSAT-2C_LTR#LTR/Gypsy', 'SC-10_LTR#LTR/Copia', 'SZ-46I#LTR/Gypsy', 'LTR-34_OS-LTR#LTR', 'OSCOPIA2_LTR#LTR/Copia', 'Gypsy-95_OS-LTR#LTR/Gypsy', 'Gypsy-40_OS-I#LTR/Gypsy', 'Gypsy-5E_OS-LTR#LTR/Gypsy', 'OSR3_I#LTR/Gypsy', 'Copia-23_OS-LTR#LTR/Copia', 'SZ-59_LTR#LTR/Gypsy', 'Copia-54_OS-LTR#LTR/Copia', 'Copia-37_OS-I#LTR/Copia', 'LTR-22_OS-LTR#LTR', 'Gypsy-110_OS-LTR#LTR/Gypsy', 'Gypsy-108_OS-LTR#LTR/Gypsy', 'Gypsy-111_OS-I#LTR/Gypsy', 'Gypsy-196_OS-I#LTR/Gypsy', 'Copia-75_OS-LTR#LTR/Copia', 'Gypsy-165_OS-I#LTR/Gypsy', 'SZ-55_LTR#LTR/Copia', 'Gypsy-91_OS-I#LTR/Gypsy', 'Gypsy-20B_OS-I#LTR/Gypsy', 'Copia-110_OS-I#LTR/Copia', 'Gypsy-144_OS-LTR#LTR/Gypsy', 'Copia-17_OS-I#LTR/Copia', 'RETRO3C_LTR#LTR/Gypsy', 'SC-8_I#LTR/Copia', 'Gypsy-64_OS-LTR#LTR/Gypsy', 'RETRO2A_LTR#LTR/Gypsy', 'Copia-123_OS-LTR#LTR/Copia', 'Gypsy-8_OS-I#LTR/Gypsy', 'Gypsy-209B_OS-I#LTR/Gypsy', 'Gypsy-218_OS-I#LTR/Gypsy', 'LTR-33_OS-LTR#LTR', 'SZ-37_I#LTR/Copia', 'Copia-136_OS-LTR#LTR/Copia', 'Gypsy-171_OS-LTR#LTR/Gypsy', 'Gypsy-92_OS-LTR#LTR/Gypsy', 'LTR-35_OS-LTR#LTR', 'Gypsy-22_OS-LTR#LTR/Gypsy', 'Copia-16_OS-I#LTR/Copia', 'Gypsy-130B_OS-LTR#LTR/Gypsy', 'RIREXI_I#LTR/Gypsy', 'Gypsy-170_OS-I#LTR/Gypsy', 'Gypsy-107_OS-LTR#LTR/Gypsy', 'RIRE8D_LTR#LTR/Gypsy', 'SZ-54_I#LTR/Gypsy', 'Gypsy-179_OS-I#LTR/Gypsy', 'Gypsy-96_OS-I#LTR/Gypsy', 'COPIA1-LTR_OS#LTR/Copia', 'SC-3_LTR#LTR/Copia', 'SZ-9LTR#LTR/Copia', 'Copia-95_OS-LTR#LTR/Copia', 'LTR-15_OS-LTR#LTR', 'CPR1_LTR#LTR/Copia', 'LTR-22_OS-I#LTR', 'GYPSY1-LTR_OS#LTR/Gypsy', 'Copia-9_OS-LTR#LTR/Copia', 'CPSC2_I#LTR/Copia', 'Gypsy-188_OS-LTR#LTR/Gypsy', 'Gypsy-22_OS-I#LTR/Gypsy', 'BAJIELTR#LTR/Gypsy', 'GYPSY1-I_OS#LTR/Gypsy', 'OSTONOR1_I#LTR/Copia', 'RN12_LTR#LTR/Gypsy', 'Copia-93_OS-I#LTR/Copia', 'Gypsy-47_OS-I#LTR/Gypsy', 'Gypsy-220_OS-I#LTR/Gypsy', 'Gypsy-146B_OS-I#LTR/Gypsy', 'Gypsy-65_OS-I#LTR/Gypsy', 'Gypsy-83_OS-LTR#LTR/Gypsy', 'Gypsy-33_OS-I#LTR/Gypsy', 'RETROSAT5_LTR#LTR/Gypsy', 'Copia-44_OS-I#LTR/Copia', 'SZ-7A_I#LTR/Gypsy', 'COPIO_LTR#LTR/Copia', 'LTR-11C_OS-I#LTR', 'Gypsy-113_OS-I#LTR/Gypsy', 'Copia-131_OS-I#LTR/Copia', 'Copia-114_OS-I#LTR/Copia', 'Gypsy-103_OS-I#LTR/Gypsy', 'Gypsy-170_OS-LTR#LTR/Gypsy', 'Gypsy-224_OS-I#LTR/Gypsy', 'Gypsy-227_OS-LTR#LTR/Gypsy', 'Gypsy-223_OS-LTR#LTR/Gypsy', 'CPSC3_I#LTR/Copia', 'ATLANTYS-LTR_OS#LTR/Gypsy', 'GYPSY-B_I#LTR/Gypsy', 'SZ-36LTR#LTR/Gypsy', 'Copia-95_OS-I#LTR/Copia', 'Copia-67_OS-I#LTR/Copia', 'RIREXD2_I#LTR/Gypsy', 'SZ-54A_I#LTR/Gypsy', 'Gypsy-207_OS-I#LTR/Gypsy', 'LTR-28_OS-LTR#LTR/Gypsy', 'Gypsy-119_OS-I#LTR/Gypsy', 'Copia-64_OS-LTR#LTR/Copia', 'SC-7_LTR#LTR/Copia', 'Gypsy-234_OS-I#LTR/Gypsy', 'SC-10_I#LTR/Copia', 'Gypsy-106_OS-LTR#LTR/Gypsy', 'RIRE3A_LTR#LTR/Gypsy', 'RETROFIT7_LTR#LTR/Copia', 'rn_179-105_LTR#LTR/Copia', 'RIRE2_I#LTR/Gypsy', 'Copia-11_OS-I#LTR/Copia', 'LTR-28_OS-I#LTR/Gypsy', 'OSR3_LTR#LTR/Gypsy', 'Gypsy-73_OS-LTR#LTR/Gypsy', 'Gypsy-25B_OS-LTR#LTR/Gypsy', 'Copia-133_OS-LTR#LTR/Copia', 'Copia-58_OS-I#LTR/Copia', 'SZ-8_LTR#LTR/Gypsy', 'RN12_I#LTR/Gypsy', 'Gypsy-192_OS-I#LTR/Gypsy', 'Copia-81_OS-I#LTR/Copia', 'LTR-38_OS-LTR#LTR', 'Copia-85_OS-LTR#LTR/Copia', 'Gypsy-20D_OS-I#LTR/Gypsy', 'SZ-b-LTR#LTR/Gypsy', 'Copia-39_OS-I#LTR/Copia', 'Gypsy-35_OS-I#LTR/Gypsy', 'Gypsy-178_OS-LTR#LTR/Gypsy', 'Copia-123_OS-I#LTR/Copia', 'RETROFIT3_I#LTR/Copia', 'LTR-34_OS-I#LTR', 'SC-1_I#LTR/Copia', 'Gypsy-21_OS-LTR#LTR/Gypsy', 'Gypsy-77_OS-LTR#LTR/Gypsy', 'Gypsy-223_OS-I#LTR/Gypsy', 'Gypsy-80_OS-I#LTR/Gypsy', 'SZ-17_LTR#LTR/Copia', 'Gypsy-166_OS-LTR#LTR/Gypsy', 'Gypsy-219B_OS-LTR#LTR/Gypsy', 'SZ-7C-LTR#LTR/Gypsy', 'Gypsy-62_OS-LTR#LTR/Gypsy', 'OSCOPIA1_I#LTR/Copia', 'LTR-23B_OS-LTR#LTR', 'Gypsy-28_OS-LTR#LTR/Gypsy', 'RIREXK_I#LTR/Gypsy', 'OSR42_LTR#LTR/Gypsy', 'RETROFIT_I#LTR/Copia', 'Gypsy-231_OS-I#LTR/Gypsy', 'SZ-56_I#LTR/Gypsy', 'Gypsy-104_OS-LTR#LTR/Gypsy', 'GYPSI_LTR#LTR/Gypsy', 'Copia-97_OS-LTR#LTR/Copia', 'Gypsy-199_OS-LTR#LTR/Gypsy', 'Gypsy-75_OS-LTR#LTR/Gypsy', 'Copia-62_OS-LTR#LTR/Copia', 'GYPSY-A_LTR#LTR/Gypsy', 'Copia-63_OS-LTR#LTR/Copia', 'SZ-66C_LTR#LTR/Gypsy', 'Gypsy-196_OS-LTR#LTR/Gypsy', 'Copia-17B_OS-LTR#LTR/Copia', 'Copia-87_OS-I#LTR/Copia', 'OSR38_LTR#LTR/Gypsy', 'SC-6_LTR#LTR/Copia', 'Gypsy-100_OS-LTR#LTR/Gypsy', 'Gypsy-40B_OS-LTR#LTR/Gypsy', 'Copia-103_OS-LTR#LTR/Copia', 'Copia-36_OS-I#LTR/Copia', 'Gypsy-62_OS-I#LTR/Gypsy', 'RETRO2A_I#LTR/Gypsy', 'LTR-36_OS-I#LTR', 'Gypsy-71_OS-I#LTR/Gypsy', 'Gypsy-234_OS-LTR#LTR/Gypsy', 'Gypsy-176_OS-I#LTR/Gypsy', 'SZ-53_I#LTR/Gypsy', 'Gypsy-174_OS-LTR#LTR/Gypsy', 'SC-4_LTR#LTR/Copia', 'Gypsy-25H_OS-LTR#LTR/Gypsy', 'SZ-52_I#LTR/Gypsy', 'Gypsy-90_OS-I#LTR/Gypsy', 'Gypsy-173_OS-LTR#LTR/Gypsy', 'Gypsy-89_OS-I#LTR/Gypsy', 'SZ-52_LTR#LTR/Gypsy', 'Copia-81_OS-LTR#LTR/Copia', 'Gypsy-208_OS-LTR#LTR/Gypsy', 'Gypsy-100_OS-I#LTR/Gypsy', 'Gypsy-284_OS-LTR#LTR/Gypsy', 'Copia-98_OS-LTR#LTR/Copia', 'Gypsy-207_OS-LTR#LTR/Gypsy', 'LTR-27_OS-LTR#LTR/Gypsy', 'Gypsy-180_OS-I#LTR/Gypsy', 'Gypsy-3_OS-LTR#LTR/Gypsy', 'Copia-14_OS-LTR#LTR/Copia', 'Gypsy-215_OS-I#LTR/Gypsy', 'RIRE7_LTR#LTR/Gypsy', 'Gypsy-285_OS-LTR#LTR/Gypsy', 'Gypsy-140_OS-LTR#LTR/Gypsy', 'Gypsy-206_OS-I#LTR/Gypsy', 'RIRE3A_I#LTR/Gypsy', 'CRMA1_I#LTR/Gypsy', 'SC-8B_LTR#LTR/Copia', 'Gypsy-7B_OS-LTR#LTR/Gypsy', 'Gypsy-127_OS-I#LTR/Gypsy', 'Gypsy-116_OS-I#LTR/Gypsy', 'Gypsy-169_OS-I#LTR/Gypsy', 'Gypsy-175_OS-LTR#LTR/Gypsy', 'Gypsy-10B_OS-I#LTR/Gypsy', 'RETROFIT2_LTR#LTR/Copia', 'LTR-35_OS-I#LTR', 'Gypsy-246_OS-LTR#LTR/Gypsy', 'Gypsy-1_OSJ-LTR#LTR/Gypsy', 'Copia-41_OS-I#LTR/Copia', 'Copia-64_OS-I#LTR/Copia', 'LTR-22B_OS-I#LTR', 'Gypsy-125_OS-I#LTR/Gypsy', 'Copia-13_OS-I#LTR/Copia', 'Copia-11B_OS-LTR#LTR/Copia', 'Gypsy-7_OS-I#LTR/Gypsy', 'Copia-120_OS-LTR#LTR/Copia', 'Gypsy-110_OS-I#LTR/Gypsy', 'Gypsy-7C_OS-I#LTR/Gypsy', 'Gypsy-93_OS-I#LTR/Gypsy', 'Gypsy-21_OS-I#LTR/Gypsy', 'Copia-73_OS-LTR#LTR/Copia', 'RETROFIT4_LTR#LTR/Copia', 'Gypsy-192_OS-LTR#LTR/Gypsy', 'Copia-40_OS-LTR#LTR/Copia', 'SZ-63_LTR#LTR/Gypsy', 'Gypsy-77_OS-I#LTR/Gypsy', 'CPSC2_LTR#LTR/Copia', 'Copia-18_OS-LTR#LTR/Copia', 'Gypsy-253_OS-LTR#LTR/Gypsy', 'Copia-61_OS-I#LTR/Copia', 'Gypsy-228_OS-LTR#LTR/Gypsy', 'Gypsy-105_OS-I#LTR/Gypsy', 'SZ-67LTR#LTR', 'SZ-4_LTR#LTR/Gypsy'}
    #
    # set2 = {'SZ-7B-I#LTR/Gypsy', 'Copia-87_OS-LTR#LTR/Copia', 'Gypsy-5_OS-I#LTR/Gypsy', 'SC-3_LTR#LTR/Copia', 'Gypsy-127_OS-I#LTR/Gypsy', 'SZ-54B_I#LTR/Gypsy', 'OSR39_LTR#LTR/Gypsy', 'Copia-41_OS-I#LTR/Copia', 'RIRE3A_LTR#LTR/Gypsy', 'OSCOPIA2_LTR#LTR/Copia', 'RETROSAT-2C_LTR#LTR/Gypsy', 'RETRO3_LTR#LTR/Gypsy', 'Gypsy-228_OS-LTR#LTR/Gypsy', 'Gypsy-234_OS-LTR#LTR/Gypsy', 'Copia-123_OS-LTR#LTR/Copia', 'Gypsy-25_OS-LTR#LTR/Gypsy', 'Copia-55_OS-I#LTR/Copia', 'RIRE5-I_OS#LTR/Copia', 'RETROSAT5_I#LTR/Gypsy', 'Copia-39_OS-I#LTR/Copia', 'LTR-15_OS-LTR#LTR', 'Gypsy-3B_OS-I#LTR/Gypsy', 'SZ-33LTR#LTR/Gypsy', 'Copia-120_OS-LTR#LTR/Copia', 'Gypsy-20D_OS-I#LTR/Gypsy', 'Copia-123_OS-I#LTR/Copia', 'Gypsy-217_OS-LTR#LTR/Gypsy', 'Gypsy-250_OS-LTR#LTR/Gypsy', 'Copia-110_OS-LTR#LTR/Copia', 'Copia-133_OS-LTR#LTR/Copia', 'Copia-108_OS-LTR#LTR/Copia', 'Gypsy-25C_OS-LTR#LTR/Gypsy', 'Copia-52_OS-I#LTR/Copia', 'LTR-36_OS-I#LTR', 'RETROSAT2_I#LTR/Gypsy', 'SC-6_I#LTR/Copia', 'OSR3_LTR#LTR/Gypsy', 'SZ-53_LTR#LTR/Gypsy', 'Gypsy-205_OS-I#LTR/Gypsy', 'LTR-35_OS-LTR#LTR', 'Copia-64_OS-I#LTR/Copia', 'Gypsy-231_OS-LTR#LTR/Gypsy', 'CRM-I_OS#LTR/Gypsy', 'Gypsy-173_OS-I#LTR/Gypsy', 'CPSC2_I#LTR/Copia', 'SZ-57_LTR#LTR', 'RETROSAT-3B_LTR#LTR/Gypsy', 'Gypsy-71_OS-LTR#LTR/Gypsy', 'LTR-22_OS-I#LTR', 'SZ-4_I#LTR/Gypsy', 'Gypsy-22_OS-LTR#LTR/Gypsy', 'SC-4_LTR#LTR/Copia', 'Copia-131_OS-I#LTR/Copia', 'Copia-11_OS-I#LTR/Copia', 'Gypsy-1_OSJ-I#LTR/Gypsy', 'Gypsy-130C_OS-LTR#LTR/Gypsy', 'LTR-26_OS-LTR#LTR', 'SZ-54B_LTR#LTR/Gypsy', 'Copia-81_OS-I#LTR/Copia', 'Gypsy-116_OS-LTR#LTR/Gypsy', 'SZ-64_I#LTR/Gypsy', 'RETROFIT4_LTR#LTR/Copia', 'Gypsy-224_OS-I#LTR/Gypsy', 'Gypsy-100_OS-LTR#LTR/Gypsy', 'Copia-114_OS-I#LTR/Copia', 'Gypsy-46_OS-I#LTR/Gypsy', 'OSTONOR-1B_I#LTR/Copia', 'Gypsy-119_OS-I#LTR/Gypsy', 'Gypsy-223_OS-LTR#LTR/Gypsy', 'SZ-55_I#LTR/Copia', 'COPIA2-I_OS#LTR/Copia', 'Gypsy-56_OS-LTR#LTR/Gypsy', 'Copia-76_OS-I#LTR/Copia', 'Gypsy-83_OS-I#LTR/Gypsy', 'Copia-13_OS-LTR#LTR/Copia', 'Gypsy-77_OS-LTR#LTR/Gypsy', 'Gypsy-178_OS-LTR#LTR/Gypsy', 'Gypsy-170_OS-I#LTR/Gypsy', 'Copia-19_OS-LTR#LTR/Copia', 'SC-4_I#LTR/Copia', 'Copia-104_OS-LTR#LTR/Copia', 'Copia-51_OS-I#LTR/Copia', 'Copia-44_OS-LTR#LTR/Copia', 'Gypsy-171_OS-I#LTR/Gypsy', 'Gypsy-145_OS-I#LTR/Gypsy', 'SZ-31B-LTR#LTR/Gypsy', 'SZ-46I#LTR/Gypsy', 'Copia-47_OS-I#LTR/Copia', 'Gypsy-5E_OS-LTR#LTR/Gypsy', 'SZ-22_I#LTR/Gypsy', 'RETRO2_LTR#LTR/Gypsy', 'SZ-54D_I#LTR/Gypsy', 'RETROSAT3_I#LTR/Gypsy', 'Copia-134_OS-LTR#LTR/Copia', 'Copia-41_OS-LTR#LTR/Copia', 'Copia-31_OS-I#LTR/Copia', 'Copia-76_OS-LTR#LTR/Copia', 'COPIA2-LTR_OS#LTR/Copia', 'LTR-33_OS-I#LTR', 'Gypsy-154_OS-I#LTR/Gypsy', 'SZ-8_LTR#LTR/Gypsy', 'SC-8_I#LTR/Copia', 'CRM-LTR_OS#LTR', 'SC-7_LTR#LTR/Copia', 'Gypsy-209_OS-LTR#LTR/Gypsy', 'Gypsy-6_OS-I#LTR/Gypsy', 'Copia-118_OS-I#LTR/Copia', 'Gypsy-248_OS-I#LTR/Gypsy', 'Gypsy-6B_OS-LTR#LTR/Gypsy', 'Copia-118_OS-LTR#LTR/Copia', 'LTR-18F_OS-LTR#LTR/Gypsy', 'Gypsy-249_OS-I#LTR/Gypsy', 'RIRE2_I#LTR/Gypsy', 'Copia-7_OS-LTR#LTR/Copia', 'Gypsy-218_OS-LTR#LTR/Gypsy', 'LTR-22_OS-LTR#LTR', 'Gypsy-212_OS-I#LTR/Gypsy', 'Gypsy-192_OS-I#LTR/Gypsy', 'Gypsy-125_OS-I#LTR/Gypsy', 'Gypsy-101_OS-LTR#LTR/Gypsy', 'CPR1_LTR#LTR/Copia', 'Gypsy-248_OS-LTR#LTR/Gypsy', 'RETROFIT6_I#LTR/Copia', 'LTR-18E_OS-LTR#LTR', 'RIREXE_I#LTR/Gypsy', 'SZ-56A_I#LTR/Gypsy', 'TRUNCATOR2_LTR#LTR/Gypsy', 'Copia-136_OS-LTR#LTR/Copia', 'Copia-101_OS-I#LTR/Copia', 'RIREX_LTR#LTR/Gypsy', 'Gypsy-175_OS-LTR#LTR/Gypsy', 'Copia-93_OS-I#LTR/Copia', 'Gypsy-138_OS-LTR#LTR/Gypsy', 'SZ-56_LTR#LTR/Gypsy', 'LTR-11C_OS-I#LTR', 'Copia-9_OS-I#LTR/Copia', 'Gypsy-76_OS-I#LTR/Gypsy', 'OSCOPIA1_I#LTR/Copia', 'Copia-96_OS-I#LTR/Copia', 'Gypsy-106_OS-I#LTR/Gypsy', 'SC-9_LTR#LTR/Copia', 'GYPSO_I#LTR/Gypsy', 'Copia-110_OS-I#LTR/Copia', 'Gypsy-110_OS-LTR#LTR/Gypsy', 'Copia-43_OS-I#LTR/Copia', 'Copia-9_OS-LTR#LTR/Copia', 'LTR-38_OS-LTR#LTR', 'Gypsy-225_OS-LTR#LTR/Gypsy', 'Copia-57_OS-LTR#LTR/Copia', 'SZ-63_I#LTR/Gypsy', 'Copia-36_OS-LTR#LTR/Copia', 'RIRE8A_LTR#LTR/Gypsy', 'GYPSY-A_I#LTR/Gypsy', 'SZ-10_I#LTR/Gypsy', 'RETROFIT4_I#LTR/Copia', 'SC-1_LTR#LTR/Copia', 'Gypsy-165_OS-LTR#LTR/Gypsy', 'Gypsy-253_OS-LTR#LTR/Gypsy', 'Gypsy-86_OS-I#LTR/Gypsy', 'GYPSY-A_LTR#LTR/Gypsy', 'RIRE3_I#LTR/Gypsy', 'Gypsy-35_OS-LTR#LTR/Gypsy', 'RIRE8A_I#LTR/Gypsy', 'GYPSY-B_I#LTR/Gypsy', 'SZ-b-LTR#LTR/Gypsy', 'RETRO3_I#LTR/Gypsy', 'RETROSAT6_LTR#LTR/Gypsy', 'RIREXH_I#LTR/Gypsy', 'Copia-81_OS-LTR#LTR/Copia', 'Copia-58_OS-LTR#LTR/Copia', 'LTR-15B_OS-I#LTR/Copia', 'Gypsy-246_OS-LTR#LTR/Gypsy', 'Gypsy-145_OS-LTR#LTR/Gypsy', 'OSR38_I#LTR/Gypsy', 'Gypsy-166_OS-LTR#LTR/Gypsy', 'RETROFIT_I#LTR/Copia', 'RETRO3C_LTR#LTR/Gypsy', 'LTR-33_OS-LTR#LTR', 'Gypsy-75_OS-I#LTR/Gypsy', 'Copia-38_OS-I#LTR/Copia', 'Gypsy-25F_OS-LTR#LTR/Gypsy', 'RETRO2_I#LTR/Gypsy', 'Gypsy-140_OS-I#LTR/Gypsy', 'Gypsy-135_OS-LTR#LTR/Gypsy', 'Gypsy-103_OS-I#LTR/Gypsy', 'Gypsy-44_OS-LTR#LTR/Gypsy', 'Copia-38_OS-LTR#LTR/Copia', 'Gypsy-91_OS-I#LTR/Gypsy', 'Gypsy-10B_OS-I#LTR/Gypsy', 'Gypsy-1_OSJ-LTR#LTR/Gypsy', 'Copia-67_OS-LTR#LTR/Copia', 'Copia-44_OS-I#LTR/Copia', 'SC-10B_LTR#LTR/Copia', 'Gypsy-107_OS-LTR#LTR/Gypsy', 'Copia-97_OS-LTR#LTR/Copia', 'COPIA3-I_OS#LTR/Copia', 'OSR38_LTR#LTR/Gypsy', 'SZLTR#LTR/Gypsy', 'Gypsy-160_OS-LTR#LTR/Gypsy', 'RIRE8D_LTR#LTR/Gypsy', 'SZ-54A_I#LTR/Gypsy', 'COPIA3-LTR_OS#LTR/Copia', 'RETROFIT7_I#LTR/Copia', 'Gypsy-56_OS-I#LTR/Gypsy', 'OSR39_I#LTR/Gypsy', 'Gypsy-223_OS-I#LTR/Gypsy', 'LTR-34_OS-I#LTR', 'LTR-18F_OS-I#LTR/Gypsy', 'Gypsy-135_OS-I#LTR/Gypsy', 'Gypsy-81_OS-I#LTR/Gypsy', 'Gypsy-284_OS-LTR#LTR/Gypsy', 'Gypsy-111_OS-LTR#LTR/Gypsy', 'Gypsy-102_OS-LTR#LTR/Gypsy', 'Gypsy-127_OS-LTR#LTR/Gypsy', 'OSR42_LTR#LTR/Gypsy', 'Gypsy-95_OS-LTR#LTR/Gypsy', 'Gypsy-5B_OS-I#LTR/Gypsy', 'RETROSAT6_I#LTR/Gypsy', 'Copia-75_OS-LTR#LTR/Copia', 'Copia-95_OS-I#LTR/Copia', 'Gypsy-180_OS-LTR#LTR/Gypsy', 'Gypsy-206_OS-I#LTR/Gypsy', 'Copia-49_OS-LTR#LTR/Copia', 'Copia-47_OS-LTR#LTR/Copia', 'Copia-13B_OS-LTR#LTR/Copia', 'Gypsy-20_OS-LTR#LTR/Gypsy', 'RETRO3D_LTR#LTR/Gypsy', 'Gypsy-255_OS-LTR#LTR/Gypsy', 'Copia-5_OS-LTR#LTR/Copia', 'Copia-72_OS-LTR#LTR/Copia', 'RETROSAT5_LTR#LTR/Gypsy', 'RIRE7_LTR#LTR/Gypsy', 'Gypsy-96_OS-LTR#LTR/Gypsy', 'Gypsy-81_OS-LTR#LTR/Gypsy', 'Copia-39_OS-LTR#LTR/Copia', 'Gypsy-101_OS-I#LTR/Gypsy', 'SZ-37_LTR#LTR/Copia', 'RETROFIT3_LTR#LTR/Copia', 'Gypsy-66_OS-I#LTR/Gypsy', 'OSCOPIA1_LTR#LTR/Copia', 'Gypsy-108_OS-LTR#LTR/Gypsy', 'RETROFIT2_I#LTR/Copia', 'Copia-14_OS-LTR#LTR/Copia', 'RETROSAT4_LTR#LTR/Gypsy', 'Gypsy-33_OS-LTR#LTR/Gypsy', 'Gypsy-93_OS-LTR#LTR/Gypsy', 'Gypsy-67_OS-I#LTR/Gypsy', 'Copia-36_OS-I#LTR/Copia', 'RETROSAT2LTRA#LTR/Gypsy', 'Gypsy-113_OS-I#LTR/Gypsy', 'Copia-43_OS-LTR#LTR/Copia', 'Gypsy-61_OS-I#LTR/Gypsy', 'GYPSY1-I_OS#LTR/Gypsy', 'Gypsy-25G_OS-LTR#LTR/Gypsy', 'SZ-52_I#LTR/Gypsy', 'SZ-56A_LTR#LTR/Gypsy', 'Gypsy-3_OS-LTR#LTR/Gypsy', 'RETROSAT2LTR#LTR/Gypsy', 'Gypsy-76_OS-LTR#LTR/Gypsy', 'Gypsy-102_OS-I#LTR/Gypsy', 'Gypsy-218_OS-I#LTR/Gypsy', 'Gypsy-73_OS-LTR#LTR/Gypsy', 'OSTONOR-1B_LTR#LTR/Copia', 'Copia-31_OS-LTR#LTR/Copia', 'Copia-95_OS-LTR#LTR/Copia', 'RN12_I#LTR/Gypsy', 'LTR-35_OS-I#LTR', 'Copia-16_OS-I#LTR/Copia', 'SC-7_I#LTR/Copia', 'Copia-17B_OS-LTR#LTR/Copia', 'Gypsy-146B_OS-I#LTR/Gypsy', 'Copia-46_OS-I#LTR/Copia', 'Gypsy-207_OS-I#LTR/Gypsy', 'SZ-59_LTR#LTR/Gypsy', 'RN12_LTR#LTR/Gypsy', 'Gypsy-118_OS-I#LTR/Gypsy', 'Copia-108_OS-I#LTR/Copia', 'RETRO2B-LTR#LTR/Gypsy', 'Copia-98_OS-LTR#LTR/Copia', 'Gypsy-188_OS-LTR#LTR/Gypsy', 'Gypsy-80_OS-LTR#LTR/Gypsy', 'Copia-23_OS-LTR#LTR/Copia', 'rn_179-105_LTR#LTR/Copia', 'Gypsy-92_OS-I#LTR/Gypsy', 'SC9A_LTR#LTR/Copia', 'Gypsy-140_OS-LTR#LTR/Gypsy', 'Gypsy-110_OS-I#LTR/Gypsy', 'Gypsy-66_OS-LTR#LTR/Gypsy', 'Gypsy-160_OS-I#LTR/Gypsy', 'SZ-54_LTR#LTR/Gypsy', 'TRUNCATOR#LTR/Gypsy', 'Gypsy-10B_OS-LTR#LTR/Gypsy', 'Copia-49_OS-I#LTR/Copia', 'Gypsy-170_OS-LTR#LTR/Gypsy', 'Copia-58_OS-I#LTR/Copia', 'Gypsy-285_OS-LTR#LTR/Gypsy', 'Gypsy-87_OS-LTR#LTR/Gypsy', 'Gypsy-126_OS-I#LTR/Gypsy', 'Gypsy-199_OS-LTR#LTR/Gypsy', 'Copia-37_OS-LTR#LTR/Copia', 'Copia-14B_OS-I#LTR/Copia', 'Gypsy-103_OS-LTR#LTR/Gypsy', 'Gypsy-246_OS-I#LTR/Gypsy', 'COPIO_LTR#LTR/Copia', 'SZ-7B-LTR#LTR/Gypsy', 'Gypsy-119_OS-LTR#LTR/Gypsy', 'COPIO_I#LTR/Copia', 'Gypsy-25D_OS-LTR#LTR/Gypsy', 'Gypsy-83_OS-LTR#LTR/Gypsy', 'SZ-53_I#LTR/Gypsy', 'Gypsy-89_OS-I#LTR/Gypsy', 'Gypsy-176_OS-I#LTR/Gypsy', 'Gypsy-192_OS-LTR#LTR/Gypsy', 'Copia-55_OS-LTR#LTR/Copia', 'Gypsy-180_OS-I#LTR/Gypsy', 'Gypsy-169_OS-I#LTR/Gypsy', 'Gypsy-224_OS-LTR#LTR/Gypsy', 'GYPSI_I#LTR/Gypsy', 'Gypsy-65_OS-I#LTR/Gypsy', 'SZ-10_LTR#LTR/Gypsy', 'CPSC4A_LTR#LTR/Copia', 'Copia-96_OS-LTR#LTR/Copia', 'Gypsy-40B_OS-LTR#LTR/Gypsy', 'Copia-73_OS-LTR#LTR/Copia', 'Gypsy-25B_OS-LTR#LTR/Gypsy', 'Gypsy-118_OS-LTR#LTR/Gypsy', 'Copia-103_OS-I#LTR/Copia', 'RIREXJ_I#LTR/Gypsy', 'Copia-131_OS-LTR#LTR/Copia', 'LTR-26_OS-I#LTR/Gypsy', 'Gypsy-62_OS-I#LTR/Gypsy', 'Gypsy-196_OS-LTR#LTR/Gypsy', 'Gypsy-247_OS-LTR#LTR/Gypsy', 'Gypsy-95_OS-I#LTR/Gypsy', 'Copia-75_OS-I#LTR/Copia', 'Gypsy-96_OS-I#LTR/Gypsy', 'SZ-37_I#LTR/Copia', 'RIRE3_LTR#LTR/Gypsy', 'Gypsy-89_OS-LTR#LTR/Gypsy', 'Copia-67_OS-I#LTR/Copia', 'SZ-56_I#LTR/Gypsy', 'ATLANTYS-I_OS#LTR/Gypsy', 'Gypsy-115_OS-LTR#LTR/Gypsy', 'COPI2_I#LTR/Copia', 'Copia-135_OS-LTR#LTR/Copia', 'Copia-122_OS-LTR#LTR/Copia', 'SC-6_LTR#LTR/Copia', 'CPSC4A_I#LTR/Copia', 'Gypsy-163_OS-I#LTR/Gypsy', 'Gypsy-21_OS-LTR#LTR/Gypsy', 'Copia-57_OS-I#LTR/Copia', 'CPSC3_LTR#LTR/Copia', 'COPI1_I#LTR/Copia', 'Gypsy-73_OS-I#LTR/Gypsy', 'Gypsy-171_OS-LTR#LTR/Gypsy', 'Copia-132_OS-LTR#LTR/Copia', 'SZ-54D_LTR#LTR/Gypsy', 'Gypsy-8B_OS-LTR#LTR/Gypsy', 'LTR-28_OS-LTR#LTR/Gypsy', 'Copia-26_OS-LTR#LTR/Copia', 'Copia-14_OS-I#LTR/Copia', 'Gypsy-80_OS-I#LTR/Gypsy', 'Gypsy-168_OS-LTR#LTR/Gypsy', 'SC-10_I#LTR/Copia', 'Copia-64_OS-LTR#LTR/Copia', 'Copia-65_OS-LTR#LTR/Copia', 'LTR-27_OS-LTR#LTR/Gypsy', 'Gypsy-105_OS-I#LTR/Gypsy', 'Gypsy-74_OS-LTR#LTR/Gypsy', 'LTR10_OS#LTR', 'OSTONOR1_LTR#LTR/Copia', 'Copia-61_OS-I#LTR/Copia', 'Gypsy-115_OS-I#LTR/Gypsy', 'Copia-97_OS-I#LTR/Copia', 'Gypsy-47_OS-I#LTR/Gypsy', 'Copia-115_OS-I#LTR/Copia', 'Gypsy-182_OS-LTR#LTR/Gypsy', 'Gypsy-21_OS-I#LTR/Gypsy', 'Gypsy-72_OS-I#LTR/Gypsy', 'Gypsy-249_OS-LTR#LTR/Gypsy', 'Gypsy-217_OS-I#LTR/Gypsy', 'Copia-107_OS-I#LTR/Copia', 'SZ-54C_I#LTR/Gypsy', 'Gypsy-62_OS-LTR#LTR/Gypsy', 'Copia-101_OS-LTR#LTR/Copia', 'Gypsy-254_OS-I#LTR/Gypsy', 'Copia-93_OS-LTR#LTR/Copia', 'LTR-23B_OS-LTR#LTR', 'SZ-17_LTR#LTR/Copia', 'COPI2_LTR#LTR/Copia', 'OSTONOR1_I#LTR/Copia', 'RETROFIT7_LTR#LTR/Copia', 'Copia-63_OS-LTR#LTR/Copia', 'RETROSOR2_I#LTR/Gypsy', 'Gypsy-228_OS-I#LTR/Gypsy', 'SZ-67LTR#LTR', 'Gypsy-104_OS-LTR#LTR/Gypsy', 'Gypsy-188_OS-I#LTR/Gypsy', 'Copia-91_OS-I#LTR/Copia', 'SZ-55C_LTR#LTR/Copia', 'Copia-46_OS-LTR#LTR/Copia', 'Gypsy-91_OS-LTR#LTR/Gypsy', 'Copia-85_OS-LTR#LTR/Copia', 'Gypsy-179B_OS-LTR#LTR/Gypsy', 'Gypsy-28_OS-LTR#LTR/Gypsy', 'RIRE7_I#LTR/Gypsy', 'Gypsy-75_OS-LTR#LTR/Gypsy', 'SZ-54_I#LTR/Gypsy', 'SZ-64_LTR#LTR/Gypsy', 'Gypsy-125_OS-LTR#LTR/Gypsy', 'Gypsy-206_OS-LTR#LTR/Gypsy', 'Gypsy-231_OS-I#LTR/Gypsy', 'Gypsy-207_OS-LTR#LTR/Gypsy', 'SZ-9LTR#LTR/Copia', 'Copia-104_OS-I#LTR/Copia', 'SZ-7C-I#LTR/Gypsy', 'SC-8B_LTR#LTR/Copia', 'GYPSI_LTR#LTR/Gypsy', 'Copia-17_OS-I#LTR/Copia', 'Gypsy-146B_OS-LTR#LTR/Gypsy', 'Copia-136_OS-I#LTR/Copia', 'Gypsy-173_OS-LTR#LTR/Gypsy', 'SZ-7A_LTR#LTR/Gypsy', 'Gypsy-33C_OS-LTR#LTR/Gypsy', 'RETROFIT3_I#LTR/Copia', 'Gypsy-22_OS-I#LTR/Gypsy', 'Copia-125_OS-LTR#LTR/Copia', 'Gypsy-47_OS-LTR#LTR/Gypsy', 'Gypsy-264_OS-I#LTR/Gypsy', 'Copia-72_OS-I#LTR/Copia', 'Gypsy-72_OS-LTR#LTR/Gypsy', 'Gypsy-61_OS-LTR#LTR/Gypsy', 'LTR-22B_OS-I#LTR', 'LTR-28_OS-I#LTR/Gypsy', 'Gypsy-251_OS-LTR#LTR/Gypsy', 'Copia-103_OS-LTR#LTR/Copia', 'Gypsy-144_OS-I#LTR/Gypsy', 'Copia-90_OS-I#LTR/Copia', 'Gypsy-40_OS-I#LTR/Gypsy', 'Gypsy-80B_OS-LTR#LTR/Gypsy', 'Gypsy-77_OS-I#LTR/Gypsy', 'SC-3_I#LTR/Copia', 'Gypsy-90_OS-LTR#LTR/Gypsy', 'CPSC2_LTR#LTR/Copia', 'COPIA1-I_OS#LTR/Copia', 'Gypsy-8_OS-LTR#LTR/Gypsy', 'GYPSY1-LTR_OS#LTR/Gypsy', 'Copia-40_OS-LTR#LTR/Copia', 'Copia-52_OS-LTR#LTR/Copia', 'OSCOPIA2_I#LTR/Copia', 'SZ-22_LTR#LTR/Gypsy', 'Gypsy-199_OS-I#LTR/Gypsy', 'Gypsy-126_OS-LTR#LTR/Gypsy', 'SZ-26_LTR#LTR/Gypsy', 'Copia-51_OS-LTR#LTR/Copia', 'Copia-107_OS-LTR#LTR/Copia', 'SZ-7C-LTR#LTR/Gypsy', 'Copia-65_OS-I#LTR/Copia', 'Gypsy-74_OS-I#LTR/Gypsy', 'Gypsy-90_OS-I#LTR/Gypsy', 'RETROSAT4_I#LTR/Gypsy', 'Copia-27_OS-I#LTR/Copia', 'Gypsy-165_OS-I#LTR/Gypsy', 'Gypsy-130B_OS-LTR#LTR/Gypsy', 'Copia-11B_OS-LTR#LTR/Copia', 'Gypsy-138_OS-I#LTR/Gypsy', 'Gypsy-46_OS-LTR#LTR/Gypsy', 'Gypsy-172B_OS-LTR#LTR/Gypsy', 'RETROSAT3_LTR#LTR/Gypsy', 'Gypsy-234_OS-I#LTR/Gypsy', 'Gypsy-212_OS-LTR#LTR/Gypsy', 'Gypsy-225_OS-I#LTR/Gypsy', 'SZ-55_LTR#LTR/Copia', 'Gypsy-111_OS-I#LTR/Gypsy', 'Gypsy-208_OS-I#LTR/Gypsy', 'Gypsy-254_OS-LTR#LTR/Gypsy', 'Copia-54_OS-LTR#LTR/Copia', 'Gypsy-219B_OS-LTR#LTR/Gypsy', 'SZ-4_LTR#LTR/Gypsy', 'Copia-61_OS-LTR#LTR/Copia', 'Gypsy-144_OS-LTR#LTR/Gypsy', 'Gypsy-7_OS-LTR#LTR/Gypsy', 'Gypsy-208_OS-LTR#LTR/Gypsy', 'Copia-62_OS-LTR#LTR/Copia', 'Copia-18_OS-I#LTR/Copia', 'RIRE5-LTR_OS#LTR/Copia', 'SZ-50_LTR#LTR/Gypsy', 'ATLANTYS-LTR_OS#LTR/Gypsy', 'LTR-36_OS-LTR#LTR', 'Gypsy-284_OS-I#LTR/Gypsy', 'Copia-54_OS-I#LTR/Copia', 'Gypsy-106_OS-LTR#LTR/Gypsy', 'Gypsy-7B_OS-LTR#LTR/Gypsy', 'Copia-90_OS-LTR#LTR/Copia', 'GYPSO_LTR#LTR/Gypsy', 'Gypsy-86_OS-LTR#LTR/Gypsy', 'Gypsy-104_OS-I#LTR/Gypsy', 'Copia-135_OS-I#LTR/Copia', 'RETRO2A_LTR#LTR/Gypsy', 'SZ-8_I#LTR/Gypsy', 'Copia-18_OS-LTR#LTR/Copia', 'Gypsy-67_OS-LTR#LTR/Gypsy', 'LTR-27_OS-I#LTR', 'Gypsy-196_OS-I#LTR/Gypsy', 'SZ-52_LTR#LTR/Gypsy', 'Gypsy-87_OS-I#LTR/Gypsy', 'Gypsy-163_OS-LTR#LTR/Gypsy', 'LTR-15_OS-I#LTR', 'Gypsy-65_OS-LTR#LTR/Gypsy', 'Gypsy-8B_OS-I#LTR/Gypsy', 'SZ-66B-LTR#LTR/Gypsy', 'Copia-16_OS-LTR#LTR/Copia', 'OSR42_I#LTR/Gypsy', 'Gypsy-40_OS-LTR#LTR/Gypsy', 'Gypsy-213_OS-LTR#LTR/Gypsy', 'SC-10_LTR#LTR/Copia', 'Gypsy-35_OS-I#LTR/Gypsy', 'Gypsy-8_OS-I#LTR/Gypsy', 'Gypsy-93_OS-I#LTR/Gypsy', 'OSR3_I#LTR/Gypsy', 'Gypsy-175_OS-I#LTR/Gypsy', 'SZ-63_LTR#LTR/Gypsy', 'Gypsy-116_OS-I#LTR/Gypsy', 'Gypsy-154_OS-LTR#LTR/Gypsy', 'Copia-37_OS-I#LTR/Copia', 'BAJIELTR#LTR/Gypsy', 'SZ-66C_LTR#LTR/Gypsy', 'Gypsy-84_OS-I#LTR/Gypsy', 'Gypsy-174_OS-LTR#LTR/Gypsy', 'SZ-55B_LTR#LTR/Copia', 'Gypsy-64_OS-LTR#LTR/Gypsy', 'Gypsy-227_OS-LTR#LTR/Gypsy', 'RETROSOR2_LTR#LTR/Gypsy', 'LTR-18K_OS-LTR#LTR', 'SZ-36LTR#LTR/Gypsy', 'Gypsy-179_OS-I#LTR/Gypsy', 'Gypsy-105_OS-LTR#LTR/Gypsy', 'SC-9_I#LTR/Copia', 'Gypsy-176_OS-LTR#LTR/Gypsy', 'LTR-11_OS-LTR#LTR', 'SZ-59_I#LTR/Gypsy', 'RETROFIT2_LTR#LTR/Copia', 'Gypsy-169_OS-LTR#LTR/Gypsy', 'SZ-54C_LTR#LTR/Gypsy', 'GYPSY-B_LTR#LTR/Gypsy', 'rn_179-105_IR#LTR/Copia', 'Copia-62_OS-I#LTR/Copia', 'Copia-114_OS-LTR#LTR/Copia', 'COPIA1-LTR_OS#LTR/Copia', 'Gypsy-25H_OS-LTR#LTR/Gypsy', 'LTR-34_OS-LTR#LTR', 'Gypsy-113_OS-LTR#LTR/Gypsy', 'SZ-54A_LTR#LTR/Gypsy', 'Gypsy-92_OS-LTR#LTR/Gypsy', 'RETROFIT6_LTR#LTR/Copia', 'Copia-73_OS-I#LTR/Copia', 'RIRE8B_I#LTR/Gypsy', 'Gypsy-100_OS-I#LTR/Gypsy', 'SZ-25LTR#LTR/Copia', 'Gypsy-84_OS-LTR#LTR/Gypsy'}
    #
    # # A - B
    # diff1 = set1.difference(set2)
    #
    # # B - A
    # diff2 = set2.difference(set1)
    #
    # # 输出结果
    # print("存在于 set1 但不存在于 set2 的元素:")
    # print(diff1)
    # print(len(diff1))
    #
    # print("\n存在于 set2 但不存在于 set1 的元素:")
    # print(diff2)
    # print(len(diff2))

    # work_dir = '/home/hukang/NeuralLTR/demo/test1'
    # left_ltr = work_dir + '/left_LTR.fa'
    # contignames, contigs = read_fasta(left_ltr)
    #
    # ltr_left_frames = work_dir + '/ltr_left_frames'
    # seq_names = set()
    # for name in os.listdir(ltr_left_frames):
    #     seq_name = name.split('.')[0]
    #     seq_names.add(seq_name)
    #
    # miss_path = work_dir + '/left_ltr_miss.fa'
    # missing_contigs = {}
    # for name in contignames:
    #     if name not in seq_names:
    #         missing_contigs[name] = contigs[name]
    # store_fasta(missing_contigs, miss_path)
    #
    # reference = '/home/hukang/NeuralLTR/demo/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fa'
    # threads = 40
    # temp_dir = work_dir + '/candidate_ltr_miss'
    # output_dir = work_dir + '/ltr_left_frames_miss'
    # split_ref_dir = work_dir + '/ref_chr'
    # generate_both_ends_frame_from_seq(miss_path, reference, threads, temp_dir, output_dir, split_ref_dir)

    # scn1 = '/home/hukang/NeuralLTR/demo/test5/confident_ltr.scn'
    # scn2 = '/home/hukang/NeuralLTR/demo/test4/confident_ltr.scn'
    # ltr_candidates1, ltr_lines1 = read_scn(scn1, None)
    # ltr_candidates2, ltr_lines2 = read_scn(scn2, None)
    # ltr_lines1 = list(ltr_lines1.values())
    # ltr_lines2 = list(ltr_lines2.values())
    # count = 0
    # for line in ltr_lines1:
    #     if line not in ltr_lines2:
    #         print(line)
    #         count += 1
    # print(count)

    # confident_ltr_internal = ''
    # deredundant_for_LTR(confident_ltr_internal, tmp_output_dir, threads)

    # align_file = '/homec/xuminghua/NeuralLTR/demo/test1/candidate_ltr/Chr448:23436931-23437224.blast.bed.fa.maf.fa'
    # output_dir = '/homec/xuminghua/NeuralLTR/demo/test1/test_dir'
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)
    # debug = 1
    # both_end_frame_path = get_both_ends_frame(query_name, cur_seq, align_file, output_dir, debug)


    # matrix_file = '/home/hukang/NeuralLTR/demo/test3/ltr_both_ends_frames/Chr451:525385-525807.matrix'
    # ltr_name = ''
    # is_TIR = is_TIR_frame(matrix_file, ltr_name)
    # print(is_TIR)


    # work_dir = '/home/hukang/NeuralLTR/demo/test3'
    # blastn2Results_path = work_dir + '/test.out'
    # full_length_out = work_dir + '/test.fl.out'
    # split_repeats_path = work_dir + '/test.fa'
    # genome_path = work_dir + '/genome.rename.fa'
    # tmp_dir = work_dir + '/test_tmp'
    # tools_dir = '/home/hukang/NeuralLTR/tools'
    # coverage_threshold = 0.95
    # category = 'LTR'
    # lines = generate_full_length_out(blastn2Results_path, full_length_out, split_repeats_path, genome_path, tmp_dir,
    #                                  tools_dir,
    #                                  coverage_threshold, category)
    # print(lines)


    # # 我们现在尝试划分 子窗口，计算每个子窗口中随机簇（指的是簇大小为1）的数量
    # threads = 40
    # redundant_ltr = '/home/hukang/NeuralLTR/demo/test15/confident_ltr.internal.fa'
    # work_dir = '/home/hukang/NeuralLTR/demo/test15'
    # deredundant_for_LTR(redundant_ltr, work_dir, threads)

    # cur_align_file = '/home/hukang/NeuralLTR/demo/test3/raw_ltr_cluster/Ninja_7/1.fa.maf.fa'
    # cons_seq = cons_from_mafft(cur_align_file)
    # print(cons_seq)

    # reference = '/home/hukang/NeuralLTR/demo/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fa'
    # scn_file = '/home/hukang/NeuralLTR/demo/test3/test.scn'
    # ref_names, ref_contigs = read_fasta(reference)
    # ltr_candidates, ltr_lines = read_scn(scn_file)
    #
    # for candidate_index in ltr_candidates.keys():
    #     (chr_name, left_ltr_start, left_ltr_end, right_ltr_start, right_ltr_end) = ltr_candidates[candidate_index]
    #
    #     ref_seq = ref_contigs[chr_name]
    #     left_ltr_name = chr_name + ':' + str(left_ltr_start) + '-' + str(left_ltr_end)
    #     left_ltr_seq = ref_seq[left_ltr_start - 1: left_ltr_end]
    #
    #     int_ltr_name = chr_name + ':' + str(left_ltr_end) + '-' + str(right_ltr_start)
    #     int_ltr_seq = ref_seq[left_ltr_end: right_ltr_start - 1]
    #
    #     print(left_ltr_name, left_ltr_seq)


    # tmp_blast_dir = '/home/hukang/NeuralLTR/demo/test/test_blastn'
    # threads = 40
    # TRsearch_dir = '/home/hukang/NeuralLTR/tools'
    # query_path = '/home/hukang/NeuralLTR/demo/test/test.fa'
    # subject_path = '/home/hukang/NeuralLTR/demo/test/genome.rename.fa'
    #
    # work_dir = '/home/hukang/NeuralLTR/demo/test'
    # std_tmp_blast_dir = work_dir + '/std_blastn'
    # standard_lib_out = work_dir + '/std_full_length.out'
    #
    # test_tmp_blast_dir = work_dir + '/test_blastn'
    # test_lib_out = work_dir + '/test_full_length.out'
    #
    # # Step 0. Obtaining the length of the genome
    # names, contigs = read_fasta(subject_path)
    # chrom_length = {}
    # for i, name in enumerate(names):
    #     chr_len = len(contigs[name])
    #     chrom_length[name] = chr_len
    #
    # coverage_threshold = 0.8
    # category = 'LTR'
    # tools_dir = TRsearch_dir
    # multi_process_align_v2(query_path, subject_path, standard_lib_out, std_tmp_blast_dir, threads, chrom_length,
    #                        coverage_threshold, category, tools_dir, is_removed_dir=True)

    # intact_dir = tmp_blast_dir + '/intact_tmp'
    # divergence_threshold = 20
    # full_length_threshold = 0.8
    # search_struct = False
    # full_length_annotations, copies_direct = get_full_length_copies_RM(query_path, subject_path, intact_dir, threads,
    #                                                                    divergence_threshold,
    #                                                                    full_length_threshold, search_struct,
    #                                                                    TRsearch_dir)
    # print(full_length_annotations)


    # cur_matrix_file = '/home/hukang/LTR_Benchmarking/LTR_libraries/Ours/zebrafish/test/BHIKHARI-5-LTR_DR.matrix'
    # # sliding_window_size = 10
    # # seq_name, is_ltr = judge_both_ends_frame(cur_matrix_file, sliding_window_size, debug=1)
    # # print(seq_name, is_ltr)
    #
    # # 是否能通过计算随机簇占整个序列的比例来过滤
    # left_lines = []
    # both_lines = []
    # with open(cur_matrix_file, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         parts = line.split('\t')
    #         left_line = parts[0]
    #         right_line = parts[1]
    #         left_lines.append(left_line)
    #         both_lines.append((left_line, right_line))
    #
    # # 我们直接对序列调用cd-hit聚类，然后统计单条序列的簇
    # cur_seq_file = '/home/hukang/LTR_Benchmarking/LTR_libraries/Ours/zebrafish/test/BHIKHARI-5-LTR_DR.fa'
    # cur_seq_cons = '/home/hukang/LTR_Benchmarking/LTR_libraries/Ours/zebrafish/test/BHIKHARI-5-LTR_DR.fa.cons'
    # contigs = {}
    # id = 0
    # for line in left_lines:
    #     name = 'seq_'+str(id)
    #     id += 1
    #     contigs[name] = line.replace('-', '')
    # store_fasta(contigs, cur_seq_file)
    # # cd_hit_command = 'cd-hit-est -aS ' + str(0.95) + ' -aL ' + str(0.95) + ' -c ' + str(0.8) \
    # #                  + ' -G 0 -g 1 -A 80 -i ' + cur_seq_file + ' -o ' + cur_seq_cons + ' -T 0 -M 0'
    # cd_hit_command = 'cd-hit-est -aS ' + str(0.5) + ' -c ' + str(0.8) \
    #                  + ' -G 0 -g 1 -A 10 -i ' + cur_seq_file + ' -o ' + cur_seq_cons + ' -T 0 -M 0'
    # os.system(cd_hit_command + ' > /dev/null 2>&1')
    # # 解析聚类文件中的单拷贝
    # cluster_file = cur_seq_cons + '.clstr'
    # cluster_idx = -1
    # clusters = {}
    # with open(cluster_file, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         if line.startswith('>'):
    #             cluster_idx = line.split(' ')[1]
    #         else:
    #             if not clusters.__contains__(cluster_idx):
    #                 clusters[cluster_idx] = []
    #             cur_cluster = clusters[cluster_idx]
    #             name = line.split(',')[1].split(' ')[1].strip()[1:]
    #             name = name[0: len(name) - 3]
    #             cur_cluster.append(name)
    # print(clusters)
    # homo_lines = []
    # random_lines = []
    # for cluster_idx in clusters.keys():
    #     cur_cluster = clusters[cluster_idx]
    #     if len(cur_cluster) == 1:
    #         seq_name = cur_cluster[0]
    #         seq_id = int(str(seq_name).replace('seq_', ''))
    #         random_line = str(both_lines[seq_id]) + '\t' + str(seq_id)
    #         random_lines.append(random_line)
    #     else:
    #         for seq_name in cur_cluster:
    #             seq_id = int(str(seq_name).replace('seq_', ''))
    #             homo_line = str(both_lines[seq_id]) + '\t' + str(seq_id)
    #             homo_lines.append(homo_line)
    # for homo_line in homo_lines:
    #     print(homo_line)
    # print('\n\n')
    # for random_line in random_lines:
    #     print(random_line)
    # # expand_left_lines = expand_frame_v1(left_lines, both_lines, expand_num=len(left_lines))
    # hamming_distance_threshold = 10
    # cluster_indexes = cluster_sequences(left_lines, hamming_distance_threshold, min_samples_threshold=1)
    # print(cluster_indexes)
    # # 将单条序列的簇放在一起
    # for cluster_index in cluster_indexes.keys():
    #     seq_ids = cluster_indexes[cluster_index]
    #     if len(seq_ids) == 1:
    #         seq_id = seq_ids[0]
    #         print(left_lines[seq_id] + '\t' + str(seq_id))

    # output_dir = '/home/hukang/left_LTR_real_dataset/raw_data/Danio_rerio/positive'
    # high_copy_output_dir = '/home/hukang/left_LTR_real_dataset/raw_data/Danio_rerio/positive_high'
    # get_high_copy_LTR(output_dir, high_copy_output_dir, copy_num_threshold=5)



    # tmp_output_dir = '/home/hukang/NeuralLTR/demo/test2'
    # # Step4.1 调用同源规则，继续过滤深度学习未能识别的假阳性序列
    # dl_output_path = tmp_output_dir + '/is_LTR_deep.txt'
    # homo_output_path = tmp_output_dir + '/is_LTR_homo.txt'
    # output_dir = tmp_output_dir + '/ltr_both_ends_frames'
    # threads = 40
    # filter_ltr_by_homo(dl_output_path, homo_output_path, output_dir, threads)

    # repbase_path = '/home/hukang/NeuralTE_dataset/Dataset1/all_repbase.ref'
    # LTR_path = '/home/hukang/NeuralTE_dataset/Dataset1/all_ltr.ref'
    # rep_names, rep_contigs = read_fasta_v1(repbase_path)
    # labels = set()
    # ltr_labels = ('Copia', 'Gypsy', 'Bel-Pao', 'Retrovirus')
    # ltr_contigs = {}
    # for name in rep_names:
    #     label = name.split('\t')[1]
    #     if label in ltr_labels:
    #         ltr_contigs[name] = rep_contigs[name]
    #     labels.add(label)
    # print(labels)
    # store_fasta(ltr_contigs, LTR_path)


    # work_dir = '/home/hukang/NeuralTE_dataset/raw_Repbase'
    #
    # # 抽取出对应物种的library
    # all_repbase = work_dir + '/all_repbase.ref'
    # cur_species = 'Homo sapiens'
    # species_lib = work_dir + '/human.lib'
    # contignames, contigs = read_fasta_v1(all_repbase)
    # new_contigNames = []
    # new_contigs = {}
    # for name in contignames:
    #     parts = name.split('\t')
    #     if len(parts) != 3:
    #         continue
    #     label = parts[1]
    #     species = parts[2]
    #     if species == cur_species:
    #         new_contigNames.append(name)
    #         new_contigs[name] = contigs[name]
    #
    # wicker2RM = {}
    # # 转换成RepeatMasker标签
    # with open(config.project_dir + '/data/Wicker2RM.info', 'r') as f_r:
    #     for line in f_r:
    #         if line.startswith('#'):
    #             continue
    #         line = line.replace('\n', '')
    #         parts = line.split('\t')
    #         Wicker_Label = parts[0]
    #         RepeatMasker_Label = parts[1]
    #         wicker2RM[Wicker_Label] = RepeatMasker_Label
    #
    # rm_contigs = {}
    # label_set = set()
    # for name in new_contigNames:
    #     parts = name.split('\t')
    #     seq_name = parts[0]
    #     label = parts[1]
    #     label_set.add(label)
    #     RepeatMasker_Label = wicker2RM[label]
    #     new_name = seq_name+'#'+RepeatMasker_Label
    #     rm_contigs[new_name] = new_contigs[name]
    # store_fasta(rm_contigs, species_lib)
    # print(rm_contigs)



    # # 过滤掉.domain文件中的test对应物种的所有记录
    # test_file = '/home/hukang/NeuralTE_experiment_bak2/Dataset3/test.ref'
    # test_names, test_contigs = read_fasta_v1(test_file)
    # test_species_set = set()
    # for name in test_names:
    #     species = name.split('\t')[2]
    #     test_species_set.add(species)
    #
    # protein_db_path = '/home/hukang/NeuralTE/data/RepeatPeps.lib'
    # # extract protein name and species
    # protein_species_dict = {}
    # p_names, p_contigs = read_fasta_v1(protein_db_path)
    # for name in p_names:
    #     protein_name = name.split(' ')[0]
    #     pattern = r'\[(.*?)\]'  # 匹配方括号内的任意字符，非贪婪模式
    #     match = re.search(pattern, name)
    #     if match:
    #         species = match.group(1)  # 获取匹配到的第一个子组
    #     else:
    #         species = 'Unknown'
    #     protein_species_dict[protein_name] = species
    #
    # domain_file = '/home/hukang/NeuralTE_experiment_bak2/Dataset3/test.ref.domain'
    # keep_lines = []
    # with open(domain_file, 'r') as f_r:
    #     for i, line in enumerate(f_r):
    #         if i <= 1:
    #             keep_lines.append(line)
    #             continue
    #         domain_name = line.split('\t')[1]
    #         domain_name = domain_name.replace(',', '')
    #         species = protein_species_dict[domain_name]
    #         if species not in test_species_set:
    #             keep_lines.append(line)
    # with open(domain_file, 'w') as f_save:
    #     for line in keep_lines:
    #         f_save.write(line)


    # # 1.取出鸟类的数据当作验证集
    # work_dir = '/home/hukang/NeuralTE/data/param_tuning/Dataset4'
    # valid_species = work_dir + '/valid_species.txt'
    # valid_species_list = set()
    # with open(valid_species, 'r') as f_r:
    #     for line in f_r:
    #         species = line.replace('\n', '').strip()
    #         valid_species_list.add(species)
    # train_ref = work_dir + '/train.ref'
    # valid_ref = work_dir + '/valid.ref'
    # train_contigs = {}
    # valid_contigs = {}
    # names, contigs = read_fasta_v1(train_ref)
    # for name in names:
    #     cur_species = name.split('\t')[2]
    #     if cur_species in valid_species_list:
    #         valid_contigs[name] = contigs[name]
    #     else:
    #         train_contigs[name] = contigs[name]
    # store_fasta(train_contigs, train_ref)
    # store_fasta(valid_contigs, valid_ref)

    # # 1. 我们抽取出哺乳动物当作测试集，剩余其他物种当作训练集
    # work_dir = '/home/hukang/NeuralTE_experiment/Dataset8'
    # test_species = work_dir + '/test_species.txt'
    # test_species_list = set()
    # with open(test_species, 'r') as f_r:
    #     for line in f_r:
    #         species = line.replace('\n', '').strip()
    #         test_species_list.add(species)
    # all_repbase = work_dir + '/all_repbase.ref'
    # train_ref = work_dir + '/train.ref'
    # test_ref = work_dir + '/test.ref'
    # train_contigs = {}
    # test_contigs = {}
    # names, contigs = read_fasta_v1(all_repbase)
    # for name in names:
    #     cur_species = name.split('\t')[2]
    #     if cur_species in test_species_list:
    #         test_contigs[name] = contigs[name]
    #     else:
    #         train_contigs[name] = contigs[name]
    # store_fasta(train_contigs, train_ref)
    # store_fasta(test_contigs, test_ref)

    # # 1. 我们抽取出开花动物当作测试集，剩余其他物种当作训练集
    # work_dir = '/home/hukang/NeuralTE_experiment/Dataset9'
    # test_species = work_dir + '/flowering_plants.txt'
    # test_species_list = set()
    # with open(test_species, 'r') as f_r:
    #     for line in f_r:
    #         species = line.replace('\n', '').strip()
    #         test_species_list.add(species)
    # all_repbase = work_dir + '/all_repbase.ref'
    # train_ref = work_dir + '/train.ref'
    # test_ref = work_dir + '/test.ref'
    # train_contigs = {}
    # test_contigs = {}
    # names, contigs = read_fasta_v1(all_repbase)
    # for name in names:
    #     cur_species = name.split('\t')[2]
    #     if cur_species in test_species_list:
    #         test_contigs[name] = contigs[name]
    #     else:
    #         train_contigs[name] = contigs[name]
    # store_fasta(train_contigs, train_ref)
    # store_fasta(test_contigs, test_ref)

    # fa = '/home/hukang/NeuralTE_experiment/Dataset7/mouse.ref'
    # names, contigs = read_fasta_v1(fa)
    # types = {}
    # for name in names:
    #     type = name.split('\t')[1]
    #     if not types.__contains__(type):
    #         types[type] = 0
    #     num = types[type]
    #     types[type] = num + 1
    # print(types)

    # protein_db_path = '/home/hukang/NeuralTE/data/RepeatPeps.lib'
    # # extract protein name and species
    # protein_species_dict = {}
    # p_names, p_contigs = read_fasta_v1(protein_db_path)
    # for name in p_names:
    #     protein_name = name.split(' ')[0]
    #     pattern = r'\[(.*?)\]'  # 匹配方括号内的任意字符，非贪婪模式
    #     match = re.search(pattern, name)
    #     if match:
    #         species = match.group(1)  # 获取匹配到的第一个子组
    #     else:
    #         species = 'Unknown'
    #     protein_species_dict[protein_name] = species
    # print(protein_species_dict['UN-Candystripe1_SB_tp#DNA/PIF-Harbinger'])


    # # 从 CDS 中选择432个CDS序列加入到原始序列中
    # all_repbase = '/home/hukang/NeuralTE_dataset/Dataset1/all_repbase.ref'
    # add_cds = '/home/hukang/NeuralTE_dataset/Dataset1/add_cds.fa'
    # CDS_dir = '/home/hukang/NeuralTE_dataset/Dataset1/CDS'
    # cds_list = ['GCF_000001635.27_GRCm39_cds_from_genomic.fna', 'GCF_000001735.4_TAIR10.1_cds_from_genomic.fna', 'GCF_000002035.6_GRCz11_cds_from_genomic.fna']
    # add_cds_contigs = {}
    # count = 0
    # for cds in cds_list:
    #     cds_file = CDS_dir + '/' + cds
    #     cds_names, cds_contigs = read_fasta_v1(cds_file)
    #     num = 144
    #     for i, name in enumerate(cds_names):
    #         if i < num:
    #             seq = cds_contigs[name]
    #             count += 1
    #             new_name = 'CDS_' + str(count) + '\tUnknown\tOryza sativa\tTSD:Unknown\tTSD_len:16\tLTR:\tTIR:'
    #             add_cds_contigs[new_name] = seq
    # store_fasta(add_cds_contigs, add_cds)

    # # 统计一下数据集中每种superfamily的比例是多少
    # all_repbase = '/home/hukang/NeuralTE_dataset/Dataset1/all_repbase.ref'
    # names, contigs = read_fasta_v1(all_repbase)
    # type_num = {}
    # count = 0
    # for name in names:
    #     superfamily = name.split('\t')[1]
    #     if superfamily != 'Unknown':
    #         count += 1
    #
    #     if not type_num.__contains__(superfamily):
    #         type_num[superfamily] = 0
    #     num = type_num[superfamily]
    #     type_num[superfamily] = num + 1
    # print(count)
    # print(len(type_num))

    # output_file = '/home/hukang/NeuralTE_dataset/all_repbase.ref'
    # names, contigs = read_fasta_v1(output_file)
    # mouse_ref = '/home/hukang/NeuralTE_dataset/mouse.ref'
    # zebrafinch_ref = '/home/hukang/NeuralTE_dataset/zebrafinch.ref'
    # chicken_ref = '/home/hukang/NeuralTE_dataset/chicken.ref'
    # ref_list = {'Mus musculus': 'mouse', 'Taeniopygia guttata': 'zebrafinch', 'Gallus gallus': 'chicken'}
    # mouse_contigs = {}
    # zebrafinch_contigs = {}
    # chicken_contigs = {}
    # for name in names:
    #     parts = name.split('\t')
    #     if len(parts) != 3:
    #         continue
    #     species = parts[2]
    #     if ref_list.__contains__(species):
    #         cur_species = ref_list[species]
    #         if cur_species == 'mouse':
    #             mouse_contigs[name] = contigs[name]
    #         elif cur_species == 'zebrafinch':
    #             zebrafinch_contigs[name] = contigs[name]
    #         elif cur_species == 'chicken':
    #             chicken_contigs[name] = contigs[name]
    # store_fasta(mouse_contigs, mouse_ref)
    # store_fasta(zebrafinch_contigs, zebrafinch_ref)
    # store_fasta(chicken_contigs, chicken_ref)


    # all_repbase_path = '/home/hukang/NeuralTE_dataset/Dataset1/all_repbase.ref'
    # names, contigs = read_fasta_v1(all_repbase_path)
    # species = set()
    # for name in names:
    #     species.add(name.split('\t')[2])
    # print(species)
    # species_txt = '/home/hukang/NeuralTE_dataset/Dataset1/species.txt'
    # with open(species_txt, 'w') as f_save:
    #     for spe in species:
    #         f_save.write(spe+'\n')
    # print('hello')
    # test_species_txt = '/home/hukang/NeuralTE_dataset/Dataset1/test_species.txt'
    # lines = []
    # with open(test_species_txt, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '').strip()
    #         print(line)
    #         lines.append(line)
    # with open(test_species_txt, 'w') as f_save:
    #     for line in lines:
    #         f_save.write(line+'\n')


    # # 重新训练DeepTE模型
    # # 1. 训练LINE模型
    # # 1.1 先提取Dataset2中的train.ref中的LINE元素对应的序列，转换成label,sequence格式
    # work_dir = '/home/hukang/NeuralTE_experiment/Dataset9/DeepTE'
    # train = work_dir + '/train.ref'
    # contigNames, contigs = read_fasta_v1(train)
    # LINE_labels = ['R2', 'RTE', 'Jockey', 'L1', 'I']
    # LINE_path = work_dir + '/ipt_shuffle_LINE_CNN_data.txt'
    # unique_labels = set()
    # with open(LINE_path, 'w') as f_save:
    #     for name in contigNames:
    #         label = name.split('\t')[1]
    #
    #         if label in LINE_labels:
    #             f_save.write(label+','+contigs[name]+'\n')
    #             unique_labels.add(label)
    # print(unique_labels)
    # # 2. 训练SINE模型
    # train = work_dir + '/train.ref'
    # contigNames, contigs = read_fasta_v1(train)
    # SINE_labels = ['tRNA', '7SL', '5S']
    # SINE_path = work_dir + '/ipt_shuffle_SINE_CNN_data.txt'
    # unique_labels = set()
    # with open(SINE_path, 'w') as f_save:
    #     for name in contigNames:
    #         label = name.split('\t')[1]
    #         if label in SINE_labels:
    #             f_save.write(label+','+contigs[name]+'\n')
    #             unique_labels.add(label)
    # print(unique_labels)
    # # 3. 训练LTR模型
    # train = work_dir + '/train.ref'
    # contigNames, contigs = read_fasta_v1(train)
    # LTR_labels = ['Copia', 'Gypsy', 'Bel-Pao', 'Retrovirus']
    # LTR_path = work_dir + '/ipt_shuffle_LTR_CNN_data.txt'
    # unique_labels = set()
    # with open(LTR_path, 'w') as f_save:
    #     for name in contigNames:
    #         label = name.split('\t')[1]
    #         if label in LTR_labels:
    #             f_save.write(label+','+contigs[name]+'\n')
    #             unique_labels.add(label)
    # print(unique_labels)
    # # # 4. 训练nLTR模型
    # train = work_dir + '/train.ref'
    # contigNames, contigs = read_fasta_v1(train)
    # nLTR_labels = {'DIRS': 'DIRS', 'Ngaro': 'DIRS', 'VIPER': 'DIRS', 'Penelope': 'PLE', 'R2': 'LINE', 'RTE': 'LINE', 'Jockey': 'LINE', 'L1': 'LINE', 'I': 'LINE', 'tRNA': 'SINE', '7SL': 'SINE', '5S': 'SINE'}
    # nLTR_path = work_dir + '/ipt_shuffle_nLTR_CNN_data.txt'
    # unique_labels = set()
    # with open(nLTR_path, 'w') as f_save:
    #     for name in contigNames:
    #         label = name.split('\t')[1]
    #         if nLTR_labels.__contains__(label):
    #             f_save.write(nLTR_labels[label] + ',' + contigs[name] + '\n')
    #             unique_labels.add(nLTR_labels[label])
    # print(unique_labels)
    # # 5. 训练ClassII模型
    # train = work_dir + '/train.ref'
    # contigNames, contigs = read_fasta_v1(train)
    # ClassII_labels = ['Tc1-Mariner', 'hAT', 'Mutator', 'Merlin', 'Transib', 'P', 'PiggyBac', 'PIF-Harbinger', 'CACTA', 'Crypton']
    # ClassII_path = work_dir + '/ipt_shuffle_ClassII_CNN_data.txt'
    # unique_labels = set()
    # with open(ClassII_path, 'w') as f_save:
    #     for name in contigNames:
    #         label = name.split('\t')[1]
    #         if label in ClassII_labels:
    #             f_save.write(label+','+contigs[name]+'\n')
    #             unique_labels.add(label)
    # print(unique_labels)
    # # 6. 训练ClassI模型
    # train = work_dir + '/train.ref'
    # contigNames, contigs = read_fasta_v1(train)
    # ClassI_labels = {'Copia': 'LTR', 'Gypsy': 'LTR', 'Bel-Pao': 'LTR', 'Retrovirus': 'LTR', 'DIRS': 'nLTR', 'Ngaro': 'nLTR', 'VIPER': 'nLTR', 'Penelope': 'nLTR', 'R2': 'nLTR', 'RTE': 'nLTR', 'Jockey': 'nLTR', 'L1': 'nLTR', 'I': 'nLTR', 'tRNA': 'nLTR', '7SL': 'nLTR', '5S': 'nLTR'}
    # ClassI_path = work_dir + '/ipt_shuffle_ClassI_CNN_data.txt'
    # unique_labels = set()
    # with open(ClassI_path, 'w') as f_save:
    #     for name in contigNames:
    #         label = name.split('\t')[1]
    #         if ClassI_labels.__contains__(label):
    #             f_save.write(ClassI_labels[label] + ',' + contigs[name] + '\n')
    #             unique_labels.add(ClassI_labels[label])
    # print(unique_labels)
    # # 7. 训练All模型
    # train = work_dir + '/train.ref'
    # contigNames, contigs = read_fasta_v1(train)
    # All_labels = {'Tc1-Mariner': 'ClassII', 'hAT': 'ClassII', 'Mutator': 'ClassII', 'Merlin': 'ClassII', 'Transib': 'ClassII', 'P': 'ClassII', 'PiggyBac': 'ClassII',
    #                 'PIF-Harbinger': 'ClassII', 'CACTA': 'ClassII', 'Crypton': 'ClassII', 'Helitron': 'ClassIII', 'Maverick': 'ClassIII', 'Copia': 'ClassI',
    #                 'Gypsy': 'ClassI', 'Bel-Pao': 'ClassI', 'Retrovirus': 'ClassI', 'DIRS': 'ClassI', 'Ngaro': 'ClassI', 'VIPER': 'ClassI',
    #                 'Penelope': 'ClassI', 'R2': 'ClassI', 'RTE': 'ClassI', 'Jockey': 'ClassI', 'L1': 'ClassI', 'I': 'ClassI', 'tRNA': 'ClassI', '7SL': 'ClassI', '5S': 'ClassI'}
    # All_path = work_dir + '/ipt_shuffle_All_CNN_data.txt'
    # unique_labels = set()
    # with open(All_path, 'w') as f_save:
    #     for name in contigNames:
    #         label = name.split('\t')[1]
    #         if All_labels.__contains__(label):
    #             f_save.write(All_labels[label] + ',' + contigs[name] + '\n')
    #             unique_labels.add(All_labels[label])
    # print(unique_labels)

    # data_dir = '/home/hukang/NeuralTE_experiment/Dataset9/DeepTE'
    # test_path = data_dir + '/test.ref'
    # predict_path = data_dir + '/results/opt_DeepTE.txt'
    # # filter random sequence
    # names, contigs = read_fasta_v1(test_path)
    # new_names = []
    # new_contigs = {}
    # random_set = set()
    # for name in names:
    #     if name.startswith('Random_'):
    #         random_set.add(contigs[name])
    #     else:
    #         new_names.append(name)
    #         new_contigs[name] = contigs[name]
    #
    # test_path = data_dir + '/test.ref.filter'
    # with open(test_path, 'w') as f_save:
    #     for name in new_names:
    #         seq = new_contigs[name]
    #         f_save.write('>'+name+'\n'+seq+'\n')
    #
    # predict_path_filter = data_dir + '/results/opt_DeepTE.txt.filter'
    # with open(predict_path_filter, 'w') as f_save:
    #     with open(predict_path, 'r') as f_r:
    #         for line in f_r:
    #             if line.split('\t')[0].__contains__('Random_'):
    #                 continue
    #             else:
    #                 f_save.write(line)
    #
    # evaluate_DeepTE(test_path, predict_path_filter)
    # # 将DeepTE的macro avg由19分类，变成13分类
    # indicators = [0.237, 0.2104, 0.1854]
    # for ind in indicators:
    #     new_ind = 24 * ind / 15
    #     print(round(new_ind, 4))

    # # 替换非ATCG字符
    # data_dir = '/home/hukang/NeuralTE_experiment/Dataset9/DeepTE'
    # train_path = data_dir + '/test.ref'
    # train_contignames, train_contigs = read_fasta_v1(train_path)
    # for name in train_contignames:
    #     seq = train_contigs[name]
    #     seq = replace_non_atcg(seq)
    #     train_contigs[name] = seq
    # store_fasta(train_contigs, train_path)

    # # 抽出非自治转座子
    # work_dir = '/home/hukang/NeuralTE_dataset/Dataset4'
    # repbase_path = '/home/hukang/NeuralTE_dataset/Dataset1/all_repbase.ref'
    # extract_non_autonomous(repbase_path, work_dir)

    # # 将Dataset2中的非ATCG字符替换成空
    # files = ['/home/hukang/NeuralTE_dataset/Dataset2/all_repbase.ref', '/home/hukang/NeuralTE_dataset/Dataset2/train.ref', '/home/hukang/NeuralTE_dataset/Dataset2/test.ref']
    # for f in files:
    #     names, contigs = read_fasta_v1(f)
    #     new_contigs = {}
    #     for name in names:
    #         seq = contigs[name]
    #         seq = replace_non_atcg(seq)
    #         new_contigs[name] = seq
    #     store_fasta(new_contigs, f)

    # pred_path = '/home/hukang/NeuralTE_experiment/Dataset9/TEsorter/test.ref.rexdb.cls.lib'
    # test_path = '/home/hukang/NeuralTE_experiment/Dataset9/TEsorter/test.ref'
    # # filter random sequence
    # names, contigs = read_fasta_v1(test_path)
    # new_names = []
    # new_contigs = {}
    # random_set = set()
    # for name in names:
    #     if name.startswith('Random_'):
    #         random_set.add(contigs[name])
    #     else:
    #         new_names.append(name)
    #         new_contigs[name] = contigs[name]
    #
    # test_path = '/home/hukang/NeuralTE_experiment/Dataset9/TEsorter/test.ref.filter'
    # with open(test_path, 'w') as f_save:
    #     for name in new_names:
    #         seq = new_contigs[name]
    #         f_save.write('>'+name+'\n'+seq+'\n')
    #
    # names, contigs = read_fasta_v1(pred_path)
    # new_names = []
    # new_contigs = {}
    # for name in names:
    #     seq = contigs[name]
    #     if seq in random_set:
    #         continue
    #     new_contigs[name] = seq
    #     new_names.append(name)
    # pred_path = '/home/hukang/NeuralTE_experiment/Dataset9/TEsorter/test.ref.rexdb.cls.lib.filter'
    # with open(pred_path, 'w') as f_save:
    #     for name in new_names:
    #         seq = new_contigs[name]
    #         f_save.write('>'+name+'\n'+seq+'\n')
    #
    #
    # evaluate_TEsorter(pred_path, test_path)
    # # 将TEsorter的macro avg由25分类，变成24分类
    # indicators = [0.5949, 0.3293, 0.3358]
    # for ind in indicators:
    #     new_ind = 4 * ind / 2
    #     print(round(new_ind, 4))

    # # 1.2 提取Dataset2中的test.ref中的LINE元素对应的序列
    # train = '/home/hukang/NeuralTE_experiment/Dataset2/DeepTE/test.ref'
    # contigNames, contigs = read_fasta_v1(train)
    # LINE_labels = ['R2', 'RTE', 'Jockey', 'L1', 'I']
    # LINE_path = '/home/hukang/NeuralTE_experiment/Dataset2/DeepTE/test.LINE.ref'
    # LINE_contigs = {}
    # unique_labels = set()
    # for name in contigNames:
    #     label = name.split('\t')[1]
    #     if label in LINE_labels:
    #         LINE_contigs[name] = contigs[name]
    #         unique_labels.add(label)
    # print(unique_labels)
    # store_fasta(LINE_contigs, LINE_path)


    # # Dataset5
    # dataset1 = '/home/hukang/NeuralTE_dataset/Repbase_raw/all_repbase.ref.raw'
    # contigNames, contigs = read_fasta_v1(dataset1)
    # species_arr = set()
    # for name in contigNames:
    #     species_arr.add(name.split('\t')[2])
    # species_arr = list(species_arr)
    # # 打乱数组顺序
    # random.shuffle(species_arr)
    # # 计算划分点
    # split_point = int(len(species_arr) * 0.8)
    # # 划分数组
    # part_1 = species_arr[:split_point]
    # part_2 = species_arr[split_point:]
    # # 输出结果
    # print("第一部分（80%）：", len(part_1))
    # print("第二部分（20%）：", len(part_2))
    # train = '/home/hukang/NeuralTE_dataset/Dataset3/train.ref'
    # test = '/home/hukang/NeuralTE_dataset/Dataset3/test.ref'
    # train_contigs = {}
    # test_contigs = {}
    # for name in contigNames:
    #     species = name.split('\t')[2]
    #     if species in part_1:
    #         train_contigs[name] = contigs[name]
    #     elif species in part_2:
    #         test_contigs[name] = contigs[name]
    # store_fasta(train_contigs, train)
    # store_fasta(test_contigs, test)
    # print(len(train_contigs), len(test_contigs))


    # feature_path = '/home/hukang/TE_Classification/ClassifyTE/data/new_features.csv.train'
    # list_path = '/home/hukang/TE_Classification/ClassifyTE/new_features/list.txt'
    # list_data_dir = '/home/hukang/TE_Classification/ClassifyTE/new_features/kanalyze-2.0.0/input_data'
    # add_ClassifyTE_classification(feature_path, list_path, list_data_dir)

    # predict_path = '/home/hukang/TE_Classification/ClassifyTE/output/predicted_out_new_features_test.csv'
    # evaluate_ClassifyTE(predict_path)

    # # 1. 调用split_train_test.py将训练集划分训练集和验证集
    # work_dir = '/home/hukang/TE_Classification/TERL/Data/DS9'
    # fasta_file = work_dir + '/train.ref'
    # split_command = 'python /home/hukang/NeuralTE/utils/split_train_test.py --data_path ' + fasta_file + ' --out_dir ' + work_dir
    # #os.system(split_command)
    # out_dir = work_dir + '/Train'
    # generate_TERL_dataset(fasta_file, out_dir)
    # fasta_file = work_dir + '/valid.ref'
    # out_dir = work_dir + '/Test'
    # generate_TERL_dataset(fasta_file, out_dir)

    # fasta_file = '/home/hukang/NeuralTE_dataset/Dataset2/all_repbase.ref'
    # generate_ClassifyTE_dataset(fasta_file)

    # # # #获取RepeatClassifier的结果评估
    # # classified_path = '/home/hukang/NeuralTE_experiment/Dataset9/RC/test.ref.classified'
    # # # 过滤包含的负例
    # # names, contigs = read_fasta_v1(classified_path)
    # # new_contigs = {}
    # # for name in names:
    # #     if name.startswith('Random_'):
    # #         continue
    # #     new_contigs[name] = contigs[name]
    # classified_path = '/home/hukang/NeuralTE_experiment/Dataset9/RC/test.ref.classified.filter'
    # # store_fasta(new_contigs, classified_path)
    # RC_name_labels = evaluate_RepeatClassifier(classified_path)
    # # 将RepeatClassifier的macro avg由25分类，变成24分类
    # indicators = [0.4418, 0.2995, 0.3027]
    # for ind in indicators:
    #     new_ind = 24 * ind / 15
    #     print(round(new_ind, 4))

    # # 将NeuralTE的macro avg由19分类，变成13分类
    # indicators = [0.4144, 0.4113, 0.4068]
    # for ind in indicators:
    #     new_ind = 25 * ind / 15
    #     print(round(new_ind, 4))


    # # 获取NeuralTE分错类的序列，看是否能够改进
    # TE_path = '/home/hukang/NeuralTE_dataset/Dataset2/test.ref'
    # contigNames, contigs = read_fasta(TE_path)
    #
    # work_dir = '/home/hukang/NeuralTE/work/dataset2'
    # info_path = work_dir + '/classified.info'
    # wrong_classified_TE = {}
    # wrong_TE = work_dir + '/wrong_TE.fa'
    # with open(info_path, 'r') as f_r:
    #     for line in f_r:
    #         line = line.replace('\n', '')
    #         if line.startswith('#'):
    #             continue
    #         parts = line.split(',')
    #         if parts[1] != parts[2]:
    #             raw_name = parts[0]
    #             wrong_classified_TE[raw_name] = contigs[raw_name]
    #             # print(line)
    # store_fasta(wrong_classified_TE, wrong_TE)


    # # 我们尝试获取小样本的比对结果，获取比对序列占query比例 > 80%，identity > 80% 的query序列，然后将query的label设置为target label
    # query_path = '/home/hukang/NeuralTE/work/dataset2/test.ref'
    # query_names, query_contigs = read_fasta(query_path)
    # target_path = '/home/hukang/NeuralTE/work/dataset2/minority/train.minority.ref'
    # target_names, target_contigs = read_fasta_v1(target_path)
    # target_labels = {}
    # target_len_dict = {}
    # for name in target_names:
    #     parts = name.split('\t')
    #     target_name = parts[0]
    #     label = parts[1]
    #     target_labels[target_name] = label
    #     target_len_dict[target_name] = len(target_contigs[name])
    #
    # test_minority_out = '/home/hukang/NeuralTE/work/dataset2/minority/test.minority.out'
    # RMOut = '/home/hukang/NeuralTE/work/dataset2/test.ref.out'
    # tools_dir = config.project_dir + '/tools'
    # transfer_RMOut2BlastnOut(RMOut, test_minority_out, tools_dir)
    #
    # query_intervals = {}
    # query_records = {}
    # with open(test_minority_out, 'r') as f_r:
    #     for line in f_r:
    #         parts = line.split('\t')
    #         query_name = parts[0]
    #         subject_name = parts[1]
    #         identity = float(parts[2])
    #         query_start = int(parts[6])
    #         query_end = int(parts[7])
    #         subject_start = int(parts[8])
    #         subject_end = int(parts[9])
    #         e_value = float(parts[10])
    #         if subject_start > subject_end:
    #             temp = subject_start
    #             subject_start = subject_end
    #             subject_end = temp
    #         if e_value > 1e-10:
    #             continue
    #         if not query_intervals.__contains__(query_name):
    #             query_intervals[query_name] = {}
    #         target_intervals = query_intervals[query_name]
    #         if not target_intervals.__contains__(subject_name):
    #             target_intervals[subject_name] = []
    #         intervals = target_intervals[subject_name]
    #         intervals.append((query_start, query_end))
    #
    #         if not query_records.__contains__(query_name):
    #             query_records[query_name] = {}
    #         target_records = query_records[query_name]
    #         if not target_records.__contains__(subject_name):
    #             target_records[subject_name] = []
    #         records = target_records[subject_name]
    #         records.append((query_start, query_end, subject_start, subject_end))
    #
    # query_labels = {}
    # for query_name in query_intervals.keys():
    #     target_intervals = query_intervals[query_name]
    #     target_records = query_records[query_name]
    #     for subject_name in target_intervals.keys():
    #         records = target_records[subject_name]
    #         target_label = target_labels[subject_name]
    #         intervals = target_intervals[subject_name]
    #         merge_intervals = merge_overlapping_intervals(intervals)
    #         # 求总共占比长度
    #         sum_len = 0
    #         for interval in merge_intervals:
    #             sum_len += abs(interval[1] - interval[0])
    #         query_len = len(query_contigs[query_name])
    #         subject_len = target_len_dict[subject_name]
    #         alignment_ratio = float(sum_len) / query_len
    #         if alignment_ratio > 0.8:
    #             if not query_labels.__contains__(query_name):
    #                 query_labels[query_name] = target_label
    #         elif target_label == 'P' or target_label == 'Merlin':
    #             # 如果target是DNA转座子，且query的终端能比对到target的终端，也算
    #             # query的比对是5'-end或者3'-end，同时subject的比对也是5'-end或者3'-end
    #             is_terminal_alignment = False
    #             for record in records:
    #                 if ((record[0] - 1) <= 5 or (query_len - record[1]) <= 5) and ((record[2] - 1) <= 5 or (subject_len - record[3]) <= 5):
    #                     is_terminal_alignment = True
    #                     query_labels[query_name] = target_label
    #                     break
    #
    # print(query_labels)
    # print(len(query_labels))
    #
    # wrong_TE = []
    # raw_NeuralTE_result = '/home/hukang/NeuralTE/work/dataset4/classified.info'
    # y_pred = []
    # y_test = []
    # with open(raw_NeuralTE_result, 'r') as f_r:
    #     for line in f_r:
    #         if line.startswith('#'):
    #             continue
    #         line = line.replace('\n', '')
    #         parts = line.split(',')
    #         seq_name = parts[0]
    #         true_label = parts[1]
    #         pred_label = parts[2]
    #         y_pred.append(pred_label)
    #         y_test.append(true_label)
    #         if true_label != pred_label:
    #             wrong_TE.append(seq_name)
    # y_test = np.array(y_test)
    # y_pred = np.array(y_pred)
    # get_metrics_by_label(y_test, y_pred)
    # print(wrong_TE)
    #
    # # 纠正
    # correct_names = []
    # wrong_TE = []
    # y_pred = []
    # y_test = []
    # with open(raw_NeuralTE_result, 'r') as f_r:
    #     for line in f_r:
    #         if line.startswith('#'):
    #             continue
    #         line = line.replace('\n', '')
    #         parts = line.split(',')
    #         seq_name = parts[0]
    #         true_label = parts[1]
    #         pred_label = parts[2]
    #         if query_labels.__contains__(seq_name):
    #             pred_label = query_labels[seq_name]
    #             correct_names.append(seq_name)
    #         y_pred.append(pred_label)
    #         y_test.append(true_label)
    #         if true_label != pred_label:
    #             wrong_TE.append(seq_name)
    # y_test = np.array(y_test)
    # y_pred = np.array(y_pred)
    # get_metrics_by_label(y_test, y_pred)
    # print(correct_names)
    # print(wrong_TE)



    # # 将训练集的小样本数据先抽出来和训练集与测试集进行同源比对
    # minority_labels = ['Crypton', '5S', '7SL', 'Merlin', 'P', 'R2']
    # work_dir = '/home/hukang/NeuralTE_experiment/Dataset2/NeuralTE'
    # minority_path = work_dir + '/minority.ref'
    # train_path = work_dir + '/train.ref'
    # test_path = work_dir + '/test.ref'
    # minority_contigs = {}
    # train_contigNames, train_contigs = read_fasta_v1(train_path)
    # test_contigNames, test_contigs = read_fasta_v1(test_path)
    # # 1. extract minority dataset
    # for name in train_contigNames:
    #     label = name.split('\t')[1]
    #     if label in minority_labels:
    #         minority_contigs[name] = train_contigs[name]
    # store_fasta(minority_contigs, minority_path)
    # # 2. 进行blastn比对
    # blastn2Results_path = work_dir + '/minority.out'
    # os.system('makeblastdb -in ' + minority_path + ' -dbtype nucl')
    # align_command = 'blastn -db ' + minority_path + ' -num_threads ' \
    #                 + str(40) + ' -query ' + train_path + ' -evalue 1e-20 -outfmt 6 > ' + blastn2Results_path
    # os.system(align_command)

    # 我们现在向Dataset2中插入Dataset1的小样本数据
    # minority_labels = ['Crypton', '5S', '7SL', 'Merlin', 'P', 'R2']
    # dataset1 = '/home/hukang/NeuralTE_dataset/Dataset1/all_repbase.ref'
    # dataset2 = '/home/hukang/NeuralTE_dataset/Dataset2/all_repbase.ref'
    # dataset3 = '/home/hukang/NeuralTE_dataset/Dataset3/all_repbase.ref'
    # DS3_contigs = {}
    # DS3_species = set()
    # DS1_contigNames, DS1_contigs = read_fasta_v1(dataset1)
    # DS2_contigNames, DS2_contigs = read_fasta_v1(dataset2)
    # DS2_raw_names = set()
    # for name in DS2_contigNames:
    #     raw_name = name.split('\t')[0]
    #     DS2_raw_names.add(raw_name)
    # for name in DS1_contigNames:
    #     parts = name.split('\t')
    #     raw_name = parts[0]
    #     label = parts[1]
    #     species = parts[2]
    #     label = config.Repbase_wicker_labels[label]
    #     if raw_name not in DS2_raw_names and label in minority_labels:
    #         new_name = parts[0] + '\t' + label + '\t' + parts[2]
    #         DS3_contigs[new_name] = DS1_contigs[name]
    #         if label == 'Merlin' or label == 'P':
    #             DS3_species.add(species)
    # store_fasta(DS3_contigs, dataset3)
    # print(DS3_species)
    # print(len(DS3_species))
    # tsd_species = ['Allomyces macrogynus', 'Ciona savignyi', 'Drosophila bifasciata', 'Eulimnadia texana',
    #                'Locusta migratoria', 'Oxytricha trifallax', 'Puccinia triticina', 'Capitella teleta',
    #                'Corbicula fluminea', 'Drosophila bocqueti', 'Folsomia candida', 'Owenia fusiformis',
    #                'Parhyale hawaiensis', 'Rhizopus arrhizus']
    # DS3_contigNames, DS3_contigs = read_fasta_v1(dataset3)
    # new_contigs = {}
    # for name in DS3_contigNames:
    #     parts = name.split('\t')
    #     species = parts[2]
    #     if species not in tsd_species:
    #         new_name = parts[0] + '\t' + parts[1] + '\t' + parts[2] + '\t' + 'TSD:\tTSD_len:0\t' + parts[5] + '\t' + parts[6]
    #         new_contigs[new_name] = DS3_contigs[name]
    #     else:
    #         new_contigs[name] = DS3_contigs[name]
    # store_fasta(new_contigs, dataset3)

    # DS3_contigNames, DS3_contigs = read_fasta_v1(dataset3)
    # species_set = set()
    # for name in DS3_contigNames:
    #     parts = name.split('\t')
    #     species = parts[2]
    #     species_set.add(species)
    # print(len(species_set), len(DS3_contigNames))

    # # 将小样本数据的TSD改为Unknown和16
    # dataset3 = '/home/hukang/NeuralTE_dataset/Dataset3/all_repbase.ref'
    # new_contigs = {}
    # contigNames, contigs = read_fasta_v1(dataset3)
    # for name in contigNames:
    #     parts = name.split('\t')
    #     new_name = parts[0] + '\t' + parts[1]+ '\t' + parts[2]+ '\t' + 'TSD:Unknown' + '\t' + 'TSD_len:16' + '\t' + parts[5] + '\t' + parts[6]
    #     new_contigs[new_name] = contigs[name]
    # store_fasta(new_contigs, dataset3)


    # # 获取Dataset1, Dataset2 和 Dfam library中少数类别的数量
    # minority_labels = ['Crypton', '5S', '7SL', 'Merlin', 'P']
    # dataset1 = '/home/hukang/NeuralTE_dataset/Dataset1/all_repbase.ref'
    # dataset2 = '/home/hukang/NeuralTE_dataset/Dataset2/all_repbase.ref'
    # contigNames, contigs = read_fasta_v1(dataset2)
    # minority_count = {}
    # for name in contigNames:
    #     label = name.split('\t')[1]
    #     # label = config.Repbase_wicker_labels[label]
    #     if label in minority_labels:
    #         if not minority_count.__contains__(label):
    #             minority_count[label] = 0
    #         prev_count = minority_count[label]
    #         minority_count[label] = prev_count + 1
    # print(minority_count)

    # minority_labels_dict = {'DNA/Crypton': 'Crypton', 'SINE/5S': '5S', 'SINE/7SL': '7SL', 'DNA/Merlin': 'Merlin', 'DNA/P': 'P'}
    # dataset3 = '/home/hukang/miniconda3/envs/HiTE/share/RepeatMasker/Libraries/RepeatMasker.lib.bak'
    # contigNames, contigs = read_fasta(dataset3)
    # minority_count = {}
    # for name in contigNames:
    #     label = name.split('#')[1]
    #     #label = minority_labels_dict[label]
    #     if label in minority_labels_dict:
    #         if not minority_count.__contains__(label):
    #             minority_count[label] = 0
    #         prev_count = minority_count[label]
    #         minority_count[label] = prev_count + 1
    # print(minority_count)


    # # 我们把train.ref中的Crypton序列都提取出来，然后和我们分错类的序列进行比对，看看什么情况
    # train = '/home/hukang/NeuralTE_dataset/Dataset2/train.ref'
    # extract_labels = ['Crypton', '5S', '7SL', 'Merlin', 'P']
    # train_contigNames, train_contigs = read_fasta_v1(train)
    # file_path = '/home/hukang/NeuralTE_dataset/Dataset2/minority.train.ref'
    # cur_contigs = {}
    # for extract_label in extract_labels:
    #     for name in train_contigNames:
    #         if name.split('\t')[1] == extract_label:
    #             cur_contigs[name] = train_contigs[name]
    # store_fasta(cur_contigs, file_path)

    # # 如果过滤掉了不具备LTR的LTR和不具备TIR的TIR序列，看剩下多少。
    # all_path = '/home/hukang/NeuralTE_dataset/Dataset2/all_repbase.ref'
    # ltr_labels = ('Copia', 'Gypsy', 'Bel-Pao', 'Retrovirus')
    # tir_labels = ('Tc1-Mariner', 'hAT', 'Mutator', 'Merlin', 'Transib', 'P', 'PiggyBac', 'PIF-Harbinger', 'CACTA')
    #
    # contigNames, contigs = read_fasta_v1(all_path)
    # print(len(contigs))
    # total_ltr_num = 0
    # delete_ltr_num = 0
    # total_tir_num = 0
    # delete_tir_num = 0
    # for name in contigNames:
    #     parts = name.split('\t')
    #     label = parts[1]
    #     ltr = parts[5]
    #     tir = parts[6]
    #     if label in ltr_labels:
    #         total_ltr_num += 1
    #         if ltr.split(':')[1] == '':
    #             delete_ltr_num += 1
    #             del contigs[name]
    #     elif label in tir_labels:
    #         total_tir_num += 1
    #         if tir.split(':')[1] == '':
    #             delete_tir_num += 1
    #             del contigs[name]
    # print(len(contigs))
    # print(total_ltr_num, delete_ltr_num)
    # print(total_tir_num, delete_tir_num)

    # # 获取水稻Repbase的拷贝，然后我们人工检查一些获得错误TSD的序列，看能否有办法获得对的TSD
    # repbase_path = '/home/hukang/NeuralTE_dataset/Dataset7/test.ref'
    # genome_path = '/home/hukang/Genome/GCF_001433935.1_IRGSP-1.0_genomic.fna'
    # flanking_len = 20
    # temp_dir = '/home/hukang/NeuralTE/work/TSD/rice_tsd'
    # threads = 40
    # batch_member_files = get_flanking_copies(repbase_path, genome_path, flanking_len, temp_dir, threads)
    #
    # species = 'rice'
    # plant = 1
    # is_expanded = 0
    # label_names, label_contigs = read_fasta_v1(repbase_path)
    # # store repbase name and label
    # repbase_labels = {}
    # for name in label_names:
    #     parts = name.split('\t')
    #     repbase_name = parts[0]
    #     classification = parts[1]
    #     species_name = parts[2]
    #     repbase_labels[repbase_name] = (classification, species_name)
    #
    # # cur_member_file = temp_dir + '/MuDR-N208E_OS.blast.bed.fa'
    # # batch_member_files = [('MuDR-N208E_OS', 'GGGTTAATTTGATCCATGCCACTGCAAATTTAGCTATTCAGAAAAATGACATTGCAATTCATCTATTCTTAACCGTGCCACTGAAATTTTGTAAAACTAAAACCGTGCCATTGACGTCACATTTTCCATCCATTCTCTTCCTTTTCCGTCTTCTTCCTTCCTTCTCCCATCTTCTTCCCGGAGTCAAGCCGGAGAGGGAGCTCGCCGGCAAGGTGAACGAACCCAACCTCGAGTGCGGTTGGCGTGGTCGGCGAATCCGGCGGTGGCGGCGTCGGACAATGGTGGCATCGGGACTCGGGCGGAACCAGCTGAGGCCTAGGCTGGGTGTCGAGCGTGATCGACGACGGTGACTCTCTTCTTCCGCGTTGCTGCTCAACCTCGGCTCCCGCTCTGGCCTCCGGGTCGGTGAGCACCTCATGCCGGCCGCTCTCCCTCGCGGCAGTGCTCTCCCCGCACTACTCATCCTTGGACCTCTCCGAGACTCCAACCGCCTCCTCGTCGCCCGCCATGAGCTCCGCAAGTAGCTGGAGCACCTCGCCGCCGTCTTCGAGTCTTCATCGCCTCTGTTGGGTCTGCTCGCGCCATGCAGCGCCAGATCTACCGTTGTTTCCATCGACGTTGGCTCCACCGCCGACACCGTCGAAGCTCGCTGCGGCCATGGATGGATGGACGACCGCCGTTGGCATCGCCGCCGCTGCTCCCGCACGAGATCTCGCCCACTTGGCTCCGGGAAGAAATGGGAGAAAGAAGGAGCCGCTGCCGCCTGCCATGGTTGGGTCGCTGGCAAGCTCCCTCTCCTCCGAGCTCGCCGGCATGACTCCGAGAAGAAATAGGAGAAAGAACCGGGAAGAAATGGGAGAAAGAAGGAAGAAGACGGAAAAGGCAGGGGATGGATGGAAAACGTGACGGCAGTGGCACGGTTCTAATTTTGTAAAATTCTAGTGGCACGGTTACGAATAGACGAATTGTAATGGCATTTTTCTTAATAGACAAATTTGCAGTGGCATAGATCAAATTAACCCTA', cur_member_file)]
    #
    # tsd_info = get_copies_TSD_info(batch_member_files, flanking_len, is_expanded, repbase_labels, threads)
    # # 将所有的序列存成文件，拷贝序列名称为在原有的名称后面加-C_{num}
    # names, contigs = read_fasta(repbase_path)
    # final_repbase_path = temp_dir + '/' + species + '.ref'
    # final_repbase_contigs = {}
    # for query_name in names:
    #     seq = contigs[query_name]
    #     label_item = repbase_labels[query_name]
    #
    #     if tsd_info.__contains__(query_name):
    #         copies_tsd_info = tsd_info[query_name]
    #     else:
    #         copies_tsd_info = [('Unknown', 16, -1)]
    #
    #     # 遍历一边，取出distance为0的TSD；
    #     new_copies_tsd_info = []
    #     for tsd_seq, tsd_len, cur_distance in copies_tsd_info:
    #         if cur_distance <= 0:
    #             new_copies_tsd_info.append((tsd_seq, tsd_len, cur_distance))
    #     copies_tsd_info = new_copies_tsd_info if len(new_copies_tsd_info) > 0 else copies_tsd_info
    #     # 将所有的拷贝对应的TSD存储起来，记录每种长度的TSD对应的出现次数和离原始边界的距离
    #     max_count_TSD = {}
    #     length_count = {}
    #     for tsd_seq, tsd_len, cur_distance in copies_tsd_info:
    #         if not length_count.__contains__(tsd_len):
    #             length_count[tsd_len] = (1, cur_distance)
    #             max_count_TSD[tsd_len] = tsd_seq
    #         else:
    #             prev_count, prev_distance = length_count[tsd_len]
    #             if cur_distance < prev_distance:
    #                 prev_distance = cur_distance
    #                 max_count_TSD[tsd_len] = tsd_seq
    #             length_count[tsd_len] = (prev_count + 1, prev_distance)
    #     # 按照(tsd_len, tsd_seq, 出现次数, 最小距离)存成数组
    #     # 取出现次数最多的TSD，如果有多个出现次数最多，取distance最小的那个，如果distance相同，取最长的那个
    #     all_tsd_set = []
    #     for tsd_len in length_count.keys():
    #         cur_count, cur_distance = length_count[tsd_len]
    #         tsd_seq = max_count_TSD[tsd_len]
    #         all_tsd_set.append((tsd_len, tsd_seq, cur_count, cur_distance))
    #     all_tsd_set = sorted(all_tsd_set, key=lambda x: (-x[2], x[3], -x[0]))
    #     final_tsd_info = all_tsd_set[0]
    #     tsd_seq = final_tsd_info[1]
    #     tsd_len = final_tsd_info[0]
    #     tsd_distance = final_tsd_info[3]
    #     if tsd_distance > 5:
    #         tsd_seq = ''
    #         tsd_len = len(tsd_seq)
    #     new_name = query_name + '\t' + label_item[0] + '\t' + label_item[1] + '\t' + 'TSD:' + str(tsd_seq) + '\t' + 'TSD_len:' + str(tsd_len)
    #     final_repbase_contigs[new_name] = seq
    # store_fasta(final_repbase_contigs, final_repbase_path)


    # # 将RepeatMasker的同源搜索库替换成train.ref，看是否能运行
    # work_dir = '/home/hukang/miniconda3/envs/HiTE/share/RepeatMasker/Libraries'
    # rm_lib = work_dir + '/RepeatMasker.lib'
    # train = '/home/hukang/NeuralTE_experiment/Dataset9/NeuralTE/train.ref'
    # test = '/home/hukang/NeuralTE_experiment/Dataset9/NeuralTE/test.ref'
    # wicker2RM = {}
    # # 转换成RepeatMasker标签
    # with open(config.project_dir + '/data/Wicker2RM.info', 'r') as f_r:
    #     for line in f_r:
    #         if line.startswith('#'):
    #             continue
    #         line = line.replace('\n', '')
    #         parts = line.split('\t')
    #         Wicker_Label = parts[0]
    #         RepeatMasker_Label = parts[1]
    #         wicker2RM[Wicker_Label] = RepeatMasker_Label
    #
    # train_contigNames, train_contigs = read_fasta_v1(train)
    # rm_contigs = {}
    # for name in train_contigNames:
    #     parts = name.split('\t')
    #     seq_name = parts[0]
    #     label = parts[1]
    #     species = parts[2].replace(' ', '_')
    #     RepeatMasker_Label = wicker2RM[label]
    #     new_name = seq_name+'#'+RepeatMasker_Label+' @'+species
    #     rm_contigs[new_name] = train_contigs[name]
    # store_fasta(rm_contigs, rm_lib)
    #
    # test_contigNames, test_contigs = read_fasta_v1(test)
    # test_species_set = set()
    # for name in test_contigNames:
    #     parts = name.split('\t')
    #     species = parts[2]
    #     test_species_set.add(species)
    #
    # # 将RepeatMasker的蛋白质库中的test物种过滤掉
    # rm_pep_lib = work_dir + '/RepeatPeps.lib.bak'
    # train_contigNames, train_contigs = read_fasta_v1(rm_pep_lib)
    # filter_pep_lib = work_dir + '/RepeatPeps.lib'
    # rm_contigs = {}
    # for name in train_contigNames:
    #     pattern = r'\[(.*?)\]'  # 匹配方括号内的任意字符，非贪婪模式
    #     match = re.search(pattern, name)
    #     if match:
    #         species = match.group(1)  # 获取匹配到的第一个子组
    #     else:
    #         species = 'Unknown'
    #     # filter content in '()'
    #     pattern = r'\([^)]*\)'
    #     species = re.sub(pattern, '', species)
    #     species = re.sub(r'\s+', ' ', species).strip()
    #     if species not in test_species_set:
    #         rm_contigs[name] = train_contigs[name]
    # store_fasta(rm_contigs, filter_pep_lib)

    # # 将RepeatMasker的蛋白质库中的特定物种过滤掉
    # rm_pep_lib = work_dir + '/RepeatPeps.lib.bak'
    # filter_species = 'Mus musculus'
    # train_contigNames, train_contigs = read_fasta_v1(rm_pep_lib)
    # filter_pep_lib = work_dir + '/RepeatPeps.lib'
    # rm_contigs = {}
    # for name in train_contigNames:
    #     if not name.__contains__(filter_species):
    #         rm_contigs[name] = train_contigs[name]
    # store_fasta(rm_contigs, filter_pep_lib)

    # # 将RepeatMasker的同源搜索Library去掉test数据集，然后测下性能。
    # work_dir = '/home/hukang/miniconda3/envs/HiTE/share/RepeatMasker/Libraries'
    # rm_lib = work_dir + '/RepeatMasker.lib.bak'
    # #test_lib = '/home/hukang/NeuralTE_experiment/Dataset2/RepeatClassifier_all-test_lib/test.ref'
    # test_lib = '/home/hukang/NeuralTE_experiment/novel_TE/rice/oryrep.ref'
    # test_contigNames, test_contigs = read_fasta(test_lib)
    # print(len(test_contigs))
    # all_test_names = set()
    # test_name_dict = {}
    # for name in test_contigNames:
    #     if name.__contains__('intactLTR'):
    #         internal_name = name.replace('intactLTR', 'I')
    #         internal_name1 = name.replace('intactLTR', 'INT')
    #         internal_name2 = name.replace('intactLTR', 'int')
    #         LTR_name = name.replace('intactLTR', 'LTR')
    #
    #         all_test_names.add(LTR_name)
    #         all_test_names.add(internal_name)
    #         all_test_names.add(internal_name1)
    #         all_test_names.add(internal_name2)
    #         test_name_dict[LTR_name] = name
    #         test_name_dict[internal_name] = name
    #         test_name_dict[internal_name1] = name
    #         test_name_dict[internal_name2] = name
    #     else:
    #         all_test_names.add(name)
    #         test_name_dict[name] = name
    #
    # filter_rm_lib = work_dir + '/RepeatMasker.lib'
    # filter_contigs = {}
    # rm_contigNames, rm_contigs = read_fasta_v1(rm_lib)
    # for name in rm_contigNames:
    #     seq_name = name.split('\t')[0].split('#')[0]
    #     if not seq_name in all_test_names:
    #         filter_contigs[name] = rm_contigs[name]
    #     else:
    #         test_name = test_name_dict[seq_name]
    #         if test_contigs.__contains__(test_name):
    #             del test_contigs[test_name]
    # store_fasta(filter_contigs, filter_rm_lib)
    #
    # # 输出Repbase有，但是RM lib中没有的序列
    # print(len(test_contigs))
    # print(test_contigs.keys())

    # # Dataset6制作
    # work_dir = '/home/hukang/NeuralTE_dataset/Dataset7'
    # total = work_dir + '/all_repbase.ref'
    # train = work_dir + '/train.ref'
    # test = work_dir + '/test.ref'
    # names, contigs = read_fasta_v1(total)
    # train_contigs = {}
    # test_contigs = {}
    # for name in names:
    #     species = name.split('\t')[2]
    #     if species == 'Zea mays':
    #         test_contigs[name] = contigs[name]
    #     elif not species.__contains__('Zea mays'):
    #         train_contigs[name] = contigs[name]
    # store_fasta(train_contigs, train)
    # store_fasta(test_contigs, test)

    # work_dir = '/home/hukang/TE_Classification/TERL'
    # test_path = work_dir + '/Data/DS9/test.ref'
    # predict_path = work_dir + '/TERL_20240326_152422_test.ref'
    # # filter random sequence
    # names, contigs = read_fasta_v1(test_path)
    # new_names = []
    # new_contigs = {}
    # random_set = set()
    # for name in names:
    #     if name.startswith('Random_'):
    #         random_set.add(contigs[name])
    #     else:
    #         new_names.append(name)
    #         new_contigs[name] = contigs[name]
    #
    # test_path = work_dir + '/Data/DS9/test.ref.filter'
    # with open(test_path, 'w') as f_save:
    #     for name in new_names:
    #         seq = new_contigs[name]
    #         f_save.write('>'+name+'\n'+seq+'\n')
    #
    # names, contigs = read_fasta_v1(predict_path)
    # new_names = []
    # new_contigs = {}
    # for name in names:
    #     seq = contigs[name]
    #     if seq in random_set:
    #         continue
    #     new_contigs[name] = seq
    #     new_names.append(name)
    # predict_path = work_dir + '/TERL_20240326_152422_test.ref.filter'
    # with open(predict_path, 'w') as f_save:
    #     for name in new_names:
    #         seq = new_contigs[name]
    #         f_save.write('>'+name+'\n'+seq+'\n')
    #
    # evaluate_TERL(test_path, predict_path)
    # # 将TERL的macro avg由19分类，变成13分类
    # indicators = [0.2386, 0.2187, 0.2032]
    # for ind in indicators:
    #     new_ind = 23 * ind / 15
    #     print(round(new_ind, 4))


    # # 获取DeepTE的结果评估
    # DeepTE的评估流程：
    # 1.将数据集转fasta格式，并根据Repbase 28.06 (Dataset1)恢复数据集的header
    # 2.使用NeuralTE提供的split_train_test.py 划分训练train_dataset和测试集test_dataset，保持类别的分布一致。
    # 3.下载DeepTE提供的Metazoans_model。
    # 4.对test_dataset 进行domain的识别。
    # 5.使用训练好的模型对test_dataset 进行预测。
    # data_dir = '/home/hukang/TE_Classification/DeepTE-master/training_example_dir/input_dir'
    # repbase_dataset = '/home/hukang/NeuralTE_dataset/Dataset1/all_repbase.ref'
    # raw_dataset = data_dir + '/ipt_shuffle_All_CNN_data.txt.bak'
    # fasta_dataset = data_dir + '/ipt_shuffle_All_CNN_data.fa'
    # transform_DeepTE_to_fasta(raw_dataset, fasta_dataset)
    # train_path = data_dir + '/train.ref'
    # test_path = data_dir +'/test.ref'
    # # transform DeepTE label to wicker label
    # train_wicker_path = data_dir + '/train.wicker.ref'
    # names, contigs = read_fasta_v1(train_path)
    # with open(train_wicker_path, 'w') as f_save:
    #     for name in names:
    #         parts = name.split('\t')
    #         seq_name = parts[0]
    #         label = parts[1]
    #         wicker_label = config.DeepTE_class[label]
    #         new_name = seq_name + '\t' + wicker_label + '\t' + 'Unknown'
    #         f_save.write('>'+new_name+'\n'+contigs[name]+'\n')
    # test_wicker_path = data_dir + '/test.wicker.ref'
    # names, contigs = read_fasta_v1(test_path)
    # with open(test_wicker_path, 'w') as f_save:
    #     for name in names:
    #         parts = name.split('\t')
    #         seq_name = parts[0]
    #         label = parts[1]
    #         wicker_label = config.DeepTE_class[label]
    #         new_name = seq_name + '\t' + wicker_label + '\t' + 'Unknown'
    #         f_save.write('>' + new_name + '\n' + contigs[name] + '\n')
    #
    # predict_path = data_dir + '/results/opt_DeepTE.txt'
    # evaluate_DeepTE(test_wicker_path, predict_path)



    # # 分析一下repbase中到底包含了多少类别
    # work_dir = '/home/hukang/RepBase28.06.fasta'
    # # 获取指定目录下的所有文件
    # files = get_all_files(work_dir)
    # unique_labels = set()
    # for file in files:
    #     names, contigs = read_fasta_v1(file)
    #     for name in names:
    #         parts = name.split('\t')
    #         if len(parts) == 3:
    #             seq_name = parts[0]
    #             if seq_name.__contains__('LTR') and not seq_name.endswith('-LTR') and not seq_name.endswith('_LTR'):
    #                 unique_labels.add(name)
    # print(unique_labels, len(unique_labels))

    # work_dir = '/home/hukang/HiTE_lib/rice_unmask'
    # test_library = work_dir + '/repbase/classified_TE.fa'
    # gold_library = work_dir + '/oryrep.RM.ref'
    # test_names, test_contigs = read_fasta(test_library)
    # gold_names, gold_contigs = read_fasta(gold_library)
    # gold_labels = {}
    # for name in gold_names:
    #     parts = name.split('#')
    #     seq_name = parts[0]
    #     label = parts[1]
    #     gold_labels[seq_name] = label
    # same_labels = {}
    # diff_labels = {}
    # diff_solo_LTR_labels = {}
    # for name in test_names:
    #     parts = name.split('#')
    #     seq_name = parts[0]
    #     label = parts[1]
    #     gold_label = gold_labels[seq_name]
    #     if label == gold_label:
    #         same_labels[seq_name] = label
    #     else:
    #         diff_labels[seq_name] = (gold_label, label)
    #
    # print(diff_labels)
    # print(len(same_labels), len(diff_labels))




    # # `>Gypsy-171_OS-I`和`>Gypsy-171_OS-LTR`
    # data_path = '/home/hukang/HiTE_lib/rice_unmask/rice-families.rename.fa'
    # # 将输入文件格式化为Repbase格式
    # names, contigs = read_fasta(data_path)
    # os.makedirs(os.path.dirname(data_path), exist_ok=True)
    # with open(data_path, 'w') as f_save:
    #     for name in names:
    #         seq = contigs[name]
    #         name = name.split('/')[0].split('#')[0]
    #         new_name = name + '\tUnknown\tUnknown'
    #         f_save.write('>' + new_name + '\n' + seq + '\n')
    #
    # config.work_dir = '/home/hukang/HiTE_lib/rice_unmask'
    # connect_LTR(data_path)

    # file_path = '/home/hukang/HiTE_lib/rice_unmask/confident_TE.cons.fa'
    # contigNames, contigs = read_fasta(file_path)
    # new_contigs = {}
    # for contigname in contigNames:
    #     seq = contigs[contigname]
    #     if contigname.endswith('_LTR') or contigname.endswith('_INT'):
    #         contigname = contigname.replace('_LTR', '-LTR').replace('_INT', '-INT')
    #     new_contigs[contigname] = seq
    # store_fasta(new_contigs, file_path)


    # work_dir = '/home/hukang/NeuralTE_dataset/Dataset7'
    # data_path = work_dir + '/test.ref'
    # wicker2RM = {}
    # # 转换成RepeatMasker标签
    # with open(config.project_dir + '/data/Wicker2RM.info', 'r') as f_r:
    #     for line in f_r:
    #         if line.startswith('#'):
    #             continue
    #         line = line.replace('\n', '')
    #         parts = line.split('\t')
    #         Wicker_Label = parts[0]
    #         RepeatMasker_Label = parts[1]
    #         wicker2RM[Wicker_Label] = RepeatMasker_Label
    #
    # filter_labels = ['SAT', 'Multicopy gene', 'Satellite', 'REP-10_OS', 'REP-1_OS', 'snRNA', 'RCH2', 'KRISPIE']
    # orig_names, orig_contigs = read_fasta_v1(data_path)
    # classified_data = work_dir + '/test.RM.ref'
    # classified_contigs = {}
    # for name in orig_names:
    #     parts = name.split('\t')
    #     map_name = parts[0]
    #     predict_label = parts[1]
    #     if predict_label in filter_labels:
    #         continue
    #     predict_label = wicker2RM[predict_label]
    #     new_name = map_name + '#' + predict_label
    #     classified_contigs[new_name] = orig_contigs[name]
    # store_fasta(classified_contigs, classified_data)

    # work_dir = '/home/hukang/Genome'
    # genome_path = work_dir + '/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa'
    # ref_names, ref_contigs = read_fasta(genome_path)
    # store_fasta(ref_contigs, genome_path)

    # # 2. 在test数据集上评估RepeatClassifier
    # # 2.2 将Dfam分类名称转成wicker格式
    # # 2.2.1 这个文件里包含了RepeatMasker类别、Repbase、wicker类别的转换
    # rmToWicker = {}
    # WickerToRM = {}
    # wicker_superfamily_set = set()
    # with open(config.project_dir + '/data/TEClasses.tsv', 'r') as f_r:
    #     for i, line in enumerate(f_r):
    #         parts = line.split('\t')
    #         rm_type = parts[5]
    #         rm_subtype = parts[6]
    #         repbase_type = parts[7]
    #         wicker_type = parts[8]
    #         wicker_type_parts = wicker_type.split('/')
    #         #print(rm_type + ',' + rm_subtype + ',' + repbase_type + ',' + wicker_type)
    #         if len(wicker_type_parts) != 3:
    #             continue
    #         wicker_superfamily_parts = wicker_type_parts[-1].strip().split(' ')
    #         if len(wicker_superfamily_parts) == 1:
    #             wicker_superfamily = wicker_superfamily_parts[0]
    #         elif len(wicker_superfamily_parts) > 1:
    #             wicker_superfamily = wicker_superfamily_parts[1].replace('(', '').replace(')', '')
    #         rm_full_type = rm_type + '/' + rm_subtype
    #         if wicker_superfamily == 'ERV':
    #             wicker_superfamily = 'Retrovirus'
    #         if wicker_superfamily == 'Viper':
    #             wicker_superfamily = 'VIPER'
    #         if wicker_superfamily == 'H':
    #             wicker_superfamily = 'Helitron'
    #         rmToWicker[rm_full_type] = (wicker_superfamily, repbase_type)
    #         WickerToRM[wicker_superfamily] = rm_full_type
    #         wicker_superfamily_set.add(wicker_superfamily)
    # # 补充一些元素
    # rmToWicker['LINE/R2'] = 'R2'
    # rmToWicker['DNA/Crypton'] = 'Crypton'
    # rmToWicker['Unknown'] = 'Unknown'
    # # 固定wicker对应的RM标签
    # WickerToRM['Retrovirus'] = 'LTR/ERV'
    # WickerToRM['DIRS'] = 'LTR/DIRS'
    # WickerToRM['R2'] = 'LINE/R2'
    # WickerToRM['RTE'] = 'LINE/RTE-RTE'
    # WickerToRM['L1'] = 'LINE/L1'
    # WickerToRM['I'] = 'LINE/I'
    # WickerToRM['tRNA'] = 'SINE/tRNA'
    # WickerToRM['7SL'] = 'SINE/7SL'
    # WickerToRM['5S'] = 'SINE/5S'
    # WickerToRM['Helitron'] = 'RC/Helitron'
    # WickerToRM['Maverick'] = 'DNA/Maverick'
    # WickerToRM['Crypton'] = 'DNA/Crypton'
    # WickerToRM['Tc1-Mariner'] = 'DNA/TcMar'
    # WickerToRM['hAT'] = 'DNA/hAT'
    # WickerToRM['Mutator'] = 'DNA/MULE'
    # WickerToRM['P'] = 'DNA/P'
    # WickerToRM['PiggyBac'] = 'DNA/PiggyBac'
    # WickerToRM['PIF-Harbinger'] = 'DNA/PIF-Harbinger'
    # WickerToRM['CACTA'] = 'DNA/CMC-EnSpm'
    # print(rmToWicker)
    # #print(len(rmToWicker))
    # #print(wicker_superfamily_set)
    # #print(len(wicker_superfamily_set))
    # ClassifySystem = config.project_dir + '/data/Wicker2RM.info'
    # with open(ClassifySystem, 'w') as f_save:
    #     f_save.write('#Wicker_Label\tRepeatMasker_Label\n')
    #     for name in WickerToRM.keys():
    #         f_save.write(name + '\t' + WickerToRM[name] + '\n')

    # # 画一个3D图
    # work_dir = '/home/hukang/NeuralTE/data/param_tuning/Dataset1'
    # Node_matrix = np.random.rand(49, 3)
    # data_path = work_dir + '/kmer_size_test.xlsx'
    # data_frame = pd.read_excel(data_path)
    # x = data_frame.iloc[:, 0].values
    # y = data_frame.iloc[:, 1].values
    # z = data_frame.iloc[:, 5].values
    # # 将k-mer映射为一个具体值
    # kmer_sizes_dict = {'[2]': 0, '[3]': 1, '[4]': 2, '[5]': 3, '[1, 2]': 4, '[1, 3]': 5, '[1, 4]': 6,
    #                             '[1, 5]': 7, '[2, 3]': 8, '[2, 4]': 9, '[2, 5]': 10, '[3, 4]': 11, '[3, 5]': 12,
    #                             '[1, 2, 3]': 13, '[1, 2, 4]': 14, '[1, 2, 5]': 15, '[1, 3, 4]': 16, '[1, 3, 5]': 17,
    #                             '[2, 3, 4]': 18, '[2, 3, 5]': 19, '[3, 4, 5]': 20}
    # new_x = []
    # new_y = []
    # new_z = []
    # for i in range(len(x)):
    #     new_x.append(kmer_sizes_dict[x[i]])
    #     new_y.append(kmer_sizes_dict[y[i]])
    #     new_z.append(z[i])
    # x = np.array(new_x)
    # y = np.array(new_y)
    # z = np.array(new_z)
    # plot_3D_param(x, y, z, work_dir)

    # # plot sequence logo
    # TE_path = '/home/hukang/NeuralTE_dataset/Dataset2/all_repbase.ref'
    # tmp_output_dir = '/home/hukang/NeuralTE_experiment/Seq_logos'
    # generate_seq_logos(TE_path, tmp_output_dir)

    # 识别HiTE中的新转座子
    # tmp_output_dir = '/home/hukang/NeuralTE_experiment/novel_TE/maize'
    # confident_tir_path = tmp_output_dir + '/maize.fa'
    # tir_repbase_path =  tmp_output_dir + '/maize.ref'
    # species = 'Zea mays'
    # names, contigs = read_fasta(confident_tir_path)
    # new_contigs = {}
    # for name in names:
    #     parts = name.split('#')
    #     raw_name = parts[0]
    #     label = parts[1]
    #     new_name = raw_name+'\t'+label+'\t'+species
    #     new_contigs[new_name] = contigs[name]
    # store_fasta(new_contigs, confident_tir_path)
    # # 连接LTRs
    # connect_LTR(confident_tir_path)
    # connect_LTR(tir_repbase_path)
    # # 去掉solo-LTR
    # names, contigs = read_fasta(confident_tir_path)
    # for name in names:
    #     if name.endswith('-LTR') or name.endswith('-INT'):
    #         del contigs[name]
    # store_fasta(contigs, confident_tir_path)
    # names, contigs = read_fasta(tir_repbase_path)
    # for name in names:
    #     if name.endswith('-LTR') or name.endswith('-INT') or name.endswith('_LTR') or name.endswith('_I') or name.endswith('-I') or name.__contains__('-I_'):
    #         del contigs[name]
    # store_fasta(contigs, tir_repbase_path)

    # identify_new_TIR(confident_tir_path, tir_repbase_path, tmp_output_dir)

    # # #识别每一种superfamily，存成文件
    # novel_tir_consensus = tmp_output_dir + '/novel_tir.fa'
    # contigNames, contigs = read_fasta(novel_tir_consensus)
    # classified_info = tmp_output_dir + '/classified.info'
    # name_label_dict = {}
    # label_nums = {}
    # with open(classified_info, 'r') as f_r:
    #     for line in f_r:
    #         if line.startswith('#'):
    #             continue
    #         parts = line.replace('\n', '').split(',')
    #         raw_name = parts[0]
    #         label = parts[2]
    #         name_label_dict[raw_name] = label
    #         if not label_nums.__contains__(label):
    #             label_nums[label] = 0
    #         num = label_nums[label]
    #         label_nums[label] = num + 1
    # print(label_nums)
    # superfamilies_contigs = {}
    # for name in contigNames:
    #     label = name_label_dict[name]
    #     if not superfamilies_contigs.__contains__(label):
    #         superfamilies_contigs[label] = {}
    #     cur_superfamily_contigs = superfamilies_contigs[label]
    #     new_name = name + '#' + label
    #     cur_superfamily_contigs[new_name] = contigs[name]
    # for label in superfamilies_contigs.keys():
    #     file_path = tmp_output_dir + '/' + label + '.fa'
    #     store_fasta(superfamilies_contigs[label], file_path)
    #
    # # 分析每一种novel TIR的插入时间
    # work_dir = tmp_output_dir
    # tir_path = tmp_output_dir + '/Tc1-Mariner.fa'
    # miu = 1.3e-8
    # type = 'Tc1-Mariner'
    # analyz_TIR_insert_time(tir_path, work_dir, miu, type)
    # TE_list = ['CACTA', 'hAT', 'Mutator', 'PIF-Harbinger', 'Tc1-Mariner']
    # colors = ['#B77072', '#F08B47', '#68A47F', '#B3A4BB', '#B6E1C5']
    # for i in range(len(TE_list)):
    #     type = TE_list[i]
    #     output_path = tmp_output_dir + '/' + type + '_insert_time.txt'
    #     output_fig = tmp_output_dir + '/' + type + '_insert_time.png'
    #     get_insert_time_dist_boxplot(output_path, output_fig, type, colors[i])


    # # 将RepeatMasker的注释转为bed文件，并过滤掉非全长TE
    # RMOut = tmp_output_dir + '/novel_tir.out'
    # out_bed = tmp_output_dir + '/novel_tir.full_length.bed'
    # consensus_path = tmp_output_dir + '/novel_tir.fa'
    # tools_dir = config.project_dir + '/tools'
    # coverage_threshold = 0.95
    # transfer_RMOut2Bed(RMOut, out_bed, consensus_path, tools_dir, coverage_threshold, name_label_dict)
    #
    # chromosomes = [f'Chr{i}' for i in range(1, 13)]
    # # 生成novel TIR的基因组注释
    # # 将水稻的gff文件染色体加上Chr，并且将TE label加到里面
    # novel_tir_gff_bak =  tmp_output_dir + '/novel_tir.gff.bak'
    # novel_tir_gff = tmp_output_dir + '/novel_tir.gff'
    # generate_predict_gff(novel_tir_gff_bak, novel_tir_gff, name_label_dict, chromosomes)
    # # total_genome_len = 374424240
    # # work_dir = tmp_output_dir
    # # analyze_class_ratio_gff(novel_tir_gff, work_dir, total_genome_len)
    #
    # # 存储基因组染色体长度
    # genome_path = '/home/hukang/Genome/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa'
    # names, contigs = read_fasta(genome_path)
    # genome_len_info = tmp_output_dir + '/Genome_len.chr'
    # with open(genome_len_info, 'w') as f_save:
    #     f_save.write('Genome Length\n')
    #     for name in names:
    #         new_name = 'Chr' + str(name)
    #         if new_name not in chromosomes:
    #             continue
    #         chr_len = len(contigs[name])
    #         f_save.write(new_name+'\t'+str(chr_len)+'\n')

    # # t-SNE可视化
    # config.is_predict = 0
    # config.use_minority = 0
    # config.work_dir = '/home/hukang/NeuralTE/work/kmer_size_search'
    # data_path = config.work_dir + "/train.ref"
    # print('data_path:', data_path)
    # # 实例化 DataProcessor 类
    # data_processor = DataProcessor()
    # # 加载数据
    # X, y, seq_names, data_path = data_processor.load_data(config.internal_kmer_sizes, config.terminal_kmer_sizes,
    #                                                       data_path)
    # X = X.reshape(X.shape[0], config.X_feature_len)
    # print(X.shape, y.shape)
    #
    # # 使用 t-SNE 进行降维
    # tsne = TSNE(n_components=2, random_state=42)
    # X_embedded = tsne.fit_transform(X)
    #
    # # 绘制降维后的数据
    # plt.figure(figsize=(8, 6))
    # for label in np.unique(y):
    #     plt.scatter(X_embedded[y == label, 0], X_embedded[y == label, 1], label=label)
    # plt.title('t-SNE Visualization of Iris Dataset')
    # plt.legend()
    # plt.show()


    # # 生成Repbase标签的GFF注释文件
    # # 将水稻的gff文件染色体加上Chr，并且将TE label加到里面
    # work_dir = '/home/hukang/NeuralTE_experiment/Dataset5/rice_annotation'
    # TE_path = work_dir + '/test.ref'
    # contigNames, contigs = read_fasta_v1(TE_path)
    # label_dict = {}
    # for name in contigNames:
    #     parts = name.split('\t')
    #     seq_name = parts[0]
    #     label = parts[1]
    #     label_dict[seq_name] = label
    # rice_gff = work_dir + '/rice.gff'
    # rice_repbase_gff = work_dir + '/rice_repbase.gff'
    # chromosomes = [f'Chr{i}' for i in range(1, 13)]
    # generate_predict_gff(rice_gff, rice_repbase_gff, label_dict, chromosomes)
    # total_genome_len = 374424240
    # repbase_work_dir = work_dir + '/repbase'
    # analyze_class_ratio_gff(rice_repbase_gff, repbase_work_dir, total_genome_len)

    # # 生成NeuralTE的标签
    # work_dir = '/home/hukang/NeuralTE_experiment/Dataset5/rice_annotation'
    # NeuralTE_pred_path = work_dir + '/classified.info'
    # label_dict = {}
    # with open(NeuralTE_pred_path, 'r') as f_r:
    #     for line in f_r:
    #         if not line.startswith('#'):
    #             line = line.replace('\n', '')
    #             parts = line.split(',')
    #             label_dict[parts[0]] = parts[2]
    # rice_gff = work_dir + '/rice.gff'
    # rice_NeuralTE_gff = work_dir + '/rice_NeuralTE.gff'
    # chromosomes = [f'Chr{i}' for i in range(1, 13)]
    # generate_predict_gff(rice_gff, rice_NeuralTE_gff, label_dict, chromosomes)
    # total_genome_len = 374424240
    # NeuralTE_work_dir = work_dir + '/NeuralTE'
    # analyze_class_ratio_gff(rice_NeuralTE_gff, NeuralTE_work_dir, total_genome_len)

    # # 生成TERL的标签
    # work_dir = '/home/hukang/NeuralTE_experiment/Dataset5/rice_annotation'
    # TERL_pred_path = '/home/hukang/TE_Classification/TERL/TERL_20240111_202829_test.ref'
    # names, contigs = read_fasta_v1(TERL_pred_path)
    # test_path = '/home/hukang/TE_Classification/TERL/Data/DS5/test.ref'
    # test_names, test_contigs = read_fasta(test_path)
    # label_dict = {}
    # names, contigs = read_fasta_v1(TERL_pred_path)
    # for i, name in enumerate(names):
    #     parts = name.split('\t')
    #     label = parts[-2]
    #     label_dict[test_names[i]] = label
    # rice_gff = work_dir + '/rice.gff'
    # rice_TERL_gff = work_dir + '/rice_TERL.gff'
    # chromosomes = [f'Chr{i}' for i in range(1, 13)]
    # generate_predict_gff(rice_gff, rice_TERL_gff, label_dict, chromosomes)
    # total_genome_len = 374424240
    # TERL_work_dir = work_dir + '/TERL'
    # analyze_class_ratio_gff(rice_TERL_gff, TERL_work_dir, total_genome_len)

    # # 生成TEsorter的标签
    # TEsorter_pred_path = '/home/hukang/NeuralTE_experiment/Dataset5/TEsorter/test.ref.rexdb.cls.lib'
    # label_dict = {'EnSpm_CACTA': 'CACTA', 'pararetrovirus': 'Retrovirus', 'LINE': 'Unknown',
    #               'MuDR_Mutator': 'Mutator', 'mixture': 'Unknown', 'Tc1_Mariner': 'Tc1-Mariner',
    #               'PIF_Harbinger': 'PIF-Harbinger'}
    # pred_names, pred_contigs = read_fasta(TEsorter_pred_path)
    # pred_label_dict = {}
    # pred_label_set = set()
    # for name in pred_names:
    #     parts = name.split('#')
    #     raw_name = parts[0]
    #     label = parts[1]
    #     pred_label_set.add(label)
    #     label_parts = label.split('/')
    #     if len(label_parts) >= 2:
    #         label = label_parts[1]
    #     if label_dict.__contains__(label):
    #         label = label_dict[label]
    #     pred_label_dict[raw_name] = label
    # print(pred_label_set)
    #
    # work_dir = '/home/hukang/NeuralTE_experiment/Dataset5/rice_annotation'
    # rice_gff = work_dir + '/rice.gff'
    # rice_TEsorter_gff = work_dir + '/rice_TEsorter.gff'
    # chromosomes = [f'Chr{i}' for i in range(1, 13)]
    # generate_predict_gff(rice_gff, rice_TEsorter_gff, pred_label_dict, chromosomes)
    # total_genome_len = 374424240
    # TEsorter_work_dir = work_dir + '/TEsorter'
    # analyze_class_ratio_gff(rice_TEsorter_gff, TEsorter_work_dir, total_genome_len)

    # # 生成DeepTE的标签
    # DeepTE_pred_path = '/home/hukang/NeuralTE_experiment/Dataset5/DeepTE/results/opt_DeepTE.fasta'
    # DeepTE_labels = {'ClassII_DNA_Mutator_unknown': 'Mutator', 'ClassII_DNA_TcMar_nMITE': 'Tc1-Mariner',
    #                  'ClassII_DNA_hAT_unknown': 'hAT', 'ClassII_DNA_P_MITE': 'P', 'ClassI_nLTR': 'Unknown',
    #                  'ClassIII_Helitron': 'Helitron', 'ClassI_LTR_Gypsy': 'Gypsy', 'ClassI_LTR': 'Unknown',
    #                  'ClassII_DNA_Mutator_MITE': 'Mutator', 'ClassI_LTR_Copia': 'Copia', 'ClassI_nLTR_LINE': 'Unknown',
    #                  'ClassII_DNA_CACTA_unknown': 'CACTA', 'ClassI_nLTR_LINE_I': 'I', 'ClassI_nLTR_DIRS': 'DIRS',
    #                  'ClassII_MITE': 'Unknown', 'unknown': 'Unknown', 'ClassII_DNA_TcMar_unknown': 'Tc1-Mariner',
    #                  'ClassII_DNA_CACTA_MITE': 'CACTA', 'ClassII_DNA_Harbinger_unknown': 'PIF-Harbinger',
    #                  'ClassII_DNA_hAT_nMITE': 'hAT', 'ClassI': 'Unknown', 'ClassI_nLTR_SINE_7SL': '7SL',
    #                  'ClassII_DNA_Harbinger_nMITE': 'PIF-Harbinger', 'ClassII_DNA_Mutator_nMITE': 'Mutator',
    #                  'ClassII_DNA_hAT_MITE': 'hAT', 'ClassII_DNA_CACTA_nMITE': 'CACTA', 'ClassI_nLTR_SINE_tRNA': 'tRNA',
    #                  'ClassII_DNA_TcMar_MITE': 'Tc1-Mariner', 'ClassII_DNA_P_nMITE': 'P', 'ClassI_nLTR_PLE': 'Penelope',
    #                  'ClassII_DNA_Harbinger_MITE': 'PIF-Harbinger', 'ClassI_nLTR_LINE_L1': 'L1',
    #                  'ClassII_nMITE': 'Unknown',
    #                  'ClassI_LTR_ERV': 'Retrovirus', 'ClassI_LTR_BEL': 'Bel-Pao', 'ClassI_nLTR_LINE_RTE': 'RTE',
    #                  'ClassI_nLTR_LINE_R2': 'R2',
    #                  'ClassII_DNA_Transib_nMITE': 'Transib', 'ClassII_DNA_PiggyBac_nMITE': 'PiggyBac',
    #                  'ClassI_nLTR_LINE_Jockey': 'Jockey',
    #                  'ClassI_nLTR_SINE_5S': '5S', 'ClassII_DNA_hAT': 'hAT', 'ClassII_DNA_Tc1-Mariner': 'Tc1-Mariner',
    #                  'ClassII_DNA_PIF-Harbinger': 'PIF-Harbinger', 'ClassII_DNA_CACTA': 'CACTA',
    #                  'ClassII_DNA_Mutator': 'Mutator',
    #                  'ClassII_DNA_P': 'P', 'ClassII_DNA_Crypton': 'Crypton', 'ClassII_DNA_Transib': 'Transib',
    #                  'ClassI_LTR_Retrovirus': 'Retrovirus',
    #                  'ClassI_LTR_Bel-Pao': 'Bel-Pao', 'ClassII_DNA_Merlin': 'Merlin'}
    # names, contigs = read_fasta(DeepTE_pred_path)
    # label_dict = {}
    # for name in names:
    #     parts = name.split('__')
    #     seq_name = parts[0]
    #     label = parts[1]
    #     wicker_label = DeepTE_labels[label]
    #     label_dict[seq_name] = wicker_label
    #
    # work_dir = '/home/hukang/NeuralTE_experiment/Dataset5/rice_annotation'
    # rice_gff = work_dir + '/rice.gff'
    # rice_DeepTE_gff = work_dir + '/rice_DeepTE.gff'
    # chromosomes = [f'Chr{i}' for i in range(1, 13)]
    # generate_predict_gff(rice_gff, rice_DeepTE_gff, label_dict, chromosomes)
    # total_genome_len = 374424240
    # DeepTE_work_dir = work_dir + '/DeepTE'
    # analyze_class_ratio_gff(rice_DeepTE_gff, DeepTE_work_dir, total_genome_len)

    # non_TE = ('tandem repeat', 'Tandem repeat',
    #           'MSAT', 'SAT', 'Satellite repetitive element', 'satellite', 'Satellite',
    #           'Simple Repeat', 'Multicopy gene', 'Pseudogene')
    # # 收集repbase中的标签
    # label_set = set()
    # startpath = '/home/hukang/RepBase28.06.fasta'
    # non_TE_path = '/home/hukang/NeuralTE/data/non_TE.fa'
    # non_TE_contigs = {}
    # for root, dirs, files in os.walk(startpath):
    #     # 遍历当前目录下的文件
    #     for filename in files:
    #         file_path = os.path.join(root, filename)
    #         names, contigs = read_fasta_v1(file_path)
    #         for name in names:
    #             parts = name.split('\t')
    #             if len(parts) == 3:
    #                 label = parts[1]
    #                 if label in non_TE:
    #                     non_TE_contigs[name] = contigs[name]
    #                     print(file_path)
    #                 label_set.add(label)
    # print(label_set)
    # store_fasta(non_TE_contigs, non_TE_path)

    # # 生成RepeatClassifier的标签
    # RepeatClassifier_pred_path = '/home/hukang/NeuralTE_experiment/novel_TE/maize/RC/novel_tir.fa.classified'
    # rmToWicker = {}
    # wicker_superfamily_set = set()
    # with open(config.project_dir + '/data/TEClasses.tsv', 'r') as f_r:
    #     for i, line in enumerate(f_r):
    #         parts = line.split('\t')
    #         rm_type = parts[5]
    #         rm_subtype = parts[6]
    #         repbase_type = parts[7]
    #         wicker_type = parts[8]
    #         wicker_type_parts = wicker_type.split('/')
    #         # print(rm_type + ',' + rm_subtype + ',' + repbase_type + ',' + wicker_type)
    #         if len(wicker_type_parts) != 3:
    #             continue
    #         wicker_superfamily_parts = wicker_type_parts[-1].strip().split(' ')
    #         if len(wicker_superfamily_parts) == 1:
    #             wicker_superfamily = wicker_superfamily_parts[0]
    #         elif len(wicker_superfamily_parts) > 1:
    #             wicker_superfamily = wicker_superfamily_parts[1].replace('(', '').replace(')', '')
    #         rm_full_type = rm_type + '/' + rm_subtype
    #         if wicker_superfamily == 'ERV':
    #             wicker_superfamily = 'Retrovirus'
    #         if wicker_superfamily == 'Viper':
    #             wicker_superfamily = 'VIPER'
    #         if wicker_superfamily == 'H':
    #             wicker_superfamily = 'Helitron'
    #         rmToWicker[rm_full_type] = wicker_superfamily
    #         wicker_superfamily_set.add(wicker_superfamily)
    # # 补充一些元素
    # rmToWicker['LINE/R2'] = 'R2'
    # rmToWicker['LINE/Tad1'] = 'I'
    # rmToWicker['LINE?/L1'] = 'L1'
    # rmToWicker['LINE/CR1'] = 'I'
    # rmToWicker['DNA/PIF'] = 'PIF-Harbinger'
    # rmToWicker['SINE/ID'] = 'tRNA'
    # rmToWicker['SINE/MIR'] = 'tRNA'
    # rmToWicker['SINE/tRNA-Deu-I'] = 'tRNA'
    # rmToWicker['DNA/CMC'] = 'CACTA'
    # rmToWicker['DNA?/hAT'] = 'hAT'
    # rmToWicker['LTR/ERVL'] = 'Retrovirus'
    # rmToWicker['LINE/R2-NeSL'] = 'R2'
    # rmToWicker['DNA/Zator'] = 'Tc1-Mariner'
    # rmToWicker['Unknown'] = 'Unknown'
    #
    # # 因为数据集里面没有这四类Ngaro，VIPER，Maverick和PiggyBac，所以我们将预测成这些类别的标记为Unknown
    # filter_labels = ('Ngaro', 'VIPER', 'Maverick', 'PiggyBac')
    #
    # ## 2.3 获取RepeatClassifier分类后的序列标签，对于未能标注到superfamily的标签或者错误的标签，我们直接标为Unknown （因为我们这里的数据集标签是直接打到了superfamily层级）
    # names, contigs = read_fasta_v1(RepeatClassifier_pred_path)
    # RC_name_labels = {}
    # label_nums = {}
    # all_unique_RM_label = set()
    # for name in names:
    #     label = name.split('#')[1].split(' ')[0]
    #     if not rmToWicker.__contains__(label):
    #         all_unique_RM_label.add(label)
    #         label = 'Unknown'
    #     else:
    #         wicker_superfamily = rmToWicker[label]
    #         label = wicker_superfamily
    #         if label in filter_labels:
    #             label = 'Unknown'
    #         all_unique_RM_label.add(label)
    #     RC_name_labels[name.split('#')[0]] = label
    #     if not label_nums.__contains__(label):
    #         label_nums[label] = 0
    #     num = label_nums[label]
    #     label_nums[label] = num + 1
    # print(label_nums)
    # print('all_unique_RM_label:' + str(all_unique_RM_label))
    #
    # #chromosomes = [f'Chr{i}' for i in range(1, 13)]
    # chromosomes = [f'chr_{i}' for i in range(0, 100)]
    # rice_gff = '/home/hukang/NeuralTE_experiment/novel_TE/maize/RC/novel_tir.gff.bak'
    # rice_RepeatClassifier_gff = '/home/hukang/NeuralTE_experiment/novel_TE/maize/RC/novel_tir.RC.gff'
    # generate_predict_gff(rice_gff, rice_RepeatClassifier_gff, RC_name_labels, chromosomes)
    # #total_genome_len = 374424240
    # total_genome_len = 2182786008
    # work_dir = '/home/hukang/NeuralTE_experiment/novel_TE/maize/RC'
    # analyze_class_ratio_gff(rice_RepeatClassifier_gff, work_dir, total_genome_len)



    # # #识别每一种superfamily，存成文件
    # work_dir = '/home/hukang/NeuralTE_experiment/novel_TE/maize'
    # consensus_path = work_dir + '/novel_tir.fa'
    # names, contigs = read_fasta(consensus_path)
    #
    # classified_info = work_dir + '/classified.info'
    # name_label_dict = {}
    # label_nums = {}
    # with open(classified_info, 'r') as f_r:
    #     for line in f_r:
    #         if line.startswith('#'):
    #             continue
    #         parts = line.replace('\n', '').split(',')
    #         raw_name = parts[0]
    #         label = parts[2]
    #         name_label_dict[raw_name] = label
    #         if not label_nums.__contains__(label):
    #             label_nums[label] = 0
    #         num = label_nums[label]
    #         label_nums[label] = num + 1
    # print(label_nums)
    # RC_name_labels = name_label_dict
    #
    # #chromosomes = [f'Chr{i}' for i in range(1, 13)]
    # chromosomes = [f'chr_{i}' for i in range(0, 100)]
    # rice_gff = work_dir + '/novel_tir.gff.bak'
    # rice_NeuralTE_gff = work_dir + '/novel_tir.NeuralTE.gff'
    # generate_predict_gff(rice_gff, rice_NeuralTE_gff, RC_name_labels, chromosomes)
    # # total_genome_len = 374424240
    # total_genome_len = 2182786008
    # NeuralTE_work_dir = work_dir + '/NeuralTE'
    # analyze_class_ratio_gff(rice_NeuralTE_gff, NeuralTE_work_dir, total_genome_len)
    #
    # superfamilies_contigs = {}
    # for name in RC_name_labels.keys():
    #     label = RC_name_labels[name]
    #     if not superfamilies_contigs.__contains__(label):
    #         superfamilies_contigs[label] = {}
    #     cur_superfamily_contigs = superfamilies_contigs[label]
    #     new_name = name + '#' + label
    #     cur_superfamily_contigs[new_name] = contigs[name]
    # for label in superfamilies_contigs.keys():
    #     file_path = work_dir + '/' + label + '.fa'
    #     store_fasta(superfamilies_contigs[label], file_path)
    #
    # # 分析每一种novel TIR的插入时间
    # miu = 1.3e-8
    # TE_list = ['CACTA', 'hAT', 'Mutator', 'PIF-Harbinger', 'Tc1-Mariner']
    # colors = ['#B77072', '#F08B47', '#68A47F', '#B3A4BB', '#B6E1C5']
    # for i in range(len(TE_list)):
    #     type = TE_list[i]
    #     tir_path = work_dir + '/'+type+'.fa'
    #     analyz_TIR_insert_time(tir_path, work_dir, miu, type, colors[i])
    #
    # # 将RepeatMasker的注释转为bed文件，并过滤掉非全长TE
    # RMOut = work_dir + '/novel_tir.out'
    # out_bed = work_dir + '/novel_tir.full_length.bed'
    # consensus_path = work_dir + '/novel_tir.fa'
    # tools_dir = config.project_dir + '/tools'
    # coverage_threshold = 0.95
    # transfer_RMOut2Bed(RMOut, out_bed, consensus_path, tools_dir, coverage_threshold, RC_name_labels)

    # # 存储基因组染色体长度
    # genome_path = '/home/hukang/Genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.rename.fna'
    # chromosomes = [f'Chr{i}' for i in range(1, 11)]
    # names, contigs = read_fasta(genome_path)
    # genome_len_info = work_dir + '/Genome_len.chr'
    # with open(genome_len_info, 'w') as f_save:
    #     f_save.write('Genome Length\n')
    #     for name in names:
    #         index = int(name.split('_')[1])
    #         new_name = 'Chr' + str(index+1)
    #         if new_name not in chromosomes:
    #             continue
    #         chr_len = len(contigs[name])
    #         f_save.write(new_name + '\t' + str(chr_len) + '\n')


    # # 存储基因组染色体长度
    # genome_path = '/home/hukang/Genome/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa'
    # chromosomes = [f'Chr{i}' for i in range(1, 13)]
    # names, contigs = read_fasta(genome_path)
    # genome_len_info = work_dir + '/Genome_len.chr'
    # with open(genome_len_info, 'w') as f_save:
    #     f_save.write('Genome Length\n')
    #     for name in names:
    #         new_name = 'Chr' + str(name)
    #         if new_name not in chromosomes:
    #             continue
    #         chr_len = len(contigs[name])
    #         f_save.write(new_name+'\t'+str(chr_len)+'\n')

    # file_path = '/home/hukang/NeuralTE_dataset/Repbase_raw/all_repbase.ref.raw'
    # names, contigs = read_fasta_v1(file_path)
    # count = 0
    # for name in names:
    #     parts = name.split('\t')
    #     seq_len = len(contigs[name])
    #     if seq_len < 80:
    #         del contigs[name]
    # store_fasta(contigs, file_path)


    # # 新增实验，在水稻上分析RepeatClassifier和NeuralTE分类标签数量交叉
    # work_dir = '/home/hukang/NeuralTE_experiment/Dataset5'
    # NeuralTE_results = work_dir + '/NeuralTE/classified.info'
    # RC_results = work_dir + '/RC/test.ref.classified'
    #
    # rmToWicker = {}
    # wicker_superfamily_set = set()
    # with open(config.project_dir + '/data/TEClasses.tsv', 'r') as f_r:
    #     for i, line in enumerate(f_r):
    #         parts = line.split('\t')
    #         rm_type = parts[5]
    #         rm_subtype = parts[6]
    #         repbase_type = parts[7]
    #         wicker_type = parts[8]
    #         wicker_type_parts = wicker_type.split('/')
    #         # print(rm_type + ',' + rm_subtype + ',' + repbase_type + ',' + wicker_type)
    #         if len(wicker_type_parts) != 3:
    #             continue
    #         wicker_superfamily_parts = wicker_type_parts[-1].strip().split(' ')
    #         if len(wicker_superfamily_parts) == 1:
    #             wicker_superfamily = wicker_superfamily_parts[0]
    #         elif len(wicker_superfamily_parts) > 1:
    #             wicker_superfamily = wicker_superfamily_parts[1].replace('(', '').replace(')', '')
    #         rm_full_type = rm_type + '/' + rm_subtype
    #         if wicker_superfamily == 'ERV':
    #             wicker_superfamily = 'Retrovirus'
    #         if wicker_superfamily == 'Viper':
    #             wicker_superfamily = 'VIPER'
    #         if wicker_superfamily == 'H':
    #             wicker_superfamily = 'Helitron'
    #         rmToWicker[rm_full_type] = wicker_superfamily
    #         wicker_superfamily_set.add(wicker_superfamily)
    # # Supplement some elements
    # rmToWicker['LINE/R2'] = 'R2'
    # rmToWicker['LINE/Tad1'] = 'I'
    # rmToWicker['LINE?/L1'] = 'L1'
    # rmToWicker['LINE/CR1'] = 'I'
    # rmToWicker['DNA/PIF'] = 'PIF-Harbinger'
    # rmToWicker['SINE/ID'] = 'tRNA'
    # rmToWicker['SINE/MIR'] = 'tRNA'
    # rmToWicker['SINE/tRNA-Deu-I'] = 'tRNA'
    # rmToWicker['DNA/CMC'] = 'CACTA'
    # rmToWicker['DNA?/hAT'] = 'hAT'
    # rmToWicker['LTR/ERVL'] = 'Retrovirus'
    # rmToWicker['LINE/R2-NeSL'] = 'R2'
    # rmToWicker['DNA/Zator'] = 'Tc1-Mariner'
    # rmToWicker['Unknown'] = 'Unknown'
    #
    # print(rmToWicker)
    # print(wicker_superfamily_set)
    # print(len(wicker_superfamily_set))
    #
    # # As the dataset lacks these four types: Ngaro, VIPER, Maverick, and PiggyBac,
    # # instances predicted as these categories are marked as Unknown
    # filter_labels = ('Ngaro', 'VIPER', 'Maverick', 'PiggyBac')
    #
    # ## 2.3 Retrieve labels classified by RepeatClassifier; for labels that didn't
    # # annotate to the superfamily or were incorrect, label them as Unknown
    # # (as our dataset labels are at the superfamily level)
    # RC_list = []
    # names, contigs = read_fasta_v1(RC_results)
    # RC_name_labels = {}
    # all_unique_RM_label = set()
    # for name in names:
    #     label = name.split('#')[1].split(' ')[0]
    #     if not rmToWicker.__contains__(label):
    #         all_unique_RM_label.add(label)
    #         label = 'Unknown'
    #     else:
    #         wicker_superfamily = rmToWicker[label]
    #         label = wicker_superfamily
    #         if label in filter_labels:
    #             label = 'Unknown'
    #         all_unique_RM_label.add(label)
    #     RC_name_labels[name.split('#')[0]] = label
    #     RC_list.append(name.split('#')[0]+'#'+label)
    # print('all_unique_RM_label:' + str(all_unique_RM_label))
    # print(RC_name_labels)
    #
    # gold_standard_list = []
    # NeuralTE_list = []
    # NeuralTE_name_labels = {}
    # with open(NeuralTE_results, 'r') as f_r:
    #     for line in f_r:
    #         if line.startswith('#'):
    #             continue
    #         line = line.replace('\n', '')
    #         parts = line.split(',')
    #         raw_name = parts[0]
    #         gold_standard = parts[1]
    #         label = parts[2]
    #         NeuralTE_name_labels[raw_name] = (gold_standard, label)
    #         gold_standard_list.append(raw_name+'#'+gold_standard)
    #         NeuralTE_list.append(raw_name + '#' + label)
    #
    # # 1.金标准、RC、NeuralTE共有
    # count1 = 0
    # # 2. 金标准与RC共有、NeuralTE不同
    # count2 = 0
    # # 3. 金标准与NeuralTE共有、RC不同
    # count3 = 0
    # # 4. RC与NeuralTE共有、金标准不同
    # count4 = 0
    # # 5. 金标准、RC、NeuralTE均不相同
    # count5 = 0
    # # 6. 金标准与RC和NeuralTE均不相同
    # count6 = 0
    # # 7. RC与金标准和NeuralTE都不同
    # count7 = 0
    # # 8. NeuralTE 与金标准和 RC 都不同
    # count8 = 0
    # for raw_name in NeuralTE_name_labels.keys():
    #     gold_standard, NeuralTE_label = NeuralTE_name_labels[raw_name]
    #     RC_label = RC_name_labels[raw_name]
    #     if gold_standard == NeuralTE_label and gold_standard == RC_label:
    #         count1 += 1
    #     elif gold_standard == RC_label and gold_standard != NeuralTE_label:
    #         count2 += 1
    #     elif gold_standard != RC_label and gold_standard == NeuralTE_label:
    #         count3 += 1
    #     elif gold_standard != RC_label and gold_standard != NeuralTE_label and NeuralTE_label == RC_label:
    #         count4 += 1
    #     elif gold_standard != RC_label and gold_standard != NeuralTE_label and NeuralTE_label != RC_label:
    #         count5 += 1
    #
    #     if gold_standard != RC_label and gold_standard != NeuralTE_label:
    #         count6 += 1
    #     if RC_label != gold_standard and RC_label != NeuralTE_label:
    #         count7 += 1
    #     if NeuralTE_label != gold_standard and NeuralTE_label != RC_label:
    #         count8 += 1
    # print('总共数量:' + str(len(NeuralTE_name_labels)))
    # print('金标准、RC、NeuralTE共有:' + str(count1))
    # print('金标准与RC共有、NeuralTE不同:' + str(count2))
    # print('金标准与NeuralTE共有、RC不同:' + str(count3))
    # print('RC与NeuralTE共有、金标准不同:' + str(count4))
    # print('金标准、RC、NeuralTE均不相同:' + str(count5))
    # print('金标准与RC和NeuralTE均不相同:' + str(count6))
    # print('RC与金标准和NeuralTE都不同:' + str(count7))
    # print('NeuralTE与金标准和RC都不同:' + str(count8))
    #
    # # import matplotlib.pyplot as plt
    # # from matplotlib_venn import venn3
    # #
    # # # 输入数据
    # # venn_labels = {'100': 118, '010': 831, '001': 198, '110': 92, '101': 725, '011': 12, '111': 1564}
    # #
    # # # 画韦恩图
    # # venn_diagram = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels=('Gold Standard', 'RepeatClassifier', 'NeuralTE'))
    # # venn_diagram.get_label_by_id('100').set_text(venn_labels['100'])
    # # venn_diagram.get_label_by_id('010').set_text(venn_labels['010'])
    # # venn_diagram.get_label_by_id('001').set_text(venn_labels['001'])
    # # venn_diagram.get_label_by_id('110').set_text(venn_labels['110'])
    # # venn_diagram.get_label_by_id('101').set_text(venn_labels['101'])
    # # venn_diagram.get_label_by_id('011').set_text(venn_labels['011'])
    # # venn_diagram.get_label_by_id('111').set_text(venn_labels['111'])
    # # output_fig = '/home/hukang/NeuralTE_experiment/Dataset5/cross_num.png'
    # # # 显示图形
    # # plt.tight_layout()
    # # # 显示图形
    # # # plt.show()
    # # plt.savefig(output_fig, format='png')
    #
    # from matplotlib import pyplot as plt
    # from matplotlib_venn import venn3
    # # 假设有两组数据
    # set1 = set(gold_standard_list)
    # set2 = set(NeuralTE_list)
    # set3 = set(RC_list)
    # # 创建一个subplot
    # fig, ax = plt.subplots()
    # # 绘制韦恩图
    # #venn_diagram = venn3([set1, set2, set3], ('Gold Standard', 'NeuralTE', 'RepeatClassifier'), set_colors=('#E66255', '#299D92', '#FFCE71'), alpha=0.9)
    # venn_diagram = venn3([set1, set2, set3], ('Gold Standard', 'NeuralTE', 'RepeatClassifier'))
    # # plt.legend(handles=[venn_diagram.get_patch_by_id('100'),  # Set 1
    # #                     venn_diagram.get_patch_by_id('010'),  # Set 2
    # #                     venn_diagram.get_patch_by_id('001')],  # Set 3
    # #            labels=['Gold Standard', 'NeuralTE', 'RepeatClassifier'], loc='lower right')
    #
    # plt.tight_layout()
    # # 设置图形标题
    # #plt.title("")
    # output_fig = '/home/hukang/NeuralTE_experiment/Dataset5/cross_num.png'
    # # 显示图形
    # plt.show()
    # #plt.savefig(output_fig, format='png')
