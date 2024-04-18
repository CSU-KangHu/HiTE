#-- coding: UTF-8 --
import os
import random
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

from src.DataProcessor import DataProcessor

current_folder = os.path.dirname(os.path.abspath(__file__))
# 添加 configs 文件夹的路径到 Python 路径
configs_folder = os.path.join(current_folder, "..")  # 需要根据实际目录结构调整
sys.path.append(configs_folder)

from configs import config
from utils.data_util import read_fasta, store_fasta, read_fasta_v1, replace_non_atcg, get_flanking_copies, \
    get_copies_TSD_info, search_TSD_regular, extract_non_autonomous, run_command, \
    transfer_RMOut2Bed, generate_random_sequences, generate_random_sequence
from utils.evaluate_util import generate_TERL_dataset, generate_ClassifyTE_dataset, evaluate_RepeatClassifier, \
    evaluate_TERL, evaluate_DeepTE, transform_DeepTE_to_fasta, add_ClassifyTE_classification, evaluate_ClassifyTE, \
    evaluate_TEsorter, merge_overlapping_intervals, get_metrics_by_label, plot_3D_param, analyze_class_ratio_gff


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

if __name__ == '__main__':
    # # 重新训练DeepTE模型
    # # 1. 训练LINE模型
    # # 1.1 先提取Dataset2中的train.ref中的LINE元素对应的序列，转换成label,sequence格式
    # work_dir = '/home/hukang/NeuralTE_experiment/exclude_rice_maize_zebrafish/DeepTE'
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

    # data_dir = '/home/hukang/NeuralTE_experiment/Dataset6/DeepTE'
    # test_path = data_dir + '/test.ref'
    # predict_path = data_dir + '/results/opt_DeepTE.txt'
    # evaluate_DeepTE(test_path, predict_path)
    # # 将DeepTE的macro avg由19分类，变成13分类
    # indicators = [0.5568, 0.5852, 0.5468]
    # for ind in indicators:
    #     new_ind = 15 * ind / 11
    #     print(round(new_ind, 4))

    # # 替换非ATCG字符
    # data_dir = '/home/hukang/NeuralTE_experiment/Dataset3/DeepTE'
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

    # pred_path = '/home/hukang/NeuralTE_experiment/Dataset6/TEsorter/test.ref.rexdb.cls.lib'
    # test_path = '/home/hukang/NeuralTE_experiment/Dataset6/TEsorter/test.ref'
    # evaluate_TEsorter(pred_path, test_path)
    # # 将TEsorter的macro avg由25分类，变成24分类
    # indicators = [0.5827, 0.3707, 0.4395]
    # for ind in indicators:
    #     new_ind = 12 * ind / 11
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

    # work_dir = '/home/hukang/TE_Classification/TERL/Data/validate_TE'
    # fasta_file = work_dir + '/train.ref'
    # outdir = work_dir + '/Train'
    # generate_TERL_dataset(fasta_file, outdir)
    # fasta_file = work_dir + '/test.ref'
    # outdir = work_dir + '/Test'
    # generate_TERL_dataset(fasta_file, outdir)

    # fasta_file = '/home/hukang/NeuralTE_dataset/Dataset2/all_repbase.ref'
    # generate_ClassifyTE_dataset(fasta_file)

    # # # #获取RepeatClassifier的结果评估
    # classified_path = '/home/hukang/NeuralTE_experiment/Dataset6/RC/test.ref.classified'
    # RC_name_labels = evaluate_RepeatClassifier(classified_path)
    # # 将RepeatClassifier的macro avg由25分类，变成24分类
    # indicators = [0.7837, 0.7324, 0.7462]
    # for ind in indicators:
    #     new_ind = 13 * ind / 11
    #     print(round(new_ind, 4))

    # # 将NeuralTE的macro avg由19分类，变成13分类
    # indicators = [0.6291, 0.6841, 0.6487]
    # for ind in indicators:
    #     new_ind = 17 * ind / 13
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
    # train = '/home/hukang/NeuralTE_experiment/Dataset6/train.ref'
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
    # test_path = work_dir + '/Data/DS6/test.ref'
    # predict_path = work_dir + '/TERL_20240111_212627_test.ref'
    # evaluate_TERL(test_path, predict_path)
    # # 将TERL的macro avg由19分类，变成13分类
    # indicators = [0.4794, 0.4648, 0.4473]
    # for ind in indicators:
    #     new_ind = 15 * ind / 11
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
    # work_dir = '/home/hukang/NeuralTE/work/kmer_size_search'
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



    # #识别每一种superfamily，存成文件
    work_dir = '/home/hukang/NeuralTE_experiment/novel_TE/maize'
    consensus_path = work_dir + '/novel_tir.fa'
    names, contigs = read_fasta(consensus_path)

    classified_info = work_dir + '/classified.info'
    name_label_dict = {}
    label_nums = {}
    with open(classified_info, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            parts = line.replace('\n', '').split(',')
            raw_name = parts[0]
            label = parts[2]
            name_label_dict[raw_name] = label
            if not label_nums.__contains__(label):
                label_nums[label] = 0
            num = label_nums[label]
            label_nums[label] = num + 1
    print(label_nums)
    RC_name_labels = name_label_dict
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
    # 将RepeatMasker的注释转为bed文件，并过滤掉非全长TE
    RMOut = work_dir + '/novel_tir.out'
    out_bed = work_dir + '/novel_tir.full_length.bed'
    consensus_path = work_dir + '/novel_tir.fa'
    tools_dir = config.project_dir + '/tools'
    coverage_threshold = 0.95
    transfer_RMOut2Bed(RMOut, out_bed, consensus_path, tools_dir, coverage_threshold, RC_name_labels)

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
    # work_dir = '/home/hukang/NeuralTE_experiment/Dataset6'
    # NeuralTE_results = work_dir + '/classified.info'
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
    # print('all_unique_RM_label:' + str(all_unique_RM_label))
    # print(RC_name_labels)
    #
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
    # import matplotlib.pyplot as plt
    # from matplotlib_venn import venn3
    #
    # # 输入数据
    # venn_labels = {'100': 4, '010': 99, '001': 20, '110': 17, '101': 96, '011': 1, '111': 1198}
    #
    # # 画韦恩图
    # venn_diagram = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels=('Gold Standard', 'RepeatClassifier', 'NeuralTE'))
    # venn_diagram.get_label_by_id('100').set_text(venn_labels['100'])
    # venn_diagram.get_label_by_id('010').set_text(venn_labels['010'])
    # venn_diagram.get_label_by_id('001').set_text(venn_labels['001'])
    # venn_diagram.get_label_by_id('110').set_text(venn_labels['110'])
    # venn_diagram.get_label_by_id('101').set_text(venn_labels['101'])
    # venn_diagram.get_label_by_id('011').set_text(venn_labels['011'])
    # venn_diagram.get_label_by_id('111').set_text(venn_labels['111'])
    # output_fig = '/home/hukang/NeuralTE_experiment/Dataset5/cross_num.png'
    # # 显示图形
    # plt.tight_layout()
    # # 显示图形
    # # plt.show()
    # plt.savefig(output_fig, format='png')
