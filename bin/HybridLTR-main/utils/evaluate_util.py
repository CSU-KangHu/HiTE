import codecs
import json
import os.path
import re

import matplotlib
#matplotlib.use('pdf')
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, \
    classification_report
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import griddata, interpolate

from configs import config
from utils.data_util import read_fasta_v1, store_fasta, read_fasta, transfer_RMOut2BlastnOut

# import tensorflow as tf

# from keras import backend as K
def focal_loss(gamma=2.):
    def focal_loss_fixed(y_true, y_pred):
        pt_1 = tf.where(tf.equal(y_true, 1), y_pred, tf.ones_like(y_pred))
        return -K.sum(K.pow(1. - pt_1, gamma) * K.log(pt_1))
    return focal_loss_fixed

def plot_confusion_matrix(y_pred, y_test):

    y_pred_set = set(y_pred)
    y_test_set = set(y_test)
    class_list = list(y_pred_set | y_test_set)
    #class_list = list(y_test_set)

    # Compute confusion matrix
    # conf_matrix = confusion_matrix(y_test, y_pred, labels=class_names)
    class_report = classification_report(y_test, y_pred)
    print(class_report)

    # # Plot the confusion matrix
    # plt.figure(figsize=(15, 15))
    # sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', xticklabels=class_names, yticklabels=class_names)
    # plt.title('Confusion Matrix')
    # plt.xlabel('Predicted Labels')
    # plt.ylabel('True Labels')
    # plt.yticks(rotation=360)
    # plt.xticks(rotation=90)
    # plt.tight_layout()
    # plt.savefig(config.work_dir + '/confusion_matrix.png', format='png')
    # #plt.show()

def get_metrics(y_pred, y_test, row_nums_test, matrix_files_test):
    y_pred = np.argmax(np.round(y_pred), axis=1)

    # for i in range(len(y_test)):
    #     if y_test[i] != y_pred[i]:
    #         if y_test[i] == 0 and y_pred[i] == 1:
    #             print('FP: ' + matrix_files_test[i])
    #         else:
    #             print('FN: ' + matrix_files_test[i])

    accuracy = 0
    precision = 0
    recall = 0
    f1 = 0

    # compute metrics
    accuracy = round(accuracy_score(y_test, y_pred), 4)
    precision = round(precision_score(y_test, y_pred, average='macro'), 4)
    recall = round(recall_score(y_test, y_pred, average='macro'), 4)
    f1 = round(f1_score(y_test, y_pred, average='macro'), 4)

    # plot confusion matrix
    plot_confusion_matrix(y_test, y_pred)
    print('accuracy, precision, recall, f1:', accuracy, precision, recall, f1)

    y_copy_test = []
    y_copy_pred = []
    for i in range(len(y_test)):
        y_copy_test += [y_test[i]] * row_nums_test[i]
        y_copy_pred += [y_pred[i]] * row_nums_test[i]

    y_copy_test = np.array(y_copy_test)
    y_copy_pred = np.array(y_copy_pred)
    copy_accuracy = 0
    copy_precision = 0
    copy_recall = 0
    copy_f1 = 0

    # compute metrics
    copy_accuracy = round(accuracy_score(y_copy_test, y_copy_pred), 4)
    copy_precision = round(precision_score(y_copy_test, y_copy_pred, average='macro'), 4)
    copy_recall = round(recall_score(y_copy_test, y_copy_pred, average='macro'), 4)
    copy_f1 = round(f1_score(y_copy_test, y_copy_pred, average='macro'), 4)

    # plot confusion matrix
    plot_confusion_matrix(y_copy_test, y_copy_pred)
    print('considering copies, accuracy, precision, recall, f1:', copy_accuracy, copy_precision, copy_recall, copy_f1)

    return accuracy, precision, recall, f1, copy_accuracy, copy_precision, copy_recall, copy_f1


def get_metrics_v1(y_pred, y_test):
    # y_pred = np.argmax(np.round(y_pred), axis=1)
    # compute metrics
    accuracy = round(accuracy_score(y_test, y_pred), 4)
    precision = round(precision_score(y_test, y_pred, average='macro'), 4)
    recall = round(recall_score(y_test, y_pred, average='macro'), 4)
    f1 = round(f1_score(y_test, y_pred, average='macro'), 4)

    # plot confusion matrix
    plot_confusion_matrix(y_test, y_pred)
    print('accuracy, precision, recall, f1:', accuracy, precision, recall, f1)

    return accuracy, precision, recall, f1

def correct_using_minority(data_path, threads):
    # Error correction for NeuralTE predictions on minority samples
    # 1. Align the test dataset with the minority sample using RepeatMasker
    target_path = config.work_dir + '/minority/train.minority.ref'
    RMOut = data_path + '.out'
    RepeatMasker_command = 'RepeatMasker -lib ' + target_path + ' -no_is -norna -nolow -pa ' + str(threads) + ' ' + data_path + ' > /dev/null 2>&1'
    print(RepeatMasker_command)
    os.system(RepeatMasker_command)
    # 2. Convert RepeatMasker output to blastn out format for easier parsing
    test_minority_out = config.work_dir + '/minority/test.minority.out'
    tools_dir = config.project_dir + '/tools'
    transfer_RMOut2BlastnOut(RMOut, test_minority_out, tools_dir)
    # 3. Obtain alignment results, get the proportion of query sequence aligned with the target;
    # for > 80% query sequence alignments, assign query label as the target label
    query_names, query_contigs = read_fasta(data_path)
    target_names, target_contigs = read_fasta_v1(target_path)
    target_labels = {}
    target_len_dict = {}
    for name in target_names:
        parts = name.split('\t')
        target_name = parts[0]
        label = parts[1]
        target_labels[target_name] = label
        target_len_dict[target_name] = len(target_contigs[name])

    query_intervals = {}
    query_records = {}
    with open(test_minority_out, 'r') as f_r:
        for line in f_r:
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            identity = float(parts[2])
            query_start = int(parts[6])
            query_end = int(parts[7])
            subject_start = int(parts[8])
            subject_end = int(parts[9])
            e_value = float(parts[10])
            if subject_start > subject_end:
                temp = subject_start
                subject_start = subject_end
                subject_end = temp
            if e_value > 1e-10:
                continue
            if not query_intervals.__contains__(query_name):
                query_intervals[query_name] = {}
            target_intervals = query_intervals[query_name]
            if not target_intervals.__contains__(subject_name):
                target_intervals[subject_name] = []
            intervals = target_intervals[subject_name]
            intervals.append((query_start, query_end))

            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            target_records = query_records[query_name]
            if not target_records.__contains__(subject_name):
                target_records[subject_name] = []
            records = target_records[subject_name]
            records.append((query_start, query_end, subject_start, subject_end))

    query_labels = {}
    for query_name in query_intervals.keys():
        target_intervals = query_intervals[query_name]
        target_records = query_records[query_name]
        for subject_name in target_intervals.keys():
            records = target_records[subject_name]
            target_label = target_labels[subject_name]
            intervals = target_intervals[subject_name]
            merge_intervals = merge_overlapping_intervals(intervals)
            # Calculate total aligned length
            sum_len = 0
            for interval in merge_intervals:
                sum_len += abs(interval[1] - interval[0])
            query_len = len(query_contigs[query_name])
            subject_len = target_len_dict[subject_name]
            alignment_ratio = float(sum_len) / query_len
            if alignment_ratio > 0.8:
                if not query_labels.__contains__(query_name):
                    query_labels[query_name] = target_label
            elif target_label == 'P' or target_label == 'Merlin':
                # If the target is a DNA transposon and the query's terminal aligns with the target's terminal, consider it
                # Query alignment is at the 5'-end or 3'-end, and subject alignment is also at the 5'-end or 3'-end
                for record in records:
                    if ((record[0] - 1) <= 5 or (query_len - record[1]) <= 5) and (
                            (record[2] - 1) <= 5 or (subject_len - record[3]) <= 5):
                        query_labels[query_name] = target_label
                        break

    raw_NeuralTE_result = config.work_dir + '/classified.info'
    # correct results
    y_pred = []
    y_test = []
    with open(raw_NeuralTE_result, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            line = line.replace('\n', '')
            parts = line.split(',')
            seq_name = parts[0]
            true_label = parts[1]
            pred_label = parts[2]
            if query_labels.__contains__(seq_name):
                pred_label = query_labels[seq_name]
            y_pred.append(pred_label)
            y_test.append(true_label)
    y_test = np.array(y_test)
    y_pred = np.array(y_pred)
    get_metrics_by_label(y_test, y_pred)


# Incorporate features into ClassifyTE training data
# Due to ClassifyTE's lack of inclusion of classification labels in feature
# extraction from training data, additional processing is required to incorporate
# classification columns when retraining the model
def add_ClassifyTE_classification(feature_path, list_path, list_data_dir):
    new_feature_path = feature_path + '.train'
    labels = ['classification']
    with open(list_path, 'r') as f_r:
        for i, line in enumerate(f_r):
            seq_name = line.replace('\n', '')
            seq_path = list_data_dir + '/' + seq_name
            names, contigs = read_fasta_v1(seq_path)
            name = names[0]
            label = name.split('\t')[1]
            label_num = config.ClassifyTE_class[label]
            labels.append(label_num)

    with open(new_feature_path, 'w') as f_save:
        with open(feature_path, 'r') as f_r:
            for i, line in enumerate(f_r):
                line = line.replace('\n', '')
                newline = line + ',' + labels[i]
                f_save.write(newline+'\n')

# Filter Repbase dataset, keeping only types included in TE_class
def filterRepbase(repbase_path, TE_class):
    filter_contigs = {}
    names, contigs = read_fasta_v1(repbase_path)
    for name in names:
        label = name.split('\t')[1]
        if label in TE_class:
            filter_contigs[name] = contigs[name]
    filter_path = repbase_path+'.filter'
    store_fasta(filter_contigs, filter_path)
    return filter_path

def get_metrics_by_label(y_test, y_pred):
    accuracy = accuracy_score(y_test, y_pred)
    print("Accuracy:", round(accuracy, 4))
    precision = precision_score(y_test, y_pred, average='macro')
    print("Precision:", round(precision, 4))
    recall = recall_score(y_test, y_pred, average='macro')
    print("Recall:", round(recall, 4))
    f1 = f1_score(y_test, y_pred, average='macro')
    print("F1:", round(f1, 4))

    plot_confusion_matrix(y_pred, y_test)

def evaluate_RepeatClassifier(classified_path):
    # 2. Evaluate RepeatClassifier on the test dataset
    # 2.2 Convert Dfam classification names to Wicker format
    # 2.2.1 This file contains conversions between RepeatMasker categories, Repbase, and Wicker categories
    rmToWicker = {}
    wicker_superfamily_set = set()
    with open(config.project_dir + '/data/TEClasses.tsv', 'r') as f_r:
        for i,line in enumerate(f_r):
            parts = line.split('\t')
            rm_type = parts[5]
            rm_subtype = parts[6]
            repbase_type = parts[7]
            wicker_type = parts[8]
            wicker_type_parts = wicker_type.split('/')
            #print(rm_type + ',' + rm_subtype + ',' + repbase_type + ',' + wicker_type)
            if len(wicker_type_parts) != 3:
                continue
            wicker_superfamily_parts = wicker_type_parts[-1].strip().split(' ')
            if len(wicker_superfamily_parts) == 1:
                wicker_superfamily = wicker_superfamily_parts[0]
            elif len(wicker_superfamily_parts) > 1:
                wicker_superfamily = wicker_superfamily_parts[1].replace('(', '').replace(')', '')
            rm_full_type = rm_type+'/'+rm_subtype
            if wicker_superfamily == 'ERV':
                wicker_superfamily = 'Retrovirus'
            if wicker_superfamily == 'Viper':
                wicker_superfamily = 'VIPER'
            if wicker_superfamily == 'H':
                wicker_superfamily = 'Helitron'
            rmToWicker[rm_full_type] = wicker_superfamily
            wicker_superfamily_set.add(wicker_superfamily)
    # Supplement some elements
    rmToWicker['LINE/R2'] = 'R2'
    rmToWicker['LINE/Tad1'] = 'I'
    rmToWicker['LINE?/L1'] = 'L1'
    rmToWicker['LINE/CR1'] = 'I'
    rmToWicker['DNA/PIF'] = 'PIF-Harbinger'
    rmToWicker['SINE/ID'] = 'tRNA'
    rmToWicker['SINE/MIR'] = 'tRNA'
    rmToWicker['SINE/tRNA-Deu-I'] = 'tRNA'
    rmToWicker['DNA/CMC'] = 'CACTA'
    rmToWicker['DNA?/hAT'] = 'hAT'
    rmToWicker['LTR/ERVL'] = 'Retrovirus'
    rmToWicker['LINE/R2-NeSL'] = 'R2'
    rmToWicker['DNA/Zator'] = 'Tc1-Mariner'
    rmToWicker['Unknown'] = 'Unknown'

    print(rmToWicker)
    print(wicker_superfamily_set)
    print(len(wicker_superfamily_set))

    # As the dataset lacks these four types: Ngaro, VIPER, Maverick, and PiggyBac,
    # instances predicted as these categories are marked as Unknown
    filter_labels = ('Ngaro', 'VIPER', 'Maverick', 'PiggyBac')

    ## 2.3 Retrieve labels classified by RepeatClassifier; for labels that didn't
    # annotate to the superfamily or were incorrect, label them as Unknown
    # (as our dataset labels are at the superfamily level)
    names, contigs = read_fasta_v1(classified_path)
    RC_name_labels = {}
    all_unique_RM_label = set()
    for name in names:
        label = name.split('#')[1].split(' ')[0]
        if not rmToWicker.__contains__(label):
            all_unique_RM_label.add(label)
            label = 'Unknown'
        else:
            wicker_superfamily = rmToWicker[label]
            label = wicker_superfamily
            if label in filter_labels:
                label = 'Unknown'
            all_unique_RM_label.add(label)
        RC_name_labels[name.split('#')[0]] = label
    print('all_unique_RM_label:' + str(all_unique_RM_label))


    # 2.4 Obtain test data names and labels, then evaluate against labels predicted by RepeatClassifier
    names, contigs = read_fasta_v1(classified_path)
    sequence_names = []
    y_labels = []
    for name in names:
        parts = name.split('\t')
        seq_name = parts[0]
        seq_name_parts = seq_name.split(' ')
        if len(seq_name_parts) > 1:
            seq_name = seq_name_parts[0].split('#')[0]
            label = seq_name_parts[1]
        else:
            label = parts[1]
        sequence_names.append(seq_name)
        y_labels.append(label)

    y_predicts = []
    for name in sequence_names:
        y_predicts.append(RC_name_labels[name])

    print(y_labels)
    print(len(y_labels))
    print(y_predicts)
    print(len(y_predicts))
    y_test = np.array(y_labels)
    y_pred = np.array(y_predicts)
    get_metrics_by_label(y_test, y_pred)
    return RC_name_labels


def transform_TERL_data(train_path, test_path, data_dir):
    # 1. First, divide the dataset into the structure required by TERL
    train_dir = data_dir + '/Train'
    test_dir = data_dir + '/Test'
    if not os.path.exists(train_dir):
        os.makedirs(train_dir)
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)
    train_names, train_contigs = read_fasta_v1(train_path)
    test_names, test_contigs = read_fasta_v1(test_path)
    train_contigs_dict = {}
    for name in train_names:
        label = name.split('\t')[1]
        if not train_contigs_dict.__contains__(label):
            train_contigs_dict[label] = {}
        cur_train_contigs = train_contigs_dict[label]
        cur_train_contigs[name] = train_contigs[name]
    for label in train_contigs_dict.keys():
        cur_path = train_dir + '/' + label + '.fa'
        store_fasta(train_contigs_dict[label], cur_path)

    test_contigs_dict = {}
    for name in test_names:
        label = name.split('\t')[1]
        if not test_contigs_dict.__contains__(label):
            test_contigs_dict[label] = {}
        cur_test_contigs = test_contigs_dict[label]
        cur_test_contigs[name] = test_contigs[name]
    for label in test_contigs_dict.keys():
        cur_path = test_dir + '/' + label + '.fa'
        store_fasta(test_contigs_dict[label], cur_path)

# Evaluate TERL performance
def evaluate_TERL(test_path, predict_path):
    # Assess TERL accuracy on the test data
    names, contigs = read_fasta_v1(test_path)
    y_labels = []
    for name in names:
        parts = name.split('\t')
        label = parts[1]
        y_labels.append(label)

    names, contigs = read_fasta_v1(predict_path)
    y_predicts = []
    for name in names:
        parts = name.split('\t')
        label = parts[-2]
        y_predicts.append(label)

    print(y_labels)
    print(len(y_labels))
    print(y_predicts)
    print(len(y_predicts))
    y_test = np.array(y_labels)
    y_pred = np.array(y_predicts)
    get_metrics_by_label(y_test, y_pred)

# Evaluate ClassifyTE performance
def evaluate_ClassifyTE(predict_path):
    y_labels = []
    y_predicts = []
    y_predicts_set = set()
    not_superfamily_labels = ('LTR', 'SubclassI', 'LINE', 'SINE')
    with open(predict_path, 'r') as f_r:
        for i, line in enumerate(f_r):
            if i == 0:
                continue
            line = line.replace('\n', '')
            parts = line.split(',')
            raw_name = parts[0]
            y_label = raw_name.split('\t')[1]
            y_predict = parts[-1]
            y_predicts_set.add(y_predict)
            if y_predict == 'gypsy':
                y_predict = 'Gypsy'
            if y_predict in not_superfamily_labels:
                y_predict = 'Unknown'
            y_labels.append(y_label)
            y_predicts.append(y_predict)
    print(y_predicts_set)
    print(y_labels)
    print(len(y_labels))
    print(y_predicts)
    print(len(y_predicts))
    y_test = np.array(y_labels)
    y_pred = np.array(y_predicts)
    get_metrics_by_label(y_test, y_pred)

def evaluate_TEsorter(pred_path, test_path):
    label_dict = {'EnSpm_CACTA': 'CACTA', 'pararetrovirus': 'Retrovirus', 'LINE': 'Unknown',
                  'MuDR_Mutator': 'Mutator', 'mixture': 'Unknown', 'Tc1_Mariner': 'Tc1-Mariner', 'PIF_Harbinger': 'PIF-Harbinger'}
    pred_names, pred_contigs = read_fasta(pred_path)
    test_names, test_contigs = read_fasta_v1(test_path)
    test_dict = {}
    for name in test_names:
        parts = name.split('\t')
        raw_name = parts[0]
        label = parts[1]
        test_dict[raw_name] = label
    pred_labels_set = set()
    y_test = []
    y_pred = []
    for name in pred_names:
        parts = name.split('#')
        raw_name = parts[0]
        label = parts[1]
        label_parts = label.split('/')
        if len(label_parts) >= 2:
            label = label_parts[1]
        if label_dict.__contains__(label):
            label = label_dict[label]
        pred_labels_set.add(label)
        y_pred.append(label)
        test_label = test_dict[raw_name]
        y_test.append(test_label)
    print(y_test)
    print(y_pred)
    print(len(y_test))
    print(len(y_pred))
    y_test = np.array(y_test)
    y_pred = np.array(y_pred)
    get_metrics_by_label(y_test, y_pred)


# Convert Repbase to input for DeepTE training model
def transform_repbase_to_DeepTE_input(repbase_path, DeepTE_input):
    names, contigs = read_fasta_v1(repbase_path)
    DeepTE_set = set()
    with open(DeepTE_input, 'w') as f_save:
        for name in names:
            label = name.split('\t')[1]
            seq = contigs[name]
            f_save.write(label+','+seq+'\n')
    print(DeepTE_set)
    print(len(DeepTE_set))

def transform_DeepTE_to_fasta_bak(raw_dataset, repbase_dataset, outpt):
    repbase_names, repbase_contigs = read_fasta_v1(repbase_dataset)
    inverted_repbase_contigs = {value: key for key, value in repbase_contigs.items()}
    filter_contigs = []
    keep_contigs = {}
    with open(raw_dataset, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split(',')
            label = parts[0]
            seq = parts[1]
            if not inverted_repbase_contigs.__contains__(seq):
                filter_contigs.append(line)
            else:
                repbase_name = inverted_repbase_contigs[seq]
                keep_contigs[repbase_name] = seq
    store_fasta(keep_contigs, outpt)
    print('filter sequences size: ' + str(len(filter_contigs)))


def transform_DeepTE_to_fasta(raw_dataset, outpt):
    node_index = 0
    new_contigs = {}
    filter_labels = ['nLTR_LINE', 'nLTR_SINE']
    with open(raw_dataset, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split(',')
            label = parts[0]
            if label in filter_labels:
                continue
            seq = parts[1]
            new_contigs['node_'+str(node_index)+'\t'+label] = seq
            node_index += 1
    store_fasta(new_contigs, outpt)

def evaluate_DeepTE(test_path, predict_path):
    names, contigs = read_fasta_v1(test_path)
    y_labels_dict = {}
    for name in names:
        parts = name.split('\t')
        seq_name = parts[0]
        label = parts[1]
        # wicker_label = config.DeepTE_class[label]
        # y_labels_dict[seq_name] = wicker_label
        y_labels_dict[seq_name] = label

    unique_predict_labels = set()
    with open(predict_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            seq_name = parts[0]
            predict = parts[1]
            unique_predict_labels.add(predict)
    print(unique_predict_labels)

    ## 4.2 Convert DeepTE classification labels to the superfamily level; if not at the superfamily level, label as unknown
    DeepTE_labels = {'ClassII_DNA_Mutator_unknown': 'Mutator', 'ClassII_DNA_TcMar_nMITE': 'Tc1-Mariner',
                     'ClassII_DNA_hAT_unknown': 'hAT', 'ClassII_DNA_P_MITE': 'P', 'ClassI_nLTR': 'Unknown',
                     'ClassIII_Helitron': 'Helitron', 'ClassI_LTR_Gypsy': 'Gypsy', 'ClassI_LTR': 'Unknown',
                     'ClassII_DNA_Mutator_MITE': 'Mutator', 'ClassI_LTR_Copia': 'Copia', 'ClassI_nLTR_LINE': 'Unknown',
                     'ClassII_DNA_CACTA_unknown': 'CACTA', 'ClassI_nLTR_LINE_I': 'I', 'ClassI_nLTR_DIRS': 'DIRS',
                     'ClassII_MITE': 'Unknown', 'unknown': 'Unknown', 'ClassII_DNA_TcMar_unknown': 'Tc1-Mariner',
                     'ClassII_DNA_CACTA_MITE': 'CACTA', 'ClassII_DNA_Harbinger_unknown': 'PIF-Harbinger',
                     'ClassII_DNA_hAT_nMITE': 'hAT', 'ClassI': 'Unknown', 'ClassI_nLTR_SINE_7SL': '7SL',
                     'ClassII_DNA_Harbinger_nMITE': 'PIF-Harbinger', 'ClassII_DNA_Mutator_nMITE': 'Mutator',
                     'ClassII_DNA_hAT_MITE': 'hAT', 'ClassII_DNA_CACTA_nMITE': 'CACTA', 'ClassI_nLTR_SINE_tRNA': 'tRNA',
                     'ClassII_DNA_TcMar_MITE': 'Tc1-Mariner', 'ClassII_DNA_P_nMITE': 'P', 'ClassI_nLTR_PLE': 'Penelope',
                     'ClassII_DNA_Harbinger_MITE': 'PIF-Harbinger', 'ClassI_nLTR_LINE_L1': 'L1', 'ClassII_nMITE': 'Unknown',
                     'ClassI_LTR_ERV': 'Retrovirus', 'ClassI_LTR_BEL': 'Bel-Pao', 'ClassI_nLTR_LINE_RTE': 'RTE', 'ClassI_nLTR_LINE_R2': 'R2',
                     'ClassII_DNA_Transib_nMITE': 'Transib', 'ClassII_DNA_PiggyBac_nMITE': 'PiggyBac', 'ClassI_nLTR_LINE_Jockey': 'Jockey',
                     'ClassI_nLTR_SINE_5S': '5S', 'ClassII_DNA_hAT': 'hAT', 'ClassII_DNA_Tc1-Mariner': 'Tc1-Mariner',
                     'ClassII_DNA_PIF-Harbinger': 'PIF-Harbinger', 'ClassII_DNA_CACTA': 'CACTA', 'ClassII_DNA_Mutator': 'Mutator',
                     'ClassII_DNA_P': 'P', 'ClassII_DNA_Crypton': 'Crypton', 'ClassII_DNA_Transib': 'Transib', 'ClassI_LTR_Retrovirus': 'Retrovirus',
                     'ClassI_LTR_Bel-Pao': 'Bel-Pao', 'ClassII_DNA_Merlin': 'Merlin'}

    y_predict_seq_names = []
    y_predicts = []
    with open(predict_path, 'r') as f_r:
        for line in f_r:
            line = line.replace('\n', '')
            parts = line.split('\t')
            seq_name = parts[0]
            predict = parts[1]
            predict = DeepTE_labels[predict]
            y_predict_seq_names.append(seq_name)
            y_predicts.append(predict)

    y_labels = []
    for seq_name in y_predict_seq_names:
        y_labels.append(y_labels_dict[seq_name])

    print(y_labels)
    print(len(y_labels))
    print(y_predicts)
    print(len(y_predicts))
    y_test = np.array(y_labels)
    y_pred = np.array(y_predicts)
    get_metrics_by_label(y_test, y_pred)

def transform_repbase_bed(repbase_bed, NeuralTE_bed, label_dict):
    new_lines = []
    with open(repbase_bed, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            line = line.replace('\n', '')
            parts = line.split('\t')
            chr_name = parts[0]
            chr_start = int(parts[1])
            chr_end = int(parts[2])
            info = parts[3]
            info_parts = info.split(';')
            te_name = info_parts[9]
            te_class = info_parts[10]
            new_label = label_dict[te_name]
            new_string = ''
            for i, item in enumerate(info_parts):
                if i == 10:
                    new_string += new_label
                else:
                    new_string += item
                if i != len(info_parts)-1:
                    new_string += ';'
            new_line = ''
            for i, item in enumerate(parts):
                if i == 3:
                    new_line += new_string
                else:
                    new_line += item
                if i != len(parts)-1:
                    new_line += '\t'
                else:
                    new_line += '\n'
            new_lines.append(new_line)
    #print(new_lines)
    with open(NeuralTE_bed, 'w') as f_save:
        for line in new_lines:
            f_save.write(line)

def evalute_genome_coverage(mazie_merge_path, genome_path, temp_dir, threads):
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    # 1. Align sequences using RepeatMasker to the genome.
    RepeatMasker_command = 'cd ' + temp_dir + ' && RepeatMasker -e ncbi -pa ' + str(threads) \
                           + ' -q -no_is -norna -nolow -div 40 -gff -lib ' + mazie_merge_path + ' -cutoff 225 ' \
                           + genome_path
    print(RepeatMasker_command)
    os.system(RepeatMasker_command)
    out_file = genome_path + '.out'
    # 2. Convert .out file to .bed file
    convert2bed_command = 'perl ' + config.project_dir + '/tools/RMout_to_bed.pl ' + out_file + ' base1'
    print(convert2bed_command)
    os.system(convert2bed_command)

    bed_file = out_file + '.bed'
    # 3. Analyze alignments in .bed files programmatically, analyze proportions
    analyze_class_ratio(bed_file, temp_dir)

def analyze_class_ratio(bed_file, work_dir):
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    # 1. Convert .out file to an analyzable format, such as .gff
    # 2. Obtain chromosomes and start/stop positions for each element,
    # which might have overlaps; merge coordinates with overlaps,
    # resulting in non-overlapping positions for each element
    # te_chrs -> {te_name: {chr_name: []}}
    te_chrs = {}
    #class_chrs = -> {class_name: {chr_name: []}}
    class_chrs = {}
    # main_class_chrs = -> {class_name: {chr_name: []}}
    # main_class_chrs = {}

    te_class_dict = {}

    # Save the number of bases for each TE
    te_bases = {}
    # Save the number of bases for each class
    class_bases = {}
    # Save the number of bases for each main_class
    main_class_bases = {}
    with open(bed_file, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            chr_name = parts[0]
            chr_start = int(parts[1])
            chr_end = int(parts[2])
            info = parts[3]
            info_parts = info.split(';')
            te_name = info_parts[9]
            te_class = info_parts[10]
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

            # if te_class.__contains__('LTR'):
            #     main_class = 'LTR'
            # elif te_class.__contains__('Helitron'):
            #     main_class = 'Helitron'
            # elif te_class.__contains__('DNA'):
            #     main_class = 'DNA'
            # elif te_class.__contains__('MITE'):
            #     main_class = 'MITE'
            #
            # if not main_class_chrs.__contains__(main_class):
            #     main_class_chrs[main_class] = {}
            # chr_dict = main_class_chrs[main_class]
            # if not chr_dict.__contains__(chr_name):
            #     chr_dict[chr_name] = []
            # chr_pos = chr_dict[chr_name]
            # chr_pos.append((chr_start, chr_end))

            te_class_dict[te_name] = te_class

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

    # # Count the top 10 DNA transposons with the highest proportions
    # te_base_list = list(te_bases.items())
    # te_base_list.sort(key=lambda x: -x[1])
    # top_dna = []
    # for item in te_base_list:
    #     te_name = item[0]
    #     class_name = te_class_dict[te_name]
    #     if not class_name.__contains__('Helitron') and (class_name.__contains__('DNA') or class_name.__contains__('MITE')):
    #         top_dna.append((te_name, class_name, item[1]))

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

    # for class_name in main_class_chrs.keys():
    #     chr_dict = main_class_chrs[class_name]
    #     total_bases = 0
    #     for chr_name in chr_dict.keys():
    #         chr_pos = chr_dict[chr_name]
    #         merged_intervals = merge_overlapping_intervals(chr_pos)
    #         for segment in merged_intervals:
    #             total_bases += abs(segment[1] - segment[0])
    #     main_class_bases[class_name] = total_bases

    # store for testing
    te_base_file = work_dir + '/te_base.json'
    with codecs.open(te_base_file, 'w', encoding='utf-8') as f:
        json.dump(te_bases, f)
    class_base_file = work_dir + '/class_base.json'
    with codecs.open(class_base_file, 'w', encoding='utf-8') as f:
        json.dump(class_bases, f)
    # main_class_base_file = work_dir + '/main_class_base.json'
    # with codecs.open(main_class_base_file, 'w', encoding='utf-8') as f:
    #     json.dump(main_class_bases, f)
    # top_dna_file = work_dir + '/top_dna.json'
    # with codecs.open(top_dna_file, 'w', encoding='utf-8') as f:
    #     json.dump(top_dna, f)


def analyze_class_ratio_gff(gff_file, work_dir, total_genome_len):
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    # 1. Convert .out file to an analyzable format, such as .gff
    # 2. Obtain chromosomes and start/stop positions for each element,
    # which might have overlaps; merge coordinates with overlaps,
    # resulting in non-overlapping positions for each element
    # te_chrs -> {te_name: {chr_name: []}}
    te_chrs = {}
    #class_chrs = -> {class_name: {chr_name: []}}
    class_chrs = {}
    # main_class_chrs = -> {class_name: {chr_name: []}}
    # main_class_chrs = {}

    te_class_dict = {}

    # Save the number of bases for each TE
    te_bases = {}
    # Save the number of bases for each class
    class_bases = {}
    # Save the number of bases for each main_class
    class_ratios = {}
    # Count the top 10 DNA transposons with the highest proportions
    main_class_bases = {}
    with open(gff_file, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            chr_name = parts[0]
            chr_start = int(parts[3])
            chr_end = int(parts[4])
            te_name = parts[8].split(' ')[1].replace('\"', '').split(':')[1]
            te_class = parts[2]
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

            # if te_class.__contains__('LTR'):
            #     main_class = 'LTR'
            # elif te_class.__contains__('Helitron'):
            #     main_class = 'Helitron'
            # elif te_class.__contains__('DNA'):
            #     main_class = 'DNA'
            # elif te_class.__contains__('MITE'):
            #     main_class = 'MITE'
            #
            # if not main_class_chrs.__contains__(main_class):
            #     main_class_chrs[main_class] = {}
            # chr_dict = main_class_chrs[main_class]
            # if not chr_dict.__contains__(chr_name):
            #     chr_dict[chr_name] = []
            # chr_pos = chr_dict[chr_name]
            # chr_pos.append((chr_start, chr_end))

            te_class_dict[te_name] = te_class

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

    # # Count the top 10 DNA transposons with the highest proportions
    # te_base_list = list(te_bases.items())
    # te_base_list.sort(key=lambda x: -x[1])
    # top_dna = []
    # for item in te_base_list:
    #     te_name = item[0]
    #     class_name = te_class_dict[te_name]
    #     if not class_name.__contains__('Helitron') and (class_name.__contains__('DNA') or class_name.__contains__('MITE')):
    #         top_dna.append((te_name, class_name, item[1]))

    for class_name in class_chrs.keys():
        chr_dict = class_chrs[class_name]
        total_bases = 0
        for chr_name in chr_dict.keys():
            chr_pos = chr_dict[chr_name]
            merged_intervals = merge_overlapping_intervals(chr_pos)
            for segment in merged_intervals:
                total_bases += abs(segment[1] - segment[0])
        class_bases[class_name] = total_bases
        class_ratios[class_name] = round(float(total_bases)/total_genome_len, 4)
    #print(class_bases)

    # for class_name in main_class_chrs.keys():
    #     chr_dict = main_class_chrs[class_name]
    #     total_bases = 0
    #     for chr_name in chr_dict.keys():
    #         chr_pos = chr_dict[chr_name]
    #         merged_intervals = merge_overlapping_intervals(chr_pos)
    #         for segment in merged_intervals:
    #             total_bases += abs(segment[1] - segment[0])
    #     main_class_bases[class_name] = total_bases

    # store for testing
    te_base_file = work_dir + '/te_base.json'
    with codecs.open(te_base_file, 'w', encoding='utf-8') as f:
        json.dump(te_bases, f)
    class_base_file = work_dir + '/class_base.json'
    with codecs.open(class_base_file, 'w', encoding='utf-8') as f:
        json.dump(class_bases, f)
    class_ratios_file = work_dir + '/class_ratios.json'
    with codecs.open(class_ratios_file, 'w', encoding='utf-8') as f:
        json.dump(class_ratios, f)
    # main_class_base_file = work_dir + '/main_class_base.json'
    # with codecs.open(main_class_base_file, 'w', encoding='utf-8') as f:
    #     json.dump(main_class_bases, f)
    # top_dna_file = work_dir + '/top_dna.json'
    # with codecs.open(top_dna_file, 'w', encoding='utf-8') as f:
    #     json.dump(top_dna, f)

def merge_overlapping_intervals(intervals):
    # Sort genomic position arrays based on the starting position
    sorted_intervals = sorted(intervals, key=lambda x: x[0])

    merged_intervals = []
    current_interval = sorted_intervals[0]

    # Traverse the sorted genomic position array
    for interval in sorted_intervals[1:]:
        # If the current genomic position overlaps with the next one
        if current_interval[1] >= interval[0]:
            # Update the end position of the current genomic position
            current_interval = (current_interval[0], max(current_interval[1], interval[1]))
        else:
            # If there's no overlap, add the current genomic position to the result array and update the current genomic position to the next one
            merged_intervals.append(current_interval)
            current_interval = interval

    # Add the last genomic position to the result array
    merged_intervals.append(current_interval)

    return merged_intervals

def get_train_except_species(total_repbase, species_name, train_path):
    names, contigs = read_fasta_v1(total_repbase)
    train_contigs = {}
    for name in names:
        parts = name.split('\t')
        cur_species_name = parts[2]
        if cur_species_name != species_name:
            train_contigs[name] = contigs[name]
    store_fasta(train_contigs, train_path)

# Create a 3D plot to visualize changes in parameters affecting the f1-score
def plot_3D_param(x, y, z, work_dir):
    print(x.shape, y.shape, z.shape)

    xi = np.linspace(min(x), max(x))
    yi = np.linspace(min(y), max(y))
    X, Y = np.meshgrid(xi, yi)

    Z = griddata((x, y), z, (X, Y), method='cubic')

    fig = plt.figure()
    ax3 = plt.axes(projection='3d')

    #plt.rcParams['font.sans-serif'] = ['FangSong']
    plt.rcParams['axes.unicode_minus'] = False

    surf=ax3.plot_surface(X, Y, Z, cmap='BuPu', linewidth=0, antialiased=False)
    ax3.set_xlabel('internal')
    ax3.set_ylabel('terminal')
    ax3.set_zlabel('F1-score')
    colorbar = fig.colorbar(surf, shrink=0.6)

    # 调整颜色条的位置
    position = ax3.get_position()
    colorbar_position = [position.x0 + position.width + 0.1, position.y0, 0.02, position.height]
    colorbar.ax.set_position(colorbar_position)

    #plt.show()
    plt.savefig(work_dir + '/output.png', dpi=500)

def generate_TERL_dataset(fasta_file, outdir):
    os.makedirs(outdir, exist_ok=True)
    contignames, contigs = read_fasta_v1(fasta_file)
    label_seqs = {}
    for name in contignames:
        label = name.split('\t')[1]
        if not label_seqs.__contains__(label):
            label_seqs[label] = {}
        cur_label_seqs = label_seqs[label]
        cur_label_seqs[name] = contigs[name]
    for label in label_seqs.keys():
        cur_label_path = outdir + '/' + label + '.fa'
        store_fasta(label_seqs[label], cur_label_path)

def generate_ClassifyTE_dataset(fasta_file):
    filter_path = filterRepbase(fasta_file, config.ClassifyTE_class)