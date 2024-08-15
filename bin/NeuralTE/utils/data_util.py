#-- coding: UTF-8 --
import os
import random
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
# import pandas as pd
from openpyxl.utils import get_column_letter
# from pandas import ExcelWriter
import numpy as np
import itertools
from configs import config
from configs import gpu_config
import tensorflow as tf


# Calculate the proportion of each category in the repbase file.
def summary_class_ratio(repbase_path):
    train_names, train_contigs = read_fasta_v1(repbase_path)
    train_class_num = {}
    train_class_ratio = {}
    train_class_ratio_reverse = {}
    class_set = set()
    species_set = set()
    for name in train_names:
        class_name = name.split('\t')[1]
        species_name = name.split('\t')[2]
        if not train_class_num.__contains__(class_name):
            train_class_num[class_name] = 0
        class_num = train_class_num[class_name]
        train_class_num[class_name] = class_num + 1
        species_set.add(species_name)
        class_set.add(class_name)
    print(train_class_num)
    for class_name in train_class_num.keys():
        class_num = train_class_num[class_name]
        ratio = round(100 * float(class_num) / len(train_names), 2)
        reverse_ratio = round(1/ratio, 2)
        train_class_ratio[class_name] = ratio
        train_class_ratio_reverse[class_name] = reverse_ratio
    print(train_class_ratio)
    #print(train_class_ratio_reverse)
    print(len(species_set))
    return train_class_num, class_set

# Divide the repbase data into plant and non-plant categories.
def split_repbase_by_plant(repbase_path, repbase_plant_path, repbase_non_plant_path, ncbi_ref_info):
    species_info = {}
    with open(ncbi_ref_info, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            line = line.replace('\n', '')
            parts = line.split('\t')
            species_name = parts[2]
            is_plant = int(parts[5])
            species_info[species_name] = is_plant

    plant_contigs = {}
    non_plant_contigs = {}
    names, contigs = read_fasta_v1(repbase_path)
    for name in names:
        species_name = name.split('\t')[2]
        if species_name not in species_info:
            continue
        else:
            is_plant = species_info[species_name]
            if is_plant:
                plant_contigs[name] = contigs[name]
            else:
                non_plant_contigs[name] = contigs[name]
    store_fasta(plant_contigs, repbase_plant_path)
    store_fasta(non_plant_contigs, repbase_non_plant_path)

def parse_embl(embl_filename):
    sequences = {}
    current_sequence = None
    sequence_section = False
    annotation_section = False
    with open(embl_filename, "r") as embl_file:
        for line in embl_file:
            line = line.strip()

            if line.startswith("ID"):
                current_sequence = line.split()[1].split(';')[0].strip()

            if current_sequence:
                if line.startswith("CC") and "RepeatMasker Annotations:" in line:
                    annotation_section = True
                    continue

                if annotation_section and line.startswith("CC"):
                    type_parts = line.split('Type:')
                    sub_type_parts = line.split('SubType:')
                    species_parts = line.split('Species:')

                    if len(sub_type_parts) == 2:
                        subtype = sub_type_parts[1].strip()
                        if len(subtype) > 0:
                            current_sequence += '/' + subtype
                    elif len(type_parts) == 2:
                        type = type_parts[1].strip()
                        if len(type) > 0:
                            current_sequence += '#' + type
                    elif len(species_parts) == 2:
                        species = species_parts[1].strip()
                        if len(species) > 0:
                            current_sequence += '\t' + species

                if annotation_section and line.startswith("XX"):
                    annotation_section = False
                    sequences[current_sequence] = ""

                if line.startswith("SQ"):
                    sequence_section = True
                elif sequence_section and line.startswith("//"):
                    sequence_section = False
                elif sequence_section and line:
                    sequence_line = "".join(line.split(" ")[:-1])
                    sequences[current_sequence] += sequence_line

    return sequences

# Extract fasta files from the Dfam embl file.
def extract_fasta_from_embl(dfam_embl, output_fasta):
    sequence_dict = parse_embl(dfam_embl)
    store_fasta(sequence_dict, output_fasta)

def generate_random_sequence(length):
    bases = ['A', 'T', 'C', 'G']
    sequence = ''.join(random.choice(bases) for _ in range(length))
    return sequence

# Function to randomly generate nucleotide sequences.
def generate_random_sequences(num_sequences):
    sequence_lengths = []

    # Generate sequences ranging from 0 to 600 in length.
    num_sequences_range1 = num_sequences // 4
    sequence_lengths.extend(random.randint(80, 600) for _ in range(num_sequences_range1))

    # Generate sequences ranging from 601 to 1800 in length.
    num_sequences_range2 = num_sequences // 4
    sequence_lengths.extend(random.randint(601, 1800) for _ in range(num_sequences_range2))

    # Generate sequences ranging from 1801 to 4000 in length.
    num_sequences_range3 = num_sequences // 4
    sequence_lengths.extend(random.randint(1801, 4000) for _ in range(num_sequences_range3))

    # Generate sequences ranging from 4000 to 20000 in length.
    num_sequences_range4 = num_sequences // 4
    sequence_lengths.extend(random.randint(4000, 20000) for _ in range(num_sequences_range4))

    return sequence_lengths

def merge_tsd_terminal_repbase(tsd_repbase, terminal_repbase, merge_repbase, total_domain, merge_domain):
    # The purpose of this function: Usually obtaining TSD (Target Site Duplication) and
    # terminal sequences are two separate actions, resulting in two separate repbase files.
    # We aim to retrieve terminal information based on the TSD sequence, specifically by
    # using seq_name, from the terminal file. If terminal information is absent, it will be
    # left empty.
    tsd_names, tsd_contigs = read_fasta_v1(tsd_repbase)
    term_names, term_contigs = read_fasta_v1(terminal_repbase)
    term_name_set = set()
    seq_name_to_term_name = {}
    merge_contigs = {}
    for term_name in term_names:
        parts = term_name.split('\t')
        seq_name = parts[0]
        LTR_info = parts[3]
        TIR_info = parts[4]
        term_name_set.add(seq_name)
        seq_name_to_term_name[seq_name] = (LTR_info, TIR_info)
    tsd_name_set = set()
    for tsd_name in tsd_names:
        seq_name = tsd_name.split('\t')[0]
        tsd_name_set.add(seq_name)
        if seq_name in term_name_set:
            LTR_info, TIR_info = seq_name_to_term_name[seq_name]
        else:
            LTR_info = 'LTR:'
            TIR_info = 'TIR:'
        common_name = tsd_name + '\t' + LTR_info + '\t' + TIR_info
        merge_contigs[common_name] = tsd_contigs[tsd_name]
    store_fasta(merge_contigs, merge_repbase)

    # Extract domain information from the overall domain file.
    domain_header = ''
    merge_domain_lines = []
    with open(total_domain, 'r') as f_r:
        for i, line in enumerate(f_r):
            if i < 2:
                domain_header += line
            else:
                seq_name = line.split('\t')[0]
                if seq_name in tsd_name_set:
                    merge_domain_lines.append(line)
    with open(merge_domain, 'w') as f_save:
        f_save.write(domain_header)
        for line in merge_domain_lines:
            f_save.write(line)



##word_seq generates eg. ['AA', 'AT', 'TC', 'CG', 'GT']
def word_seq(seq, k, stride=1):
    i = 0
    words_list = []
    while i <= len(seq) - k:
        words_list.append(seq[i: i + k])
        i += stride
    return (words_list)

def generate_kmer_dic(repeat_num):
    ##initiate a dic to store the kmer dic
    ##kmer_dic = {'ATC':0,'TTC':1,...}
    kmer_dic = {}
    bases = ['A','G','C','T']
    kmer_list = list(itertools.product(bases, repeat=int(repeat_num)))
    for eachitem in kmer_list:
        #print(eachitem)
        each_kmer = ''.join(eachitem)
        kmer_dic[each_kmer] = 0

    return (kmer_dic)

def generate_mat(words_list,kmer_dic):
    for eachword in words_list:
        kmer_dic[eachword] += 1
    num_list = []  ##this dic stores num_dic = [0,1,1,0,3,4,5,8,2...]
    for eachkmer in kmer_dic:
        num_list.append(kmer_dic[eachkmer])
    return (num_list)

def get_batch_kmer_freq_v1(grouped_x, internal_kmer_sizes, terminal_kmer_sizes, minority_labels_class, all_wicker_class):
    # Build vocabulary and reverse vocabulary for nucleotide sequences.
    base_vocab = {'A': 1, 'T': 2, 'C': 3, 'G': 4}

    group_dict = {}
    for x in grouped_x:
        # Split sequences into internal_Seq, LTR, and TIR segments.
        seq_name = x[0]
        seq = x[1]
        TSD_seq = x[2]
        TSD_len = x[3]
        LTR_pos = x[4]
        TIR_pos = x[5]
        domain_label_set = x[6]
        # minority_label_set = x[7]
        # LTR_pos = x[2]
        # TIR_pos = x[3]
        # domain_label_set = x[4]
        internal_seq = ''
        LTR_seq = ''
        TIR_seq = ''
        LTR_pos_str = str(LTR_pos.split(':')[1]).strip()
        TIR_pos_str = str(TIR_pos.split(':')[1]).strip()
        if LTR_pos_str == '' and TIR_pos_str == '':
            internal_seq = seq
        if TIR_pos_str != '':
            TIR_parts = TIR_pos_str.split(',')
            left_TIR_start = int(TIR_parts[0].split('-')[0])
            left_TIR_end = int(TIR_parts[0].split('-')[1])
            right_TIR_start = int(TIR_parts[1].split('-')[0])
            right_TIR_end = int(TIR_parts[1].split('-')[1])
            TIR_seq = seq[left_TIR_start-1: left_TIR_end]
            internal_seq = seq[left_TIR_end: right_TIR_start-1]

            # Attempt to extract only the first 60 bp of the TIR sequence.
            TIR_seq = TIR_seq[0: 60]
        if LTR_pos_str != '':
            LTR_parts = LTR_pos_str.split(',')
            left_LTR_start = int(LTR_parts[0].split('-')[0])
            left_LTR_end = int(LTR_parts[0].split('-')[1])
            right_LTR_start = int(LTR_parts[1].split('-')[0])
            right_LTR_end = int(LTR_parts[1].split('-')[1])
            LTR_seq = seq[left_LTR_start-1: left_LTR_end]
            internal_seq = seq[left_LTR_end: right_LTR_start-1]

        # Represent internal_seq and LTR as k-mer frequencies and positional information.
        connected_num_list = []
        if config.use_kmers:
            for kmer_size in internal_kmer_sizes:
                if config.use_terminal:
                    # Represent internal_seq as k-mer frequencies.
                    words_list = word_seq(internal_seq, kmer_size, stride=1)
                    kmer_dic = generate_kmer_dic(kmer_size)
                    num_list = generate_mat(words_list, kmer_dic)
                    connected_num_list += num_list
                else:
                    # Represent seq as k-mer frequencies.
                    words_list = word_seq(seq, kmer_size, stride=1)
                    kmer_dic = generate_kmer_dic(kmer_size)
                    num_list = generate_mat(words_list, kmer_dic)
                    connected_num_list += num_list

            for kmer_size in terminal_kmer_sizes:
                if config.use_terminal:
                    # Represent LTR and TIR as k-mer frequencies.
                    words_list = word_seq(LTR_seq, kmer_size, stride=1)
                    kmer_dic = generate_kmer_dic(kmer_size)
                    num_list = generate_mat(words_list, kmer_dic)
                    connected_num_list += num_list

                    words_list = word_seq(TIR_seq, kmer_size, stride=1)
                    kmer_dic = generate_kmer_dic(kmer_size)
                    num_list = generate_mat(words_list, kmer_dic)
                    connected_num_list += num_list

        if config.use_TSD:
            # Convert TSD into one-hot encoding.
            encoder = np.eye(4, dtype=np.int8)
            max_length = config.max_tsd_length
            encoded_TSD = np.zeros((max_length, 4), dtype=np.int8)
            if TSD_seq == 'Unknown' or 'N' in TSD_seq:
                TSD_seq = ''
                padding_length = max_length - len(TSD_seq)
                # TSD_len
                connected_num_list.append(max_length+1)
                for i in range(padding_length):
                    encoded_TSD[len(TSD_seq) + i] = np.ones(4)
            else:
                # TSD_len
                connected_num_list.append(TSD_len)
                padding_length = max_length - len(TSD_seq)
                for i, base in enumerate(TSD_seq):
                    if base == 'A':
                        encoded_TSD[i] = encoder[0]
                    elif base == 'T':
                        encoded_TSD[i] = encoder[1]
                    elif base == 'C':
                        encoded_TSD[i] = encoder[2]
                    elif base == 'G':
                        encoded_TSD[i] = encoder[3]
                for i in range(padding_length):
                    encoded_TSD[len(TSD_seq) + i] = np.zeros(4)
            onehot_encoded_flat = encoded_TSD.reshape(-1)
            connected_num_list = np.concatenate((connected_num_list, onehot_encoded_flat))

        if config.use_ends:
            # Take the first and last 5 bp of a sequence to form a 10bp vector.
            end_seq = seq[:5] + seq[-5:]
            # Convert end_seq into one-hot encoding.
            encoder = np.eye(4, dtype=np.int8)
            max_length = 10
            encoded_end_seq = np.zeros((max_length, 4), dtype=np.int8)
            for i, base in enumerate(end_seq):
                if base == 'A':
                    encoded_end_seq[i] = encoder[0]
                elif base == 'T':
                    encoded_end_seq[i] = encoder[1]
                elif base == 'C':
                    encoded_end_seq[i] = encoder[2]
                elif base == 'G':
                    encoded_end_seq[i] = encoder[3]
            onehot_encoded_flat = encoded_end_seq.reshape(-1)
            connected_num_list = np.concatenate((connected_num_list, onehot_encoded_flat))

        if config.use_domain:
            # Convert domain set into one-hot encoding.
            encoder = [0] * len(all_wicker_class)
            for domain_label in domain_label_set:
                domain_label_num = all_wicker_class[domain_label]
                encoder[domain_label_num] = 1
            connected_num_list = np.concatenate((connected_num_list, encoder))

        group_dict[seq_name] = connected_num_list
    return group_dict

def get_gpu_config(start_gpu_num, use_gpu_num):
    gpus = tf.config.experimental.list_physical_devices('GPU')
    # For GPU memory growth
    for device in gpus:
        tf.config.experimental.set_memory_growth(device, True)
    use_devices = gpu_config.all_devices[start_gpu_num: start_gpu_num + use_gpu_num]
    tf.config.experimental.set_visible_devices(gpus[start_gpu_num: start_gpu_num + use_gpu_num], 'GPU')
    # Create a MirroredStrategy to use multiple GPUs
    strategy = tf.distribute.MirroredStrategy(devices=use_devices)
    # print('Number of devices: {}'.format(strategy.num_replicas_in_sync))

def get_feature_len():
    # Obtain CNN input dimensions
    X_feature_len = 0
    # Dimensions of TE terminal and internal sequences
    if config.use_kmers != 0:
        for kmer_size in config.internal_kmer_sizes:
            X_feature_len += pow(4, kmer_size)
        if config.use_terminal != 0:
            for i in range(2):
                for kmer_size in config.terminal_kmer_sizes:
                    X_feature_len += pow(4, kmer_size)
    if config.use_TSD != 0:
        X_feature_len += config.max_tsd_length * 4 + 1
    # if config.use_minority != 0:
    #     X_feature_len += len(config.minority_labels_class)
    if config.use_domain != 0:
        X_feature_len += len(config.all_wicker_class)
    if config.use_ends != 0:
        X_feature_len += 10 * 4
    return X_feature_len

def split_list_into_groups(lst, group_size):
    return [lst[i:i+group_size] for i in range(0, len(lst), group_size)]

def generate_feature_mats(X, Y, seq_names, minority_labels_class, all_wicker_class, internal_kmer_sizes, terminal_kmer_sizes, ex):
    seq_mats = {}
    jobs = []
    grouped_X = split_list_into_groups(X, 100)

    for grouped_x in grouped_X:
        job = ex.submit(get_batch_kmer_freq_v1, grouped_x, internal_kmer_sizes, terminal_kmer_sizes, minority_labels_class, all_wicker_class)
        jobs.append(job)
    ex.shutdown(wait=True)

    for job in as_completed(jobs):
        cur_group_dict = job.result()
        seq_mats.update(cur_group_dict)

    final_X = []
    final_Y = []
    for item in seq_names:
        seq_name = item[0]
        x = seq_mats[seq_name]
        final_X.append(x)
        label = Y[seq_name]
        label_num = all_wicker_class[label]
        final_Y.append(label_num)
    return np.array(final_X), np.array(final_Y)

def replace_non_atcg(sequence):
    return re.sub("[^ATCG]", "", sequence)

def getRMToWicker(RM_Wicker_struct):
    # 3.2 Convert Dfam classification names into Wicker format.
    ## 3.2.1 This file contains the conversion between RepeatMasker category, Repbase, and Wicker category.
    rmToWicker = {}
    wicker_superfamily_set = set()
    with open(RM_Wicker_struct, 'r') as f_r:
        for i, line in enumerate(f_r):
            parts = line.split('\t')
            rm_type = parts[5]
            rm_subtype = parts[6]
            repbase_type = parts[7]
            wicker_type = parts[8]
            wicker_type_parts = wicker_type.split('/')
            # print(rm_type + ',' + rm_subtype + ',' + repbase_type + ',' + wicker_type)
            # if len(wicker_type_parts) != 3:
            #     continue
            wicker_superfamily_parts = wicker_type_parts[-1].strip().split(' ')
            if len(wicker_superfamily_parts) == 1:
                wicker_superfamily = wicker_superfamily_parts[0]
            elif len(wicker_superfamily_parts) > 1:
                wicker_superfamily = wicker_superfamily_parts[1].replace('(', '').replace(')', '')
            rm_full_type = rm_type + '/' + rm_subtype
            if wicker_superfamily == 'ERV':
                wicker_superfamily = 'Retrovirus'
            rmToWicker[rm_full_type] = wicker_superfamily
            wicker_superfamily_set.add(wicker_superfamily)
    # Supplement some elements.
    rmToWicker['LINE/R2'] = 'R2'
    rmToWicker['LINE/RTE'] = 'RTE'
    rmToWicker['LTR/ERVL'] = 'Retrovirus'
    rmToWicker['LTR/Ngaro'] = 'DIRS'
    return rmToWicker

def transfer_RMOut2Bed(RMOut, out_bed, consensus_path, tools_dir, coverage_threshold, name_label_dict):
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
            # query_name = 'Chr'+parts[0]
            query_name = 'Chr'+str(int(parts[0].split('_')[1]) + 1)
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
                new_line = query_name+'\t'+q_start+'\t'+q_end+'\t'+name_label_dict[subject_name]+'\n'
                lines.append(new_line)
    with open(out_bed, 'w') as f_save:
        for line in lines:
            f_save.write(line)

def transfer_RMOut2BlastnOut(RMOut, BlastnOut, tools_dir):
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
            new_line = query_name+'\t'+subject_name+'\t'+'-1'+'\t'+'-1'+'\t'+'-1'+'\t'+'-1'+'\t'+q_start+'\t'+q_end+'\t'+s_start+'\t'+s_end+'\t'+'-1'+'\t'+'-1'+'\n'
            lines.append(new_line)
    with open(BlastnOut, 'w') as f_save:
        for line in lines:
            f_save.write(line)

def load_repbase_with_TSD(path, domain_path, minority_train_path, minority_out, all_wicker_class, RM_Wicker_struct):
    rmToWicker = getRMToWicker(RM_Wicker_struct)
    domain_name_labels = {}
    if config.use_domain == 1 and os.path.exists(domain_path):
        # Load the domain file and read the TE-contained domain labels.
        with open(domain_path, 'r') as f_r:
            for i, line in enumerate(f_r):
                if i < 2:
                    continue
                parts = line.split('\t')
                TE_name = parts[0]
                label = parts[1].split('#')[1]
                if not rmToWicker.__contains__(label):
                    label = 'Unknown'
                else:
                    wicker_superfamily = rmToWicker[label]
                    label = wicker_superfamily
                    if not all_wicker_class.__contains__(label):
                        label = 'Unknown'
                if not domain_name_labels.__contains__(TE_name):
                    domain_name_labels[TE_name] = set()
                label_set = domain_name_labels[TE_name]
                label_set.add(label)

    names, contigs = read_fasta_v1(path)
    X = []
    Y = {}
    seq_names = []
    for name in names:
        feature_info = {}
        parts = name.split('\t')
        seq_name = parts[0]
        label = parts[1]
        species_name = parts[2]
        for p_name in parts:
            if 'TSD:' in p_name:
                TSD_seq = p_name.split(':')[1]
                feature_info['TSD_seq'] = TSD_seq
            elif 'TSD_len:' in p_name:
                tsd_len_str = p_name.split(':')[1]
                if tsd_len_str == '':
                    TSD_len = 0
                else:
                    TSD_len = int(tsd_len_str)
                feature_info['TSD_len'] = TSD_len
            elif 'LTR:' in p_name:
                LTR_info = p_name
                feature_info['LTR_info'] = LTR_info
            elif 'TIR:' in p_name:
                TIR_info = p_name
                feature_info['TIR_info'] = TIR_info
        if config.use_TSD:
            TSD_seq = feature_info['TSD_seq']
            TSD_len = feature_info['TSD_len']
        else:
            TSD_seq = ''
            TSD_len = 0

        if config.use_terminal:
            LTR_info = feature_info['LTR_info']
            TIR_info = feature_info['TIR_info']
        else:
            LTR_info = 'LTR:'
            TIR_info = 'TIR:'

        if seq_name.endswith('-RC'):
            raw_seq_name = seq_name[:-3]
        else:
            raw_seq_name = seq_name
        if domain_name_labels.__contains__(raw_seq_name):
            domain_label_set = domain_name_labels[raw_seq_name]
        else:
            domain_label_set = {'Unknown'}

        seq = contigs[name]
        seq = replace_non_atcg(seq)  # undetermined nucleotides in splice
        x_feature = (seq_name, seq, TSD_seq, TSD_len, LTR_info, TIR_info, domain_label_set)
        X.append(x_feature)
        Y[seq_name] = label
        seq_names.append((seq_name, label))
    return X, Y, seq_names

def generate_non_autonomous_data(total_repbase, out_path):
    pattern = r'[^-_]+-\d*N\d*[a-zA-Z]*_[^-_]+$'
    names, contigs = read_fasta_v1(total_repbase)
    DNA_labels = ['Tc1-Mariner', 'hAT', 'Mutator', 'Merlin', 'Transib', 'P', 'PiggyBac', 'PIF-Harbinger', 'CACTA']
    match_names = []
    out_contigs = {}
    for name in names:
        parts = name.split('\t')
        seq_name = parts[0]
        label = parts[1]
        if label in DNA_labels and re.match(pattern, seq_name):
            match_names.append(name)
            out_contigs[name] = contigs[name]
    print(len(match_names))
    store_fasta(out_contigs, out_path)

def generate_only_dna_data(total_repbase, out_path):
    names, contigs = read_fasta_v1(total_repbase)
    DNA_labels = ['Tc1-Mariner', 'hAT', 'Mutator', 'Merlin', 'Transib', 'P', 'PiggyBac', 'PIF-Harbinger', 'CACTA']
    match_names = []
    out_contigs = {}
    for name in names:
        parts = name.split('\t')
        seq_name = parts[0]
        label = parts[1]
        if label in DNA_labels:
            match_names.append(name)
            out_contigs[name] = contigs[name]
    print(len(match_names))
    store_fasta(out_contigs, out_path)

def split_fasta(cur_path, output_dir, num_chunks):
    split_files = []

    if os.path.exists(output_dir):
        os.system('rm -rf ' + output_dir)
    os.makedirs(output_dir)

    names, contigs = read_fasta_v1(cur_path)
    num_names = len(names)
    chunk_size = num_names // num_chunks

    for i in range(num_chunks):
        chunk_start = i * chunk_size
        chunk_end = chunk_start + chunk_size if i < num_chunks - 1 else num_names
        chunk = names[chunk_start:chunk_end]
        output_path = output_dir + '/out_' + str(i) + '.fa'
        with open(output_path, 'w') as out_file:
            for name in chunk:
                seq = contigs[name]
                out_file.write('>'+name+'\n'+seq+'\n')
        split_files.append(output_path)
    return split_files

def run_command(command):
    subprocess.run(command, check=True, shell=True)

def identify_terminals(split_file, output_dir, tool_dir):
    try:
        ltrsearch_command = 'cd ' + output_dir + ' && ' + tool_dir + '/ltrsearch -l 50 ' + split_file + ' > /dev/null 2>&1'
        itrsearch_command = 'cd ' + output_dir + ' && ' + tool_dir + '/itrsearch -i 0.7 -l 7 ' + split_file+ ' > /dev/null 2>&1'
        run_command(ltrsearch_command)
        run_command(itrsearch_command)
        # os.system(ltrsearch_command)
        # os.system(itrsearch_command)
        ltr_file = split_file + '.ltr'
        tir_file = split_file + '.itr'

        # Read ltr and itr files to get the start and end positions of ltr and itr.
        ltr_names, ltr_contigs = read_fasta_v1(ltr_file)
        tir_names, tir_contigs = read_fasta_v1(tir_file)
        LTR_info = {}
        for ltr_name in ltr_names:
            parts = ltr_name.split('\t')
            orig_name = parts[0]
            terminal_info = parts[-1]
            LTR_info_parts = terminal_info.split('LTR')[1].split(' ')[0].replace('(', '').replace(')', '').split('..')
            LTR_left_pos_parts = LTR_info_parts[0].split(',')
            LTR_right_pos_parts = LTR_info_parts[1].split(',')
            lLTR_start = int(LTR_left_pos_parts[0])
            lLTR_end = int(LTR_left_pos_parts[1])
            rLTR_start = int(LTR_right_pos_parts[1])
            rLTR_end = int(LTR_right_pos_parts[0])
            LTR_info[orig_name] = (lLTR_start, lLTR_end, rLTR_start, rLTR_end)
        TIR_info = {}
        for tir_name in tir_names:
            parts = tir_name.split('\t')
            orig_name = parts[0]
            terminal_info = parts[-1]
            TIR_info_parts = terminal_info.split('ITR')[1].split(' ')[0].replace('(', '').replace(')', '').split('..')
            TIR_left_pos_parts = TIR_info_parts[0].split(',')
            TIR_right_pos_parts = TIR_info_parts[1].split(',')
            lTIR_start = int(TIR_left_pos_parts[0])
            lTIR_end = int(TIR_left_pos_parts[1])
            rTIR_start = int(TIR_right_pos_parts[1])
            rTIR_end = int(TIR_right_pos_parts[0])
            TIR_info[orig_name] = (lTIR_start, lTIR_end, rTIR_start, rTIR_end)

        # Update the header of the split_file, adding two columns LTR:1-206,4552-4757 TIR:1-33,3869-3836.
        update_split_file = split_file + '.updated'
        update_contigs = {}
        names, contigs = read_fasta_v1(split_file)
        for name in names:
            orig_name = name.split('\t')[0]
            LTR_str = 'LTR:'
            if LTR_info.__contains__(orig_name):
                lLTR_start, lLTR_end, rLTR_start, rLTR_end = LTR_info[orig_name]
                LTR_str += str(lLTR_start) + '-' + str(lLTR_end) + ',' + str(rLTR_start) + '-' + str(rLTR_end)
            TIR_str = 'TIR:'
            if TIR_info.__contains__(orig_name):
                lTIR_start, lTIR_end, rTIR_start, rTIR_end = TIR_info[orig_name]
                TIR_str += str(lTIR_start) + '-' + str(lTIR_end) + ',' + str(rTIR_start) + '-' + str(rTIR_end)
            update_name = name + '\t' + LTR_str + '\t' + TIR_str
            update_contigs[update_name] = contigs[name]
        store_fasta(update_contigs, update_split_file)
        return update_split_file
    except Exception as e:
        return e



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

    # Get all LTR sequences.
    LTR_names = set()
    for name in raw_names:
        if name.__contains__('Non_LTR'):
            continue
        # Identify LTR terminal sequences and obtain corresponding internal sequence names.
        pattern = r'\b(\w+(-|_)?)LTR((-|_)?\w*)\b'
        matches = re.findall(pattern, name)
        if matches:
            replacement = r'\1I\3'
            internal_name1 = re.sub(pattern, replacement, name)
            replacement = r'\1INT\3'
            internal_name2 = re.sub(pattern, replacement, name)
            LTR_names.add(name)
            LTR_names.add(internal_name1)
            LTR_names.add(internal_name2)

    # Store segmented LTR-to-complete LTR correspondences.
    SegLTR2intactLTR = {}
    new_names = []
    new_contigs = {}
    for name in raw_names:
        if name in LTR_names:
            # If the current sequence is LTR, check if the internal sequence exists.
            pattern = r'\b(\w+(-|_)?)LTR((-|_)?\w*)\b'
            matches = re.findall(pattern, name)
            if matches:
                ltr_name = name
                replacement = r'\1I\3'
                internal_name1 = re.sub(pattern, replacement, name)
                replacement = r'\1INT\3'
                internal_name2 = re.sub(pattern, replacement, name)

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
                        replacement = r'\1intactLTR\3'
                        intact_ltr_name = re.sub(pattern, replacement, name)
                        intact_ltr_seq = raw_contigs[ltr_name] + internal_seq + raw_contigs[ltr_name]
                        new_names.append(intact_ltr_name)
                        new_contigs[intact_ltr_name] = intact_ltr_seq
                        repbase_labels[intact_ltr_name] = repbase_labels[ltr_name]
                        SegLTR2intactLTR[ltr_name] = intact_ltr_name
                        SegLTR2intactLTR[internal_name] = intact_ltr_name
        else:
            # If the current sequence is INT, discard it directly because an INT with an LTR will surely be recognized,
            # while an INT without an LTR should be treated as an incomplete LTR and discarded.
            pattern = r'\b(\w+(-|_)?)INT((-|_)?\w*)\b'
            matches = re.findall(pattern, name)
            if not matches:
                # Retain other types of transposons.
                new_names.append(name)
                new_contigs[name] = raw_contigs[name]

    # Step4. store Repbase sequence with classification, species_name, and TSD sequence
    # get all classification
    final_repbase_contigs = {}
    for query_name in new_names:
        label_item = repbase_labels[query_name]
        new_name = query_name + '\t' + label_item[0] + '\t' + label_item[1]
        final_repbase_contigs[new_name] = new_contigs[query_name]
    store_fasta(final_repbase_contigs, repbase_path)

    # Store segmented LTR-to-complete LTR correspondences.
    SegLTR2intactLTRMap = config.work_dir + '/segLTR2intactLTR.map'
    with open(SegLTR2intactLTRMap, 'a+') as f_save:
        for name in SegLTR2intactLTR.keys():
            intact_ltr_name = SegLTR2intactLTR[name]
            f_save.write(name + '\t' + intact_ltr_name + '\n')
    return repbase_path, repbase_labels


def generate_terminal_info(data_path, work_dir, tool_dir, threads):
    output_dir = work_dir + '/temp'
    # Split the file into threads blocks.
    split_files = split_fasta(data_path, output_dir, threads)

    # Parallelize the identification of LTR and TIR.
    cur_update_path = data_path + '.update'
    os.system('rm -f ' + cur_update_path)
    with ProcessPoolExecutor(threads) as executor:
        futures = []
        for split_file in split_files:
            future = executor.submit(identify_terminals, split_file, output_dir, tool_dir)
            futures.append(future)
        executor.shutdown(wait=True)

        is_exit = False
        for future in as_completed(futures):
            update_split_file = future.result()
            if isinstance(update_split_file, str):
                os.system('cat ' + update_split_file + ' >> ' + cur_update_path)
            else:
                print(f"An error occurred: {update_split_file}")
                is_exit = True
                break
        if is_exit:
            print('Error occur, exit...')
            exit(1)
        else:
            os.system('mv ' + cur_update_path + ' ' + data_path)

    return data_path

def generate_domain_info(input_path, domain_path, work_dir, threads):
    output_table = input_path + '.domain'
    temp_dir = work_dir + '/domain'
    get_domain_info(input_path, domain_path, output_table, threads, temp_dir)

def generate_minority_info(train_path, minority_train_path, minority_out, threads, is_train):
    if is_train:
        minority_labels_class = config.minority_labels_class
        minority_contigs = {}
        train_contigNames, train_contigs = read_fasta_v1(train_path)
        # 1. extract minority dataset
        for name in train_contigNames:
            label = name.split('\t')[1]
            if minority_labels_class.__contains__(label):
                minority_contigs[name] = train_contigs[name]
        store_fasta(minority_contigs, minority_train_path)
    # elif not os.path.exists(minority_train_path):
    #     print('We are currently in the model prediction step, attempting to use the minority feature. '
    #           'However, the minority data from the training set at: ' + minority_train_path + ' cannot be found. '
    #                                                                                     'Please verify if this data exists or consider setting the parameter `--use_minority 0`.')
    #     sys.exit(-1)
    #
    # # 2. conduct blastn alignment
    # blastn2Results_path = minority_out
    # os.system('makeblastdb -in ' + minority_train_path + ' -dbtype nucl')
    # align_command = 'blastn -db ' + minority_train_path + ' -num_threads ' \
    #                 + str(threads) + ' -query ' + train_path + ' -evalue 1e-20 -outfmt 6 > ' + blastn2Results_path
    # os.system(align_command)
    # return blastn2Results_path

def store2file(data_partition, cur_consensus_path):
    if len(data_partition) > 0:
        with open(cur_consensus_path, 'w') as f_save:
            for item in data_partition:
                f_save.write('>'+item[0]+'\n'+item[1]+'\n')
        f_save.close()

def PET(seq_item, partitions):
    # sort contigs by length
    original = seq_item
    original = sorted(original, key=lambda x: len(x[1]), reverse=True)
    return divided_array(original, partitions)

def divided_array(original_array, partitions):
    final_partitions = [[] for _ in range(partitions)]
    node_index = 0

    read_from_start = True
    read_from_end = False
    i = 0
    j = len(original_array) - 1
    while i <= j:
        # read from file start
        if read_from_start:
            final_partitions[node_index % partitions].append(original_array[i])
            i += 1
        if read_from_end:
            final_partitions[node_index % partitions].append(original_array[j])
            j -= 1
        node_index += 1
        if node_index % partitions == 0:
            # reverse
            read_from_end = bool(1 - read_from_end)
            read_from_start = bool(1 - read_from_start)
    return final_partitions

def get_domain_info(cons, lib, output_table, threads, temp_dir):
    if os.path.exists(temp_dir):
        os.system('rm -rf ' + temp_dir)
    os.makedirs(temp_dir)

    blast_db_command = 'makeblastdb -dbtype prot -in ' + lib
    os.system(blast_db_command)
    # 1. Divide the cons, and for each block, use blastx -num_threads 1 -evalue 1e-20 to compare cons with domain.
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
        cur_table = temp_dir + '/' + str(partition_index) + '.tbl'
        cur_file = (cur_consensus_path, lib, cur_output, cur_table)
        job = ex.submit(multiple_alignment_blastx_v1, cur_file, merge_distance)
        jobs.append(job)
    ex.shutdown(wait=True)

    # 2. Generate a table of the best matches between query and domain.
    os.system("echo 'TE_name\tdomain_name\tTE_start\tTE_end\tdomain_start\tdomain_end\n' > " + output_table)
    is_exit = False
    for job in as_completed(jobs):
        cur_table = job.result()
        if isinstance(cur_table, str):
            os.system('cat ' + cur_table + ' >> ' + output_table)
        else:
            print(f"An error occurred: {cur_table}")
            is_exit = True
            break
    if is_exit:
        print('Error occur, exit...')
        exit(1)


def multiple_alignment_blastx_v1(repeats_path, merge_distance):
    try:
        split_repeats_path = repeats_path[0]
        protein_db_path = repeats_path[1]
        blastx2Results_path = repeats_path[2]
        cur_table = repeats_path[3]
        align_command = 'blastx -db ' + protein_db_path + ' -num_threads ' \
                        + str(1) + ' -evalue 1e-20 -query ' + split_repeats_path + ' -outfmt 6 > ' + blastx2Results_path
        #os.system(align_command)
        run_command(align_command)

        fixed_extend_base_threshold = merge_distance
        # Merge segmented blastx alignments.
        query_names, query_contigs = read_fasta(split_repeats_path)

        # parse blastn output, determine the repeat boundary
        # query_records = {query_name: {subject_name: [(q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end), (q_start, q_end, s_start, s_end)] }}
        query_records = {}
        with open(blastx2Results_path, 'r') as f_r:
            for idx, line in enumerate(f_r):
                # print('current line idx: %d' % (idx))
                parts = line.split('\t')
                query_name = parts[0]
                subject_name = parts[1]
                identity = float(parts[2])
                alignment_len = int(parts[3])
                q_start = int(parts[6])
                q_end = int(parts[7])
                s_start = int(parts[8])
                s_end = int(parts[9])
                if not query_records.__contains__(query_name):
                    query_records[query_name] = {}
                subject_dict = query_records[query_name]

                if not subject_dict.__contains__(subject_name):
                    subject_dict[subject_name] = []
                subject_pos = subject_dict[subject_name]
                subject_pos.append((q_start, q_end, s_start, s_end))
        f_r.close()

        keep_longest_query = {}
        longest_repeats = {}
        for idx, query_name in enumerate(query_records.keys()):
            query_len = len(query_contigs[query_name])
            # print('total query size: %d, current query name: %s, idx: %d' % (len(query_records), query_name, idx))

            subject_dict = query_records[query_name]

            # if there are more than one longest query overlap with the final longest query over 90%,
            # then it probably the true TE
            longest_queries = []
            for subject_name in subject_dict.keys():
                subject_pos = subject_dict[subject_name]
                # subject_pos.sort(key=lambda x: (x[2], x[3]))

                # cluster all closed fragments, split forward and reverse records
                forward_pos = []
                reverse_pos = []
                for pos_item in subject_pos:
                    if pos_item[0] > pos_item[1]:
                        reverse_pos.append(pos_item)
                    else:
                        forward_pos.append(pos_item)
                forward_pos.sort(key=lambda x: (x[2], x[3]))
                reverse_pos.sort(key=lambda x: (-x[0], -x[1]))

                clusters = {}
                cluster_index = 0
                for k, frag in enumerate(forward_pos):
                    if not clusters.__contains__(cluster_index):
                        clusters[cluster_index] = []
                    cur_cluster = clusters[cluster_index]
                    if k == 0:
                        cur_cluster.append(frag)
                    else:
                        is_closed = False
                        for exist_frag in reversed(cur_cluster):
                            if (frag[0] - exist_frag[1] < fixed_extend_base_threshold):
                                is_closed = True
                                break
                        if is_closed:
                            cur_cluster.append(frag)
                        else:
                            cluster_index += 1
                            if not clusters.__contains__(cluster_index):
                                clusters[cluster_index] = []
                            cur_cluster = clusters[cluster_index]
                            cur_cluster.append(frag)

                cluster_index += 1
                for k, frag in enumerate(reverse_pos):
                    if not clusters.__contains__(cluster_index):
                        clusters[cluster_index] = []
                    cur_cluster = clusters[cluster_index]
                    if k == 0:
                        cur_cluster.append(frag)
                    else:
                        is_closed = False
                        for exist_frag in reversed(cur_cluster):
                            if (exist_frag[1] - frag[0] < fixed_extend_base_threshold):
                                is_closed = True
                                break
                        if is_closed:
                            cur_cluster.append(frag)
                        else:
                            cluster_index += 1
                            if not clusters.__contains__(cluster_index):
                                clusters[cluster_index] = []
                            cur_cluster = clusters[cluster_index]
                            cur_cluster.append(frag)

                for cluster_index in clusters.keys():
                    cur_cluster = clusters[cluster_index]
                    cur_cluster.sort(key=lambda x: (x[2], x[3]))

                    cluster_longest_query_start = -1
                    cluster_longest_query_end = -1
                    cluster_longest_query_len = -1

                    cluster_longest_subject_start = -1
                    cluster_longest_subject_end = -1
                    cluster_longest_subject_len = -1

                    cluster_extend_num = 0

                    # print('subject pos size: %d' %(len(cur_cluster)))
                    # record visited fragments
                    visited_frag = {}
                    for i in range(len(cur_cluster)):
                        # keep a longest query start from each fragment
                        origin_frag = cur_cluster[i]
                        if visited_frag.__contains__(origin_frag):
                            continue
                        cur_frag_len = abs(origin_frag[1] - origin_frag[0])
                        cur_longest_query_len = cur_frag_len
                        longest_query_start = origin_frag[0]
                        longest_query_end = origin_frag[1]
                        longest_subject_start = origin_frag[2]
                        longest_subject_end = origin_frag[3]

                        cur_extend_num = 0

                        visited_frag[origin_frag] = 1
                        # try to extend query
                        for j in range(i + 1, len(cur_cluster)):
                            ext_frag = cur_cluster[j]
                            if visited_frag.__contains__(ext_frag):
                                continue

                            # could extend
                            # extend right
                            if ext_frag[3] > longest_subject_end:
                                # judge query direction
                                if longest_query_start < longest_query_end and ext_frag[0] < ext_frag[1]:
                                    # +
                                    if ext_frag[1] > longest_query_end:
                                        # forward extend
                                        if ext_frag[0] - longest_query_end < fixed_extend_base_threshold and ext_frag[
                                            2] - longest_subject_end < fixed_extend_base_threshold / 3:
                                            # update the longest path
                                            longest_query_start = longest_query_start
                                            longest_query_end = ext_frag[1]
                                            longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                             ext_frag[
                                                                                                 2] else ext_frag[2]
                                            longest_subject_end = ext_frag[3]
                                            cur_longest_query_len = longest_query_end - longest_query_start
                                            cur_extend_num += 1
                                            visited_frag[ext_frag] = 1
                                        elif ext_frag[0] - longest_query_end >= fixed_extend_base_threshold:
                                            break
                                elif longest_query_start > longest_query_end and ext_frag[0] > ext_frag[1]:
                                    # reverse
                                    if ext_frag[1] < longest_query_end:
                                        # reverse extend
                                        if longest_query_end - ext_frag[0] < fixed_extend_base_threshold and ext_frag[
                                            2] - longest_subject_end < fixed_extend_base_threshold / 3:
                                            # update the longest path
                                            longest_query_start = longest_query_start
                                            longest_query_end = ext_frag[1]
                                            longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                             ext_frag[
                                                                                                 2] else ext_frag[2]
                                            longest_subject_end = ext_frag[3]
                                            cur_longest_query_len = longest_query_start - longest_query_end
                                            cur_extend_num += 1
                                            visited_frag[ext_frag] = 1
                                        elif longest_query_end - ext_frag[0] >= fixed_extend_base_threshold:
                                            break
                        if cur_longest_query_len > cluster_longest_query_len:
                            cluster_longest_query_start = longest_query_start
                            cluster_longest_query_end = longest_query_end
                            cluster_longest_query_len = cur_longest_query_len

                            cluster_longest_subject_start = longest_subject_start
                            cluster_longest_subject_end = longest_subject_end
                            cluster_longest_subject_len = longest_subject_end - longest_subject_start

                            cluster_extend_num = cur_extend_num
                    # keep this longest query
                    if cluster_longest_query_len != -1:
                        longest_queries.append((cluster_longest_query_start, cluster_longest_query_end,
                                                cluster_longest_query_len, cluster_longest_subject_start,
                                                cluster_longest_subject_end, cluster_longest_subject_len, subject_name,
                                                cluster_extend_num))

            # we now consider, we should take some sequences from longest_queries to represent this query sequence.
            # we take the longest sequence by length, if the latter sequence overlap with the former sequence largely (50%),
            # continue find next sequence until the ratio of query sequence over 90% or no more sequences.
            longest_queries.sort(key=lambda x: -x[2])
            keep_longest_query[query_name] = longest_queries
        # print(keep_longest_query)
        # Save the record as a table, removing redundant records
        # (more than 50% of the latter sequence's area is within the former).
        with open(cur_table, 'w') as f_save:
            for query_name in keep_longest_query.keys():
                domain_array = keep_longest_query[query_name]
                # for domain_info in domain_array:
                #     f_save.write(query_name+'\t'+str(domain_info[6])+'\t'+str(domain_info[0])+'\t'+str(domain_info[1])+'\t'+str(domain_info[3])+'\t'+str(domain_info[4])+'\n')
                merge_domains = []
                # Merge domain_array.
                domain_array.sort(key=lambda x: -x[2])
                for domain_info in domain_array:
                    if len(merge_domains) == 0:
                        merge_domains.append(domain_info)
                    else:
                        is_new_domain = True
                        for pre_domain in merge_domains:
                            pre_start = pre_domain[0]
                            pre_end = pre_domain[1]
                            # Calculate overlap.
                            if pre_start > pre_end:
                                tmp = pre_start
                                pre_start = pre_end
                                pre_end = tmp
                            cur_start = domain_info[0]
                            cur_end = domain_info[1]
                            if cur_start > cur_end:
                                tmp = cur_start
                                cur_start = cur_end
                                cur_end = tmp
                            if cur_end >= pre_start and cur_end <= pre_end:
                                if cur_start <= pre_start:
                                    overlap = cur_end - pre_start
                                else:
                                    overlap = cur_end - cur_start
                            elif cur_end > pre_end:
                                if cur_start >= pre_start and cur_start <= pre_end:
                                    overlap = pre_end - cur_start
                                else:
                                    overlap = 0
                            else:
                                overlap = 0

                            if float(overlap / domain_info[2]) > 0.5:
                                is_new_domain = False
                        if is_new_domain:
                            merge_domains.append(domain_info)

                for domain_info in merge_domains:
                    f_save.write(query_name + '\t' + str(domain_info[6]) + '\t' + str(domain_info[0]) + '\t' + str(
                        domain_info[1]) + '\t' + str(domain_info[3]) + '\t' + str(domain_info[4]) + '\n')

        f_save.close()
        return cur_table
    except Exception as e:
        return e

def generate_TSD_info(all_repbase_path, ncbi_ref_info, work_dir, is_expanded, keep_raw, threads):
    print('Preprocess: Start get TSD information')
    # Step1. Count sequences by species.
    names, contigs = read_fasta_v1(all_repbase_path)
    speices_TE_contigs = {}
    species_TE_summary = {}
    for name in names:
        parts = name.split('\t')
        species_name = parts[2]
        label = parts[1]
        TE_name = parts[0]
        if not species_TE_summary.__contains__(species_name):
            species_TE_summary[species_name] = {}
        label_TE_summary = species_TE_summary[species_name]
        if not label_TE_summary.__contains__(label):
            label_TE_summary[label] = 0
        label_count = label_TE_summary[label]
        label_count += 1
        label_TE_summary[label] = label_count

        if not speices_TE_contigs.__contains__(species_name):
            speices_TE_contigs[species_name] = {}
        TE_contigs = speices_TE_contigs[species_name]
        TE_contigs[name] = contigs[name]

    all_sequence_count = 0
    species_count = []
    species_count_dict = {}
    for species_name in species_TE_summary:
        label_TE_summary = species_TE_summary[species_name]
        total_num = 0
        for label in label_TE_summary.keys():
            total_num += label_TE_summary[label]
        all_sequence_count += total_num
        species_count.append((species_name, total_num))
        species_count_dict[species_name] = total_num
    species_count.sort(key=lambda x: -x[1])

    # Step2. Store the number of TE sequences corresponding to species in the training data.
    data = {}
    species_names = []
    TE_sequences_nums = []
    for item in species_count:
        species_names.append(item[0])
        TE_sequences_nums.append(item[1])
    data['Species Name'] = species_names
    data['Total TE number'] = TE_sequences_nums

    temp_dir = work_dir + '/temp'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    df = pd.DataFrame(data)
    # Save the DataFrame to an Excel file.
    with pd.ExcelWriter(temp_dir + '/species_TE_num.xlsx', engine="openpyxl") as writer:
        to_excel_auto_column_weight(df, writer, f'Species TE number information')

    # Check how many species we actually have genomes for in the top 100 species.
    keep_species_names = set()
    species_genome = {}
    with open(ncbi_ref_info, 'r') as f_r:
        for line in f_r:
            if line.startswith('#'):
                continue
            species_name = line.split('\t')[0]
            genome = line.split('\t')[1]
            is_plant = line.split('\t')[2]
            species_genome[species_name] = (genome, is_plant)
            keep_species_names.add(species_name)

    top_num = 100
    top_seq_count = 0
    not_keep_species_names = []
    for i, species_name in enumerate(species_names):
        if i > top_num:
            break
        if species_name not in keep_species_names:
            not_keep_species_names.append((species_name, species_count_dict[species_name]))
        top_seq_count += species_count_dict[species_name]
    print('top num:' + str(top_num) + ', all sequence count:' + str(top_seq_count))

    # Step3. ① Store repbase data by species; ② Align sequences to the corresponding genomes to obtain TSD information.
    # Store TE sequences as separate files according to species name.
    species_dir = work_dir + '/species'
    processed_species_dir = work_dir + '/species_processed'
    if os.path.exists(species_dir):
        os.system('rm -rf ' + species_dir)
    if os.path.exists(processed_species_dir):
        os.system('rm -rf ' + processed_species_dir)
    os.makedirs(species_dir)
    os.makedirs(processed_species_dir)
    species_TE_files = {}
    for species_name in speices_TE_contigs.keys():
        TE_contigs = speices_TE_contigs[species_name]
        species = species_name.replace(' ', '_')
        species_TE_files[species_name] = species_dir + '/' + species + '.ref'
        store_fasta(TE_contigs, species_dir + '/' + species + '.ref')
    for i, species_name in enumerate(species_names):
        if species_genome.__contains__(species_name):
            genome_path = species_genome[species_name][0]
            is_plant = int(species_genome[species_name][1])
            print('species name:' + species_name + ', genome_path:' + genome_path + ', is_plant:' + str(is_plant))
            species = species_name.replace(' ', '_')
            repbase_path = species_TE_files[species_name]
            print(repbase_path)

            flanking_len = 20
            final_repbase_path = expandRepBase(repbase_path, genome_path, temp_dir, threads, flanking_len, is_plant, species, is_expanded=is_expanded)
            os.system('mv ' + final_repbase_path + ' ' + processed_species_dir)
    # Merge all files with TSD.
    tsd_file = work_dir + '/repbase.tsd_merge.ref'
    with open(tsd_file, 'w') as out:
        for filename in os.listdir(processed_species_dir):
            file_path = os.path.join(processed_species_dir, filename)
            if os.path.isfile(file_path):
                with open(file_path, 'r') as input_file:
                    content = input_file.read()
                    out.write(content)
    print('Preprocess: End get TSD information')

    if keep_raw == 1:
        keep_file = work_dir + '/repbase.tsd_processed.ref'
        tsd_names, tsd_contigs = read_fasta_v1(tsd_file)
        seq_name_to_tsd_header = {}
        for tsd_name in tsd_names:
            parts = tsd_name.split('\t')
            seq_name = parts[0]
            seq_name_to_tsd_header[seq_name] = tsd_name
        with open(keep_file, 'w') as f_save:
            for name in names:
                seq_name = name.split('\t')[0]
                if seq_name_to_tsd_header.__contains__(seq_name):
                    new_name = seq_name_to_tsd_header[seq_name]
                else:
                    new_name = name + '\t' + 'TSD:' + '\t' + 'TSD_len:0'
                f_save.write('>'+new_name+'\n'+contigs[name]+'\n')
        os.system('mv ' + keep_file + ' ' + all_repbase_path)
    else:
        os.system('mv ' + tsd_file + ' ' + all_repbase_path)
    return all_repbase_path

def search_TSD_regular(motif, sequence):
    motif_length = len(motif)
    pattern = ''

    # Build a regular expression pattern based on motif length.
    if motif_length >= 8:
        for i in range(motif_length):
            pattern += f"{motif[:i]}[ACGT]{motif[i + 1:]}" if i < motif_length - 1 else motif[:i] + "[ACGT]"
            if i < motif_length - 1:
                pattern += "|"
    else:
        pattern = motif

    matches = re.finditer(pattern, sequence)

    found = False
    pos = None
    for match in matches:
        #print(f"Found motif at position {match.start()}: {match.group()}")
        found = True
        pos = match.start()
        break
    return found, pos

def search_confident_tsd(orig_seq, raw_tir_start, raw_tir_end, tsd_search_distance):
    # Change all coordinates to start from 0.
    raw_tir_start -= 1
    raw_tir_end -= 1

    orig_seq_len = len(orig_seq)
    # 1. First, take 2 * tsd_search_distance sequences near the start and end positions
    left_start = raw_tir_start - tsd_search_distance
    if left_start < 0:
        left_start = 0
    # We don’t search inwards here because we consider Repbase boundaries to be correct.
    # If we consider the boundaries to be incorrect, many abnormal TSDs may meet the requirements.
    # For simplicity, we assume that Repbase boundaries are correct.
    left_end = raw_tir_start
    left_round_seq = orig_seq[left_start: left_end]
    # Obtain the position offset of left_round_seq relative to the entire sequence to correct the subsequent TSD boundary positions.
    left_offset = left_start
    right_start = raw_tir_end + 1
    if right_start < 0:
        right_start = 0
    right_end = raw_tir_end + tsd_search_distance + 1
    right_round_seq = orig_seq[right_start: right_end]
    # Obtain the position offset of right_round_seq relative to the entire sequence to correct the subsequent TSD boundary positions.
    right_offset = right_start

    # 2. Split the left sequence into k-mers from large to small, then search for the right sequence with k-mers.
    # If found, record as a candidate TSD, and finally select the closest one to the original boundary as the TSD.
    TIR_TSDs = [15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2]
    # Record the position nearest to the original boundary.
    is_found = False
    tsd_set = []
    for k_num in TIR_TSDs:
        for i in range(len(left_round_seq) - k_num, -1, -1):
            left_kmer = left_round_seq[i: i + k_num]
            left_pos = left_offset + i + k_num
            if left_pos < 0 or left_pos > orig_seq_len-1:
                continue
            found_tsd, right_pos = search_TSD_regular(left_kmer, right_round_seq)
            if found_tsd and not left_kmer.__contains__('N'):
                right_pos = right_offset + right_pos - 1
                is_found = True
                # Calculate the distance from the original boundary.
                left_distance = abs(left_pos - raw_tir_start)
                right_distance = abs(right_pos - raw_tir_end)
                distance = left_distance + right_distance
                TSD_seq = left_kmer
                TSD_len = len(TSD_seq)
                tsd_set.append((distance, TSD_len, TSD_seq))
    tsd_set = sorted(tsd_set, key=lambda x: (x[0], -x[1]))

    if not is_found:
        TSD_seq = 'Unknown'
        TSD_len = 16
        min_distance = -1
    else:
        TSD_seq = tsd_set[0][2]
        TSD_len = tsd_set[0][1]
        min_distance = tsd_set[0][0]
    return TSD_seq, TSD_len, min_distance

def getReverseSequence(sequence):
    base_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    res = ''
    length = len(sequence)
    i = length - 1
    while i >= 0:
        base = sequence[i]
        if base not in base_map.keys():
            base = 'N'
        else:
            base = base_map[base]
        res += base
        i -= 1
    return res

def get_full_length_copies(query_path, reference):
    blastn2Results_path = query_path + '.blast.out'
    repeats_path = (query_path, reference, blastn2Results_path)
    all_copies = multiple_alignment_blast_and_get_copies(repeats_path)
    return all_copies, blastn2Results_path

def multiple_alignment_blast_and_get_copies(repeats_path):
    split_repeats_path = repeats_path[0]
    ref_db = repeats_path[1]
    blastn2Results_path = repeats_path[2]
    os.system('rm -f ' + blastn2Results_path)
    all_copies = None
    repeat_names, repeat_contigs = read_fasta(split_repeats_path)
    if len(repeat_contigs) > 0:
        align_command = 'blastn -db ' + ref_db + ' -num_threads ' \
                        + str(1) + ' -query ' + split_repeats_path + ' -evalue 1e-20 -outfmt 6 >> ' + blastn2Results_path
        os.system(align_command)
        all_copies = get_copies_v1(blastn2Results_path, split_repeats_path, '')
    return all_copies

def get_query_copies(cur_segments, query_contigs, subject_path, query_coverage, subject_coverage, query_fixed_extend_base_threshold=1000, subject_fixed_extend_base_threshold=1000, max_copy_num=100):
    all_copies = {}

    if subject_coverage > 0:
        subject_names, subject_contigs = read_fasta(subject_path)

    for item in cur_segments:
        query_name = item[0]
        subject_dict = item[1]

        longest_queries = []
        for subject_name in subject_dict.keys():
            subject_pos = subject_dict[subject_name]

            # cluster all closed fragments, split forward and reverse records
            forward_pos = []
            reverse_pos = []
            for pos_item in subject_pos:
                if pos_item[2] > pos_item[3]:
                    reverse_pos.append(pos_item)
                else:
                    forward_pos.append(pos_item)
            forward_pos.sort(key=lambda x: (x[2], x[3]))
            reverse_pos.sort(key=lambda x: (-x[2], -x[3]))

            clusters = {}
            cluster_index = 0
            for k, frag in enumerate(forward_pos):
                if not clusters.__contains__(cluster_index):
                    clusters[cluster_index] = []
                cur_cluster = clusters[cluster_index]
                if k == 0:
                    cur_cluster.append(frag)
                else:
                    is_closed = False
                    for exist_frag in reversed(cur_cluster):
                        if (frag[2] - exist_frag[3] < subject_fixed_extend_base_threshold and frag[1] > exist_frag[1]):
                            is_closed = True
                            break
                    if is_closed:
                        cur_cluster.append(frag)
                    else:
                        cluster_index += 1
                        if not clusters.__contains__(cluster_index):
                            clusters[cluster_index] = []
                        cur_cluster = clusters[cluster_index]
                        cur_cluster.append(frag)

            cluster_index += 1
            for k, frag in enumerate(reverse_pos):
                if not clusters.__contains__(cluster_index):
                    clusters[cluster_index] = []
                cur_cluster = clusters[cluster_index]
                if k == 0:
                    cur_cluster.append(frag)
                else:
                    is_closed = False
                    for exist_frag in reversed(cur_cluster):
                        if (exist_frag[3] - frag[2] < subject_fixed_extend_base_threshold and frag[1] > exist_frag[1]):
                            is_closed = True
                            break
                    if is_closed:
                        cur_cluster.append(frag)
                    else:
                        cluster_index += 1
                        if not clusters.__contains__(cluster_index):
                            clusters[cluster_index] = []
                        cur_cluster = clusters[cluster_index]
                        cur_cluster.append(frag)

            for cluster_index in clusters.keys():
                cur_cluster = clusters[cluster_index]
                cur_cluster.sort(key=lambda x: (x[0], x[1]))

                cluster_longest_query_start = -1
                cluster_longest_query_end = -1
                cluster_longest_query_len = -1

                cluster_longest_subject_start = -1
                cluster_longest_subject_end = -1
                cluster_longest_subject_len = -1

                cluster_identity = 0
                cluster_extend_num = 0

                # record visited fragments
                visited_frag = {}
                for i in range(len(cur_cluster)):
                    # keep a longest query start from each fragment
                    origin_frag = cur_cluster[i]
                    if visited_frag.__contains__(origin_frag):
                        continue
                    cur_frag_len = origin_frag[1] - origin_frag[0] + 1
                    cur_longest_query_len = cur_frag_len
                    longest_query_start = origin_frag[0]
                    longest_query_end = origin_frag[1]
                    longest_subject_start = origin_frag[2]
                    longest_subject_end = origin_frag[3]

                    cur_identity = origin_frag[4]
                    cur_extend_num = 0

                    visited_frag[origin_frag] = 1
                    # try to extend query
                    for j in range(i + 1, len(cur_cluster)):
                        ext_frag = cur_cluster[j]
                        if visited_frag.__contains__(ext_frag):
                            continue

                        # could extend
                        # extend right
                        if ext_frag[1] > longest_query_end:
                            # judge subject direction
                            if longest_subject_start < longest_subject_end and ext_frag[2] < ext_frag[3]:
                                # +
                                if ext_frag[3] > longest_subject_end:
                                    # forward extend
                                    if ext_frag[0] - longest_query_end < query_fixed_extend_base_threshold and ext_frag[
                                        2] - longest_subject_end < subject_fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start < \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start

                                        cur_identity += ext_frag[4]
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= query_fixed_extend_base_threshold:
                                        break
                            elif longest_subject_start > longest_subject_end and ext_frag[2] > ext_frag[3]:
                                # reverse
                                if ext_frag[3] < longest_subject_end:
                                    # reverse extend
                                    if ext_frag[
                                        0] - longest_query_end < query_fixed_extend_base_threshold and longest_subject_end - \
                                            ext_frag[2] < subject_fixed_extend_base_threshold:
                                        # update the longest path
                                        longest_query_start = longest_query_start
                                        longest_query_end = ext_frag[1]
                                        longest_subject_start = longest_subject_start if longest_subject_start > \
                                                                                         ext_frag[
                                                                                             2] else ext_frag[2]
                                        longest_subject_end = ext_frag[3]
                                        cur_longest_query_len = longest_query_end - longest_query_start

                                        cur_identity += ext_frag[4]
                                        cur_extend_num += 1
                                        visited_frag[ext_frag] = 1
                                    elif ext_frag[0] - longest_query_end >= query_fixed_extend_base_threshold:
                                        break
                    if cur_longest_query_len > cluster_longest_query_len:
                        cluster_longest_query_start = longest_query_start
                        cluster_longest_query_end = longest_query_end
                        cluster_longest_query_len = cur_longest_query_len

                        cluster_longest_subject_start = longest_subject_start
                        cluster_longest_subject_end = longest_subject_end
                        cluster_longest_subject_len = abs(longest_subject_end - longest_subject_start) + 1

                        cluster_identity = cur_identity
                        cluster_extend_num = cur_extend_num
                # keep this longest query
                if cluster_longest_query_len != -1:
                    longest_queries.append((cluster_longest_query_start, cluster_longest_query_end,
                                            cluster_longest_query_len, cluster_longest_subject_start,
                                            cluster_longest_subject_end, cluster_longest_subject_len, subject_name,
                                            cluster_extend_num, cluster_identity))

        longest_queries.sort(key=lambda x: -x[2])
        query_len = len(query_contigs[query_name])
        # query_len = int(query_name.split('-')[1].split('_')[1])
        copies = []
        keeped_copies = set()
        for query in longest_queries:
            if len(copies) > max_copy_num:
                break
            subject_name = query[6]
            subject_start = query[3]
            subject_end = query[4]
            direct = '+'
            if subject_start > subject_end:
                tmp = subject_start
                subject_start = subject_end
                subject_end = tmp
                direct = '-'
            item = (subject_name, subject_start, subject_end)
            if subject_coverage > 0:
                subject_len = len(subject_contigs[subject_name])
                cur_subject_coverage = float(query[5])/subject_len
                if float(query[2])/query_len >= query_coverage and cur_subject_coverage >= subject_coverage and item not in keeped_copies:
                    copies.append((subject_name, subject_start, subject_end, query[2], direct))
                    keeped_copies.add(item)
            else:
                if float(query[2]) / query_len >= query_coverage and item not in keeped_copies:
                    copies.append((subject_name, subject_start, subject_end, query[2], direct))
                    keeped_copies.add(item)
        #copies.sort(key=lambda x: abs(x[3]-(x[2]-x[1]+1)))
        all_copies[query_name] = copies
    return all_copies

def get_copies_v1(blastnResults_path, query_path, subject_path, query_coverage=0.95, subject_coverage=0):
    query_records = {}
    with open(blastnResults_path, 'r') as f_r:
        for idx, line in enumerate(f_r):
            parts = line.split('\t')
            query_name = parts[0]
            subject_name = parts[1]
            identity = float(parts[2])
            alignment_len = int(parts[3])
            q_start = int(parts[6])
            q_end = int(parts[7])
            s_start = int(parts[8])
            s_end = int(parts[9])
            if query_name == subject_name:
                continue
            if not query_records.__contains__(query_name):
                query_records[query_name] = {}
            subject_dict = query_records[query_name]

            if not subject_dict.__contains__(subject_name):
                subject_dict[subject_name] = []
            subject_pos = subject_dict[subject_name]
            subject_pos.append((q_start, q_end, s_start, s_end, identity))
    f_r.close()

    query_names, query_contigs = read_fasta(query_path)
    cur_segments = list(query_records.items())
    all_copies = get_query_copies(cur_segments, query_contigs, subject_path, query_coverage, subject_coverage)

    return all_copies

def read_fasta(fasta_path):
    contignames = []
    contigs = {}
    if os.path.exists(fasta_path):
        with open(fasta_path, 'r') as rf:
            contigname = ''
            contigseq = ''
            for line in rf:
                if line.startswith('>'):
                    if contigname != '' and contigseq != '':
                        contigs[contigname] = contigseq
                        contignames.append(contigname)
                    contigname = line.strip()[1:].split(" ")[0].split('\t')[0]
                    contigseq = ''
                else:
                    contigseq += line.strip().upper()
            if contigname != '' and contigseq != '':
                contigs[contigname] = contigseq
                contignames.append(contigname)
        rf.close()
    return contignames, contigs

def read_fasta_v1(fasta_path):
    contignames = []
    contigs = {}
    if os.path.exists(fasta_path):
        with open(fasta_path, 'r') as rf:
            contigname = ''
            contigseq = ''
            for line in rf:
                if line.startswith('>'):
                    if contigname != '' and contigseq != '':
                        contigs[contigname] = contigseq
                        contignames.append(contigname)
                    contigname = line.strip()[1:]
                    contigseq = ''
                else:
                    contigseq += line.strip().upper()
            if contigname != '' and contigseq != '':
                contigs[contigname] = contigseq
                contignames.append(contigname)
        rf.close()
    return contignames, contigs


def store_fasta(contigs, file_path):
    with open(file_path, 'w') as f_save:
        for name in contigs.keys():
            seq = contigs[name]
            f_save.write('>'+name+'\n'+seq+'\n')
    f_save.close()


def get_flanking_copies(processed_TE_path, genome_path, flanking_len, temp_dir, threads):
    names, contigs = read_fasta(processed_TE_path)
    # Step1. align the TEs to genome, and get copies
    os.system('makeblastdb -in ' + genome_path + ' -dbtype nucl')
    batch_size = 10
    batch_id = 0
    total_names = set(names)
    split_files = []
    cur_contigs = {}
    for i, name in enumerate(names):
        cur_file = temp_dir + '/' + str(batch_id) + '.fa'
        cur_contigs[name] = contigs[name]
        if len(cur_contigs) == batch_size:
            store_fasta(cur_contigs, cur_file)
            split_files.append(cur_file)
            cur_contigs = {}
            batch_id += 1
    if len(cur_contigs) > 0:
        cur_file = temp_dir + '/' + str(batch_id) + '.fa'
        store_fasta(cur_contigs, cur_file)
        split_files.append(cur_file)
        batch_id += 1

    blastn2Results_path = temp_dir + '/all.out'
    os.system('rm -f ' + blastn2Results_path)
    ex = ProcessPoolExecutor(threads)
    jobs = []
    for cur_split_files in split_files:
        job = ex.submit(get_full_length_copies, cur_split_files, genome_path)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_copies = {}
    for job in as_completed(jobs):
        cur_all_copies, cur_blastn2Results_path = job.result()
        all_copies.update(cur_all_copies)
        os.system('cat ' + cur_blastn2Results_path + ' >> ' + blastn2Results_path)

    # Step2. flank all copies to obtain TSDs.
    ref_names, ref_contigs = read_fasta(genome_path)
    batch_member_files = []
    new_all_copies = {}
    for query_name in all_copies.keys():
        copies = all_copies[query_name]
        for copy in copies:
            ref_name = copy[0]
            copy_ref_start = int(copy[1])
            copy_ref_end = int(copy[2])
            direct = copy[4]
            copy_len = copy_ref_end - copy_ref_start + 1
            if copy_ref_start - 1 - flanking_len < 0 or copy_ref_end + flanking_len > len(ref_contigs[ref_name]):
                continue
            copy_seq = ref_contigs[ref_name][copy_ref_start - 1 - flanking_len: copy_ref_end + flanking_len]
            if direct == '-':
                copy_seq = getReverseSequence(copy_seq)
            if len(copy_seq) < 100:
                continue
            new_name = ref_name + ':' + str(copy_ref_start) + '-' + str(copy_ref_end) + '(' + direct + ')'
            if not new_all_copies.__contains__(query_name):
                new_all_copies[query_name] = {}
            copy_contigs = new_all_copies[query_name]
            copy_contigs[new_name] = copy_seq
            new_all_copies[query_name] = copy_contigs
            if len(cur_contigs) >= 100:
                break
    for query_name in new_all_copies.keys():
        copy_contigs = new_all_copies[query_name]
        cur_member_file = temp_dir + '/' + query_name + '.blast.bed.fa'
        store_fasta(copy_contigs, cur_member_file)
        query_seq = contigs[query_name]
        batch_member_files.append((query_name, query_seq, cur_member_file))
    return batch_member_files

def get_seq_TSDs(batch_member_file, flanking_len, is_expanded, expandClassNum, repbase_labels):
    tsd_info = {}
    cur_query_name = batch_member_file[0]
    tsd_info[cur_query_name] = []
    copies_tsd_info = tsd_info[cur_query_name]

    if is_expanded == 0:
        max_copy_num = 1
    else:
        cur_label = repbase_labels[cur_query_name][0]
        if expandClassNum.__contains__(cur_label):
            max_copy_num = expandClassNum[cur_label]
        else:
            max_copy_num = 1

    cur_member_file = batch_member_file[2]
    cur_contignames, cur_contigs = read_fasta(cur_member_file)
    # Get all copies corresponding to TSD.
    for i, query_name in enumerate(cur_contignames):
        seq = cur_contigs[query_name]
        tir_start = flanking_len + 1
        tir_end = len(seq) - flanking_len
        tsd_search_distance = flanking_len
        cur_tsd, cur_tsd_len, min_distance = search_confident_tsd(seq, tir_start, tir_end, tsd_search_distance)
        copies_tsd_info.append((cur_tsd, cur_tsd_len, min_distance))
    return tsd_info

def get_copies_TSD_info(batch_member_files, flanking_len, is_expanded, repbase_labels, threads):
    # search TSDs in all flanked copies.
    expandClassNum = config.expandClassNum

    ex = ProcessPoolExecutor(threads)
    jobs = []
    tsd_info = {}
    for batch_member_file in batch_member_files:
        job = ex.submit(get_seq_TSDs, batch_member_file, flanking_len, is_expanded, expandClassNum, repbase_labels)
        jobs.append(job)
    ex.shutdown(wait=True)
    all_copies = {}
    for job in as_completed(jobs):
        cur_tsd_info = job.result()
        tsd_info.update(cur_tsd_info)

    return tsd_info

def get_copies_TSD_info_bak(batch_member_files, flanking_len, plant, is_expanded, repbase_labels):
    # search TSDs in all flanked copies.
    expandClassNum = config.expandClassNum

    tsd_info = {}
    for batch_member_file in batch_member_files:
        cur_query_name = batch_member_file[0]
        tsd_info[cur_query_name] = []
        copies_tsd_info = tsd_info[cur_query_name]

        if is_expanded == 0:
            max_copy_num = 1
        else:
            cur_label = repbase_labels[cur_query_name][0]
            if expandClassNum.__contains__(cur_label):
                max_copy_num = expandClassNum[cur_label]
            else:
                max_copy_num = 1

        cur_member_file = batch_member_file[2]
        cur_contignames, cur_contigs = read_fasta(cur_member_file)
        # summarize the count of tsd_len, get the one occurs most times
        for i, query_name in enumerate(cur_contignames):
            if i >= max_copy_num:
                break
            seq = cur_contigs[query_name]
            tir_start = flanking_len + 1
            tir_end = len(seq) - flanking_len
            tsd_search_distance = flanking_len
            cur_tsd, cur_tsd_len, cur_TE = search_confident_tsd(seq, tir_start, tir_end, tsd_search_distance, query_name, plant)
            copies_tsd_info.append((cur_tsd, cur_tsd_len, cur_TE))
    return tsd_info

def expandRepBase(repbase_path, genome_path, temp_dir, threads, flanking_len, plant, species, is_expanded=0):
    # The purpose of this function is to expand the relatively less abundant categories in the Repbase
    # database to balance the entire Repbase dataset. The expansion method is as follows:
    # 1. Align all data to the genome to obtain their copies.
    # 2. Get the number of copies needed for each type of sequence expansion.
    # For example, Merlin: 20 times. The number of expansions for each type varies.
    # Refer to the expandClassNum parameter for details.
    # By default, no expansion is performed, that is, only the TSD of the current sequence is obtained.
    # Expansion is only carried out when we want to obtain a balanced dataset.

    os.system('rm -rf ' + temp_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    label_names, label_contigs = read_fasta_v1(repbase_path)
    # store repbase name and label
    repbase_labels = {}
    for name in label_names:
        parts = name.split('\t')
        repbase_name = parts[0]
        classification = parts[1]
        species_name = parts[2]
        repbase_labels[repbase_name] = (classification, species_name)
    # # Concatenate LTR elements in Repbase.
    # processed_TE_path, repbase_labels = connect_LTR(repbase_path)
    # Get all copies of elements with side sequence. Returns an array,
    # where each element is a tuple (query_name, query_seq, cur_member_file),
    # and cur_member_file is the collection file corresponding to this element
    batch_member_files = get_flanking_copies(repbase_path, genome_path, flanking_len, temp_dir, threads)
    # Obtain TSD information corresponding to each copy of the original sequence
    # and the number of copies for each classification, according to expandClassNum.
    tsd_info = get_copies_TSD_info(batch_member_files, flanking_len, is_expanded, repbase_labels, threads)
    # Save all sequences as files, with the copy sequence name added with -C_{num}.
    names, contigs = read_fasta(repbase_path)
    final_repbase_path = temp_dir + '/' + species + '.ref'
    final_repbase_contigs = {}
    for query_name in names:
        seq = contigs[query_name]
        label_item = repbase_labels[query_name]

        if tsd_info.__contains__(query_name):
            copies_tsd_info = tsd_info[query_name]
        else:
            copies_tsd_info = [('Unknown', 16, -1)]

        # Traverse and extract TSDs with a distance of 0;
        new_copies_tsd_info = []
        for tsd_seq, tsd_len, cur_distance in copies_tsd_info:
            if cur_distance <= 0:
                new_copies_tsd_info.append((tsd_seq, tsd_len, cur_distance))
        copies_tsd_info = new_copies_tsd_info if len(new_copies_tsd_info) > 0 else copies_tsd_info
        # Store all copies corresponding to TSDs, recording the number of occurrences and
        # distances from the original boundary for each length of TSD.
        max_count_TSD = {}
        length_count = {}
        for tsd_seq, tsd_len, cur_distance in copies_tsd_info:
            if not length_count.__contains__(tsd_len):
                length_count[tsd_len] = (1, cur_distance)
                max_count_TSD[tsd_len] = tsd_seq
            else:
                prev_count, prev_distance = length_count[tsd_len]
                if cur_distance < prev_distance:
                    prev_distance = cur_distance
                    max_count_TSD[tsd_len] = tsd_seq
                length_count[tsd_len] = (prev_count + 1, prev_distance)
        # Store as an array based on (tsd_len, tsd_seq, occurrence count, minimum distance).
        # Take the most frequently occurring TSD. If multiple occurrences are the same,
        # take the one with the smallest distance from the original boundary and the longest length.
        all_tsd_set = []
        for tsd_len in length_count.keys():
            cur_count, cur_distance = length_count[tsd_len]
            tsd_seq = max_count_TSD[tsd_len]
            all_tsd_set.append((tsd_len, tsd_seq, cur_count, cur_distance))
        all_tsd_set = sorted(all_tsd_set, key=lambda x: (-x[2], x[3], -x[0]))
        final_tsd_info = all_tsd_set[0]
        tsd_seq = final_tsd_info[1]
        tsd_len = final_tsd_info[0]
        tsd_distance = final_tsd_info[3]
        if tsd_distance > 5:
            tsd_seq = ''
            tsd_len = len(tsd_seq)
        new_name = query_name + '\t' + label_item[0] + '\t' + label_item[1] + '\t' + 'TSD:' + str(tsd_seq) + '\t' + 'TSD_len:' + str(tsd_len)
        final_repbase_contigs[new_name] = seq
    store_fasta(final_repbase_contigs, final_repbase_path)
    return final_repbase_path

def to_excel_auto_column_weight(df: pd.DataFrame, writer: ExcelWriter, sheet_name="Shee1"):
    df.to_excel(writer, sheet_name=sheet_name, index=False)
    column_widths = (
        df.columns.to_series().apply(lambda x: len(str(x).encode('gbk'))).values
    )
    max_widths = (
        df.astype(str).applymap(lambda x: len(str(x).encode('gbk'))).agg(max).values
    )
    widths = np.max([column_widths, max_widths], axis=0)
    worksheet = writer.sheets[sheet_name]
    for i, width in enumerate(widths, 1):
        worksheet.column_dimensions[get_column_letter(i)].width = width + 2

def split2TrainTest(total_path, train_path, test_path):
    # 1. Split the Repbase data into a 8-2 ratio for training and testing sets.
    names, contigs = read_fasta_v1(total_path)
    # Shuffle the list randomly.
    random.shuffle(names)
    # Calculate the split index positions.
    split_index = int(0.8 * len(names))
    # Split into two lists of 80% and 20%.
    train_list = names[:split_index]
    test_list = names[split_index:]
    train_contigs = {}
    test_contigs = {}
    for name in train_list:
        train_contigs[name] = contigs[name]
    for name in test_list:
        test_contigs[name] = contigs[name]
    store_fasta(train_contigs, train_path)
    store_fasta(test_contigs, test_path)

def get_species_TE(total_repbase, mazie_path, mazie_merge_path, species_name):
    names, contigs = read_fasta_v1(total_repbase)
    merge_contigs = {}
    new_contigs = {}
    for name in names:
        parts = name.split('\t')
        seq_name = parts[0]
        label = parts[1]
        if parts[2] == species_name:
            merge_contigs[seq_name+'#'+label] = contigs[name]
            new_contigs[name] = contigs[name]
    store_fasta(merge_contigs, mazie_merge_path)
    store_fasta(new_contigs, mazie_path)

def extract_60_species(total_repbase, train_species_path, test_species_path):
    names, contigs = read_fasta_v1(total_repbase)
    species_names = set()
    for name in names:
        parts = name.split('\t')
        species_name = parts[2]
        species_names.add(species_name)
    species_names = list(species_names)
    train_contigs = {}
    test_contigs = {}
    random.shuffle(species_names)
    split_index = int(0.9 * len(species_names))
    train_list = species_names[:split_index]
    test_list = species_names[split_index:]
    for name in names:
        parts = name.split('\t')
        seq_name = parts[0]
        label = parts[1]
        species_name = parts[2]
        if species_name in train_list:
            train_contigs[name] = contigs[name]
        if species_name in test_list:
            test_contigs[name] = contigs[name]
    store_fasta(train_contigs, train_species_path)
    store_fasta(test_contigs, test_species_path)
    print(len(train_list))
    print(len(test_list))

def get_other_species_from_raw_repbase(raw_repbase, total_repbase, train_species_path, test_species_path):
    raw_names, raw_contigs = read_fasta_v1(raw_repbase)
    names, contigs = read_fasta_v1(total_repbase)
    test_contigs = {}

    name_set = set()
    for name in names:
        parts = name.split('\t')
        seq_name = parts[0]
        name_set.add(seq_name)

    current_name_set = set()
    for raw_name in raw_names:
        parts = raw_name.split('\t')
        seq_name = parts[0]
        if seq_name in name_set:
            current_name_set.add(seq_name)
        else:
            test_contigs[raw_name] = raw_contigs[raw_name]
    store_fasta(contigs, train_species_path)
    store_fasta(test_contigs, test_species_path)
    print(len(test_contigs))
    print(len(name_set))
    print(len(current_name_set))

def extract_non_autonomous(repbase_path, out_path):
    contigNames, contigs = read_fasta_v1(repbase_path)
    # Extract all non-autonomous transposons from Repbase sequences with TSDs, excluding DIRS, Helitron, and Crypton.
    # Define a regular expression to identify all non-autonomous transposons.
    pattern = r'\b\w+[-|_]\d*N\d*[-|_]\w+\b'
    filter_labels = ('Helitron', 'Crypton', 'DIRS', 'Penelope', 'RTE', 'Gypsy', 'L1', 'Retrovirus', 'Jockey', 'I', 'Bel-Pao', 'Copia','R2', 'Unknown')
    non_auto_TEs = {}
    non_auto_repbase_path = out_path
    for name in contigNames:
        parts = name.split('\t')
        seq_name = parts[0]
        label = parts[1]
        matches = re.findall(pattern, seq_name)
        if matches:
            if label not in filter_labels:
                non_auto_TEs[name] = contigs[name]
    print('non-autonomous DNA number: ' + str(len(non_auto_TEs)))
    store_fasta(non_auto_TEs, non_auto_repbase_path)

    contignames, contigs = read_fasta_v1(non_auto_repbase_path)
    type_count = {}
    for name in contignames:
        label = name.split('\t')[1]
        if not type_count.__contains__(label):
            type_count[label] = 0
        count = type_count[label]
        type_count[label] = count + 1
    print('non-autonomous DNA type:')
    print(type_count)
