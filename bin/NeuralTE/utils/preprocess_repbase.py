import argparse
import os
import random
import re
import sys

current_folder = os.path.dirname(os.path.abspath(__file__))
# Add the path to the 'configs' folder to the Python path
configs_folder = os.path.join(current_folder, "..")
sys.path.append(configs_folder)

from configs import config
from utils.data_util import read_fasta_v1, read_fasta, store_fasta, generate_random_sequences, generate_random_sequence


def split_fasta(cur_path, output_dir, num_chunks):
    split_files = []
    os.system('rm -rf ' + output_dir)
    if not os.path.exists(output_dir):
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


def identify_terminals(split_file, output_dir, tool_dir):
    ltrsearch_command = 'cd ' + output_dir + ' && ' + tool_dir + '/ltrsearch ' + split_file
    itrsearch_command = 'cd ' + output_dir + ' && ' + tool_dir + '/itrsearch -i 0.7 -l 7 ' + split_file
    os.system(ltrsearch_command)
    os.system(itrsearch_command)
    ltr_file = split_file + '.ltr'
    tir_file = split_file + '.itr'

    # Read LTR and ITR files to obtain start and end positions of LTR and ITR
    ltr_names, ltr_contigs = read_fasta_v1(ltr_file)
    tir_names, tir_contigs = read_fasta_v1(tir_file)
    LTR_info = {}
    for ltr_name in ltr_names:
        orig_name = ltr_name.split('\t')[0]
        terminal_info = ltr_name.split('\t')[2]
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
        orig_name = tir_name.split('\t')[0]
        terminal_info = tir_name.split('\t')[2]
        TIR_info_parts = terminal_info.split('ITR')[1].split(' ')[0].replace('(', '').replace(')', '').split('..')
        TIR_left_pos_parts = TIR_info_parts[0].split(',')
        TIR_right_pos_parts = TIR_info_parts[1].split(',')
        lTIR_start = int(TIR_left_pos_parts[0])
        lTIR_end = int(TIR_left_pos_parts[1])
        rTIR_start = int(TIR_right_pos_parts[1])
        rTIR_end = int(TIR_right_pos_parts[0])
        TIR_info[orig_name] = (lTIR_start, lTIR_end, rTIR_start, rTIR_end)

    # Update split_file's header, adding two columns LTR:1-206,4552-4757 TIR:1-33,3869-3836
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

    # Get all LTR sequences
    LTR_names = set()
    for name in raw_names:
        # Identify LTR terminal sequences and obtain corresponding internal sequence names
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

    # Store segmented LTR and complete LTR associations
    SegLTR2intactLTR = {}
    new_names = []
    new_contigs = {}
    for name in raw_names:
        if name in LTR_names:
            # If the current sequence is LTR, check for the existence of internal sequences
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
            # If the current sequence is INT, discard it directly, as INTs with LTRs will surely be identified,
            # while INTs without LTRs should be discarded as incomplete LTRs
            pattern = r'\b(\w+(-|_)?)INT((-|_)?\w*)\b'
            matches = re.findall(pattern, name)
            if not matches:
                # Retain other types of transposons
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

    # Store segmented LTR and complete LTR associations
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

def split_dataset(sequences, train_output_file, test_output_file, split_ratio=0.8):
    # Randomly partition into training and test sets
    all_ids = list(sequences.keys())
    random.shuffle(all_ids)
    train_ids = all_ids[:int(split_ratio * len(all_ids))]
    test_ids = all_ids[int(split_ratio * len(all_ids)):]

    # Write to the training set fasta file
    with open(train_output_file, 'w') as f:
        for seq_id in train_ids:
            f.write('>' + seq_id + '\n')
            f.write(sequences[seq_id] + '\n')

    # Write to the test set fasta file
    with open(test_output_file, 'w') as f:
        for seq_id in test_ids:
            f.write('>' + seq_id + '\n')
            f.write(sequences[seq_id] + '\n')


def print_dataset_info(repbase_path):
    repbase_names, repbase_contigs = read_fasta_v1(repbase_path)
    # Count the number of sequences and species
    unique_species = set()
    for name in repbase_names:
        unique_species.add(name.split('\t')[2])
    print('pre-processed Repbase database sequence size: ' + str(len(repbase_names)) + ', total species num: ' + str(
        len(unique_species)))


def main():
    # 1.parse args
    describe_info = '########################## NeuralTE-preprocess_repbase, version ' + str(config.version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--repbase_dir', metavar='repbase_dir', help='Input Repbase directory')
    parser.add_argument('--out_dir', metavar='out_dir', help='Output directory')

    args = parser.parse_args()

    repbase_dir = args.repbase_dir
    out_dir = args.out_dir

    repbase_dir = os.path.realpath(repbase_dir)
    out_dir = os.path.realpath(out_dir)
    config.work_dir = out_dir

    repbase_path = config.work_dir + '/all_repbase.ref'

    # 1. Merge all Repbase files under the Repbase directory, retaining only sequences with headers in the format seq_name\tlabel\tspecies_name
    # 2. Retain sequences that can be converted to Wicker superfamily labels; it's difficult to determine the superfamily category for other sequences
    non_TE = ('tandem repeat', 'Tandem repeat',
              'MSAT', 'SAT', 'Satellite repetitive element', 'satellite', 'Satellite',
              'Simple Repeat', 'Multicopy gene', 'Pseudogene')
    non_TE_count = 0
    files = get_all_files(repbase_dir)
    all_repbase_contigs = {}
    for file in files:
        names, contigs = read_fasta_v1(file)
        for name in names:
            parts = name.split('\t')
            if len(parts) == 3:
                if config.Repbase_wicker_labels.__contains__(parts[1]):
                    all_repbase_contigs[name] = contigs[name]
                elif parts[1] in non_TE:
                    new_name = parts[0] + '\t' + 'Unknown' + '\t' + parts[2]
                    non_TE_count += 1
                    all_repbase_contigs[new_name] = contigs[name]
    # generate random sequences
    num_sequences = 1000
    sequences_lengths = generate_random_sequences(num_sequences)
    for i, l in enumerate(sequences_lengths):
        seq = generate_random_sequence(l)
        all_repbase_contigs['Random_' + str(i + 1) + '\t' + 'Unknown' + '\t' + 'Oryza sativa'] = seq
    print('non_TE_count: ' + str(non_TE_count))
    store_fasta(all_repbase_contigs, repbase_path)
    # 3. Concatenate Repbase sequences' LTRs and Internals, filtering out incomplete LTR sequences
    repbase_path, repbase_labels = connect_LTR(repbase_path)
    print_dataset_info(repbase_path)


if __name__ == '__main__':
    main()